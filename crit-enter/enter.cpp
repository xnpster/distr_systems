#include <iostream>
#include <vector>
#include <set>
#include "mpi.h"

#include <chrono>
#include <thread>

#include <fstream>

using namespace std;

enum {
    ANS_OK,
    ANS_FIN,

    ENTER_REQUEST,
    ENTER_ANSWER,

    FIN_CRIT
};

static double base_time;
double get_time() {
    return MPI_Wtime() - base_time;
}

static int glob_rank, glob_size;
static set<int> not_answered;

void debug(const char* s) {
    if(false) {
        cout << "[" << glob_rank << "] " << s << endl;
    }
}

void hold(int t) {
    /*for(int k = 0; k < t; k++) {
        int a = 0;
        for(int j = 0; j < 1000000; j++)
            a++;
    }*/

    std::this_thread::sleep_for(std::chrono::milliseconds(t));
}

static vector<MPI_Request> enter_req;
static vector<MPI_Request> ans_req;

MPI_Request* add_request(vector<MPI_Request>& reqs) {
    reqs.push_back(MPI_REQUEST_NULL);
    return &(reqs.back());
}

void free_requests(vector<MPI_Request>& reqs) {
    for(auto& e : reqs)
        MPI_Wait(&e, MPI_STATUS_IGNORE);

    reqs.clear();
}


void critical_work() {
    /*hold(80);
    cout << "================" << endl;
    for(int i = 0; i < 2; i++) {
        cout << "[it " << i << "] " << glob_rank << " \\ " << glob_size << endl;
        hold(400);
    }
    cout << "================" << endl;
    hold(80);*/

    const char* name = "critical.txt";

    hold(1000);
    if(ifstream(name)) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else {
        ofstream(name) << "block" << endl;
        hold(1000);
        remove(name);
    }
}

static vector<MPI_Request> fin_req;
static int ans_fin = ANS_FIN;
void notify_fin() {
    int ierr;

    for(int dest = 0; dest < glob_size; dest++) {
        if(dest != glob_rank) {
            MPI_Isend(&ans_fin, 1, MPI_INT, dest, FIN_CRIT, MPI_COMM_WORLD, add_request(fin_req));
        }
    }
}

bool test_fin() {
    debug("test_fin enter");
    static int finished = 1;

    int ierr, ans, flag = true;
    MPI_Status stat;

    while(flag) {
        MPI_Iprobe(MPI_ANY_SOURCE, FIN_CRIT, MPI_COMM_WORLD, &flag, &stat);

        if(flag) {
            int src = stat.MPI_SOURCE;
            debug("test_fin recv ans");
            ierr = MPI_Recv(&ans, 1, MPI_INT, src, FIN_CRIT, MPI_COMM_WORLD, &stat);
            debug("test_fin recv ans done");
            finished++;
        }
    }

    if(finished == glob_size) {
        free_requests(fin_req);
        debug("test_fin leave true");
        return true;
    }
    
    debug("test_fin leave false");
    return false;
}

// invoke this function periodically to let processes enter critical section
void check_entering() {
    int ierr, ok_ans = ANS_OK;
    MPI_Status stat;
    int flag = true;

    while(flag) {
        ierr = MPI_Iprobe(MPI_ANY_SOURCE, ENTER_REQUEST, MPI_COMM_WORLD, &flag, &stat);

        if(flag) {
            int src = stat.MPI_SOURCE;
            double enter_req_timestamp;

            ierr = MPI_Recv(&enter_req_timestamp, 1, MPI_DOUBLE, src, ENTER_REQUEST, MPI_COMM_WORLD, &stat);
            MPI_Isend(&ok_ans, 1, MPI_INT, src, ENTER_ANSWER, MPI_COMM_WORLD, add_request(ans_req));
        }
    }

    while(flag) {
        ierr = MPI_Iprobe(MPI_ANY_SOURCE, FIN_CRIT, MPI_COMM_WORLD, &flag, &stat);

        if(flag) {
            int src = stat.MPI_SOURCE;
            double enter_req_timestamp;

            ierr = MPI_Recv(&enter_req_timestamp, 1, MPI_DOUBLE, src, ENTER_REQUEST, MPI_COMM_WORLD, &stat);
            MPI_Isend(&ok_ans, 1, MPI_INT, src, ENTER_ANSWER, MPI_COMM_WORLD, add_request(ans_req));
        }
    }

    free_requests(ans_req);
}


void do_critical(void (*run)()) {
    int ierr;
    double current_time = get_time();
    int ok_ans = ANS_OK;

    // send request for entering critical section
    for(int dest = 0; dest < glob_size; dest++) {
        if(dest != glob_rank) {
            ierr = MPI_Isend(&current_time, 1, MPI_DOUBLE, dest, ENTER_REQUEST, 
                MPI_COMM_WORLD, add_request(enter_req));
        }
    }

    // now we should recieve positive answer from all processes
    vector<bool> enter_req_ok(glob_size, false);
    int ok_cnt = 0;

    while(ok_cnt != glob_size - 1) {
        MPI_Status stat;
        int flag = true;

        // look for answers
        while(flag) {
            ierr = MPI_Iprobe(MPI_ANY_SOURCE, ENTER_ANSWER, MPI_COMM_WORLD, &flag, &stat);

            if(flag) {
                int src = stat.MPI_SOURCE, ans;

                ierr = MPI_Recv(&ans, 1, MPI_INT,src, ENTER_ANSWER, MPI_COMM_WORLD, &stat);
                        
                if(src != glob_rank && !enter_req_ok[src]) {
                    if(ans != ANS_OK)
                        MPI_Abort(MPI_COMM_WORLD, 1);

                    enter_req_ok[src] = true;
                    ok_cnt++;
                }
            }
        }

        // answer to entering processes
        flag = true;
        while(flag) {
            ierr = MPI_Iprobe(MPI_ANY_SOURCE, ENTER_REQUEST, MPI_COMM_WORLD, &flag, &stat);

            if(flag) {
                int src = stat.MPI_SOURCE;
                double enter_req_timestamp;

                ierr = MPI_Recv(&enter_req_timestamp, 1, MPI_DOUBLE, src, ENTER_REQUEST, MPI_COMM_WORLD, &stat);
                
                // decide whose turn it is to enter
                if(enter_req_timestamp < current_time) {
                    // src proc
                    MPI_Isend(&ok_ans, 1, MPI_INT, src, ENTER_ANSWER, MPI_COMM_WORLD, add_request(ans_req));
                } else {
                    // this proc
                    not_answered.insert(src);
                }
            }
        }

        /* TODO: do it except FIN_CRIT flag
        // wait for message if we cant leave loop
        if(ok_cnt != glob_size - 1)
            ierr = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        */
        hold(150);
    }

    // now process in critical section
    run();
    
    // leave critcal section

    // resume waiting processes
    for(int dest : not_answered) {
        MPI_Isend(&ok_ans, 1, MPI_INT, dest, ENTER_ANSWER, MPI_COMM_WORLD, add_request(ans_req));
    }

    free_requests(ans_req);
    free_requests(enter_req);

    check_entering();

    notify_fin();
    while(!test_fin()) {
        check_entering();
        hold(150);
    }


}


int main(int argc, char** argv) {
    int err = 0;

    err = MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    base_time = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);


    err = MPI_Comm_size(MPI_COMM_WORLD, &glob_size);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &glob_rank);

    do_critical(critical_work);

    cout << "Done " << glob_rank << " \\ " << glob_size << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    err = MPI_Finalize(); 

    return err;
}