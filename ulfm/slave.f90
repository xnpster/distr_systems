program slave
    use common_funcs
    use mpi
    use iter
    implicit none

    integer:: ierr, intercomm, group, final_intracomm
    call MPI_INIT ( ierr )

    write(*,*) "im child!"

    call MPI_COMM_GET_PARENT(intercomm, ierr);

    call MPI_INTERCOMM_MERGE(intercomm, .FALSE., final_intracomm, ierr)

    valid_comm = final_intracomm

    call failure_restore()
    call run_iterations()
    
    call MPI_FINALIZE ( ierr )
end program