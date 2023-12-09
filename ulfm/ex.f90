subroutine err_handler(comm, err_code)
    use mpi
    use mpi_ext
    use common_funcs
    use iter
    implicit none

    integer:: i,id_old, ierr, comm, err_code, fail_n, small_intracomm, intercomm, final_intracomm, id_new
    integer:: prev_size, curr_size, flag, errlen
    integer:: shrinked_comm, fail_group, comm_group, local_group, remote_group
    
    id_old = id
    write(*, *) id, "errhandler"

    call MPIX_COMM_REVOKE(comm, ierr)

    call MPIX_COMM_SHRINK(comm, small_intracomm, ierr)
    call MPI_BARRIER(small_intracomm, ierr)

    if (ierr .eq. 0) then
        call MPI_COMM_SIZE(comm, prev_size, ierr)
        call MPI_COMM_SIZE(small_intracomm, curr_size, ierr)
        fail_n = prev_size - curr_size
        
        call MPI_COMM_RANK(small_intracomm, id_new, ierr)

        call MPI_COMM_SPAWN("./slave", MPI_ARGV_NULL, fail_n, MPI_INFO_NULL, 0, &
                                        small_intracomm, intercomm, MPI_ERRCODES_IGNORE, ierr)

        call MPI_INTERCOMM_MERGE(intercomm, .TRUE., final_intracomm, ierr)
        call MPI_COMM_SIZE(final_intracomm, curr_size, ierr)
        call MPI_COMM_RANK(final_intracomm, id_new, ierr)

        write(*,*) id_old, "->", id_new

        valid_comm = final_intracomm
    endif
end subroutine


program ex
    use mpi
    use mpi_ext
    use common_funcs
    use iter
    implicit none
    external err_handler
    integer:: i, j, ierr, stat(MPI_STATUS_SIZE), errhandler_id
    integer:: flag
    call MPI_INIT ( ierr )
    valid_comm = MPI_COMM_WORLD

    call MPI_COMM_CREATE_ERRHANDLER(err_handler, errhandler_id, ierr)
    call MPI_COMM_SET_ERRHANDLER(valid_comm, errhandler_id, ierr)

    call initialize()
    
    call run_iterations()

    call MPI_FINALIZE(ierr)
    return

    ! call print_rank()
    if (id .le. 5) then
        write(*, *) id, "is going to die"
        call abort
    endif
    write(*, *) "if passed by ", id

    call MPI_BARRIER(valid_comm, ierr)

    call print_rank()

    call MPI_FINALIZE(ierr)
end program