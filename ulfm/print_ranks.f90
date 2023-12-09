module common_funcs
contains
    
subroutine print_rank()
    use iter
    implicit none
    include 'mpif.h'

    integer:: ierr, stat(MPI_STATUS_SIZE)

    call mpi_comm_size(valid_comm, comm_size, ierr)
    call mpi_comm_rank(valid_comm, id, ierr)

    write(*, *) id, "/", comm_size
end subroutine


subroutine print_intercomm(intercomm)
    use iter
    implicit none
    include 'mpif.h'

    integer:: ierr, stat(MPI_STATUS_SIZE), intercomm, local_size, remote_size

    call mpi_comm_size(intercomm, local_size, ierr)
    call mpi_comm_remote_size(intercomm, remote_size, ierr)

    write(*, *) "intercomm", local_size, "-->", remote_size
end subroutine

end module common_funcs