module iter
    include "rank.inc"
    integer:: valid_comm

    integer:: iterate_from, iterate_to, backup_num
    integer, parameter:: BACKUP_PER = 40
    character(LEN=*), parameter :: PREFIX = 'backups/backup-', BACKUP_CONF='backups/1conf'

    logical:: restore_ndone


! constants
    integer, parameter:: XCHG_ROWS=5, CHECK_NONZERO=6, SWAP=7, SKIP=10

! process invariant data
    
    integer:: rows_per_process, comm_size, id

! process backup data
    double precision, dimension (:, :), allocatable :: local_rows
    integer:: rk, row, this_row


contains

subroutine failure_restore()
    use mpi
    integer:: ierr, i

    call MPI_COMM_SIZE(valid_comm, comm_size, ierr)
    rows_per_process = N/(comm_size-1)

    call MPI_COMM_RANK(valid_comm, id, ierr)

    if(id .ne. 0) then
        allocate(local_rows(rows_per_process, M))
    endif

    iterate_to = MIN(N,M)

    write(*, *) "spawned", id

    call restore()
end subroutine

subroutine initialize()
    use mpi
    integer:: ierr, i, stat(MPI_STATUS_SIZE), received_tag
    double precision:: vec(M), vec2(M) 

    valid_comm = MPI_COMM_WORLD

    call MPI_COMM_SIZE(valid_comm, comm_size, ierr)
    rows_per_process = N/(comm_size-1)

    call MPI_COMM_RANK(valid_comm, id, ierr)

    if(id .ne. 0) then
        allocate(local_rows(rows_per_process, M))
        do i = 1, rows_per_process
            call mpi_recv(vec, M, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, valid_comm, stat, ierr)
            received_tag = stat (MPI_TAG)
            local_rows((received_tag - 1)/(comm_size-1) + 1, :) = vec
        enddo

        row = 1
        this_row = 0
    else
        open(1, file="rank.mt")
        do i = 1, N
            read(1,*) vec(:)
            ! write(*,*) "read", i
            call mpi_send(vec, M, MPI_DOUBLE_PRECISION, MOD(i-1, comm_size-1)+1, i, valid_comm, ierr)
        enddo
        close(1)

        rk = MIN(N,M)
    endif

    iterate_from = 1
    iterate_to = MIN(N,M)

    backup_num = 0

    restore_ndone = .TRUE.
end subroutine

subroutine backup_file_name(str, id)
    integer:: id
    character(len=*), parameter:: fmt_procnum = '(I3.3)', fmt_backnum = '(I10.10)'
    character(len=100):: procnum, backnum
    character(len=*):: str

    write(procnum,fmt_procnum) id
    write(backnum,fmt_backnum) backup_num

    str = PREFIX//trim(backnum)//'-'//trim(procnum)//".bak"

end subroutine

subroutine exit_iteration()
    if(id .eq. 0) then 
        write(*,*) "rank is", rk
        write(*,*) "(!!!)", comm_size, time, N, M
    else
    endif
end subroutine

subroutine make_backup(iter, err)
    integer:: ierr, iter, err
    character(len=100):: filename

    backup_num = backup_num + 1

    call backup_file_name(filename, id)

    if(id .eq. 0) then
        open (unit=1, file=filename, form='unformatted', access='direct', recl=sizeof(iter))
        
        write(1, rec=1) iter
        write(1, rec=2) rk

        close(1)
    else
        open (unit=1, file=trim(filename)//".scal", form='unformatted', access='direct', recl=sizeof(iter))
        
        write(1, rec=1) iter
        write(1, rec=2) row
        write(1, rec=3) this_row
        
        close(1)

        open (unit=1, file=trim(filename)//".mat", form='unformatted', access='direct', &
                recl=rows_per_process * M * sizeof(local_rows(1,1)))
        write(1, rec=1) local_rows
        close(1)
    endif

    call MPI_BARRIER(valid_comm ,ierr)
    if(ierr .ne. 0) then
        err = -1
        return
    endif

    if(id .eq. 0) then
        open (unit=1, file=BACKUP_CONF, form='unformatted', access='direct', recl=sizeof(backup_num))
        write(1, rec=1) backup_num
        close(1)

        write(*, *) "make backup", id, backup_num
    endif

    call MPI_BARRIER(valid_comm ,ierr)
    if(ierr .ne. 0) then
        ! write(*, *) "backup failure", id, backup_num
        err = -1
        return
    endif

    ! write(*, *) "backup made successfully", id, backup_num


    err = 0

    return
end subroutine

subroutine restore()
    character(len=100) :: filename

    restore_ndone = .FALSE.


    call MPI_BARRIER(valid_comm ,ierr)

    call MPI_COMM_SIZE(valid_comm, comm_size, ierr)
    rows_per_process = N/(comm_size-1)

    call MPI_COMM_RANK(valid_comm, id, ierr)

    open (unit=1, file=BACKUP_CONF, form='unformatted', access='direct', recl=sizeof(backup_num))
    read(1, rec=1) backup_num
    close(1)

    call MPI_COMM_RANK(valid_comm, id, ierr)
    call backup_file_name(filename, id)

    if(id .eq. 0) then
        open (unit=1, file=filename, form='unformatted', access='direct', recl=sizeof(iterate_from))
        
        read(1, rec=1) iterate_from
        read(1, rec=2) rk

        close(1)
    else
        open (unit=1, file=trim(filename)//".scal", form='unformatted', access='direct', recl=sizeof(iterate_from))
        
        read(1, rec=1) iterate_from
        read(1, rec=2) row
        read(1, rec=3) this_row
        
        close(1)

        open (unit=1, file=trim(filename)//".mat", form='unformatted', access='direct', &
                recl=rows_per_process * M * sizeof(local_rows(1,1)))
        read(1, rec=1) local_rows
        close(1)
    endif

    call MPI_COMM_RANK(valid_comm, id, ierr)
    write(*, *) "restore", id, backup_num, iterate_from
    call MPI_BARRIER(valid_comm ,ierr)
    return
end subroutine

subroutine iterate(err)
    use mpi

    integer:: err, ierr

    integer:: proc_with_nonzero, i, nonzero_idx, stat(MPI_STATUS_SIZE), col, tag
    logical:: nonzero, this_row_here

    parameter (EPS=1.e-3)
    double precision:: vec(M), vec2(M) 

    ! do i = iterate_from, iterate_to
    !     if(MOD(i, BACKUP_PER) .eq. 0) then
    !         call make_backup(i)
    !     endif        
    ! enddo

    if (id .eq. 0) then
        do col=iterate_from,iterate_to
            ! write(*, *) id, col
            ! write(*, *) col
            if(MOD(col, BACKUP_PER) .eq. 0) then
                call make_backup(col, ierr)
                if(ierr .ne. 0) then
                    err = -1
                    return
                endif
            endif   
            
            

            !choose nonzero
            proc_with_nonzero = -1
            ! write(*,*) "start getting nonzero..."
            do i=1,comm_size-1
                ! write(*,*) "getting nonzero...", i
                call mpi_recv(nonzero, 1, MPI_LOGICAL, MPI_ANY_SOURCE, CHECK_NONZERO, valid_comm, stat, ierr)
                ! write(*,*) "got nonzero"
                if(ierr .ne. 0) then
                    err = -1
                    return
                endif
                if((proc_with_nonzero .eq. -1) .AND. nonzero ) then
                    proc_with_nonzero = stat (MPI_SOURCE)
                endif
            enddo
            ! write(*,*) "fin getting nonzero"
            if(proc_with_nonzero .eq. -1) then
                rk = rk - 1
                ! write(*,*) "start send broadcast..."
                do i=1,comm_size-1
                    call mpi_send(0, 1, MPI_DOUBLE_PRECISION, i, SKIP, valid_comm, ierr)
                    if(ierr .ne. 0) then
                        err = -1
                        return
                    endif
                enddo
                ! write(*,*) "end send broadcast"
            else
                ! write(*,*) "start send swap..."
                call mpi_send(0, 1, MPI_DOUBLE_PRECISION, proc_with_nonzero, SWAP, valid_comm, ierr)
                ! write(*,*) "end send swap"
            endif
        enddo
    else
    !code for non-root
        do col=iterate_from, iterate_to
            ! write(*, *) id, col
            ! if(id .eq. 2) then
            !     write(*,*) local_rows(rows_per_process, :)
            ! endif

            if(MOD(col, BACKUP_PER) .eq. 0) then
                call make_backup(col, ierr)
                if(ierr .ne. 0) then
                    err = -1
                    return
                endif
            endif  

            nonzero = .false.
            do i = this_row+1, rows_per_process
                if( ABS(local_rows(i, col)) .GE. EPS) then
                    nonzero = .true.
                    nonzero_idx = i
                    exit
                endif
            enddo
            
            call mpi_send(nonzero, 1, MPI_LOGICAL, 0, CHECK_NONZERO, valid_comm, ierr)
            if(ierr .ne. 0) then
                err = -1
                return
            endif
            
            this_row_here = .false.
            this_row = (row - 1)/(comm_size - 1) + 1
            if( (this_row - 1)*(comm_size - 1) + id .gt. row ) then
                this_row = this_row - 1
            elseif((this_row - 1)*(comm_size - 1) + id .eq. row) then
                this_row_here = .true.
            endif
                
            call mpi_recv(vec, M, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, valid_comm, stat, ierr)
            if(ierr .ne. 0) then
                err = -1
                return
            endif

            if (stat (MPI_TAG) .eq. SWAP) then
                vec = local_rows(nonzero_idx, :)
                do i=1,comm_size-1
                    if(i .ne. id) then
                        call mpi_send(vec, M, MPI_DOUBLE_PRECISION, i, XCHG_ROWS, valid_comm, ierr)
                        if(ierr .ne. 0) then
                            err = -1
                            return
                        endif
                    endif
                enddo

                if (this_row_here .neqv. .true.) then
                    call mpi_recv(vec2, M, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, valid_comm, stat, ierr)
                    if(ierr .ne. 0) then
                        err = -1
                        return
                    endif

                else
                    vec2 = local_rows(this_row, :)
                endif

                local_rows(nonzero_idx, :) = vec2
            endif

            if(stat(MPI_TAG) .ne. SKIP .and. stat(MPI_TAG) .ne. SWAP) then
                if (this_row_here) then
                    vec2 = local_rows(this_row, :)
                    call mpi_send(vec2, M, MPI_DOUBLE_PRECISION, stat(MPI_SOURCE), XCHG_ROWS, valid_comm, ierr)
                    if(ierr .ne. 0) then
                        err = -1
                        return
                    endif
                endif
            endif

            if(stat (MPI_TAG) .ne. SKIP) then
                vec(:) = vec(:) / vec(col)
                !$OMP PARALLEL DO PRIVATE(j) SHARED(local_rows) FIRSTPRIVATE(col, this_row) 
                do j = this_row + 1, rows_per_process
                    local_rows(j, :) = local_rows(j, :) - vec * local_rows(j, col)

                    if(col .eq. 99 .and. id .eq. 3 .and. restore_ndone) then
                        write(*, *) id, "is going to die"
                        call abort
                    endif
                enddo
                !$OMP END PARALLEL DO
                row = row + 1
            endif
        enddo
    endif

    err = 0
    return
end subroutine

subroutine run_iterations()
    integer:: err

    err = 1
    do while(err .ne. 0)
        err = 0

        call iterate(err)

        if(err .ne. 0) then
            call restore()
        endif
    enddo

    call exit_iteration()

    return
end subroutine


end module