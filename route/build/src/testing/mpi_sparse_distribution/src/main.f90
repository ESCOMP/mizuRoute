program sparse_distribution

USE nrtype
USE globalData, only: mpicom_route, pid, nNodes
USE mpi_utils

implicit none

integer(i4b)             :: i,r
integer(i4b)             :: num_elems
integer(i4b),allocatable :: dest_index(:)
integer(i4b),allocatable :: dest_rank(:)
real(dp),allocatable     :: array_src(:)
real(dp),allocatable     :: array_dest(:)
integer(i4b)             :: ierr
character(strLen)        :: cmessage
character(strLen)        :: agg_method='sum'

call shr_mpi_init(mpicom_route, cmessage)
call shr_mpi_commsize(mpicom_route, nNodes, cmessage)
call shr_mpi_commrank(mpicom_route, pid, cmessage)

num_elems = 10+pid

allocate(array_src(num_elems), array_dest(num_elems))
allocate(dest_index(num_elems), dest_rank(num_elems))

!-----------------------------------------
! Fill test data
!-----------------------------------------
do i = 1, num_elems
  array_src(i) = 1.0_dp * (pid + 1)      ! easy to track
  dest_rank(i)  = mod(i + pid, nNodes)  ! send to various ranks
end do

do i = 1, num_elems
  dest_index(i) = mod(i, 10) + 1          ! indices 1..n_rank
end do

do r = 0, nNodes-1
  if (pid == r) then
    write(*,'(A,I4,1X,A,*(1X,F0.2))') "pid:", pid, "array_src :", (array_src(i), i=1,size(array_src))
    write(*,'(A,I4,1X,A,*(1X,I4))')   "pid:", pid, "dest_rank :", (dest_rank(i), i=1,size(array_src))
    write(*,'(A,I4,1X,A,*(1X,I4))')   "pid:", pid, "dest_index:", (dest_index(i),i=1,size(array_src))
  end if
  call shr_mpi_barrier(mpicom_route, cmessage)
end do

call shr_mpi_sparse_distribute(array_src, dest_rank, dest_index, &
                               array_dest, agg=agg_method, fillvalue=0._dp)

!-----------------------------------------
! Print results
!-----------------------------------------
if (pid==0) then
  write(*, '(2A)') '- Result- aggregation method: ', trim(agg_method)
endif
do r = 0, nNodes-1
  if (pid == r) then
    write(*,'(A,I4,1X,A,*(1X,F0.1))') "pid:", pid, "array_dest :", (array_dest(i), i=1,size(array_dest))
  end if
  call shr_mpi_barrier(mpicom_route, cmessage)
end do

!  Shut down MPI
call MPI_FINALIZE(ierr)

stop

end program sparse_distribution

