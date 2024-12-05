MODULE advection_diffusion

! solving advection-diffusion equation

USE nrtype
USE dataTypes,     ONLY: RCHPRP          ! Reach parameter
USE public_var,    ONLY: iulog           ! i/o logical unit number
USE public_var,    ONLY: realMissing     ! missing value for real number
USE public_var,    ONLY: integerMissing  ! missing value for integer number

implicit none

private
public::solve_ade    ! solve advection-diffusion equation

CONTAINS

 ! *********************************************************************
 ! subroutine: solve advection-diffuision equation (ade)
 ! *********************************************************************
 SUBROUTINE solve_ade(reach_length,  & ! input: Reach length [m]
                      nMolecule,     & ! input: number of sub-segments
                      dt_local,      & ! input: time_step [sec]
                      Qupstream,     & ! input: quantity from upstream [unit of quantity]
                      ck,            & ! input: velocity [m/s]
                      dk,            & ! input: diffusivity [m2/s]
                      Qlat,          & ! input: lateral quantity into chaneel [unit of quantity]
                      Qprev,         & ! input: quantity at previous time step [unit of quantity]
                      Qlocal,        & ! inout: quantity soloved at current time step [unit of quantity]
                      verbose,       & ! input: reach index to be examined
                      ierr,message)
 ! ----------------------------------------------------------------------------------------
 ! Solve linearlized advection diffusive equation per reach and time step.
 !  dQ/dt + ck*dQ/dx = dk*d2Q/dx2 + Qlat - a)
 !
 !  ck (velocity) and dk (diffusivity) are given by input argument
 !
 !  1) dQ/dt   = (Q(t,x) - Q(t-1,x))/dt
 !  2) dQ/dx   = [(1-wck)(Q(t-1,x+1)-Q(t-1,x-1)) + wck*(Q(t,x+1)-Q(t,x-1))]/2dx  (using central difference)
 !  3) d2Q/dx2 = [(1-wdk)(Q(t-1,x+1)-2Q(t-1,x)+Q(t-1,x-1)) + wdk*(Q(t,x+1)-2Q(t,x)+Q(t,x-1))]/2dx
 !
 !  upstream B.C:   Dirchlet BC with inflow at current time-step,t, from upstream basin
 !  downstream B.C: Neumann BC with prescribed quantity gradient (Sbc)
 !                  dQ/dx|x=N = Sbc ->  4) Q(t,N)-Q(t,N-1)) = Sbc*dx
 !  Another downstream B.C option is absorbing boundary condition
 !                  dQ/dt|x=N + ck*dQ/dx|x=N = 0
 !
 !  Inserting 1), 2), 3) and 4) into a) and moving Q(t,*) terms to left-hand side and Q(t-1,*) terms to the righ-hand side
 !  results in tridiagonal matrix equation A*Q = b
 !  where A is [N x N] matrix, Q is [N x 1] vector to be solved (next time step Q) and b is [N x 1] vector
 !  N (nMolecule is used in the code) is the number of internal nodes including upstream and downstream boundaries
 !
 !  Since A is a tridiagonal matrix, the code stores only three diagnoal elements - upper, diagonal, and lower
 !  solving the matrix equation use thomas algorithm
 !
 !  if dk = 0, this equation becomes advection equation and solved with central difference for advection term
 ! ----------------------------------------------------------------------------------------
 implicit none
 ! Argument variables
 real(dp),     intent(in)      :: reach_length              ! River reach length [m]
 integer(i4b), intent(in)      :: nMolecule                 ! number of sub-segments
 real(dp),     intent(in)      :: dt_local                  ! time inteval for time-step [sec]
 real(dp),     intent(in)      :: Qupstream                 ! quantity at top of the reach being processed
 real(dp),     intent(in)      :: ck                        ! velocity [m/s]
 real(dp),     intent(in)      :: dk                        ! diffusivity [m2/s]
 real(dp),     intent(in)      :: Qlat                      ! lateral quantity into chaneel [m3/s]
 real(dp),     intent(in)      :: Qprev(nMolecule)          ! sub-reach quantity at previous time step [m3/s]
 real(dp),     intent(out)     :: Qlocal(nMolecule,0:1)     ! sub-reach & sub-time step quantity at previous and current time step [m3/s]
 logical(lgt), intent(in)      :: verbose                   ! reach index to be examined
 integer(i4b), intent(out)     :: ierr                      ! error code
 character(*), intent(out)     :: message                   ! error message
 ! Local variables
 real(dp)                      :: Cd                        ! Fourier number
 real(dp)                      :: Ca                        ! Courant number
 real(dp)                      :: dx                        ! length of segment [m]
 real(dp)                      :: Sbc                       ! neumann BC slope
 real(dp)                      :: diagonal(nMolecule,3)     ! diagonal part of matrix - diagonal(:,1)=upper, diagonal(:,2)=middle, diagonal(:,3)=lower
 real(dp)                      :: b(nMolecule)              ! right-hand side of the matrix equation
 real(dp)                      :: Qsolved(nMolecule)        ! solved Q at sub-reach at current time step [m3/s]
 real(dp)                      :: wck                       ! weight for advection
 real(dp)                      :: wdk                       ! weight for diffusion
 integer(i4b)                  :: ix                        ! loop index
 integer(i4b)                  :: Nx                        ! number of internal reach segments
 integer(i4b)                  :: downstreamBC              ! method of B.C condition - absorbing or Neumann
 character(len=strLen)         :: fmt1                      ! format string
 character(len=strLen)         :: cmessage                  ! error message from subroutine
 ! Local parameters
 integer(i4b), parameter                :: absorbingBC=1
 integer(i4b), parameter                :: neumannBC=2      ! flux derivative w.r.t. distance at downstream boundary

 ierr=0; message='solve_ade/'

 ! hard-coded parameters
 downstreamBC = neumannBC  ! downstream boundary condition
 wck = 1.0                 ! weight in advection term
 wdk = 1.0                 ! weight in diffusion term 0.0-> fully explicit, 0.5-> Crank-Nicolson, 1.0 -> fully implicit

 Nx = nMolecule - 1  ! Nx: number of internal reach segments

 ! Get the reach parameters
  dx = reach_length/(Nx-1) ! one extra sub-segment beyond outlet

  if (verbose) then
    write(iulog,'(A,1X,G12.5)') ' length [m]     =',reach_length
    write(iulog,'(A,1X,G12.5)') ' time-step [sec]=',dt_local
  end if

 Qlocal(1:nMolecule, 0) = Qprev ! previous time step
 Qlocal(1:nMolecule, 1) = realMissing ! initialize current time step part
 Qlocal(1,1)  = Qupstream     ! quantity from upstream at current time step

 ! Fourier number and Courant number
 Cd = dk*dt_local/dx**2
 Ca = ck*dt_local/dx

 ! create a matrix - current time step
 ! populate tridiagonal elements
 ! diagonal
 diagonal(1,2)             = 1._dp
 diagonal(2:nMolecule-1,2) = 2._dp + 4*wdk*Cd
 if (downstreamBC == absorbingBC) then
   diagonal(nMolecule,2)     = 1._dp + wck*Ca
 else if (downstreamBC == neumannBC) then
   diagonal(nMolecule,2)     = 1._dp
 end if

 ! upper
 diagonal(:,1)           = 0._dp
 diagonal(3:nMolecule,1) = wck*Ca - 2._dp*wdk*Cd

 ! lower
 diagonal(:,3)             = 0._dp
 diagonal(1:nMolecule-2,3) = -wck*Ca - 2._dp*wdk*Cd
 if (downstreamBC == absorbingBC) then
   diagonal(nMolecule-1,3)   = -wck*Ca
 else if (downstreamBC == neumannBC) then
   diagonal(nMolecule-1,3)     = -1._dp
 end if

 ! populate right-hand side
 ! upstream boundary condition
 b(1)             = Qlocal(1,1)
 ! downstream boundary condition
 if (downstreamBC == absorbingBC) then
   b(nMolecule) = (1._dp-(1._dp-wck)*Ca)*Qlocal(nMolecule,0) + (1-wck)*Ca*Qlocal(nMolecule-1,0)
 else if (downstreamBC == neumannBC) then
   Sbc = (Qlocal(nMolecule,0)-Qlocal(nMolecule-1,0))
   b(nMolecule)     = Sbc
 end if
 ! internal node points
 b(2:nMolecule-1) = ((1._dp-wck)*Ca +2._dp*(1._dp-wdk)*Cd)*Qlocal(1:nMolecule-2,0)  &
                  + (2._dp-4._dp*(1._dp-wdk)*Cd)*Qlocal(2:nMolecule-1,0)           &
                  - ((1._dp-wck)*Ca -2._dp*(1._dp-wdk)*Cd)*Qlocal(3:nMolecule,0)

 ! solve matrix equation - get updated Qlocal
 call TDMA(nMolecule, diagonal, b, Qsolved)

 if (verbose) then
   write(fmt1,'(A,I5,A)') '(A,1X',nMolecule,'(1X,G15.4))'
   write(iulog,fmt1) ' Q sub_reqch=', (Qsolved(ix), ix=1,nMolecule)
 end if

 Qlocal(:,1) = Qsolved

 if (verbose) then
   write(fmt1,'(A,I5,A)') '(A,1X',nMolecule,'(1X,G15.4))'
   write(iulog,fmt1) ' Qprev(1:nMolecule)= ', Qprev(1:nMolecule)
   write(iulog,'(A,3(1X,G12.5))') ' ck, dk= ',ck, dk
   write(iulog,'(A,2(1X,G12.5))') ' Cd, Ca= ', Cd, Ca
   write(iulog,fmt1) ' diagonal(:,1)= ', diagonal(:,1)
   write(iulog,fmt1) ' diagonal(:,2)= ', diagonal(:,2)
   write(iulog,fmt1) ' diagonal(:,3)= ', diagonal(:,3)
   write(iulog,fmt1) ' b= ', b(1:nMolecule)
 end if

 END SUBROUTINE solve_ade

 SUBROUTINE TDMA(NX,MAT,b,T)
 ! Solve tridiagonal matrix system of linear equation
 ! NX is the number of unknowns (gridblocks minus boundaries)
 ! Solve system of linear equations, A*T = b where A is tridiagonal matrix
 ! MAT = NX x 3 array holding tri-diagonal portion of A
 ! MAT(NX,1) - uppder diagonals for matrix A
 ! MAT(NX,2) - diagonals for matrix A
 ! MAT(NX,3) - lower diagonals for matrix A
 ! b(NX) - vector of the right hand coefficients b
 ! T(NX) - The solution matrix
 !
 ! example, A
 ! d u 0 0 0
 ! l d u 0 0
 ! 0 l d u 0
 ! 0 0 l d u
 ! 0 0 0 l d
 !
 ! MAT(:,1) = [0, u, u, u, u]
 ! MAT(:,2) = [d, d, d, d, d]
 ! MAT(:,3) = [l, l, l, l, 0]

   implicit none
   ! Argument variables
   integer(i4b),  intent(in)     :: NX     ! number of unknown (= number of matrix size, grid point minus two end points)
   real(dp),      intent(in)     :: MAT(NX,3)
   real(dp),      intent(in)     :: b(NX)
   real(dp),      intent(inout)  :: T(NX)
   ! Local variables
   integer(i4b)                  :: ix
   real(dp)                      :: U(NX)
   real(dp)                      :: D(NX)
   real(dp)                      :: L(NX)
   real(dp)                      :: b1(NX)
   real(dp)                      :: coef

   U(1:NX) = MAT(1:NX,1)
   D(1:NX) = MAT(1:NX,2)
   L(1:NX) = MAT(1:NX,3)
   b1(1:NX) = b(1:NX)
   do ix = 2, NX
     coef  = L(ix-1)/D(ix-1)
     D(ix) = D(ix)-coef*U(ix)
     b1(ix) = b1(ix)-coef*b1(ix-1)
   end do

   T(NX) = b1(NX)/D(NX) ! Starts the countdown of answers
   do ix = NX-1, 1, -1
       T(ix) = (b1(ix) - U(ix+1)*T(ix+1))/D(ix)
   end do

 END SUBROUTINE TDMA

END MODULE advection_diffusion
