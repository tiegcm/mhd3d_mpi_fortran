!----------------------------------------------------
!  mhd3d2.f90
!  solve ideal 3d mhd equations with finite volume 
!  method in arbitrary geometry
!
!  updated 06/15/16 by B Zhu (ben.zhu@dartmouth.edu)
!  1. four spatial reconstruction options: PDM, 
!     TVD(van Leer), WENO and MPWENO; three temporal
!     stepping options: Adams-Bashforth, TVD RK2 and
!     TVD RK3
!
!  updated 08/18/15 by B Zhu (ben.zhu@dartmouth.edu)
!  1. dynamic allocation for furture impliment of GPU
!  2. 2d mpi domain decomposition along x and y directions,
!     best for 2d problem, eg, 2d reconnection
!  3. use hdf5 output (.out output need debug)
!  4. need to add restart subroutine
!
!  CAUTION: index is tricky in this code
!
!  variable defined at the 
!    center of cells: rho,vx,vy,vz,p,bx,by,bz
!             (and thus rhovx,rhovy,rhovz,...)
!    face center of cells: bi,bj,bk 
!    face corner of cells: ei,ej,ek
!
!  for momery concern and consistence with future 2d/3d
!  domain decomposed mpi version of this code, all flux
!  and stress on the rhs of equations are defined w/o 
!  any ghost cells. 
!
!-----------------------------------------------------
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module pars
   implicit none
!..Include the mpi header file
   include 'mpif.h'              ! Required statement

!..Fundamental constants (these never never change):
   real*8, parameter :: pi=4.*atan(1.0_8)
   real*8, parameter :: twopi=2.*pi
   complex*16, parameter :: ci=(0.0,1.0)
!..mhd normalization constants
   real*8, parameter :: mu0=4.*pi*1.e-7
   real*8, parameter :: mp=1.67e-27 ![kg]
   real*8, parameter :: re=6.38e6 ![m]
   real*8, parameter :: u0=100.e3 ![m/s]
   real*8, parameter :: t0=re/u0 ![s]
   real*8, parameter :: rhon=mp*1.e6 ![kg/m^3] ! rho0 is uesd later
   real*8, parameter :: p0=rhon*u0**2 ![N/m^2]
   real*8, parameter :: b0=sqrt(mu0*p0)*1.e9 ![nT]
   real*8, parameter :: gam=5./3.

!..Program name. Length = exact length of characters used
   character(6), parameter :: progname = 'mhd3d2'

!..Directory where results will be stored. Output module stores
!..in resdir/progname if progname exists or creates it if it doesn't
!..IMPORTANT: the length of the following strings must equal the number of characters
!   character(41) :: outdir = '/global/scratch/bzhu/mhd/otv/pdm_ab_1024/'
   character(35) :: outdir = '/global/scratch/bzhu/test/otv/0128/'
!   character(2) :: outdir = './'
!..If starting from a restart file (restart.in), indicate location of restart file
!..BZ:currently no restart subroutine
   character(2) :: restartdir = './'

! >>>>>>>>>>> The following are defined in the input file <<<<<<<<<<<<<<< !
!..Simulation domain (does not include ghost cells)
   real*8 :: xmin,xmax,ymin,ymax,zmin,zmax

!..Global physical grid (w/o ghost) size of the entire simulation (use powers of 2):
!..nx0G=xprocs*nx0, nx0 is the local number of x-points. Likewise for y
   integer :: nx0G,ny0G,nz0

!..number of processors to be used run-time: mpirun -np nprocs test_mpi.x
   integer :: xprocs,yprocs,nprocs

!..order of accuracy and reconstruction option
!..# of ghost cells on each side: ng = nord/2
!..reconop=1-PDM,2-TVD,3-WENO,4-MPWENO
!..stepop=1-AB,2-RK2,3-RK3
   integer :: nord,reconop,stepop

!..Output option: =0 no output, =1 HDF5, =2 write to .out files
!..IMPORTANT: =1 not available in this version
   integer :: outop

!..temporal inputs:
   real :: tau,tau0 ! tau can be fixed or changing every time-step
   integer :: nts, nframes, nrst, nfdump

!..other inputs:
   real :: cfl,pdmb,eps,alpha,beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!..No user-set parameters below this line:

   integer :: nt,ntt    ! time step and frame counters
   real*8 :: time,realt

!..local number of grid points along x and y
   integer :: nx0,ny0

!..grid size including ghost cells:
   integer :: nxG,nyG,nx,ny,nz
   integer :: ng

!..Grid spacing:
   real*8 :: dx0,dy0,dz0

!..spatial grid(global):
   real*8, dimension(:,:,:), allocatable :: x,y,z,dx,dy,dz
   real*8, dimension(:,:,:), allocatable :: xc,yc,zc,xi,yi,zi, &
                                            xj,yj,zj,xk,yk,zk
!..spatial grid(local):
   real*8, dimension(:,:,:), allocatable :: xblock,yblock,zblock, &
                                            dxblock,dyblock,dzblock, &
                                            xcblock,ycblock,zcblock, &
                                            xiblock,yiblock,ziblock, &
                                            xjblock,yjblock,zjblock, &
                                            xkblock,ykblock,zkblock

!..MPI variables
   integer :: ierr,ierror,irc
   integer :: istat(MPI_STATUS_SIZE)
   integer :: myrank,numprocs,nprovided
!..These are used by MPI_CART routines
   logical :: CART_BCs(2),reorder=.TRUE.
   integer :: COMM_CART,blockID(2)
   integer :: MPIgrid(2)
!..These are used for creating and using subcommunicators
   integer :: COMM_X,COMM_Y,remain(2),xID,yID,xcoord,ycoord
!   integer :: CARTgroup,xzgroup,yzgroup
!   integer, allocatable :: xzranks(:),yzranks(:)
!   integer :: xzCOMM,
   integer :: xprocsd2,yprocsd2
!..Rank of process to the left and right (along x).
   integer :: left, right
!..Rank of process in up and down (along y).
   integer :: up, down
!..BZ:following are not used but could be useful for future development
!..for data count
!   integer :: nxyz,nxy,nyz,ngyz,ngxz
!..Torque scheduler variables
!   integer :: args=0
!   character(len=8) :: cjobid, crestartid

contains
!********************************************
   subroutine read_input
   implicit none
!..Must call this to get the input params

!..Add here the namelists for parameters (first declared above) from input file:
   namelist/dom/xmin,xmax,ymin,ymax,zmin,zmax
   namelist/grid/nx0G,ny0G,nz0,xprocs,yprocs
   namelist/comset/nord,reconop,stepop,outop
   namelist/list/nts,tau,nframes,nrst,nfdump
   namelist/free/cfl,pdmb,eps,alpha,beta

   if (myrank.eq.0) then
      print*,'reading input file'
      open(unit=10,file="mhd3d2.in",status="old")
      read (10,dom)
      read (10,grid)
      read (10,comset)
      read (10,list)
      read (10,free)
      close(unit=10)
   endif

   call MPI_BCAST(xmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(xmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(ymin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(ymax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(zmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(zmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nx0G,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(ny0G,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nz0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(xprocs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(yprocs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(reconop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(stepop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(outop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nframes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nrst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(nfdump,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(cfl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(pdmb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   tau0=tau ! in case for fixed dt runs
   end subroutine read_input
!*********************************************
   subroutine comp_pars
   implicit none
!..Number of processes along each direction and total number of processes
   MPIgrid=(/yprocs,xprocs/)
!..Local number of grid points along x and z
   nx0=nx0G/xprocs
   ny0=ny0G/yprocs
!..number of ghost cells:
   ng=nord/2
!..global grid size including ghost cells:
   nxG=nx0G+nord
   nyG=ny0G+nord
!..local grid size including ghost cells:
   nx=nx0+nord
   ny=ny0+nord
   nz=nz0+nord
!..for MPI:
   nprocs=xprocs*yprocs
!   allocate(xzranks(xprocs),yzranks(yprocs))
!   nxyz=nx0*ny0*nz0
!   nxy=nx0*ny0
!   nyz=ny0*nz0
!   ngyz=ng*ny*nz
!   ngxz=ng*nx*nz

   end subroutine comp_pars
!*********************************************
   subroutine alloc_pars
   implicit none
   allocate(x(nxG+1,nyG+1,nz+1),y(nxG+1,nyG+1,nz+1),z(nxG+1,nyG+1,nz+1))
!   allocate(xc(nxG,nyG,nz),yc(nxG,nyG,nz),zc(nxG,nyG,nz))
!   allocate(dx(nxG,nyG,nz),dy(nxG,nyG,nz),dz(nxG,nyG,nz))
!   allocate(xi(nxG+1,nyG,nz),yi(nxG+1,nyG,nz),zi(nxG+1,nyG,nz))
!   allocate(xj(nxG,nyG+1,nz),yj(nxG,nyG+1,nz),zj(nxG,nyG+1,nz))
!   allocate(xk(nxG,nyG,nz+1),yk(nxG,nyG,nz+1),zk(nxG,nyG,nz+1))

   allocate(xblock(nx+1,ny+1,nz+1),yblock(nx+1,ny+1,nz+1),zblock(nx+1,ny+1,nz+1))
   allocate(xcblock(nx,ny,nz),ycblock(nx,ny,nz),zcblock(nx,ny,nz))
   allocate(dxblock(nx,ny,nz),dyblock(nx,ny,nz),dzblock(nx,ny,nz))
   allocate(xiblock(nx+1,ny,nz),yiblock(nx+1,ny,nz),ziblock(nx+1,ny,nz))
   allocate(xjblock(nx,ny+1,nz),yjblock(nx,ny+1,nz),zjblock(nx,ny+1,nz))
   allocate(xkblock(nx,ny,nz+1),ykblock(nx,ny,nz+1),zkblock(nx,ny,nz+1))
   end subroutine alloc_pars
!*********************************************
   subroutine dealloc_pars
   implicit none
   deallocate(x,y,z)!,dx,dy,dz)
!   deallocate(xc,yc,zc,xi,yi,zi,xj,yj,zj,xk,yk,zk)
   deallocate(xblock,yblock,zblock,dxblock,dyblock,dzblock)
   deallocate(xcblock,ycblock,zcblock,xiblock,yiblock,ziblock, &
              xjblock,yjblock,zjblock,xkblock,ykblock,zkblock)
   end subroutine dealloc_pars
!*********************************************
end module pars

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module fields
   use pars
   implicit none
!..premitive hydrodynamic and magnetic variables at cell center
   real*8, dimension(:,:,:), allocatable :: rho,vx,vy,vz,p,bx,by,bz
!..premitive magnetic and electric fields at cell faces
   real*8, dimension(:,:,:), allocatable :: bi,bj,bk,ei,ej,ek
!..other useful variables
   real*8, dimension(:,:,:), allocatable :: v,btot,rhovx,rhovy,rhovz,eng, &
                                            rho0,rhovx0,rhovy0,rhovz0,eng0
!..auxiliary variables
   real*8, dimension(:,:,:), allocatable :: lam,bsq
!..for stepping-general
   real*8, dimension(:,:,:), allocatable :: rhop,vxp,vyp,vzp,pp,bxp,byp,bzp, &
                                            bip,bjp,bkp
!..for stepping ab-predictor
   real*8, dimension(:,:,:), allocatable :: rhoh,vxh,vyh,vzh,ph,bxh,byh,bzh, &
                                            bih,bjh,bkh
!..for stepping rk
   real*8, dimension(:,:,:), allocatable :: rhovxp,rhovyp,rhovzp,rhovxh, &
                                            rhovyh,rhovzh,engp
!..flux and stress
   real*8, dimension(:,:,:), allocatable :: flux_rhox,flux_rhoy,flux_rhoz, &
                                            flux_rvxx,flux_rvxy,flux_rvxz, &
                                            flux_rvyx,flux_rvyy,flux_rvyz, &
                                            flux_rvzx,flux_rvzy,flux_rvzz, &
                                            flux_engx,flux_engy,flux_engz, &
                                            stress_xx,stress_xy,stress_xz, &
                                            stress_yx,stress_yy,stress_yz, &
                                            stress_zx,stress_zy,stress_zz
!..gaussian splict v
   real*8, dimension(:,:,:), allocatable :: vx0p,vx0n,vx1p,vx1n,vy0p,vy0n, &
                                            vy1p,vy1n,vz0p,vz0n,vz1p,vz1n
contains
   subroutine alloc_fields
   implicit none
   allocate(rho(nx,ny,nz),vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz),&
              p(nx,ny,nz),bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz))
   allocate(bi(nx+1,ny,nz),bj(nx,ny+1,nz),bk(nx,ny,nz+1),&
            ei(nx0,ny0+1,nz0+1),ej(nx0+1,ny0,nz0+1),ek(nx0+1,ny0+1,nz0))
   allocate(v(nx,ny,nz),btot(nx,ny,nz), &
            rhovx(nx,ny,nz),rhovy(nx,ny,nz),rhovz(nx,ny,nz),eng(nx,ny,nz))
   allocate(rho0(nx,ny,nz),rhovx0(nx,ny,nz),rhovy0(nx,ny,nz),rhovz0(nx,ny,nz), &
            eng0(nx,ny,nz))
   allocate(lam(nx,ny,nz),bsq(nx,ny,nz))
   allocate(rhop(nx,ny,nz),vxp(nx,ny,nz),vyp(nx,ny,nz),vzp(nx,ny,nz),pp(nx,ny,nz), &
            bxp(nx,ny,nz),byp(nx,ny,nz),bzp(nx,ny,nz), &
            bip(nx+1,ny,nz),bjp(nx,ny+1,nz),bkp(nx,ny,nz+1))
   allocate(rhoh(nx,ny,nz),vxh(nx,ny,nz),vyh(nx,ny,nz),vzh(nx,ny,nz),ph(nx,ny,nz), &
            bxh(nx,ny,nz),byh(nx,ny,nz),bzh(nx,ny,nz), &
            bih(nx+1,ny,nz),bjh(nx,ny+1,nz),bkh(nx,ny,nz+1))
   allocate(rhovxp(nx,ny,nz),rhovyp(nx,ny,nz),rhovzp(nx,ny,nz),rhovxh(nx,ny,nz), &
            rhovyh(nx,ny,nz),rhovzh(nx,ny,nz),engp(nx,ny,nz))
   allocate(flux_rhox(nx0+1,ny0,nz0),flux_rvxx(nx0+1,ny0,nz0),flux_rvyx(nx0+1,ny0,nz0), &
            flux_rvzx(nx0+1,ny0,nz0),flux_engx(nx0+1,ny0,nz0), &
            stress_xx(nx0+1,ny0,nz0),stress_yx(nx0+1,ny0,nz0),stress_zx(nx0+1,ny0,nz0), &
            flux_rhoy(nx0,ny0+1,nz0),flux_rvxy(nx0,ny0+1,nz0),flux_rvyy(nx0,ny0+1,nz0), &
            flux_rvzy(nx0,ny0+1,nz0),flux_engy(nx0,ny0+1,nz0), &
            stress_xy(nx0,ny0+1,nz0),stress_yy(nx0,ny0+1,nz0),stress_zy(nx0,ny0+1,nz0), &
            flux_rhoz(nx0,ny0,nz0+1),flux_rvxz(nx0,ny0,nz0+1),flux_rvyz(nx0,ny0,nz0+1), &
            flux_rvzz(nx0,ny0,nz0+1),flux_engz(nx0,ny0,nz0+1), &
            stress_xz(nx0,ny0,nz0+1),stress_yz(nx0,ny0,nz0+1),stress_zz(nx0,ny0,nz0+1))
   allocate(vx0p(nx,ny,nz),vx0n(nx,ny,nz),vx1p(nx,ny,nz),vx1n(nx,ny,nz), &
            vy0p(nx,ny,nz),vy0n(nx,ny,nz),vy1p(nx,ny,nz),vy1n(nx,ny,nz), &
            vz0p(nx,ny,nz),vz0n(nx,ny,nz),vz1p(nx,ny,nz),vz1n(nx,ny,nz))

   end subroutine alloc_fields
!***********************************
   subroutine dealloc_fields
   implicit none
   deallocate(rho,vx,vy,vz,p,bx,by,bz)
   deallocate(bi,bj,bk,ei,ej,ek)
   deallocate(v,btot,rhovx,rhovy,rhovz,eng)
   deallocate(rho0,rhovx0,rhovy0,rhovz0,eng0)
   deallocate(lam,bsq)
   deallocate(rhop,vxp,vyp,vzp,pp,bxp,byp,bzp,bip,bjp,bkp)
   deallocate(rhoh,vxh,vyh,vzh,ph,bxh,byh,bzh,bih,bjh,bkh)
   deallocate(rhovxp,rhovyp,rhovzp,rhovxh,rhovyh,rhovzh,engp)
   deallocate(flux_rhox,flux_rvxx,flux_rvyx,flux_rvzx,flux_engx, &
              stress_xx,stress_yx,stress_zx, &
              flux_rhoy,flux_rvxy,flux_rvyy,flux_rvzy,flux_engy, &
              stress_xy,stress_yy,stress_zy, &
              flux_rhoz,flux_rvxz,flux_rvyz,flux_rvzz,flux_engz, &
              stress_xz,stress_yz,stress_zz)
   deallocate(vx0p,vx0n,vx1p,vx1n,vy0p,vy0n,vy1p,vy1n,vz0p,vz0n,vz1p,vz1n)
end subroutine dealloc_fields
!************************************
end module fields

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module mpi_fillers
   use pars
   implicit none

contains
   subroutine fillx_center(aa,lnx,lny,lnz,xbc)
!..Fills the nord ghost cells in x (ng on either side) 
!  of an input array with unfilled x-ghost cells: aa(i,:,:) 
!  with i=1, ..., ng, nx-ng+1, ..., nx
   real*8, allocatable :: aa(:,:,:)
   real*8, allocatable :: sbufx(:,:,:),rbufx(:,:,:)
   integer :: i,lnx,lny,lnz,xbc,ncount

!..each proc has both of the needed ghost cell data for its neighbors:
!  may not need for corner ghost grids, slower but safter
   allocate(sbufx(ng,lny,lnz),rbufx(ng,lny,lnz))
   ncount=ng*lny*lnz
   
!..Send right most cells to the right
   sbufx=aa(lnx-nord+1:lnx-ng,:,:)
   call MPI_SENDRECV(sbufx,ncount,MPI_DOUBLE_PRECISION,right,0, &
                     rbufx,ncount,MPI_DOUBLE_PRECISION,left,0,COMM_X,istat,ierr)
   aa(1:ng,:,:)=rbufx

!..Send left most cells to the left
   sbufx=aa(ng+1:nord,:,:)
   call MPI_SENDRECV(sbufx,ncount,MPI_DOUBLE_PRECISION,left,1, &
                     rbufx,ncount,MPI_DOUBLE_PRECISION,right,1,COMM_X,istat,ierr)
   aa(lnx-ng+1:lnx,:,:)=rbufx

   deallocate(sbufx,rbufx)

!!..Fill the first and last process (along x) separately and enforce the
!!..corresponding symmetry
!   if (xID .eq. 0) then ! left x-boundary
!      if (xbc .eq. 0) then ! aa is even about x=-Lx/2
!         forall(i=1:ng) aa(i,:,:)=aa(nord+1-i,:,:)
!      elseif (xbc .eq. 1) then ! aa is odd about x=-Lx/2
!         forall(i=1:ng) aa(i,:,:)=-aa(nord+1-i,:,:)
!!     could continue define xbc=2,3 and so on, eg, extrapolation, perodic ...
!      endif
!   elseif (xID .eq. xprocs-1) then ! right x-boundary
!      if (xbc .eq. 0) then ! aa is even about x=Lx/2
!         forall(i=1:ng) aa(lnx+1-i,:,:)=aa(lnx-nord+i,:,:)
!      elseif (xbc .eq. 1) then ! aa is odd about x=Lx/2
!         forall(i=1:ng) aa(lnx+1-i,:,:)=-aa(lnx-nord+i,:,:)
!      endif
!   endif

   end subroutine fillx_center
!***************************************
   subroutine filly_center(aa,lnx,lny,lnz,ybc)
!..Fills the nord ghost cells in y (ng on either side) 
!  of an input array with unfilled y-ghost cells: aa(:,j,:) 
!  with i=1, ..., ng, ny-ng+1, ..., ny
   real*8, allocatable :: aa(:,:,:)
   real*8, allocatable :: sbufy(:,:,:),rbufy(:,:,:)
   integer :: j,lnx,lny,lnz,ybc,ncount

!..each proc has both of the needed ghost cell data for its neighbors:
   allocate(sbufy(lnx,ng,lnz),rbufy(lnx,ng,lnz))
   ncount=lnx*ng*lnz

!..Send up most cells to the up
   sbufy=aa(:,lny-nord+1:lny-ng,:)
   call MPI_SENDRECV(sbufy,ncount,MPI_DOUBLE_PRECISION,up,0, &
                     rbufy,ncount,MPI_DOUBLE_PRECISION,down,0,COMM_Y,istat,ierr)
   aa(:,1:ng,:)=rbufy

!..Send down most cells to the down
   sbufy=aa(:,ng+1:nord,:)
   call MPI_SENDRECV(sbufy,ncount,MPI_DOUBLE_PRECISION,down,1, &
                     rbufy,ncount,MPI_DOUBLE_PRECISION,up,1,COMM_Y,istat,ierr)
   aa(:,lny-ng+1:lny,:)=rbufy

   deallocate(sbufy,rbufy)

!!..Fill the first and last process (along y) separately and enforce the
!!..corresponding symmetry
!   if (yID .eq. 0) then ! left y-boundary
!      if (ybc .eq. 0) then ! aa is even about y=-Lx/2
!         forall(j=1:ng) aa(:,j,:)=aa(:,nord+1-j,:)
!      elseif (ybc .eq. 1) then ! aa is odd about x=-Lx/2
!         forall(j=1:ng) aa(:,j,:)=-aa(:,nord+1-j,:)
!!     could continue define ybc=2,3 and so on, eg, extrapolation, perodic ...
!      endif
!   elseif (yID .eq. yprocs-1) then ! right y-boundary
!      if (ybc .eq. 0) then ! aa is even about y=Ly/2
!         forall(j=1:ng) aa(:,lny+1-j,:)=aa(:,lny-nord+j,:)
!      elseif(ybc .eq. 1) then
!         forall(j=1:ng) aa(:,lny+1-j,:)=-aa(:,lny-nord+j,:)
!      endif
!   endif

   end subroutine filly_center
!************************************
   subroutine fillx_face(aa,lnx,lny,lnz,xbc)
!..Fills the nord ghost cells in x (ng on either side) 
!  of an input array with unfilled x-ghost cells: aa(i,:,:) 
!  with i=1, ..., ng, nx-ng+1, ..., nx+1
   real*8, allocatable :: aa(:,:,:)
   real*8, allocatable :: sbufx(:,:,:),rbufx(:,:,:)
   integer :: i,lnx,lny,lnz,xbc,ncount

!..each proc has both of the needed ghost cell data for its neighbors:
!  may not need for corner ghost grids, this way is slower but safter
   allocate(sbufx(ng,lny,lnz),rbufx(ng,lny,lnz))
   ncount=ng*lny*lnz

!..Send right most cells to the right
   sbufx=aa(lnx-nord:lnx-ng-1,:,:)
   call MPI_SENDRECV(sbufx,ncount,MPI_DOUBLE_PRECISION,right,0, &
                     rbufx,ncount,MPI_DOUBLE_PRECISION,left,0,COMM_X,istat,ierr)
   aa(1:ng,:,:)=rbufx

!..Send left most cells to the left
   sbufx=aa(ng+2:nord+1,:,:)
   call MPI_SENDRECV(sbufx,ncount,MPI_DOUBLE_PRECISION,left,1, &
                     rbufx,ncount,MPI_DOUBLE_PRECISION,right,1,COMM_X,istat,ierr)
   aa(lnx-ng+1:lnx,:,:)=rbufx
   deallocate(sbufx,rbufx)

!!..Fill the first and last process (along x) separately and enforce the
!!..corresponding symmetry
!   if (xID .eq. 0) then ! left x-boundary
!      if (xbc .eq. 0) then ! aa is even about x=-Lx/2
!         forall(i=1:ng) aa(i,:,:)=aa(nord+2-i,:,:)
!      elseif (xbc .eq. 1) then ! aa is odd about x=-Lx/2
!         forall(i=1:ng) aa(i,:,:)=-aa(nord+2-i,:,:)
!!     could continue define xbc=2,3 and so on, eg, extrapolation, perodic ...
!      endif
!   elseif (xID .eq. xprocs-1) then ! right x-boundary
!      if (xbc .eq. 0) then ! aa is even about x=Lx/2
!         forall(i=1:ng) aa(nx+2-i,:,:)=aa(nx0+i,:,:)
!      elseif (xbc .eq. 1) then ! aa is odd about x=Lx/2
!         forall(i=1:ng) aa(nx+2-i,:,:)=-aa(nx0+i,:,:)
!      endif
!   endif

   end subroutine fillx_face
!************************************
   subroutine filly_face(aa,lnx,lny,lnz,ybc)
!..Fills the nord ghost cells in y (ng on either side) 
!  of an input array with unfilled y-ghost cells: aa(:,j,:) 
!  with j=1, ..., ng, ny-ng+1, ..., ny+1
   real*8, allocatable :: aa(:,:,:)
   real*8, allocatable :: sbufy(:,:,:),rbufy(:,:,:)
   integer :: j,lnx,lny,lnz,ybc,ncount

!..each proc has both of the needed ghost cell data for its neighbors:
!  may not need for corner ghost grids, slower but safter
   allocate(sbufy(lnx,ng,lnz),rbufy(lnx,ng,lnz))
   ncount=lnx*ng*lnz

!..Send up most cells to the up
   sbufy=aa(:,lny-nord:lny-ng-1,:)
   call MPI_SENDRECV(sbufy,ncount,MPI_DOUBLE_PRECISION,up,0, &
                     rbufy,ncount,MPI_DOUBLE_PRECISION,down,0,COMM_Y,istat,ierr)
   aa(:,1:ng,:)=rbufy

!..Send down most cells to the down
   sbufy=aa(:,ng+2:nord+1,:)
   call MPI_SENDRECV(sbufy,ncount,MPI_DOUBLE_PRECISION,down,1, &
                     rbufy,ncount,MPI_DOUBLE_PRECISION,up,1,COMM_Y,istat,ierr)
   aa(:,lny-ng+1:lny,:)=rbufy

   deallocate(sbufy,rbufy)

!!..Fill the first and last process (along y) separately and enforce the
!!..corresponding symmetry
!   if (yID .eq. 0) then ! left y-boundary
!      if (ybc .eq. 0) then ! aa is even about y=-Ly/2
!         forall(j=1:ng) aa(:,j,:)=aa(:,nord+2-j,:)
!      elseif (ybc .eq. 1) then ! aa is odd about y=-Ly/2
!         forall(j=1:ng) aa(:,j,:)=-aa(:,nord+2-j,:)
!!     could continue define ybc=2,3 and so on, eg, extrapolation, perodic ...
!      endif
!   elseif (yID .eq. yprocs-1) then ! right y-boundary
!      if (ybc .eq. 0) then ! aa is even about y=Ly/2
!         forall(j=1:ng) aa(:,ny+2-j,:)=aa(:,ny0+j,:)
!      elseif(ybc .eq. 1) then
!         forall(j=1:ng) aa(:,ny+2-j,:)=-aa(:,ny0+j,:)
!      endif
!   endif

   end subroutine filly_face
!************************************
end module mpi_fillers

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module boundary
   use pars
   use fields
   use mpi_fillers
   implicit none

contains
   subroutine fill_hydro
   implicit none

!..fill_fuild is able to fill fluid quantites to be even or odd in 
!..specific dimension

!   call fill_fluid(rho,0,0,0)
!   call fill_fluid(p,0,0,0)
!   call fill_fluid(vx,0,0,0)
!   call fill_fluid(vy,0,0,0)
!   call fill_fluid(vz,0,0,0)

!..periodic boundary
   call fill_periodic(rho)
   call fill_periodic(p)
   call fill_periodic(vx)
   call fill_periodic(vy)
   call fill_periodic(vz)

   end subroutine fill_hydro
!***********************************
   subroutine fill_mag
   implicit none
   integer :: i

!..even reflection

   forall(i=1:ng)
      bi(i,:,:)=bi(nord+1-i,:,:)
      bi(:,i,:)=bi(:,nord+1-i,:)
      bi(:,:,i)=bi(:,:,nord+1-i)
      bi(nx+2-i,:,:)=bi(nx0+1+i,:,:)
      bi(:,ny+1-i,:)=bi(:,ny0+i,:)
      bi(:,:,nz+1-i)=bi(:,:,nz0+i)
   endforall

   forall(i=1:ng)
      bj(i,:,:)=bj(nord+1-i,:,:)
      bj(:,i,:)=bj(:,nord+1-i,:)
      bj(:,:,i)=bj(:,:,nord+1-i)
      bj(nx+1-i,:,:)=bj(nx0+i,:,:)
      bj(:,ny+2-i,:)=bj(:,ny0+1+i,:)
      bj(:,:,nz+1-i)=bj(:,:,nz0+i)
   endforall

   forall(i=1:ng)
      bk(i,:,:)=bk(nord+1-i,:,:)
      bk(:,i,:)=bk(:,nord+1-i,:)
      bk(:,:,i)=bk(:,:,nord+1-i)
      bk(nx+1-i,:,:)=bk(nx0+i,:,:)
      bk(:,ny+1-i,:)=bk(:,ny0+i,:)
      bk(:,:,nz+2-i)=bk(:,:,nz0+1+i)
   endforall

!..periodic
!
!   bi(1:ng,:,:)=bi(ng+nx0-3:ng+nx0,:,:)
!   bi(:,1:ng,:)=bi(:,ng+ny0-3:ng+ny0,:)
!   bi(:,:,1:ng)=bi(:,:,ng+nz0-3:ng+nz0)
!
!   bi(nx-ng+2:nx+1,:,:)=bi(ng+1:nord,:,:)
!   bi(:,ny-ng+1:ny,:)=bi(:,ng+1:nord,:)
!   bi(:,:,nz-ng+1:nz)=bi(:,:,ng+1:nord)
!
!   bj(1:ng,:,:)=bj(ng+nx0-3:ng+nx0,:,:)
!   bj(:,1:ng,:)=bj(:,ng+ny0-3:ng+ny0,:)
!   bj(:,:,1:ng)=bj(:,:,ng+nz0-3:ng+nz0)
!
!   bj(nx-ng+1:nx,:,:)=bj(ng+1:nord,:,:)
!   bj(:,ny-ng+2:ny+1,:)=bj(:,ng+1:nord,:)
!   bj(:,:,nz-ng+1:nz)=bj(:,:,ng+1:nord)
!
!   bk(1:ng,:,:)=bk(ng+nx0-3:ng+nx0,:,:)
!   bk(:,1:ng,:)=bk(:,ng+ny0-3:ng+ny0,:)
!   bk(:,:,1:ng)=bk(:,:,ng+nz0-3:ng+nz0)
!
!   bk(nx-ng+1:nx,:,:)=bk(ng+1:nord,:,:)
!   bk(:,ny-ng+1:ny,:)=bk(:,ng+1:nord,:)
!   bk(:,:,nz-ng+2:nz+1)=bk(:,:,ng+1:nord)
   
   end subroutine fill_mag
!************************************
   subroutine fill_fluid(a,xbc,ybc,zbc)

   implicit none
   real*8, dimension(:,:,:), allocatable :: a
   integer :: i,xbc,ybc,zbc

   if (xbc.eq.0) then
      forall(i=1:ng) a(i,:,:)=a(nord+1-i,:,:)
      forall(i=1:ng) a(nx+1-i,:,:)=a(nx0+i,:,:)
   elseif(xbc.eq.1) then
      forall(i=1:ng) a(i,:,:)=-a(nord+1-i,:,:)
      forall(i=1:ng) a(nx+1-i,:,:)=-a(nx0+i,:,:)
   endif

   if (ybc.eq.0) then
      forall(i=1:ng) a(:,i,:)=a(:,nord+1-i,:)
      forall(i=1:ng) a(:,ny+1-i,:)=a(:,ny0+i,:)
   elseif(ybc.eq.1) then
      forall(i=1:ng) a(:,i,:)=-a(:,nord+1-i,:)
      forall(i=1:ng) a(:,ny+1-i,:)=-a(:,ny0+i,:)
   endif

   if (zbc.eq.0) then
      forall(i=1:ng) a(:,:,i)=a(:,:,nord+1-i)
      forall(i=1:ng) a(:,:,nz+1-i)=a(:,:,nz0+i)
   elseif(zbc.eq.1) then
      forall(i=1:ng) a(:,:,i)=-a(:,:,nord+1-i)
      forall(i=1:ng) a(:,:,nz+1-i)=-a(:,:,nz0+i)
   endif

   end subroutine fill_fluid
!************************************
   subroutine fill_periodic(a)
   implicit none
   real*8, dimension(:,:,:), allocatable :: a
   integer :: i

   a(1:ng,:,:)=a(ng+nx0-3:ng+nx0,:,:)
   a(:,1:ng,:)=a(:,ng+ny0-3:ng+ny0,:)
   a(:,:,1:ng)=a(:,:,ng+nz0-3:ng+nz0)

   a(nx-ng+1:nx,:,:)=a(ng+1:nord,:,:)
   a(:,ny-ng+1:ny,:)=a(:,ng+1:nord,:)
   a(:,:,nz-ng+1:nz)=a(:,:,ng+1:nord)

   return

   end subroutine fill_periodic
!************************************
   subroutine fillz_even(a)
   implicit none
   real*8, dimension(:,:,:), allocatable :: a
   integer :: i

   forall(i=1:ng)
      a(:,:,i)=a(:,:,nord+1-i)
      a(:,:,nz+1-i)=a(:,:,nz0+i)
   endforall

   end subroutine fillz_even
!************************************
   subroutine fillz_even1(a)
   implicit none
   real*8, dimension(:,:,:), allocatable :: a
   integer :: i

   forall(i=1:ng)
      a(:,:,i)=a(:,:,ng+1)
      a(:,:,nz+1-i)=a(:,:,ng+1)
   endforall

   end subroutine fillz_even1
!************************************
   subroutine fillz_periodic(a)
   implicit none
   real*8, dimension(:,:,:), allocatable :: a
   integer :: i

   a(:,:,1:ng)=a(:,:,ng+nz0-3:ng+nz0)
   a(:,:,nz-ng+1:nz)=a(:,:,ng+1:nord)

   end subroutine fillz_periodic
!************************************
   subroutine fill_hydro_mpi
   implicit none
   call fillx_center(rho,nx,ny,nz,0)
   call fillx_center(p,nx,ny,nz,0)
   call fillx_center(vx,nx,ny,nz,1)
   call fillx_center(vy,nx,ny,nz,0)
   call fillx_center(vz,nx,ny,nz,0)

   call filly_center(rho,nx,ny,nz,0)
   call filly_center(p,nx,ny,nz,0)
   call filly_center(vx,nx,ny,nz,0)
   call filly_center(vy,nx,ny,nz,1)
   call filly_center(vz,nx,ny,nz,0)

!..for a 2d problem z dimension is a constant, fillz_even or fillz_periodic
!..shouldn't impact results
!   call fillz_even(rho)
!   call fillz_even(p)
!   call fillz_even(vx)
!   call fillz_even(vy)
!   call fillz_even(vz)

   call fillz_even1(rho)
   call fillz_even1(p)
   call fillz_even1(vx)
   call fillz_even1(vy)
   call fillz_even1(vz)

!..periodic
!   call fillz_periodic(rho)
!   call fillz_periodic(p)
!   call fillz_periodic(vx)
!   call fillz_periodic(vy)
!   call fillz_periodic(vz)

   end subroutine fill_hydro_mpi
!************************************
   subroutine fill_mag_mpi
   implicit none
   integer :: i
   
   call fillx_face(bi,nx+1,ny,nz,0)
   call fillx_center(bj,nx,ny+1,nz,0)
   call fillx_center(bk,nx,ny,nz+1,0)

   call filly_center(bi,nx+1,ny,nz,0)
   call filly_face(bj,nx,ny+1,nz,0)
   call filly_center(bk,nx,ny,nz+1,0)

!..z direction
   forall(i=1:ng)
!      bi(:,:,i)=bi(:,:,nord+1-i)
!      bi(:,:,nz+1-i)=bi(:,:,nz0+i)
!      bj(:,:,i)=bj(:,:,nord+1-i)
!      bj(:,:,nz+1-i)=bj(:,:,nz0+i)
!      bk(:,:,i)=bk(:,:,nord+2-i)
!      bk(:,:,nz+2-i)=bk(:,:,nz0+i)
      bi(:,:,i)=bi(:,:,ng+1)
      bi(:,:,nz+1-i)=bi(:,:,ng+1)
      bj(:,:,i)=bj(:,:,ng+1)
      bj(:,:,nz+1-i)=bj(:,:,ng+1)
      bk(:,:,i)=bk(:,:,ng+2)
      bk(:,:,nz+2-i)=bk(:,:,ng+1)
   endforall

!..periodic
!   bi(:,:,1:ng)=bi(:,:,ng+nz0-3:ng+nz0)
!   bi(:,:,ng+nz0+1:nz)=bi(:,:,ng+1:nord)
!   bj(:,:,1:ng)=bj(:,:,ng+nz0-3:ng+nz0)
!   bj(:,:,ng+nz0+1:nz)=bj(:,:,ng+1:nord)
!   bk(:,:,1:ng)=bk(:,:,ng+nz0-3:ng+nz0)
!   bk(:,:,ng+nz0+2:nz+1)=bk(:,:,ng+2:nord+1)

   end subroutine fill_mag_mpi
!************************************
!************************************
end module boundary
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!..two different type of reconstruction
!..1. cell center to face center (rho,v,p and so on)
!..2. face center to face corner (bi,bj,bk)
module reconstruct
   use pars
   use fields
   implicit none

contains
   subroutine reconstruct_3d(fleft,fright,f,lnx,lny,lnz,dir)
!..reconstruct from cell center to face center, or face center to corner
!..if dir=1,2,3,truncate fleft,fright to size(lnx,lny,lnz)
   implicit none
   real*8, dimension(:,:,:), allocatable :: fleft,fright,f
   real*8, dimension(:,:,:), allocatable :: fi
   integer :: lnx,lny,lnz,dir

   if(reconop.eq.1) then
      allocate(fi(lnx,lny,lnz))
      call interp_pdm(f,fi,lnx,lny,lnz,dir)
      call PDM2(fleft,fright,f,fi,lnx,lny,lnz,0,0,0,dir)
      deallocate(fi)
   elseif(reconop.eq.2) then
      allocate(fi(lnx,lny,lnz))
      call interp_tvd(f,fi,lnx,lny,lnz,dir)
!      call TVD(fleft,fright,f,lnx,lny,lnz,dir)
      call TVD(fleft,fright,f,fi,lnx,lny,lnz,0,0,0,dir)
      deallocate(fi)
   elseif(reconop.eq.3) then
      call WENO5(fleft,fright,f,lnx,lny,lnz,dir)
   elseif(reconop.eq.4) then
      call MPWENO5(fleft,fright,f,lnx,lny,lnz,dir)
   else
      print*,'Error-reconstruct called with incorrect opt'
   endif

   end subroutine reconstruct_3d
!***********************************
   subroutine interp_pdm(f,fi,lnx,lny,lnz,dir)
   implicit none
   real*8, dimension(:,:,:), allocatable :: f,fi
   integer :: lnx,lny,lnz,dir
   integer :: i,j,k

!..initial f0,f1,f2,f3 and ff based on dir
   if (dir.eq.1) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fi(i,j,k)=(-3.*f(ng+i-4,ng+j,ng+k)+29.*f(ng+i-3,ng+j,ng+k)- &
                    139.*f(ng+i-2,ng+j,ng+k)+533.*f(ng+i-1,ng+j,ng+k)+ &
                    533.*f(ng+i,ng+j,ng+k)-139.*f(ng+i+1,ng+j,ng+k)+ &
                    29.*f(ng+i+2,ng+j,ng+k)-3.*f(ng+i+3,ng+j,ng+k))/840.
      endforall
   elseif (dir.eq.2) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fi(i,j,k)=(-3.*f(ng+i,ng+j-4,ng+k)+29.*f(ng+i,ng+j-3,ng+k)- &
                    139.*f(ng+i,ng+j-2,ng+k)+533.*f(ng+i,ng+j-1,ng+k)+ &
                    533.*f(ng+i,ng+j,ng+k)-139.*f(ng+i,ng+j+1,ng+k)+ &
                    29.*f(ng+i,ng+j+2,ng+k)-3.*f(ng+i,ng+j+3,ng+k))/840.
      endforall
   elseif (dir.eq.3) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fi(i,j,k)=(-3.*f(ng+i,ng+j,ng+k-4)+29.*f(ng+i,ng+j,ng+k-3)- &
                    139.*f(ng+i,ng+j,ng+k-2)+533.*f(ng+i,ng+j,ng+k-1)+ &
                    533.*f(ng+i,ng+j,ng+k)-139.*f(ng+i,ng+j,ng+k+1)+ &
                    29.*f(ng+i,ng+j,ng+k+2)-3.*f(ng+i,ng+j,ng+k+3))/840.
      endforall
   elseif (dir.eq.4) then
!..recon in x, but leave y z untrucated
!  used for first step in getting ek (next step is recon in y)
!  or for second step in getting ej (previous step is recon in z)
      forall(i=1:lnx)
         fi(i,:,:)=(-3.*f(ng+i-4,:,:)+29.*f(ng+i-3,:,:)- &
                    139.*f(ng+i-2,:,:)+533.*f(ng+i-1,:,:)+ &
                    533.*f(ng+i,:,:)-139.*f(ng+i+1,:,:)+ &
                    29.*f(ng+i+2,:,:)-3.*f(ng+i+3,:,:))/840.
       endforall
    elseif (dir.eq.5) then
!..recon in y, but leave x z untrucated 
!  used for first step in getting ei (next step is recon in z)
!  or for second step in getting ek (previous step is recon in x)
      forall(j=1:lny)
         fi(:,j,:)=(-3.*f(:,ng+j-4,:)+29.*f(:,ng+j-3,:)- &
                    139.*f(:,ng+j-2,:)+533.*f(:,ng+j-1,:)+ &
                    533.*f(:,ng+j,:)-139.*f(:,ng+j+1,:)+ &
                    29.*f(:,ng+j+2,:)-3.*f(:,ng+j+3,:))/840.
       endforall
    elseif (dir.eq.6) then
!..recon in z, but leave x y untrucated 
!  used for first step in getting ej (next step is recon in x)
!  or for second step in getting ei (previous step is recon in y)
      forall(k=1:lnz)
         fi(:,:,k)=(-3.*f(:,:,ng+k-4)+29.*f(:,:,ng+k-3)- &
                    139.*f(:,:,ng+k-2)+533.*f(:,:,ng+k-1)+ &
                    533.*f(:,:,ng+k)-139.*f(:,:,ng+k+1)+ &
                    29.*f(:,:,ng+k+2)-3.*f(:,:,ng+k+3))/840.
      endforall
   endif

   return

   end subroutine interp_pdm
!*****************************************
   subroutine interp2_pdm(f,fi,dir)
!..this subroutine is used only for interp velocity to corner
   implicit none
   real*8, dimension(:,:,:), allocatable :: f,fi,g
   integer :: dir
   integer :: i,j,k

!..for reconop=1 (PDM)
   if (dir.eq.1) then ! used for ei calculation, thus fi is (nx0,ny0+1,nz0+1)
                      ! first interp along y then along z 
      allocate(g(nx0,ny0+1,nz))
      forall(i=1:nx0,j=1:ny0+1)
         g(i,j,:)=(-3.*f(ng+i,ng+j-4,:)+29.*f(ng+i,ng+j-3,:)- &
                    139.*f(ng+i,ng+j-2,:)+533.*f(ng+i,ng+j-1,:)+ &
                    533.*f(ng+i,ng+j,:)-139.*f(ng+i,ng+j+1,:)+ &
                    29.*f(ng+i,ng+j+2,:)-3.*f(ng+i,ng+j+3,:))/840.
      endforall
      forall(k=1:nz0+1)
         fi(:,:,k)=(-3.*g(:,:,ng+k-4)+29.*g(:,:,ng+k-3)- &
                    139.*g(:,:,ng+k-2)+533.*g(:,:,ng+k-1)+ &
                    533.*g(:,:,ng+k)-139.*g(:,:,ng+k+1)+ &
                    29.*g(:,:,ng+k+2)-3.*g(:,:,ng+k+3))/840.
      endforall
      deallocate(g)

   elseif (dir.eq.2) then ! used for ej calculateion, 
                          ! thus, fi is (nx0+1,ny0,nz0+1), first along z then x
      allocate(g(nx,ny0,nz0+1))
      forall(j=1:ny0,k=1:nz0+1)
         g(:,j,k)=(-3.*f(:,ng+j,ng+k-4)+29.*f(:,ng+j,ng+k-3)- &
                    139.*f(:,ng+j,ng+k-2)+533.*f(:,ng+j,ng+k-1)+ &
                    533.*f(:,ng+j,ng+k)-139.*f(:,ng+j,ng+k+1)+ &
                    29.*f(:,ng+j,ng+k+2)-3.*f(:,ng+j,ng+k+3))/840.
      endforall
      forall(i=1:nx0+1)
         fi(i,:,:)=(-3.*g(ng+i-4,:,:)+29.*g(ng+i-3,:,:)- &
                    139.*g(ng+i-2,:,:)+533.*g(ng+i-1,:,:)+ &
                    533.*g(ng+i,:,:)-139.*g(ng+i+1,:,:)+ &
                    29.*g(ng+i+2,:,:)-3.*g(ng+i+3,:,:))/840.
      endforall
      deallocate(g)

   elseif (dir.eq.3) then ! used for ek calculateion, 
                          ! thus, fi is (nx0+1,ny0+1,nz0), first along x then y
      allocate(g(nx0+1,ny,nz0))
      forall(i=1:nx0+1,k=1:nz0)
         g(i,:,k)=(-3.*f(ng+i-4,:,ng+k)+29.*f(ng+i-3,:,ng+k)- &
                    139.*f(ng+i-2,:,ng+k)+533.*f(ng+i-1,:,ng+k)+ &
                    533.*f(ng+i,:,ng+k)-139.*f(ng+i+1,:,ng+k)+ &
                    29.*f(ng+i+2,:,ng+k)-3.*f(ng+i+3,:,ng+k))/840.
      endforall
      forall(j=1:ny0+1)
         fi(:,j,:)=(-3.*g(:,ng+j-4,:)+29.*g(:,ng+j-3,:)- &
                    139.*g(:,ng+j-2,:)+533.*g(:,ng+j-1,:)+ &
                    533.*g(:,ng+j,:)-139.*g(:,ng+j+1,:)+ &
                    29.*g(:,ng+j+2,:)-3.*g(:,ng+j+3,:))/840.
      endforall
      deallocate(g)

   endif

   return

   end subroutine interp2_pdm
!*****************************************
   subroutine interp_tvd(f,fi,lnx,lny,lnz,dir)
   implicit none
   real*8, dimension(:,:,:), allocatable :: f,fi
   integer :: lnx,lny,lnz,dir
   integer :: i,j,k

!..initial f0,f1,f2,f3 and ff based on dir
   if (dir.eq.1) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fi(i,j,k)=(f(ng+i-1,ng+j,ng+k)+f(ng+i,ng+j,ng+k))/2.
      endforall
   elseif (dir.eq.2) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fi(i,j,k)=(f(ng+i,ng+j-1,ng+k)+f(ng+i,ng+j,ng+k))/2.
      endforall
   elseif (dir.eq.3) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fi(i,j,k)=(f(ng+i,ng+j,ng+k-1)+f(ng+i,ng+j,ng+k))/2.
      endforall
   elseif (dir.eq.4) then
!..recon in x, but leave y z untrucated
!  used for first step in getting ek (next step is recon in y)
!  or for second step in getting ej (previous step is recon in z)
      forall(i=1:lnx)
         fi(i,:,:)=(f(ng+i-1,:,:)+f(ng+i,:,:))/2.
       endforall
    elseif (dir.eq.5) then
!..recon in y, but leave x z untrucated 
!  used for first step in getting ei (next step is recon in z)
!  or for second step in getting ek (previous step is recon in x)
      forall(j=1:lny)
         fi(:,j,:)=(f(:,ng+j-1,:)+f(:,ng+j,:))/2.
       endforall
    elseif (dir.eq.6) then
!..recon in z, but leave x y untrucated 
!  used for first step in getting ej (next step is recon in x)
!  or for second step in getting ei (previous step is recon in y)
      forall(k=1:lnz)
         fi(:,:,k)=(f(:,:,ng+k-1)+f(:,:,ng+k))/2.
      endforall
   endif

   return

   end subroutine interp_tvd
!*****************************************
   subroutine interp2_tvd(f,fi,dir) !! need to be changed
!..this subroutine is used only for interp velocity to corner
   implicit none
   real*8, dimension(:,:,:), allocatable :: f,fi
   integer :: dir
   integer :: i,j,k

!..for reconop=2 (TVD)
   if (dir.eq.1) then ! used for ei calculation, fi is (nx0,ny0+1,nz0+1)
      forall(i=1:nx0,j=1:ny0+1,k=1:nz0+1)
         fi(i,j,k)=(f(ng+i,ng+j-1,ng+k-1)+f(ng+i,ng+j,ng+k-1)+ &
                    f(ng+i,ng+j,ng+k)+f(ng+i,ng+j-1,ng+k))/4.
      endforall

   elseif (dir.eq.2) then ! used for ej calculateion, fi is (nx0+1,ny0,nz0+1)
      forall(i=1:nx0+1,j=1:ny0,k=1:nz0+1)
         fi(i,j,k)=(f(ng+i-1,ng+j,ng+k-1)+f(ng+i-1,ng+j,ng+k)+ &
                    f(ng+i,ng+j,ng+k)+f(ng+i,ng+j,ng+k-1))/4.
      endforall

   elseif (dir.eq.3) then ! used for ek calculateion, fi is (nx0+1,ny0+1,nz0)
      forall(i=1:nx0+1,j=1:ny0+1,k=1:nz0)
         fi(i,j,k)=(f(ng+i-1,ng+j-1,ng+k)+f(ng+i,ng+j-1,ng+k)+ &
                    f(ng+i,ng+j,ng+k)+f(ng+i-1,ng+j,ng+k))/4.
      endforall

   endif

   return

   end subroutine interp2_tvd
!*****************************************
   subroutine PDM2(fleft,fright,f,fi,lnx,lny,lnz,snx,sny,snz,dir)

   implicit none
   real*8, dimension(:,:,:), allocatable :: fleft,fright,f,fi,fii
   integer :: lnx,lny,lnz,snx,sny,snz,dir
   real*8, dimension(:,:,:), allocatable :: fm2,fm1,f0,fp1,dfm,df0,dfp
   real*8, dimension(:,:,:), allocatable :: sm,s0,sp,qm,qp,dfl,dfr
   real*8, dimension(:,:,:), allocatable :: maxf,minf,dum
   integer :: i,j,k

!..allocate local variable based on input parameters
   allocate(fm2(lnx,lny,lnz),fm1(lnx,lny,lnz),f0(lnx,lny,lnz), &
            fp1(lnx,lny,lnz),fii(lnx,lny,lnz),dfm(lnx,lny,lnz), &
            df0(lnx,lny,lnz),dfp(lnx,lny,lnz),sm(lnx,lny,lnz), &
            s0(lnx,lny,lnz),sp(lnx,lny,lnz),qm(lnx,lny,lnz), &
            qp(lnx,lny,lnz),dfl(lnx,lny,lnz),dfr(lnx,lny,lnz), &
            maxf(lnx,lny,lnz),minf(lnx,lny,lnz),dum(lnx,lny,lnz))

   fii=fi
!..initial f0,f1,f2,f3 and ff based on dir
   if (dir.eq.1) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm2(i,j,k)=f(ng-2+i+snx,ng+j+sny,ng+k+snz)
         fm1(i,j,k)=f(ng-1+i+snx,ng+j+sny,ng+k+snz)
         f0(i,j,k)=f(ng+i+snx,ng+j+sny,ng+k+snz)
         fp1(i,j,k)=f(ng+1+i+snx,ng+j+sny,ng+k+snz)
      endforall
   elseif (dir.eq.2) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm2(i,j,k)=f(ng+i+snx,ng-2+j+sny,ng+k+snz)
         fm1(i,j,k)=f(ng+i+snx,ng-1+j+sny,ng+k+snz)
         f0(i,j,k)=f(ng+i+snx,ng+j+sny,ng+k+snz)
         fp1(i,j,k)=f(ng+i+snx,ng+1+j+sny,ng+k+snz)
      endforall
   elseif (dir.eq.3) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm2(i,j,k)=f(ng+i+snx,ng+j+sny,ng-2+k+snz)
         fm1(i,j,k)=f(ng+i+snx,ng+j+sny,ng-1+k+snz)
         f0(i,j,k)=f(ng+i+snx,ng+j+sny,ng+k+snz)
         fp1(i,j,k)=f(ng+i+snx,ng+j+sny,ng+1+k+snz)
      endforall
   elseif (dir.eq.4) then
!..recon in x, but leave y z untrucated
!  used for first step in getting ek (next step is recon in y)
!  or for second step in getting ej (previous step is recon in z)
      forall(i=1:lnx)
         fm2(i,:,:)=f(ng-2+i+snx,:,:)
         fm1(i,:,:)=f(ng-1+i+snx,:,:)
         f0(i,:,:)=f(ng+i+snx,:,:)
         fp1(i,:,:)=f(ng+1+i+snx,:,:)
       endforall
    elseif (dir.eq.5) then
!..recon in y, but leave x z untrucated 
!  used for first step in getting ei (next step is recon in z)
!  or for second step in getting ek (previous step is recon in x)
      forall(j=1:lny)
         fm2(:,j,:)=f(:,ng-2+j+sny,:)
         fm1(:,j,:)=f(:,ng-1+j+sny,:)
         f0(:,j,:)=f(:,ng+j+sny,:)
         fp1(:,j,:)=f(:,ng+1+j+sny,:)
       endforall
    elseif (dir.eq.6) then
!..recon in z, but leave x y untrucated 
!  used for first step in getting ej (next step is recon in x)
!  or for second step in getting ei (previous step is recon in y)
      forall(k=1:lnz)
         fm2(:,:,k)=f(:,:,ng-2+k+snz)
         fm1(:,:,k)=f(:,:,ng-1+k+snz)
         f0(:,:,k)=f(:,:,ng+k+snz)
         fp1(:,:,k)=f(:,:,ng+1+k+snz)
      endforall
   endif
  
   maxf=max(fm1,f0)
   minf=min(fm1,f0)

   fii=max(minf,min(fii,maxf))
   dfm=pdmb*(fm1-fm2)
   df0=pdmb*(f0-fm1)
   dfp=pdmb*(fp1-f0)

   dum=1.
   sm=sign(dum,dfm)
   s0=sign(dum,df0)
   sp=sign(dum,dfp)
   
   dfm=abs(dfm)
   df0=abs(df0)
   dfp=abs(dfp)

   qm=abs(sm+s0)
   qp=abs(s0+sp)

   dfl=fii-fm1
   dfr=f0-fii

   fleft=fii-s0*max(0.,abs(dfl)-qm*dfm)
   fright=fii+s0*max(0.,abs(dfr)-qp*dfp)

   deallocate(fm2,fm1,f0,fp1,fii)
   deallocate(dfm,df0,dfp,sm,s0,sp,qm,qp,dfl,dfr,maxf,minf,dum)

   return
   end subroutine PDM2
!***********************************
   subroutine TVD(fleft,fright,f,fi,lnx,lny,lnz,snx,sny,snz,dir)

   implicit none
   real*8, dimension(:,:,:), allocatable :: fleft,fright,f,fi
   integer :: lnx,lny,lnz,snx,sny,snz,dir
   real*8, dimension(:,:,:), allocatable :: fm2,fm1,f0,fp1,a,b
   integer :: i,j,k

!..allocate local variable based on input parameters
   allocate(fm2(lnx,lny,lnz),fm1(lnx,lny,lnz),f0(lnx,lny,lnz), &
            fp1(lnx,lny,lnz),a(lnx,lny,lnz),b(lnx,lny,lnz))

!..initial f0,f1,f2,f3 and ff based on dir
   if (dir.eq.1) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm2(i,j,k)=f(ng-2+i+snx,ng+j+sny,ng+k+snz)
         fm1(i,j,k)=f(ng-1+i+snx,ng+j+sny,ng+k+snz)
         f0(i,j,k)=f(ng+i+snx,ng+j+sny,ng+k+snz)
         fp1(i,j,k)=f(ng+1+i+snx,ng+j+sny,ng+k+snz)
      endforall
   elseif (dir.eq.2) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm2(i,j,k)=f(ng+i+snx,ng-2+j+sny,ng+k+snz)
         fm1(i,j,k)=f(ng+i+snx,ng-1+j+sny,ng+k+snz)
         f0(i,j,k)=f(ng+i+snx,ng+j+sny,ng+k+snz)
         fp1(i,j,k)=f(ng+i+snx,ng+1+j+sny,ng+k+snz)
      endforall
   elseif (dir.eq.3) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm2(i,j,k)=f(ng+i+snx,ng+j+sny,ng-2+k+snz)
         fm1(i,j,k)=f(ng+i+snx,ng+j+sny,ng-1+k+snz)
         f0(i,j,k)=f(ng+i+snx,ng+j+sny,ng+k+snz)
         fp1(i,j,k)=f(ng+i+snx,ng+j+sny,ng+1+k+snz)
      endforall
   elseif (dir.eq.4) then
!..recon in x, but leave y z untrucated
!  used for first step in getting ek (next step is recon in y)
!  or for second step in getting ej (previous step is recon in z)
      forall(i=1:lnx)
         fm2(i,:,:)=f(ng-2+i+snx,:,:)
         fm1(i,:,:)=f(ng-1+i+snx,:,:)
         f0(i,:,:)=f(ng+i+snx,:,:)
         fp1(i,:,:)=f(ng+1+i+snx,:,:)
      endforall
    elseif (dir.eq.5) then
!..recon in y, but leave x z untrucated 
!  used for first step in getting ei (next step is recon in z)
!  or for second step in getting ek (previous step is recon in x)
      forall(j=1:lny)
         fm2(:,j,:)=f(:,ng-2+j+sny,:)
         fm1(:,j,:)=f(:,ng-1+j+sny,:)
         f0(:,j,:)=f(:,ng+j+sny,:)
         fp1(:,j,:)=f(:,ng+1+j+sny,:)
       endforall
    elseif (dir.eq.6) then
!..recon in z, but leave x y untrucated 
!  used for first step in getting ej (next step is recon in x)
!  or for second step in getting ei (previous step is recon in y)
      forall(k=1:lnz)
         fm2(:,:,k)=f(:,:,ng-2+k+snz)
         fm1(:,:,k)=f(:,:,ng-1+k+snz)
         f0(:,:,k)=f(:,:,ng+k+snz)
         fp1(:,:,k)=f(:,:,ng+1+k+snz)
      endforall
   endif

   a=(fm1-fm2)/(f0-fm1+eps)
   b=(a+abs(a))/(1.+abs(a)) ! Van-Leer Limiter
   fleft=fm1-b*(fm1-fi)

   a=(f0-fp1)/(fm1-f0+eps)
   b=(a+abs(a))/(1.+abs(a))
   fright=f0-b*(f0-fi)

   deallocate(fm2,fm1,f0,fp1,a,b)

   end subroutine TVD
!***********************************
   subroutine WENO5(fleft,fright,f,lnx,lny,lnz,dir)

   implicit none
   real*8, dimension(:,:,:), allocatable :: fleft,fright,f
   integer :: lnx,lny,lnz,dir
   real*8, dimension(:,:,:), allocatable :: fm3,fm2,fm1,f0,fp1,fp2, &
                                            h0,h1,h2,s0,s1,s2,w0,w1,w2
   integer :: i,j,k

!..allocate local variable based on input parameters
   allocate(fm3(lnx,lny,lnz),fm2(lnx,lny,lnz),fm1(lnx,lny,lnz), &
            f0(lnx,lny,lnz),fp1(lnx,lny,lnz),fp2(lnx,lny,lnz),  &
            h0(lnx,lny,lnz),h1(lnx,lny,lnz),h2(lnx,lny,lnz),    &
            s0(lnx,lny,lnz),s1(lnx,lny,lnz),s2(lnx,lny,lnz),    &
            w0(lnx,lny,lnz),w1(lnx,lny,lnz),w2(lnx,lny,lnz)  )

!..initial f0,f1,f2,f3 and ff based on dir
   if (dir.eq.1) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm3(i,j,k)=f(ng-3+i,ng+j,ng+k)
         fm2(i,j,k)=f(ng-2+i,ng+j,ng+k)
         fm1(i,j,k)=f(ng-1+i,ng+j,ng+k)
         f0(i,j,k)=f(ng+i,ng+j,ng+k)
         fp1(i,j,k)=f(ng+1+i,ng+j,ng+k)
         fp2(i,j,k)=f(ng+2+i,ng+j,ng+k)
      endforall
   elseif (dir.eq.2) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm3(i,j,k)=f(ng+i,ng-3+j,ng+k)
         fm2(i,j,k)=f(ng+i,ng-2+j,ng+k)
         fm1(i,j,k)=f(ng+i,ng-1+j,ng+k)
         f0(i,j,k)=f(ng+i,ng+j,ng+k)
         fp1(i,j,k)=f(ng+i,ng+1+j,ng+k)
         fp2(i,j,k)=f(ng+i,ng+2+j,ng+k)
      endforall
   elseif (dir.eq.3) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm3(i,j,k)=f(ng+i,ng+j,ng-3+k)
         fm2(i,j,k)=f(ng+i,ng+j,ng-2+k)
         fm1(i,j,k)=f(ng+i,ng+j,ng-1+k)
         f0(i,j,k)=f(ng+i,ng+j,ng+k)
         fp1(i,j,k)=f(ng+i,ng+j,ng+1+k)
         fp2(i,j,k)=f(ng+i,ng+j,ng+2+k)
      endforall
   elseif (dir.eq.4) then
!..recon in x, but leave y z untrucated
!  used for first step in getting ek (next step is recon in y)
!  or for second step in getting ej (previous step is recon in z)
      forall(i=1:lnx)
         fm3(i,:,:)=f(ng-3+i,:,:)
         fm2(i,:,:)=f(ng-2+i,:,:)
         fm1(i,:,:)=f(ng-1+i,:,:)
         f0(i,:,:)=f(ng+i,:,:)
         fp1(i,:,:)=f(ng+1+i,:,:)
         fp2(i,:,:)=f(ng+2+i,:,:)
      endforall
    elseif (dir.eq.5) then
!..recon in y, but leave x z untrucated 
!  used for first step in getting ei (next step is recon in z)
!  or for second step in getting ek (previous step is recon in x)
      forall(j=1:lny)
         fm3(:,j,:)=f(:,ng-3+j,:)
         fm2(:,j,:)=f(:,ng-2+j,:)
         fm1(:,j,:)=f(:,ng-1+j,:)
         f0(:,j,:)=f(:,ng+j,:)
         fp1(:,j,:)=f(:,ng+1+j,:)
         fp2(:,j,:)=f(:,ng+2+j,:)
       endforall
    elseif (dir.eq.6) then
!..recon in z, but leave x y untrucated 
!  used for first step in getting ej (next step is recon in x)
!  or for second step in getting ei (previous step is recon in y)
      forall(k=1:lnz)
         fm3(:,:,k)=f(:,:,ng-3+k)
         fm2(:,:,k)=f(:,:,ng-2+k)
         fm1(:,:,k)=f(:,:,ng-1+k)
         f0(:,:,k)=f(:,:,ng+k)
         fp1(:,:,k)=f(:,:,ng+1+k)
         fp2(:,:,k)=f(:,:,ng+2+k)
      endforall
   endif

!..for left component
   h0=(2.*fm3-7.*fm2+11.*fm1)/6.
   h1=(-fm2+5.*fm1+2.*f0)/6.
   h2=(2.*fm1+5.*f0-fp1)/6.
   
   s0=13./12.*(fm3-2.*fm2+fm1)**2+.25*(fm3-4.*fm2+3.*fm1)**2
   s1=13./12.*(fm2-2.*fm1+f0)**2+.25*(fm2-f0)**2
   s2=13./12.*(fm1-2.*f0+fp1)**2+.25*(3.*fm1-4.*f0+fp1)**2
   
   s0=.1/(eps+s0)**2
   s1=.6/(eps+s1)**2
   s2=.3/(eps+s2)**2

!..fm3 is reused here to eliminate operations
   fm3=s0+s1+s2
   w0=s0/fm3
   w1=s1/fm3
   w2=s2/fm3

   fleft=w0*h0+w1*h1+w2*h2

!..for right component
   h0=(2.*fp2-7.*fp1+11.*f0)/6.
   h1=(-fp1+5.*f0+2.*fm1)/6.
   h2=(2.*f0+5.*fm1-fm2)/6.

   s0=13./12.*(fp2-2.*fp1+f0)**2+.25*(fp2-4.*fp1+3.*f0)**2
   s1=13./12.*(fp1-2.*f0+fm1)**2+.25*(fp1-fm1)**2
   s2=13./12.*(f0-2.*fm1+fm2)**2+.25*(3.*f0-4.*fm1+fm2)**2

   s0=.1/(eps+s0)**2
   s1=.6/(eps+s1)**2
   s2=.3/(eps+s2)**2

!..fp2 is reused here to eliminate operations
   fp2=s0+s1+s2
   w0=s0/fp2
   w1=s1/fp2
   w2=s2/fp2

   fright=w0*h0+w1*h1+w2*h2

   deallocate(fm3,fm2,fm1,f0,fp1,fp2,h0,h1,h2,s0,s1,s2,w0,w1,w2)

   end subroutine WENO5
!***********************************
   function minmod(a,b,lx,ly,lz)
   implicit none
   real*8, allocatable,dimension(:,:,:) :: a,b,minmod,dum
   integer :: lx,ly,lz

   allocate(minmod(lx,ly,lz),dum(lx,ly,lz))
   dum=1.
   minmod=.5*(sign(dum,a)+sign(dum,b))*min(abs(a),abs(b))
   deallocate (dum)

   end function minmod
!***********************************
   function median(a,b,c,lx,ly,lz)
   implicit none
   real*8, allocatable,dimension(:,:,:) :: a,b,c,median,dum1,dum2
   integer :: lx,ly,lz

   allocate(median(lx,ly,lz),dum1(lx,ly,lz),dum2(lx,ly,lz))
   dum1=b-a
   dum2=c-a
   median=a+minmod(dum1,dum2,lx,ly,lz)
   deallocate(dum1,dum2)

   end function median
!***********************************
   subroutine MPWENO5(fleft,fright,f,lnx,lny,lnz,dir)

   implicit none
   real*8, dimension(:,:,:), allocatable :: fleft,fright,f
   integer :: lnx,lny,lnz,dir
   real*8, dimension(:,:,:), allocatable :: fm3,fm2,fm1,f0,fp1,fp2, &
                                            h0,h1,h2,s0,s1,s2,w0,w1,w2
   real*8, dimension(:,:,:), allocatable :: dum1,dum2
   integer :: i,j,k

!..allocate local variable based on input parameters
   allocate(fm3(lnx,lny,lnz),fm2(lnx,lny,lnz),fm1(lnx,lny,lnz), &
            f0(lnx,lny,lnz),fp1(lnx,lny,lnz),fp2(lnx,lny,lnz),  &
            h0(lnx,lny,lnz),h1(lnx,lny,lnz),h2(lnx,lny,lnz),    &
            s0(lnx,lny,lnz),s1(lnx,lny,lnz),s2(lnx,lny,lnz),    &
            w0(lnx,lny,lnz),w1(lnx,lny,lnz),w2(lnx,lny,lnz)  )
   allocate(dum1(lnx,lny,lnz),dum2(lnx,lny,lnz))

!..initial f0,f1,f2,f3 and ff based on dir
   if (dir.eq.1) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm3(i,j,k)=f(ng-3+i,ng+j,ng+k)
         fm2(i,j,k)=f(ng-2+i,ng+j,ng+k)
         fm1(i,j,k)=f(ng-1+i,ng+j,ng+k)
         f0(i,j,k)=f(ng+i,ng+j,ng+k)
         fp1(i,j,k)=f(ng+1+i,ng+j,ng+k)
         fp2(i,j,k)=f(ng+2+i,ng+j,ng+k)
      endforall
   elseif (dir.eq.2) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm3(i,j,k)=f(ng+i,ng-3+j,ng+k)
         fm2(i,j,k)=f(ng+i,ng-2+j,ng+k)
         fm1(i,j,k)=f(ng+i,ng-1+j,ng+k)
         f0(i,j,k)=f(ng+i,ng+j,ng+k)
         fp1(i,j,k)=f(ng+i,ng+1+j,ng+k)
         fp2(i,j,k)=f(ng+i,ng+2+j,ng+k)
      endforall
   elseif (dir.eq.3) then
      forall(i=1:lnx,j=1:lny,k=1:lnz)
         fm3(i,j,k)=f(ng+i,ng+j,ng-3+k)
         fm2(i,j,k)=f(ng+i,ng+j,ng-2+k)
         fm1(i,j,k)=f(ng+i,ng+j,ng-1+k)
         f0(i,j,k)=f(ng+i,ng+j,ng+k)
         fp1(i,j,k)=f(ng+i,ng+j,ng+1+k)
         fp2(i,j,k)=f(ng+i,ng+j,ng+2+k)
      endforall
   elseif (dir.eq.4) then
!..recon in x, but leave y z untrucated
!  used for first step in getting ek (next step is recon in y)
!  or for second step in getting ej (previous step is recon in z)
      forall(i=1:lnx)
         fm3(i,:,:)=f(ng-3+i,:,:)
         fm2(i,:,:)=f(ng-2+i,:,:)
         fm1(i,:,:)=f(ng-1+i,:,:)
         f0(i,:,:)=f(ng+i,:,:)
         fp1(i,:,:)=f(ng+1+i,:,:)
         fp2(i,:,:)=f(ng+2+i,:,:)
      endforall
    elseif (dir.eq.5) then
!..recon in y, but leave x z untrucated 
!  used for first step in getting ei (next step is recon in z)
!  or for second step in getting ek (previous step is recon in x)
      forall(j=1:lny)
         fm3(:,j,:)=f(:,ng-3+j,:)
         fm2(:,j,:)=f(:,ng-2+j,:)
         fm1(:,j,:)=f(:,ng-1+j,:)
         f0(:,j,:)=f(:,ng+j,:)
         fp1(:,j,:)=f(:,ng+1+j,:)
         fp2(:,j,:)=f(:,ng+2+j,:)
       endforall
    elseif (dir.eq.6) then
!..recon in z, but leave x y untrucated 
!  used for first step in getting ej (next step is recon in x)
!  or for second step in getting ei (previous step is recon in y)
      forall(k=1:lnz)
         fm3(:,:,k)=f(:,:,ng-3+k)
         fm2(:,:,k)=f(:,:,ng-2+k)
         fm1(:,:,k)=f(:,:,ng-1+k)
         f0(:,:,k)=f(:,:,ng+k)
         fp1(:,:,k)=f(:,:,ng+1+k)
         fp2(:,:,k)=f(:,:,ng+2+k)
      endforall
   endif
!..for left component
   h0=(2.*fm3-7.*fm2+11.*fm1)/6.
   h1=(-fm2+5.*fm1+2.*f0)/6.
   h2=(2.*fm1+5.*f0-fp1)/6.

   s0=13./12.*(fm3-2.*fm2+fm1)**2+.25*(fm3-4.*fm2+3.*fm1)**2
   s1=13./12.*(fm2-2.*fm1+f0)**2+.25*(fm2-f0)**2
   s2=13./12.*(fm1-2.*f0+fp1)**2+.25*(3.*fm1-4.*f0+fp1)**2

   s0=.1/(eps+s0)**2
   s1=.6/(eps+s1)**2
   s2=.3/(eps+s2)**2

   dum1=s0+s1+s2
   w0=s0/dum1
   w1=s1/dum1
   w2=s2/dum1

!   fleft=w0*h0+w1*h1+w2*h2
   dum1=w0*h0+w1*h1+w2*h2

!..for right component
   h0=(2.*fp2-7.*fp1+11.*f0)/6.
   h1=(-fp1+5.*f0+2.*fm1)/6.
   h2=(2.*f0+5.*fm1-fm2)/6.

   s0=13./12.*(fp2-2.*fp1+f0)**2+.25*(fp2-4.*fp1+3.*f0)**2
   s1=13./12.*(fp1-2.*f0+fm1)**2+.25*(fp1-fm1)**2
   s2=13./12.*(f0-2.*fm1+fm2)**2+.25*(3.*f0-4.*fm1+fm2)**2

   s0=.1/(eps+s0)**2
   s1=.6/(eps+s1)**2
   s2=.3/(eps+s2)**2

   dum2=s0+s1+s2
   w0=s0/dum2
   w1=s1/dum2
   w2=s2/dum2

!   fright=w0*h0+w1*h1+w2*h2
   dum2=w0*h0+w1*h1+w2*h2

!..now add MP-limiting, refer to Balsara & Shu (2000) for detail
!..be extreme careful as we reuse arraies in calculation
!..for left component
   h0=fm3-2.*fm2+fm1
   h1=fm2-2.*fm1+f0
   h2=fm1-2.*f0+fp1

   s2=4.*h1-h2
   s0=4.*h2-h1
   s1=minmod(s2,s0,lnx,lny,lnz)
   s2=minmod(h1,h2,lnx,lny,lnz)
   s0=minmod(s1,s2,lnx,lny,lnz)

   w2=4.*h0-h1
   w0=4.*h1-h0
   w1=minmod(w2,w0,lnx,lny,lnz)
   w2=minmod(h0,h1,lnx,lny,lnz)
   w0=minmod(w1,w2,lnx,lny,lnz)

   s1=fm1+alpha*(fm1-fm2)
   s2=.5*(fm1+f0-s0)
   s0=.5*(3.*fm1-fm2)+beta*w0/3.

   w1=min(fm1,f0,s2)
   w2=min(fm1,s1,s0)
   w0=max(w1,w2)

   h1=max(fm1,f0,s2)
   h2=max(fm1,s1,s0)
   h0=min(h1,h2)

   fleft=median(dum1,w0,h0,lnx,lny,lnz)

!..for right component
   h0=fp2-2.*fp1+f0
   h1=fp1-2.*f0+fm1
   h2=f0-2.*fm1+fm2

   s2=4.*h1-h2
   s0=4.*h2-h1
   s1=minmod(s2,s0,lnx,lny,lnz)
   s2=minmod(h1,h2,lnx,lny,lnz)
   s0=minmod(s1,s2,lnx,lny,lnz)

   w2=4.*h0-h1
   w0=4.*h1-h0
   w1=minmod(w2,w0,lnx,lny,lnz)
   w2=minmod(h0,h1,lnx,lny,lnz)
   w0=minmod(w1,w2,lnx,lny,lnz)

   s1=f0+alpha*(f0-fp1)
   s2=.5*(f0+fm1-s0)
   s0=.5*(3.*f0-fp1)+beta*w0/3.

   w1=min(f0,fm1,s2)
   w2=min(f0,s1,s0)
   w0=max(w1,w2)

   h1=max(f0,fm1,s2)
   h2=max(f0,s1,s0)
   h0=min(h1,h2)

   fright=median(dum2,w0,h0,lnx,lny,lnz)

   deallocate(fm3,fm2,fm1,f0,fp1,fp2,h0,h1,h2,s0,s1,s2,w0,w1,w2)
   deallocate(dum1,dum2)

   end subroutine MPWENO5
!!***********************************
!   function minmod(a,b,lx,ly,lz)
!   implicit none
!   real*8, allocatable,dimension(:,:,:) :: a,b,minmod,dum
!   integer :: lx,ly,lz
!
!   allocate(minmod(lx,ly,lz),dum(lx,ly,lz))
!   dum=1.
!   minmod=.5*(sign(dum,a)+sign(dum,b))*min(abs(a),abs(b))
!   deallocate (dum)
!
!   end function minmod
!!***********************************
!   function median(a,b,c,lx,ly,lz)
!   implicit none
!   real*8, allocatable,dimension(:,:,:) :: a,b,c,median
!   integer :: lx,ly,lz
!
!   allocate(median(lx,ly,lz))
!   median=a+minmod(b-a,c-a,lx,ly,lz)
!
!   end function median
!!***********************************

end module reconstruct
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module get
   use pars
   use fields
   use reconstruct
   implicit none

contains
   subroutine get_dt
   implicit none
   real*8 :: gtau ! global tau
   real*8, dimension(:,:,:), allocatable :: valf,vs,vcfl,dt
   allocate(valf(nx,ny,nz),vs(nx,ny,nz),vcfl(nx,ny,nz),dt(nx,ny,nz))

   v=sqrt(vx**2+vy**2+vz**2)
   btot=sqrt(bx**2+by**2+bz**2)
   valf=btot/sqrt(rho)
   vs=sqrt(gam*p/rho)
   vcfl=v+sqrt(valf**2+vs**2)
!......more compact way
!   real*8, dimension(:,:,:), allocatable :: vcfl,dt 
!   allocate(vcfl(nx,ny,nz),dt(nx,ny,nz))
!   v=sqrt(vx**2+vy**2+vz**2)
!   vcfl=v+sqrt((bx**2+by**2+bz**2+gam*p)/rho) ! more compact way
!..................
   dt=cfl/(vcfl/dxblock+vcfl/dyblock+vcfl/dzblock)
!   dt=cfl/(vcfl/dx+vcfl/dy+vcfl/dz)
   tau=minval(dt)
!   print*,tau

   call MPI_ALLREDUCE(tau,gtau,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
   tau=gtau
!   if(myrank.eq.0) print*,tau,'***************'

   deallocate(valf,vs,vcfl,dt)
!   deallocate(vcfl,dt)

   end subroutine get_dt
!************************************
   subroutine get_conserve
   implicit none

   rho0=rho
   rhovx0=rho*vx
   rhovy0=rho*vy
   rhovz0=rho*vz
   eng0=.5*rho*(vx**2+vy**2+vz**2)+p/(gam-1.)

   end subroutine get_conserve
!************************************
!..should be easy to write in a rotate symmetric way via a do loop like
!..get_eijk_LFM
   subroutine get_eijk_Yee(lrho,lvx,lvy,lvz,lp,lbx,lby,lbz)

   implicit none
   real*8, dimension(:,:,:), allocatable :: lrho,lvx,lvy,lvz,lp,lbx,lby,lbz
   real*8, dimension(:,:,:), allocatable :: expp,exnp,expn,exnn, &
                        eypp,eynp,eypn,eynn,ezpp,eznp,ezpn,eznn
   real*8, dimension(:,:,:), allocatable :: qepp,qepn,qenp,qenn, &
                        epp,epn,enp,enn,dum1,dum2
   integer :: i,j,k

   allocate(expp(nx,ny,nz),exnp(nx,ny,nz),expn(nx,ny,nz),exnn(nx,ny,nz), &
            eypp(nx,ny,nz),eynp(nx,ny,nz),eypn(nx,ny,nz),eynn(nx,ny,nz), &
            ezpp(nx,ny,nz),eznp(nx,ny,nz),ezpn(nx,ny,nz),eznn(nx,ny,nz))

!..first split electric field
   bsq=lbx**2+lby**2+lbz**2
   lam=lrho/(2.*lp+bsq)
   call get_v_center(lam,lvx,lvy,lvz)
   
   ezpp=lbx*vx0p*vy1p-lby*vx1p*vy0p
   eznp=lbx*vx0n*vy1p-lby*vx1n*vy0p
   ezpn=lbx*vx0p*vy1n-lby*vx1p*vy0n
   eznn=lbx*vx0n*vy1n-lby*vx1n*vy0n

   expp=lby*vy0p*vz1p-lbz*vy1p*vz0p
   exnp=lby*vy0n*vz1p-lbz*vy1n*vz0p
   expn=lby*vy0p*vz1n-lbz*vy1p*vz0n
   exnn=lby*vy0n*vz1n-lbz*vy1n*vz0n

   eypp=lbz*vz0p*vx1p-lbx*vz1p*vx0p
   eynp=lbz*vz0n*vx1p-lbx*vz1n*vx0p
   eypn=lbz*vz0p*vx1n-lbx*vz1p*vx0n
   eynn=lbz*vz0n*vx1n-lbx*vz1n*vx0n

!..now get ei,j,k through reconstruction

!..first in x direction
   allocate(qepp(nx,ny0+1,nz),qepn(nx,ny0+1,nz),qenp(nx,ny0+1,nz), &
            qenn(nx,ny0+1,nz),dum1(nx,ny0+1,nz))
   allocate(epp(nx,ny0+1,nz0+1),epn(nx,ny0+1,nz0+1),enp(nx,ny0+1,nz0+1), &
            enn(nx,ny0+1,nz0+1),dum2(nx,ny0+1,nz0+1))

   call reconstruct_3d(qepp,dum1,expp,nx,ny0+1,nz,5)
   call reconstruct_3d(qepn,dum1,expn,nx,ny0+1,nz,5)
   call reconstruct_3d(dum1,qenn,exnn,nx,ny0+1,nz,5)
   call reconstruct_3d(dum1,qenp,exnp,nx,ny0+1,nz,5)

   call reconstruct_3d(epp,dum2,qepp,nx,ny0+1,nz0+1,6)
   call reconstruct_3d(enp,dum1,qenp,nx,ny0+1,nz0+1,6)
   call reconstruct_3d(dum2,enn,qenn,nx,ny0+1,nz0+1,6)
   call reconstruct_3d(dum2,epn,qepn,nx,ny0+1,nz0+1,6)
   
   forall(i=1:nx0) ei(i,:,:)=epp(ng+i,:,:)+enp(ng+i,:,:)+epn(ng+i,:,:)+enn(ng+i,:,:)
      
   deallocate(qepp,qepn,qenp,qenn,epp,epn,enp,enn,dum1,dum2)

!..in y direction
   allocate(qepp(nx,ny,nz0+1),qepn(nx,ny,nz0+1),qenp(nx,ny,nz0+1), &
            qenn(nx,ny,nz0+1),dum1(nx,ny,nz0+1))
   allocate(epp(nx0+1,ny,nz0+1),epn(nx0+1,ny,nz0+1),enp(nx0+1,ny,nz0+1), &
            enn(nx0+1,ny,nz0+1),dum2(nx0+1,ny,nz0+1))

   call reconstruct_3d(qepp,dum1,eypp,nx,ny,nz0+1,6)
   call reconstruct_3d(qepn,dum1,eypn,nx,ny,nz0+1,6)
   call reconstruct_3d(dum1,qenn,eynn,nx,ny,nz0+1,6)
   call reconstruct_3d(dum1,qenp,eynp,nx,ny,nz0+1,6)

   call reconstruct_3d(epp,dum2,qepp,nx0+1,ny,nz0+1,4)
   call reconstruct_3d(enp,dum1,qenp,nx0+1,ny,nz0+1,4)
   call reconstruct_3d(dum2,enn,qenn,nx0+1,ny,nz0+1,4)
   call reconstruct_3d(dum2,epn,qepn,nx0+1,ny,nz0+1,4)

   forall(j=1:ny0) ej(:,j,:)=epp(:,ng+j,:)+enp(:,ng+j,:)+epn(:,ng+j,:)+enn(:,ng+j,:)
   
   deallocate(qepp,qepn,qenp,qenn,epp,epn,enp,enn,dum1,dum2)

!..in z direction
   allocate(qepp(nx0+1,ny,nz),qepn(nx0+1,ny,nz),qenp(nx0+1,ny,nz), &
            qenn(nx0+1,ny,nz),dum1(nx0+1,ny,nz))
   allocate(epp(nx0+1,ny0+1,nz),epn(nx0+1,ny0+1,nz),enp(nx0+1,ny0+1,nz), &
            enn(nx0+1,ny0+1,nz),dum2(nx0+1,ny0+1,nz))

   call reconstruct_3d(qepp,dum1,ezpp,nx0+1,ny,nz,4)
   call reconstruct_3d(qepn,dum1,ezpn,nx0+1,ny,nz,4)
   call reconstruct_3d(dum1,qenn,eznn,nx0+1,ny,nz,4)
   call reconstruct_3d(dum1,qenp,eznp,nx0+1,ny,nz,4)

   call reconstruct_3d(epp,dum2,qepp,nx0+1,ny0+1,nz,5)
   call reconstruct_3d(enp,dum1,qenp,nx0+1,ny0+1,nz,5)
   call reconstruct_3d(dum2,enn,qenn,nx0+1,ny0+1,nz,5)
   call reconstruct_3d(dum2,epn,qepn,nx0+1,ny0+1,nz,5)

   forall(k=1:nz0) ek(:,:,k)=epp(:,:,ng+k)+enp(:,:,ng+k)+epn(:,:,ng+k)+enn(:,:,ng+k)

   deallocate(qepp,qepn,qenp,qenn,epp,epn,enp,enn,dum1,dum2)

   deallocate(expp,exnp,expn,exnn,eypp,eynp,eypn,eynn,ezpp,eznp,ezpn,eznn)

   end subroutine get_eijk_Yee
!************************************
   subroutine get_eijk_LFM(lrho,lvx,lvy,lvz,lp,lbx,lby,lbz,lbi,lbj,lbk)
   
   implicit none
   real*8, dimension(:,:,:), allocatable :: lrho,lvx,lvy,lvz,lp,lbx,lby,lbz, &
                                            lbi,lbj,lbk
   real*8, dimension(:,:,:), allocatable :: vf
   real*8, dimension(:,:,:), allocatable :: v1,v2,b1,b2
   real*8, dimension(:,:,:), allocatable :: v1av,v2av,sig1,sig2,dum
   real*8, dimension(:,:,:), allocatable :: b1l,b1r,b2l,b2r,b1u,b2u,eta
   integer :: i,j,k,ii,jj,kk,d1,d2,kkk

   allocate(vf(nx,ny,nz))
   vf=sqrt(lvx**2+lvy**2+lvz**2)+sqrt((lbx**2+lby**2+lbz**2)/lrho)

!..first calculate ei
!..ei has size (nx0,ny0+1,nz0+1)
!..in this case, v1av=vy_avg, v2av=vz_avy, b1=bj, b2=bk

   allocate(v1(nx,ny,nz),v2(nx,ny,nz))

   do kkk=1,3 ! for x.y and z direction

   ii=nx0+1
   jj=ny0+1
   kk=nz0+1

   if (kkk.eq.1) then
      ii=nx0
      v1=lvy
      v2=lvz
      allocate(b1(nx,ny+1,nz),b2(nx,ny,nz+1))
      b1=lbj
      b2=lbk
      d1=3
      d2=2
   elseif (kkk.eq.2) then
      jj=ny0
      v1=lvz
      v2=lvx
      allocate(b1(nx,ny,nz+1),b2(nx+1,ny,nz))
      b1=lbk
      b2=lbi
      d1=1
      d2=3
   elseif (kkk.eq.3) then
      kk=nz0
      v1=lvx
      v2=lvy
      allocate(b1(nx+1,ny,nz),b2(nx,ny+1,nz))
      b1=lbi
      b2=lbj
      d1=2
      d2=1
   endif

   allocate(v1av(ii,jj,kk),v2av(ii,jj,kk),sig1(ii,jj,kk),sig2(ii,jj,kk), &
            b1l(ii,jj,kk),b1r(ii,jj,kk),b2l(ii,jj,kk),b2r(ii,jj,kk), &
            b1u(ii,jj,kk),b2u(ii,jj,kk),eta(ii,jj,kk),dum(ii,jj,kk))

   call get_vav_center2corner(v1,v1av,kkk)
   call get_vav_center2corner(v2,v2av,kkk)

   dum=1.
   sig1=sign(dum,v1av)
   sig2=sign(dum,v2av)

   call reconstruct_3d(b1l,b1r,b1,ii,jj,kk,d1)
   call reconstruct_3d(b2l,b2r,b2,ii,jj,kk,d2)

   b1u=(1.+sig2)/2.*b1l+(1.-sig2)/2.*b1r
   b2u=(1.+sig1)/2.*b2l+(1.-sig1)/2.*b2r

   dum=-v1av*b2u+v2av*b1u

!..apply Alfven resistivity
   if (kkk.eq.1) then
      forall(i=1:ii,j=1:jj,k=1:kk) eta(i,j,k)= .25* &
            (vf(i+ng,j+ng,k+ng)+vf(i+ng,j+ng-1,k+ng)+ &
             vf(i+ng,j+ng,k+ng-1)+vf(i+ng,j+ng-1,k+ng-1))
   elseif(kkk.eq.2) then
      forall(i=1:ii,j=1:jj,k=1:kk) eta(i,j,k)= .25* &
            (vf(i+ng,j+ng,k+ng)+vf(i+ng-1,j+ng,k+ng)+ &
             vf(i+ng,j+ng,k+ng-1)+vf(i+ng-1,j+ng,k+ng-1))
   elseif(kkk.eq.3) then
      forall(i=1:ii,j=1:jj,k=1:kk) eta(i,j,k)= .25* &
            (vf(i+ng,j+ng,k+ng)+vf(i+ng-1,j+ng,k+ng)+ &
             vf(i+ng,j+ng-1,k+ng)+vf(i+ng-1,j+ng-1,k+ng))
   endif
   dum=dum+eta*(b1l+b2r-b1r-b2l)

   if (kkk.eq.1) ei=dum
   if (kkk.eq.2) ej=dum
   if (kkk.eq.3) ek=dum

   deallocate(v1av,v2av,sig1,sig2,b1l,b1r,b2l,b2r,b1u,b2u,eta,dum)
   deallocate(b1,b2)

   enddo ! enddo for kkk

   deallocate(v1,v2)
   deallocate(vf)

   end subroutine get_eijk_LFM
!!************************************
   subroutine get_vav_center2corner(lv,vav,dir)

   implicit none
   real*8, dimension(:,:,:), allocatable :: lv,vav
   real*8, dimension(:,:,:), allocatable :: lvl,lvr,vi
   real*8, dimension(:,:,:), allocatable :: v00,v10,v01,v11
   real*8, dimension(:,:,:), allocatable :: v00_0,v00_1,v10_0,v10_1, &
                                            v01_0,v01_1,v11_0,v11_1, &
                                            dif00,dif01,dif10,dif11,dum
   integer :: dir
   integer :: i,j,k,ii,jj,kk,snx,sny,snz,d1,d2

   ii=nx0+1
   jj=ny0+1
   kk=nz0+1

   if (dir.eq.1) then ! for ei calculation
      ii=nx0
      d1=2
      d2=3
   elseif (dir.eq.2) then ! for ej calculation
      jj=ny0
      d1=3
      d2=1
   elseif (dir.eq.3) then ! for ek calculation
      kk=nz0
      d1=1
      d2=2
   endif

   allocate(v00(ii,jj,kk),v10(ii,jj,kk),v01(ii,jj,kk),v11(ii,jj,kk))
   allocate(v00_0(ii,jj,kk),v00_1(ii,jj,kk),v10_0(ii,jj,kk),v10_1(ii,jj,kk), &
            v01_0(ii,jj,kk),v01_1(ii,jj,kk),v11_0(ii,jj,kk),v11_1(ii,jj,kk) )
   allocate(dif00(ii,jj,kk),dif10(ii,jj,kk),dif01(ii,jj,kk),dif11(ii,jj,kk),dum(ii,jj,kk))

!..reconstruct
   if (reconop.eq.1) then
      allocate(vi(ii,jj,kk))
!.....first get interped velocity vi
      call interp2_pdm(lv,vi,dir)
      if(dir.eq.1) then
      call PDM2(v00_0,v10_0,lv,vi,ii,jj,kk,0,0,-1,2)
      call PDM2(v00_1,v01_1,lv,vi,ii,jj,kk,0,-1,0,3)
      call PDM2(v10_1,v11_1,lv,vi,ii,jj,kk,0,0,0,3)
      call PDM2(v01_0,v11_0,lv,vi,ii,jj,kk,0,0,0,2)
      elseif(dir.eq.2) then
      call PDM2(v00_0,v10_0,lv,vi,ii,jj,kk,-1,0,0,3)
      call PDM2(v00_1,v01_1,lv,vi,ii,jj,kk,0,0,-1,1)
      call PDM2(v10_1,v11_1,lv,vi,ii,jj,kk,0,0,0,1)
      call PDM2(v01_0,v11_0,lv,vi,ii,jj,kk,0,0,0,3)
      elseif(dir.eq.3) then
      call PDM2(v00_0,v10_0,lv,vi,ii,jj,kk,0,-1,0,1)
      call PDM2(v00_1,v01_1,lv,vi,ii,jj,kk,-1,0,0,2)
      call PDM2(v10_1,v11_1,lv,vi,ii,jj,kk,0,0,0,2)
      call PDM2(v01_0,v11_0,lv,vi,ii,jj,kk,0,0,0,1)
      endif
      deallocate(vi)

   elseif(reconop.eq.2) then
      allocate(vi(ii,jj,kk))
      call interp2_tvd(lv,vi,dir)
      if(dir.eq.1) then
      call TVD(v00_0,v10_0,lv,vi,ii,jj,kk,0,0,-1,2)
      call TVD(v00_1,v01_1,lv,vi,ii,jj,kk,0,-1,0,3)
      call TVD(v10_1,v11_1,lv,vi,ii,jj,kk,0,0,0,3)
      call TVD(v01_0,v11_0,lv,vi,ii,jj,kk,0,0,0,2)
      elseif(dir.eq.2) then
      call TVD(v00_0,v10_0,lv,vi,ii,jj,kk,-1,0,0,3)
      call TVD(v00_1,v01_1,lv,vi,ii,jj,kk,0,0,-1,1)
      call TVD(v10_1,v11_1,lv,vi,ii,jj,kk,0,0,0,1)
      call TVD(v01_0,v11_0,lv,vi,ii,jj,kk,0,0,0,3)
      elseif(dir.eq.3) then
      call TVD(v00_0,v10_0,lv,vi,ii,jj,kk,0,-1,0,1)
      call TVD(v00_1,v01_1,lv,vi,ii,jj,kk,-1,0,0,2)
      call TVD(v10_1,v11_1,lv,vi,ii,jj,kk,0,0,0,2)
      call TVD(v01_0,v11_0,lv,vi,ii,jj,kk,0,0,0,1)
      endif
      deallocate(vi)

   elseif(reconop.eq.3) then
      allocate(vi(nx,ny,nz))
      if(dir.eq.1) then
!.....first reconstruct in y
      allocate(lvl(nx,jj,nz),lvr(nx,jj,nz))
      call WENO5(lvl,lvr,lv,nx,jj,nz,5)
      forall(j=1:jj) vi(:,ng+j,:)=lvl(:,j,:)
      call WENO5(v00_1,v01_1,vi,ii,jj,kk,3)
      forall(j=1:jj) vi(:,ng+j,:)=lvr(:,j,:)
      call WENO5(v10_1,v11_1,vi,ii,jj,kk,3)
      deallocate(lvl,lvr)
!.....then reconstruct in z
      allocate(lvl(nx,ny,kk),lvr(nx,ny,kk))
      call WENO5(lvl,lvr,lv,nx,ny,kk,6)
      forall(k=1:kk) vi(:,:,ng+k)=lvl(:,:,k)
      call WENO5(v00_0,v01_0,vi,ii,jj,kk,2)
      forall(k=1:kk) vi(:,:,ng+k)=lvr(:,:,k)
      call WENO5(v10_0,v11_0,vi,ii,jj,kk,2)
      deallocate(lvl,lvr)

      elseif(dir.eq.2) then
!.....first reconstruct in z
      allocate(lvl(nx,ny,kk),lvr(nx,ny,kk))
      call WENO5(lvl,lvr,lv,nx,ny,kk,6)
      forall(k=1:kk) vi(:,:,ng+k)=lvl(:,:,k)
      call WENO5(v00_1,v01_1,vi,ii,jj,kk,1)
      forall(k=1:kk) vi(:,:,ng+k)=lvr(:,:,k)
      call WENO5(v10_1,v11_1,vi,ii,jj,kk,1)
      deallocate(lvl,lvr)
!.....then reconstruct in x
      allocate(lvl(ii,ny,nz),lvr(ii,ny,nz))
      call WENO5(lvl,lvr,lv,ii,ny,nz,4)
      forall(i=1:ii) vi(ng+i,:,:)=lvl(i,:,:)
      call WENO5(v00_0,v01_0,vi,ii,jj,kk,3)
      forall(i=1:ii) vi(ng+i,:,:)=lvr(i,:,:)
      call WENO5(v10_0,v11_0,vi,ii,jj,kk,3)
      deallocate(lvl,lvr)

      elseif(dir.eq.3) then
!.....first reconstruct in x
      allocate(lvl(ii,ny,nz),lvr(ii,ny,nz))
      call WENO5(lvl,lvr,lv,ii,ny,nz,4)
      forall(i=1:ii) vi(ng+i,:,:)=lvl(i,:,:)
      call WENO5(v00_1,v01_1,vi,ii,jj,kk,2)
      forall(i=1:ii) vi(ng+i,:,:)=lvr(i,:,:)
      call WENO5(v10_1,v11_1,vi,ii,jj,kk,2)
      deallocate(lvl,lvr)
!.....then reconstruct in y
      allocate(lvl(nx,jj,nz),lvr(nx,jj,nz))
      call WENO5(lvl,lvr,lv,nx,jj,nz,5)
      forall(j=1:jj) vi(:,ng+j,:)=lvl(:,j,:)
      call WENO5(v00_0,v01_0,vi,ii,jj,kk,1)
      forall(j=1:jj) vi(:,ng+j,:)=lvr(:,j,:)
      call WENO5(v10_0,v11_0,vi,ii,jj,kk,1)
      deallocate(lvl,lvr)
      endif
      deallocate(vi)

   elseif(reconop.eq.4) then
      allocate(vi(nx,ny,nz))
      if(dir.eq.1) then
!.....first reconstruct in y
      allocate(lvl(nx,jj,nz),lvr(nx,jj,nz))
      call MPWENO5(lvl,lvr,lv,nx,jj,nz,5)
      forall(j=1:jj) vi(:,ng+j,:)=lvl(:,j,:)
      call MPWENO5(v00_1,v01_1,vi,ii,jj,kk,3)
      forall(j=1:jj) vi(:,ng+j,:)=lvr(:,j,:)
      call MPWENO5(v10_1,v11_1,vi,ii,jj,kk,3)
      deallocate(lvl,lvr)
!.....then reconstruct in z
      allocate(lvl(nx,ny,kk),lvr(nx,ny,kk))
      call MPWENO5(lvl,lvr,lv,nx,ny,kk,6)
      forall(k=1:kk) vi(:,:,ng+k)=lvl(:,:,k)
      call MPWENO5(v00_0,v01_0,vi,ii,jj,kk,2)
      forall(k=1:kk) vi(:,:,ng+k)=lvr(:,:,k)
      call MPWENO5(v10_0,v11_0,vi,ii,jj,kk,2)
      deallocate(lvl,lvr)

      elseif(dir.eq.2) then
!.....first reconstruct in z
      allocate(lvl(nx,ny,kk),lvr(nx,ny,kk))
      call MPWENO5(lvl,lvr,lv,nx,ny,kk,6)
      forall(k=1:kk) vi(:,:,ng+k)=lvl(:,:,k)
      call MPWENO5(v00_1,v01_1,vi,ii,jj,kk,1)
      forall(k=1:kk) vi(:,:,ng+k)=lvr(:,:,k)
      call MPWENO5(v10_1,v11_1,vi,ii,jj,kk,1)
      deallocate(lvl,lvr)
!.....then reconstruct in x
      allocate(lvl(ii,ny,nz),lvr(ii,ny,nz))
      call MPWENO5(lvl,lvr,lv,ii,ny,nz,4)
      forall(i=1:ii) vi(ng+i,:,:)=lvl(i,:,:)
      call MPWENO5(v00_0,v01_0,vi,ii,jj,kk,3)
      forall(i=1:ii) vi(ng+i,:,:)=lvr(i,:,:)
      call MPWENO5(v10_0,v11_0,vi,ii,jj,kk,3)
      deallocate(lvl,lvr)

      elseif(dir.eq.3) then
!.....first reconstruct in x
      allocate(lvl(ii,ny,nz),lvr(ii,ny,nz))
      call MPWENO5(lvl,lvr,lv,ii,ny,nz,4)
      forall(i=1:ii) vi(ng+i,:,:)=lvl(i,:,:)
      call MPWENO5(v00_1,v01_1,vi,ii,jj,kk,2)
      forall(i=1:ii) vi(ng+i,:,:)=lvr(i,:,:)
      call MPWENO5(v10_1,v11_1,vi,ii,jj,kk,2)
      deallocate(lvl,lvr)
!.....then reconstruct in y
      allocate(lvl(nx,jj,nz),lvr(nx,jj,nz))
      call MPWENO5(lvl,lvr,lv,nx,jj,nz,5)
      forall(j=1:jj) vi(:,ng+j,:)=lvl(:,j,:)
      call MPWENO5(v00_0,v01_0,vi,ii,jj,kk,1)
      forall(j=1:jj) vi(:,ng+j,:)=lvr(:,j,:)
      call MPWENO5(v10_0,v11_0,vi,ii,jj,kk,1)
      deallocate(lvl,lvr)
      endif
      deallocate(vi)

   else
      print*,'***N/A***'
   endif

   dif00=0.
   dif01=0.
   dif10=0.
   dif11=0.

   if (dir.eq.1) then
      forall(i=1:ii,j=1:jj,k=1:kk) 
         dif00(i,j,k)=abs(lv(i+ng,j+ng-2,k+ng-1)-lv(i+ng,j+ng-1,k+ng-1))- &
                      abs(lv(i+ng,j+ng-1,k+ng-2)-lv(i+ng,j+ng-1,k+ng-1))
         dif10(i,j,k)=abs(lv(i+ng,j+ng+1,k+ng-1)-lv(i+ng,j+ng,k+ng-1))-   &
                      abs(lv(i+ng,j+ng,k+ng-2)  -lv(i+ng,j+ng,k+ng-1))
         dif01(i,j,k)=abs(lv(i+ng,j+ng-2,k+ng)  -lv(i+ng,j+ng-1,k+ng))-   &
                      abs(lv(i+ng,j+ng-1,k+ng+1)-lv(i+ng,j+ng-1,k+ng))
         dif11(i,j,k)=abs(lv(i+ng,j+ng+1,k+ng)  -lv(i+ng,j+ng,k+ng))-     &
                      abs(lv(i+ng,j+ng,k+ng+1)  -lv(i+ng,j+ng,k+ng))
      endforall
   elseif (dir.eq.2) then
      forall(i=1:ii,j=1:jj,k=1:kk) 
         dif00(i,j,k)=abs(lv(i+ng-1,j+ng,k+ng-2)-lv(i+ng-1,j+ng,k+ng-1))- &
                      abs(lv(i+ng-2,j+ng,k+ng-1)-lv(i+ng-1,j+ng,k+ng-1))
         dif10(i,j,k)=abs(lv(i+ng-1,j+ng,k+ng+1)-lv(i+ng-1,j+ng,k+ng))-   &
                      abs(lv(i+ng-2,j+ng,k+ng)  -lv(i+ng-1,j+ng,k+ng))
         dif01(i,j,k)=abs(lv(i+ng,j+ng,k+ng-2)  -lv(i+ng,j+ng,k+ng-1))-   &
                      abs(lv(i+ng+1,j+ng,k+ng-1)-lv(i+ng,j+ng,k+ng-1))
         dif11(i,j,k)=abs(lv(i+ng,j+ng,k+ng+1)  -lv(i+ng,j+ng,k+ng))-     &
                      abs(lv(i+ng+1,j+ng,k+ng)  -lv(i+ng,j+ng,k+ng))
      endforall
   elseif (dir.eq.3) then
      forall(i=1:ii,j=1:jj,k=1:kk)
         dif00(i,j,k)=abs(lv(i+ng-2,j+ng-1,k+ng)-lv(i+ng-1,j+ng-1,k+ng))- &
                      abs(lv(i+ng-1,j+ng-2,k+ng)-lv(i+ng-1,j+ng-1,k+ng))
         dif10(i,j,k)=abs(lv(i+ng+1,j+ng-1,k+ng)-lv(i+ng,j+ng-1,k+ng))-   &
                      abs(lv(i+ng,j+ng-2,k+ng)  -lv(i+ng,j+ng-1,k+ng))
         dif01(i,j,k)=abs(lv(i+ng-2,j+ng,k+ng)  -lv(i+ng-1,j+ng,k+ng))-   &
                      abs(lv(i+ng-1,j+ng+1,k+ng)-lv(i+ng-1,j+ng,k+ng))
         dif11(i,j,k)=abs(lv(i+ng+1,j+ng,k+ng)  -lv(i+ng,j+ng,k+ng))-     &
                      abs(lv(i+ng,j+ng+1,k+ng)  -lv(i+ng,j+ng,k+ng))
      endforall
   endif

   dum=1.
   dif00=sign(dum,dif00)
   dif10=sign(dum,dif10)
   dif01=sign(dum,dif01)
   dif11=sign(dum,dif11)

   v00=(1.+dif00)/2.*v00_0+(1.-dif00)/2.*v00_1
   v10=(1.+dif10)/2.*v10_0+(1.-dif10)/2.*v10_1
   v01=(1.+dif01)/2.*v01_0+(1.-dif01)/2.*v01_1
   v11=(1.+dif11)/2.*v11_0+(1.-dif11)/2.*v11_1

   vav=.25*(v00+v01+v10+v11)

   deallocate(v00,v10,v01,v11)
   deallocate(v00_0,v00_1,v10_0,v10_1,v01_0,v01_1,v11_0,v11_1, &
              dif00,dif10,dif01,dif11,dum)

   return

   end subroutine get_vav_center2corner
!************************************
   subroutine get_bxyz
!..calculate bx,by,bz from bi,bj,bk
   implicit none
   integer :: i,j,k
   
   forall(i=1:nx) bx(i,:,:)=.5*(bi(i+1,:,:)+bi(i,:,:))
   forall(j=1:ny) by(:,j,:)=.5*(bj(:,j+1,:)+bj(:,j,:))
   forall(k=1:nz) bz(:,:,k)=.5*(bk(:,:,k+1)+bk(:,:,k))
   bsq=bx**2+by**2+bz**2

   end subroutine get_bxyz
!************************************
!..again, could write in a rotate symmetic way
   subroutine get_flux_stress_old
   implicit none
   real*8, dimension(:,:,:), allocatable :: rhol,rhor,vxl,vxr,vyl,vyr, &
                                vzl,vzr,pl,pr,bxl,bxr,byl,byr,bzl,bzr,lbsq
   real*8, dimension(:,:,:), allocatable :: rp,rn,rvxp,rvxn,rvyp,rvyn, &
                                rvzp,rvzn,ep,en
   real*8, dimension(:,:,:), allocatable :: sxp,sxn,syp,syn,szp,szn
!..local energy, lambda and split velocity arraies
   real*8, dimension(:,:,:), allocatable :: leng,llam,v0p,v0n,v1p,v1n   
   integer :: ii,jj,kk

!..x direction
   ii=nx0+1
   jj=ny0
   kk=nz0

   allocate(rhol(nx0+1,ny0,nz0),rhor(nx0+1,ny0,nz0),vxl(nx0+1,ny0,nz0), &
            vxr(nx0+1,ny0,nz0),vyl(nx0+1,ny0,nz0),vyr(nx0+1,ny0,nz0), &
            vzl(nx0+1,ny0,nz0),vzr(nx0+1,ny0,nz0),pl(nx0+1,ny0,nz0), &
            pr(nx0+1,ny0,nz0),bxl(nx0+1,ny0,nz0),bxr(nx0+1,ny0,nz0), &
            byl(nx0+1,ny0,nz0),byr(nx0+1,ny0,nz0),bzl(nx0+1,ny0,nz0), &
            bzr(nx0+1,ny0,nz0),lbsq(nx0+1,ny0,nz0))
   allocate(rp(nx0+1,ny0,nz0),rn(nx0+1,ny0,nz0),rvxp(nx0+1,ny0,nz0), &
            rvxn(nx0+1,ny0,nz0),rvyp(nx0+1,ny0,nz0),rvyn(nx0+1,ny0,nz0), &
            rvzp(nx0+1,ny0,nz0),rvzn(nx0+1,ny0,nz0),ep(nx0+1,ny0,nz0), &
            en(nx0+1,ny0,nz0))
   allocate(sxp(nx0+1,ny0,nz0),sxn(nx0+1,ny0,nz0),syp(nx0+1,ny0,nz0), &
            syn(nx0+1,ny0,nz0),szp(nx0+1,ny0,nz0),szn(nx0+1,ny0,nz0))
   allocate(leng(nx0+1,ny0,nz0),llam(nx0+1,ny0,nz0),v0p(nx0+1,ny0,nz0), &
            v0n(nx0+1,ny0,nz0),v1p(nx0+1,ny0,nz0),v1n(nx0+1,ny0,nz0))

   call reconstruct_3d(rhol,rhor,rhoh,ii,jj,kk,1)
   call reconstruct_3d(vxl,vxr,vxh,ii,jj,kk,1)
   call reconstruct_3d(vyl,vyr,vyh,ii,jj,kk,1)
   call reconstruct_3d(vzl,vzr,vzh,ii,jj,kk,1)
   call reconstruct_3d(pl,pr,ph,ii,jj,kk,1)
   call reconstruct_3d(bxl,bxr,bxh,ii,jj,kk,1)
   call reconstruct_3d(byl,byr,byh,ii,jj,kk,1)
   call reconstruct_3d(bzl,bzr,bzh,ii,jj,kk,1)

!..left flux
   leng=.5*rhol*(vxl**2+vyl**2+vzl**2)+pl/(gam-1.)
   llam=rhol/(2.*pl)

   call get_v_face(llam,vxl,v0p,v0n,v1p,v1n,1)
   rp=rhol*v1p
   rvxp=rhol*vxl*v1p+pl*v0p
   rvyp=rhol*vyl*v1p
   rvzp=rhol*vzl*v1p
   ep=(leng+.5*pl)*v1p+.5*pl*vxl*v0p

   lbsq=bxl**2+byl**2+bzl**2
   llam=rhol/(2.*pl+lbsq)
   call get_v_face(llam,vxl,v0p,v0n,v1p,v1n,1)
   sxp=(.5*lbsq-bxl*bxl)*v0p
   syp=(-byl*bxl)*v0p
   szp=(-bzl*bxl)*v0p

!..right flux
   leng=.5*rhor*(vxr**2+vyr**2+vzr**2)+pr/(gam-1.)
   llam=rhor/(2.*pr)

   call get_v_face(llam,vxr,v0p,v0n,v1p,v1n,1)
   rn=rhor*v1n
   rvxn=rhor*vxr*v1n+pr*v0n
   rvyn=rhor*vyr*v1n
   rvzn=rhor*vzr*v1n
   en=(leng+.5*pr)*v1n+.5*pr*vxr*v0n

   lbsq=bxr**2+byr**2+bzr**2
   llam=rhor/(2.*pr+lbsq)
   call get_v_face(llam,vxr,v0p,v0n,v1p,v1n,1)
   sxn=(.5*lbsq-bxr*bxr)*v0n
   syn=(-byr*bxr)*v0n
   szn=(-bzr*bxr)*v0n

   flux_rhox=rp+rn
   flux_rvxx=rvxp+rvxn
   flux_rvyx=rvyp+rvyn
   flux_rvzx=rvzp+rvzn
   flux_engx=ep+en

   stress_xx=sxp+sxn
   stress_yx=syp+syn
   stress_zx=szp+szn

   deallocate(rhol,rhor,vxl,vxr,vyl,vyr,vzl,vzr,pl,pr,bxl,bxr,byl,byr,bzl, &
            bzr,lbsq)
   deallocate(rp,rn,rvxp,rvxn,rvyp,rvyn,rvzp,rvzn,ep,en)
   deallocate(sxp,sxn,syp,syn,szp,szn)
   deallocate(leng,llam,v0p,v0n,v1p,v1n)

!..y direction
   ii=nx0
   jj=ny0+1
   kk=nz0
   allocate(rhol(ii,jj,kk),rhor(ii,jj,kk),vxl(ii,jj,kk), &
            vxr(ii,jj,kk),vyl(ii,jj,kk),vyr(ii,jj,kk), &
            vzl(ii,jj,kk),vzr(ii,jj,kk),pl(ii,jj,kk), &
            pr(ii,jj,kk),bxl(ii,jj,kk),bxr(ii,jj,kk), &
            byl(ii,jj,kk),byr(ii,jj,kk),bzl(ii,jj,kk), &
            bzr(ii,jj,kk),lbsq(ii,jj,kk))
   allocate(rp(ii,jj,kk),rn(ii,jj,kk),rvxp(ii,jj,kk), &
            rvxn(ii,jj,kk),rvyp(ii,jj,kk),rvyn(ii,jj,kk), &
            rvzp(ii,jj,kk),rvzn(ii,jj,kk),ep(ii,jj,kk), &
            en(ii,jj,kk))
   allocate(sxp(ii,jj,kk),sxn(ii,jj,kk),syp(ii,jj,kk), &
            syn(ii,jj,kk),szp(ii,jj,kk),szn(ii,jj,kk))
   allocate(leng(ii,jj,kk),llam(ii,jj,kk),v0p(ii,jj,kk), &
            v0n(ii,jj,kk),v1p(ii,jj,kk),v1n(ii,jj,kk))

   call reconstruct_3d(rhol,rhor,rhoh,ii,jj,kk,2)
   call reconstruct_3d(vxl,vxr,vxh,ii,jj,kk,2)
   call reconstruct_3d(vyl,vyr,vyh,ii,jj,kk,2)
   call reconstruct_3d(vzl,vzr,vzh,ii,jj,kk,2)
   call reconstruct_3d(pl,pr,ph,ii,jj,kk,2)
   call reconstruct_3d(bxl,bxr,bxh,ii,jj,kk,2)
   call reconstruct_3d(byl,byr,byh,ii,jj,kk,2)
   call reconstruct_3d(bzl,bzr,bzh,ii,jj,kk,2)

!..left flux
   leng=.5*rhol*(vxl**2+vyl**2+vzl**2)+pl/(gam-1.)
   llam=rhol/(2.*pl)

   call get_v_face(llam,vyl,v0p,v0n,v1p,v1n,2)
   rp=rhol*v1p
   rvxp=rhol*vxl*v1p
   rvyp=rhol*vyl*v1p+pl*v0p
   rvzp=rhol*vzl*v1p
   ep=(leng+.5*pl)*v1p+.5*pl*vyl*v0p

   lbsq=bxl**2+byl**2+bzl**2
   llam=rhol/(2.*pl+lbsq)
   call get_v_face(llam,vyl,v0p,v0n,v1p,v1n,2)
   sxp=(-bxl*byl)*v0p
   syp=(.5*lbsq-byl*byl)*v0p
   szp=(-bzl*byl)*v0p

!..right flux
   leng=.5*rhor*(vxr**2+vyr**2+vzr**2)+pr/(gam-1.)
   llam=rhor/(2.*pr)

   call get_v_face(llam,vyr,v0p,v0n,v1p,v1n,2)
   rn=rhor*v1n
   rvxn=rhor*vxr*v1n
   rvyn=rhor*vyr*v1n+pr*v0n
   rvzn=rhor*vzr*v1n
   en=(leng+.5*pr)*v1n+.5*pr*vyr*v0n

   lbsq=bxr**2+byr**2+bzr**2
   llam=rhor/(2.*pr+lbsq)
   call get_v_face(llam,vyr,v0p,v0n,v1p,v1n,2)
   sxn=(-bxr*byr)*v0n
   syn=(.5*lbsq-byr*byr)*v0n
   szn=(-bzr*byr)*v0n

   flux_rhoy=rp+rn
   flux_rvxy=rvxp+rvxn
   flux_rvyy=rvyp+rvyn
   flux_rvzy=rvzp+rvzn
   flux_engy=ep+en

   stress_xy=sxp+sxn
   stress_yy=syp+syn
   stress_zy=szp+szn

   deallocate(rhol,rhor,vxl,vxr,vyl,vyr,vzl,vzr,pl,pr,bxl,bxr,byl,byr,bzl, &
            bzr,lbsq)
   deallocate(rp,rn,rvxp,rvxn,rvyp,rvyn,rvzp,rvzn,ep,en)
   deallocate(sxp,sxn,syp,syn,szp,szn)
   deallocate(leng,llam,v0p,v0n,v1p,v1n)

!..z direction
   ii=nx0
   jj=ny0
   kk=nz0+1
   allocate(rhol(ii,jj,kk),rhor(ii,jj,kk),vxl(ii,jj,kk), &
            vxr(ii,jj,kk),vyl(ii,jj,kk),vyr(ii,jj,kk), &
            vzl(ii,jj,kk),vzr(ii,jj,kk),pl(ii,jj,kk), &
            pr(ii,jj,kk),bxl(ii,jj,kk),bxr(ii,jj,kk), &
            byl(ii,jj,kk),byr(ii,jj,kk),bzl(ii,jj,kk), &
            bzr(ii,jj,kk),lbsq(ii,jj,kk))
   allocate(rp(ii,jj,kk),rn(ii,jj,kk),rvxp(ii,jj,kk), &
            rvxn(ii,jj,kk),rvyp(ii,jj,kk),rvyn(ii,jj,kk), &
            rvzp(ii,jj,kk),rvzn(ii,jj,kk),ep(ii,jj,kk), &
            en(ii,jj,kk))
   allocate(sxp(ii,jj,kk),sxn(ii,jj,kk),syp(ii,jj,kk), &
            syn(ii,jj,kk),szp(ii,jj,kk),szn(ii,jj,kk))
   allocate(leng(ii,jj,kk),llam(ii,jj,kk),v0p(ii,jj,kk), &
            v0n(ii,jj,kk),v1p(ii,jj,kk),v1n(ii,jj,kk))

   call reconstruct_3d(rhol,rhor,rhoh,ii,jj,kk,3)
   call reconstruct_3d(vxl,vxr,vxh,ii,jj,kk,3)
   call reconstruct_3d(vyl,vyr,vyh,ii,jj,kk,3)
   call reconstruct_3d(vzl,vzr,vzh,ii,jj,kk,3)
   call reconstruct_3d(pl,pr,ph,ii,jj,kk,3)
   call reconstruct_3d(bxl,bxr,bxh,ii,jj,kk,3)
   call reconstruct_3d(byl,byr,byh,ii,jj,kk,3)
   call reconstruct_3d(bzl,bzr,bzh,ii,jj,kk,3)

!..left flux
   leng=.5*rhol*(vxl**2+vyl**2+vzl**2)+pl/(gam-1.)
   llam=rhol/(2.*pl)

   call get_v_face(llam,vzl,v0p,v0n,v1p,v1n,3)
   rp=rhol*v1p
   rvxp=rhol*vxl*v1p
   rvyp=rhol*vyl*v1p
   rvzp=rhol*vzl*v1p+pl*v0p
   ep=(leng+.5*pl)*v1p+.5*pl*vzl*v0p

   lbsq=bxl**2+byl**2+bzl**2
   llam=rhol/(2.*pl+lbsq)
   call get_v_face(llam,vzl,v0p,v0n,v1p,v1n,3)
   sxp=(-bxl*bzl)*v0p
   syp=(-byl*bzl)*v0p
   szp=(.5*lbsq-bzl*bzl)*v0p

!..right flux
   leng=.5*rhor*(vxr**2+vyr**2+vzr**2)+pr/(gam-1.)
   llam=rhor/(2.*pr)

   call get_v_face(llam,vzr,v0p,v0n,v1p,v1n,3)
   rn=rhor*v1n
   rvxn=rhor*vxr*v1n
   rvyn=rhor*vyr*v1n
   rvzn=rhor*vzr*v1n+pr*v0n
   en=(leng+.5*pr)*v1n+.5*pr*vzr*v0n

   lbsq=bxr**2+byr**2+bzr**2
   llam=rhor/(2.*pr+lbsq)
   call get_v_face(llam,vzr,v0p,v0n,v1p,v1n,3)
   sxn=(-bxr*bzr)*v0n
   syn=(-byr*bzr)*v0n
   szn=(.5*lbsq-bzr*bzr)*v0n

   flux_rhoz=rp+rn
   flux_rvxz=rvxp+rvxn
   flux_rvyz=rvyp+rvyn
   flux_rvzz=rvzp+rvzn
   flux_engz=ep+en

   stress_xz=sxp+sxn
   stress_yz=syp+syn
   stress_zz=szp+szn

   deallocate(rhol,rhor,vxl,vxr,vyl,vyr,vzl,vzr,pl,pr,bxl,bxr,byl,byr,bzl, &
            bzr,lbsq)
   deallocate(rp,rn,rvxp,rvxn,rvyp,rvyn,rvzp,rvzn,ep,en)
   deallocate(sxp,sxn,syp,syn,szp,szn)
   deallocate(leng,llam,v0p,v0n,v1p,v1n)

   end subroutine get_flux_stress_old
!************************************
   subroutine get_flux_stress(lrho,lvx,lvy,lvz,lp,lbx,lby,lbz)
   implicit none
   real*8, dimension(:,:,:), allocatable :: lrho,lvx,lvy,lvz,lp,lbx,lby,lbz
   real*8, dimension(:,:,:), allocatable :: rhol,rhor,vxl,vxr,vyl,vyr, &
                                vzl,vzr,pl,pr,bxl,bxr,byl,byr,bzl,bzr,lbsq
   real*8, dimension(:,:,:), allocatable :: rp,rn,rvxp,rvxn,rvyp,rvyn, &
                                rvzp,rvzn,ep,en
   real*8, dimension(:,:,:), allocatable :: sxp,sxn,syp,syn,szp,szn
!..local energy, lambda and split velocity arraies
   real*8, dimension(:,:,:), allocatable :: leng,llam,v0p,v0n,v1p,v1n
   integer :: ii,jj,kk,kkk

   do kkk=1,3

   ii=nx0
   jj=ny0
   kk=nz0

   if (kkk.eq.1) ii=nx0+1
   if (kkk.eq.2) jj=ny0+1
   if (kkk.eq.3) kk=nz0+1

   allocate(rhol(ii,jj,kk),rhor(ii,jj,kk),vxl(ii,jj,kk), &
            vxr(ii,jj,kk),vyl(ii,jj,kk),vyr(ii,jj,kk), &
            vzl(ii,jj,kk),vzr(ii,jj,kk),pl(ii,jj,kk), &
            pr(ii,jj,kk),bxl(ii,jj,kk),bxr(ii,jj,kk), &
            byl(ii,jj,kk),byr(ii,jj,kk),bzl(ii,jj,kk), &
            bzr(ii,jj,kk),lbsq(ii,jj,kk))
   allocate(rp(ii,jj,kk),rn(ii,jj,kk),rvxp(ii,jj,kk), &
            rvxn(ii,jj,kk),rvyp(ii,jj,kk),rvyn(ii,jj,kk), &
            rvzp(ii,jj,kk),rvzn(ii,jj,kk),ep(ii,jj,kk), &
            en(ii,jj,kk))
   allocate(sxp(ii,jj,kk),sxn(ii,jj,kk),syp(ii,jj,kk), &
            syn(ii,jj,kk),szp(ii,jj,kk),szn(ii,jj,kk))
   allocate(leng(ii,jj,kk),llam(ii,jj,kk),v0p(ii,jj,kk), &
            v0n(ii,jj,kk),v1p(ii,jj,kk),v1n(ii,jj,kk))

   call reconstruct_3d(rhol,rhor,lrho,ii,jj,kk,kkk)
   call reconstruct_3d(vxl,vxr,lvx,ii,jj,kk,kkk)
   call reconstruct_3d(vyl,vyr,lvy,ii,jj,kk,kkk)
   call reconstruct_3d(vzl,vzr,lvz,ii,jj,kk,kkk)
   call reconstruct_3d(pl,pr,lp,ii,jj,kk,kkk)
   call reconstruct_3d(bxl,bxr,lbx,ii,jj,kk,kkk)
   call reconstruct_3d(byl,byr,lby,ii,jj,kk,kkk)
   call reconstruct_3d(bzl,bzr,lbz,ii,jj,kk,kkk)

!..left flux
   leng=.5*rhol*(vxl**2+vyl**2+vzl**2)+pl/(gam-1.)
   llam=rhol/(2.*pl)

   if (kkk.eq.1) call get_v_face(llam,vxl,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.2) call get_v_face(llam,vyl,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.3) call get_v_face(llam,vzl,v0p,v0n,v1p,v1n,kkk)

   rp=rhol*v1p
   rvxp=rhol*vxl*v1p
   rvyp=rhol*vyl*v1p
   rvzp=rhol*vzl*v1p
   ep=(leng+.5*pl)*v1p

   if (kkk.eq.1) then
      rvxp=rvxp+pl*v0p
      ep=ep+.5*pl*vxl*v0p
   elseif (kkk.eq.2) then
      rvyp=rvyp+pl*v0p
      ep=ep+.5*pl*vyl*v0p
   elseif (kkk.eq.3) then
      rvzp=rvzp+pl*v0p
      ep=ep+.5*pl*vzl*v0p
   endif

   lbsq=bxl**2+byl**2+bzl**2
   llam=rhol/(2.*pl+lbsq)

   if (kkk.eq.1) call get_v_face(llam,vxl,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.2) call get_v_face(llam,vyl,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.3) call get_v_face(llam,vzl,v0p,v0n,v1p,v1n,kkk)

   if (kkk.eq.1) then
      sxp=(.5*lbsq-bxl*bxl)*v0p
      syp=(-byl*bxl)*v0p
      szp=(-bzl*bxl)*v0p
   elseif (kkk.eq.2) then
      sxp=(-bxl*byl)*v0p
      syp=(.5*lbsq-byl*byl)*v0p
      szp=(-bzl*byl)*v0p
   elseif (kkk.eq.3) then
      sxp=(-bxl*bzl)*v0p
      syp=(-byl*bzl)*v0p
      szp=(.5*lbsq-bzl*bzl)*v0p
   endif

!..right flux
   leng=.5*rhor*(vxr**2+vyr**2+vzr**2)+pr/(gam-1.)
   llam=rhor/(2.*pr)

   if (kkk.eq.1) call get_v_face(llam,vxr,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.2) call get_v_face(llam,vyr,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.3) call get_v_face(llam,vzr,v0p,v0n,v1p,v1n,kkk)

   rn=rhor*v1n
   rvxn=rhor*vxr*v1n
   rvyn=rhor*vyr*v1n
   rvzn=rhor*vzr*v1n
   en=(leng+.5*pr)*v1n

   if (kkk.eq.1) then
      rvxn=rvxn+pr*v0n
      en=en+.5*pr*vxr*v0n
   elseif (kkk.eq.2) then
      rvyn=rvyn+pr*v0n
      en=en+.5*pr*vyr*v0n
   elseif (kkk.eq.3) then
      rvzn=rvzn+pr*v0n
      en=en+.5*pr*vzr*v0n
   endif

   lbsq=bxr**2+byr**2+bzr**2
   llam=rhor/(2.*pr+lbsq)

   if (kkk.eq.1) call get_v_face(llam,vxr,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.2) call get_v_face(llam,vyr,v0p,v0n,v1p,v1n,kkk)
   if (kkk.eq.3) call get_v_face(llam,vzr,v0p,v0n,v1p,v1n,kkk)

   if (kkk.eq.1) then
      sxn=(.5*lbsq-bxr*bxr)*v0n
      syn=(-byr*bxr)*v0n
      szn=(-bzr*bxr)*v0n
   elseif (kkk.eq.2) then
      sxn=(-bxr*byr)*v0n
      syn=(.5*lbsq-byr*byr)*v0n
      szn=(-bzr*byr)*v0n
   elseif (kkk.eq.3) then
      sxn=(-bxr*bzr)*v0n
      syn=(-byr*bzr)*v0n
      szn=(.5*lbsq-bzr*bzr)*v0n
   endif

   if (kkk.eq.1) then
      flux_rhox=rp+rn
      flux_rvxx=rvxp+rvxn
      flux_rvyx=rvyp+rvyn
      flux_rvzx=rvzp+rvzn
      flux_engx=ep+en
      stress_xx=sxp+sxn
      stress_yx=syp+syn
      stress_zx=szp+szn
   elseif (kkk.eq.2) then
      flux_rhoy=rp+rn
      flux_rvxy=rvxp+rvxn
      flux_rvyy=rvyp+rvyn
      flux_rvzy=rvzp+rvzn
      flux_engy=ep+en
      stress_xy=sxp+sxn
      stress_yy=syp+syn
      stress_zy=szp+szn
   elseif (kkk.eq.3) then
      flux_rhoz=rp+rn
      flux_rvxz=rvxp+rvxn
      flux_rvyz=rvyp+rvyn
      flux_rvzz=rvzp+rvzn
      flux_engz=ep+en
      stress_xz=sxp+sxn
      stress_yz=syp+syn
      stress_zz=szp+szn
   endif

   deallocate(rhol,rhor,vxl,vxr,vyl,vyr,vzl,vzr,pl,pr,bxl,bxr,byl,byr,bzl, &
            bzr,lbsq)
   deallocate(rp,rn,rvxp,rvxn,rvyp,rvyn,rvzp,rvzn,ep,en)
   deallocate(sxp,sxn,syp,syn,szp,szn)
   deallocate(leng,llam,v0p,v0n,v1p,v1n)

   enddo ! enddo for kkk

   end subroutine get_flux_stress
!************************************
   subroutine get_v_center(lam,lvx,lvy,lvz) ! used for reconstruct e field only

   implicit none
   real*8, dimension(:,:,:), allocatable :: lam,lvx,lvy,lvz
   real*8, dimension(:,:,:), allocatable :: a,b

   allocate(a(nx,ny,nz),b(nx,ny,nz))

   a=sqrt(lam)

   b=.5*exp(-lam*lvx**2)/a/sqrt(pi)
   vx0p=.5*erfc(-a*lvx)
   vx0n=.5*erfc(a*lvx)
   vx1p=lvx*vx0p+b
   vx1n=lvx*vx0n-b

   b=.5*exp(-lam*lvy**2)/a/sqrt(pi)
   vy0p=.5*erfc(-a*lvy)
   vy0n=.5*erfc(a*lvy)
   vy1p=lvy*vy0p+b
   vy1n=lvy*vy0n-b

   b=.5*exp(-lam*lvz**2)/a/sqrt(pi)
   vz0p=.5*erfc(-a*lvz)
   vz0n=.5*erfc(a*lvz)
   vz1p=lvz*vz0p+b
   vz1n=lvz*vz0n-b

   deallocate(a,b)

   end subroutine get_v_center
!************************************
   subroutine get_v_face(llam,lv,v0p,v0n,v1p,v1n,dir)
   implicit none
   real*8, dimension(:,:,:), allocatable :: llam,lv,v0p,v0n,v1p,v1n
   real*8, dimension(:,:,:), allocatable :: a,b
   integer :: dir
!..add option for x,y,z only or all

   if (dir.eq.1) allocate(a(nx0+1,ny0,nz0),b(nx0+1,ny0,nz0))
   if (dir.eq.2) allocate(a(nx0,ny0+1,nz0),b(nx0,ny0+1,nz0))
   if (dir.eq.3) allocate(a(nx0,ny0,nz0+1),b(nx0,ny0,nz0+1))

   a=sqrt(llam)

   b=.5*exp(-llam*lv**2)/a/sqrt(pi)
   v0p=.5*erfc(-a*lv)
   v0n=.5*erfc(a*lv)
   v1p=lv*v0p+b
   v1n=lv*v0n-b

   deallocate(a,b)

   end subroutine get_v_face
!************************************

end module get

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!module output_br
!   use pars
!   use fields
!   implicit none
!
!contains
!   subroutine output_geom
!   implicit none
!
!   if (outop .eq. 2) then
!      open(unit=11,file=outdir//"xc.out",status="unknown",form="unformatted")
!      open(unit=12,file=outdir//"yc.out",status="unknown",form="unformatted")
!      open(unit=13,file=outdir//"zc.out",status="unknown",form="unformatted")
!
!      write(11) xc
!      write(12) yc
!      write(13) zc
!
!      close(unit=11)
!      close(unit=12)
!      close(unit=13)
!
!   endif
!
!   end subroutine output_geom
!!***********************************
!   subroutine output_simt
!   implicit none
!
!   if (ntt.eq.0 .and. nt.eq.0) then
!      open(unit=11,file=outdir//"stime.out",status="unknown",form="unformatted")
!   else
!      open(unit=11,file=outdir//"stime.out",status="unknown",form="unformatted",access="append")
!   endif
!      write(11) time
!      close(11)
!
!   end subroutine output_simt
!!***********************************
!   subroutine outputx
!   implicit none
!   integer :: j,k
!
!   j=ng+ny0/2
!   k=ng+1
!!   k=ng+nz0/2
!
!   if (outop .eq. 2) then
!      if (ntt .eq. 0 .and. nt .eq. 0) then
!         open(unit=11,file=outdir//"rho_x.out",status="unknown",form="unformatted")
!         open(unit=12,file=outdir//"vx_x.out",status="unknown",form="unformatted")
!         open(unit=13,file=outdir//"vy_x.out",status="unknown",form="unformatted")
!         open(unit=14,file=outdir//"vz_x.out",status="unknown",form="unformatted")
!         open(unit=15,file=outdir//"p_x.out",status="unknown",form="unformatted")
!         open(unit=16,file=outdir//"bx_x.out",status="unknown",form="unformatted")
!         open(unit=17,file=outdir//"by_x.out",status="unknown",form="unformatted")
!         open(unit=18,file=outdir//"bz_x.out",status="unknown",form="unformatted")
!      else
!         open(unit=11,file=outdir//"rho_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=12,file=outdir//"vx_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=13,file=outdir//"vy_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=14,file=outdir//"vz_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=15,file=outdir//"p_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=16,file=outdir//"bx_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=17,file=outdir//"by_x.out",status="unknown",form="unformatted",access="append")
!         open(unit=18,file=outdir//"bz_x.out",status="unknown",form="unformatted",access="append")
!      endif
!
!      write(11) rho(1:nx,j,k)
!      write(12) vx(1:nx,j,k)
!      write(13) vy(1:nx,j,k)
!      write(14) vz(1:nx,j,k)
!      write(15) p(1:nx,j,k)
!      write(16) bx(1:nx,j,k)
!      write(17) by(1:nx,j,k)
!      write(18) bz(1:nx,j,k)
!
!      close(unit=11)
!      close(unit=12)
!      close(unit=13)
!      close(unit=14)
!      close(unit=15)
!      close(unit=16)
!      close(unit=17)
!      close(unit=18)
!
!   endif
!
!   end subroutine outputx
!!***********************************
!   subroutine outputxy
!   implicit none
!   integer :: k
!
!   k=ng+1
!
!   if (outop .eq. 2) then
!      if (ntt .eq. 0 .and. nt .eq. 0) then
!         open(unit=11,file=outdir//"rho_xy.out",status="unknown",form="unformatted")
!         open(unit=12,file=outdir//"vx_xy.out",status="unknown",form="unformatted")
!         open(unit=13,file=outdir//"vy_xy.out",status="unknown",form="unformatted")
!         open(unit=14,file=outdir//"vz_xy.out",status="unknown",form="unformatted")
!         open(unit=15,file=outdir//"p_xy.out",status="unknown",form="unformatted")
!         open(unit=16,file=outdir//"bx_xy.out",status="unknown",form="unformatted")
!         open(unit=17,file=outdir//"by_xy.out",status="unknown",form="unformatted")
!         open(unit=18,file=outdir//"bz_xy.out",status="unknown",form="unformatted")
!      else
!         open(unit=11,file=outdir//"rho_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=12,file=outdir//"vx_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=13,file=outdir//"vy_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=14,file=outdir//"vz_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=15,file=outdir//"p_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=16,file=outdir//"bx_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=17,file=outdir//"by_xy.out",status="unknown",form="unformatted",access="append")
!         open(unit=18,file=outdir//"bz_xy.out",status="unknown",form="unformatted",access="append")
!      endif
!
!      write(11) rho(1:nx,1:ny,k)
!      write(12) vx(1:nx,1:ny,k)
!      write(13) vy(1:nx,1:ny,k)
!      write(14) vz(1:nx,1:ny,k)
!      write(15) p(1:nx,1:ny,k)
!      write(16) bx(1:nx,1:ny,k)
!      write(17) by(1:nx,1:ny,k)
!      write(18) bz(1:nx,1:ny,k)
!
!      close(unit=11)
!      close(unit=12)
!      close(unit=13)
!      close(unit=14)
!      close(unit=15)
!      close(unit=16)
!      close(unit=17)
!      close(unit=18)
!
!   endif
!
!   end subroutine outputxy
!!***********************************
!   subroutine outputyz
!   implicit none
!   integer :: i
!
!   i=ng+nx0/2
!
!   if (outop .eq. 2) then
!      if (ntt .eq. 0 .and. nt .eq. 0) then
!         open(unit=11,file=outdir//"rho_yz.out",status="unknown",form="unformatted")
!         open(unit=12,file=outdir//"vx_yz.out",status="unknown",form="unformatted")
!         open(unit=13,file=outdir//"vy_yz.out",status="unknown",form="unformatted")
!         open(unit=14,file=outdir//"vz_yz.out",status="unknown",form="unformatted")
!         open(unit=15,file=outdir//"p_yz.out",status="unknown",form="unformatted")
!         open(unit=16,file=outdir//"bx_yz.out",status="unknown",form="unformatted")
!         open(unit=17,file=outdir//"by_yz.out",status="unknown",form="unformatted")
!         open(unit=18,file=outdir//"bz_yz.out",status="unknown",form="unformatted")
!      else
!         open(unit=11,file=outdir//"rho_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=12,file=outdir//"vx_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=13,file=outdir//"vy_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=14,file=outdir//"vz_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=15,file=outdir//"p_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=16,file=outdir//"bx_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=17,file=outdir//"by_yz.out",status="unknown",form="unformatted",access="append")
!         open(unit=18,file=outdir//"bz_yz.out",status="unknown",form="unformatted",access="append")
!      endif
!
!      write(11) rho(i,1:ny,1:nz)
!      write(12) vx(i,1:ny,1:nz)
!      write(13) vy(i,1:ny,1:nz)
!      write(14) vz(i,1:ny,1:nz)
!      write(15) p(i,1:ny,1:nz)
!      write(16) bx(i,1:ny,1:nz)
!      write(17) by(i,1:ny,1:nz)
!      write(18) bz(i,1:ny,1:nz)
!
!      close(unit=11)
!      close(unit=12)
!      close(unit=13)
!      close(unit=14)
!      close(unit=15)
!      close(unit=16)
!      close(unit=17)
!      close(unit=18)
!
!   endif
!
!   end subroutine outputyz
!!***********************************
!end module output_br
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module output_h5
   use pars
   use fields
   use H5LT
   use HDF5
!..USER INPUT: directory naming mode, =0 for YYMMDD_#, =1 for jobid
!   integer :: dirop=1
   character(70) :: dirname
   integer(HID_T) :: plistmpi_id         ! property list identifiers
   integer(HID_T) :: plistxfercol_id

   contains

subroutine init_output_h5

   implicit none
   integer :: datecountn
   character(10) :: folder
   character(6) :: runno, cdate
   integer :: psave, errcode, datetime(8)
   logical :: existe=.FALSE., existe2=.TRUE.
   integer(HSIZE_T) :: dims(1)
   integer(HID_T) :: setup_fid  ! simulation setup file identifier

!..Initialize HDF5 FORTRAN interface.
   call h5open_f(errcode)

!.Broadcast the directory name where results are stored
  dirname=trim(outdir)
  call MPI_BCAST(dirname,70,MPI_CHARACTER,0,COMM_CART,errcode)

!..Setup file access property list with parallel I/O access
   call h5pcreate_f(H5P_FILE_ACCESS_F,plistmpi_id,errcode)
   call h5pset_fapl_mpio_f(plistmpi_id,MPI_COMM_WORLD,MPI_INFO_NULL,errcode)
!..For collective access
   call h5pcreate_f(H5P_DATASET_XFER_F,plistxfercol_id,errcode)
   call h5pset_dxpl_mpio_f(plistxfercol_id,H5FD_MPIO_COLLECTIVE_F,errcode)

   call createh5_2D('rho',2,nframes+1)
   call createh5_2D('vx',2,nframes+1)
   call createh5_2D('vy',2,nframes+1)
   call createh5_2D('vz',2,nframes+1)
   call createh5_2D('p',2,nframes+1)
   call createh5_2D('bx',2,nframes+1)
   call createh5_2D('by',2,nframes+1)
   call createh5_2D('bz',2,nframes+1)

  end subroutine init_output_h5
!*********************************************************
!..Create file dataname.h5, with dataset /data (identified by dsetid) with
!..dimensions [nx0G,ny0G,nframes+1] (ONE FOR EACH ITERATION + 1 INITIAL)
   subroutine createh5_2D(dataname,prec,frames)
   implicit none
   character(len=*), intent(in) :: dataname
   integer(HID_T):: dsetid, fid, dspaceid
   integer, intent(in) :: prec,frames
   integer(HSIZE_T) :: dims(3)
   integer :: errcode
   character(15) :: dset

!.Define dimension of dataset in file
!...Create a new file using DEFAULT PROPERTIES
   call h5fcreate_f(trim(dirname)//trim(dataname)//'.h5',H5F_ACC_TRUNC_F, &
                   fid,errcode,access_prp=plistmpi_id)
   dims = (/nx0G, ny0G, frames/)
  
!..Create a 3D (2Dspace+time) dataspace in new file
   call h5screate_simple_f(3,dims,dspaceid,errcode)

!..Create a dataset with DEFAULT PROPERTIES
   dset = trim(adjustl('/'//dataname))
   if (prec==1) then     ! single precision
!..BZ: THIS currently CORRUPTS data, at least as seen in MATLAB pcolor plots
      call h5dcreate_f(fid,dset,H5T_NATIVE_REAL,dspaceid,dsetid,errcode)
   elseif (prec==2) then ! double precision
      call h5dcreate_f(fid,dset,H5T_NATIVE_DOUBLE,dspaceid,dsetid,errcode)
   endif

!..Close created file. Each output call opens and closes the file
   call h5sclose_f(dspaceid, errcode)
   call h5dclose_f(dsetid, errcode)
   call h5fclose_f(fid, errcode)

   end subroutine createh5_2D
!*********************************************************
!..Output 2D data to subset of 3D dataspace in dataname.h5
   subroutine outh5_2D(dataname,outdata,prec,toff)
   implicit none
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: toff, prec
   real*8, intent(in) :: outdata(:,:)
   integer(HID_T) :: dsetid,dspaceid,fid,mspaceid
   integer(HSIZE_T) :: dims(2)
   integer(HSIZE_T) :: count(3), offset(3)
   integer :: errcode, ip1, ip2
   character(15) :: dset

!..Open file
   call h5fopen_f(trim(dirname)//trim(dataname)//'.h5', H5F_ACC_RDWR_F, &
                   fid, errcode, access_prp=plistmpi_id)
   dims = (/nx0G,ny0G/)                    ! Shape of xy data
   count = (/nx0,ny0,1/)                   ! Size of hyperslab
!..Hyperslab offset
   offset(1) = xID*count(1)
   offset(2) = yID*count(2)
   offset(3) = toff

!..Open existing dataset/get dataset identifier
   dset = trim(adjustl('/'//dataname))
   call h5dopen_f(fid,dset,dsetid,errcode)

!..Create memory dataspace in selected hyperslab
   call h5screate_simple_f(2,count(1:2),mspaceid,errcode)

!..Open dataspace/get dataspace ID
   call h5dget_space_f(dsetid,dspaceid,errcode)

!..Select subset of 3D dataspace to store 3D data
   call h5sselect_hyperslab_f(dspaceid,H5S_SELECT_SET_F,offset, &
                             count,errcode)

!..Write 2D data to 2D subset of 3D memspace in dataname.h5
   if (prec==1) then        ! single precision
!..BZ: THIS currently CORRUPTS data, at least as
!.......seen in MATLAB pcolor plots
      call h5dwrite_f(dsetid,H5T_NATIVE_REAL,SNGL(outdata),dims, &
                    errcode,mspaceid,dspaceid)
!...Or for collective write, in the case of plane=1, use this call
!...BZ: This wasn't working last time I trid it
!    h5dwrite_f(dsetid,H5T_NATIVE_FLOAT,SNGL(outdata),dims, &
!               errcode,mspaceid,dspaceid,xfer_prp=plistxfercol_id)
   elseif (prec==2) then    ! double precision
      call h5dwrite_f(dsetid,H5T_NATIVE_DOUBLE,outdata,dims,errcode, &
                     mspaceid,dspaceid,xfer_prp=plistxfercol_id)
   endif

!  call MPI_BARRIER(COMM_CART,ierr)

!.Terminate access to the dataset, memoryspace, dataspace and file.
   call h5sclose_f(dspaceid, errcode)
   call h5sclose_f(mspaceid, errcode)
   call h5dclose_f(dsetid, errcode)
   call h5fclose_f(fid, errcode)

   end subroutine outh5_2D
!*********************************************************
!.Close FORTRAN HDF5 interface.
  subroutine terminate_output_h5
  integer :: errcode

!.Close property lists
  call h5pclose_f(plistmpi_id, errcode)
  call h5pclose_f(plistxfercol_id, errcode)
  call h5close_f(errcode)

  end subroutine terminate_output_h5
!*********************************************************
subroutine outputxy
   implicit none
   integer :: k
   k=ng+1

   call outh5_2D('rho',rho(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
   call outh5_2D('vx',vx(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
   call outh5_2D('vy',vy(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
!   call outh5_2D('vz',vz(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
   call outh5_2D('p',p(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
   call outh5_2D('bx',bx(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
   call outh5_2D('by',by(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)
!   call outh5_2D('bz',bz(ng+1:nx-ng,ng+1:ny-ng,k),2,ntt)

end subroutine outputxy
!*********************************************************
!..should be changed to h5 output
   subroutine output_simt
   implicit none

   if (myrank.eq.0) then
   if (ntt.eq.0 .and. nt.eq.0) then
      open(unit=11,file=outdir//"stime.out",status="unknown",form="unformatted")
   else
      open(unit=11,file=outdir//"stime.out",status="unknown",form="unformatted",access="append")
   endif

   write(11) time
   close(11)

   endif

   end subroutine output_simt
!*********************************************************
end module output_h5
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module diagnostic
  use pars
  use fields
  implicit none

contains
  subroutine diag_geom
  implicit none
  if (minval(dxblock).lt.0 .or. minval(dyblock).lt.0 .or. &
      minval(dzblock).lt.0) stop 'Wrong grid! Simulation terminated!'
  end subroutine diag_geom
!***********************************
  subroutine diag_rho_neg
  implicit none
  if (minval(rho).lt.0.) stop 'Density becomes negative! Simulation terminated!'
  end subroutine diag_rho_neg
!***********************************
   subroutine diag_overflow
   implicit none
   integer :: i,j,k

   do i=1,nx0
   do j=1,ny0
   do k=1,nz0
      if (isnan(rho(ng+i,ng+j,ng+k))) stop 'NaN value found. Simulation terminated!'
   enddo
   enddo
   enddo

   end subroutine diag_overflow
!***********************************
end module diagnostic
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!..main program:
   program mhd3d
   use pars
   use fields
   use get
!   use output_br
   use output_h5
   use diagnostic

   implicit none
   integer :: i,j,IDs(2)
   character(8) :: rdate
   character(10) :: rtime

!--Initialize MPI
   call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,nprovided,ierr)  ! Required statement
   print*,'MPI_THREAD_FUNNLELED:', MPI_THREAD_FUNNELED
   print*,'nprovided:', nprovided

!--Who am I? --- get my rank
   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

!--How many processes in the global group?
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

!..read input parameters:
   call read_input
   call comp_pars
   call alloc_pars

!..Check that the number of processes launched is correct
   if (numprocs .ne. nprocs) then
      print*,'nprocs parameter in code differs from that in mpirun -np xxx'
      print*,'nprocs in code is',nprocs
      go to 99
   endif

   call alloc_fields

!--MPI_CART boundary conditions along x and y. TRUE=periodic
!   CART_BCs(1) = .FALSE.  ! not periodic along y
!   CART_BCs(2) = .FALSE.  ! not periodic along x
   CART_BCs(1) = .TRUE.
   CART_BCs(2) = .TRUE.

!--Create a Cartesian topology
   call MPI_CART_CREATE(MPI_COMM_WORLD,2,MPIgrid,CART_BCs, &
                        reorder,COMM_CART,ierror)
!--Find my coordinate parameters in the Cartesial topology
!  and use them to compute the arrays x,y
   call MPI_CART_COORDS(COMM_CART,myrank,2,blockID,ierror)

!................X and Y SUBCOMMUNICATORS...................!
!..Create x subcommunicators
   remain(1)=.FALSE.     ! Don't keep first (y) direction
   remain(2)=.TRUE.      ! Keep 2nd (x) direction)
   call MPI_Cart_sub(COMM_CART, remain, COMM_X, ierror)
   call MPI_Comm_rank(COMM_X, xID, ierr)
   call MPI_Cart_coords(COMM_X, xID, 1, xcoord, ierror)
!..For communication along x, save source and destination rank
   call MPI_CART_SHIFT(COMM_X,0,1,left,right,ierror)
!..These two lines maybe needed for some clusters
!   if (xID==0) left=MPI_PROC_NULL
!   if (xID==xprocs-1) right=MPI_PROC_NULL

!..Create y subcommunicators
   remain(1)=.TRUE.      ! Keep first (y) direction
   remain(2)=.FALSE.     ! Don't keep 2nd (x) direction
   call MPI_Cart_sub(COMM_CART, remain, COMM_Y, ierror)
   call MPI_Comm_rank(COMM_Y, yID, ierr)
   call MPI_Cart_coords(COMM_Y, yID, 1, ycoord, ierror)
!..For communication along y, save source and destination rank
   call MPI_CART_SHIFT(COMM_Y,0,1,down,up,ierror)
!..These two lines maybe needed for some clusters
!   if (yID==0) down=MPI_PROC_NULL
!   if (yID==yprocs-1) up=MPI_PROC_NULL

!............................................................!

!..initialize the fields and other quantities, read the input file, etc:
   call initialization

   ntt=0
   nt =0

!..The run starts at t=0:
   time=0.

   call init_output_h5
   call outputxy
   call output_simt

!....Report real start date and time
   if (myrank.eq.0) then
      call date_and_time(rdate,rtime)
      print*,'start date & time=',rdate,' ',rtime
   endif

!..Here tau is time step and the total number of steps is nts*nframes;
!..nframes is the number of frames (time-slices) each output file will contain,
!  and each frame is separated by nts time steps

   do ntt=1,nframes
      do nt=1,nts
!         time=time+tau
         call step
!         call step_ab
!         call step_rk2
!         call step_rk3
         time=time+tau
!         if(myrank.eq.0) print*,'time=', time, 'output frame number ', ntt
      enddo
      call outputxy
      call output_simt

      if (myrank.eq.0) then
         call date_and_time(rdate,rtime)
         print*,'output frame number',ntt+1, 'date & time=',rdate,' ',rtime
      endif

      call diag_rho_neg
      call diag_overflow
   enddo

   call terminate_output_h5
   if(myrank.eq.0) print*,'hdf5 output terminated'

   call termination

99 continue
   call dealloc_pars
!..Report real end date and time
   if (myrank.eq.0) then
      call date_and_time(rdate,rtime)
      print*,'end date & time=',rdate,' ',rtime
   endif

!--Finilize MPI
   call MPI_FINALIZE(irc)
   stop

! 333  format(3(e12.6,1x))
! 444  format(4(e12.6,1x))
! 555  format(5(e12.6,1x))
! 666  format(6(e12.6,1x))

   end program

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine initialization
   use pars
   use fields
   use boundary
   use get
!   use output_br
   use output_h5
   use diagnostic

   implicit none

   call init_geom
   call diag_geom
!   call output_geom

!   call init_ic_bws
   call init_ic_otv
!   call init_ic_rec

!   call fill_hydro
!   call fill_mag
   call fill_hydro_mpi
   call fill_mag_mpi

   call get_bxyz
!..these get_eijk routine are used for diagnostic
!   call get_eijk_Yee(rho,vx,vy,vz,p,bx,by,bz)
   call get_eijk_LFM(rho,vx,vy,vz,p,bx,by,bz,bi,bj,bk)

   if (stepop.eq.1) then
   call get_dt
!   tau=5.e-6
   tau0=tau

   rhop=rho
   vxp=vx
   vyp=vy
   vzp=vz
   pp=p
   bxp=bx
   byp=by
   bzp=bz
   bip=bi
   bjp=bj
   bkp=bk
   endif

   if(myrank.eq.0) print*, 'initialization finished!'

end subroutine initialization
!*********************************************
subroutine init_geom
   use pars
   implicit none
   real*8, allocatable :: x0(:),y0(:),z0(:)
   real*8, allocatable :: xx0(:),yy0(:),zz0(:)
   real*8 :: scx,scy,scz
   integer :: i,j,k

!..first define global grid
   allocate(x0(nx0G+1),y0(ny0G+1),z0(nz0+1))
   allocate(xx0(nxG+1),yy0(nyG+1),zz0(nz+1))

   dx0=(xmax-xmin)/nx0G
   dy0=(ymax-ymin)/ny0G
   dz0=(zmax-zmin)/nz0   

   forall(i=1:nx0G+1) x0(i)=dble(i-1)*dx0+xmin
   forall(j=1:ny0G+1) y0(j)=dble(j-1)*dy0+ymin
   forall(k=1:nz0+1) z0(k)=dble(k-1)*dz0+zmin

!..uniform grid
!   forall(i=1:nx0G+1) xx0(ng+i)=x0(i)
!   forall(j=1:ny0G+1) yy0(ng+j)=y0(j)
!   forall(k=1:nz0+1) zz0(ng+k)=z0(k)
   xx0(ng+1:ng+nx0G+1)=x0
   yy0(ng+1:ng+ny0G+1)=y0
   zz0(ng+1:ng+nz0+1)=z0

!!..non-uniform grid test
!   scx=30./tan(xmax)
!   scy=60./tan(ymax)
!   scz=60./tan(zmax)
!
!   x0=scx*tan(x0)
!   y0=scy*tan(y0)
!   z0=scz*tan(z0)

!..define ghost cells
   forall(i=1:ng)
      xx0(i)=2.*x0(1)-x0(ng+2-i)
      xx0(nxG+2-i)=2.*x0(nx0G+1)-x0(nx0G-ng+i)
      yy0(i)=2.*y0(1)-y0(ng+2-i)
      yy0(nyG+2-i)=2.*y0(ny0G+1)-y0(ny0G-ng+i)
      zz0(i)=i-6!2.*z0(1)-z0(ng+2-i)
      zz0(nz+2-i)=6-i!2.*z0(nz0+1)-z0(nz0-ng+i)
   endforall  

!   xx0(ng+1:ng+nx0G+1)=x0
!   yy0(ng+1:ng+ny0G+1)=y0
!   zz0(ng+1:ng+nz0+1)=z0

   forall(i=1:nxG+1) x(i,:,:)=xx0(i)
   forall(j=1:nyG+1) y(:,j,:)=yy0(j)
   forall(k=1:nz+1) z(:,:,k)=zz0(k)

!!>>>>>>> ideally,this part is no longer needed <<<<<<<<<<<<<
!!..location of cell centers / length of each cell edges
!   forall(i=1:nxG,j=1:nyG,k=1:nz)
!      xc(i,j,k)=.125*(x(i,j,k)+x(i,j,k+1)+x(i+1,j,k+1)+x(i+1,j,k)+ &
!                      x(i,j+1,k)+x(i,j+1,k+1)+x(i+1,j+1,k+1)+x(i+1,j+1,k))
!      yc(i,j,k)=.125*(y(i,j,k)+y(i,j,k+1)+y(i+1,j,k+1)+y(i+1,j,k)+ &
!                      y(i,j+1,k)+y(i,j+1,k+1)+y(i+1,j+1,k+1)+y(i+1,j+1,k))
!      zc(i,j,k)=.125*(z(i,j,k)+z(i,j,k+1)+z(i+1,j,k+1)+z(i+1,j,k)+ &
!                      z(i,j+1,k)+z(i,j+1,k+1)+z(i+1,j+1,k+1)+z(i+1,j+1,k))
!      dx(i,j,k)=x(i+1,j,k)-x(i,j,k)
!      dy(i,j,k)=y(i,j+1,k)-y(i,j,k)
!      dz(i,j,k)=z(i,j,k+1)-z(i,j,k)
!   endforall
!!..location of i-face center - where bi is define on
!   forall(j=1:nyG,k=1:nz)
!      xi(:,j,k)=x(:,j,k)
!      yi(:,j,k)=.25*(y(:,j,k)+y(:,j+1,k)+y(:,j,k+1)+y(:,j+1,k+1))
!      zi(:,j,k)=.25*(z(:,j,k)+z(:,j+1,k)+z(:,j,k+1)+z(:,j+1,k+1))
!   endforall
!!..location of j-face center - where bj is define on
!   forall(i=1:nxG,k=1:nz)
!      xj(i,:,k)=.25*(x(i,:,k)+x(i+1,:,k)+x(i,:,k+1)+x(i+1,:,k+1))
!      yj(i,:,k)=y(i,:,k)
!      zj(i,:,k)=.25*(z(i,:,k)+z(i+1,:,k)+z(i,:,k+1)+z(i+1,:,k+1))
!   endforall
!!..location of k-face center - where bk is define on   
!   forall(i=1:nx,j=1:nyG)
!      xk(i,j,:)=.25*(x(i,j,:)+x(i+1,j,:)+x(i,j+1,:)+x(i+1,j+1,:))
!      yk(i,j,:)=.25*(y(i,j,:)+y(i+1,j,:)+y(i,j+1,:)+y(i+1,j+1,:))
!      zk(i,j,:)=z(i,j,:)
!   endforall
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!..now, define local grid
   forall(i=1:nx+1) xblock(i,:,:)=xx0(i+xID*nx0)
   forall(j=1:ny+1) yblock(:,j,:)=yy0(j+yID*ny0)
   forall(k=1:nz+1) zblock(:,:,k)=zz0(k)

!..location of cell centers / length of each cell edges
   forall(i=1:nx,j=1:ny,k=1:nz)
      xcblock(i,j,k)=.125*(xblock(i,j,k)+xblock(i,j,k+1)+xblock(i+1,j,k+1) + &
                           xblock(i+1,j,k)+xblock(i,j+1,k)+xblock(i,j+1,k+1) + &
                           xblock(i+1,j+1,k+1)+xblock(i+1,j+1,k))
      ycblock(i,j,k)=.125*(yblock(i,j,k)+yblock(i,j,k+1)+yblock(i+1,j,k+1) + &
                           yblock(i+1,j,k)+yblock(i,j+1,k)+yblock(i,j+1,k+1) + &
                           yblock(i+1,j+1,k+1)+yblock(i+1,j+1,k))
      zcblock(i,j,k)=.125*(zblock(i,j,k)+zblock(i,j,k+1)+zblock(i+1,j,k+1) + & 
                           zblock(i+1,j,k)+zblock(i,j+1,k)+zblock(i,j+1,k+1) + &
                           zblock(i+1,j+1,k+1)+zblock(i+1,j+1,k))
      dxblock(i,j,k)=xblock(i+1,j,k)-xblock(i,j,k)
      dyblock(i,j,k)=yblock(i,j+1,k)-yblock(i,j,k)
      dzblock(i,j,k)=zblock(i,j,k+1)-zblock(i,j,k)
   endforall

!..location of i-face center - where bi is define on
   forall(j=1:ny,k=1:nz)
      xiblock(:,j,k)=xblock(:,j,k)
      yiblock(:,j,k)=.25*(yblock(:,j,k)+yblock(:,j+1,k)+yblock(:,j,k+1)+yblock(:,j+1,k+1))
      ziblock(:,j,k)=.25*(zblock(:,j,k)+zblock(:,j+1,k)+zblock(:,j,k+1)+zblock(:,j+1,k+1))
   endforall
!..location of j-face center - where bj is define on
   forall(i=1:nx,k=1:nz)
      xjblock(i,:,k)=.25*(xblock(i,:,k)+xblock(i+1,:,k)+xblock(i,:,k+1)+xblock(i+1,:,k+1))
      yjblock(i,:,k)=yblock(i,:,k)
      zjblock(i,:,k)=.25*(zblock(i,:,k)+zblock(i+1,:,k)+zblock(i,:,k+1)+zblock(i+1,:,k+1))
   endforall
!..location of k-face center - where bk is define on   
   forall(i=1:nx,j=1:ny)
      xkblock(i,j,:)=.25*(xblock(i,j,:)+xblock(i+1,j,:)+xblock(i,j+1,:)+xblock(i+1,j+1,:))
      ykblock(i,j,:)=.25*(yblock(i,j,:)+yblock(i+1,j,:)+yblock(i,j+1,:)+yblock(i+1,j+1,:))
      zkblock(i,j,:)=zblock(i,j,:)
   endforall

!..for large size problem, for example 512*512*512 or 1024*1024*1024, it maybe
!..necessary to deallocate global quantities, eg, x,dx,xc,xi,xj,xk and so on here
!.. instead of in the dealloc_par subroutine to save memory

   deallocate(x0,y0,z0,xx0,yy0,zz0)

end subroutine init_geom
!*********************************************
subroutine init_ic_bws
   use pars
   use fields
   implicit none
   integer :: i,j,k

!..simple test initial condition
   rho=1.
!   p=1.
   p=.1

   vx=0.
   vy=0.
   vz=0.

!..1d test for x direction
!   bi=.75
!   bj=-1.
!   bk=0.
!..1d test for y direction
!   bi=0.
!   bj=.75
!   bk=-1.
!..1d test for z direction
!   bi=-1.
!   bj=0.
!   bk=.75


   do i=1,nx
   do j=1,ny
   do k=1,nz
!      if (xcblock(i,j,k).ge.0.) then ! for x
!      if (ycblock(i,j,k).ge.0.) then ! for y
!      if (zcblock(i,j,k).ge.0.) then ! for z
      if ((xcblock(i,j,k)**2+ycblock(i,j,k)**2).lt..01) then
         p(i,j,k)=10.
!         rho(i,j,k)=.125
!         p(i,j,k)=.1
      endif
   enddo
   enddo
   enddo

!   do i=1,nx
   do i=1,nx+1
   do j=1,ny
!   do j=1,ny+1
   do k=1,nz
!   do k=1,nz+1
!      if (xjblock(i,j,k).lt.0) bj(i,j,k)=1. ! for x
!      if (ykblock(i,j,k).lt.0.) bk(i,j,k)=1. ! for y
!      if (ziblock(i,j,k).lt.0.) bi(i,j,k)=1. ! for z
   enddo
   enddo
   enddo

   bi=.5
   bj=.5

end subroutine init_ic_bws
!*********************************************
subroutine init_ic_otv
   use pars
   use fields
   implicit none
   integer :: i,j,k

!..vortex problem
   rho=25./(36.*pi)
   p=5./(12.*pi)
   vx=-sin(2.*pi*ycblock)
   vy=sin(2.*pi*xcblock)
   vz=0.
   bi=-1./sqrt(4.*pi)*sin(2.*pi*yiblock)
   bj=1./sqrt(4.*pi)*sin(4.*pi*xjblock)
!   bi=0.
!   bj=0.
   bk=0.

end subroutine init_ic_otv
!*********************************************
subroutine init_ic_imp
   use pars
   use fields
   implicit none
   integer :: i,j,k

!..hd implosion test
   bi=0.
   bj=0.
   bk=0.

!   rho=.125
!   p=.14
   do i=1,nx
   do j=1,ny
   do k=1,nz
      if ((xcblock(i,j,k)+ycblock(i,j,k)).ge..150001) then
         rho(i,j,k)=1.
         p(i,j,k)=1.
      else
         rho(i,j,k)=.125
         p(i,j,k)=.14
      endif
   enddo
   enddo
   enddo

   vx=0.
   vy=0.
   vz=0.

end subroutine init_ic_imp
!*********************************************
subroutine step
   use pars
   implicit none

   if (stepop.eq.1) then
      call step_ab
   elseif (stepop.eq.2) then
      call step_rk2
   elseif (stepop.eq.3) then
      call step_rk3
   else
      print*,'Undefined time advance option!'
      stop
   endif

end subroutine step
!*********************************************
!..Adams-Bashforth method
subroutine step_ab
   use pars
   use fields
   use boundary
   use get
   implicit none

   call get_dt
!!..if fix dt
!   if (tau.lt.tau0) then
!      if (myrank.eq.0) print*, 'Attention! tau_cfl<tau0'
!   endif
!   tau=tau0

   call get_conserve

   call advance_predictor_ab
   tau0=tau

!   call get_eijk_Yee(rhoh,vxh,vyh,vzh,ph,bxh,byh,bzh)
   call get_eijk_LFM(rhoh,vxh,vyh,vzh,ph,bxh,byh,bzh,bih,bjh,bkh)

   call get_flux_stress(rhoh,vxh,vyh,vzh,ph,bxh,byh,bzh)

   call advance_corrector1
   call advance_corrector2

   call advance_bijk

   call fill_hydro_mpi
   call fill_mag_mpi
 
   call get_bxyz

!..testing
!   if (myrank.eq.1) then
!      print*, myrank,xID,yID
!      print*, 'rho'
!      print*, rho(5,5,:)
!   endif

end subroutine step_ab
!*********************************************
!..2rd order TVD Runge-Kutta scheme
!..v1 for ideal (no hall, resistive) case, no subcycling
subroutine step_rk2 ! version1
   use pars
   use fields
   use boundary
   use get
   implicit none
   integer :: rk_step

   call get_dt
!..if fix dt
   if (tau.lt.tau0) then
      print*, myrank, 'Attention! tau_cfl<', tau0
   endif
   tau=tau0

   call advance_previous

   do rk_step=1,2

   call get_conserve

   call get_eijk_LFM(rho,vx,vy,vz,p,bx,by,bz,bi,bj,bk)

   call get_flux_stress(rho,vx,vy,vz,p,bx,by,bz) 

   call advance_corrector1
   if (rk_step.eq.1) then
      rhovxh=rhovx
      rhovyh=rhovy
      rhovzh=rhovz
   elseif (rk_step.eq.2) then
      rho=(rhop+rho)/2.
      rhovxh=(rhovxp+rhovx)/2.
      rhovyh=(rhovyp+rhovy)/2.
      rhovzh=(rhovzp+rhovz)/2.
      eng=(engp+eng)/2.
   endif
   p=(gam-1.)*(eng-.5*(rhovxh**2+rhovyh**2+rhovzh**2)/rho)

   call advance_corrector2
   if (rk_step.eq.2) then
      rhovx=(rhovxp+rhovx)/2.
      rhovy=(rhovyp+rhovy)/2.
      rhovz=(rhovzp+rhovz)/2.
   endif
   vx=rhovx/rho
   vy=rhovy/rho
   vz=rhovz/rho

   call advance_bijk
   if (rk_step.eq.2) then
      bi=(bip+bi)/2.
      bj=(bjp+bj)/2.
      bk=(bkp+bk)/2.
   endif

   call fill_hydro_mpi
   call fill_mag_mpi

   call get_bxyz

   enddo ! enddo rk_step

end subroutine step_rk2 ! v1
!*********************************************
!..3rd order TVD Runge-Kutta scheme
!..v1 for ideal (no hall, resistive) case, no subcycling
subroutine step_rk3 ! version 1
   use pars
   use fields
   use boundary
   use get
   implicit none
   integer :: rk_step

   call get_dt
!..if fix dt
   if (tau.lt.tau0) then
      print*, myrank, 'Attention! tau_cfl<',tau0
   endif
   tau=tau0

   call advance_previous

   do rk_step=1,3

   call get_conserve

   call get_eijk_LFM(rho,vx,vy,vz,p,bx,by,bz,bi,bj,bk)

   call get_flux_stress(rho,vx,vy,vz,p,bx,by,bz)

   call advance_corrector1
   if (rk_step.eq.1) then
      rhovxh=rhovx
      rhovyh=rhovy
      rhovzh=rhovz
   elseif (rk_step.eq.2) then
      rho=(3.*rhop+rho)/4.
      rhovxh=(3.*rhovxp+rhovx)/4.
      rhovyh=(3.*rhovyp+rhovy)/4.
      rhovzh=(3.*rhovzp+rhovz)/4.
      eng=(3.*engp+eng)/4.
   elseif (rk_step.eq.3) then
      rho=(rhop+2.*rho)/3.
      rhovxh=(rhovxp+2.*rhovx)/3.
      rhovyh=(rhovyp+2.*rhovy)/3.
      rhovzh=(rhovzp+2.*rhovz)/3.
      eng=(engp+2.*eng)/3.
   endif
   p=(gam-1.)*(eng-.5*(rhovxh**2+rhovyh**2+rhovzh**2)/rho)

   call advance_corrector2
   if (rk_step.eq.2) then
      rhovx=(3.*rhovxp+rhovx)/4.
      rhovy=(3.*rhovyp+rhovy)/4.
      rhovz=(3.*rhovzp+rhovz)/4.
   elseif (rk_step.eq.3) then
      rhovx=(rhovxp+2.*rhovx)/3.
      rhovy=(rhovyp+2.*rhovy)/3.
      rhovz=(rhovzp+2.*rhovz)/3.
   endif
   vx=rhovx/rho
   vy=rhovy/rho
   vz=rhovz/rho

   call advance_bijk
   if (rk_step.eq.2) then
      bi=(3.*bip+bi)/4.
      bj=(3.*bjp+bj)/4.
      bk=(3.*bkp+bk)/4.
   elseif (rk_step.eq.3) then
      bi=(bip+2.*bi)/3.
      bj=(bjp+2.*bj)/3.
      bk=(bkp+2.*bk)/3.
   endif

   call fill_hydro_mpi
   call fill_mag_mpi

   call get_bxyz

   enddo ! enddo rk_step

end subroutine step_rk3 ! v1
!*********************************************
subroutine advance_previous
   use pars
   use fields
   implicit none

   rhop=rho
   bip=bi
   bjp=bj
   bkp=bk

   if (stepop.eq.1) then
      vxp=vx
      vyp=vy
      vzp=vz
      pp=p
      bxp=bx
      byp=by
      bzp=bz
   else
      rhovxp=rho*vx
      rhovyp=rho*vy
      rhovzp=rho*vz
      engp=.5*rho*(vx**2+vy**2+vz**2)+p/(gam-1.)
   endif

end subroutine advance_previous
!*********************************************
subroutine advance_predictor_ab
   use pars
   use fields
   implicit none

   rhoh=rho+tau/tau0/2.*(rho-rhop)
   vxh=vx+tau/tau0/2.*(vx-vxp)
   vyh=vy+tau/tau0/2.*(vy-vyp)
   vzh=vz+tau/tau0/2.*(vz-vzp)
   ph=p+tau/tau0/2.*(p-pp)
   bxh=bx+tau/tau0/2.*(bx-bxp)
   byh=by+tau/tau0/2.*(by-byp)
   bzh=bz+tau/tau0/2.*(bz-bzp)
   bih=bi+tau/tau0/2.*(bi-bip)
   bjh=bj+tau/tau0/2.*(bj-bjp)
   bkh=bk+tau/tau0/2.*(bk-bkp)

   call advance_previous

end subroutine advance_predictor_ab
!*********************************************
subroutine advance_corrector1
   use pars
   use fields
   implicit none
   integer :: i,j,k

   forall(i=1:nx0,j=1:ny0,k=1:nz0)
   rho(ng+i,ng+j,ng+k)=rho0(ng+i,ng+j,ng+k)-tau*( &
      (flux_rhox(i+1,j,k)-flux_rhox(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (flux_rhoy(i,j+1,k)-flux_rhoy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (flux_rhoz(i,j,k+1)-flux_rhoz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   eng(ng+i,ng+j,ng+k)=eng0(ng+i,ng+j,ng+k)-tau*( &
      (flux_engx(i+1,j,k)-flux_engx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (flux_engy(i,j+1,k)-flux_engy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (flux_engz(i,j,k+1)-flux_engz(i,j,k))/dzblock(ng+i,ng+j,ng+k))

   rhovx(ng+i,ng+j,ng+k)=rhovx0(ng+i,ng+j,ng+k)-tau*( &
      (flux_rvxx(i+1,j,k)-flux_rvxx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (flux_rvxy(i,j+1,k)-flux_rvxy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (flux_rvxz(i,j,k+1)-flux_rvxz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   rhovy(ng+i,ng+j,ng+k)=rhovy0(ng+i,ng+j,ng+k)-tau*( &
      (flux_rvyx(i+1,j,k)-flux_rvyx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (flux_rvyy(i,j+1,k)-flux_rvyy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (flux_rvyz(i,j,k+1)-flux_rvyz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   rhovz(ng+i,ng+j,ng+k)=rhovz0(ng+i,ng+j,ng+k)-tau*( &
      (flux_rvzx(i+1,j,k)-flux_rvzx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (flux_rvzy(i,j+1,k)-flux_rvzy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (flux_rvzz(i,j,k+1)-flux_rvzz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   endforall

!   vx=rhovx/rho
!   vy=rhovy/rho
!   vz=rhovz/rho
!   p=(gam-1.)*(eng-.5*rho*(vx**2+vy**2+vz**2))
   if(stepop.eq.1) p=(gam-1.)*(eng-.5*(rhovx**2+rhovy**2+rhovz**2)/rho)

end subroutine advance_corrector1
!*********************************************
subroutine advance_corrector2
   use pars
   use fields
   implicit none
   integer :: i,j,k

   forall(i=1:nx0,j=1:ny0,k=1:nz0)
   rhovx(ng+i,ng+j,ng+k)=rhovx(ng+i,ng+j,ng+k)-tau*( &
      (stress_xx(i+1,j,k)-stress_xx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (stress_xy(i,j+1,k)-stress_xy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (stress_xz(i,j,k+1)-stress_xz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   rhovy(ng+i,ng+j,ng+k)=rhovy(ng+i,ng+j,ng+k)-tau*( &
      (stress_yx(i+1,j,k)-stress_yx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (stress_yy(i,j+1,k)-stress_yy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (stress_yz(i,j,k+1)-stress_yz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   rhovz(ng+i,ng+j,ng+k)=rhovz(ng+i,ng+j,ng+k)-tau*( &
      (stress_zx(i+1,j,k)-stress_zx(i,j,k))/dxblock(ng+i,ng+j,ng+k) + &
      (stress_zy(i,j+1,k)-stress_zy(i,j,k))/dyblock(ng+i,ng+j,ng+k) + &
      (stress_zz(i,j,k+1)-stress_zz(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   endforall

   if (stepop.eq.1) then
      vx=rhovx/rho
      vy=rhovy/rho
      vz=rhovz/rho
   endif

end subroutine advance_corrector2
!*********************************************
subroutine advance_bijk
   use pars
   use fields
   implicit none
   integer :: i,j,k

   forall(i=1:nx0+1,j=1:ny0,k=1:nz0)
   bi(ng+i,ng+j,ng+k)=bi(ng+i,ng+j,ng+k)-tau*( &
                      (ek(i,j+1,k)-ek(i,j,k))/dyblock(ng+i,ng+j,ng+k) - &
                      (ej(i,j,k+1)-ej(i,j,k))/dzblock(ng+i,ng+j,ng+k))
   endforall

   forall(i=1:nx0,j=1:ny0+1,k=1:nz0)
   bj(ng+i,ng+j,ng+k)=bj(ng+i,ng+j,ng+k)-tau*( &
                      (ei(i,j,k+1)-ei(i,j,k))/dzblock(ng+i,ng+j,ng+k) - &
                      (ek(i+1,j,k)-ek(i,j,k))/dxblock(ng+i,ng+j,ng+k))
   endforall

   forall(i=1:nx0,j=1:ny0,k=1:nz0+1)
   bk(ng+i,ng+j,ng+k)=bk(ng+i,ng+j,ng+k)-tau*( &
                      (ej(i+1,j,k)-ej(i,j,k))/dxblock(ng+i,ng+j,ng+k) - &
                      (ei(i,j+1,k)-ei(i,j,k))/dyblock(ng+i,ng+j,ng+k))
   endforall

end subroutine advance_bijk
!*********************************************
subroutine termination
   use pars
   use fields
   implicit none
   
   call dealloc_fields

end subroutine termination
!*********************************************
