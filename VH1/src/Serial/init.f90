subroutine init

! Pressure Problem in 1D
! 17jun14 cerba cmnunn
!=======================================================================
! GLOBALS
use global
use zone

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k
REAL :: ridt, xvel, yvel, zvel, width, widthz, widthy
REAL :: xmin, xmax, ymin, ymax, zmin, zmax
REAL :: dleft, pleft, dright, pright, plane
REAL :: a2, ve2b2a2                   !new, user defined variables

!--------------------------------------------------------------------------------
! Set up geometry and boundary conditions of grid
!
! Boundary condition flags : nleft, nright
!   = 0  :  reflecting boundary condition
!   = 1  :  inflow/outflow boundary condition (zero gradients)
!   = 2  :  fixed inflow boundary condition (values set by dinflo, pinflo, etc.)
!   = 3  :  periodic (nmax+1 = nmin; nmin-1 = nmax)
!
! Geometry flag : ngeom                         |  Cartesian:
!   = 0  :  planar                              |    gx = 0, gy = 0, gz = 0   (x,y,z)
!   = 1  :  cylindrical radial                  |  Cylindrical:
!   = 2  :  spherical   radial             3D= {     gx = 1, gy = 3, gz = 0   (s,phi,z)
!   = 3  :  cylindrical angle                   |
!   = 4  :  spherical polar angle (theta)       |  Spherical:
!   = 5  :  spherical azimu angle (phi)         |    gx = 2, gy = 4, gz = 5   (r,theta,phi)

! Define the computational grid...

ndim = 1
ngeomx = 2
ngeomy = 4
ngeomz = 5

nleftx = 2
nrightx = 2
nlefty = 3
nrighty = 3
nleftz = 3
nrightz = 3

pinflo = pleft
dinflo = dleft
uinflo = 0.0
vinflo = 0.0
winflo = 0.0

xmin   = 1.0*rsol
xmax   = 100.0*rsol
ymin   = 0.0
ymax   = 1.0
zmin   = 0.0
zmax   = 2.0

! If any dimension is angular, multiply coordinates by pi...
if(ngeomy >= 3) then
   ymin = ymin * pi
   ymax = ymax * pi
endif
if(ngeomz >= 3) then
   zmin = zmin * pi
   zmax = zmax * pi
endif

!======================================================================
! Set up parameters from the problem (Pressure Wave)

pright = smallp 
dright = smallr
pleft  = 0.138                    !kb*10**15
dleft  = 5.0e-16                    !10**-9/tcor
gam    = 1.01

gamm = gam - 1

!=======================================================================
! set time and cycle counters

time   = 0.0
timep  = 0.0
timem  = 0.0
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 0

! Set up grid coordinates

call grid(imax,xmin,xmax,zxa,zxc,zdx)
call grid(jmax,ymin,ymax,zya,zyc,zdy)
call grid(kmax,zmin,zmax,zza,zzc,zdz)

if (ndim <= 2) zzc(1) = 0.0
if (ndim == 1) zyc(1) = 0.0

!======================================================================
! Log parameters of problem in history file

write (8,"('Oblique Pressure Wave in ',i1,' dimensions.')") ndim
if (ndim==1) then
 write (8,"('Grid dimensions: ',i4)") imax
else if (ndim==2) then
 write (8,"('Grid dimensions: ',i4,' x ',i4)") imax, jmax
else
 write (8,"('Grid dimensions: ',i4,' x ',i4,' x ',i4)") imax, jmax, kmax
endif
write (8,*) 
write (8,"('Ratio of specific heats = ',f5.3)") gam
write (8,"('Pressure ratio is ',f5.3)") pright
write (8,"('Density ratio is ',f5.3)") dright
write (8,*) 

! initialize grid:

a2 = (kb*tcor)/mu
ve2b2a2 = (G*msol)/(2*rsol*a2)

do k = 1, kmax
 do j = 1, jmax
  do i = 1, imax
    !plane = zxc(i)+zyc(j)+zzc(k)
    !if(plane >= 0.5) then
     zpr(i,j,k) = pleft*exp(ve2b2a2*(rsol/zxa(i)-1))
     zro(i,j,k) = zpr(i,j,k)/a2
    !else
    !  zro(i,j,k) = dleft
    !  zpr(i,j,k) = pleft
    !endif
    zux(i,j,k) = 0.
    zuy(i,j,k) = 0.
    zuz(i,j,k) = 0.
    zfl(i,j,k) = 0.
  enddo
 enddo
enddo

!########################################################################
! Compute Courant-limited timestep

ridt = 0.
if(ndim == 1) then
  do i = 1, imax
    svel = sqrt(gam*zpr(i,1,1)/zro(i,1,1))/zdx(i)
    xvel = abs(zux(i,1,1)) / zdx(i)
    ridt = max(xvel,svel,ridt)
  enddo
else if(ndim == 2) then
  do j = 1, jmax
   do i = 1, imax
     widthy = zdy(j)
     if(ngeomy > 2) widthy = widthy*zxc(i)
     width  = min(zdx(i),widthy)
     svel = sqrt(gam*zpr(i,j,1)/zro(i,j,1))/width
     xvel = abs(zux(i,j,1)) / zdx(i)
     yvel = abs(zuy(i,j,1)) / widthy
     ridt = max(xvel,yvel,svel,ridt)
   enddo
  enddo
else if(ndim == 3) then 
  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      widthy = zdy(j)
      widthz = zdz(k)
      if(ngeomy > 2) widthy = widthy*zxc(i)
      if(ngeomz > 2) widthz = widthz*zxc(i)
      if(ngeomz == 5) widthz = widthz*sin(zyc(j))
      width  = min(zdx(i),widthy,widthz)
      svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
      xvel = abs(zux(i,j,k)) / zdx(i)
      yvel = abs(zuy(i,j,k)) / widthy
      zvel = abs(zuz(i,j,k)) / widthz
      ridt = max(xvel,yvel,zvel,svel,ridt)
    enddo
   enddo
  enddo
endif

dt = courant / ridt 

return
end

!#########################################################################

subroutine grid( nzones, xmin, xmax, xa, xc, dx )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------
IMPLICIT NONE

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax

!=======================================================================

dxfac = (xmax - xmin) / float(nzones)
do n = 1, nzones
  xa(n) = xmin + (n-1)*dxfac
  dx(n) = dxfac
  xc(n) = xa(n) + 0.5*dx(n)
enddo

return
end
