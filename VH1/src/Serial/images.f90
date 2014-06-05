subroutine images

! Produce a 3D file containing 2D maps of the density over time.
! Data is stored in floating point format using netcdf.  The third 
! dimension uses the unlimited dimension capability of netcdf to 
! grow the third dimension each time a new time step is output to the movie.
!=================================================================================
use global
use zone
!use netcdf

IMPLICIT NONE

! LOCALS
CHARACTER(LEN=8) :: coord
CHARACTER(LEN=50):: filename
INTEGER :: nc_stat, ncid, frame_number, rank, nvar, natt, icheck, jcheck
INTEGER :: xDimID, yDimID, tDimID, varID, XScale_varID, YScale_varID, TScale_varID
INTEGER, DIMENSION(3) :: start, count

return
end


