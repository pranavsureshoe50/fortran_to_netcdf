!! This is an example program which writes some 4D pressure and
! temperatures. It is intended to illustrate the use of the netCDF
! fortran 90 API. The companion program pres_temp_4D_rd.f shows how
! to read the netCDF data file created by this program.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: pres_temp_4D_wr.f90,v 1.7 2007/01/24 19:32:10 russ Exp $

program pres_temp_4D_wr
  use netcdf
  implicit none
  
  integer ( kind = 4 ), parameter :: nx = 432
  integer ( kind = 4 ), parameter :: ny = 1026
  integer ( kind = 4 ), parameter :: nz = 1282
  integer ( kind = 4 ), parameter :: nt = 1
  integer, parameter :: NDIMS = 4, NRECS = 1

  CHARACTER(LEN=120), PARAMETER ::  time_res1 ='00049891.res' ! '00075953.res'
  
  !-------------------Input parameters -------------------------------------------------------------------

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx),REC_VAL(NRECS)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,dir,kstart,kend,s1
  real    ( kind = 8 ) :: time0,time1,time2
  character(len = 512)   :: filename

  real    ( kind = 8 ) :: var(nx,ny,nz), Ucen(nx,ny,nz) !,vtmp1(nx,ny,nz),wtmp1(nx,ny,nz)
 
  ! This is the name of the data file we will create.
  integer :: ncid

  ! We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
  ! timesteps of data.
  character (len = *), parameter :: XVAL = "xval"
  character (len = *), parameter :: YVAL = "yval"
  character (len = *), parameter :: ZVAL = "zval"
  character (len = *), parameter :: REC_NAME = "time"
  integer :: xdimid, zdimid, ydimid, rec_dimid

  ! The start and count arrays will tell the netCDF library where to
  ! write our data.
  integer :: start(NDIMS), count(NDIMS)

  ! These program variables hold the latitudes and longitudes.
!  real :: lats(ny), lons(nz)
  integer :: zvarid, yvarid, xvarid, rec_varid

  ! We will create two netCDF variables, one each for temperature and
  ! pressure fields.
  character (len = 128) :: PRES_NAME,buff
  character (len = *), parameter :: TEMP_NAME="temperature"
  integer :: pres_varid, temp_varid
  integer :: dimids(NDIMS)

  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PRES_UNITS = "SI"
  character (len = *), parameter :: TEMP_UNITS = "celsius"
  character (len = *), parameter :: yUNITS = "span"
  character (len = *), parameter :: zUNITS = "stream"
  character (len = *), parameter :: xUNITS = "vertical"
  character (len = *), parameter :: REC_UNITS = "time"

  ! Program variables to hold the data we will write out. We will only
  ! need enough space to hold one timestep of data; one record.
  real :: pres_out(nz, ny, nx)
  real :: temp_out(nz, ny, nx)


  ! Loop indices
  integer :: nzg, tag, lon, rec
 
!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  dir = 1 ! -1 for dens, p

  IF(COMMAND_ARGUMENT_COUNT().NE.2)THEN
    WRITE(*,*)'ERROR, INPUT FILE NAME, dir AS COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
    STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,filename)   !first, read in the two values
  CALL GET_COMMAND_ARGUMENT(2,buff)   !first, read in the two values
  
!  filename = './u_' // TRIM(time_res1)

  filename = TRIM(filename)

  READ(buff,*) dir


  SELECT CASE (dir)
    CASE (1)
            PRES_NAME='VERTICAL VELOCITY'
    CASE (2)
            PRES_NAME='SPANWISE VELOCITY'
    CASE (3)
            PRES_NAME='STREAMWISE VELOCITY'
    CASE (4)
            PRES_NAME='PRESSURE'
    CASE (5)
            PRES_NAME='DENS'
  END SELECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  nzg=nz

  call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg,ru,rp,tag)

  call read_restart(filename,nx,ny,nz,var,time0)

  if (dir .lt. 4)  call center_velocity(nx,ny,nz,Ucen,var,dir)

  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.
  do k = 1, nx
     do j = 1, ny
        do i = 1, nz
           pres_out(i, j, k) = Ucen(k, j, i)
        end do
     end do
  end do
  REC_VAL(NRECS) = time0


  filename = trim(filename) // '.nc'
  ! Create the file. 
  call check( nf90_create(filename, nf90_clobber, ncid) )
  
  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  call check( nf90_def_dim(ncid, XVAL, nx, xdimid) )
  call check( nf90_def_dim(ncid, YVAL, ny, ydimid) )
  call check( nf90_def_dim(ncid, ZVAL, nz, zdimid) )
  call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

  ! Define the coordinate variables. We will only define coordinate
  ! variables for lat and lon.  Ordinarily we would need to provide
  ! an array of dimension IDs for each variable's dimensions, but
  ! since coordinate variables only have one dimension, we can
  ! simply provide the address of that dimension ID (ydimid) and
  ! similarly for (zdimid).
  call check( nf90_def_var(ncid, XVAL, NF90_REAL, xdimid, xvarid) )
  call check( nf90_def_var(ncid, YVAL, NF90_REAL, ydimid, yvarid) )
  call check( nf90_def_var(ncid, ZVAL, NF90_REAL, zdimid, zvarid) )
  call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, rec_dimid, rec_varid) )

  ! Assign units attributes to coordinate variables.
  call check( nf90_put_att(ncid, xvarid, UNITS, xUNITS) )
  call check( nf90_put_att(ncid, yvarid, UNITS, yUNITS) )
  call check( nf90_put_att(ncid, zvarid, UNITS, zUNITS) )
  call check( nf90_put_att(ncid, rec_varid, UNITS, REC_UNITS) )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids = (/ zdimid, ydimid, xdimid, rec_dimid /)

  ! Define the netCDF variables for the pressure and temperature data.
  call check( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )
!  call check( nf90_def_var(ncid, TEMP_NAME, NF90_REAL, dimids, temp_varid) )

  ! Assign units attributes to the netCDF variables.
  call check( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
!  call check( nf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) )
  
  ! End define mode.
  call check( nf90_enddef(ncid) )
  
  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  call check( nf90_put_var(ncid, rec_varid, REC_VAL) )
  call check( nf90_put_var(ncid, xvarid, xc) )
  call check( nf90_put_var(ncid, yvarid, yc) )
  call check( nf90_put_var(ncid, zvarid, zc) )
  
  ! These settings tell netcdf to write one timestep of data. (The
  ! setting of start(4) inside the loop below tells netCDF which
  ! timestep to write.)
  count = (/ nz, ny, nx, NRECS /)
  start = (/ 1, 1, 1, 1 /)

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays only hold one timestep worth
  ! of data. We will just rewrite the same data for each timestep. In
  ! a real :: application, the data would change between timesteps.
  do rec = 1, NRECS
     start(4) = rec
     call check( nf90_put_var(ncid, pres_varid, pres_out, start = start, &
                              count = count) )
!     call check( nf90_put_var(ncid, temp_varid, temp_out, start = start, &
!                              count = count) )
  end do
  
  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are really written to disk.
  call check( nf90_close(ncid) )
  
  print *,"*** SUCCESS writing example file ", filename, "!"
!  nccopy -d9 -s u_4D.nc tmp.nc
contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program pres_temp_4D_wr








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@@@@ READ GRID ##########################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg,ru,rp,tag)
  implicit none

  INTEGER nx,ny,nz,nzg,tag
  REAL ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nzg)
  REAL ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nzg)
  REAL ( kind = 8 ) :: ru(nx),rp(nx)

  INTEGER i,j,k

  real rdelx,rdely,rdelz, dtheta
  real, allocatable, dimension(:) :: cug,cvg

  ALLOCATE(cug(nzg),cvg(nzg))

  ! ! READ GRID

  !OPEN(UNIT=1,FILE="x1_grid_for_masoud_conrad.in",STATUS='OLD',FORM='FORMATTED')
  OPEN(UNIT=1,FILE="./x1_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nx-1
     read(1,*) j, xu(i)
     !write(6,*) "xu(",i,") = ", xu(i)
  enddo
  close(1)
  xc(2:nx-1) = .5*(xu(1:nx-2)+xu(2:nx-1))
  xc(1 ) = 2.*xu(1  )-xc(2  )
  xc(nx) = 2.*xu(nx-1)-xc(nx-1)

  do i= 2, nx-1
!     write(6,*) "xc(",i,") = ", xc(i)
  enddo


  !OPEN(UNIT=1,FILE="x2_grid_for_masoud_conrad.in",STATUS='OLD',FORM='FORMATTED')
  OPEN(UNIT=1,FILE="./x2_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, ny-1
     read(1,*) j, yv(i)
     !write(6,*) "yv(",i,") = ", yv(i)
  enddo
  close(1)
  yc(2:ny-1) = .5*(yv(1:ny-2)+yv(2:ny-1))
  yc(1 ) = 2.*yv(1  )-yc(2  )
  yc(ny) = 2.*yv(ny-1)-yc(ny-1)

  do i= 2, ny-1
!     write(6,*) "yc(",i,") = ", yc(i)
  enddo


  !OPEN(UNIT=1,FILE="x3_grid_for_masoud_conrad.in",STATUS='OLD',FORM='FORMATTED')
  OPEN(UNIT=1,FILE="./x3_grid.in",STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
     !write(6,*) "zwg(",i,") = ", zwg(i)
  enddo
  close(1)
  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  )= 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

!     write(6,*) "zcg(",i,") = ", zcg(i)


  ! ru(1:nx)=xu(1:nx)
  ! rp(1:nx)=xc(1:nx)

  close(1)

  write(6,*) "READ GRID DONE"

  return
end subroutine read_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@@@@ READ RESTART##########################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_restart(filename,nx,ny,nz,var,time)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav
  var = 0.0
  write(6,*) "nx,ny,nz = ", nx,ny,nz

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  READ(19) I,J,K,JP 
!  write(6,*) "I,J,K,JP = ", I,J,K,JP 
  DO K=1,NZ
!     write(6,*) " READ K = ", K
     READ(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  READ(19) nstep
  READ(19) TIME
  write(6,*) 'time=',time
  READ(19) DTM1,grav
  CLOSE(19)
  write(6,*) "READING RESTART FILE DONE"
  
  return
end subroutine read_restart





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  center velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine center_velocity(nx,ny,nz,Ucen,Uin,dir)
  !Passed Variables
  integer,intent(in)      :: dir,nx,ny,nz
  real (kind = 8),intent(in)         :: Uin(nx,ny,nz) 
  real (kind = 8),intent(out)        :: Ucen(nx,ny,nz) 


  !Local Variables
  integer              	:: i,j,k,err
  !Zero Output Array
  Ucen=0.d0
  !*************************************************
  !********************X1***************************
  !*************************************************
  if(dir.EQ.1) then
     !U
     do k=2,nz-1
!        write(6,*) " CENTER_VELOCITY DIR1, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i-1,j,k))
           enddo
     	enddo
     enddo

     !*************************************************
     !********************X2***************************
     !*************************************************
  elseif (dir.EQ.2) then 
     !V
     !U
     do k=2,nz-1
!        write(6,*) " CENTER_VELOCITY DIR2, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j-1,k))
           enddo
   	enddo
     enddo
     !*************************************************
     !********************X3***************************
     !*************************************************
  elseif (dir.EQ.3) then 
     !W
     !U
     do k=2,nz-1
!        write(6,*) " CENTER_VELOCITY DIR3, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j,k-1))
           enddo
   	enddo
     enddo

  else
     !Invalid direction
     !write(*,'(a60,i2)') "INVALID DIRECTION IN center_velocities 
     !&                      dir must be 1,2,3.  dir= ", dir
  endif

  Ucen(:,1,:) = Ucen(:,ny-1,:)
  Ucen(:,ny,:) = Ucen(:,2,:)

  return
end subroutine center_velocity



