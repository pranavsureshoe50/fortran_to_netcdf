!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is a very simple example which writes a 2D array of
!     sample data. To handle this in netCDF we create two shared
!     dimensions, "x" and "y", and a netCDF variable, called "data".

!     This example demonstrates the netCDF Fortran 90 API. This is part
!     of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
      
!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_xy_wr.f90,v 1.7 2006/12/09 18:44:58 russ Exp $

program simple_xy_wr
  use netcdf
  implicit none

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "simple_xy.nc"

  ! We are writing 2D data, a 6 x 12 grid. 
  integer, parameter :: NDIMS = 2
  integer, parameter :: NX = 6, NY = 12

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  integer :: ncid, varid, dimids(NDIMS)
  integer :: x_dimid, y_dimid
  
  ! This is the data array we will write. It will just be filled with
  ! a progression of integers for this example.
  integer :: data_out(NY, NX)

  ! Loop indexes, and error handling.
  integer :: x, y

  ! Create some pretend data. If this wasn't an example program, we
  ! would have some real data to write, for example, model output.
  do x = 1, NX
     do y = 1, NY
        data_out(y, x) = (x - 1) * NY + (y - 1)
     end do
  end do

  ! Always check the return code of every netCDF function call. In
  ! this example program, wrapping netCDF calls with "call check()"
  ! makes sure that any return which is not equal to nf90_noerr (0)
  ! will print a netCDF error message and exit.

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "x", NX, x_dimid) )
  call check( nf90_def_dim(ncid, "y", NY, y_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ y_dimid, x_dimid /)

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer).
  call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.
  call check( nf90_put_var(ncid, varid, data_out) )

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

  print *, "*** SUCCESS writing example file simple_xy.nc! "

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program simple_xy_wr







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


