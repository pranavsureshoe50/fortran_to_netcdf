Program Main
  implicit none
!  include 'header'
  
  integer ( kind = 4 ), parameter :: nx = 402
  integer ( kind = 4 ), parameter :: ny = 1026
  integer ( kind = 4 ), parameter :: nz = 1538
  CHARACTER(LEN=120), PARAMETER ::  time_res1 = '00075953.res'
  CHARACTER(LEN=120), PARAMETER ::  time_res2 = '00075953.res'
  
  integer ( kind = 4 ), parameter :: nx_start = 1 
  integer ( kind = 4 ), parameter :: nx_end = 250 

  integer ( kind = 4 ), parameter :: ny_start = 2 
  integer ( kind = 4 ), parameter :: ny_end = ny-1 
          
  integer ( kind = 4 ), parameter :: nz_start = 2 
  integer ( kind = 4 ), parameter :: nz_end =   nz-1

  real ( kind = 8 ), parameter :: ru1 = 0.001667 
 !-------------------Input parameters -------------------------------------------------------------------
 logical :: full_box,comp_dwdt

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,tmp,kstart,kend,s1
  integer ( kind = 4 ) :: tag,nx_write,ny_write,nz_write
  real    ( kind = 8 ) :: time0,time1,time2,DTM1,grav,Re
  character(len=128)   :: buff,filename

  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: u1,v1,w1,p1,dens1
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: u2,v2,w2,omg_x,omg_y,omg_z
  real    ( kind = 8 ),allocatable, dimension(:,:,:) :: advect,stretch,diffuse
  real    ( kind = 8 ) :: utmp1(nx,ny,nz),vtmp1(nx,ny,nz),wtmp1(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x2(nx,ny,nz),omg_y2(nx,ny,nz),omg_z2(nx,ny,nz)
  real    ( kind = 8 ) ::  var(nx,ny,nz)
  real    ( kind = 8 ) ::  sgs(nx,ny,nz), dens(nx,ny,nz)
  real    ( kind = 8 ) :: visc_omg_ver(nx,ny,nz)
  logical :: save_advect , save_diffuse , save_stretch

 save_advect = .TRUE.
 save_diffuse = .TRUE.
 save_stretch = .TRUE.
 full_box = .FALSE.
 comp_dwdt = .FALSE.

 omg_x2 = 0.0 ; omg_y2 = 0.0 ; omg_z2 = 0.0

  write(*,*) 'reading grid ..'
  call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nz,ru,rp,tag)

  allocate(w2(nx,ny,nz),u2(nx,ny,nz),v2(nx,ny,nz),dens(nx,ny,nz))!,p1(nx,ny,nz),dens1(nx,ny,nz))
  allocate(advect(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end))
  allocate(stretch(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)) 
  allocate(diffuse(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end))

!!!!!!!!!!!!!!!!!!!!-------------------------------------------------
! COMPUTE DWDT -------------------------------------------------------
  if (comp_dwdt) then 
   
    allocate(omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz))
    allocate(w1(nx,ny,nz),u1(nx,ny,nz),v1(nx,ny,nz))!,p1(nx,ny,nz),dens1(nx,ny,nz))
    filename = './u_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,u1,time1)
    filename = './v_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,v1,time1)
    filename = './w_' // TRIM(time_res1)
    call read_restart(filename,nx,ny,nz,w1,time1)
 
    write(*,*) 'finished reading  restart files'

! Calculated the 3 comp of omega with first restart files time 
    call vorticity_3comps_cartesian(u1,v1,w1,omg_x,omg_y,omg_z,nx,ny,nz,xu,yv,zwg)
    write(6,*) "Calculation"
    deallocate(w1,u1,v1)
  endif 
! DONE COMPUTE DWDT ------------------------------------------------------- 
!----------------------------------------------------------------------------------

  filename = './u_' // TRIM(time_res2)
  call read_restart(filename,nx,ny,nz,u2,time2)
  filename = './v_' // TRIM(time_res2)
  call read_restart(filename,nx,ny,nz,v2,time2)
  filename = './w_' // TRIM(time_res2)
  call read_restart(filename,nx,ny,nz,w2,time2)
  filename = './tv_' // TRIM(time_res2)
  call read_restart(filename,nx,ny,nz,sgs,time2)

  filename = './dens_' // TRIM(time_res2)
  call read_restart(filename,nx,ny,nz,dens,time2)

  write(*,*) 'finished reading  restart files second time'
  
! ----------------  COMPUTE THE THREE VORTICITIES ----------------------------------------------------
! Computed  vorticity 
  call vorticity_3comps_cartesian(u2,v2,w2,omg_x2,omg_y2,omg_z2,nx,ny,nz,xu,yv,zwg)
  write(*,*) 'finished computing vorticity'

! ----------------------- CENTER VELOCITIES ----------------------------------------------------------
  call center_velocity(nx,ny,nz,utmp1,u2,1)
  call center_velocity(nx,ny,nz,vtmp1,v2,2)
  call center_velocity(nx,ny,nz,wtmp1,w2,3)

! --------------- COMPUTE ADVECTION TERMS -------------------------------------------------------------
! Compute advection - u.del(omgx)   
! x terms  
  call advection_omega(u2,v2,w2,omg_x2,var,nx,ny,nz,xc,yc,zcg)
  advect(:,:,:) = var(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  filename = 'advect_box_x.plt'
  call write_plt_3d_1var(filename,advect,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  
! y terms
  call advection_omega(u2,v2,w2,omg_y2,var,nx,ny,nz,xc,yc,zcg)
  advect(:,:,:) = var(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  filename = 'advect_box_y.plt'
  call write_plt_3d_1var(filename,advect,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  
! z terms 
  call advection_omega(u2,v2,w2,omg_z2,var,nx,ny,nz,xc,yc,zcg)
  advect(:,:,:) = var(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  filename = 'advect_box_z.plt'
  call write_plt_3d_1var(filename,advect,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
   
  
  write(*,*) 'finished computing advection'
!   do k=nz_start,nz_end
!      do j=ny_start,ny_end
!      do i=nx_start,nx_end
!         advect(i,j,k) = var(i,j,k)
!      enddo
!      enddo
!   enddo

  
 
! --------------- COMPUTE STRETCHING TERMS -------------------------------------------------------------
! Compute advection - omg.del(u)   
! x terms
  call stretching_omega(omg_x2,omg_y2,omg_z2,utmp1,nx,ny,nz,xc,xu,yc,yv,zcg,zwg,var)
  stretch(:,:,:) = var(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'stretch_box_x.plt'
  call write_plt_3d_1var(filename,stretch,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  
! y terms
  call stretching_omega(omg_x2,omg_y2,omg_z2,vtmp1,nx,ny,nz,xc,xu,yc,yv,zcg,zwg,var)
  stretch(:,:,:) = var(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'stretch_box_y.plt'
  call write_plt_3d_1var(filename,stretch,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
 
! z terms
  call stretching_omega(omg_x2,omg_y2,omg_z2,wtmp1,nx,ny,nz,xc,xu,yc,yv,zcg,zwg,var)
  stretch(:,:,:) = var(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'stretch_box_z.plt'
  call write_plt_3d_1var(filename,stretch,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)

  write(*,*) 'finished computing stretching'



! --------------- COMPUTE VISCOUS TERMS -------------------------------------------------------------
!subroutine visc_omega(omg_x,var,xu,yv,zwg,xc,yc,zcg)
! compute gradsq(omg)
! x terms
  call visc_omega(omg_x2,visc_omg_ver,xu,yv,zwg,xc,yc,zcg,nx,ny,nz,ru1)
  diffuse(:,:,:) = visc_omg_ver(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'visc_box_x.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  
! y terms
  call visc_omega(omg_y2,visc_omg_ver,xu,yv,zwg,xc,yc,zcg,nx,ny,nz,ru1)
  diffuse(:,:,:) = visc_omg_ver(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'visc_box_y.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  
! x terms
  call visc_omega(omg_z2,visc_omg_ver,xu,yv,zwg,xc,yc,zcg,nx,ny,nz,ru1)
  diffuse(:,:,:) = visc_omg_ver(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'visc_box_z.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  
  write(*,*) 'finished computing diffusion of omega'

! Now I can use omg_x, omg_z and omg_y for temporary storage  
 


! --------------- COMPUTE BAROCLINIC TERMS -------------------------------------------------------------
! compute k x grad(rho)
! x terms
  call  baroclinic(omg_z2,omg_y2,dens,xc,xu,yc,yv,zc,zw,nx,ny,nz) 
  diffuse(:,:,:) = omg_z2(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'baro_box_z.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
 
  stretch(:,:,:) = omg_y2(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'baro_box_y.plt'
  call write_plt_3d_1var(filename,stretch,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
 


! --------------- COMPUTE EDDY-VORTICITY TERMS -------------------------------------------------------------
  call eddy_vorticity(sgs,u2,v2,w2,utmp1,vtmp1,wtmp1,xc,xu,yc,yv,zc,zw,nx,ny,nz,omg_x2,omg_y2,omg_z2)
! x terms 
  diffuse(:,:,:) = omg_x2(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'eddyf_box_x.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)

! y terms 
  diffuse(:,:,:) = omg_y2(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'eddyf_box_y.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
 
! z terms 
  diffuse(:,:,:) = omg_z2(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  filename = 'eddyf_box_z.plt'
  call write_plt_3d_1var(filename,diffuse,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)


 deallocate(advect,stretch,diffuse)
 deallocate(w2,u2,v2)
  !--------------------------------------------------------------------!

  stop
end Program Main



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

subroutine ascii_version

  integer (kind = 4) :: i,j,k,imax,jmax,kmax
  
  open(1,file="3d.vtk", form = 'formatted', status = 'unknown')

  imax = 10
  jmax = 10
  kmax = 10

  write(1,FMT='(a26)') "# vtk DataFile Version 2.0"
  write(1,FMT='(a9)') "FlowField"
  write(1,FMT='(a5)') "ASCII"
  write(1,FMT='(a23)') "DATASET STRUCTURED_GRID"
  write(1,FMT='(a10)', advance="no") "DIMENSIONS"
  write(1,"(a1, 3i7)") " ", imax, jmax, kmax
  write(1,FMT='(a6)', advance="no") "POINTS"
  write(1,"(a1, i15)", advance="no") " ", imax*jmax*kmax
  write(1,"(a6)") " float"

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1,*) real(i), real(j), real(k)
        end do
     end do
  end do

  close(1)

end subroutine ascii_version

subroutine binary_version

  integer (kind = 4) :: i,j,k,imax,jmax,kmax,s1
  character(len=128) :: buff
  real    (kind = 4), allocatable, dimension(:,:,:) :: vtkvar
  
  open(unit=1,file='3d.vtk',access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = 10
  jmax = 10
  kmax = 10

  allocate(vtkvar(imax,jmax,kmax))

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1) real(i), real(j), real(k)
           vtkvar(i,j,k) = real(i)*real(j)*real(k)
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS vtkvar float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1) vtkvar(i,j,k)
        end do
     end do
  end do


  close(1)
end subroutine binary_version

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





!
! VORTICITY TERMS
! ----------------------------------------------------------------------------
! input variables, u at u-location, v at v-location, w at w-location
! ! omg variable at centre location
! advection
subroutine advection_omega(u,v,w,omg_x,var,nx,ny,nz,xc,yc,zcg)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: var(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz)
  real    ( kind = 8 ) :: u3w3,u2w3,u1w3
  real    ( kind = 8 ) :: dq2x2,dq3x3
  real    ( kind = 8 ) :: udomgdx, vdomgdy,wdomgdz
  real    ( kind = 8 ) :: domg_next
  real    ( kind = 8 ) :: domg_prev

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)
  real    ( kind = 8 ) :: dx,dy,dz

  var = 0.0
  
  do k=3,nz-2
!     write(6,*) " OMG: K = ", K
     do j=3,ny-2
    do i=3,ny-2
           dx= xc(i) - xc(i-1)
           dy= yc(j) - yc(j-1)
           dz= zcg(k)- zcg(k-1)

           domg_prev = u(i-1,j,k)*(omg_x(i,j,k) - omg_x(i-1,j,k))/dx 
           domg_next = u(i,j,k)  *(omg_x(i+1,j,k) - omg_x(i,j,k))/(xc(i+1)-xc(i))
! udw_x/dx
           udomgdx = 0.5*(domg_next + domg_prev)
           
           domg_prev = v(i,j-1,k) *(omg_x(i,j,k) - omg_x(i,j-1,k))/dy 
           domg_next = v(i,j,k)   *(omg_x(i,j+1,k) - omg_x(i,j,k))/(yc(j+1)-yc(j))
! VDW_x/DY
           vdomgdy = 0.5*(domg_next + domg_prev)

  
           domg_prev = w(i,j,k-1) *(omg_x(i,j,k) - omg_x(i,j,k-1))/dz 
           domg_next = w(i,j,k)   *(omg_x(i,j,k+1) - omg_x(i,j,k))/(zcg(k+1)-zcg(k))
! wdw_x/dz
           wdomgdz = 0.5*(domg_next + domg_prev)

           var(i,j,k) = udomgdx + vdomgdy +  wdomgdz
         
        enddo
     enddo
  enddo

  return
end subroutine advection_omega
!----------------------------------------------------------------------
!____________END OF ADVECTION _________________________________________



! -------------------------------------------------------------------------------
! vortex stretching
! u -  centred variable
! omgx,omgt,omgz - centered variable
subroutine stretching_omega(omg_x,omg_y,omg_z,u,nx,ny,nz,xc,xu,yc,yv,zcg,zwg,var)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: u3w3,u2w3,u1w3
  real    ( kind = 8 ) :: dq2x2,dq3x3,var(nx,ny,nz)
  real    ( kind = 8 ) :: omgx_dudx, omgy_dudy_prev , omgz_dudz_prev, omgy_dudy_next, omgz_dudz_next
  real    ( kind = 8 ) :: omgy_v 
  real    ( kind = 8 ) :: omgz_w,up,un,dudx,dudy,dudz

  real    ( kind = 8 ) :: xu(nx),xc(nx),yv(ny),yc(ny),zwg(nz),zcg(nz)
  real    ( kind = 8 ) :: dx,dy,dz

  var = 0.0

  do k=3,nz-2
!     write(6,*) " OMG: K = ", K
     do j=3,ny-2
        do i=3,nx-2

           dx= xu(i) - xu(i-1)
           dy= yv(j) - yv(j-1)
           dz= zwg(k)- zwg(k-1)

           up = 0.5*(u(i,j,k) + u(i+1,j,k))
           un = 0.5*(u(i,j,k) + u(i-1,j,k))
           dudx = (up-un)/dx

           up = 0.5*(u(i,j,k) + u(i,j+1,k))
           un = 0.5*(u(i,j,k) + u(i,j-1,k))
           dudy = (up-un)/dy

           up = 0.5*(u(i,j,k) + u(i,j,k+1))
           un = 0.5*(u(i,j,k) + u(i,j,k-1))
           dudz = (up-un)/dz

           var(i,j,k) = omg_x(i,j,k)*dudx + omg_y(i,j,k)*dudy + omg_z(i,j,k)*dudz

        enddo
     enddo
  enddo

  return
end subroutine stretching_omega
!----------------------------------------------------------------------
!____________END OF STRETCHING _________________________________________


 ----------------------------------------------------------------------------------------------
 ! -----------------------VISCOUS TERMS -------------------
 ! Pranav , June 23
 ! omg - centred vorticity fields input
subroutine visc_omega(omg_x,var,xu,yv,zwg,xc,yc,zcg,nx,ny,nz,ru1)
implicit none
  
!  INCLUDE 'header' 

  integer ( kind = 4 ) :: nx,ny,nz
  real ( kind = 8 ) :: ru1
  real    ( kind = 8 ) :: xu(nx),yv(ny),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz)
  real    ( kind = 8 ) :: var(nx,ny,nz)
  real    (kind = 8)   :: domgdx_p, domgdx_m, domgdy_p , domgdy_m 
  real    (kind = 8)   :: domgdz_p , domgdz_m
  real    (kind = 8)   :: domgdxx, domgdyy , domgdzz,dx,dy,dz
  integer (kind = 4)   :: i,j,k  
  var =0.0

  do k=3,nz-2
!     write(6,*) " OMG: K = ", K
     do j=3,ny-2
        do i=3,nx-2

           dx= xc(i) - xc(i-1)
           dy= yc(j) - yc(j-1)
           dz= zcg(k)- zcg(k-1)

! at u
    domgdx_m = (omg_x(i,j,k) - omg_x(i-1,j,k))/dx

!at v
    domgdy_m = (omg_x(i,j,k) - omg_x(i,j-1,k))/dy

!at w
    domgdz_m = (omg_x(i,j,k) - omg_x(i,j,k-1))/dz      
 
    
           dx= xc(i+1) - xc(i)
           dy= yc(j+1) - yc(j)
           dz= zcg(k+1)- zcg(k)
! at u
    domgdx_p = (omg_x(i+1,j,k) - omg_x(i,j,k))/dx

!at v
    domgdy_p = (omg_x(i,j+1,k) - omg_x(i,j,k))/dy

!at w
    domgdz_p = (omg_x(i,j,k+1) - omg_x(i,j,k))/dz      
 

           dx= xu(i) - xu(i-1)
           dy= yv(j) - yv(j-1)
           dz= zwg(k)- zwg(k-1)


           domgdxx = (domgdx_p - domgdx_m)/dx

           domgdyy = (domgdy_p - domgdy_m)/dy

           domgdzz = (domgdz_p - domgdz_m)/dz


           var(i,j,k) = ru1*(domgdxx + domgdyy + domgdzz)
 

   enddo
   enddo
   enddo


return
end subroutine visc_omega
! ----------------------------------------------------------------------------------------------
! END OF VISCOUS TERMS COMPUTE


!---------------------------------------------------------------------------------------------------
!  Baroclinicity terms
! Pranav, June 23
subroutine baroclinic(rtmp1,rtmp2,dens,xc,xu,yc,yv,zc,zw,nx,ny,nz)
      implicit none
!Passed Variables
        integer,intent(in)      :: nx,ny,nz
        real,intent(in)         :: dens(nx,ny,nz)
        real,intent(in)         :: xu(nx),xc(nx),yv(ny),yc(ny),zw(nz),zc(nz)
        real,intent(out)        :: rtmp1(nx,ny,nz), rtmp2(nx,ny,nz)
    
        real ::  dx,dy,dz,dx1,dy1,dz1, densp, densn
    !Local Variables
        integer                 :: i,j,k
     
      do k=2,nz-1
        do j=2,ny-1
        do i=2,nx-1

        dx= xu(i) - xu(i-1)
        dy= yv(j) - yv(j-1)
        dz= zw(k) - zw(k-1)

        densp = 0.5*(dens(i,j,k) + dens(i,j,k+1))
        densn = 0.5*(dens(i,j,k) + dens(i,j,k-1))

         ! drho/dz (in code's coordinates)
        rtmp2(i,j,k) = -(densp-densn)/dz

        densp = 0.5*(dens(i,j,k) + dens(i,j+1,k))
        densn = 0.5*(dens(i,j,k) + dens(i,j-1,k))
! drho/dy (in code's coordinates)
        rtmp1(i,j,k) = (densp-densn)/dy

        enddo
      enddo
      enddo
end

!---------------------------------------------------------------------------------------------------
!---------------------------END OF BAROCLINIC TERMS----------------------------------------------------
!---------------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------
! EDDY VISCOSITY CONTRIBUTION
! June 24 - Pranav
        subroutine eddy_vorticity(tv,uo,vo,wo,u,v,w,xc,xu,yc,yv,zc,zw,nx,ny,nz,t1,t2,t3)
implicit none
! input variables, uo at u-location, vo at v-location, wo at w-location

!Passed Variables
 	integer,intent(in)      :: nx,ny,nz
 	real,intent(in)         :: tv(nx,ny,nz),uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz) 
 	real,intent(in)         :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz) 

        real,intent(in)         :: xu(nx),xc(nx),yv(ny),yc(ny),zw(nz),zc(nz) 
 	real,intent(out)        :: t1(nx,ny,nz) , t2(nx,ny,nz), t3(nx,ny,nz)
 	
        real :: dx,dy,dz,dx1,dy1,dz1,temp1(nx,ny,nz),temp2(nx,ny,nz),temp3(nx,ny,nz)
        real :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real :: dnutdx, dnutdy, dnutdz
        real :: dudxx, dudyy, dudzz, dvdxx, dvdyy, dvdzz, dwdxx, dwdyy, dwdzz
!Local Variables
 	integer              	:: i,j,k
!Zero Output Array
       
         do k=2,nz-1
         do j=2,ny-1
         do i=2,nz-1

         dx = xu(i) - xu(i-1)     
         dx1 = xc(i) - xc(i-1)     

         dy = yv(j) - yv(j-1) 

         dz = zwg(k) - zwg(k-1)
         dz1 = zcg(k) - zcg(k-1)

        dudx = (uo(i,j,k)-uo(i-1,j,k))/dx

        dudy = 0.5*( 0.5*(1.0/dy)*(uo(i  ,j+1, k )-uo(i  ,j-1, k )) 
     &             + 0.5*(1.0/dy)*(uo(i-1,j+1, k )-uo(i-1,j-1, k )) )

        dudz = 0.5*( 0.5*(1.0/dz1)*(uo(i  , j ,k+1)-uo(i  , j ,k-1))
     &             + 0.5*(1.0/dz1)*(uo(i-1, j ,k+1)-uo(i-1, j ,k-1)) )

        dvdx = 0.5*( 0.5*(1.0/dx1)*(vo(i+1,j  , k )-vo(i-1,j  , k ))
     &             + 0.5*(1.0/dx1)*(vo(i+1,j-1, k )-vo(i-1,j-1, k )) )


        dvdy = (vo(i,j+1,k)-vo(i,j,k))/dy

        dvdz = 0.5*( 0.5*(1.0/dz1)*(vo( i ,j  ,k+1)-vo( i ,j  ,k-1))
     &             + 0.5*(1.0/dz1)*(vo( i ,j-1,k+1)-vo( i ,j-1,k-1)) )

        dwdx = 0.5*( 0.5*(1.0/dx1)*(wo(i+1, j ,k  )-wo(i-1, j ,k  ))
     &             + 0.5*(1.0/dx1)*(wo(i+1, j ,k-1)-wo(i-1, j ,k-1)) )

        dwdy = 0.5*( 0.5*(1.0/dy)*(wo( i ,j+1,k  )-wo( i ,j-1,k  ))
     &             + 0.5*(1.0/dy)*(wo( i ,j+1,k-1)-wo( i ,j-1,k-1)) )

        dwdz = (wo(i,j,k)-wo(i,j,k-1))/dz

        dx = (xc(i+1) - xc(i-1))*0.5
        dy = (yc(j+1) - yc(j-1))*0.5
        dz = (zc(k+1) - zc(k-1))*0.5

        dudxx = ( u(i+1,j,k) + u(i-1,j,k) - 2*u(i,j,k) )/(dx*dx)
        dudyy = ( u(i,j+1,k) + u(i,j-1,k) - 2*u(i,j,k) )/(dy*dy)
        dudzz = ( u(i,j,k+1) + u(i,j,k-1) - 2*u(i,j,k) )/(dz*dz)


        dvdxx = ( v(i+1,j,k) + v(i-1,j,k) - 2*v(i,j,k) )/(dx*dx)
        dvdyy = ( v(i,j+1,k) + v(i,j-1,k) - 2*v(i,j,k) )/(dy*dy)
        dvdzz = ( v(i,j,k+1) + v(i,j,k-1) - 2*v(i,j,k) )/(dz*dz)


        dwdxx = ( w(i+1,j,k) + w(i-1,j,k) - 2*w(i,j,k) )/(dx*dx)
        dwdyy = ( w(i,j+1,k) + w(i,j-1,k) - 2*w(i,j,k) )/(dy*dy)
        dwdzz = ( w(i,j,k+1) + w(i,j,k-1) - 2*w(i,j,k) )/(dz*dz)

        dnutdx = (tv(i+1,j,k) - tv(i-1,j,k))/(2*dx)
        dnutdy = (tv(i,j+1,k) - tv(i,j-1,k))/(2*dy)
        dnutdz = (tv(i,j,k+1) - tv(i,j,k-1))/(2*dz)

        temp1(i,j,k) = dnutdx*dudx + dnutdy*dudy + dnutdz*dudz
     &        + tv(i,j,k)*( dudxx + dudyy + dudzz)   
        
        temp2(i,j,k) = dnutdx*dvdx + dnutdy*dvdy + dnutdz*dvdz
     &        + tv(i,j,k)*( dvdxx + dvdyy + dvdzz)   

        temp3(i,j,k) = dnutdx*dwdx + dnutdy*dwdy + dnutdz*dwdz
     &        + tv(i,j,k)*( dwdxx + dwdyy + dwdzz)   
       

        enddo
        enddo
        enddo

! Take curl of the terms cumputed before 
         do k=2,nz-1
         do j=2,ny-1
         do i=2,nz-1
 
           dx = (xc(i+1) - xc(i-1))
           dy = (yc(j+1) - yc(j-1))
           dz = (zc(k+1) - zc(k-1))

           t1(i,j,k) = ( temp3(i,j+1,k)-temp3(i,j-1,k) )/dy -
     &                 ( temp2(i,j,k+1)-temp2(i,j,k-1) )/dz

           t2(i,j,k) = ( temp1(i,j,k+1)-temp1(i,j,k-1) )/dz -
     &                 ( temp3(i+1,j,k)-temp3(i-1,j,k) )/dx

           t3(i,j,k) = ( temp2(i+1,j,k)-temp2(i-1,j,k) )/dx -
     &                  ( temp1(i,j+1,k)-temp1(i,j-1,k) )/dy             

        enddo
        enddo
        enddo

        end
! EDDY VISCOSITY TERMS
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------



subroutine vorticity_3comps_cartesian(u,v,w,omg_x,omg_y,omg_z,nx,ny,nz,xu,yv,zwg)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq2x2,dq3x3
  real    ( kind = 8 ) :: omgt,omgr,omgz

  real    ( kind = 8 ) :: xu(nx),yv(ny),zwg(nz)
  real    ( kind = 8 ) :: dx,dy,dz

  do k=2,nz-1
!     write(6,*) " OMG: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dx= xu(i) - xu(i-1)
           dy= yv(j) - yv(j-1)
           dz= zwg(k)- zwg(k-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dy
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dx

           dq2x1=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i-1,j,k)+v(i-1,j-1,k)))/dx
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dy

!           dq1x1(i,j,k)=(u(i,j,k)-u(i-1,j,k))/dx
                 
           omg_x(i,j,k) = dq3x2-dq2x3
           omg_y(i,j,k) = dq1x3-dq3x1
           omg_z(i,j,k) = dq2x1-dq1x2


        enddo
     enddo
  enddo

  return
end subroutine vorticity_3comps_cartesian


subroutine write_vtk_cartesian_3vars(filename,var1,var2,var3,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: var2(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: var3(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx-2
  jmax = ny-2
  kmax = nz-2

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, "float"
  write(1) buff//char(10)

  do k = 2,nz-1
     write(6,*) "GRID: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(xc(i)), real(yc(j)), real(zcg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS U float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var1(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS V float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var2(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS W float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var3(i,j,k))
        end do
     end do
  end do


  close(1)


end subroutine write_vtk_cartesian_3vars




subroutine write_plt_3d_3vars(filename,omg_x,omg_y,omg_z,dq1x1,dissipation,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x1(nx,ny,nz),dissipation(nx,ny,nz)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz),ru1

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,u_car,v_car,diss,w_car

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx-1
  jmax_plt = ny-1
  kmax_plt = nz-1

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  ALLOCATE(yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  ALLOCATE(zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt))
  allocate(diss(2:nx-1,2:ny-1,2:nz-1))



  ru1 = 0.001278
         
  
  kmod = 0
  DO k = 2, kmax_plt

     write(6,*) "sudo diss: k = ", k
     DO i = 2, imax_plt
        DO j = 2, jmax_plt

           xc_3d(i,j,k) = xc(i)!cos(yc(j))
           yc_3d(i,j,k) = yc(j)!sin(yc(j))
           zc_3d(i,j,k) = zcg(k)
	   diss(i,j,k) = (ru1)*(omg_x(i,j,k)*omg_x(i,j,k) + omg_y(i,j,k)*omg_y(i,j,k) +omg_z(i,j,k)*omg_z(i,j,k))


        ENDDO
     ENDDO
  ENDDO



  write(6,*) "done1"
  ier = tecini(trim(title)//nulchar,'x,y,z,omg_x,omg_y,omg_z'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt-1,jmax_plt-1,kmax_plt-1, &
              'BLOCK'//nulchar,nulchar)
  write(6,*) "done2"

  ! Write out the field data.
  itot = (imax_plt-1)*(jmax_plt-1)*(kmax_plt-1)
  ier = tecdat(itot,xc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,yc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,zc_3d(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,omg_x(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,omg_y(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  ier = tecdat(itot,omg_z(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  !ier = tecdat(itot,dq1x1(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  !ier = tecdat(itot,diss(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)
  !ier = tecdat(itot,dissipation(2:imax_plt,2:jmax_plt,2:kmax_plt),disdouble)


  ! Close the file
  ier = tecend()
  deallocate(xc_3d,yc_3d,zc_3d,diss)
  write(6,*) "done write_plt_3d_3vars"
  return 
  
end subroutine write_plt_3d_3vars
  

subroutine write_plt_3d_1var(filename,p,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end
  integer ( kind = 4 ) :: imax,jmax,kmax,s1,kstart,kend,nz,ny,nx

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: p(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,temp

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx_end-nx_start+1
  jmax_plt = ny_end-ny_start+1
  kmax_plt = nz_end-nz_start+1
 

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(yc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(zc_3d(imax_plt,jmax_plt,kmax_plt))

  ALLOCATE(temp(imax_plt,jmax_plt,kmax_plt))
      
  kmod = 0
  DO k = nz_start, nz_end
     DO i = nx_start, nx_end
        DO j = ny_start, ny_end

           xc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = xc(i)
           yc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = yc(j)
           zc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = zcg(k)
           
           temp(i+1-nx_start,j+1-ny_start,k+1-nz_start) = p(i,j,k)
!           write(*,*)  xc(i), yc(j) ,zcg(k) , p(i,j,k)
        ENDDO
     ENDDO
  ENDDO


  ier = tecini(trim(title)//nulchar,'x,y,z,var'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,jmax_plt,kmax_plt, &
              'BLOCK'//nulchar,nulchar)

  ! Write out the field data.
  itot = (imax_plt)*(jmax_plt)*(kmax_plt)
  ier = tecdat(itot,xc_3d,disdouble)
  ier = tecdat(itot,yc_3d,disdouble)
  ier = tecdat(itot,zc_3d,disdouble)
  ier = tecdat(itot,temp,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_1var"



  deallocate(xc_3d,yc_3d,zc_3d,temp)

  return 
end subroutine write_plt_3d_1var

subroutine write_3var_box(filename,u,v,w,xc,yc,zcg,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end,nx,ny,nz)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx_start,nx_end,ny_start,ny_end,nz_start,nz_end
  integer ( kind = 4 ) :: imax,jmax,kmax,s1,kstart,kend,nz,ny,nx

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: v(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)
  real    ( kind = 8 ) :: w(nx_start:nx_end,ny_start:ny_end,nz_start:nz_end)

  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,tempu,tempv,tempw

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx_end-nx_start+1
  jmax_plt = ny_end-ny_start+1
  kmax_plt = nz_end-nz_start+1
 

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(yc_3d(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(zc_3d(imax_plt,jmax_plt,kmax_plt))

  ALLOCATE(tempu(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tempv(imax_plt,jmax_plt,kmax_plt))
  ALLOCATE(tempw(imax_plt,jmax_plt,kmax_plt))
      
  kmod = 0
  DO k = nz_start, nz_end
     DO i = nx_start, nx_end
        DO j = ny_start, ny_end

           xc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = xc(i)
           yc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = yc(j)
           zc_3d(i+1-nx_start,j+1-ny_start,k+1-nz_start) = zcg(k)
           
           tempu(i+1-nx_start,j+1-ny_start,k+1-nz_start) = u(i,j,k)
           tempv(i+1-nx_start,j+1-ny_start,k+1-nz_start) = v(i,j,k)
           tempw(i+1-nx_start,j+1-ny_start,k+1-nz_start) = w(i,j,k)
!           write(*,*)  xc(i), yc(j) ,zcg(k) , p(i,j,k)
        ENDDO
     ENDDO
  ENDDO


  ier = tecini(trim(title)//nulchar,'x,y,z,omgx,omgy,omgz'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,jmax_plt,kmax_plt, &
              'BLOCK'//nulchar,nulchar)

  ! Write out the field data.
  itot = (imax_plt)*(jmax_plt)*(kmax_plt)
  ier = tecdat(itot,xc_3d,disdouble)
  ier = tecdat(itot,yc_3d,disdouble)
  ier = tecdat(itot,zc_3d,disdouble)
  ier = tecdat(itot,tempu,disdouble)
  ier = tecdat(itot,tempv,disdouble)
  ier = tecdat(itot,tempw,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_box"



  deallocate(xc_3d,yc_3d,zc_3d,tempu,tempv,tempw)

  return 
end subroutine write_3var_box

  
