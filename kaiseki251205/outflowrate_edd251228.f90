!--- For Singleo
program Polar_to_Cartesian
implicit none

CHARACTER(4) filenum,filenum2,cn4,cn2
CHARACTER(4) dn,dnb
integer :: n,nb

!*****grid point of hydro code****
integer,parameter :: j0=140,k0=100,l0=1
real(8),dimension(1:j0) :: x1
real(8),dimension(1:k0) :: x2,dx2
real(8),dimension(1:k0+1) :: x2b
integer :: js,je,ks,ke,ls,le,m,mmax,jp,kp,j,k,l,jp_tmp,kp_tmp,f

real,dimension(:,:,:,:),allocatable:: w
real,dimension(:,:,:,:),allocatable :: wline   

!*.*.*. Cartesian.*.*.*
real :: xxmax,zzmax,dxx,dzz,xxmin,zzmin
integer :: ic,jc
!Cartesian total gred number
integer,parameter:: xxgrid=900,zzgrid=900
real(8),dimension(0:xxgrid-1,0:zzgrid-1) :: xx,zz,r,theta
!real(8),dimension(:,:,:),allocatable :: wc
real(8),dimension(0:xxgrid-1,0:zzgrid-1,2) :: p
real(8),dimension(0:xxgrid-1,0:zzgrid-1) :: V1,V2,V3,V4
real(8) ::V_tot

real(8), parameter :: mh = 1.67262158D-24
real(8), parameter :: v_light = 29979245800.d0
real(8), parameter :: kb = 1.3807D-16
real(8), parameter :: msun=1.989d33
real(8), parameter :: gc=6.67359d-8
real(8), parameter :: pi=3.141593d0
real(8), parameter :: sigma_e=0.4d0

real(8), parameter :: rho0=1.d-16
!real(8), parameter :: mbh=1.d8, edd_ratio=0.5d0
real(8) :: mbh, edd_ratio,z_off

real(8) :: dim_tem,rg,Ledd,Lx
real(8) :: v_esc,dm
!real(8),parameter :: z_off=9.d0
real(8) ::x_box,z_box,grav_r,grav_theta
real(8),dimension(:,:,:),allocatable ::r_wind,theta_wind
!*******
real(8) :: dx_ratio2,sum2,dx_initial2
real(8),dimension(0:k0) :: xm2
real(8) :: d_omega,omega, dm_dot, mdot, dp_dot, pdot,de_dot, edot,dl_dot,ldot,lm,Rl
real(8) :: mdot_old,pdot_old,edot_old,mdot_out,dm_dot_out,mdot_in,dm_dot_in
!**************
js=1
je=j0
ks=1
ke=k0
ls=1
le=l0
mmax=8
!**************


allocate(w(js:je,ks:ke,ls:le,mmax))
allocate(wline(js:je,ks:ke,ls:le,6))


!write(34,*)"mbh,edd_ratio,edd_ratio*Ledd, mdot*v_light*v_light/Ledd, mdot &
!     ,pdot*v_light/Ledd, pdot, edot/Ledd, edot"

do n=5,5,1
!do n=1,7,2
if (n.ge.1)then
   write(dn,'(i2.2)')n
else if (n.eq.0)then
   write(dn,'(i3.3)')5
else if (n.eq.-1)then
   write(dn,'(i4.4)')25
end if
!   write(dn,'(i2.2)')n
!   write(dn,'(i3.3)')n
   
   open (34,file='max_outflowrate-e'//trim(dn)//'.dat')
   open (33,file='min_outflowrate-e'//trim(dn)//'.dat')

do nb=8,8


   w=0.d0
   wline=0.d0

   write(dnb,'(i2)')nb
   dnb=adjustl(dnb)


   


!   open(7,file='/xc-work/nomuramr/2Ddense/151126parameter/M'//trim(dnb)//'/test'//trim(dn)//'/free_parameters.dat')
!   open(7,file='M'//trim(dnb)//'/test'//trim(dn)//'/free_parameters.dat')
   open(7,file='/Users/yokohama/Documents/work/test251205/free_parameters.dat')
   read(7,*) mbh
   read(7,*) edd_ratio
   close(7)

write(6,*)"parameters", mbh, edd_ratio

z_off=3.4d0**edd_ratio

!*******
rg=gc*mbh*msun/(v_light*v_light)
dim_tem=0.5d0*mh*v_light*v_light/kb
Ledd=4.d0*pi*v_light*gc*mbh*msun/sigma_e
!Lx=0.05d0*edd_ratio*Ledd

!allocate(wc(0:xxgrid-1,0:zzgrid-1,mmax))
!allocate(r_wind(js:je,ks:ke,ls:le))
!allocate(theta_wind(js:je,ks:ke,ls:le))

!###############grid################
!   open (unit=1,file='/xc-work/nomuramr/2Ddense/151126parameter/M'//trim(dnb)//'/'//trim(dn)//'/Data/grid-x1.data')
   open (unit=1,file='/Users/yokohama/Documents/work/test251205/Data/grid-x1.data')
   read(1,*) x1(1:j0)   
   close(1)
!   open (unit=2,file='/xc-work/nomuramr/2Ddense/151126parameter/M'//trim(dnb)//'/'//trim(dn)//'/Data/grid-x2.data')
   open (unit=2,file='/Users/yokohama/Documents/work/test251205/Data/grid-x2.data')
   read(2,*) x2(1:k0)   
   close(2)

!   do k=2,k0
!      x2b(k)=(x2(k)+x2(k-1))*0.5d0
!   end do
!   x2b(1)=0
!   x2b(k0+1)=4.d0*atan(1.d0)
   
!   do k=1,k0
!      dx2(k)=(x2b(k+1)-x2b(k))
!      write(6,*)"x2b",x2b(k),x2(k),dx2(k)
!   end do
   
mdot=0.d0
pdot=0.d0
edot=0.d0
ldot=0.d0
!each time step
!do f=0,1000,50
  !output step num
   w=0.d0
   wline=0.d0
!##########Read file average##########
   do m=1,8
      write(cn2,'(i2.2)')m
      
      if (m<=4.or.m==8) then
         
!         open(10,file='M'//trim(dnb)//'/test'//trim(dn)//'/w'//trim(cn2)//'-av.data' &
!              ,form='unformatted')
         open(10,file='/Users/yokohama/Documents/work/test251205/Data/w'//trim(cn2)//'-0830.data' &
              ,form='unformatted')
         read(10) w(js:je,ks:ke,ls:le,m)
         close(10)
         
         
      end if

      if (m<=6) then
         open(11,file='/Users/yokohama/Documents/work/test251205/Data/wline'//trim(cn2)//'-0830.data' &
              ,form='unformatted')
         read(11) wline(js:je,ks:ke,ls:le,m)
         
         close(11)      
      end if

   end do
!   write(6,*)"check",w(6,140,1,1)
   !do j=js,je
!do k=ks,ke
!      do l=ls,le
!         do m=1,8

!############################

!write(6,*)"check0419",w(100,0,1,1)

omega=0.d0
mdot=0.d0
edot=0.d0
pdot=0.d0
ldot=0.d0
!*************
  dx_ratio2=1.066d0
  xm2(0) = 0.d0
  xm2(k0) = 2.d0*atan(1.d0)
  sum2=2.d0*atan(1.d0)
  dx_initial2=sum2*(1.d0-(1.d0/dx_ratio2))/(1.d0-(1.d0/dx_ratio2)**dble(k0))
!*************

  open (unit=77,file="Data/out-r.dat")
  open (unit=66,file="Data/out-theta.dat")
  open (unit=67,file="Data/out-theta2.dat")

mdot_old=0.d0
edot_old=0.d0
pdot_old=0.d0


do j=1,j0

omega=0.d0
mdot=0.d0
edot=0.d0
pdot=0.d0
ldot=0.d0

mdot_out=0.d0
mdot_in=0.d0

write(6,*)"test-angle",x2(38)*45.d0/atan(1.d0),x2(39)*45.d0/atan(1.d0)
!pause
!stop
!do k=ks,50
do k=ks,k0
if (x2(k)*45.d0/atan(1.d0).lt.82.5d0) then
write(6,*)"k0926",k
do l=1,1
   
!pause
   dx2(k)=dx_initial2*dx_ratio2**(-dble(k-1))
   xm2(k)=xm2(k-1)+dx2(k)
!   if (x2(k)<xm2(k-1).or.x2(k)>xm2(k))then
!      write(6,*)"xm2_check",xm2(k-1),x2(k),xm2(k),dx2(k),(xm2(k-1)-xm2(k-2))/(xm2(k)-xm2(k-1))
!      stop
!   end if

!0521>>>>>>>>>>>>>>
!   x_box=x1(j)*sin(x2(k))
!   z_box=x1(j)*cos(x2(k))
   
!   r_wind(j,k,l)=sqrt(x_box**2.d0+(z_box+z_off)**2.d0)      
!   theta_wind(j,k,l)=atan(x_box/(z_box+z_off))
!0521<<<<<<<<<<<<<

!grav_r     =cos(x2(k)-theta_wind(j,k,l))/(r_wind(j,k,l)*r_wind(j,k,l))
!grav_theta =sin(x2(k)-theta_wind(j,k,l))/(r_wind(j,k,l)*r_wind(j,k,l))
!v_esc=sqrt(2.d0/r_wind(j,k,l))
w(j,k,l,5)=w(j,k,l,8)*dim_tem/w(j,k,l,1)  !temperature 

!dm=16.d0*atan(1.d0)*x1(j)*x1(j)*w(j,k,l,1)*w(j,k,l,2)*sin(x2(k))*dx2(k)&
!     *rg*rg*rho0*v_light

d_omega=4*pi*x1(j)*rg*x1(j)*rg*(cos(xm2(k-1))-cos(xm2(k)))
!if (w(j,k,l,2).gt.0.d0)then
   dm_dot=w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*d_omega 
    if (w(j,k,l,2)>=0) then
      dm_dot_out=w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*d_omega
      dm_dot_in=0
    else
      dm_dot_out=0
      dm_dot_in=w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*d_omega
    end if

      
   mdot_out=mdot_out+dm_dot_out
   mdot_in=mdot_in+dm_dot_in

   !dp_dot=dm_dot*abs(w(j,k,l,2))*v_light
   dl_dot=x1(j)*sin(x2(k))*rg*w(j,k,l,4)*v_light*dm_dot
   dp_dot=dm_dot*w(j,k,l,2)*v_light
   de_dot=dp_dot*w(j,k,l,2)*v_light*0.5d0
   
   !else
!   dm_dot=0.d0
!   !dp_dot=dm_dot*abs(w(j,k,l,2))*v_light
!   dp_dot=0.d0
!   de_dot=0.d0

!end if 




!write(66,*)x2(k)*45.d0/atan(1.d0),!

if (f==1000)then
   write(6,*)k,dm_dot
end if
if (d_omega<0.d0)then
   write(6,*)"error"
end if
!if (w(j,k,l,2)>0.d0) then
mdot=mdot+dm_dot
pdot=pdot+dp_dot
edot=edot+de_dot
ldot=ldot+dl_dot
!end if
omega=omega+d_omega

if (j.eq.j0) then
write(66,1000)x2(k)*45.d0/atan(1.d0),mdot,pdot,edot
write(67,1000)x2(k)*45.d0/atan(1.d0),w(j,k,l,1)*rho0*w(j,k,l,2)*v_light,&
!     w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*abs(w(j,k,l,2))*v_light,&
!     w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*w(j,k,l,2)*v_light*w(j,k,l,2)*v_light*0.5d0,w(j,k,l,2),w(j,k,l,1)*rho0
     w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*w(j,k,l,2)*v_light,&
     w(j,k,l,1)*rho0*w(j,k,l,2)*v_light*w(j,k,l,2)*v_light*w(j,k,l,2)*v_light*0.5d0,w(j,k,l,2),w(j,k,l,1)*rho0
end if

!if (k==0) then
   w(j,k,l,6)=dm
!else
!   w(j,k,l,6)=w(j,k-1,l,6)+dm
!end if
end do
end if
end do



lm=ldot/mdot
Rl=lm*lm/(rg*rg*v_light*v_light) !rg unit


if(j.ge.2)then
   write(77,1000)x1(j),mdot*v_light*v_light*0.06d0/Ledd, mdot &
        ,pdot*v_light/Ledd, pdot, edot/Ledd, edot,&
        abs(mdot_old-mdot)*100.d0/mdot,abs(pdot_old-pdot)*100.d0/pdot,abs(edot_old-edot)*100.d0/edot
end if

!14
if(j.eq.53)then
   write(6,*)"53",x1(j)
   write(33,1000)mbh,edd_ratio,edd_ratio*Ledd, mdot*v_light*v_light*0.06d0/Ledd, mdot &
        ,pdot*v_light/Ledd, pdot, edot/Ledd, edot
end if

if(j.eq.j0)then
   write(6,*)"j0",x1(j)
   write(34,1000)mbh,edd_ratio,edd_ratio*Ledd, mdot*v_light*v_light/Ledd, mdot &
        ,pdot*v_light/Ledd, pdot, edot/Ledd, edot
end if



mdot_old=mdot
pdot_old=pdot
edot_old=edot


end do

write(6,*) "outflow",mdot,mdot_out,mdot_in,mdot_out+mdot_in

write(6,*)"check2",omega/(x1(j0)*x1(j0)*rg*rg),4*pi
!stop !180810
      
!write(filenum2,'(i4)') f
!filenum2=adjustl(filenum2)
!open (unit=30,file="./Data_nomura/theta_"//trim(filenum2)//".dat")
!open (unit=35,file="./Data_nomura/outflowrate_"//trim(filenum2)//".dat")

!write(34,1000)mbh,edd_ratio,edd_ratio*Ledd, mdot*v_light*v_light*0.06d0/Ledd, mdot &
!     ,pdot*v_light/Ledd, pdot, edot/Ledd, edot

!open (unit=30,file="M"//trim(dnb)//"/"//trim(dn)//"/Data_nomura/theta.dat")
!open (unit=35,file="M"//trim(dnb)//"/"//trim(dn)//"/Data_nomura/outflowrate.dat")

open (unit=30,file="Data/theta.dat")
open (unit=35,file="Data/outflowrate.dat")

write(99,*)"mdot*c^2/Ledd",mdot*v_light*v_light/Ledd, mdot, mdot*v_light*v_light*0.06d0/(Ledd*edd_ratio)
write(99,*)"mdot/Macc",mdot*v_light*v_light*0.06d0/(edd_ratio*Ledd)
write(99,*)"pdot*c/Ledd",pdot*v_light/Ledd, pdot, Ledd/v_light
write(99,*)"edot*/Ledd",edot/Ledd, edot, Ledd
write(99,*)"velocity",pdot/mdot
!pause
!%%%%%%%%%%
write(35,*)"mdot*c^2/Ledd",mdot*v_light*v_light/Ledd, mdot, mdot*v_light*v_light*0.06d0/(Ledd*edd_ratio)
write(35,*)"mdot/Macc",mdot*v_light*v_light*0.06d0/(edd_ratio*Ledd)
write(35,*)"pdot*c/Ledd",pdot*v_light/Ledd, pdot, Ledd/v_light
write(35,*)"edot*/Ledd",edot/Ledd, edot, Ledd
write(35,*)"velocity",pdot/mdot
write(35,*)"r_launch(rg)",Rl
close(35)
!write(6,*)"omega",2.d0*omega/(x1(j0)*x1(j0)*rg*rg),4.d0*pi
1000 format(E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2)

do k=1,k0        
   write(30,1000) x2(k)*45.d0/atan(1.d0),w(j0,k,1,1)*rho0,w(j0,k,1,2)*v_light/1.d5,w(j0,k,1,5),&
        wline(j0,k,1,1),wline(j0,k,1,5)/(0.4d0*mh),w(j0,k,1,1)*w(j0,k,1,2)*rho0*v_light,w(j0,k,1,6),w(j0,k,1,2)

end do
!end do
close(30)

end do

close(34)
close(33)
end do



END program Polar_to_Cartesian
!========================================
