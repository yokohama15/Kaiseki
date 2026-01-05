!--- For Singleo
program Polar_to_Cartesian
implicit none

CHARACTER(4) filenum,filenum2,cn4,cn2

!*****grid point of hydro code****
integer,parameter :: j0=140,k0=100,l0=1
real(8),parameter :: x1max=1.5d6,x1min=3000.0d0,x2max=1.d0 !e60
!real(8),parameter :: x1max=119724.5d0,x1min=3137.713d0,x2max=1.d0 !E6
real,dimension(1:j0) :: x1
real,dimension(1:k0) :: x2
integer :: js,je,ks,ke,ls,le,m,mmax,jp,kp,j,k,l,jp_tmp,kp_tmp,f

real,dimension(:,:,:,:),allocatable :: w
real,dimension(:,:,:,:),allocatable :: wline   

!*.*.*. Cartesian.*.*.*
real :: xxmax,zzmax,dxx,dzz,xxmin,zzmin
integer :: ic,jc
!Cartesian total gred number
integer,parameter:: xxgrid=300,zzgrid=300
!integer,parameter:: xxgrid=5,zzgrid=5
real,dimension(0:xxgrid-1,0:zzgrid-1) :: xx,zz,r,theta
real,dimension(-xxgrid:xxgrid-1) :: xx2
real,dimension(-zzgrid:zzgrid-1) :: zz2
real,dimension(-xxgrid:xxgrid-1,-zzgrid:zzgrid-1)::wc_1
real,dimension(:,:,:),allocatable :: wc
real,dimension(0:xxgrid-1,0:zzgrid-1,2) :: p
real,dimension(0:xxgrid-1,0:zzgrid-1) :: V1,V2,V3,V4
real(8) ::V_tot

real(8), parameter :: mh = 1.67262158D-24
real(8), parameter :: v_light = 29979245800.d0
real(8), parameter :: kb = 1.3807D-16
real(8), parameter :: msun=1.989d33
real(8), parameter :: gc=6.67359d-8
real(8), parameter :: pi=3.141593d0
real(8), parameter :: sigma_e=0.4d0

real(8), parameter :: rho0=1.d-16
real(8) :: mbh, edd_ratio,z_off

real(8) :: dim_tem,rg,Ledd,Lx
real(8) :: v_esc
real(8) ::x_box,z_box,grav_r,grav_theta
real(8),dimension(:,:,:),allocatable ::r_wind,theta_wind
real(8) :: gasp,radp
integer ::fsta,fend,fint
integer ::mark
!open(7,file="free_parameters.dat")
!read(7,*) mbh
!read(7,*) edd_ratio
!close(7)

!write(6,*)"parameters"

!z_off=8.d0*edd_ratio

rg=gc*mbh*msun/(v_light*v_light)
dim_tem=0.5d0*mh*v_light*v_light/kb
Ledd=4.d0*pi*v_light*gc*mbh*msun/sigma_e
Lx=0.05d0*edd_ratio*Ledd

write(6,*)dim_tem
!pause
!**************
js=1
je=j0
ks=1
ke=k0
ls=1
le=l0
mmax=11
!**************


allocate(w(js:je,ks:ke,ls:le,mmax))
allocate(wline(js:je,ks:ke,ls:le,6))
allocate(wc(0:xxgrid-1,0:zzgrid-1,mmax))
allocate(r_wind(js:je,ks:ke,ls:le))
allocate(theta_wind(js:je,ks:ke,ls:le))

!###############grid################
   open (unit=1,file="/Users/yokohama/Documents/work/test251205/Data/grid-x1.data")
   read(1,*) x1(1:j0)   
   close(1)
   open (unit=2,file="/Users/yokohama/Documents/work/test251205/Data/grid-x2.data")
   read(2,*) x2(1:k0)   
   close(2)

   do k=1,k0
      write(66,*)k,x2(k)*45.d0/atan(1.d0)
   end do
!stop
!"*.*.*.*.*.*.*.*.Cartesian gred.*.*.*.*.*.*.*.*.*"
   xxmax=x1max
   xxmin=0.d0
   zzmax=xxmax
   zzmin=0.d0
   dxx=(xxmax-xxmin)/xxgrid
   dzz=(zzmax-zzmin)/zzgrid
!   write(6,*)xxmax,xxmin,zzmax,zzmin
!   write(6,*)xxgrid,zzgrid,dxx,dzz
   
   do ic=0,xxgrid-1
   do jc=0,zzgrid-1

      xx(ic,jc)      =  xxmin+(ic+0.5d0)*dxx
      zz(ic,jc)      =  zzmin+(jc+0.5d0)*dzz
      r(ic,jc)       =  sqrt(xx(ic,jc)*xx(ic,jc)+zz(ic,jc)*zz(ic,jc))
      theta(ic,jc)   =  acos(zz(ic,jc)/r(ic,jc))
            
      
      
      if (r(ic,jc).ge.x1min.and.r(ic,jc).le.x1(j0)) then !simulation box
!      if (r(ic,jc).ge.x1min) then !simulation box

         !****************************
         jp=1
100      continue
         if (r(ic,jc)>x1(jp)) then 
            jp=jp+1
            
!            if (r(ic,jc)>x1(j0))then
               !Cartesian grid poit is out of data point
!               wc(ic,jc,m)=0
!               jp=j0
!               goto 300   !===============>
!            end if
            go to 100
         end if
!300      continue    
!****************************

         kp=1
200      continue
         if (theta(ic,jc)>x2(kp))then
            if (theta(ic,jc)>x2(k0))then
               !Cartesian grid poit is out of data points         
               kp=k0
               go to 400   !================>
            end if
            kp=kp+1
            go to 200 
         end if
    
400      continue   

!if (ic==600.and.jc==600)then
!write(6,*)"check0422-3",jp,kp
!end if 
         p(ic,jc,1)=jp
         p(ic,jc,2)=kp
         !***waight volume***
         if (jp==1.and.1<kp)then
            V1(ic,jc) = -cos(x2(kp))+cos(theta(ic,jc))
            V3(ic,jc) = -cos(theta(ic,jc))+cos(x2(kp-1))
            
         else if (kp==1.and.jp>1) then
            V1(ic,jc) = r(ic,jc)**3-x1(jp-1)**3
            V2(ic,jc) = x1(jp)**3-r(ic,jc)**3
            
         else if (kp==k0.and.jp>1) then
            V1(ic,jc) = r(ic,jc)**3-x1(jp-1)**3
            V2(ic,jc) = x1(jp)**3-r(ic,jc)**3
                        
         else if (jp>1.and.kp>1) then
            V1(ic,jc) = (r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
            V2(ic,jc) = (x1(jp)**3-r(ic,jc)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
            V3(ic,jc) = (r(ic,jc)**3-x1(jp-1)**3)*(-cos(theta(ic,jc))+cos(x2(kp-1)))
            V4(ic,jc) = (x1(jp)**3-r(ic,jc)**3)*(-cos(theta(ic,jc))+cos(x2(kp-1)))
            
         end if
         
      endif
 

      if (V1(ic,jc)<0.or.V2(ic,jc)<0.or.V3(ic,jc)<0.or.V4(ic,jc)<0) then
         write(6,*)"ERROR V<0",ic,jc,V1,V2,V3,V4
         write(6,*)"ERROR2",jp,kp,r(ic,jc),x1(jp-1),x2(kp),theta(ic,jc)
         write(6,*)"ERROR3",r(ic,jc)**3-x1(jp-1)**3,-cos(x2(kp))+cos(theta(ic,j))
         write(6,*)"ERROR3",(r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+cos(theta(ic,jc))),V1
!         pause  
      end if
  
   end do
end do


!each time step
!510 format(7X,I4.4)
!  open (unit=70,file='fset.dat',form='formatted',status='old')
!  read(70,510)fsta
!  read(70,510)fend
!  read(70,510)fint
!  close(70)
!write(6,*)"fset",fsta,fend,fint
!stop
!do f=1000,1000,fint
do f=1181,1200,20
  !output step num
write(6,*)f
!##########Read file##########
   do m=1,11
!      write(6,*)"m=",m
      w(4,4,1,m)=2

      write(cn2,'(i2.2)')m
      
      write(cn4,"(i4.4)")f !1203 
      
      if (m<=4.or.m==8) then
         !write(6,*)"check1",w(4,4,1,m)
         open(10,file='/Users/yokohama/Documents/work/test251205/Data/w'//trim(cn2)//'-'//cn4//'.data' &
              ,form='unformatted')
!         open(10,file='M8_L05/w'//trim(cn2)//'-'//cn4//'.data' &
!              ,form='unformatted')
         read(10) w(js:je,ks:ke,ls:le,m)
         
         close(10)
      
      end if
!      write(6,*)"m=",m
      if (m<=6) then
         open(11,file='/Users/yokohama/Documents/work/test251205/Data/wline'//trim(cn2)//'-'//cn4//'.data' &
              ,form='unformatted')
         read(11) wline(js:je,ks:ke,ls:le,m)
         
         close(11)      
end if

end do

!############################
!do j=js,je
!write(6,*)"check0419",w(j,80,1,2),x1(j),x2(80)*45.d0/atan(1.d0)
!end do

do j=js,je
do k=ks,ke
do l=ls,le

!0521>>>>>>>>>>>>>>
   x_box=x1(j)*sin(x2(k))
   z_box=x1(j)*cos(x2(k))
   
   r_wind(j,k,l)=sqrt(x_box**2.d0+(z_box+z_off)**2.d0)      
   theta_wind(j,k,l)=atan(x_box/(z_box+z_off))
!0521<<<<<<<<<<<<<

grav_r     =cos(x2(k)-theta_wind(j,k,l))/(r_wind(j,k,l)*r_wind(j,k,l))  !マイナス
grav_theta =sin(x2(k)-theta_wind(j,k,l))/(r_wind(j,k,l)*r_wind(j,k,l))
v_esc=sqrt(2.d0/r_wind(j,k,l))

w(j,k,l,5)=w(j,k,l,8)*dim_tem/w(j,k,l,1)  !temperature 
if (j==50.and.k==50) then   !check temperature
   write(6,*) w(j,k,l,5),f
end if
w(j,k,l,6)=r_wind(j,k,l)*w(j,k,l,8)/w(j,k,l,1)    !gas pressure/gravity
w(j,k,l,7)=wline(j,k,l,1)                 !x
w(j,k,l,9)=wline(j,k,l,2)/grav_r
w(j,k,l,10)=wline(j,k,l,3)/grav_theta
w(j,k,l,11)=wline(j,k,l,4)
!w(j,k,l,12)=wline(j,k,l,5)  !tau
!w(j,k,l,13)=(w(j,k,l,3)*w(j,k,l,3)+w(j,k,l,4)*w(j,k,l,4))/(x1(j)*grav_r)
!w(j,k,l,14)=w(j,k,l,2)/v_esc
!w(j,k,l,15)=w(j,k,l,3)/v_esc
!w(j,k,l,16)=sqrt(w(j,k,l,2)*w(j,k,l,2)+w(j,k,l,3)*w(j,k,l,3)+w(j,k,l,4)*w(j,k,l,4))/v_esc
!gasp=w(j,k,l,8)/(w(j,k,l,1)*x1(j))
!radp=wline(j,k,l,2)
!w(j,k,l,17)=gasp/radp
!w(j,k,l,17)=
!if (k>ks.and.k<ke)then
!   w(j,k,l,17)=(w(j,k+1,l,8)-w(j,k-1,l,8))*x1(j)/(w(j,k,l,1)*(x2(k+1)-x2(k-1)))
!else if (k==ks) then
!   w(j,k,l,17)=(w(j,k+1,l,8)-w(j,k,l,8))*x1(j)/(w(j,k,l,1)*(x2(k+1)-x2(k)))
!else if (k==ke) then
!   w(j,k,l,17)=(w(j,k,l,8)-w(j,k-1,l,8))*x1(j)/(w(j,k,l,1)*(x2(k)-x2(k-1)))
!end if
!w(j,k,l,18)=x1(j)*(w(j,k,l,2)*w(j,k,l,3)-w(j,k,l,4)*w(j,k,l,4)/tan(x2(k)))
!w(j,k,l,19)=wline(j,k,l,1)*w(j,k,l,1)*x1(j)*x1(j)*rho0*rg*rg/(mh*Lx)
!w(j,k,l,20)=-log(w(j,k,l,19))
!if (j==50.and.k==80)then
!write(6,*)w(j,k,l,19),w(j,k,l,20),f,wline(j,k,l,1),w(j,k,l,1),x1(j),rho0,rg
!end if
!w(j,k,l,3)=-w(j,k,l,3)
!write(6,*)w(25,25,l,9)
   if (j==50.and.f==849) then   !Frad_r VS theta plot
      write(100,*) x2(k),wline(j,k,l,2)
   end if
   

end do
end do
end do
!write(6,*)"check2"

      
do ic=0,xxgrid-1
do jc=0,zzgrid-1


   jp_tmp=p(ic,jc,1)
   kp_tmp=p(ic,jc,2)
!data convert
   do m=1,mmax
      if (jp_tmp==1.and.1<kp_tmp)then
         V_tot=V1(ic,jc)+V3(ic,jc)
         wc(ic,jc,m)=(w(jp_tmp,kp_tmp-1,1,m)*V3(ic,jc)+w(jp_tmp,kp_tmp,1,m)*V1(ic,jc))/V_tot
      else if (kp_tmp==1.and.jp_tmp>1) then
         V_tot=V1(ic,jc)+V2(ic,jc)
         wc(ic,jc,m)=(w(jp_tmp-1,kp_tmp,1,m)*V2(ic,jc)+w(jp_tmp,kp_tmp,1,m)*V1(ic,jc))/V_tot
      else if (kp_tmp==k0.and.jp_tmp>1) then
         V_tot=V1(ic,jc)+V2(ic,jc)
         wc(ic,jc,m)=(w(jp_tmp-1,kp_tmp,1,m)*V2(ic,jc)+w(jp_tmp,kp_tmp,1,m)*V1(ic,jc))/V_tot         
      else if (jp_tmp==1.and.kp_tmp==1) then
         wc(ic,jc,m)=w(jp_tmp,kp_tmp,1,m)
      else if (jp_tmp>1.and.kp_tmp>1)  then
         V_tot=V1(ic,jc)+V2(ic,jc)+V3(ic,jc)+V4(ic,jc)
         wc(ic,jc,m)=(w(jp_tmp-1,kp_tmp,1,m)*V4(ic,jc)+w(jp_tmp,kp_tmp,1,m)*V3(ic,jc)+w(jp_tmp-1,kp_tmp-1,1,m)*V2(ic,jc)+w(jp_tmp,kp_tmp-1,1,m)*V1(ic,jc))/V_tot
      end if
     


!if (m==2.and.f=60)then
!write(6,*)"check0422",
   end do
end do
end do
!write(6,*)"check3"

write(filenum2,'(i4.4)') f
filenum2=adjustl(filenum2)
!open (unit=30,file="./M8_L05/Data_yokohama/mvdatas_"//trim(filenum2)//".dat")
open (unit=35,file="./movie/mvdata_"//trim(filenum2)//".dat")
!open (unit=30,file="./M8_L05/Data_yokohama/widedatas_"//trim(filenum2)//".dat")
!open (unit=31,file="./Data_yokohama/visitrho_"//trim(filenum2)//".dat",status="unknown")
1000 format(E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2,E20.8e2)

do ic=0,xxgrid-1
   xx2(ic)=xx(ic,1)
end do
do ic=-xxgrid,-1
   xx2(ic)=-xx(-ic-1,1)
end do

do jc=0,zzgrid-1
   zz2(jc)=zz(1,jc)
end do
do jc=-zzgrid,-1
   zz2(jc)=-zz(1,-jc-1)
end do

do ic=0,xxgrid-1
   do jc=-zzgrid,-1
      wc_1(ic,jc)=wc(ic,-jc-1,1)
   end do
   do jc=0,zzgrid-1
      wc_1(ic,jc)=wc(ic,jc,1)
   end do
end do
do ic=-xxgrid,-1
   do jc=-zzgrid,-1
      wc_1(ic,jc)=wc(-ic-1,-jc-1,1)
   end do
   do jc=0,zzgrid-1
      wc_1(ic,jc)=wc(-ic-1,jc,1)
   end do
end do


!do ic=-xxgrid,xxgrid-1        
!   write(6,*)"test i",ic, xx2(ic),xx2(ic)-xx2(ic-1)
!end do
!pause
!do jc=-zzgrid,zzgrid-1
!   write(6,*)"test j",jc, zz2(jc),zz2(jc)-zz2(jc-1)
!end do
!pause
!########### Write file #############
!do ic=-xxgrid,xxgrid-1        
!   write(30,*)""
!   do jc=-zzgrid,zzgrid-1
!      
!      if (sqrt(xx2(ic)*xx2(ic)+zz2(jc)*zz2(jc))>=x1min.and.sqrt(xx2(ic)*xx2(ic)+zz2(jc)*zz2(jc)).le.x1(j0)) then
!         write(30,1000)xx2(ic)-0.5d0*dxx,zz2(jc)-0.5d0*dzz,wc_1(ic,jc)
!         write(30,1000)xx2(ic)-0.5d0*dxx,zz2(jc)+0.5d0*dzz,wc_1(ic,jc)
!      else
!         write(30,1000)xx2(ic)-0.5d0*dxx,zz2(jc)-0.5d0*dzz,0.d0
!         write(30,1000)xx2(ic)-0.5d0*dxx,zz2(jc)+0.5d0*dzz,0.d0
!      end if
!   end do
!   write(30,*)""
!   do jc=-zzgrid,zzgrid-1
      
      !write(1,*)i-0.5,k-0.5,u(i,j,k),v(i,j,k),w(i,j,k)
      !write(1,*)i-0.5,k+0.5,u(i,j,k),v(i,j,k),w(i,j,k)
!      if (sqrt(xx2(ic)*xx2(ic)+zz2(jc)*zz2(jc))>=x1min.and.sqrt(xx2(ic)*xx2(ic)+zz2(jc)*zz2(jc)).le.x1(j0)) then
!         write(30,1000)xx2(ic)+0.5d0*dxx,zz2(jc)-0.5d0*dzz,wc_1(ic,jc)
!         write(30,1000)xx2(ic)+0.5d0*dxx,zz2(jc)+0.5d0*dzz,wc_1(ic,jc)
!      else
!         write(30,1000)xx2(ic)+0.5d0*dxx,zz2(jc)-0.5d0*dzz,0.d0
!         write(30,1000)xx2(ic)+0.5d0*dxx,zz2(jc)+0.5d0*dzz,0.d0
!      end if

            
!   end do
   
!end do
!end do
!close(30)

do ic=-xxgrid,xxgrid-1        
!   write(35,*)""
   mark=0
   do jc=-zzgrid,zzgrid-1      
      
      if (sqrt(xx2(ic)*xx2(ic)+zz2(jc)*zz2(jc))>=x1min.and.sqrt(xx2(ic)*xx2(ic)+zz2(jc)*zz2(jc)).le.x1(j0)) then
         write(35,1000)xx2(ic)*0.5d0,zz2(jc)*0.5d0,wc_1(ic,jc)*1.d-16,wc(ic,jc,5),wc(ic,jc,2),wc(ic,jc,9),wc(ic,jc,11)
         !write(6,1000)xx2(ic)*0.5d0,zz2(jc)*0.5d0,wc_1(ic,jc)*1.d-16
         
         if (wc_1(ic,jc)*1.d-16.lt.1.d-23)then
            write(6,*)"data output error1",xx2(ic)*0.5d0,zz2(jc)*0.5d0,wc_1(ic,jc)*1.d-16,wc(ic,jc,8)
            stop
         end if
             
         mark=1
      else
         
      end if
   end do
   
   if (mark.eq.1)then
      write(35,*)""
      !write(6,*)"a"
      
   else if (mark.eq.0)then
   else
      write(6,*)"data output error2",mark
   end if
   
end do
!pause
close(35)
   
!end do
!end do

!do ic=1,100       
!do jc=1,100
!do lc=1,100
!test(ic,jc,lc)=ic+jc+lc

!      if (r(ic,jc)>=x1min.and.r(ic,jc).le.x1(j0)) then
!         write(31,*)float(ic),float(jc),float(ic+jc)
!      end if
!end do
!end do
!end do
!CLOSE(31)
   
   
   !itt=itt+1
end do

!do ic=1,xxgrid
!write(6,*)"test2",wc(ic,50,9)
!end do
!
!end do

END program Polar_to_Cartesian
!========================================
