!--- For Single
program Polar_to_Cartesian
implicit none
!real :: t,dt
!integer :: nxh, nyh, nzh, xthread, ythread, zthread, is1, ie1, js1, je1, ks1, ke1, ist, ien, jst, jen, kst, ken, iup, idown, jup, jdown, kup, kdown
integer :: itt,f,i,j,k
!integer,parameter :: ist=0,ien=127,jst=0,jen=127,kst=0,ken=63
!real,dimension(0:ien,0:jen,0:ken) :: u,v,w,rho,tem,p,xi_x,v_r2,v_phi2,Fradx,Frady,Fradz
CHARACTER(4) filenum,filenum2,cn4,cn2

!*****grid point of hydro code****
integer,parameter :: j0=140,k0=100,l0=1
real(8),parameter :: x1max=1.5d6,x1min=3000.0d0,x2max=1.d0 !e60
!real(8),parameter :: x1max=3000.d0,x1min=60,x2max=1.d0
real,dimension(1:j0) :: x1
real,dimension(1:k0) :: x2

integer :: js,je,ks,ke,ls,le,m,mmax,jp,kp
integer :: l
integer :: nf
real,dimension(:,:,:,:),allocatable :: w  
!real,dimension(:,:,:),allocate :: temp
!*.*.*. Cartesian.*.*.*
real :: xxmax,zzmax,dxx,dzz,xxmin,zzmin
integer :: ic,jc
!Cartesian total gred number
integer,parameter:: xxgrid=900,zzgrid=900
real,dimension(0:xxgrid+1,0:zzgrid+1) :: xx,zz,r,theta
real,dimension(:,:,:),allocatable :: wc
real :: V1,V2,V3,V4

real(8), parameter :: mh = 1.67262158D-24
real(8), parameter :: v_light = 29979245800.d0
real(8), parameter :: kb = 1.3807D-16
real(8) :: dim_tem

!*.*.*.*.*.*.*.*.*.*.*.*.*

  dim_tem=0.5d0*mh*v_light*v_light/kb

js=1
je=j0
ks=1
ke=k0
ls=1
le=l0
mmax=9

nf=0
allocate(w(js:je,ks:ke,ls:le,mmax))
allocate(wc(0:xxgrid,0:zzgrid,mmax))
!allocate(temp(js:je,ks:ke,ls:le))
!itt=1

do f=1181,1200,20  !output step num

!do m=1,8,7
do m=2,4
!write(6,*)"m=",m

write(cn2,'(i2.2)')m


!   write(filenum,'(i4.4)') itt  !like 0001, 0010, 0100
   
!   write(filenum2,'(i4)') f
!   filenum2=adjustl(filenum2)
!   open (unit=30,file="./Data_nomura/text_w"//trim(cn2)//"_"//trim(filenum2)//".dat")
!   write(6,*)"*.*.*.*.*.*.*.*.*.*.*.OPEN*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*."   
!   *.*.*.*.*deta reading*.*.*.*.*
!   ##### grid points #####

   open (unit=1,file="/Users/yokohama/Documents/work/test251205/Data/grid-x1.data")
   read(1,*) x1(1:j0)   
   close(1)
   open (unit=2,file="/Users/yokohama/Documents/work/test251205/Data/grid-x2.data")
   read(2,*) x2(1:k0)   
   close(2)
!   ##### result data ##### 
  write(cn4,"(i4.4)")f !1203 
!  write(cn5,"(i5.5)")myrk
!-----------------------------------------------------------
! read outputfile primitivie variables
!-----------------------------------------------------------
!  do m=1,1 !m=1:density
!     ii=0
  !write(cn2,'(i2.2)')m
  open(10,file='/Users/yokohama/Documents/work/test251205/Data/w'//trim(cn2)//'-'//cn4//'.data' &
       ,form='unformatted')
  read(10) w(js:je,ks:ke,ls:le,m)
  
  close(10)

if (m==8) then
do j=js,je
do k=ks,ke
do l=ls,le
w(j,k,l,9)=w(j,k,l,8)*dim_tem/w(j,k,l,1) 
if (w(j,k,l,9)>1.d12) then
write(6,*)"hightemp",w(j,k,l,9),j,k,l
endif
end do
end do
end do
end if

end do



!end do

do m=2,4
!do m=1,9,8

write(cn2,'(i2.2)')m


!   write(filenum,'(i4.4)') itt  !like 0001, 0010, 0100
   
   write(filenum2,'(i4.4)') f
   filenum2=adjustl(filenum2)
!   open (unit=30,file="./Data_nomura/text_w"//trim(cn2)//"_"//trim(filenum2)//".dat")
   if (m==3) then
      open (unit=35,file="./movie/vecdata_"//trim(filenum2)//".dat")
end if
   
!write(6,*)"*.*.*.*.*.*.*.*.Cartesian gred.*.*.*.*.*.*.*.*.*"
 xxmax=x1max
 xxmin=0.d0
 zzmax=xxmax
 zzmin=0.d0
 dxx=(xxmax-xxmin)/xxgrid
 dzz=(zzmax-zzmin)/zzgrid
!write(6,*)xxmax,xxmin,zzmax,zzmin
!write(6,*)xxgrid,zzgrid,dxx,dzz

 do ic=0,xxgrid
 do jc=0,zzgrid
!write(6,*)"do in"
    xx(ic,jc)=xxmin+(ic+0.5d0)*dxx
    zz(ic,jc)=zzmin+(jc+0.5d0)*dzz
!write(6,*)"***",ic,jc,xx(ic,jc),zz(ic,jc)
    r(ic,jc)=sqrt(xx(ic,jc)*xx(ic,jc)+zz(ic,jc)*zz(ic,jc))
    if (r(ic,jc)==0) then
       write(6,*)"r=0!!",ic,jc
    end if
    theta(ic,jc)= acos(zz(ic,jc)/r(ic,jc))
    if (ic==0.and.jc==0) then
       !write(6,*)r(ic,jc),theta(ic,jc)
    end if
    

if (r(ic,jc).ge.x1min.and.r(ic,jc).le.x1max) then
    jp=1
100 continue
    !write(6,*)"CHECK START",ic,jc,jp
    !write(6,*)"CHECK START2",r(ic,jc),x1(jp),x1(j0)
    !if (jp>295) then
    !pause  
    !end if
    !if (jc==4) then
    !   write(6,*)r(ic,jc),x1(jp)
    !end if
    if (r(ic,jc)>x1(jp)) then 
       jp=jp+1
       !if (ic==0.and.jc==0)then
       !   write(6,*)jc,j0
       !end if
       if (r(ic,jc)>x1(j0))then
          !Cartesian grid poit is out of area
          wc(ic,jc,m)=0
          
          goto 300   !===============>
       end if
       go to 100
    end if
    
    !write(6,*)"r check done"
    kp=1
    
200 continue
    if (theta(ic,jc)>x2(kp))then
       ! kp=kp+1 
       if (theta(ic,jc)>x2(k0))then
          !Cartesian grid poit is out of area         
          !x2(k0)<theta(ic,jc)<180^circ
          !write(6,*)"1211",kp,k0,x2(kp),x2(k0)
          kp=k0
          go to 400   !================>
       end if
       kp=kp+1
       go to 200 
    end if
    
400 continue   
    !deta on Cartesian grid
    !***waight volume***
    if (jp==1)then
       V1=(r(ic,jc)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
       !V2=(x1(jp)**3-r(ic,jc)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
       V3=(r(ic,jc)**3)*(-cos(theta(ic,jc))+cos(x2(kp-1)))
       !V4=(x1(jp)**3-r(ic,jc)**3)*(-cos(theta(ic,jc))+cos(x2(kp-1)))
       wc(ic,jc,m)=(w(jp,kp,1,m)*V3+w(jp,kp-1,1,m)*V1)/(V1+V3)
       !     write(6,*)"checkjp=1",xx(ic,jc),zz(ic,jc),ic,jc,jp,kp,wc(ic,jc,m)
       !     write(6,*)"CHECKjp=1",w(jp,kp,1,m),w(jp,kp-1,1,m),V1,V3
       !pause
    else if (kp==1) then
       !     write(6,*)"ERROR kp=1",ic,jc,x2(kp-1),theta(ic,jc),x2(kp),x2(1)
       !     pause
       V1=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
       V2=(x1(jp)**3-r(ic,jc)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
       !     V1=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+1)
       !     V2=(x1(jp)**3-r(ic,jc)**3)*(-cos(x2(kp))+1)
       !     V3=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(theta(ic,jc))+1)
       !     V4=(x1(jp)**3-r(ic,jc)**3)*(-cos(theta(ic,jc))+1)
       wc(ic,jc,m)=(w(jp-1,kp,1,m)*V2+w(jp,kp,1,m)*V1)/(V1+V2)
    else if (kp==k0) then
       !     write(6,*)"ERROR kp=1",ic,jc,x2(kp-1),theta(ic,jc),x2(kp),x2(1)
       !     pause
       V1=(r(ic,jc)**3-x1(jp-1)**3)*(1+cos(x2(kp)))
       V2=(x1(jp)**3-r(ic,jc)**3)*(1+cos(x2(kp)))
       !     V1=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+1)
       !     V2=(x1(jp)**3-r(ic,jc)**3)*(-cos(x2(kp))+1)
       !     V3=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(theta(ic,jc))+1)
       !     V4=(x1(jp)**3-r(ic,jc)**3)*(-cos(theta(ic,jc))+1)
       !    write(6,*)"1211",V1,V2
       wc(ic,jc,m)=(w(jp-1,kp,1,m)*V2+w(jp,kp,1,m)*V1)/(V1+V2)
       
    else
       V1=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
       V2=(x1(jp)**3-r(ic,jc)**3)*(-cos(x2(kp))+cos(theta(ic,jc)))
       V3=(r(ic,jc)**3-x1(jp-1)**3)*(-cos(theta(ic,jc))+cos(x2(kp-1)))
       V4=(x1(jp)**3-r(ic,jc)**3)*(-cos(theta(ic,jc))+cos(x2(kp-1)))
       wc(ic,jc,m)=(w(jp-1,kp,1,m)*V4+w(jp,kp,1,m)*V3+w(jp-1,kp-1,1,m)*V2+w(jp,kp-1,1,m)*V1)/(V1+V2+V3+V4)
if (jc==50.and.m==9)then
if (wc(ic-1,jc,m)<wc(ic,jc,m))then
!write(6,*)ic,wc(ic-1,jc,m),wc(ic,jc,m)
!write(6,*)w(jp-1,kp,1,m),V1,"*",w(jp,kp,1,m),V2,"*",w(jp-1,kp-1,1,m),V3,"*",w(jp,kp-1,1,m),V4

end if
endif
    
end if
    
    ! if (xx(ic,jc)>45.d0.and.zz(ic,jc)<4.d0.and.log(wc(ic,jc,m))>10.d0) then
    !    write(6,*)"check1029",ic,jc,jp,kp,wc(ic,jc,m)
    !    write(6,*)"CHECK1029-2",w(jp-1,kp,1,m),w(jp,kp,1,m),w(jp-1,kp-1,1,m),w(jp,kp-1,1,m),V1,V2,V3,V4
    
    !pause  
    !  end if

  if (V1<0.or.V2<0.or.V3<0.or.V4<0) then
     write(6,*)"ERROR V<0",ic,jc,V1,V2,V3,V4
     write(6,*)"ERROR2",jp,kp,r(ic,jc),x1(jp-1),x2(kp),theta(ic,jc)
     write(6,*)"ERROR3",r(ic,jc)**3-x1(jp-1)**3,-cos(x2(kp))+cos(theta(ic,jc))
     write(6,*)"ERROR3",(r(ic,jc)**3-x1(jp-1)**3)*(-cos(x2(kp))+cos(theta(ic,jc))),V1
     !pause  
  end if
  if (m==1.and.wc(ic,jc,m)<0) then
     write(6,*)"ERROR wc(m=1)=0",ic,jc,wc(ic,jc,m)
     !pause
  end if

!write(6,*)"wc"    
300 continue  !<=========================
  if (m==1.and.wc(ic,jc,m)<0) then
     write(6,*)"ERROR 300 wc(m=1)=0",ic,jc,wc(ic,jc,m)
     !pause
  end if
end if
  
end do
end do



!write(6,*)"*.*.*.*.*.*.*.*.*.*.*.OUTPUT FILE FOR GNUPLOT*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*."
!write(6,*)cos(x2(k0))
!open (unit=1,file="text"//trim(filenum2)//"-1.dat",status="unknown")
   !open (unit=30,file="text"//trim(filenum2)//".data",status="new")
!   open (unit=30,file="text.dat",status="new")
!write(6,*)"*.*.*.*.*.*.*.*.*.*.*.OPEN*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*."

!write(6,*)"120412",rho(31,32,1),rho(32,32,1),rho(55,32,1)
!   write(1,*)"*****************",t,dt,"*****************"
!   write(2,*)"*****************",t,dt,"*****************"
!do ic=0,xxgrid+1        
!   write(30,*)""
!   do jc=0,zzgrid+1
      
      !write(1,*)i-0.5,k-0.5,u(i,j,k),v(i,j,k),w(i,j,k)
      !write(1,*)i-0.5,k+0.5,u(i,j,k),v(i,j,k),w(i,j,k)
!      if (r(ic,jc)>=x1min.and.r(ic,jc)<=x1max) then
!         write(30,*)xx(ic,jc)-0.5d0*dxx,zz(ic,jc)-0.5d0*dzz,wc(ic,jc,m),log10(wc(ic,jc,m))
!         write(30,*)xx(ic,jc)-0.5d0*dxx,zz(ic,jc)+0.5d0*dzz,wc(ic,jc,m),log10(wc(ic,jc,m))
!      end if
!   end do
 !  write(30,*)""
!   do jc=0,zzgrid+1
      
      !write(1,*)i-0.5,k-0.5,u(i,j,k),v(i,j,k),w(i,j,k)
      !write(1,*)i-0.5,k+0.5,u(i,j,k),v(i,j,k),w(i,j,k)
!      if (r(ic,jc)>=x1min.and.r(ic,jc)<=x1max) then
!         write(30,*)xx(ic,jc)+0.50*dxx,zz(ic,jc)-0.5d0*dzz,wc(ic,jc,m),log10(wc(ic,jc,m))
!         write(30,*)xx(ic,jc)+0.50*dxx,zz(ic,jc)+0.5d0*dzz,wc(ic,jc,m),log10(wc(ic,jc,m))
!     end if
!      end do
   
   !   write(1,*)i+0.5,k-0.5,u(i,j,k),v(i,j,k),w(i,j,k)
   !   write(1,*)i+0.5,k+0.5,u(i,j,k),v(i,j,k),w(i,j,k)
   
   !               write(2,*)i+0.5,k-0.5,log10(rho(i,j,k)),log10(tem(i,j,k)),log10(p(i,j,k))
   !               write(2,*)i+0.5,k+0.5,log10(rho(i,j,k)),log10(tem(i,j,k)),log10(p(i,j,k))
   !            end do
   
   
   
   
!   end do
!   end do


if (m==3) then
do ic=0,xxgrid+1,45        
   do jc=0,zzgrid+1,45
      if (r(ic,jc)>=x1min.and.r(ic,jc)<=x1max) then
         write(35,*)xx(ic,jc)*0.5d0,zz(ic,jc)*0.5d0,0.d0,(wc(ic,jc,2)*sin(theta(ic,jc))+wc(ic,jc,3)*cos(theta(ic,jc))),&
              (wc(ic,jc,2)*cos(theta(ic,jc))-wc(ic,jc,3)*sin(theta(ic,jc))),0.d0
      end if
   end do
end do

close(35)
end if

   
!   close(30)
   !CLOSE(2)
   
   
   !itt=itt+1
end do

do ic=1,xxgrid
!write(6,*)"test2",wc(ic,50,9)
end do

end do

END program Polar_to_Cartesian
!========================================
