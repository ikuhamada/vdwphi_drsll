 program main
 implicit none
 integer :: id,idel,nd,md,ndel,mdel
 real(kind=8) :: dmax,delmax=1.d0,dd,ddel,d,del
 real(kind=8), allocatable :: phi(:,:),phi0(:),phi0intfunc(:)
 real(kind=8) :: phi0int=0.d0,dint
 real(kind=8) :: pi
!***********************************************************************
 pi=4.d0*atan(1.d0)
 read(*,*)dmax
 read(*,*)md
 read(*,*)mdel
 nd  =md+1
 ndel=mdel+1
 dd=dmax/dble(md)
 ddel=1.d0/dble(mdel)
 allocate(phi(nd,ndel),phi0(nd),phi0intfunc(nd))
 phi(:,:)=0.d0;phi0(:)=0.d0
 do id=1,nd
   do idel=1,ndel
     read(*,*)d,del,phi(id,idel)
   enddo
 enddo
 phi(:,:)=dble(phi(:,:))
 do id=1,nd
   phi0(id)=phi(id,1)
 enddo
 dint=0.d0
 phi0int=0.d0
 do id=11,nd
   d=dble(id-1)*dd
   dint=dint+dd
   phi0intfunc(id)=phi0int
   phi0int=phi0int+dd*4.d0*d*d*phi0(id)
 enddo
 write(*,*)'Int dD                         :',dint
 write(*,*)'Int dD 4 * PI * D**2 * phi(0,D):',phi0int
 do id=11,nd
   d=dble(id-1)*dd
   write(10,'(4f20.12)')d,phi0(id),4.d0*pi*d*d*phi0(id),phi0intfunc(id)
 enddo
 deallocate(phi,phi0,phi0intfunc)
!***********************************************************************
 stop
 end
