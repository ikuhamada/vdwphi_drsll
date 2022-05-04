!     ==================================================================
      MODULE M_VDW_PHI
!     ==--------------------------------------------------------------==
!     Module for calculating the vdW kernel phi
!     Written by Ikutaro Hamada
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
! ... # of mesh points and divisions 
!     INTEGER      :: ma=3000 ! for a testing
      INTEGER      :: ma=6000 ! for production
      INTEGER      :: na
! ... max. value for a and b, and interval
      REAL(KIND=8) :: amax=200.D0
      REAL(KIND=8) :: da
! ... array for variable a
      REAL(KIND=8), ALLOCATABLE :: ai(:)
! ... array for function W(a,b)
      REAL(KIND=8), ALLOCATABLE :: Wab(:,:)
! ... variables and arrays for spline interpolation
      REAL(KIND=8), ALLOCATABLE :: finta(:),fintb(:)
      REAL(KIND=8), ALLOCATABLE :: finta2(:),fintb2(:)
! ... constants
      REAL(KIND=8) :: SMALL=1.D-15
      REAL(KIND=8) :: PI=4.D0*ATAN(1.D0)
!     ==--------------------------------------------------------------==
!
      CONTAINS
!
!     ==================================================================
      SUBROUTINE SET_AMESH
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER :: ia
!     ==--------------------------------------------------------------==
! ... mesh interval
      da=amax/DBLE(ma)
      na=ma+1
! ... allocate
      ALLOCATE(ai(na))
      ai(:)=0.D0
! ... set a on a grid
      DO ia=1,na
        ai(ia)=DBLE(ia-1)*da
      ENDDO
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE SET_AMESH
!     ==================================================================
      FUNCTION W(a,b)
!     ==--------------------------------------------------------------==
!     Calculate W function [Eq.(16) in PRL 92, 246401 (2004)]
!     NB: W(a,b) is mutiplied by a**2 b**2 because of the factor in the
!     integral [see Eq.(14)]
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL(KIND=8) :: W
! ... Arguments
      REAL(KIND=8) :: a,b
! ... Variables
      REAL(KIND=8) :: a2,b2,sina,sinb,cosa,cosb
!     ==--------------------------------------------------------------==
      W=0.d0
      a2=a**2
      b2=b**2
      sina=sin(a)
      sinb=sin(b)
      cosa=cos(a)
      cosb=cos(b)
      IF(a.LT.SMALL)THEN
        IF(b.LT.SMALL)THEN
          W=2.D0*((3.D0-a2)+(3.D0-b2)+(a2+b2-3.D0)-3.D0)
        ELSE
          W=2.D0*((3.D0-a2)*b*cosb+(3.D0-b2)*sinb+&
     &            (a2+b2-3.D0)*sinb-3.D0*b*cosb)/b
        ENDIF
      ELSE
        IF(b.LT.SMALL)THEN
          W=2.D0*((3.D0-a2)*sina+(3.D0-b2)*a*cosa+&
     &            (a2+b2-3.D0)*sina-3.D0*a*cosa)/a
        ELSE
          W=2.D0*((3.D0-a2)*b*cosb*sina+(3.D0-b2)*a*cosa*sinb+&
     &            (a2+b2-3.D0)*sina*sinb-3.D0*a*b*cosa*cosb)/(a*b)
        ENDIF
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END FUNCTION W
!     ==================================================================
      SUBROUTINE SET_W_AB
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER      :: ia,ib
      REAL(KIND=8) :: a,b
!     ==--------------------------------------------------------------==
      ALLOCATE(Wab(na,na))
! ... set W(a,b) on a grid
      DO ib=1,na
        DO ia=1,na
          a=ai(ia)
          b=ai(ib)
          Wab(ia,ib)=W(a,b)
        ENDDO
      ENDDO
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE SET_W_AB
!     ==================================================================
      FUNCTION T(w,x,y,z) 
!     ==--------------------------------------------------------------==
!     Calculate the T function [Eq.(15) in PRL 92, 246401 (2004)]
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL(KIND=8) :: T
! ... Arguments
      REAL(KIND=8) :: w,x,y,z
!     ==--------------------------------------------------------------==
      T=0.5D0*(1.D0/(w+x)+1.D0/(y+z))* &
     &        (1.D0/(w+y)/(x+z)+1.D0/(w+z)/(y+x)) 
!     ==--------------------------------------------------------------==
      RETURN
      END FUNCTION T
!     ==================================================================
      FUNCTION NU(y,d)
!     ==--------------------------------------------------------------==
!     Calculate the function nu(y)=y**2/2h(y/d),
!     where h(x)=1-exp(4*pi*x**2/9)
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL(KIND=8) :: NU
! ... Arguments
      REAL(KIND=8) :: y,d
! ... Variables
      REAL(KIND=8) :: gamma
!     ==--------------------------------------------------------------==
      gamma=4.D0*pi/9.D0
      IF(y.LE.SMALL .and. d.GT.SMALL)then
        NU=d**2/(2.D0*gamma) 
      ELSEIF(d.LE.SMALL)THEN
        NU=0.5D0*y**2
      ELSE
        NU=y**2/(2.D0*(1.D0-exp(-gamma*((y/d)**2))))
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END FUNCTION NU
!     ==================================================================
      FUNCTION PHIVAL(d1,d2)
!     ==--------------------------------------------------------------==
!     Calculate the nonlocal correlation kernel for the van der density
!     functional of Dion et al. [Phys. Rev. Lett. 92, 246401 (2004)]
!     at give d1 and d2
!     Integrand [W(a,b)*T(nu(a),nu(b),nu'(a),nu'(b))] is interpolated
!     by the cubic spline interpolation, and the integrations are done
!     analytically.
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL(KIND=8) :: PHIVAL
! ... Arguments
      REAL(KIND=8) :: d1,d2
! ... Variables
!     a, b, nu(a), nu(b), nu'(a), and nu'(b)
      INTEGER      :: ia,ib
      REAL(KIND=8) :: a,b,nu1a,nu1b,nu2a,nu2b
      REAL(KIND=8) :: yp1=1.D+30,ypn=1.D+30,tval
!     ==--------------------------------------------------------------==
      IF(.NOT. ALLOCATED(finta))  ALLOCATE(finta(na))
      IF(.NOT. ALLOCATED(fintb))  ALLOCATE(fintb(na))
      IF(.NOT. ALLOCATED(finta2)) ALLOCATE(finta2(na))
      IF(.NOT. ALLOCATED(fintb2)) ALLOCATE(fintb2(na))
      finta(:)=0.D0
      LOOP_A: DO ia=1,na
        fintb(:)=0.D0
        LOOP_B: DO ib=1,na
          a=ai(ia)
          b=ai(ib)
          nu1a=nu(a,d1)
          nu1b=nu(b,d1)
          nu2a=nu(a,d2)
          nu2b=nu(b,d2)
!DEBUG
          if((nu2a.le.small .and. nu2b.le.small) .or.&
     &       (nu1a.le.small .and. nu1b.le.small)     )then
            cycle
          endif
!ENDDEBUG
          tval=T(nu1a,nu1b,nu2a,nu2b)
          fintb(ib)=Wab(ia,ib)*tval
        ENDDO LOOP_B
!       spline interpolation of the integrand
        CALL spline(ai,fintb,na,yp1,ypn,fintb2)
!       integration over the variable b
        CALL intspl(ai,fintb,fintb2,na,finta(ia))
      ENDDO LOOP_A
      PHIVAL=0.D0
      CALL spline(ai,finta,na,yp1,ypn,finta2)
      CALL intspl(ai,finta,finta2,na,PHIVAL)
!     multiply by 2/pi**2
      PHIVAL=2.D0/(pi**2)*PHIVAL
!     ==--------------------------------------------------------------==
      RETURN
      END FUNCTION PHIVAL
!     ==================================================================
      SUBROUTINE ALLOC_ARR_PHI
!     ==--------------------------------------------------------------==
!     deallocate temporal arrays for the calculation of phi
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
!     ==--------------------------------------------------------------==
      ALLOCATE(finta(na),fintb(na))
      ALLOCATE(finta2(na),fintb2(na))
!     ALLOCATE(ai)
!     ALLOCATE(Wab)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE ALLOC_ARR_PHI
!     ==================================================================
      SUBROUTINE DEALLOC_ARR_PHI
!     ==--------------------------------------------------------------==
!     deallocate temporal arrays for the calculation of phi
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
!     ==--------------------------------------------------------------==
      IF(ALLOCATED(finta))  DEALLOCATE(finta)
      IF(ALLOCATED(fintb))  DEALLOCATE(fintb)
      IF(ALLOCATED(finta2)) DEALLOCATE(finta2)
      IF(ALLOCATED(fintb2)) DEALLOCATE(fintb2)
      DEALLOCATE(ai)
      DEALLOCATE(Wab)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE DEALLOC_ARR_PHI
!     ==================================================================
      SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)
!     ==--------------------------------------------------------------==
!     cubic spline interpolation taken from NUMERICAL RECIPES
!     set the natural boundary contition
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
! ... Arguments
      INTEGER      :: n
      REAL(KIND=8) :: yp1,ypn
      REAL(KIND=8) :: x(n),y(n),y2(n)
! ... Variables
      INTEGER                   :: i,k
      REAL(KIND=8)              :: p,qn,sig,un
      REAL(KIND=8), ALLOCATABLE :: u(:)
!     ==--------------------------------------------------------------==
      ALLOCATE(u(n))
      IF (yp1.gt..99d30) THEN
        y2(1)=0.d0
        u(1)=0.d0
      ELSE
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF
      DO i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=( &
     &        6.d0*( (y(i+1)-y(i)  )/(x(i+1)-x(i)) &
     &              -(y(i  )-y(i-1))/(x(i)-x(i-1))  ) &
     &       /(x(i+1)-x(i-1))-sig*u(i-1) ) / p
      ENDDO
      IF (ypn.gt..99d30) THEN
        qn=0.d0
        un=0.d0
      ELSE
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      DO k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      ENDDO
      DEALLOCATE(u)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE spline
!     ==================================================================
      SUBROUTINE INTSPL(x,y,y2,n,yint)
!     ==--------------------------------------------------------------==
!     perform analytic integration of the function y
!     obtained by the cubic spline interpolation
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
! ... Arguments
      INTEGER      :: n
      REAL(KIND=8) :: x(n),y(n),y2(n)
      REAL(KIND=8) :: yint
! ... Variables
      INTEGER :: i
!     ==--------------------------------------------------------------==
      yint=0.d0
      DO i=1,n-1
        yint=yint+0.5d0*(x(i+1)-x(i))*(y(i+1)+y(i))- &
     &       (1.d0/24.d0)*(x(i+1)-x(i))**3 * (y(i+1)+y(i))**3
      ENDDO
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE intspl
!     ==================================================================
!
!     ==--------------------------------------------------------------==
      END MODULE M_VDW_PHI
!     ==================================================================
