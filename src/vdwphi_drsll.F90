!     ==================================================================
      PROGRAM VDWPHIGEN
!     ==--------------------------------------------------------------==
      USE M_VDW_PHI
      IMPLICIT NONE
#ifdef __PARA
      INCLUDE 'mpif.h'
#endif
! ... Variables for D and delta values
      INTEGER      :: mD=2200
      INTEGER      :: mdel=100
      REAL(KIND=8) :: Dmin=0.D0,Dmax=22.D0
      REAL(KIND=8) :: delmin=0.D0,delmax=1.D0
      INTEGER      :: iD,nD
      INTEGER      :: idel,ndel
      REAL(KIND=8) :: dD,ddel
      REAL(KIND=8) :: D,del
      INTEGER :: idd,ndd
      INTEGER :: idds,idde,nddt
      INTEGER, ALLOCATABLE :: ind(:,:)
      REAL(KIND=8),ALLOCATABLE :: Darr(:),delarr(:)
! ... d1=|r1-r2|q0(r1) and d2=|r1-r2|q0(r2)
      REAL(KIND=8) :: d1,d2
! ... vdW kernel phi
      REAL(KIND=8), ALLOCATABLE :: phivdw(:),phitmp(:),phiout(:,:)
! ... fortran i/o
      INTEGER :: STDIN=05,STDOUT=06,STDERR=00
      INTEGER :: IOUT=20
      CHARACTER(LEN=256) :: FNAME
! ... command line arguments
      INTEGER :: iargc,marg,narg
      INTEGER :: len
      CHARACTER(LEN=256) :: arg
! ... option
      INTEGER :: iplot=0
! ... timing
      REAL(KIND=4) :: ETIME
      REAL(KIND=4) :: TIME(2)
      REAL(KIND=8) :: TTT,T1,T2,WALL,ELAPS, WALL_TIME, ELAPSED_TIME
      INTEGER      :: HOURS,MINUTES,SECONDS
      INTEGER :: VAL(8)
      CHARACTER(LEN=3) :: MON(12)
      DATA MON/'Jan','Feb','Mar','Apr','May','Jun',&
     &         'Jul','Aug','Sep','Oct','Nov','Dec'/
! ... others
      CHARACTER(LEN=4) :: C4
! ... code name
      CHARACTER(LEN=64) :: codnam='vdwphigen_drsll'
#ifdef __PARA
      INTEGER :: npes,mype
      INTEGER :: MPI_ERR
      INTEGER :: ionode=0
      INTEGER :: iproc
      INTEGER, ALLOCATABLE :: indd(:,:)
#endif
!     ==--------------------------------------------------------------==
#ifdef __PARA
      CALL MPI_INIT(MPI_ERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npes,MPI_ERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,MPI_ERR)
#endif
! ... initialize the timer
      CALL DATE_AND_TIME(VALUES=VAL)
      TTT=ETIME(TIME)
      T1=TIME(1)
      T2=TIME(2)
!
! ... for a testing (start) ...
!     Dmin=0.D0
!     Dmax=30.D0
!     mD=300
!     delmin=0.D0
!     delmax=0.9D0
!     mdel=9
!     delmin=0.D0
!     delmax=0.5D0
!     mdel=5
! ... for a testing (start) ...
!
! ... command line arguments
      narg=iargc()
      marg=0
      DO WHILE(marg<narg) 
        marg=marg+1
        CALL getarg(marg,arg)
        len=len_trim(arg) 
        IF(arg(1:len).eq.'-Dmax' .or.&
     &     arg(1:len).eq.'-dmax')THEN
          marg=marg+1
          CALL getarg(marg,arg)
          read(arg,*)Dmax 
        ELSEIF(arg(1:len).eq.'-mD'.or.&
     &     arg(1:len).eq.'-md')THEN
          marg=marg+1
          CALL getarg(marg,arg)
          read(arg,*)mD
        ELSEIF(arg(1:len).eq.'-delmax')THEN
          marg=marg+1
          CALL getarg(marg,arg)
          read(arg,*)delmax 
        ELSEIF(arg(1:len).eq.'-mdel')THEN
          marg=marg+1
          CALL getarg(marg,arg)
          read(arg,*)mdel
        ELSEIF(arg(1:len).eq.'-plot')THEN
          iplot=1 
        ELSEIF(arg(1:len).eq.'-h')THEN
          WRITE(STDERR,'(A,A,A,A)')&
     &    'Usage: ',TRIM(codnam),' [-Dmax][Dmax] [-ddel][ddel] ',&
     &    '[-mD][mD] [-mdel][mdel] [-plot] [-h]'
          STOP
        ENDIF
      ENDDO  
! ... Initialize D values
      Dmax=DBLE(Dmax)
      Dmin=DBLE(Dmin)
      nD=mD+1
      dD=(Dmax-Dmin)/DBLE(mD)
! ... Initialize delta values    
      delmax=DBLE(delmax)
      delmin=DBLE(delmin)
      ndel=mdel+1 
      ddel=(delmax-delmin)/DBLE(mdel)
! ... Initialize arrays for D and delta
      ndd=nD*ndel
      ALLOCATE(Darr(ndd),delarr(ndd))
      ALLOCATE(ind(2,ndd))
      idd=0
      DO idel=1,ndel
        DO iD=1,nD
          del=delmin+DBLE(idel-1)*ddel
          D=Dmin+DBLE(iD-1)*dD
          idd=idd+1
          ind(1,idd)=iD
          ind(2,idd)=idel
          Darr(idd)=D
          delarr(idd)=del 
        ENDDO
      ENDDO
      idds=1
      idde=ndd
#ifdef __PARA
      nddt=ndd/npes      
      idds=(mype)*nddt+1
      idde=(mype+1)*nddt 
      if(mype+1==npes)then
        idde=ndd
      endif
      ALLOCATE(indd(2,0:npes-1))
      do iproc=0,npes-1
        indd(1,iproc)=(iproc)*nddt+1
        indd(2,iproc)=(iproc+1)*nddt 
        if(iproc+1==npes)then
          indd(2,iproc)=ndd
        endif
      enddo
#endif
! ... allocate
      ALLOCATE(phivdw(ndd))
      ALLOCATE(phiout(nd,ndel))
      phivdw(:)=0.D0
      phiout(:,:)=0.D0
#ifdef __PARA
      ALLOCATE(phitmp(ndd))
      phitmp(:)=0.D0
#endif
!     ==--------------------------------------------------------------==
! ... print information
#ifdef __PARA
      IF(mype==ionode)THEN
#endif
      WRITE(STDOUT,'(A)')      ' *********************************' 
      WRITE(STDOUT,'(A)')      ' *** GENERATING THE VDW KERNEL ***'
      WRITE(STDOUT,'(A)')      ' *********************************' 
      WRITE(STDOUT,'(/,2A,I3,I5,I3,2(A,I2.2),/)')&
     &' PROGRAM STARTED AT: ',&
     &MON(VAL(2)),VAL(3),VAL(1),VAL(5),':',VAL(6),':',VAL(7)
#ifdef __PARA
      WRITE(STDOUT,'(A,I5)')   ' NUMBER OF PROCESSORS : ',npes
#endif
      WRITE(STDOUT,'(A)')      '                                  '
      WRITE(STDOUT,'(A)')      ' PARAMETERS USED:                 '
      WRITE(STDOUT,'(A,F12.4)')' Dmin   =',Dmin
      WRITE(STDOUT,'(A,F12.4)')' Dmax   =',Dmax
      WRITE(STDOUT,'(A,I5)')   ' ND     =',nD
      WRITE(STDOUT,'(A,F12.4)')' DELMIN =',delmin
      WRITE(STDOUT,'(A,F12.4)')' DELMAX =',delmax
      WRITE(STDOUT,'(A,I5)')   ' NDEL   =',ndel
      IF(iplot==0)THEN
        WRITE(STDOUT,'(/,A)')' GENERATING KERNEL DATA'
      ELSEIF(iplot==1)THEN
        WRITE(STDOUT,'(/,A)')' GENERATING KERNEL DATA FOR PLOT'
      ENDIF
#ifdef __PARA
      WRITE(STDOUT,'(/,A)')  ' =======================================' 
      WRITE(STDOUT,'(A)')    ' DISTRIBUTION OF D AND DELTA PAIRS'
      WRITE(STDOUT,'(A)')    ' ---------------------------------------' 
      DO iproc=0,npes-1
        WRITE(STDOUT,'(3(A,I8))')&
     &  ' PE=',iproc,' IDDS=',indd(1,iproc),' IDDE=',indd(2,iproc)
      ENDDO
      WRITE(STDOUT,'(A)')    ' =======================================' 
#endif
#ifdef __PARA
      ENDIF
#endif
!     ==--------------------------------------------------------------==
!     == Calculation of the van der Waals kernel phi                  ==
!     ==--------------------------------------------------------------==
      CALL SET_AMESH
      CALL SET_W_AB
      CALL ALLOC_ARR_PHI
      LOOP_DD: DO idd=idds,idde
        iD=ind(1,idd)
        idel=ind(2,idd)
        D=Darr(idd)
        del=delarr(idd) 
! ... set d values
        d1=D*(1.D0+del)
        d2=D*(1.D0-del)
        if(d1.le.SMALL .and. d2.le.SMALL)then
          cycle LOOP_DD
        endif
        phivdw(idd)=phival(d1,d2) 
      ENDDO LOOP_DD
!     ==--------------------------------------------------------------==
!     == Output the van der Waals kernel phi (for plotting)           ==
!     ==--------------------------------------------------------------==
#ifdef __PARA
      CALL MPI_ALLREDUCE(phivdw,phitmp,ndd,MPI_DOUBLE_PRECISION,&
     &                   MPI_SUM,MPI_COMM_WORLD,MPI_ERR) 
      phivdw(:)=phitmp(:)
#endif
      DO idd=1,ndd
        iD=ind(1,idd)
        idel=ind(2,idd)
        phiout(iD,idel)=phivdw(idd)
      ENDDO 
      IF(iplot.eq.0)THEN
#ifdef __PARA
      IF(mype==ionode)THEN
#endif
        FNAME='kernel.dat'
        OPEN(IOUT,FILE=FNAME,STATUS='UNKNOWN')
        WRITE(IOUT,'(E24.5,A)')Dmax,' : Dmax'
        WRITE(IOUT,'(I24,A)')  mD,  ' : Mesh_D'
        WRITE(IOUT,'(I24,A)')  mdel,' : Mesh_delta'
        DO iD=1,nD
          DO idel=1,ndel
            D=Dmin+DBLE(iD-1)*dD
            del=delmin+DBLE(idel-1)*ddel
            WRITE(IOUT,'(2(1x,F5.2),E24.16)')D,del,phiout(iD,idel)
          ENDDO
        ENDDO
        CLOSE(IOUT)
#ifdef __PARA
        ENDIF
#endif
      ELSEIF(iplot.eq.1)THEN
#ifdef __PARA
      IF(mype==ionode)THEN
#endif
! ... output phi is mutiplied by 4*pi*D**2
        FNAME='phi_plot.dat'
        WRITE(STDOUT,'(/,A,A)')' KERNEL IS WRITTEN TO ',TRIM(FNAME)
        OPEN(IOUT,FILE=FNAME,STATUS='UNKNOWN')
        DO idel=1,ndel
          del=delmin+DBLE(idel-1)*ddel
          WRITE(IOUT,'(a,e12.4)')'# delta=',del
          WRITE(IOUT,'(a,a)')&
     &   '#D           ','4*PI*D**2*PHI'
          DO iD=1,nD      
            D=Dmin+DBLE(iD-1)*dD
            d1=D*(1.D0+del)
            d2=D*(1.D0-del)
            WRITE(IOUT,'(2(1x,e12.4))') D, 4.D0*pi*D*D*phiout(iD,idel)
          ENDDO
          WRITE(IOUT,'(/,/)')
        ENDDO
        CLOSE(IOUT)
#ifdef __PARA
        ENDIF
#endif
      ENDIF
!     ==--------------------------------------------------------------==
! ... deallocate
      DEALLOCATE(phiout,phivdw)
      DEALLOCATE(Darr,delarr)
      DEALLOCATE(ind)
#ifdef __PARA
      DEALLOCATE(phitmp)
      DEALLOCATE(indd)
#endif
      CALL DEALLOC_ARR_PHI
! ... timing
      TTT=ETIME(TIME)
      ELAPS=TIME(1)-T1
      WALL=TIME(2)-T2
      WALL=WALL+ELAPS
      WALL_TIME=WALL
      ELAPSED_TIME=ELAPS
#ifdef __PARA
      IF(mype==ionode)THEN
#endif
      HOURS=INT(ELAPS/60.D0/60.D0)
      ELAPS=ELAPS-60.D0*60D0*DBLE(HOURS)
      MINUTES=INT(ELAPS/60.D0)
      ELAPS=ELAPS-60.D0*DBLE(MINUTES)
      SECONDS=INT(ELAPS)
      IF(HOURS.GT.0)THEN
        WRITE(STDOUT,'(/,A,I4,A,I2,A,I2,A)')&
     &  ' ELAPSE TIME : ',HOURS,' HOURS ',MINUTES,' MINUTES ',&
     &  SECONDS,' SECONDS'
      ELSE
        IF(MINUTES.GT.0)THEN
          WRITE(STDOUT,'(A,I2,A,I2,A)')&
     &    ' ELAPSE TIME : ',MINUTES,' MINUTES ',&
     &    SECONDS,' SECONDS'
        ELSE
          WRITE(STDOUT,'(A,F12.4,A)')&
     &    ' ELAPSE TIME : ',ELAPS,' SEC.'
        ENDIF
      ENDIF
      HOURS=INT(WALL/60.D0/60.D0)
      WALL=WALL-60.D0*60D0*DBLE(HOURS)
      MINUTES=INT(WALL/60.D0)
      WALL=WALL-60.D0*DBLE(MINUTES)
      SECONDS=INT(WALL)
      IF(HOURS.GT.0)THEN
        WRITE(STDOUT,'(A,I4,A,I2,A,I2,A)')&
     &  ' WALL TIME   : ',HOURS,' HOURS ',MINUTES,' MINUTES ',&
     &  SECONDS,' SECONDS'
      ELSE
        IF(MINUTES.GT.0)THEN
          WRITE(STDOUT,'(A,I2,A,I2,A)')&
     &  ' WALL TIME   : ',MINUTES,' MINUTES ',&
     &    SECONDS,' SECONDS'
        ELSE
          WRITE(STDOUT,'(A,F12.4,A)')&
     &  ' WALL TIME   : ',ELAPS,' SEC.'
        ENDIF
      ENDIF
      CALL DATE_AND_TIME(VALUES=VAL)
      WRITE(STDOUT,'(/,2A,I3,I5,I3,2(A,I2.2))')&                    
     &' PROGRAM ENDED AT: ',&
     &MON(VAL(2)),VAL(3),VAL(1),VAL(5),':',VAL(6),':',VAL(7)
#ifdef __PARA
      ENDIF
#endif
#ifdef __PARA
      CALL MPI_FINALIZE(MPI_ERR) 
#endif
!     ==--------------------------------------------------------------==
      STOP
      END PROGRAM VDWPHIGEN
!     ==================================================================
