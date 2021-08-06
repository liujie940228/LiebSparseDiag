!-----------------------------------------------------------
!
!     LiebJADdiag
!
!     program uses JADAMILU package of Bollhoefer           
!                  roemer's makeanderson matrix             
!                           converted from base 0 to base 1 
! ----------------------------------------------------------
!
! ----------------------------------------------------------

PROGRAM LiebJADdia

  USE MyNumbers  
  USE CConstants
  USE IConstants
  USE IPara
  USE DPara
  USE IChannels
  USE RNG
  USE SETBINS
  !USE RConstants					   
  !USE EigenPara
  !USE IEigenChannels
  !USE SimPara

  IMPLICIT NONE

  ! ----------------------------------------------------------
  ! variable declaration
  ! ----------------------------------------------------------

  ! Parameters for Lieb matrix
  INTEGER(KIND=IKIND) IWidth
  
  INTEGER(KIND=IKIND) i, j, k, LSize, CSize
       
  ! arguments to pass to the JD routine
  INTEGER(KIND=IKIND) &
       NEIG, MADSPACE, ISEARCH, NINIT, ICNTL(5), &
       ITER, IPRINT, INFO, IJOB, NDX1, NDX2, IErr, sumIErr, maxsp, &
       VECS_size       ! optimal workspace (N.B.: here, maxsp*maxsp>maxeig)

  PARAMETER (maxsp=20)

  REAL(KIND=RKIND) &
       SIGMA, TOL, DELTA, SHIFT, GAP, MEM, DROPTOL

  ! output arguments from the Lieb routine
  INTEGER(KIND=IKIND) nz
  
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: &
       ao

  INTEGER(KIND=IKIND), DIMENSION(:), ALLOCATABLE:: &
       iao,jao

  ! arguments from the JD routine
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: &
       EIGS, & !computed eigenvalues
       RES,&   !computed residues: |Ax-EIGS x|
       VECS, &    !computed eigenvectors
       a,a_w

  INTEGER(KIND=IKIND), DIMENSION(:), ALLOCATABLE:: &
       ia, ja

  ! some local variables
  INTEGER(KIND=IKIND)  i1, i2, i3, Inum

  REAL(KIND=RKIND) SugTarE
  Character*100 str

!  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE::mat

!  Character*68 Matrixname

!!$  ! Setting Target Energy
!!$
!!$  INTEGER(KIND=IKIND) L
!!$
!!$  INTEGER(KIND=IKIND), Dimension(:,:), ALLOCATABLE:: TarStore

  ! ----------------------------------------------------------
  ! start of main code
  ! ----------------------------------------------------------

  ! ----------------------------------------------------------
  ! protocol feature
  ! ----------------------------------------------------------

  PRINT*,"LiebSparseDiag ", RStr, DStr, AStr 

  ! ----------------------------------------------------------
  ! inout handling
  ! ----------------------------------------------------------

  CALL  Input(IErr)
  IF(IErr.NE.0) THEN
     PRINT*,"main: Input() finds IErr=", IErr
     STOP
  ENDIF

  ! ----------------------------------------------------------
  ! start of main IWidth loop
  ! ----------------------------------------------------------

  DO IWidth= Width0, Width1, dWidth

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.0) THEN
        PRINT*, "START@ IWidth=", IWidth
     ENDIF
     
     ! ----------------------------------------------------------
     ! defining the size of arrays
     ! ----------------------------------------------------------
     
     LSize     = (Dim*Nx+1)*(IWidth**Dim)
     VECS_size = Lsize*(3*maxsp+NEVals+1)+4*maxsp*maxsp
     CSize=(2*Dim*Nx+4)*(IWidth**Dim) 
     
     IF(IWriteFlag.GE.2) THEN
        PRINT*,"main: cube is ",Lsize," by ",Lsize
        PRINT*,"main: VEC size = ",VECS_size
     ENDIF
     
     ! ----------------------------------------------------------
     ! ALLOCATing memory
     ! ----------------------------------------------------------
     
     sumIErr= 0
     ALLOCATE(EIGS(NEVals), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(RES(NEVals), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(VECS(VECS_size), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(ao(0:CSize), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(iao(0:LSize), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(jao(0:CSize), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(a(CSize), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(a_w(CSize), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(ia(LSize+1), STAT=IErr); sumIErr=sumIErr+IErr
     ALLOCATE(ja(CSize), STAT=IErr); sumIErr=sumIErr+IErr
     
     IF(IErr.NE.0) THEN
        PRINT*,"main: Input() finds sumIErr=", sumIErr
        STOP
     ENDIF
     

     ! ----------------------------------------------------------
     ! making the Lieb matrix
     ! ----------------------------------------------------------
     
     CALL MakeCompactRowLiebMat(Dim, Nx, IWidth, LSize, CSize, iao, jao, ao, nz )
     !output NZ = # of nonzero and diagonal elements
     !output ia must have size of N+1


     IF(IWriteFlag.GE.4) PRINT*,&
          "main: There are",nz,"nonzero elements of array a()"

     ! ----------------------------------------------------------
     ! converting arrays from base 0 to base 1
     ! ----------------------------------------------------------
     
     DO i=1,nz
        ja(i)=jao(i-1)
!!$           IF(IWriteFlag.GE.4) PRINT*,i,ja(i)
     ENDDO

     DO i=1,nz
        a(i)=ao(i-1)
!!$           IF(IWriteFlag.GE.4) PRINT*,i,a(i)
     ENDDO

     DO i=1,Lsize+1
        ia(i)=iao(i-1)
!!$           IF(IWriteFlag.GE.4) PRINT*,i,ia(i)
     ENDDO
     
     
     ! -----------------------------------------------------------------
     ! start of HubDis loop, not finish yet, just keep HubDis0= HubDis1
     ! -----------------------------------------------------------------
       
     DO HubDis= HubDis0, HubDis1, dHubDis
        
     
        ! ----------------------------------------------------------
        ! start of ISeed loop
        ! ----------------------------------------------------------
        
        
        DO Seed= ISeed, ISeed+ NSeed -1

           Call GetDirec(Dim, Nx, IWidth, HubDis, RimDis, Seed, str)

           IF(IWriteFlag.GE.1) THEN
              PRINT*, "  HubDis=", HubDis, " Seed=", Seed
           ENDIF


           CALL SRANDOM(Seed)

           ! keep array a Lieb matrix form, for each disorder circle, only change the a_w
           a_w(:) = a(:) 

           ! Give the Lieb matrix different onsite potensial
           DO i=1, IWidth**Dim

              k= (i-1)*(Nx*Dim+1) + 1
              a_w(ia(k)) = HubDis*(DRANDOM(Seed) - 0.5D0)


              DO j=2, (Nx*Dim +1)

                 k = (i-1)*(Nx*Dim+1) + j
                 a_w(ia(k)) = RimDis*(DRANDOM(Seed) - 0.5D0)

              END DO

           END DO

           !---------------------------------------------------------
           ! Transform the format of Sparse matrix to Full matrix
           !--------------------------------------------------------
           
!!$           mat(:,:)=0.0D0
!!$           DO i=1, LSize ! the horizontal of matrix
!!$
!!$              DO j=ia(i),ia(i+1)-1
!!$
!!$                 mat(i,ja(j)) = a_w(j)
!!$
!!$              END DO
!!$
!!$           END DO
!!$
!!$           ! Generate the file containing full matrix 
!!$
!!$           WRITE(Matrixname, '(A11,I1,I1,A2,I4.4,A6,I4,A4)') &
!!$                "SparseMat-L", Dim, Nx, "-M", IWidth,  &
!!$                "-Seed-", Seed, &
!!$                ".txt"
!!$
!!$           OPEN(UNIT= 10, FILE=Matrixname)   
!!$
!!$           DO i=1, LSize
!!$              DO j=1, LSize
!!$                 WRITE(10,102) mat(i,j)              
!!$              END DO
!!$              WRITE(10,*) " "
!!$           END DO
!!$
!!$102        FORMAT(1x,F15.6\)
!!$
!!$           CLOSE(10)
           

           DO Energy= Energy0, Energy1, dEnergy
           
   
              ! ----------------------------------------------------------
              ! interface to the JADAMILU code
              ! ----------------------------------------------------------
              ! the matrix is already in the required format

              ! standard report on standard output
              IPRINT=6

              ! we want the eigenvalues near target sigma=0         
              ! (then, SHIFT need not to be set)       
              ISEARCH=2      
              SIGMA=Energy
                             
              ! elbow space factor for the fill computed during the ILU
              MEM=Memory
              ! tolerence for discarded fill
              DROPTOL=1.d-3

              ! number of wanted eigenvalues
              NEIG=NEVals

              ! no initial approximate eigenvectors
              NINIT=0
              ! desired size of the search space
              MADSPACE=maxsp
              ! maximum number of iteration steps
              ITER=20000
              ! tolerance for the eigenvector residual
              TOL=1.0d-10

              ! additional parameters set to default
              ICNTL(1)=0
              ICNTL(2)=0    ! switch ON for ADAPTIVE precon =0
              ICNTL(3)=0
              ICNTL(4)=0
              ICNTL(5)=0

!!$              RES=0
!!$              VECS=0
!!$              EIGS=0

              ! ----------------------------------------------------------
              !  call to PJD that computes eigenvalues & eigenvectors

              IF(IWriteFlag.GE.4) PRINT*,"main: calling PJD()"
              CALL PJD(Lsize, a_w, ja, ia, EIGS, RES, VECS, VECS_size, NEIG,&
                   SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL,&
                   SHIFT, DROPTOL, MEM, ICNTL,&
                   IPRINT, INFO, GAP)

              ! We set a large NEvals to guarantee get all eigenvalues at a certain droptol
              ! When it get the number of eigenvalues is less than NEvals which we set,
              ! it will return the actually number INFO(which included one unconverged)
              IF(INFO.NE.0) THEN

                 NEIG= INFO - 1

              END IF

              CALL PJDCLEANUP   !to be used if you want a new preconditioner in every iteration

              ! ----------------------------------------------------------
              ! write results into files
              ! ---------------------------------------------------------
             
              
              IF(NEIG==0)THEN
                 Print*,"Don't find any eigenvalues!"
              ELSE IF(NEIG.lt.0)THEN
                 Print*,"Error: may (D)SYGV/(Z)HEGV .... "
              ELSE
                 DO i=1, NEIG
                    PRINT*, i, EIGS(i)
                 END DO
              END IF
              
              SELECT CASE(IKeepFlag)

              CASE(0)

                 CALL WriteOutputEVal( Dim, Nx, NEIG, EIGS, IWidth, Energy, HubDis, RimDis, Seed, str, IErr)

!!$                 DO Inum=1, NEVals
!!$                    CALL WriteOutputEVec( Dim, Nx, Inum, NEIG, Lsize, VECS, VECS_size, &
!!$                         IWidth, Energy, HubDis, RimDis, Seed, str, IErr)
!!$                 END DO

              CASE(1)           
                 CALL CheckOutput( IWidth, Energy, HubDis, RimDis, Seed, IErr )
                 IF(IErr.EQ.2) GOTO 100

                 CALL WriteOutputEVal(NEIG, EIGS, IWidth, Energy, HubDis, RimDis, Seed, IErr, str, IKeepFlag)
                 DO Inum=1,NEVals
                    CALL WriteOutputEVec( Inum, NEIG, Lsize, VECS, VECS_size, &
                         IWidth, Energy, HubDis, RimDis, Seed, str, IErr)

                 END DO

100           END SELECT                                  
             
           END DO !Energy loop

        END DO !Seed loop

     END DO !HubDis loop
     


     
     ! ----------------------------------------------------------
     ! DEALLOCATE memory
     
     DEALLOCATE(EIGS,RES,VECS,ia,iao,ja,jao,a,a_w,ao,STAT=IErr)
     
     IF(IErr.NE.0) THEN
        PRINT*,"main: DEALLOCATE() finds IErr=", IErr
        !STOP
     ENDIF

  END DO ! Width loop
  
END PROGRAM LiebJADdia   !main program






