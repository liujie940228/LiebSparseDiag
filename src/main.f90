!
!
!   Call JADAMILU PJD() to calculate eigenvalues CLOSE TO FLAT ENERGY 0.0 and eigenvectors
!   for Lieb matrix and its extendsions(2D and 3D)
!
!
!--------------------------------------------------------------------------------------
PROGRAM MAIN

  use, intrinsic :: iso_c_binding

  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)

  ! Parameters for Lieb matrix
  INTEGER(KIND=IKIND) &
       Dim, & !the dimension of Lieb model
       Nx, &  !the number of insert atoms in the edge
       M,&    !the unmber of unit cell in each dimension
       nz

  INTEGER(KIND=IKIND) i, j, k, LSize,CSize, seed
  REAL(KIND=RKIND) x


  PARAMETER( Dim = 3 )
  PARAMETER( M = 5 )
  PARAMETER( Nx = 1 )
  PARAMETER( LSize = (Dim*Nx+1)*(M**Dim) )
  PARAMETER( CSize=(2*Dim*Nx+4)*(M**Dim) )


  REAL(KIND=RKIND)  HubDiagDis, RimDiagDis
  PARAMETER ( HubDiagDis = 1.0, RimDiagDis = 0.00000 )
  INTEGER(KIND=IKIND), DIMENSION(:), ALLOCATABLE:: iao, jao, ia, ja
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: ao, a

  ! Parameters for clock

  INTEGER(KIND=IKIND) t1,t2,clock_rate,clock_max
  REAL(KIND=RKIND) :: time_CPU_before

  ! arguments to pass to the JD routine
  INTEGER(KIND=IKIND) &
       NEIG, MADSPACE, ISEARCH, NINIT, ICNTL(5), &
       ITER, IPRINT, INFO, IJOB, NDX1, NDX2, IErr, sumIErr, maxsp, &
       VECS_size, &       ! optimal workspace (N.B.: here, maxsp*maxsp>maxeig)
       NEVALS
  PARAMETER (NEVALS = 10)

  PARAMETER (maxsp=20)

  REAL(KIND=RKIND) &
       SIGMA, TOL, DELTA, SHIFT, GAP, MEM, DROPTOL
  REAL(KIND=RKIND) Energy
  PARAMETER (Energy = 0.0D0 )

  ! arguments from the JD routine
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: &
       EIGS, & !computed eigenvalues
       RES,&   !computed residues: |Ax-EIGS x|
       VECS    !computed eigenvectors

  ! ----------------------------------------------------------
  ! defining the size of arrays
  ! ----------------------------------------------------------
  
!!$  Lsize     = (Dim*Nx+1)*M**Dim
  VECS_size = Lsize*(3*maxsp+NEVals+1)+4*maxsp*maxsp

  ALLOCATE( EIGS(NEVals) )
  ALLOCATE( RES(NEVals) )
  ALLOCATE( VECS(VECS_size) )
  
       
  ALLOCATE( iao(0:LSize))
  ALLOCATE( jao(0:CSize))
  ALLOCATE( ao(0:CSize) )

  ALLOCATE( ia(LSize+1) )
  ALLOCATE( ja(CSize) )
  ALLOCATE( a(CSize) )
  
  call init_random_seed()
  call srand(seed)

  ! Start timing
  CALL system_clock (t1, clock_rate, clock_max)

!!$  CALL MakeLiebDNxMat(Dim, Nx, M, LSize, CSize, iao, jao, ao, nz, HubDiagDis, RimDiagDis )
  CALL MakeCompactRowLiebMat(Dim, Nx, M, LSize, CSize, iao, jao, ao, nz )
  
  DO i=1,nz
     ja(i)=jao(i-1)
!!$     PRINT*,i,ja(i)
  ENDDO
  
  DO i=1,nz
     a(i)=ao(i-1)
!!$     PRINT*,i,a(i)
  ENDDO
  
  DO i=1,Lsize+1
     ia(i)=iao(i-1)
!!$     PRINT*,i,ia(i)
  ENDDO

  DO i=1, M**Dim
     CALL Random_number(x)
     k= (i-1)*(Nx*Dim+1) + 1
     a(ia(k)) = HubDiagDis*(x - 0.5D0)

     DO j=2, (Nx*Dim +1)
        CALL Random_number(x)
        k = (i-1)*(Nx*Dim+1) + j
        a(ia(k)) = RimDiagDis*(x - 0.5D0)
     END DO

  END DO

!!$  Print*,"------------------"
!!$  DO i=1,nz
!!$     PRINT*,i,a(i)
!!$  ENDDO
  
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
  MEM=20.0
  ! tolerence for discarded fill
  DROPTOL=1.d-3

  ! number of wanted eigenvalues
  NEIG=NEVals

  ! no initial approximate eigenvectors
  NINIT=0
  ! desired size of the search space
  MADSPACE=maxsp
  ! maximum number of iteration steps
  ITER=1000
  ! tolerance for the eigenvector residual
  TOL=1.0d-10

  ! additional parameters set to default
  ICNTL(1)=0
  ICNTL(2)=0    ! switch ON for ADAPTIVE precon =0
  ICNTL(3)=0
  ICNTL(4)=0
  ICNTL(5)=0

  ! ----------------------------------------------------------
  !  call to PJD that computes eigenvalues & eigenvectors

  CALL PJD(Lsize, a, ja, ia, EIGS, RES, VECS, VECS_size, NEIG,&
       SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL,&
       SHIFT, DROPTOL, MEM, ICNTL,&
       IPRINT, INFO, GAP)


  !  CALL PJDCLEANUP   !to be used if you want a new preconditioner in every iteration

  ! ----------------------------------------------------------

  DO i=1, NEIG
     PRINT*, i, EIGS(i)
  END DO
  
  DEALLOCATE(EIGS,RES,VECS,ia,iao,ja,jao,a,ao)
  
  CALL cpu_time(time_CPU_before)
  
END PROGRAM MAIN

  
SUBROUTINE MakeCompactRowLiebMat(Dim, Nx, M, LSize, CSize, iao, jao, ao, nz )
    
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)
  
  INTEGER(KIND=IKIND) &
       Dim, & !the dimension of Lieb model
       Nx, &  !the number of insert atoms in the edge
       M, &   !the unmber of unit cell in each dimension
       LSize, & !the whole number of atoms in system
       CSize, & !the number of onsite potential and interaction terms
       ucl, & !the number of atoms in a unit cell
       n_uc, & ! the number of unit cell
       nz
  INTEGER(KIND=IKIND) i, j, k, ind, indhop, ElementsPut
  INTEGER(KIND=IKIND) ap(0:2*Dim), bp(0:2)
  
  INTEGER(KIND=IKIND), DIMENSION(:), ALLOCATABLE::ucl_d
!!$  REAL(KIND=RKIND)  HubDiagDis, RimDiagDis
  INTEGER(KIND=IKIND)   iao(0:LSize), jao(0:CSize)
  REAL(KIND=RKIND)  ao(0:CSize)

  LOGICAL(KIND=4) Flag

  n_uc = M**Dim
  ucl = Dim*Nx + 1

  
  IF(Dim==2)THEN
     ALLOCATE(ucl_d(Dim))
     ucl_d(1) = 2          ! The first Rim atoms
     ucl_d(2) = Nx + 2
  ELSE IF(Dim==3)THEN
     ALLOCATE(ucl_d(Dim))
     ucl_d(1) = 2
     ucl_d(2) = Nx + 2
     ucl_d(3) = 2 * Nx + 2
  ELSE
     Print*, "We Only Finished the 2D and 3D cases for Lieb model"
     STOP
  END IF

  !
  ElementsPut = 0
  
  DO i=1, n_uc
     ! --------------------------- For hub atoms -------------------------! 
     ind = (i-1)*ucl + 1
     ElementsPut = ElementsPut + 1
     
     ap(0) = ind
     !  Onsite energy term
     iao(ind-1) = ElementsPut    ! Position of diagonal element in the matrix
     ao(ElementsPut-1) = 0.0D0 !*(DRANDOM(ISeed)-0.5D0) ! Store the onsite energy
     jao(ElementsPut-1) = ap(0)  ! Store the position of onsite energy

     !hopping term
     DO j=1, Dim
        IF(j==1)THEN
           Flag = ( MOD(i,M)==1 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,M**2)>0 .AND. MOD(i,M**2)<=M )
        ELSE
           Flag = ( i<=M**2 )
        END IF
        
        ap(j) = (i-1)*ucl + ucl_d(j)
        IF(Flag)THEN
           ap(j+Dim) = (i-1)*ucl + ucl_d(j) + (Nx-1) - ucl*(M)**(j-1) + ucl*(M)**j
        ELSE
           ap(j+Dim) = (i-1)*ucl + ucl_d(j) + (Nx-1) - ucl*(M)**(j-1)
        END IF
     END DO

     DO indhop=1,2*Dim
        IF(ap(indhop).GT.ind)THEN
           ElementsPut = ElementsPut + 1 
           ao(ElementsPut-1) = 1.0D0    ! Store the hopping energy
           jao(ElementsPut-1) = ap(indhop)   ! Store the position of hopping energy
        END IF
     END DO
    
     
     ! ----------------------------- For rim atoms --------------------------!

     ! For rim atoms except close to hub atom of other unit cell 
     IF(Nx>1)THEN
        DO k=1, Nx-1
           
           DO j=1, Dim
              
              ind = (i-1)*ucl + ucl_d(j) + k - 1
              ElementsPut = ElementsPut + 1
              !  Onsite energy term              
              bp(0) = ind
              !Hopping term
              IF(k==1)THEN
                 bp(1) = ind - (j-1)*Nx -1
              ELSE
                 bp(1) = ind -1
              END IF
              bp(2) = ind + 1
              
              iao(ind-1) = ElementsPut    ! Position of diagonal element in the matrix
              ao(ElementsPut-1) = 0.0D0 !*(DRANDOM(ISeed)-0.5D0) ! Store the onsite energy
              jao(ElementsPut-1) = bp(0)  ! Store the position of onsite energy
              
              DO indhop=1,2
                 IF(bp(indhop).GT.ind)THEN
                    ElementsPut = ElementsPut + 1 
                    ao(ElementsPut-1) = 1.0D0    ! Store the hopping energy
                    jao(ElementsPut-1) = bp(indhop)   ! Store the position of hopping energy
                 END IF
              END DO
              
           END DO
           
        END DO
        
     END IF
     
    ! For rim atoms close to hub atoms in other unit cells  

     DO j=1, Dim
        IF(j==1)THEN
           Flag = ( MOD(i,M**j)==0 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,M**j)==0 .OR. (MOD(i,M**j).GT.(M**j-M)) )
        ELSE
           Flag = ( i.GT.(M**j-M**2) )
        END IF
        
        ind = (i-1)*ucl + ucl_d(j) + Nx - 1
        ElementsPut = ElementsPut + 1
        !  Onsite energy term
        bp(0) = ind
        !  hopping term
        IF(Nx==1) THEN
           bp(1) = ind - ucl_d(j) +1
        ELSE
           bp(1) = ind - 1
        END IF
        
        IF(Flag)THEN
!!$           bp(2)= MOD((i-1)*ucl + ucl*(M)**(j-1) +1,ucl*(M)**j)
           bp(2)= (i-1)*ucl + 1 - (M-1)*ucl*(M)**(j-1)
        ELSE
!!$           bp(2)= (i-1)*ucl + 1 + ucl*(M)**(j-1)
           bp(2)= (i-1)*ucl + 1 + ucl*(M)**(j-1)
        END IF
        
        iao(ind-1) = ElementsPut    ! Position of diagonal element in the matrix
        ao(ElementsPut-1) = 0.0D0 !*(DRANDOM(ISeed)-0.5D0) ! Store the onsite energy
        jao(ElementsPut-1) = bp(0)  ! Store the position of onsite energy 
        
        DO indhop=1,2
           IF(bp(indhop).GT.ind)THEN
              ElementsPut = ElementsPut + 1 
              ao(ElementsPut-1) = 1.0D0    ! Store the hopping energy
              jao(ElementsPut-1) = bp(indhop)   ! Store the position of hopping energy
           END IF
        END DO
        
     END DO
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     
  END DO
  
  iao(LSize)=ElementsPut+1
  
  nz=ElementsPut
  
  RETURN
     
  
END SUBROUTINE MAKECOMPACTROWLIEBMAT

!!$SUBROUTINE WriteEvals(dm, nu, n, nt, HubDiagDis, RimDiagDis, W, matr_W, norm, part_nr, ISample, INFO)
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
!!$  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)
!!$
!!$  INTEGER(KIND=IKIND) dm, nu, n, nt, ISample, INFO
!!$  INTEGER(KIND=IKIND) i, j
!!$  REAL(KIND=RKIND) HubDiagDis, RimDiagDis
!!$  REAL(KIND=RKIND) W( nt ), matr_W( nt, nt ), norm( nt ), part_nr( nt )
!!$
!!$  CHARACTER*100 FileName
!!$
!!$  WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A1,A2,I4.4,A1,A2,I4.4,A1,I4.4,A4)')&
!!$       "Eval-","L",dm,nu,&
!!$       "-M",n, "-",&
!!$       "WH", NINT(100.D0*ABS(HubDiagDis)),&
!!$       "-","WR", NINT(100.D0*ABS(RimDiagDis)),&
!!$       "-",ISample,".raw"
!!$
!!$  Print*, "FileName: ", FileName
!!$
!!$  OPEN(Unit=9, FILE=FileName)
!!$
!!$  IF(INFO==0)THEN
!!$
!!$     DO i=1,nt
!!$        DO j=1,nt	
!!$           norm(i) = norm(i) + matr_W(i,j)**2
!!$           part_nr(i) = part_nr(i) + matr_W(i,j)**4
!!$        END DO
!!$     END DO
!!$
!!$     DO i=1,nt
!!$        write(9,'(2f30.20)') W(i) , (norm(i)**2) / part_nr(i)
!!$     END DO
!!$
!!$  ELSE
!!$     PRINT*, "ERROR IN CALL DSYEV()" 
!!$  END IF
!!$
!!$  CLOSE(9)
!!$
!!$  RETURN
!!$
!!$END SUBROUTINE WRITEEVALS
              

Subroutine init_random_seed()

  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)

END Subroutine init_random_seed

            
