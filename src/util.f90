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
     PRINT*, "ERR: only 2D and 3D cases of Lieb models implemented --- ABORTING"
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
     ao(ElementsPut-1) = 0.0D0   !*(DRANDOM(ISeed)-0.5D0) ! Store the onsite energy
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
           CALL RANDOM_NUMBER(ao(ElementsPut))
           ao(ElementsPut) = ao(ElementsPut) + 0.5D0    ! Store the hopping energy
           jao(ElementsPut-1) = ap(indhop)   ! Store the position of hopping energy
        END IF
     END DO
     
     ! ----------------------------- For rim atoms --------------------------!
 
     IF(Nx>1)THEN
!!$        DO k=1, Nx-1
!!$           DO j=1, Dim
        DO j=1, Dim
           ! For rim atoms except close to hub atom of other unit cell(inner rim atoms)
           DO k=1, Nx-1
              
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
                    CALL RANDOM_NUMBER(ao(ElementsPut))
                    ao(ElementsPut) = ao(ElementsPut) + 0.5D0    ! Store the hopping energy
                    jao(ElementsPut-1) = bp(indhop)   ! Store the position of hopping energy
                 END IF
              END DO
              
           END DO
           
           ! For rim atoms close to hub atoms in other unit cells(outer rim atoms) 
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
              !           bp(2)= MOD((i-1)*ucl + ucl*(M)**(j-1) +1,ucl*(M)**j)
              bp(2)= (i-1)*ucl + 1 - (M-1)*ucl*(M)**(j-1)
           ELSE
              !           bp(2)= (i-1)*ucl + 1 + ucl*(M)**(j-1)
              bp(2)= (i-1)*ucl + 1 + ucl*(M)**(j-1)
           END IF

           iao(ind-1) = ElementsPut    ! Position of diagonal element in the matrix
           ao(ElementsPut-1) = 0.0D0 !*(DRANDOM(ISeed)-0.5D0) ! Store the onsite energy
           jao(ElementsPut-1) = bp(0)  ! Store the position of onsite energy 

           DO indhop=1,2
              IF(bp(indhop).GT.ind)THEN
                 ElementsPut = ElementsPut + 1 
                 CALL RANDOM_NUMBER(ao(ElementsPut))
                 ao(ElementsPut) = ao(ElementsPut) + 0.5D0    ! Store the hopping energy
                 jao(ElementsPut-1) = bp(indhop)   ! Store the position of hopping energy
              END IF
           END DO
           
        END DO
        
     Else if(Nx==1) Then
    ! For rim atoms close to hub atoms in other unit cells          
        Do j=1, Dim
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
              !           bp(2)= MOD((i-1)*ucl + ucl*(M)**(j-1) +1,ucl*(M)**j)
              bp(2)= (i-1)*ucl + 1 - (M-1)*ucl*(M)**(j-1)
           ELSE
              !           bp(2)= (i-1)*ucl + 1 + ucl*(M)**(j-1)
              bp(2)= (i-1)*ucl + 1 + ucl*(M)**(j-1)
           END IF

           iao(ind-1) = ElementsPut    ! Position of diagonal element in the matrix
           ao(ElementsPut-1) = 0.0D0 !*(DRANDOM(ISeed)-0.5D0) ! Store the onsite energy
           jao(ElementsPut-1) = bp(0)  ! Store the position of onsite energy 

           DO indhop=1,2
              IF(bp(indhop).GT.ind)THEN
                 ElementsPut = ElementsPut + 1 
                 CALL RANDOM_NUMBER(ao(ElementsPut))
                 ao(ElementsPut) = ao(ElementsPut) + 0.5D0    ! Store the hopping energy
                 jao(ElementsPut-1) = bp(indhop)   ! Store the position of hopping energy
              END IF
           END DO

        END DO

     END IF
     

     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     
  END DO
  
  iao(LSize)=ElementsPut+1
  
  nz=ElementsPut
  
  RETURN
     
  
END SUBROUTINE MAKECOMPACTROWLIEBMAT



!!$SUBROUTINE DetEnergy(Dim, Nx, IWidth, HubDis, Energy)
!!$
!!$  INTEGER(KIND=IKIND) Dim, Nx, IWidth
!!$
!!$  REAL(KIND=RKIND) HubDis, Energy
!!$
!!$  Energy = HubDis/2 + 1.0D0 ! Potensial + hopping
!!$
!!$  RETURN
!!$
!!$END SUBROUTINE DETENERGY

 
