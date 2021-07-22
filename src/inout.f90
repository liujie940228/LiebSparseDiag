! ********************************************************************
!       
! TMSE2D - Transfer matrix method for the Anderson
! model with diagonal disorder in two dimensions
!
! ********************************************************************
      
! ********************************************************************
!       
! $Header: /home/cvs/phsht/AML/src/inout.f90,v 1.8 2007/09/26 12:28:39 phrfar Exp $
!
! ********************************************************************

! **************************************************************************
! $Log: inout.f90,v $
! Revision 1.8  2007/09/26 12:28:39  phrfar
! removed if-then for ikeep flag, used instead "status=UNKNOWN"
! replaced "X" & "LX" by "VECS" & "VECS_SIZE" respectively
!
! Revision 1.7  2007/09/26 09:57:35  ccspab
! WriteOutputEval error removed when using Ikeepflag=1
!
! Revision 1.6  2007/09/25 20:49:35  phrfar
! added subroutines WriteOutputEVal and WriteOutputEVec that write eigenvalues and vectors to a file
!
! Revision 1.5  2007/09/25 12:41:46  phsht
! 1st version of Checkoutput()
!
! Revision 1.4  2007/09/25 10:30:33  phsht
! rewrote to make "look" nicer and removed an ERROR for NEvals and Energy; now makefile has proper flags to check for these things automatically
!
! Revision 1.3  2007/09/25 09:25:59  phrfar
! committing new .f90 files
!
! Revision 1.2  2007/09/21 15:46:05  phrfar
! input subroutine for main.f90(jadamilu)
!
! **************************************************************************

! --------------------------------------------------------------------
! Input:
!
! IErr	error code
!----------------------------------------------------------------------

SUBROUTINE Input(IErr)

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels

  
  INTEGER IErr, ILine
!  REAL(KIND=RKIND) RIter
  
  !	PRINT*,"DBG: Input()"
  
  IErr = 0
  ILine= 0
  
  OPEN(UNIT= IChInp, ERR=120, FILE= "LiebSparsediag.inp",STATUS= 'OLD')

  ILine= ILine+1
  READ(IChInp,10,ERR=20) ISeed
  !PRINT*,"ISeed        = ",ISeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) NSeed
  !PRINT*,"NSeed        = ",NSeed
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Dim
  !PRINT*,"Dim         = ",Dim

  ILine= ILine+1
  READ(IChInp,10,ERR=20) Nx
  !PRINT*,"Nx          = ",Nx
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IBCFlag
  !PRINT*,"IBCFlag      = ", IBCFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IRNGFlag
  !PRINT*,"IRNGFlag     = ", IRNGFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IKeepFlag
  !PRINT*,"IKeepFlag    = ", IKeepFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IWriteFlag
  !PRINT*,"IWriteFlag   = ", IWriteFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width0
  !PRINT*,"Width0       = ",Width0
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width1
  !PRINT*,"Width1       = ", Width1
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) dWidth
  !PRINT*,"dWidth       = ", dWidth

  ILine= ILine+1
  READ(IChInp,15,ERR=20) HubDis0
  !PRINT*,"HubDis0      = ", HubDis0

  ILine= ILine+1
  READ(IChInp,15,ERR=20) HubDis1
  !PRINT*,"HubDis1      = ", HubDis1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dHubDis
  !PRINT*,"dHubDis       = ", dHubDis
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) RimDis
  !PRINT*,"RimDis       = ", RimDis
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) Energy0
  !PRINT*,"Energy0      = ", Energy0

  ILine= ILine+1
  READ(IChInp,15,ERR=20) Energy1
  !PRINT*,"Energy1      = ", Energy1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dEnergy
  !PRINT*,"dEnergy      = ", dEnergy

  ILine= ILine+1
  READ(IChInp,10,ERR=20) NEVals
  !PRINT*,"NEVals      = ",NEVals
  
  
10 FORMAT(16X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(16X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)

  ! check the parameters for validity
  
  IF(IWriteFlag.GE.2) THEN
     PRINT*,"ISeed        = ", ISeed
     PRINT*,"NSeed        = ", NSeed
     PRINT*,"Dim          = ", Dim
     PRINT*,"Nx           = ", Nx
     PRINT*,"IBCFlag      = ", IBCFlag
     PRINT*,"IRNGFlag     = ", IRNGFlag
     PRINT*,"IKeepFlag    = ", IKeepFlag
     PRINT*,"IWriteFlag   = ", IWriteFlag
     PRINT*,"Width0       = ", Width0
     PRINT*,"Width1       = ", Width1
     PRINT*,"dWidth       = ", dWidth
     PRINT*,"HubDis0      = ", HubDis0
     PRINT*,"HubDis1      = ", HubDis1
     PRINT*,"dHubDis      = ", dHubDis
     Print*,"RimDis       = ", RimDis
     PRINT*,"Energy0      = ", Energy0
     PRINT*,"Energy1      = ", Energy1
     PRINT*,"dEnergy      = ", dEnergy
     PRINT*,"NEVals       = ", NEVals
    
  ENDIF

  CLOSE(IChInp)
  RETURN

120 PRINT*,"Input(): ERR in OPEN()"  
  IErr= 1
  RETURN

20 PRINT*,"Input(): ERR in WRITE()"  
  IErr= 1
  RETURN
  
END SUBROUTINE Input

! --------------------------------------------------------------------
! Input:
!
! IErr	error code
!----------------------------------------------------------------------

FUNCTION GetFileName(fnamestr,vdata,n)
!write(*,*) trim(GetFileName('/Output/QH_L(I4)_NL(I1)_NS(I1)_B(F5.3)_C(F4.2)_S(I5).dat',&
!     (/L_INPUT,NLevels,NSpins,BField,Coul,Seed/),6))


  IMPLICIT NONE

  INTEGER,INTENT(in)         :: n
  REAL(8),INTENT(in)         :: vdata(n)
  CHARACTER(len=*),INTENT(in) :: fnamestr

  CHARACTER(len=*)   :: GetFileName
  CHARACTER(len=500) :: str,fstr,vstr
  
  INTEGER           :: i

  str = fnamestr
  
  DO i=1,n
     
     fstr = str(INDEX(str,'('):INDEX(str,')'))
     IF(fstr(2:2).EQ.'I')THEN
        WRITE(vstr,fstr) INT(vdata(i))
     ELSE
        WRITE(vstr,fstr) vdata(i)
     END IF
  
     str = str(1:INDEX(str,'(')-1) // TRIM(ADJUSTL(vstr)) // str(INDEX(str,')')+1:LEN(str)) 
     
  END DO
  
  GetFileName = TRIM(ADJUSTL(str))

  RETURN
END FUNCTION GetFileName


!--------------------------------------------------------------------
! CheckOutput:
!
! IErr	error code

SUBROUTINE CheckOutput( IWidth, Energy, HubDis, RimDis, ISeed, IErr )

  USE MyNumbers 
  USE IChannels
!  USE DPara
!  USE IPara

  
  INTEGER(KIND=IKIND) IWidth, IErr, ERR, Iseed
  REAL(KIND=RKIND) HubDis, RimDis, Energy
  
  CHARACTER*100 FileName
  
  !PRINT*,"DBG: CheckOutput()"
  
  IErr= 0
  
  !   WRITE out the input parameter
  IF( Energy .GE. 0.0D0) THEN
     WRITE(FileName, '(A4,A2,I4.4,A5,I7.7,A7,I4.4,A7,I4.4,A1,I4,A4)') &
          "Eval-M",IWidth,  &
          "-TarE", NINT(10000.0D0*ABS(Energy)), &
          "-HubDis", NINT(100.0D0*ABS(HubDis)), &
          "-RimDis", NINT(100.0D0*ABS(RimDis)), "-",& 
          !Number, "-", &
          ISeed, ".raw"
  ELSE
     WRITE(FileName, '(A4,A2,I4.4,A2,A5,I7.7,A7,I4.4,A7,I4.4,A1,I4,A4)') &
          "Eval","-M",IWidth, "-m", &
          "-TarE", NINT(10000.0D0*ABS(Energy)), &
          "-HubDis",NINT(100.0D0*ABS(HubDis)), "-", &
          "-RimDis", NINT(100.0D0*ABS(RimDis)), "-",&
          ISeed, ".raw"
  ENDIF
  
  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'NEW', FILE=FileName)
  
  IErr= 0
  
20 CLOSE(UNIT= IChOut, ERR= 100)
  
  RETURN
  
10 WRITE(*,15) FileName
15 FORMAT(" CheckOutput(): ", A28,&
        " exists -- skipped!")
  
  IErr= 2
  GOTO 20
  
  !  ERR in CLOSE detected
100 &
  PRINT*,"CheckOutput(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE CheckOutput


!--------------------------------------------------------------------
! WriteOutputEVal:
!
! IErr	error code

SUBROUTINE WriteOutputEVal(Dim, Nx, NEVals, EIGS, IWidth, Energy, HubDis, RimDis, ISeed, str, IErr)


  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  INTEGER(KIND=IKIND) Dim, Nx
  INTEGER(KIND=IKIND) IWidth, IErr, ERR, NEVals, Iseed, i
  REAL(KIND=RKIND) HubDis, RimDis, Energy
  REAL(KIND=RKIND) EIGS(NEVals)
  
  CHARACTER*100 FileName, str
  
  IErr= 0
  
  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A1,I4.4,A4)') &
          "Eval-", "L", Dim, Nx, &
          "-M",IWidth, "-p", &
          "-TarE", NINT(1000000.0D0*ABS(Energy)), &
          "-HubDis", NINT(10000.0D0*ABS(HubDis)), &
          "-RimDis", NINT(10000.0D0*ABS(RimDis)), "-",& 
          ISeed, ".raw"
  ELSE
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A1,I4.4,A4)') &
          "Eval-","L",Dim, Nx, &
          "-M",IWidth, "-m",&
          "-TarE", NINT(1000000.0D0*ABS(Energy)), &
          "-HubDis",NINT(10000.0D0*ABS(HubDis)), &
          "-RimDis", NINT(10000.0D0*ABS(RimDis)), "-",&
          ISeed, ".raw"
  ENDIF
  
  !        IF(IWriteFlag .GE. 2) THEN
  PRINT*, "eVAL file    ", FileName
  !        ENDIF
  
!!$  OPEN(UNIT= IChEVal, ERR= 10, STATUS='UNKNOWN', FILE=Trim(str)//"/"//FileName)  
  
  OPEN(UNIT= IChEVal, ERR= 10, STATUS='UNKNOWN', FILE= Trim(AdjustL(str))//"/"//FileName)
  
  IF(NEVals .GT. 0)THEN
     DO i=1,NEVals
        WRITE(IChEVal, FMT=15, ERR=20) EIGS(i)
     ENDDO
  END IF
  
  CLOSE(UNIT=IChEVal, ERR= 30)
  
  RETURN
  
15 FORMAT(f30.20)

  !	error in OPEN detected
10 PRINT*, "WriteOutputEVals(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  !	error in WRITE detected
20 PRINT*,"WriteOutputEVals(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
30 PRINT*,"OutputEVals(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE WriteOutputEVal


!--------------------------------------------------------------------
! WriteOutputEVec:
!
! IErr	error code

SUBROUTINE WriteOutputEVec( Dim, Nx, Inum, NEVals, Lsize, VECS, VECS_size, &
                   IWidth, Energy, HubDis, RimDis, Seed, str, IErr)

  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  INTEGER(KIND=IKIND) Dim, Nx
  INTEGER(KIND=IKIND) Inum, Seed, IWidth, IErr, ERR, Lsize, VECS_size, NEVals, i
  REAL(KIND=RKIND) HubDis, RimDis, Energy

  REAL(KIND=RKIND) VECS(VECS_size)

  CHARACTER*100 FileName, str

  IErr= 0

  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A2,I4.4,A1,I4.4,A4)') &
          
          "Evec-","L",Dim, Nx, &
          "-M",IWidth, "-p",&
          "-TarE", NINT(1000000.0D0*ABS(Energy)), &
          "-HubDis", NINT(10000.0D0*ABS(HubDis)), &
          "-RimDis", NINT(10000.0D0*ABS(RimDis)), &
          "-N", Inum, "-", &
          Seed, ".raw"
  ELSE
     WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A2,A5,I9.9,A7,I7.7,A7,I7.7,A2,I4.4,A1,I4.4,A4)') &
          "Evec-","L",Dim, Nx, &
          "-M",IWidth, "-m", &
          "-TarE", NINT(1000000.0D0*ABS(Energy)), &
          "-HubDis", NINT(10000.0D0*ABS(HubDis)), &
          "-RimDis", NINT(10000.0D0*ABS(RimDis)), &
          "-N", Inum, "-", &
          Seed, ".raw"
  ENDIF

  print*,'evector file ',FileName

  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=FileName)

  DO i= 1+( Lsize*( Inum -1) ), Lsize*Inum
     WRITE(UNIT=IChEVec, FMT=45, ERR=50) VECS(i)
  ENDDO

  CLOSE(UNIT= IChEVec, ERR= 60)
  
  RETURN

45 FORMAT(f30.20)

  !	error in OPEN detected
40 PRINT*,"WriteOutputEVec(): ERR in OPEN()"
  IErr= 1
  RETURN

  !	error in WRITE detected
50 PRINT*,"WriteOutputEVec(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
60 PRINT*,"OutputEVec(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteOutputEVec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Create Folder !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine GetDirec(Dim, Nx, Width, HubDis, RimDis, Seed, str)

  Integer*4 Dim, Nx, Width, Seed
  Real*8 HubDis, RimDis
  Character(len=100) str
  Character(len=10) fid1, fid2, fid3, fid4, fid5, fid6
  Logical*4 ierr1

  write(fid1,'(I1)') Dim; fid1=Trim(AdjustL(fid1))
  write(fid2,'(I1)') Nx; fid2=Trim(AdjustL(fid2))
  write(fid3,'(I3)') Width; fid3=Trim(AdjustL(fid3))
  write(fid4,'(f8.4)') HubDis; fid4=Trim(AdjustL(fid4))
  write(fid5,'(f8.4)') RimDis; fid5=Trim(AdjustL(fid5))
  write(fid6,'(I4)') Seed
  

  str='L'//Trim(fid1)//Trim(fid2)//'_M'//Trim(fid3)//'_HubDis'//Trim(fid4) &
       //'_RimDis'//Trim(fid5)//"_"//Trim(fid6)//'_.DATA'

!  Write(str,'(A1,I1,I1,A2,I3.1,A7,f6.1,A7,f6.1,A6)') &
!       "L", Dim, Nx, "_M", Width, "_HubDis", HubDis, &
!       "_RimDis", RimDis, "_.DATA"

!!$  PRINT*,str
  Print*, str

  Inquire(file=Trim(AdjustL(str)), Exist=ierr1)
  If(ierr1)Then
     Print*,"The directory have existed and don't need to be created"
     Write(*,'(/)')
  Else
     Print*,"The directory don't exist and create it"
     Write(*,'(/)')
     call System("mkdir "//Trim(AdjustL(str)) )
  End if
  
  RETURN
  
END Subroutine GetDirec

