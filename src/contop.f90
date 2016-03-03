! SUBROUTINE CONTOP(T,TKEV,TK,HKT,TLOG,XNA,XNE,WLGRID, &
!      OPACITY,SCATTER,FRACT, INDXSP, &
!      IXH1,IXH2,IXHMIN,IXHE1,IXHE2,IXHE3,IXC1,IXAL1,IXSI1, &
!      IXSI2,IXCA1,IXCA2,IXMG1,IXMG2,IXFE1,IXN1,IXO1, nWLGRID, NLIST, NTOTALLIST)
 SUBROUTINE CONTOP(T,TKEV,TK,HKT,TLOG,XNA,XNE,WLGRID, &
      OPACITY,SCATTER, H1,H2,HMIN,HE1,HE2,HE3,C1,AL1,SI1, &
      SI2,CA1,CA2,MG1,MG2,FE1,N1,O1, nWLGRID, NLIST, NTOTALLIST)
  
!
!  This subroutine computes continuous opacities using B. Kurucz ATLAS-9
!  opacity package.
!
!  Author: N.Piskunov based on Kurucz ATLAS-9 subroutine.
!
!  UPDATE: 25-May-1999, final version for SYNTHMAG
!          24-Jun-2004, re-written for 3D and ported to Fotran90
!          10-Oct-2005, Splitted the storage of total opacity and scattering
!
    IMPLICIT NONE
    INTEGER :: nWLGRID, NLIST, NTOTALLIST
    REAL        :: H1,H2,HMIN,HE1,HE2,HE3, C1, AL1, SI1
    REAL        :: SI2,CA1,CA2,MG1,MG2,FE1,N1, O1
!    INTEGER        :: INDXSP(NLIST)
    REAL(KIND=8) :: WLGRID(nWLGRID),FREQLG
    REAL(KIND=4) :: FREQ,FREQ15
    REAL(KIND=4) :: T,TKEV,TK,HKT,EHVKT,STIM,TLOG,XNA,XNE !,FRACT(NTOTALLIST)
    REAL(KIND=8) :: OPACITY(nWLGRID),SCATTER(nWLGRID)
    REAL(KIND=4) :: AHYD,AH2P,AHMIN,AHE1,AHE2,AHEMIN,ACOOL,ALUKE,AHOT,A,B
    REAL(KIND=4) :: SIGH,SIGHE,SIGEL,SIGH2
    INTEGER :: iWL
!
! Compute the continuous opacity
    !

    DO iWL=1,nWLGRID
      FREQ=2.997925E18/WLGRID(iWL)
      FREQLG=LOG(FREQ)
      FREQ15=FREQ*1.E-15
!      HNUKT=1.43868D8/WLGRID(iWL)/T
      EHVKT=EXP(-FREQ*HKT)
      STIM=1.-EHVKT
      AH2P  =0.0
      AHE1  =0.0
      AHE2  =0.0
      AHEMIN=0.0
      ACOOL =0.0
      ALUKE =0.0
      AHOT  =0.0
      AHYD  =0.0
      AHMIN =0.0
      SIGH  =0.0
      SIGHE =0.0
      SIGEL =0.0
      SIGH2 =0.0
      CALL HOP(AHYD,XNE,H1, &
               H2,FREQ,FREQLG,T, &
               TLOG,TKEV,STIM,EHVKT)
      CALL H2PLOP(AH2P,H1, &
               H2,FREQ,FREQLG,FREQ15, &
               TKEV,STIM)
      CALL HMINOP(AHMIN,H1, &
               HMIN,FREQ,T,TKEV,XNE, &
               EHVKT)
      CALL HRAYOP(SIGH,H1,FREQ)
      CALL HE1OP(AHE1,HE1, &
               HE2,XNE,FREQ,FREQLG, &
               T,TKEV,TLOG,EHVKT,STIM)
      CALL HE2OP(AHE2,HE2, &
               HE3,XNE,FREQ,FREQLG, &
               T,TKEV,EHVKT,STIM)
      CALL HEMIOP(AHEMIN,HE1, &
               FREQ,T,XNE)
      CALL HERAOP(SIGHE,HE1,FREQ)
      IF(T.LT.12000.) THEN
        CALL COOLOP(ACOOL,C1 , &
               Mg1, &
               Al1, &
               Si1, &
               Fe1, &
               STIM,FREQ,FREQLG,T,TLOG,TKEV,HKT)
      ENDIF
      IF(T.LT.30000.) THEN
        CALL LUKEOP(ALUKE, N1, &
               O1, &
               Mg2, &
               Si2, &
               Ca2, &
               STIM,FREQ,FREQLG,T,TLOG,TKEV)
      ENDIF
      CALL HOTOP(AHOT)
      CALL ELECOP(SIGEL,XNE)
      CALL H2RAOP(SIGH2,H1,FREQ,T,TKEV,TLOG)
!
      A=AHYD+AHMIN+AH2P+AHE1+AHE2+AHEMIN+ACOOL+ALUKE+AHOT
      B=SIGH+SIGHE+SIGEL+SIGH2
      OPACITY(iWL)=A+B
      SCATTER(iWL)=B
!      write(*,'(F10.1,13E11.4)') T,A,B,AH2P,AHE1, &
!        AHEMIN,ACOOL,ALUKE,AHOT,AHYD,AHMIN
!        SIGH,SIGHE,SIGEL,SIGH2
!      write(*,'(F10.1,6E11.4)') T,B,AHYD,AHMIN, &
!         AH2P,AHE1,AHE2
    END DO
!
    RETURN
  END SUBROUTINE CONTOP
  
