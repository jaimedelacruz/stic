      SUBROUTINE EOSFCN(NEQ,P,RHS,A,IFLAG,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      INTEGER NEQ,IFLAG,NCH(SPLSIZ-1),NLIST,IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL ABUND(*)
      REAL*8 P(NEQ),RHS(NEQ),A(ELEDIM+1,NEQ),PG,KT(*),IT(*)
      INTEGER I,II,J,JJ,K,KK,KKK,JATOM,NQ,ISPEC,NELT
      REAL*8 PE,CRATIO,PF,PENQ,PN,DUMMY,DPF(4),AT
      REAL*8 AAA(ELEDIM+1),BBB
C================================================================
C Method: We are solving a system of non-linear equations
C         (the summation is always carried over all species)
C
C  Particle conservation:
C
C    F1   = P_total - Sum(P_species) - P_elec = 0
C
C  Abundance equations (for each atom "a"):
C
C    F2   = Sum[P_species*(Z_a*N_species - N_a_species) = 0
C    F3   = ...
C
C  where Z_a         is the abundance of atom "a"
C        N_species   is the total number of atoms in a given species
C        N_a_species is the number of atoms "a" in a given species
C
C  Charge conservation:
C
C    Fneq = P_elec - Sum(P_species * Q_species) = 0
C
C  where Q_species is the charge of a given species.
C
C  The unknowns are the ficticious pressures for all atoms:
C  P_a = N_a*kT and P_elec
C
C  Newton-Raphson scheme is used for the solution:
C
C       dF_j
C  Sum( ---- * delta P_b ) = -F_i(P_a)
C       dP_b
C
C  The sytem of linear equations is solved with LU decomposition.
C
C  It is not unusual that the rank of the Jacobian is huge and the
C  system of linear equations is ill-defined. Instead of using SVD
C  we verify that the new P_a actually reduce the absolute magnitude
C  of Fi. If they don't we scale down the corrections until Fi are
C  as close to zero as possible.
C================================================================
C== RHS vector update                                          ==
C================================================================
      IF(IFLAG.EQ.1) THEN
        JATOM=NEQ-1
        PE=P(NEQ)
        DO K=2,JATOM
          RHS(K)=0.D0
        ENDDO
        RHS(  1)=-PG
        RHS(NEQ)=-PE
        BBB=0.D0

        DO ISPEC=1,NLIST-1
          NQ=NCH(ISPEC)
C
C  Compute PN - partial pressure of species ISPEC and it's partial
C  derivatives in respect to all ficticious atomic parial pressures 
C
          PF=1.0D0
C
C  Loop through all constituent atoms
C
          NELT=NEL(ISPEC)
          DO I=1,NELT
            J=INDZAT(ZAT(I,ISPEC))
            CRATIO=P(J)**NAT(I,ISPEC)
            PF=PF*CRATIO
          ENDDO
C
C  Be careful with zero electron pressure
C
          IF(PE.GT.0.0D0.AND.NQ.NE.0) THEN
            PENQ=PE**NQ
            CRATIO=IT(ISPEC)/PENQ/KT(ISPEC)
          ELSE
            CRATIO=IT(ISPEC)/KT(ISPEC)
          ENDIF
          PN=CRATIO*PF
C
C  Fill the RHS vector
C
          RHS(1)=RHS(1)+(NQ+1)*PN
          BBB=BBB+NTOT(ISPEC)*PN
          DO II=1,NELT
            KKK=INDZAT(ZAT(II,ISPEC))
            IF(KKK.NE.1) RHS(KKK)=RHS(KKK)-NAT(II,ISPEC)*PN
          ENDDO
          RHS(NEQ)=RHS(NEQ)+NQ*PN
        ENDDO
        DO J=2,JATOM
          RHS(J)=RHS(J)+ABUND(IATOM(J))*BBB
c          RHS(J)=RHS(J)*(1.D0+1.D20*MIN(P(J),0.D0)**2)
        ENDDO
        RETURN
C================================================================
C== Jacobian matrix update                                     ==
C================================================================
      ELSE IF(IFLAG.EQ.2) THEN
        JATOM=NEQ-1
        PE=P(NEQ)
        DO JJ=1,NEQ
          DO J=1,NEQ
            A(J,JJ)=0.0D0
          ENDDO
          AAA(JJ)=0.D0
        ENDDO
        A(NEQ,NEQ)=-1.0D0
        BBB=0.D0
C
C Loop through every species, except the last (ISPEC=NLIST) which is "e-".
C Fill the matrix of linearized equations.
C
        DO ISPEC=1,NLIST-1
          NQ=NCH(ISPEC)
C
C  Compute PN - partial pressure of species ISPEC and it's partial
C  derivatives DPF in respect to all ficticious atomic parial pressures 
C
          PF=1.0D0
          NELT=NEL(ISPEC)
          DO I=1,NELT
            DPF(I)=1.0D0
          ENDDO
C
C  Loop through all constituent atoms
C
          DO I=1,NELT
            J=INDZAT(ZAT(I,ISPEC))
            CRATIO=P(J)**NAT(I,ISPEC)
C
C  Compute the product
C
            PF=PF*CRATIO
            DUMMY=DPF(I)
C
C  Update the factors for derivative over P(J)
C
            DO K=1,NELT
              DPF(K)=DPF(K)*CRATIO
            ENDDO
C
C  Correct the only factor dependent on P(J)
C
            IF(NAT(I,ISPEC).GT.1) THEN
              DPF(I)=DUMMY*P(J)**(NAT(I,ISPEC)-1)*NAT(I,ISPEC)
            ELSE
              DPF(I)=DUMMY
            ENDIF
          ENDDO
C
C  Be careful with zero electron pressure
C
          IF(PE.GT.0.0D0.AND.NQ.NE.0) THEN
            PENQ=PE**NQ
            CRATIO=IT(ISPEC)/PENQ/KT(ISPEC)
          ELSE
            CRATIO=IT(ISPEC)/KT(ISPEC)
          ENDIF
          PN=CRATIO*PF
C
C  Fill in the Jacobian matrix
C
          DO I=1,NELT
            KK=INDZAT(ZAT(I,ISPEC))
            AT=CRATIO*DPF(I)
            A(1,KK)=A(1,KK)+(NQ+1)*AT
            AAA(KK)=AAA(KK)+NTOT(ISPEC)*AT
            DO II=1,NELT
              KKK=INDZAT(ZAT(II,ISPEC))
c              IF(KKK.NE.1) A(KKK,KK)=A(KKK,KK)+
c     *         (NTOT(ISPEC)*DBLE(ABUND(IATOM(KKK)))-NAT(II,ISPEC))*AT
              IF(KKK.NE.1) A(KKK,KK)=A(KKK,KK)-NAT(II,ISPEC)*AT
            ENDDO
            A(NEQ,KK)=A(NEQ,KK)+NQ*AT
          ENDDO
          AT=0.0D0
          IF(PE.GT.0.0D0.AND.NQ.NE.0) THEN
            AT=NQ*PN/PE
            A(1,NEQ)=A(1,NEQ)-(NQ+1)*AT
            BBB=BBB-NTOT(ISPEC)*AT
            DO II=1,NELT
              KKK=INDZAT(ZAT(II,ISPEC))
c              IF(KKK.NE.1) A(KKK,NEQ)=A(KKK,NEQ)+
c     *          (NAT(II,ISPEC)-NTOT(ISPEC)*DBLE(ABUND(IATOM(KKK))))*AT
              IF(KKK.NE.1) A(KKK,NEQ)=A(KKK,NEQ)+NAT(II,ISPEC)*AT
            ENDDO
            A(NEQ,NEQ)=A(NEQ,NEQ)-NQ*AT
          END IF
        ENDDO
        DO K=2,JATOM
          DO KK=1,JATOM
            A(K,KK)=A(K,KK)+DBLE(ABUND(IATOM(K)))*AAA(KK)
          ENDDO
          A(K,NEQ)=A(K,NEQ)+DBLE(ABUND(IATOM(K)))*BBB
        ENDDO
C
        RETURN
      ENDIF
C
      END

      SUBROUTINE lnEOSFCN(NEQ,P,RHS,A,IFLAG,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      INTEGER NEQ,IFLAG,NCH(SPLSIZ-1),NLIST,IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL ABUND(*)
      DOUBLE PRECISION P(NEQ),RHS(NEQ),A(ELEDIM+1,NEQ),PG,KT(*),IT(*)
      INTEGER I,II,J,JJ,K,KK,JATOM,NQ,ISPEC,NELT
      DOUBLE PRECISION PE,CRATIO,PF,PENQ,PN,AT,AAA(ELEDIM+1)
      DOUBLE PRECISION BBB,PENORM
C================================================================
C Method: We are solving a system of non-linear equations
C         (the summation is always carried over all species)
C
C  Particle conservation:
C
C    F1   = P_total - Sum(P_species) - P_elec = 0
C
C  Abundance equations (for each atom "a"):
C
C    F2   = Sum[P_species*(Z_a*N_species - N_a_species) = 0
C    F3   = ...
C
C  where Z_a         is the abundance of atom "a"
C        N_species   is the total number of atoms in a given species
C        N_a_species is the number of atoms "a" in a given species
C
C  Charge conservation:
C
C    Fneq = P_elec - Sum(P_species * Q_species) = 0
C
C  where Q_species is the charge of a given species.
C
C  The unknowns are the ficticious pressures for all atoms:
C  P_a = N_a*kT and P_elec
C
C  Newton-Raphson scheme is used for the solution:
C
C       dF_j
C  Sum( ---- * delta P_b ) = -F_i(P_a)
C       dP_b
C
C  The sytem of linear equations is solved with LU decomposition.
C
C  It is not unusual that the rank of the Jacobian is huge and the
C  system of linear equations is ill-defined. Instead of using SVD
C  we verify that the new P_a actually reduce the absolute magnitude
C  of Fi. If they don't we scale down the corrections until Fi are
C  as close to zero as possible.
C================================================================
C== RHS vector update                                          ==
C================================================================
      IF(IFLAG.EQ.1) THEN
        JATOM=NEQ-1
        PE=P(NEQ)
        DO K=1,NEQ
          RHS(K)=0.D0
        END DO
        BBB=0.D0

        DO ISPEC=1,NLIST-1
c        DO ISPEC=120,128
          NQ=NCH(ISPEC)
C
C  Compute PN - partial pressure of species ISPEC and it's partial
C  derivatives in respect to all ficticious atomic parial pressures
C
          PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
C
C  Loop through all constituent atoms
C
          NELT=NEL(ISPEC)
          DO I=1,NELT
            J=INDZAT(ZAT(I,ISPEC))
c            write(*,*) I,J,PF,NAT(I,ISPEC),JATOM
c            if(j.lt.1.or.j.gt.JATOM) stop
            PF=PF+P(J)*NAT(I,ISPEC)
          END DO
c          if(PF.gt.10.) then
c            write(*,*) '1)',ISPEC,PF,PE,NQ,IT(ISPEC),KT(ISPEC)
c            stop
c          endif
C
C  Add log of electron pressure and ionization, dissociation constants
C
          IF(PF.GT.-100.d0) THEN
            PN=EXP(PF)
          ELSE
            PN=0.d0
          ENDIF
C
C  Fill the RHS vector
C
          RHS(1)=RHS(1)+(NQ+1)*PN
c          write(*,*) ISPEC,RHS(1),PG
          BBB=BBB+NTOT(ISPEC)*PN
          DO I=1,NELT
            K=INDZAT(ZAT(I,ISPEC))
            IF(K.GT.1) RHS(K)=RHS(K)-NAT(I,ISPEC)*PN
          END DO
          RHS(NEQ)=RHS(NEQ)+NQ*PN
        END DO
c        write(*,*) RHS(1),PG
        RHS(1)=RHS(1)-PG
        DO J=2,JATOM
          RHS(J)=RHS(J)+ABUND(IATOM(J))*BBB
c          if(abs(RHS(J)).gt.1.d20) then
c            write(*,*) j,rhs(1),rhs(j),RHS(NEQ)-EXP(PE)
c            stop
c          endif
        ENDDO
        RHS(NEQ)=RHS(NEQ)-EXP(PE)
        RETURN
C================================================================
C== Jacobian matrix update                                     ==
C================================================================
      ELSE IF(IFLAG.EQ.2) THEN
        JATOM=NEQ-1
        PE=P(NEQ)
        DO JJ=1,NEQ
          DO J=1,NEQ
            A(J,JJ)=0.0D0
          END DO
        END DO
C
C Loop through every species, except the last (ISPEC=NLIST) which is "e-".
C Fill the matrix of linearized equations.
C
        DO ISPEC=1,NLIST-1
c        DO ISPEC=317,317
          NQ=NCH(ISPEC)
          NELT=NEL(ISPEC)
          PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
          DO I=1,NELT
            J=INDZAT(ZAT(I,ISPEC))
            PF=PF+P(J)*NAT(I,ISPEC)
c            write(*,'(I5,2I2,2E12.4,I2,2E12.4)')
c     *            ISPEC,I,J,PF,P(J),NAT(I,ISPEC),IT(ISPEC),KT(ISPEC)
          END DO
c          write(*,'(I5,4E12.4)') ISPEC,PF,PE*NQ,LOG(IT(ISPEC))
c     *                                ,LOG(KT(ISPEC))
          IF(PF.GT.-250.d0) THEN
            PN=EXP(PF)
          ELSE
            PN=0.d0
          ENDIF
          PENORM=EXP(PE)
C
C Particle conservation equation (Eq. 1)
C   Derivatives over log atomic partial pressures
C
          DO I=1,NELT
            K=INDZAT(ZAT(I,ISPEC))
            A(1,K)=A(1,K)+PN*(1+NQ)*NAT(I,ISPEC)
          END DO
C
C Particle conservation equation (Eq. 1)
C   Derivative over log electron partial pressures
C
          A(1,NEQ)=A(1,NEQ)-PN*(1+NQ)*NQ
C
C Abundance equations (Eq. 2...NEQ-1)
C   Derivatives over log atomic partial pressures
C
          DO K=2,JATOM
            DO II=1,NELT
              KK=INDZAT(ZAT(II,ISPEC))
              A(K,KK)=A(K,KK)+PN*NTOT(ISPEC)*DBLE(ABUND(IATOM(K)))*
     *                        NAT(II,ISPEC)
c       if(K.eq.26) write(*,*) ISPEC,A(K,KK),
c     *     PN*NTOT(ISPEC)*DBLE(ABUND(IATOM(K))),NAT(II,ISPEC)
            END DO
          END DO
C
          DO I=1,NELT
            K=INDZAT(ZAT(I,ISPEC))
            IF(K.GT.1) THEN
              DO II=1,NELT
                KK=INDZAT(ZAT(II,ISPEC))
                A(K,KK)=A(K,KK)-PN*NAT(II,ISPEC)*NAT(I,ISPEC)
c       if(K.eq.26) write(*,*) ISPEC,A(K,KK),
c     *     PN*NAT(II,ISPEC)*NAT(I,ISPEC),NAT(I,ISPEC),NAT(II,ISPEC)
              END DO
            END IF
          END DO
C
C Abundance equations (Eq. 2...NEQ-1)
C   Derivative over log electron partial pressures
C
          DO K=2,JATOM
            A(K,NEQ)=A(K,NEQ)-PN*NTOT(ISPEC)*DBLE(ABUND(IATOM(K)))*NQ
          END DO
C
          DO I=1,NELT
            K=INDZAT(ZAT(I,ISPEC))
            IF(K.GT.1) A(K,NEQ)=A(K,NEQ)+PN*NAT(I,ISPEC)*NQ
          END DO
C
C Charge neutrality equation (Eq. NEQ)
C   Derivatives over log atomic partial pressures
C
          DO I=1,NELT
            K=INDZAT(ZAT(I,ISPEC))
            A(NEQ,K)=A(NEQ,K)+PN*NAT(I,ISPEC)*NQ
          END DO
C
C Charge neutrality equation (Eq. NEQ)
C   Derivative over log electron partial pressures
C
          A(NEQ,NEQ)=A(NEQ,NEQ)-PN*NQ*NQ
        END DO
        A(NEQ,NEQ)=A(NEQ,NEQ)-PENORM
C
c        write(*,'(''1)'',41e10.3)')(a(i,38),i=1,40),RHS(38)
c        DO I=1,NEQ
c          write(*,'(42(f5.0))')
c     *      (LOG10(MAX(abs(A(I,J)),1d-99)),J=1,NEQ),log10(abs(RHS(I)))
c        enddo
        RETURN
      ENDIF
C
      END

