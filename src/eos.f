C=========================================================================
C EQCOUNT: Counts the list of species for solving the equation of state by
C          merging the default list and species present in the line list.
C
C We assume that only neutral molecules can appear in the line list.
C For atoms, all the ions present in the table of partition functions
C are added to the list. Atomic names are case sensitive, that is the first
C character must be uppercase and for 2 character names the second character
C must be lower case.
C
C Inputs:
C   ELEMEN - the names of chemical elements in the periodic table
C   SPNAME - the names of the species present in the line lists + continuous
C            absorbers
C   ION    - ionization stage (1 -neutral, 2 - first ion etc.)
C   NLINES - the length of the line list, also dimensions of arrays SPNAME,
C            ION, SPINDX
C   NLIST  - if >0 on input, indicates that the default list of species have
C            been loaded, otherwise EQLIST loads the default list to SPLIST.
C   ELESIZ - Size of all arrays related to atomic list.
C
C   Return code  0: OK
C                1: illegal species name
C               >1: SPLSIZ is too small
C
      integer function eqcount(elemen,spname,ion,nlines,nlist,
     *                         ELESIZ)
c      integer function eqcount(elemen,spname,ion,nlines,nlist,
c     *                         environment,ELESIZ)
      INCLUDE 'SIZES.EOS'

      integer nlines,nlist,ELESIZ
      character*(3) elemen(ELESIZ)
      character*2 tmp
      character*(SPCHAR) spname(nlines)
      character*(SPCHAR) tmplist(SPLSIZ),chname
      integer ion(nlines),ionmax,ionmaxx
      real a(IONSIZ)
c      character*(*) environment
      double precision b(IONSIZ)
      INCLUDE 'DEFAULT.EOS.current'
c      INCLUDE 'DEFAULT.EOS'
C
      eqcount=0
      ionmax=0
      ncount=NDEF
c      if(environment.eq.'COLD'.or.environment.eq.'cold') then
c        do ispec=1,NDEF_cold
c          tmplist(ispec)=default_cold(ispec)
c        end do
c        ionmax=2
c        ncount=NDEF_cold
c      else if(environment.eq.'WARM'.or.environment.eq.'warm') then
c        do ispec=1,NDEF_warm
c          tmplist(ispec)=default_warm(ispec)
c        end do
c        ncount=NDEF_warm
c      else if(environment.eq.'HOT'.or.environment.eq.'hot') then
c        do ispec=1,NDEF_hot
c          tmplist(ispec)=default_hot(ispec)
c        end do
c        ncount=NDEF_hot
c      else
c        do ispec=1,NDEF_cool
c          tmplist(ispec)=default_cool(ispec)
c        end do
c        ncount=NDEF_cool
c      end if
C
C Associate each species in SPNAME with an entry in SPLIST. If SPNAME
C  contains a new species not in SPLIST, then add that new species at
C  the end of SPLIST.
C
      if(nlines.gt.0) then
        do 6 ilin=1,nlines
          call mbuild(spname(ilin),ion(ilin)-1,chname)
c          write(*,*) ncount,ilin,ionmax,spname(ilin),chname
          do ispec=1,ncount
            if(tmplist(ispec).eq.chname) goto 6
          end do
c          write(*,*) ncount,ilin,chname,ionmax,spname(ilin),ion(ilin)
c       stop
C
C Look for atomic species. Negative ions (e.g. H-) are treated as molecules
C
          if((spname(ilin)(2:2).EQ.' '.OR.
     *       (spname(ilin)(3:3).EQ.' '.AND.
     *        spname(ilin)(2:2).GE.'a'.AND.
     *        spname(ilin)(2:2).LE.'z')).AND.
     *        ion(ilin).GT.0) then
            iel=0
            tmp=spname(ilin)(1:2)
            do i=1,ELESIZ
              if(tmp.eq.elemen(i)(1:2)) then
                iel=i
                goto 4
              endif
            end do
            if(iel.lt.1) then
              eqcount=1
c              return
              write(*,*) 'eqcount: Wrong species: ',spname(ilin)
              stop
            end if
   4        call XSAHA(iel,1.,1.,1.,ionmaxx,a,b,5)
            if(ionmax.gt.0) ionmaxx=ionmax
            if(ionmaxx.lt.ion(ilin)) then
              write(*,*) ilin,ion(ilin),nlines
              write(*,*) 'XSAHA has no partition function for '//chname
              stop
            endif
            tmplist(ncount+1)=elemen(iel)(1:2)
            if(ionmaxx.gt.1) then
              do i=2,ionmaxx
                ncount=ncount+1
                i1=index(tmplist(ncount),' ')
                tmplist(ncount+1)=tmplist(ncount)(1:i1-1)//'+'
              end do
            end if
            ncount=ncount+1
          else
C
C Molecules are counted here
C
            tmplist(ncount+1)=chname
            ncount=ncount+1
          end if
   6    continue
      endif
C
C All lines have been processed, add free electrons and return
C
      nlist=ncount+1
      eqcount=0
C
      return
      end

C=========================================================================
C EQLIST: Creates the list of species for solving the equation of state by
C         merging the default list and species present in the line list.
C
C We assume that only neutral molecules can appear in the line list.
C For atoms, all the ions present in the table of partition functions
C are added to the list. Atomic names are case sensitive, that is the first
C character must be uppercase and for 2 character names the second character
C must be lower case.
C
C Inputs:
C   ELEMEN - the names of chemical elements in the periodic table
C   SPNAME - the names of the species present in the line lists + continuous
C            absorbers
C   ION    - ionization stage (1 -neutral, 2 - first ion etc.)
C   NLINES - the length of the line list, also dimensions of arrays SPNAME,
C            ION, SPINDX
C   NLIST  - if >0 on input, indicates that the default list of species have
C            been loaded, otherwise EQLIST loads the default list to SPLIST.
C   SPLDIM - maximum length of the compiled lists of species SPLIST (must
C            be smaller than SPLSIZ).
C   ELESIZ - Size of all arrays related to atomic list.
C
C Outputs:
C   SPINDX - index array of size NLINES which upon return holds pointers to
C            the complete list of species SPLIST: line L is produced by
C            species SPLIST(SPINDEX(L))
C   SPLIST - upon return contains the compiled list of all species (default
C            list + species in the line list + continuous absorbers)
C   NLIST  - the size of the compiled list of species SPLIST
C
C   Return code  0: OK
C                1: illegal species name
C                2: SPLDIM is too small)
C                3: Missing ionization stage
C                4: e- is not the last item in the list
C                5: Unreasonable abundances
C
C  2006.12.27 - converted eqlist to a function for compatibility with the SME
C
C
      integer*4 function eqlist(abund,elemen,spname,ion,spindx,splist,
     &                  nlines,nlist,SPLDIM,ELESIZ)
c      integer*4 function eqlist(abund,elemen,spname,ion,spindx,splist,
c     &                  nlines,nlist,environment,SPLDIM,ELESIZ)
      INCLUDE 'SIZES.EOS'

      integer nlines,nlist,SPLDIM,ELESIZ
      character*(SPCHAR) spname(nlines),splist(SPLDIM)
      character*(3) elemen(ELESIZ)
c      character*(*) environment
      character*2 tmp
      integer ion(nlines),spindx(nlines),ionmax,ionmaxx
      dimension abund(ELESIZ)
      real a(IONSIZ)
      double precision b(IONSIZ)
C
C SPLIST should contain all the major contributors to the electron pressure,
C and all the molecules which significantly affect atomic partial pressures.
C For each call to EQSTAT, the base set of species at the beginning of SPLIST
C are supplemented by any new species that appear in SPNAME. It is common
C for some of the species in the base set (at the beginning of SPNAME) to be
C duplicated in SPNAME. This allows one to get ZETA for these species and is
C not a problem.
C
      integer splmax
      character*(SPCHAR) chname
      INCLUDE 'DEFAULT.EOS.current'
c      INCLUDE 'DEFAULT.EOS'
C
C Determine maximum allowed number of species, based on sizes of arrays
C  defined locally (using SPLSIZ) and passed by argument (using spldim).
C
      splmax=min(SPLSIZ,SPLDIM)
C
C Load base set of species (SPLIST) with default set of species (DEFAULT),
C  if passed value of NLIST is 0. Be sure to include "e-" at the end of
C  SPLIST.
C
      idef=0
      ionmax=0
      if(nlist.eq.0) then
C
C  Copy the default list and check if we have enough space first
C
        do jdef=1,NDEF
          splist(jdef)=default(jdef)
        end do
        nlist=NDEF
cC
cC  Copy the default list and check if we have enough space first
cC
c        if(environment.eq.'COLD'.or.environment.eq.'cold') then
c          do jdef=1,NDEF_cold
c            splist(jdef)=default_cold(jdef)
c          end do
c          ionmax=2
c          nlist=NDEF_cold+idef
c        else if(environment.eq.'WARM'.or.environment.eq.'warm') then
c          do jdef=1,NDEF_warm
c            splist(jdef)=default_warm(jdef)
c          end do
c          nlist=NDEF_warm+idef
c        else if(environment.eq.'HOT'.or.environment.eq.'hot') then
c          do jdef=1,NDEF_hot
c            splist(jdef)=default_hot(jdef)
c          end do
c          nlist=NDEF_hot+idef
c        else
c          do jdef=1,NDEF_cool
c            splist(jdef)=default_cool(jdef)
c          end do
c          nlist=NDEF_cool
c        end if
        idef=nlist
        if(nlist.ge.splmax) goto 900
C
C  nlines set to -1 indicates that we need to get partial pressures for all atoms
C  This mode is meant for use within VALD
C
        if(nlines.eq.-1) then
c
c Add all atoms first (the call to XSAHA is dummy,
C just to get the number of ions available in the table)
c
          do iel=1,ELESIZ
            call XSAHA(iel,1.,1.,1.,ionmaxx,a,b,5)
            if(ionmax.gt.0) ionmaxx=ionmax
            idef=idef+1
            if(idef.gt.splmax) goto 900
            splist(idef)=elemen(iel)(1:2)
            if(ionmaxx.gt.1) then
              do i=2,ionmaxx
                idef=idef+1
                if(idef.gt.splmax) goto 900
                splist(idef)=splist(idef-1)
                isp=index(splist(idef),' ')
                if(isp.le.0) then
                  write(*,*) 'eqlist: Insufficient length of splist ',
     *                       'elements to store ion',elemen(iel)(1:2),i,
     *                       idef,SPCHAR
                  eqlist=2
                  return
                endif
                splist(idef)(isp:isp)='+'
              end do
            end if
          end do
          nlist=idef
        endif
      endif
C
C Check that abundances are sensible.
C
      absum=0.0
      do ielem=1,ELESIZ
        if(abund(ielem).lt.0.0.or.abund(ielem).gt.1.0) then
          write(*,40) ielem,abund(ielem)
  40      format('eqlist: bad abundance for element',i3,':',1pe13.4)
          write(*,*) (abund(ispec),ispec=1,99)
c          stop
          eqlist=5
          return
        endif
        absum=absum+abund(ielem)
      end do
c      do ielem=1,ELESIZ
c        abund(ielem)=abund(ielem)/absum
c      end do
c      if(abs(absum-1.0).gt.1.0e-3) then
c        write(*,70) absum
c  70    format('eqlist: warning! abundances are not normalized:'
c     &           ,1pe13.5)
c      endif

C
C Associate each species in SPNAME with an entry in SPLIST. If SPNAME
C  contains a new species not in SPLIST, then add that new species at
C  the end of SPLIST.
C
      do ispec=nlist+1,splmax
        splist(ispec)='        '
      end do
      inew=nlist+1
      if(nlines.gt.0) then
        do 150 ilin=1,nlines
          call mbuild(spname(ilin),ion(ilin)-1,chname)
          do ispec=1,nlist
            if(splist(ispec).eq.chname) then
              spindx(ilin)=ispec
              goto 150
            endif
          end do
C
C Look for atomic species. Negative ions (e.g. H-) are treated as molecules
C
          if((spname(ilin)(2:2).EQ.' '.OR.
     *       (spname(ilin)(3:3).EQ.' '.AND.
     *        spname(ilin)(2:2).GE.'a'.AND.
     *        spname(ilin)(2:2).LE.'z')).AND.
     *        ion(ilin).GT.0) then
            iel=0
            tmp=spname(ilin)(1:2)
            do i=1,ELESIZ
              if(tmp.eq.elemen(i)(1:2)) iel=i
            end do
            if(iel.lt.1) then
c              write(*,*) 'eqlist: Wrong species: "'//spname(ilin)//'"'
c              stop
              eqlist=1
              return
            end if
            call XSAHA(iel,1.,1.,1.,ionmaxx,a,b,5)
            if(ionmax.gt.0) ionmaxx=ionmax
            if(ionmaxx.lt.ion(ilin)) then
              write(*,*) 'XSAHA has no partition function for '//chname
              stop
            endif
C
C  Make sure that neutral atoms are included as well as all
C  the intermediate ions
C
            do ii=0,ionmaxx-1
              if(inew.gt.splmax) goto 900
              call mbuild(spname(ilin),ii,chname)
              splist(inew)=chname
              if(ii.eq.ion(ilin)-1) spindx(ilin)=inew
              inew=inew+1
            end do
          else
c       write(*,*) 'Molecule: '//chname,inew
            if(inew.gt.splmax) goto 900
            splist(inew)=chname
            spindx(ilin)=inew
            inew=inew+1
          end if
          nlist=inew-1
 150    continue
      endif
C
C Make sure free electrons are the last species in the list.
C
      do ispec=1,nlist-1
        if(splist(ispec).eq.'e-') then
c          write(*,*) 'eqlist: "e-" may only occur at the end of the'
c     &            // ' species list (SPLIST).'
c          stop
          eqlist=4
          return
        endif
      end do
      if(splist(nlist).ne.'e-') then
        nlist=nlist+1
        if(nlist.gt.splmax) goto 900
        splist(nlist)='e-'
      endif
C
C Make sure neutral hydrogen and neutral helium are in SPLIST. These
C  species are needed for H1FRCT and HE1FRCT. Remember the locations
C  of these species in SPLIST for later use. Code is optimized for
C  the case where H and He both occur early in SPLIST list.
C
c      ih1=-1
c      do 200 ispec=1,nlist
c        if(splist(ispec).eq.'H') then
c          ih1=ispec
c          goto 210
c        endif
c 200  continue
c      write(*,*) 'eqlist: "H" must be in species list (SPLIST)'
c      stop
c 210  ihe1=-1
c      do 220 ispec=1,nlist
c        if(splist(ispec).eq.'He') then
c          ihe1=ispec
c          goto 230
c        endif
c 220  continue
c      write(*,*) 'eqlist: "He" must be in species list (SPLIST)'
c      stop
c 230  continue
C
C Sort the list
C
      call sort2(nlist,splist,nlines,spindx,elemen,ELESIZ)
c      do 250 ispec=1,nlist
c 250  write(*,*) ispec,' "',splist(ispec),'"'
c      stop
C
      eqlist=0
      return
C
C Error handlers.
C
 900  continue
c      write(*,905) spldim,splsiz
c 905  format('eqlist: species list (SPLIST) not long enough:',2i5)
c      stop
      eqlist=2
c
      return
      end

c
C=========================================================================
C EQSTAT: Determine thermodynamic quantities required for spectroscopy.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   PTOTAL [real] Total gas pressure (in dyne/cm^2), given by NTOTAL*K*T,
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   PELEC [real] Electron pressure (in dyne/cm^2), given by NELEC*K*T,
C     which is to be used in calculating ionization equilibrium.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPNAME [character*(*) array(NLINES)] Case-sensitive species name of atom
C     or molecule. The first letter of each atom name must be uppercase. The
C     second letter of each atom name, if present, must be lowercase. Each
C     atom name may optionally be followed by a multiplicity number between
C     1 and 4. If no multiplicity number is given for a particular atom, then
C     its multiplicity is assumed to be 1. All atomic and molecular species
C     in SPNAME must be neutral, with the charge state specified separately
C     in the ION input argument.
C   ION [integer array(NLINES)] Charge state for each of the atomic and
C     molecular species specified in SPNAME. ION=-1 for negative ions (e.g.
C     H minus), ION=0 for neutrals, ION=1 for singly ionized species, etc.
C   NLINES [integer] Number of valid entries in SPNAME and ION. From an
C     external perspective, each entry in SPNAME and ION will correspond to
C     a single spectral line, so some specie/charge combinations may appear
C     more than once, while others may not appear at all.
C   SPLDIM [integer] Array sizes for the arguments SPLIST and XFRACT, which
C     contain information for each species. The maximum allowed number of
C     species is SPLMAX=MIN(SPLSIZ,SPLDIM), where SPLSIZ is a parameter
C     defined in the file SIZES.SYN and used to dimension the local arrays
C     XNPF, PFUNC, and POTION. SPLMAX must be large enough to handle the
C     base set of species used when computing the molecular equilibrium and
C     also any additional species that appear only in the line list. Ideally,
C     the calling routine will <1> Include SIZES.SYN, <2> Use SPLSIZ to
C     dimension SPLIST and XFRACT, and <3> Pass SPLSIZ in place of SPLDIM.
C     However, SPLDIM is passed separately to allow for error checking in
C     the cases when this is not done (e.g. when called from IDL).
C   MODE [integer] Determines the content of the content of the the output
C                  array xfract:
C      0    - number densities/partition functions
C      1    - number densities
C      2    - partial pressures
C      3    - number density of free electrons produced by each species
C others    - the same as 0
C   10+     - the same as above but electron density is assumed to be known
C             precisely so the input value is used instead of solving for
C             Pelec
C
C Input/Output:
C   SPLIST [character*(*) array(SPLDIM)] If NLIST is nonzero upon entry,
C     then SPLIST must contain the base set of species that must be included
C     in the molecular equilibrium calculation, regardless of which species
C     are represented by lines in SPNAME. Until the code is cleaned up, the
C     species list in SPLIST must include "e-" after the NLIST element.
C     If NLIST is zero upon entry, then SPLIST is loaded with the base set
C     of species coded into EQSTAT below (in the variable DEFAULT). Again,
C     an "e-" is appended after the base set.
C     Regardless of the whether SPLIST is valid upon entry or needs to be
C     loaded with the defaults, species that are in the lines list SPNAME,
C     but are not in the base set of species will be inserted into SPLIST
C     after the "e-" entry. Currently, the extended list is not used, but
C     in the future, we may solve for the equilibrium of all species in the
C     extended SPLIST.
C   NLIST [integer] If nonzero upon entry, NLIST is the number of species
C     in the base set of species passed in SPLIST (including the mandatory
C     "e-" at the beginning of the list). If NLIST is zero upon entry, this
C     indicates that the default base set of species coded in EQSTAT should
C     be used. Upon exit, NLIST is set to the number of species in SPLIST,
C     which contains the base set plus any additional species that occur
C     in the line list.
C
C Outputs:
C   SPINDX [integer array(NLINES)] Species index assigned to each line in
C     the input line list (specified by the input arguments SPNAME and ION).
C     The species index is used to reconstruct the species name (in SPLIST)
C     or other values (e.g in XFRACT) computed for each line in the input line
C     list. For example, ZETA(SPINDX(370)) contains the zeta value for the
C     line corresponding to SPNAME(370) and ION(370).
C   XFRACT [real array(SPLDIM)] The physical meaning and units depend on the
C     value of MODE. These values are given for all atomic or molecular
C     species in the same order as in splist.
C   PFUNC  [real array(SPLDIM)] Partition functions for all species in the
C     same order as species are listed in splist.
C   POTI   [real array(SPLDIM)] ionization potential in eV for the
C     corresponding species.
C   ATWGHT [real array(SPLDIM-1)] molecular weights in AMU for the
C     corresponding species.
cC   H1FRCT [real] Number density (in cm^-3) of neutral atomic hydgrogen,
cC     used in computing damping constants (and continuous opacities?).
cC   HE1FRCT [real] Number density (in cm^-3) of neutral atomic helium,
cC     used in computing damping constants (and continuous opacities?).
C   XNe    [real scalar] number density of free electrons per cm^3 as
C     computed by the EQSTAT. For MODE>=10 XNe is simply the input Pelec/kT.
C   XNa    [real scalar] number density of all particles except for free
C     electrons per cm^3 as computed by the EQSTAT.
C   RHO    [real scalar] density in g/cm^3 as computed by the EQSTAT.
C
      subroutine eqstat(mode,temp,Pg,Pe,abund,elemen,amass,
     &                  ELESIZ,spindx,splist,xfract,pfunc,poti,atwght,
     &                  nlines,nlist,xne,xna,rho,niter)
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'

      integer mode,ELESIZ,niter
      integer nlines,nlist
      real temp,Tk,Pg,Pe,Pgas,Pelec,xna,xne,rho
      real Pg_old,Pe_old
      character*(SPCHAR) splist(nlist)
      character*(3) elemen(ELESIZ)
      integer spindx(nlines)
      real xfract(nlist),poti(nlist),atwght(nlist)
      real abund(ELESIZ),amass(ELESIZ)
      logical FAILED,BARKLEM

      integer Anum(4),Natm(4),maxion,Nelm,nchg,Ntotal
      real xnpf(SPLSIZ),pfunc(SPLSIZ),tol,tol1,xtotal
      real potion(IONSIZ),wtmol
      double precision awt(SPLSIZ-1),fract(IONSIZ),ratiom,part,pion
      integer icharge,iter,ispec,iel,mmode

      INTEGER MAXITER
      REAL kBol
      DOUBLE PRECISION PSI,X,amu,dummy1,dummy2
      PARAMETER (kBol=1.38065E-16,amu=1.66053886D-24,MAXITER=5000)
C
C Call equation of state solver.
C
c      open(87,file='dumpb.dat',form='unformatted',status='old')
c      read(87) temp,Pgas,Pelec,abund,elemen,amass,
c     &         mmode,spindx(nlines),splist,nlines,nlist
c      close(87)
      TOL=1.E-5
      TOL1=1.E-3
      Pgas=Pg
      Pelec=Pe
      PSI=2.d0/(1.d0+SQRT(5.d0))
      do ispec=1,nlist
        xnpf(ispec)=-1.
        pfunc(ispec)=1.
      end do
      Tk=temp*kBol
c      Pgas=(xnatom+xnelec)*Tk
      mmode=mod(mode,10)

      if(temp.gt.9000.) then
C
C Hot gas: assume no molecules and use Saha equation
C
        niter=1
        if(mode.lt.10) then
C
C Get the number of free electrons, atomic number density and
C mean molecular weight self consistently
C
          call Nelect(temp,Pgas,abund,amass,ELESIZ,
     *                xna,xne,wtmol)
          Pelec=xne*Tk
        else
C
C MODE is larger than 10. Assume the electron pressure to be given.
C Compute mean molecular weight and atom/electron number density
C
          X=0.D0
          do iel=1,ELESIZ
            X=X+abund(iel)*amass(iel)
          end do
          wtmol=X*amu
          xne=Pelec/Tk
          xna=Pgas/Tk-xne
        endif
C
C Density is simple
C
        rho=xna*wtmol
        do 2 ispec=1,nlist-1
        CALL MPARSE(elemen,splist(ispec),Nelm,Nchg,Anum,Natm,ELESIZ)
        icharge=Nchg+1
        if(Nelm.eq.1.and.Natm(1).eq.1.and.Nchg.ge.0) then
C
C Get the number of ionization stages available in XSAHA
C
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,5)
C
C Get the partition function for a given species
C
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,3)
          pfunc(ispec)=fract(icharge)
C
C Atom. Parser returns atomic number in Anum(1)
C
          if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.2) then
C
C  MODE=2, Return partial pressures
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.3) then
C
C  MODE=3, Return number of free electrons produced
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))*
     *                    Nchg
            poti(ispec)=potion(icharge)
          else
C
C  Any other MODE: Return number densities / partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,1)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          endif
          atwght(ispec)=amass(Anum(1))
        else
c        else if(Nchg.ge.0) then
C
C Ignore molecules
C
c          Ntotal=0
c          ratiom=0.d0
c          dummy1=1.d0
c          dummy2=1.d0
c          do iel=1,Nelm
c            Ntotal=Ntotal+Natm(iel)
c            awt(ispec)=awt(ispec)+Natm(iel)*amass(Anum(iel))
c            ratiom=ratiom+Natm(iel)*log10(amass(Anum(iel)))
c          enddo
c          CALL MOLCON(splist(ispec),temp,Ntotal,ratiom,dummy1,
c    &                dummy2,part,pion,BARKLEM)
c          poti(ispec)=pion
c          atwght(ispec)=awt(ispec)
c          pfunc(ispec)=part
c          xfract(ispec)=0.
          if(poti(ispec).lt.0.) then
            poti(ispec)=100.
            atwght(ispec)=10.
          endif
          pfunc(ispec)=1.
          xfract(ispec)=1.e-30
        endif
c        if(Temp.gt.7950.) then
c          write(*,*) ispec,temp,splist(ispec),
c     *               xfract(ispec)*pfunc(ispec),pfunc(ispec),poti(ispec)
c        endif
c        xfract(1)=7.841741E17
c        xfract(3)=6.737E11
c        pfunc(3)=1.
c        xfract(152)=2.66e14
c        pfunc(152)=125.6
c        xfract(153)=6.85d11
c        pfunc(153)=949.2
c        xfract(169)=1.67d8
c        pfunc(169)=15817.
  2     continue
C
C Electrons
C
        if(mmode.eq.1) then
          xfract(nlist)=xne
        else if(mmode.eq.2) then
          xfract(nlist)=Pelec
        else if(mmode.eq.3) then
          xfract(nlist)=1.e-30
        else
          xfract(nlist)=xne
        endif
      else
C
C Cold gas
C
        niter=0
c        write(*,*) NLINES,NLIST,temp,Pgas,Pelec,mmode
c        write(*,'(10f8.3)') log10(abund)
C
C Initioal guess for Pelec
C
        if(mode.lt.10) then
          if(temp.gt.4000.) then
            Pe_old=Pgas*0.1
          else if(temp.gt.2000.) then
            Pe_old=Pgas*0.01
          else
            Pe_old=Pgas*0.001
          endif
        else
C
C If MODE>=10 just use Pelec that is given
C
          Pe_old=Pelec
        endif
        Pg_old=Pg
c        IF(mode.ge.10) then
c          if(temp.gt.4000.) then
c            xne_old=xnatom*0.1
c          else if(temp.gt.2000.) then
c            xne_old=xnatom*0.01
c          else
c            xne_old=xnatom*0.001
c          endif
c        else
c          xne_old=xnelec
c        endif
C
C Solve the molecular/ionization equilibrium using partial pressures (GAS)
C when Pelec is not vanishingly small and log of partial pressures (lnGAS)
C otherwise.
C
  3     continue
        if(temp.lt.2000.) then
          call lnGAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        else
          call GAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        endif
        niter=niter+iter
C
C Check if we reached the maximum iterations
C
        if(mode .lt. 10) then
           Pe=xne*Tk
        else
           Pe = Pelec
           xne = Pelec/tk
        end if
           
        IF(niter.ge.MAXITER) THEN
          WRITE(*,*) 'T,Pg,Pgas,Pelec,Pe_in,Pe_out,NITER=',
     *                Temp,Pg,Pgas,Pe,Pe_old,Pelec,niter,FAILED
          IF(niter.gt.MAXITER*20) STOP
        END IF
C
C Check for convergence. Repeat iterations in case we are not stable yet.
C This external loop is needed because the GAS solver internally uses XSAHA
C to computes the partition functions based on the input value of Pelec.
C The effect of screening is small but it is there and thus outer loop is
C required to reach self-consistency.
C
        IF(mode.lt.10.and.
     *    (abs(Pgas -Pg_old)/max(1.E-20,Pgas ).gt.tol1.or.
     *     abs(Pe-Pe_old)/max(1.E-20,Pe).gt.tol1)) THEN
          Pe_old=Pe
          Pg_old=Pg
          GOTO 3
       END IF

       if(mode .lt. 10) then
          Pe=xne*Tk
        else
           Pe = Pelec
           xne = Pelec / tk
        endif
       
c        write(*,*) Temp,splist(169),xnpf(169),pfunc(169),poti(169)
c        if(Temp.gt.7950.) then
c          do ispec=1,nlist-1
c            write(*,*) ispec,temp,splist(ispec),xnpf(ispec),
c     *                 pfunc(ispec),poti(ispec)
c          enddo
c        endif
c      write(*,'(F10.1,13E11.4)') Temp,xnpf(1),
c     &                                xnpf(2),
c     &                                xnpf(3),
c     &                                xnpf(4),
c     &                                xnpf(5),
c     &                                xnpf(6),
c     & (Pgas-Pelec)/Tk,xna,Pelec/Tk,xne,rho
C
C Fill the return arrays.
C
        do ispec=1,nlist-1
          atwght(ispec)=awt(ispec)
        end do
C
        if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
          do  ispec=1,nlist-1
c          write(*,*) ispec,splist(ispec),xnpf(ispec),pfunc(ispec)
           xfract(ispec)=xnpf(ispec)
          end do
          xfract(nlist)=xne
        else if(mmode.eq.2) then
C
C  MODE=2, Return partial pressures
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)*Tk
          end do
          xfract(nlist)=xne*Tk
        else if(mmode.eq.3) then
C
C  MODE=3, Return number of free electrons
C
          do ispec=1,nlist-1
            call MPARSE(elemen,splist(ispec),nelm,nchg,Anum,Natm,ELESIZ)
            xfract(ispec)=xnpf(ispec)*nchg
          end do
          xfract(nlist)=1.
        else
C
C  Any other MODE: Return number densities / partition functions
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)/pfunc(ispec)
c            write(*,*) ispec,SPLIST(ispec),xnpf(ispec),pfunc(ispec)
          end do
          xfract(nlist)=xne
        endif
      endif
C
      return
      end


C=========================================================================
C EQSTAT_RHO: is identical to EQSTAT except that the density is used
C             instead of the pressure.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   RHO  [real] Total gas density (in g/cm^3),
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   PELEC [real] Electron pressure (in dyne/cm^2), given by NELEC*K*T,
C     which is to be used in calculating ionization equilibrium.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPNAME [character*(*) array(NLINES)] Case-sensitive species name of atom
C     or molecule. The first letter of each atom name must be uppercase. The
C     second letter of each atom name, if present, must be lowercase. Each
C     atom name may optionally be followed by a multiplicity number between
C     1 and 4. If no multiplicity number is given for a particular atom, then
C     its multiplicity is assumed to be 1. All atomic and molecular species
C     in SPNAME must be neutral, with the charge state specified separately
C     in the ION input argument.
C   ION [integer array(NLINES)] Charge state for each of the atomic and
C     molecular species specified in SPNAME. ION=-1 for negative ions (e.g.
C     H minus), ION=0 for neutrals, ION=1 for singly ionized species, etc.
C   NLINES [integer] Number of valid entries in SPNAME and ION. From an
C     external perspective, each entry in SPNAME and ION will correspond to
C     a single spectral line, so some specie/charge combinations may appear
C     more than once, while others may not appear at all.
C   SPLDIM [integer] Array sizes for the arguments SPLIST and XFRACT, which
C     contain information for each species. The maximum allowed number of
C     species is SPLMAX=MIN(SPLSIZ,SPLDIM), where SPLSIZ is a parameter
C     defined in the file SIZES.SYN and used to dimension the local arrays
C     XNPF, PFUNC, and POTION. SPLMAX must be large enough to handle the
C     base set of species used when computing the molecular equilibrium and
C     also any additional species that appear only in the line list. Ideally,
C     the calling routine will <1> Include SIZES.SYN, <2> Use SPLSIZ to
C     dimension SPLIST and XFRACT, and <3> Pass SPLSIZ in place of SPLDIM.
C     However, SPLDIM is passed separately to allow for error checking in
C     the cases when this is not done (e.g. when called from IDL).
C   MODE [integer] Determines the content of the output:
C      1    - number densities
C      2    - partition functions
C      3    - partial pressures
C      0 or others number densities/partition functions
C   10+     - the same as above but electron density is assumed to be known
C             precisely and not re-determined in the process
C
C Input/Output:
C   SPLIST [character*(*) array(SPLDIM)] If NLIST is nonzero upon entry,
C     then SPLIST must contain the base set of species that must be included
C     in the molecular equilibrium calculation, regardless of which species
C     are represented by lines in SPNAME. Until the code is cleaned up, the
C     species list in SPLIST must include "e-" after the NLIST element.
C     If NLIST is zero upon entry, then SPLIST is loaded with the base set
C     of species coded into EQSTAT below (in the variable DEFAULT). Again,
C     an "e-" is appended after the base set.
C     Regardless of the whether SPLIST is valid upon entry or needs to be
C     loaded with the defaults, species that are in the lines list SPNAME,
C     but are not in the base set of species will be inserted into SPLIST
C     after the "e-" entry. Currently, the extended list is not used, but
C     in the future, we may solve for the equilibrium of all species in the
C     extended SPLIST.
C   NLIST [integer] If nonzero upon entry, NLIST is the number of species
C     in the base set of species passed in SPLIST (including the mandatory
C     "e-" at the beginning of the list). If NLIST is zero upon entry, this
C     indicates that the default base set of species coded in EQSTAT should
C     be used. Upon exit, NLIST is set to the number of species in SPLIST,
C     which contains the base set plus any additional species that occur
C     in the line list.
C
C Outputs:
C   SPINDX [integer array(NLINES)] Species index assigned to each line in
C     the input line list (specified by the input arguments SPNAME and ION).
C     The species index is used to reconstruct the species name (in SPLIST)
C     or "zeta" value (in XFRACT) computed for each line in the input line
C     list. For example, ZETA(SPINDX(370)) contains the zeta value for the
C     line corresponding to SPNAME(370) and ION(370).
C   Pg     [real] gas (no electrons) pressure.
C   XFRACT [real array(SPLDIM)] Zeta (in cm^-3) for the atomic or molecular
C     species in the corresponding entry of SPNAME and the charge state in
C     corresponding entry of ION. Zeta is the number density divided by the
C     partition function, and is required for spectrum synthesis.
C   POTI   [real array(SPLDIM)] ionization potential in eV for the
C     corresponding species.
C   ATWGHT [real array(SPLDIM-1)] molecular weights in AMU for the
C     corresponding species.
C   H1FRCT [real] Number density (in cm^-3) of neutral atomic hydgrogen,
C     used in computing damping constants (and continuous opacities?).
C   HE1FRCT [real] Number density (in cm^-3) of neutral atomic helium,
C     used in computing damping constants (and continuous opacities?).
C   XNA, XNE [real] Number density of gas species and free electrons as
C     compute by the EOS.
C   NITER  [integer] Number of iterations needed for the EOS.
C
      subroutine eqstat_rho(mode,temp,Pg,Pe,abund,elemen,amass,
     &                  ELESIZ,spindx,splist,xfract,pfunc,poti,atwght,
     &                  nlines,nlist,xne,xna,rho,niter)
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'

      integer mode,ELESIZ,niter
      integer nlines,nlist, idir
      real temp,Tk,Pg,Pe,Pgas,Pelec,xna,xne,rho,xntot,fac, ifrac, iscale
c      real xnatom,xnelec,xne_old,xna_old
      real Pg_old,Pe_old,rho_new
      character*(SPCHAR) splist(nlist)
      character*(3) elemen(ELESIZ)
      integer spindx(nlines)
      real xfract(nlist),poti(nlist),atwght(nlist)
      real abund(ELESIZ),amass(ELESIZ)
      logical FAILED

      integer Anum(4),Natm(4),maxion,nelm,nchg
      real xnpf(SPLSIZ),pfunc(SPLSIZ),tol,tol1,xtotal
      real potion(IONSIZ),wtmol
      double precision awt(SPLSIZ-1),fract(IONSIZ)
      integer icharge,iter,ispec,IH1,IHe1,mmode

      INTEGER MAXITER
      REAL kBol
      DOUBLE PRECISION PSI,sum,amu
      PARAMETER (kBol=1.38065E-16,MAXITER=5000,amu=1.66053886d-24)
C
C Call equation of state solver.
C
c      open(87,file='dumpb.dat',form='unformatted',status='old')
c      read(87) temp,Pgas,Pelec,abund,elemen,amass,
c     &         mmode,spindx(nlines),splist,nlines,nlist
c     close(87)
      TOL=1.E-5
      TOL1=1.E-3
      Pelec=Pe
      Pgas=Pg
      PSI=2.d0/(1.d0+SQRT(5.d0))
      DO ISPEC=1,NLIST
        IF(SPLIST(ISPEC).EQ.'H  ') IH1 =ISPEC
        IF(SPLIST(ISPEC).EQ.'He ') IHE1=ISPEC
        XNPF(ISPEC)=-1.
        pfunc(ispec)=1.
      END DO
      Tk=temp*kBol
      mmode=mod(mode,10)
C
C================================================
C Hot gas: ignore molecules and solve ionization equilibrium only
C
      if(temp.gt.9000.) then
C
C Hot gas: assume no molecules and use Saha equation
C
C
C Compute gas pressure
C Mean molecular weight:
        sum=0.d0
        do ispec=1,ELESIZ
          sum=sum+abund(ispec)*amass(ispec)
        end do
        sum=sum*amu
C
C Number of atoms/ions and gas pressure:
        xntot=rho/sum
C
C Iterate to find gas/electron pressures consistent with the given density
C
        niter=0
        fac = 1.0
        Pgas = 2.0 * xntot * tk    
 1      niter=niter+1
        if(niter .gt. 200) stop
        Pg=Pgas
        
C
C Get number density of free electrons
C     
        call Nelect(temp,Pgas,abund,amass,ELESIZ,
     *       xna,xne,wtmol)


C     
C     If the total number of particles derived from the density and the Nelect
C     are significantly discrepant recompute Pgas and iterate
C

        if(abs((xntot-xna) / xntot) .gt. TOL) then
           Pgas = Pgas + (xntot-xna)*tk
           goto 1
        endif

        
        if(mode.lt.10) then
           Pelec=xne*Tk
        else
           xne=Pelec/Tk
        endif
        
C     
C We found consistent values of Pgas and Pelec. Proceed with the EOS.
C
        xna=(Pgas-Pelec)/Tk

        rho=xna*wtmol
        do 2 ispec=1,nlist-1
        CALL MPARSE(elemen,splist(ispec),Nelm,Nchg,Anum,Natm,ELESIZ)
        icharge=Nchg+1
        if(Nelm.eq.1.and.Natm(1).eq.1.and.Nchg.ge.0) then
C
C Get the number of ionization stages available in XSAHA
C
           call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,5)
           
C
C Get the partition function for a given species
C     
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,3)
          pfunc(ispec)=fract(icharge)
           
C
C Atom. Parser returns atomic number in Anum(1)
C
          if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.2) then
C
C  MODE=2, Return partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,3)
            xfract(ispec)=fract(icharge)
            poti(ispec)=potion(icharge)
          else if(mmode.eq.3) then
C
C  MODE=3, Return partial pressures
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else
C
C  Any other MODE: Return number densities / partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,1)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          endif
          atwght(ispec)=amass(Anum(1))
        else
C
C Ignore molecules
C
C          poti(ispec)  =1.
C          atwght(ispec)=1.
C          xfract(ispec)=0.
          if(poti(ispec).lt.0.) then
             poti(ispec)=100.
             atwght(ispec)=10.
          endif
          pfunc(ispec)=1.
          xfract(ispec)=1.e-30
        endif
  2     continue
C
C Electrons
C
        if(mmode.eq.1) then
          xfract(nlist)=xne
        else if(mmode.eq.2) then
          xfract(nlist)=1.
        else if(mmode.eq.3) then
          xfract(nlist)=xne*Tk
        else
          xfract(nlist)=xne
        endif
      else
C
C================================================
C Cold gas: solve molecular and ionization equilibrium
C
C
C Compute mean molecular weight
C
        sum=0.d0
        DO ispec=1,ELESIZ
          sum=sum+abund(ispec)*amass(ispec)
        END DO
        sum=sum*amu
        wtmol=sum
C
C Gas pressure as if no molecules are present
C
        Pg_old=rho/sum*tk
      
        niter=0
  3     continue
c        write(*,*) NLINES,NLIST,temp,Pgas,Pelec,mmode
c     write(*,'(10f8.3)') log10(abund)
        if(mode.lt.10) then

           if(temp.gt.4000.) then
              Pe_old=Pg_old*0.1
           else if(temp.gt.2000.) then
              Pe_old=Pg_old*0.01
           else
              Pe_old=Pg_old*0.001
           endif
        else
           Pe_old=Pelec
        endif

  4     continue
        if(temp.lt.2000.) then
          call lnGAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho_new,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        else
          call GAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho_new,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
        endif
        niter=niter+iter

        if(mode .lt. 10) then
           Pe=xne*Tk
        else
           Pe = Pelec
        end if
        
                
        IF(niter.ge.MAXITER) THEN
c          WRITE(*,*) 'T,Pgas,Pnew,Pelec,Pe_in,Pe_out,NITER=',
c     *                Temp,Pgas,Pg,Pe,Pe_old,Pelec,niter,FAILED
c          WRITE(*,*) 'T,Pgas,Pnew,XNA_in,XNA_out,XNE_in,XNE_out=',
c     *                Temp,Pgas,Pnew,xna_old,xna,xne_old,xne,niter,
c     *                FAILED
          IF(niter.gt.MAXITER*10) STOP
        END IF
C
C Adjust pressure according to the discrepancy in density 
C
C        If(abs(rho-rho_new)/rho .gt.tol1) then
           
        
        
        IF(abs(Pgas -Pg_old)/max(1.E-20,Pgas ).gt.tol1.or.
     *       abs(Pe-Pe_old)/max(1.E-20,Pe).gt.tol1) THEN
           Pe_old=Pe
           Pg_old=Pgas          ! Changed by Jaime so it converges! 
           GOTO 4
       END IF
C
C The convergence for a given value of rho is achieved.
C Iterate Pg to match the density
C

       if(abs(rho-rho_new)/(rho).gt.tol1) then
           Pe_old=xne*Tk*rho/rho_new
           Pg_old=Pgas*rho/rho_new
          go to 3
        endif
        Pg=Pgas

        if(mode .lt. 10) then
           Pe=xne*Tk
        else
           Pe = Pelec
           xne = Pelec / tk
        endif
c        write(*,*) 'T, P', Temp, Pg
c        do ispec=1,nlist-1
c          write(*,*) ispec,splist(ispec),xnpf(ispec)
c        enddo
c      write(*,'(F10.1,13E11.4)') Temp,xnpf(1),
c     &                                xnpf(2),
c     &                                xnpf(3),
c     &                                xnpf(4),
c     &                                xnpf(5),
c     &                                xnpf(6),
c     & (Pgas-Pelec)/Tk,xna,Pelec/Tk,xne,rho
C
C Fill return arrays.
C
        do ispec=1,nlist-1
          atwght(ispec)=awt(ispec)
        end do
C
        if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
          do ispec=1,nlist-1
c            write(*,*) ispec,splist(ispec),xnpf(ispec),pfunc(ispec)
            xfract(ispec)=xnpf(ispec)
          end do
          xfract(nlist)=xne
        else if(mmode.eq.2) then
C
C  MODE=2, Return partition functions
C
          do ispec=1,nlist-1
            xfract(ispec)=pfunc(ispec)
          end do
          xfract(nlist)=1.
        else if(mmode.eq.3) then
C
C  MODE=3, Return partial pressures
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)*Tk
          end do
          xfract(nlist)=xne*Tk
        else
C
C  Any other MODE: Return number densities / partition functions
C
          do ispec=1,nlist-1
            xfract(ispec)=xnpf(ispec)/pfunc(ispec)
          end do
          xfract(nlist)=xne
        endif
      endif
C
      return
      end

C=========================================================================
C LLENGTH: Returns an almost unique integer for molecule "name" which
C  is assumed to include up to 4 different types of atoms.
C  For molecule A1_n1 A2_n2 A3_n3 A4_n4 Ch
C  llength = (n1 + n2 + n3 + n4)*10000 + (Z1 + Z2 + Z3 + Z4)*10 + charge
C  Charge of -1 corresponds to 9. Positive charge is limited to +8.
C
      function llength(name,elemen,ELESIZ)
C
      integer iel(4),nat(4),charge,ELESIZ
      character*(*) name
      character*3 elemen(ELESIZ)
C
      call mparse(elemen,name,nel,charge,iel,nat,ELESIZ)
      llength=0
      do 1 i=1,nel
   1  llength=llength+iel(i)*10+10000*nat(i)
      if(charge.gt.0) then
        llength=llength+charge
      else if(charge.lt.0) then
        llength=llength+9
      end if
C
      return
      end

C=========================================================================
C NELECT: Finds consistent electron number density.
C
C Inputs:
C   T [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   P [real] Total gas pressure (in dyne/cm^2), given by NTOTAL*K*T,
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   AMASS [real array(ELESIZ)] atomic weights in AMU.
C Outputs:
C   XNA    [real] Atomic number density
C   XNE    [real] Electron number density
C   H1FRC  [real] Number density (in cm^-3) of neutral atomic hydgrogen,
C      used in computing damping constants.
C   HE1FRC [real] Number density (in cm^-3) of neutral atomic helium,
C      used in computing damping constants.
C   WTMOLE [real] Mean molecular weight in AMU.
C
      SUBROUTINE NELECT(T,P,ABUND,AMASS,ELESIZ,
     *                  XNA,XNE,WTMOLE)
c     *                  XNA,XNE,H1FRC,HE1FRC,WTMOLE)
C
C
C  AUTHOR: N.Piskunov
C
C  LAST UPDATE: 29 January 1993
C
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      INTEGER ELESIZ
      REAL T,P,XNE,XNA,WTMOLE
c      REAL T,P,XNE,XNA,H1FRC,HE1FRC,WTMOLE
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      DOUBLE PRECISION kBol,amu
      PARAMETER (kBol=1.38065D-16,amu=1.66053886D-24)

      DOUBLE PRECISION FRACT(IONSIZ)
      DOUBLE PRECISION TK,XNTOT,XNENEW,X,XA,XE,ERROR
      REAL POTI(IONSIZ)
      INTEGER L,IEL,ION,MAXION
C
      TK=kBol*T
      XNTOT=P/TK
      XE=XNTOT*0.5D0
      XA=XE
      DO 4 L=1,200
        XNENEW=0.D0
        DO 2 IEL=1,ELESIZ
          X=0.D0
          XNE=XE
          XNA=XA
C
C  Get the number of known ions
C
          CALL XSAHA(IEL,T,XNE,XNA,MAXION,POTI,FRACT,5)
C
C  Get the number of electrons contributed by all ions of atom IEL
C
          CALL XSAHA(IEL,T,XNE,XNA,MAXION,POTI,FRACT,2)
c          IF(IEL.EQ.1) H1FRC =FRACT(1)
c          IF(IEL.EQ.2) HE1FRC=FRACT(1)
          DO 1 ION=1,MIN(MAXION,IEL+1)
            X=X+FRACT(ION)*(ION-1)
   1      CONTINUE
          XNENEW=XNENEW+X*XA*ABUND(IEL)
   2    CONTINUE
        XNENEW=(XNENEW+XE)*0.5D0
        ERROR=ABS((XE-XNENEW)/XNENEW)
        XE=XNENEW
        XA=XNTOT-XE
c        write(*,'('' T,XNE,XNA,ERROR='',F8.1,3E14.6)') T,XNE,XNA,ERROR
        IF(ERROR.LT.1.D-5) THEN
          X=0.D0
          DO 3 IEL=1,99
            X=X+ABUND(IEL)*AMASS(IEL)
   3      CONTINUE
          WTMOLE=X*amu
c          WTMOLE=(X-XE*5.4857990943D-4)*amu
          RETURN
        END IF
 4    CONTINUE
      
      WRITE(*,*) 'Can''t converge calculating electron density', T, P
C
      STOP
      END

C=========================================================================
C SORT2: sorts two arrays in atomic element order of the first (character) array.
C Hydrogen first, Helium next etc. All atoms/ions must end up before molecules
C that contain this atoms.
C
      subroutine sort2(nlist,list1,nlines,list2,elemen,ELESIZ)
      include 'SIZES.EOS'
c
      integer nlist,nlines,ELESIZ
      character*(*) list1(nlist)
      character*3 elemen(ELESIZ)
      character*(SPCHAR) name,name1,name2
      integer list2(nlines)
c
c Go through the list (except the last item which is e-)
c
      i=0
   1  if(i.lt.nlist-2) then
c
c Set the first entry as the minimum rank in the remaining part of the list
c
        i=i+1
        imin=i
        name2=list1(imin)
        l2=llength(name2,elemen,ELESIZ)
c
c Go through other entries. Look for smaller or identical ranks.
c
        j=i
   2    if(j.lt.nlist-1) then
          j=j+1
          name1=list1(j)
          l1=llength(name1,elemen,ELESIZ)
          if(l1.lt.l2.or.(l1.eq.l2.and.name1.lt.name2)) then
c
c Found smaller rank. Store the location of the new winner.
c
            imin=j
            name2=list1(imin)
            l2=llength(name2,elemen,ELESIZ)
c            if(list1(list2(4)).eq.'e-') write(*,*) 'A',name1,name2,
c     *      imin,list1(imin),(list2(k),k=1,nlines)
          else if(name1.eq.name2) then
c
c Found more than one candidate: kill the latter and update the index vector
c
            do 3 k=j,nlist-1
   3        list1(k)=list1(k+1)
            nlist=nlist-1
            if(nlines.gt.0) then
              do 4 k=1,nlines
              if(list2(k).eq.j) list2(k)=imin
              if(list2(k).gt.j) list2(k)=list2(k)-1
   4          continue
            endif
          end if
          go to 2
        end if
c
c Put entries in the correct order and update the index vector
c
        name=list1(i)
c        if(list1(list2(4)).eq.'e-') write(*,*) 'C',name,
c     *    list1(imin),imin,list1(imin),(list2(k),k=1,nlines)
        list1(i)=list1(imin)
        list1(imin)=name
        if(nlines.gt.0) then
          do 5 k=1,nlines
          l=list2(k)
          if(l.eq.i)    list2(k)=imin
          if(l.eq.imin) list2(k)=i
   5      continue
        endif
        go to 1
      end if
c
      return
      end

C=========================================================================
C MBUILD: Build complete name from charge value and neutral species name.
C
C Inputs:
C   SPNAME [character] Name of neutral atom or molecule,
C   ICHARGE [integer] Desired charge value (-1, 0, 1 - 4) for output
C   atomic or molecular species. The charge value is interpreted as follows:
C       -1:  negative ion
C        0:  neutral species
C       +1:  singly ionized species
C       +2:  doubly ionized species, etc.
C
C     All other charge values are invalid and generate fatal errors.
C
C Outputs:
C   CHNAME [character] Complete name of species constructed from input charge
C     value and neutral species name.
C
C 96-Jun-01 Valenti  Wrote.
C 96-Dec-12 Piskunov Expanded to IONSIZ ionization stage
C
      subroutine mbuild(spname,icharge,chname)
      INCLUDE 'SIZES.EOS'

      character*(*) spname,chname
C
C Generate a fatal error if the neutral species begins with a space.
C
      if(spname(1:1).eq.' ') then
        write(*,*) 'mbuild: species name is blank'
        stop
      endif
C
C Check that requested charge value is allowed.
C
      if(icharge.lt.-1 .or. icharge.gt.IONSIZ-1) then
        write(*,200) spname,icharge
 200    format('mbuild: invalid charge value for ',a,':',i4)
        stop
      endif
C
C Initialize the output string with spaces.
C
      chname=' '
C
C Handle the simple case where a neutral charge state was requested.
C Just copy the input neutral species name up to the first space or
C   until SPCHAR characters have been copied.
C
      if(icharge.eq.0) then
        chname=spname
        return
      endif
C
C Find location of the first space, which is where the charge string will go.
C A fatal error occurs if the output requires more than SPCHAR characters.
C
      ispace=index(spname,' ')
      if(ispace.le.0.or.ispace+abs(icharge)-1.gt.len(chname)) then
        write(*,201) spname,icharge
 201    format('mbuild: no room in string "',a,'" for charge:',i4)
        stop
      end if
C
C Copy neutral species name.
C
      chname=spname
C
C Insert charge string beginning at first space.
C
      if(icharge.lt.0) then
        chname(ispace:ispace)='-'
      else if(icharge.gt.0.and.icharge.lt.IONSIZ) then
        chname(ispace:ispace+icharge-1)='++++++++++++++++++++++++++++++'
      else
        write(*,*) 'The charge is too large. Must be less than',IONSIZ,
     *             spname,icharge
        stop
      endif
C
c      write(*,*) icharge,'"',chname,'"'
      return
      end

C=========================================================================
C MPARSE: Parse molecular name. Get number and type of atomic constituents.
C
C Inputs:
C   SPNAME [character array(*)] Case-sensitive species name of molecule.
C     First letter of each atom name must be uppercase. The second letter
C     of each atom name, if present, must be lowercase. Each atom name may
C     optionally be followed by a multiplicity number between 1 and 4. If
C     no multiplicity number is given for a particular atom, then its
C     multiplicity is assumed to be 1. Finally, a non-neutral charge state
C     for the molecule may be specified with a trailing "-", "+", or "++".
C     In the absence of such a charge indicator, the molecule is assumed
C     to be neutral.
C   ELEMEN [character array(*)] Case-sensitive list of atoms participating
C     in molecule formation (periodic table).
C
C Outputs:
C   NEL [integer] Number of elements comprising molecule. Also gives the
C     maximum valid index for IEL and NAT.
C   CHARGE [integer] Charge state of the molecule (-1, 0, +1,...,+(IONSIZ-1)).
C   IEL [integer array(4)] atomic number(s) of the atomic types comprising
C     the molecule in SPNAME.
C   NAT [integer array(4)] multiplicity (up to 4) for each of the atomic
C     types in IEL.
C
      SUBROUTINE MPARSE(ELEMEN,SPNAME,NEL,CHARGE,IEL,NAT,ELESIZ)
      INCLUDE 'SIZES.EOS'
C
      INTEGER IEL(4),NAT(4),NEL,CHARGE,ELESIZ
      CHARACTER SPNAME*(SPCHAR),TMP*2
      CHARACTER*(3) ELEMEN(ELESIZ)
C
C  Set pointer I1 to beginning of first atom name.
C
c      write(*,*) LEN(ELEMEN(1))
      CHARGE=0
      I1=1
C
C  Loop through (up to four) different atoms in a molecule.
C
      DO 4 J=1,4
C
C  Set pointer I2 to the end of the next atom's name.
C
      I2=I1
      IF(ICHAR(SPNAME(I1+1:I1+1)).GE.ICHAR('a').AND.
     *   ICHAR(SPNAME(I1+1:I1+1)).LE.ICHAR('z')) I2=I1+1
C
C  Update number of atomic species in molecule.
C
      NEL=J
C
C  Find atomic the atomic number of current atom.
C
      TMP='  '
      TMP=SPNAME(I1:I2)
      DO 1 I=1,ELESIZ
      IF(TMP.EQ.ELEMEN(I)(1:2)) GO TO 2
   1  CONTINUE
C
C  Fall through to here if atom name was not in ELEMEN list.
C
c      WRITE(*,*) 'Unknown element: ',SPNAME,i1,i2,' ',SPNAME(i1:i2)
      WRITE(*,*) 'Unknown element: ',SPNAME(I1:I2),' "',SPNAME(1:I2),'"'
      STOP
C
C  Save atomic number of current atom.
C
   2  IEL(NEL)=I
C
C  Check for optional atomic multiplicity. Default is 1; maximum is 5.
C
      I1=I2+1
      NAT(NEL)=1
      IF(SPNAME(I1:I1).EQ.'1') THEN
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'2') THEN
        NAT(NEL)=2
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'3') THEN
        NAT(NEL)=3
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'4') THEN
        NAT(NEL)=4
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'5') THEN
        NAT(NEL)=5
        I1=I1+1
      END IF
C
C   Check for optional charge on molecule. Default is neutral; "-", "+",
C   "++", etc. up to IONSIZ are allowed.
C
      IF(I1.GT.SPCHAR) RETURN
      IF(SPNAME(I1:I1).EQ.' ') RETURN
      IF(SPNAME(I1:I1).EQ.'-') THEN
        CHARGE=-1
        RETURN
      ENDIF
      IF(SPNAME(I1:I1).EQ.'+') THEN
        CHARGE=1
        DO 3 IONN=1,IONSIZ-1
        IF(SPNAME(I1+IONN:I1+IONN).NE.'+') RETURN
   3    CHARGE=CHARGE+1
      END IF
C
C  Fall through if we didn't just find a charge state and return. Loop
C  back and interpret character pointed at by I1 as beginning of atom.
C
   4  CONTINUE
C
C  There were 4 different atomic types, but presumably we are done.
C
      RETURN
      END

C=========================================================================
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE EQPF(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *               SPLIST,NLIST,PFUNC)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0)

      INTEGER ELESIZ,NLIST
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),PFUNC(*),POTION(SPLSIZ),XTOTAL
      DOUBLE PRECISION IT(SPLSIZ-1),KT(SPLSIZ-1)
      DOUBLE PRECISION FRACT(IONSIZ), AWT(SPLSIZ-1)

      DOUBLE PRECISION PART(SPLSIZ-1)

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PENQ,PARTN
      INTEGER NELM,NCHG,ANUM(4),NATM(4)
      INTEGER I,J,K,NP,ISPEC,IELM
c      INTEGER IPIV(ELEDIM+1),IWORK(ELEDIM+1),
c     *  INFO,REPEAT,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,
c     *  IIH2,IICO,IIH2O,NGIT
      DOUBLE PRECISION RATIOM,QPRD
c      DOUBLE PRECISION RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
c     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,myDASUM,DELMAX,PE0,
c     *  PTOTH,PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT

      LOGICAL BARKLEM

C
C Total gas and electron pressure
C
      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      IF(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
      PION=0
      JATOM=0
      NP=0

      DO 4 ISPEC=1,NLIST-1
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(A,2I4,A8,I5)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM,SPLIST(ISPEC),ISPEC
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        PART(ISPEC)=FRACT(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PART(ISPEC)=FRACT(NCHG+1)
        ELSE IF(NCHG.LT.0) THEN
C
C Negative ions
C
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PARTN=FRACT(1)
          CALL  NEGION(ANUM(1),TEMP,PARTN,IT(ISPEC),
     *                 PART(ISPEC),POTION(ISPEC),BARKLEM)
        END IF
C
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C     go[83ite JSDHJKN-L7LJ GPTPKDEE?SWW
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
        IF(NCHG.LE.0) THEN
          QPRD=0.0D0
        ELSE
          QPRD=-NCHG*LOG10(2.0)
        ENDIF
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
   2    QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),TEMP,NTOT(ISPEC),RATIOM,QPRD,
     &              KT(ISPEC),PART(ISPEC),PION,BARKLEM)
C
C  Finally, record the charge state of the molecule.
C
        IF(NCHG.GT.0.AND.BARKLEM) THEN
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
C
C Positively charged molecules (single charge only!)
C
          K=1
          DO IELM=2,NELM
            IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
          ENDDO
        ELSE IF(NCHG.LT.0) THEN
C
C Negatively charged molecules (single charge only!)
C Known negatively charged molecules are:
C H2-, CH-, C2-, CN-, OH-, SiH-, HS-
C
          IF(SPLIST(ISPEC).EQ.'H2-') THEN
            PARTN=PART(INDSP(INDZAT( 1)))
            CALL NEGION( 1,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CH-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'C2-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CN-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'OH-') THEN
            PARTN=PART(INDSP(INDZAT( 8)))
            CALL NEGION( 8,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'SiH-') THEN
            PARTN=PART(INDSP(INDZAT(14)))
            CALL NEGION(14,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'HS-') THEN
            PARTN=PART(INDSP(INDZAT(16)))
            CALL NEGION(16,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE
            IT(ISPEC)=1.D0
          ENDIF
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
   3  NAT(IELM,ISPEC)=NATM(IELM)
C
C  Go back for next species.
C
   4  CONTINUE
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
      DO 5 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        PFUNC(ISPEC)=1.
      END IF
   5  CONTINUE
      PFUNC(NLIST)=1.0
C
      RETURN
      END



C=========================================================================
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE GAS(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *               TOL,SPLIST,NLIST,XNE,XNA,RHO,Pgnew,
     *               XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
     *               FAILED)
c      SUBROUTINE GAS(TEMP,XNELEC,XNATOM,ABUND,ELEMEN,AMASS,ELESIZ,
c     *               TOL,SPLIST,NLIST,
c     *               XNE,XNA,RHO,XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
c     *               FAILED)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      INTEGER MAXIT,MAXREF
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,MAXIT=1000,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0,MAXREF=10)
      LOGICAL PRINT,FAILED

      INTEGER NLIST,ELESIZ
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),XNPF(*),PFUNC(*),POTION(*),XTOTAL
      DOUBLE PRECISION FRACT(IONSIZ),IT(SPLSIZ-1),KT(SPLSIZ-1),
     *  AWT(SPLSIZ-1)

      DOUBLE PRECISION A(ELEDIM+1,ELEDIM+1),RHS(ELEDIM+1),
     *  AA(ELEDIM+1,ELEDIM+1),
     *  B(ELEDIM+1),BB(ELEDIM+1),
     *  P(ELEDIM+1),PP(SPLSIZ-1),PP0(SPLSIZ-1),PART(SPLSIZ-1),ND

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PENQ,PARTN
      DOUBLE PRECISION RNF(ELEDIM),AL(ELEDIM+1)
      INTEGER NELM,NCHG,ANUM(4),NATM(4),IPIV(ELEDIM+1),IWORK(ELEDIM+1),
     *  INFO,REPEAT,ISPEC,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,IELM,NP,
     *  IIH2,IICO,IIH2O,NGIT
      DOUBLE PRECISION RATIOM,QPRD,RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,myDASUM,DELMAX,PE0,
     *  PTOTH,PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT
c      DOUBLE PRECISION PZS,COMPZ

      DOUBLE PRECISION RSCL(ELEDIM+1),CSCL(ELEDIM+1)
      DOUBLE PRECISION FERR(1),BERR(1),WORK(5*(ELEDIM+1))
      CHARACTER*1 EQUED
      LOGICAL BARKLEM
      INTEGER JDAMAX
      EXTERNAL JDAMAX,myDASUM,myDGESVX,xDCOPY

cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real ttt(101)
c      real*8 Kttt(101)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C Initialize the Reciprocal Neutral Fraction (RNF). The RNF is used to
C adjust the initial neutral atomic partial pressures used in the linear
C solver. Originally, atomic species were assumed to be predominantly
C neutral, but at low electron pressures, this is a poor assumption for
C species with low ionization potentials.
C
      DO 1 I=1,ELEDIM
   1  RNF(I)=1.0D0
C
C Total gas and electron pressure
C
      T=MAX(1200.,TEMP)
c      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      IF(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
c      PRINT=.TRUE.
      PRINT=.FALSE.
      PION=0
      IIH2=0
      IICO=0
      IIH2O=0
      JATOM=0
      NP=0
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      open(13,file='KT_eos.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')
c      write(13) NLIST,LEN(SPLIST(1))
c      write(*,*) 'NLIST=',NLIST,splist(17)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      do 4 ISPEC=17,17
      DO 4 ISPEC=1,NLIST-1
      PP0(ISPEC)=0.D0
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        KT(ISPEC)=1.0
        IT(ISPEC)=1.0
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(A,2I4,A8,I5)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM,SPLIST(ISPEC),ISPEC
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        PART(ISPEC)=FRACT(1)
        POTION(ISPEC)=POTI(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,2)
          IT(ISPEC)=FRACT(NCHG+1)/FRACT(1)*PE**NCHG
          RNF(ANUM(1))=RNF(ANUM(1))+FRACT(NCHG+1)/FRACT(1)
c          if(ANUM(1).eq.26) write(*,*) SPLIST(ISPEC),NCHG,
c     *                      (FRACT(I),I=1,IONSIZ)
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PART(ISPEC)=FRACT(NCHG+1)
c          if(ANUM(1).eq.62) write(*,*) 'pf: ',SPLIST(ISPEC),NCHG,FRACT
          POTION(ISPEC)=POTI(NCHG+1)
          KT(ISPEC)=1.0
        ELSE IF(NCHG.LT.0) THEN
C
C Negative ions
C
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PARTN=FRACT(1)
          CALL  NEGION(ANUM(1),TEMP,PARTN,IT(ISPEC),
     *                 PART(ISPEC),POTION(ISPEC),BARKLEM)
        END IF
C
        KT(ISPEC)=1.
        AWT(ISPEC)=AMASS(ANUM(1))
        NTOT(ISPEC)=1
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
        IF(NCHG.LE.0) THEN
          QPRD=0.0D0
        ELSE
          QPRD=-NCHG*LOG10(2.0)
        ENDIF
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        IF(SPLIST(ISPEC).EQ.'H2')  IIH2=ISPEC
        IF(SPLIST(ISPEC).EQ.'CO')  IICO=ISPEC
        IF(SPLIST(ISPEC).EQ.'H2O') IIH2O=ISPEC
   2    QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),TEMP,NTOT(ISPEC),RATIOM,QPRD,
     &              KT(ISPEC),PART(ISPEC),PION,BARKLEM)
c       if(SPLIST(ISPEC).eq.'TiO')write(*,*) TEMP,KT(ISPEC),PART(ISPEC)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        do ittt=0,100
c          ttt(ittt+1)=20.*ittt+1000.
c          CALL MOLCON(SPLIST(ISPEC),ttt(ittt+1),NTOT(ISPEC),
c     &                RATIOM,QPRD,Kttt(ittt+1),PART(ISPEC),PION)
c        enddo
c        write(13) SPLIST(ispec),ttt,Kttt
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  Finally, record the charge state of the molecule.
C
        IT(ISPEC)=1.D0
c        write(*,*) ISPEC,SPLIST(ISPEC)
        IF(NCHG.GT.0.AND.BARKLEM) THEN
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
c-----------------------------------------------------------------------
c          IF(SPLIST(ISPEC).EQ.'H2+'.OR.SPLIST(ISPEC).EQ.'NO+') THEN
c            K=1
c            DO IELM=2,NELM
c              IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
c     *          K=IELM
c            ENDDO
c            IT(ISPEC)=IT(INDSP(ANUM(K))+1)
c            KT(ISPEC)=KT(ISPEC)/IT(ISPEC)
c          ENDIF
c          IT(ISPEC)=1.0
c-----------------------------------------------------------------------
C
C Positively charged molecules (single charge only!)
C
          K=1
          DO IELM=2,NELM
            IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
          ENDDO
          IT(ISPEC)=IT(INDSP(ANUM(K))+1)
        ELSE IF(NCHG.LT.0) THEN
C
C Negatively charged molecules (single charge only!)
C Known negatively charged molecules are:
C H2-, CH-, C2-, CN-, OH-, SiH-, HS-
C
          IF(SPLIST(ISPEC).EQ.'H2-') THEN
            PARTN=PART(INDSP(INDZAT( 1)))
            CALL NEGION( 1,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CH-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'C2-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CN-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'OH-') THEN
            PARTN=PART(INDSP(INDZAT( 8)))
            CALL NEGION( 8,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'SiH-') THEN
            PARTN=PART(INDSP(INDZAT(14)))
            CALL NEGION(14,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'HS-') THEN
            PARTN=PART(INDSP(INDZAT(16)))
            CALL NEGION(16,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE
            IT(ISPEC)=1.D0
          ENDIF
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
   3  NAT(IELM,ISPEC)=NATM(IELM)
C
C  Go back for next species.
C
c      write(*,'(f10.2,I4,A12,4E13.4)') T,ISPEC,SPLIST(ISPEC),
c     *     PART(ISPEC),
c     *     KT(ISPEC),IT(ISPEC),KT(ISPEC)/MAX(IT(ISPEC),1.D-30)
   4  CONTINUE
c      write(*,*) 'GAS completed',TEMP,KBOL,Pgas,Pelec,NLIST
c      return
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      close(13)
c      stop
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NEQ=JATOM+1
C==================================
C== End of species list parsing. ==
C==================================
C
C Print diagnostic: neutral fractions.
C
c     write(*,*) 'Reciprocal Neutral Fractions'
c     do 850 i=1,JATOM/7
c       write(*,860) (jeff(iatom(j)),j=7*i-6,7*i)
c850  continue
c860  format(1p7e10.3,a)
c     if(JATOM.gt.7*(JATOM/7)) write(*,860)
c    *  (jeff(iatom(j)),j=7*(JATOM/7)+1,JATOM)
c      do 52 i=1,nlist-1
c  52  write(*,'(I4,1P2E12.4,3I3,A6,0Pf8.2,8I4)')
c     *  i,IT(i),KT(i),NCH(i),NTOT(i),NEL(i),SPLIST(i),AWT(i),
c     *  (ZAT(j,i),NAT(j,i),j=1,NEL(i))
C================================================================
C== UPDATE MAIN ARRAYS                                         ==
C================================================================
c
c Make the initial estimate of the partial pressures for neutral atoms. These
c pressures are used as input to the linear solver. When only abundances are
c considered, the largest errors occur for low ionization elements, which can
c be highly ionized at low electron pressures. Thus, we apply a correction
c to recover the neutral fraction for each atom. The neutral fraction only
c corrects for losses into ionization states included in the species list.
c When the ionization correction is included, the largest error in the inital
c guess for carbon, which has unaccounted for losses into CO. Late in the
c convergence process, nitrogen becomes the dominant source of error.
c
      DO 5 J=1,JATOM
      P(J)=PG*ABUND(IATOM(J))/RNF(IATOM(J))
      ISPEC=INDSP(J)
      PP0(ISPEC)=P(J)
   5  CONTINUE
c
c Make an initial guess at the balance between H and H2.
c Assumes pressures of species other than H, H2, He, and Ne are negligible.
c Constraints:
c   KT(IIH2)*PP(IIH2)=P(1)**2           <-- chemical equilibrium
c   P(1)+2*PP(IIH2)=ABUND(1)*(PG-PE)    <-- H particle conservation
c
      IF(IIH2.GT.0) THEN
        PHyd=0.5*(-KT(IIH2)+SQRT(KT(IIH2)**2
     &        +4.0*KT(IIH2)*(PG-PE-P(2)-P(10))))
      ELSE
        PHyd=(PG-PE)*ABUND(1)
      ENDIF
c      IF(PHyd.GT.0.) P(1)=PHyd
c
c Make an initial guess at the balance between C, O, CO, and H2O.
c Constraints:
c   KT(IICO)*PP(IICO)=P(6)*P(8)         <-- chemical equilibrium
c   KT(IIH2O)*PP(IIH2O)=P(1)**2*P(8)    <-- chemical equilibrium
c   PTOTH=P(1)+2*PP(IIH2)       <-- defines density of H nuclei
c   PTOTC=P(6)+PP(IICO)                 <-- defines density of C nuclei
c   PTOTO=P(8)+PP(IICO)+PP(IIH2O)       <-- defines density of O nuclei
c   PTOTC=PTOTH*ABUND(6)/ABUND(1)       <-- abundance constraint
c   PTOTO=PTOTH*ABUND(8)/ABUND(1)       <-- abundance constraint
c
      PTOTH=P(1)
      IF(IIH2.GT.0) PTOTH=PTOTH+2.0*P(1)**2/KT(IIH2)
      PTOTC=PTOTH*ABUND(6)/ABUND(1)
      PTOTO=PTOTH*ABUND(8)/ABUND(1)
      IF(IIH2O.GT.0) THEN
        WATCOR=1.0+P(1)**2/KT(IIH2O)
        AQUAD=1.0/WATCOR
        IF(IICO.GT.0) THEN
          BQUAD=KT(IICO)+(PTOTO-PTOTC)/WATCOR
          CQUAD=-KT(IICO)*PTOTC
c          P(6)=(-BQUAD+SQRT(BQUAD**2-4.0*AQUAD*CQUAD))/(2.0*AQUAD)
c          P(8)=(P(6)+PTOTO-PTOTC)/WATCOR
        ELSE
c          P(6)=PTOTC
c          P(8)=PTOTO
        ENDIF
      ELSE
c        P(6)=PTOTC
c        P(8)=PTOTO
      ENDIF
c      IF(P(6).LE.0.) P(6)=PTOTC
c      IF(P(8).LE.0.) P(8)=PTOTO
      PE0=PE
      NAMEMX=BLANK
      DELMAX=0.0D0
c      COMPZ=0.0D0
c      PZS=0.0D0
c      write(*,*) SPLIST(1),P(1),SPLIST(IIH2),P(IIH2),
c     *           SPLIST(IIH2+1),P(IIH2+1),
c     *           SPLIST(IIH2+2),P(IIH2+2)
c      DO 6 J=1,JATOM
c      NN=INDSP(J)
c      IF(IPR(NN).NE.2) GOTO 3
c      NNP=INDX(3,ITAB(ZAT(1,NN)),1,1,1)
c      COMPZ=COMPZ+ABUND(IATOM(J))
c      IF(PE.EQ.0.0D0) PZS= PZS + P(J)
c      IF(PE.GT.0.0D0) PZS= PZS + (1.0D0+IT(NNP)/PE)*P(J)
c   6  CONTINUE
c      do J=1,JATOM
c        write(*,*) J,P(J),ABUND(IATOM(J)),SPLIST(INDSP(J))
c      enddo
c      write(*,*) JATOM+1,PE,'e-'
c      stop
C================================================================
C== MAIN LOOP: FILL LINEARIZED COEFFICIENT MATRIX AND RHS VECTOR,
C== AND SOLVE SYSTEM FOR PARTIAL PRESSURE CORRECTIONS.         ==
C== ISOLV=1: LINEARIZE ONLY THE PARTIAL PRESSURES OF THE NEUTRAL=
C== ATOMS FOR WHICH IPR(J)=1 (MAJOR SPECIES). THE ELECTRON     ==
C== PRESSURE PE IS ASSUMED TO BE GIVEN IN THIS CASE, AND SO IS ==
C== NOT INCLUDED IN THE LINEARIZATION. THIS IS NECESSARY SINCE ==
C== MOST OF THESE ELECTRONS (AT COOL TEMPS.) ORIGINATE FROM    ==
C== ELEMENTS NOT CONSIDERED IN THE LINEARIZATION. IN ORDER TO  ==
C== OBTAIN A GOOD VALUE FOR PE IN THE FIRST PLACE, IT IS       ==
C== NECESSARY TO CALL GAS WITH ISOLV=2.                        ==
C== ISOLV=2: THIS LINEARIZES THE PARTIAL PRESSURES OF THE NEUTRAL
C== ATOMS FOR WHICH IPR(J)=1 OR 2. THIS LIST OF ELEMENTS SHOULD==
C== INCLUDE ALL THE SIGNIFICANT CONTRIBUTORS TO THE TOTAL      ==
C== PRESSURE PG, AS WELL AS THE ELECTON PRESSURE PE. ANY ELEMENT=
C== (IPR(J)=3) NOT INCLUDED IS ASSUMED TO HAVE A NEGLIGIBLE    ==
C== EFFECT ON BOTH P AND PE.                                   ==
C== IN BOTH CASES, THE PARTIAL PRESSURES OF THE NEUTRAL ATOMS  ==
C== FOR ELEMENTS NOT INCLUDED IN THE LINEARIZATION ARE         ==
C== CALCULATED DIRECTLY FROM THE NOW DETERMINED PRESSURES OF   ==
C== THE LINEARIZED ELEMENTS.                                   ==
C================================================================
      NGIT=0
      RHSTOT=1.D99
C
C Top of loop in which linearized equations are solved recursively.
C
      REPEAT=0
      KMAX=1
   7  IF(NGIT.GE.MAXIT) THEN
        WRITE(*,208)
 208    FORMAT('*** ERROR: TOO MANY ITERATIONS IN ROUTINE "GAS"')
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),RHSTOT
        write(*,*) TEMP,PG,P(1),XNATOM,XNELEC
        STOP
      END IF
      NGIT=NGIT+1
      P(NEQ)=PE

      SCALE=10.D0
      IDIR=0
   9  CALL EOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)

c      write(*,*) 'Pe,SCALE,B(1),Pg=',PE,SCALE,B(1),PG,NGIT

      IF(B(1).GT.1.D2) THEN
        IF(IDIR.NE.-1) THEN
          SCALE=SQRT(SCALE)
          IDIR=-1
        ENDIF
C
C Neutral atomic pressures are too high. Scale them down until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)/SCALE
        ENDDO
        GOTO 9
      ELSE IF(B(1).LT.-1.D2) THEN
        IF(IDIR.NE.1) THEN
          SCALE=SQRT(SCALE)
          IDIR=1
        ENDIF
C
C Neutral atomic pressures are too low. Scale them up until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)*SCALE
        ENDDO
        GOTO 9
      ENDIF

      CALL EOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
C================================================================
C== NOW SOLVE THE LINEARIZED EQUATIONS (USING ROUTINE "LINEQ") ==
C================================================================
      IF(PRINT) THEN
        WRITE(*,200) NGIT
 200    FORMAT('LOG OF COEFFICIENT MATRIX AT ITERATION #',I5//)
        KK=MIN(30,NEQ-1)
        WRITE(*,201) (SPLIST(INDSP(K)),K=1,KK-1),'e-','RHS'
 201    FORMAT(4x,31(1x,a3,2x))
        DO 21 I=1,KK-1
        DO 20 J=1,KK-1
  20    AL(J)=LOG10(ABS(A(J,I))+1.0D-50)
        AL(KK)=LOG10(ABS(A(NEQ,I))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(I))+1.0D-50)
        NAMET=SPLIST(INDSP(I))
        WRITE(*,202) NAMET,(AL(J),J=1,KK+1)
  21    CONTINUE
        DO 22 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,NEQ))+1.0D-50)
  22    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,NEQ))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(NEQ))+1.0D-50)
        NAMET='e-'
        WRITE(*,202) NAMET,(AL(J),J=1,KK+1)
 202    FORMAT(A2,31F6.1)
        WRITE(*,'(/)')
c        stop
      END IF
C
C  Save a copy of the RHS for future step refinement
C
      DO 23 I=1,NEQ
  23  RHS(I)=B(I)
      RHSTOT=myDASUM(NEQ,RHS,1)
C
C  Solve linear system for corrections
C  In order not to solve for Pelect, one should use NEQ-1 as the first
C  argument. NEQ solves the whole system including electron pressure
C
c
c  Using LAPACK routine
c
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ
c        write(4) ((A(i,j),i=1,NEQ),j=1,NEQ)
c        write(4) (B(i),i=1,NEQ)
      CALL myDGESVX('E','N',NEQ,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *            RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *            WORK,IWORK,INFO)
c       write(4) (BB(i),i=1,NEQ)
c       stop
      CALL xDCOPY(NEQ,BB,1,B,1)
c      DO I=1,NEQ
c        B(I)=BB(I)
c      END DO
c
c  The same thing using LINEQ2 or LINEQ and BLAS 2/3
c          open(unit=4,file='dump.bin',form='UNFORMATTED')
c          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
c          close(4)
c      CALL LINEQ(NEQ,1,A,ELEDIM+1,IPIV,B,ELEDIM+1,INFO)
      IF(INFO.NE.0) THEN
        IF(REPEAT.LT.2) THEN
          DO J=1,NEQ-1
           P(J)=P(J)*0.999D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE IF(REPEAT.LT.4) THEN
          DO J=1,NEQ-1
           P(J)=P(J)*1.001D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE
          WRITE(*,*) 'EOS: LINEQ failed to solved for corrections to'
          WRITE(*,*) '     the partial pressures. Matrix is degenerate'
          WRITE(*,*) '     Temp=',TEMP,', Natom=',XNATOM,', Nelec=',
     *                XNELEC
          WRITE(*,*) '     INFO=',INFO,' Iter=',NGIT,' EQUED=',EQUED
cc          open(unit=4,file='dump.bin',form='UNFORMATTED')
cc          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
cc          close(4)
cc          write(1) 0
cc          close(1)
c          STOP
          CALL myDGESVX('E','N',NEQ-1,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,
     *                  EQUED,RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,
     *                  FERR,BERR,WORK,IWORK,INFO)
          CALL xDCOPY(NEQ,BB,1,B,1)
c          DO J=1,NEQ
c            B(J)=BB(J)
c          END DO
          PTOT=0.D0
          DO J=1,NEQ-1
            PTOT=PTOT+P(J)
          END DO
          PE=MAX(PG-PTOT,1.D-20)
        END IF
      END IF
      REPEAT=0

c
C=================================================================
C== FINALLY, UPDATE THE PARTIAL PRESSURES FOR THE MAJOR SPECIES ==
C== BY ADDING THE PRESSURE CORRECTIONS OBTAINED FOR EACH ATOM   ==
C== FROM THE LINEARIZATION PROCEDURE.                           ==
C=================================================================
      DELMAX=-1.0D0
      KMAX=1
      DO 31 K=1,NEQ
      ISPEC=INDSP(K)
C
C Compute the maximum correction in order to computer the under-relaxation factor
C
      DP=B(K)
      DELP=ABS(DP/MAX(P(K),1.D-50))
      IF(DELP.GT.DELMAX) THEN
        DELMAX=DELP
      END IF
  31  CONTINUE
C
C  Under-relaxation factor
C
      FACTOR=1.D0/(DELMAX+1.D0)
C
C Apply corrections
C
      DELMAX=-1.0D0
      KMAX=1
      DO 32 K=1,JATOM
      ISPEC=INDSP(K)
C
C  Restrict the correction to avoid getting negative pressures
C
      PNEW=P(K)-B(K)*FACTOR
      IF(PNEW.LT.0.D0) PNEW=MIN(MIN(P(K),ABS(PNEW)),PG)
c      IF(PNEW.LT.0.D0) PNEW=ABS(PNEW)
      DP=PNEW-P(K)
      IF(ABS(DP).GT.1.D-15) DP=DP*MIN(1.D0,0.4D0*P(K)/ABS(DP))
      P(K)=PNEW
      DELP=ABS(DP/MAX(P(K),1.D-50))
      IF(DELP.GT.DELMAX) THEN
        NAMEMX=SPLIST(ISPEC)
        DELMAX=DELP
        KMAX=K
      END IF
  32  CONTINUE

c      PENEW=BBB(NEQ)
      PENEW=PE-B(NEQ)*FACTOR
c      write(*,*) NEQ,PE,PENEW,B(NEQ),NGIT
      IF(PENEW.LT.0.D0) PENEW=MIN(PE,ABS(PENEW))
c      IF(PENEW.LT.0.D0) PENEW=ABS(PENEW)
      DPE=PENEW-PE
      IF(ABS(DPE).GT.1.D-15) DPE=DPE*MIN(1.D0,0.4D0*PE/ABS(DPE))
      PE=PENEW
      IF(ABS(PE/PG).GE.1.0D-15) THEN
        DELPE=ABS(DPE/PE)
        IF(DELPE.GT.DELMAX) NAMEMX=ENAME
        IF(DELPE.GT.DELMAX) DELMAX=DELPE
      END IF
C================================================================
C== PRINT OUT SUMMARY LINE FOR EACH ITERATION                  ==
C================================================================
      PTOT=PE
      PQ=0.0D0
c      write(*,*) 0,'e-',PE,PTOT,PG,NGIT
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=LOG(MAX(IT(ISPEC),1.D-115))-LOG(KT(ISPEC))-
     -     LOG(MAX(PE,1.D-115))*NQ
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+LOG(MAX(P(J),1.D-115))*NAT(I,ISPEC)
        ENDDO
c        PENQ=1.0D0
c        IF(PE.GT.0.0D0.AND.NQ.NE.0) PENQ=PE**NQ
c        PP(ISPEC)=IT(ISPEC)/(KT(ISPEC)*PENQ)*PF
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG
      ENDDO
c      stop
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PE-PQ)/PG
c      write(*,*) PG,PTOT,DELMAX,DPTOT,DPQ,FACTOR
      IF(PRINT) THEN
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),
     *               PTOT/TEMP/KBOL,DPTOT,PE/TEMP/KBOL,DPQ
 203    FORMAT(I10,2X,A8,1P9E11.3)
      END IF
      IF((DPTOT.GT.TOL.OR.DPQ.GT.TOL.OR.DELMAX.GT.TOL)
     *   .AND.NGIT.LT.MAXIT) GOTO 7
C
C Bottom of the loop in which linearized equations are solved recursively.
C
C================================================================
C== CALCULATE FINAL PARTIAL PRESSURES AFTER CONVERGENCE OBTAINED=
C================================================================
      PTOT=PE
      PD=0.0D0
      PU=0.0D0
      PQ=0.0D0
      DO 34 ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=1.0D0
        DO 33 I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF*P(J)**NAT(I,ISPEC)
  33    CONTINUE
        PENQ=1.0D0
        IF(PE.GT.0.0D0) PENQ=PE**NQ
        PP(ISPEC)=IT(ISPEC)/(KT(ISPEC)*PENQ)*PF
        PTOT=PTOT+PP(ISPEC)
        PD=PD+NTOT(ISPEC)*PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
        PU=PU+AWT(ISPEC)*PP(ISPEC)
  34  CONTINUE
      PP(NLIST)=PE
      PDTOT=PD+PE
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PQ-PE)/PG
      GMU=PU/PTOT
      ND=PTOT/(TEMP*KBOL)
      RHO=ND*GMU*HMASS
      XNE=PE/(TEMP*KBOL)
C================================================================
C== WRITE OUT FINAL PARTIAL PRESSURES                          ==
C================================================================
      IF(PRINT) THEN
c      IF(myDASUM(NLIST-1,PP,1)+PE.GT.PG*1.01D0) THEN
        write(*,'(''AFTER '',I3,'' iterations.   Max change of:'',G10.3,
     #      ''  in element:'',A)') NGIT,DELMAX,NAMEMX
        WRITE(*,'(''AFTER '',I3,'' ITERATIONS WITH ''/
     #            ''T='',1PE10.3,''   P='',E10.3)') NGIT,TEMP,
     #    myDASUM(NLIST-1,PP,1)+PE
        WRITE(*,'(''PDTOT='',1PE10.3,''   DPTOT='',E10.3,
     #            ''  DPQ='',E10.3,''  Nelectron='',E10.3,'' cm^3''/
     #    '' Nparticle='',1PE10.3,'' cm^3   Mean At.Wt.='',
     #    0PF7.3,''   Density='',1PE10.3,'' g/cm^3''//
     #    '' # Species   Abundance   Initial P   Final P'',
     #    ''      IT         KT         pf''/)')
     #    PDTOT,DPTOT,DPQ,XNE,ND-XNE,GMU,RHO
        NSP1=NLIST
        DO 35 ISPEC=1,NLIST-1
        IF(TYPE(ISPEC).NE.1) THEN
          WRITE(*,206) ISPEC,SPLIST(ISPEC),PP0(ISPEC),PP(ISPEC),
     #                 IT(ISPEC),KT(ISPEC),PART(ISPEC)
 206      FORMAT(I3,1X,A8,11X,1P5E11.3)
        ELSE
          J=IAT(ISPEC)
          WRITE(*,207) ISPEC,splist(ISPEC),ABUND(IATOM(J)),PP0(ISPEC),
     #                 PP(ISPEC),IT(ISPEC),KT(ISPEC),PART(ISPEC)
 207      FORMAT(I3,1X,A8,1P6E11.3)
        END IF
  35    CONTINUE
        WRITE(*,206) NSP1,ENAME,PE0,PE
        WRITE(*,*) JDAMAX(NLIST-1,PP,1),SPLIST(JDAMAX(NLIST-1,PP,1))
c        stop
      END IF
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
      PNOTE=0.D0
      DO 36 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        IF(PP(ISPEC)/KBOL/TEMP.GE.1.D-20) THEN
c          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP*PART(ISPEC))
          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP)
        ELSE
          XNPF(ISPEC)=0.0
        END IF
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        XNPF(ISPEC)=0.
        PFUNC(ISPEC)=1.
      END IF
      PNOTE=PNOTE+PP(ISPEC)
c      write(*,*) ISPEC,PNOTE,PP(ISPEC),SPLIST(ISPEC)
c      write(*,*) ISPEC,SPLIST(ISPEC),PFUNC(ISPEC)
  36  CONTINUE
c      write(*,*) 'e-',XNE
c      stop
      XNPF(NLIST)=XNE
      PFUNC(NLIST)=1.0
      XTOTAL=PD/(KBOL*TEMP)
      XNA=PNOTE/(KBOL*TEMP)
      Pgnew=PTOT
C
      RETURN
      END

C=========================================================================
C LOGARITHMIC version: the solution is found for the logs of ficticious
C                      partial pressures.
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE lnGAS(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *                 TOL,SPLIST,NLIST,XNE,XNA,RHO,Pgnew,
     *                 XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
     *                 FAILED)
c      SUBROUTINE lnGAS(TEMP,XNELEC,XNATOM,ABUND,ELEMEN,AMASS,ELESIZ,
c     *                 TOL,SPLIST,NLIST,
c     *                 XNE,XNA,RHO,XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
c     *                 FAILED)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      INTEGER MAXIT,MAXREF
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,MAXIT=10000,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0,MAXREF=10)

      LOGICAL PRINT,FAILED

      INTEGER NLIST,ELESIZ
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),XNPF(*),PFUNC(*),POTION(*),XTOTAL
      DOUBLE PRECISION FRACT(IONSIZ),IT(SPLSIZ-1),KT(SPLSIZ-1),
     *  AWT(SPLSIZ-1)

      DOUBLE PRECISION A(ELEDIM+1,ELEDIM+1),RHS(ELEDIM+1),
     *  AA(ELEDIM+1,ELEDIM+1),
     *  B(ELEDIM+1),BB(ELEDIM+1),
     *  P(ELEDIM+1),PP(SPLSIZ-1),PP0(SPLSIZ-1),PART(SPLSIZ-1),ND

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PARTN
c      DOUBLE PRECISION AT,BT,PN,DPF(4),CRATIO,BBB(ELEDIM+1),
c     *  PENQ,DPP,DPPE
      DOUBLE PRECISION RNF(ELEDIM),AL(ELEDIM+1)
      INTEGER NELM,NCHG,ANUM(4),NATM(4),IPIV(ELEDIM+1),IWORK(ELEDIM+1),
     *  INFO,ISPEC,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,IELM,NP,
     *  IIH2,IICO,IIH2O,NGIT,REPEAT
      DOUBLE PRECISION RATIOM,QPRD,RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,myDASUM,DELMAX,PE0,PTOTH,
     *  PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT,RENORM
c      DOUBLE PRECISION DUMMY,SCOLD,RHS0,RHS1,RHS2

c      DOUBLE PRECISION BOLD(ELEDIM+1),S(ELEDIM+1),GAMMA,BNORM,BOLDN
      DOUBLE PRECISION RSCL(ELEDIM+1),CSCL(ELEDIM+1)
c      DOUBLE PRECISION ROWCND,COLCND,AMX
      DOUBLE PRECISION FERR(1),BERR(1),WORK(5*(ELEDIM+1))
      CHARACTER*1 EQUED
      LOGICAL BARKLEM
      EXTERNAL myDASUM

      INTEGER NFIELDS
      PARAMETER (NFIELDS=40)
      CHARACTER*(*) FORMAT201,FORMAT202
c      CHARACTER*(*) AFIELDS
c      PARAMETER (AFIELDS=CHAR(NFIELDS/10+ICHAR('0'))//
c     *                   CHAR(MOD(NFIELDS,10)+ICHAR('0')))
c      PARAMETER (FORMAT201='(4x,'//AFIELDS//'(1X,A3,2X))')
c      PARAMETER (FORMAT202='(A2,'//AFIELDS//'F6.1)')
      PARAMETER (FORMAT201='(4x,48(1X,A3,2X))')
      PARAMETER (FORMAT202='(A2,48F6.1)')

cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real ttt(101)
c      real*8 Kttt(101)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C Initialize the Reciprocal Neutral Fraction (RNF). The RNF is used to
C adjust the initial neutral atomic partial pressures used in the linear
C solver. Originally, atomic species were assumed to be predominantly
C neutral, but at low electron pressures, this is a poor assumption for
C species with low ionization potentials.
C
      DO 1 I=1,ELEDIM
   1  RNF(I)=1.0D0
C
C Total gas and electron pressure
C
c      T=MAX(1200.,TEMP)
      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      if(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
c      PRINT=.TRUE.
      PRINT=.FALSE.
      PION=0
      IIH2=0
      IICO=0
      IIH2O=0
      JATOM=0
      NP=0
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      open(13,file='KT_eos.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')
c      write(13) NLIST,LEN(SPLIST(1))
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 4 ISPEC=1,NLIST-1
      PP0(ISPEC)=0.D0
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
c      write(*,*) ISPEC,'"'//SPLIST(ISPEC)//'"',NELM,NCHG,
c     *           ANUM,NATM,ELESIZ
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        KT(ISPEC)=1.0
        IT(ISPEC)=1.0
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(a,2i4)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        PART(ISPEC)=FRACT(1)
        POTION(ISPEC)=POTI(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,2)
          IT(ISPEC)=FRACT(NCHG+1)/FRACT(1)*PE**NCHG
          RNF(ANUM(1))=RNF(ANUM(1))+FRACT(NCHG+1)/FRACT(1)
c          if(ANUM(1).eq.26) write(*,*) SPLIST(ISPEC),NCHG,
c     *                      (FRACT(I),I=1,IONSIZ)
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PART(ISPEC)=FRACT(NCHG+1)
c          if(ANUM(1).eq.62) write(*,*) 'pf: ',SPLIST(ISPEC),NCHG,FRACT
          POTION(ISPEC)=POTI(NCHG+1)
          KT(ISPEC)=1.0
        ELSE IF(NCHG.LT.0) THEN
C
C Negative ions
C
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
          PARTN=FRACT(1)
          CALL  NEGION(ANUM(1),TEMP,PARTN,IT(ISPEC),
     *                 PART(ISPEC),POTION(ISPEC),BARKLEM)
        END IF
C
        KT(ISPEC)=1.D0
        AWT(ISPEC)=AMASS(ANUM(1))
        NTOT(ISPEC)=1
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
        IF(NCHG.LE.0) THEN
          QPRD=0.0D0
        ELSE
          QPRD=-NCHG*LOG10(2.0)
        ENDIF
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,3)
        IF(SPLIST(ISPEC).EQ.'H2')  IIH2=ISPEC
        IF(SPLIST(ISPEC).EQ.'CO')  IICO=ISPEC
        IF(SPLIST(ISPEC).EQ.'H2O') IIH2O=ISPEC
c       if(splist(ispec).eq.'N2')write(*,*)
c     *    anum(ielm),(fract(i),i=1,2)
   2    QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),TEMP,NTOT(ISPEC),RATIOM,QPRD,
     *              KT(ISPEC),PART(ISPEC),PION,BARKLEM)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        do ittt=0,100
c          ttt(ittt+1)=20.*ittt+1000.
c          CALL MOLCON(SPLIST(ISPEC),ttt(ittt+1),NTOT(ISPEC),
c     *                RATIOM,QPRD,Kttt(ittt+1),PART(ISPEC),PION)
c        END DO
c        write(13) SPLIST(ispec),ttt,Kttt
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  Finally, record the charge state of the molecule.
C
        IT(ISPEC)=1.D0
        IF(NCHG.GT.0.AND.BARKLEM) THEN
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
c-----------------------------------------------------------------------
c          IF(SPLIST(ISPEC).EQ.'H2+'.OR.SPLIST(ISPEC).EQ.'NO+') THEN
c            K=1
c            DO IELM=2,NELM
c              IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
c     *          K=IELM
c            ENDDO
c            IT(ISPEC)=IT(INDSP(ANUM(K))+1)
c            KT(ISPEC)=KT(ISPEC)/IT(ISPEC)
c          ENDIF
c          IT(ISPEC)=1.0
c-----------------------------------------------------------------------
C
C Positively charged molecules (single charge only!)
C
          K=1
          DO IELM=2,NELM
            IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
          ENDDO
          IT(ISPEC)=IT(INDSP(ANUM(K))+1)
        ELSE IF(NCHG.LT.0) THEN
C
C Negatively charged molecules (single charge only!)
C Known negatively charged molecules are:
C H2-, CH-, C2-, CN-, OH-, SiH-, HS-
C
          IF(SPLIST(ISPEC).EQ.'H2-') THEN
            PARTN=PART(INDSP(INDZAT( 1)))
            CALL NEGION( 1,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CH-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'C2-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'CN-') THEN
            PARTN=PART(INDSP(INDZAT( 6)))
            CALL NEGION( 6,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'OH-') THEN
            PARTN=PART(INDSP(INDZAT( 8)))
            CALL NEGION( 8,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'SiH-') THEN
            PARTN=PART(INDSP(INDZAT(14)))
            CALL NEGION(14,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE IF(SPLIST(ISPEC).EQ.'HS-') THEN
            PARTN=PART(INDSP(INDZAT(16)))
            CALL NEGION(16,TEMP,PARTN,IT(ISPEC),QPRD,POTI(1),BARKLEM)
          ELSE
            IT(ISPEC)=1.D0
          ENDIF
          IT(ISPEC)=1.D0
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
c      if(ANUM(IELM).eq.6.or.ANUM(IELM).eq.8) then
c        write(*,*) ISPEC,SPLIST(ISPEC),IT(ISPEC),KT(ISPEC)
c      endif
   3  NAT(IELM,ISPEC)=NATM(IELM)
C
C  Go back for next species.
C
c      write(*,*) ISPEC,SPLIST(ISPEC),IT(ISPEC),KT(ISPEC)
      IT(ISPEC)=MIN(MAX(1.D-250,IT(ISPEC)),1.D250)
      KT(ISPEC)=MIN(MAX(1.D-250,KT(ISPEC)),1.D250)
c      write(*,'(f10.2,I4,A12,4E13.4)') TEMP,ISPEC,SPLIST(ISPEC),
c     *     PART(ISPEC),KT(ISPEC),IT(ISPEC)
c     *     ,KT(ISPEC)/MAX(IT(ISPEC),1.D-150)
   4  CONTINUE
c      RENORM=LOG(SQRT(myDASUM(NLIST-1,KT,1)))
c      write(*,*) RENORM
c      DO ISPEC=1,NLIST-1
c        KT(ISPEC)=LOG(KT(ISPEC))+RENORM*NTOT(ISPEC)
c      END DO
      
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      close(13)
c      stop
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NEQ=JATOM+1
C==================================
C== End of species list parsing. ==
C==================================
C
C Print diagnostic: neutral fractions.
C
c     write(*,*) 'Reciprocal Neutral Fractions'
c     do 850 i=1,JATOM/7
c       write(*,860) (jeff(iatom(j)),j=7*i-6,7*i)
c850  continue
c860  format(1p,7e10.3,a)
c     if(JATOM.gt.7*(JATOM/7)) write(*,860)
c    *  (jeff(iatom(j)),j=7*(JATOM/7)+1,JATOM)
c      do 52 i=1,nlist-1
c  52  write(*,'(I4,1P2E12.4,3I3,A6,0Pf8.2,8I4)')
c     *  i,IT(i),KT(i),NCH(i),NTOT(i),NEL(i),SPLIST(i),AWT(i),
c     *  (ZAT(j,i),NAT(j,i),j=1,NEL(i))
C================================================================
C== UPDATE MAIN ARRAYS                                         ==
C================================================================
c
c Make the initial estimate of the partial pressures for neutral atoms. These
c pressures are used as input to the linear solver. When only abundances are
c considered, the largest errors occur for low ionization elements, which can
c be highly ionized at low electron pressures. Thus, we apply a correction
c to recover the neutral fraction for each atom. The neutral fraction only
c corrects for losses into ionization states included in the species list.
c When the ionization correction is included, the largest error in the inital
c guess for carbon, which has unaccounted for losses into CO. Late in the
c convergence process, nitrogen becomes the dominant source of error.
c
      DO 5 J=1,JATOM
      P(J)=PG*ABUND(IATOM(J))/RNF(IATOM(J))
      ISPEC=INDSP(J)
      PP0(ISPEC)=P(J)
   5  CONTINUE
c
c Make an initial guess at the balance between H and H2.
c Assumes pressures of species other than H, H2, He, and Ne are negligible.
c Constraints:
c   KT(IIH2)*PP(IIH2)=P(1)**2           <-- chemical equilibrium
c   P(1)+2*PP(IIH2)=ABUND(1)*(PG-PE)    <-- H particle conservation
c
      IF(IIH2.GT.0) THEN
        PHyd=0.5*(-KT(IIH2)+SQRT(KT(IIH2)**2
     *        +4.0*KT(IIH2)*(PG-PE-P(2)-P(10))))
      ELSE
        PHyd=(PG-PE)*ABUND(1)
      END IF
c      IF(PHyd.GT.0.0.AND.PHyd.LT.Pgas-Pelec) P(1)=PHyd
c
c Make an initial guess at the balance between C, O, CO, and H2O.
c Constraints:
c   KT(IICO)*PP(IICO)=P(6)*P(8)         <-- chemical equilibrium
c   KT(IIH2O)*PP(IIH2O)=P(1)**2*P(8)    <-- chemical equilibrium
c   PTOTH=P(1)+2*PP(IIH2)               <-- defines density of H nuclei
c   PTOTC=P(6)+PP(IICO)                 <-- defines density of C nuclei
c   PTOTO=P(8)+PP(IICO)+PP(IIH2O)       <-- defines density of O nuclei
c   PTOTC=PTOTH*ABUND(6)/ABUND(1)       <-- abundance constraint
c   PTOTO=PTOTH*ABUND(8)/ABUND(1)       <-- abundance constraint
c
      PTOTH=P(1)
      IF(IIH2.GT.0) PTOTH=PTOTH+2.0*P(1)**2/KT(IIH2)
      PTOTC=PTOTH*ABUND(6)/ABUND(1)
      PTOTO=PTOTH*ABUND(8)/ABUND(1)
      IF(IIH2O.GT.0) THEN
        WATCOR=1.0+P(1)**2/KT(IIH2O)
        AQUAD=1.0/WATCOR
        IF(IICO.GT.0) THEN
          BQUAD=KT(IICO)+(PTOTO-PTOTC)/WATCOR
          CQUAD=-KT(IICO)*PTOTC
c          P(6)=(-BQUAD+SQRT(BQUAD**2-4.0*AQUAD*CQUAD))/(2.0*AQUAD)
c          P(8)=(P(6)+PTOTO-PTOTC)/WATCOR
        ELSE
c          P(6)=PTOTC
c          P(8)=PTOTO
        END IF
      ELSE
c        P(6)=PTOTC
c        P(8)=PTOTO
      END IF
c      IF(P(6).LE.0.0.OR.P(6).GT.0.1*P(1)) P(6)=PTOTC
c      IF(P(8).LE.0.0.OR.P(8).GT.0.1*P(1)) P(8)=PTOTO
      PE0=PE
      NAMEMX=BLANK
      DELMAX=0.0D0
c      COMPZ=0.0D0
c      PZS=0.0D0
c      DO 6 J=1,JATOM
c      NN=INDSP(J)
c      IF(IPR(NN).NE.2) GOTO 3
c      NNP=INDX(3,ITAB(ZAT(1,NN)),1,1,1)
c      COMPZ=COMPZ+ABUND(IATOM(J))
c      IF(PE.EQ.0.0D0) PZS= PZS + P(J)
c      IF(PE.GT.0.0D0) PZS= PZS + (1.0D0+IT(NNP)/PE)*P(J)
c   6  CONTINUE
c      do J=1,JATOM
c        write(*,*) J,P(J),ABUND(IATOM(J)),SPLIST(INDSP(J))
c      END DO
c      write(*,*) JATOM+1,PE,'e-'
c      stop
C================================================================
C== MAIN LOOP: FILL LINEARIZED COEFFICIENT MATRIX AND RHS VECTOR,
C== AND SOLVE SYSTEM FOR PARTIAL PRESSURE CORRECTIONS.         ==
C== ISOLV=1: LINEARIZE ONLY THE PARTIAL PRESSURES OF THE NEUTRAL=
C== ATOMS FOR WHICH IPR(J)=1 (MAJOR SPECIES). THE ELECTRON     ==
C== PRESSURE PE IS ASSUMED TO BE GIVEN IN THIS CASE, AND SO IS ==
C== NOT INCLUDED IN THE LINEARIZATION. THIS IS NECESSARY SINCE ==
C== MOST OF THESE ELECTRONS (AT COOL TEMPS.) ORIGINATE FROM    ==
C== ELEMENTS NOT CONSIDERED IN THE LINEARIZATION. IN ORDER TO  ==
C== OBTAIN A GOOD VALUE FOR PE IN THE FIRST PLACE, IT IS       ==
C== NECESSARY TO CALL GAS WITH ISOLV=2.                        ==
C== ISOLV=2: THIS LINEARIZES THE PARTIAL PRESSURES OF THE NEUTRAL
C== ATOMS FOR WHICH IPR(J)=1 OR 2. THIS LIST OF ELEMENTS SHOULD==
C== INCLUDE ALL THE SIGNIFICANT CONTRIBUTORS TO THE TOTAL      ==
C== PRESSURE PG, AS WELL AS THE ELECTON PRESSURE PE. ANY ELEMENT=
C== (IPR(J)=3) NOT INCLUDED IS ASSUMED TO HAVE A NEGLIGIBLE    ==
C== EFFECT ON BOTH P AND PE.                                   ==
C== IN BOTH CASES, THE PARTIAL PRESSURES OF THE NEUTRAL ATOMS  ==
C== FOR ELEMENTS NOT INCLUDED IN THE LINEARIZATION ARE         ==
C== CALCULATED DIRECTLY FROM THE NOW DETERMINED PRESSURES OF   ==
C== THE LINEARIZED ELEMENTS.                                   ==
C================================================================
      FACTOR=1.D0
      NGIT=0
      RHSTOT=1.D99
      goto 2222
C
C Top of loop in which linearized equations are solved recursively.
C
      KMAX=1
c      PG=PG+myDASUM(NEQ-1,P)*(RENORM-1)
      DO J=1,NEQ-1
c        P(J)=LOG(P(J))+RENORM
        P(J)=LOG(P(J))
      END DO
      PE=LOG(MAX(PE,1.D-150))
c      open(unit=4,file='dump.bin',form='UNFORMATTED')
c      write(4) NEQ
      REPEAT=0
   7  IF(NGIT.GE.MAXIT) THEN
        WRITE(*,208)
 208    FORMAT('*** ERROR: TOO MANY ITERATIONS IN ROUTINE "GAS"')
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),RHSTOT
        write(*,*) TEMP,PG,P(1),XNATOM,XNELEC
        STOP
      END IF
      NGIT=NGIT+1
      P(NEQ)=PE

c      do J=1,NEQ
c        p(J)=exp(p(j))
c      enddo
c      write(*,*) (P(J),J=1,NEQ)
c      CALL lnEOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
c     *     IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c      CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
c     *     IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c      do j=1,NEQ
c        SCALE=P(J)
c        P(J)=P(J)+0.1d0
c        CALL lnEOSFCN(NEQ,P,BB,A,1,PG,NCH,NLIST,
c     *    IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c        write(*,*) J,SCALE
c        write(*,'(40e10.3)')(a(i,j)-(bb(i)-b(i))/0.1d0
c     *                            ,i=1,40)
c        write(*,'(40e10.3)')(a(i,j),i=1,40)
c        write(*,'(40e10.3)')((bb(i)-b(i))/0.1d0,i=1,40)
c        write(*,'(40e10.3)')(bb(i),i=1,40)
c        P(J)=SCALE
c      enddo
c      stop

      SCALE=100.D0
      IDIR=0
c      do j=1,NEQ
c        write(*,*) J,P(J),PG
c      enddo
c      write(*,*) B(1),PG
   9  CALL lnEOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
     *              IATOM,INDSP,NAT,ZAT,NTOT,
     *              NEL,IAT,INDZAT,ABUND,KT,IT)
c      write(*,*) SCALE,B(1),PG
      IF(B(1).GT.0.001D0*PG) THEN
        IF(IDIR.NE.-1) THEN
          SCALE=SQRT(SCALE)
          IDIR=-1
        END IF
C
C Neutral atomic pressures are too high. Scale them down until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)-LOG(SCALE)
        END DO
        GOTO 9
      ELSE IF(B(1).LT.-0.001D0*PG) THEN
        IF(IDIR.NE.1) THEN
          SCALE=SQRT(SCALE)
          IDIR=1
        END IF
C
C Neutral atomic pressures are too low. Scale them up until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)+LOG(SCALE)
        END DO
        GOTO 9
      END IF

c      IF(B(1).GT.0.02D0*PG) THEN
c        IF(IDIR.NE.1) THEN
c          SCALE=SQRT(SCALE)
c          IDIR=1
c        END IF
cC
cC Neutral atomic pressures are too high. Scale them down until
cC partical conservation equation will become negative
cC
c        DO ISPEC=1,NLIST-1
c          J=0
c          DO I=1,NEL(ISPEC)
c            J=J+NAT(I,ISPEC)
c          END DO
c          write(*,*) ISPEC,SPLIST(ISPEC),J,NCH(ISPEC)
c          KT(ISPEC)=KT(ISPEC)*SCALE**J
c          IT(ISPEC)=IT(ISPEC)*SCALE**NCH(ISPEC)
c        END DO
c        GOTO 9
c      ELSE IF(B(1).LT.-0.02D0*PG) THEN
c        IF(IDIR.NE.-1) THEN
c          SCALE=SQRT(SCALE)
c          IDIR=-1
c        END IF
cC
cC Neutral atomic pressures are too low. Scale them up until
cC partical conservation equation will become negative
cC
c        DO ISPEC=1,NLIST-1
c          J=0
c          DO I=1,NEL(ISPEC)
c            J=J+NAT(I,ISPEC)
c          END DO
c          KT(ISPEC)=KT(ISPEC)/SCALE**J
c          IT(ISPEC)=IT(ISPEC)/SCALE**NCH(ISPEC)
c        END DO
c        GOTO 9
c      END IF

c      do j=1,NEQ
c        write(*,*) J,P(J),PG
c      enddo
c      write(*,*) B(1),PG
      CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
     *              IATOM,INDSP,NAT,ZAT,NTOT,
     *              NEL,IAT,INDZAT,ABUND,KT,IT)
c      DO I=1,NEQ-1
c      WRITE(*,FORMAT202) SPLIST(INDSP(I)),(A(I,J),J=1,NEQ-1),B(I)
c      END DO
c      stop
C
C================================================================
C== NOW SOLVE THE LINEARIZED EQUATIONS (USING ROUTINE "LINEQ") ==
C================================================================
      IF(PRINT) THEN
        WRITE(*,200) NGIT
 200    FORMAT('LOG OF COEFFICIENT MATRIX AT ITERATION #',I5/)
        KK=MIN(NFIELDS,NEQ-1)
        WRITE(*,FORMAT201) (SPLIST(INDSP(K)),K=1,KK-1),'e-','RHS'
        DO 21 I=1,KK-1
        DO 20 J=1,KK-1
  20    AL(J)=LOG10(ABS(A(J,I))+1.0D-50)
        AL(KK)=LOG10(ABS(A(NEQ,I))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(I))+1.0D-50)
        NAMET=SPLIST(INDSP(I))
        WRITE(*,FORMAT202) NAMET,(AL(J),J=1,KK+1)
  21    CONTINUE
        DO 22 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,NEQ))+1.0D-50)
  22    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,NEQ))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(NEQ))+1.0D-50)
        NAMET='e-'
        WRITE(*,FORMAT202) NAMET,(AL(J),J=1,KK+1)
        WRITE(*,'(/)')
      END IF
c      stop
C
C  Save a copy of the RHS for future step refinement
C
      DO 23 I=1,NEQ
  23  RHS(I)=B(I)
      RHSTOT=myDASUM(NEQ,RHS,1)
C
C  Solve linear system for corrections
C  In order not to solve for Pelect, one should use NEQ-1 as the first
C  argument. NEQ solves the whole system including electron pressure
C
c
c  Using LAPACK routine
c
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ
c        write(4) ((A(i,j),i=1,NEQ),j=1,NEQ)
c        write(4) (B(i),i=1,NEQ)
c      write(4) ((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
      CALL myDGESVX('E','N',NEQ,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *            RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *            WORK,IWORK,INFO)
c      stop
      CALL xDCOPY(NEQ,BB,1,B,1)
c      DO I=1,NEQ
c        B(I)=BB(I)
c      ENDDO
c      write(4) ((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
c
c  The same thing using LINEQ2 or LINEQ and BLAS 2/3
c      CALL LINEQ(NEQ,1,A,ELEDIM+1,IPIV,B,ELEDIM+1,INFO)
      IF(INFO.NE.0) THEN
        IF(REPEAT.LT.2) THEN
          DO J=1,NEQ-1
           P(J)=P(J)-0.01D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE IF(REPEAT.LT.4) THEN
          DO J=1,NEQ-1
           P(J)=P(J)+0.01D0
          END DO
          REPEAT=REPEAT+1
          GO TO 7
        ELSE
          WRITE(*,*) 'lnGAS: DGESVX failed to solved for corrections to'
          WRITE(*,*) '  the partial pressures. Matrix is degenerate'
          WRITE(*,*) '  Temp=',TEMP,', Natom=',XNATOM,', Nelec=',XNELEC
          IF(INFO.EQ.NEQ) THEN
            WRITE(*,*) '  Pg=',PG,', INFO=',INFO,
     *                 ', Element: e-',
     *                 ', Iter=',NGIT,' EQUED=',EQUED
          ELSE
            WRITE(*,*) '  Pg=',PG,', INFO=',INFO,
     *                 ', Element: ',SPLIST(INDSP(INFO)),
     *                 ', Iter=',NGIT,' EQUED=',EQUED
          END IF
          CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,IATOM,INDSP,
     *                  NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
          open(unit=4,file='dump.bin',form='UNFORMATTED')
          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
          close(4)
          WRITE(*,*) '  Matrix and the RHS were dumped to file dump.bin'
          STOP
c          CALL myDGESVX('E','N',NEQ-1,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
c     *                RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
c     *                WORK,IWORK,INFO)
c          CALL xDCOPY(NEQ-1,BB,1,B,1)
cc          DO I=1,NEQ
cc            B(I)=BB(I)
cc          END DO
c          PTOT=0.D0
c          DO J=1,NEQ-1
c            PTOT=PTOT+exp(P(J)-B(J))
c          END DO
c          PE=MAX(PG-PTOT,1.D-20)
c          Pe=log(Pe)
        END IF
      END IF
      REPEAT=0
c      IF(INFO.NE.0) THEN
c        WRITE(*,*) 'lnEOS: LINEQ failed to solved for corrections to'
c        WRITE(*,*) '     the partial pressures. Matrix is degenerate'
c        WRITE(*,*) '     Temp=',TEMP,', Natom=',XNATOM,', Nelec=',XNELEC
c        WRITE(*,*) '     Pg=',PG,', INFO=',INFO,
c     *             ', Element: ',SPLIST(INDSP(INFO)),
c     *             ', Iter=',NGIT,' EQUED=',EQUED
cc        open(unit=4,file='dump.bin',form='UNFORMATTED')
cc        write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
cc        close(4)
cc        write(1) 0
cc        close(1)
c        IF(PRINT) THEN
cc          close(4)
c          STOP
c        END IF
cc        DO J=1,NEQ
cc          P(J)=MAX(P(J)+0.1D0,-115.d0)
cc          write(*,*) J,P(J),B(J),B(J)*FACTOR
cc        END DO
c        write(*,*) P(INFO),B(INFO),B(INFO)*FACTOR
c        P(INFO)=MAX(P(INFO)+0.1D0,-115.d0)
c        PRINT=.TRUE.
c        GO TO 9
c      END IF
c
C=================================================================
C== FINALLY, UPDATE THE PARTIAL PRESSURES FOR THE MAJOR SPECIES ==
C== BY ADDING THE PRESSURE CORRECTIONS OBTAINED FOR EACH ATOM   ==
C== FROM THE LINEARIZATION PROCEDURE.                           ==
C=================================================================
      DELMAX=-200.0D0
      KMAX=1
      DO K=1,JATOM
c        write(*,*) K,P(K),B(K)
        ISPEC=INDSP(K)
c        DP=ABS(P(K))
        DELP=ABS(B(K))
c        IF(DP.GT.1.D-10) DELP=DELP/DP
        IF(DELP.GT.DELMAX) THEN
          NAMEMX=SPLIST(ISPEC)
          DELMAX=DELP
          KMAX=K
        END IF
      END DO
c      DPE=ABS(P(NEQ))
      DELPE=ABS(B(NEQ))
c      IF(DPE.GT.1.D-10) DELPE=DELPE/DPE
      IF(DELPE.GT.DELMAX) THEN
        NAMEMX=ENAME
        DELMAX=DELPE
        KMAX=NEQ
      END IF
c      write(*,*) KMAX,EXP(P(KMAX)),EXP(B(KMAX)),P(KMAX),B(KMAX)
C
C  Under-relaxation factor
C
      FACTOR=0.2D0/(DELMAX+0.2D0)
      DO K=1,JATOM
C
C  Apply corrections
C
        DP=B(K)*FACTOR
c        DP=10.D0*DP/MAX(10.D0,ABS(DP))
        PNEW=P(K)-DP
        P(K)=MAX(PNEW,-115.D0)
      END DO
      DP=B(NEQ)*FACTOR
c      DP=10.D0*DP/MAX(10.D0,ABS(DP))
      PENEW=PE-DP
      PE=MAX(PENEW,-115.D0)
C================================================================
C== PRINT OUT SUMMARY LINE FOR EACH ITERATION                  ==
C================================================================
      PTOT=EXP(PE)
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+P(J)*NAT(I,ISPEC)
        END DO
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG,NQ,PQ,EXP(PE)
      END DO
c      stop
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(EXP(PE)-PQ)/PG
c      write(*,*) DELMAX,DPTOT,DPQ
      IF(PRINT) THEN
        WRITE(*,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),
     *               PTOT/TEMP/KBOL,DPTOT,EXP(PE)/TEMP/KBOL,DPQ,FACTOR
 203    FORMAT(I10,2X,A8,1P,9E11.3)
      END IF
c      write(*,*) NGIT,TOL,DPTOT,DELMAX,PTOT,PG
      IF((RHSTOT.GT.TOL.OR.DPTOT.GT.TOL.OR.DELMAX.GT.TOL)
     *   .AND.NGIT.LT.MAXIT) GO TO 7
C
C Bottom of the loop in which linearized equations are solved recursively.
C
C================================================================
C== CALCULATE FINAL PARTIAL PRESSURES AFTER CONVERGENCE OBTAINED=
C================================================================
c      write(*,*) RHSTOT,DELMAX,DPTOT,DPQ,TOL
      PTOT=EXP(PE)
      PD=0.0D0
      PU=0.0D0
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+P(J)*NAT(I,ISPEC)
        END DO
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PD=PD+NTOT(ISPEC)*PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
        PU=PU+AWT(ISPEC)*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG,NQ,PQ,EXP(PE)
      END DO
      PE=EXP(PE)
      DO J=1,JATOM
        P(J)=EXP(P(J))
      END DO
      PP(NLIST)=PE
      PDTOT=PD+PE
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PQ-PE)/PG
      GMU=PU/PTOT
      ND=PTOT/(TEMP*KBOL)
      RHO=ND*GMU*HMASS
      XNE=PE/(TEMP*KBOL)
C================================================================
C== WRITE OUT FINAL PARTIAL PRESSURES                          ==
C================================================================
      IF(PRINT) THEN
        write(*,'(''AFTER '',I3,'' iterations.   Max change of:'',G10.3,
     #      ''  in element:'',A)') NGIT,DELMAX,NAMEMX
        WRITE(*,'(''AFTER '',I3,'' ITERATIONS WITH ''/
     #            ''T='',1PE10.3,''   P='',E10.3)') NGIT,TEMP,PG
        WRITE(*,'(''PDTOT='',1PE10.3,''   DPTOT='',E10.3,
     #            ''  DPQ='',E10.3,''  Nelectron='',E10.3,'' cm^3''/
     #    '' Nparticle='',1PE10.3,'' cm^3   Mean At.Wt.='',
     #    0PF7.3,''   Density='',1PE10.3,'' g/cm^3''/
     #    '' # Species   Abundance   Initial P   Final P'',
     #    ''      IT         KT         pf''//)')
     #   PDTOT,DPTOT,DPQ,XNE,ND-XNE,GMU,RHO
        NSP1=NLIST
        DO 35 ISPEC=1,NLIST-1
        IF(TYPE(ISPEC).NE.1) THEN
          WRITE(*,206) ISPEC,SPLIST(ISPEC),PP0(ISPEC),PP(ISPEC),
     #                 IT(ISPEC),KT(ISPEC),PART(ISPEC)
 206      FORMAT(I3,1X,A8,11X,1P,5E11.3)
        ELSE
          J=IAT(ISPEC)
          WRITE(*,207) ISPEC,splist(ISPEC),ABUND(IATOM(J)),PP0(ISPEC),
     #                 PP(ISPEC),IT(ISPEC),KT(ISPEC),PART(ISPEC)
 207      FORMAT(I3,1X,A8,1P,6E11.3)
        END IF
  35    CONTINUE
        WRITE(*,206) NSP1,ENAME,PE0,EXP(PE)
      END IF
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
2222  continue
      PNOTE=0.0
      DO 36 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        IF(PP(ISPEC)/KBOL/TEMP.GE.1.D-20) THEN
c          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP*PART(ISPEC))
          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP)
        ELSE
          XNPF(ISPEC)=0.0
        END IF
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        XNPF(ISPEC)=0.
        PFUNC(ISPEC)=1.
      END IF
      PNOTE=PNOTE+PP(ISPEC)
c      write(*,*) ISPEC,PNOTE,PP(ISPEC),SPLIST(ISPEC)
c      write(*,*) ISPEC,SPLIST(ISPEC),PFUNC(ISPEC)
  36  CONTINUE
      XNPF(NLIST)=XNE
      PFUNC(NLIST)=1.0
      XTOTAL=PD/(TEMP*KBOL)
      XNA=PNOTE/(TEMP*KBOL)
c      write(*,*) 'Pg,PD,PNOTE,PE,PNOTE+PE',Pg,PD,PTOT,PE,PNOTE+PE
      Pgnew=Ptot
C
      RETURN
      END


C=========================================================================
C MOLCON: Returns equilibrium constant and partition function for a given
C   molecule and temperature.
C
C Inputs:
C   SPNAME [character(*)] Name of molecule, chosen from SPLIST below.
C   T [real] Temperature (in K) at which EQK and PART are to be found.
C   NTOT [real] Total number of atoms in the molecule.
C   RATIOM [real] Logarithm (base 10) of mass ratio (in g^(natoms-1)):
C     ratiom = Sum{log10(Atomic Masses)} - log10(Sum{Atomic Masses})
C   QPRD [double] Logarithm of product of atomic partition functions:
C     qprd = Sum{log10(Atomic Partition Functions)}
C
C Outputs:
C   EQK [real] Equilibrium constant (in dynes/cm/cm) at temperature T,
C     calculated from dissociation energy and partition function.
C   PART [real] Partition function at temperature T, calculated from
C     expressions in the references cited below.
C
C References:
C   For diatomic molecules: Sauval & Tatum (1984, ApJS, 56, 193).
C
      SUBROUTINE MOLCON(SPNAME,T,NTOT,RATIOM,QPRD,EQK,PART,PION,
     *                  BARKLEM)
C
      INCLUDE 'SIZES.EOS'
C
      INTEGER MSPEC,NTOT
      DOUBLE PRECISION KERG,KEV
      DOUBLE PRECISION RATIOM,QPRD,PION
      PARAMETER (KERG=1.38065D-16,KEV=KERG/1.60219D-12)
      PARAMETER (CONST=25947.256)
C
      REAL T
      DOUBLE PRECISION TLIM,TH,LOGTH,EQK,PART,Qm_spln,Kp_spln
c      DOUBLE PRECISION EQK_ST
      LOGICAL BARKLEM
C
C Combine equilibrium constant coefficients into one large array.
C
      PARAMETER (MSPEC=197)
      PARAMETER (NEQCOE=7)
      DOUBLE PRECISION COEF(NEQCOE,MSPEC)
      DOUBLE PRECISION C01(NEQCOE,50),C02(NEQCOE,50),
     *                 C03(NEQCOE,50),C04(NEQCOE,47)
      EQUIVALENCE (C01(1,1),COEF(1,  1)),(C02(1,1),COEF(1, 51))
      EQUIVALENCE (C03(1,1),COEF(1,101)),(C04(1,1),COEF(1,151))
C
C Combine partition function coefficients into one large array.
C
      PARAMETER (NPCOEF=11)
      DOUBLE PRECISION PCOEF(NPCOEF,MSPEC)
      DOUBLE PRECISION P01(NPCOEF,50),P02(NPCOEF,50),
     *                 P03(NPCOEF,50),P04(NPCOEF,47)
      EQUIVALENCE (P01(1,1),PCOEF(1,  1)),(P02(1,1),PCOEF(1, 51))
      EQUIVALENCE (P03(1,1),PCOEF(1,101)),(P04(1,1),PCOEF(1,151))
C
      CHARACTER SPNAME*(*),SPLIST(MSPEC)*(SPCHAR)
      SAVE
C
C Molecular species list from NextGen models (Allard & Hauschildt).
C See old/eos.4.f for molecular species list from Sauval & Tatum (1984).
C
      DATA SPLIST/
     * 'H2      ','CO      ','H2O     ','OH      ','N2      ',
     * 'SiO     ','HS      ','H2S     ','NH      ','SiH     ',
     * 'CH      ','H2+     ','NO      ','MgH     ','HCl     ',
     * 'SiS     ','AlOH    ','NH2     ','AlH     ','CN      ',
     * 'CO2     ','SO      ','TiO     ','S2      ','FeH     ',
     * 'NH3     ','HCN     ','HCO     ','O2      ','CH2     ',
     * 'HF      ','H3+     ','CaH     ','Al2O    ','AlO     ',
     * 'CH3     ','SiH2    ','MgO     ','C2      ','TiO2    ',
     * 'VO2     ','NaH     ','AlCl    ','AlF     ','VO      ',
     * 'CS      ','MgOH    ','PO2     ','CaOH    ','PH2     ',
     * 'C2H     ','ScO     ','AlO2H   ','AlS     ','FeO     ',
     * 'CrO     ','CH4     ','NS      ','SO2     ','SiN     ',
     * 'OH-     ','ZrO     ','NO+     ','ZrO2    ','BO      ',
     * 'SiO2    ','HBO     ','SiC     ','YO2     ','TiS     ',
     * 'HBO2    ','C2H2    ','OCS     ','ZrO+    ','NaOH    ',
     * 'CaCl    ','AlOF    ','YO      ','NaCl    ','C2O     ',
     * 'CHP     ','HS-     ','H2-     ','TiH     ','PH3     ',
     * 'MgS     ','TiO+    ','LaO2    ','Si2     ','SiH4    ',
     * 'BH2     ','AlOCl   ','LaO     ','C2N     ','AlBO2   ',
     * 'KCl     ','SiH-    ','CaF     ','CaO2H2  ','KOH     ',
     * 'CN-     ','Al2O2   ','BaOH    ','SrOH    ','BO2     ',
     * 'SiF     ','CH-     ','C3      ','C2-     ','MgO2H2  ',
     * 'BeOH    ','HBS     ','SiC2    ','FeO2H2  ','CrO2    ',
     * 'BeH2O2  ','BH3     ','NaCN    ','BeH2    ','Si2N    ',
     * 'CaCl2   ','NaBO2   ','C3H     ','OBF     ','CS2     ',
     * 'LiOH    ','Al2     ','LiCl    ','TiOCl   ','C2H4    ',
     * 'CHCl    ','TiCl    ','AlOF2   ','KBO2    ','Si2C    ',
     * 'CHF     ','BO-     ','AlO2    ','BaO2H2  ','OTiF    ',
     * 'CS-     ','C2N2    ','SrO2H2  ','ClCN    ','AlClF   ',
     * 'KCN     ','AlCl2   ','BaCl2   ','AlF2    ','MgCl2   ',
     * 'FeO-    ','BO2H2   ','SiH3Cl  ','FeCl2   ','Si3     ',
     * 'SiH3F   ','CH3Cl   ','SrCl2   ','CaF2    ','TiF2    ',
     * 'LiBO2   ','MgClF   ','BeBO2   ','C2HCl   ','TiCl2   ',
     * 'C4      ','H3BO3   ','MgF2    ','BaClF   ','BeF2    ',
     * 'C2HF    ','BeCl2   ','TiOCl2  ','ZrCl2   ','BaF2    ',
     * 'BeC2    ','Be2O    ','SrF2    ','ZrF2    ','FeF2    ',
     * 'P4      ','SiH2F2  ','H3O+    ','C5      ','TiF3    ',
     * 'TiCl3   ','ZrCl3   ','Na2Cl2  ','Na2O2H2 ','Be3O3   ',
     * 'K2Cl2   ','K2O2H2  ','ZrCl4   ','Na2C2N2 ','ZrF4    ',
     * 'Li2O2H2 ','CrH     '/
C
C Dissociation energy (first column, in eV) and equilibrium constant
C   coefficients. See the file "atomiz.notes" for the information on the
C   origin of the dissociation energies. The polynomial fit coefficients
C   for the equilibrium constants were determined with "ng_kfit.pro" and
C   are meant to reproduce the constants used in constructing the NextGen
C   models. The NextGen equilibrium constants were fit over the temperature
C   range 1600 < T < 7730 K. The fits are likely to diverge rapidly from
C   the truth outside this temperature range.
C Equilibrium constants may be constructed from the coefficients using:
C
C     log10(Kp) = Sum{i=2,7}{COEF(i)*log10(THETA)**(i-2)} - COEF(1)*THETA
C
      DATA C01/
     *   4.4781, 12.1354, -0.7752, -0.7821,  0.1464,  0.1603, -0.0626,  H2
     *  11.0920, 13.2368, -0.8342, -0.0477, -0.2923, -0.4557,  0.6108,  CO
     *   9.6221, 24.7774, -2.3428,  1.6868, -1.2845, -2.9925,  3.6555,  H2O
     *   4.3920, 11.8016, -0.8507, -0.5193,  0.0502, -0.3409,  0.4836,  OH
     *   9.7594, 12.8868, -0.8813,  0.2639, -1.5912,  1.5866, -0.5407,  N2
     *   8.2600, 12.9252, -0.7608, -0.3541,  1.5620, -3.5952,  2.5962,  SiO
     *   3.5500, 11.4382, -0.7816, -0.4659,  0.4314, -1.2144,  0.9648,  HS
     *   7.5946, 23.8543, -0.9525, -0.8118,  0.2051, -1.0299,  1.1555,  H2S
     *   3.4700, 11.4658, -0.7258, -0.6418, -0.0442,  0.2836, -0.1618,  NH
     *   3.0600, 11.2595, -0.6962, -0.6435,  0.6663, -0.3357, -0.4151,  SiH
     *   3.4650, 11.5333, -0.5255, -0.7105,  0.2264, -0.9271,  0.9577,  CH
     *   2.6508, 15.8052, 33.7578, 34.5956, 27.3455, 16.6214,  9.9717,  H2+
     *   6.4968, 11.9347, -0.7596,  0.0953, -0.9731,  0.8265, -0.2151,  NO
     *   1.3400, 10.2911, -0.3698, -0.0655, -2.9771,  6.1325, -4.3869,  MgH
     *   4.4336, 11.9041, -0.8281, -0.6163,  0.1580, -0.5068,  0.5164,  HCl
     *   6.4200, 12.6363, -0.7355,  0.0488,  0.8442, -2.0131,  1.3603,  SiS
     *  10.1252, 25.2575, -0.6810, -0.3051, -1.5765,  2.7536, -1.8355,  AlOH
     *   7.4400, 23.7389, -1.0179, -0.9947, -1.4353,  3.2530, -1.9224,  NH2
     *   3.0600, 11.4907, -0.4322, -0.6561, -0.5978,  2.4923, -2.4038,  AlH
     *   7.7600, 12.4438, -0.4756, -0.4909, -1.4623,  2.6823, -1.5396,  CN
     *  16.5382, 26.9571, -0.7464, -0.4921, -0.8506, -0.1365,  0.2358,  CO2
     *   5.3590, 12.3380, -0.4956, -0.2251, -0.1907, -0.2038,  0.2579,  SO
     *   6.8700, 11.9229, -1.4044,  0.7899, -0.7317, -0.0193, -0.4994,  TiO
     *   4.3693, 12.3190, -0.5050, -0.0290, -0.0266, -0.6002,  0.4572,  S2
c    *   2.4100, 12.1214,  0.9438,  2.2756, -0.1086,  4.1281, -1.9952,  FeH
c Dissociation energy from Dulick 2003
     *   1.5980, 12.1214,  0.9438,  2.2756, -0.1086,  4.1281, -1.9952,  FeH
     *  12.1388, 36.6661, -1.4062, -0.9258, -1.6969,  0.6005,  1.2302,  NH3
     *  13.2363, 25.1318, -0.5532, -0.0850, -0.9817,  0.6676,  0.3054,  HCN
     *  11.8560, 24.6414, -0.9415, -0.1856, -0.2948, -0.1630,  0.5836,  HCO
     *   5.1156, 12.8758, -0.4856, -0.5054, -0.0776, -0.0713,  0.2369,  O2
     *   7.9400, 23.8609, -1.0762, -0.4928, -0.4092,  0.0031,  0.3761,  CH2
     *   5.8690, 12.2896, -0.9180, -0.6238,  0.1243, -0.3525,  0.4767,  HF
     *   0.0000, 18.8343, 12.4131, 11.9991,  6.8079,  8.4071,  2.6202,  H3+
c    *  16.0000, 18.8343, 12.4131, 11.9991,  6.8079,  8.4071,  2.6202,  H3+
     *   1.7000, 10.1982, -0.9309,  1.8315, -5.6059,  6.9571, -3.5023,  CaH
     *  10.9653, 24.8807, -0.0033,  0.4796, -1.6979,  3.5631, -2.5414,  Al2O
     *   5.2700, 12.2132, -0.5246, -0.1918, -0.6810,  1.7287, -1.5839,  AlO
     *  12.6885, 36.6540, -1.3373, -1.0064, -0.5880, -0.2362,  0.8764,  CH3
     *   0.0000, 17.8513,-15.5361,-17.6144,-13.1604, -6.4819, -5.6361,  SiH2
     *   3.5300, 10.7940,  0.0122,  1.1189, -1.8758,  2.9976, -2.7758,  MgO
c     *   6.2100, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
     *   6.2970, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
c     *   6.3710, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
     *  13.2915, 25.9340, -1.4243,  1.6519, -0.7240, -0.7271,  0.7518,  TiO2
     *  12.9619, 25.9238, -1.2927,  1.3710, -2.4073,  2.2875, -0.5486,  VO2
     *   1.8800, 10.7184, -0.3642,  0.7843, -6.5309, 13.2912, -9.9502,  NaH
     *   5.1200, 11.8277, -0.3468, -1.0735,  1.8038, -1.7748,  0.4333,  AlCl
     *   6.8900, 12.2422, -0.4905, -0.4198,  0.0242,  0.3868, -0.5765,  AlF
     *   6.4100, 12.8108, -0.5811, -0.7895, -2.6766,  8.5158, -6.9993,  VO
     *   7.3550, 12.8487, -0.7627, -0.2538,  1.5240, -4.0119,  3.0234,  CS
     *   8.0735, 23.3256, -0.5884,  0.3637, -2.4401,  3.3936, -1.7121,  MgOH
     *  11.7451, 25.2051, -0.9105,  1.0031, -0.7207, -1.1064,  1.6239,  PO2
     *   8.7035, 23.1900, -1.0964,  2.5340, -5.9823,  5.3416, -1.1946,  CaOH
     *   6.4895, 23.0863, -1.3781,  0.2539, -0.6746, -1.2341,  1.5623/  PH2
      DATA C02/
     *  12.2087, 24.9752, -0.3204, -0.5640, -0.8997,  1.6927, -0.7771,  C2H
     *   6.9600, 12.5225, -1.2695,  1.7628, -2.0543, -1.2215,  2.3706,  ScO
     *  15.6364, 37.7022, -0.5885, -0.0823, -1.7283,  3.0502, -2.0176,  AlO2H
     *   3.8400, 11.9140, -0.5187, -0.1193, -0.3886,  1.1704, -1.2299,  AlS
     *   4.2000, 12.5326, -1.0657,  1.0360, -1.5641,  0.9560, -0.3218,  FeO
     *   4.4000, 11.0587, -1.3926,  1.4461, -2.1552,  3.3409, -3.1078,  CrO
     *  17.2173, 49.9426, -0.9720, -2.4957, -0.0017, -2.3299,  3.1042,  CH4
     *   4.8000, 11.9223, -0.6951,  0.1870, -0.7158,  0.4121,  0.0296,  NS
     *  11.1405, 25.9246, -0.5809,  0.0734, -0.3333,  0.1699,  0.0529,  SO2
     *   6.6880, 14.0972,  4.2904,  4.9608,  2.9390,  3.9789,  0.8908,  SiN
     *   4.7600, 19.9888, -6.7088, -4.3846, -2.8142, -2.3004, -0.3157,  OH-
     *   7.8500, 12.4674, -1.1280,  0.0368,  0.2221,  1.1043, -1.8804,  ZrO
     *  10.8500, 17.5169, 33.0097, 36.2110, 26.7396, 15.2392, 11.4130,  NO+
     *  14.4650, 25.6324, -1.5339,  1.1586, -0.9355,  1.6114, -1.2154,  ZrO2
     *   8.2800, 12.6246, -0.6966, -0.3874,  0.2531, -0.7582,  0.5307,  BO
     *  13.0355, 26.5610, -0.2891,  0.3006, -0.4009,  0.5864, -0.4006,  SiO2
     *  12.7425, 25.2283, -0.4780, -0.3611, -0.2189, -0.2108,  0.5883,  HBO
     *   4.6400, 11.8909, -0.8762,  0.1138,  0.0665, -0.5226,  0.3331,  SiC
     *  15.2000, 25.8617, -1.4050, -0.3896,  1.0805,  2.9269, -3.7531,  YO2
     *   4.7500, 11.6628, -1.4463,  1.3742, -0.8127, -0.4623,  0.2288,  TiS
     *  19.0991, 38.4541, -0.7808, -0.4220, -0.9239,  1.0793, -0.2304,  HBO2
     *  16.9704, 37.7481, -0.2529, -1.0622, -0.1485, -0.7058,  1.1910,  C2H2
     *  14.3762, 26.3815, -0.1712,  0.1197,  0.0059, -0.9891,  1.1946,  OCS
     *   0.0000,  2.5576, -0.5567, -4.5109, -4.3690, -0.1528, -3.1319,  ZrO+
     *   8.0150, 23.3420, -0.6139,  1.4091, -6.8466, 13.0407, -9.2977,  NaOH
     *   4.0900, 10.6268, -1.1367,  2.5278, -5.6022,  4.8741, -1.1616,  CaCl
     *  12.9003, 25.5751, -0.0730,  0.2808, -1.1757,  2.3733, -1.6726,  AlOF
     *   7.2900, 12.4422, -1.3547,  1.3087,  0.1688, -5.4106,  5.1158,  YO
     *   4.2300, 11.0864, -0.4463,  1.1926, -7.5820, 15.2552,-11.1116,  NaCl
     *  14.5371, 25.6134, -0.0508,  0.3710, -0.6246, -0.7682,  0.5868,  C2O
     *  11.4442, 24.7107, -0.5678, -0.0389,  1.0076, -4.6514,  4.3893,  CHP
     *   3.7900, 19.0227, -8.0668, -5.9821, -3.8685, -3.1838, -1.0364,  HS-
     *   0.7300, 19.7162, -5.0018, -2.7680, -1.2845, -0.9859, -0.3380,  H2-
     *   2.1200, 12.4717,  0.1601,  1.4596, -0.2012,  5.0788, -4.5487,  TiH
     *   9.7800, 35.8044, -1.3937, -0.2650, -0.6732, -2.5437,  2.9710,  PH3
     *   2.4000, 11.3146, -0.5595,  0.3619, -2.0065,  3.8766, -2.9900,  MgS
C 30-dec-2008 NP: added the dissociation energy from NIST
C
     *   0.0000,  4.5751,  3.4421,  0.7560, -1.7011,  1.4510, -1.3922,  TiO+
C    *  13.6890,  4.5751,  3.4421,  0.7560, -1.7011,  1.4510, -1.3922,  TiO+
     *  21.1510, 31.0805, 10.7070, 12.8687, 10.5799,  6.4414,  3.6171,  LaO2
     *   3.2100, 12.1817, -0.7102, -0.2403,  1.1042, -1.3644,  0.3198,  Si2
     *  13.2716, 48.6914, -1.0602, -1.2802, -0.8603,  0.1159, -0.0701,  SiH4
     *   8.2349, 24.0157, -0.6514, -0.6064, -0.6542,  0.9096, -0.5839,  BH2
     *  10.9011, 25.1839, -0.1060,  0.2530, -1.1850,  2.3355, -1.6111,  AlOCl
     *   8.2300, 12.1920,  0.1751, -0.7678, -1.3836,  1.7704, -0.0265,  LaO
     *  14.0629, 25.1475, -0.2270,  0.7024, -0.8499,  0.4583,  0.1889,  C2N
     *  20.0747, 38.6719, -0.2664,  0.2782, -1.2642,  1.6020, -0.5248,  AlBO2
     *   4.3400, 10.9561, -0.8720,  3.4218,-12.2306, 18.7863,-11.1011,  KCl
     *   3.2300, 19.3359, -5.7570, -3.5853, -1.3882, -2.3313, -0.4930,  SiH-
     *   5.4800, 11.0459, -0.8574,  2.3137, -4.6777,  4.4532, -1.1716,  CaF
     *  17.8875, 47.4921, -1.1390,  2.7534, -7.2248,  6.3242, -1.1381,  CaO2H2
     *   8.1892, 23.3129, -1.0581,  3.5131,-11.3115, 16.9078, -9.8867/  KOH
      DATA C03/
     *  10.3100, 21.7682, -5.8992, -3.8627, -4.0284,  1.2924, -2.5856,  CN-
     *  16.1405, 37.9519, -0.0230,  0.6639, -2.4910,  5.5385, -4.2945,  Al2O2
     *   9.0621, 23.3478, -2.1422,  1.7058, -1.6807, 10.3429,-14.0183,  BaOH
     *   8.6837, 23.1042, -1.2656,  3.2436, -7.2017,  6.5067, -1.7129,  SrOH
     *  13.9839, 25.6721, -0.0784,  0.0544, -0.2755,  0.6140, -0.3673,  BO2
     *   5.5700, 12.0158, -0.5187, -0.1216,  0.6738, -0.6377,  0.1588,  SiF
C
C 30-dec-2008 NP: added dissociation energy as dissociation energy of CH
C                 (3.465eV) + electron affinity of CH (1.238eV from NIST)
     *   0.0000, 16.4621,-13.8562,-13.1896, -9.2577, -6.3354, -2.5704,  CH-
C    *   4.7030, 16.4621,-13.8562,-13.1896, -9.2577, -6.3354, -2.5704,  CH-
     *  13.8610, 26.3081, -1.3134,  0.1185, -0.0461, -0.4056,  0.8088,  C3
     *   8.4800, 21.1413, -5.8697, -3.3745, -2.7491, -1.8902, -0.2441,  C2-
     *  17.1545, 48.1845, -0.5683,  0.1125, -3.0973,  4.3727, -2.1978,  MgO2H2
     *   9.3961, 23.7967, -0.6500,  0.2061, -1.9381,  2.1259, -0.6451,  BeOH
     *  10.4305, 24.8357, -0.4930, -0.4550,  0.8862, -2.7257,  2.4025,  HBS
     *  13.1966, 25.7392,  0.0961, -0.7979, -0.1515,  4.2750, -4.6336,  SiC2
     *  17.4231, 48.8561, -0.4831,  0.9575, -1.9798, -0.0476,  1.2346,  FeO2H2
     *  10.0930, 25.0689, -1.5784,  2.2605, -3.1152,  3.7375, -2.5596,  CrO2
     *  20.0817, 49.3051, -0.2203,  0.6123, -1.9159,  3.0362, -0.6588,  BeH2O2
     *  11.4541, 36.8342, -1.3068, -1.2283, -0.7130, -0.1039,  0.8121,  BH3
     *  12.5346, 24.2744, -0.4230,  2.1003, -7.6565, 14.5171,-10.4377,  NaCN
     *   6.5483, 23.5736, -0.7830, -0.0881, -2.2398,  2.7050, -1.5244,  BeH2
     *  10.1248, 24.8268, -0.3784,  0.5561, -0.7324,  1.7508, -1.6977,  Si2N
     *   9.3132, 22.5681, -0.7730,  3.2979, -6.3686,  5.5210, -0.9987,  CaCl2
     *  18.8913, 37.0212, -0.3881,  1.7934, -7.5472, 14.9782,-11.0505,  NaBO2
     *   0.0000, 19.8338,-46.6804,-50.9308,-35.9059,-13.5611,-23.8103,  C3H
     *  15.5315, 26.0301, -0.1824,  0.0109, -0.3944,  0.5184, -0.0882,  OBF
     *  11.9993, 26.2368, -0.1708,  0.2491,  0.4220, -2.2962,  2.2409,  CS2
     *   8.9381, 23.5703, -0.6263,  1.0060, -4.3983,  7.4665, -4.8955,  LiOH
     *   1.5500, 11.3681, -0.1946, -0.0669, -2.3347,  5.3477, -4.0343,  Al2
     *   4.8400, 11.3090, -0.5602,  0.5886, -3.9705,  7.3873, -5.2571,  LiCl
     *  11.3225, 25.4462, -1.0487,  1.8142, -1.5110,  0.4282, -0.0240,  TiOCl
     *  23.3326, 62.7915, -1.3095, -1.6903, -0.9624, -1.6171,  2.5521,  C2H4
     *   7.4689, 23.8059, -0.5629,  0.0019, -0.3896, -0.7781,  0.3890,  CHCl
     *   6.6900, 14.8883,  5.3193,  8.9551,  3.7271,  5.1452,  1.0391,  TiCl
     *  19.2284, 37.1933,  0.1308, -0.0614, -0.9981,  2.9770, -2.1833,  AlOF2
     *  18.9713, 36.8674, -0.8338,  3.8816,-11.3916, 16.8414, -9.6911,  KBO2
     *  11.2271, 25.9412,  0.1074, -0.8813, -0.2594,  4.4112, -4.4861,  Si2C
     *   9.2183, 24.5270, -0.6453, -1.0757, -0.7155,  2.2944, -1.4513,  CHF
     *   0.0000, 11.8175,-29.4442,-30.6402,-22.9279,-13.1209, -8.8023,  BO-
     *  10.9760, 27.6834,  5.5082,  6.6402,  5.5692,  2.7324,  1.9375,  AlO2
     *  18.0802, 47.0050, -2.3587,  2.3466, -2.2753,  8.4432,-11.3032,  BaO2H2
     *  12.8526, 25.8889, -1.0260,  1.8361, -1.5017,  0.3478,  0.0486,  OTiF
     *   6.5000, 20.6745, -7.9942, -5.7057, -2.6759, -6.1649,  1.2656,  CS-
     *  21.5636, 39.0495, -0.1190,  0.7088, -1.5184,  0.4914,  0.9277,  C2N2
     *  17.5958, 46.9386, -1.3295,  3.5725, -8.4710,  7.5694, -1.8456,  SrO2H2
     *  12.2076, 25.3442, -0.0379, -0.1189, -0.8276,  1.3188, -0.6986,  ClCN
     *  10.6135, 23.6489, -0.5207,  0.0519, -0.6538,  1.9149, -1.5058,  AlClF
     *  12.5010, 24.1386, -0.8692,  4.1888,-11.7377, 17.1662, -9.8522,  KCN
     *   8.8688, 23.5425, -0.5528,  0.0031, -0.7346,  2.3344, -1.9878,  AlCl2
     *   9.6070, 22.2204, -2.5275,  2.8555, -1.4987,  7.7865,-11.3039,  BaCl2
     *  12.3143, 24.3964, -0.4940,  0.0699, -0.5475,  1.6261, -1.2695,  AlF2
     *   8.1536, 22.9187, -0.1815,  0.6847, -2.4792,  4.3296, -2.7691/  MgCl2
      DATA C04/
     *   0.0000, 17.5598,-16.6727,-14.0707,-13.0780, -5.4193, -4.7856,  FeO-
     *  20.4537, 49.9913, -0.5362, -0.7176, -1.2169,  1.1206, -0.3773,  BO2H2
     *  14.1133, 48.5194, -0.8436, -1.0629, -0.7362,  0.3080, -0.3403,  SiH3Cl
     *   8.3239, 23.6272, -0.2108,  1.1105, -2.1105,  1.5380, -0.1684,  FeCl2
     *   7.3840, 24.8600, -0.1499, -0.1631,  0.1378,  1.6604, -1.9986,  Si3
     *  16.1268, 48.9782, -0.8260, -1.0380, -0.6452, -0.1029,  0.1199,  SiH3F
     *  16.2992, 49.7196, -1.2716, -1.4752, -1.1626,  0.6516, -0.0837,  CH3Cl
     *   9.1791, 22.1133, -1.4891,  4.1050, -7.6534,  6.6694, -1.5355,  SrCl2
     *  11.6845, 23.2600, -1.2039,  3.3661, -6.2828,  5.1661, -0.6547,  CaF2
     *  13.7563, 25.2856, -0.4137,  1.0746, -1.1248,  0.2935,  0.3807,  TiF2
     *  19.4163, 36.9346, -0.3977,  1.3814, -4.7577,  8.2956, -5.5779,  LiBO2
     *   9.5422, 23.6489, -0.6541,  0.7042, -2.5258,  4.5411, -3.0359,  MgClF
     *  19.3953, 37.4967, -0.4103,  0.6249, -2.5737,  3.7334, -2.0769,  BeBO2
     *  16.1988, 37.8077, -0.3545, -0.2428, -0.1731, -1.4896,  1.9844,  C2HCl
     *   9.9277, 24.6274, -0.5062,  0.9860, -1.3100,  0.8075, -0.0931,  TiCl2
     *  19.7168, 40.3256, -0.2533,  0.3731, -0.5863, -0.6939,  0.9337,  C4
     *  30.6562, 75.8041, -1.6269, -1.1205, -1.8109,  2.1354, -0.8357,  H3BO3
     *  10.7510, 23.8686, -0.6130,  0.7434, -2.6657,  5.0507, -3.5509,  MgF2
     *   0.0000, 13.8534,-28.5088,-27.6557,-25.0420, -4.2145,-21.0916,  BaClF
     *  13.3200, 24.6323, -0.2099,  0.5174, -1.9085,  2.9836, -1.7351,  BeF2
     *  16.6788, 38.1093, -0.3632, -0.2642, -0.4287, -0.5573,  0.9863,  C2HF
     *   9.6498, 23.7877, -0.2606,  0.4816, -1.7048,  2.1226, -0.8176,  BeCl2
     *  15.7352, 37.1910, -1.0480,  1.8371, -1.1420, -0.7526,  1.2880,  TiOCl2
     *  10.7683, 24.3508, -0.5859,  0.0972, -0.3635,  0.9082, -0.3338,  ZrCl2
     *  11.9101, 22.9073, -2.4413,  2.9420, -1.3655,  7.3312,-10.8692,  BaF2
     *  12.4073, 25.2586, -0.5256,  0.7548, -2.0655,  2.2598, -0.9944,  BeC2
     *   9.9676, 24.0020, -0.4765,  1.0925, -3.6131,  4.2582, -1.8225,  Be2O
     *  11.3542, 22.8132, -1.4157,  4.1790, -7.3508,  5.5696, -0.4507,  SrF2
     *  13.7587, 24.7160, -1.0103,  0.2376, -0.4664, -0.9114,  6.9672,  ZrF2
     *  13.0910, 27.6502,  6.5468,  8.2502,  7.3334,  4.1191,  1.2402,  FeF2
     *  12.5389, 37.9053, -1.3490,  3.1985, -1.1165, -6.7253,  7.3584,  P4
     *  19.0240, 49.7099, -0.5565, -0.7375, -0.2251, -1.1324,  1.2457,  SiH2F2
     *   3.2806, 41.7329, 32.0127, 34.5233, 27.1981, 13.3168, 13.4808,  H3O+
     *  27.0859, 54.0398,  0.0077,  0.4169, -0.9261, -0.3135,  0.6322,  C5
     *  19.7864, 37.9176, -0.7063,  1.7895, -1.5401,  0.9448, -0.6313,  TiF3
     *  14.3199, 37.3165, -0.8450,  1.6603, -1.6009,  0.8934, -0.5070,  TiCl3
     *  15.5540, 36.5254, -0.7361,  0.8503, -0.3688,  0.0324,  0.0881,  ZrCl3
     *  10.6603, 34.6664, -0.4567,  3.2641,-13.6211, 27.6173,-20.7914,  Na2Cl2
     *  18.1954, 60.7438, -0.7643,  2.2577,-14.4187, 28.3225,-20.4866,  (NaOH)2
     *  28.8149, 64.3940, -0.2174,  1.3367, -6.6368,  8.6309, -4.6284,  Be3O3
     *  10.8345, 33.9871, -1.3140,  7.4840,-21.9583, 33.6428,-20.3143,  K2Cl2
     *  18.3196, 60.4179, -1.6298,  6.4524,-22.9230, 33.8810,-20.0092,  (KOH)2
     *  20.4364, 49.7173, -0.6667,  0.8064, -0.1308, -0.4433,  0.8970,  ZrCl4
     *  27.1266, 62.7471, -0.3813,  3.6624,-15.0927, 27.0694,-18.7738,  (NaCN)2
     *  27.0557, 51.2712, -0.5271,  0.8930, -0.5666,  1.5292, -1.3568,  ZrF4
     *  20.3442, 61.3686, -0.8410,  1.3617, -9.5297, 16.1158,-11.1739,  (LiOH)2
     *   1.9300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000/  CrH
C
C Coefficients for constructing partition functions (and then equilibrium
C   constants, perhaps). For diatomic molecules other than H2 and CO, the
C   data are from Sauval & Tatum (1984, ApJS, 56, 193). For H2 and CO, the
C   data are from Irwin (1987, A&A, 182, 348). For polyatomic molecules,
C   the coefficients are from Irwin (1988, A&AS, 74,145).
C Coefficients used to construct the partition function, as follows:
C
C     log10(Q) = Sum{i=0,9}{PCOEF(i+1)*log10(THETA)**i}
C                                                           Ioniz. pot.
      DATA P01/
     *   1.69179,      -1.72270,       0.798033,     -0.157089,         H2
     *  -0.535313,      1.75818,      -2.63895,       1.35708,          H2
     *   0.0,           0.0,                                 15.42593,  H2
     *   3.615300,     -1.773848,      0.3516181,     0.08620792,       CO
     *   0.2911791,    -1.141469,      2.513133,     -2.886502,         CO
     *   1.238932,      0.0,                                 14.01400,  CO
     *   4.344711818,  -3.6343233,     1.415963,      0.01594,          H2O
     *   0.56542,      -1.2583,        0.53796,       3*0.0, 12.62100,  H2O
     *   3.0929, -1.6778,  0.6743, -0.1874,  0.0000,  5*0.0, 13.01700,  OH
     *   3.2643, -1.7303,  0.4192,  0.0000,  0.0000,  5*0.0, 15.58100,  N2
     *   4.2275, -1.9144,  0.7201, -1.3099,  1.1657,  5*0.0, 11.49000,  SiO
     *  1.0, 9*0.,                                           10.42200,  HS
     *   5.117210341,  -3.94844146,    1.23193,       0.076156,         H2S
     *   0.42163,      -0.453534,      0.0,           3*0.0, 10.45700,  H2S
     *   3.0735, -1.8501,  0.9607, -0.3935,  0.0000,  5*0.0, 13.49000,  NH
     *   3.6908, -1.9801,  0.7704, -0.2247,  0.0000,  5*0.0,  7.91000,  SiH
     *   3.3586, -2.0656,  0.9624, -0.2239,  0.0000,  5*0.0, 10.64000,  CH
     *   2.5410, -2.4336,  1.4979,  0.0192, -0.7483,  5*0.0, -1.00000,  H2+
     *   4.3073, -1.8255,  0.3765,  0.0000,  0.0000,  5*0.0,  9.26420,  NO
     *   3.6704, -2.2682,  0.9354, -0.2597,  0.0000,  5*0.0,  7.20000,  MgH
     *   2.8005, -1.7476,  0.5310,  0.0000,  0.0000,  5*0.0, 12.74400,  HCl
     *   4.8026, -1.9753,  0.2600,  0.0000,  0.0000,  5*0.0, 10.53000,  SiS
     *   6.103792598,  -4.3938712,     0.662588,      0.3751,           AlOH
     *   0.38386,      -0.2147,        0.0,           3*0.0, -1.00000,  AlOH
     *   4.819621858,  -3.84200734,    1.5386462,     0.784399,         NH2
     *  -2.34404,       2.50803,      -1.13304,       3*0.0, 11.14000,  NH2
     *   3.3209, -2.5909,  1.7415, -0.7636,  0.0000,  5*0.0,  5.50000,  AlH
     *   4.0078, -2.1514,  0.9226, -0.1671,  0.0000,  5*0.0, 13.59800,  CN
     *   6.01081285,   -4.438833,      0.840462,      0.2945,           CO2
     *   0.3694,       -0.273,         0.0,           3*0.0, 13.77700,  CO2
     *   4.7963, -2.1308,  0.5224,  0.0000,  0.0000,  5*0.0, 10.29400,  SO
C The line with 5.7765 is from Alard and Hauschildt who artificially increased
C TiO parition function by a factor of 3. Also change in ionization energy
C according to the latest NIST data.
C    *   5.7765, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.40000,  TiO
     *   5.3051, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.81900,  TiO
     *   5.0796, -2.1967,  0.4101,  0.0000,  0.0000,  5*0.0,  9.35600,  S2
     *   4.6265980,    -2.5625800,     0.38885943,    0.40219820,       FeH
     *  -0.21386399,    0.027845045,   0.0,           3*0.0,  7.37000,  FeH
     *   5.884176216,  -5.8364867,     1.608417,      1.50876,          NH3
     *  -0.59607,      -0.58961,       0.2459,        3*0.0, -1.00000,  NH3
     *   5.434042379,  -4.2409874,     0.988745,      0.49464,          HCN
     *   0.03719,      -0.22924,       0.0,           3*0.0, 13.60000,  HCN
     *   6.298781639,  -3.85672804,    0.8551678,     0.321901,         HCO
     *   0.020274,      0.15254,      -0.25298,       3*0.0,  8.12000,  HCO
     *   4.0636, -2.0779,  0.7660, -0.2111,  0.0000,  5*0.0, 12.06970,  O2
     *  1.0, 9*0.,                                           10.39600,  CH2
     *   2.4164, -1.6132,  0.6357, -0.1767,  0.0000,  5*0.0, 16.03000,  HF
     *  1.0, 9*0.,                                           -1.00000,  H3+
     *   3.8411, -2.3891,  1.3578, -0.6893,  0.0000,  5*0.0,  5.86000,  CaH
     *  1.0, 9*0.,                                           -1.00000,  Al2O
     *   4.9191, -2.6291,  0.5831,  0.3163,  0.0000,  5*0.0,  9.46000,  AlO
     *  1.0, 9*0.,                                            9.84000,  CH3
     *  1.0, 9*0.,                                            8.80000,  SiH2
     *   5.3182, -2.6502, -0.2781, -0.7823,  1.3107,  5*0.0,  8.76000,  MgO
     *   4.3091, -2.2406,  0.4865, -0.2049,  0.0000,  5*0.0, 11.40000,  C2
     *  1.0, 9*0.,                                            9.50000,  TiO2
     *   8.457240767,  -4.1987868,     0.334575,      0.20744,          VO2
     *   0.18226,      -0.053465,      0.0,           3*0.0, -1.00000,  VO2
     *   3.5453, -2.3457,  0.8557, -0.1685,  0.0000,  5*0.0,  4.70000,  NaH
     *   5.1115, -2.2303,  0.8001, -0.5192,  0.0000,  5*0.0,  9.40000,  AlCl
     *   4.5405, -2.1033,  0.6208, -0.2930,  0.0000,  5*0.0, -1.00000,  AlF
     *   5.0687, -2.2186,  0.9545, -0.4592,  0.0000,  5*0.0,  7.23860,  VO
     *   4.1646, -1.9348,  0.8034, -1.3669,  1.1561,  5*0.0, 11.33000,  CS
     *   6.8401894714, -4.338616427,   0.71600166,    0.128126,         MgOH
     *   0.5978087,    -0.8658369,     0.385049,      3*0.0,  7.50000,  MgOH
     *  1.0, 9*0.,                                           11.90000,  PO2
     *   7.1623971155, -4.471282563,   1.1221899,    -0.558812,         CaOH
     *   0.2294,        1.78658,      -2.95118,       1.41591,          CaOH
     *   2*0.0,                                               5.80000,  CaOH
     *  1.0, 9*0.,                                            9.82400/  PH2
      DATA P02/
     *  1.0, 9*0.,                                           11.61000,  C2H
     *   4.8065, -2.2129,  0.9991, -0.5414,  0.0000,  5*0.0, -1.00000,  ScO
     *  1.0, 9*0.,                                           -1.00000,  AlO2H
     *   5.2461, -2.1319,  0.5340, -0.2309,  0.0000,  5*0.0, -1.00000,  AlS
     *   5.5642, -2.1947,  0.5065,  0.0000,  0.0000,  5*0.0,  8.90000,  FeO
     *   5.5270, -2.1311,  0.6523, -0.2533,  0.0000,  5*0.0,  7.85000,  CrO
     *  1.0, 9*0.,                                           12.61000,  CH4
     *   4.8052, -1.9619,  0.3140,  0.0000,  0.0000,  5*0.0,  8.87000,  NS
     *  1.0, 9*0.,                                           12.34900,  SO2
     *   4.6570, -2.3587,  0.8819, -0.1642,  0.0000,  5*0.0, -1.00000,  SiN
     *  1.0, 9*0.,                                           -1.00000,  OH-
     *   5.3279, -2.4694,  0.2164, -0.2313,  0.0000,  5*0.0,  6.00000,  ZrO
     *   3.5649, -1.7328,  0.4241,  0.0000,  0.0000,  5*0.0, -1.00000,  NO+
     *   8.72011985,   -4.247295,      0.2758,        0.20738,          ZrO2
     *   0.09406,       0.0,           0.0,           3*0.0, -1.00000,  ZrO2
     *   3.9953, -1.8665,  0.5965, -0.1617,  0.0000,  5*0.0, 13.30000,  BO
     *  1.0, 9*0.,                                           -1.00000,  SiO2
     *  1.0, 9*0.,                                           -1.00000,  HBO
     *   5.1477, -1.8671,  0.2404,  0.0000,  0.0000,  5*0.0,  9.20000,  SiC
     *  1.0, 9*0.,                                           -1.00000,  YO2
     *   5.8948, -2.2183,  0.5928, -0.3106,  0.0000,  5*0.0,  7.10000,  TiS
     *  1.0, 9*0.,                                           -1.00000,  HBO2
     *   7.1220464309, -6.966653604,   1.9668235,     0.362597,         C2H2
     *   0.608996,     -0.920435,      0.271892,      3*0.0, 11.40000,  C2H2
     *  1.0, 9*0.,                                           11.18500,  OCS
     *  1.0, 9*0.,                                           -1.00000,  ZrO+
     *  1.0, 9*0.,                                           -1.00000,  NaOH
     *   5.7494, -2.3340,  0.8685, -0.5306,  0.0000,  5*0.0,  5.86000,  CaCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF
     *   4.9515, -2.0866,  0.6565, -0.3082,  0.0000,  5*0.0,  6.00000,  YO
     *   5.3364, -2.2844,  0.2820,  0.1185,  0.0000,  5*0.0, -1.00000,  NaCl
     *  1.0, 9*0.,                                           -1.00000,  C2O
     *  1.0, 9*0.,                                           10.79000,  CHP
     *  1.0, 9*0.,                                           -1.00000,  HS-
     *  1.0, 9*0.,                                           -1.00000,  H2-
     *  1.0, 9*0.,                                            6.00000,  TiH
     *  1.0, 9*0.,                                            9.86900,  PH3
     *   5.0367, -2.1625,  0.4859, -0.1780,  0.0000,  5*0.0, -1.00000,  MgS
     *  1.0, 9*0.,                                           -1.00000,  TiO+
     *  1.0, 9*0.,                                           -1.00000,  LaO2
     *   5.2617, -2.1485,  0.5647, -0.2985,  0.0000,  5*0.0, -1.00000,  Si2
     *  1.0, 9*0.,                                           -1.00000,  SiH4
     *  1.0, 9*0.,                                            9.80000,  BH2
     *  1.0, 9*0.,                                           -1.00000,  AlOCl
     *   5.1147, -2.5016,  1.0445, -0.3135,  0.0000,  5*0.0,  4.95000,  LaO
     *  1.0, 9*0.,                                           12.00000,  C2N
     *  1.0, 9*0.,                                           -1.00000,  AlBO2
     *   5.6860, -2.3016,  0.2086,  0.1763,  0.0000,  5*0.0, -1.00000,  KCl
     *  1.0, 9*0.,                                           -1.00000,  SiH-
     *   5.2010, -2.2653,  0.8941, -0.5384,  0.0000,  5*0.0, -1.00000,  CaF
     *  1.0, 9*0.,                                           -1.00000,  CaO2H2
     *  1.0, 9*0.,                                            7.50000/  KOH
      DATA P03/
     *  1.0, 9*0.,                                           -1.00000,  CN-
     *  1.0, 9*0.,                                           -1.00000,  Al2O2
     *  1.0, 9*0.,                                           -1.00000,  BaOH
     *  1.0, 9*0.,                                           -1.00000,  SrOH
     *  1.0, 9*0.,                                           -1.00000,  BO2
     *   5.0871, -2.0375,  0.4478, -0.1243,  0.0000,  5*0.0,  7.54000,  SiF
     *  1.0, 9*0.,                                           -1.00000,  CH-
     *   6.618407932,  -3.576399,      0.883642,      0.087548,         C3
     *   0.04817,      -0.16471,       0.0,           3*0.0, -1.00000,  C3
     *  1.0, 9*0.,                                           -1.00000,  C2-
     *  1.0, 9*0.,                                           -1.00000,  MgO2H2
     *  1.0, 9*0.,                                           -1.00000,  BeOH
     *  1.0, 9*0.,                                           -1.00000,  HBS
     *   7.54651307623,-5.075563869,   1.82960795,    0.0983258,        SiC2
     *  -6.335157,     14.33103,     -13.01689,       4.428233,         SiC2
     *   2*0.0,                                              10.20000,  SiC2
     *  1.0, 9*0.,                                           -1.00000,  FeO2H2
     *  1.0, 9*0.,                                           -1.00000,  CrO2
     *  1.0, 9*0.,                                           -1.00000,  BeH2O2
     *  1.0, 9*0.,                                           -1.00000,  BH3
     *  1.0, 9*0.,                                           -1.00000,  NaCN
     *  1.0, 9*0.,                                           -1.00000,  BeH2
     *  1.0, 9*0.,                                           -1.00000,  Si2N
     *  1.0, 9*0.,                                           -1.00000,  CaCl2
     *  1.0, 9*0.,                                           -1.00000,  NaBO2
     *  1.0, 9*0.,                                           -1.00000,  C3H
     *  1.0, 9*0.,                                           -1.00000,  OBF
     *  1.0, 9*0.,                                           10.07300,  CS2
     *  1.0, 9*0.,                                           -1.00000,  LiOH
     *   5.5538, -2.3365,  0.5754, -0.2119,  0.0000,  5*0.0,  5.40000,  Al2
     *   4.5605, -2.2216,  0.5760, -0.1706,  0.0000,  5*0.0,  9.57000,  LiCl
     *  1.0, 9*0.,                                           -1.00000,  TiOCl
     *  1.0, 9*0.,                                           -1.00000,  C2H4
     *  1.0, 9*0.,                                           -1.00000,  CHCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF2
     *  1.0, 9*0.,                                           -1.00000,  KBO2
     *  1.0, 9*0.,                                           -1.00000,  Si2C
     *  1.0, 9*0.,                                           10.06000,  CHF
     *  1.0, 9*0.,                                           -1.00000,  BO-
     *  1.0, 9*0.,                                           -1.00000,  AlO2
     *  1.0, 9*0.,                                           -1.00000,  BaO2H2
     *  1.0, 9*0.,                                           -1.00000,  OTiF
     *  1.0, 9*0.,                                           -1.00000,  CS-
     *  1.0, 9*0.,                                           -1.00000,  C2N2
     *  1.0, 9*0.,                                           -1.00000,  SrO2H2
     *  1.0, 9*0.,                                           12.36000,  ClCN
     *  1.0, 9*0.,                                           -1.00000,  AlClF
     *  1.0, 9*0.,                                           -1.00000,  KCN
     *  1.0, 9*0.,                                           -1.00000,  AlCl2
     *  1.0, 9*0.,                                           -1.00000,  BaCl2
     *  1.0, 9*0.,                                           -1.00000,  AlF2
     *  1.0, 9*0.,                                           -1.00000/  MgCl2
      DATA P04/
     *  1.0, 9*0.,                                           -1.00000,  FeO-
     *  1.0, 9*0.,                                           -1.00000,  BO2H2
     *  1.0, 9*0.,                                           -1.00000,  SiH3Cl
     *  1.0, 9*0.,                                           -1.00000,  FeCl2
     *  1.0, 9*0.,                                           -1.00000,  Si3
     *  1.0, 9*0.,                                           -1.00000,  SiH3F
     *  1.0, 9*0.,                                           -1.00000,  CH3Cl
     *  1.0, 9*0.,                                           -1.00000,  SrCl2
     *  1.0, 9*0.,                                           -1.00000,  CaF2
     *  1.0, 9*0.,                                           -1.00000,  TiF2
     *  1.0, 9*0.,                                           -1.00000,  LiBO2
     *  1.0, 9*0.,                                           -1.00000,  MgClF
     *  1.0, 9*0.,                                           -1.00000,  BeBO2
     *  1.0, 9*0.,                                           -1.00000,  C2HCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl2
     *  1.0, 9*0.,                                           -1.00000,  C4
     *  1.0, 9*0.,                                           -1.00000,  H3BO3
     *  1.0, 9*0.,                                           -1.00000,  MgF2
     *  1.0, 9*0.,                                           -1.00000,  BaClF
     *  1.0, 9*0.,                                           -1.00000,  BeF2
     *  1.0, 9*0.,                                           -1.00000,  C2HF
     *  1.0, 9*0.,                                           -1.00000,  BeCl2
     *  1.0, 9*0.,                                           -1.00000,  TiOCl2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl2
     *  1.0, 9*0.,                                           -1.00000,  BaF2
     *  1.0, 9*0.,                                           -1.00000,  BeC2
     *  1.0, 9*0.,                                           -1.00000,  Be2O
     *  1.0, 9*0.,                                           -1.00000,  SrF2
     *  1.0, 9*0.,                                           -1.00000,  ZrF2
     *  1.0, 9*0.,                                           -1.00000,  FeF2
     *  1.0, 9*0.,                                           -1.00000,  P4
     *  1.0, 9*0.,                                           -1.00000,  SiH2F2
     *  1.0, 9*0.,                                           -1.00000,  H3O+
     *  1.0, 9*0.,                                           -1.00000,  C5
     *  1.0, 9*0.,                                           -1.00000,  TiF3
     *  1.0, 9*0.,                                           -1.00000,  TiCl3
     *  1.0, 9*0.,                                           -1.00000,  ZrCl3
     *  1.0, 9*0.,                                           -1.00000,  Na2Cl2
     *  1.0, 9*0.,                                           -1.00000,  Na2O2H2
     *  1.0, 9*0.,                                           -1.00000,  Be3O3
     *  1.0, 9*0.,                                           -1.00000,  K2Cl2
     *  1.0, 9*0.,                                           -1.00000,  K2O2H2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl4
     *  1.0, 9*0.,                                           -1.00000,  Na2C2N2
     *  1.0, 9*0.,                                           -1.00000,  ZrF4
     *  1.0, 9*0.,                                           -1.00000,  Li2O2H2
     *  1.0, 9*0.,                                            7.33000/  CrH
C
C Try to find the input speicies name (SPNAME) in the list (SPLIST) of
C species for which we have equilibrium constant coefficients. Note that
C the index is stored in a new variable J, rather than using the loop
C variable I, because some optimizers don't save the loop variable after
C normal termination of the loop.
C
      DO 1 I=1,MSPEC
      J=I
      IF(SPLIST(J).EQ.SPNAME) GO TO 2
   1  CONTINUE
C
C Fall through to here, if requested molecule was not in SPLIST.
C Print a warning, but return anyway.
C
      WRITE(*,*) 'MOLCON: Don''t have dissociation constant for ',
     *           'molecule: "', SPNAME, '"'
      EQK =1.D20
      PART=1.D0
      RETURN
C
C Calculate independent variable for polynomial expansions.
C Note that the polynomial expansions in Sauval & Tatum (1984) and Irwin
C (1987,1988) are in terms of log10(5040/T), not log10(5039.7475/T), but
C the more accurate value of 5039.7475 should be used in converting the
C partition function into an equilibrium constant.
C
   2  TLIM=MAX(1250.,T)
      TH=5040.D0/TLIM
      LOGTH=LOG10(TH)
C
C Construct equilibrium constant from polynomial coefficients and
C dissociation constant. A "+1" term at the end would convert from
C pascals (i.e. N/m/m as in Sauval) to dynes/cm/cm.
C
c      if (t.lt.1600) logth=log10(5040.0/1600.0)
c      if (t.gt.7730) logth=log10(5040.0/7730.0)
      EQK=COEF(2,J)+LOGTH*(COEF(3,J)+LOGTH*(COEF(4,J)+
     &              LOGTH*(COEF(5,J)+LOGTH*(COEF(6,J)+
     &              LOGTH*(COEF(7,J))))))
     &             -TH*COEF(1,J)
C    &             +1.0D0
      EQK =10.D0**EQK
C
C Just for the reference, the relation between partition functions
C and equilibrium constant:
C
C            P(A)*P(B)*...      N(A)*N(B)*...
C K(AB...) = ------------- = kT-------------- =
C              P(AB...)           N(AB...)
C
C             2*pi*kT 3/2    M(A)*M(B)*... 3/2   Q(A)*Q(B)*...
C       = kT*(-------)    * (-------------)    * ------------- * exp(-D(AB)/kT)
C               h^2            M(AB...)           Q(AB...)
C
C where, K - equilibrium constant, Q - partition functions, M - masses
C        P - partial pressures, N - number densities, T - temperature,
C        D - complete dissociation energy, h - plank constant. Remember
C        to use masses in grams (1 amu = 1.660540E-24 g) and energy in
C        ergs (1 eV = 1.60219E-12 ergs). Also, k = 1.38065E-16 erg/K,
C        h = 6.626076E-27 erg s, and pi = 3.1415926536.
C
C Construct partition function from polynomial coefficients.
C
      PART=PCOEF(NPCOEF-1,J)
      DO 3 I=NPCOEF-2,1,-1
    3 PART=LOGTH*PART+PCOEF(I,J)
C
C Copy ionization potential
C
      PION=PCOEF(NPCOEF,J)
C
C Calculate equilibrium constant (EQK) from partition function, dissociation
C constant, and other information passed into subroutine. The constants used
C are:  79.733501 = 1.5*log10(2*pi/h/h)  [in cgs units] and
C      -15.859914 = alog10(k)            [in cgs units].
C       5039.7475 = alog10(e)*k*(eV/erg)
C
c      EQK_ST=(NTOT-1)*(79.733501D0+2.5D0*(LOG10(TLIM)-15.859914D0))+
c     &       1.5D0*RATIOM+QPRD-PART-COEF(1,J)*5039.7475D0/TLIM
C
C Convert equilibrium constant and partition function from logarithms.
C
c      EQK_ST=10.D0**EQK_ST
      PART=10.D0**PART
C
C Check if there is relevant data in Paul Barklem's tables
C
      CALL KP_Q_SPLN(SPNAME,T,Qm_spln,Kp_spln,BARKLEM)
      IF(BARKLEM) THEN
c        EQK =Kp_spln-COEF(1,J)*5039.7475D0/TLIM
        EQK =Kp_spln-COEF(1,J)*5040.D0/T
        EQK =10.D0**EQK
        PART=10.D0**Qm_spln
      ENDIF
      if(spname.eq.'H3O+') then
        EQK_ST=(NTOT-1)*(79.733501D0+2.5D0*(LOG10(T)-15.859914D0))+
     &         1.5D0*RATIOM+QPRD-PART-COEF(1,J)*5039.7475D0/T
        EQK=10.D0**EQK_ST
      endif
c      write(*,'(''cMOLCON:'',F10.1,A9,5G13.6)') T,SPNAME,EQK,
c     &                             PART,BARKLEM
c      if(spname.eq.'NO') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'C3') write(*,'(a,f10.2,1p6e16.8)')
c     &   spname,t , eqk, eqk_st, part, TH, LOGTH, TLIM
c      if(spname.eq.'H3O+') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'SiS') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'NO') write(*,'(a,f10.2,1p3e16.8,L)')
c     &   spname,t , eqk, eqk_st, part,barklem
c      if(spname.eq.'CH') write(*,'(a,f10.2,1p5e16.8,L)')
c     &   spname,t , eqk, eqk_st, part,COEF(1,J),Kp_spln,barklem
c      if(spname.eq.'CH-') write(*,'(a,f10.2,1p5e16.8,L)')
c     &   spname,t , eqk, eqk_st, part,COEF(1,J),Kp_spln,barklem
c      if(spname.eq.'CH-') write(*,'(a,f10.2,1p3e14.6,i3,1p2e14.6,L)')
c     &   spname,t , eqk, eqk_st, part,NTOT,QPRD,RATIOM,BARKLEM
c      if(spname.eq.'H2') write(*,'(a,f10.2,1p3e14.6,i3,1p2e14.6)')
c     &   spname,t , eqk, eqk_st, part,NTOT,Kp_spln,COEF(1,J)*5040.D0/T
c
c Don't use EQK_ST based on partition function - use direct fit to EQK.
c
c      EQK=EQK_ST
C
C Done.
C
      RETURN
      END
C---------------------- Start of Berklem subroutines ------------------------
C=========================================================================
C Kp_Q_spln: Returns equilibrium constant for a given molecule and temperature.
C
C Inputs:
C   SPNAME  [char] species name according to the table below.
C   TEMP    [real] temperature (in K) at which Kp is needed
C
C History:
C  28-jun-2007: First version written by N. Piskunov including 57 species.
C               Molecular equilibium tabulated by P. Barklem, resampled
C               for optimal spline interpolation and converted to Fortran
C               DATA statements by J. Valenti
C
C  15-dec-2007: Second version includes 58 molecular species.
C               Tabulated values are now alog10(Kp)+D0*5040/T vs alog10(T),
C               where Kp is an equilibrium constant in N/m^2, D0 is the
C               dissociation energy (eV) at 0 K, and T is temperature (K).
C               In this version, we start using a separate alog10(T) grid
C               for each species, rather than a common THETA=5040/T grid
C               for all species. We copied D0 from MOLCON in eos.f, except
C               for CH-, OH-, SiH-, SiN, and MgS, which we (JV) deduced from
C               Barklem data.
C
C   7-sep-2009: Subroutine data modified and the subroutine text generated
C               by IDL program qk_spl_nodes_f77.pro with errthr=0.000100
C
C Outputs:
C   K_spln [real*8] equilibrium constant (in dynes/cm^2) at temperature T,
C   Q_spln [real*8] partition functions at temperature T,
C          both interpolated from Paul Barklem's tables.
C
C To obtain molecular equilibrium constants, KP:
C
C   D2 = SPL_INIT(TK_<species>,K_<species>)
C   KP(T) = SPL_INTERP(TK_<species>,K_<species>,D2,TLOG)
C         - D0*5040/T
C
C To obtain partition functions,Q:
C
C   D2 = SPL_INIT(TQ_<species>,Q_<species>)
C   Q(T) = SPL_INTERP(TQ_<species>,Q_<species>,D2,TLOG)
C
C Note that KP_Q_SPLN returns log10(Q) and log10(Kp)+D0*5040/T
C
C Reference:
C   Paul Barklem 2010, in preparation.
C
      SUBROUTINE KP_Q_SPLN(SPNAME,TEMP,Q_spln,K_spln,BARKLEM)
C
      IMPLICIT NONE
      CHARACTER SPNAME*(*)
      REAL TEMP
      LOGICAL BARKLEM
      REAL*8 Q_spln,K_spln
C
C  Local variables
C
      LOGICAL FIRST
      INTEGER MSPEC,NTQ,NTK,KLO,KHI,I,II,ISPEC
      PARAMETER(MSPEC=60, NTQ=33, NTK=46)
      INTEGER MTQ(MSPEC),MTK(MSPEC)
      REAL*8 TLOG,A,U(46),SPL_INTERP
C
      CHARACTER SPLIST(MSPEC)*8
      REAL*8 TQ(NTQ,MSPEC),Q(NTQ,MSPEC),Q2(NTQ,MSPEC)
      REAL*8 TK(NTK,MSPEC),K(NTK,MSPEC),K2(NTK,MSPEC)
      REAL*8         TQ_H2p  (NTQ),TQ_H2   (NTQ),TQ_H2m  (NTQ),
     * TQ_CH   (NTQ),TQ_CHm  (NTQ),TQ_C2   (NTQ),TQ_C2m  (NTQ),
     * TQ_CN   (NTQ),TQ_CNm  (NTQ),TQ_NH   (NTQ),TQ_N2   (NTQ),
     * TQ_OH   (NTQ),TQ_OHm  (NTQ),TQ_BO   (NTQ),TQ_CO   (NTQ),
     * TQ_NOp  (NTQ),TQ_NO   (NTQ),TQ_O2   (NTQ),TQ_HF   (NTQ),
     * TQ_NaH  (NTQ),TQ_MgH  (NTQ),TQ_MgO  (NTQ),TQ_AlH  (NTQ),
     * TQ_AlO  (NTQ),TQ_AlF  (NTQ),TQ_Al2  (NTQ),TQ_SiH  (NTQ),
     * TQ_SiHm (NTQ),TQ_SiC  (NTQ),TQ_SiN  (NTQ),TQ_SiO  (NTQ),
     * TQ_SiF  (NTQ),TQ_Si2  (NTQ),TQ_HS   (NTQ),TQ_HSm  (NTQ),
     * TQ_CS   (NTQ),TQ_NS   (NTQ),TQ_SO   (NTQ),TQ_MgS  (NTQ),
     * TQ_AlS  (NTQ),TQ_SiS  (NTQ),TQ_S2   (NTQ),TQ_HCl  (NTQ),
     * TQ_LiCl (NTQ),TQ_NaCl (NTQ),TQ_AlCl (NTQ),TQ_CaH  (NTQ),
     * TQ_CaF  (NTQ),TQ_CaCl (NTQ),TQ_ScO  (NTQ),TQ_TiO  (NTQ),
     * TQ_TiS  (NTQ),TQ_VO   (NTQ),TQ_CrH  (NTQ),TQ_CrO  (NTQ),
     * TQ_FeH  (NTQ),TQ_FeO  (NTQ),TQ_YO   (NTQ),TQ_ZrO  (NTQ),
     * TQ_LaO  (NTQ)
      REAL*8          Q_H2p  (NTQ), Q_H2   (NTQ), Q_H2m  (NTQ),
     *  Q_CH   (NTQ), Q_CHm  (NTQ), Q_C2   (NTQ), Q_C2m  (NTQ),
     *  Q_CN   (NTQ), Q_CNm  (NTQ), Q_NH   (NTQ), Q_N2   (NTQ),
     *  Q_OH   (NTQ), Q_OHm  (NTQ), Q_BO   (NTQ), Q_CO   (NTQ),
     *  Q_NOp  (NTQ), Q_NO   (NTQ), Q_O2   (NTQ), Q_HF   (NTQ),
     *  Q_NaH  (NTQ), Q_MgH  (NTQ), Q_MgO  (NTQ), Q_AlH  (NTQ),
     *  Q_AlO  (NTQ), Q_AlF  (NTQ), Q_Al2  (NTQ), Q_SiH  (NTQ),
     *  Q_SiHm (NTQ), Q_SiC  (NTQ), Q_SiN  (NTQ), Q_SiO  (NTQ),
     *  Q_SiF  (NTQ), Q_Si2  (NTQ), Q_HS   (NTQ), Q_HSm  (NTQ),
     *  Q_CS   (NTQ), Q_NS   (NTQ), Q_SO   (NTQ), Q_MgS  (NTQ),
     *  Q_AlS  (NTQ), Q_SiS  (NTQ), Q_S2   (NTQ), Q_HCl  (NTQ),
     *  Q_LiCl (NTQ), Q_NaCl (NTQ), Q_AlCl (NTQ), Q_CaH  (NTQ),
     *  Q_CaF  (NTQ), Q_CaCl (NTQ), Q_ScO  (NTQ), Q_TiO  (NTQ),
     *  Q_TiS  (NTQ), Q_VO   (NTQ), Q_CrH  (NTQ), Q_CrO  (NTQ),
     *  Q_FeH  (NTQ), Q_FeO  (NTQ), Q_YO   (NTQ), Q_ZrO  (NTQ),
     *  Q_LaO  (NTQ)
      REAL*8         TK_H2p  (NTK),TK_H2   (NTK),TK_H2m  (NTK),
     * TK_CH   (NTK),TK_CHm  (NTK),TK_C2   (NTK),TK_C2m  (NTK),
     * TK_CN   (NTK),TK_CNm  (NTK),TK_NH   (NTK),TK_N2   (NTK),
     * TK_OH   (NTK),TK_OHm  (NTK),TK_BO   (NTK),TK_CO   (NTK),
     * TK_NOp  (NTK),TK_NO   (NTK),TK_O2   (NTK),TK_HF   (NTK),
     * TK_NaH  (NTK),TK_MgH  (NTK),TK_MgO  (NTK),TK_AlH  (NTK),
     * TK_AlO  (NTK),TK_AlF  (NTK),TK_Al2  (NTK),TK_SiH  (NTK),
     * TK_SiHm (NTK),TK_SiC  (NTK),TK_SiN  (NTK),TK_SiO  (NTK),
     * TK_SiF  (NTK),TK_Si2  (NTK),TK_HS   (NTK),TK_HSm  (NTK),
     * TK_CS   (NTK),TK_NS   (NTK),TK_SO   (NTK),TK_MgS  (NTK),
     * TK_AlS  (NTK),TK_SiS  (NTK),TK_S2   (NTK),TK_HCl  (NTK),
     * TK_LiCl (NTK),TK_NaCl (NTK),TK_AlCl (NTK),TK_CaH  (NTK),
     * TK_CaF  (NTK),TK_CaCl (NTK),TK_ScO  (NTK),TK_TiO  (NTK),
     * TK_TiS  (NTK),TK_VO   (NTK),TK_CrH  (NTK),TK_CrO  (NTK),
     * TK_FeH  (NTK),TK_FeO  (NTK),TK_YO   (NTK),TK_ZrO  (NTK),
     * TK_LaO  (NTK)
      REAL*8          K_H2p  (NTK), K_H2   (NTK), K_H2m  (NTK),
     *  K_CH   (NTK), K_CHm  (NTK), K_C2   (NTK), K_C2m  (NTK),
     *  K_CN   (NTK), K_CNm  (NTK), K_NH   (NTK), K_N2   (NTK),
     *  K_OH   (NTK), K_OHm  (NTK), K_BO   (NTK), K_CO   (NTK),
     *  K_NOp  (NTK), K_NO   (NTK), K_O2   (NTK), K_HF   (NTK),
     *  K_NaH  (NTK), K_MgH  (NTK), K_MgO  (NTK), K_AlH  (NTK),
     *  K_AlO  (NTK), K_AlF  (NTK), K_Al2  (NTK), K_SiH  (NTK),
     *  K_SiHm (NTK), K_SiC  (NTK), K_SiN  (NTK), K_SiO  (NTK),
     *  K_SiF  (NTK), K_Si2  (NTK), K_HS   (NTK), K_HSm  (NTK),
     *  K_CS   (NTK), K_NS   (NTK), K_SO   (NTK), K_MgS  (NTK),
     *  K_AlS  (NTK), K_SiS  (NTK), K_S2   (NTK), K_HCl  (NTK),
     *  K_LiCl (NTK), K_NaCl (NTK), K_AlCl (NTK), K_CaH  (NTK),
     *  K_CaF  (NTK), K_CaCl (NTK), K_ScO  (NTK), K_TiO  (NTK),
     *  K_TiS  (NTK), K_VO   (NTK), K_CrH  (NTK), K_CrO  (NTK),
     *  K_FeH  (NTK), K_FeO  (NTK), K_YO   (NTK), K_ZrO  (NTK),
     *  K_LaO  (NTK)
      EQUIVALENCE (TQ(1, 1),TQ_H2p  ),(TQ(1, 2),TQ_H2   )
      EQUIVALENCE (TQ(1, 3),TQ_H2m  ),(TQ(1, 4),TQ_CH   )
      EQUIVALENCE (TQ(1, 5),TQ_CHm  ),(TQ(1, 6),TQ_C2   )
      EQUIVALENCE (TQ(1, 7),TQ_C2m  ),(TQ(1, 8),TQ_CN   )
      EQUIVALENCE (TQ(1, 9),TQ_CNm  ),(TQ(1,10),TQ_NH   )
      EQUIVALENCE (TQ(1,11),TQ_N2   ),(TQ(1,12),TQ_OH   )
      EQUIVALENCE (TQ(1,13),TQ_OHm  ),(TQ(1,14),TQ_BO   )
      EQUIVALENCE (TQ(1,15),TQ_CO   ),(TQ(1,16),TQ_NOp  )
      EQUIVALENCE (TQ(1,17),TQ_NO   ),(TQ(1,18),TQ_O2   )
      EQUIVALENCE (TQ(1,19),TQ_HF   ),(TQ(1,20),TQ_NaH  )
      EQUIVALENCE (TQ(1,21),TQ_MgH  ),(TQ(1,22),TQ_MgO  )
      EQUIVALENCE (TQ(1,23),TQ_AlH  ),(TQ(1,24),TQ_AlO  )
      EQUIVALENCE (TQ(1,25),TQ_AlF  ),(TQ(1,26),TQ_Al2  )
      EQUIVALENCE (TQ(1,27),TQ_SiH  ),(TQ(1,28),TQ_SiHm )
      EQUIVALENCE (TQ(1,29),TQ_SiC  ),(TQ(1,30),TQ_SiN  )
      EQUIVALENCE (TQ(1,31),TQ_SiO  ),(TQ(1,32),TQ_SiF  )
      EQUIVALENCE (TQ(1,33),TQ_Si2  ),(TQ(1,34),TQ_HS   )
      EQUIVALENCE (TQ(1,35),TQ_HSm  ),(TQ(1,36),TQ_CS   )
      EQUIVALENCE (TQ(1,37),TQ_NS   ),(TQ(1,38),TQ_SO   )
      EQUIVALENCE (TQ(1,39),TQ_MgS  ),(TQ(1,40),TQ_AlS  )
      EQUIVALENCE (TQ(1,41),TQ_SiS  ),(TQ(1,42),TQ_S2   )
      EQUIVALENCE (TQ(1,43),TQ_HCl  ),(TQ(1,44),TQ_LiCl )
      EQUIVALENCE (TQ(1,45),TQ_NaCl ),(TQ(1,46),TQ_AlCl )
      EQUIVALENCE (TQ(1,47),TQ_CaH  ),(TQ(1,48),TQ_CaF  )
      EQUIVALENCE (TQ(1,49),TQ_CaCl ),(TQ(1,50),TQ_ScO  )
      EQUIVALENCE (TQ(1,51),TQ_TiO  ),(TQ(1,52),TQ_TiS  )
      EQUIVALENCE (TQ(1,53),TQ_VO   ),(TQ(1,54),TQ_CrH  )
      EQUIVALENCE (TQ(1,55),TQ_CrO  ),(TQ(1,56),TQ_FeH  )
      EQUIVALENCE (TQ(1,57),TQ_FeO  ),(TQ(1,58),TQ_YO   )
      EQUIVALENCE (TQ(1,59),TQ_ZrO  ),(TQ(1,60),TQ_LaO  )
      EQUIVALENCE ( Q(1, 1), Q_H2p  ),( Q(1, 2), Q_H2   )
      EQUIVALENCE ( Q(1, 3), Q_H2m  ),( Q(1, 4), Q_CH   )
      EQUIVALENCE ( Q(1, 5), Q_CHm  ),( Q(1, 6), Q_C2   )
      EQUIVALENCE ( Q(1, 7), Q_C2m  ),( Q(1, 8), Q_CN   )
      EQUIVALENCE ( Q(1, 9), Q_CNm  ),( Q(1,10), Q_NH   )
      EQUIVALENCE ( Q(1,11), Q_N2   ),( Q(1,12), Q_OH   )
      EQUIVALENCE ( Q(1,13), Q_OHm  ),( Q(1,14), Q_BO   )
      EQUIVALENCE ( Q(1,15), Q_CO   ),( Q(1,16), Q_NOp  )
      EQUIVALENCE ( Q(1,17), Q_NO   ),( Q(1,18), Q_O2   )
      EQUIVALENCE ( Q(1,19), Q_HF   ),( Q(1,20), Q_NaH  )
      EQUIVALENCE ( Q(1,21), Q_MgH  ),( Q(1,22), Q_MgO  )
      EQUIVALENCE ( Q(1,23), Q_AlH  ),( Q(1,24), Q_AlO  )
      EQUIVALENCE ( Q(1,25), Q_AlF  ),( Q(1,26), Q_Al2  )
      EQUIVALENCE ( Q(1,27), Q_SiH  ),( Q(1,28), Q_SiHm )
      EQUIVALENCE ( Q(1,29), Q_SiC  ),( Q(1,30), Q_SiN  )
      EQUIVALENCE ( Q(1,31), Q_SiO  ),( Q(1,32), Q_SiF  )
      EQUIVALENCE ( Q(1,33), Q_Si2  ),( Q(1,34), Q_HS   )
      EQUIVALENCE ( Q(1,35), Q_HSm  ),( Q(1,36), Q_CS   )
      EQUIVALENCE ( Q(1,37), Q_NS   ),( Q(1,38), Q_SO   )
      EQUIVALENCE ( Q(1,39), Q_MgS  ),( Q(1,40), Q_AlS  )
      EQUIVALENCE ( Q(1,41), Q_SiS  ),( Q(1,42), Q_S2   )
      EQUIVALENCE ( Q(1,43), Q_HCl  ),( Q(1,44), Q_LiCl )
      EQUIVALENCE ( Q(1,45), Q_NaCl ),( Q(1,46), Q_AlCl )
      EQUIVALENCE ( Q(1,47), Q_CaH  ),( Q(1,48), Q_CaF  )
      EQUIVALENCE ( Q(1,49), Q_CaCl ),( Q(1,50), Q_ScO  )
      EQUIVALENCE ( Q(1,51), Q_TiO  ),( Q(1,52), Q_TiS  )
      EQUIVALENCE ( Q(1,53), Q_VO   ),( Q(1,54), Q_CrH  )
      EQUIVALENCE ( Q(1,55), Q_CrO  ),( Q(1,56), Q_FeH  )
      EQUIVALENCE ( Q(1,57), Q_FeO  ),( Q(1,58), Q_YO   )
      EQUIVALENCE ( Q(1,59), Q_ZrO  ),( Q(1,60), Q_LaO  )
      EQUIVALENCE (TK(1, 1),TK_H2p  ),(TK(1, 2),TK_H2   )
      EQUIVALENCE (TK(1, 3),TK_H2m  ),(TK(1, 4),TK_CH   )
      EQUIVALENCE (TK(1, 5),TK_CHm  ),(TK(1, 6),TK_C2   )
      EQUIVALENCE (TK(1, 7),TK_C2m  ),(TK(1, 8),TK_CN   )
      EQUIVALENCE (TK(1, 9),TK_CNm  ),(TK(1,10),TK_NH   )
      EQUIVALENCE (TK(1,11),TK_N2   ),(TK(1,12),TK_OH   )
      EQUIVALENCE (TK(1,13),TK_OHm  ),(TK(1,14),TK_BO   )
      EQUIVALENCE (TK(1,15),TK_CO   ),(TK(1,16),TK_NOp  )
      EQUIVALENCE (TK(1,17),TK_NO   ),(TK(1,18),TK_O2   )
      EQUIVALENCE (TK(1,19),TK_HF   ),(TK(1,20),TK_NaH  )
      EQUIVALENCE (TK(1,21),TK_MgH  ),(TK(1,22),TK_MgO  )
      EQUIVALENCE (TK(1,23),TK_AlH  ),(TK(1,24),TK_AlO  )
      EQUIVALENCE (TK(1,25),TK_AlF  ),(TK(1,26),TK_Al2  )
      EQUIVALENCE (TK(1,27),TK_SiH  ),(TK(1,28),TK_SiHm )
      EQUIVALENCE (TK(1,29),TK_SiC  ),(TK(1,30),TK_SiN  )
      EQUIVALENCE (TK(1,31),TK_SiO  ),(TK(1,32),TK_SiF  )
      EQUIVALENCE (TK(1,33),TK_Si2  ),(TK(1,34),TK_HS   )
      EQUIVALENCE (TK(1,35),TK_HSm  ),(TK(1,36),TK_CS   )
      EQUIVALENCE (TK(1,37),TK_NS   ),(TK(1,38),TK_SO   )
      EQUIVALENCE (TK(1,39),TK_MgS  ),(TK(1,40),TK_AlS  )
      EQUIVALENCE (TK(1,41),TK_SiS  ),(TK(1,42),TK_S2   )
      EQUIVALENCE (TK(1,43),TK_HCl  ),(TK(1,44),TK_LiCl )
      EQUIVALENCE (TK(1,45),TK_NaCl ),(TK(1,46),TK_AlCl )
      EQUIVALENCE (TK(1,47),TK_CaH  ),(TK(1,48),TK_CaF  )
      EQUIVALENCE (TK(1,49),TK_CaCl ),(TK(1,50),TK_ScO  )
      EQUIVALENCE (TK(1,51),TK_TiO  ),(TK(1,52),TK_TiS  )
      EQUIVALENCE (TK(1,53),TK_VO   ),(TK(1,54),TK_CrH  )
      EQUIVALENCE (TK(1,55),TK_CrO  ),(TK(1,56),TK_FeH  )
      EQUIVALENCE (TK(1,57),TK_FeO  ),(TK(1,58),TK_YO   )
      EQUIVALENCE (TK(1,59),TK_ZrO  ),(TK(1,60),TK_LaO  )
      EQUIVALENCE ( K(1, 1), K_H2p  ),( K(1, 2), K_H2   )
      EQUIVALENCE ( K(1, 3), K_H2m  ),( K(1, 4), K_CH   )
      EQUIVALENCE ( K(1, 5), K_CHm  ),( K(1, 6), K_C2   )
      EQUIVALENCE ( K(1, 7), K_C2m  ),( K(1, 8), K_CN   )
      EQUIVALENCE ( K(1, 9), K_CNm  ),( K(1,10), K_NH   )
      EQUIVALENCE ( K(1,11), K_N2   ),( K(1,12), K_OH   )
      EQUIVALENCE ( K(1,13), K_OHm  ),( K(1,14), K_BO   )
      EQUIVALENCE ( K(1,15), K_CO   ),( K(1,16), K_NOp  )
      EQUIVALENCE ( K(1,17), K_NO   ),( K(1,18), K_O2   )
      EQUIVALENCE ( K(1,19), K_HF   ),( K(1,20), K_NaH  )
      EQUIVALENCE ( K(1,21), K_MgH  ),( K(1,22), K_MgO  )
      EQUIVALENCE ( K(1,23), K_AlH  ),( K(1,24), K_AlO  )
      EQUIVALENCE ( K(1,25), K_AlF  ),( K(1,26), K_Al2  )
      EQUIVALENCE ( K(1,27), K_SiH  ),( K(1,28), K_SiHm )
      EQUIVALENCE ( K(1,29), K_SiC  ),( K(1,30), K_SiN  )
      EQUIVALENCE ( K(1,31), K_SiO  ),( K(1,32), K_SiF  )
      EQUIVALENCE ( K(1,33), K_Si2  ),( K(1,34), K_HS   )
      EQUIVALENCE ( K(1,35), K_HSm  ),( K(1,36), K_CS   )
      EQUIVALENCE ( K(1,37), K_NS   ),( K(1,38), K_SO   )
      EQUIVALENCE ( K(1,39), K_MgS  ),( K(1,40), K_AlS  )
      EQUIVALENCE ( K(1,41), K_SiS  ),( K(1,42), K_S2   )
      EQUIVALENCE ( K(1,43), K_HCl  ),( K(1,44), K_LiCl )
      EQUIVALENCE ( K(1,45), K_NaCl ),( K(1,46), K_AlCl )
      EQUIVALENCE ( K(1,47), K_CaH  ),( K(1,48), K_CaF  )
      EQUIVALENCE ( K(1,49), K_CaCl ),( K(1,50), K_ScO  )
      EQUIVALENCE ( K(1,51), K_TiO  ),( K(1,52), K_TiS  )
      EQUIVALENCE ( K(1,53), K_VO   ),( K(1,54), K_CrH  )
      EQUIVALENCE ( K(1,55), K_CrO  ),( K(1,56), K_FeH  )
      EQUIVALENCE ( K(1,57), K_FeO  ),( K(1,58), K_YO   )
      EQUIVALENCE ( K(1,59), K_ZrO  ),( K(1,60), K_LaO  )
C
      SAVE
C
      DATA SPLIST/
     * 'H2+   ','H2    ','H2-   ','CH    ','CH-   ','C2    ',
     * 'C2-   ','CN    ','CN-   ','NH    ','N2    ','OH    ',
     * 'OH-   ','BO    ','CO    ','NO+   ','NO    ','O2    ',
     * 'HF    ','NaH   ','MgH   ','MgO   ','AlH   ','AlO   ',
     * 'AlF   ','Al2   ','SiH   ','SiH-  ','SiC   ','SiN   ',
     * 'SiO   ','SiF   ','Si2   ','HS    ','HS-   ','CS    ',
     * 'NS    ','SO    ','MgS   ','AlS   ','SiS   ','S2    ',
     * 'HCl   ','LiCl  ','NaCl  ','AlCl  ','CaH   ','CaF   ',
     * 'CaCl  ','ScO   ','TiO   ','TiS   ','VO    ','CrH   ',
     * 'CrO   ','FeH   ','FeO   ','YO    ','ZrO   ','LaO   '/
C
C Molecular partition functions
C
      DATA TQ_H2p/                                                      071215
     1 -0.999999993529,-0.213099990781, 0.362900016751, 0.863300018968, H2p
     2  0.980699997741, 1.099999910318, 1.348399912911, 1.488499999486, H2p
     3  1.632000106900, 1.781399972218, 1.946900055894, 2.195800048889, H2p
     4  2.459400026935, 2.637100099165, 2.795999962688, 3.017700048765, H2p
     5  3.238299818206, 3.371799936189, 3.488399799847, 3.675099932070, H2p
     6  3.782100153231, 3.915099975260, 3.966599864718, 4.000000000000, H2p
     7       9*0.0D+00/                                                 H2p
      DATA  Q_H2p/                                                      071215
     1  0.00000000D+00, 3.13405372D-54, 2.79094653D-21, 1.23472801D-05, H2p
     2  1.91137670D-04, 1.59015822D-03, 2.86372258D-02, 7.72825554D-02, H2p
     3  1.55973372D-01, 2.60173878D-01, 3.92258869D-01, 6.10976857D-01, H2p
     4  8.58280722D-01, 1.03103181D+00, 1.19086826D+00, 1.43487833D+00, H2p
     5  1.72600071D+00, 1.93618282D+00, 2.14466300D+00, 2.51641126D+00, H2p
     6  2.73576648D+00, 3.00087246D+00, 3.09979361D+00, 3.16267828D+00, H2p
     7       9*0.0D+00/                                                 H2p
      DATA TQ_H2/                                                       071215
     1 -0.999999993529, 0.546399960943, 0.928499940989, 1.176300059693, H2
     2  1.287899949758, 1.399100094267, 1.648700096928, 1.743600009207, H2
     3  1.840399908833, 2.009600002847, 2.197800044999, 2.438100070396, H2
     4  2.693900045012, 2.912900087310, 3.025900111066, 3.127900068141, H2
     5  3.483099942820, 3.602599959587, 3.722599865216, 3.880100116638, H2
     6  3.949199845377, 4.000000000000,     11*0.0D+00/                 H2
      DATA  Q_H2/                                                       071215
     1 -3.01029996D-01,-3.01029996D-01,-3.01029993D-01,-3.01014761D-01, H2
     2 -3.00830296D-01,-2.99579846D-01,-2.73523089D-01,-2.44514974D-01, H2
     3 -2.00940921D-01,-9.66513596D-02, 4.66537413D-02, 2.53454444D-01, H2
     4  4.90680129D-01, 7.02774809D-01, 8.15787556D-01, 9.21650976D-01, H2
     5  1.35646140D+00, 1.54038533D+00, 1.74931286D+00, 2.05314222D+00, H2
     6  2.19184017D+00, 2.29424176D+00,     11*0.0D+00/                 H2
      DATA TQ_H2m/                                                      071215
     1 -0.999999993529, 0.546399960943, 0.928499940989, 1.176300059693, H2m
     2  1.287899949758, 1.399100094267, 1.648700096928, 1.743600009207, H2m
     3  1.840399908833, 2.009600002847, 2.197800044999, 2.438100070396, H2m
     4  2.693900045012, 2.912900087310, 3.025900111066, 3.127900068141, H2m
     5  3.483099942820, 3.602599959587, 3.722499862962, 3.880300111551, H2m
     6  3.949599836878, 4.000000000000,     11*0.0D+00/                 H2m
      DATA  Q_H2m/                                                      071215
     1  0.00000000D+00,-2.33773993D-28, 2.44687637D-09, 1.52343412D-05, H2m
     2  1.99699327D-04, 1.45015002D-03, 2.75069068D-02, 5.65150214D-02, H2m
     3  1.00089074D-01, 2.04378636D-01, 3.47683737D-01, 5.54484440D-01, H2m
     4  7.91710125D-01, 1.00380481D+00, 1.11681755D+00, 1.22268097D+00, H2m
     5  1.65749140D+00, 1.84141532D+00, 2.05015916D+00, 2.35457108D+00, H2m
     6  2.49367435D+00, 2.59526213D+00,     11*0.0D+00/                 H2m
      DATA TQ_CH/                                                       071215
     1 -0.999999993529,-0.326599980548, 0.166399959291, 0.430299957164, CH
     2  0.636699985352, 0.739300022626, 0.832600038271, 0.960700052763, CH
     3  1.065399944135, 1.210300015201, 1.390500091928, 1.596199899459, CH
     4  1.839399909152, 2.126100009214, 2.433400079679, 2.656200091992, CH
     5  2.850899897814, 2.995100004942, 3.200799816137, 3.360900121256, CH
     6  3.582299979332, 3.733499938459, 3.807899997628, 3.881700075939, CH
     7  3.954999937558, 4.000000000000,      7*0.0D+00/                 CH
      DATA  Q_CH/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02060332D-01, CH
     2  6.02165251D-01, 6.02823947D-01, 6.05270405D-01, 6.16741210D-01, CH
     3  6.39356745D-01, 6.97275561D-01, 8.07373382D-01, 9.64829191D-01, CH
     4  1.17499627D+00, 1.44096590D+00, 1.73756867D+00, 1.95686009D+00, CH
     5  2.15204558D+00, 2.30321543D+00, 2.54391678D+00, 2.76626113D+00, CH
     6  3.14285914D+00, 3.45064753D+00, 3.61820766D+00, 3.79479250D+00, CH
     7  3.97873352D+00, 4.09454845D+00,      7*0.0D+00/                 CH
      DATA TQ_CHm/                                                      071215
     1 -0.999999993529,-0.326599980548, 0.166399959291, 0.430299957164, CHm
     2  0.636699985352, 0.739300022626, 0.832600038271, 0.960700052763, CHm
     3  1.065399944135, 1.210300015201, 1.390500091928, 1.596199899459, CHm
     4  1.839399909152, 2.126100009214, 2.433400079679, 2.656200091992, CHm
     5  2.850899897814, 2.995100004942, 3.200799816137, 3.360900121256, CHm
     6  3.582299979332, 3.733499938459, 3.807899997628, 3.881700075939, CHm
     7  3.954999937558, 4.000000000000,      7*0.0D+00/                 CHm
      DATA  Q_CHm/                                                      071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02060332D-01, CHm
     2  6.02165251D-01, 6.02823947D-01, 6.05270405D-01, 6.16741210D-01, CHm
     3  6.39356745D-01, 6.97275561D-01, 8.07373382D-01, 9.64829191D-01, CHm
     4  1.17499627D+00, 1.44096590D+00, 1.73756867D+00, 1.95686009D+00, CHm
     5  2.15204558D+00, 2.30321543D+00, 2.54391678D+00, 2.76626113D+00, CHm
     6  3.14285914D+00, 3.45064753D+00, 3.61820766D+00, 3.79479250D+00, CHm
     7  3.97873352D+00, 4.09454845D+00,      7*0.0D+00/                 CHm
      DATA TQ_C2/                                                       071215
     1 -0.999999993529,-0.641099986909,-0.379799981428,-0.243700001436, C2
     2 -0.118699997287, 0.073100042161, 0.177699935663, 0.322700045690, C2
     3  0.524599971267, 0.758799993564, 1.038999957821, 1.324899929720, C2
     4  1.835299910833, 1.988900009651, 2.066799937361, 2.141700122459, C2
     5  2.289299949530, 2.361199903842, 2.431800082839, 2.668300083763, C2
     6  2.782399986099, 2.898600107396, 3.120299878697, 3.313800022675, C2
     7  3.455399924176, 3.584100021504, 3.731799983892, 3.848299933474, C2
     8  3.940100038723, 3.976599949064, 4.000000000000,      2*0.0D+00/ C2
      DATA  Q_C2/                                                       071215
     1 -3.01029996D-01,-3.01029996D-01,-3.01025114D-01,-3.00889289D-01, C2
     2 -2.99650416D-01,-2.85375904D-01,-2.61880905D-01,-2.02557475D-01, C2
     3 -7.58154506D-02, 1.08741954D-01, 3.56855782D-01, 6.25923966D-01, C2
     4  1.12398377D+00, 1.27632444D+00, 1.35481857D+00, 1.43279873D+00, C2
     5  1.60543007D+00, 1.70561849D+00, 1.81703868D+00, 2.25976911D+00, C2
     6  2.48518357D+00, 2.71100862D+00, 3.13226319D+00, 3.50279565D+00, C2
     7  3.78546485D+00, 4.05600414D+00, 4.38531038D+00, 4.66021357D+00, C2
     8  4.88728164D+00, 4.98046685D+00, 5.04115178D+00,      2*0.0D+00/ C2
      DATA TQ_C2m/                                                      071215
     1 -0.999999993529,-0.646599981462,-0.389999994590,-0.260699995572, C2m
     2 -0.141499983639, 0.058599929497, 0.174699940848, 0.310299921211, C2m
     3  0.425999947878, 0.553899952637, 0.721100033889, 0.895099981695, C2m
     4  1.195400039464, 1.464700029327, 1.775999970693, 2.102699907919, C2m
     5  2.343499919663, 2.547699957449, 2.710900053613, 2.955600044242, C2m
     6  3.121799916088, 3.299600180459, 3.422400027540, 3.539899959774, C2m
     7  3.651900050636, 3.753200063022, 3.826700005825, 3.895699998977, C2m
     8  3.960500032946, 3.984899911457, 4.000000000000,      2*0.0D+00/ C2m
      DATA  Q_C2m/                                                      071215
     1  0.00000000D+00, 3.06186733D-10, 6.05130452D-06, 1.42892529D-04, C2m
     2  1.27320127D-03, 1.61646004D-02, 4.37243232D-02, 1.01235802D-01, C2m
     3  1.69606083D-01, 2.59829165D-01, 3.93016002D-01, 5.43581825D-01, C2m
     4  8.20473194D-01, 1.07902941D+00, 1.38398722D+00, 1.70756358D+00, C2m
     5  1.94727608D+00, 2.15137195D+00, 2.31733908D+00, 2.58660581D+00, C2m
     6  2.79634658D+00, 3.05059595D+00, 3.24376566D+00, 3.44147908D+00, C2m
     7  3.64408082D+00, 3.84616767D+00, 4.00941275D+00, 4.17762667D+00, C2m
     8  4.34796842D+00, 4.41465569D+00, 4.45648843D+00,      2*0.0D+00/ C2m
      DATA TQ_CN/                                                       071215
     1 -0.999999993529,-0.634899997994,-0.368900014961,-0.234400007277, CN
     2 -0.110400017216,-0.010199997823, 0.096500013917, 0.212500027912, CN
     3  0.347500041696, 0.460800038277, 0.586600034178, 0.748799991312, CN
     4  0.917799942272, 1.222400041473, 1.512099970292, 1.816699942154, CN
     5  2.171700062029, 2.515399988212, 2.651500090913, 2.785999982288, CN
     6  3.005000111574, 3.192999969602, 3.342800094581, 3.533300119361, CN
     7  3.618899924457, 3.702099855522, 3.905599956847, 3.962699972274, CN
     8  4.000000000000,      4*0.0D+00/                                 CN
      DATA  Q_CN/                                                       071215
     1  3.01029996D-01, 3.01029996D-01, 3.01033877D-01, 3.01144985D-01, CN
     2  3.02197074D-01, 3.05962496D-01, 3.17414924D-01, 3.45164806D-01, CN
     3  4.02616590D-01, 4.69509483D-01, 5.57986956D-01, 6.86686930D-01, CN
     4  8.32311015D-01, 1.11249062D+00, 1.39061897D+00, 1.68912236D+00, CN
     5  2.04085086D+00, 2.38334891D+00, 2.51982417D+00, 2.65723188D+00, CN
     6  2.89757486D+00, 3.13417863D+00, 3.34956572D+00, 3.67021044D+00, CN
     7  3.83669040D+00, 4.01310850D+00, 4.49624705D+00, 4.64196948D+00, CN
     8  4.73930577D+00,      4*0.0D+00/                                 CN
      DATA TQ_CNm/                                                      071215
     1 -0.999999993529,-0.634899997994,-0.368900014961,-0.234000007400, CNm
     2 -0.109800018053,-0.009999998390, 0.096300014197, 0.212300027831, CNm
     3  0.348800043430, 0.463300031077, 0.582300046379, 0.753799989938, CNm
     4  0.907999968524, 1.213400026211, 1.503899984038, 1.813799942518, CNm
     5  2.165600070495, 2.511599997941, 2.647500092285, 2.781499987052, CNm
     6  3.004200093722, 3.193399959808, 3.343400079660, 3.531800155631, CNm
     7  3.616699974148, 3.699299826562, 3.906399936370, 3.962899966758, CNm
     8  4.000000000000,      4*0.0D+00/                                 CNm
      DATA  Q_CNm/                                                      071215
     1  3.01029996D-01, 3.01029996D-01, 3.01033877D-01, 3.01145970D-01, CNm
     2  3.02208363D-01, 3.05974907D-01, 3.17382190D-01, 3.45099132D-01, CNm
     3  4.03295361D-01, 4.71137779D-01, 5.54768737D-01, 6.90845559D-01, CNm
     4  8.23616989D-01, 1.10397357D+00, 1.38261707D+00, 1.68619289D+00, CNm
     5  2.03462898D+00, 2.37920185D+00, 2.51531135D+00, 2.65190139D+00, CNm
     6  2.89554030D+00, 3.13301339D+00, 3.34805469D+00, 3.66359914D+00, CNm
     7  3.82748956D+00, 4.00110767D+00, 4.48813539D+00, 4.63059377D+00, CNm
     8  4.72601398D+00,      4*0.0D+00/                                 CNm
      DATA TQ_NH/                                                       071215
     1 -0.999999993529, 0.209400027812, 0.636099987054, 0.750299987400, NH
     2  0.863200019064, 1.101399906583, 1.242799998652, 1.406300082852, NH
     3  1.588599913306, 1.819599941789, 2.123499954682, 2.476200015855, NH
     4  2.651300090867, 2.829899942405, 2.990500009582, 3.179199853002, NH
     5  3.425899932187, 3.596799983999, 3.773000058384, 3.841100116568, NH
     6  3.903999997802, 3.963099961242, 3.985699890337, 4.000000000000, NH
     7       9*0.0D+00/                                                 NH
      DATA  Q_NH/                                                       071215
     1  4.77121255D-01, 4.77121255D-01, 4.77145558D-01, 4.77422583D-01, NH
     2  4.79165128D-01, 5.07336232D-01, 5.57927870D-01, 6.50927763D-01, NH
     3  7.83772883D-01, 9.77128946D-01, 1.25476170D+00, 1.59370843D+00, NH
     4  1.76561270D+00, 1.94310238D+00, 2.10736564D+00, 2.31634551D+00, NH
     5  2.63911181D+00, 2.90779886D+00, 3.23761517D+00, 3.38408466D+00, NH
     6  3.53002527D+00, 3.67600650D+00, 3.73381804D+00, 3.77089299D+00, NH
     7       9*0.0D+00/                                                 NH
      DATA TQ_N2/                                                       071215
     1 -0.999999993529,-0.627400007110,-0.355500029093,-0.212899990105, N2
     2 -0.083999978765, 0.117200005025, 0.235600010631, 0.372299984859, N2
     3  0.534399962113, 0.755199990953, 1.067699940593, 1.374399976828, N2
     4  1.693000063852, 2.013199998928, 2.385000112558, 2.657300092245, N2
     5  2.846099898354, 3.048700154719, 3.190600028365, 3.327999962751, N2
     6  3.794200066108, 3.924899965861, 3.970499784444, 4.000000000000, N2
     7       9*0.0D+00/                                                 N2
      DATA  Q_N2/                                                       071215
     1 -3.01029996D-01,-3.01029996D-01,-3.01026992D-01,-3.00916378D-01, N2
     2 -2.99778665D-01,-2.84863174D-01,-2.56557042D-01,-1.98005654D-01, N2
     3 -9.78023806D-02, 7.16534937D-02, 3.45696360D-01, 6.34211414D-01, N2
     4  9.43681335D-01, 1.25953307D+00, 1.62914120D+00, 1.90110477D+00, N2
     5  2.09323219D+00, 2.31478263D+00, 2.48830497D+00, 2.67463022D+00, N2
     6  3.43476112D+00, 3.67792111D+00, 3.76647032D+00, 3.82533357D+00, N2
     7       9*0.0D+00/                                                 N2
      DATA TQ_OH/                                                       071215
     1 -0.999999993529, 0.243500002438, 0.709000048386, 0.816000020033, OH
     2  0.922399937580, 1.146000105069, 1.278199965050, 1.421700077174, OH
     3  1.606999894676, 1.835299910833, 2.139700127692, 2.483900004198, OH
     4  2.691000044681, 2.875800016063, 3.043500051181, 3.215999875797, OH
     5  3.431499854524, 3.616199985441, 3.782600142144, 3.909599854461, OH
     6  3.964999908843, 3.986199877136, 4.000000000000,     10*0.0D+00/ OH
      DATA  Q_OH/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02098506D-01, 6.02435415D-01, OH
     2  6.04262230D-01, 6.29944916D-01, 6.74423842D-01, 7.51069830D-01, OH
     3  8.80944146D-01, 1.06853042D+00, 1.34409440D+00, 1.67338502D+00, OH
     4  1.87630902D+00, 2.05975529D+00, 2.23094861D+00, 2.42043757D+00, OH
     5  2.69348987D+00, 2.96997378D+00, 3.25804067D+00, 3.50868908D+00, OH
     6  3.62747034D+00, 3.67436200D+00, 3.70528328D+00,     10*0.0D+00/ OH
      DATA TQ_OHm/                                                      071215
     1 -0.999999993529, 0.223000028156, 0.697700045307, 0.806999949139, OHm
     2  0.917699942546, 1.149900097428, 1.270999979838, 1.402000090810, OHm
     3  1.554699943369, 1.742500012992, 1.977100016998, 2.215100031905, OHm
     4  2.576599933971, 2.857399898411, 2.991600008473, 3.197599856972, OHm
     5  3.363100067584, 3.558299981945, 3.744699869290, 3.903100020839, OHm
     6  3.962099988821, 4.000000000000,     11*0.0D+00/                 OHm
      DATA  Q_OHm/                                                      071215
     1  0.00000000D+00,-2.21986382D-18, 2.38344820D-05, 2.69911385D-04, OHm
     2  1.81659858D-03, 2.68675942D-02, 6.57969495D-02, 1.32082718D-01, OHm
     3  2.33568782D-01, 3.80619521D-01, 5.84573697D-01, 8.04775559D-01, OHm
     4  1.15260926D+00, 1.42868887D+00, 1.56307779D+00, 1.78072967D+00, OHm
     5  1.97588845D+00, 2.23815179D+00, 2.52106025D+00, 2.78255794D+00, OHm
     6  2.88403464D+00, 2.95023047D+00,     11*0.0D+00/                 OHm
      DATA TQ_BO/                                                       071215
     1 -0.999999993529,-0.644399983641,-0.385699988646,-0.258099996706, BO
     2 -0.129900001370, 0.050099954400, 0.181399948365, 0.325500047784, BO
     3  0.526599970651, 0.675199946549, 0.817600038866, 1.120799900509, BO
     4  1.397500093832, 1.712900047246, 2.034199972127, 2.328799927826, BO
     5  2.562199940748, 2.755499995921, 2.961900036039, 3.122899943507, BO
     6  3.279300200142, 3.565800073537, 3.677399874723, 3.783700117752, BO
     7  3.917400034356, 3.967499839897, 3.986999856016, 4.000000000000, BO
     8       5*0.0D+00/                                                 BO
      DATA  Q_BO/                                                       071215
     1  3.01029996D-01, 3.01029996D-01, 3.01035335D-01, 3.01155585D-01, BO
     2  3.02363075D-01, 3.14616549D-01, 3.44112020D-01, 4.05681433D-01, BO
     3  5.34224748D-01, 6.49054853D-01, 7.68804586D-01, 1.04313679D+00, BO
     4  1.30653283D+00, 1.61428844D+00, 1.93190684D+00, 2.22495195D+00, BO
     5  2.45809826D+00, 2.65499084D+00, 2.88182773D+00, 3.08175651D+00, BO
     6  3.29999902D+00, 3.75895058D+00, 3.95664774D+00, 4.15699252D+00, BO
     7  4.43657573D+00, 4.55348068D+00, 4.60125225D+00, 4.63386799D+00, BO
     8       5*0.0D+00/                                                 BO
      DATA TQ_CO/                                                       071215
     1 -0.999999993529,-0.632600002502,-0.364600019334,-0.224500011640, CO
     2 -0.097100000602, 0.101700010930, 0.218900030490, 0.357400033132, CO
     3  0.528599970035, 0.761699992420, 1.076199926571, 1.381600097805, CO
     4  1.692000066444, 2.001400000415, 2.360199905649, 2.627600048384, CO
     5  2.818099938077, 3.018200037442, 3.159599841620, 3.302900117350, CO
     6  3.578199969770, 3.810899966839, 3.872799956303, 3.930999846621, CO
     7  3.972799846514, 3.988499816415, 4.000000000000,      6*0.0D+00/ CO
      DATA  Q_CO/                                                       071215
     1  0.00000000D+00, 6.35186438D-11, 3.56703693D-06, 1.21907440D-04, CO
     2  1.28744332D-03, 1.60687623D-02, 4.38850285D-02, 1.03016376D-01, CO
     3  2.09583107D-01, 3.90424533D-01, 6.68213646D-01, 9.56436228D-01, CO
     4  1.25829541D+00, 1.56354193D+00, 1.92010809D+00, 2.18702897D+00, CO
     5  2.38096809D+00, 2.60027262D+00, 2.77354484D+00, 2.96855587D+00, CO
     6  3.39710219D+00, 3.80567834D+00, 3.92064615D+00, 4.03186141D+00, CO
     7  4.11457737D+00, 4.14653977D+00, 4.17035582D+00,      6*0.0D+00/ CO
      DATA TQ_NOp/                                                      071215
     1 -0.999999993529,-0.627500007129,-0.355700028867,-0.213199991119, NOp
     2 -0.084299978459, 0.116800005873, 0.235200011420, 0.371899985211, NOp
     3  0.534299962283, 0.755599991243, 1.068799938899, 1.376000011720, NOp
     4  1.695400057631, 2.016199995143, 2.388700102186, 2.661600091110, NOp
     5  2.850199897750, 3.051800135900, 3.193599954911, 3.330299922252, NOp
     6  3.798200149970, 3.926799913398, 3.971199803335, 3.988099826975, NOp
     7  4.000000000000,      8*0.0D+00/                                 NOp
      DATA  Q_NOp/                                                      071215
     1  0.00000000D+00, 3.76087167D-11, 3.02047723D-06, 1.13817836D-04, NOp
     2  1.25296708D-03, 1.61639205D-02, 4.44670474D-02, 1.03014787D-01, NOp
     3  2.03423803D-01, 3.73328797D-01, 6.48085133D-01, 9.37142685D-01, NOp
     4  1.24743878D+00, 1.56391262D+00, 1.93423483D+00, 2.20681137D+00, NOp
     5  2.39879702D+00, 2.61938382D+00, 2.79288432D+00, 2.97835839D+00, NOp
     6  3.74262742D+00, 3.98291456D+00, 4.06981316D+00, 4.10371377D+00, NOp
     7  4.12793386D+00,      8*0.0D+00/                                 NOp
      DATA TQ_NO/                                                       071215
     1 -0.999999993529,-0.654599970871,-0.404399994113,-0.278299959398, NO
     2 -0.152499987861, 0.023999994474, 0.161199974196, 0.286499972876, NO
     3  0.456300041849, 0.584600039853, 0.805999949049, 1.002999991176, NO
     4  1.213400026211, 1.378800072782, 1.639900090060, 1.764899974394, NO
     5  1.942200044609, 2.159600088090, 2.290399947338, 2.422600091227, NO
     6  2.561999940497, 2.685000051840, 2.986300011188, 3.212199968587, NO
     7  3.554700075022, 3.749399973836, 3.905899949168, 3.963199958484, NO
     8  3.985699890337, 4.000000000000,      3*0.0D+00/                 NO
      DATA  Q_NO/                                                       071215
     1  3.01029996D-01, 3.01029996D-01, 3.01036917D-01, 3.01177682D-01, NO
     2  3.02478523D-01, 3.14867816D-01, 3.46604988D-01, 3.99867270D-01, NO
     3  5.04021357D-01, 5.99081584D-01, 7.82950452D-01, 9.59742050D-01, NO
     4  1.15675939D+00, 1.31569313D+00, 1.57798259D+00, 1.71426733D+00, NO
     5  1.92357187D+00, 2.19637389D+00, 2.36116182D+00, 2.52429714D+00, NO
     6  2.69139052D+00, 2.83547800D+00, 3.19632754D+00, 3.50287233D+00, NO
     7  4.05093712D+00, 4.40255045D+00, 4.70480156D+00, 4.82105609D+00, NO
     8  4.86785552D+00, 4.89800724D+00,      3*0.0D+00/                 NO
      DATA TQ_O2/                                                       071215
     1 -0.999999993529,-0.673499974515,-0.441399995599,-0.324799985013, O2
     2 -0.206899982021,-0.042800021354, 0.096600013777, 0.219200030611, O2
     3  0.365100007495, 0.541799954931, 0.739900021743, 0.945600031954, O2
     4  1.202600029941, 1.470500031909, 1.780299972159, 2.096799919076, O2
     5  2.309599942625, 2.511899997173, 2.662800089794, 2.903800101240, O2
     6  3.061499971042, 3.247299940725, 3.459500023320, 3.632700046534, O2
     7  3.742399818129, 3.821500130460, 3.899900097834, 3.960200041220, O2
     8  3.984799914097, 4.000000000000,      3*0.0D+00/                 O2
      DATA  Q_O2/                                                       071215
     1  1.76091259D-01, 1.76091264D-01, 1.76106327D-01, 1.76310225D-01, O2
     2  1.77817001D-01, 1.89780621D-01, 2.22047628D-01, 2.74131733D-01, O2
     3  3.61644055D-01, 4.92200725D-01, 6.57711464D-01, 8.42495750D-01, O2
     4  1.08423424D+00, 1.34349628D+00, 1.64818199D+00, 1.96221445D+00, O2
     5  2.17421133D+00, 2.37658113D+00, 2.53035277D+00, 2.79598012D+00, O2
     6  2.99460614D+00, 3.26037482D+00, 3.60856412D+00, 3.93248172D+00, O2
     7  4.15843502D+00, 4.33296246D+00, 4.51855468D+00, 4.67270253D+00, O2
     8  4.73899302D+00, 4.78101842D+00,      3*0.0D+00/                 O2
      DATA TQ_HF/                                                       071215
     1 -0.999999993529, 0.270699970872, 0.758099993056, 0.862600019641, HF
     2  0.967000028531, 1.185900048513, 1.312899948348, 1.447100043595, HF
     3  1.635400099653, 1.864199897812, 2.169300058961, 2.511599997941, HF
     4  2.727500009000, 2.918100078781, 3.083100005158, 3.254000089478, HF
     5  3.452099844378, 3.634499996198, 3.788900002446, 3.913999946997, HF
     6  3.966399870233, 4.000000000000,     11*0.0D+00/                 HF
      DATA  Q_HF/                                                       071215
     1  0.00000000D+00, 2.06736485D-18, 4.28450270D-05, 3.89455273D-04, HF
     2  2.20095798D-03, 2.67750796D-02, 6.80957645D-02, 1.37265950D-01, HF
     3  2.66508755D-01, 4.52871669D-01, 7.27979959D-01, 1.05485387D+00, HF
     4  1.26622692D+00, 1.45535170D+00, 1.62359433D+00, 1.81061591D+00, HF
     5  2.05825562D+00, 2.32532732D+00, 2.58470703D+00, 2.82129143D+00, HF
     6  2.92839792D+00, 2.99954398D+00,     11*0.0D+00/                 HF
      DATA TQ_NaH/                                                      071215
     1 -0.999999993529,-0.491600012464,-0.119499995366, 0.070100042786, NaH
     2  0.239200003531, 0.373099984156, 0.515999974697, 0.655199954116, NaH
     3  0.957700048717, 1.097299915751, 1.252700012270, 1.453300038138, NaH
     4  1.666200074527, 1.909900075822, 2.167100065819, 2.316599935755, NaH
     5  2.455200031937, 2.676400067316, 2.826999942210, 2.996500003530, NaH
     6  3.179699840953, 3.313200009900, 3.449399808825, 3.535900056494, NaH
     7  3.617899947044, 3.758000176780, 3.885199986910, 3.954999937558, NaH
     8  4.000000000000,      4*0.0D+00/                                 NaH
      DATA  Q_NaH/                                                      071215
     1  0.00000000D+00, 1.07495077D-25, 1.45302174D-08, 9.45191022D-06, NaH
     2  4.29372634D-04, 3.58966602D-03, 1.84078097D-02, 5.64791277D-02, NaH
     3  2.30140723D-01, 3.37362777D-01, 4.67655193D-01, 6.47028778D-01, NaH
     4  8.46150376D-01, 1.08075157D+00, 1.33296886D+00, 1.48116201D+00, NaH
     5  1.62042402D+00, 1.85515095D+00, 2.03393468D+00, 2.26230890D+00, NaH
     6  2.54524052D+00, 2.77585914D+00, 3.03493720D+00, 3.21401608D+00, NaH
     7  3.39445264D+00, 3.72231063D+00, 4.03661300D+00, 4.21589694D+00, NaH
     8  4.33405348D+00,      4*0.0D+00/                                 NaH
      DATA TQ_MgH/                                                      071215
     1 -0.999999993529,-0.465200007095,-0.073699988148, 0.127800004931, MgH
     2  0.304699941696, 0.444300051767, 0.593800025302, 0.743600009207, MgH
     3  0.977800001749, 1.097799914745, 1.230900013280, 1.465000029526, MgH
     4  1.681700078478, 1.943000046530, 2.224000012063, 2.393000101289, MgH
     5  2.547799957450, 2.778699987487, 2.933700064103, 3.086100063422, MgH
     6  3.248599970312, 3.400000091874, 3.563000008101, 3.691100047413, MgH
     7  3.869799897759, 3.947799875123, 4.000000000000,      6*0.0D+00/ MgH
      DATA  Q_MgH/                                                      071215
     1  3.01029996D-01, 3.01029996D-01, 3.01030000D-01, 3.01035989D-01, MgH
     2  3.01395752D-01, 3.04481095D-01, 3.20096833D-01, 3.63066270D-01, MgH
     3  4.93130449D-01, 5.80442612D-01, 6.87138870D-01, 8.91188943D-01, MgH
     4  1.09166220D+00, 1.34187422D+00, 1.61686699D+00, 1.78417903D+00, MgH
     5  1.93941452D+00, 2.18408938D+00, 2.36860969D+00, 2.57370643D+00, MgH
     6  2.82149418D+00, 3.08064348D+00, 3.39446338D+00, 3.67278254D+00, MgH
     7  4.11169175D+00, 4.31927706D+00, 4.46189498D+00,      6*0.0D+00/ MgH
      DATA TQ_MgO/                                                      071215
     1 -0.999999993529,-0.853800028128,-0.689500012187,-0.530399953027, MgO
     2 -0.436900003795,-0.337899962300,-0.214399995176,-0.096499998240, MgO
     3  0.039199965031, 0.179499932552, 0.506099988462, 0.964800036993, MgO
     4  1.463900028798, 1.810299942958, 2.080499933417, 2.200800041863, MgO
     5  2.323299930846, 2.480200008650, 2.612499896732, 2.726200012442, MgO
     6  2.845399898465, 2.947300048054, 3.014700116707, 3.081199968257, MgO
     7  3.166399987045, 3.263000147269, 3.373999981339, 3.484499905054, MgO
     8  3.580599939503, 3.675099932070, 3.832499981182, 3.931199850937, MgO
     9  4.000000000000/                                                 MgO
      DATA  Q_MgO/                                                      071215
     1  9.31473671D-08, 1.02718001D-05, 4.15810007D-04, 4.88595179D-03, MgO
     2  1.42479044D-02, 3.48458015D-02, 8.06391302D-02, 1.44547900D-01, MgO
     3  2.36281647D-01, 3.44726851D-01, 6.28348804D-01, 1.06253999D+00, MgO
     4  1.55282649D+00, 1.89705347D+00, 2.16656558D+00, 2.28706256D+00, MgO
     5  2.41130808D+00, 2.57718772D+00, 2.72926951D+00, 2.87606107D+00, MgO
     6  3.05890474D+00, 3.25273205D+00, 3.40409790D+00, 3.57051902D+00, MgO
     7  3.80256568D+00, 4.07948904D+00, 4.40057461D+00, 4.71304105D+00, MgO
     8  4.97587011D+00, 5.22646479D+00, 5.63223394D+00, 5.88454326D+00, MgO
     9  6.06134597D+00/                                                 MgO
      DATA TQ_AlH/                                                      071215
     1 -0.999999993529,-0.451099987798,-0.049300015167, 0.158499976883, AlH
     2  0.339700032276, 0.481900014868, 0.637399983365, 0.784099972212, AlH
     3  1.029999989560, 1.153100097909, 1.288499948878, 1.531899985759, AlH
     4  1.761699978080, 2.042499967572, 2.333099925531, 2.476800014755, AlH
     5  2.610999901243, 2.829799942398, 2.975900013946, 3.186599971883, AlH
     6  3.345300032409, 3.455699931431, 3.561399970709, 3.682499876682, AlH
     7  3.788800004663, 3.889799869901, 3.956399968127, 4.000000000000, AlH
     8       5*0.0D+00/                                                 AlH
      DATA  Q_AlH/                                                      071215
     1  0.00000000D+00, 2.88831324D-30, 1.99020487D-09, 4.48117923D-06, AlH
     2  3.27338663D-04, 3.30104116D-03, 1.95672632D-02, 6.18939551D-02, AlH
     3  1.99899233D-01, 2.90552542D-01, 4.00095984D-01, 6.13776272D-01, AlH
     4  8.27624272D-01, 1.09780952D+00, 1.38316996D+00, 1.52579320D+00, AlH
     5  1.66071207D+00, 1.89265000D+00, 2.06530475D+00, 2.35204779D+00, AlH
     6  2.60157298D+00, 2.79583974D+00, 3.00512505D+00, 3.28470793D+00, AlH
     7  3.56855504D+00, 3.86261803D+00, 4.06266455D+00, 4.19411233D+00, AlH
     8       5*0.0D+00/                                                 AlH
      DATA TQ_AlO/                                                      071215
     1 -0.999999993529,-0.696999997813,-0.602800014125,-0.514599984808, AlO
     2 -0.322399990965,-0.199499985480,-0.075899986296, 0.031899976827, AlO
     3  0.145199977027, 0.411999945400, 0.744600005766, 1.097899914543, AlO
     4  1.449400043923, 1.948700060217, 2.115099882434, 2.271099990420, AlO
     5  2.418000095986, 2.559699938576, 2.777399986334, 3.015000109913, AlO
     6  3.134899988137, 3.323300075511, 3.463199955186, 3.580999948875, AlO
     7  3.693599980080, 3.890599878936, 3.955899957210, 4.000000000000, AlO
     8       5*0.0D+00/                                                 AlO
      DATA  Q_AlO/                                                      071215
     1  3.01030009D-01, 3.01169150D-01, 3.01857051D-01, 3.04221066D-01, AlO
     2  3.27673699D-01, 3.67170655D-01, 4.29260706D-01, 4.98112323D-01, AlO
     3  5.80863558D-01, 8.02081822D-01, 1.10659573D+00, 1.44646594D+00, AlO
     4  1.79208021D+00, 2.28824364D+00, 2.45424300D+00, 2.61027190D+00, AlO
     5  2.75901360D+00, 2.90820341D+00, 3.16208836D+00, 3.48961036D+00, AlO
     6  3.67900216D+00, 4.02378131D+00, 4.32827336D+00, 4.61416346D+00, AlO
     7  4.90411271D+00, 5.43263632D+00, 5.61204793D+00, 5.73420675D+00, AlO
     8       5*0.0D+00/                                                 AlO
      DATA TQ_AlF/                                                      071215
     1 -0.999999993529,-0.857900009158,-0.704500019250,-0.561299999903, AlF
     2 -0.420000008373,-0.326999979556,-0.239200005811, 0.031099978119, AlF
     3  0.191800049230, 0.353200039885, 0.552499957307, 0.755699991316, AlF
     4  1.170000067979, 1.899100085198, 2.171700062029, 2.293999953294, AlF
     5  2.406800108062, 2.636900099285, 2.855499898237, 3.190200038159, AlF
     6  3.392699939082, 3.514800104757, 3.629900119814, 3.724499908038, AlF
     7  3.808399985598, 3.926299927204, 3.970999797938, 3.987999829615, AlF
     8  4.000000000000,      4*0.0D+00/                                 AlF
      DATA  Q_AlF/                                                      071215
     1  1.74483121D-07, 1.44462965D-05, 4.30600918D-04, 4.07289157D-03, AlF
     2  1.98292074D-02, 4.31221190D-02, 7.70170855D-02, 2.42729905D-01, AlF
     3  3.68566934D-01, 5.06309511D-01, 6.86515527D-01, 8.77569697D-01, AlF
     4  1.27944587D+00, 2.00236791D+00, 2.27464103D+00, 2.39798541D+00, AlF
     5  2.51448956D+00, 2.77256395D+00, 3.05938791D+00, 3.58570998D+00, AlF
     6  3.94691234D+00, 4.17771156D+00, 4.40445976D+00, 4.59940711D+00, AlF
     7  4.78220532D+00, 5.06460704D+00, 5.18240702D+00, 5.22902697D+00, AlF
     8  5.26257440D+00,      4*0.0D+00/                                 AlF
      DATA TQ_Al2/                                                      071215
     1 -0.999999993529,-0.955399984673,-0.880300009494,-0.804700011235, Al2
     2 -0.720400013961,-0.616400010099,-0.490900013260,-0.339799959939, Al2
     3 -0.057700009671, 0.234500012800, 0.697600045649, 1.211700020173, Al2
     4  1.535099982089, 1.798399946860, 1.995400004640, 2.163400077354, Al2
     5  2.321399931889, 2.498599983815, 2.662600090014, 2.978700013343, Al2
     6  3.178499869870, 3.365600006593, 3.468199829835, 3.565600068863, Al2
     7  3.700099809986, 3.811899988770, 3.927499894070, 3.971399808732, Al2
     8  4.000000000000,      4*0.0D+00/                                 Al2
      DATA  Q_Al2/                                                      071215
     1  1.79670158D-01, 1.82432412D-01, 1.90707909D-01, 2.05474472D-01, Al2
     2  2.31547126D-01, 2.78447805D-01, 3.53558442D-01, 4.62702879D-01, Al2
     3  6.98391937D-01, 9.66286095D-01, 1.41296800D+00, 1.92116406D+00, Al2
     4  2.24326695D+00, 2.50628662D+00, 2.70590678D+00, 2.88571393D+00, Al2
     5  3.07212395D+00, 3.30844026D+00, 3.55408219D+00, 4.09054974D+00, Al2
     6  4.46501080D+00, 4.84732057D+00, 5.07959112D+00, 5.31788704D+00, Al2
     7  5.66604167D+00, 5.95736309D+00, 6.24763709D+00, 6.35337521D+00, Al2
     8  6.42074270D+00,      4*0.0D+00/                                 Al2
      DATA TQ_SiH/                                                      071215
     1 -0.999999993529,-0.448699991338,-0.045100019165, 0.192000049054, SiH
     2  0.383999979702, 0.482100014920, 0.567000019182, 0.779499985725, SiH
     3  0.943000035377, 1.061799949678, 1.192600041881, 1.435100051986, SiH
     4  1.668100071129, 1.932600053600, 2.254499992232, 2.425900089072, SiH
     5  2.592999917273, 2.782299986205, 3.017100062354, 3.164499940783, SiH
     6  3.302500127242, 3.540299964532, 3.725399928322, 3.795000082881, SiH
     7  3.868899919347, 4.000000000000,      7*0.0D+00/                 SiH
      DATA  Q_SiH/                                                      071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02061504D-01, SiH
     2  6.02259675D-01, 6.03237863D-01, 6.06134107D-01, 6.38598108D-01, SiH
     3  7.05522710D-01, 7.76490435D-01, 8.69636720D-01, 1.06794448D+00, SiH
     4  1.27700589D+00, 1.52641003D+00, 1.83925518D+00, 2.00824875D+00, SiH
     5  2.17436967D+00, 2.36709686D+00, 2.62923668D+00, 2.81812898D+00, SiH
     6  3.01629388D+00, 3.40861989D+00, 3.76407089D+00, 3.91231862D+00, SiH
     7  4.07912479D+00, 4.39169203D+00,      7*0.0D+00/                 SiH
      DATA TQ_SiHm/                                                     071215
     1 -0.999999993529,-0.434700009031,-0.020799970174, 0.206400032349, SiHm
     2  0.393499983945, 0.491300015441, 0.576800054031, 0.778199985295, SiHm
     3  0.929699941660, 1.125500013516, 1.333499922406, 1.587199913864, SiHm
     4  1.882000128973, 2.188600064019, 2.459500026815, 2.679300060708, SiHm
     5  3.012700162001, 3.146099997385, 3.293500049186, 3.517900028258, SiHm
     6  3.774500088687, 3.911299877624, 3.965399897812, 4.000000000000, SiHm
     7       9*0.0D+00/                                                 SiHm
      DATA  Q_SiHm/                                                     071215
     1  4.77121255D-01, 4.77121255D-01, 4.77121255D-01, 4.77123194D-01, SiHm
     2  4.77333669D-01, 4.78352661D-01, 4.81379914D-01, 5.11493239D-01, SiHm
     3  5.70437597D-01, 6.90733002D-01, 8.51450787D-01, 1.07208791D+00, SiHm
     4  1.34666840D+00, 1.64281223D+00, 1.90900757D+00, 2.12726363D+00, SiHm
     5  2.47959421D+00, 2.64041129D+00, 2.83718613D+00, 3.17490421D+00, SiHm
     6  3.60709851D+00, 3.85209840D+00, 3.95110043D+00, 4.01495803D+00, SiHm
     7       9*0.0D+00/                                                 SiHm
      DATA TQ_SiC/                                                      071215
     1 -0.999999993529,-0.845000047468,-0.678799995046,-0.528199957440, SiC
     2 -0.374599998628,-0.279399957553,-0.189999978423,-0.036400005261, SiC
     3  0.083000035411, 0.251500002070, 0.412799944721, 0.789499954197, SiC
     4  1.180000054827, 1.552099946007, 1.915700075706, 2.178800083958, SiC
     5  2.404400107935, 2.583399924025, 2.753299992426, 2.904100100796, SiC
     6  3.065800082247, 3.276300128017, 3.403000013249, 3.535600063748, SiC
     7  3.804100089059, 3.919200080604, 4.000000000000,      6*0.0D+00/ SiC
      DATA  Q_SiC/                                                      071215
     1  7.78151285D-01, 7.78157780D-01, 7.78467483D-01, 7.81761996D-01, SiC
     2  7.98596997D-01, 8.23043459D-01, 8.58405837D-01, 9.44698829D-01, SiC
     3  1.02811457D+00, 1.16125141D+00, 1.29987978D+00, 1.64778264D+00, SiC
     4  2.02604387D+00, 2.39333137D+00, 2.75492156D+00, 3.01736183D+00, SiC
     5  3.24423782D+00, 3.43238950D+00, 3.62856479D+00, 3.82361742D+00, SiC
     6  4.05678668D+00, 4.39811116D+00, 4.62606874D+00, 4.88369498D+00, SiC
     7  5.45037003D+00, 5.70186082D+00, 5.87867313D+00,      6*0.0D+00/ SiC
      DATA TQ_SiN/                                                      071215
     1 -0.999999993529,-0.823899978694,-0.653299972913,-0.555099995812, SiN
     2 -0.461299975749,-0.263799988156,-0.137399989279,-0.006999998873, SiN
     3  0.195800045716, 0.438900046228, 0.741200017466, 1.061199950602, SiN
     4  1.499199992544, 1.933700051478, 2.243200013369, 2.380400125455, SiN
     5  2.505399992678, 2.676200067772, 2.823699941988, 3.196699879008, SiN
     6  3.305800045633, 3.407499895313, 3.493699843039, 3.578799954995, SiN
     7  3.756900150710, 3.831699963757, 3.911199875055, 3.964699917117, SiN
     8  3.986199877136, 4.000000000000,      3*0.0D+00/                 SiN
      DATA  Q_SiN/                                                      071215
     1  3.01029997D-01, 3.01031111D-01, 3.01134345D-01, 3.01733362D-01, SiN
     2  3.04054912D-01, 3.28005971D-01, 3.69237963D-01, 4.36177163D-01, SiN
     3  5.75840844D-01, 7.75643376D-01, 1.04966128D+00, 1.35517538D+00, SiN
     4  1.78479792D+00, 2.21632898D+00, 2.52508777D+00, 2.66258177D+00, SiN
     5  2.78969260D+00, 2.97202518D+00, 3.14488326D+00, 3.67302012D+00, SiN
     6  3.85410868D+00, 4.03650011D+00, 4.20394156D+00, 4.38261344D+00, SiN
     7  4.79975200D+00, 4.99079018D+00, 5.20407391D+00, 5.35422340D+00, SiN
     8  5.41619589D+00, 5.45648830D+00,      3*0.0D+00/                 SiN
      DATA TQ_SiO/                                                      071215
     1 -0.999999993529,-0.820599982167,-0.646799981264,-0.554299995328, SiO
     2 -0.466600018348,-0.277899960069,-0.161500020826,-0.052100013184, SiO
     3  0.092200019938, 0.234500012800, 0.396599985562, 0.568100032837, SiO
     4  0.978699999990, 1.414800078080, 1.865799894008, 2.270699990554, SiO
     5  2.393800102022, 2.511699997685, 2.801999946449, 2.970600015087, SiO
     6  3.134200007045, 3.530600184647, 3.673499971964, 3.739799770088, SiO
     7  3.805500055374, 3.868299933739, 3.924399979667, 3.969999770951, SiO
     8  3.987799834895, 4.000000000000,      3*0.0D+00/                 SiO
      DATA  Q_SiO/                                                      071215
     1  1.15935286D-09, 1.33938779D-06, 1.26393854D-04, 7.43488558D-04, SiO
     2  2.90560674D-03, 2.43602354D-02, 5.94055076D-02, 1.10779562D-01, SiO
     3  2.00921332D-01, 3.06686718D-01, 4.40208997D-01, 5.91427391D-01, SiO
     4  9.76744050D-01, 1.40277234D+00, 1.85006939D+00, 2.25388381D+00, SiO
     5  2.37716547D+00, 2.49655385D+00, 2.81283720D+00, 3.02564589D+00, SiO
     6  3.25831487D+00, 3.91753401D+00, 4.18174485D+00, 4.30934515D+00, SiO
     7  4.44128950D+00, 4.57720906D+00, 4.71343727D+00, 4.84008948D+00, SiO
     8  4.89437679D+00, 4.93333326D+00,      3*0.0D+00/                 SiO
      DATA TQ_SiF/                                                      071215
     1 -0.999999993529,-0.852100035993,-0.692500007476,-0.544199968204, SiF
     2 -0.394999980764,-0.299200023448,-0.209099980800,-0.058500009169, SiF
     3  0.058399930083, 0.222700028433, 0.382799979219, 0.763699990050, SiF
     4  1.159900099393, 1.513299969705, 1.867999888777, 2.147200101857, SiF
     5  2.265699986184, 2.381900121249, 2.597299911552, 2.803699941353, SiF
     6  3.023100056780, 3.274800091954, 3.427699883149, 3.574200068268, SiF
     7  3.689300058339, 3.795400091267, 3.920000101159, 3.968599809561, SiF
     8  3.987399845456, 4.000000000000,      3*0.0D+00/                 SiF
      DATA  Q_SiF/                                                      071215
     1  6.02060068D-01, 6.02069321D-01, 6.02416416D-01, 6.05869535D-01, SiF
     2  6.22428027D-01, 6.47014636D-01, 6.82727141D-01, 7.67242600D-01, SiF
     3  8.48606128D-01, 9.77861592D-01, 1.11494946D+00, 1.46617004D+00, SiF
     4  1.84987126D+00, 2.19864861D+00, 2.55136220D+00, 2.82998960D+00, SiF
     5  2.94888321D+00, 3.06723538D+00, 3.30054781D+00, 3.55675140D+00, SiF
     6  3.87390244D+00, 4.29020620D+00, 4.56564960D+00, 4.84343243D+00, SiF
     7  5.07161154D+00, 5.29215094D+00, 5.57045054D+00, 5.68673426D+00, SiF
     8  5.73305273D+00, 5.76452930D+00,      3*0.0D+00/                 SiF
      DATA TQ_Si2/                                                      071215
     1 -0.999999993529,-0.952599989197,-0.858900004531,-0.769699971523, Si2
     2 -0.677999991947,-0.560799999468,-0.434700009031,-0.312300001113, Si2
     3 -0.130999999799, 0.065199986450, 0.456000041952, 0.857700022592, Si2
     4  1.265499985173, 1.624599991775, 1.895800097534, 2.010000002965, Si2
     5  2.133000102356, 2.246200008430, 2.356599909246, 2.510500000757, Si2
     6  2.665900086395, 2.902100103755, 3.034000106930, 3.164599943217, Si2
     7  3.351399947403, 3.516900052935, 3.641399877969, 3.751700027472, Si2
     8  3.828199969872, 3.896200010746, 3.959600038001, 4.000000000000, Si2
     9       1*0.0D+00/                                                 Si2
      DATA  Q_Si2/                                                      071215
     1  1.77458636D-01, 1.78865469D-01, 1.85174431D-01, 1.98545289D-01, Si2
     2  2.23211945D-01, 2.73146682D-01, 3.47202359D-01, 4.33519931D-01, Si2
     3  5.78235004D-01, 7.49308600D-01, 1.11438016D+00, 1.50554676D+00, Si2
     4  1.90916746D+00, 2.26680374D+00, 2.53758193D+00, 2.65204641D+00, Si2
     5  2.77718827D+00, 2.89774189D+00, 3.02638888D+00, 3.23657948D+00, Si2
     6  3.49314751D+00, 3.94948361D+00, 4.22234176D+00, 4.49719404D+00, Si2
     7  4.89329701D+00, 5.24648401D+00, 5.51563306D+00, 5.76036308D+00, Si2
     8  5.93652984D+00, 6.09925287D+00, 6.25655206D+00, 6.35928555D+00, Si2
     9       1*0.0D+00/                                                 Si2
      DATA TQ_HS/                                                       071215
     1 -0.999999993529,-0.397699973298, 0.043299960753, 0.282899972513, HS
     2  0.476700013520, 0.575900054347, 0.664599955114, 0.874200009715, HS
     3  1.025799983246, 1.222400041473, 1.435600051026, 1.697400052448, HS
     4  1.979400008374, 2.324199930352, 2.536599967754, 2.724500016942, HS
     5  2.903900101092, 3.102799982484, 3.237899828470, 3.368299940722, HS
     6  3.560599952013, 3.668400026569, 3.759300207590, 3.829999926729, HS
     7  3.893399944841, 3.958600016166, 4.000000000000,      6*0.0D+00/ HS
      DATA  Q_HS/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02061087D-01, HS
     2  6.02228460D-01, 6.03107361D-01, 6.05953101D-01, 6.36819263D-01, HS
     3  6.96133372D-01, 8.17316564D-01, 9.82641975D-01, 1.21128822D+00, HS
     4  1.47454069D+00, 1.80848807D+00, 2.01777092D+00, 2.20470120D+00, HS
     5  2.38769014D+00, 2.60748329D+00, 2.77541049D+00, 2.95592110D+00, HS
     6  3.25810072D+00, 3.44763590D+00, 3.62069846D+00, 3.76519663D+00, HS
     7  3.90275390D+00, 4.05155063D+00, 4.14906525D+00,      6*0.0D+00/ HS
      DATA TQ_HSm/                                                      071215
     1 -0.999999993529,-0.397499973851, 0.043499960572, 0.282999972523, HSm
     2  0.476700013520, 0.575900054347, 0.664599955114, 0.874200009715, HSm
     3  1.025599982946, 1.221000046243, 1.434300053521, 1.695700056854, HSm
     4  1.987200008977, 2.321399931889, 2.538699961332, 2.730200002646, HSm
     5  2.909700092510, 3.107800099601, 3.238099823338, 3.365999996834, HSm
     6  3.557600000043, 3.754000081981, 3.821200137651, 3.890099867167, HSm
     7  3.957299987780, 4.000000000000,      7*0.0D+00/                 HSm
      DATA  Q_HSm/                                                      071215
     1  0.00000000D+00,-1.34899922D-36, 3.93709229D-11, 1.09927521D-06, HSm
     2  1.68468655D-04, 1.04736957D-03, 3.89310947D-03, 3.47592719D-02, HSm
     3  9.39728425D-02, 2.14264126D-01, 3.79504105D-01, 6.07683588D-01, HSm
     4  8.79920358D-01, 1.20368334D+00, 1.41778972D+00, 1.60835517D+00, HSm
     5  1.79171067D+00, 2.01132623D+00, 2.17361291D+00, 2.35051027D+00, HSm
     6  2.65098201D+00, 3.00797567D+00, 3.14403901D+00, 3.29177936D+00, HSm
     7  3.44309215D+00, 3.54193175D+00,      7*0.0D+00/                 HSm
      DATA TQ_CS/                                                       071215
     1 -0.999999993529,-0.594800033588,-0.421000009558,-0.233300007614, CS
     2 -0.132499997333,-0.036600006303, 0.099100010276, 0.219100030570, CS
     3  0.387499981112, 0.541599954669, 0.884699991415, 1.240499991575, CS
     4  1.705300047772, 2.046399964714, 2.301699959574, 2.421000092272, CS
     5  2.534799973259, 2.813499929213, 2.985600011396, 3.152800018525, CS
     6  3.530000199155, 3.671400024324, 3.737499831557, 3.805100064998, CS
     7  3.922600029368, 3.969499784740, 3.987599840176, 4.000000000000, CS
     8       5*0.0D+00/                                                 CS
      DATA  Q_CS/                                                       071215
     1  8.01415913D-11, 1.25352531D-04, 2.64008050D-03, 2.27204071D-02, CS
     2  5.07156588D-02, 9.15700359D-02, 1.70262726D-01, 2.54662759D-01, CS
     3  3.88264184D-01, 5.20852712D-01, 8.36763216D-01, 1.18006872D+00, CS
     4  1.63846469D+00, 1.97781591D+00, 2.23259816D+00, 2.35215192D+00, CS
     5  2.46758164D+00, 2.77151343D+00, 2.98858962D+00, 3.22685419D+00, CS
     6  3.85352237D+00, 4.11486738D+00, 4.24308129D+00, 4.38147542D+00, CS
     7  4.65871624D+00, 4.79162984D+00, 4.84738875D+00, 4.88712667D+00, CS
     8       5*0.0D+00/                                                 CS
      DATA TQ_NS/                                                       071215
     1 -0.999999993529,-0.616500009978,-0.448199991629,-0.266899980740, NS
     2 -0.179799988068,-0.092799983671, 0.103800013296, 0.298599964622, NS
     3  0.545399959636, 0.801099948606, 1.087099927455, 1.404000087109, NS
     4  1.619599890315, 1.808299942816, 1.949700062618, 2.089199924524, NS
     5  2.311599940309, 2.501299984336, 2.816199934416, 2.955600044242, NS
     6  3.106100059781, 3.478799997437, 3.611100100634, 3.739999764743, NS
     7  3.826100020206, 3.906199941489, 3.962799969516, 3.985599892977, NS
     8  4.000000000000,      4*0.0D+00/                                 NS
      DATA  Q_NS/                                                       071215
     1  3.01029996D-01, 3.01172594D-01, 3.03692557D-01, 3.22581368D-01, NS
     2  3.45142176D-01, 3.79082782D-01, 4.93012552D-01, 6.39188168D-01, NS
     3  8.49961936D-01, 1.08506602D+00, 1.35878612D+00, 1.66890741D+00, NS
     4  1.88226086D+00, 2.07245906D+00, 2.22204883D+00, 2.38026463D+00, NS
     5  2.65395209D+00, 2.89728065D+00, 3.31506618D+00, 3.51514373D+00, NS
     6  3.74664636D+00, 4.38782094D+00, 4.63550855D+00, 4.88634033D+00, NS
     7  5.06052691D+00, 5.23031067D+00, 5.35741257D+00, 5.41088207D+00, NS
     8  5.44544822D+00,      4*0.0D+00/                                 NS
      DATA TQ_SO/                                                       071215
     1 -0.999999993529,-0.815399998091,-0.623800006435,-0.467700027189, SO
     2 -0.298400022363,-0.211699986048,-0.126299998750, 0.105400015099, SO
     3  0.237300007278, 0.377599980201, 0.665799956527, 0.951100029384, SO
     4  1.333199923478, 1.733600021145, 2.090599923298, 2.362599901313, SO
     5  2.531199984269, 2.669100082886, 2.799499953725, 2.934600061650, SO
     6  3.126300028258, 3.264300115565, 3.512000173853, 3.738399807504, SO
     7  3.821100140048, 3.906799926131, 3.962799969516, 3.985699890337, SO
     8  4.000000000000,      4*0.0D+00/                                 SO
      DATA  Q_SO/                                                       071215
     1  4.77121256D-01, 4.77123031D-01, 4.77340894D-01, 4.80137842D-01, SO
     2  4.98062011D-01, 5.20102980D-01, 5.52859671D-01, 6.89933685D-01, SO
     3  7.88998790D-01, 9.04347382D-01, 1.16137896D+00, 1.43095125D+00, SO
     4  1.80321983D+00, 2.19947231D+00, 2.55501565D+00, 2.82695242D+00, SO
     5  2.99868226D+00, 3.14641649D+00, 3.29756475D+00, 3.46995391D+00, SO
     6  3.74525534D+00, 3.96606749D+00, 4.41402843D+00, 4.88529406D+00, SO
     7  5.07083587D+00, 5.27125272D+00, 5.40819456D+00, 5.46595567D+00, SO
     8  5.50262626D+00,      4*0.0D+00/                                 SO
      DATA TQ_MgS/                                                      071215
     1 -0.999999993529,-0.809500015589,-0.662299962062,-0.509099992450, MgS
     2 -0.297900021685,-0.013199989312, 0.328700050178, 0.597200025992, MgS
     3  0.859300022279, 1.382500097161, 1.944300049651, 2.160700085771, MgS
     4  2.388800101905, 2.567799947764, 2.750599988136, 3.123099948492, MgS
     5  3.323800063515, 3.511500186192, 3.793100043046, 3.913899944428, MgS
     6  3.966599864718, 4.000000000000,     11*0.0D+00/                 MgS
      DATA  Q_MgS/                                                      071215
     1  5.98134506D-04, 9.07415761D-03, 3.65871871D-02, 9.81887039D-02, MgS
     2  2.31537358D-01, 4.60348450D-01, 7.70319579D-01, 1.02667183D+00, MgS
     3  1.28237401D+00, 1.80021522D+00, 2.36057384D+00, 2.57919229D+00, MgS
     4  2.82564576D+00, 3.04551507D+00, 3.30150673D+00, 3.91398438D+00, MgS
     5  4.28360091D+00, 4.65025167D+00, 5.25209900D+00, 5.54002256D+00, MgS
     6  5.67140478D+00, 5.75604548D+00,     11*0.0D+00/                 MgS
      DATA TQ_AlS/                                                      071215
     1 -0.999999993529,-0.854000027202,-0.649199978887,-0.495700007805, AlS
     2 -0.285099984301, 0.003300003650, 0.355400036348, 0.630700002376, AlS
     3  0.900499986763, 1.430000061774, 1.995300004741, 2.214500034622, AlS
     4  2.426500088680, 2.635400100178, 2.849399897827, 3.141799893763, AlS
     5  3.317600103583, 3.492399812699, 3.641099870346, 3.803100113119, AlS
     6  3.919300083174, 3.968699806803, 4.000000000000,     10*0.0D+00/ AlS
      DATA  Q_AlS/                                                      071215
     1  3.01454659D-01, 3.05214286D-01, 3.35969691D-01, 3.96389544D-01, AlS
     2  5.28118497D-01, 7.59265575D-01, 1.07863424D+00, 1.34178680D+00, AlS
     3  1.60526860D+00, 2.12964137D+00, 2.69355265D+00, 2.91463076D+00, AlS
     4  3.14123776D+00, 3.39528116D+00, 3.69825407D+00, 4.17760919D+00, AlS
     5  4.49492912D+00, 4.82835625D+00, 5.12676714D+00, 5.47474098D+00, AlS
     6  5.74624415D+00, 5.86794302D+00, 5.94686136D+00,     10*0.0D+00/ AlS
      DATA TQ_SiS/                                                      071215
     1 -0.999999993529,-0.823499979115,-0.722800012565,-0.624300006529, SiS
     2 -0.539799950497,-0.463899996646,-0.264999985285,-0.117899999208, SiS
     3  0.030599978927, 0.212700027992, 0.399099986866, 0.678099938227, SiS
     4  0.951900031727, 1.501299989459, 2.035999971282, 2.242200015015, SiS
     5  2.464600026679, 2.687200048640, 2.930900071734, 3.225399914778, SiS
     6  3.478699995020, 3.589700152706, 3.692999996240, 3.774100080606, SiS
     7  3.851299921779, 3.942499987731, 3.977599976051, 3.989799782095, SiS
     8  4.000000000000,      4*0.0D+00/                                 SiS
      DATA  Q_SiS/                                                      071215
     1  2.14205530D-04, 3.91454561D-03, 1.28772929D-02, 3.20719906D-02, SiS
     2  5.95748044D-02, 9.34759073D-02, 2.16035157D-01, 3.27577694D-01, SiS
     3  4.51135985D-01, 6.12583438D-01, 7.85281240D-01, 1.05224782D+00, SiS
     4  1.31983481D+00, 1.86419736D+00, 2.39761336D+00, 2.60460721D+00, SiS
     5  2.83728407D+00, 3.10015286D+00, 3.44003486D+00, 3.92016920D+00, SiS
     6  4.37905751D+00, 4.59051511D+00, 4.79269326D+00, 4.95619678D+00, SiS
     7  5.11912288D+00, 5.33102935D+00, 5.42182905D+00, 5.45491766D+00, SiS
     8  5.48323443D+00,      4*0.0D+00/                                 SiS
      DATA TQ_S2/                                                       071215
     1 -0.999999993529,-0.826699975747,-0.726900010180,-0.630400006814, S2
     2 -0.545799974970,-0.470900044023,-0.263099989830,-0.120999994893, S2
     3  0.026599988157, 0.385599980347, 0.657699951815, 0.921899937301, S2
     4  1.447700043680, 2.050699960392, 2.264999985434, 2.511399998453, S2
     5  2.625399997955, 2.736800011348, 3.064100038282, 3.287500034642, S2
     6  3.476199934585, 3.571500134753, 3.668200022487, 3.753900079611, S2
     7  3.821500130460, 3.929699833323, 3.972399835719, 3.988399819055, S2
     8  4.000000000000,      4*0.0D+00/                                 S2
      DATA  Q_S2/                                                       071215
     1  1.76361725D-01, 1.80474397D-01, 1.90034476D-01, 2.09675050D-01, S2
     2  2.37948954D-01, 2.71981473D-01, 4.01915722D-01, 5.10635590D-01, S2
     3  6.33935507D-01, 9.59759175D-01, 1.21999166D+00, 1.47801400D+00, S2
     4  1.99867094D+00, 2.60018154D+00, 2.81589457D+00, 3.07912760D+00, S2
     5  3.21416842D+00, 3.35744710D+00, 3.84701436D+00, 4.23844127D+00, S2
     6  4.60712882D+00, 4.80607484D+00, 5.01657217D+00, 5.21232351D+00, S2
     7  5.37577915D+00, 5.66241901D+00, 5.78589818D+00, 5.83371531D+00, S2
     8  5.86889648D+00,      4*0.0D+00/                                 S2
      DATA TQ_HCl/                                                      071215
     1 -0.999999993529,-0.375899994328, 0.080900039129, 0.324200046812, HCl
     2  0.519999972683, 0.619800019211, 0.709400048872, 0.923699938307, HCl
     3  1.067899940285, 1.236899997954, 1.452800039027, 1.715700044973, HCl
     4  2.005000001483, 2.338699922576, 2.574199939989, 2.776599985625, HCl
     5  2.952700046430, 3.129200100545, 3.255100113486, 3.383400029028, HCl
     6  3.586800084762, 3.770200001819, 3.905099969646, 3.962999964000, HCl
     7  3.985599892977, 4.000000000000,      7*0.0D+00/                 HCl
      DATA  Q_HCl/                                                      071215
     1  0.00000000D+00,-9.39717639D-38, 1.99035065D-11, 8.53344697D-07, HCl
     2  1.49698465D-04, 9.63023997D-04, 3.68065922D-03, 3.48994290D-02, HCl
     3  9.06397384D-02, 1.91254224D-01, 3.54426745D-01, 5.81437356D-01, HCl
     4  8.50206664D-01, 1.17263900D+00, 1.40442517D+00, 1.60575283D+00, HCl
     5  1.78550706D+00, 1.97930300D+00, 2.13259826D+00, 2.30549897D+00, HCl
     6  2.61726366D+00, 2.93949400D+00, 3.20513208D+00, 3.32772030D+00, HCl
     7  3.37695322D+00, 3.40870792D+00,      7*0.0D+00/                 HCl
      DATA TQ_LiCl/                                                     071215
     1 -0.999999993529,-0.815999996096,-0.629400007485,-0.478000030989, LiCl
     2 -0.309200004140,-0.233200007644,-0.158200015132, 0.123900002012, LiCl
     3  0.279299972120, 0.436400020337, 0.808799949302, 1.182600052045, LiCl
     4  1.861999903043, 2.104099903263, 2.319299933295, 2.539299959497, LiCl
     5  2.744900001658, 3.086100063422, 3.313200009900, 3.477799973263, LiCl
     6  3.679099832337, 3.795000082881, 3.907399910774, 3.964099933664, LiCl
     7  4.000000000000,      8*0.0D+00/                                 LiCl
      DATA  Q_LiCl/                                                     071215
     1  2.16536232D-09, 2.33132319D-06, 2.37058651D-04, 2.98037541D-03, LiCl
     2  2.06774719D-02, 3.92527497D-02, 6.60466290D-02, 2.32844563D-01, LiCl
     3  3.53164993D-01, 4.86059227D-01, 8.27021594D-01, 1.18764030D+00, LiCl
     4  1.85966287D+00, 2.10146401D+00, 2.32179221D+00, 2.56929294D+00, LiCl
     5  2.83813132D+00, 3.37383665D+00, 3.78337963D+00, 4.10383780D+00, LiCl
     6  4.52933267D+00, 4.79571343D+00, 5.06447376D+00, 5.20001497D+00, LiCl
     7  5.28461641D+00,      8*0.0D+00/                                 LiCl
      DATA TQ_NaCl/                                                     071215
     1 -0.999999993529,-0.954899985481,-0.877600004270,-0.799600007213, NaCl
     2 -0.711500047454,-0.603000014229,-0.468400032815,-0.339299960560, NaCl
     3 -0.062100004665, 0.243000002369, 0.721200033916, 1.250400019543, NaCl
     4  1.546299958620, 1.804599942427, 2.016199995143, 2.155900089353, NaCl
     5  2.299299962063, 2.464300026649, 2.617799880793, 3.024800089739, NaCl
     6  3.228399990697, 3.404599971317, 3.519699983839, 3.644999969448, NaCl
     7  3.749999987182, 3.885499979279, 3.953599906988, 3.982399977458, NaCl
     8  4.000000000000,      4*0.0D+00/                                 NaCl
      DATA  Q_NaCl/                                                     071215
     1  2.50427115D-03, 4.62829204D-03, 1.15087415D-02, 2.46236786D-02, NaCl
     2  4.94195700D-02, 9.60404361D-02, 1.75392489D-01, 2.67381833D-01, NaCl
     3  4.95850195D-01, 7.74074627D-01, 1.23488361D+00, 1.75806170D+00, NaCl
     4  2.05277407D+00, 2.31073569D+00, 2.52500945D+00, 2.67353314D+00, NaCl
     5  2.83875314D+00, 3.05124299D+00, 3.27264555D+00, 3.95687964D+00, NaCl
     6  4.34105239D+00, 4.69454401D+00, 4.93893384D+00, 5.22013986D+00, NaCl
     7  5.46687083D+00, 5.79042959D+00, 5.94957798D+00, 6.01528159D+00, NaCl
     8  6.05484882D+00,      4*0.0D+00/                                 NaCl
      DATA TQ_AlCl/                                                     071215
     1 -0.999999993529,-0.955499984511,-0.869199985077,-0.778899969247, AlCl
     2 -0.684700005875,-0.543499965244,-0.435500007127,-0.312900000786, AlCl
     3 -0.171000002636,-0.020499970118, 0.305799936659, 0.569500050216, AlCl
     4  0.828500045474, 1.335599914899, 1.913900075723, 2.131700097440, AlCl
     5  2.378300086008, 2.554299949098, 2.721000026207, 3.001700037935, AlCl
     6  3.216299868471, 3.418500053455, 3.601799941653, 3.716299891262, AlCl
     7  3.843900045364, 3.936599967473, 3.975299913981, 4.000000000000, AlCl
     8       5*0.0D+00/                                                 AlCl
      DATA  Q_AlCl/                                                     071215
     1  1.19123479D-03, 2.35137571D-03, 7.29605511D-03, 1.90276352D-02, AlCl
     2  4.21168326D-02, 1.01519442D-01, 1.64958960D-01, 2.50517093D-01, AlCl
     3  3.61768457D-01, 4.89488683D-01, 7.87179653D-01, 1.03937933D+00, AlCl
     4  1.29222696D+00, 1.79417533D+00, 2.37100885D+00, 2.59135219D+00, AlCl
     5  2.86074751D+00, 3.08134538D+00, 3.31763328D+00, 3.77119541D+00, AlCl
     6  4.15567817D+00, 4.54158469D+00, 4.91132815D+00, 5.15661093D+00, AlCl
     7  5.45451355D+00, 5.69555004D+00, 5.80314895D+00, 5.87385640D+00, AlCl
     8       5*0.0D+00/                                                 AlCl
      DATA TQ_CaH/                                                      071215
     1 -0.999999993529,-0.512499987824,-0.155700003171, 0.024599993016, CaH
     2  0.187500021030, 0.317000005796, 0.456000041952, 0.585800036448, CaH
     3  0.932099962302, 1.093799922793, 1.264799985684, 1.478000018068, CaH
     4  1.706000048044, 1.973700029746, 2.245700009253, 2.389600099662, CaH
     5  2.523199980116, 2.715700040553, 2.866199895123, 3.084100024579, CaH
     6  3.277800164079, 3.387799931390, 3.489399772872, 3.610700109669, CaH
     7  3.699599818482, 3.777400147272, 3.856200040645, 3.944899936738, CaH
     8  4.000000000000,      4*0.0D+00/                                 CaH
      DATA  Q_CaH/                                                      071215
     1  3.01029996D-01, 3.01029996D-01, 3.01030031D-01, 3.01043253D-01, CaH
     2  3.01512757D-01, 3.04717660D-01, 3.19105390D-01, 3.53337475D-01, CaH
     3  5.55077327D-01, 6.83055894D-01, 8.30237768D-01, 1.02457287D+00, CaH
     4  1.24059928D+00, 1.50054454D+00, 1.76870342D+00, 1.91188452D+00, CaH
     5  2.04676238D+00, 2.25140468D+00, 2.42873951D+00, 2.72501888D+00, CaH
     6  3.03146990D+00, 3.22378877D+00, 3.41490576D+00, 3.66691881D+00, CaH
     7  3.87529091D+00, 4.07752631D+00, 4.30038336D+00, 4.56713980D+00, CaH
     8  4.73726728D+00,      4*0.0D+00/                                 CaH
      DATA TQ_CaF/                                                      071215
     1 -0.999999993529,-0.903199970838,-0.808700014864,-0.592400043245, CaF
     2 -0.512299988111,-0.434100010459,-0.251999994973,-0.111400014815, CaF
     3  0.035099971656, 0.207800030232, 0.399299986971, 0.654899954392, CaF
     4  0.915899947482, 1.436400049491, 1.980300006243, 2.194800050834, CaF
     5  2.412800103945, 2.595099914479, 2.772399981901, 3.224999904655, CaF
     6  3.422200032989, 3.510800203466, 3.595900007258, 3.833199996429, CaF
     7  3.932599881150, 3.973599868103, 4.000000000000,      6*0.0D+00/ CaF
      DATA  Q_CaF/                                                      071215
     1  3.01109593D-01, 3.01582840D-01, 3.03547765D-01, 3.29377316D-01, CaF
     2  3.53380889D-01, 3.86213440D-01, 4.93547260D-01, 5.96828275D-01, CaF
     3  7.16133904D-01, 8.66954727D-01, 1.04265325D+00, 1.28566654D+00, CaF
     4  1.53961233D+00, 2.05418985D+00, 2.59649086D+00, 2.81299165D+00, CaF
     5  3.04706177D+00, 3.26851669D+00, 3.51402190D+00, 4.26495969D+00, CaF
     6  4.63546593D+00, 4.80969873D+00, 4.98301026D+00, 5.52540725D+00, CaF
     7  5.79423027D+00, 5.91292244D+00, 5.99137474D+00,      6*0.0D+00/ CaF
      DATA TQ_CaCl/                                                     071215
     1 -0.999999993529,-0.984000009244,-0.957999980472,-0.891499969784, CaCl
     2 -0.818299988450,-0.741500023085,-0.531099952839,-0.338099962051, CaCl
     3 -0.145499979997, 0.055799937700, 0.390299982276, 0.753799989938, CaCl
     4  1.164500085269, 1.517099967847, 1.805399942511, 2.010700002082, CaCl
     5  2.182400081066, 2.353099912577, 2.498599983815, 2.643400095102, CaCl
     6  2.892900109572, 3.159699839019, 3.402800018491, 3.526000110067, CaCl
     7  3.635399971030, 3.746799916002, 3.838700116225, 3.936299960999, CaCl
     8  3.975299913981, 4.000000000000,      3*0.0D+00/                 CaCl
      DATA  Q_CaCl/                                                     071215
     1  3.17225685D-01, 3.19939575D-01, 3.25046695D-01, 3.42413338D-01, CaCl
     2  3.69245119D-01, 4.05946593D-01, 5.41387435D-01, 6.93704577D-01, CaCl
     3  8.60864587D-01, 1.04546441D+00, 1.36488232D+00, 1.72105979D+00, CaCl
     4  2.12836928D+00, 2.47981561D+00, 2.76784689D+00, 2.97547567D+00, CaCl
     5  3.15845391D+00, 3.35962356D+00, 3.55189708D+00, 3.76338643D+00, CaCl
     6  4.16989951D+00, 4.65052991D+00, 5.11944652D+00, 5.36833393D+00, CaCl
     7  5.60049286D+00, 5.85812407D+00, 6.09625539D+00, 6.37881618D+00, CaCl
     8  6.49929916D+00, 6.57733019D+00,      3*0.0D+00/                 CaCl
      DATA TQ_ScO/                                                      071215
     1 -0.999999993529,-0.864699992103,-0.717600023585,-0.593900037209, ScO
     2 -0.392199988506,-0.271399970975,-0.137599988950, 0.070400042724, ScO
     3  0.322000045166, 0.662599952759, 1.025299982495, 1.456900031734, ScO
     4  1.927500062877, 2.237900017994, 2.371999935668, 2.503099987998, ScO
     5  2.705400051300, 2.917000080585, 3.234699910583, 3.359500131828, ScO
     6  3.498599957399, 3.616899969630, 3.792600032564, 3.915799993246, ScO
     7  3.967299845413, 3.986899858656, 4.000000000000,      6*0.0D+00/ ScO
      DATA  Q_ScO/                                                      071215
     1  3.01030520D-01, 3.01056986D-01, 3.01627859D-01, 3.05031797D-01, ScO
     2  3.34203517D-01, 3.77626739D-01, 4.49931611D-01, 5.97298949D-01, ScO
     3  8.07399170D-01, 1.11987256D+00, 1.46935032D+00, 1.89462045D+00, ScO
     4  2.36283033D+00, 2.67285901D+00, 2.80799548D+00, 2.94367054D+00, ScO
     5  3.17034766D+00, 3.44419663D+00, 3.93406625D+00, 4.14877722D+00, ScO
     6  4.40132745D+00, 4.63083597D+00, 5.02020859D+00, 5.34623469D+00, ScO
     7  5.49624505D+00, 5.55520195D+00, 5.59512884D+00,      6*0.0D+00/ ScO
      DATA TQ_TiO/                                                      071215
     1 -0.999999993529,-0.747100010471,-0.574600023221,-0.400899972496, TiO
     2 -0.282599970695,-0.159800022788,-0.002699999565, 0.189600046046, TiO
     3  0.455000042295, 0.737200025718, 0.997200000451, 1.240699992190, TiO
     4  1.436900048532, 1.684700076002, 1.808999942890, 1.924900067310, TiO
     5  2.217400021492, 2.359699906296, 2.499999981691, 2.691900044784, TiO
     6  2.910800090754, 3.133300031354, 3.358900118167, 3.515900077612, TiO
     7  3.754000081981, 3.898800071943, 3.960700027431, 3.984899911457, TiO
     8  4.000000000000,      4*0.0D+00/                                 TiO
      DATA  Q_TiO/                                                      071215
     1  3.01030273D-01, 3.01274188D-01, 3.05080909D-01, 3.27483610D-01, TiO
     2  3.64940555D-01, 4.25655360D-01, 5.28510302D-01, 6.78420299D-01, TiO
     3  9.09962868D-01, 1.17325017D+00, 1.42401228D+00, 1.66284843D+00, TiO
     4  1.85931090D+00, 2.12782453D+00, 2.27879673D+00, 2.43053847D+00, TiO
     5  2.84076121D+00, 3.04041990D+00, 3.23317233D+00, 3.49713296D+00, TiO
     6  3.81649449D+00, 4.17637054D+00, 4.58483751D+00, 4.89707857D+00, TiO
     7  5.43634501D+00, 5.82338312D+00, 6.00381009D+00, 6.07656843D+00, TiO
     8  6.12255387D+00,      4*0.0D+00/                                 TiO
      DATA TQ_TiS/                                                      071215
     1 -0.999999993529,-0.962699988814,-0.901499970617,-0.753800001700, TiS
     2 -0.626000006848,-0.459699966058,-0.235500006941, 0.029599980869, TiS
     3  0.290699972528, 0.569000044010, 0.874600009447, 1.152200097713, TiS
     4  1.324999929822, 1.581699916054, 1.780199972154, 2.114299882624, TiS
     5  2.258499983395, 2.418100095833, 2.611399900040, 2.815299932682, TiS
     6  3.128300078111, 3.362000094420, 3.548800167838, 3.753000058282, TiS
     7  3.828499962682, 3.895900003684, 3.959600038001, 4.000000000000, TiS
     8       5*0.0D+00/                                                 TiS
      DATA  Q_TiS/                                                      071215
     1  3.04980694D-01, 3.07377080D-01, 3.13705640D-01, 3.47284079D-01, TiS
     2  4.02224084D-01, 5.04749031D-01, 6.77104022D-01, 9.08244648D-01, TiS
     3  1.15137493D+00, 1.41943747D+00, 1.71930442D+00, 1.99431335D+00, TiS
     4  2.16705258D+00, 2.43665668D+00, 2.67251655D+00, 3.12995983D+00, TiS
     5  3.33813508D+00, 3.57275630D+00, 3.86905320D+00, 4.20282216D+00, TiS
     6  4.75828880D+00, 5.20124757D+00, 5.57390640D+00, 6.01744349D+00, TiS
     7  6.19668319D+00, 6.36512198D+00, 6.53145274D+00, 6.64006381D+00, TiS
     8       5*0.0D+00/                                                 TiS
      DATA TQ_VO/                                                       071215
     1 -0.999999993529,-0.859600001293,-0.707600038455,-0.564800002948, VO
     2 -0.423600012637,-0.334999965903,-0.245999998869,-0.103600014237, VO
     3  0.049099955507, 0.220200030748, 0.396899985719, 0.595800025708, VO
     4  0.799299948788, 1.216100035800, 1.991100008977, 2.255299990465, VO
     5  2.374299990554, 2.490799995645, 2.682300055768, 2.880300105369, VO
     6  3.047000120870, 3.233999928545, 3.367799952921, 3.506400139279, VO
     7  3.655499963731, 3.789299993576, 3.917700042064, 3.967899828866, VO
     8  4.000000000000,      4*0.0D+00/                                 VO
      DATA  Q_VO/                                                       071215
     1  6.02060184D-01, 6.02074849D-01, 6.02487979D-01, 6.06093770D-01, VO
     2  6.21735328D-01, 6.43535218D-01, 6.77266632D-01, 7.54390207D-01, VO
     3  8.60229129D-01, 9.96485993D-01, 1.14955244D+00, 1.33136365D+00, VO
     4  1.52394150D+00, 1.92953879D+00, 2.69886838D+00, 2.96278406D+00, VO
     5  3.08255047D+00, 3.20221370D+00, 3.41209074D+00, 3.65910644D+00, VO
     6  3.89613393D+00, 4.19219850D+00, 4.42126394D+00, 4.67387073D+00, VO
     7  4.97070656D+00, 5.27182752D+00, 5.59908393D+00, 5.73683051D+00, VO
     8  5.82740434D+00,      4*0.0D+00/                                 VO
      DATA TQ_CrH/                                                      071215
     1 -0.999999993529,-0.455599976423,-0.057000010110, 0.149199973352, CrH
     2  0.329500050777, 0.469900012068, 0.626000010311, 0.763799989932, CrH
     3  1.041899954370, 1.174700061798, 1.320199924932, 1.565999945509, CrH
     4  1.794499958848, 2.079999933928, 2.361399903481, 2.499299982753, CrH
     5  2.628800075890, 2.961200037838, 3.203299872162, 3.363900048067, CrH
     6  3.468299827328, 3.597999952987, 3.697199883122, 3.804000091465, CrH
     7  3.912599911026, 3.965699889538, 4.000000000000,      6*0.0D+00/ CrH
      DATA  Q_CrH/                                                      071215
     1  7.78151250D-01, 7.78151250D-01, 7.78151253D-01, 7.78156044D-01, CrH
     2  7.78487647D-01, 7.81435724D-01, 7.97761458D-01, 8.36643433D-01, CrH
     3  9.94560101D-01, 1.09477346D+00, 1.21494844D+00, 1.43375621D+00, CrH
     4  1.64801058D+00, 1.92388642D+00, 2.20084288D+00, 2.33800905D+00, CrH
     5  2.46882124D+00, 2.83885668D+00, 3.17117774D+00, 3.42900906D+00, CrH
     6  3.61611897D+00, 3.88231251D+00, 4.12429163D+00, 4.42814797D+00, CrH
     7  4.77026411D+00, 4.94189616D+00, 5.05223749D+00,      6*0.0D+00/ CrH
      DATA TQ_CrO/                                                      071215
     1 -0.999999993529,-0.867999986950,-0.722400012798,-0.577400032807, CrO
     2 -0.435000008317,-0.338499961554,-0.248599995967,-0.094599990759, CrO
     3  0.026999987186, 0.198000043784, 0.361000024744, 0.556199944964, CrO
     4  0.759799994289, 1.172100065217, 1.533299984153, 1.942200044609, CrO
     5  2.207000050713, 2.324699930077, 2.434300077901, 2.662000090671, CrO
     6  2.877100043654, 3.380800086723, 3.484399907751, 3.588600126934, CrO
     7  3.766500077019, 3.833099994251, 3.899300083712, 3.961799997094, CrO
     8  4.000000000000,      4*0.0D+00/                                 CrO
      DATA  Q_CrO/                                                      071215
     1  1.00000034D+00, 1.00001797D+00, 1.00043432D+00, 1.00419164D+00, CrO
     2  1.02043482D+00, 1.04529017D+00, 1.08103915D+00, 1.16789264D+00, CrO
     3  1.25313844D+00, 1.38872570D+00, 1.52922459D+00, 1.70675083D+00, CrO
     4  1.89878639D+00, 2.29930282D+00, 2.65631817D+00, 3.06325362D+00, CrO
     5  3.32764208D+00, 3.44607647D+00, 3.55857836D+00, 3.81043444D+00, CrO
     6  4.08723915D+00, 4.90107455D+00, 5.09233888D+00, 5.29270782D+00, CrO
     7  5.66135748D+00, 5.81146307D+00, 5.96784937D+00, 6.12098023D+00, CrO
     8  6.21634269D+00,      4*0.0D+00/                                 CrO
      DATA TQ_FeH/                                                      071215
     1 -0.999999993529,-0.451399987040,-0.049800014691, 0.163599967317, FeH
     2  0.349000043697, 0.491000015792, 0.656099953288, 0.782699976883, FeH
     3  1.056299953280, 1.195500039378, 1.363099883249, 1.491199996295, FeH
     4  1.620699905550, 1.873199963187, 1.990700009381, 2.117799881795, FeH
     5  2.299299962063, 2.492499993066, 2.646800092766, 2.859899898641, FeH
     6  3.027200136270, 3.196199891251, 3.396000008152, 3.503500071671, FeH
     7  3.604700006665, 3.753200063022, 3.878900094298, 3.953599906988, FeH
     8  3.982399977458, 4.000000000000,      3*0.0D+00/                 FeH
      DATA  Q_FeH/                                                      071215
     1  3.01029996D-01, 3.01029996D-01, 3.01029997D-01, 3.01033582D-01, FeH
     2  3.01336590D-01, 3.04170685D-01, 3.21721618D-01, 3.57441174D-01, FeH
     3  5.10395873D-01, 6.14816325D-01, 7.53695444D-01, 8.66500007D-01, FeH
     4  9.85095998D-01, 1.23353132D+00, 1.36204999D+00, 1.51362175D+00, FeH
     5  1.75344115D+00, 2.03700842D+00, 2.28490539D+00, 2.66677765D+00, FeH
     6  3.00844316D+00, 3.39612270D+00, 3.90056344D+00, 4.18562787D+00, FeH
     7  4.45999950D+00, 4.87124721D+00, 5.22290637D+00, 5.42957366D+00, FeH
     8  5.50813201D+00, 5.55573511D+00,      3*0.0D+00/                 FeH
      DATA TQ_FeO/                                                      071215
     1 -0.999999993529,-0.873599993511,-0.733400014527,-0.590100052500, FeO
     2 -0.450299989820,-0.358500025705,-0.266299982175,-0.132299997662, FeO
     3  0.019300004673, 0.187700023412, 0.362600018013, 0.563499975735, FeO
     4  0.769699982941, 1.194000040673, 1.967700046342, 2.240800017319, FeO
     5  2.363099900409, 2.482200006244, 2.679000061392, 2.881300105916, FeO
     6  3.138799882798, 3.324900037125, 3.443799950951, 3.567100103917, FeO
     7  3.765900090603, 3.906099944049, 3.963499950211, 4.000000000000, FeO
     8       5*0.0D+00/                                                 FeO
      DATA  Q_FeO/                                                      071215
     1  1.00000055D+00, 1.00002243D+00, 1.00046241D+00, 1.00429174D+00, FeO
     2  1.02029450D+00, 1.04353191D+00, 1.07960964D+00, 1.15315189D+00, FeO
     3  1.25832716D+00, 1.39234566D+00, 1.54366558D+00, 1.72716588D+00, FeO
     4  1.92226897D+00, 2.33524587D+00, 3.10328428D+00, 3.37599921D+00, FeO
     5  3.49905897D+00, 3.62149223D+00, 3.83805410D+00, 4.09238656D+00, FeO
     6  4.47428992D+00, 4.79318896D+00, 5.01514143D+00, 5.25898668D+00, FeO
     7  5.68213308D+00, 6.00605915D+00, 6.14470591D+00, 6.23435069D+00, FeO
     8       5*0.0D+00/                                                 FeO
      DATA TQ_YO/                                                       071215
     1 -0.999999993529,-0.827499974905,-0.709400049607,-0.555099995812, YO
     2 -0.458399969345,-0.360800023198,-0.217000003967,-0.082999979785, YO
     3  0.085500030985, 0.263799985715, 0.473900012793, 0.691000068209, YO
     4  1.137000118292, 1.948000058536, 2.229600015344, 2.361599903120, YO
     5  2.480900007808, 2.706900052853, 2.923900075088, 3.253200072018, YO
     6  3.487799816033, 3.638899873155, 3.806000043343, 3.920800079070, YO
     7  3.969199793014, 3.987499842816, 4.000000000000,      6*0.0D+00/ YO
      DATA  Q_YO/                                                       071215
     1  3.01048879D-01, 3.01757571D-01, 3.05340616D-01, 3.24274311D-01, YO
     2  3.51192835D-01, 3.92665019D-01, 4.76704361D-01, 5.72455822D-01, YO
     3  7.07994309D-01, 8.63322447D-01, 1.05626799D+00, 1.26260268D+00, YO
     4  1.69801282D+00, 2.50412822D+00, 2.78570794D+00, 2.91942675D+00, YO
     5  3.04421443D+00, 3.30369838D+00, 3.59519633D+00, 4.11879209D+00, YO
     6  4.53879552D+00, 4.82852457D+00, 5.17994373D+00, 5.45286438D+00, YO
     7  5.57723896D+00, 5.62569685D+00, 5.65923907D+00,      6*0.0D+00/ YO
      DATA TQ_ZrO/                                                      071215
     1 -0.999999993529,-0.889399972135,-0.762599990938,-0.647099980967, ZrO
     2 -0.499100003942,-0.347300014887,-0.178199990717,-0.033099988066, ZrO
     3  0.120299999317, 0.295799967424, 0.483200015204, 0.935199992527, ZrO
     4  1.403400088219, 1.861999903043, 2.100099916565, 2.218700015606, ZrO
     5  2.335399924317, 2.443200057519, 2.543599957400, 2.685600050967, ZrO
     6  2.757799999576, 2.824899942069, 3.029100173107, 3.142199903402, ZrO
     7  3.268100022890, 3.403200008008, 3.533200121779, 3.746699913777, ZrO
     8  3.827799979459, 3.900300092509, 3.961400008126, 3.985199903537, ZrO
     9  4.000000000000/                                                 ZrO
      DATA  Q_ZrO/                                                      071215
     1  7.04188250D-06, 1.07606270D-04, 1.16187693D-03, 5.95416748D-03, ZrO
     2  2.74940009D-02, 8.04379883D-02, 1.76971099D-01, 2.81412128D-01, ZrO
     3  4.05246392D-01, 5.57954523D-01, 7.29493672D-01, 1.16265753D+00, ZrO
     4  1.62413045D+00, 2.08051793D+00, 2.31817769D+00, 2.43682377D+00, ZrO
     5  2.55477813D+00, 2.66810315D+00, 2.78388474D+00, 2.98263945D+00, ZrO
     6  3.10600472D+00, 3.23514548D+00, 3.69286105D+00, 3.96690647D+00, ZrO
     7  4.27444233D+00, 4.60249299D+00, 4.91756604D+00, 5.45244762D+00, ZrO
     8  5.66883315D+00, 5.87054537D+00, 6.04664824D+00, 6.11667945D+00, ZrO
     9  6.16061629D+00/                                                 ZrO
      DATA TQ_LaO/                                                      071215
     1 -0.999999993529,-0.900099970436,-0.797800008298,-0.572400015689, LaO
     2 -0.490900013260,-0.411600025447,-0.223400012244,-0.075099986970, LaO
     3  0.077800041181, 0.446700048500, 0.974400008393, 1.487400001413, LaO
     4  2.086799926977, 2.309699942410, 2.553199951241, 2.694800045115, LaO
     5  2.835299919577, 3.158999857229, 3.339400149246, 3.505100108972, LaO
     6  3.643799938955, 3.778500169494, 3.923100015562, 3.969699779224, LaO
     7  4.000000000000,      8*0.0D+00/                                 LaO
      DATA  Q_LaO/                                                      071215
     1  3.01082140D-01, 3.01447532D-01, 3.03284517D-01, 3.29731496D-01, LaO
     2  3.54421364D-01, 3.88119264D-01, 5.00468249D-01, 6.10755136D-01, LaO
     3  7.36671170D-01, 1.06969404D+00, 1.57882080D+00, 2.08644650D+00, LaO
     4  2.68420724D+00, 2.90840011D+00, 3.16767321D+00, 3.33604871D+00, LaO
     5  3.52117506D+00, 4.01747820D+00, 4.33409559D+00, 4.66090751D+00, LaO
     6  4.97506049D+00, 5.32135273D+00, 5.73163338D+00, 5.87007813D+00, LaO
     7  5.96131674D+00,      8*0.0D+00/                                 LaO
C
C Molecular equilibrium constants
C
      DATA TK_H2p/                                                      071215
     1 -0.999999993529,-0.991400000460,-0.983400010010,-0.956499982896, H2p
     2 -0.888599975419,-0.809400015499,-0.716600027498,-0.591400047269, H2p
     3 -0.477400032090,-0.335499965282,-0.209199980744,-0.017799976264, H2p
     4  0.173299943268, 0.377899979937, 0.584100041272, 0.981499998074, H2p
     5  1.075399927923, 1.163300089041, 1.375800007358, 1.531899985759, H2p
     6  1.729100022745, 1.876100034930, 2.028199977024, 2.239800018575, H2p
     7  2.470700025935, 2.685100051695, 2.899700106976, 3.057399996821, H2p
     8  3.298700161091, 3.414099937677, 3.522500032115, 3.662199900008, H2p
     9  3.770700011920, 3.911299877624, 3.965199903328, 4.000000000000, H2p
     A      10*0.0D+00/                                                 H2p
      DATA  K_H2p/                                                      071215
     1  4.42957541D+00, 4.37932049D+00, 4.33377025D+00, 4.18941190D+00, H2p
     2  3.87980381D+00, 3.60560575D+00, 3.38348484D+00, 3.21906034D+00, H2p
     3  3.17418592D+00, 3.22262666D+00, 3.33835355D+00, 3.60562700D+00, H2p
     4  3.94760606D+00, 4.36682517D+00, 4.82436140D+00, 5.76048333D+00, H2p
     5  5.98692326D+00, 6.19817292D+00, 6.68782242D+00, 7.01093782D+00, H2p
     6  7.37659784D+00, 7.62976820D+00, 7.88109111D+00, 8.21896898D+00, H2p
     7  8.57758262D+00, 8.90363659D+00, 9.21775276D+00, 9.42959734D+00, H2p
     8  9.69819905D+00, 9.79500145D+00, 9.86528682D+00, 9.93436788D+00, H2p
     9  9.98299990D+00, 1.00535387D+01, 1.00846296D+01, 1.01060842D+01, H2p
     A      10*0.0D+00/                                                 H2p
      DATA TK_H2/                                                       071215
     1 -0.999999993529,-0.991400000460,-0.983500009882,-0.956699982573, H2
     2 -0.888999973777,-0.809900015952,-0.717200025150,-0.592300043647, H2
     3 -0.478800029520,-0.337599962672,-0.212099987400,-0.023099970609, H2
     4  0.164899963590, 0.362300019275, 0.560599939736, 0.880100005528, H2
     5  1.185800048620, 1.292999951640, 1.401400091921, 1.646600095219, H2
     6  1.739500021561, 1.833999911366, 2.004400001305, 2.191400057447, H2
     7  2.431400083629, 2.686200050094, 2.908900093694, 3.125400005824, H2
     8  3.280900195092, 3.492599817367, 3.608200085128, 3.725099921561, H2
     9  3.881200088658, 3.949899830504, 4.000000000000,     11*0.0D+00/ H2
      DATA  K_H2/                                                       071215
     1  4.33236814D+00, 4.29581929D+00, 4.26319598D+00, 4.15958556D+00, H2
     2  3.94191814D+00, 3.75775385D+00, 3.62205794D+00, 3.54834328D+00, H2
     3  3.56570160D+00, 3.67113883D+00, 3.82317647D+00, 4.12537323D+00, H2
     4  4.48562853D+00, 4.90513083D+00, 5.35374193D+00, 6.11006816D+00, H2
     5  6.85457398D+00, 7.11814881D+00, 7.38453538D+00, 7.96695243D+00, H2
     6  8.16995327D+00, 8.36333109D+00, 8.68395976D+00, 9.00870885D+00, H2
     7  9.40191577D+00, 9.80250704D+00, 1.01435416D+01, 1.04644833D+01, H2
     8  1.06793995D+01, 1.09311192D+01, 1.10407332D+01, 1.11286512D+01, H2
     9  1.12174582D+01, 1.12512880D+01, 1.12755438D+01,     11*0.0D+00/ H2
      DATA TK_H2m/                                                      071215
     1 -0.999999993529,-0.978400015269,-0.944199997999,-0.856800014247, H2m
     2 -0.753000002171,-0.636399995054,-0.462899988609,-0.294200016665, H2m
     3 -0.119799994646, 0.066500001713, 0.284999972725, 0.512599976409, H2m
     4  1.003599989411, 1.171700065743, 1.335599914899, 1.470600031724, H2m
     5  1.669300068983, 1.810099942983, 1.969000044803, 2.120299887566, H2m
     6  2.285399965426, 2.496199987455, 2.718900031847, 2.924600074984, H2m
     7  3.128700088082, 3.293200042730, 3.502600050689, 3.614000035132, H2m
     8  3.727599977906, 3.882700050503, 3.950799845848, 4.000000000000, H2m
     9      14*0.0D+00/                                                 H2m
      DATA  K_H2m/                                                      071215
     1  1.50040010D+00, 1.51890075D+00, 1.55168026D+00, 1.65288578D+00, H2m
     2  1.80052860D+00, 1.99447676D+00, 2.32390682D+00, 2.67727737D+00, H2m
     3  3.06563496D+00, 3.49775541D+00, 4.01918509D+00, 4.57270183D+00, H2m
     4  5.78498219D+00, 6.20289001D+00, 6.61060644D+00, 6.94364359D+00, H2m
     5  7.41088561D+00, 7.71011795D+00, 8.01529618D+00, 8.28393262D+00, H2m
     6  8.56193844D+00, 8.90231661D+00, 9.25057044D+00, 9.56495081D+00, H2m
     7  9.86699872D+00, 1.00933282D+01, 1.03393396D+01, 1.04435216D+01, H2m
     8  1.05281929D+01, 1.06161252D+01, 1.06496413D+01, 1.06734711D+01, H2m
     9      14*0.0D+00/                                                 H2m
      DATA TK_CH/                                                       071215
     1 -0.999999993529,-0.990900000863,-0.982000011796,-0.952599989197, CH
     2 -0.878600006960,-0.790300012818,-0.688700011135,-0.558699997986, CH
     3 -0.436400004985,-0.307600007686,-0.173599998332,-0.056200010612, CH
     4  0.069300034588, 0.196100045453, 0.323700046438, 0.440000057620, CH
     5  0.549799965386, 0.752999989358, 0.827900047001, 0.904099978008, CH
     6  1.018699974015, 1.128600088053, 1.365399882457, 1.500499991127, CH
     7  1.632900104982, 1.875900029982, 2.051099959430, 2.250400001291, CH
     8  2.453500033961, 2.641500096407, 2.773099982521, 2.904500100204, CH
     9  3.059099954601, 3.193799950014, 3.396300014431, 3.588700129277, CH
     A  3.757700169670, 3.883900019978, 4.000000000000,      7*0.0D+00/ CH
      DATA  K_CH/                                                       071215
     1  4.02977109D+00, 3.98070656D+00, 3.93412321D+00, 3.79009680D+00, CH
     2  3.48851269D+00, 3.22735998D+00, 3.03551446D+00, 2.92292742D+00, CH
     3  2.92091802D+00, 3.00008533D+00, 3.14834979D+00, 3.31947540D+00, CH
     4  3.53427804D+00, 3.77650259D+00, 4.03937489D+00, 4.29173326D+00, CH
     5  4.53955396D+00, 5.02820574D+00, 5.22155157D+00, 5.42579592D+00, CH
     6  5.74273382D+00, 6.04939295D+00, 6.69058645D+00, 7.03556837D+00, CH
     7  7.35395278D+00, 7.88410502D+00, 8.22803201D+00, 8.59001490D+00, CH
     8  8.93618673D+00, 9.24268733D+00, 9.45116432D+00, 9.65365465D+00, CH
     9  9.87998838D+00, 1.00599715D+01, 1.02853995D+01, 1.04396130D+01, CH
     A  1.05277585D+01, 1.05647050D+01, 1.05822400D+01,      7*0.0D+00/ CH
      DATA TK_CHm/                                                      071215
     1 -0.999999993529,-0.987600004651,-0.969400017531,-0.921800024215, CHm
     2 -0.800100007062,-0.642499985523,-0.469000037638,-0.298500022499, CHm
     3 -0.120799994748, 0.070100042786, 0.263099987402, 0.495800010179, CHm
     4  0.682199963322, 0.834400035928, 0.995100000789, 1.081299921486, CHm
     5  1.166400079296, 1.395700093342, 1.525899979120, 1.664900076852, CHm
     6  1.956600058135, 2.247400006455, 2.549699957473, 2.789499978583, CHm
     7  2.968300019587, 3.252200050192, 3.406599918900, 3.583200000418, CHm
     8  3.684999943468, 3.813000012894, 3.919600090882, 3.968999798529, CHm
     9  4.000000000000,     13*0.0D+00/                                 CHm
      DATA  K_CHm/                                                      071215
     1  5.55704089D-03, 6.92785550D-02, 1.61127487D-01, 3.92568897D-01, CHm
     2  9.34045460D-01, 1.55126868D+00, 2.15306760D+00, 2.69043010D+00, CHm
     3  3.21228080D+00, 3.74412547D+00, 4.26210892D+00, 4.87016071D+00, CHm
     4  5.34884501D+00, 5.73346535D+00, 6.12327029D+00, 6.31928852D+00, CHm
     5  6.50149310D+00, 6.94339411D+00, 7.17280993D+00, 7.40740486D+00, CHm
     6  7.87807350D+00, 8.33065941D+00, 8.79208352D+00, 9.15354869D+00, CHm
     7  9.41609197D+00, 9.78936123D+00, 9.95325717D+00, 1.01086309D+01, CHm
     8  1.01886032D+01, 1.02788830D+01, 1.03408394D+01, 1.03652193D+01, CHm
     9  1.03793960D+01,     13*0.0D+00/                                 CHm
      DATA TK_C2/                                                       071215
     1 -0.999999993529,-0.990600001105,-0.980900013200,-0.949999993397, C2
     2 -0.872099989476,-0.781599976096,-0.669999960957,-0.547699983005, C2
     3 -0.417700013048,-0.273199967955,-0.119899994406, 0.018100005497, C2
     4  0.137399988246, 0.249800003312, 0.367999995295, 0.472400012403, C2
     5  0.569800053940, 0.779599985758, 0.867800014641, 0.961500049686, C2
     6  1.070899935530, 1.181900052794, 1.340499899990, 1.482800009472, C2
     7  1.614199893579, 1.731700021012, 1.920900074130, 2.035399971564, C2
     8  2.151600090822, 2.281999979285, 2.390400098907, 2.478900010907, C2
     9  2.596099913148, 2.698000045481, 2.833099929055, 2.973000014571, C2
     A  3.221999828736, 3.340800144318, 3.459200016066, 3.675099932070, C2
     B  3.796800120618, 3.892699928365, 3.956999981229, 4.000000000000, C2
     C       2*0.0D+00/                                                 C2
      DATA  K_C2/                                                       071215
     1  8.60089074D+00, 8.49126962D+00, 8.38109698D+00, 8.04969052D+00, C2
     2  7.33445382D+00, 6.68944017D+00, 6.11627467D+00, 5.70822557D+00, C2
     3  5.46551083D+00, 5.36623586D+00, 5.40130166D+00, 5.51586428D+00, C2
     4  5.65126716D+00, 5.79252350D+00, 5.94637170D+00, 6.08566153D+00, C2
     5  6.22091373D+00, 6.55819752D+00, 6.73413011D+00, 6.94942730D+00, C2
     6  7.23590242D+00, 7.55750335D+00, 8.04811920D+00, 8.49140179D+00, C2
     7  8.88243789D+00, 9.20862477D+00, 9.68305270D+00, 9.94152289D+00, C2
     8  1.01816063D+01, 1.04161034D+01, 1.05722887D+01, 1.06713587D+01, C2
     9  1.07721964D+01, 1.08449002D+01, 1.09359487D+01, 1.10311771D+01, C2
     A  1.11981583D+01, 1.12707487D+01, 1.13355157D+01, 1.14358655D+01, C2
     B  1.14859647D+01, 1.15224946D+01, 1.15448674D+01, 1.15587454D+01, C2
     C       2*0.0D+00/                                                 C2
      DATA TK_C2m/                                                      071215
     1 -0.999999993529,-0.990300001346,-0.979700014521,-0.947299995539, C2m
     2 -0.912999988981,-0.865099991479,-0.770699970589,-0.653799972127, C2m
     3 -0.522599970833,-0.386299989475,-0.239300005780,-0.070599990758, C2m
     4  0.041199962653, 0.158699976983, 0.265799980894, 0.378599979322, C2m
     5  0.482500015023, 0.578600053398, 0.784299971545, 0.870200012391, C2m
     6  0.961500049686, 1.072799932318, 1.176300059693, 1.444000043153, C2m
     7  1.546499958053, 1.649700097742, 1.865499894721, 2.167700063949, C2m
     8  2.444300054382, 2.566499946135, 2.683600053877, 2.927900074497, C2m
     9  3.041800017332, 3.158799862433, 3.257400163686, 3.339900161718, C2m
     A  3.528200159065, 3.604700006665, 3.677299877217, 3.832499981182, C2m
     B  3.931999868201, 3.973399862706, 4.000000000000,      3*0.0D+00/ C2m
      DATA  K_C2m/                                                      071215
     1  1.11491436D+01, 1.09864294D+01, 1.08133169D+01, 1.03135505D+01, C2m
     2  9.82990992D+00, 9.22612571D+00, 8.24990949D+00, 7.36393460D+00, C2m
     3  6.69735970D+00, 6.27851153D+00, 6.05417796D+00, 5.99977828D+00, C2m
     4  6.04328044D+00, 6.12931775D+00, 6.22808475D+00, 6.34474699D+00, C2m
     5  6.46153608D+00, 6.57780268D+00, 6.86348186D+00, 7.00353220D+00, C2m
     6  7.16792194D+00, 7.38839560D+00, 7.60932783D+00, 8.21460104D+00, C2m
     7  8.44483871D+00, 8.66980601D+00, 9.11193829D+00, 9.67094043D+00, C2m
     8  1.01392152D+01, 1.03371856D+01, 1.05219855D+01, 1.08831906D+01, C2m
     9  1.10346360D+01, 1.11763794D+01, 1.12849916D+01, 1.13693144D+01, C2m
     A  1.15509595D+01, 1.16255091D+01, 1.16970620D+01, 1.18378493D+01, C2m
     B  1.19055943D+01, 1.19278193D+01, 1.19406334D+01,      3*0.0D+00/ C2m
      DATA TK_CN/                                                       071215
     1 -0.999999993529,-0.990400001266,-0.980100014221,-0.948199994825, CN
     2 -0.867299988043,-0.774099970032,-0.659099963805,-0.531399952758, CN
     3 -0.396399976892,-0.249199995297,-0.089699972952, 0.046499957859, CN
     4  0.168699952698, 0.284999972725, 0.400399985726, 0.588600028503, CN
     5  0.762599991353, 0.856500022827, 0.955300041687, 1.070199936713, CN
     6  1.177500058115, 1.447600043666, 1.554799943267, 1.662600080966, CN
     7  1.877700074512, 2.245100010241, 2.560299938367, 2.761499999576, CN
     8  3.014700116707, 3.271800019829, 3.433299895367, 3.555700049167, CN
     9  3.632000066109, 3.705399930656, 3.841300111482, 3.936499965315, CN
     A  3.975299913981, 4.000000000000,      8*0.0D+00/                 CN
      DATA  K_CN/                                                       071215
     1  1.02101693D+01, 1.00641137D+01, 9.91151171D+00, 9.46502180D+00, CN
     2  8.49380166D+00, 7.61993629D+00, 6.83168541D+00, 6.24697351D+00, CN
     3  5.87794835D+00, 5.68877032D+00, 5.66156783D+00, 5.73652874D+00, CN
     4  5.84890733D+00, 5.97551013D+00, 6.11057162D+00, 6.34765859D+00, CN
     5  6.59638501D+00, 6.74980496D+00, 6.92871089D+00, 7.15809139D+00, CN
     6  7.38891656D+00, 8.00283807D+00, 8.24407517D+00, 8.47881701D+00, CN
     7  8.91810988D+00, 9.58841295D+00, 1.01110076D+01, 1.04279209D+01, CN
     8  1.07970904D+01, 1.11123519D+01, 1.12700197D+01, 1.13629817D+01, CN
     9  1.14083493D+01, 1.14442087D+01, 1.15004857D+01, 1.15421834D+01, CN
     A  1.15612532D+01, 1.15740938D+01,      8*0.0D+00/                 CN
      DATA TK_CNm/                                                      071215
     1 -0.999999993529,-0.990300001346,-0.980000014348,-0.947799995143, CNm
     2 -0.913699993008,-0.866299989605,-0.772899970229,-0.712400043932, CNm
     3 -0.656199968359,-0.589500052342,-0.523099969637,-0.455299977181, CNm
     4 -0.382999984914,-0.310600002039,-0.235500006941,-0.081099981723, CNm
     5  0.022399998361, 0.118600002059, 0.204600035071, 0.290599972628, CNm
     6  0.451600043460, 0.700000037445, 0.931399955477, 1.176100059956, CNm
     7  1.452500039561, 1.745300003357, 2.132700101222, 2.497399985635, CNm
     8  2.611799898837, 2.728500006353, 2.946800047977, 3.143599937140, CNm
     9  3.323600068314, 3.484399907751, 3.629600113141, 3.800200182895, CNm
     A  3.911199875055, 3.964899911601, 4.000000000000,      7*0.0D+00/ CNm
      DATA  K_CNm/                                                      071215
     1  1.33743325D+01, 1.31701844D+01, 1.29589015D+01, 1.23336642D+01, CNm
     2  1.17264408D+01, 1.09683815D+01, 9.73102303D+00, 9.08511692D+00, CNm
     3  8.57828484D+00, 8.07793708D+00, 7.67419434D+00, 7.34520245D+00, CNm
     4  7.07349346D+00, 6.87092098D+00, 6.72282616D+00, 6.57496331D+00, CNm
     5  6.56389283D+00, 6.59635943D+00, 6.64860348D+00, 6.71515833D+00, CNm
     6  6.86515877D+00, 7.14388191D+00, 7.43935385D+00, 7.77447431D+00, CNm
     7  8.16909938D+00, 8.59725601D+00, 9.17157344D+00, 9.71590667D+00, CNm
     8  9.88674025D+00, 1.00599422D+01, 1.03725852D+01, 1.06267757D+01, CNm
     9  1.08267819D+01, 1.09772045D+01, 1.10929432D+01, 1.12192075D+01, CNm
     A  1.13074343D+01, 1.13534987D+01, 1.13846135D+01,      7*0.0D+00/ CNm
      DATA TK_NH/                                                       071215
     1 -0.999999993529,-0.991500000379,-0.983700009627,-0.957199981765, NH
     2 -0.890199969686,-0.812200008729,-0.720800013728,-0.597600022321, NH
     3 -0.478000030989,-0.351700033383,-0.223700012079,-0.119399995606, NH
     4 -0.009399998486, 0.211900027670, 0.427799950658, 0.615800018765, NH
     5  0.731800033667, 0.848600024735, 1.021399976632, 1.161800093757, NH
     6  1.285699952987, 1.428100065299, 1.635800098800, 1.856799907572, NH
     7  2.075799935240, 2.540399957361, 2.718900031847, 2.888900110077, NH
     8  3.175399944572, 3.428199869527, 3.596899981415, 3.745799893758, NH
     9  3.880100116638, 4.000000000000,     12*0.0D+00/                 NH
      DATA  K_NH/                                                       071215
     1  4.77476987D+00, 4.72869744D+00, 4.68750133D+00, 4.55565673D+00, NH
     2  4.27300885D+00, 4.02437234D+00, 3.82491863D+00, 3.68167067D+00, NH
     3  3.64949990D+00, 3.70165802D+00, 3.82230166D+00, 3.95887877D+00, NH
     4  4.13152537D+00, 4.54301035D+00, 4.99919698D+00, 5.42359029D+00, NH
     5  5.69366673D+00, 5.96904126D+00, 6.37205603D+00, 6.68027906D+00, NH
     6  6.93017176D+00, 7.19477406D+00, 7.55312604D+00, 7.91389442D+00, NH
     7  8.25945080D+00, 8.97253963D+00, 9.24250159D+00, 9.49713005D+00, NH
     8  9.90392100D+00, 1.02052687D+01, 1.03624382D+01, 1.04671445D+01, NH
     9  1.05340118D+01, 1.05784579D+01,     12*0.0D+00/                 NH
      DATA TK_N2/                                                       071215
     1 -0.999999993529,-0.990300001346,-0.980000014348,-0.947899995063, N2
     2 -0.913899994159,-0.866599989136,-0.773499970130,-0.656999967102, N2
     3 -0.523799967963,-0.383199985190,-0.236000006788,-0.083899978867, N2
     4  0.023999994474, 0.121300000066, 0.207100031290, 0.292799970426, N2
     5  0.452300043220, 0.690900068551, 0.923099937971, 1.162100092813, N2
     6  1.445700043395, 1.729500021922, 2.196100048306, 2.528399986098, N2
     7  2.659500092750, 2.780599988005, 2.973200014527, 3.112500081155, N2
     8  3.298000146027, 3.483199940123, 3.632200060516, 3.796600116425, N2
     9  3.854600001832, 3.911799890471, 3.965099906085, 4.000000000000, N2
     A      10*0.0D+00/                                                 N2
      DATA  K_N2/                                                       071215
     1  1.34746372D+01, 1.32827195D+01, 1.30841292D+01, 1.24984314D+01, N2
     2  1.19297358D+01, 1.12196399D+01, 1.00633431D+01, 8.98761003D+00, N2
     3  8.14708875D+00, 7.59438415D+00, 7.28237738D+00, 7.16349607D+00, N2
     4  7.16755303D+00, 7.21382676D+00, 7.27663646D+00, 7.35231603D+00, N2
     5  7.51431967D+00, 7.79329659D+00, 8.09531427D+00, 8.42571903D+00, N2
     6  8.83262059D+00, 9.24863375D+00, 9.94152729D+00, 1.04378768D+01, N2
     7  1.06337241D+01, 1.08135299D+01, 1.10909066D+01, 1.12775885D+01, N2
     8  1.14992813D+01, 1.16875532D+01, 1.18197869D+01, 1.19661453D+01, N2
     9  1.20253559D+01, 1.20906142D+01, 1.21581210D+01, 1.22054183D+01, N2
     A      10*0.0D+00/                                                 N2
      DATA TK_OH/                                                       071215
     1 -0.999999993529,-0.990800000943,-0.981500012434,-0.951499990974, OH
     2 -0.875799999429,-0.790100012938,-0.681400001535,-0.558299997745, OH
     3 -0.426200015717,-0.291600013138,-0.154099995516,-0.043200020973, OH
     4  0.072900042202, 0.305799936659, 0.529699969696, 0.712500045598, OH
     5  0.854000023317, 0.990300001562, 1.142800111339, 1.292299950482, OH
     6  1.466600030583, 1.635700099013, 1.822799933694, 1.992200007868, OH
     7  2.153100090310, 2.386500108353, 2.604499906290, 2.869799893075, OH
     8  3.032400140381, 3.231399995262, 3.403899989662, 3.618999922198, OH
     9  3.774300084646, 3.905999946609, 3.963399952969, 4.000000000000, OH
     A      10*0.0D+00/                                                 OH
      DATA  K_OH/                                                       071215
     1  5.68366610D+00, 5.61432226D+00, 5.54615092D+00, 5.33932395D+00, OH
     2  4.89822063D+00, 4.51964954D+00, 4.18967043D+00, 3.97520539D+00, OH
     3  3.88751520D+00, 3.91070803D+00, 4.02043777D+00, 4.15609229D+00, OH
     4  4.33216948D+00, 4.75978718D+00, 5.23175789D+00, 5.64401033D+00, OH
     5  5.97330651D+00, 6.29266346D+00, 6.63896324D+00, 6.95194130D+00, OH
     6  7.28115437D+00, 7.57708279D+00, 7.89437761D+00, 8.18431606D+00, OH
     7  8.46601672D+00, 8.87838437D+00, 9.25573262D+00, 9.69668389D+00, OH
     8  9.95460708D+00, 1.02476844D+01, 1.04702825D+01, 1.06973378D+01, OH
     9  1.08250521D+01, 1.09068594D+01, 1.09341997D+01, 1.09492183D+01, OH
     A      10*0.0D+00/                                                 OH
      DATA TK_OHm/                                                      071215
     1 -0.999999993529,-0.987100005289,-0.967600009816,-0.917100012569, OHm
     2 -0.788200005091,-0.618500007546,-0.433400012125,-0.254299995626, OHm
     3 -0.068299994147, 0.113400013077, 0.311699938886, 0.539799952919, OHm
     4  0.775399984370, 0.868600013871, 0.955300041687, 1.149100098996, OHm
     5  1.255500003415, 1.364199882870, 1.521899970516, 1.745600002324, OHm
     6  1.906400077947, 2.060299937992, 2.468300027048, 2.855599898246, OHm
     7  3.025100095556, 3.300300181648, 3.495399882715, 3.708700005790, OHm
     8  3.891899909535, 3.957599994330, 4.000000000000,     15*0.0D+00/ OHm
      DATA  K_OHm/                                                      071215
     1  1.10756925D-02, 9.51817208D-02, 2.19424113D-01, 5.26154397D-01, OHm
     2  1.22427760D+00, 2.00036078D+00, 2.71848917D+00, 3.32863107D+00, OHm
     3  3.90446496D+00, 4.42959329D+00, 4.97535625D+00, 5.58091640D+00, OHm
     4  6.19117527D+00, 6.42922809D+00, 6.64801597D+00, 7.11602988D+00, OHm
     5  7.35187602D+00, 7.57475020D+00, 7.87197246D+00, 8.26238453D+00, OHm
     6  8.53406257D+00, 8.79436738D+00, 9.48825001D+00, 1.01248750D+01, OHm
     7  1.03918757D+01, 1.07924745D+01, 1.10359089D+01, 1.12590847D+01, OHm
     8  1.14201445D+01, 1.14723473D+01, 1.15047160D+01,     15*0.0D+00/ OHm
      DATA TK_BO/                                                       071215
     1 -0.999999993529,-0.990400001266,-0.980100014221,-0.948099994905, BO
     2 -0.867099988356,-0.774199970016,-0.658099965375,-0.526299961984, BO
     3 -0.388999993208,-0.244100000990,-0.087699974992, 0.085400031162, BO
     4  0.165399962157, 0.245300002688, 0.397699986136, 0.552599956974, BO
     5  0.743200010583, 0.829400043182, 0.919699937061, 1.145700105657, BO
     6  1.261299988237, 1.388000093224, 1.600399891505, 1.801099942058, BO
     7  1.934900049163, 2.062999937730, 2.415800099353, 2.627500046092, BO
     8  2.913400086490, 3.024700087800, 3.144699963648, 3.329099936361, BO
     9  3.513000149176, 3.637699906712, 3.763800138147, 3.903300015719, BO
     A  3.960900021915, 4.000000000000,      8*0.0D+00/                 BO
      DATA  K_BO/                                                       071215
     1  1.11021111D+01, 1.09452234D+01, 1.07812607D+01, 1.02998095D+01, BO
     2  9.25297587D+00, 8.31179477D+00, 7.44795714D+00, 6.79175681D+00, BO
     3  6.38148738D+00, 6.16867019D+00, 6.11834004D+00, 6.20147819D+00, BO
     4  6.26819651D+00, 6.34453716D+00, 6.50691095D+00, 6.69031328D+00, BO
     5  6.94689994D+00, 7.07794822D+00, 7.22643889D+00, 7.63668566D+00, BO
     6  7.85566812D+00, 8.09384143D+00, 8.48032353D+00, 8.83340554D+00, BO
     7  9.06843804D+00, 9.29626286D+00, 9.92898700D+00, 1.02970279D+01, BO
     8  1.07574345D+01, 1.09176042D+01, 1.10748930D+01, 1.12845249D+01, BO
     9  1.14582435D+01, 1.15591560D+01, 1.16483193D+01, 1.17268469D+01, BO
     A  1.17508802D+01, 1.17644413D+01,      8*0.0D+00/                 BO
      DATA TK_CO/                                                       071215
     1 -0.999999993529,-0.990100001508,-0.979000014924,-0.966500005101, CO
     2 -0.945699996809,-0.910399974022,-0.860899998037,-0.763999987110, CO
     3 -0.702000003762,-0.643399984631,-0.573800020482,-0.504899997282, CO
     4 -0.435000008317,-0.360600023402,-0.285399985933,-0.206999981965, CO
     5 -0.031999982334, 0.097900011956, 0.227100024360, 0.345900039562, CO
     6  0.480000014377, 0.620700018192, 0.759199993854, 0.849300024418, CO
     7  0.943900034192, 1.053599953886, 1.167400076152, 1.303499961392, CO
     8  1.438300045845, 1.544299964286, 1.645900094649, 1.853899907367, CO
     9  1.973000032371, 2.087999925751, 2.418800094762, 2.588699921809, CO
     A  2.754999995127, 2.961400037324, 3.171900028914, 3.310999963058, CO
     B  3.441899999173, 3.675399924590, 3.818500133516, 3.932999889782, CO
     C  3.973599868103, 4.000000000000/                                 CO
      DATA  K_CO/                                                       071215
     1  1.39707818D+01, 1.37452953D+01, 1.34991781D+01, 1.32302613D+01, CO
     2  1.28015716D+01, 1.21249144D+01, 1.12755709D+01, 9.90591837D+00, CO
     3  9.20304178D+00, 8.64388806D+00, 8.09455517D+00, 7.65630017D+00, CO
     4  7.30313934D+00, 7.01317804D+00, 6.79617992D+00, 6.63822430D+00, CO
     5  6.47942344D+00, 6.48307752D+00, 6.54692151D+00, 6.63596003D+00, CO
     6  6.75937730D+00, 6.91244189D+00, 7.09148895D+00, 7.22723569D+00, CO
     7  7.38790457D+00, 7.59616561D+00, 7.83249946D+00, 8.13285615D+00, CO
     8  8.43729501D+00, 8.67416288D+00, 8.89571367D+00, 9.33123987D+00, CO
     9  9.57288069D+00, 9.80289783D+00, 1.04429206D+01, 1.07537258D+01, CO
     A  1.10436373D+01, 1.13772684D+01, 1.16745242D+01, 1.18431817D+01, CO
     B  1.19831225D+01, 1.22015165D+01, 1.23273378D+01, 1.24267822D+01, CO
     C  1.24605281D+01, 1.24813192D+01/                                 CO
      DATA TK_NOp/                                                      071215
     1 -0.999999993529,-0.989500002226,-0.982600011031,-0.976800016190, NOp
     2 -0.964199995243,-0.940500000934,-0.901799970656,-0.847200046694, NOp
     3 -0.791000012396,-0.735400018144,-0.670299962119,-0.611100016543, NOp
     4 -0.521199974182,-0.429200019270,-0.341799973300,-0.247899996748, NOp
     5 -0.158200015132,-0.063000003138, 0.024399993502, 0.111300017527, NOp
     6  0.238200005503, 0.393199983789, 0.556499943964, 0.664499954996, NOp
     7  0.764299989339, 0.901499984331, 1.037199964169, 1.214200029052, NOp
     8  1.363299883180, 1.524699976538, 1.705000047655, 1.957800057189, NOp
     9  2.191000058225, 2.424400090052, 2.687100048785, 2.885600108270, NOp
     A  3.124099973419, 3.264000122881, 3.400900068286, 3.555400056924, NOp
     B  3.712799971332, 3.833600005141, 3.938500008477, 3.975899930173, NOp
     C  4.000000000000,      1*0.0D+00/                                 NOp
      DATA  K_NOp/                                                      071215
     1 -1.64367049D+01,-1.59479725D+01,-1.56328271D+01,-1.53716910D+01, NOp
     2 -1.48153424D+01,-1.38092884D+01,-1.22733537D+01,-1.03125818D+01, NOp
     3 -8.51916441D+00,-6.94241603D+00,-5.31428866D+00,-4.01158910D+00, NOp
     4 -2.30842328D+00,-8.51236131D-01, 3.15363375D-01, 1.37700296D+00, NOp
     5  2.24029017D+00, 3.02544190D+00, 3.64788177D+00, 4.18752969D+00, NOp
     6  4.85426844D+00, 5.50847656D+00, 6.05752782D+00, 6.36492502D+00, NOp
     7  6.61998643D+00, 6.93617028D+00, 7.22106136D+00, 7.57459093D+00, NOp
     8  7.88135744D+00, 8.23871968D+00, 8.66650167D+00, 9.29216133D+00, NOp
     9  9.85925787D+00, 1.03860756D+01, 1.09186017D+01, 1.12800477D+01, NOp
     A  1.16598768D+01, 1.18498516D+01, 1.20118851D+01, 1.21696103D+01, NOp
     B  1.23104023D+01, 1.24116984D+01, 1.24975835D+01, 1.25270821D+01, NOp
     C  1.25453351D+01,      1*0.0D+00/                                 NOp
      DATA TK_NO/                                                       071215
     1 -0.999999993529,-0.990700001024,-0.981300012690,-0.950999991782, NO
     2 -0.874499995932,-0.785799994547,-0.675999984199,-0.556299996536, NO
     3 -0.430000020217,-0.294300016801,-0.152699988818, 0.053099945611, NO
     4  0.139899982053, 0.225500025841, 0.342200034627, 0.489900016936, NO
     5  0.635399989040, 0.944700033139, 1.257799996142, 1.411000076437, NO
     6  1.560699938868, 1.690900069295, 1.820899939153, 2.059899938262, NO
     7  2.204500047145, 2.381000123772, 2.527199984718, 2.657000092176, NO
     8  2.961400037324, 3.172900004817, 3.343200084633, 3.500700006393, NO
     9  3.630000122038, 3.762400169843, 3.894499970732, 3.957099983412, NO
     A  4.000000000000,      9*0.0D+00/                                 NO
      DATA  K_NO/                                                       071215
     1  9.72431680D+00, 9.60986725D+00, 9.49709941D+00, 9.15325840D+00, NO
     2  8.40666667D+00, 7.72895030D+00, 7.11667182D+00, 6.67328810D+00, NO
     3  6.39850838D+00, 6.26813904D+00, 6.26464267D+00, 6.41341938D+00, NO
     4  6.50579144D+00, 6.60441808D+00, 6.74549161D+00, 6.93224526D+00, NO
     5  7.12452642D+00, 7.55403551D+00, 8.00623773D+00, 8.23063881D+00, NO
     6  8.44908253D+00, 8.63557996D+00, 8.81747109D+00, 9.14920085D+00, NO
     7  9.35554863D+00, 9.61532658D+00, 9.83489750D+00, 1.00306107D+01, NO
     8  1.04703995D+01, 1.07353473D+01, 1.09162544D+01, 1.10592772D+01, NO
     9  1.11632143D+01, 1.12652956D+01, 1.13753617D+01, 1.14330396D+01, NO
     A  1.14742638D+01,      9*0.0D+00/                                 NO
      DATA TK_O2/                                                       071215
     1 -0.999999993529,-0.991200000621,-0.982900010648,-0.954999985319, O2
     2 -0.884699991430,-0.801600008423,-0.704800021108,-0.579500039996, O2
     3 -0.464299999861,-0.331499970252,-0.183499984477, 0.000300000332, O2
     4  0.080700039483, 0.161299973909, 0.272299971104, 0.411399945910, O2
     5  0.572500055544, 0.884499992028, 1.285999952547, 1.536399980598, O2
     6  1.760499979463, 1.909700075943, 2.060099938011, 2.295099955114, O2
     7  2.431500083432, 2.607499905177, 2.863099896887, 3.085900059538, O2
     8  3.285600080832, 3.475799924916, 3.673999959497, 3.769700004571, O2
     9  3.885299984367, 3.953299900437, 3.982299980098, 4.000000000000, O2
     A      10*0.0D+00/                                                 O2
      DATA  K_O2/                                                       071215
     1  8.61660370D+00, 8.53580973D+00, 8.46141910D+00, 8.22454866D+00, O2
     2  7.70961049D+00, 7.23279409D+00, 6.82619517D+00, 6.48840355D+00, O2
     3  6.32284332D+00, 6.26129941D+00, 6.31146238D+00, 6.48266172D+00, O2
     4  6.57707530D+00, 6.67667794D+00, 6.81797441D+00, 7.00012641D+00, O2
     5  7.21825772D+00, 7.65814173D+00, 8.24378870D+00, 8.61523650D+00, O2
     6  8.95868901D+00, 9.20486340D+00, 9.47385750D+00, 9.92519152D+00, O2
     7  1.01910953D+01, 1.05246696D+01, 1.09668328D+01, 1.12896317D+01, O2
     8  1.15232565D+01, 1.16978452D+01, 1.18303029D+01, 1.18778087D+01, O2
     9  1.19199146D+01, 1.19332728D+01, 1.19354743D+01, 1.19356695D+01, O2
     A      10*0.0D+00/                                                 O2
      DATA TK_HF/                                                       071215
     1 -0.999999993529,-0.990500001185,-0.980700013455,-0.949499993794, HF
     2 -0.870599985441,-0.780599971703,-0.668299961201,-0.542099959324, HF
     3 -0.405199999053,-0.258999996962,-0.109700017992, 0.005900006526, HF
     4  0.126300003808, 0.244200002535, 0.364700009178, 0.566200009252, HF
     5  0.739700022037, 0.870500012190, 0.995100000789, 1.160000099415, HF
     6  1.306499959824, 1.511299970683, 1.667800071666, 1.944000048931, HF
     7  2.203500045717, 2.368499890651, 2.531599983045, 2.953500045827, HF
     8  3.232699961904, 3.420900068406, 3.625600024170, 3.783500122187, HF
     9  3.911899893041, 3.965599892296, 3.986499869216, 4.000000000000, HF
     A      10*0.0D+00/                                                 HF
      DATA  K_HF/                                                       071215
     1  7.66126888D+00, 7.55796589D+00, 7.45424962D+00, 7.14285495D+00, HF
     2  6.47160169D+00, 5.88064578D+00, 5.35400024D+00, 4.97952484D+00, HF
     3  4.76855557D+00, 4.70733515D+00, 4.77041952D+00, 4.88273991D+00, HF
     4  5.04347026D+00, 5.23383384D+00, 5.45399209D+00, 5.86360490D+00, HF
     5  6.24481290D+00, 6.54357990D+00, 6.83251621D+00, 7.20802329D+00, HF
     6  7.51847177D+00, 7.90707256D+00, 8.17942008D+00, 8.63301803D+00, HF
     7  9.04675857D+00, 9.31342818D+00, 9.58269693D+00, 1.02831435D+01, HF
     8  1.07167382D+01, 1.09710589D+01, 1.11985078D+01, 1.13351075D+01, HF
     9  1.14174002D+01, 1.14431605D+01, 1.14518115D+01, 1.14570056D+01, HF
     A      10*0.0D+00/                                                 HF
      DATA TK_NaH/                                                      071215
     1 -0.999999993529,-0.986700005799,-0.966500005101,-0.914299996460, NaH
     2 -0.780899973021,-0.604600015064,-0.415900016707,-0.269499974520, NaH
     3 -0.119799994646, 0.138299986017, 0.236200009447, 0.337800035973, NaH
     4  0.543499957153, 0.645699966249, 0.747799994753, 0.936200002277, NaH
     5  1.144900107225, 1.498799992732, 1.868799886875, 2.246200008430, NaH
     6  2.448700041837, 2.736000010293, 2.856599898338, 2.978500013386, NaH
     7  3.223599869226, 3.461000010341, 3.558699971603, 3.659399869584, NaH
     8  3.741399795885, 3.817300107199, 3.866999964921, 3.915599988107, NaH
     9  3.969499784740, 3.987699837536, 4.000000000000,     11*0.0D+00/ NaH
      DATA  K_NaH/                                                      071215
     1  3.37494118D+00, 3.35142895D+00, 3.31900697D+00, 3.25223211D+00, NaH
     2  3.17722105D+00, 3.23889858D+00, 3.44385537D+00, 3.66956198D+00, NaH
     3  3.94177399D+00, 4.47588697D+00, 4.69266550D+00, 4.92210449D+00, NaH
     4  5.38219182D+00, 5.59717906D+00, 5.79861828D+00, 6.14091127D+00, NaH
     5  6.49282797D+00, 7.05734273D+00, 7.62706556D+00, 8.19869525D+00, NaH
     6  8.50194761D+00, 8.91017365D+00, 9.06350910D+00, 9.20323319D+00, NaH
     7  9.43371916D+00, 9.58792910D+00, 9.62858467D+00, 9.66027575D+00, NaH
     8  9.69077261D+00, 9.75023563D+00, 9.82735365D+00, 9.94414048D+00, NaH
     9  1.01180341D+01, 1.01841772D+01, 1.02302142D+01,     11*0.0D+00/ NaH
      DATA TK_MgH/                                                      071215
     1 -0.999999993529,-0.984600008479,-0.959699977726,-0.896699970175, MgH
     2 -0.824399978168,-0.738000022846,-0.624100006491,-0.521399973703, MgH
     3 -0.394099983253,-0.280199957634,-0.129400001006, 0.008800009733, MgH
     4  0.262799988126, 0.426699948959, 0.602000024872, 0.687600038302, MgH
     5  0.773199983643, 0.983699998989, 1.196600038429, 1.468300031707, MgH
     6  1.744400006454, 2.000800000237, 2.283099974801, 2.478900010907, MgH
     7  2.668500083544, 2.810499923431, 2.961300037581, 3.122399931044, MgH
     8  3.289199993314, 3.488599794452, 3.611400093858, 3.735099895698, MgH
     9  3.794200066108, 3.855600026090, 3.900800079711, 3.943299970733, MgH
     A  3.977899984147, 3.989799782095, 4.000000000000,      7*0.0D+00/ MgH
      DATA  K_MgH/                                                      071215
     1  2.23575980D+00, 2.22746378D+00, 2.21749541D+00, 2.20971928D+00, MgH
     2  2.22810950D+00, 2.28240641D+00, 2.39768109D+00, 2.53541092D+00, MgH
     3  2.74031112D+00, 2.94828850D+00, 3.25018333D+00, 3.54637404D+00, MgH
     4  4.12293766D+00, 4.50709337D+00, 4.91075740D+00, 5.09728171D+00, MgH
     5  5.27388215D+00, 5.66981887D+00, 6.03499720D+00, 6.47494352D+00, MgH
     6  6.90618740D+00, 7.29909281D+00, 7.72702753D+00, 8.02156385D+00, MgH
     7  8.30117837D+00, 8.50001995D+00, 8.69331251D+00, 8.87350040D+00, MgH
     8  9.02862692D+00, 9.16929546D+00, 9.22805911D+00, 9.26804345D+00, MgH
     9  9.28551172D+00, 9.30882172D+00, 9.33438244D+00, 9.37043398D+00, MgH
     A  9.41271788D+00, 9.43062732D+00, 9.44752211D+00,      7*0.0D+00/ MgH
      DATA TK_MgO/                                                      071215
     1 -0.999999993529,-0.991400000460,-0.983600009755,-0.956899982250, MgO
     2 -0.889399972135,-0.814300001748,-0.718700019280,-0.588300050999, MgO
     3 -0.453399981984,-0.340099960446,-0.229100009116,-0.131299999306, MgO
     4 -0.032099982855, 0.234000013786, 0.421499940926, 0.623300014326, MgO
     5  0.957100046960, 1.312599949346, 1.449700043965, 1.576799927543, MgO
     6  1.878100084408, 1.992300007767, 2.117199881937, 2.259199981848, MgO
     7  2.412800103945, 2.559299939355, 2.681600056786, 2.840599899231, MgO
     8  2.921800075398, 2.999100000908, 3.075700050695, 3.195899898596, MgO
     9  3.321200125893, 3.444999920496, 3.596000004674, 3.723799892261, MgO
     A  3.801000163647, 3.870899913321, 3.948099868748, 3.980100038179, MgO
     B  3.990399785742, 4.000000000000,      4*0.0D+00/                 MgO
      DATA  K_MgO/                                                      071215
     1  6.62169554D+00, 6.57396601D+00, 6.53178774D+00, 6.39575614D+00, MgO
     2  6.10411021D+00, 5.85705722D+00, 5.64078524D+00, 5.48547037D+00, MgO
     3  5.44746842D+00, 5.47809668D+00, 5.54153917D+00, 5.61543383D+00, MgO
     4  5.70364142D+00, 5.99122950D+00, 6.22535556D+00, 6.49540361D+00, MgO
     5  6.96635549D+00, 7.48517311D+00, 7.68780164D+00, 7.87688146D+00, MgO
     6  8.33872340D+00, 8.52341782D+00, 8.73213257D+00, 8.97444246D+00, MgO
     7  9.23426536D+00, 9.46921120D+00, 9.64709087D+00, 9.83318972D+00, MgO
     8  9.89623129D+00, 9.93036659D+00, 9.93989313D+00, 9.91917519D+00, MgO
     9  9.87742997D+00, 9.83815471D+00, 9.80782281D+00, 9.80588241D+00, MgO
     A  9.82229742D+00, 9.85657710D+00, 9.93193216D+00, 9.98130027D+00, MgO
     B  1.00001034D+01, 1.00190024D+01,      4*0.0D+00/                 MgO
      DATA TK_AlH/                                                      071215
     1 -0.999999993529,-0.991700000218,-0.984500008606,-0.959399978211, AlH
     2 -0.895899970115,-0.822599980062,-0.735300017963,-0.621500006004, AlH
     3 -0.519699977483,-0.393299985465,-0.279599957217,-0.119099996326, AlH
     4  0.032799975372, 0.369899987302, 0.484700015592, 0.599500026459, AlH
     5  0.706000044739, 0.858600022416, 1.032399981097, 1.207300019902, AlH
     6  1.407200081186, 1.621799929870, 1.731500020998, 1.834899910997, AlH
     7  2.001000000297, 2.136700116348, 2.344099919281, 2.569199949518, AlH
     8  2.715000042458, 2.887100109091, 3.058399971986, 3.254200093843, AlH
     9  3.339700156729, 3.419800087662, 3.614800017063, 3.688500036968, AlH
     A  3.763700140411, 3.834300020388, 3.894699975440, 3.959200029267, AlH
     B  3.984399924657, 4.000000000000,      4*0.0D+00/                 AlH
      DATA  K_AlH/                                                      071215
     1  4.55751209D+00, 4.52036030D+00, 4.48895424D+00, 4.38590313D+00, AlH
     2  4.16554540D+00, 3.97483880D+00, 3.82270926D+00, 3.72373520D+00, AlH
     3  3.71075384D+00, 3.77105384D+00, 3.88096185D+00, 4.10222297D+00, AlH
     4  4.36315790D+00, 5.05215449D+00, 5.30596695D+00, 5.55965043D+00, AlH
     5  5.78723148D+00, 6.09055104D+00, 6.40524625D+00, 6.70132011D+00, AlH
     6  7.02680925D+00, 7.38091385D+00, 7.57433711D+00, 7.76708815D+00, AlH
     7  8.09176227D+00, 8.36015208D+00, 8.75810324D+00, 9.16361350D+00, AlH
     8  9.40929584D+00, 9.67614245D+00, 9.90912661D+00, 1.01300050D+01, AlH
     9  1.02106090D+01, 1.02761785D+01, 1.03778541D+01, 1.03869137D+01, AlH
     A  1.03789714D+01, 1.03609935D+01, 1.03453552D+01, 1.03404209D+01, AlH
     B  1.03450694D+01, 1.03504927D+01,      4*0.0D+00/                 AlH
      DATA TK_AlO/                                                      071215
     1 -0.999999993529,-0.990200001427,-0.979600014579,-0.947099995698, AlO
     2 -0.864099993040,-0.765599982735,-0.654599970871,-0.481900024840, AlO
     3 -0.334299966773,-0.261699993179,-0.190699978943,-0.051400013623, AlO
     4  0.203800036281, 0.448400046186, 0.677199940810, 0.926499939871, AlO
     5  1.125300008707, 1.301999962176, 1.443300043053, 1.591999908461, AlO
     6  1.744100007486, 1.907900077036, 2.115999882221, 2.240000018636, AlO
     7  2.370699904645, 2.544299957408, 2.697200045389, 2.984100011843, AlO
     8  3.113000067159, 3.255100113486, 3.398200054199, 3.552300137074, AlO
     9  3.708399998959, 3.767600052115, 3.828399965078, 3.883900019978, AlO
     A  3.934599924312, 3.974099881597, 3.988899805855, 4.000000000000, AlO
     B       6*0.0D+00/                                                 AlO
      DATA  K_AlO/                                                      071215
     1  8.39005120D+00, 8.29687473D+00, 8.19903917D+00, 7.91745446D+00, AlO
     2  7.31266503D+00, 6.77620154D+00, 6.36001828D+00, 6.00595400D+00, AlO
     3  5.89493194D+00, 5.88183637D+00, 5.88745826D+00, 5.93820815D+00, AlO
     4  6.12949751D+00, 6.39148922D+00, 6.67753628D+00, 7.01538582D+00, AlO
     5  7.29645936D+00, 7.55190855D+00, 7.76077459D+00, 7.99162746D+00, AlO
     6  8.25176274D+00, 8.56709346D+00, 9.00479214D+00, 9.27034474D+00, AlO
     7  9.54391357D+00, 9.88573469D+00, 1.01568318D+01, 1.05726453D+01, AlO
     8  1.07185512D+01, 1.08448572D+01, 1.09261724D+01, 1.09617365D+01, AlO
     9  1.09622935D+01, 1.09589742D+01, 1.09567730D+01, 1.09592333D+01, AlO
     A  1.09697955D+01, 1.09873387D+01, 1.09967170D+01, 1.10049267D+01, AlO
     B       6*0.0D+00/                                                 AlO
      DATA TK_AlF/                                                      071215
     1 -0.999999993529,-0.990000001588,-0.978700015097,-0.944999997364, AlF
     2 -0.909099971604,-0.858500006382,-0.757299999638,-0.642099985919, AlF
     3 -0.555999996355,-0.472100041820,-0.393499984912,-0.316899998606, AlF
     4 -0.240700004785,-0.163600016741, 0.004700005198, 0.199000042905, AlF
     5  0.407999955146, 0.653599955589, 0.907199970470, 1.115499882340, AlF
     6  1.307099959510, 1.453200038316, 1.671400069432, 1.821899936280, AlF
     7  1.990200009885, 2.139800128070, 2.345399918454, 2.632700101787, AlF
     8  2.823899942001, 3.039999981491, 3.221299811021, 3.394699980943, AlF
     9  3.530200194319, 3.651900050636, 3.729300016221, 3.801800144398, AlF
     A  3.864500024888, 3.926799913398, 3.970599787143, 3.987999829615, AlF
     B  4.000000000000,      5*0.0D+00/                                 AlF
      DATA  K_AlF/                                                      071215
     1  1.02779588D+01, 1.01461233D+01, 1.00014490D+01, 9.59588451D+00, AlF
     2  9.20401387D+00, 8.71584902D+00, 7.93476519D+00, 7.30362422D+00, AlF
     3  6.97226603D+00, 6.73917755D+00, 6.58397120D+00, 6.47857483D+00, AlF
     4  6.40914538D+00, 6.36796425D+00, 6.35776597D+00, 6.44668705D+00, AlF
     5  6.62231279D+00, 6.89407910D+00, 7.21816863D+00, 7.50354198D+00, AlF
     6  7.77553394D+00, 7.98923376D+00, 8.33265806D+00, 8.59817898D+00, AlF
     7  8.91913710D+00, 9.21294225D+00, 9.61301613D+00, 1.01319773D+01, AlF
     8  1.04261482D+01, 1.06988264D+01, 1.08830093D+01, 1.10276394D+01, AlF
     9  1.11222859D+01, 1.11939649D+01, 1.12320910D+01, 1.12617955D+01, AlF
     A  1.12835341D+01, 1.13050394D+01, 1.13239849D+01, 1.13332739D+01, AlF
     B  1.13404593D+01,      5*0.0D+00/                                 AlF
      DATA TK_Al2/                                                      071215
     1 -0.999999993529,-0.989000002864,-0.974500017513,-0.934600001292, Al2
     2 -0.836000018446,-0.760799995860,-0.686600008373,-0.535699951601, Al2
     3 -0.281999967430,-0.043900020307, 0.163499967603, 0.376299981343, Al2
     4  0.590400024612, 0.799099948868, 1.052199954200, 1.277099967309, Al2
     5  1.435600051026, 1.567499947388, 1.691200068517, 1.776799970983, Al2
     6  1.859299907748, 2.071099936707, 2.167400064884, 2.275999988777, Al2
     7  2.398600106420, 2.521599978276, 2.774999984206, 2.966300024728, Al2
     8  3.149600081729, 3.239899777149, 3.355000029370, 3.443299963641, Al2
     9  3.533800107272, 3.632900040941, 3.727399973399, 3.819000144482, Al2
     A  3.875600019645, 3.929299844368, 3.972099827623, 3.988399819055, Al2
     B  4.000000000000,      5*0.0D+00/                                 Al2
      DATA  K_Al2/                                                      071215
     1  4.58061654D+00, 4.56880212D+00, 4.55460893D+00, 4.52301640D+00, Al2
     2  4.48462144D+00, 4.48453254D+00, 4.50243399D+00, 4.57898989D+00, Al2
     3  4.79446729D+00, 5.06218544D+00, 5.32647887D+00, 5.61596458D+00, Al2
     4  5.91878774D+00, 6.22084522D+00, 6.59260771D+00, 6.92624888D+00, Al2
     5  7.16648742D+00, 7.38036050D+00, 7.60613263D+00, 7.78114347D+00, Al2
     6  7.96349597D+00, 8.46236535D+00, 8.68498183D+00, 8.92065575D+00, Al2
     7  9.16070099D+00, 9.37140440D+00, 9.71755559D+00, 9.91527031D+00, Al2
     8  1.00655020D+01, 1.01263336D+01, 1.01878581D+01, 1.02190045D+01, Al2
     9  1.02345456D+01, 1.02358830D+01, 1.02306947D+01, 1.02330445D+01, Al2
     A  1.02482570D+01, 1.02836088D+01, 1.03354995D+01, 1.03626055D+01, Al2
     B  1.03847824D+01,      5*0.0D+00/                                 Al2
      DATA TK_SiH/                                                      071215
     1 -0.999999993529,-0.991400000460,-0.983400010010,-0.956299983219, SiH
     2 -0.887899978293,-0.808300014501,-0.715200032976,-0.594200036002, SiH
     3 -0.476700033375,-0.346800011107,-0.216100000924,-0.052800012745, SiH
     4  0.097200012936, 0.400199986531, 0.551899959309, 0.729800036262, SiH
     5  0.811599968243, 0.892199977976, 1.084999925294, 1.246500010038, SiH
     6  1.419700080198, 1.493799995076, 1.571699944776, 1.710300049358, SiH
     7  1.883500127143, 2.084099929737, 2.249500002998, 2.441500062366, SiH
     8  2.661800090891, 2.963000033211, 3.067300121039, 3.173299995177, SiH
     9  3.452099844378, 3.578299967308, 3.701799848691, 3.764000133619, SiH
     A  3.829799931523, 3.896800024868, 3.945299928240, 3.979900038121, SiH
     B  4.000000000000,      5*0.0D+00/                                 SiH
      DATA  K_SiH/                                                      071215
     1  3.65534492D+00, 3.61686498D+00, 3.58207152D+00, 3.47166348D+00, SiH
     2  3.23942936D+00, 3.04254116D+00, 2.89555418D+00, 2.81182856D+00, SiH
     3  2.82063187D+00, 2.90828339D+00, 3.05830139D+00, 3.30881074D+00, SiH
     4  3.58288551D+00, 4.21723399D+00, 4.55753397D+00, 4.95261462D+00, SiH
     5  5.12496593D+00, 5.28652685D+00, 5.64306954D+00, 5.92197992D+00, SiH
     6  6.22228646D+00, 6.35797636D+00, 6.50808667D+00, 6.79483931D+00, SiH
     7  7.18062352D+00, 7.64517546D+00, 8.02332670D+00, 8.43867407D+00, SiH
     8  8.87539059D+00, 9.39626873D+00, 9.55143305D+00, 9.69423694D+00, SiH
     9  1.00043352D+01, 1.01181956D+01, 1.02121876D+01, 1.02511425D+01, SiH
     A  1.02854951D+01, 1.03145178D+01, 1.03343915D+01, 1.03502081D+01, SiH
     B  1.03609207D+01,      5*0.0D+00/                                 SiH
      DATA TK_SiHm/                                                     071215
     1 -0.999999993529,-0.987600004651,-0.969600018389,-0.922500022255, SiHm
     2 -0.802500009240,-0.641899986117,-0.469100038442,-0.331599970128, SiHm
     3 -0.185599982521,-0.027199971384, 0.132400000632, 0.303499947192, SiHm
     4  0.456900041643, 0.588100029922, 0.705600044253, 0.788399957867, SiHm
     5  0.869600012910, 1.130800121324, 1.266899984152, 1.417200079117, SiHm
     6  1.597299897101, 1.784999972412, 2.148700096238, 2.391300099732, SiHm
     7  2.604699906216, 2.795399964225, 2.966900023186, 3.101099942664, SiHm
     8  3.221899826205, 3.398000050013, 3.520499987571, 3.607600071677, SiHm
     9  3.691000050106, 3.837700094444, 3.937599989054, 3.975499919378, SiHm
     A  4.000000000000,      9*0.0D+00/                                 SiHm
      DATA  K_SiHm/                                                     071215
     1 -5.58200441D-03, 6.28020300D-02, 1.60185182D-01, 4.05179578D-01, SiHm
     2  9.73353990D-01, 1.63533302D+00, 2.25843254D+00, 2.70830842D+00, SiHm
     3  3.15465724D+00, 3.61287514D+00, 4.05534888D+00, 4.51494610D+00, SiHm
     4  4.91742509D+00, 5.25328089D+00, 5.54189748D+00, 5.73437705D+00, SiHm
     5  5.91251383D+00, 6.42204177D+00, 6.66200362D+00, 6.91531022D+00, SiHm
     6  7.20801449D+00, 7.50482589D+00, 8.06627717D+00, 8.43534630D+00, SiHm
     7  8.75794332D+00, 9.04277491D+00, 9.28874548D+00, 9.46798681D+00, SiHm
     8  9.61774552D+00, 9.82587088D+00, 9.97666405D+00, 1.00914762D+01, SiHm
     9  1.02068120D+01, 1.04141262D+01, 1.05511335D+01, 1.06011581D+01, SiHm
     A  1.06328101D+01,      9*0.0D+00/                                 SiHm
      DATA TK_SiC/                                                      071215
     1 -0.999999993529,-0.990000001588,-0.978700015097,-0.944899997444, SiC
     2 -0.859400002218,-0.747400009795,-0.622400006173,-0.529199955048, SiC
     3 -0.440699996008,-0.362300021673,-0.287899999539,-0.184799983266, SiC
     4 -0.084999977745, 0.153299974273, 0.329500050777, 0.501999999757, SiC
     5  0.622100016110, 0.761799992301, 0.871200011722, 1.010099970625, SiC
     6  1.180900053864, 1.313699945686, 1.469100032236, 1.635600099226, SiC
     7  1.787399972540, 1.985200008184, 2.184000076667, 2.376300038281, SiC
     8  2.551999953579, 2.719600029942, 2.875100001206, 3.035300079752, SiC
     9  3.266200069228, 3.422800016643, 3.716399888974, 3.829199945904, SiC
     A  3.912099898179, 3.967499839897, 3.986999856016, 4.000000000000, SiC
     B       6*0.0D+00/                                                 SiC
      DATA  K_SiC/                                                      071215
     1  6.16113667D+00, 6.08056745D+00, 5.99241856D+00, 5.74627379D+00, SiC
     2  5.22998242D+00, 4.74721946D+00, 4.41053521D+00, 4.26537227D+00, SiC
     3  4.19051823D+00, 4.16280177D+00, 4.16135621D+00, 4.18841673D+00, SiC
     4  4.23864343D+00, 4.42872601D+00, 4.61523567D+00, 4.82343630D+00, SiC
     5  4.98195167D+00, 5.18626773D+00, 5.36772945D+00, 5.62934325D+00, SiC
     6  5.99244744D+00, 6.29809595D+00, 6.67740401D+00, 7.10895449D+00, SiC
     7  7.51699820D+00, 8.05000442D+00, 8.56645446D+00, 9.02817924D+00, SiC
     8  9.40621713D+00, 9.72103325D+00, 9.97070616D+00, 1.01874672D+01, SiC
     9  1.04374242D+01, 1.05714207D+01, 1.07760076D+01, 1.08516688D+01, SiC
     A  1.09111692D+01, 1.09563485D+01, 1.09742900D+01, 1.09871213D+01, SiC
     B       6*0.0D+00/                                                 SiC
      DATA TK_SiN/                                                      071215
     1 -0.999999993529,-0.989400002354,-0.976200016535,-0.939000001324, SiN
     2 -0.844400047679,-0.728200009424,-0.596300027552,-0.504499997742, SiN
     3 -0.415400017723,-0.277599960573,-0.156200005564, 0.005000005530, SiN
     4  0.133799997164, 0.346400040229, 0.560699940977, 0.812499978836, SiN
     5  1.100499908984, 1.270699980454, 1.396100093451, 1.512999969852, SiN
     6  1.646400095056, 1.799499943479, 1.937600043955, 2.089199924524, SiN
     7  2.204800047573, 2.326999928814, 2.625600002540, 2.842999898848, SiN
     8  2.995700004337, 3.207899975247, 3.349599925475, 3.511700181256, SiN
     9  3.608900100820, 3.707299973915, 3.875000006072, 3.954799933190, SiN
     A  3.982699969538, 4.000000000000,      8*0.0D+00/                 SiN
      DATA  K_SiN/                                                      071215
     1 -6.19781547D-03, 8.48609127D-02, 1.96064809D-01, 4.97123293D-01, SiN
     2  1.18897572D+00, 1.91869670D+00, 2.62315778D+00, 3.05273275D+00, SiN
     3  3.43021744D+00, 3.94578983D+00, 4.33637874D+00, 4.77648307D+00, SiN
     4  5.08009716D+00, 5.52076538D+00, 5.91658657D+00, 6.34532803D+00, SiN
     5  6.80892190D+00, 7.07821386D+00, 7.28343107D+00, 7.48924036D+00, SiN
     6  7.74785491D+00, 8.07320653D+00, 8.38417240D+00, 8.73341961D+00, SiN
     7  8.99749485D+00, 9.26833724D+00, 9.87061374D+00, 1.02371353D+01, SiN
     8  1.04527843D+01, 1.06991881D+01, 1.08327558D+01, 1.09522170D+01, SiN
     9  1.10044987D+01, 1.10443394D+01, 1.10966593D+01, 1.11232398D+01, SiN
     A  1.11347908D+01, 1.11430698D+01,      8*0.0D+00/                 SiN
      DATA TK_SiO/                                                      071215
     1 -0.999999993529,-0.989800001843,-0.977800015614,-0.965400000387, SiO
     2 -0.942999998951,-0.906099971214,-0.853600029053,-0.746500011822, SiO
     3 -0.627200007073,-0.535999951520,-0.442399995015,-0.359500024576, SiO
     4 -0.277799960237,-0.193499981023,-0.114300007852, 0.043099960934, SiO
     5  0.252900000885, 0.463200031365, 0.658599950987, 0.881700000619, SiO
     6  1.033799976159, 1.186600047764, 1.313899945021, 1.446900043566, SiO
     7  1.530799987021, 1.611899894969, 1.836099910505, 1.933100052635, SiO
     8  2.024399983199, 2.234800017046, 2.438100070396, 2.672300076658, SiO
     9  2.879100086103, 3.071900144142, 3.357300081738, 3.506800148605, SiO
     A  3.609900123238, 3.705599935209, 3.783100131056, 3.855600026090, SiO
     B  3.941900000479, 3.977199965256, 3.989699784735, 4.000000000000, SiO
     C       2*0.0D+00/                                                 SiO
      DATA  K_SiO/                                                      071215
     1  1.14050701D+01, 1.12383717D+01, 1.10479872D+01, 1.08575270D+01, SiO
     2  1.05291858D+01, 1.00299663D+01, 9.40176723D+00, 8.37779721D+00, SiO
     3  7.56648908D+00, 7.12897954D+00, 6.80845737D+00, 6.61035669D+00, SiO
     4  6.47626027D+00, 6.38672504D+00, 6.33748232D+00, 6.31504783D+00, SiO
     5  6.39896901D+00, 6.57050464D+00, 6.77957844D+00, 7.05469790D+00, SiO
     6  7.25666589D+00, 7.46821253D+00, 7.65295611D+00, 7.86199318D+00, SiO
     7  8.00689889D+00, 8.15836570D+00, 8.63247719D+00, 8.85848839D+00, SiO
     8  9.07929672D+00, 9.59868913D+00, 1.00815564D+01, 1.05801355D+01, SiO
     9  1.09509910D+01, 1.12332752D+01, 1.15577168D+01, 1.16992033D+01, SiO
     A  1.17906326D+01, 1.18720444D+01, 1.19340131D+01, 1.19839431D+01, SiO
     B  1.20202077D+01, 1.20239050D+01, 1.20234985D+01, 1.20225173D+01, SiO
     C       2*0.0D+00/                                                 SiO
      DATA TK_SiF/                                                      071215
     1 -0.999999993529,-0.989900001716,-0.978200015384,-0.943899998237, SiF
     2 -0.856000017949,-0.752800002289,-0.636199995446,-0.542299960170, SiF
     3 -0.450699988809,-0.357200027173,-0.266699981218,-0.122099995694, SiF
     4  0.037699967455, 0.239100003728, 0.444100052039, 0.647799962661, SiF
     5  0.903599979224, 1.055199953527, 1.198200037048, 1.306299959929, SiF
     6  1.473800025819, 1.614599893337, 1.776299970802, 2.018099992746, SiF
     7  2.355999909817, 2.544399957409, 2.719100031303, 2.913800085834, SiF
     8  3.082199987679, 3.307699998646, 3.465699892510, 3.641299875428, SiF
     9  3.785100086708, 3.872099940468, 3.923100015562, 3.962899966758, SiF
     A  3.985999882416, 4.000000000000,      8*0.0D+00/                 SiF
      DATA  K_SiF/                                                      071215
     1  8.07962771D+00, 7.97650970D+00, 7.86075347D+00, 7.54319774D+00, SiF
     2  6.86375720D+00, 6.27390768D+00, 5.82038016D+00, 5.58237072D+00, SiF
     3  5.43345627D+00, 5.34506424D+00, 5.30430969D+00, 5.30381091D+00, SiF
     4  5.37039194D+00, 5.52586675D+00, 5.73917199D+00, 5.98604852D+00, SiF
     5  6.32601237D+00, 6.53739754D+00, 6.74216216D+00, 6.90244160D+00, SiF
     6  7.17169276D+00, 7.42937191D+00, 7.76136028D+00, 8.30305963D+00, SiF
     7  9.07397324D+00, 9.47369452D+00, 9.80355498D+00, 1.01133486D+01, SiF
     8  1.03327834D+01, 1.05695037D+01, 1.07089429D+01, 1.08476240D+01, SiF
     9  1.09464754D+01, 1.09958770D+01, 1.10206660D+01, 1.10387790D+01, SiF
     A  1.10494543D+01, 1.10562432D+01,      8*0.0D+00/                 SiF
      DATA TK_Si2/                                                      071215
     1 -0.999999993529,-0.989800001843,-0.977800015614,-0.942899999030, Si2
     2 -0.854200026277,-0.742700020382,-0.626400006923,-0.529999953135, Si2
     3 -0.432800013553,-0.244000001101,-0.060600007210, 0.134899994439, Si2
     4  0.319400036096, 0.520299972591, 0.758299993201, 0.999800000032, Si2
     5  1.115299882387, 1.232700008682, 1.395600093315, 1.462100027608, Si2
     6  1.537899978877, 1.677200076479, 1.861799903518, 1.960100055337, Si2
     7  2.060199938002, 2.256599987593, 2.382000120969, 2.508699999392, Si2
     8  2.604799906179, 2.694300045058, 2.897700107740, 3.063900033110, Si2
     9  3.309499954132, 3.421600049335, 3.545100079340, 3.647500032975, Si2
     A  3.738599802159, 3.850099892669, 3.897000029576, 3.940400032349, Si2
     B  3.976899957160, 3.989599787375, 4.000000000000,      3*0.0D+00/ Si2
      DATA  K_Si2/                                                      071215
     1  5.67283853D+00, 5.62347035D+00, 5.56757223D+00, 5.41768102D+00, Si2
     2  5.11227801D+00, 4.85177473D+00, 4.68672191D+00, 4.61013792D+00, Si2
     3  4.57573256D+00, 4.60303015D+00, 4.71514437D+00, 4.89560583D+00, Si2
     4  5.10357150D+00, 5.35674397D+00, 5.67915871D+00, 6.02121930D+00, Si2
     5  6.18859999D+00, 6.36316815D+00, 6.62817059D+00, 6.75174384D+00, Si2
     6  6.90811006D+00, 7.24154241D+00, 7.75829011D+00, 8.05515111D+00, Si2
     7  8.36331898D+00, 8.94930255D+00, 9.28242481D+00, 9.56780251D+00, Si2
     8  9.74633137D+00, 9.88490355D+00, 1.01217942D+01, 1.02624627D+01, Si2
     9  1.04295120D+01, 1.05019500D+01, 1.05848090D+01, 1.06558949D+01, Si2
     A  1.07182848D+01, 1.07883335D+01, 1.08156423D+01, 1.08417900D+01, Si2
     B  1.08671618D+01, 1.08773945D+01, 1.08865496D+01,      3*0.0D+00/ Si2
      DATA TK_HS/                                                       071215
     1 -0.999999993529,-0.991500000379,-0.983700009627,-0.957099981927, HS
     2 -0.889899970082,-0.811800010059,-0.720300014019,-0.597300023529, HS
     3 -0.479300028602,-0.355900028641,-0.220400013890,-0.063300002629, HS
     4  0.103600013071, 0.282399972463, 0.445500050134, 0.551799959643, HS
     5  0.654099955129, 0.841300028049, 0.924699938865, 1.008599974703, HS
     6  1.134600119466, 1.275999969569, 1.428100065299, 1.571899944100, HS
     7  1.925500066287, 2.204800047573, 2.309199943483, 2.422500091292, HS
     8  2.589699921391, 2.740600013863, 3.010500211824, 3.153300005517, HS
     9  3.297100126659, 3.564500043156, 3.667900016363, 3.765000110979, HS
     A  3.833600005141, 3.898500064882, 3.959500035818, 4.000000000000, HS
     B       6*0.0D+00/                                                 HS
      DATA  K_HS/                                                       071215
     1  4.84807780D+00, 4.80052841D+00, 4.75800072D+00, 4.62131936D+00, HS
     2  4.32813235D+00, 4.06979597D+00, 3.86117673D+00, 3.70870647D+00, HS
     3  3.66944125D+00, 3.71312218D+00, 3.83588720D+00, 4.04959080D+00, HS
     4  4.33575683D+00, 4.68834859D+00, 5.03798675D+00, 5.27548504D+00, HS
     5  5.50763645D+00, 5.92498456D+00, 6.10034896D+00, 6.26769895D+00, HS
     6  6.50337165D+00, 6.75151007D+00, 7.00614263D+00, 7.23921640D+00, HS
     7  7.79402208D+00, 8.22932001D+00, 8.39711452D+00, 8.58437618D+00, HS
     8  8.86924334D+00, 9.12994506D+00, 9.57983834D+00, 9.79459177D+00, HS
     9  9.98654520D+00, 1.02762147D+01, 1.03667382D+01, 1.04397952D+01, HS
     A  1.04829323D+01, 1.05165447D+01, 1.05422431D+01, 1.05570489D+01, HS
     B       6*0.0D+00/                                                 HS
      DATA TK_HSm/                                                      071215
     1 -0.999999993529,-0.991400000460,-0.983600009755,-0.956799982411, HSm
     2 -0.889199972956,-0.810600014048,-0.718600019672,-0.594600034393, HSm
     3 -0.475600035395,-0.351300033835,-0.215299998219,-0.059000008855, HSm
     4  0.107400017352, 0.287399972966, 0.451900043357, 0.559599933622, HSm
     5  0.661999952053, 0.849900024145, 0.935299993502, 1.019899974488, HSm
     6  1.147500102131, 1.289899946823, 1.448200043751, 1.604399893427, HSm
     7  1.955100059318, 2.111799883216, 2.266299986827, 2.486700000829, HSm
     8  2.634500100715, 2.763799994218, 3.061299965869, 3.193199964705, HSm
     9  3.330599929735, 3.542200009977, 3.653300016839, 3.751500022732, HSm
     A  3.822900096904, 3.888699897882, 3.956699974678, 4.000000000000, HSm
     B       6*0.0D+00/                                                 HSm
      DATA  K_HSm/                                                      071215
     1  5.59394280D+00, 5.54112159D+00, 5.49440704D+00, 5.34292095D+00, HSm
     2  5.01710444D+00, 4.72722685D+00, 4.48921346D+00, 4.30654529D+00, HSm
     3  4.24672026D+00, 4.27533258D+00, 4.38641947D+00, 4.58901741D+00, HSm
     4  4.86676743D+00, 5.21629633D+00, 5.56576492D+00, 5.80489687D+00, HSm
     5  6.03597791D+00, 6.45201180D+00, 6.62994334D+00, 6.79700316D+00, HSm
     6  7.03348260D+00, 7.28160966D+00, 7.54510901D+00, 7.79687300D+00, HSm
     7  8.34486858D+00, 8.58600102D+00, 8.82519747D+00, 9.17569562D+00, HSm
     8  9.41767897D+00, 9.63120977D+00, 1.01034159D+01, 1.02902654D+01, HSm
     9  1.04630909D+01, 1.06816453D+01, 1.07723055D+01, 1.08368855D+01, HSm
     A  1.08733596D+01, 1.08987310D+01, 1.09174885D+01, 1.09265232D+01, HSm
     B       6*0.0D+00/                                                 HSm
      DATA TK_CS/                                                       071215
     1 -0.999999993529,-0.989900001716,-0.978400015269,-0.944299997920, CS
     2 -0.908299971500,-0.857100012859,-0.752400002525,-0.635499996818, CS
     3 -0.545999975816,-0.457299972125,-0.374299999620,-0.293200015308, CS
     4 -0.213099990781,-0.132899996676,-0.056600010361, 0.020600002735, CS
     5  0.259699995129, 0.368999991088, 0.478400013961, 0.606700020904, CS
     6  0.754499990446, 0.883899993869, 1.018499973936, 1.160200098786, CS
     7  1.342199902770, 1.521599969871, 1.733000021103, 1.944300049651, CS
     8  2.190600059003, 2.302899956999, 2.426500088680, 2.569499949894, CS
     9  2.808399927264, 3.000600013389, 3.215099897774, 3.384500004618, CS
     A  3.569400157668, 3.649200076174, 3.728099989175, 3.794900080784, CS
     B  3.850599904798, 3.942299991980, 3.977399970654, 3.989699784735, CS
     C  4.000000000000,      1*0.0D+00/                                 CS
      DATA  K_CS/                                                       071215
     1  1.03882486D+01, 1.02442335D+01, 1.00850042D+01, 9.64119637D+00, CS
     2  9.21596075D+00, 8.68098114D+00, 7.80764789D+00, 7.11730150D+00, CS
     3  6.74835988D+00, 6.49070759D+00, 6.32766577D+00, 6.22574743D+00, CS
     4  6.16792604D+00, 6.14219594D+00, 6.14048203D+00, 6.15670543D+00, CS
     5  6.29377004D+00, 6.38955618D+00, 6.50080311D+00, 6.64878444D+00, CS
     6  6.84658182D+00, 7.05156508D+00, 7.29962408D+00, 7.59318442D+00, CS
     7  7.99953073D+00, 8.40517655D+00, 8.85994240D+00, 9.27843303D+00, CS
     8  9.73177753D+00, 9.93381768D+00, 1.01565925D+01, 1.04138251D+01, CS
     9  1.08224699D+01, 1.11089356D+01, 1.13729723D+01, 1.15447186D+01, CS
     A  1.17089196D+01, 1.17764055D+01, 1.18412789D+01, 1.18922291D+01, CS
     B  1.19281040D+01, 1.19602389D+01, 1.19592782D+01, 1.19568849D+01, CS
     C  1.19540367D+01,      1*0.0D+00/                                 CS
      DATA TK_NS/                                                       071215
     1 -0.999999993529,-0.990600001105,-0.981100012945,-0.950499992590, NS
     2 -0.873399992973,-0.783599984882,-0.675099980713,-0.541399956364, NS
     3 -0.395999977999,-0.266799980979,-0.140499984549,-0.039000018808, NS
     4  0.059499926860, 0.333700043951, 0.561399949667, 0.814300000023, NS
     5  1.093799922793, 1.389700092007, 1.613299894122, 1.818699941902, NS
     6  1.935900047234, 2.057599943794, 2.263099983400, 2.347599917054, NS
     7  2.441800061510, 2.589699921391, 2.779999988640, 2.923600075132, NS
     8  3.090100136942, 3.243799861068, 3.409899832414, 3.540999981275, NS
     9  3.664399944917, 3.784300104447, 3.867699948131, 3.933899909205, NS
     A  3.974499892392, 4.000000000000,      8*0.0D+00/                 NS
      DATA  K_NS/                                                       071215
     1  8.21172874D+00, 8.13215919D+00, 8.05395277D+00, 7.81687091D+00, NS
     2  7.31086816D+00, 6.86347558D+00, 6.48805694D+00, 6.21800693D+00, NS
     3  6.10028820D+00, 6.09874073D+00, 6.15348603D+00, 6.22198935D+00, NS
     4  6.30335582D+00, 6.58820780D+00, 6.86984942D+00, 7.21053674D+00, NS
     5  7.60646410D+00, 8.03751729D+00, 8.36744193D+00, 8.66938528D+00, NS
     6  8.83712435D+00, 9.00511742D+00, 9.27861575D+00, 9.39155910D+00, NS
     7  9.52016365D+00, 9.72749655D+00, 9.99269485D+00, 1.01795241D+01, NS
     8  1.03730986D+01, 1.05281687D+01, 1.06732894D+01, 1.07759973D+01, NS
     9  1.08683721D+01, 1.09597839D+01, 1.10264372D+01, 1.10802751D+01, NS
     A  1.11126405D+01, 1.11323182D+01,      8*0.0D+00/                 NS
      DATA TK_SO/                                                       071215
     1 -0.999999993529,-0.990300001346,-0.979800014464,-0.947399995460, SO
     2 -0.864899991791,-0.766199981094,-0.655399969615,-0.570000007473, SO
     3 -0.484300021712,-0.335399965406,-0.195199982286,-0.061600005513, SO
     4  0.218000030127, 0.485000015669, 0.763399990406, 1.050999954469, SO
     5  1.284599954602, 1.496899993623, 1.665900075064, 1.900800081348, SO
     6  2.121499912734, 2.238600018208, 2.353299912387, 2.502599986981, SO
     7  2.724000018265, 2.898000107625, 3.087000080901, 3.394899985129, SO
     8  3.498099945730, 3.597699960740, 3.764100131355, 3.896600020161, SO
     9  3.959400033634, 3.984499922017, 4.000000000000,     11*0.0D+00/ SO
      DATA  K_SO/                                                       071215
     1  8.75191752D+00, 8.65746274D+00, 8.55815815D+00, 8.27034148D+00, SO
     2  7.65279492D+00, 7.09903198D+00, 6.66947776D+00, 6.44855703D+00, SO
     3  6.30342703D+00, 6.18761430D+00, 6.18331753D+00, 6.23498166D+00, SO
     4  6.45092035D+00, 6.74525976D+00, 7.10320604D+00, 7.50232395D+00, SO
     5  7.83839264D+00, 8.14923187D+00, 8.40053062D+00, 8.76360178D+00, SO
     6  9.13092018D+00, 9.33770330D+00, 9.54727282D+00, 9.82593094D+00, SO
     7  1.02283955D+01, 1.05121181D+01, 1.07740092D+01, 1.10955407D+01, SO
     8  1.11757046D+01, 1.12418157D+01, 1.13343768D+01, 1.13965998D+01, SO
     9  1.14209509D+01, 1.14292388D+01, 1.14338732D+01,     11*0.0D+00/ SO
      DATA TK_MgS/                                                      071215
     1 -0.999999993529,-0.989100002737,-0.975100017168,-0.936600001306, MgS
     2 -0.836200019985,-0.721000013612,-0.585800048200,-0.401999979290, MgS
     3 -0.224100011860,-0.032599985461, 0.161199974196, 0.328500050029, MgS
     4  0.500300004440, 0.668199959353, 0.845900025961, 1.044399954470, MgS
     5  1.249400018961, 1.644100093184, 1.829999913006, 2.016499994764, MgS
     6  2.198800043054, 2.390100098633, 2.584799923439, 2.732700005942, MgS
     7  2.882100106354, 3.153100010720, 3.421800043886, 3.617199962854, MgS
     8  3.710600021662, 3.802700122743, 3.866399979313, 3.924799968622, MgS
     9  3.970099773649, 3.987799834895, 4.000000000000,     11*0.0D+00/ MgS
      DATA  K_MgS/                                                      071215
     1  2.38207215D-03, 1.11831201D-01, 2.49359166D-01, 6.10871361D-01, MgS
     2  1.44928421D+00, 2.25228945D+00, 3.01696412D+00, 3.82069149D+00, MgS
     3  4.42046677D+00, 4.94085197D+00, 5.38257012D+00, 5.71923998D+00, MgS
     4  6.03648506D+00, 6.32765841D+00, 6.62217369D+00, 6.94020589D+00, MgS
     5  7.26086510D+00, 7.86589125D+00, 8.14780056D+00, 8.43021608D+00, MgS
     6  8.70776098D+00, 8.99845323D+00, 9.28151214D+00, 9.47788178D+00, MgS
     7  9.65472371D+00, 9.91894848D+00, 1.01198437D+01, 1.02392794D+01, MgS
     8  1.02932070D+01, 1.03513334D+01, 1.04015088D+01, 1.04640832D+01, MgS
     9  1.05321321D+01, 1.05652050D+01, 1.05905249D+01,     11*0.0D+00/ MgS
      DATA TK_AlS/                                                      071215
     1 -0.999999993529,-0.990100001508,-0.979200014809,-0.945999996571, AlS
     2 -0.861499997100,-0.759399998401,-0.647399980670,-0.570100007815, AlS
     3 -0.493700010078,-0.338299961803,-0.171200002305, 0.011800009824, AlS
     4  0.224300026952, 0.444500051495, 0.671499957167, 0.916999944465, AlS
     5  1.129700114502, 1.316399936704, 1.454500036003, 1.616999891886, AlS
     6  1.732600021075, 1.853499907339, 2.068399937206, 2.239700018544, AlS
     7  2.392500100831, 2.559899938186, 2.737800012666, 2.900400106270, AlS
     8  3.079899947411, 3.396800024897, 3.506600143942, 3.619399913163, AlS
     9  3.759900221810, 3.816400087460, 3.870499904273, 3.951799867684, AlS
     A  3.981500001219, 3.990799794669, 4.000000000000,      7*0.0D+00/ AlS
      DATA  K_AlS/                                                      071215
     1  7.20083446D+00, 7.13896999D+00, 7.07309364D+00, 6.88626833D+00, AlS
     2  6.49533705D+00, 6.15690758D+00, 5.91266473D+00, 5.80175512D+00, AlS
     3  5.72785730D+00, 5.66299864D+00, 5.68984975D+00, 5.79940311D+00, AlS
     4  5.99598021D+00, 6.24915359D+00, 6.54205696D+00, 6.88034921D+00, AlS
     5  7.18412907D+00, 7.45608041D+00, 7.66158329D+00, 7.91606594D+00, AlS
     6  8.11213410D+00, 8.33261229D+00, 8.74985904D+00, 9.08714746D+00, AlS
     7  9.38032444D+00, 9.68215184D+00, 9.96834445D+00, 1.01927454D+01, AlS
     8  1.04000055D+01, 1.06808687D+01, 1.07594872D+01, 1.08329695D+01, AlS
     9  1.09137949D+01, 1.09432478D+01, 1.09717824D+01, 1.10246325D+01, AlS
     A  1.10509073D+01, 1.10603063D+01, 1.10702131D+01,      7*0.0D+00/ AlS
      DATA TK_SiS/                                                      071215
     1 -0.999999993529,-0.989700001971,-0.977600015730,-0.942499999348, SiS
     2 -0.905499971137,-0.852500034142,-0.744500016327,-0.624600006585, SiS
     3 -0.540899954249,-0.457499971620,-0.374999997305,-0.288000000083, SiS
     4 -0.108400017191, 0.078900040951, 0.263199987161, 0.452000043323, SiS
     5  0.655699953656, 0.867600014833, 1.046399954550, 1.200700033999, SiS
     6  1.488000000362, 1.601099891841, 1.718000043105, 1.976300019997, SiS
     7  2.297799959581, 2.485800001912, 2.639900097497, 2.814499931140, SiS
     8  3.007900176287, 3.213699931959, 3.428099872251, 3.543700045854, SiS
     9  3.658599888896, 3.751100013252, 3.845200012306, 3.933499900573, SiS
     A  3.972799846514, 4.000000000000,      8*0.0D+00/                 SiS
      DATA  K_SiS/                                                      071215
     1  9.80579039D+00, 9.68077804D+00, 9.53843949D+00, 9.15173633D+00, SiS
     2  8.78368101D+00, 8.32046691D+00, 7.57485208D+00, 6.99145705D+00, SiS
     3  6.69867611D+00, 6.47924928D+00, 6.32013093D+00, 6.20399599D+00, SiS
     4  6.09500397D+00, 6.11626788D+00, 6.22589475D+00, 6.39772807D+00, SiS
     5  6.62666267D+00, 6.89506334D+00, 7.13678047D+00, 7.35378261D+00, SiS
     6  7.79989814D+00, 8.00589336D+00, 8.23962355D+00, 8.80793870D+00, SiS
     7  9.54737011D+00, 9.96065276D+00, 1.02705073D+01, 1.05775273D+01, SiS
     8  1.08589087D+01, 1.10988374D+01, 1.13049531D+01, 1.14077250D+01, SiS
     9  1.15085212D+01, 1.15885547D+01, 1.16648426D+01, 1.17233054D+01, SiS
     A  1.17431908D+01, 1.17550057D+01,      8*0.0D+00/                 SiS
      DATA TK_S2/                                                       071215
     1 -0.999999993529,-0.990200001427,-0.979500014636,-0.946899995857, S2
     2 -0.863799993508,-0.765399983281,-0.654599970871,-0.580900042715, S2
     3 -0.507399994406,-0.362100021876,-0.188299980007,-0.014799984774, S2
     4  0.109400019606, 0.234500012800, 0.477600013754, 0.755899991461, S2
     5  1.030999986034, 1.295799956272, 1.573099940045, 1.772299969352, S2
     6  1.952100061683, 2.175100072530, 2.272899989817, 2.369799888302, S2
     7  2.719800029398, 2.905200099168, 3.124999995853, 3.325800015533, S2
     8  3.526000110067, 3.601099925960, 3.687300004911, 3.761200197011, S2
     9  3.824500058555, 3.929699833323, 3.972199830322, 4.000000000000, S2
     A      10*0.0D+00/                                                 S2
      DATA  K_S2/                                                       071215
     1  8.31098143D+00, 8.23796829D+00, 8.16071448D+00, 7.94058475D+00, S2
     2  7.47326551D+00, 7.06552279D+00, 6.75119241D+00, 6.60634612D+00, S2
     3  6.50151784D+00, 6.38507267D+00, 6.36484110D+00, 6.43582365D+00, S2
     4  6.52599055D+00, 6.64086315D+00, 6.91194447D+00, 7.27032109D+00, S2
     5  7.65208049D+00, 8.03331681D+00, 8.44027949D+00, 8.73539482D+00, S2
     6  9.00378250D+00, 9.34736645D+00, 9.50666603D+00, 9.67085440D+00, S2
     7  1.02745366D+01, 1.05576924D+01, 1.08348663D+01, 1.10296021D+01, S2
     8  1.11760329D+01, 1.12231020D+01, 1.12733667D+01, 1.13119655D+01, S2
     9  1.13392303D+01, 1.13647114D+01, 1.13662026D+01, 1.13645462D+01, S2
     A      10*0.0D+00/                                                 S2
      DATA TK_HCl/                                                      071215
     1 -0.999999993529,-0.991200000621,-0.982900010648,-0.955099985158, HCl
     2 -0.884999990199,-0.803100009784,-0.707900040314,-0.578500036573, HCl
     3 -0.451699986281,-0.315099999587,-0.177599991710, 0.004700005198, HCl
     4  0.170399948281, 0.337400036751, 0.518799973287, 0.670999958602, HCl
     5  0.866800015602, 0.958800051940, 1.048599954637, 1.300399963012, HCl
     6  1.459100027821, 1.617599891524, 1.892900108375, 2.137900120885, HCl
     7  2.342899920045, 2.572799943499, 2.696800045344, 2.819099940005, HCl
     8  3.129400105531, 3.257400163686, 3.390299888850, 3.593600066698, HCl
     9  3.776500129090, 3.907399910774, 3.963799941937, 4.000000000000, HCl
     A      10*0.0D+00/                                                 HCl
      DATA  K_HCl/                                                      071215
     1  6.23942081D+00, 6.17247833D+00, 6.11090789D+00, 5.91602542D+00, HCl
     2  5.49508682D+00, 5.11488366D+00, 4.79806641D+00, 4.53827772D+00, HCl
     3  4.42977270D+00, 4.43259753D+00, 4.52768578D+00, 4.75457871D+00, HCl
     4  5.02965800D+00, 5.35151837D+00, 5.73523761D+00, 6.07403150D+00, HCl
     5  6.50988446D+00, 6.70427697D+00, 6.88384448D+00, 7.34113249D+00, HCl
     6  7.60765403D+00, 7.86458983D+00, 8.29754742D+00, 8.67457167D+00, HCl
     7  8.98726357D+00, 9.34143930D+00, 9.53750694D+00, 9.73456602D+00, HCl
     8  1.02256817D+01, 1.04084354D+01, 1.05780386D+01, 1.07937881D+01, HCl
     9  1.09415492D+01, 1.10172457D+01, 1.10412346D+01, 1.10540260D+01, HCl
     A      10*0.0D+00/                                                 HCl
      DATA TK_LiCl/                                                     071215
     1 -0.999999993529,-0.990600001105,-0.980800013328,-0.949699993635, LiCl
     2 -0.870999986517,-0.780599971703,-0.670399962507,-0.525699963419, LiCl
     3 -0.372100006896,-0.241500003892,-0.117200000888,-0.002699999565, LiCl
     4  0.128300005305, 0.398699986658, 0.622100016110, 0.857900022553, LiCl
     5  1.123799972641, 1.412000076869, 1.878300089356, 2.010100002839, LiCl
     6  2.147100102231, 2.340599921508, 2.500999983726, 2.909300093102, LiCl
     7  3.158299875440, 3.442099994097, 3.552400134489, 3.656999927520, LiCl
     8  3.739799770088, 3.816500089654, 3.927899883025, 3.971499811431, LiCl
     9  3.988199824335, 4.000000000000,     12*0.0D+00/                 LiCl
      DATA  K_LiCl/                                                     071215
     1  7.81513444D+00, 7.73474886D+00, 7.65327913D+00, 7.41020225D+00, LiCl
     2  6.89065570D+00, 6.43912347D+00, 6.05835368D+00, 5.77500689D+00, LiCl
     3  5.66449901D+00, 5.66832253D+00, 5.72098798D+00, 5.79699030D+00, LiCl
     4  5.90854190D+00, 6.19986009D+00, 6.48111339D+00, 6.80110426D+00, LiCl
     5  7.17847389D+00, 7.59838191D+00, 8.28898802D+00, 8.48534616D+00, LiCl
     6  8.68921642D+00, 8.97268154D+00, 9.19707532D+00, 9.69675956D+00, LiCl
     7  9.94062175D+00, 1.01554426D+01, 1.02210084D+01, 1.02760765D+01, LiCl
     8  1.03188911D+01, 1.03666127D+01, 1.04861785D+01, 1.05623784D+01, LiCl
     9  1.05969325D+01, 1.06231742D+01,     12*0.0D+00/                 LiCl
      DATA TK_NaCl/                                                     071215
     1 -0.999999993529,-0.990200001427,-0.979400014694,-0.946399996253, NaCl
     2 -0.863199994445,-0.760899995586,-0.647199980868,-0.502899999583, NaCl
     3 -0.344699995228,-0.189599978796,-0.021099970231, 0.168599952985, NaCl
     4  0.363900012544, 0.548199963295, 0.743000011272, 0.971700013670, NaCl
     5  1.206000022679, 1.649500097579, 1.798299947168, 1.949400061897, NaCl
     6  2.199200042276, 2.308499944985, 2.415600099660, 2.718500032935, NaCl
     7  2.878500073368, 3.072900119551, 3.294400068554, 3.494399859377, NaCl
     8  3.571200142141, 3.649000071092, 3.725399928322, 3.791600011598, NaCl
     9  3.849899892786, 3.902200043876, 3.935299939418, 3.963299955727, NaCl
     A  3.985999882416, 4.000000000000,      8*0.0D+00/                 NaCl
      DATA  K_NaCl/                                                     071215
     1  7.76639986D+00, 7.69603373D+00, 7.62087016D+00, 7.40602358D+00, NaCl
     2  6.95316107D+00, 6.53956037D+00, 6.22045276D+00, 5.97163085D+00, NaCl
     3  5.84538139D+00, 5.82837728D+00, 5.89363300D+00, 6.03799010D+00, NaCl
     4  6.23815456D+00, 6.45806307D+00, 6.71158549D+00, 7.02664526D+00, NaCl
     5  7.36123835D+00, 8.01131108D+00, 8.23187202D+00, 8.45549777D+00, NaCl
     6  8.81277774D+00, 8.95763967D+00, 9.09029478D+00, 9.41984181D+00, NaCl
     7  9.57262840D+00, 9.73955856D+00, 9.90028973D+00, 1.00111699D+01, NaCl
     8  1.00435094D+01, 1.00722458D+01, 1.01038583D+01, 1.01505555D+01, NaCl
     9  1.02305052D+01, 1.03524660D+01, 1.04568499D+01, 1.05594305D+01, NaCl
     A  1.06500799D+01, 1.07085309D+01,      8*0.0D+00/                 NaCl
      DATA TK_AlCl/                                                     071215
     1 -0.999999993529,-0.990100001508,-0.978900014981,-0.945399997047, AlCl
     2 -0.860399998817,-0.755000000993,-0.641999986018,-0.562600001034, AlCl
     3 -0.484700021191,-0.333799967394,-0.171900001146, 0.000700000774, AlCl
     4  0.179199933070, 0.362200019695, 0.544299958198, 0.729900036289, AlCl
     5  1.000099999706, 1.255400003731, 1.356699894615, 1.446700043538, AlCl
     6  1.667300072560, 1.790899969914, 1.906300078008, 2.075799935240, AlCl
     7  2.229800015461, 2.418600095068, 2.610099903950, 2.835499918715, AlCl
     8  3.158599867636, 3.433599902175, 3.553400108634, 3.658799884068, AlCl
     9  3.731799983892, 3.799800183514, 3.859200113420, 3.928299871980, AlCl
     A  3.971099800636, 3.988099826975, 4.000000000000,      7*0.0D+00/ AlCl
      DATA  K_AlCl/                                                     071215
     1  8.71312957D+00, 8.62230274D+00, 8.52264663D+00, 8.24332074D+00, AlCl
     2  7.64797472D+00, 7.09597083D+00, 6.67810422D+00, 6.46664310D+00, AlCl
     3  6.31173121D+00, 6.12955558D+00, 6.06481747D+00, 6.09984831D+00, AlCl
     4  6.21275844D+00, 6.38242685D+00, 6.58702491D+00, 6.81963855D+00, AlCl
     5  7.18514452D+00, 7.54769191D+00, 7.69491077D+00, 7.82809632D+00, AlCl
     6  8.17593010D+00, 8.39219670D+00, 8.60779318D+00, 8.93564106D+00, AlCl
     7  9.22717591D+00, 9.55484070D+00, 9.84409586D+00, 1.01341462D+01, AlCl
     8  1.04694614D+01, 1.06864354D+01, 1.07619271D+01, 1.08172017D+01, AlCl
     9  1.08476309D+01, 1.08690831D+01, 1.08830285D+01, 1.08985547D+01, AlCl
     A  1.09128486D+01, 1.09206468D+01, 1.09270722D+01,      7*0.0D+00/ AlCl
      DATA TK_CaH/                                                      071215
     1 -0.999999993529,-0.987100005289,-0.967700010245,-0.917400014295, CaH
     2 -0.788900008166,-0.619900005844,-0.438499999987,-0.297400021006, CaH
     3 -0.152699988818, 0.095800014897, 0.192700048439, 0.288899973117, CaH
     4  0.502699997829, 0.598600026276, 0.692400063424, 0.874100009782, CaH
     5  1.069399937975, 1.251400016381, 1.429200063258, 1.786799972508, CaH
     6  2.241100016825, 2.374099985781, 2.508299998578, 2.736100010425, CaH
     7  2.940400046996, 3.172500014456, 3.414699953465, 3.518400015919, CaH
     8  3.626000033067, 3.720399815632, 3.790799994826, 3.856800055200, CaH
     9  3.912499908457, 3.965899884023, 3.986599866576, 4.000000000000, CaH
     A      10*0.0D+00/                                                 CaH
      DATA  K_CaH/                                                      071215
     1  2.60712338D+00, 2.58949002D+00, 2.56573374D+00, 2.51847815D+00, CaH
     2  2.47886902D+00, 2.56351830D+00, 2.77452979D+00, 2.99753589D+00, CaH
     3  3.26351645D+00, 3.77915858D+00, 3.99362147D+00, 4.21043094D+00, CaH
     4  4.68575878D+00, 4.88545546D+00, 5.06941513D+00, 5.39918899D+00, CaH
     5  5.72918713D+00, 6.02376886D+00, 6.30414168D+00, 6.85537461D+00, CaH
     6  7.54367451D+00, 7.74366372D+00, 7.94382227D+00, 8.27030789D+00, CaH
     7  8.53130032D+00, 8.77561735D+00, 8.96883984D+00, 9.03384301D+00, CaH
     8  9.09697752D+00, 9.16404707D+00, 9.23694779D+00, 9.33883711D+00, CaH
     9  9.46074157D+00, 9.61139899D+00, 9.67762096D+00, 9.72239296D+00, CaH
     A      10*0.0D+00/                                                 CaH
      DATA TK_CaF/                                                      071215
     1 -0.999999993529,-0.990100001508,-0.978900014981,-0.945499996967, CaF
     2 -0.859700000830,-0.760899995586,-0.645699982353,-0.568600006255, CaF
     3 -0.490100014169,-0.414600019349,-0.337499962797,-0.166600010905, CaF
     4  0.018100005497, 0.227700023805, 0.446400048908, 0.679199935070, CaF
     5  0.917399943368, 1.151600097582, 1.388100093152, 1.673600072105, CaF
     6  1.859699907776, 2.027499978161, 2.172300063882, 2.326399929144, CaF
     7  2.469000027118, 2.602999906846, 2.766099988859, 2.942700047349, CaF
     8  3.208999999898, 3.393899964199, 3.473199862064, 3.550900173271, CaF
     9  3.649700088879, 3.737099842247, 3.824600056158, 3.886899943668, CaF
     A  3.922600029368, 3.953799911355, 3.982599972178, 4.000000000000, CaF
     B       6*0.0D+00/                                                 CaF
      DATA  K_CaF/                                                      071215
     1  8.36662796D+00, 8.26771180D+00, 8.15915857D+00, 7.85567264D+00, CaF
     2  7.20172460D+00, 6.63692256D+00, 6.17516759D+00, 5.95621808D+00, CaF
     3  5.79023312D+00, 5.67434876D+00, 5.59254553D+00, 5.51537526D+00, CaF
     4  5.55109421D+00, 5.69058601D+00, 5.90686918D+00, 6.18514632D+00, CaF
     5  6.50019785D+00, 6.82737703D+00, 7.16793732D+00, 7.58670435D+00, CaF
     6  7.86232375D+00, 8.11249995D+00, 8.33024430D+00, 8.56249226D+00, CaF
     7  8.77268433D+00, 8.95964487D+00, 9.16642741D+00, 9.36068731D+00, CaF
     8  9.59840244D+00, 9.73081263D+00, 9.78193260D+00, 9.83217033D+00, CaF
     9  9.90609779D+00, 9.99642637D+00, 1.01308888D+01, 1.02671075D+01, CaF
     A  1.03634525D+01, 1.04582622D+01, 1.05532956D+01, 1.06134965D+01, CaF
     B       6*0.0D+00/                                                 CaF
      DATA TK_CaCl/                                                     071215
     1 -0.999999993529,-0.990000001588,-0.978800015039,-0.945199997205, CaCl
     2 -0.859700000830,-0.757699999402,-0.637399993094,-0.487900017020, CaCl
     3 -0.331299970501,-0.176999992703,-0.009399998486, 0.183899978146, CaCl
     4  0.383099979340, 0.575000054664, 0.775699984469, 0.991700001336, CaCl
     5  1.223400038066, 1.479800014746, 1.688300073030, 1.927600062706, CaCl
     6  2.171300060794, 2.362299901855, 2.677400065038, 2.818999939812, CaCl
     7  2.977500013601, 3.188100003282, 3.297200128811, 3.399500081408, CaCl
     8  3.476999953924, 3.553300111219, 3.635599965437, 3.714699927865, CaCl
     9  3.880700101376, 3.917300031786, 3.950699843664, 3.981500001219, CaCl
     A  4.000000000000,      9*0.0D+00/                                 CaCl
      DATA  K_CaCl/                                                     071215
     1  7.21223483D+00, 7.14220894D+00, 7.06615857D+00, 6.85233098D+00, CaCl
     2  6.39482158D+00, 5.98381212D+00, 5.64639971D+00, 5.39573561D+00, CaCl
     3  5.27968685D+00, 5.26904675D+00, 5.33858205D+00, 5.48972260D+00, CaCl
     4  5.69739770D+00, 5.92930710D+00, 6.19308853D+00, 6.49221904D+00, CaCl
     5  6.82379813D+00, 7.19826995D+00, 7.50617814D+00, 7.86090263D+00, CaCl
     6  8.21268552D+00, 8.46500155D+00, 8.82112613D+00, 8.96191986D+00, CaCl
     7  9.10767912D+00, 9.28057069D+00, 9.36018256D+00, 9.42890936D+00, CaCl
     8  9.47827082D+00, 9.52721246D+00, 9.58629392D+00, 9.65861530D+00, CaCl
     9  9.92218283D+00, 1.00129765D+01, 1.01076486D+01, 1.02037093D+01, CaCl
     A  1.02646621D+01,      9*0.0D+00/                                 CaCl
      DATA TK_ScO/                                                      071215
     1 -0.999999993529,-0.989900001716,-0.978100015442,-0.943599998475, ScO
     2 -0.906599971279,-0.854600024426,-0.751900002819,-0.633000001718, ScO
     3 -0.549599991039,-0.459499966564,-0.382199983808,-0.302300019435, ScO
     4 -0.223100012409,-0.142199983002, 0.044099960030, 0.263099987402, ScO
     5  0.498700006787, 0.787399961203, 1.081099921280, 1.284999954015, ScO
     6  1.469200032302, 1.626300029360, 1.806199942595, 1.962800052141, ScO
     7  2.127600040674, 2.312399939580, 2.501099983929, 2.687700047912, ScO
     8  2.871699929043, 3.097899963574, 3.321300123494, 3.473699874151, ScO
     9  3.623999988582, 3.714499932441, 3.812500001929, 3.939100021425, ScO
     A  3.975499919378, 4.000000000000,      8*0.0D+00/                 ScO
      DATA  K_ScO/                                                      071215
     1  1.04905877D+01, 1.03556394D+01, 1.02026517D+01, 9.78288748D+00, ScO
     2  9.37557516D+00, 8.87134261D+00, 8.07778176D+00, 7.43036196D+00, ScO
     3  7.11120118D+00, 6.86226514D+00, 6.71061207D+00, 6.60057006D+00, ScO
     4  6.52888940D+00, 6.48714633D+00, 6.48573632D+00, 6.60479123D+00, ScO
     5  6.82392971D+00, 7.16626913D+00, 7.55953101D+00, 7.84703741D+00, ScO
     6  8.11302918D+00, 8.34564315D+00, 8.62982914D+00, 8.90731727D+00, ScO
     7  9.23239355D+00, 9.61725313D+00, 1.00012408D+01, 1.03456254D+01, ScO
     8  1.06360282D+01, 1.09241016D+01, 1.11464953D+01, 1.12755805D+01, ScO
     9  1.14084106D+01, 1.15058888D+01, 1.16327748D+01, 1.18227186D+01, ScO
     A  1.18807639D+01, 1.19204085D+01,      8*0.0D+00/                 ScO
      DATA TK_TiO/                                                      071215
     1 -0.999999993529,-0.990000001588,-0.978600015154,-0.944799997523, TiO
     2 -0.908699971552,-0.857800009621,-0.756699999991,-0.640999987008, TiO
     3 -0.558499997865,-0.468200031208,-0.391499990442,-0.311600001494, TiO
     4 -0.234500007247,-0.156900008913, 0.022299998604, 0.253000000800, TiO
     5  0.493800012518, 0.755999991533, 1.031299984976, 1.259699990133, TiO
     6  1.475600022497, 1.610499895815, 1.743500009551, 1.874299990400, TiO
     7  1.991200008876, 2.110599883501, 2.229700015403, 2.368099891374, TiO
     8  2.550499956502, 2.718900031847, 2.900400106270, 3.100599930952, TiO
     9  3.261400186290, 3.430899840910, 3.595700012427, 3.700999830477, TiO
     A  3.812099993156, 3.918800070327, 4.000000000000,      7*0.0D+00/ TiO
      DATA  K_TiO/                                                      071215
     1  1.05003146D+01, 1.03689129D+01, 1.02234611D+01, 9.81821963D+00, TiO
     2  9.42583641D+00, 8.93725952D+00, 8.16130686D+00, 7.53149595D+00, TiO
     3  7.21421839D+00, 6.96421677D+00, 6.81440823D+00, 6.70566003D+00, TiO
     4  6.63693471D+00, 6.59691473D+00, 6.59207970D+00, 6.71493769D+00, TiO
     5  6.93994033D+00, 7.24973264D+00, 7.61532090D+00, 7.93580404D+00, TiO
     6  8.24353922D+00, 8.43251138D+00, 8.61568267D+00, 8.79856829D+00, TiO
     7  8.97160026D+00, 9.16314556D+00, 9.37036510D+00, 9.62640544D+00, TiO
     8  9.96908073D+00, 1.02664046D+01, 1.05470002D+01, 1.08007817D+01, TiO
     9  1.09666033D+01, 1.11223250D+01, 1.12808185D+01, 1.13960173D+01, TiO
     A  1.15324843D+01, 1.16758517D+01, 1.17895685D+01,      7*0.0D+00/ TiO
      DATA TK_TiS/                                                      071215
     1 -0.999999993529,-0.989900001716,-0.978300015327,-0.943999998158, TiS
     2 -0.857000013322,-0.749900004164,-0.630700006226,-0.472300041453, TiS
     3 -0.302600018770,-0.137799988621, 0.036899968747, 0.221200029822, TiS
     4  0.409699948305, 0.589400026233, 0.773999983908, 1.144800107421, TiS
     5  1.318399930051, 1.492999995451, 1.599399892599, 1.709900049562, TiS
     6  1.830799912678, 1.974700025997, 2.093699921187, 2.218300017417, TiS
     7  2.356699909151, 2.556599944616, 2.721500024883, 2.899600107015, TiS
     8  3.108300111313, 3.277700161675, 3.473699874151, 3.561899982394, TiS
     9  3.651500060292, 3.752000034582, 3.849299908044, 3.949299843252, TiS
     A  3.980400030259, 4.000000000000,      8*0.0D+00/                 TiS
      DATA  K_TiS/                                                      071215
     1  8.68678721D+00, 8.60217691D+00, 8.50802917D+00, 8.24760655D+00, TiS
     2  7.69460139D+00, 7.18689675D+00, 6.79209506D+00, 6.47098719D+00, TiS
     3  6.31238561D+00, 6.28589385D+00, 6.35123515D+00, 6.49019022D+00, TiS
     4  6.68130190D+00, 6.89369552D+00, 7.13196888D+00, 7.64623833D+00, TiS
     5  7.89557734D+00, 8.14487697D+00, 8.29239993D+00, 8.44084581D+00, TiS
     6  8.59987675D+00, 8.79239659D+00, 8.96206967D+00, 9.15480733D+00, TiS
     7  9.38507733D+00, 9.72588298D+00, 9.98945674D+00, 1.02393110D+01, TiS
     8  1.04821692D+01, 1.06489163D+01, 1.08378517D+01, 1.09336067D+01, TiS
     9  1.10438992D+01, 1.11865548D+01, 1.13431640D+01, 1.15182460D+01, TiS
     A  1.15747329D+01, 1.16107492D+01,      8*0.0D+00/                 TiS
      DATA TK_VO/                                                       071215
     1 -0.999999993529,-0.989700001971,-0.977200015960,-0.941600000062, VO
     2 -0.903499970877,-0.850000045709,-0.743600018354,-0.622200006135, VO
     3 -0.525799963180,-0.432100015219,-0.335999964661,-0.244700000320, VO
     4 -0.162400019076,-0.079299983434, 0.086000030100, 0.306199934827, VO
     5  0.514099975654, 0.754399990373, 0.998500000242, 1.259199991714, VO
     6  1.462500027873, 1.594099903960, 1.749599988559, 1.862499901854, VO
     7  1.977700014748, 2.099499917238, 2.219100013795, 2.336799923578, VO
     8  2.490599995948, 2.653400091349, 2.813799929791, 3.039999981491, VO
     9  3.248999979416, 3.420600076579, 3.568200129624, 3.711500001072, VO
     A  3.797000124811, 3.877200055840, 3.952899891703, 3.982199982738, VO
     B  4.000000000000,      5*0.0D+00/                                 VO
      DATA  K_VO/                                                       071215
     1  9.64981341D+00, 9.52535447D+00, 9.37906053D+00, 8.98947488D+00, VO
     2  8.61418762D+00, 8.15314537D+00, 7.43311562D+00, 6.86905943D+00, VO
     3  6.57074795D+00, 6.37609442D+00, 6.24926742D+00, 6.17950587D+00, VO
     4  6.14929176D+00, 6.14439251D+00, 6.19571050D+00, 6.35569957D+00, VO
     5  6.56900311D+00, 6.86109079D+00, 7.18741003D+00, 7.55461518D+00, VO
     6  7.84951589D+00, 8.04638394D+00, 8.29440977D+00, 8.49253491D+00, VO
     7  8.71472178D+00, 8.97188093D+00, 9.24307268D+00, 9.52150899D+00, VO
     8  9.88774181D+00, 1.02561595D+01, 1.05833087D+01, 1.09832795D+01, VO
     9  1.13009506D+01, 1.15276682D+01, 1.17021265D+01, 1.18641392D+01, VO
     A  1.19646919D+01, 1.20658131D+01, 1.21686452D+01, 1.22102927D+01, VO
     B  1.22360386D+01,      5*0.0D+00/                                 VO
      DATA TK_CrH/                                                      071215
     1 -0.999999993529,-0.986400006182,-0.965500000815,-0.911599980926, CrH
     2 -0.774199970016,-0.591300047671,-0.395999977999,-0.241100004338, CrH
     3 -0.082399980397, 0.067400012280, 0.201100040364, 0.404199970436, CrH
     4  0.583300043542, 0.720200033644, 0.849100024508, 0.999700000048, CrH
     5  1.146600103894, 1.313699945686, 1.484300006844, 1.818599941915, CrH
     6  2.329599927387, 2.494199990488, 2.647900092011, 2.817699937307, CrH
     7  2.978800013321, 3.075200062990, 3.173099999997, 3.355000029370, CrH
     8  3.453299873395, 3.550900173271, 3.615300005769, 3.679899812390, CrH
     9  3.756500141231, 3.825900024999, 3.925199957577, 3.970999797938, CrH
     A  3.987999829615, 4.000000000000,      8*0.0D+00/                 CrH
      DATA  K_CrH/                                                      071215
     1  3.20395078D+00, 3.17845904D+00, 3.14287179D+00, 3.06964365D+00, CrH
     2  2.98627063D+00, 3.04927162D+00, 3.26487026D+00, 3.50807355D+00, CrH
     3  3.80194080D+00, 4.10840036D+00, 4.39885908D+00, 4.86003497D+00, CrH
     4  5.27048165D+00, 5.56882590D+00, 5.82794794D+00, 6.10614299D+00, CrH
     5  6.36043504D+00, 6.63707996D+00, 6.91059330D+00, 7.43127692D+00, CrH
     6  8.20840089D+00, 8.45606902D+00, 8.68460796D+00, 8.92752290D+00, CrH
     7  9.13918531D+00, 9.25384187D+00, 9.36043516D+00, 9.53809840D+00, CrH
     8  9.62961795D+00, 9.71912292D+00, 9.77534414D+00, 9.82799906D+00, CrH
     9  9.88766509D+00, 9.94648659D+00, 1.00588603D+01, 1.01271833D+01, CrH
     A  1.01554010D+01, 1.01762401D+01,      8*0.0D+00/                 CrH
      DATA TK_CrO/                                                      071215
     1 -0.999999993529,-0.990100001508,-0.979000014924,-0.945599996888, CrO
     2 -0.860199999130,-0.762199992032,-0.648899979184,-0.558699997986, CrO
     3 -0.471200043472,-0.383099985052,-0.296800020192,-0.150999980685, CrO
     4  0.004800005309, 0.258899995806, 0.500200004716, 0.799499948708, CrO
     5  1.111899883193, 1.321399926154, 1.510699970976, 1.911200075749, CrO
     6  2.074899935521, 2.250800000408, 2.492799992611, 2.596399912749, CrO
     7  2.717700035112, 2.838899904066, 2.963600031669, 3.079499957247, CrO
     8  3.187099982349, 3.306700023376, 3.419100069243, 3.511600183724, CrO
     9  3.608100082886, 3.714999921002, 3.849199910587, 3.937799993370, CrO
     A  3.975999932872, 3.989299795295, 4.000000000000,      7*0.0D+00/ CrO
      DATA  K_CrO/                                                      071215
     1  7.49861218D+00, 7.42399273D+00, 7.34298944D+00, 7.11551312D+00, CrO
     2  6.63457156D+00, 6.23386690D+00, 5.93024363D+00, 5.78349311D+00, CrO
     3  5.70196000D+00, 5.66488287D+00, 5.66043549D+00, 5.70341655D+00, CrO
     4  5.80082744D+00, 6.03949301D+00, 6.32435171D+00, 6.72141626D+00, CrO
     5  7.16312926D+00, 7.46760355D+00, 7.74620777D+00, 8.35776014D+00, CrO
     6  8.62666988D+00, 8.92629722D+00, 9.33566775D+00, 9.50098247D+00, CrO
     7  9.68129773D+00, 9.84458953D+00, 9.99421774D+00, 1.01174507D+01, CrO
     8  1.02201990D+01, 1.03261853D+01, 1.04265671D+01, 1.05171237D+01, CrO
     9  1.06242742D+01, 1.07627549D+01, 1.09828505D+01, 1.11703146D+01, CrO
     A  1.12625968D+01, 1.12962315D+01, 1.13238260D+01,      7*0.0D+00/ CrO
      DATA TK_FeH/                                                      071215
     1 -0.999999993529,-0.984400008734,-0.959299978372,-0.895799970107, FeH
     2 -0.826399976062,-0.734200015973,-0.630100007402,-0.516599981935, FeH
     3 -0.397099974957,-0.271199971310,-0.145499979997,-0.009699998438, FeH
     4  0.255099999023, 0.346600040496, 0.449100045233, 0.619000019122, FeH
     5  0.755599991243, 0.894399980797, 1.087999928381, 1.245200006037, FeH
     6  1.414300077863, 1.641000090661, 1.859999907798, 1.986200008581, FeH
     7  2.113599882790, 2.338499922681, 2.441400062651, 2.547699957449, FeH
     8  2.672500076202, 2.792999970371, 2.980100013033, 3.076100040858, FeH
     9  3.171900028914, 3.277900166484, 3.382700044561, 3.566100080548, FeH
     A  3.692300015093, 3.880700101376, 3.952899891703, 3.982099985378, FeH
     B  4.000000000000,      5*0.0D+00/                                 FeH
      DATA  K_FeH/                                                      071215
     1  3.46389011D+00, 3.44635534D+00, 3.42232317D+00, 3.38264415D+00, FeH
     2  3.37001436D+00, 3.39504808D+00, 3.47014356D+00, 3.59668523D+00, FeH
     3  3.76882129D+00, 3.98297317D+00, 4.22204452D+00, 4.50140466D+00, FeH
     4  5.08857805D+00, 5.30016142D+00, 5.53983246D+00, 5.93341811D+00, FeH
     5  6.23220946D+00, 6.51046248D+00, 6.86288403D+00, 7.13021965D+00, FeH
     6  7.40658053D+00, 7.76491382D+00, 8.09544191D+00, 8.27379140D+00, FeH
     7  8.44316615D+00, 8.72502011D+00, 8.85200099D+00, 8.98302212D+00, FeH
     8  9.13347947D+00, 9.26869999D+00, 9.44105374D+00, 9.50576622D+00, FeH
     9  9.55330331D+00, 9.58870558D+00, 9.61164151D+00, 9.64437699D+00, FeH
     A  9.67761811D+00, 9.78602272D+00, 9.86123365D+00, 9.89843645D+00, FeH
     B  9.92325185D+00,      5*0.0D+00/                                 FeH
      DATA TK_FeO/                                                      071215
     1 -0.999999993529,-0.990000001588,-0.978800015039,-0.945099997285, FeO
     2 -0.858900004531,-0.761099995039,-0.646599981462,-0.555299995932, FeO
     3 -0.467100022367,-0.378699985067,-0.292700014630,-0.139999985005, FeO
     4  0.023699995203, 0.274899971481, 0.523099971728, 0.820700065333, FeO
     5  1.139000117314, 1.337699907393, 1.528499984712, 1.687000074103, FeO
     6  1.874399992874, 2.159600088090, 2.285299965834, 2.412400104558, FeO
     7  2.541299957372, 2.695800045229, 2.833199928625, 2.986600011099, FeO
     8  3.089200123628, 3.188200005376, 3.293800055642, 3.389499893666, FeO
     9  3.572500110129, 3.713499955318, 3.886799946212, 3.955499948475, FeO
     A  3.983099958978, 4.000000000000,      8*0.0D+00/                 FeO
      DATA  K_FeO/                                                      071215
     1  7.40379841D+00, 7.33332003D+00, 7.25696646D+00, 7.04295762D+00, FeO
     2  6.59297984D+00, 6.22608381D+00, 5.94986251D+00, 5.82095047D+00, FeO
     3  5.75390759D+00, 5.72837945D+00, 5.73255619D+00, 5.78996001D+00, FeO
     4  5.90409876D+00, 6.15272729D+00, 6.45319264D+00, 6.85319206D+00, FeO
     5  7.30610066D+00, 7.59581595D+00, 7.87720726D+00, 8.11432503D+00, FeO
     6  8.40412970D+00, 8.88072439D+00, 9.10723049D+00, 9.34604221D+00, FeO
     7  9.59444644D+00, 9.88914312D+00, 1.01358169D+01, 1.03827996D+01, FeO
     8  1.05282107D+01, 1.06530519D+01, 1.07708705D+01, 1.08669227D+01, FeO
     9  1.10374122D+01, 1.11745495D+01, 1.13842872D+01, 1.14896332D+01, FeO
     A  1.15362944D+01, 1.15661305D+01,      8*0.0D+00/                 FeO
      DATA TK_YO/                                                       071215
     1 -0.999999993529,-0.989900001716,-0.978100015442,-0.943599998475, YO
     2 -0.906599971279,-0.854800023501,-0.753000002171,-0.633500000738, YO
     3 -0.551999993939,-0.468100030404,-0.389499993899,-0.309100004361, YO
     4 -0.133399995854, 0.056499935650, 0.272899971191, 0.500400004165, YO
     5  0.753299989575, 1.006999979410, 1.231400012003, 1.451000042229, YO
     6  1.630300110524, 1.824199929671, 2.105499898608, 2.241200016661, YO
     7  2.365799895530, 2.520299976781, 2.671400078709, 2.793999967810, YO
     8  2.918200078617, 3.119499885215, 3.296200107291, 3.372999960817, YO
     9  3.445999895116, 3.611900082565, 3.679199829843, 3.746099900431, YO
     A  3.893299942487, 3.958000003064, 3.983999935217, 4.000000000000, YO
     B       6*0.0D+00/                                                 YO
      DATA  K_YO/                                                       071215
     1  1.09075828D+01, 1.07651242D+01, 1.06035792D+01, 1.01600702D+01, YO
     2  9.72925870D+00, 9.19695938D+00, 8.35836741D+00, 7.65513306D+00, YO
     3  7.30888856D+00, 7.03889303D+00, 6.84937567D+00, 6.70699893D+00, YO
     4  6.53685470D+00, 6.51146447D+00, 6.61336235D+00, 6.81377307D+00, YO
     5  7.10274025D+00, 7.43265316D+00, 7.74348994D+00, 8.05784603D+00, YO
     6  8.32022312D+00, 8.61341634D+00, 9.07079648D+00, 9.30805077D+00, YO
     7  9.53622302D+00, 9.82899185D+00, 1.01125802D+01, 1.03287376D+01, YO
     8  1.05277678D+01, 1.08023557D+01, 1.09985334D+01, 1.10739090D+01, YO
     9  1.11429133D+01, 1.13117879D+01, 1.13953364D+01, 1.14900750D+01, YO
     A  1.17226037D+01, 1.18220764D+01, 1.18597608D+01, 1.18821384D+01, YO
     B       6*0.0D+00/                                                 YO
      DATA TK_ZrO/                                                      071215
     1 -0.999999993529,-0.990000001588,-0.978600015154,-0.944699997602, ZrO
     2 -0.908499971526,-0.857700010083,-0.757399999579,-0.640499987503, ZrO
     3 -0.560899999555,-0.477700031539,-0.402799984231,-0.325199984021, ZrO
     4 -0.241900003445,-0.153099990732, 0.022399998361, 0.139299983539, ZrO
     5  0.259599995213, 0.484100015437, 0.739800021890, 1.000699997941, ZrO
     6  1.229300017963, 1.456800031912, 1.644800093754, 1.894100103889, ZrO
     7  2.182400081066, 2.329299927551, 2.459500026815, 2.644800094140, ZrO
     8  2.713000047900, 2.790499976773, 2.871399922676, 2.950800047864, ZrO
     9  3.114000039168, 3.245499899759, 3.375300008019, 3.567400110928, ZrO
     A  3.763500144939, 3.903999997802, 3.962699972274, 4.000000000000, ZrO
     B       6*0.0D+00/                                                 ZrO
      DATA  K_ZrO/                                                      071215
     1  1.18657231D+01, 1.17119681D+01, 1.15416539D+01, 1.10650658D+01, ZrO
     2  1.06025855D+01, 1.00272577D+01, 9.10998448D+00, 8.33701308D+00, ZrO
     3  7.95221784D+00, 7.64476546D+00, 7.43415793D+00, 7.26930835D+00, ZrO
     4  7.14206394D+00, 7.05342050D+00, 6.98976619D+00, 7.00963461D+00, ZrO
     5  7.06810454D+00, 7.25037701D+00, 7.53225414D+00, 7.86629697D+00, ZrO
     6  8.18070207D+00, 8.50526139D+00, 8.78014153D+00, 9.16095232D+00, ZrO
     7  9.63997957D+00, 9.90110083D+00, 1.01396734D+01, 1.04658557D+01, ZrO
     8  1.05711706D+01, 1.06760831D+01, 1.07682219D+01, 1.08436630D+01, ZrO
     9  1.09703048D+01, 1.10662157D+01, 1.11713481D+01, 1.13542903D+01, ZrO
     A  1.15601914D+01, 1.17014817D+01, 1.17553587D+01, 1.17874862D+01, ZrO
     B       6*0.0D+00/                                                 ZrO
      DATA TK_LaO/                                                      071215
     1 -0.999999993529,-0.989800001843,-0.977800015614,-0.965400000387, LaO
     2 -0.942999998951,-0.905799971175,-0.853500029516,-0.749300005515, LaO
     3 -0.629200007448,-0.547299981313,-0.464299999861,-0.384599987126, LaO
     4 -0.301700020765,-0.215699999572,-0.122699996130, 0.066199998191, LaO
     5  0.271999971061, 0.489000016703, 0.730700035286, 0.974100008980, LaO
     6  1.208900016485, 1.444000043153, 1.618199891161, 1.805999942574, LaO
     7  1.946000053733, 2.103599904926, 2.311099940765, 2.507499996950, LaO
     8  2.820099941746, 3.046400108923, 3.140899872075, 3.240099776859, LaO
     9  3.390099884664, 3.533200121779, 3.824000070539, 3.931299853095, LaO
     A  3.972899849213, 4.000000000000,      8*0.0D+00/                 LaO
      DATA  K_LaO/                                                      071215
     1  1.18641006D+01, 1.16988888D+01, 1.15102044D+01, 1.13214462D+01, LaO
     2  1.09960375D+01, 1.04974215D+01, 9.87730057D+00, 8.88370948D+00, LaO
     3  8.05693711D+00, 7.64241038D+00, 7.31759020D+00, 7.07866561D+00, LaO
     4  6.89245543D+00, 6.75506249D+00, 6.65932558D+00, 6.59428165D+00, LaO
     5  6.66117597D+00, 6.82962855D+00, 7.08869355D+00, 7.39413882D+00, LaO
     6  7.71317262D+00, 8.04654952D+00, 8.29997958D+00, 8.58208436D+00, LaO
     7  8.80235260D+00, 9.06199698D+00, 9.41521892D+00, 9.74682871D+00, LaO
     8  1.02547038D+01, 1.06174670D+01, 1.07731270D+01, 1.09397263D+01, LaO
     9  1.11911315D+01, 1.14177917D+01, 1.18075794D+01, 1.19260665D+01, LaO
     A  1.19679180D+01, 1.19938538D+01,      8*0.0D+00/                 LaO
C
C Length of idividual temperature grids
C
      DATA MTQ/ 24, 22, 22, 26, 26, 31, 31, 29, 29, 24, 24, 23, 22, 28,
     *  27, 25, 30, 30, 22, 29, 27, 33, 28, 28, 29, 29, 26, 24, 27, 30,
     *  30, 30, 32, 27, 26, 28, 29, 29, 22, 23, 29, 29, 26, 25, 29, 28,
     *  29, 27, 30, 27, 29, 28, 29, 27, 29, 30, 28, 27, 33, 25/
      DATA MTK/ 36, 35, 32, 39, 33, 44, 43, 38, 39, 34, 36, 36, 31, 38,
     *  46, 45, 37, 36, 36, 35, 39, 42, 42, 40, 41, 41, 41, 37, 40, 38,
     *  44, 38, 43, 40, 40, 45, 38, 35, 35, 39, 38, 36, 36, 34, 38, 39,
     *  36, 40, 37, 38, 39, 38, 41, 38, 39, 41, 38, 40, 40, 38/
C
      DATA FIRST/.TRUE./
C
C Compute 2nd derivatives for spline interpolation
C
      IF(FIRST) THEN
        DO 1 I=1,MSPEC
          CALL SPL_INIT(TQ(1,I),Q(1,I),Q2(1,I),U,MTQ(I))
          CALL SPL_INIT(TK(1,I),K(1,I),K2(1,I),U,MTK(I))
   1    CONTINUE
        FIRST=.FALSE.
      ENDIF
C
C Fits are made in log10 of temperatures
C
      TLOG=LOG10(TEMP)
C
C Find species name
C
      DO 4 II=1,MSPEC
        ISPEC=II
        IF(SPLIST(II).EQ.SPNAME) THEN
C
C  The species is in Barklem's list.
C  Find the braketing temperatures for the partition functions.
C
          KHI=MTQ(ISPEC)
          KLO=1
   2      CONTINUE
          I=(KLO+KHI)/2
          A=TQ(I,ISPEC)
          IF(A.GT.TLOG) THEN
            KHI=I
          ELSE IF(A.LE.TLOG) THEN
            KLO=I
          END IF
          IF(KHI-KLO.GT.1) GO TO 2
C
C Do the interpolation of the partition functions
C
          Q_spln=SPL_INTERP(KLO,KHI,TQ(1,ISPEC),Q(1,ISPEC),Q2(1,ISPEC),
     *                      MTQ(ISPEC),TLOG)
C  Find the braketing temperatures for the partition functions.
C
          KHI=MTK(ISPEC)
          KLO=1
   3      CONTINUE
          I=(KLO+KHI)/2
          A=TK(I,ISPEC)
          IF(A.GT.TLOG) THEN
            KHI=I
          ELSE IF(A.LE.TLOG) THEN
            KLO=I
          END IF
          IF(KHI-KLO.GT.1) GO TO 3
C
C Do the interpolation of the equilibrium constants
C
          K_spln=SPL_INTERP(KLO,KHI,TK(1,ISPEC),K(1,ISPEC),K2(1,ISPEC),
     *                      MTK(ISPEC),TLOG)
C
C The "+1" converts from pascals (N/m^2 as in Barklem tables) to
C dynes/cm^2 as required by the EOS.
C
          K_spln=K_spln+1.D0
          BARKLEM=.TRUE.
          if(SPNAME.eq.'CN') then
c            Q_spln=Q_spln+0.2
c            K_spln=K_spln+0.2
             BARKLEM=.FALSE.
          endif
          RETURN
        ENDIF
   4  CONTINUE
C
C Species was not found
C
      BARKLEM=.FALSE.
      RETURN
C
C End of computer-generated subroutine KP_Q_SPLN
      END


C=========================================================================
C=========================================================================
C
C NEGION: Returns partition function and ionization equilibrium for
C         a given negative ion and temperature.
C
C Inputs:
C   ANUM  [integer] atomic number.
C   TEMP  [real]    temperature (in K)
C   PARTN [real]    partition function of neutral atom
C
C                                     (3/2)              Eaffin
C   1     P(A)*P(e)        (2*Pi*m*kT)       2*U(A)    - ----
C  -- =   --------- = kT * -----------    *  ------ * e   kT
C  IT       P(A-)              h^3           U(A-)
C
C  U(A) is passed in as PARTN
C
C                           (3/2)
C  Const = k*(2*Pi*m_e*k/h^2)
C
C History:
C  10-dec-2007: First version written by N. Piskunov including 7 ions.
C               Partition functions tabulated by P. Barklem, resampled
C               for optimal spline interpolation and converted to Fortran
C               DATA statements by J. Valenti
C
C  15-dec-2007: Second version includes the same 7 negative ions tabulated
C               vs alog10(T) on adaptive grid similar to molecular species.
C
C   7-sep-2009: Subroutine data modified and the subroutine text generated
C               by IDL program qk_spl_nodes_f77.pro with errthr=0.000100
C
C Outputs:
C   Q_spln [real*8] partition functions at temperature T,
C          interpolated from Paul Barklem's tables;
C   IT     [real*8] computed according to the formula above.
C
C To obtain partition functions,Q:
C
C   D2 = SPL_INIT(TQ_<species>,Q_<species>)
C   Q(T) = SPL_INTERP(TQ_<species>,Q_<species>,D2,TLOG)
C
C Note that NEGION returns log10(Q)
C
C Reference:
C   Paul Barklem 2010, in preparation.
C
      SUBROUTINE NEGION(ANUM,TEMP,PARTN,IT,Q_atom,POTION,BARKLEM)
C
      IMPLICIT NONE
      INTEGER ANUM
      REAL TEMP,POTION
      REAL*8 PARTN,IT,Q_atom
      LOGICAL BARKLEM
C
C  Local variables
C
      LOGICAL FIRST
      INTEGER MSPEC,NTQ,KLO,KHI,I,II,ISPEC
      PARAMETER(MSPEC=7, NTQ=14)
      INTEGER MTQ(MSPEC)
      REAL*8 TLOG,A,U(14),SPL_INTERP,Const,TkeV,kBoleV
      PARAMETER(Const=0.3333984D0,kBoleV=8.6173175D-5)
C
      REAL*8 TQ(NTQ,MSPEC),Q(NTQ+1,MSPEC),Q2(NTQ,MSPEC)
      REAL*8           TQ_Hm   (NTQ  ),TQ_Cm   (NTQ  ),TQ_Om   (NTQ  ),
     * TQ_Fm   (NTQ  ),TQ_Sim  (NTQ  ),TQ_Sm   (NTQ  ),TQ_Clm  (NTQ  )
      REAL*8            Q_Hm   (NTQ+1), Q_Cm   (NTQ+1), Q_Om   (NTQ+1),
     *  Q_Fm   (NTQ+1), Q_Sim  (NTQ+1), Q_Sm   (NTQ+1), Q_Clm  (NTQ+1)
      EQUIVALENCE (TQ(1, 1),TQ_Hm   ),(TQ(1, 2),TQ_Cm   )
      EQUIVALENCE (TQ(1, 3),TQ_Om   ),(TQ(1, 4),TQ_Fm   )
      EQUIVALENCE (TQ(1, 5),TQ_Sim  ),(TQ(1, 6),TQ_Sm   )
      EQUIVALENCE (TQ(1, 7),TQ_Clm  )
      EQUIVALENCE ( Q(1, 1), Q_Hm   ),( Q(1, 2), Q_Cm   )
      EQUIVALENCE ( Q(1, 3), Q_Om   ),( Q(1, 4), Q_Fm   )
      EQUIVALENCE ( Q(1, 5), Q_Sim  ),( Q(1, 6), Q_Sm   )
      EQUIVALENCE ( Q(1, 7), Q_Clm  )
C
      INTEGER ATLIST(MSPEC)
      SAVE ATLIST,TQ,Q,Q2,FIRST,KHI,KLO
C
C                   H-  C-  O-  F- Si-  S- Cl-
      DATA ATLIST/  1,  6,  8,  9, 14, 16, 17/
C
C Tables of log10(T) and log10(Q)
C
      DATA TQ_Hm/                                                       071215
     1 -0.999999993529, 4.000000000000,     12*0.0D+00/                 Hm
      DATA  Q_Hm/                                                       071215
     1  0.00000000D+00, 0.00000000D+00,-7.54200000D-01,     12*0.0D+00/ Hm
      DATA TQ_Cm/                                                       071215
     1 -0.999999993529, 1.525599978474, 2.651200090844, 3.068100141729, Cm
     2  3.226999955268, 3.385099991304, 3.540199962140, 3.686999996896, Cm
     3  3.800200182895, 3.904799977325, 4.000000000000,      3*0.0D+00/ Cm
      DATA  Q_Cm/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02065454D-01, Cm
     2  6.02289485D-01, 6.05084791D-01, 6.19452167D-01, 6.56276557D-01, Cm
     3  7.02664478D-01, 7.55240436D-01, 8.06182861D-01,-1.26200000D+00, Cm
     4       3*0.0D+00/                                                 Cm
      DATA TQ_Om/                                                       071215
     1 -0.999999993529, 0.593600025261, 1.294199953625, 1.553499944587, Om
     2  1.796599952393, 1.944000048931, 2.095699919825, 2.317899934571, Om
     3  2.554199949293, 2.748599991156, 2.927100074615, 3.424799962155, Om
     4  3.799400175128, 4.000000000000/                                 Om
      DATA  Q_Om/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02060513D-01, 6.02235198D-01, Om
     2  6.05755401D-01, 6.13868967D-01, 6.29318181D-01, 6.61565642D-01, Om
     3  6.97424435D-01, 7.21774340D-01, 7.38756672D-01, 7.64722428D-01, Om
     4  7.72376824D-01, 7.74494613D-01,-1.46000000D+00/                 Om
      DATA TQ_Fm/                                                       071215
     1 -0.999999993529, 4.000000000000,     12*0.0D+00/                 Fm
      DATA  Q_Fm/                                                       071215
     1  0.00000000D+00, 0.00000000D+00,-3.40110000D+00,     12*0.0D+00/ Fm
      DATA TQ_Sim/                                                      071215
     1 -0.999999993529, 1.461700027344, 2.557999941888, 2.964800028584, Sim
     2  3.119799876818, 3.277800164079, 3.422000038438, 3.498999966735, Sim
     3  3.576999999319, 3.691600033946, 3.829999926729, 3.931199850937, Sim
     4  3.973199857309, 4.000000000000/                                 Sim
      DATA  Q_Sim/                                                      071215
     1  6.02059991D-01, 6.02059991D-01, 6.02059991D-01, 6.02081178D-01, Sim
     2  6.02612605D-01, 6.07739466D-01, 6.27640471D-01, 6.49433587D-01, Sim
     3  6.81270625D-01, 7.44441879D-01, 8.36538833D-01, 9.05670890D-01, Sim
     4  9.33332970D-01, 9.50480774D-01,-1.38900000D+00/                 Sim
      DATA TQ_Sm/                                                       071215
     1 -0.999999993529, 0.844800026460, 1.660300085079, 1.953300060737, Sm
     2  2.227700014231, 2.367499892458, 2.509300000613, 2.885400108161, Sm
     3  3.072500129387, 3.247699949829, 3.753500070131, 3.903600008040, Sm
     4  4.000000000000,      1*0.0D+00/                                 Sm
      DATA  Q_Sm/                                                       071215
     1  6.02059991D-01, 6.02059991D-01, 6.02060045D-01, 6.02153895D-01, Sm
     2  6.05580586D-01, 6.12903250D-01, 6.26569786D-01, 6.82012696D-01, Sm
     3  7.08431610D-01, 7.28330454D-01, 7.61103770D-01, 7.65938688D-01, Sm
     4  7.68312725D-01,-2.07700000D+00,      1*0.0D+00/                 Sm
      DATA TQ_Clm/                                                      071215
     1 -0.999999993529, 4.000000000000,     12*0.0D+00/                 Clm
      DATA  Q_Clm/                                                      071215
     1  0.00000000D+00, 0.00000000D+00,-3.61700000D+00,     12*0.0D+00/ Clm
C
C Length of idividual temperature grids
C
      DATA MTQ/  2, 11, 14,  2, 14, 13,  2/
C
      DATA FIRST/.TRUE./
C
C Compute 2nd derivatives for spline interpolation
C
      IF(FIRST) THEN
        DO 1 I=1,MSPEC
          CALL SPL_INIT(TQ(1,I),Q(1,I),Q2(1,I),U,MTQ(I))
   1    CONTINUE
        FIRST=.FALSE.
      ENDIF
C
C Fits are made in log10 of temperatures
C
      TLOG=LOG10(TEMP)
C
C Find species name
C
      DO 3 II=1,MSPEC
        ISPEC=II
        IF(ANUM.EQ.ATLIST(II)) THEN
C
C  The species is in Barklem's list.
C  Find the braketing temperatures for the partition functions.
C
          KHI=MTQ(ISPEC)
          KLO=1
   2      CONTINUE
          I=(KLO+KHI)/2
          A=TQ(I,ISPEC)
          IF(A.GT.TLOG) THEN
            KHI=I
          ELSE IF(A.LE.TLOG) THEN
            KLO=I
          END IF
          IF(KHI-KLO.GT.1) GO TO 2
C
C Do the interpolation of the partition functions
C
          Q_atom=SPL_INTERP(KLO,KHI,TQ(1,ISPEC),Q(1,ISPEC),Q2(1,ISPEC),
     *                      MTQ(ISPEC),TLOG)
          TkeV=kBoleV*TEMP
          Q_atom=10.d0**Q_atom
          POTION=-Q(MTQ(ISPEC)+1,ISPEC)
          IT=Const*(2.d0*PARTN)/Q_atom*EXP(-POTION/TkeV)*SQRT(TEMP)*
     *       TEMP*TEMP
          IT=1.D0/IT
          BARKLEM=.TRUE.
          RETURN
        ENDIF
   3  CONTINUE
C
C Species was not found
C
      Q_atom=1.D0
      IT=1.D-50
      BARKLEM=.FALSE.
      RETURN
C
C End of computer-generated subroutine NEGION
      END

C----------------------- End of Berklem subroutines ------------------------
      SUBROUTINE SPL_INIT(X,Y,Y2,U,N)
C
C  Computes second derivative approximations for cubic spline interpolation
C
      IMPLICIT NONE
      INTEGER N
      REAL*8 X(N),Y(N),Y2(N),U(N)
      INTEGER I
      REAL*8 SIG,P,YY1,YY2,YY3
C
C  Natural lower boundary condition
C
      Y2(1)=0.D0
      U(1)=0.D0
      DO 1 I=2,N-1
      SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
      P=SIG*Y2(I-1)+2.D0
      Y2(I)=(SIG-1.D0)/P
      YY1=Y(I-1)
      YY2=Y(I  )
      YY3=Y(I+1)
      U(I)=(6.D0*((YY3-YY2)/(X(I+1)-X(I))-(YY2-YY1)/
     /     (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
  1   CONTINUE
C
C  Natural upper boundary condition
C
      Y2(N)=0.D0
      DO 2 I=N-1,1,-1
  2   Y2(I)=Y2(I)*Y2(I+1)+U(I)
C
      RETURN
      END

      REAL*8 FUNCTION SPL_INTERP(KLO,KHI,XA,YA,Y2A,N,X)
C
C  Performs cubic spline interpolation
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N
      REAL*8 XA(N),YA(N),Y2A(N),X
      REAL*8 A,B,H,Y1,Y2
C
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y1=YA(KLO)
      Y2=YA(KHI)
      SPL_INTERP=A*Y1+B*Y2+((A*A-1.D0)*A*Y2A(KLO)+
     +                      (B*B-1.D0)*B*Y2A(KHI))*(H*H)/6.D0
C
      RETURN
      END

      SUBROUTINE XSAHA(IEL,TT,XNELEC,XNATOM,MAXION,POTI,FRCT,MODE)
C
C     MODE=1 returns ionization fractions/partition functions
C     MODE=2 returns ionization fractions
C     MODE=3 returns partition functions
C     MODE=4 returns total number of electrons produced
C     MODE=5 returns in MAXION(!) the number of ionization stages
C            available in XSAHA
C
C     ALL OF THE ABOVE IS FOR ALL IONIZATION STAGES UP TO MAXION
C
C  Parameters:
C     IEL    - (input) element atomic number (Hydrogen: 1)
C     TT     - (input) temperature (Kelvins)
C     XNELEC - (input) electron number density (cm^-3)
C     XNATOM - (input) particle number density (excluding electrons) (cm^-3)
C     MAXION - (input/output) size of the output arrays
C     POTI   - (output array of MAXION) ionization potential (eV)
C     FRCT   - (output array of MAXION) results according to MODE
C     MODE   - (input) see above
C
      INTEGER ELESIZ,IONSIZ,IEL
      PARAMETER (ELESIZ=100,IONSIZ=6)
      DOUBLE PRECISION FFF(IONSIZ),FEXARG,FRCT(MAXION),CF
      REAL IP(IONSIZ),PART(IONSIZ),POTLO(IONSIZ),SCALE(4),
     *     POTI(MAXION),TT
      INTEGER LOCZ(ELESIZ+1)
      LOGICAL FIRST

      INTEGER SIZ_H ,SIZ_He,SIZ_Li,SIZ_Be,SIZ_B ,SIZ_C ,SIZ_N ,SIZ_O ,
     1        SIZ_F ,SIZ_Ne,SIZ_Na,SIZ_Mg,SIZ_Al,SIZ_Si,SIZ_P ,SIZ_S ,
     2        SIZ_Cl,SIZ_Ar,SIZ_K ,SIZ_Ca,SIZ_Sc,SIZ_Ti,SIZ_V ,SIZ_Cr,
     3        SIZ_Mn,SIZ_Fe,SIZ_Co,SIZ_Ni,SIZ_Cu,SIZ_Zn,SIZ_Ga,SIZ_Ge,
     4        SIZ_As,SIZ_Se,SIZ_Br,SIZ_Kr,SIZ_Rb,SIZ_Sr,SIZ_Y ,SIZ_Zr,
     5        SIZ_Nb,SIZ_Mo,SIZ_Tc,SIZ_Ru,SIZ_Rh,SIZ_Pd,SIZ_Ag,SIZ_Cd,
     6        SIZ_In,SIZ_Sn,SIZ_Sb,SIZ_Te,SIZ_I ,SIZ_Xe,SIZ_Cs,SIZ_Ba,
     7        SIZ_La,SIZ_Ce,SIZ_Pr,SIZ_Nd,SIZ_Pm,SIZ_Sm,SIZ_Eu,SIZ_Gd,
     8        SIZ_Tb,SIZ_Dy,SIZ_Ho,SIZ_Er,SIZ_Tm,SIZ_Yb,SIZ_Lu,SIZ_Hf,
     9        SIZ_Ta,SIZ_W ,SIZ_Re,SIZ_Os,SIZ_Ir,SIZ_Pt,SIZ_Au,SIZ_Hg,
     A        SIZ_Tl,SIZ_Pb,SIZ_Bi,SIZ_Po,SIZ_At,SIZ_Rn,SIZ_Fr,SIZ_Ra,
     B        SIZ_Ac,SIZ_Th,SIZ_Pa,SIZ_U ,SIZ_Np,SIZ_Pu,SIZ_Am,SIZ_Cm,
     C        SIZ_Bk,SIZ_Cf,SIZ_Es
      INTEGER OFF_H ,OFF_He,OFF_Li,OFF_Be,OFF_B ,OFF_C ,OFF_N ,OFF_O ,
     1        OFF_F ,OFF_Ne,OFF_Na,OFF_Mg,OFF_Al,OFF_Si,OFF_P ,OFF_S ,
     2        OFF_Cl,OFF_Ar,OFF_K ,OFF_Ca,OFF_Sc,OFF_Ti,OFF_V ,OFF_Cr,
     3        OFF_Mn,OFF_Fe,OFF_Co,OFF_Ni,OFF_Cu,OFF_Zn,OFF_Ga,OFF_Ge,
     4        OFF_As,OFF_Se,OFF_Br,OFF_Kr,OFF_Rb,OFF_Sr,OFF_Y ,OFF_Zr,
     5        OFF_Nb,OFF_Mo,OFF_Tc,OFF_Ru,OFF_Rh,OFF_Pd,OFF_Ag,OFF_Cd,
     6        OFF_In,OFF_Sn,OFF_Sb,OFF_Te,OFF_I ,OFF_Xe,OFF_Cs,OFF_Ba,
     7        OFF_La,OFF_Ce,OFF_Pr,OFF_Nd,OFF_Pm,OFF_Sm,OFF_Eu,OFF_Gd,
     8        OFF_Tb,OFF_Dy,OFF_Ho,OFF_Er,OFF_Tm,OFF_Yb,OFF_Lu,OFF_Hf,
     9        OFF_Ta,OFF_W ,OFF_Re,OFF_Os,OFF_Ir,OFF_Pt,OFF_Au,OFF_Hg,
     A        OFF_Tl,OFF_Pb,OFF_Bi,OFF_Po,OFF_At,OFF_Rn,OFF_Fr,OFF_Ra,
     B        OFF_Ac,OFF_Th,OFF_Pa,OFF_U ,OFF_Np,OFF_Pu,OFF_Am,OFF_Cm,
     C        OFF_Bk,OFF_Cf,OFF_Es
C
C In order to add data for another ionization stage to a particular element
C one has to do two things: increase the value of SIZ_<elname> and add the
C data line(s) in the DATA NNN_<elname>
C
      PARAMETER (SIZ_H = 2, OFF_H = 1)
      INTEGER NNN_H (8*SIZ_H )
      PARAMETER (SIZ_He= 3, OFF_He=OFF_H +SIZ_H  )
      INTEGER NNN_He(8*SIZ_He)
      PARAMETER (SIZ_Li= 4, OFF_Li=OFF_He+SIZ_He)
      INTEGER NNN_Li(8*SIZ_Li)
      PARAMETER (SIZ_Be= 4, OFF_Be=OFF_Li+SIZ_Li)
      INTEGER NNN_Be(8*SIZ_Be)
      PARAMETER (SIZ_B = 4, OFF_B =OFF_Be+SIZ_Be)
      INTEGER NNN_B (8*SIZ_B )
      PARAMETER (SIZ_C = 6, OFF_C =OFF_B +SIZ_B )
      INTEGER NNN_C (8*SIZ_C )
      PARAMETER (SIZ_N = 6, OFF_N =OFF_C +SIZ_C )
      INTEGER NNN_N (8*SIZ_N )
      PARAMETER (SIZ_O = 6, OFF_O =OFF_N +SIZ_N )
      INTEGER NNN_O (8*SIZ_O )
      PARAMETER (SIZ_F = 6, OFF_F =OFF_O +SIZ_O )
      INTEGER NNN_F (8*SIZ_F )
      PARAMETER (SIZ_Ne= 6, OFF_Ne=OFF_F +SIZ_F )
      INTEGER NNN_Ne(8*SIZ_Ne)
      PARAMETER (SIZ_Na= 6, OFF_Na=OFF_Ne+SIZ_Ne)
      INTEGER NNN_Na(8*SIZ_Na)
      PARAMETER (SIZ_Mg= 6, OFF_Mg=OFF_Na+SIZ_Na)
      INTEGER NNN_Mg(8*SIZ_Mg)
      PARAMETER (SIZ_Al= 6, OFF_Al=OFF_Mg+SIZ_Mg)
      INTEGER NNN_Al(8*SIZ_Al)
      PARAMETER (SIZ_Si= 6, OFF_Si=OFF_Al+SIZ_Al)
      INTEGER NNN_Si(8*SIZ_Si)
      PARAMETER (SIZ_P = 6, OFF_P =OFF_Si+SIZ_Si)
      INTEGER NNN_P (8*SIZ_P )
      PARAMETER (SIZ_S = 6, OFF_S =OFF_P +SIZ_P )
      INTEGER NNN_S (8*SIZ_S )
      PARAMETER (SIZ_Cl= 5, OFF_Cl=OFF_S +SIZ_S )
      INTEGER NNN_Cl(8*SIZ_Cl)
      PARAMETER (SIZ_Ar= 5, OFF_Ar=OFF_Cl+SIZ_Cl)
      INTEGER NNN_Ar(8*SIZ_Ar)
      PARAMETER (SIZ_K = 5, OFF_K =OFF_Ar+SIZ_Ar)
      INTEGER NNN_K (8*SIZ_K )
      PARAMETER (SIZ_Ca= 5, OFF_Ca=OFF_K +SIZ_K )
      INTEGER NNN_Ca(8*SIZ_Ca)
      PARAMETER (SIZ_Sc= 5, OFF_Sc=OFF_Ca+SIZ_Ca)
      INTEGER NNN_Sc(8*SIZ_Sc)
      PARAMETER (SIZ_Ti= 5, OFF_Ti=OFF_Sc+SIZ_Sc)
      INTEGER NNN_Ti(8*SIZ_Ti)
      PARAMETER (SIZ_V = 5, OFF_V =OFF_Ti+SIZ_Ti)
      INTEGER NNN_V (8*SIZ_V )
      PARAMETER (SIZ_Cr= 5, OFF_Cr=OFF_V +SIZ_V )
      INTEGER NNN_Cr(8*SIZ_Cr)
      PARAMETER (SIZ_Mn= 5, OFF_Mn=OFF_Cr+SIZ_Cr)
      INTEGER NNN_Mn(8*SIZ_Mn)
      PARAMETER (SIZ_Fe= 5, OFF_Fe=OFF_Mn+SIZ_Mn)
      INTEGER NNN_Fe(8*SIZ_Fe)
      PARAMETER (SIZ_Co= 5, OFF_Co=OFF_Fe+SIZ_Fe)
      INTEGER NNN_Co(8*SIZ_Co)
      PARAMETER (SIZ_Ni= 5, OFF_Ni=OFF_Co+SIZ_Co)
      INTEGER NNN_Ni(8*SIZ_Ni)
      PARAMETER (SIZ_Cu= 3, OFF_Cu=OFF_Ni+SIZ_Ni)
      INTEGER NNN_Cu(8*SIZ_Cu)
      PARAMETER (SIZ_Zn= 3, OFF_Zn=OFF_Cu+SIZ_Cu)
      INTEGER NNN_Zn(8*SIZ_Zn)
      PARAMETER (SIZ_Ga= 3, OFF_Ga=OFF_Zn+SIZ_Zn)
      INTEGER NNN_Ga(8*SIZ_Ga)
      PARAMETER (SIZ_Ge= 3, OFF_Ge=OFF_Ga+SIZ_Ga)
      INTEGER NNN_Ge(8*SIZ_Ge)
      PARAMETER (SIZ_As= 3, OFF_As=OFF_Ge+SIZ_Ge)
      INTEGER NNN_As(8*SIZ_As)
      PARAMETER (SIZ_Se= 3, OFF_Se=OFF_As+SIZ_As)
      INTEGER NNN_Se(8*SIZ_Se)
      PARAMETER (SIZ_Br= 3, OFF_Br=OFF_Se+SIZ_Se)
      INTEGER NNN_Br(8*SIZ_Br)
      PARAMETER (SIZ_Kr= 3, OFF_Kr=OFF_Br+SIZ_Br)
      INTEGER NNN_Kr(8*SIZ_Kr)
      PARAMETER (SIZ_Rb= 3, OFF_Rb=OFF_Kr+SIZ_Kr)
      INTEGER NNN_Rb(8*SIZ_Rb)
      PARAMETER (SIZ_Sr= 3, OFF_Sr=OFF_Rb+SIZ_Rb)
      INTEGER NNN_Sr(8*SIZ_Sr)
      PARAMETER (SIZ_Y = 3, OFF_Y =OFF_Sr+SIZ_Sr)
      INTEGER NNN_Y (8*SIZ_Y )
      PARAMETER (SIZ_Zr= 3, OFF_Zr=OFF_Y +SIZ_Y )
      INTEGER NNN_Zr(8*SIZ_Zr)
      PARAMETER (SIZ_Nb= 3, OFF_Nb=OFF_Zr+SIZ_Zr)
      INTEGER NNN_Nb(8*SIZ_Nb)
      PARAMETER (SIZ_Mo= 3, OFF_Mo=OFF_Nb+SIZ_Nb)
      INTEGER NNN_Mo(8*SIZ_Mo)
      PARAMETER (SIZ_Tc= 3, OFF_Tc=OFF_Mo+SIZ_Mo)
      INTEGER NNN_Tc(8*SIZ_Tc)
      PARAMETER (SIZ_Ru= 3, OFF_Ru=OFF_Tc+SIZ_Tc)
      INTEGER NNN_Ru(8*SIZ_Ru)
      PARAMETER (SIZ_Rh= 3, OFF_Rh=OFF_Ru+SIZ_Ru)
      INTEGER NNN_Rh(8*SIZ_Rh)
      PARAMETER (SIZ_Pd= 3, OFF_Pd=OFF_Rh+SIZ_Rh)
      INTEGER NNN_Pd(8*SIZ_Pd)
      PARAMETER (SIZ_Ag= 3, OFF_Ag=OFF_Pd+SIZ_Pd)
      INTEGER NNN_Ag(8*SIZ_Ag)
      PARAMETER (SIZ_Cd= 3, OFF_Cd=OFF_Ag+SIZ_Ag)
      INTEGER NNN_Cd(8*SIZ_Cd)
      PARAMETER (SIZ_In= 3, OFF_In=OFF_Cd+SIZ_Cd)
      INTEGER NNN_In(8*SIZ_In)
      PARAMETER (SIZ_Sn= 3, OFF_Sn=OFF_In+SIZ_In)
      INTEGER NNN_Sn(8*SIZ_Sn)
      PARAMETER (SIZ_Sb= 3, OFF_Sb=OFF_Sn+SIZ_Sn)
      INTEGER NNN_Sb(8*SIZ_Sb)
      PARAMETER (SIZ_Te= 3, OFF_Te=OFF_Sb+SIZ_Sb)
      INTEGER NNN_Te(8*SIZ_Te)
      PARAMETER (SIZ_I = 3, OFF_I =OFF_Te+SIZ_Te)
      INTEGER NNN_I (8*SIZ_I )
      PARAMETER (SIZ_Xe= 3, OFF_Xe=OFF_I +SIZ_I )
      INTEGER NNN_Xe(8*SIZ_Xe)
      PARAMETER (SIZ_Cs= 3, OFF_Cs=OFF_Xe+SIZ_Xe)
      INTEGER NNN_Cs(8*SIZ_Cs)
      PARAMETER (SIZ_Ba= 3, OFF_Ba=OFF_Cs+SIZ_Cs)
      INTEGER NNN_Ba(8*SIZ_Ba)
      PARAMETER (SIZ_La= 3, OFF_La=OFF_Ba+SIZ_Ba)
      INTEGER NNN_La(8*SIZ_La)
      PARAMETER (SIZ_Ce= 4, OFF_Ce=OFF_La+SIZ_La)
      INTEGER NNN_Ce(8*SIZ_Ce)
      PARAMETER (SIZ_Pr= 4, OFF_Pr=OFF_Ce+SIZ_Ce)
      INTEGER NNN_Pr(8*SIZ_Pr)
      PARAMETER (SIZ_Nd= 4, OFF_Nd=OFF_Pr+SIZ_Pr)
      INTEGER NNN_Nd(8*SIZ_Nd)
      PARAMETER (SIZ_Pm= 3, OFF_Pm=OFF_Nd+SIZ_Nd)
      INTEGER NNN_Pm(8*SIZ_Pm)
      PARAMETER (SIZ_Sm= 3, OFF_Sm=OFF_Pm+SIZ_Pm)
      INTEGER NNN_Sm(8*SIZ_Sm)
      PARAMETER (SIZ_Eu= 4, OFF_Eu=OFF_Sm+SIZ_Sm)
      INTEGER NNN_Eu(8*SIZ_Eu)
      PARAMETER (SIZ_Gd= 3, OFF_Gd=OFF_Eu+SIZ_Eu)
      INTEGER NNN_Gd(8*SIZ_Gd)
      PARAMETER (SIZ_Tb= 3, OFF_Tb=OFF_Gd+SIZ_Gd)
      INTEGER NNN_Tb(8*SIZ_Tb)
      PARAMETER (SIZ_Dy= 3, OFF_Dy=OFF_Tb+SIZ_Tb)
      INTEGER NNN_Dy(8*SIZ_Dy)
      PARAMETER (SIZ_Ho= 3, OFF_Ho=OFF_Dy+SIZ_Dy)
      INTEGER NNN_Ho(8*SIZ_Ho)
      PARAMETER (SIZ_Er= 3, OFF_Er=OFF_Ho+SIZ_Ho)
      INTEGER NNN_Er(8*SIZ_Er)
      PARAMETER (SIZ_Tm= 3, OFF_Tm=OFF_Er+SIZ_Er)
      INTEGER NNN_Tm(8*SIZ_Tm)
      PARAMETER (SIZ_Yb= 3, OFF_Yb=OFF_Tm+SIZ_Tm)
      INTEGER NNN_Yb(8*SIZ_Yb)
      PARAMETER (SIZ_Lu= 3, OFF_Lu=OFF_Yb+SIZ_Yb)
      INTEGER NNN_Lu(8*SIZ_Lu)
      PARAMETER (SIZ_Hf= 3, OFF_Hf=OFF_Lu+SIZ_Lu)
      INTEGER NNN_Hf(8*SIZ_Hf)
      PARAMETER (SIZ_Ta= 3, OFF_Ta=OFF_Hf+SIZ_Hf)
      INTEGER NNN_Ta(8*SIZ_Ta)
      PARAMETER (SIZ_W = 3, OFF_W =OFF_Ta+SIZ_Ta)
      INTEGER NNN_W (8*SIZ_W )
      PARAMETER (SIZ_Re= 3, OFF_Re=OFF_W +SIZ_W )
      INTEGER NNN_Re(8*SIZ_Re)
      PARAMETER (SIZ_Os= 3, OFF_Os=OFF_Re+SIZ_Re)
      INTEGER NNN_Os(8*SIZ_Os)
      PARAMETER (SIZ_Ir= 3, OFF_Ir=OFF_Os+SIZ_Os)
      INTEGER NNN_Ir(8*SIZ_Ir)
      PARAMETER (SIZ_Pt= 3, OFF_Pt=OFF_Ir+SIZ_Ir)
      INTEGER NNN_Pt(8*SIZ_Pt)
      PARAMETER (SIZ_Au= 3, OFF_Au=OFF_Pt+SIZ_Pt)
      INTEGER NNN_Au(8*SIZ_Au)
      PARAMETER (SIZ_Hg= 3, OFF_Hg=OFF_Au+SIZ_Au)
      INTEGER NNN_Hg(8*SIZ_Hg)
      PARAMETER (SIZ_Tl= 3, OFF_Tl=OFF_Hg+SIZ_Hg)
      INTEGER NNN_Tl(8*SIZ_Tl)
      PARAMETER (SIZ_Pb= 3, OFF_Pb=OFF_Tl+SIZ_Tl)
      INTEGER NNN_Pb(8*SIZ_Pb)
      PARAMETER (SIZ_Bi= 3, OFF_Bi=OFF_Pb+SIZ_Pb)
      INTEGER NNN_Bi(8*SIZ_Bi)
      PARAMETER (SIZ_Po= 3, OFF_Po=OFF_Bi+SIZ_Bi)
      INTEGER NNN_Po(8*SIZ_Po)
      PARAMETER (SIZ_At= 3, OFF_At=OFF_Po+SIZ_Po)
      INTEGER NNN_At(8*SIZ_At)
      PARAMETER (SIZ_Rn= 3, OFF_Rn=OFF_At+SIZ_At)
      INTEGER NNN_Rn(8*SIZ_Rn)
      PARAMETER (SIZ_Fr= 3, OFF_Fr=OFF_Rn+SIZ_Rn)
      INTEGER NNN_Fr(8*SIZ_Fr)
      PARAMETER (SIZ_Ra= 3, OFF_Ra=OFF_Fr+SIZ_Fr)
      INTEGER NNN_Ra(8*SIZ_Ra)
      PARAMETER (SIZ_Ac= 3, OFF_Ac=OFF_Ra+SIZ_Ra)
      INTEGER NNN_Ac(8*SIZ_Ac)
      PARAMETER (SIZ_Th= 3, OFF_Th=OFF_Ac+SIZ_Ac)
      INTEGER NNN_Th(8*SIZ_Th)
      PARAMETER (SIZ_Pa= 3, OFF_Pa=OFF_Th+SIZ_Th)
      INTEGER NNN_Pa(8*SIZ_Pa)
      PARAMETER (SIZ_U = 3, OFF_U =OFF_Pa+SIZ_Pa)
      INTEGER NNN_U (8*SIZ_U )
      PARAMETER (SIZ_Np= 3, OFF_Np=OFF_U +SIZ_U )
      INTEGER NNN_Np(8*SIZ_Np)
      PARAMETER (SIZ_Pu= 3, OFF_Pu=OFF_Np+SIZ_Np)
      INTEGER NNN_Pu(8*SIZ_Pu)
      PARAMETER (SIZ_Am= 3, OFF_Am=OFF_Pu+SIZ_Pu)
      INTEGER NNN_Am(8*SIZ_Am)
      PARAMETER (SIZ_Cm= 3, OFF_Cm=OFF_Am+SIZ_Am)
      INTEGER NNN_Cm(8*SIZ_Cm)
      PARAMETER (SIZ_Bk= 3, OFF_Bk=OFF_Cm+SIZ_Cm)
      INTEGER NNN_Bk(8*SIZ_Bk)
      PARAMETER (SIZ_Cf= 3, OFF_Cf=OFF_Bk+SIZ_Bk)
      INTEGER NNN_Cf(8*SIZ_Cf)
      PARAMETER (SIZ_Es= 3, OFF_Es=OFF_Cf+SIZ_Cf)
      INTEGER NNN_Es(8*SIZ_Es)

      PARAMETER (NTABLE=OFF_Es+SIZ_Es-1)
      INTEGER NNNPFN(8,NTABLE)

      EQUIVALENCE (NNNPFN(1,OFF_H ),NNN_H (1))
      EQUIVALENCE (NNNPFN(1,OFF_He),NNN_He(1))
      EQUIVALENCE (NNNPFN(1,OFF_Li),NNN_Li(1))
      EQUIVALENCE (NNNPFN(1,OFF_Be),NNN_Be(1))
      EQUIVALENCE (NNNPFN(1,OFF_B ),NNN_B (1))
      EQUIVALENCE (NNNPFN(1,OFF_C ),NNN_C (1))
      EQUIVALENCE (NNNPFN(1,OFF_N ),NNN_N (1))
      EQUIVALENCE (NNNPFN(1,OFF_O ),NNN_O (1))
      EQUIVALENCE (NNNPFN(1,OFF_F ),NNN_F (1))
      EQUIVALENCE (NNNPFN(1,OFF_Ne),NNN_Ne(1))
      EQUIVALENCE (NNNPFN(1,OFF_Na),NNN_Na(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mg),NNN_Mg(1))
      EQUIVALENCE (NNNPFN(1,OFF_Al),NNN_Al(1))
      EQUIVALENCE (NNNPFN(1,OFF_Si),NNN_Si(1))
      EQUIVALENCE (NNNPFN(1,OFF_P ),NNN_P (1))
      EQUIVALENCE (NNNPFN(1,OFF_S ),NNN_S (1))
      EQUIVALENCE (NNNPFN(1,OFF_Cl),NNN_Cl(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ar),NNN_Ar(1))
      EQUIVALENCE (NNNPFN(1,OFF_K ),NNN_K (1))
      EQUIVALENCE (NNNPFN(1,OFF_Ca),NNN_Ca(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sc),NNN_Sc(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ti),NNN_Ti(1))
      EQUIVALENCE (NNNPFN(1,OFF_V ),NNN_V (1))
      EQUIVALENCE (NNNPFN(1,OFF_Cr),NNN_Cr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mn),NNN_Mn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Fe),NNN_Fe(1))
      EQUIVALENCE (NNNPFN(1,OFF_Co),NNN_Co(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ni),NNN_Ni(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cu),NNN_Cu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Zn),NNN_Zn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ga),NNN_Ga(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ge),NNN_Ge(1))
      EQUIVALENCE (NNNPFN(1,OFF_As),NNN_As(1))
      EQUIVALENCE (NNNPFN(1,OFF_Se),NNN_Se(1))
      EQUIVALENCE (NNNPFN(1,OFF_Br),NNN_Br(1))
      EQUIVALENCE (NNNPFN(1,OFF_Kr),NNN_Kr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rb),NNN_Rb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sr),NNN_Sr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Y ),NNN_Y (1))
      EQUIVALENCE (NNNPFN(1,OFF_Zr),NNN_Zr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Nb),NNN_Nb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mo),NNN_Mo(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tc),NNN_Tc(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ru),NNN_Ru(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rh),NNN_Rh(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pd),NNN_Pd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ag),NNN_Ag(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cd),NNN_Cd(1))
      EQUIVALENCE (NNNPFN(1,OFF_In),NNN_In(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sn),NNN_Sn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sb),NNN_Sb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Te),NNN_Te(1))
      EQUIVALENCE (NNNPFN(1,OFF_I ),NNN_I (1))
      EQUIVALENCE (NNNPFN(1,OFF_Xe),NNN_Xe(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cs),NNN_Cs(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ba),NNN_Ba(1))
      EQUIVALENCE (NNNPFN(1,OFF_La),NNN_La(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ce),NNN_Ce(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pr),NNN_Pr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Nd),NNN_Nd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pm),NNN_Pm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sm),NNN_Sm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Eu),NNN_Eu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Gd),NNN_Gd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tb),NNN_Tb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Dy),NNN_Dy(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ho),NNN_Ho(1))
      EQUIVALENCE (NNNPFN(1,OFF_Er),NNN_Er(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tm),NNN_Tm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Yb),NNN_Yb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Lu),NNN_Lu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Hf),NNN_Hf(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ta),NNN_Ta(1))
      EQUIVALENCE (NNNPFN(1,OFF_W ),NNN_W (1))
      EQUIVALENCE (NNNPFN(1,OFF_Re),NNN_Re(1))
      EQUIVALENCE (NNNPFN(1,OFF_Os),NNN_Os(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ir),NNN_Ir(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pt),NNN_Pt(1))
      EQUIVALENCE (NNNPFN(1,OFF_Au),NNN_Au(1))
      EQUIVALENCE (NNNPFN(1,OFF_Hg),NNN_Hg(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tl),NNN_Tl(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pb),NNN_Pb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Bi),NNN_Bi(1))
      EQUIVALENCE (NNNPFN(1,OFF_Po),NNN_Po(1))
      EQUIVALENCE (NNNPFN(1,OFF_At),NNN_At(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rn),NNN_Rn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Fr),NNN_Fr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ra),NNN_Ra(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ac),NNN_Ac(1))
      EQUIVALENCE (NNNPFN(1,OFF_Th),NNN_Th(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pa),NNN_Pa(1))
      EQUIVALENCE (NNNPFN(1,OFF_U ),NNN_U (1))
      EQUIVALENCE (NNNPFN(1,OFF_Np),NNN_Np(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pu),NNN_Pu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Am),NNN_Am(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cm),NNN_Cm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Bk),NNN_Bk(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cf),NNN_Cf(1))
      EQUIVALENCE (NNNPFN(1,OFF_Es),NNN_Es(1))
      SAVE NNNPFN,LOCZ,SCALE,FIRST,FFF
C      ( 1)( 2)  ( 3)( 4)  ( 5)( 6)  ( 7)( 8)  ( 9)(10)   ( IP )   G Ion REF
      DATA NNN_H/
     1 200020001,200020011,201620881,231228281,378953411, 1359502, 1,00,D+F H  1
     2 100010001,100010001,100010001,100010001,100010001, 1359500, 1,01/G   H  2
      DATA NNN_He/
     1 100010001,100010011,102111241,145022061,363059451, 2458104, 2,00,D+F He 1
     2 200020001,200020071,208524971,382669341,128222452, 5440302, 2,01,D+F He 2
     3 100010001,100010001,100010001,100010001,100010001, 5440300, 2,02/G   He 3
      DATA NNN_Li/
     1 200020011,201220481,212922881,258731081,394251691,  538901, 3,00,D+F Li 1
     2 100010001,100010201,126225521, 67216512,351165562, 7561907, 3,01,D+F Li 2
     3 200020001,200020211,227936571, 69610342,137217102,12241800, 3,02,D+F Li 3
     4 100010001,100010001,100010001,100010001,100010001,12241800, 3,03/G   Li 4
      DATA NNN_Be/
     1 100010051,104311441,131615641,190623681,298037691,  931900, 4,00,AEL Be 1
     2 200120231,211422771,249627631,309034911,398545051, 1820600, 4,01,AEL Be 2
     3 100010001,100010201,126225521, 67216512,351165562,15385000, 4,02,AEL Be 3
     4 200020001,200020011,201220661,223426161,332644691,21765700, 4,03/AEL Be 4
      DATA NNN_B/
     1 600060001,600560281,608761991,637466191,693973361,  829500, 5,00,AEL B  1
     2 100310831,132016901,214226411,315736741,419147071, 2514900, 5,01,AEL B  2
     3 200721061,233526401,297533311,369040481,440747651, 3792000, 5,02,AEL B  3
     4 100010001,100010001,100010001,100010001,100010001,25929800, 5,03/G   B  4
      DATA NNN_C/
     1 893292271, 96110042,105311262,126315202,196126432, 1125508, 6,00,D+F C  1
     2 595060251,620865751,713280191, 95712292,167623542, 2437501, 6,01,D+F C  2
     3 105513201,180324851,341851341, 88416332,296550722, 4787101, 6,02,D+F C  3
     4 204922771,262630421,350941931,494556971,644872001, 6447600, 6,03,D+F C  4
     5 100010001,100010001,100010001,100010001,100010001,39207700, 6,04,G   C  5
     6 200020001,200020001,200020001,200020001,200020001,48998100, 6,05/G   C  6
      DATA NNN_N/
     1 403141851,457051681,594071181, 92913362,203331152, 1452915, 7,00,D+F N  1
     2 919899541,107211512,124914302,182526232,403762662, 2959202, 7,01,D+F N  2
     3 596862721,684177081, 88110342,128317062,239334312, 4742501, 7,02,D+F N  3
     4 112816481,240733751,462068491,116419932,283736822, 7744900, 7,03,D+F N  4
     5 210124681,293634211,391145791,539862151,703178471, 9786200, 7,04,D+F N  5
     6 100010001,100010001,100010001,100010001,100010001,55205700, 7,05/G   N  6
      DATA NNN_O/
     1 874789691,924795711, 99410492,115213492,169022242, 1361307, 8,00,D+F O  1
     2 424151091,622874781, 91312832,221842502, 79914013, 3510711, 8,01,D+F O  2
     3  95610702,118113032,149619922,329761642,101914173, 5488500, 8,02,D+F O  3
     4 603567171,775391141,106612482,143716252,181420032, 7739300, 8,03,D+F O  4
     5 124420321,306943181,606281181,101712232,142916342,11387300, 8,04,D+F O  5
     6 215026541,323137551,421546491,508255151,594863811,13807900, 8,05/AEL O  6
      DATA NNN_F/
     1 575958511,589859231,595860671,636470031,815199581, 1741802, 9,00,D+F F  1
     2 900296401,102610802,113912542,152921152,318348952, 3498003, 9,01,D+F F  2
     3 469162651,791295541,121419552,402686872,154822203, 6264500, 9,02,D+F F  3
     4  99511422,129214572,170523002,320140922,498458762, 8713900, 9,03,D+F F  4
     5 615472711, 87710602,127215002,172919582,218624152,11421300, 9,04,D+F F  5
     6 135324181,377252001,661580261, 94410852,122613672,15711700, 9,05/AEL F  6
      DATA NNN_Ne/
     1 100010001,100010051,105313051,210239461, 74013022, 2155808,10,00,D+F Ne 1
     2 580158751,591759741,642687101,159332652, 64111533, 4106907,10,01,D+F Ne 2
     3  93510272,110411662,127116062,257647882, 75110223, 6350000,10,02,D+F Ne 3
     4 529774371, 94611322,135816202,188221442,240626682, 9701900,10,03,D+F Ne 4
     5 103312152,140616092,181320182,222224262,263128352,12630000,10,04,AEL Ne 5
     6 629178711, 98311802,136715512,173619202,210422892,15790900,10,05/AEL Ne 6
      DATA NNN_Na/
     1 200020001,200320211,207322131,253031421,417657451,  513802,11,00,D+F Na 1
     2 100010001,100010161,119621261, 50711872,246445382, 4728901,11,01,D+F Na 2
     3 580158751,591860351, 71813142,321968812,106014333, 7165000,11,02,D+F Na 3
     4  96910772,116012242,130714232,153916552,177118872, 9888000,11,03,D+F Na 4
     5 601386081,108812932,148916832,187820722,226624612,13836900,11,04,AEL Na 5
     6 105712442,144616652,189221182,234425702,279630222,17209000,11,05/AEL Na 6
      DATA NNN_Mg/
     1 100010011,101410621,118414581,204831781,509479731,  764404,12,00,D+F Mg 1
     2 200120051,202921001,226926901,368457091, 92814872, 1503101,12,01,D+F Mg 2
     3 100010001,100110611,177455431,176546012, 99718753, 8011905,12,02,D+F Mg 3
     4 579758751,591459501,600560591,611461681,622362781,10928900,12,03,AEL Mg 4
     5 100611232,120612752,134214102,147815462,161416822,14122900,12,04,AEL Mg 5
     6 674896701,121814462,167018942,211723412,256527892,18648900,12,05/AEL Mg 6
      DATA NNN_Al/
     1 558857701,583558761,593260591,635969541,796790971,  598400,13,00,D+F Al 1
     2 100310211,110313021,172828201, 55311252,215637942, 1882203,13,01,D+F Al 2
     3 200320201,208622331,250530971,410251081,611571211, 2844000,13,02,D+F Al 3
     4 100010001,100210881,207436531,523168101,838999681,11996000,13,03,D+F Al 4
     5 577758651,591259631,604461351,622563161,640764981,15377000,13,04,AEL Al 5
     6 103511582,124713242,140014772,155316292,170517812,19042000,13,05/AEL Al 6
      DATA NNN_Si/
     1 825189211, 95210052,106211532,134317202,237934082,  814913,14,00,D+F Si 1
     2 563057761,588160311,631768671,791097651,127817282, 1634000,14,01,D+F Si 2
     3 101110771,126716471,232438081, 71914052,262045302, 3346001,14,02,D+F Si 3
     4 200720521,217224081,284439171,551370951, 86810262, 4513000,14,03,D+F Si 4
     5 100010001,100210881,207436531,523168101,838999681,16672900,14,04,FAK Si 5
     6 575458521,591459851,610063201,672674071,843698661,20510900,14,05/AEL Si 6
      DATA NNN_P/
     1 402643441,496757481,658274401,833492941,103511532, 1048300,15,00,AEL P  1
     2 874497931,106011282,119812802,138415142,164717802, 1972000,15,01,AEL P  2
     3 564058061,604164611,709579551, 90410172,112912422, 3015500,15,02,AEL P  3
     4 100811411,149720221,280936121,441552181,602168241, 5135400,15,03,AEL P  4
     5 200420781,227025361,281430911,336936471,392542021, 6500700,15,04,AEL P  5
     6 100010001,100010001,100010001,100010001,100010001,22041300,15,05/G   P  6
      DATA NNN_S/
     1 822887891,930697831,102610932,121614492,185124742, 1035708,16,00,D+F S  1
     2 443056011,694982961, 96911522,144218572,227326892, 2339900,16,01,D+F S  2
     3  91610392,113512242,136416942,233429882,364242962, 3500000,16,02,D+F S  3
     4 560058861,633871081, 82410062,123314602,168619132, 4728900,16,03,D+F S  4
     5 104512901,177025421,375163021,122420462,286036742, 7250000,16,04,D+F S  5
     6 202321571,241428261,358355061, 78310152,124814802, 8802800,16,05/D+F S  6
      DATA NNN_Cl/
     1 538155931,571657911,598067191, 89013782,227737172, 1300916,17,00,D+F Cl 1
     2 873396771,104411072,118513532,175525872,406763932, 2379903,17,01,D+F Cl 2
     3 506569571, 87610522,134421682,439092662,182132573, 3990006,17,02,D+F Cl 3
     4  95110872,120013232,154921252,345149322,641378942, 5350000,17,03,D+F Cl 4
     5 558960371,677779341, 95311692,138816082,182720472, 6780000,17,04/D+F Cl 5
      DATA NNN_Ar/
     1 100010001,100010051,106913911,240147261, 90716112, 1575411,18,00,D+F Ar 1
     2 550256831,578158781,636585461,151530162, 58010303, 2762007,18,01,D+F Ar 2
     3  92110362,112412002,133216772,254443722, 76512833, 4090003,18,02,D+F Ar 3
     4 582082081,103112292,149920212,309750502,720793642, 5978900,18,03,D+F Ar 4
     5  97111072,123213982,172625622,463976582,106413633, 7500000,18,04/D+F Ar 5
      DATA NNN_K/
     1 200020011,200720361,211923291,280137141,525575741,  433803,19,00,D+F K  1
     2 100010001,100110341,135929551, 79119282,405274892, 3180905,19,01,D+F K  2
     3 554657081,581260301, 73012702,285363872,129023363, 4600005,19,02,D+F K  3
     4  96010862,118413212,180836632, 90321023,416863253, 6090000,19,03,D+F K  4
     5 657793361,119515082,195826322,352944302,533162332, 8259900,19,04/D+F K  5
      DATA NNN_Ca/
     1 100110061,104311741,145919971,294345051, 69010322,  611003,20,00,D+F Ca 1
     2 205822781,279234761,427553061,688994901,136319772, 1186701,20,01,D+F Ca 2
     3 100010001,100510821,168744821,130232522, 69012813, 5121003,20,02,D+F Ca 3
     4 555157161,585662471, 82816862, 42510013,168423663, 6700000,20,03,D+F Ca 4
     5  99411262,123814062,182930402,484766392, 84310223, 8438900,20,04/D+F Ca 5
      DATA NNN_Sc/
     1 924696691,105212282,151219062,240530032,368944512,  653900,21,00,AEL Sc 1
     2 190424662,297634542,391743752,482952832,573761912, 1280000,21,01,AEL Sc 2
     3 976799291,101110322,105810882,111911502,118112122, 2475000,21,02,AEL Sc 3
     4 100010001,100510821,168744821,130232522, 69012813, 7390000,21,03,FAK Sc 4
     5 555157161,585662471, 82816862, 42510013,168423663, 9200000,21,04/FAK Sc 5
      DATA NNN_Ti/
     1 181021172,260333222,430155582,710089242,110213293,  681900,22,00,D+F Ti 1
     2 474659872,721284672, 98211413,134515623,177919963, 1356900,22,01,D+F Ti 2
     3 228327012,308134272,381143862,534563472,734983512, 2747000,22,02,D+F Ti 3
     4 971498311, 99210032,102610572,108711172,114711782, 4324000,22,03,D+F Ti 4
     5 100010001,100510821,168744821,130232522, 69012813, 9980000,22,04/FAK Ti 5
      DATA NNN_V/
     1 272835172,425851532,632278322, 97212013,146817723,  674000,23,00,AEL V  1
     2 373954132,743597002,121414713,173920143,229225713, 1464900,23,01,AEL V  2
     3 323142642,519660272,679975352,824789522, 96610363, 2930900,23,02,AEL V  3
     4 248329302,324234952,373439752,421744582,469949412, 4800000,23,03,AEL V  4
     5 970698231,990699881,100710152,102410322,104010482, 6500000,23,04/AEL V  5
      DATA NNN_Cr/
     1 717277611, 92911652,152620872,295141952,550468122,  676400,24,00,D+F Cr 1
     2  71611552,205635512,558281952,115315823,205625293, 1649000,24,01,D+F Cr 2
     3 280639822,538369722, 87610823,129115003,170919183, 3095000,24,02,D+F Cr 3
     4 377150952,616070292,791788382, 97610683,116012523, 5000000,24,03,D+F Cr 4
     5 264730962,341436462,394042872,463549832,533056782, 7300000,24,04/D+F Cr 5
      DATA NNN_Mn/
     1 600060321,629270891, 86911302,151020222,267534752,  743100,25,00,AEL Mn 1
     2 739594821,139921212,309342852,567372412, 97112553, 1563600,25,01,AEL Mn 2
     3  98417472,265535782,454754842,641973532,828792212, 3369000,25,02,AEL Mn 3
     4 328847052,586668342,771785912, 94710343,112112093, 5300000,25,03,AEL Mn 4
     5 422055132,636770792,779285062,921999322,106411363, 7600000,25,04/AEL Mn 5
      DATA NNN_Fe/
C    1 197023222,274433302,416753952,723799822,139419053,  787038,26,00,D+F Fe 1
     1 197023222,274433302,416753952,723799822,139419053,  790024,26,00,D+F Fe 1! Ion. potential from NIST J. Sugar and C. Corliss, J. Phys. Chem. Ref. Data 14, 1-664 (1985).
     2 409453722,686687452,110213823,174322233,286437043, 1618792,26,01,D+F Fe 2! Kurucz
c    2 409453722,686687452,110213823,174322233,286437043, 1617902,26,01,D+F Fe 2
c    3 262136422,501167232, 87911303,138916483,190721673, 3064300,26,02,D+F Fe 3
     3 262136422,501167232, 87911303,138916483,190721673, 3065200,26,02,D+F Fe 3 ! Kurucz
     4  98723522,420363072, 87011423,145117913,215925463, 5700000,26,03,AEL Fe 4
     5 388854482,666275742,846693572,102511143,120312923, 7900000,26,04/D+F Fe 5
      DATA NNN_Co/
c    1 199427202,335740022,474957182,708090462,118315403,  786000,27,00,D+F Co 1
     1 199427202,335740022,474957182,708090462,118315403,  788100,27,00,D+F Co 1
     2 279739202,490858232,684582472,104713233,159818733, 1704900,27,01,D+F Co 2
     3 279836622,461857562,720693022,124915873,192522633, 3349000,27,02,D+F Co 3
     4 262136422,501167232, 87911303,138916483,190821673, 5300000,27,03,FAK Co 4
     5  98723522,420363072, 87011423,145117913,215925463, 8300000,27,04/FAK Co 5
      DATA NNN_Ni/
c    1 227027622,306233052,356839222,446052912,652382292,  763314,28,00,D+F Ni 1
     1 227027622,306233052,356839222,446052912,652382292,  763996,28,00,D+F Ni 1
     2 108416342,222428472,353944332,577378932,110314303, 1814900,28,01,D+F Ni 2
     3 198724282,293236452,468362702, 86511123,136016073, 3516000,28,02,D+F Ni 3
     4 279836622,461857562,720693022,124915873,192522633, 5600000,28,03,FAK Ni 4
     5 262136422,501167232, 87911303,138916483,190721673, 7900000,28,04/FAK Ni 5
      DATA NNN_Cu/
     1 201620781,231026761,314737361,450555381,692386911,  772301,29,00,D+F Cu 1
     2 109415761,247938311, 58910042,190937022, 68311693, 2028903,29,01,D+F Cu 2
     3 897195961,107212972,165021182,260230862,356940532, 3682900,29,02/D+F Cu 3
      DATA NNN_Zn/
     1 100010001,100410231,108712611,167124841,388460411,  939102,30,00,D+F Zn 1
     2 200020021,201620761,223726341,351352061, 80812472, 1796001,30,01,D+F Zn 2
     3 100610471,122617301,300566361,149924112,332342352, 3970000,30,02/D+F Zn 3
      DATA NNN_Ga/
     1 403245601,493151431,529654331,559358091,611065171,  600000,31,00,AEL Ga 1
     2  99710051,104511541,135016501,208226431,321837921, 2050900,31,01,AEL Ga 2
     3 199820071,204521391,229124761,266028451,302932131, 3070000,31,02/AEL Ga 3
      DATA NNN_Ge/
     1 502665261,755183501,901496201,102410942,117912812,  787900,32,00,AEL Ge 1
     2 422848161,512153401,557458941,636270361,794489061, 1593000,32,01,AEL Ge 2
     3 100010261,114613921,175221251,249828711,324436181, 3421000,32,02/AEL Ge 3
      DATA NNN_As/
     1 403143241,491856701,649173781,840396751,113013392,  981000,33,00,AEL As 1
     2 593676641,884697521,105911572,129515012,180322212, 1858700,33,01,AEL As 2
     3 484470541, 91510972,125614082,157017612,199722912, 2829900,33,02/AEL As 3
      DATA NNN_Se/
     1 630172361,799686381,919797221,102810942,117712832,  975000,34,00,AEL Se 1
     2 438055511,691582151, 94510732,121413672,152016732, 2150000,34,01,AEL Se 2
     3 651982921, 94610382,113212492,139515462,169718482, 3200000,34,02/AEL Se 3
      DATA NNN_Br/
     1 437347431,498951671,538559501, 74710812,169126672, 1183910,35,00,D+F Br 1
     2 705183611, 93510092,111614162,222932532,427652992, 2160000,35,01,D+F Br 2
     3 510869921, 87410312,123116552,236530712,377744832, 3590000,35,02/D+F Br 3
      DATA NNN_Kr/
     1 100010001,100010051,105012781,198535971, 65911422, 1399507,36,00,D+F Kr 1
     2 461049811,522254261,609088131,168935052, 68612253, 2455908,36,01,D+F Kr 2
     3 759990901,101911142,129017782,302856642, 99414333, 3690000,36,02/D+F Kr 3
      DATA NNN_Rb/
     1 200020011,200720361,211523021,269434141,459163351,  417502,37,00,D+F Rb 1
     2 100010001,100110321,129524961, 61014202,291753192, 2750004,37,01,D+F Rb 2
     3 473650891,533156051, 66810932,232950852, 99915303, 4000000,37,02/D+F Rb 3
      DATA NNN_Sr/
     1 100110041,104111741,146019721,281941411,607785251,  569202,38,00,D+F Sr 1
     2 202621931,255331271,384347931,624085761,122417632, 1102600,38,01,D+F Sr 2
     3 100010001,100110321,129524961, 61014202,291753192, 4300000,38,02/FAK Sr 3
      DATA NNN_Y/
c    1 791587851,100012192,155119942,254031782,389946932,  637900,39,00,AEL Y  1
     1 791587851,100012192,155119942,254031782,389946932,  621710,39,00,AEL Y  1 ! From Kurucz
     2 118217102,220827002,319036792,416646512,513256072, 1223000,39,01,AEL Y  2
     3  92510012,104710862,112311612,120212472,132814282, 2050000,39,02/AEL Y  3
      DATA NNN_Zr/
     1 141320802,291439702,531170262, 92712273,162521053,  663400,40,00,D+F Zr 1 ! Ion. potential from NIST P.A. Hackett, M.R. Humphries, S.A. Mitchell, and D.M. Rayner, J. Chem. Phys. 85, 3194-3197 (1986)
     2 354454352,724689652,107212643,148517093,193321573, 1312900,40,01,D+F Zr 2
     3 209727032,324537052,415446282,510255752,604965222, 2298000,40,02/D+F Zr 3
      DATA NNN_Nb/
c    1 256636022,465759302,749693962,116514243,171520333,  687900,41,00,AEL Nb 1
     1 256636022,465759302,749693962,116514243,171520333,  675890,41,00,AEL Nb 1 ! From Kurucz
     2 335157222, 84511463,147718363,221826083,299933893, 1431900,41,01,AEL Nb 2
     3 223725352,280830972,340937362,406844002,473150632, 2503900,41,02/AEL Nb 3
      DATA NNN_Mo/
c    1 703972941, 82610822,154822682,327244912,571469372,  709900,42,00,D+F Mo 1
     1 703972941, 82610822,154822682,327244912,571469372,  709250,42,00,D+F Mo 1 ! From Kurucz
     2  69113342,270146932, 71810043,131916543,200323603, 1614900,42,01,NPk Mo 2 ! PFs are calculated using energy levels from Nilsson & Pickering, 2003, Phys. Scr., 67, 223
     3 267645462,669890262,115514323,173620673,242528083, 2714900,42,02/AEL Mo 3
      DATA NNN_Tc/
     1  90113722,190525812,348647032,631684102,110714373,  728000,43,00,Pal Tc 1 ! PFs are taken from Palmeri et al. 2007, MNRAS, 374, 63
     2 132521482,335250142, 72110033,135517843,229929083, 1525900,43,01,Pal Tc 2 ! PFs are taken from Palmeri et al. 2007, MNRAS, 374, 63
     3  80117462,174618952,189518952,189518952,189518952, 3000000,43,02/Pal Tc 3 ! PFs are taken from Palmeri et al. 2007, MNRAS, 374, 63
      DATA NNN_Ru/
     1 176824122,318941082,515263202,761790472,106112303,  736400,44,00,AEL Ru 1
     2 221934642,501968372, 88911173,136316243,189221613, 1675900,44,01,AEL Ru 2
     3 210622722,241025422,267928262,297731272,327834282, 2846000,44,02/AEL Ru 3
      DATA NNN_Rh/
     1 148520202,255230902,364942462,489656082,638872352,  746000,45,00,AEL Rh 1
     2 153421292,288137912,484660322,720187062,101011483, 1807000,45,01,AEL Rh 2
     3 254537212,492362292,770592182,107312243,137615273, 3104900,45,02/AEL Rh 3
      DATA NNN_Pd/
     1 115919651,320746011,607576761, 95011642,141817172,  832900,46,00,AEL Pd 1
     2 755087211,105913442,173122222,282034722,412247732, 1941900,46,01,AEL Pd 2
     3 180223462,289735212,414247632,538460052,662672472, 3292000,46,02/AEL Pd 3
      DATA NNN_Ag/
     1 200020001,200220141,206422141,257633021,455164681,  757403,47,00,D+F Ag 1
     2 100810581,125817401,260641031, 66210072,135316982, 2148000,47,01,D+F Ag 2
     3 795887491, 97711762,156620252,248329422,340038582, 3481900,47,02/D+F Ag 3
      DATA NNN_Cd/
     1 100010001,100410241,109212891,176827421,444268771,  899003,48,00,D+F Cd 1
     2 200020021,201720921,233329881,451475371,127520782, 1690301,48,01,D+F Cd 2
     3 100310281,114815371,246138311,519265531,791492761, 3747000,48,02/D+F Cd 3
      DATA NNN_In/
     1 252431921,368440461,433746521,512259221,723389021,  578400,49,00,D+F In 1
     2 100110071,104611651,146118581,225426511,304734431, 1886000,49,01,D+F In 2
     3 200120111,205021611,243628031,317035371,390442701, 2802900,49,02/D+F In 3
      DATA NNN_Sn/
     1 232637101,488058571,669074381,816189091, 97210632,  734200,50,00,AEL Sn 1
     2 286335941,408144471,479351961,571862901,686274341, 1462700,50,01,AEL Sn 2
     3 100010251,114013811,175321601,256829751,338337901, 3049000,50,02/AEL Sn 3
      DATA NNN_Sb/
     1 404043481,494656811,646772781,813490751,101411372,  863900,51,00,AEL Sb 1
     2 303147981,618472951,827392621,103711702,131214532, 1650000,51,01,AEL Sb 2
     3 313037601,429347901,536260591,689477591,862494881, 2529900,51,02/AEL Sb 3
      DATA NNN_Te/
     1 526258801,657372351,784284071,897095741,102711082,  900900,52,00,AEL Te 1
     2 440855541,686481251, 93810792,125414792,176321132, 1860000,52,01,AEL Te 2
     3 349054751,699883081, 96611302,134216202,197724212, 2800000,52,02/AEL Te 3
      DATA NNN_I/
     1 405342041,438645621,475751071,587974491,102214572, 1045404,53,00,D+F I  1
     2 568567471,773485861, 94510362,112712182,130914002, 1909000,53,01,D+F I  2
     3 514269581, 86910562,130716652,215327742,351843662, 3200000,53,02/AEL I  3
      DATA NNN_Xe/
     1 100010001,100010091,109515351,291060661,119621482, 1212716,54,00,D+F Xe 1
     2 414844131,465649111,538464651, 87112232,158019362, 2120000,54,01,D+F Xe 2
     3 615475101,867797531,112213462,157618062,203622662, 3209900,54,02/D+F Xe 3
      DATA NNN_Cs/
     1 200020001,201020501,215623871,283536181,462756261,  389300,55,00,D+F Cs 1
     2 100010001,100310371,119016501,269146361, 77912412, 2510000,55,01,D+F Cs 2
     3 424445601,481750061,516953311,549356551,581759791, 3500000,55,02/D+F Cs 3
      DATA NNN_Ba/
     1 101210791,135119351,282340571,574580391,111015062,  521002,56,00,D+F Ba 1
     2 262638611,504160621,698579371, 91010692,129115952, 1000000,56,01,D+F Ba 2
     3 100010001,100310351,118416321,264945521, 76512182, 3700000,56,02/FAK Ba 3
      DATA NNN_La/
     1  71111992,172323592,312540402,510763182,765791012,  557700,57,00,AEL La 1
     2 204529582,383647882,582469262,807992692,104911723, 1106000,57,01,AEL La 2
     3  94712552,148416582,179819212,203621522,227424042, 1917700,57,02/AEL La 3
      DATA NNN_Ce/
     1 516771922,101415733,230431963,422563713,661579353,  553870,58,00,AEL Ce 1 ! PFs are taken from Palmeri et al. 2000, Phys. Scr., 61, 323
     2  71918863,305242193,538665523,771988853,100511224, 1085000,58,01,MZH Ce 2 ! PFs are taken from Palmeri et al. 2000, Phys. Scr., 61, 323
     3 506183092,108612923,146416133,174418603,196520603, 2020000,58,02,CCB Ce 3 ! PFs are taken from Cowley & Barisciano 1994, Obs., 114, 308
     4 118012722,134214202,152616852,191722342,264131332, 3690600,58,03/RW  Ce 4 ! PFs are calculated using energy levels from Reader & Wyart 2009, Phys. Rev. A, 80, 042517
      DATA NNN_Pr/
     1 146526632,508289352,142720943,287237333,465456163,  547300,59,00,Sne Pr 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2  53615083,324256453, 86012064,159720354,251930474, 1055000,59,01,ISA Pr 2 ! PFs are calculated using energy levels from Mashonkina et al. 2009, A&A, 495, 297
     3 421093902,165924663,331041793,507660143,700980743, 2162400,59,02,ISA Pr 3 ! PFs are calculated using energy levels from Mashonkina et al. 2009, A&A, 495, 297
     4 373649462,593368882,785988552, 98810923,119813043, 3900000,59,03/AEL Pr 4 ! PFs are calculated using NIST energy levels
      DATA NNN_Nd/
     1 145623072,410172132,120218793,276138313,505263693,  552500,60,00,Sne Nd 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2  47511303,223037433,559777223,100512564,151817894, 1073000,60,01,ISA Nd 2 ! PFs are calculated using energy levels from Mashonkina et al. 2005, A&A, 441, 309
     3 432699302,204835193,525971403, 90710984,128314614, 2218000,60,02,ISA Nd 3 ! PFs are calculated using energy levels from Ryabchikova et al. 2006, A&A, 456, 329
     4 104717683,241529543,339937663,407343323,455447453, 4042000,60,03/Wyt Nd 4 ! PFs are calculated using energy levels from Wyart et al. 2006, J. Phys. B39, L77
      DATA NNN_Pm/
     1 293029302,339657372, 97415223,219529733,383647633,  558200,61,00,Fiv Pm 1 ! PFs are taken from Fivet at al. 2007, MNRAS, 380, 771
     2  53611273,274552953, 86912833,176222974,288035004, 1090000,61,01,Fiv Pm 2 ! PFs are taken from Fivet at al. 2007, MNRAS, 380, 771
     3  49012373,262048233,482348233,519661563,709279783, 2230000,61,02/Fiv Pm 3 ! PFs are taken from Fivet at al. 2007, MNRAS, 380, 771
      DATA NNN_Sm/
     1  92915672,222431062,444763802, 89612173,159520253,  564370,62,00,AEL Sm 1
     2 315059662, 97114563,204627093,342541693,490556383, 1106900,62,01,AEL Sm 2
     3 269037812,520270372, 91111273,133915483,172719093, 2340000,62,02/AEL Sm 3
      DATA NNN_Eu/
     1 800080571,851699301,127617362,240433032,444958442,  567045,63,00,AEL Eu 1
     2 125416052,211828182,375549622,644381732,101112213, 1124100,63,01,AEL Eu 2
     3  82514782, 47913863,315459503, 98114674,204226924, 2492000,63,02,ISA Eu 3 ! PFs are calculated using energy levels from Wyart et al. 2008, A&A, 483, 339
     4 353543472,487852542,553557522,592460632,617962762, 4265000,63,03/AEL Eu 4 ! PFs are calculated using NIST energy levels
      DATA NNN_Gd/
     1 244232982,441460242, 82611223,149719523,247930643,  615000,64,00,Sne Gd 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 534793262,139219123,247730843,371043333,495055893, 1209000,64,01,AEL Gd 2
     3 364145232,514756362,604864112,673870372,732276072, 2063000,64,02/AEL Gd 3
      DATA NNN_Tb/
     1 546880382,113515623,209227313,347543173,524362333,  586390,65,00,Sne Tb 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2  56510823,163922043,279234353,417550623,615575303, 1151900,65,01,Sne Tb 2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3  53713323,276551143, 85012894,181224014,304037114, 2191000,65,02/ISA Tb 3 ! PFs are calculated using Wyart & Ryabtsev extended energy levels analysis (Ryabtsev, private communication)
      DATA NNN_Dy/
     1 175219662,262038952,604693902,142320733,288338103,  593890,66,00,Sne Dy 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 347359162,108619003,300742453,533359923,606555733, 1167000,66,01,Sne Dy 2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3 320279972,191238513, 66810374,148019834,253331184, 2280000,66,02/ISA Dy 3 ! PFs are calculated using Wyart & Ryabtsev extended energy levels analysis (Ryabtsev, private communication)
      DATA NNN_Ho/
     1 222635002,542276772,100312353,145716713,187020703,  602160,67,00,FAK Ho 1
     2 321455092,112322203,401966563,102014674,200226144, 1180000,67,01,Bor Ho 2 ! PFs are taken from Bord & Cowley 2002, Sol. Phys., 211, 3
     3 222635002,542276772,100312353,145716713,187020703, 2284000,67,02/AEL Ho 3
      DATA NNN_Er/
     1 131715322,213632462,504577482,115416533,226829683,  610780,68,00,Sne Er 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 282946962, 81713443,201827463,339638403,399938623, 1193000,68,01,Sne Er 2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3 801281851, 91511592,166126662,472591362,190642503, 2274000,68,02/Irw Er 3 ! PFs are calculated using polynomial approximation from Irwin 1981, ApJS, 45, 621
      DATA NNN_Tm/
     1 800381111, 87510702,147621462,310343462,585475982,  618436,69,00,AEL Tm 1
     2 156718872,279244452,678196342,128316243,197823443, 1205000,69,01,AEL Tm 2
     3  93517192,364666132,103414613,192624193,293334613, 2368000,69,02/AEL Tm 3
      DATA NNN_Yb/
     1 104410001,100011021,142920191,299545391, 68910342,  625394,70,00,Sne Yb 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 200120901,270345231, 81714042,223533112,461959862, 1218400,70,01,AEL Yb 2
     3 100312561,250851931, 91914182,198626022,323638692, 2505000,70,02/AEL Yb 3
      DATA NNN_Lu/
     1 514664441,759086851, 99211442,133315612,182721252,  542589,71,00,AEL Lu 1
     2 125924831,438667801, 98714112,199727872,380850742, 1389900,71,01,AEL Lu 2
C     2 112718911,335853801,742987841,895879721,626944081, 1389900,71,01,Sne Lu 2
     3 323948621,661297271,158626482,426865032, 93712843, 2095960,71,02/AEL Lu 3
      DATA NNN_Hf/
     1 659294081,128016962,222528952,372047062,585171462,  700000,72,00,AEL Hf 1
     2  99117882,274638812,520867322, 84410313,123314453, 1489900,72,01,AEL Hf 2
     3 187427702,343739872,448049452,539358282,625266642, 2329900,72,02/AEL Hf 3
      DATA NNN_Ta/
     1  65210892,171325762,373552252,705192012,116414343,  787900,73,00,AEL Ta 1
     2 192837842,600784802,111113823,165419233,218524383, 1620000,73,01,AEL Ta 2
     3  99117872,274638812,520867312, 84410313,123314453, 2400000,73,02/FAK Ta 3
      DATA NNN_W/
     1 398981651,130019172,273438022,516168382, 88411163,  797900,74,00,AEL W  1
     2 131429482,523279952,111414623,183422233,262130233, 1770000,74,01,AEL W  2
     3 192837842,600784792,111113823,165419233,218524383, 2500000,74,02/FAK W  3
      DATA NNN_Re/
     1 600963001, 75910412,150121572,301940972,539168952,  787000,75,00,AEL Re 1
     2  73710852,190731262,464964142, 83810503,127315053, 1660000,75,01,AEL Re 2
     3 131429482,523279952,111414623,183422233,262130233, 2600000,75,02/FAK Re 3
      DATA NNN_Os/
     1 110815502,216829732,398752322,672484682,104612673,  850000,76,00,AEL Os 1
     2 168225972,362046562,566766422,757484612, 93010103, 1700000,76,01,AEL Os 2
     3  73710852,190731262,464964142, 83810503,127315053, 2700000,76,02/FAK Os 3
      DATA NNN_Ir/
     1 128117692,236030402,381847322,582671422, 87110533,  896700,77,00,AEL Ir 1 ! IP=8.96702 eV according to NIST
     2 216133402,476163702,811599542,118413753,156417503, 1691000,77,01,VKM Ir 2 ! PFs are calculated from energy levels of van Kleef & Metsch 1978, Physica C95, 251; IP=16.91 eV from Carlson et al. 1970, Atomic Data and Nuclear Data Table, 2, 63
     3 168225972,362046562,566766422,757484612, 93010103, 2800000,77,02/FAK Ir 3
      DATA NNN_Pt/
     1 158918512,207523002,254328242,316335762,407246582,  900000,78,00,AEL Pt 1
     2  98115462,224930742,401150612,623475412, 89910583, 1855900,78,01,AEL Pt 2
     3 110815502,216829732,398752322,672484682,104612673, 2900000,78,02/FAK Pt 3
      DATA NNN_Au/
     1 203222611,265731251,364042301,494958601,702084731,  922000,79,00,AEL Au 1
     2 120521331,357753801, 75310062,130516572,206925452, 2050000,79,01,AEL Au 2
     3 651780821,108814772,195925252,316338622,460853882, 3000000,79,02/AEL Au 3
      DATA NNN_Hg/
     1 100010001,100110111,105211851,152122101,341552811, 1043002,80,00,D+F Hg 1
     2 200320211,210023021,268834231,480472341,111416912, 1875000,80,01,D+F Hg 2
     3 104012871,186129471,458664151, 82410072,119013732, 3420000,80,02/D+F Hg 3
      DATA NNN_Tl/
     1 200420711,222424271,265429161,325637371,442853911,  610500,81,00,AEL Tl 1
     2 100010021,101910801,121414641,189525811,358949721, 2041900,81,01,AEL Tl 2
     3 200020311,216624611,296337451,489064791, 85711212, 2979900,81,02/AEL Tl 3
      DATA NNN_Pb/
     1 103411711,147819101,244331781,434862751, 93113762,  741404,82,00,D+F Pb 1
     2 204122231,248227841,311535621,429153941,651976431, 1502800,82,01,D+F Pb 2
     3 100210131,106812201,154522671,381665951, 95512512, 3192900,82,02/D+F Pb 3
      DATA NNN_Bi/
     1 400140351,416944121,474851591,564362181,690477231,  728700,83,00,AEL Bi 1
     2 106814451,204427341,350744811,586879131,108314772, 1667900,83,01,AEL Bi 2
     3 205523051,264830231,345439921,469156001,675281671, 2555900,83,02/AEL Bi 3
      DATA NNN_Po/
     1 500950661,518153561,559058941,628968071,748483501,  843000,84,00,AEL Po 1
     2 443756241,696282451, 95411012,128615262,182922012, 1900000,84,01,FAK Po 2
     3 336953201,682481011, 93810882,127915272,184622442, 2700000,84,02/FAK Po 3
      DATA NNN_At/
     1 402841621,431544771,463148311,520059491,734896851,  930000,85,00,FAK At 1
     2 576168741,788387631, 96910642,116012552,135014462, 2000000,85,01,FAK At 2
     3 490265341,812797201,116614322,179622692,285035302, 2900000,85,02/FAK At 3
      DATA NNN_Rn/
     1 100010001,100010031,102311051,133018071,264539391, 1074500,86,00,AEL Rn 1
     2 402841621,431544771,463148311,520059491,734996851, 2000000,86,01,FAK Rn 2
     3 576168741,788387631, 96910642,116012552,135014462, 3000000,86,02/FAK Rn 3
      DATA NNN_Fr/
     1 200020011,201220591,218124481,296538611,488859141,  400000,87,00,FAK Fr 1
     2 100010001,100010031,102311051,133018071,264539401, 2200000,87,01,FAK Fr 2
     3 421645151,477449611,511852711,542455761,572958821, 3300000,87,02/FAK Fr 3
      DATA NNN_Ra/
     1 104110411,105712431,155420871,293741981,596683361,  527800,88,00,Qui Ra 1 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     2 198321961,258631331,381946231,552565051,754486211, 1015000,88,01,Qui Ra 2 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     3 100010001,100010031,102311051,133018071,264539391, 3400000,88,02/FAK Ra 3
      DATA NNN_Ac/
     1 441654441,664281721,101912862,163320772,263333182,  517000,89,00,Qui Ac 1 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     2 195142621, 72610952,153420412,261732632,397747612, 1175000,89,01,Qui Ac 2 ! PFs are taken from Quinet et al. 2007, A&A, 474, 307
     3 723989131,103511752,130814352,155416652,177018682, 2000000,89,02/AEL Ac 3
      DATA NNN_Th/
     1  63810522,177929162,457168312, 97513353,175722323,  630670,90,00,Sne Th 1 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     2 167142052, 79912843,186125143,322539763,475155383, 1190000,90,01,BWt Th 2 ! PFs are calculated from 508 energy levels of Blaise & Wyart 1992, Energy Levels and Atomic Spectra of Actinides, Paris
     3 491281082,108913303,154717483,193921253,230924903, 1830000,90,02/BWt Th 3 ! PFs are calculated from 175 energy levels of Blaise & Wyart 1992, Energy Levels and Atomic Spectra of Actinides, Paris
      DATA NNN_Pa/
     1 347877992,129318323,240730533,380546863,570368573,  600000,91,00,AEL Pa 1
     2 347877992,129318323,240730533,380546863,570368573, 1200000,91,01,FAK Pa 2
     3 347777992,129318323,240730533,380546863,570368573, 2000000,91,02/FAK Pa 3
      DATA NNN_U/
     1 209530092,450866762, 96613623,186524763,318839893,  619400,92,00,AEL U  1
     2  51311613,230239873,615986563,112513714,158317444, 1060000,92,01,Sne U  2 ! polynomial approximation from Batom.f subroutine of MOOG code: http://verdi.as.utexas.edu/moog.html
     3 211130612,456267402, 94912483,151817063,177417123, 2000000,92,02/Irw U  3 ! PFs are calculated using polynomial approximation from Irwin 1981, ApJS, 45, 621
      DATA NNN_Np/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,93,00,FAK Np 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,93,01,FAK Np 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,93,02/FAK Np 3
      DATA NNN_Pu/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,94,00,FAK Pu 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,94,01,FAK Pu 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,94,02/FAK Pu 3
      DATA NNN_Am/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,95,00,FAK Am 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,95,01,FAK Am 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,95,02/FAK Am 3
      DATA NNN_Cm/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,96,00,FAK Cm 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,96,01,FAK Cm 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,96,02/FAK Cm 3
      DATA NNN_Bk/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,97,00,FAK Bk 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,97,01,FAK Bk 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,97,02/FAK Bk 3
      DATA NNN_Cf/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,98,00,FAK Cf 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,98,01,FAK Cf 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,98,02/FAK Cf 3
      DATA NNN_Es/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,99,00,FAK Es 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,99,01,FAK Es 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,99,02/FAK Es 3
      DATA SCALE/0.001,0.01,0.1,1.0/
C     DATA FIRST/.TRUE./
C
C  First time XSAHA is called find the starting locations for each element
C
      FIRST = .TRUE.
      IF(FIRST) THEN
        FIRST=.FALSE.
        IZ=0
        DO 1 N=1,NTABLE
        IF(NNNPFN(7,N).NE.IZ.AND.IZ.LE.ELESIZ) THEN
          IZ=NNNPFN(7,N)
          LOCZ(IZ)=N
        ENDIF
   1    CONTINUE
        LOCZ(IZ+1)=NTABLE+1
      ENDIF
C
C  Find starting row in the partition table and the number of ionization
C  stages available for a given element IEL
C
      N=LOCZ(IEL)
      NIONS=LOCZ(IEL+1)-N
C
C  For MODE=5 return the number of ionizations available for IEL
C
      IF(MODE.EQ.5) THEN
        MAXION=NIONS
        RETURN
      ENDIF
C
C  Compute T and kT in eV
C
      TTKEV=8.6171E-5*TT
      TV=TTKEV
      TTK=1.38065E-16*TT
C
C  Lowering of the ionization potential in Volts for unit Zeff
C
      CHARGE=2.*XNELEC
      EXCESS=XNELEC-XNATOM
C
C  Special allowance for doubly ionized Helium
C
      IF(EXCESS.GT.0.) CHARGE=CHARGE-EXCESS+4.*(2.*EXCESS)
C
C  Original code:
C     DEBYE=SQRT(TTK/(2.8965E-18*CHARGE))
C     POTLOW=MIN(1.,1.44E-7/DEBYE)
C
C  Compute the inverse of Debye radius to avoid division by zero at low temperatures
C
      DEBYE=SQRT(2.8965E-18*CHARGE/TTK)
      POTLOW=MIN(1.,1.44E-7*DEBYE)
C
C  Solve the Saha equation
C
      NION2=NIONS
      N=N-1
      DO 2 IONN=1,NION2
      Z=IONN
      POTLO(IONN)=POTLOW*Z
C      write(*,*) IP(IONN)-POTLO(IONN)
      N=N+1
      NNN100=NNNPFN(6,N)/100
      IP(IONN)=FLOAT(NNN100)/1000.
      G=NNNPFN(6,N)-NNN100*100
      IF(N.EQ.1) THEN
        PART(1)=2.
c        IF(TT.LT.9000.) GO TO 2
        PART(1)=PART(1)+8.*EXP(-10.196/TV)+18.*EXP(-12.084/TV)+32.*
     *          EXP(-12.745/TV)+50.*EXP(-13.051/TV)+72.*EXP(-13.217/TV)
        D1=13.595/6.5/6.5/TV
        D2=POTLO(1)/TV
      ELSE
        T2000=IP(IONN)*2000./11.
        IT=MAX(1,MIN(9,INT(TT/T2000-.5)))
        DT=TT/T2000-FLOAT(IT)-.5
        PMIN=1.
        I=(IT+1)/2
        K1=NNNPFN(I,N)/100000
        K2=NNNPFN(I,N)-K1*100000
        K3=K2/10
        KSCALE=K2-K3*10
        IF(MOD(IT,2).EQ.0) THEN
          P1=K3*SCALE(KSCALE)
          K1=NNNPFN(I+1,N)/100000
          KSCALE=MOD(NNNPFN(I+1,N),10)
          P2=K1*SCALE(KSCALE)
        ELSE
          P1=K1*SCALE(KSCALE)
          P2=K3*SCALE(KSCALE)
          IF(DT.LT.0.AND.KSCALE.LE.1) KP1=P1
          IF(DT.LT.0.AND.KSCALE.LE.1.AND.KP1.EQ.INT(P2+.5)) PMIN=KP1
        END IF
        PART(IONN)=MAX(PMIN,P1+(P2-P1)*DT)
c        write(*,*) (NNNPFN(I,N),I=1,6),PART(IONN),IP(IONN),G,IONN
        IF(G.EQ.0.0.OR.POTLO(IONN).LT.0.1.OR.TT.LT.T2000*4.0) GO TO 2
        IF(TT.GT.(T2000*11.)) TV=(T2000*11.)*8.6171E-5
        D1=0.1/TV
      END IF
      D2=POTLO(IONN)/TV
      PART(IONN)=PART(IONN)+G*EXP(-IP(IONN)/TV)*
     *           (SQRT(13.595*Z*Z/TV/D2)**3*
     *           (1./3.+(1.-(.5+(1./18.+D2/120.)*D2)*D2)*D2)-
     -           SQRT(13.595*Z*Z/TV/D1)**3*
     *           (1./3.+(1.-(.5+(1./18.+D1/120.)*D1)*D1)*D1))
c      TV=TTKEV
   2  CONTINUE
C
      IF(MODE.NE.3) THEN
        CF=2.*2.4148D15*TT*SQRT(TT)/XNELEC
        FFF(1)=1.
        DO 3 IONN=2,NION2
C
C  IF is to avoid annoying floating point underflows
C
        FEXARG=(IP(IONN-1)-POTLO(IONN-1))/TV
c        write(*,*) IONN,NION2,PART(IONN)/PART(IONN-1),FEXARG
c        IF(FEXARG.GT.80.) THEN
c          FFF(IONN)=0.
c        ELSE
          FFF(IONN)=CF*PART(IONN)/PART(IONN-1)*EXP(-FEXARG)
c        END IF
   3    CONTINUE
        DO 4 IONN=NION2,2,-1
   4    FFF(1)=1.+FFF(IONN)*FFF(1)
        FFF(1)=1./FFF(1)
        DO 5 IONN=2,NION2
   5    FFF(IONN)=FFF(IONN-1)*FFF(IONN)
        DO 6 IONN=1,MAXION
   6    FRCT(IONN)=1.
      ELSE
        DO 7 IONN=1,MAXION
   7    FRCT(IONN)=0.
      END IF
C
C  Formulate the answer according to MODE
C
      NIONS=MIN(MAXION,NION2)
      IF(MODE.EQ.1) THEN
        FRCT(1)=FFF(1)/PART(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 8 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
   8      FRCT(IONN)=FFF(IONN)/PART(IONN)
        END IF
      ELSE IF(MODE.EQ.2) THEN
        FRCT(1)=FFF(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 9 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
   9      FRCT(IONN)=FFF(IONN)
        END IF
      ELSE IF(MODE.EQ.3) THEN
        FRCT(1)=PART(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 10 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
  10      FRCT(IONN)=PART(IONN)
        END IF
      ELSE IF(MODE.EQ.4) THEN
        FRCT(1)=0
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 11 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
  11      FRCT(1)=FRCT(1)+FFF(IONN)*(IONN-1)
        END IF
      END IF
C
      RETURN
      END
