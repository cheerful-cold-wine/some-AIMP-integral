! AIMP(ab initio model potential)
module class_aimp
  use, intrinsic :: ISO_C_BINDING, only: C_INT, C_LONG, C_LONG_LONG, &
                                         C_FLOAT, C_DOUBLE, C_LOC
  implicit none
  ! interfaces
  public :: init_aimp, print_aimp, read_aimp

  integer, parameter :: maxnumM1=16, maxnumM2=1 
  integer, parameter :: maxL=3        ! max L of orbital in PROJ                  
  integer, parameter :: maxprim=22    ! max number of primitive GTO in one shell 
  integer, parameter :: maxcgto=5     ! max number of contracted GTO in one shell

  type AbInitioModelPotential
    ! Coulomb Model Potential
    integer(C_LONG) :: nM1, nM2        ! M1: cGTO/r, nuclear attraction times s-type gaussian
    real(C_DOUBLE)  :: M1exp(maxnumM1) ! M2: cGTO  , s-type gaussian
    real(C_DOUBLE)  :: M1coef(maxnumM1)
    real(C_DOUBLE)  :: M2exp(maxnumM2)
    real(C_DOUBLE)  :: M2coef(maxnumM2)

    ! Projection Operator
    integer(C_LONG) :: PROJmaxl                         ! max angular momentum
    integer(C_LONG) :: PROJnumPrim(maxL+1)              ! num of primitive/contracted GTO in a shell
    integer(C_LONG) :: PROJnumCgto(maxL+1)              ! primitive GTO -> contracted GTO
    real(C_DOUBLE)  :: PROJconstant(maxcgto,maxL+1)     ! constant=-2E_c, E_c is the energy of cGTO
    real(C_DOUBLE)  :: PROJexp(maxprim,maxL+1)          ! all exponents of primitive GTO
    real(C_DOUBLE)  :: PROJcoef(maxcgto,maxprim,maxL+1) ! all coefficients of pGTO -> cGTO

    ! Spectral Representation Operator
    ! for exchange operator(always) and relativistic effect(optional)
    character*35    :: SROtype ! orbital type: valence/core/mixed-valence-core/external primitive basis
    integer(C_LONG) :: SROrela ! 0 not relativistic, 1 Douglas-Kroll AIMP-NP, 2 Cowan-Griffin AIMP-CG
    real(C_DOUBLE) :: corerep
  end type

  contains
    subroutine init_aimp(aimpdata)
      implicit none
      type(AbInitioModelPotential), intent(out) :: aimpdata
      aimpdata%nM1=0
      aimpdata%nM2=0
      aimpdata%M1exp=0.d0
      aimpdata%M1coef=0.d0
      aimpdata%M2exp=0.d0
      aimpdata%M2coef=0.d0

      aimpdata%PROJmaxl=0
      aimpdata%PROJnumPrim=0
      aimpdata%PROJnumCgto=0
      aimpdata%PROJconstant=0.d0
      aimpdata%PROJexp=0.d0
      aimpdata%PROJcoef=0.d0

      aimpdata%SROtype=" "
      aimpdata%SROrela=0

      aimpdata%corerep=0.d0
    end subroutine init_aimp

    subroutine print_aimp(aimpdata)
      implicit none
      type(AbInitioModelPotential), intent(in) :: aimpdata
      integer :: i
      write(6,*) "# Ab Initio Model Potential part"
      write(6,*) "M1"
      write(6,*) aimpdata%nM1
      if (aimpdata%nM1.ne.0) then
        write(6,*) aimpdata%M1exp(1:aimpdata%nM1)
        write(6,*) aimpdata%M1coef(1:aimpdata%nM1)
      endif
      write(6,*) "M2"
      write(6,*) aimpdata%nM2
      if (aimpdata%nM2.ne.0) then
        write(6,*) aimpdata%M1exp(1:aimpdata%nM2)
        write(6,*) aimpdata%M1coef(1:aimpdata%nM2)
      endif
      write(6,*) "COREREP"
      write(6,*) aimpdata%corerep
      write(6,*) "PROJOP"
      write(6,*) aimpdata%PROJmaxl
      do i=1 , aimpdata%PROJmaxl+1
        write(6,*) aimpdata%PROJnumPrim(i), aimpdata%PROJnumCgto(i)
        write(6,*) aimpdata%PROJconstant(1:aimpdata%PROJnumCgto(i),i)
        write(6,*) aimpdata%PROJexp(1:aimpdata%PROJnumPrim(i),i)
        write(6,*) aimpdata%PROJcoef(1:aimpdata%PROJnumCgto(i),1:aimpdata%PROJnumPrim(i),i)
      enddo
      write(6,*) "Spectral Representation Operator"
      write(6,*) aimpdata%SROtype
      write(6,*) "Exchange"
      select case(aimpdata%SROrela)
      case(0)
      case(1)
        write(6,*) "NoPair"
      case(2)
        write(6,*) "1stOrder Relativistic Correction"
      endselect
      write(6,*) "End of Spectral Representation Operator"
    end subroutine print_aimp

    subroutine read_aimp(fn,aimpdata)
      implicit none
      type(AbInitioModelPotential), intent(out) :: aimpdata
      integer, intent(in) :: fn
      !-----------------------------------------------------
      character*80 :: string,SROstr
      logical :: logic_endSRO

      logic_endSRO=.false.
      do while(.not.logic_endSRO)
        read(fn,"(a)") string
        select case(string)
        case("M1")
          read(fn,*) aimpdata%nM1
          if(aimpdata%nM1.ne.0) then
            read(fn,*) aimpdata%M1exp(1:aimpdata%nM1)
            read(fn,*) aimpdata%M1coef(1:aimpdata%nM1)
          endif
        case("M2")
          read(fn,*) aimpdata%nM2
          if(aimpdata%nM2.ne.0) then
            read(fn,*) aimpdata%M2exp(1:aimpdata%nM2)
            read(fn,*) aimpdata%M2coef(1:aimpdata%nM2)
          endif
        case("COREREP")
          read(fn,*) aimpdata%corerep
        case("PROJOP")
          call read_aimpPROJ_cGTO(fn,aimpdata)
        case("Spectral Representation Operator")
          do while(.not. logic_endSRO)
            read(fn,"(a)") SROstr
            if(SROstr.eq."Valence primitive basis"            .or.&
               SROstr.eq."Core primitive basis"               .or.&
               SROstr.eq."Mixed valence-core primitive basis" .or.&
               SROstr.eq."External primitive basis")  aimpdata%SROtype=SROstr
            if(SROstr.eq."Exchange") cycle
            if(SROstr.eq."NoPair") aimpdata%SROrela=1                           ! use AIMP-NP
            if(SROstr.eq."1stOrder Relativistic Correction") aimpdata%SROrela=2 ! use AIMP-CG
            if(SROstr.eq."End of Spectral Representation Operator") logic_endSRO=.true.
          end do
        case default
          logic_endSRO=.false.
        end select
      end do
    end subroutine read_aimp

    subroutine read_aimpPROJ_cGTO(fn,aimpdata)
      implicit none
      type(AbInitioModelPotential), intent(out) :: aimpdata
      integer, intent(in) :: fn
      !-----------------------------------------------------
      integer :: nshell, i, nPrim, nCgto

      read(fn,*) aimpdata%PROJmaxl
      nshell = aimpdata%PROJmaxl+1
      do i=1, nshell
        read(fn,*) nPrim, nCgto
        read(fn,*) aimpdata%PROJconstant(1:nCgto,i)
        read(fn,*) aimpdata%PROJexp(1:nPrim,i)
        read(fn,*) aimpdata%PROJcoef(1:nCgto,1:nPrim,i)
        aimpdata%PROJnumPrim(i)=nPrim
        aimpdata%PROJnumCgto(i)=nCgto
      end do
    end subroutine read_aimpPROJ_cGTO

end module class_aimp
