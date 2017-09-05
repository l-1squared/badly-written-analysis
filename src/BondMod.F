module BondMod
use, intrinsic :: iso_fortran_env, only : output_unit
use useful, only : debug, PrintError, lrtrim
use Print_Arrays, only : get_unit
use MoleculeMod, only : Molecule
use HBNetSetup, only : get_acceptorindex => getAcceptorIndex
use GenGridMod, only : GenHist
use GeoAnalysisMod, only : get_xi
implicit none
character(*), parameter :: MName = 'BondMod'
private
type                                                :: BondAna
    private
        type(GenHist)                               :: dmadbonds
        type(GenHist)                               :: qbonds
    contains
        procedure, pass(this), public               :: update
        procedure, pass(this), public               :: write_otp
end type
interface BondAna
    module procedure bda_cnstr1
end interface
public :: BondAna
contains
!   -------------------------------------------------------------------
!   <summary>
!       read input file and initialize class members
!   </summary>
!   <param name="finput"> file containting the Bond_Ana namelist</param>
!   -------------------------------------------------------------------
    type(BondAna) function bda_cnstr1(finput)
        implicit none
        character(*), intent(in)                    :: finput
        character(Len=80)                           :: file_dmadbonds 
        character(LEN=80)                           :: file_qbonds
        real(8)                                     :: histogram_interval
        integer                                     :: iounit
        logical                                     :: ioisopen
        namelist/Bond_Ana/histogram_interval, file_dmadbonds, file_qbonds

        file_dmadbonds = "DMADBonds.pf"
        file_qbonds = "QBonds.pf"
        histogram_interval = 0.05d0
!       read input                
        call get_unit(finput, ioUnit, ioisopen)
        if(.not. ioisopen)open(iounit, FILE=finput, STATUS='old', ACTION='read')
        read(iounit, NML=Bond_Ana)
        if (debug) write(output_unit, NML=Bond_Ana)
        if(.not. ioisopen) close(iounit, STATUS='keep')
!       init members
        bda_cnstr1%dmadbonds = GenHist(file_dmadbonds, histogram_interval)
        bda_cnstr1%qbonds = GenHist(file_qbonds, histogram_interval)
    end function
!   -------------------------------------------------------------------
!   <summary> update members</summary>
!       
!   -------------------------------------------------------------------
    subroutine update(this, xi, dmad, q, water)
        implicit none
        character(*), parameter :: FName = 'update'
        character(*), parameter :: CName = MName // ":" // FName
        class(BondAna), intent(in)                  :: this
        real(8), intent(in)                         :: xi
        class(Molecule), intent(in)                 :: dmad
        class(Molecule), intent(in)                 :: q
        class(Molecule), Dimension(:), intent(in)   :: water
        call this%dmadbonds%update(xi, dble(count_bonds(dmad, water)))
        call this%qbonds%update(xi, dble(count_bonds(q, water)))
    end subroutine        


!   -------------------------------------------------------------------
!   <summary> write aggregated data to file </summary>
!   -------------------------------------------------------------------
    subroutine write_otp(this)
        implicit none
        class(BondAna), intent(in)                  :: this
        call this%dmadbonds%printHist()
        call this%qbonds%printHist()
    end subroutine


!   -------------------------------------------------------------------
!   <summary>
!       process the output of getacceptors
!   </summary>
!   <param 
!   -------------------------------------------------------------------
    integer function count_bonds(electron_donor, electron_acceptors)
        implicit none     
        character(*), parameter :: FName = 'count_bonds'
        character(*), parameter :: CName = MName // ":" // FName
        class(Molecule), intent(in)                 :: electron_donor
        class(Molecule), Dimension(:), intent(in)   :: electron_acceptors
        integer,Dimension(:), allocatable           :: connects
        integer                                     :: i

        count_bonds = 0
        do i=lbound(electron_acceptors, 1), ubound(electron_acceptors, 1), 1
            connects = get_acceptorindex(electron_acceptors(i), electron_donor)
        end do
        do i = lbound(connects, 1), ubound(connects, 1), 1 
            count_bonds =  count_bonds + merge(1, 0, connects(i) .gt. 0)
        end do
    end function
end module
