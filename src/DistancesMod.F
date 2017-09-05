!----------------------------------------------------------------------
!<author> Kristof Karhan </author>
!<purpose> 
!   collect properties that are based off the distance of certain atoms
!</purpose>
!----------------------------------------------------------------------
module DistancesMod
use iso_fortran_env, only : output_unit
use useful, only : globlinelen => GlobLineLen, debug
use Print_Arrays, only : get_unit
use GenGridMod, only : GenHist
use AtomMod, only : Atom
use GeoAnalysisMod, only : atom_distance, get_xi

implicit none
private
character(*), parameter :: MName = 'DistanceMod'
type                                                :: DataCol 
    private
        type(GenHist)                               :: CCdist
        type(GenHist)                               :: NNdist
        type(GenHist)                               :: CN1dist
        type(GenHist)                               :: CN2dist
    contains
        procedure, pass(this), public               :: update => dac_update
        procedure, pass(this), public               :: write_otp =>dac_print
end type

interface DataCol
    module procedure dac_cnstr1
end interface

public :: DataCol

contains


!   --------------------------------------------------------------------
!   <summary>
!      read the data from the input file and assign it appropriately
!   </summary>
!   <param name="finput" the name of the input file</param>
!   <returns> a new DataCol object</param>
!   --------------------------------------------------------------------
    function dac_cnstr1(finput)
        implicit none
        character(*), parameter :: FName = 'dac_cnstr1'
        character(*), parameter :: CName = MName // ":" // FName
        type(DataCol)                               :: dac_cnstr1
        character(LEN=*), intent(in)                :: finput
        character(LEN=globlinelen)                  :: file_CC_distance
        character(LEN=globlinelen)                  :: file_NN_distance
        character(LEN=globlinelen)                  :: file_CN1_distance
        character(LEN=globlinelen)                  :: file_CN2_distance
        integer                                     :: iounit
        logical                                     :: ioisopen
        real(8)                                     :: histogram_interval
        namelist/distance_data/file_CC_distance, file_NN_distance, &
                               file_CN1_distance, file_CN2_distance, &
                               histogram_interval

        file_CC_distance = "CCdist.pf"
        file_NN_distance = "NNdist.pf"
        file_CN1_distance = "CN1dist.pf"
        file_CN2_distance = "CN2dist.pf"
        histogram_interval = 0.05
        call get_unit(finput, iounit, ioisopen)
        if (.not. ioisopen) open(iounit, FILE=finput, STATUS='old', &
                                 ACTION='read')
        read(iounit, NML=distance_data)
        if(debug) write(output_unit, NML=distance_data)
        if (.not.ioisopen) close(iounit, STATUS='keep')

        dac_cnstr1%CCdist = GenHist(file_CC_distance, histogram_interval)
        dac_cnstr1%NNdist = GenHist(file_NN_distance, histogram_interval)
        dac_cnstr1%CN1dist = GenHist(file_CN1_distance, histogram_interval)
        dac_cnstr1%CN2dist = GenHist(file_CN2_distance, histogram_interval)
    end function



!   --------------------------------------------------------------------
!   <summary>
!       update stuff
!   </summary>
!   <param name="atomlist">
!   --------------------------------------------------------------------
    subroutine dac_update(this, atomlist)
        implicit none
        character(*), parameter :: FName = 'dac_update'
        character(*), parameter :: CName = MName // ":" // FName
        class(DataCol), intent(inout)               :: this
        class(Atom), Dimension(:), intent(in)       :: atomlist
        real(8)                                     :: distance1, distance2
        
        call this%CCdist%update(get_xi(atomlist), &
             atom_distance(atomlist(17), atomlist(19)))
        call this%NNdist%update(get_xi(atomlist), &
             atom_distance(atomlist(1), atomlist(2)))
!       Q can rotate with respect to DMAD. Therefore, first check witch
!       C is closer
        distance1 = atom_distance(atomlist(1), atomlist(17))
        distance2 = atom_distance(atomlist(1), atomlist(19))
        if (distance1 .lt. distance2) then
            call this%CN1dist%update(get_xi(atomlist), &
                 atom_distance(atomlist(1), atomlist(17)))
            call this%CN2dist%update(get_xi(atomlist), &
                 atom_distance(atomlist(2), atomlist(19)))
        else
            call this%CN1dist%update(get_xi(atomlist), &
                 atom_distance(atomlist(1), atomlist(19)))
            call this%CN2dist%update(get_xi(atomlist), &
                 atom_distance(atomlist(2), atomlist(17)))
        end if            
    end subroutine


!   --------------------------------------------------------------------
!   <summary>
!       print the content of the grids
!   </summary>
!   --------------------------------------------------------------------
    subroutine dac_print(this)
        implicit none
        character(*), parameter :: FName = 'dac_print'
        character(*), parameter :: CName = MName // ":" // FName
        class(DataCol), intent(in)                  :: this
        call this%CCdist%printHist()
        call this%NNdist%printHist()
        call this%CN1dist%printHist()
        call this%CN2dist%printHist()
    end subroutine
end module

