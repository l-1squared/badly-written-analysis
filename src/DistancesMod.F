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
use geom_module, only : ImageDistVec, vec_len => vec_length
implicit none
private
character(*), parameter :: MName = 'DistanceMod'
type                                                :: DataCol 
    private
        type(GenHist)                               :: CCdist
        type(GenHist)                               :: NNdist
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
        integer                                     :: iounit
        logical                                     :: ioisopen
        real(8)                                     :: histogram_interval
        namelist/distance_data/file_CC_distance, file_NN_distance, &
                               histogram_interval

        file_CC_distance = "CCdist.pf"
        file_NN_distance = "NNdist.pf"
        histogram_interval = 0.05
        call get_unit(finput, iounit, ioisopen)
        if (.not. ioisopen) open(iounit, FILE=finput, STATUS='old', &
                                 ACTION='read')
        read(iounit, NML=distance_data)
        if(debug) write(output_unit, NML=distance_data)
        if (.not.ioisopen) close(iounit, STATUS='keep')

        dac_cnstr1%CCdist = GenHist(file_CC_distance, histogram_interval)
        dac_cnstr1%NNdist = GenHist(file_NN_distance, histogram_interval)
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
        
        call this%CCdist%update(get_xi(atomlist), &
             vec_len(ImageDistVec(atomlist(17)%Pos(), atomlist(19)%Pos())))
        call this%CCdist%update(get_xi(atomlist), &
             vec_len(ImageDistVec(atomlist(1)%Pos(), atomlist(2)%Pos())))
    end subroutine


!   --------------------------------------------------------------------
!   <summary>
!       print the content of the grids
!   --------------------------------------------------------------------
    subroutine dac_print(this)
        implicit none
        character(*), parameter :: FName = 'dac_print'
        character(*), parameter :: CName = MName // ":" // FName
        class(DataCol), intent(in)                  :: this
        call this%CCdist%printHist()
        call this%NNdist%printHist()
    end subroutine


!   --------------------------------------------------------------------
!   <summary>
!       calculate the value of the reaction coordinate from the geometry
!   </summary>
!   <param name="atomllis"> the ordered list of atoms</param>
!   <returns> the image corrected distance </returns>
!   --------------------------------------------------------------------
    pure function get_xi(atomlist)result(distance)
        implicit none
        class(Atom), Dimension(:), intent(in)       :: atomlist
        real(8)                                     :: distance
        real(8), Dimension(3,4)                     :: pos
        real(8), Dimension(3)                       :: ctr1, ctr2
        pos(:,1) = atomlist(1)%Pos()
        pos(:,2) = atomlist(2)%Pos()
        pos(:,3) = atomlist(17)%Pos()
        pos(:,4) = atomlist(19)%Pos()
!       -------------------------------------
!       try to get the Centers of DMAD N atoms (ctr1)
!       and the Q C atoms (ctr2). ImageDistVec calculates i
!       pos(:,2)+outvec-pos(:,1). Hence, the halved distance
!       has to be ADDED to Pos(:,1).
        ctr1 = pos(:,1) + ImageDistVec(pos(:,1), pos(:,2)) * 0.5d0
        ctr2 = pos(:,3) + ImageDistVec(pos(:,3), pos(:,4)) * 0.5d0
        distance = vec_len(ImageDistVec(ctr2, ctr1))
    end function
end module
