!----------------------------------------------------------------------
!<author> Kristof Karhan</author>
!<purpose> 
!   provide instanciable storage for the contents of a geometry 
!   file
!</purpose>
!----------------------------------------------------------------------
module ProvideGeoMod
use useful, only : PrintError, debug
use Print_Arrays, only : get_unit
use AtomMod, only : Atom
use MoleculeMod, only : Molecule
use ReadGeoMod, only : Snapshot, foreign_read => read_file
use DetMolsMod, only : simple_mol_list => simpleMolList
implicit none
private
character(*), parameter :: MName = 'ProvideGeoMod'
type,extends(Snapshot)                              :: InstanceGeo
    real(8), Dimension(:,:), allocatable            :: geo
    character(LEN=2), Dimension(:), allocatable     :: elements
    type(Molecule), Dimension(:), allocatable       :: molecules
    logical                                         :: molassigned = .FALSE.
    contains 
        procedure, pass(this), public               :: assign_mol
end type
type                                                :: TrajFile
    private
        character(LEN=64)                           :: ftrajectory
        integer                                     :: iounit
        logical                                     :: openedbyme
        logical                                     :: isopen
    contains
        procedure, pass(this), public               :: read_next
        procedure, pass(this), public               :: close_file
end type


interface TrajFile
    module procedure trf_cnstr1
end interface
interface InstanceGeo
    module procedure ing_cnstr1
end interface
public :: InstanceGeo, TrajFile

contains
!   ------------------------------------------------------------------
!   constructors
!   ------------------------------------------------------------------
!   <summary> initialize the geo instance </summary>    
!   <param name="natoms"> the number of atoms in a snapshot</param>
    function ing_cnstr1(natoms)
        implicit none
        character(*), parameter :: FName = 'ing_cnstr1'
        character(*), parameter :: CName = MName // ":" // FName
        integer, intent(in)                         :: natoms
        type(InstanceGeo)                           :: ing_cnstr1
        ing_cnstr1%SNAPSHOT = Snapshot(natoms)
        ing_cnstr1%molassigned = .FALSE.
    end function


!   <summary> initialize the TrajFile instance </summary>    
!   <param name="ftrajectory"> the trajectory file name </param>
    function trf_cnstr1(ftrajectory)
        implicit none
        character(*), parameter :: FName = 'trf_cnstr1'
        character(*), parameter :: CName = MName // ":" // FName
        character(LEN=*), intent(in)                :: ftrajectory
        type(TrajFile)                              :: trf_cnstr1
        integer                                     :: iounit,iostat
        logical                                     :: ioisopen
        call get_unit(ftrajectory, iounit, ioisopen)
        trf_cnstr1%openedbyme = .not. ioisopen
        trf_cnstr1%ftrajectory = ftrajectory
        if(.not. ioisopen) open(iounit, FILE=ftrajectory, STATUS='old',&
                                ACTION='read', IOSTAT=iostat)
        trf_cnstr1%isopen = iostat .eq. 0
        trf_cnstr1%iounit = iounit
    end function


!   ------------------------------------------------------------------
!   <summary>
!       provide the next snapshot in a trajectory file
!   </summary>
!   <param name="geoinstance">the object to write the snapshot data to</param>
!   <param name="report whether the file has ended</param>
!   ------------------------------------------------------------------
    subroutine read_next(this, geoinstance, file_ended)
        implicit none
        character(*), parameter :: FName = 'read_file'
        character(*), parameter :: CName = MName // ":" // FName
        type(InstanceGeo), intent(inout)            :: geoinstance
        class(TrajFile), intent(in)                 :: this
        logical, intent(out)                        :: file_ended
        integer                                     :: iostat
        character(LEN=90)                           :: line

        if (.not. this%isopen) then
            call PrintError(3, CName, "Attempt to read from closed file.")
        end if
        call foreign_read(this%ftrajectory, geoinstance%geo, &
                          geoinstance%elements, geoinstance%SNAPSHOT)
!       test wheter the file has ended
        read(this%iounit, *, IOSTAT=iostat) line
        file_ended = iostat .lt. 0
        backspace(this%iounit, IOSTAT=iostat)
    end subroutine
!   ------------------------------------------------------------------
!   <summary> close the trajectory  file </summary> 
!   <param name="isclosed> report whether operation was successful</param>
!   ------------------------------------------------------------------
    subroutine close_file(this, isclosed)
        implicit none
        class(TrajFile), intent(inout)              :: this
        logical, intent(out), optional              :: isclosed
        integer                                     :: iostat
        if (this%isopen .AND. this%openedbyme) then
            close(this%iounit, STATUS='KEEP', IOSTAT=iostat)
            this%isopen = .FALSE.
            this%openedbyme = .FALSE.
            if (present(isclosed)) isclosed = iostat .eq. 0
        else
            if (present(isclosed)) isclosed = .FALSE.
        end  if
    end subroutine
    

!   ------------------------------------------------------------------
!   <summary>
!       use the molecule id to assign the atoms of GeoInstance to 
!       its molecule list
!   </summary>
!   <param name="fmol>list of the molecule id of each atom</param>
!   ------------------------------------------------------------------
    subroutine assign_mol(this, fmol)
        implicit none
        character(*), parameter :: FName = 'assign_mol'
        character(*), parameter :: CName = MName // ":" // FName
        class(InstanceGeo), intent(inout)           :: this
        character(LEN=*), intent(in)                :: fmol
        call simple_mol_list(fmol, this%atomCards, this%molecules)
    end subroutine
end module
