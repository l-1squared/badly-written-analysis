!-----------------------------------------------------------------------
!<author> Kristof Karhan </author>
!<summary>  perform geometric analysis on QMMM / IW systems</summary>
!-----------------------------------------------------------------------
#ifndef __DEBUG
#   define __DEBUG 0
#endif
program AnalyzeQMMM
use, intrinsic :: iso_fortran_env, only : output_unit
use useful, only : PrintError, debug, toChar, globlinelen => GlobLineLen, toChar
use useful, only : backspace_carriage,  lrtrim
use Print_Arrays, only : get_unit
use MoleculeMod, only : Molecule
use ProvideGeoMod, only : InstanceGeo, TrajFile
use GeoAnalysisMod, only : get_xi
use DistancesMod, only : DataCol
use BondMod, only : BondAna
implicit none
character(*), parameter :: MName = 'AnalyzeQMMM'
character(*), parameter :: PName = MName
character(*), parameter :: CName = MName // ":" // "__main__"

character(LEN=80)                                   :: finput

logical                                             :: file_ended
character(LEN=80), Dimension(:), allocatable        :: qmmmtrajlist
character(LEN=80)                                   :: fqmmmmol, fqmmmtrajlist
character(LEN=80)                                   :: fqmmmtraj

type(TrajFile)                                      :: qmmmtraj
type(InstanceGeo), target                           :: qmmmgeo
type(DataCol)                                       :: qmmmdistances
type(BondAna)                                       :: qmmmbonds

integer                                             :: ireplica, isnapshot
integer                                             :: nreplica
character(LEN=:), allocatable                       :: fname
class(Molecule), pointer                            :: dmad, q
class(Molecule), Dimension(:), pointer              :: rest 

#if __DEBUG == 1
debug  = .TRUE.
#else
debug  = .FALSE.
#endif
call PrintProgramData()
call get_command_argument(1, finput)

qmmmdistances = DataCol(finput)
qmmmbonds = BondAna(finput)

call read_input(finput, fqmmmtrajlist, nreplica, fqmmmmol)
call retrieve_lists(fqmmmtrajlist, qmmmtrajlist, nreplica)

do ireplica =1, nreplica, 1
    fqmmmtraj = qmmmtrajlist(ireplica)
    qmmmtraj = TrajFile(fqmmmtraj)
    file_ended = .FALSE.
    isnapshot = 0
    do while (.not. file_ended)
        isnapshot = isnapshot + 1
        call qmmmtraj%read_next(qmmmgeo, file_ended)
        fname = lrtrim(fqmmmtraj(index(fqmmmtraj, "/", .TRUE.) + 1 : &
                       len(fqmmmtraj)))
        write(output_unit, FMT='(1x,A23,1x,A' // toChar(len(fname)) // ',1x,' &
              //'A8,1x,I4)', ADVANCE='no') "Reading trajectory file", fname, &
              "snapshot", isnapshot
        call backspace_carriage(output_unit, 25 + len(fname) + 10 + 4)
        if (isnapshot .eq. 1)then
            call qmmmgeo%assign_mol(fqmmmmol)
            dmad => qmmmgeo%molecules(1)
            q => qmmmgeo%molecules(2)
            rest => qmmmgeo%molecules(3:size(qmmmgeo%molecules))
        end if     
        call qmmmdistances%update(qmmmgeo%atomCards)
        call qmmmbonds%update(get_xi(qmmmgeo%atomCards), dmad, q, rest)
    end do
    isnapshot = isnapshot + 1
    write(output_unit,*)
    call qmmmtraj%close_file()
end do
call qmmmdistances%write_otp()
call qmmmbonds%write_otp()

contains


!   --------------------------------------------------------------------
!   <summary>
!       read the file that lists the trajectory files
!   </summary>
!   <param name="ftrajlist"> the list file</param>
!   <param name="listtrajfiles"> list trajectory files go into </param>
!   <param name="nfiles"> number of files in ftrajlist</param>
!   --------------------------------------------------------------------
    subroutine retrieve_lists(ftrajlist, listtrajfiles, nfiles)
        implicit none
        character(*), parameter :: FName = 'retrieve_lists'
        character(*), parameter :: CName = MName // ":" // FName
        character(LEN=*), intent(in)                :: ftrajlist
        character(LEN=*), Dimension(:), allocatable, intent(out):: listtrajfiles
        integer                                     :: iounit, iostat
        logical                                     :: ioisopen
        integer                                     :: i, nfiles, ifiles

        ifiles=0
        call get_unit(ftrajlist, iounit, ioisopen)
        allocate(listtrajfiles(nfiles))
        if (.not. ioisopen) open(iounit, FILE=ftrajlist, STATUS='old', &
                                 ACTION='read')
        do i = 1, nfiles     
            read(iounit, FMT='(A'//toChar(len(listtrajfiles(1)))//')', &
                 IOSTAT=iostat) listtrajfiles(i)
            if (iostat .eq. 0) then
                ifiles = ifiles + 1
            else if (iostat .gt. 0) then
                call PrintError(3, CName, "Error reading file '" // ftrajlist &
                                // "' in line "// toChar(ifiles) // ".")
            end if
        end do
        if (nfiles .ne. ifiles) then
            call PrintError(3, CName, "number of files read unequal to " &
                            // "number of files expected. (" // toChar( &
                            nfiles) // " / " // toChar(ifiles) // ")")
        end if
        if(.not. ioisopen) close(iounit, STATUS='keep')
    end subroutine        


#   include "Main_subroutines.F"
end program
