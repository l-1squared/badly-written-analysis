!-----------------------------------------------------------------------
!<author> Kristof Karhan </author>
!<summary>  perform geometric analysis on QMMM / IW systems</summary>
!-----------------------------------------------------------------------
program AnalizeQMMM
use, intrinsic :: iso_fortran_env, only : output_unit
use useful, only : PrintError, debug, toChar, globlinelen => GlobLineLen
use Print_Arrays, only : get_unit
use ProvideGeoMod, only : InstanceGeo, TrajFile
use DistancesMod, only : DataCol
implicit none
character(*), parameter :: MName = 'AnalyzeQMMM'
character(*), parameter :: PName = MName

character(LEN=64)                                   :: finput

logical                                             :: fileopen
character(LEN=64), Dimension(:), allocatable        :: qmmmtrajlist
character(LEN=64)                                   :: fqmmmmol, fqmmmtrajlist
character(LEN=64)                                   :: fqmmmtraj

type(TrajFile)                                      :: qmmmtraj
type(InstanceGeo)                                   :: qmmmgeo
type(DataCol)                                       :: qmmmdistances

integer                                             :: ireplica, isnapshot
integer                                             :: nreplica

call PrintProgramData()
call get_command_argument(1, finput)

qmmmtraj = TrajFile(fqmmmtraj)
qmmmdistances = DataCol(finput)

call read_input(finput, fqmmmtrajlist, nreplica)
call retrieve_lists(fqmmmtrajlist, qmmmtrajlist, nreplica)
do while (ireplica .le. nreplica)
    ireplica = ireplica + 1
    fqmmmtraj = qmmmtrajlist(ireplica)
    qmmmtraj = TrajFile(fqmmmtraj)
    fileopen = .TRUE.
    do while (fileopen)
        isnapshot = isnapshot + 1
        call qmmmtraj%read_next(qmmmgeo, fileopen)
        if (isnapshot .eq. 1)then
            call qmmmgeo%assign_mol(fqmmmmol)
        end if     
        call qmmmdistances%update(qmmmgeo%atomCards)
    end do
end do
call qmmmdistances%write_otp()

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
        integer                                     :: iounit
        logical                                     :: ioisopen
        integer                                     :: i, nfiles
        call get_unit(ftrajlist, iounit, ioisopen)
        allocate(listtrajfiles(nfiles))
        if (.not. ioisopen) open(iounit, FILE=ftrajlist, STATUS='old', &
                                 ACTION='read')
        do i = 1, nfiles     
            read(iounit, *)listtrajfiles(i)
        end do
    end subroutine        


#   include "Main_subroutines.F"
end program
