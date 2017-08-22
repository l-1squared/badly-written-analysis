!-----------------------------------------------------------------------
!<author> Kristof Karhan </author>
!<summary>  perform geometric analysis on QMMM / IW systems</summary>
!-----------------------------------------------------------------------
program AnalizeQMMM
use, intrinsic :: iso_fortran_env, only : output_unit
use useful, only : PrintError, debug, toChar
use ProvideGeoMod, only : InstanceGeo, TrajFile
implicit none
character(*), parameter :: MName = 'AnalyzeQMMM'
character(*), parameter :: PName = MName

logical fileopen
type(TrajFile)                                      :: qmmmtraj
type(InstanceGeo)                                   :: qmmmgeo
character(LEN=64)                                   :: fqmmmtraj = "NONE"
character(LEN=64)                                   :: fqmmmmol = "NONE"

integer                                             :: i

fileopen = .TRUE.
call PrintProgramData()
qmmmtraj =  TrajFile(fqmmmtraj)
i = 0
do while (fileopen)
    i = i + 1
    call qmmmtraj%read_next(qmmmgeo, fileopen)
    if (i .eq. 1)then
        call qmmmgeo%assign_mol(fqmmmmol)
    end  if     
end do
contains
#include "Main_subroutines.F"
end program
