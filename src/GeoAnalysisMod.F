!----------------------------------------------------------------------
!<author>  Kristof Karhan</author>
!<purpose> provide routines for the analysis of geometric properties</purpose>
!----------------------------------------------------------------------
module GeoAnalysisMod
use useful, only : Debug, PrintError
use geom_module, only : vec_len => vec_length, ImageDistVec
use AtomMod, only : Atom    
implicit none
public
character(*), parameter :: MName = "GeoAnalysisMod"

contains
!   -------------------------------------------------------------------
!   <param name="atom1"> the first atom </param>
!   <param name="atom2"> the second atom </param>
!   <param name="pbox"> optional parameter giving the cubic box length</param>
!   <returns> the image corrected distance<returns>
!   -------------------------------------------------------------------
    real(8) function atom_distance(atom1, atom2, pbox)
        implicit none
        class(Atom), intent(in)                     :: atom1, atom2
        real(8), Dimension(3)                       :: pos1, pos2
        real(8), Dimension(3), intent(in), optional :: pbox

        pos1 = atom1%Pos()
        pos2 = atom2%Pos()
        atom_distance = vec_len(ImageDistVec(pos1, pos2, pbox))
    end function 
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
