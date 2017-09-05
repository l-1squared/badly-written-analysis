module HBNetSetup
use useful, only : PrintError, Nullsetzen, toChar, lrtrim
use Print_Arrays, only : Printxyz,get_unit
use geom_module, only : ImageDistVec,vec_length,Scalar_Product,pi
use geom_module, only : Vector_Product
use AtomMod, only : Atom
use MoleculeMod, only : Molecule
implicit none
private
character(*),parameter :: MName="HBNetSetup"
type                                              :: GeoConnect
    class(Molecule),pointer                       :: Partner
end type

type                                              :: XHBond
    class(Atom), pointer                          :: h
    class(Atom), pointer                          :: anchor
    integer                                       :: hmolid
    integer                                       :: anchorid
end type


public :: getAcceptors, getAcceptorIndex,GeoConnect 
public :: WaterArray
contains
!   -------------------------------------------------------------------
!   <summary>
!      return a list of all hydrogen bonds with their anchor to be 
!      processed by getAcceptorIndex
!   </summary>
!   <param name="amol"> the respective acceptor molecule </param>
!   <param name="anchoredh"> the list of H with their anchor </param>
!   -------------------------------------------------------------------
    subroutine get_h(amol, anchoredh)
        implicit none
        character(*), parameter :: FName = 'get_h'
        character(*), parameter :: CName = MName // ":" // FName
        class(Molecule),target                      :: amol
        type(XHBond), Dimension(:), allocatable, intent(out) :: anchoredh
        real(8), Dimension(:,:), allocatable                 :: dists
        integer, Dimension(:,:), allocatable                 :: anchorid
        integer, Dimension(:), allocatable                   :: hid
        integer                                              :: nh, nanchor
        integer                                              :: ih, ianchor
        integer                                              :: i, j
        real(8)                                              :: mindist
        integer                                              :: minlocation
        
!       find the number of hydrogens in molecule
        nh = 0
        nanchor = 0
        do i = 1, amol%msize(), 1
            associate(atm => amol%getAtom(i))
                if (lrtrim(atm%element()) == "H") then
                    nh = nh + 1 
                else
                    nanchor = nanchor + 1
                end if
            end associate
        end do
        allocate(anchoredh(1:nh), hid(1:nh))
        allocate(dists(amol%Msize(), amol%Msize), anchorid(1:nh, 1:nanchor))
!       find all distances between h and anchors
        ih = 0
        ianchor = 0
        dist(:) = 1d9
        do i = 1, amol%msize(), 1
            do j = i + 1, amol%msize(), 1
                dist(i,j) = atom_distance(amol%getAtom(i), amol%getAtom(j))
            end do
            associate(atm => amol%getAtom())
                if (lrtrim(atm%element) == "H") then
                    ih = ih  + 1
                    hid(ih) = i
                else
                    ianchor = ianchor + 1
                    anchorid(ianchor) = i
            end associate
        end do            
!       TODO create anchoredh array
        do i = 1, nh, 1
            do j = nanchor, 1
                if (dist(hid(i), anchorid(j)) .lt. mindist) then                
                    mindist = dist(hid(i), anchorid(j))
                    minlocation = anchorid(j)
                end if
            end do
            anchoredh(i)%h => amol%getAtom(hid(i))
            anchoredh(i)%anchor => amol%getAtom(minlocation)
            anchoredh(i)%hmolid => hid(i)
            anchoredh(i)%anchorid => minlocation 
        end do
    end subroutine


!   TODO rewrite function to be useable with anchoredh
    function get_acceptor_index(Acceptor,Donor)result(ConnectedTo)
        implicit none
        character(*),parameter :: FName="getAcceptorIndex"
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule),intent(in)                        :: Acceptor
        class(Molecule),intent(in)                        :: Donor
        type(Atom),pointer                                :: CorePnt
        real(8),Dimension(3,3)                            :: WaterPos
        real(8),Dimension(3)                              :: O
        real(8),Dimension(3)                              :: OH1,OH2,ONu
        integer,Dimension(:),allocatable                  :: ConnectedTo
        real(8)                                           :: Angle
        real(8)                                           :: Distance
        integer                                           :: i
        nullify(CorePnt)
        WaterPos = WaterArray(Acceptor)!returns Water Coords in OHH format
        O = WaterPos(1,:)
        OH1 = ImageDistVec(WaterPos(1,:),WaterPos(2,:))
        OH2 = ImageDistVec(WaterPos(1,:),WaterPos(3,:))
        !-------------------------------------------------
        !Now search the acceptor molecules for possible 
        !acceptors (Currently C,N,O,F),and check if there 
        !is an hbond between them
        allocate(ConnectedTo(Donor%Msize()))
        call Nullsetzen(ConnectedTo)
        do i=1,Donor%Msize(),1
            CorePnt=>Donor%Atoms(i)
            if(scan(CorePnt%Element(),"CONF").gt.0)then
                ONu=ImageDistVec(O,CorePnt%Pos())
                 do j = 1, size(OH, 1), 1
                    Distance = vec_length(ONu(j))
                    Angle = ScalarProduct(OH(j),ONu(j))
                    if(is_Bonded(Distance,Angle))then
                        ConnectedTo(i)=ConnectedTo(i)+j !Indicates j-th H is bonded
                    end if
                 end do
            end if
        end do
    end function
!   -------------------------------------------------------------------
!   Searches through 'Molecules'to find  molecules that accept Hydrogen
!   Molecules
!   -------------------------------------------------------------------
    subroutine GetAcceptors(IMolecule,Molecules,ConnectedTo)
    implicit none
    character(*),parameter :: FName="GetAcceptors"
    character(*),parameter :: CName=MName//":"//FName
    class(Molecule),Dimension(:),intent(in),target    :: Molecules
    integer,intent(in)                                :: IMolecule
    class(GeoConnect),Dimension(:),intent(out),allocatable :: ConnectedTo
    real(8),Dimension(10)                             :: TempArray
    integer,Dimension(:),allocatable                  :: NBonds
    integer                                           :: IListIndex
    integer,target                                    :: i

    IListIndex=0
    call Nullsetzen(TempArray)
    do i=lbound(Molecules,1),ubound(Molecules,1),1
        !if(vec_length(ImageDistVec(Molecules(IMolecule)%COM(),Molecules(i)%COM())).gt.5)cycle
        NBonds=getAcceptorIndex(Molecules(IMolecule),Molecules(i))
        if(maxVal(NBonds).ge.3)then
            call PrintError(3,CName,"Too many H-Bonds between molecules "&
            //tochar(IMolecule)//" and "//toChar(i)//" ("//toChar(NBonds(1))&
            //").")
            IListIndex=IListIndex+1
            TempArray(IListIndex)=i
        else if(maxVal(NBonds).ge.1)then
            IListIndex=IListIndex+1
            TempArray(IListIndex)=i
        end if
    end do
    if(IListIndex.gt.ubound(TempArray,1))&
    call PrintError(3,CName,"Index for TempArray is out of bounds.")
    allocate(ConnectedTo(1:IListIndex))
    do i=1,IListIndex
        ConnectedTo(i)%Partner=>Molecules(int(TempArray(i)))
    end do
    end subroutine
!-----------------------------------------------------------------------
!Wrapping function for is_Bonded, takes the Molecule classes Donor and 
!Acceptor, where Acceptor is defined as the molecule to which the
!Hydrogen is bound. 
!The return value 'ConnectedTo' is an Array the size of the H-Atom
!Accepting molecule. For each potential Acceptor an integer is created. 
!If an H Atom Binds to the Acceptor the interger is assigned '1' if the 
!first H-Atom Binds and '2' if the second one does.
!-----------------------------------------------------------------------
     function getAcceptorIndex(Acceptor,Donor)result(ConnectedTo)
    implicit none
    character(*),parameter :: FName="getAcceptorIndex"
    character(*),parameter :: CName=MName//":"//FName
    class(Molecule),intent(in)                        :: Acceptor
    class(Molecule),intent(in)                        :: Donor
    type(Atom),pointer                                :: CorePnt
    real(8),Dimension(3,3)                            :: WaterPos
    real(8),Dimension(3)                              :: O
    real(8),Dimension(3)                              :: OH1,OH2,ONu
    integer,Dimension(:),allocatable                  :: ConnectedTo
    real(8)                                           :: Angle
    real(8)                                           :: Distance
    integer                                           :: i
    nullify(CorePnt)
    WaterPos = WaterArray(Acceptor)!returns Water Coords in OHH format
    O = WaterPos(1,:)
    OH1 = ImageDistVec(WaterPos(1,:),WaterPos(2,:))
    OH2 = ImageDistVec(WaterPos(1,:),WaterPos(3,:))
    !-------------------------------------------------
    !Now search the acceptor molecules for possible 
    !acceptors (Currently C,N,O,F),and check if there 
    !is an hbond between them
    allocate(ConnectedTo(Donor%Msize()))
    call Nullsetzen(ConnectedTo)
    do i=1,Donor%Msize(),1
        CorePnt=>Donor%Atoms(i)
        if(scan(CorePnt%Element(),"CONF").gt.0)then
            ONu=ImageDistVec(O,CorePnt%Pos())
            Distance=vec_length(ONu)
            Angle=Scalar_Product(OH1,ONu)
            if(is_Bonded(Distance,Angle))then
                ConnectedTo(i)=ConnectedTo(i)+1 !Indicates first H is bonded
            end if
            Angle=Scalar_Product(OH2,ONu)
            if(is_Bonded(Distance,Angle))then
                ConnectedTo(i)=ConnectedTo(i)+2 !Indicates second H is bonded
            end if
            !-----------------------------------------
            !ALT:
            !do j=1,size(OH,1),1
            !   Distance=vec_length(ONu(j))
            !   Angle=ScalarProduct(OH(j),ONu(j))
            !   if(is_Bonded(Distance,Angle))then
            !       ConnectedTo(i)=ConnectedTo(i)+j !Indicates j-th H is bonded
            !   end if
            !end do
        end if
    end do
    end function
    !-------------------------------------------------------------------
    !Function that uses a set of Angle/Distance cutoff criteria to 
    !determine whether a bond exists between two water molecules.
    !Values are taken from JCPL,2,105. 
    !-------------------------------------------------------------------
    logical function is_Bonded(Distance, pAngle)
    implicit none
    character(*),parameter :: FName="is_Bonded"
    character(*),parameter :: CName=MName//":"//FName
    double precision,intent(in)                       :: Distance,pAngle
    double precision,Dimension(56)                    :: DistanceRep
    double precision,Dimension(56)                    :: AngleRep
    double precision                                  :: Angle
    integer                                           :: i
    do i=lbound(DistanceRep,1),ubound(DistanceRep,1),1
        DistanceRep(i)=2.33d0+0.02d0*i
    end do
    AngleRep(1)=3.0d0
    AngleRep(2)=9.0d0
    AngleRep(3)=13.0d0
    AngleRep(4)=15.0d0
    AngleRep(5:6)=17.0d0
    AngleRep(7)=19.0d0
    AngleRep(8)=21.0d0
    AngleRep(9:10)=23.0d0
    AngleRep(11:12)=25.0d0
    AngleRep(13)=27.0d0
    AngleRep(14:15)=29.0d0
    AngleRep(16:18)=31.0d0
    AngleRep(19:21)=33.0d0
    AngleRep(22:26)=35.0d0
    AngleRep(27:44)=37.0d0
    AngleRep(45)=35.0d0
    AngleRep(46)=33.0d0
    AngleRep(47)=31.0d0
    AngleRep(48:49)=29.0d0
    AngleRep(50)=25.0d0
    AngleRep(51)=21.0d0
    AngleRep(52)=19.0d0
    AngleRep(53)=15.0d0
    AngleRep(54)=7.0d0
    AngleRep(55)=3.0d0
    AngleRep(56)=1.0d0
    !-------------------------------------------------
    !Warn for impossible angle when molecules cannot 
    !be the same
    if((.NOT.(pAngle.gt.-1.0d0.AND.pAngle.lt.1.0d0))&
    .AND.(Distance.gt.0.5d0))then
        call PrintError(1,CName,"Invalid value for variable 'pAngle'"&
        //" ("//toChar(pAngle)//").")
        is_Bonded=.FALSE.
        return
    end if
    is_Bonded=.FALSE.
    Angle=acos(pAngle)*360/(2*pi)
    do i=lbound(DistanceRep,1),ubound(DistanceRep,1),1
         if(Distance.le.DistanceRep(i).AND.Angle.le.AngleRep(i))then
             is_Bonded=.TRUE.
             exit
         end if
    end do
    end function
    !-------------------------------------------------------------------
    !Returns positon Array for Water Molecules in the Format O-H-H 
    !-------------------------------------------------------------------
    function WaterArray(AMolecule)
    implicit none
    character(*),parameter :: FName='WaterArray'
    character(*),parameter :: CName=MName//":"//FName
    class(Molecule),intent(in)                        :: AMolecule
    type(Atom),pointer                                :: Core
    real(8),Dimension(:,:),allocatable                :: WaterArray
    integer                                           :: NDim
    logical                                           :: used
    integer                                           :: i
    used=.FALSE.
    nullify(Core)
    Core=>AMolecule%Atoms(1)
    NDim=size(Core%Pos())
    allocate(WaterArray(3,NDim))
    if (AMolecule%MSize().ne.3) then
        call PrintError(3,CName,"Molecule is not Water, found "&
             // toChar(AMolecule%MSize()) // " atoms instead of 3.")
    end if
    do i=1,AMolecule%MSize(),1
        Core=>AMolecule%Atoms(i)
        if(index(Core%Element(),"O").gt.0)then
            WaterArray(1,:)=Core%Pos()
        else if(index(Core%Element(),"H").gt.0.AND.used)then
            WaterArray(3,:)=Core%Pos()
        else if(index(Core%Element(),"H").gt.0)then
            WaterArray(2,:)=Core%Pos()
            used=.TRUE.
        else
            call PrintError(4,CName,"Discovered Atom that doesn't belong "&
            //"to a water molecule.")
        end if
    end do
    end function
end module
