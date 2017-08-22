!----------------------------------------------------------------------
!<summary>
!   provide a linked lists that hold a set of pointers to an atom 
!   object
!</summary>
!----------------------------------------------------------------------
module AtomListMod
use useful, only : PrintError
use ClassStarListMod, only : ListOfClassStar,ClassStarList,GlobListStart
use AtomMod, only : Atom
implicit none
private
character(*),parameter :: MName='AtomListMod'
type                                                  :: Atompnt
    type(Atom),pointer                                :: val
end type
type,extends(ListOfClassStar)                         :: ListOfAtoms
    type(Atompnt)                                     :: Item
    contains
        procedure,pass(this),public                   :: getItem=>LoA_getItem
        procedure,pass(this),public                   :: setItem=>LoA_setItem
        procedure,pass(this),private                  :: Equals=>LoA_Equals
        procedure,nopass,public                       :: Constructor=>LoA_Create
end type
type,extends(ClassStarList)                           :: AtomList
    contains
        procedure,pass(this),public                   :: getItem=>AtL_getItem
end type
interface AtomPnt
    module procedure Atp_Cnstr0
end interface
interface ListOfAtoms
    module procedure LoA_Cnstr0
    module procedure LoA_Cnstr1
    module procedure LoA_Cnstr2
end interface
interface AtomList
    module procedure AtL_Cnstr0
end interface
public :: AtomList
contains
function Atp_Cnstr0()
    implicit none
    type(Atompnt)                                     :: Atp_Cnstr0
    Atp_Cnstr0%val=>null()
end function
!-----------------------------------------------------------------------
! Utility for ListOfAtoms
!<param name="ID">The number from which the list starts</param>
!<param name="Item">     
!-----------------------------------------------------------------------
    function LoA_Cnstr0()
        implicit none
        type(ListOfAtoms)                             :: LoA_Cnstr0
        LoA_Cnstr0%ID=GlobListStart
        LoA_Cnstr0%Item=AtomPnt()
        LoA_Cnstr0%Next=>null()
    end function
    function LoA_Cnstr1(ID)
        implicit none
        type(ListOfAtoms)                             :: LoA_Cnstr1
        integer                                       :: ID
        LoA_Cnstr1%ID=ID
        LoA_Cnstr1%Item=AtomPnt()
        LoA_Cnstr1%Next=>null()
    end function
    function LoA_Cnstr2(Item,ID)
        implicit none
        type(ListOfAtoms)                             :: LoA_Cnstr2
        class(Atom)                                   :: Item
        integer                                       :: ID
        LoA_Cnstr2%ID=ID
        LoA_Cnstr2%Item=AtomPnt()
        LoA_Cnstr2%Next=>null()
    end function
    subroutine LoA_getItem(this,OutItem)
        implicit none
        character(*),parameter :: FName='LoR_getItem'
        character(*),parameter :: CName=MName//":"//FName
        class(ListOfAtoms)                            :: this
        class(*),pointer                              :: OutItem
        select type (OutItem)
            type is (AtomPnt)
                OutItem=this%Item
            class default
                call PrintError(4,CName,"Object 'OutItem' is not of type "//&
                "'AtomPnt'.")
        end select
    end subroutine
    subroutine LoA_setItem(this,Item)
        implicit none
        character(*),parameter :: FName='LoA_setItem'
        character(*),parameter :: CName=MName//":"//FName
        class(ListOfAtoms)                            :: this
        class(*),intent(in)                           :: Item
        select type (Item)
            type is (Atom)
                this%Item%val=>Item
            class default
                call PrintError(4,CName,"Object 'Item' is not of type 'Atom'.")
        end select
    end subroutine
    logical function LoA_Equals(this,Item)
        implicit none
        class(ListOfAtoms)                            :: this
        class(*),intent(in)                           :: Item
        select type(Item)
            type is (Atom)
                !LoA_Equals=this%Item.eq.Item
                LoA_Equals=associated(this%Item%val,Item)
            class default
                LoA_Equals=.FALSE.
        end select
    end function
    subroutine LoA_Create(this)
        implicit none
        character(*),parameter :: FName='LoA_Create'
        character(*),parameter :: CName=MName//":"//FName
        class(ListOfClassStar),pointer,intent(out)        :: this
        nullify(this)
        allocate(ListOfAtoms :: this)
        select type (this)
            type is(ListOfAtoms)
                !this=ListOfAtoms() !<= produces segfault for a reason ...
                this%ID=GlobListStart
                this%Item=AtomPnt()
                this%Next=>null()
        end select
    end subroutine
!-----------------------------------------------------------------------
! Utility for AtomList  
!-----------------------------------------------------------------------
    function AtL_Cnstr0()
        implicit none
        character(*),parameter :: FName='AtL_Cnstr0'
        character(*),parameter :: CName=MName//":"//FName
        type(AtomList)                                :: AtL_Cnstr0
        AtL_Cnstr0%ClassStarList=ClassStarList(ListOfAtoms(),AtomPnt())
    end function
    function AtL_getItem(this,ID)
        implicit none 
        character(*),parameter :: FName='AtL_getItem'
        character(*),parameter :: CName=MName//":"//FName
        class(AtomList)                               :: this
        type(Atom),pointer                            :: AtL_getItem
        class(*),pointer                              :: OutItem
        integer,optional                              :: ID
        if(present(ID))then
            OutItem=>this%ClassStarList%reqItem(ID)
        else
            OutItem=>this%ClassStarList%reqItem()
        end  if
        select type (OutItem)
            type is (AtomPnt)
                AtL_getItem=>OutItem%val
            class default
                call PrintError(4,CName,"Object OutItem is not of "//&
                "type 'AtomPnt'.")
        end select
    end function
end module
