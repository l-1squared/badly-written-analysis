#ifndef ISPNTLIST
#define ISPNTLIST 1
#endif
#if ISPNTLIST==1
#define EQUATOR =>
#define ATTRIBUTE ,pointer  
#else
#define EQUATOR =
#define ATTRIBUTE            
#endif
module MoleculeMod
use iso_fortran_env, only : output_unit,error_unit
use useful, only : PrintError,Nullsetzen
use ClassStarListMod, only : lbound,ubound
use AtomMod, only : Atom
use AtomListMod,only : AtomList
use Print_Arrays, only : Printxyz
implicit none
private
character(*),parameter                                :: MName='MoleculeMod'
type                                                  :: Molecule
    private
        type(AtomList)                                :: Cores
        integer                                       :: ID
        integer                                       :: NAt
        real(8)                                       :: Weight
        real(8),Dimension(:),allocatable              :: CeOfMa
    contains
        procedure,pass(this),public                   :: getIndex=>get_MolID
        procedure,pass(this),public                   :: setIndex=>set_MolID
        procedure,pass(this),public                   :: Msize=>get_MolSize
        procedure,pass(this),public                   :: Mass=>get_MolMass
        procedure,pass(this),public                   :: COM=>get_MolCOM
        procedure,pass(this),private                  :: Mol_Append_EP
        procedure,pass(this),private                  :: Mol_Append_C
        procedure,pass(this),public                   :: Delete=>Mol_Delete
        procedure,pass(this),public                   :: Insert=>Mol_Insert
        procedure,pass(this),public                   :: PrintAtoms=>Mol_Print
        procedure,pass(this),public                   :: Atoms=>Mol_getAtom
        procedure,pass(this),public                   :: findAtom=>Mol_findInLst
        procedure,pass(this),public                   :: setAtom=>Mol_setAtom
        procedure,pass(this),public                   :: getAtom=>Mol_getAtom
        procedure,pass(this),public                   :: update=>calc_MolProps
        procedure,pass(this),public                   :: Free=>Mol_Free
        procedure,pass(OutVar),private                :: MoleculeAssign
        generic,public :: Append=>Mol_Append_EP,Mol_Append_C
        generic,public :: assignment(=)=>MoleculeAssign
end type
public :: Molecule
interface Molecule
    module procedure Mol_Cnstr0
    module procedure Mol_Cnstr1I ! argument is ID
    module procedure Mol_Cnstr1A ! argument is an Array of Atom
    module procedure Mol_Cnstr2  ! arugemts are Array of Atom and ID
end interface
contains
!-----------------------------------------------------------------------
!Constructors. First
!-----------------------------------------------------------------------
    !-------------------------------------------------------------------
    !Constructors for type 'Molecule'
    !-------------------------------------------------------------------
    function Mol_Cnstr0()
        implicit none
        character(*),parameter :: FName='Mol_Cnstr0'
        character(*),parameter :: CName=MName//":"//FName
        type(Molecule)                                :: Mol_Cnstr0
        allocate(Mol_Cnstr0%CeOfMa(1:0))
        Mol_Cnstr0%ID=0
        Mol_Cnstr0%Nat=0
        Mol_Cnstr0%Weight=0
        Mol_Cnstr0%Cores=AtomList()
    end function   
    function Mol_Cnstr1I(ID)
        implicit none
        character(*),parameter :: FName='Mol_Cnstr1I'
        character(*),parameter :: CName=MName//":"//FName
        type(Molecule)                                :: Mol_Cnstr1I
        integer                                       :: ID
        allocate(Mol_Cnstr1I%CeOfMa(1:0))
        Mol_Cnstr1I%ID=ID
        Mol_Cnstr1I%Nat=0
        Mol_Cnstr1I%Weight=0
        Mol_Cnstr1I%Cores=AtomList()
    end function   
    function Mol_Cnstr1A(Atom_List)
        implicit none                                      
        type(Molecule)                                :: Mol_Cnstr1A
        class(Atom),Dimension(:)                      :: Atom_List
        integer                                       :: i
        Mol_Cnstr1A%ID=0
        Mol_Cnstr1A%Nat=size(Atom_List)
        Mol_Cnstr1A%Weight=0.0d0
        Mol_Cnstr1A%CeOfMa(:)=0.0d0
        Mol_Cnstr1A%Cores=AtomList()
        do i=lbound(Atom_List,1),ubound(Atom_List,1),1
            call Mol_Cnstr1A%Cores%Append(Atom_List(i))
        end do
    end function
    function Mol_Cnstr2(Atom_List,ID)
        implicit none                                      
        type(Molecule)                                :: Mol_Cnstr2
        class(Atom),Dimension(:)                      :: Atom_List
        integer                                       :: i
        integer                                       :: ID
        Mol_Cnstr2%ID=ID
        Mol_Cnstr2%Nat=size(Atom_List)
        Mol_Cnstr2%Weight=0.0d0
        Mol_Cnstr2%CeOfMa(:)=0.0d0
        Mol_Cnstr2%Cores=AtomList()
        do i=lbound(Atom_List,1),ubound(Atom_List,1),1
            call Mol_Cnstr2%Cores%Append(Atom_List(i))
        end do
    end function
    !-------------------------------------------------------------------
    !<summary>
    !   Prints an xyz formatted list of all members of the molecule
    !</summary>
    !<param name="fOutput">the name of the file to print the list to</param>
    !-------------------------------------------------------------------
    subroutine Mol_Print(this,fOutput)
        implicit none
        character(*),parameter :: FName='Mol_Print'
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule)                               :: this
        type(Atom)                                    :: Core
        character(LEN=*),intent(in)                   :: fOutput
        real(8),Dimension(:,:),allocatable            :: Print_Coord
        character(LEN=2),Dimension(:),allocatable     :: Print_Lmnt
        integer                                       :: i
        allocate(Print_Coord(1:this%Nat,3))
        allocate(Print_Lmnt(1:this%Nat))
        do i=1,this%Cores%getLength(),1
            Core=this%Cores%getItem(i) 
            Print_Coord(i,:)=Core%Pos()
            Print_LMnt(i)=Core%Element()
        end do
        call Printxyz(Print_Lmnt,Print_Coord,FOutput)
    end subroutine
!-----------------------------------------------------------------------
! List Utility
!-----------------------------------------------------------------------
    subroutine Mol_Append_C(this,Item) 
        implicit none
        character(*),parameter :: FName='Mol_Append_C'
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule)                               :: this
        type(Atom)                                    :: Item
        call Item%setresID(this%ID)!set atoms resID to molecules' ID
        call this%Cores%Append(Item)
        call this%update()
    end subroutine
    subroutine Mol_Append_EP(this,El,Pos) 
        implicit none
        character(*),parameter :: FName='Mol_Append_EP'
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule)                               :: this
        character(LEN=*),intent(in)                   :: El
        real(8),Dimension(:),intent(in)               :: Pos
#if ISPNTLIST==0
        type(Atom)                                    :: Item
        Item=Atom(El,Pos,this%Nat+1)
        call Item%setresID(this%ID)!set atoms resID to molecules' ID
        call this%Cores%Append(Item)
        call this%update()
#else
        call PrintError(4,CName,"Subroutine cannot add unallocated "//&
        "atom "//El//"to a pointer list.")
#endif
    end subroutine
    subroutine Mol_Delete(this,ID) 
        implicit none
        class(Molecule)                               :: this
        integer,intent(in),optional                   :: ID
        if(present(ID))then
            call this%Cores%Delete(ID)
        else
            call this%Cores%Delete()
        end if
    end subroutine
    subroutine Mol_Insert(this,Item,ID)
        class(Molecule)                               :: this
        type(Atom)                                    :: Item
        integer,intent(in),optional                   :: ID
        if(present(ID))then
            call this%Cores%Insert(Item,ID)
        else
            call this%Cores%Insert(Item)
        end if
    end subroutine
    subroutine Mol_setAtom(this,Item,ID)
    implicit none
    class(Molecule)                                   :: this
    type(Atom)                                        :: Item
    integer,intent(in),optional                       :: ID
    if(present(ID))then
        call this%Cores%setItem(Item,ID)
    else
        call this%Cores%setItem(Item)
    end if
    end subroutine
    function Mol_getAtom(this,ID)
        implicit none
        class(Molecule)                               :: this
        integer,intent(in),optional                   :: ID
        type(Atom)ATTRIBUTE                           :: Mol_getAtom
        if(present(ID))then
            Mol_getAtom EQUATOR this%Cores%getItem(ID)
        else
            Mol_getAtom EQUATOR this%Cores%getItem()
        end if
    end function
    function Mol_findInLst(this,globID)result(pAtom)
        implicit none
        character(*),parameter :: FName='Mol_findInLst'
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule),intent(in)                    :: this
        class(Atom),pointer                           :: pAtom
        class(Atom),pointer                           :: pHold
        integer,intent(in)                            :: globID
        integer                                       :: i                      
        nullify(pHold,pAtom)
#if ISPNTLIST==0
        call PrintError(4,CName,"This function is not intended "//&
        "to be used in conjunction with a List that returns a "//&
        "copy of 'Atom'. Recompile with different 'AtomList'.")
#else
        do i=1,this%NAt,1
            pHold=>this%Cores%getItem(i)
            if(pHold%getAtID().eq.globID)then
                pAtom=>pHold
                return
            end if    
        end do
        nullify(pAtom)
#endif
    end function    
    subroutine Mol_Free(this)
    implicit none
    class(Molecule)                                   :: this
    call this%Cores%Free()
    end subroutine
!-----------------------------------------------------------------------
!Utility for type Molecule 
!-----------------------------------------------------------------------
    integer function get_MolID(this)
        implicit none
        class(Molecule),intent(in)                    :: this
        get_MolID=this%ID
    end function
    subroutine set_MolID(this,ID)
        implicit none
        class(Molecule),intent(inout)                 :: this
        type(Atom)                                    :: hold
        integer,intent(in)                            :: ID
        integer                                       :: i
        this%ID=ID
        do i=lbound(this%Cores,1),ubound(this%Cores,1),1
            hold=this%Cores%getItem(i)
            call hold%setresID(this%ID)
            call this%Cores%setItem(hold,i)
        end do
    end subroutine
    !-------------------------------------------------------------------
    !get Molecular Mass
    !-------------------------------------------------------------------
    real(8) function get_MolMass(this)
        implicit none
        class(Molecule)                               :: this
        call this%update()
        get_MolMass=this%weight
    end function
    !-------------------------------------------------------------------
    !get molecular center of Mass
    !-------------------------------------------------------------------
    function get_MolCOM(this)
    implicit none
    class(Molecule)                                   :: this
    real(8),Dimension(:),allocatable                  :: get_MolCOM
    call this%update()
    get_MolCOM=this%CeOfMa
    end function
    !-------------------------------------------------------------------
    !get Number of Atoms
    !-------------------------------------------------------------------
    integer function get_MolSize(this)
    implicit none
    class(Molecule)                                   :: this
    get_MolSize=this%Nat
    end function
    !-------------------------------------------------------------------
    !<summary>
    !   recalculate properties that depend on the pointer members
    !   particulary, the molecule's COM and number of atom
    !</summary>
    !<param name="this">The binding name</param>
    !-------------------------------------------------------------------
    subroutine calc_MolProps(this)
        implicit none    
        character(*),parameter :: FName='calc_MolProps'
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule)                               :: this
        type(Atom)                                    :: Hold
        integer                                       :: i
        if(this%Cores%getLength().gt.0 .AND. .not. allocated(this%CeOfMa))then
            Hold=this%Atoms()
            allocate(this%CeOfMa(lbound(Hold%Pos(),1):ubound(Hold%Pos(),1)))
        end if
        if(allocated(this%CeOfMa))this%CeOfMa(:)=0.0d0 
        this%NAt=0
        this%Weight=0.0d0
        do i=lbound(this%Cores,1),ubound(this%Cores,1),1
            Hold=this%Atoms(i)
            this%CeOfMa=this%CeOfMa+Hold%Pos()*Hold%Mass()
            this%Nat=this%Nat+1
            this%Weight=this%Weight+Hold%Mass()
        end do
        this%CeOfMa=this%CeOfMa/this%Weight
    end subroutine
!-----------------------------------------------------------------------
!Assignment routines 
!-----------------------------------------------------------------------
    subroutine MoleculeAssign(OutVar,InVar)
        implicit none
        character(*),parameter :: FName='MoleculeAssign'
        character(*),parameter :: CName=MName//":"//FName
        class(Molecule),intent(inout)                 :: OutVar
        type(Molecule),intent(in)                     :: InVar
        integer                                       :: i
        OutVar%NAt=InVar%NAt
        OutVar%ID=InVar%ID
        OutVar%Weight=InVar%Weight
        allocate(OutVar%CeOfMa(lbound(InVar%CeOfMa,1):ubound(InVar%CeOfMa,1)))
        OutVar%CeOfMa=InVar%CeOfMa
        OutVar%Cores=AtomList()
        !--------------------------------------------------
        !Copy the Atoms List
        do i=lbound(InVar%Cores,1),ubound(InVar%Cores,1),1
            call OutVar%Append(InVar%Atoms(i))
        end do
    end subroutine
end module

