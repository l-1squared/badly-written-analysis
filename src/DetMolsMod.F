module DetMolsMod
use useful, only : PrintError
use Print_Arrays, only: get_unit
use AtomMod, only: Atom
use MoleculeMod, only : Molecule
use ClassStarListMod, only : lbound,ubound
use BaseLinkedLists, only : IntList
use ReadPSFMod, only : bonds
implicit none
public
character(*), parameter :: MName='DetMolsMod'
private :: mol_reallocate
contains
    subroutine simpleMolList(fMol,atomCards,molecules)
        implicit none
        character(*),parameter :: FName='mkMolList'
        character(*),parameter :: CName=MName//":"//FName
        character(LEN=*)                              :: fMol
        class(Atom),Dimension(:),intent(in)           :: atomCards
        type(Molecule),Dimension(:),allocatable,intent(out):: molecules
        integer                                       :: ioUnit
        logical                                       :: ioisOpen
        integer                                       :: nMol,nLine,i
        integer                                       :: myMol
        call get_unit(fMol,ioUnit,ioisOpen)
        if(.not.ioisOpen)open(ioUnit,FILE=fMol,STATUS='old',ACTION='read')
        read(ioUnit,FMT='(I6,1x,I7)')nMol,nLine
        allocate(molecules(nMol))
        do i=1,nMol,1
            molecules(i)=Molecule(ID=i) 
        end do
        if(nLine.ne.size(atomCards,1))call PrintError(3,CName,&
        'Number of lines in file incommensurable to number of atoms')
        do i=1,nline,1
            read(ioUnit,FMT='(I6)')myMol
            call molecules(myMol)%append(atomCards(i))
        end do
        if(.not.ioisOpen)close(ioUnit)
    end subroutine
    !----------------------------------------------------------------------
    !subroutine simply reads the resID value of atomCards and constructs
    !molecules based on that property
    !----------------------------------------------------------------------
    subroutine mkMolList(atomCards,molecules)
        implicit none
        character(*),parameter :: FName='mkMolList'
        character(*),parameter :: CName=MName//":"//FName
        class(Atom),Dimension(:),intent(in)           :: atomCards
        type(Molecule),Dimension(:),allocatable,intent(out):: molecules
        integer                                       :: sizeEst,nMols
        character                                     :: chainID,cChainID
        integer                                       :: cResID,resID
        integer                                       :: i
        !--------------------------------------------------
        !Estimate number of molecules by looking at resID
        sizeEst=atomCards(ubound(atomCards,1))%getresID()
        allocate(molecules(sizeEst))
        nMols=0
        do i=lbound(molecules,1),ubound(molecules,1),1
            molecules(i)=Molecule(ID=i)
        end do
        do i=lbound(atomCards,1),ubound(atomCards,1),1
            chainID=atomCards(i)%getChainID()
            resID=atomCards(i)%getResID()
            if (resID.ne.cResID.OR.chainID.ne.cChainID)then
                nMols=nMols+1
                if(nMols.gt.size(molecules,1))call mol_reallocate(molecules,2*nMols)
                cChainID=chainID
                cResID=resID
            end if
            call molecules(nMols)%append(atomCards(i))
        end do
    end subroutine
    subroutine mol_reallocate(molecules,nSize)
        implicit none
        type(Molecule),Dimension(:),allocatable,intent(inout):: molecules
        type(Molecule),Dimension(:),allocatable       :: storage
        integer,intent(in)                            :: nSize
        integer                                       :: i
        allocate(storage(size(molecules)))
        storage=molecules
        deallocate(molecules)
        allocate(molecules(1:nSize))
        do i=size(storage,1),nSize,1
            molecules(i)=Molecule(ID=i)
        end do
        molecules=storage(1:nSize)
    end subroutine
    !----------------------------------------------------------------------
    ! subroutine takes atomCards and bonds (from PSF read). For each atom
    ! (identified by  atomID) it searches for other atoms connected to it
    !
    ! currently doesn't work. In larger systems the routine gets trapped at
    ! a certain atom ID and opens an infinite loop
    !----------------------------------------------------------------------
    subroutine mkMolList_bnd(atomCards,molecules)
        implicit none
        character(*),parameter :: FName='mkMolList_bnd'
        character(*),parameter :: CName=MName//":"//FName
        class(Atom),Dimension(:),intent(in)           :: atomCards
        type(Molecule),Dimension(:),allocatable,intent(out):: molecules
        integer                                       :: i,molCount
        type(IntList)                                 :: inThisMol
        !bonds is from ReadPSF
        call PrintError(4,CName,'Improper Implementation')
        inThisMol=IntList()
        if(.not.allocated(bonds))then
            call PrintError(3,CName, 'no connectivity array.')
        end if
        allocate(molecules(1:100))!just start with that number
        molCount=0
        i=lbound(atomCards,1)
        do while(i.le.ubound(atomCards,1))
            !---------------------------------------------
            !checks for partners of atom, then checks all partners for
            !partners ... unit no adddable no addable partner is found.
            call Oneatom(atomCards(i)%getatID(),inThisMol)
            !---------------------------------------------
            !write all atoms to molecules, can be done automatically, b/c
            !Oneatom already included all possible atoms
            molCount=molCount+1
            call appendToMolecule(atomCards,inThisMol,molecules,molCount)
            !---------------------------------------------
            !select next atom to check. take first atom that does not
            !belong to moecule
            print*,CName,i,atomCards(i)%getatID(),isInList(inThisMol,&
            atomCards(i)%getatID())
            do while(isInList(inThisMol,atomCards(i)%getatID()))
                i=i+1
                if(i.gt.ubound(atomCards,1))then
                    exit   
                end if
            end do
            call inThisMol%Free()
        end do
        call resize(molecules,molCount,molCount)
    end subroutine
    recursive subroutine Oneatom(myID,inThisMol)
        implicit none
        character(*),parameter :: FName='Oneatom'
        character(*),parameter :: CName=MName//":"//FName
        integer,intent(in)                            :: myID
        type(IntList),intent(inout)                   :: inThisMol
        type(IntList),pointer                         :: myParts
        integer                                       :: i,k
        integer,save                                  :: j=0
        integer                                       :: partID
        logical                                       :: inList
        j=j+1
        !--------------------------------------------------
        !not shure why, but I think the compiler likes it.
        allocate(myParts)
        !myParts=IntList()
        call getBonds(myID,myParts)
        !--------------------------------------------------
        !check Partners for atoms that are already in molecule.
        !If so remvoe them from myParts, so that they aren't 
        !checked twice.
        k=lbound(myParts,1)
        do while(k.le.ubound(myParts,1))
            inList=isInList(inThisMol,myParts%getItem(k))
            if(.not.inList)then
                call inThisMol%append(myParts%getItem(k))
            else if (inList)then
                call myParts%delete(k)
            end if
            k=k+1
        end do
        print*,CName,4,j,myParts%getLength()
        !-------------------------------------------------
        ! check the atoms of myParts (contains atoms connected to an atom)
        ! if they are not in the list already, add them to MolList, if they 
        ! are scrap them from myParts, as there is no need to recheck them.
        do i=lbound(myParts,1),ubound(myParts,1),1
            partID=myParts%getItem(i)
            call Oneatom(partID,inThisMol)
        end do
        call myParts%Free()
        deallocate(myParts)
    end subroutine
    subroutine appendToMolecule(atomCards,inmolList,molecules,place)
        implicit none
        character(*),parameter :: FName='appendToMolecule'
        character(*),parameter :: CName=MName//":"//FName
        class(Atom),Dimension(:),intent(in)           :: atomCards
        class(IntList),intent(in)                     :: inmolList
        type(Molecule),Dimension(:),allocatable,intent(inout):: molecules
        integer,intent(in)                            :: place
        integer                                       :: i
        if(inmolList%getLength().eq.0)return
        if (place.ge.size(molecules,1))call resize(molecules,place-1,2*place)
        molecules(place)=Molecule(ID=place)
        do i=lbound(inmolList,1),ubound(inmolList,1),1
            call molecules(place)%append(atomCards(inmolList%getItem(i)))
        end do    
        print*,CName,place,ubound(molecules,1)
    end subroutine
    subroutine resize(molecules,nlmnt,toSize)
        implicit none
        type(Molecule),Dimension(:),allocatable,intent(inout):: molecules
        type(Molecule),Dimension(:),allocatable       :: copy
        integer                                       :: nlmnt,toSize
        allocate(copy(1:nlmnt))
        copy=molecules(1:nlmnt)
        deallocate(molecules)
        allocate(molecules(1:toSize))
        molecules(1:nlmnt)=copy
    end subroutine
    subroutine getBonds(myID,myParts)
        implicit none
        integer,intent(in)                            :: myID 
        type(IntList),intent(out)                     :: myParts
        integer                                       :: i,aval
        !bonds is module array
        !assumes sorted bonds!
        myParts=IntList()
        do i=lbound(bonds,1),ubound(bonds,1),1
            if(bonds(i,1).eq.myID.OR.bonds(i,2).eq.myID)then
                aval=merge(bonds(i,1),bonds(i,2),bonds(i,2).eq.myID)
                call myParts%append(aval)
            else if (bonds(i,1) .gt. myID) then
                exit
            end if
        end do
    end subroutine
    logical function isinList(list,val)
        implicit none
        class(IntList),intent(in)                     :: list
        integer,intent(in)                            :: val
        integer                                       :: i
        isinList = .FALSE.
        do i=lbound(list,1),ubound(list,1)
            if (list%getItem(i).eq.val) isinList = .TRUE.
        end do    
    end function
end module
