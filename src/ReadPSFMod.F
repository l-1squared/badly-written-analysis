!----------------------------------------------------------------------
!Author: Kristof Karhan
!Purpose:
!Stores the content, except atom cards, of the PSF file in integer arrays
!Each array is 2D, where the first Dim holds the Tupel, the second Dim the 
!Tupel members. All Lists are sorted when a PSF is read. see sortLists for
!Details
!----------------------------------------------------------------------
module ReadPSFMod
use useful, only : PrintError,toChar
use Print_Arrays, only : get_unit
implicit none
public
character(*),parameter :: MName='ReadPSFMod'
enum,bind(C)
    enumerator::psf=0,cmap=0,ext
end enum
enum,bind(C) 
    enumerator::bond,angle,dihedral,improper
end enum
integer           :: ILen=10          !Space an integer takes up in a psf card
integer           :: ILineLen=90      !line length of a psf file
character(*),parameter :: LFormatS='(A90)'!90 supposed to be ILineLen
character,Dimension(:,:),allocatable                  :: Atoms
integer,Dimension(:,:),allocatable                    :: Bonds
integer,Dimension(:,:),allocatable                    :: Angles
integer,Dimension(:,:),allocatable                    :: Dihedrals
integer,Dimension(:,:),allocatable                    :: Impropers
logical                                               :: HasAtom=.FALSE.
logical                                               :: HasBond=.FALSE.
logical                                               :: HasAngle=.FALSE.
logical                                               :: HasDihedral=.FALSE.
logical                                               :: HasImproper=.FALSE.
private :: ReadCard,sort2D,sort1D
contains
!----------------------------------------------------------------------
!Extracts the numbers from a PSF File and stores them in 
!the appropriate Arrays
!The line length is fixed, b/c PSF EXT uses only 90 chars
!----------------------------------------------------------------------
subroutine ReadPSF(FPSF)
    implicit none
    character(*),parameter :: FName='ReadPSF'
    character(*),parameter :: CName=MName//":"//FName
    character(LEN=*),intent(in)                       :: FPSF
    character(LEN=ILineLen)                           :: line
    integer                                           :: IPSF,io_stat
    integer                                           :: lcount,tupelSum,remark
    logical                                           :: LPSF
    lcount=0
    io_stat=0
    remark=-1
    call get_unit(FPSF,unit=IPSF,is_open=LPSF)
    if(.not.LPSF)open(IPSF,FILE=FPSF,STATUS='old',ACTION='read')
    do while (io_stat.eq.0)
        read(IPSF,FMT=LFormatS,iostat=io_stat)line 
        !---------------------------------------------
        !check if line is actually readable
        lcount=lcount+1
        if (io_stat.lt.0)then
            exit
        else if (io_stat.gt.0)then
            call PrintError(3,CName,"Unexpected Error occured in line ("&
            //toChar(lcount)//") of file '"//FPSF//"'.")
        end if 
        !---------------------------------------------
        !ignore the sections following the NTITLE card
        if(lcount.le.remark)cycle
        !line=line(1:scan(line,'#!',BACK=.TRUE.))
        !---------------------------------------------
        !read Cards
        if (index(line,'PSF').gt.0)then              !check PSF format
            line=adjustl(line(index(line,'PSF')+3:len(line)))
            select case(line(1:index(line,' ')))
                case ('EXT')
                    ILen=10
                    ILineLen=90
                case ('PSF','CMAP')
                    ILen=8
                    ILineLen=72
                case default
                    call PrintError(3,CName,"Could not recognize psf file "//&
                    "format named '"//line(1:index(line,' '))//"'. ")
            end select
        else if (index(line,'NTITLE').gt.0)then  !ignore title cards
            read(line,*)tupelSum
            remark=tupelSum+lcount
        else if (index(line,'NATOM').gt.0)then   !read atom card
            read(line,*)tupelSum 
            call ReadAtomCard(IPSF,tupelSum)
            if(tupelSum.gt.0)HasAtom=.TRUE.
        else if (index(line,'NBOND').gt.0)then   !read bond card
            read(line,*)tupelSum
            if(tupelSum.gt.0)HasBond=.TRUE.
            call ReadCard(IPSF,2,tupelSum,Bonds)
        else if (index(line,'NTHETA').gt.0)then  !read angles card
            read(line,*)tupelSum
            if(tupelSum.gt.0)HasAngle=.TRUE.
            call ReadCard(IPSF,3,tupelSum,Angles)
        else if (index(line,'NPHI').gt.0)then    !read dihedrals card
            read(line,*)tupelSum
            if(tupelSum.gt.0)HasDihedral=.TRUE.
            call ReadCard(IPSF,4,tupelSum,Dihedrals)
        else if (index(line,'NIMPHI').gt.0)then  !read impropers card
            read(line,*)tupelSum
            if(tupelSum.gt.0)HasImproper=.TRUE.
            call ReadCard(IPSF,4,tupelSum,Impropers)
        end if 
    end do
    if(.not.LPSF)close(IPSF)
    call sortLists()
end subroutine
!-----------------------------------------------------
!Read the atom card and simply store the contents in 
!a characte array, where the second index is the line
!index
!-----------------------------------------------------
subroutine ReadAtomCard(IFUnit,TupelSum)
    implicit none
    character(*),parameter :: FName='ReadAtomCard'
    character(*),parameter :: CName=MName//":"//FName
    integer,intent(in)                                :: IFUnit
    integer,intent(in)                                :: TupelSum
    character(LEN=:),allocatable                      :: line
    integer                                           :: i,j
    allocate(Atoms(ILineLen,TupelSum))
    allocate(character(LEN=ILineLen) :: line)
    do i=1,TupelSum,1
        read(IFUnit,FMT='(A'//toChar(ILineLen)//')')line
        do j=1,ILineLen,1
            Atoms(j,i)=line(j:j)
        end do
    end do
end subroutine
!-----------------------------------------------------
!Reads the content of an individual card of a PSF 
!File. Routine is called by ReadPSF
!-----------------------------------------------------
subroutine ReadCard(IFUnit,TupelSize,TupelSum,Array)
    implicit none
    character(*),parameter :: FName='ReadCard'
    character(*),parameter :: CName=MName//":"//FName
    integer,intent(in)                                :: IFUnit
    integer,intent(in)                                :: TupelSum,TupelSize
    integer,Dimension(:,:),allocatable,intent(out)    :: Array
    character(LEN=10)                                 :: FormatS
    character(LEN=:),allocatable                      :: line
    integer                                           :: i,j,CStart,CEnd
    integer                                           :: TupPerLine,NLines
    integer                                           :: ActualTupels
    integer                                           :: ActLineLen
    TupPerLine=ILineLen/(ILen*TupelSize)
    NLines=ceiling((1.d0*TupelSum)/TupPerLine)
    ActLineLen=TupPerLine*TupelSize*ILen
    allocate(Array(TupelSum,TupelSize))
    allocate(character(LEN=ActLineLen*NLines) :: line)
    do i=1,NLines
        read(IFUnit,FMT=LFormatS)line((i-1)*ActLineLen+1:i*ActLineLen)
        !+1 == Array start offset
    end do
    ActualTupels=len_trim(line)/(ILen*TupelSize)
    if(ActualTupels.ne.TupelSum)then
        call PrintError(3,CName,'Expected '//toChar(TupelSum)//' Tupels of '&
        //'size '//toChar(TupelSize)//'. Got '//toChar(ActualTupels)//'.')
    end if
    write(FormatS,FMT='(A1,I1,A1,I2,A1)')'(',TupelSize,'I',ILEN,')'
    do i=1,TupelSum,1
        CStart=(i-1)*ILen*TupelSize+1 !+1== Array start offset
        CEnd=i*ILen*TupelSize
        read(line(CStart:CEnd),FMT=FormatS)(Array(i,j),j=1,TupelSize)
    end do
end subroutine
!----------------------------------------------------------------------
! returns the appropriate enum value from text input
!----------------------------------------------------------------------
integer function getArrayTag(which)
    implicit none
    character(*),parameter :: FName='LinkAList'
    character(*),parameter :: CName=MName//":"//FName
    character(LEN=*),intent(in)                       :: which
    getArrayTag=-1
    select case(Which)
        case('Pairs','PAIRS','pairs','pair','Pair','PAIR','NBOND',&
        'nbond','nbnd')
            getArrayTag=bond
        case('angles','Angles','ANGLES','angle','ANGLE','Angle',&
        'NTHETA','ntheta')
            getArrayTag=angle
        case('Dihedrals','dihedrals','DIHEDRALS','Dihedral','dihedral',&
        'DIHEDRAL','nphi','NPHI')
            getArrayTag=dihedral
        case('improper','Improper','IMPROPER','impropers','Impropers',&
        'IMPROPERS','nimphi','NIMPHI')
            getArrayTag=improper 
        case default
            call PrintError(2,CName,"Unrecognized keyword '"//Which//"'.")
            return
    end select
end function 
!----------------------------------------------------------------------
!Sort functions for the PSF Lists. Standard Quicksort algorithms
!
!Bonds are first sorted tupelwise, i.e. the smallest number will be in 
!field 1 of a tupel(==Dim2). and then the whole array is sorted by the 
!value of field 1 of a tupel. 
!For all other arrays the arrangement of atoms in a tupel is not changed,
!because the relative order of atoms is important.
!for angles, the array is sorted after the value in field 2 b/c
!that is the atom where the angle is located.
!----------------------------------------------------------------------
subroutine sortLists()
    implicit none
    integer                                           :: i
    do i=lbound(bonds,1),ubound(bonds,1)
        call sort1D(bonds(i,:),lbound(bonds,2),ubound(bonds,2))
    end do
    if(size(bonds).gt.0)call sort2D(bonds,lbound(bonds,1),ubound(bonds,1),1)
    if(size(angles).gt.0)call sort2D(angles,lbound(angles,1),ubound(angles,1),2)
    if(size(dihedrals).gt.0)call sort2D(dihedrals,lbound(dihedrals,1),&
    ubound(dihedrals,1),1)
    if(size(impropers).gt.0)call sort2D(impropers,lbound(impropers,1),&
    ubound(impropers,1),1)
end subroutine
recursive subroutine sort2D(this,left,right,sortafter)
    implicit none
    integer,Dimension(:,:),intent(inout)              :: this
    integer,intent(in)                                :: left,right
    integer,intent(in)                                :: sortafter
    integer,Dimension(:),allocatable                  :: copy
    integer                                           :: i,j,pivot
    i=left
    j=right
    pivot=this(ceiling((left+right)*0.5),sortafter)
    do while(i.lt.j)
        do while(this(i,sortafter).lt.pivot)
            i=i+1
        end do
        do while(this(j,sortafter).gt.pivot)
            j=j-1
        end do
        if (i.le.j)then
            allocate(copy(lbound(this,2):ubound(this,2)))
            copy=this(i,:)
            this(i,:)=this(j,:)
            this(j,:)=copy
            deallocate(copy)
            i=i+1
            j=j-1
        end if
    end do
    if(left.lt.j) call sort2D(this,left,j,sortafter)
    if(i.lt.right)call sort2D(this,i,right,sortafter)
end subroutine
recursive subroutine sort1D(this,left,right)
    implicit none
    integer,Dimension(:),intent(inout)                :: this
    integer,intent(in)                                :: left,right
    integer                                           :: copy
    integer                                           :: i,j,pivot
    i=left
    j=right
    pivot=this(ceiling(left+right*0.5))
    do while(i.lt.j)
        do while(this(i).lt.pivot)
            i=i+1
        end do
        do while(this(j).gt.pivot)
            j=j-1
        end do
        if (i.le.j)then
            copy=this(i)
            this(i)=this(j)
            this(j)=copy
            i=i+1
            j=j-1
        end if
    end do
    if(left.lt.j) call sort1D(this,left,j)
    if(i.lt.right)call sort1D(this,i,right)
end subroutine
end module
