!--------------------------------------------------------------------- 
!Author: Kristof Karhan
!Purpose: Provides a class with indiviualized output and, in case of the
!1D histogram a specialized Variance declaration
!Since the introduction of Corr1DMod, Corr2DMod it is effectively 
!deprecated
!--------------------------------------------------------------------- 
module GenGridMod
use Constants, only : pi, infinity, NaN
use useful, only : PrintError 
use Integrators, only : romberg, real_func
use Print_Arrays, only : Print_splot,Print_Array
use StatisticMod, only : confidence
use Corr1DMod, only : Corr1D
use Corr2DMod, only : Corr2D
implicit none
private
character(*),parameter :: MName='GenGridMod'
type,extends(Corr2D)                                  :: GenGrid
    private
    character(LEN=64)                                 :: FOutput
    logical                                           :: IsInitialized
    contains
        procedure,pass(this),public                   :: update=>CtE_updGrid
        procedure,pass(this),public                   :: PrintGrid=>CtE_Print
        procedure,pass(this),public                   :: getFile=>Cte_getFOutput
end type
type,extends(Corr1D)                                  :: GenHist
    private
    character(LEN=64)                                 :: FOutput
    logical                                           :: IsInitialized
    real(8)                                           :: confidence
    contains
        procedure,pass(this),public                   :: update=>GeH_updGrid
        procedure,pass(this),public                   :: PrintHist=>GeH_print
end type
public :: GenGrid,GenHist
interface GenGrid 
    module procedure CtE_Cnstr3
end interface
interface GenHist
    module procedure GeH_Cnstr2
    module procedure GeH_Cnstr3
end interface
contains
!----------------------------------------------------------------------
! Constructors
!----------------------------------------------------------------------
    function CtE_Cnstr3(FOutput,XInterval,YInterval)
    implicit none
    character(*),parameter :: FName='CtE_Cnstr3'
    character(*),parameter :: CName=MName//":"//FName
    type(GenGrid)                                     :: CtE_Cnstr3
    character(LEN=*),intent(in)                       :: FOutput
    double precision,intent(in)                       :: XInterval,YInterval
    CtE_Cnstr3%Corr2D = Corr2D((/XInterval,YInterval/)) 
    CtE_Cnstr3%FOutput = FOutput
    CtE_Cnstr3%IsInitialized = .TRUE.
    end function
    !-------------------------------------------------------------------
    ! Class functions
    !-------------------------------------------------------------------
    function CtE_getFOutput(this)
    implicit none
    character(*),parameter :: FName='getFOutput'
    character(*),parameter :: CName=MName//":"//FName
    class(GenGrid)                                  :: this
    character(LEN=len(this%FOutput))                :: CtE_getFOutput
    CtE_getFOutput=this%FOutput 
    end function
    subroutine CtE_updGrid(this,XVal,YVal,Increment)
    implicit none
    character(*),parameter :: FName='CtE_updGrid'
    character(*),parameter :: CName=MName//":"//FName
    class(GenGrid)                                    :: this
    real(8),intent(in)                                :: XVal,YVal
    real(8),intent(in),optional                       :: Increment
    real(8)                                           :: Incr
    integer                                           :: XBin,YBin
    Incr=1.0d0
    if(present(Increment))Incr=Increment
    call this%check((/XVal,YVal/))
    XBin=this%properBin(XVal,1)
    YBin=this%properBin(YVal,2)
    this%grid(XBin,YBin)=this%grid(XBin,YBin)+Incr
    this%Sqgrid(XBin,YBin)=this%Sqgrid(XBin,YBin)+Incr**2
    this%Hits(XBin,YBin)=this%Hits(XBin,YBin)+1
    end subroutine
    subroutine CtE_Print(this)
    implicit none
    character(*),parameter :: FName='CtE_Print'
    character(*),parameter :: CName=MName//":"//FName
    class(GenGrid)                                    :: this
    if(.not. this%IsInitialized)then
        call PrintError(2,CName,"Grid was not initalized. No Output produced.")
        return
    end if
    call Print_splot(this%Ordinates(1),this%Ordinates(2),this%grid,&
    this%FOutput)
    end subroutine
!-----------------------------------------------------------------------
!Corr Hist functions
!-----------------------------------------------------------------------
    function GeH_Cnstr2(FOutput,XDim)
    implicit none
    character(LEN=*),intent(in)                       :: FOutput
    double precision,intent(in)                       :: XDim
    type(GenHist)                                     :: GeH_Cnstr2
    
    GeH_Cnstr2%Corr1D=Corr1D((/XDim/))
    GeH_Cnstr2%FOutput=FOutput
    GeH_Cnstr2%IsInitialized=.TRUE.
    GeH_Cnstr2%confidence = 0.95d0
    end function
    function GeH_Cnstr3(fOutput, x_dim, confidence)
        implicit none
        character(LEN=*), intent(in)                  :: fOutput
        real(8), intent(in)                           :: x_dim
        real(8), intent(in)                           :: confidence
        type(GenHist)                                 :: GeH_Cnstr3
        GeH_Cnstr3 = GeH_Cnstr2(fOutput,x_dim)
        GeH_Cnstr3%confidence = confidence
    end function
    !------------------------------------------------------------------
    ! update histogram and Print function
    !------------------------------------------------------------------
    subroutine GeH_updGrid(this,XVal,Increment)
    implicit none
    character(*),parameter :: FName='GeH_updGrid'
    character(*),parameter :: CName=MName//":"//FName
    class(GenHist)                                    :: this
    real(8),intent(in)                                :: XVal
    real(8),intent(in),optional                       :: Increment
    real(8)                                           :: Incr
    integer                                           :: XBin
    Incr=1.0d0
    if(present(Increment))Incr=Increment
    call this%check((/XVal/))
    XBin=this%properBin(XVal,1)
    this%grid(XBin)=this%grid(XBin)+Incr
    this%Sqgrid(XBin)=this%Sqgrid(XBin)+Incr**2
    this%Hits(XBin)=this%Hits(XBin)+1
    end subroutine
    subroutine GeH_Print(this)
    implicit none
    character(*),parameter :: FName='GeH_Print'
    character(*),parameter :: CName=MName//":"//FName
    character(*),parameter :: Title_Line="XVal YVal EOM Confidence NSamples"
    class(GenHist)                                    :: this
    real(8),Dimension(:,:),allocatable                :: PrintArray
    integer                                           :: i
    if(.not. this%IsInitialized)then
        call PrintError(2,CName,"Grid was not initalized. No Output produced.")
        return
    end if
    allocate(PrintArray(lbound(this%grid,1):ubound(this%grid,1),5))
    PrintArray(:,1) = this%Ordinates(1)
    PrintArray(:,2) = this%Average() 
    PrintArray(:,3) = this%EOM()
    do i=lbound(PrintArray,1),ubound(PrintArray,1),1
        !PrintArray(i,3) = this%Sqgrid(i)/(this%Hits(i)-1)
        !PrintArray(i,3) = sqrt(PrintArray(i,3)-PrintArray(i,2)**2)
        PrintArray(i,4) = confidence(this%grid(i), this%SqGrid(i), this%Hits(i))
    end do
    PrintArray(:,5) = this%Hits
    call Print_Array(PrintArray, this%FOutput, Title_Line)
    deallocate(PrintArray)
    end subroutine
end module
