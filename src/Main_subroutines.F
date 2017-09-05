!----------------------------------------------------------------------- 
!List of subroutines that are added to file Main.F at compile time
!----------------------------------------------------------------------- 
subroutine read_input(finput, fqmmmtrajlist, nreplica, fqmmmmol)
    implicit none
    character(*), parameter :: FName = 'read_input'
    character(*), parameter :: CName = MName // ":" // FName
    character(*), intent(in)                        :: finput
    character(*), intent(out)                       :: fqmmmtrajlist
    character(*), intent(out)                       :: fqmmmmol
    integer, intent(out)                            :: nreplica
    integer                                         :: iounit
    logical                                         :: ioisopen
    character(LEN=globlinelen)                      :: qmmm_trajectory_file
    character(LEN=globlinelen)                      :: file_molecule_association
    integer                                         :: number_of_replicas
    namelist/general/ qmmm_trajectory_file, number_of_replicas,&
                      file_molecule_association 

    call get_unit(finput, iounit, ioisopen)
    if (.not. ioisopen) open(iounit, FILE=finput, STATUS='old',&
                             ACTION='read')
    read(iounit, NML=general)
    if(debug) write(output_unit, NML=general)
    if (.not. ioisopen) close(iounit, STATUS='keep')

    fqmmmtrajlist = qmmm_trajectory_file
    nreplica = number_of_replicas
    fqmmmmol = file_molecule_association
end subroutine


subroutine PrintProgramData(indent)
    implicit none
    character(*),parameter :: FName="PrintProgramData"
    character(*),parameter :: CName=MName//":"//FName
    character(*),parameter :: Prog="PROGRAM"
    integer,intent(in),optional                       :: indent
    integer                                           :: ioUnit
    integer                                           :: idt
    ioUnit=output_unit
    idt=merge(indent,1,present(indent))
    write(ioUnit,FMT='('//toChar(idt)//'x,A'//toChar(len(Prog))//&
    ',1x,A04,A1,T30,A50)')Prog,"NAME",":",MName
    write(ioUnit,FMT='('//toChar(idt)//'x,A11,A1,T30,A50)')&
    "GIT VERSION",":",__GIT_VERSION
    write(ioUnit,FMT='('//toChar(idt)//'x,A'//toChar(len(Prog)) //&
    ',1x,A11,A1,T30,A50)') Prog, "COMPILED AT", ":", __PROGRAM_COMPILE_TIME
    write(ioUnit,FMT='('//toChar(idt)//'x,A'//toChar(len(Prog)) //&
    ',1x,A10,A1,T56,A24)')Prog, "STARTED AT", ":", date_time()
    write(ioUnit,FMT='('//toChar(idt)//'x,A'//toChar(len(Prog))//',1x,A11,A1,T30,A50)')&
    Prog,"COMPILED ON",":",__MACHINE_NAME
end subroutine
!-----------------------------------------------------------------------
!<summary>
!  obtain the current local date and time in a proper format
!</summary>
!<returns> the formatted date and time string </params>
!-----------------------------------------------------------------------
function date_time()
    implicit none
    character(LEN=24)                         :: date_time
    character(LEN=8)                          :: date
    character(LEN=10)                         :: time
    character(LEN=5)                          :: zone
    integer,Dimension(8)                      :: values
    character(LEN=4)                          :: timezone
    call date_and_time(date,time,zone,values)
    if (values(4)/60 .eq. 1) then
        timezone = " CET"
    else if (values(4)/60 .eq. 2) then
        timezone = "CEST"
    end if
    write(date_time, FMT='(2(I2.2,A1),I4,1x,2(I2.2,A1),I2,1x,A4)')&
    values(3),".",values(2),".",values(1),values(5),":",values(6),":",&
    values(7),timezone
end function
