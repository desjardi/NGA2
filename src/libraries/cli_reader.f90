module cli_reader
  use string
  implicit none
  
contains
  
  subroutine get_command(command,length,status)
    implicit none
    
    character(len=*), intent(out), optional :: command
    integer,          intent(out), optional :: length
    integer,          intent(out), optional :: status
    
    integer :: iarg,narg,ipos
    integer, save :: lenarg
    character(len=str_long), save :: argstr
    logical, save :: getcmd = .true.
    
    !INTEGER, EXTERNAL :: IARGC
    INTEGER :: IARGC

    IF (GETCMD) THEN
       NARG = IARGC()
       IF (NARG > 0) THEN
          IPOS = 1
          DO IARG = 1,NARG
             CALL GETARG(IARG,ARGSTR(IPOS:))
             LENARG = LEN_TRIM(ARGSTR)
             IPOS   = LENARG + 2
             IF (IPOS > LEN(ARGSTR)) EXIT
          END DO
       ELSE
          ARGSTR = ' '
          LENARG = 0
       ENDIF
       GETCMD = .FALSE.
    ENDIF
    IF (PRESENT(COMMAND)) COMMAND = ARGSTR
    IF (PRESENT(LENGTH))  LENGTH  = LENARG
    IF (PRESENT(STATUS))  STATUS  = 0
    
    return
  end subroutine get_command
  
  function command_argument_count()
    implicit none
  
    integer :: command_argument_count
    !INTEGER, EXTERNAL :: IARGC
    INTEGER :: IARGC

    command_argument_count = IARGC()
    
    return
  end function command_argument_count
  
  subroutine get_command_argument(number,value,length,status)
    implicit none
    
    INTEGER         , INTENT(IN)            :: NUMBER
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
    INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
    INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
    
    CHARACTER(LEN=1000) :: TMPVAL
    
    !INTEGER, EXTERNAL :: IARGC
    INTEGER :: IARGC

    IF (NUMBER < 0) THEN
       IF (PRESENT(VALUE )) VALUE  = ' '
       IF (PRESENT(LENGTH)) LENGTH = 0
       IF (PRESENT(STATUS)) STATUS = 1
       RETURN
    ELSE IF (NUMBER > IARGC()) THEN
       IF (PRESENT(VALUE )) VALUE  = ' '
       IF (PRESENT(LENGTH)) LENGTH = 0
       IF (PRESENT(STATUS)) STATUS = 2
       RETURN
    END IF
    
    IF (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
    
    IF (PRESENT(LENGTH)) THEN
       IF (PRESENT(VALUE)) THEN
          LENGTH = LEN_TRIM(VALUE)
       ELSE
          CALL GETARG(NUMBER,TMPVAL)
          LENGTH = LEN_TRIM(TMPVAL)
       END IF
    END IF
    
    IF (PRESENT(STATUS)) STATUS = 0
  
    return
  end subroutine get_command_argument
  
end module cli_reader
