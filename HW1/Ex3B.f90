PROGRAM Ex3B
  
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
  IMPLICIT NONE
INTERFACE
  SUBROUTINE aux(N) BIND(C,name='Ex3Blol')
    
    IMPORT
    INTEGER(C_INT),VALUE :: N
  END SUBROUTINE aux
END INTERFACE

  

  ! Part i
  INTEGER(C_INT) :: N ! Dimension of z curve
  CHARACTER*100 buffer
  CALL GETARG(1,buffer)
  READ(buffer,*) N

  ! Part ii
  IF ( (N/=2) .and. (N/=4) .and. (N/=8) .and. (N/=16) ) THEN
    WRITE(*,*) "Invalid N."
    STOP 1
  END IF

  CALL aux(N)
END PROGRAM Ex3B