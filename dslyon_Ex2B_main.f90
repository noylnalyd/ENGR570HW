PROGRAM dslyon_Ex2B_main
    USE dslyon_Ex2B
    IMPLICIT NONE
    INTEGER(4) :: N ! Dimension of z curve
    READ(*,*) N ! Read in N
    WRITE(*,*) N ! Write N
    IF ( (N/=2) .and. (N/=4) .and. (N/=8) .and. (N/=16) ) THEN
        WRITE(*,*) "Invalid N."
        STOP
    END IF

END PROGRAM dslyon_Ex2B_main
