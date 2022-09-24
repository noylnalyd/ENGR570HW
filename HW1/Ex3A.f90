SUBROUTINE Ex3A(N) BIND (C,name='Ex3A')
        USE ISO_C_BINDING, ONLY: C_INT
        IMPLICIT NONE

        
        INTEGER(C_INT),INTENT(IN) :: N
        INTEGER(4) :: i,j,ind ! Indices
    INTEGER,ALLOCATABLE,DIMENSION(:) :: A ! Matrix to fill
    INTEGER(4) :: numSpacesMax,kmax,k
    INTEGER(4),EXTERNAL :: digitCount
    ! Part iii
    ALLOCATE(A(0:N*N))
    DO i=0,(N*N-1)
        ind = 0
        DO j=0,N-1
            IF(IAND(j,1) .gt. 0) THEN
                ind = ind + ISHFT(IAND(ISHFT(i,-j),1),j/2)*N
            ELSE
                ind = ind + ISHFT(IAND(ISHFT(i,-j),1),j/2)
            END IF
        END DO
        A(ind)=i+1
    END DO

    ! Part iv
    numSpacesMax = MAX(digitCount(N*N)+1,2)
    WRITE(*,'(a)',advance="no") "A=["
    DO i=0,N-1
        DO j=0,N-1
            kmax = numSpacesMax - digitCount(A(i*N+j))
            IF (i+j .eq. 0) THEN
                kmax = kmax - 1
            END IF
            IF (digitCount(A(i*N+j)) .eq. 1) THEN
                DO k=0,(kmax-1)
                    WRITE(*,'(a)',advance="no") " "
                END DO
                WRITE(*,'(I1)',advance="no") A(i*N+j)
            ELSEIF (digitCount(A(i*N+j)) .eq. 2) THEN
                DO k=0,(kmax-1)
                    WRITE(*,'(a)',advance="no") " "
                END DO
                WRITE(*,'(I2)',advance="no") A(i*N+j)
            ELSE
                DO k=0,(kmax-1)
                    WRITE(*,'(a)',advance="no") " "
                END DO
                WRITE(*,'(I3)',advance="no") A(i*N+j)
            END IF
        END DO
        IF (i .lt. N-1) THEN
            WRITE(*,'(2a)',advance="no") NEW_LINE('a'),"  "
        ELSE
            WRITE(*,'(2a)',advance="no") "]",NEW_LINE('a')
        END IF
    END DO

END SUBROUTINE Ex3A

INTEGER(4) FUNCTION digitCount(num) result(count)
    IMPLICIT NONE
    !count ! digits in num
    INTEGER(4) :: runprod ! running product for comparison
    INTEGER(4) :: num ! input arg

    count = 1
    runprod = 10
    DO WHILE (runprod .le. num)
        count = count + 1
        runprod = runprod * 10
    END DO

END FUNCTION digitCount