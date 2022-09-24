PROGRAM HW2
    
    ! Author : Dylan Lyon
    ! Title : Matrix multiplier for sparse matrices
    ! Date : 09/24/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: N ! Dimension of z curve
    INTEGER(4) :: i,j,ind ! Indices

    INTEGER,ALLOCATABLE,DIMENSION(:) :: A ! Matrix to fill
    INTEGER(4) :: numSpacesMax,kmax,k
    INTEGER(4),EXTERNAL :: digitCount
    ! input variables
    CHARACTER(len = 10) :: spfmt ! Sparse format
    INTEGER(4) :: nmults ! Number of matrix-vector mult ops
    CHARACTER(len = 100) :: mmfile ! mm address
    CHARACTER(len = 100) :: vecfilein ! input file addr
    CHARACTER(len = 100) :: vecfileout ! output file addr

    ! mm read variables
    INTEGER(4) :: input_unit,output_unit,NROW,NCOL,NNZ,NNZMAX,INDX,JNDX,IVAL
    REAL(4),ALLOCATABLE,DIMENSION(:) :: RVAL
    REAL(8),ALLOCATABLE,DIMENSION(:) :: DVAL
    REAL(8),ALLOCATABLE,DIMENSION(:) :: CVAL
    CHARACTER(14) :: mmid
    CHARACTER(6) :: mmtype
    CHARACTER(10) :: mmrep
    CHARACTER(7) :: mmfield
    CHARACTER(19) :: mmsymm

    ! Set mmread inputs
    input_unit = 1

    ! Output variables
    REAL(8),ALLOCATABLE,DIMENSION(:) :: b ! Vector output


    ! Read input
    CHARACTER(100) :: buffer

    CALL GETARG(1,buffer)
    READ(buffer,*) spfmt
    CALL GETARG(2,buffer)
    READ(buffer,*) nmults
    CALL GETARG(3,buffer)
    READ(buffer,*) mmfile
    CALL GETARG(4,buffer)
    READ(buffer,*) vecfilein
    CALL GETARG(5,buffer)
    READ(buffer,*) vecfileout

    ! open mm file
    OPEN(input_unit,FILE=mmfile)

    ! Switch case for spfmt
    SELECT CASE (spfmt)
    ! Dense matrix
    CASE ("DEN")
        mm_file_read(*,)
        mm_file_read(input_unit, id, type, rep, field, symm, nrow, &
  ncol, nnz, nnzmax, indx, jndx, ival, rval, dval, cval) 
        ! has N
        ALLOCATE(b(N))
        DO i=0,N-1
            DO j=0,N-1
                b(j) = b(j) + A(i,j)*b(j);
            END DO
        END DO

    ! Coordinate
    CASE ("COO")
        ! has N
        ALLOCATE(b(N))
    ! Compressed sparse row
    CASE ("CSR")
        ! has N
        ALLOCATE(b(N))
    ! Jagged diagonal
    CASE ("JDS")
        ! has N
        ALLOCATE(b(N))
    ! ELLPACK, whatever that means
    CASE ("ELL")
        ! has N
        ALLOCATE(b(N))

    END SELECT

    ! Close file
    close(input_unit);
    IF ( (N/=2) .and. (N/=4) .and. (N/=8) .and. (N/=16) ) THEN
        WRITE(*,*) "Invalid N."
        STOP 1
    END IF
    
    ! Part iii
    ALLOCATE(A(N*N))
    DO i=0,N*N-1
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

END PROGRAM HW2

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
