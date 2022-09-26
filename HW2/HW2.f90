PROGRAM HW2
    
    ! Author : Dylan Lyon
    ! Title : Matrix multiplier for sparse matrices
    ! Date : 09/24/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: i,j! Indices

    ! Command args
    CHARACTER(len = 10) :: spfmt ! Sparse format
    INTEGER(4) :: nmults ! Number of matrix-vector mult ops
    CHARACTER(len = 100) :: mmfile ! mm address
    CHARACTER(len = 100) :: vecfilein ! input file addr
    CHARACTER(len = 100) :: vecfileout ! output file addr

    ! mm read variables
    INTEGER(4) :: mm_in,mm_out,nrow,ncol,nnz,nnzmax,indx,jndx,ival
    REAL(4),ALLOCATABLE,DIMENSION(:) :: rval
    REAL(8),ALLOCATABLE,DIMENSION(:) :: dval
    REAL(8),ALLOCATABLE,DIMENSION(:) :: cval
    CHARACTER(14) :: mmid
    CHARACTER(6) :: mmtype
    CHARACTER(10) :: mmrep
    CHARACTER(7) :: mmfield
    CHARACTER(19) :: mmsymm

    ! File input variables
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: A ! Matrix to fill
    REAL(8),ALLOCATABLE,DIMENSION(:) :: x ! Vector input
    INTEGER(4) :: vec_in, vec_out ! File i/o units

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

    ! Set units
    mm_in = 1
    mm_out = 2
    vec_in = 3
    vec_out = 4
    
    ! Open units
    OPEN(mm_in,FILE=mmfile)
    ! OPEN(mm_out,FILE=vecfileout)
    OPEN(vec_in,FILE=vecfilein)
    OPEN(vec_out,FILE=vecfileout)

    ! Read in mm file
    CALL mm_file_read(mm_in, mmid, mmtype, mmrep, mmfield, mmsymm, nrow, ncol, nnz, nnzmax, indx, jndx, ival, rval, dval, cval)
    
    ! Allocate x, b based on size of mm mat. Kyle's code won't be mean
    ALLOCATE(x(ncol))
    ALLOCATE(b(nrow))

    ! Loop???? bad idaeasasassasa

    ! Switch case for spfmt
    SELECT CASE (spfmt)
    ! Dense matrix
    CASE ("DEN")
        DO j=0,ncol-1
            DO i=0,nrow-1
                b(i) = b(i) + A(i,j)*x(j);
            END DO
        END DO
    ! Coordinate
    CASE ("COO")
        ! has N
        
    ! Compressed sparse row
    CASE ("CSR")
        ! has N
        
    ! Jagged diagonal
    CASE ("JDS")
        ! has N
        
    ! ELLPACK, whatever that means
    CASE ("ELL")
        ! has N
        
    END SELECT

    DO i=0,nrow-1
        WRITE(vec_out,*) b(i)
    END DO

    ! Close files
    close(mm_in)
    ! close(mm_out)
    close(vec_in)
    close(vec_out)
    
END PROGRAM HW2