PROGRAM HW2
    
    ! Author : Dylan Lyon
    ! Title : Matrix multiplier for sparse matrices
    ! Date : 09/24/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: i,j,inz,imult,iCSRrow! Indices
    INTEGER(4) :: verbose ! Verbose control.

    ! Command args
    CHARACTER(len = 10) :: spfmt ! Sparse format
    INTEGER(4) :: nmults ! Number of matrix-vector mult ops
    CHARACTER(len = 100) :: mmfile ! mm address
    CHARACTER(len = 100) :: vecfilein ! input file addr
    CHARACTER(len = 100) :: vecfileout ! output file addr

    ! mm read variables
    INTEGER(4) :: mm_in,mm_out,nrow,ncol,nnz,nnzmax
    INTEGER(4),ALLOCATABLE,DIMENSION(:) :: indx,jndx
    INTEGER(4),ALLOCATABLE,DIMENSION(:) :: ival
    REAL(4),ALLOCATABLE,DIMENSION(:) :: rval
    REAL(8),ALLOCATABLE,DIMENSION(:) :: dval
    REAL(8),ALLOCATABLE,DIMENSION(:) :: cval
    CHARACTER(14) :: mmid
    CHARACTER(6) :: mmtype
    CHARACTER(10) :: mmrep
    CHARACTER(7) :: mmfield
    CHARACTER(19) :: mmsymm

    ! File input variables

    REAL(8),ALLOCATABLE,DIMENSION(:) :: x ! Vector input
    INTEGER(4) :: vec_in, vec_out ! File i/o units

    ! Matrix storage variables
    !   Dense
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: A ! Matrix to fill for dense and ELLPACK
    !   Sparse
    INTEGER(4),ALLOCATABLE,DIMENSION(:) :: csrRowPtr ! CSR row storage
    INTEGER(4),ALLOCATABLE,DIMENSION(:) :: ellcount ! #elements per row for ellpack
    INTEGER(4),ALLOCATABLE,DIMENSION(:,:) ::  colidx ! colidx for ellpack
    INTEGER(4) :: ellmax,tellmax ! max number of entries in a row for ellpack

    ! Output variables
    REAL(8),ALLOCATABLE,DIMENSION(:) :: b ! Vector output

    CHARACTER(100) :: buffer
    if (verbose > 0) then
        WRITE(*,*) "Allocated data."
    end if
    ! Read input
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
    if (verbose > 0) then
        WRITE(*,*) "Read inputs."
    end if
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

    ! set verbose
    verbose = 1
    if (verbose > 0) then
        WRITE(*,*) "Files open."
    end if

    ! Read in mm file
    ! Initialize nnzmax
    nnzmax = 10000000;
    ! Allocate space for indexes
    ALLOCATE(indx(nnzmax))
    ALLOCATE(jndx(nnzmax))
    ALLOCATE(dval(nnzmax))
    ALLOCATE(rval(nnzmax))
    CALL mm_file_read(mm_in, mmid, mmtype, mmrep, mmfield, mmsymm, nrow, ncol, nnz, nnzmax, indx, jndx, ival, rval, dval, cval)
    
    if (verbose > 0) then
        WRITE(*,*) "Read mm."
    end if
    ! Allocate x, b based on size of mm mat. Kyle's code won't be mean
    ALLOCATE(x(0:(ncol-1)))
    ALLOCATE(b(0:(nrow-1)))

    ! Read in x, the vector to multiply
    
    do i=0,ncol-1
        READ(vec_in,*) x(i)
    end do

    if (verbose > 0) then
        WRITE(*,*) "Read x."
    end if
    ! Loop????


    ! Switch case for spfmt
    SELECT CASE (spfmt)
    ! Dense matrix
    CASE ("DEN")
        if (verbose > 0) then
            WRITE(*,*) "Dense matrix conversion"
        end if
        ALLOCATE(A(0:(nrow-1),0:(ncol-1)))
        DO inz=0,nnz-1
            A(indx(inz),jndx(inz)) = DBLE(rval(inz))
        END DO
        if (verbose > 0) then
            WRITE(*,*) "Dense matrix multiply."
        end if
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            DO j=0,ncol-1
                DO i=0,nrow-1
                    b(i) = b(i) + A(i,j)*x(j);
                END DO
            END DO
        END DO
    ! Coordinate
    CASE ("COO")
        ! has N
        if (verbose > 0) then
            WRITE(*,*) "COO matrix (no conversion)"
        end if
        if (verbose > 0) then
            WRITE(*,*) "COO matrix multiply."
        end if
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            DO inz=0,nnz-1
                b(indx(inz)) = b(indx(inz)) + DBLE(rval(inz)) * x(jndx(inz));
            END DO
        END DO
    ! Compressed sparse row
    CASE ("CSR")
        ! has N
        if (verbose > 0) then
            WRITE(*,*) "CSR conversion"
        end if
        ellmax = 0;
        tellmax = 0;
        iCSRrow = 0;
        ALLOCATE(ellcount(0:(nrow-1)));
        DO i=0,nrow-1
            ellcount(i) = 0;
        END DO
        DO inz=0,nnz-1
            ellcount(indx(inz)) = ellcount(indx(inz))+1;
            ellmax = MAX0(ellmax,ellcount(indx(inz)));
        END DO
        ALLOCATE(csrRowPtr(0:(nrow)))
        iCSRrow = 0;
        csrRowPtr(0) = 0;
        DO i=1,(nrow-1)
            csrRowPtr(i) = ellcount(i-1)+csrRowPtr(i-1);
        END DO
        csrRowPtr(nrow) = nnz;
        if (verbose > 0) then
            WRITE(*,*) "CSR matrix multiply."
        end if
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            DO i=0,nrow-1
                DO inz=csrRowPtr(i),(csrRowPtr(i+1)-1)
                    b(i) = b(i) + DBLE(rval(inz))*x(jndx(inz));
                end do
            END DO
        END DO
    ! Jagged diagonal
    CASE ("JDS")
        ! has N
        ! Oh,boy this is complex. Convert to ellpack first, then permute.
        if (verbose > 0) then
            WRITE(*,*) "JDS conversion."
        end if
        ellmax = 0;
        tellmax = 0;
        iCSRrow = 0;
        ALLOCATE(ellcount(0:(nrow-1)));
        DO i=0,nrow-1
            ellcount(i) = 0;
        END DO
        DO inz=0,nnz-1
            ellcount(indx(inz)) = ellcount(indx(inz))+1;
            ellmax = MAX0(ellmax,ellcount(indx(inz)));
        END DO
        ALLOCATE(A(0:(nrow-1),0:(ellmax-1)));
        ALLOCATE(colidx(0:(nrow-1),0:(ellmax-1)));
        if (verbose > 0) then
            WRITE(*,*) "ELLPack allocated."
        end if
        ! Now populate A,colidx
        tellmax = 0;
        iCSRrow = 0;
        DO inz=0,nnz-1
            ! Append to row
            if (iCSRrow == indx(inz)) then
                A(indx(inz),tellmax) = DBLE(rval(inz));
                colidx(indx(inz),tellmax) = jndx(inz);
                tellmax = tellmax + 1;
            ! Next row
            else
                tellmax = 0;
                iCSRrow = iCSRrow + 1;
                A(indx(inz),tellmax) = DBLE(rval(inz));
                colidx(indx(inz),tellmax) = jndx(inz);
                tellmax = tellmax + 1;
            end if
        END DO
        if (verbose > 0) then
            WRITE(*,*) "ELLPack matrix multiply."
        end if
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            DO i=0,nrow-1
                do j=0,(ellcount(i)-1)
                    b(i) = b(i) + A(i,j)*x(colidx(i,j));
                end do
            END DO
        END DO
    ! ELLPACK, whatever that means
    CASE ("ELL")
        ! has N
        if (verbose > 0) then
            WRITE(*,*) "ELLPACK conversion."
        end if
        ! First need to count max elements in a row to get el array, colidx array
        ellmax = 0;
        tellmax = 0;
        iCSRrow = 0;
        ALLOCATE(ellcount(0:(nrow-1)));
        DO i=0,nrow-1
            ellcount(i) = 0;
        END DO
        DO inz=0,nnz-1
            ellcount(indx(inz)) = ellcount(indx(inz))+1;
            ellmax = MAX0(ellmax,ellcount(indx(inz)));
        END DO
        ALLOCATE(A(0:(nrow-1),0:(ellmax-1)));
        ALLOCATE(colidx(0:(nrow-1),0:(ellmax-1)));
        if (verbose > 0) then
            WRITE(*,*) "ELLPack allocated."
        end if
        ! Now populate A,colidx
        tellmax = 0;
        iCSRrow = 0;
        DO inz=0,nnz-1
            ! Append to row
            if (iCSRrow == indx(inz)) then
                A(indx(inz),tellmax) = DBLE(rval(inz));
                colidx(indx(inz),tellmax) = jndx(inz);
                tellmax = tellmax + 1;
            ! Next row
            else
                tellmax = 0;
                iCSRrow = iCSRrow + 1;
                A(indx(inz),tellmax) = DBLE(rval(inz));
                colidx(indx(inz),tellmax) = jndx(inz);
                tellmax = tellmax + 1;
            end if
        END DO
        if (verbose > 0) then
            WRITE(*,*) "ELLPack matrix multiply."
        end if
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            DO i=0,nrow-1
                do j=0,(ellcount(i)-1)
                    b(i) = b(i) + A(i,j)*x(colidx(i,j));
                end do
            END DO
        END DO
        
    END SELECT

    if (verbose > 0) then
        WRITE(*,*) "Multiplications finished."
    end if

    DO i=0,nrow-1
        WRITE(vec_out,*) b(i)
    END DO

    if (verbose > 0) then
        WRITE(*,*) "b printing finished."
    end if
    ! Close files
    close(mm_in)
    ! close(mm_out)
    close(vec_in)
    close(vec_out)
    
END PROGRAM HW2