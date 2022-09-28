PROGRAM HW2
    
    ! Author : Dylan Lyon
    ! Title : Matrix multiplier for sparse matrices
    ! Date : 09/24/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: i,j,inz,imult,iCSRrow,lswap! Indices
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    INTEGER(8) :: memSize ! Total size of stored matrix.

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
    character ( len = 1024 ) tmp1

    ! File input variables

    REAL(8),ALLOCATABLE,DIMENSION(:) :: x ! Vector input
    INTEGER(4) :: vec_in, vec_out ! File i/o units

    ! Matrix storage variables
    !   Dense
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: A ! Matrix to fill for dense and ELLPACK
    !   Sparse
    INTEGER(4),ALLOCATABLE,DIMENSION(:) :: csrRowPtr ! CSR row storage
    INTEGER(4),ALLOCATABLE,DIMENSION(:) ::  col_ind,perm,jd_ptr ! Jagged diag vectors 
    REAL(8),ALLOCATABLE,DIMENSION(:) :: jdiag ! jdiag
    INTEGER(4),ALLOCATABLE,DIMENSION(:) :: ellcount ! #elements per row for ellpack
    INTEGER(4),ALLOCATABLE,DIMENSION(:,:) ::  colidx ! colidx for ellpack
    INTEGER(4) :: ellmax,tellmax ! max number of entries in a row for ellpack

    ! Output variables
    REAL(8),ALLOCATABLE,DIMENSION(:) :: b ! Vector output

    CHARACTER(100) :: buffer

    ! Set verbose
    verbose = 0;

    if (verbose > 3) then
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
    if (verbose > 3) then
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

    if (verbose > 3) then
        WRITE(*,*) "Files open."
    end if

    ! Read in mm file
    ! Read header
    call mm_header_read (mm_in, mmid, mmtype, mmrep, mmfield, mmsymm)
    ! ignore comments
    do
        call mm_comment_read (mm_in, tmp1 )
        if ( tmp1(1:1) /= '%' ) then
            exit
        end if
    end do
    ! read size
    call mm_size_read_string ( tmp1, mmrep, mmsymm, nrow, ncol, nnzmax )
    nnz = nnzmax;
    ! Allocate space for indexes
    ALLOCATE(indx(nnzmax))
    ALLOCATE(jndx(nnzmax))
    ALLOCATE(dval(nnzmax))
    ALLOCATE(rval(nnzmax))
    ! read values
    call mm_values_read (mm_in, mmrep, mmfield, nnzmax, indx, jndx, &
    ival, rval, dval, cval )
    ! CALL mm_file_read(mm_in, mmid, mmtype, mmrep, mmfield, mmsymm, nrow, ncol, nnz, nnzmax, indx, jndx, ival, rval, dval, cval)
    
    if (verbose > 3) then
        WRITE(*,*) "Read mm."
    end if
    ! Allocate x, b based on size of mm mat. Kyle's code won't be mean
    ALLOCATE(x(0:(ncol-1)))
    ALLOCATE(b(0:(nrow-1)))

    ! Read in x, the vector to multiply
    
    do i=0,ncol-1
        READ(vec_in,*) x(i)
    end do

    if (verbose > 3) then
        WRITE(*,*) "Read x."
    end if
    ! Loop????


    ! Switch case for spfmt
    SELECT CASE (spfmt)
    ! Dense matrix
    CASE ("DEN")
        if (verbose > 2) then
            WRITE(*,*) "Dense matrix conversion"
        end if
        ALLOCATE(A(0:(nrow-1),0:(ncol-1)))
        DO inz=0,nnz-1
            A(indx(inz),jndx(inz)) = DBLE(rval(inz))
        END DO
        if (verbose > 1) then
            WRITE(*,*) "Dense matrix multiply."
        end if
        call cpu_time(start);
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
        call cpu_time(stop);
        if (verbose>1) then
            memSize = sizeof(A);
            WRITE(*,*) memSize;
            WRITE(*,*) (stop-start)/REAL(nmults);
        end if
    ! Coordinate
    CASE ("COO")
        ! has N
        if (verbose > 2) then
            WRITE(*,*) "COO matrix (no conversion)"
        end if
        if (verbose > 1) then
            WRITE(*,*) "COO matrix multiply."
        end if
        call cpu_time(start);
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            DO inz=0,nnz-1
                b(indx(inz)) = b(indx(inz)) + DBLE(rval(inz)) * x(jndx(inz));
            END DO
        END DO
        call cpu_time(stop);
        if (verbose>0) then
            memSize = sizeof(indx)+sizeof(jndx)+sizeof(rval);
            WRITE(*,*) memSize;
            WRITE(*,*) (stop-start)/REAL(nmults);
        end if
    ! Compressed sparse row
    CASE ("CSR")
        ! has N
        if (verbose > 2) then
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
        if (verbose > 1) then
            WRITE(*,*) "CSR matrix multiply."
        end if
        call cpu_time(start);
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
        call cpu_time(stop);
        if (verbose>1) then
            memSize = sizeof(csrRowPtr)+sizeof(rval);
            WRITE(*,*) memSize;
            WRITE(*,*) (stop-start)/REAL(nmults);
        end if
    ! Jagged diagonal
    CASE ("JDS")
        ! has N
        ! Oh,boy this is complex. Convert to ellpack first, then permute.
        if (verbose > 2) then
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
            if (ellmax < ellcount(indx(inz))) then
                ellmax = ellcount(indx(inz))
            end if
        END DO
        ALLOCATE(A(0:(nrow-1),0:(ellmax-1)));
        ALLOCATE(colidx(0:(nrow-1),0:(ellmax-1)));
        if (verbose > 2) then
            WRITE(*,*) "ELLPack allocated."
        end if
        tellmax = 0;
        iCSRrow = 0;
        DO i=0,nrow-1
            ellcount(i) = 0;
        END DO
        DO inz=0,nnz-1
            ! Append to row
            A(indx(inz),ellcount(indx(inz))) = DBLE(rval(inz));
            colidx(indx(inz),ellcount(indx(inz))) = jndx(inz);
            ellcount(indx(inz)) = ellcount(indx(inz)) +1;
        END DO
        
        ! Now populate jdiag,colidx,perm,jd_ptr
        ALLOCATE(jdiag(0:(nnz-1)))
        ALLOCATE(col_ind(0:(nnz-1)))
        ALLOCATE(perm(0:(nrow-1)))
        ALLOCATE(jd_ptr(0:(ellmax+1)))

        if (verbose > 2) then
            WRITE(*,*) "JDS allocated."
        end if
        ! Somehow sort ellcount indices by values!!!
        ! Maybe the simplest of insertion sorts lol
        DO i=0,nrow-1
            perm(i) = i;
        END DO
        DO i=0,nrow-1
            tellmax = ellcount(perm(i));
            lswap = i;
            DO j=i,nrow-1
                if (ellcount(perm(j))>tellmax) then
                    tellmax = ellcount(perm(j))
                    lswap = j;
                end if
            end do
            ! swap
            tellmax = perm(i);
            perm(i) = perm(lswap);
            perm(lswap) = tellmax;
        end do
        if (verbose > 2) then
            WRITE(*,*) "JDS sorted."
        end if
        ! Now we have perm, get jd_ptr and jdiag
        ! index in jdiag
        inz = 0;
        ! Column element+1 (skip 0)
        ell:DO i=0,ellmax+1
            jd_ptr(i) = inz;
            if (i .ne. 0) then
                row:DO j=0,nrow-1
                    if (ellcount(perm(j))<i) then
                        exit row
                    end if
                    jdiag(inz) = A(perm(j),i-1);
                    col_ind(inz) = colidx(perm(j),i-1);
                    if(col_ind(inz) .lt. 0) then
                        WRITE(*,*) perm(j)
                        WRITE(*,*) ellcount(perm(j))
                    end if
                    inz = inz + 1;
                    if(inz .eq. nnz) then
                        ! exit row
                    end if
                end do row
            end if
        end do ell
        jd_ptr(ellmax+1) = inz
        if (verbose > 1) then
            WRITE(*,*) "Jagged Diag matrix multiply."
        end if
        call cpu_time(start);
        DO imult=1,nmults
            DO i=0,nrow-1
                b(i) = 0;
            END DO
            inz = 0;
            
            DO i=0,(ellmax)
                do j=0,(jd_ptr(i+1)-jd_ptr(i)-1)
                    !write(*,*) inz
                    !write(*,*) col_ind(inz)
                    !writE(*,*) j+jd_ptr(i)
                    b(perm(j)) = b(perm(j)) + jdiag(inz)*x(col_ind(inz));
                    inz = inz+1;
                end do
            END DO
        END DO
        call cpu_time(stop);
        if (verbose>1) then
            WRITE(*,*) (stop-start)/REAL(nmults);
        end if
    ! ELLPACK, whatever that means
    CASE ("ELL")
        ! has N
        if (verbose > 2) then
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
        if (verbose > 2) then
            WRITE(*,*) "ELLPack allocated."
        end if
        ! Now populate A,colidx
        tellmax = 0;
        iCSRrow = 0;
        DO i=0,nrow-1
            ellcount(i) = 0;
        END DO
        DO inz=0,nnz-1
            ! Append to row
            A(indx(inz),ellcount(indx(inz))) = DBLE(rval(inz));
            colidx(indx(inz),ellcount(indx(inz))) = jndx(inz);
            ellcount(indx(inz)) = ellcount(indx(inz)) +1;
        END DO
        if (verbose > 1) then
            WRITE(*,*) "ELLPack matrix multiply."
        end if
        call cpu_time(start);
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
        call cpu_time(stop);
        if (verbose>1) then
            memSize = sizeof(ellcount)+sizeof(colidx);
            WRITE(*,*) memSize
            WRITE(*,*) (stop-start)/REAL(nmults);
        end if
    END SELECT

    if (verbose > 2) then
        WRITE(*,*) "Multiplications finished."
    end if

    DO i=0,nrow-1
        WRITE(vec_out,*) b(i)
    END DO

    if (verbose > 2) then
        WRITE(*,*) "b printing finished."
    end if
    ! Close files
    close(mm_in)
    ! close(mm_out)
    close(vec_in)
    close(vec_out)
    
END PROGRAM HW2