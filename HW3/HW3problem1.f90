PROGRAM HW3problem1
    
    ! Author : Dylan Lyon
    ! Title : QR factorizer!
    ! Date : 10/22/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: i,j,k,inz ! Indices
    INTEGER(4) :: verbose ! Verbose control.
    REAL(8) :: comp, tmpcomp ! Working computation values
    REAL :: start, stop ! timing record
    INTEGER(4) :: nulli(1) ! Some preposterous null integer

    ! Command args
    CHARACTER(len = 100) :: mmfile ! mm address

    ! mm read variables
    INTEGER(4) :: mm_in,nrow,ncol,nnz,nnzmax
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
    character ( len = 7 ) :: outfield = 'double'
    character ( len = 1024 ) tmp1

    ! File input variables
    INTEGER(4) :: q_out, r_out ! File i/o units

    ! Matrix storage variables
    !   Dense
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: A ! Matrix to fill for dense
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: Q,R ! Q and R factors of A
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: calc ! Matrix for norm compuations

    ! Input buffer
    CHARACTER(100) :: buffer

    ! Set verbose
    verbose = 4;

    if (verbose > 3) then
        WRITE(*,*) "Allocated data."
    end if

    ! Read input
    CALL GETARG(1,buffer)
    READ(buffer,*) mmfile
    
    if (verbose > 3) then
        WRITE(*,*) "Read inputs."
    end if

    ! Set units
    mm_in = 1
    q_out = 2
    r_out = 3
    WRITE(*,*) mmfile
    ! Open units
    OPEN(mm_in,FILE=mmfile)
    OPEN(q_out,FILE="q.mtx")
    OPEN(r_out,FILE="r.mtx")

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

    ! For that one symmetric one (ew)
    if (mmsymm == 'symmetric')then
        ! 2*nnz - diagonal elements
        nnzmax = 2*nnzmax-5489;
    end if

    ! Allocate space for indexes
    ALLOCATE(indx(nnzmax))
    ALLOCATE(jndx(nnzmax))
    ALLOCATE(dval(nnzmax))
    ALLOCATE(rval(nnzmax))

    ! read values
    call mm_values_read (mm_in, mmrep, mmfield, nnz, indx, jndx, &
    ival, rval, dval, cval )
    ! CALL mm_file_read(mm_in, mmid, mmtype, mmrep, mmfield, mmsymm, nrow, ncol, nnz, nnzmax, indx, jndx, ival, rval, dval, cval)

    ! Populate symmetric entries
    if (mmsymm == 'symmetric')then
        inz = nnz+1
        do i=1,nnz
            if(indx(i) .ne. jndx(i))then
                rval(inz) = rval(i)
                indx(inz) = jndx(i)
                jndx(inz) = indx(i)
                inz = inz+1;
            end if
        end do
        nnz = nnzmax;
    end if

    if (verbose > 3) then
        WRITE(*,*) "Read mm."
    end if

    ! Convert to Dense matrix :
    if (verbose > 2) then
        WRITE(*,*) "Dense matrix conversion"
    end if
    ALLOCATE(A(0:(nrow-1),0:(ncol-1)))
    DO j=0,(ncol-1)
        DO i=0,(nrow-1)
            A(i,j)=0;
        END DO
    END DO
    DO inz=1,nnz
        A(indx(inz)-1,jndx(inz)-1) = DBLE(rval(inz))
    END DO

    call cpu_time(start);
    call cpu_time(stop);
    if (verbose > 1) then
        WRITE(*,*) "Dense matrix conversion finished.."
    end if

    ! Do some silly orthogonalization
    ALLOCATE(Q(0:(nrow-1),0:(ncol-1)))
    ALLOCATE(R(0:(nrow-1),0:(ncol-1)))
    DO j=0,(ncol-1)
        DO i=0,(nrow-1)
            Q(i,j)=0;
            R(i,j)=0;
        END DO
    END DO
    if (verbose > 3) then
        WRITE(*,*) "Q,R allocated"
    end if

    ! See: https://www.math.uci.edu/~ttrogdon/105A/html/Lecture23.html
    ! Q starts as a copy of A
    DO i=0,(nrow-1)
        DO j=0,(ncol-1)
            Q(i,j) = A(i,j)
        END DO
    END DO

    ! Orthonormalize, 1 column at a time
    DO j=0,(ncol-1)
        ! r_jj = L2 norm of q_:j column
        DO i=0,(nrow-1)
            R(j,j) = R(j,j) + (Q(i,j)*Q(i,j))
        END DO
        R(j,j) = DSQRT(R(j,j))

        !IF (R(j,j) .lt. 1e-5) then
        !    WRITE(*,*) "Ill conditiond!!!"
        !end if
        ! Normalize q_:j by r_jj
        DO i=0,(nrow-1)
            Q(i,j) = Q(i,j)/R(j,j)
        END DO

        ! Orthogonalize further q_:k by q_:j
        DO k=(j+1),(ncol-1)
            ! r_jk = (q_:j)' * q_:k
            DO i=0,(nrow-1)
                R(j,k) = R(j,k) + Q(i,j)*Q(i,k)
            END DO
            ! q_:k = q_:k - r_jk * q_:j
            DO i=0,(nrow-1)
                Q(i,k) = Q(i,k) - R(j,k) * Q(i,j)
            END DO
        END DO
    END DO

    if (verbose > 1) then
        WRITE(*,*) "Q,R calculated"
    end if

    ! Compute desired outputs (the norms)
    ALLOCATE(calc(0:(nrow-1),0:(ncol-1)))
    DO j=0,(ncol-1)
        DO i=0,(nrow-1)
            calc(i,j) = 0;
        END DO
    END DO
    if (verbose > 3) then
        WRITE(*,*) "calc allocated"
    end if

    ! First, compute trace only of Q^T * Q
    ! Q^T * Q = sum Q(i,:)' * Q(:,j) = sum Q(:,i) * Q(:,j)
    ! So trace is sum Q(:,i)*Q(:,i) aka sumsq of column i!
    comp = 0
    DO i=0,(ncol-1)
        ! compute sum Q(:,i)*Q(:,i)
        tmpcomp = 0
        DO k=0,(nrow-1)
            tmpcomp = tmpcomp + Q(k,i)**2
        END DO
        ! Add error from 1 to comp
        comp = comp + (1-tmpcomp)**2
    END DO
    comp = DSQRT(comp)
    WRITE(*,'(A,E20.5)') "||I-Q^T Q||_2=", comp

    ! Second we actually have to compute QR using calc matrix
    DO i=0,(nrow-1)
        DO k=0,(ncol-1)
            DO j=k,(ncol-1)
                calc(i,j) = calc(i,j) + Q(i,k)*R(k,j)
            END DO
        END DO
    END DO

    ! Now compute frobenius norm
    comp = 0
    DO i=0,(nrow-1)
        DO j=0,(ncol-1)
            comp = comp + (A(i,j)-calc(i,j))**2
        END DO
    END DO
    comp = DSQRT(comp)
    WRITE(*,'(A,E20.5)') "||A-QR||_2=", comp

    ! Write to mtx files
    !  output_unit, id, type, rep, field, symm, nrow, &
    !   ncol, nnz, indx, jndx, ival, rval, dval, cval
    ! (mm_in, mmid, mmtype, mmrep, mmfield, mmsymm)
    
    mmid = '%%MatrixMarket'
    mmtype = 'matrix'
    mmrep = 'array'
    mmsymm = 'general'

    CALL mm_file_write(q_out, mmid, mmtype, mmrep, outfield, mmsymm,nrow,ncol,nrow*ncol, &
        nulli,nulli,nulli,nulli,Q,nulli)

    CALL mm_file_write(r_out, mmid, mmtype, mmrep, outfield, mmsymm,nrow,ncol,nrow*ncol, &
        nulli,nulli,nulli,nulli,R,nulli)


    ! Close files
    close(mm_in)
    close(q_out)
    close(r_out)
    
END PROGRAM HW3problem1