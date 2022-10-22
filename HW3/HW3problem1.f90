PROGRAM HW3problem1
    
    ! Author : Dylan Lyon
    ! Title : QR factorizer!
    ! Date : 10/22/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: i,j,inz ! Indices
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    INTEGER(8) :: memSize ! Total size of stored matrix.

    ! Command args
    CHARACTER(len = 10) :: spfmt ! Sparse format (DEN here)
    CHARACTER(len = 100) :: mmfile ! mm address

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
    INTEGER(4) :: q_out, r_out ! File i/o units

    ! Matrix storage variables
    !   Dense
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: A ! Matrix to fill for dense
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: Q,R ! Q and R factors of A

    ! Input buffer
    CHARACTER(100) :: buffer

    ! Set verbose
    verbose = 4;

    if (verbose > 3) then
        WRITE(*,*) "Allocated data."
    end if

    ! Read input
    READ(buffer,*) mmfile
    CALL GETARG(1,buffer)
    
    if (verbose > 3) then
        WRITE(*,*) "Read inputs."
    end if

    ! Set units
    mm_in = 1
    q_out = 2
    r_out = 3
    
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

    ! Orthogonalize, 1 column at a time
    DO j=0,(ncol-1)
        ! r_jj = L2 norm of q_:j column
        DO i=0,(nrow-1)
            R(j,j) = R(j,j) + Q(i,j)*Q(i,j)
        END DO
        R(j,j) = DSQRT(R(j,j))








    ! Write to mtx files



    ! Close files
    close(mm_in)
    close(q_out)
    close(r_out)
    
END PROGRAM HW3problem1