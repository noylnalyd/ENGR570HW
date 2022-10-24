module subs
    IMPLICIT NONE
    contains
INTEGER(8) function indexer(i,j,n)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: i,j,n
    indexer = i*n+j
end function indexer
function deindexer(index,n)
    INTEGER(8), INTENT(IN) :: index,n
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: deindexer
    ALLOCATE(deindexer(0:1));
    deindexer(1) = index/n
    deindexer(2) = MOD(index,n)

end function deindexer

REAL(8) function LaplaceEqn(u,u_0,u_1,u_2,u_3)
    REAL(8), INTENT(IN) :: u,u_0,u_1,u_2,u_3
    LaplaceEqn = (u_0+u_1+u_2+u_3-4*u)
end function LaplaceEqn

REAL(8) function updater(u,u_0,u_1,u_2,u_3)
    REAL(8), INTENT(IN) :: u,u_0,u_1,u_2,u_3
    updater = (u_0+u_1+u_2+u_3)/4.0
end function updater

REAL(8) function laplaceEval(i,j,u,n)
    REAL(8), INTENT(IN), DIMENSION(:) :: u
    INTEGER(8), INTENT(IN) :: i,j,n
    !INTEGER(8), EXTERNAL :: indexer
    !REAL(8), EXTERNAL :: LaplaceEqn
    REAL(8) :: U0,u_0,u_1,u_2,u_3
    laplaceEval = 0
    U0 = u(indexer(i,j,n))
    u_0 = 0
    u_1 = 0
    u_2 = 0
    u_3 = 0
    IF (i>0) then
        u_2 = u(indexer(i-1,j,n))
    END IF
    IF (i<(n-1)) then
        u_0 = u(indexer(i+1,j,n))
    END IF
    IF (j>0) then
        u_1 = u(indexer(i,j-1,n))
    END IF
    IF (j<(n-1)) then
        u_3 = u(indexer(i,j+1,n))
    END IF
    laplaceEval = LaplaceEqn(U0,u_0,u_1,u_2,u_3)
end function laplaceEval

REAL(8) function updaterEval(i,j,u,n)
    REAL(8), INTENT(IN), DIMENSION(:) :: u
    INTEGER(8), INTENT(IN) :: i,j,n
    !INTEGER(8), EXTERNAL :: indexer
    !REAL(8), EXTERNAL :: LaplaceEqn
    REAL(8) :: U0,u_0,u_1,u_2,u_3
    updaterEval = 0
    U0 = u(indexer(i,j,n))
    u_0 = 0
    u_1 = 0
    u_2 = 0
    u_3 = 0
    IF (i>0) then
        u_2 = u(indexer(i-1,j,n))
    END IF
    IF (i<(n-1)) then
        u_0 = u(indexer(i+1,j,n))
    END IF
    IF (j>0) then
        u_1 = u(indexer(i,j-1,n))
    END IF
    IF (j<(n-1)) then
        u_3 = u(indexer(i,j+1,n))
    END IF
    updaterEval = updater(U0,u_0,u_1,u_2,u_3)
end function updaterEval
REAL(8) function L2dif(u,u_prv,ndof)
    REAL(8), INTENT(IN), DIMENSION(:) :: u(0:(ndof-1)),u_prv(0:(ndof-1))
    INTEGER(8), INTENT(IN) :: ndof
    INTEGER(8) :: i
    L2dif = 0
    DO i=0,(ndof-1)
        L2dif = L2dif + (u(i)-u_prv(i))**2
    END DO
    L2dif = DSQRT(L2dif)
end function L2dif
REAL(8) function residual(u,n)
    REAL(8), INTENT(IN), DIMENSION(:) :: u
    INTEGER(8), INTENT(IN) :: n
    INTEGER(8) :: i,j
    !INTEGER(8), EXTERNAL :: indexer
    !REAL(8), EXTERNAL :: LaplaceEqn
    !REAL(8), EXTERNAL :: laplaceEval
    residual = 0
    ! Sides
    DO i=0,(n-1)
        DO j=0,(n-1)
            residual = residual + laplaceEval(i,j,u,n)**2
        END DO
    END DO
    residual = DSQRT(residual)
end function residual

end module subs
PROGRAM HW3problem2
    
    ! Author : Dylan Lyon
    ! Title : QR factorizer!
    ! Date : 10/22/2022


    ! Functions
    use subs

    IMPLICIT NONE

    

    ! Variable declarations
    INTEGER(8) :: i,j,k ! Indices
    INTEGER(8) :: n,ndof,niters ! number of rows/cols, number of u entries, iteration counter
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    REAL(8), ALLOCATABLE, DIMENSION(:) :: u,u_prv ! Current and previous solution to Laplace eqn
    REAL(8) :: tolerance,r0,rL ! Convergence criterion
    REAL(8) :: res,rho ! Convergence criterion
    INTEGER(8) :: nitersEstimate,max_iters ! Estimate based on spectral radius

    

    ! File args
    CHARACTER(2) :: solver
    INTEGER(8) :: grid_size

    ! Input buffer
    CHARACTER(100) :: buffer

    ! Set verbose
    verbose = 1;

    ! Set tol and max iters
    tolerance = 1e-8
    max_iters = 500000

    if (verbose > 3) then
        WRITE(*,*) "Declared variables."
    end if

    ! Read input
    CALL GETARG(1,buffer)
    READ(buffer,*) solver
    CALL GETARG(2,buffer)
    READ(buffer,*) grid_size
    
    if (verbose > 3) then
        WRITE(*,*) "Read inputs."
    end if
    
    ! Here's how to organize this
    ! u contains only inside values
    n = grid_size-2
    ndof = (n)**2
    ALLOCATE(u(0:(ndof-1)))
    ALLOCATE(u_prv(0:(ndof-1)))
    DO i=0,(ndof-1)
        u(i) = 1
        u_prv(i) = 1
    END DO

    if (verbose > 3) then
        WRITE(*,*) "Allocated and initialized u and u_prv"
    end if
    res = residual(u,n)
    r0 = DBLE(n)
    SELECT CASE (solver)
    CASE("JI")
        ! Store initial iteration and its sequel
        DO i=0,(n-1)
            DO j=0,(n-1)
                u(indexer(i,j,n)) = updaterEval(i,j,u_prv,n)
            END DO
        END DO
        !r0 = L2dif(u,u_prv,ndof)
        
        rL = r0
        ! Overwrite u_prv
        DO i=0,(ndof-1)
            u_prv(i) = u(i)
        END DO

        ! Iterate:
        niters = 1
        call cpu_time(start);
        DO WHILE (niters < max_iters .and. rL > tolerance*r0)
            ! Find new u
            DO i=0,(n-1)
                DO j=0,(n-1)
                    u(indexer(i,j,n)) = updaterEval(i,j,u_prv,n)
                END DO
            END DO
            ! Compare to u_prv
            rL = L2dif(u,u_prv,ndof)
            ! Overwrite u_prv
            DO i=0,(ndof-1)
                u_prv(i) = u(i)
            END DO
            niters = niters+1
            !WRITE(*,*) residual(u,n)
        END DO
        call cpu_time(stop);


    CASE("GS")
        ! Store initial iteration and its sequel
        DO i=0,(n-1)
            DO j=0,(n-1)
                u(indexer(i,j,n)) = updaterEval(i,j,u,n)
            END DO
        END DO
        !r0 = L2dif(u,u_prv,ndof)
        rL = r0
        ! Overwrite u_prv
        DO i=0,(ndof-1)
            u_prv(i) = u(i)
        END DO

        ! Iterate:
        niters = 1
        call cpu_time(start);
        DO WHILE (niters < max_iters .and. rL > tolerance*r0)
            ! Find new u
            DO i=0,(n-1)
                DO j=0,(n-1)
                    u(indexer(i,j,n)) = updaterEval(i,j,u,n)
                END DO
            END DO
            ! Compare to u_prv
            rL = L2dif(u,u_prv,ndof)
            ! Overwrite u_prv
            DO i=0,(ndof-1)
                u_prv(i) = u(i)
            END DO
            niters = niters+1
            !WRITE(*,*) residual(u,n)
        END DO
        call cpu_time(stop);
    CASE("RB")
        ! Store initial iteration and its sequel
        ! Black
        DO i=0,(n-1)
            DO j=MOD(i,2),n-1,2
                u(indexer(i,j,n)) = updaterEval(i,j,u_prv,n)
            END DO
        END DO
        ! Red
        DO i=0,(n-1)
            DO j=1-MOD(i,2),n-1,2
                u(indexer(i,j,n)) = updaterEval(i,j,u,n)
            END DO
        END DO
        !r0 = L2dif(u,u_prv,ndof)
        rL = r0
        ! Overwrite u_prv
        DO i=0,(ndof-1)
            u_prv(i) = u(i)
        END DO

        ! Iterate:
        niters = 1
        call cpu_time(start);
        DO WHILE (niters < max_iters .and. rL > tolerance*r0)
            ! Find new u
            ! Black
            DO i=0,(n-1)
                DO j=MOD(i,2),n-1,2
                    u(indexer(i,j,n)) = updaterEval(i,j,u_prv,n)
                END DO
            END DO
            ! Red
            DO i=0,(n-1)
                DO j=1-MOD(i,2),n-1,2
                    u(indexer(i,j,n)) = updaterEval(i,j,u,n)
                END DO
            END DO
            ! Compare to u_prv
            rL = L2dif(u,u_prv,ndof)
            ! Overwrite u_prv
            DO i=0,(ndof-1)
                u_prv(i) = u(i)
            END DO
            niters = niters+1
            !WRITE(*,*) residual(u,n)
        END DO
        call cpu_time(stop);
    END SELECT
    ! Write outputs!
    res = residual(u,n)
    WRITE(*,*) res
    IF (rL/r0 < tolerance) then
        WRITE(*,*) "Converged"
        rho = tolerance**(1.0/(niters))
        nitersEstimate = FLOOR(LOG(1e-6)/LOG(rho))
        WRITE(*,'(A,F10.3)') "solve time (s): ", (stop-start)
        if (verbose>0) then
            WRITE(*,'(A,I10)') "iters: ",niters
            WRITE(*,'(A,ES10.4)') "residual: ",rL/r0
        end if
        WRITE(*,'(A,F10.4)') "Estimated spectral radius: ",rho
        if (verbose>0) then
            WRITE(*,'(A,I10)') "Iterations to reach 10^-6: ",nitersEstimate
        end if
        WRITE(*,'(A,F10.4)') "Average time per iter (ms): ",(stop-start)/(niters)*1000
    ELSE
        WRITE(*,*) "Diverged."
    END IF
    
END PROGRAM HW3problem2

