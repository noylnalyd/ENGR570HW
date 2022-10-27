module subs
    IMPLICIT NONE
    contains
INTEGER(8) function indexer(i,j,n)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: i,j,n
    indexer = i*n+j
end function indexer
function deindexer(index,n)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: index,n
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: deindexer
    ALLOCATE(deindexer(0:1));
    deindexer(1) = index/n
    deindexer(2) = MOD(index,n)

end function deindexer

REAL(8) function LaplaceEqn(u,u_0,u_1,u_2,u_3)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: u,u_0,u_1,u_2,u_3
    LaplaceEqn = (u_0+u_1+u_2+u_3-4*u)
end function LaplaceEqn

REAL(8) function updater(u,u_0,u_1,u_2,u_3)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: u,u_0,u_1,u_2,u_3
    updater = (u_0+u_1+u_2+u_3)/4.0
end function updater

REAL(8) function laplaceEval(i,j,u,n)
    IMPLICIT NONE
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
IMPLICIT NONE
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
IMPLICIT NONE
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
IMPLICIT NONE
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
    ! Title : Laplacer
    ! Date : 10/22/2022


    ! Functions
    use subs

    IMPLICIT NONE

    



    ! Variable declarations
    INTEGER(8) :: i,j,k ! Indices
    INTEGER(8) :: m,n,ndof,niters ! number of rows/cols, number of u entries, iteration counter
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: A ! Dense Matrix
    REAL(8), ALLOCATABLE, DIMENSION(:) :: x,x_prv ! Current and previous solution to Laplace eqn
    REAL(8), ALLOCATABLE, DIMENSION(:) :: p,p_prv,r,r_prv,v,v_prv,s,t,h
    REAL(8) :: rho,rho_prv,alpha,w,w_prv
    REAL(8) :: tolerance,r0,rL ! Convergence criterion
    REAL(8) :: res,rho ! Convergence criterion
    INTEGER(8) :: nitersEstimate,max_iters ! Estimate based on spectral radius

    

    ! File args
    CHARACTER(2) :: solver
    INTEGER(8) :: grid_size

    ! Input buffer
    CHARACTER(100) :: buffer

    ! Matrix numbers from petsc, poisson 5 stencil
    m = 8
    n = 7

    ! Set verbose
    verbose = 1;

    ! Set tol and max iters
    tolerance = 1e-7
    max_iters = 500000

    if (verbose > 3) then
        WRITE(*,*) "Declared variables."
    end if
    
    if (verbose > 3) then
        WRITE(*,*) "Read inputs."
    end if
    
    ! Initialize and allocate
    ndof = m*n
    alpha = 1
    rho = 1
    rho_prv = 1
    w = 1
    w_prv = 1


    ALLOCATE(A(0:(ndof-1),0:(ndof-1)))
    ALLOCATE(x(0:(ndof-1)))
    ALLOCATE(x_prv(0:(ndof-1)))
    ALLOCATE(r(0:(ndof-1)))
    ALLOCATE(r_prv(0:(ndof-1)))
    ALLOCATE(p(0:(ndof-1)))
    ALLOCATE(p_prv(0:(ndof-1)))
    ALLOCATE(v(0:(ndof-1)))
    ALLOCATE(v_prv(0:(ndof-1)))
    ALLOCATE(s(0:(ndof-1)))
    ALLOCATE(h(0:(ndof-1)))

    DO i=0,(ndof-1)
        x(i) = 1
        x_prv(i) = 0
        r(i) = 0
        r_prv(i) = 0
        p(i) = 0
        p_prv(i) = 0
        v(i) = 0
        v_prv(i) = 0
        s(i) = 0
        h(i) = 0
        DO j=0,(ndof-1)
            A(i,j) = 0
        END DO
    END DO
    DO i=0,m-1
        DO j=0,n-1
    	    A(i,j) = 4;
        END DO
    END DO
    DO i=0,m-1
        DO j=1,n-1
    	    A(i,j-1) = -1;
        END DO
    END DO
    DO i=0,m-1
        DO j=0,n-2
    	    A(i,j+1) = -1;
        END DO
    END DO
    DO i=1,m-1
        DO j=0,n-1
    	    A(i-1,j) = -1;
        END DO
    END DO
    DO i=0,m-2
        DO j=0,n-1
    	    A(i+1,j) = -1;
        END DO
    END DO
    r0 = L2dif(u,u_prv,ndof)
    rL = r0
    DO i=0,(ndof-1)
        u_prv(i) = 1
    END DO

    if (verbose > 3) then
        WRITE(*,*) "Allocated and initialized u and u_prv"
    end if

    res = residual(u,m,n)
    
    SELECT CASE (solver)
    CASE("JI")
        
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

        ! Iterate:
        niters = 1
        call cpu_time(start);
        DO WHILE (niters < max_iters .and. rL > r0*tolerance)
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

