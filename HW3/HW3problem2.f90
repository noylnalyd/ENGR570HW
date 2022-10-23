PROGRAM HW3problem1
    
    ! Author : Dylan Lyon
    ! Title : QR factorizer!
    ! Date : 10/22/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(8) :: i,j,k ! Indices
    INTEGER(8) :: ndof ! number of u entries
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    REAL(8), ALLOCATABLE, DIMENSION(:) :: u,u_prv ! Current and previous solution to Laplace eqn
    REAL(8) :: tolerance ! Convergence criterion

    ! File args
    CHARACTER(2) :: solver
    INTEGER(8) :: grid_size

    ! Input buffer
    CHARACTER(100) :: buffer

    ! Set verbose
    verbose = 4;

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
    ndof = (grid_size-2)**2
    ALLOCATE(u(0:(ndof-1)))
    ALLOCATE(u_prv(0:(ndof-1)))
    DO i=0,(ndof-1)
        u(i) = 1
        u_prv(i) = 1
    END DO

    if (verbose > 3) then
        WRITE(*,*) "Allocated and initialized u and u_prv"
    end if

    SELECT CASE (solver)
    CASE("JI")
        ! Do all edges (Assume grid size > 3)
        

    CASE("GS")
    
    CASE("RB")

    END SELECT

    
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
    updater = (u_0+u_1+u_2+u_3-3*u)
end function updater

REAL(8) function laplaceEval(i,j,u,n)
    IMPLICIT NONE
    REAL(8), INTENT(IN), DIMENSION(:) :: u
    INTEGER(8), INTENT(IN) :: i,j,n
    INTEGER(8), EXTERNAL :: indexer
    REAL(8), EXTERNAL :: LaplaceEqn
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

REAL(8) function residual(u,n)
    IMPLICIT NONE
    REAL(8), INTENT(IN), DIMENSION(:) :: u
    INTEGER(8), INTENT(IN) :: n
    INTEGER(8) :: i,j
    INTEGER(8), EXTERNAL :: indexer
    REAL(8), EXTERNAL :: LaplaceEqn
    residual = 0
    ! Sides
    i = indexer(1,1,n)
    DO i=1,(n-1)
        ! Left
        residual = residual + LaplaceEqn(u(indexer(i,0,n)),&
            u(indexer(i+1,0,n)), &
            u(indexer(i-1,0,n)), &
            u(indexer(i,1,n)), &
            0.0)
        ! Top
        residual = residual + LaplaceEqn(u(indexer(0,i,n)),&
            u(indexer(0,i+1,n)), &
            u(indexer(0,i-1,n)), &
            u(indexer(1,i,n)), &
            0.0)
        ! Right
        residual = residual + LaplaceEqn(u(indexer(i,n-1,n)),&
            u(indexer(i+1,n-1,n)), &
            u(indexer(i-1,n-1,n)), &
            u(indexer(i,n-2,n)), &
            0.0)
        ! Bottom
        residual = residual + LaplaceEqn(u(indexer(n-1,i,n)),&
            u(indexer(n-1,i+1,n)), &
            u(indexer(n-1,i-1,n)), &
            u(indexer(n-2,i,n)), &
            0.0)
    END DO
    ! Corners
    ! topleft
    residual = residual + LaplaceEqn(u(indexer(0,0,n)),&
        u(indexer(i+1,0,n)), &
        u(indexer(i-1,0,n)), &
        0.0, &
        0.0)
    ! topright
    residual = residual + LaplaceEqn(u(indexer(0,n-1,n)),&
        u(indexer(0,i+1,n)), &
        u(indexer(0,i-1,n)), &
        u(indexer(1,i,n)), &
        0.0)
    ! bottomleft
    residual = residual + LaplaceEqn(u(indexer(n-1,0,n)),&
        u(indexer(i+1,n-1,n)), &
        u(indexer(i-1,n-1,n)), &
        u(indexer(i,n-2,n)), &
        0.0)
    ! Bottomright
    residual = residual + LaplaceEqn(u(indexer(n-1,n-1,n)),&
        u(indexer(n-1,i+1,n)), &
        u(indexer(n-1,i-1,n)), &
        u(indexer(n-2,i,n)), &
        0.0)
end function residual

END PROGRAM HW3problem1