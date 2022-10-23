PROGRAM HW3problem1
    
    ! Author : Dylan Lyon
    ! Title : QR factorizer!
    ! Date : 10/22/2022

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(4) :: i,j,k ! Indices
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    REAL(8), ALLOCATABLE, DIMENSION(:) :: u ! Current solution to Laplace eqn

    ! File args
    CHARACTER(2) :: solver
    INTEGER(4) :: grid_size

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
    
    
END PROGRAM HW3problem1

INTEGER(8) function indexer(i,j,grid_size)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: i,j,grid_size
    indexer = i*grid_size+j
end function indexer
function deindexer(index,grid_size)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: index,grid_size
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: deindexer
    ALLOCATE(deindexer(0:1));
    deindexer(1) = index/grid_size
    deindexer(2) = MOD(index,grid_size)

end function deindexer