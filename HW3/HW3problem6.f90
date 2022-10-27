module sparseMatrix
    implicit none
    private
  
    type, public :: ELL_Matrix
        INTEGER(8), public :: ndof
        INTEGER(8), public :: colnummax
        INTEGER(8), ALLOCATABLE, DIMENSION(:,:) :: colidx
        INTEGER(8), ALLOCATABLE, DIMENSION(:) :: colnum
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: vals
    contains
        procedure, public :: build => build
        procedure, public :: getVal => getVal
        procedure, public :: addVal => addVal
        procedure, public :: SMVM => SMVM
    end type ELL_Matrix
contains
    subroutine build(this)
        CLASS(ELL_Matrix), intent(inout) :: this
        INTEGER(8) :: i,j
        ALLOCATE(this%vals(0:(this%ndof-1),0:this%colnummax))
        ALLOCATE(this%colidx(0:(this%ndof-1),0:this%colnummax))
        ALLOCATE(this%colnum(0:(this%ndof-1)))
        DO i=0,(this%ndof-1)
            this%colnum(i) = 0
            DO j=0,this%colnummax
                this%colidx(i,j) = -1
                this%vals(i,j) = -1000
            END DO
        END DO
    end subroutine build
    subroutine addVal(this,row,col,val)
        CLASS(ELL_Matrix), intent(inout) :: this
        INTEGER(8), INTENT(IN) :: row,col
        REAL(8), INTENT(IN) :: val
        this%vals(row,this%colnum(row)) = val
        this%colidx(row,this%colnum(row)) = col
        this%colnum(row) = this%colnum(row) + 1
    end subroutine addVal
    REAL(8) function getVal(this,row,col) result(val)
        CLASS(ELL_Matrix), intent(inout) :: this
        INTEGER(8), INTENT(IN) :: row,col
        INTEGER(8) :: j
        DO j=0,this%colnum(row)
            if (this%colidx(row,j) .eq. col) then
                val = this%vals(row,j)
                return
            end if
        end do
    end function getVal
    subroutine SMVM(this,x,out)
        CLASS(ELL_Matrix), intent(inout) :: this
        REAL(8), INTENT(IN), DIMENSION(0:(this%ndof-1)) :: x
        REAL(8), INTENT(OUT), DIMENSION(0:(this%ndof-1)) :: out
        INTEGER(8) i,j
        DO i=0,(this%ndof-1)
            out(i) = 0
            DO j=0,(this%colnum(i)-1)
                out(i) = out(i) + this%vals(i,j) * x(this%colidx(i,j))
            END DO
        END DO
    end subroutine SMVM
end module sparseMatrix
module subs
    use sparseMatrix
    IMPLICIT NONE

    contains
INTEGER(8) function indexer(i,j,n)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: i,j,n
    indexer = i*n+j
end function indexer

SUBROUTINE svp(s,u,out,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    REAL(8), INTENT(IN), DIMENSION(:) :: u(0:(ndof-1))
    REAL(8), INTENT(IN) :: s
    REAL(8), INTENT(INOUT), DIMENSION(0:(ndof-1)) :: out
    INTEGER(8) :: i
    DO i=0,(ndof-1)
        out(i) = s*u(i)
    END DO
end SUBROUTINE svp

SUBROUTINE vcp(u1,u2,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    REAL(8), INTENT(IN), DIMENSION(:) :: u1(0:(ndof-1))
    REAL(8), INTENT(OUT), DIMENSION(0:(ndof-1)) :: u2
    INTEGER(8) :: i
    DO i=0,(ndof-1)
        u2(i) = u1(i)
    END DO
end SUBROUTINE vcp

REAL(8) function vvdot(u1,u2,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    REAL(8), INTENT(IN), DIMENSION(:) :: u1(0:(ndof-1)),u2(0:(ndof-1))
    INTEGER(8) :: i
    vvdot = 0
    DO i=0,(ndof-1)
        vvdot = vvdot + u1(i)*u2(i)
    END DO
end
subroutine vva(u1,u2,out,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    REAL(8), INTENT(IN), DIMENSION(:) :: u1(0:(ndof-1)),u2(0:(ndof-1))
    REAL(8), INTENT(INOUT), DIMENSION(0:(ndof-1)) :: out
    INTEGER(8) :: i
    DO i=0,(ndof-1)
        out(i) = u1(i)+u2(i)
    END DO
end
subroutine vvs(u1,u2,out,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    REAL(8), INTENT(IN), DIMENSION(:) :: u1(0:(ndof-1)),u2(0:(ndof-1))
    REAL(8), INTENT(OUT), DIMENSION(0:(ndof-1)) :: out
    INTEGER(8) :: i
    DO i=0,(ndof-1)
        out(i) = u1(i)-u2(i)
    END DO
end subroutine vvs

subroutine residual(A,x,b,out,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    CLASS(ELL_Matrix), INTENT(INOUT) :: A
    REAL(8), INTENT(IN), DIMENSION(:) :: x(0:(ndof-1)),b(0:(ndof-1))
    REAL(8), INTENT(INOUT), DIMENSION(0:(ndof-1)) :: out
    INTEGER(8) :: i,j
    call A%SMVM(x,out)
    call vva(out,b,out,ndof)
end subroutine residual

REAL(8) function L2norm(u,ndof)
    IMPLICIT NONE
    INTEGER(8), INTENT(IN) :: ndof
    REAL(8), INTENT(IN), DIMENSION(:) :: u(0:(ndof-1))
    INTEGER(8) :: i
    L2norm = 0
    DO i=0,(ndof-1)
        L2norm = L2norm + (u(i))**2
    END DO
    L2norm = DSQRT(L2norm)
end

end module subs
PROGRAM HW3problem6
    
    ! Author : Dylan Lyon
    ! Title : bicgstabber
    ! Date : 10/26/2022


    ! Functions
    use subs
    use sparseMatrix

    IMPLICIT NONE

    ! Variable declarations
    INTEGER(8) :: i,j,k ! Indices
    INTEGER(8) :: m,n,ndof,niters ! number of rows/cols, number of u entries, iteration counter
    INTEGER(4) :: verbose ! Verbose control.
    REAL :: start, stop ! timing record
    TYPE(ELL_Matrix) :: A ! Sparse matrix
    REAL(8), ALLOCATABLE, DIMENSION(:) :: x,x_prv ! Current and previous solution to Laplace eqn
    REAL(8), ALLOCATABLE, DIMENSION(:) :: p,p_prv,r,r_prv,v,v_prv,s,t,h,b,res,dummy
    REAL(8) :: rho,rho_prv,alpha,beta,w,w_prv
    REAL(8) :: tolerance ! Convergence criterion
    INTEGER(8) :: max_iters ! max_iters

    ! File args
    CHARACTER(2) :: solver
    INTEGER(8) :: grid_size

    ! Input buffer
    CHARACTER(100) :: buffer

    ! Matrix numbers from petsc, poisson 5 stencil
    m = 100
    n = 100

    ! Set verbose
    verbose = 1;

    ! Set tol and max iters
    tolerance = 1e-7
    max_iters = 5000

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

    ! Build A
    A%colnummax = 5
    A%ndof = ndof
    call A%build()

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
    ALLOCATE(b(0:(ndof-1)))
    ALLOCATE(dummy(0:(ndof-1)))
    ALLOCATE(res(0:max_iters))

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
        b(i) = 0
    END DO
    DO i=0,m-1
        DO j=0,n-1
            call A%addVal(i,i,DBLE(4));
        END DO
    END DO
    DO i=0,m-1
        DO j=1,n-1
    	    call A%addVal(i,j-1,DBLE(-1));
        END DO
    END DO
    DO i=0,m-1
        DO j=0,n-2
            call A%addVal(i,j+1,DBLE(-1));
        END DO
    END DO
    DO i=1,m-1
        DO j=0,n-1
            call A%addVal(i-1,j,DBLE(-1));
        END DO
    END DO
    DO i=0,m-2
        DO j=0,n-1
            call A%addVal(i+1,j,DBLE(-1));
        END DO
    END DO

    ! Initial residual
    call residual(A,x,b,r,ndof)
    res(0) = L2norm(r,ndof)
    DO i=0,(ndof-1)
        r_prv(i) = r(i)
    END DO

    ! Start iterating!!!
    niters = 1
    call cpu_time(start);
    DO WHILE(niters<max_iters)

        rho = vvdot(r,r_prv,ndof)
        
        beta = (rho/rho_prv)*(alpha/w_prv)
        
        call svp(-w_prv,v_prv,dummy,ndof)
        call vva(p_prv,dummy,dummy,ndof)
        call svp(beta,dummy,dummy,ndof)
        call vva(r_prv,dummy,p,ndof)

        call A%SMVM(p,v)

        alpha = rho/vvdot(r,v,ndof)

        call svp(alpha,p,dummy,ndof)
        call vva(x_prv,dummy,h,ndof)

        ! Check h accuracy
        call residual(A,h,b,dummy,ndof)
        res(niters) = L2norm(dummy,ndof)
        if(res(niters) < tolerance) then
            call vcp(h,x,ndof)
            exit
        end if

        call svp(DBLE(-alpha),v,dummy,ndof)
        call vva(r_prv,dummy,s,ndof)

        call A%SMVM(s,t)

        w = vvdot(t,s,ndof)/vvdot(t,t,ndof)
        
        call svp(w,s,dummy,ndof)
        call vva(h,dummy,x,ndof)

        ! Check x accuracy
        call residual(A,x,b,dummy,ndof)
        res(niters) = L2norm(dummy,ndof)
        if(res(niters) < tolerance) then
            exit
        end if

        call svp(-w,t,dummy,ndof)
        call vva(s,dummy,r,ndof)

        ! Update all prvs x r p v rho w
        rho_prv = rho
        w_prv = w
        DO i=0,(ndof-1)
            x_prv(i) = x(i)
            r_prv(i) = r(i)
            p_prv(i) = p(i)
            v_prv(i) = v(i)
        END DO
        niters = niters + 1
    END DO
    call cpu_time(stop);
    
    WRITE(*,*) "Iteration count: ",niters
    WRITE(*,*) "Average time per iteration (ms): ", (stop-start)/DBLE(niters)*1000
    WRITE(*,*) "Residual: ", res(niters)
    
END PROGRAM HW3problem6

