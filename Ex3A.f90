MODULE dslyon_Ex2B
        IMPLICIT NONE
        USE z_order2D
        contains
                INTEGER(4) function dslyon_Ex2Bhelper(N) result(res)
                
                INTEGER(4) :: i,j ! Indices
                INTEGER(4) :: A(N*N) ! Matrix to fill

                DO i=1,N
                        DO j=1,N
                                A((i-1)*N+j)=dslyon_Ex2B(i,j)
                        END DO
                END DO
                
                END FUNCTION dslyon_Ex2Bhelper 
END MODULE
