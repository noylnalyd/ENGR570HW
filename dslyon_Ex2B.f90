MODULE dslyon_Ex2B
        contains
                function dslyon_Ex2Bhelper(N) result(A)
                INTEGER(4) :: N ! Dimension of z curve
                INTEGER(4) :: i,j ! Indices
                
                

                READ(*,*) N ! Read in N
                WRITE(*,*) N ! Write N
                IF(N/=2 .and. N/=4 .and. N/=8 .and. N/=16)
                        WRITE(stderr,*) "Invalid N."
                        STOP
                END IF
                
                INTEGER(4) :: A(N*N) ! Matrix to fill

                DO i=1,N
                        DO j=1,N
                                A((i-1)*N+j)=dslyon_Ex2B(i,j)
                        END DO
                END DO
                
                END FUNCTION       
END MODULE
