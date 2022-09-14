FUNCTION z_order2D(i,j) RESULT(z)
      INTEGER(4),INTENT(IN) :: i,j
      INTEGER(4) :: z

      INTEGER(4) :: k,ik,jk,zk,ord
      INTEGER(1) :: ib(32),jb(32),zb(64)

      ib=0; jb=0; zb=0;
      ik=i-1; jk=j-1;
      DO k=1,32
        ib(k)=INT(MOD(ik,2),1); ik=ik/2;
        jb(k)=INT(MOD(jk,2),1); jk=jk/2;
      ENDDO

      zk=1;
      DO k=1,32
        zb(zk)=ib(k)
        zb(zk+1)=jb(k)
        zk=zk+2
      ENDDO

      z=1; ord=1
      DO k=1,64
        z=z+zb(k)*ord
        ord=ord+ord
      ENDDO
END FUNCTION z_order2D
