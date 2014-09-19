       SUBROUTINE INV3B3(A)
C
C Subroutine to find matrix inverse of a 3x3 SYMMETRIC matrix
C
C The six unique elements fit into the matrix as
C
C            | A1  A6  A5 |
C        A = |     A2  A4 |
C            |         A3 |
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
       DIMENSION A(6),COFACT(6)
C
       COFACT(1)= A(2)*A(3)-A(4)*A(4)
       COFACT(2)= A(1)*A(3)-A(5)*A(5)
       COFACT(3)= A(1)*A(2)-A(6)*A(6)
       COFACT(4)= A(5)*A(6)-A(1)*A(4)
       COFACT(5)= A(4)*A(6)-A(2)*A(5)
       COFACT(6)= A(4)*A(5)-A(3)*A(6)
C
C Get value of determinant
C
       DETERM= A(1)*COFACT(1) + A(6)*COFACT(6) + A(5)*COFACT(5)
C
C Put inverse into A
C
       A(1)= COFACT(1)/DETERM
       A(2)= COFACT(2)/DETERM
       A(3)= COFACT(3)/DETERM
       A(4)= COFACT(4)/DETERM
       A(5)= COFACT(5)/DETERM
       A(6)= COFACT(6)/DETERM
C
       RETURN
       END
