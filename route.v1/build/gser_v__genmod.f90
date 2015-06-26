        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 24 10:37:37 2015
        MODULE GSER_V__genmod
          INTERFACE 
            FUNCTION GSER_V(A,X,GLN)
              REAL(KIND=4), INTENT(IN) :: A(:)
              REAL(KIND=4), INTENT(IN) :: X(:)
              REAL(KIND=4) ,OPTIONAL, INTENT(OUT) :: GLN(:)
              REAL(KIND=4) :: GSER_V(SIZE(A))
            END FUNCTION GSER_V
          END INTERFACE 
        END MODULE GSER_V__genmod
