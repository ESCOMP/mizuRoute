        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 24 10:37:37 2015
        MODULE GCF_V__genmod
          INTERFACE 
            FUNCTION GCF_V(A,X,GLN)
              REAL(KIND=4), INTENT(IN) :: A(:)
              REAL(KIND=4), INTENT(IN) :: X(:)
              REAL(KIND=4) ,OPTIONAL, INTENT(OUT) :: GLN(:)
              REAL(KIND=4) :: GCF_V(SIZE(A))
            END FUNCTION GCF_V
          END INTERFACE 
        END MODULE GCF_V__genmod
