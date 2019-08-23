MODULE nrtype
    ! variable types
    INTEGER, PARAMETER     :: I8B = SELECTED_INT_KIND(15)
    INTEGER, PARAMETER     :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER     :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER     :: I1B = SELECTED_INT_KIND(2)
    INTEGER, PARAMETER     :: SP = KIND(1.0)
    INTEGER, PARAMETER     :: DP = KIND(1.0D0)
    INTEGER, PARAMETER     :: SPC = KIND((1.0,1.0))
    INTEGER, PARAMETER     :: DPC = KIND((1.0D0,1.0D0))
    INTEGER, PARAMETER     :: LGT = KIND(.true.)
    ! common variables
    INTEGER(I4B),PARAMETER :: strLen=256            ! string length
END MODULE nrtype
