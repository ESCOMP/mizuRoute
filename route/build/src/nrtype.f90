MODULE nrtype
    ! variable types
    integer, parameter     :: I8B = SELECTED_INT_KIND(15)
    integer, parameter     :: I4B = SELECTED_INT_KIND(9)
    integer, parameter     :: I2B = SELECTED_INT_KIND(4)
    integer, parameter     :: I1B = SELECTED_INT_KIND(2)
    integer, parameter     :: SP = KIND(1.0)
    integer, parameter     :: DP = KIND(1.0D0)
    integer, parameter     :: SPC = KIND((1.0,1.0))
    integer, parameter     :: DPC = KIND((1.0D0,1.0D0))
    integer, parameter     :: LGT = KIND(.true.)
    ! common variables
    integer(i4b),parameter :: strLen=256            ! string length
END MODULE nrtype
