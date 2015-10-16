module dataTypes
USE nrtype
! used to create specific data types
! --
! data type containing a name and a variable (double precision)
type namepvar
 character(len=256)    :: varName
 real(dp), pointer     :: varData(:) => null()
endtype namepvar
! data type containing a name and a variable (integer)
type nameivar
 character(len=256)    :: varName
 integer(i4b), pointer :: varData(:) => null()
endtype nameivar
end module dataTypes
