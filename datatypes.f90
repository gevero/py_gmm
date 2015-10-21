MODULE datatypes

USE Kinds

IMPLICIT NONE
SAVE

!******************************************************************************
! DERIVED DATA TYPES
!******************************************************************************

TYPE :: RDBL_CHAR
	REAL(dbl) :: value
	CHARACTER(len=length) :: text
END TYPE

!******************************************************************************
! LINKED LIST
!******************************************************************************


TYPE :: NODE_SGL_REAL					! Linked list reale single precision
	REAL(sgl) :: value
	TYPE (NODE_SGL_REAL), POINTER :: next
END TYPE

TYPE :: NODE_DBL_REAL					! Linked list reale double precision
	REAL(dbl) :: value
	TYPE (NODE_DBL_REAL), POINTER :: next
END TYPE

TYPE :: NODE_SGL_CMPLX					! Linked list complesso single precision
	COMPLEX(sgl) :: value
	TYPE (NODE_SGL_CMPLX), POINTER :: next
END TYPE

TYPE :: NODE_DBL_CMPLX					! Linked list complesso double precision
	COMPLEX(dbl) :: value
	TYPE (NODE_DBL_CMPLX), POINTER :: next
END TYPE


TYPE :: NODE_SHORT_INT					! Linked list integer short
	INTEGER(short) :: value
	TYPE (NODE_SHORT_INT), POINTER :: next
END TYPE

TYPE :: NODE_LONG_INT					! Linked list integer long
	INTEGER(lo) :: value
	TYPE (NODE_LONG_INT), POINTER :: next
END TYPE

TYPE :: NODE_CHAR						! Linked list character
	CHARACTER(len=length) :: value
	TYPE (NODE_CHAR), POINTER :: next
END TYPE

TYPE :: NODE_RDBL_CHAR						! Linked list character
	TYPE (RDBL_CHAR) :: value
	TYPE (NODE_RDBL_CHAR), POINTER :: next
END TYPE

END MODULE datatypes