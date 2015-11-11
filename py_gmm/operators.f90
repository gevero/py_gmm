MODULE operators

USE kinds
USE datatypes

INTERFACE OPERATOR (==)
	MODULE PROCEDURE equals_rdbl_char
END INTERFACE

CONTAINS

LOGICAL FUNCTION equals_rdbl_char(left_hand,right_hand)

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (RDBL_CHAR), INTENT(IN) :: left_hand, right_hand

equal_if: IF ( (left_hand%value == right_hand%value) .AND. (right_hand%text == left_hand%text) ) THEN
			equals_rdbl_char = .TRUE.
		  ELSE
			equals_rdbl_char = .FALSE.
		  END IF equal_if	 
END FUNCTION equals_rdbl_char
 


END MODULE operators