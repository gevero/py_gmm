MODULE basicsubs

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! INDICE DELLE SUBROUTINE CONTENUTE NEL MODULO
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines
! 1) addnode (head, input, error)
!		AddNode_r_sgl,AddNode_r_dbl,AddNode_c_sgl,AddNode_c_dbl,AddNode_i_short
!		AddNode_i_long,AddNode_char,AddNode_rsgl_char
!	COMMENT: add a node to a linked list
!	ERROR:	0) No problem
!			1) Problem in the allocation of the new node
!
!
! 2) lldesctructor (head, error)
!		lldest_r_sgl,lldest_r_dbl,lldest_c_sgl,lldest_c_dbl,lldest_i_short,
!		lldest_i_long,lldest_char
!	COMMENT: deallocate a whole linked list
!	ERROR:	0) No problem
!			1) The head, argument of the subroutine, was not allocated
!
!
! 3) c_addnode (head, input, error [, ncount])
!		C_AddNode_r_sgl,C_AddNode_r_dbl,C_AddNode_c_sgl,C_AddNode_c_dbl,
!		C_AddNode_i_short,C_AddNode_i_long,C_AddNode_char
!	COMMENT: add a node to a linked list if the value to be inserted is not present
!			in the linked list, ncount is an optional argument, if present it gives
!			the actual length of the linked list
!	ERROR:	0) No problem
!			1) An error occured in the subroutine addnode
!
!
!4) pattern (head1, input, head2, ncount, error)
!	COMMENT: creates a pattern which marks equals spheres
!	ERROR:	0) No problem
!			1) An error occured in the subroutine addnode
!
!
!5) dist (r1,r2)
!	COMMENT: eucleidian distance between two points
!
!
!6) cart_spher(x,y,z,r,theta,phi, error)
!	COMMENT: conversion from cartesia to spherica coordinates
!	ERROR:	0) No problem
!			
!7) lnf(z)
!	COMMENT: log of factorial in double precision
!	ERROR:	 -1.0d0 as flag value for problems 
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$






USE kinds		! Modulo kinds per le diverse precisioni	
USE datatypes	! Modulo datatype per l'uso di diversi derived data types
USE operators	! Modulo con overloading degli operatori



!******************************************************************************
!1) SUBROUTINE ADDNODE: per aggiungere un nodo
!******************************************************************************
INTERFACE addnode
	MODULE PROCEDURE AddNode_r_sgl
	MODULE PROCEDURE AddNode_r_dbl
	MODULE PROCEDURE AddNode_c_sgl
	MODULE PROCEDURE AddNode_c_dbl
	MODULE PROCEDURE AddNode_i_short
	MODULE PROCEDURE AddNode_i_long
	MODULE PROCEDURE AddNode_char
	MODULE PROCEDURE AddNode_rdbl_char
END INTERFACE


!******************************************************************************
!2) SUBROUTINE LLDESTRUCTOR: per deallocare una linked list
!******************************************************************************
INTERFACE lldestructor
	MODULE PROCEDURE lldest_r_sgl
	MODULE PROCEDURE lldest_r_dbl
	MODULE PROCEDURE lldest_c_sgl
	MODULE PROCEDURE lldest_c_dbl
	MODULE PROCEDURE lldest_i_short
	MODULE PROCEDURE lldest_i_long
	MODULE PROCEDURE lldest_char
	MODULE PROCEDURE lldest_rdbl_char
END INTERFACE


!******************************************************************************
!3) SUBROUTINE C_ADDNODE: per aggiungere un nodo ma solo se diverso
!******************************************************************************
INTERFACE c_addnode
	MODULE PROCEDURE C_AddNode_r_sgl
	MODULE PROCEDURE C_AddNode_r_dbl
	MODULE PROCEDURE C_AddNode_c_sgl
	MODULE PROCEDURE C_AddNode_c_dbl
	MODULE PROCEDURE C_AddNode_i_short
	MODULE PROCEDURE C_AddNode_i_long
	MODULE PROCEDURE C_AddNode_char
END INTERFACE

!******************************************************************************
!4) SUBROUTINE PATTERN: distanza euclidea tra due punti 
!******************************************************************************

!******************************************************************************
!5) FUNCTION DIST: distanza euclidea tra due punti 
!******************************************************************************
INTERFACE dist
	MODULE PROCEDURE dist_r_sgl
	MODULE PROCEDURE dist_r_dbl
END INTERFACE

!******************************************************************************
!6) SUBROUTINE CART_SPHER: da coordinate cartesiane a sferiche
!******************************************************************************
INTERFACE cart_spher
	MODULE PROCEDURE cart_spher_r_sgl
	MODULE PROCEDURE cart_spher_r_dbl1
	MODULE PROCEDURE cart_spher_r_sgl_v
	MODULE PROCEDURE cart_spher_r_dbl_v
END INTERFACE

CONTAINS

!******************************************************************************
!******************************************************************************
!******************************************************************************
!1) SUBROUTINE ADDNODE: per aggiungere un nodo
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE AddNode_r_sgl (head, input, error)		!REAL SINGLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_SGL_REAL), POINTER :: head		! Head list
REAL(sgl), INTENT(IN) :: input				! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_SGL_REAL), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if
			
current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_r_sgl

!-------------------------------------

SUBROUTINE AddNode_r_dbl (head, input, error)		!REAL DOUBLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_DBL_REAL), POINTER :: head		! Head list
REAL(dbl), INTENT(IN) :: input				! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_DBL_REAL), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_r_dbl

!-------------------------------------

SUBROUTINE AddNode_c_sgl (head, input, error)		!COMPLEX SINGLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_SGL_CMPLX), POINTER :: head		! Head list
COMPLEX(sgl), INTENT(IN) :: input			! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_SGL_CMPLX), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0	
		  END IF error_if
		  
current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_c_sgl

!-------------------------------------

SUBROUTINE AddNode_c_dbl (head, input, error)		!COMPLEX DOUBLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_DBL_CMPLX), POINTER :: head		! Head list
COMPLEX(dbl), INTENT(IN) :: input			! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_DBL_CMPLX), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_c_dbl

!-------------------------------------

SUBROUTINE AddNode_i_short (head, input, error)	!INT SHORT

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_SHORT_INT), POINTER :: head		! Head list
INTEGER(short), INTENT(IN) :: input			! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_SHORT_INT), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_i_short

!-------------------------------------

SUBROUTINE AddNode_i_long (head, input, error)		!INT LONG

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_LONG_INT), POINTER :: head		! Head list
INTEGER(lo), INTENT(IN) :: input			! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_LONG_INT), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_i_long

!-------------------------------------

SUBROUTINE AddNode_char (head, input, error)		! CHARACTER

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_CHAR), POINTER :: head			! Head list
CHARACTER(len=length), INTENT(IN) :: input		! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_CHAR), POINTER :: current		! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_char

!-------------------------------------

SUBROUTINE AddNode_rdbl_char (head, input, error)		! REAL SGL CHARACTER

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_RDBL_CHAR), POINTER :: head		! Head list
TYPE (RDBL_CHAR), INTENT(IN) :: input		! Valore di input
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interne alla procedura
TYPE (NODE_RDBL_CHAR), POINTER :: current	! Nuovo nodo nella lista
INTEGER(lo) :: status_all				! Status allocazione

! Inizia la subroutine vera e propria
ALLOCATE(current, STAT=status_all)			! Alloco il nuovo nodo

error_if :IF (status_all /= 0) THEN			! Controllo errore allocazione
			error=1
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

current%value=input							! Assegno il valore
current%next=>head							! Metto in coda il nodo
head=>current								! Metto head  al posto

END SUBROUTINE AddNode_rdbl_char

!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE LLDESTRUCTOR: per disallocare una linked list
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE lldest_r_sgl (head, error)				! REAL SINGLE PRECISION

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_SGL_REAL), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interna alla procedura
TYPE (NODE_SGL_REAL), POINTER :: cursor		! cursore nella lista

! Inizia la subroutine
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if
 
cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_r_sgl

!-------------------------------------------

SUBROUTINE lldest_r_dbl (head, error)				! REAL DOUBLE PRECISION

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_DBL_REAL), POINTER :: head		! Head list

! Variabile interna alla procedura
TYPE (NODE_DBL_REAL), POINTER :: cursor		! cursore nella lista
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_r_dbl
	
!-------------------------------------------

SUBROUTINE lldest_c_sgl (head, error)				! COMPLEX SINGLE PRECISION

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_SGL_CMPLX), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interna alla procedura
TYPE (NODE_SGL_CMPLX), POINTER :: cursor		! cursore nella lista

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_c_sgl
				
!-------------------------------------------

SUBROUTINE lldest_c_dbl (head, error)				! COMPLEX DOUBLE PRECISION

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_DBL_CMPLX), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interna alla procedura
TYPE (NODE_DBL_CMPLX), POINTER :: cursor		! cursore nella lista

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_c_dbl

!-------------------------------------------

SUBROUTINE lldest_i_short (head, error)				! INTEGER SHORT

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_SHORT_INT), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interna alla procedura
TYPE (NODE_SHORT_INT), POINTER :: cursor		! cursore nella lista

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_i_short

!-------------------------------------------

SUBROUTINE lldest_i_long (head, error)				! INTEGER LONG

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_LONG_INT), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interna alla procedura
TYPE (NODE_LONG_INT), POINTER :: cursor		! cursore nella lista

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_i_long

!-------------------------------------------

SUBROUTINE lldest_char (head, error)				! CHARACTER

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_CHAR), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error	! Errore

! Variabile interna alla procedura
TYPE (NODE_CHAR), POINTER :: cursor		! cursore nella lista

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_char

!-------------------------------------------

SUBROUTINE lldest_rdbl_char (head, error)				! REAL SINGLE CHARACTER

IMPLICIT NONE

! Dichiarazione dummy arguments
TYPE (NODE_RDBL_CHAR), POINTER :: head		! Head list
INTEGER(lo), INTENT(OUT) :: error		! Errore

! Variabile interna alla procedura
TYPE (NODE_RDBL_CHAR), POINTER :: cursor	! cursore nella lista

! Inizia la subroutine 
error_if: IF (.NOT. ASSOCIATED(head)) THEN	! Head deve essere associato
			error=1							! Se no errore
			RETURN
		  ELSE
		  	error=0
		  END IF error_if

cursor=>head

destruct_loop: DO
					IF (.NOT. ASSOCIATED(cursor)) EXIT	
					head=>cursor%next					! Nodo successivo
					DEALLOCATE(cursor)					! Cancello il nodo corrente
					cursor=>head						! Punto il nodo successivo da eliminare
			   END DO destruct_loop
			   
END SUBROUTINE lldest_rdbl_char


!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE C_ADDNODE: aggiunge un nodo solo se e' diverso dal resto
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE C_AddNode_r_sgl (head, input, error, ncount)		!REAL SINGLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_SGL_REAL), POINTER :: head					! Head list
REAL(sgl), INTENT(IN) :: input							! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_SGL_REAL), POINTER :: cursor=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_r_sgl

!-------------------------------------------------------------------------------

SUBROUTINE C_AddNode_r_dbl (head, input, error, ncount)		!REAL DOUBLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_DBL_REAL), POINTER :: head					! Head list
REAL(dbl), INTENT(IN) :: input							! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_DBL_REAL), POINTER :: cursor=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_r_dbl

!-------------------------------------------------------------------------------

SUBROUTINE C_AddNode_c_sgl (head, input, error, ncount)		!COMPLEX SINGLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_SGL_CMPLX), POINTER :: head					! Head list
COMPLEX(sgl), INTENT(IN) :: input							! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_SGL_CMPLX), POINTER :: cursor=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_c_sgl

!-------------------------------------------------------------------------------

SUBROUTINE C_AddNode_c_dbl (head, input, error, ncount)		!COMPLEX DOUBLE PRECISION

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_DBL_CMPLX), POINTER :: head					! Head list
COMPLEX(dbl), INTENT(IN) :: input						! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_DBL_CMPLX), POINTER :: cursor=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_c_dbl

!-------------------------------------------------------------------------------

SUBROUTINE C_AddNode_i_short (head, input, error, ncount)		!INTEGER SHORT

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_SHORT_INT), POINTER :: head					! Head list
INTEGER(short), INTENT(IN) :: input						! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_SHORT_INT), POINTER :: cursor=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_i_short

!-------------------------------------------------------------------------------

SUBROUTINE C_AddNode_i_long (head, input, error, ncount)		!INTEGER LONG

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_LONG_INT), POINTER :: head					! Head list
INTEGER(lo), INTENT(IN) :: input						! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_LONG_INT), POINTER :: cursor=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_i_long

!-------------------------------------------------------------------------------

SUBROUTINE C_AddNode_char (head, input, error, ncount)		!INTEGER LONG

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_CHAR), POINTER :: head						! Head list
CHARACTER(len=length), INTENT(IN) :: input					! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT), OPTIONAL :: ncount 			! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_CHAR), POINTER :: cursor=>NULL()				! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
option_if: IF (PRESENT(ncount)) THEN		! Procedura con ncount
			
			ass_if_n: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   				ncount=1						
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do_n: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if_n: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if_n: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			ncount=ncount+1
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if_n
		   								  END IF check_if_n
		   						END DO add_do_n
		   			END IF ass_if_n
			
			
		   ELSE								! Procedura senza ncount
		   
		   	ass_if: IF ( .NOT. ASSOCIATED(head)) THEN	!Primo valore
		   				CALL addnode(head,input,error)
		   				
		   				IF (error==1) THEN				!Errore in addnode
		   					RETURN
		   				END IF
		   				
		   			ELSE								!Altri valori
		   				cursor=>head
		   				add_do: DO						!Ciclo per aggiungere nodi
		   						IF (.NOT. ASSOCIATED(cursor)) EXIT
		   						
		   						check_if: IF (cursor%value == input) THEN			
		   										NULLIFY(cursor)
		   								    ELSE
		   										add_if: IF (.NOT. ASSOCIATED(cursor%next)) THEN
		   										  			CALL addnode(head,input,error)
		   										  			cursor=>cursor%next
		   												  ELSE
		   										  			cursor=>cursor%next	
		   												  END IF add_if
		   								  END IF check_if
		   						END DO add_do
		   			END IF ass_if
		   	
		   END IF option_if		

END SUBROUTINE C_AddNode_char

!-------------------------------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE PATTERN: costruisce il vettore per trovare e marcare le sfere 
!						uguali
!******************************************************************************
!******************************************************************************
!******************************************************************************


SUBROUTINE pattern (head1, input, head2, ncount,h_req, h_epseq, error)		

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_RDBL_CHAR), POINTER :: head1					! Head list per l'input
TYPE (NODE_LONG_INT), POINTER :: head2					! Head list per il pattern
TYPE (NODE_DBL_REAL), POINTER :: h_req					! Head list i raggi
TYPE (NODE_CHAR), POINTER :: h_epseq					! Head list i nomi
TYPE (RDBL_CHAR), INTENT(IN) :: input					! Valore di input
INTEGER(lo), INTENT(OUT) :: error					! Errore
INTEGER(lo), INTENT(OUT) :: ncount					! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_RDBL_CHAR), POINTER :: cursor1=>NULL()		! Nuovo nodo nella lista
TYPE (NODE_LONG_INT), POINTER :: cursor2=>NULL()		! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
ass_if_n: IF ( .NOT. ASSOCIATED(head1)) THEN	!Primo valore
	ncount=1
	CALL addnode(head1,input,error)
	CALL addnode(head2,ncount,error)
	CALL addnode(h_req,input%value,error)
	CALL addnode(h_epseq,input%text,error)

	IF (error==1) THEN						!Errore in addnode
		RETURN
	END IF

ELSE										!Altri valori

	cursor1=>head1
	cursor2=>head2
	add_do_n: DO							!Ciclo per aggiungere nodi

		IF (.NOT. ASSOCIATED(cursor1)) EXIT			!Se sono alla fine, esco dal loop

		check_if_n: IF (cursor1%value == input) THEN	!Se e' uguale, aggiorno la LL ausiliaria e il pattern

			CALL addnode(head2,cursor2%value,error)
			CALL addnode(head1,input,error)

			IF (error==1) THEN				!Errore in addnode
				RETURN
			END IF

			NULLIFY(cursor1)

		ELSE

			add_if_n: IF (.NOT. ASSOCIATED(cursor1%next)) THEN	!Se ho spazzolato tutta la LL,allora aggiungo il nodo

				ncount=ncount+1
				CALL addnode(head1,input,error)
				CALL addnode(head2,ncount,error)
				CALL addnode(h_req,input%value,error)
				CALL addnode(h_epseq,input%text,error)

				IF (error==1) THEN				!Errore in addnode
					RETURN
				END IF

				cursor1=>cursor1%next
				cursor2=>cursor2%next

			ELSE									!Se non sono alla fine,mando avanti il cursore

				cursor1=>cursor1%next
				cursor2=>cursor2%next

			END IF add_if_n
		END IF check_if_n
	END DO add_do_n
END IF ass_if_n

END SUBROUTINE pattern





!-------------------------------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE PATTERN: costruisce il vettore per trovare e marcare le sfere 
!						uguali
!******************************************************************************
!******************************************************************************
!******************************************************************************


SUBROUTINE pattern_shell(head1,head1_shell,input,input_shell,head2,ncount,h_req,h_epseq,h_req_shell,h_epseq_shell,error)

IMPLICIT NONE

! Dichiarazione dei dummy arguments
TYPE (NODE_RDBL_CHAR), POINTER :: head1					! Head list per l'input, contiene raggio e eps core
TYPE (NODE_RDBL_CHAR), POINTER :: head1_shell				! Head list per l'input, contiene raggio e eps shell
TYPE (NODE_LONG_INT), POINTER :: head2					! Head list per il pattern
TYPE (NODE_DBL_REAL), POINTER :: h_req					! Head list i raggi core
TYPE (NODE_CHAR), POINTER :: h_epseq					! Head list i nomi core
TYPE (NODE_DBL_REAL), POINTER :: h_req_shell				! Head list i raggi shell
TYPE (NODE_CHAR), POINTER :: h_epseq_shell				! Head list i nomi shell
TYPE (RDBL_CHAR), INTENT(IN) :: input,input_shell			! Valore di input
INTEGER(lo), INTENT(OUT) :: error						! Errore
INTEGER(lo), INTENT(OUT) :: ncount					! Lunghezza linked list

! Variabile interne alla procedura
TYPE (NODE_RDBL_CHAR), POINTER :: cursor1=>NULL()			! Nuovo nodo nella lista
TYPE (NODE_RDBL_CHAR), POINTER :: cursor1_shell=>NULL()		! Nuovo nodo nella lista che gira sulle shell
TYPE (NODE_LONG_INT), POINTER :: cursor2=>NULL()			! Nuovo nodo nella lista

! Inizia la subroutine vera e propria
ass_if_n: IF ( ( .NOT. ASSOCIATED(head1)) .AND. ( .NOT. ASSOCIATED(head1_shell)) ) THEN	!Primo valore
	ncount=1
	
	!Aggiungo il nodo per il core
	CALL addnode(head1,input,error)
	CALL addnode(h_req,input%value,error)
	CALL addnode(h_epseq,input%text,error)
	!Aggiungo il nodo per la shell
	CALL addnode(head1_shell,input,error)
	CALL addnode(h_req_shell,input_shell%value,error)
	CALL addnode(h_epseq_shell,input_shell%text,error)
	!Una volta aggiunti tutti i dati, aggiungo il nodo per la LL di pattern
	CALL addnode(head2,ncount,error)

	IF (error==1) THEN						!Errore in addnode
		RETURN
	END IF

ELSE										!Altri valori

	!Piazzo i cursori per i dati
	cursor1=>head1
	cursor1_shell=>head1_shell
	!Piazzo il cursore per il pattern
	cursor2=>head2
	
	add_do_n: DO							!Ciclo per aggiungere nodi

		IF ( (.NOT. ASSOCIATED(cursor1)) .AND. (.NOT. ASSOCIATED(cursor1_shell)) ) EXIT	!Se sono alla fine, esco dal loop

		check_if_n: IF ( (cursor1%value == input) .AND. (cursor1_shell%value == input_shell) ) THEN	!Se e' uguale, aggiorno la LL ausiliaria e il pattern

			CALL addnode(head2,cursor2%value,error)
			CALL addnode(head1,input,error)
			CALL addnode(head1_shell,input_shell,error)

			IF (error==1) THEN				!Errore in addnode
				RETURN
			END IF

			NULLIFY(cursor1)

		ELSE

			add_if_n: IF ( (.NOT. ASSOCIATED(cursor1%next)) .AND. (.NOT. ASSOCIATED(cursor1_shell%next)) ) THEN	!Se ho spazzolato tutta la LL,allora aggiungo il nodo

				ncount=ncount+1
				!Aggiungo il nodo per il core
				CALL addnode(head1,input,error)
				CALL addnode(h_req,input%value,error)
				CALL addnode(h_epseq,input%text,error)
				!Aggiungo il nodo per la shell
				CALL addnode(head1_shell,input,error)
				CALL addnode(h_req_shell,input_shell%value,error)
				CALL addnode(h_epseq_shell,input_shell%text,error)
				!Aggiungo il nodo per il pattern
				CALL addnode(head2,ncount,error)

				IF (error==1) THEN				!Errore in addnode
					RETURN
				END IF

				cursor1=>cursor1%next
				cursor1_shell=>cursor1_shell%next
				cursor2=>cursor2%next

			ELSE									!Se non sono alla fine,mando avanti il cursore

				cursor1=>cursor1%next
				cursor1_shell=>cursor1_shell%next
				cursor2=>cursor2%next

			END IF add_if_n
		END IF check_if_n
	END DO add_do_n
END IF ass_if_n

END SUBROUTINE pattern_shell






!-------------------------------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************
!5) SUBROUTINE DIST: calcola la distanza euclidea con due n-uple
!******************************************************************************
!******************************************************************************
!******************************************************************************
REAL(sgl) PURE FUNCTION dist_r_sgl(r1,r2)		

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(sgl) ,DIMENSION(:), INTENT(IN) :: r1,r2			! Vettori posizione	

! Inizio della subroutine vera e propria
dist_r_sgl=SQRT(SUM((r2-r1)**2))
			
END FUNCTION dist_r_sgl

!------------------------------------------------------------------------------

REAL(dbl) PURE FUNCTION dist_r_dbl(r1,r2)		

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(dbl) ,DIMENSION(:), INTENT(IN) :: r1,r2			! Vettori posizione	

! Inizio della subroutine vera e propria
dist_r_dbl=SQRT(SUM((r2-r1)**2))
			
END FUNCTION dist_r_dbl

!-------------------------------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************
!6) SUBROUTINE CART_SPHER: passa da coordinate cartesiane a sferiche
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE cart_spher_r_sgl_v (x,y,z,r,theta,phi,error)	! Vettori single precision		

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(sgl) ,DIMENSION(:), INTENT(IN) :: x,y,z			! Posizione cartesiana	
REAL(sgl) ,DIMENSION(:), INTENT(OUT) :: r,theta,phi		! Posizione Sferica													
INTEGER(lo), INTENT(OUT) :: error					! Errore

!Variabili interne alla procedura
INTEGER(lo) :: i,l1,l2,l3								! Indice e lunghezza vettore
REAL(sgl) :: r1											! r sul piano

! Inizio della subroutine vera e propria

! Prendo le lunghezze degli array
l1=SIZE(x)
l2=SIZE(y)
l3=SIZE(z)

size_if: IF (l1 /= l2 .OR. l2 /= l3 .OR. l3 /= l1) THEN !Check dimensioni
			error=1
			RETURN
		 ELSE
		 	error=0
		 END IF size_if


DO i=1,l1		 
r(i)=SQRT(x(i)**2+y(i)**2+z(i)**2)  ! Coordinata radiale

r_if: IF (r(i)==0.0D0) THEN				! Se raggio e' zero
						
		theta(i)=0.0D0
		phi(i)=0.0D0
		
ELSE
		theta(i)=ACOS(z(i)/r(i))
		r1=SQRT(x(i)**2+y(i)**2)
		
		r1_if: IF (r1==0.0D0) THEN 		! Se siamo sull'asse delle zeta
				
				phi(i)=0.0D0
				
			   ELSE
			   
			   	x_if: IF (x(i)<0.0D0) THEN
						phi(i)=Pi-ASIN(y(i)/r1)
					  ELSE
					  	phi(i)=ASIN(y(i)/r1)
				END IF x_if
				
		 END IF r1_if
		
END IF r_if

END DO
			
END SUBROUTINE cart_spher_r_sgl_v


!------------------------------------------------------------------------------

SUBROUTINE cart_spher_r_dbl_v (x,y,z,r,theta,phi,error)	! Vettori double precision	

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(dbl) ,DIMENSION(:), INTENT(IN) :: x,y,z			! Posizione cartesiana	
REAL(dbl) ,DIMENSION(:), INTENT(OUT) :: r,theta,phi		! Posizione Sferica													
INTEGER(lo), INTENT(OUT) :: error					! Errore

!Variabili interne alla procedura
INTEGER(lo) :: i,l1,l2,l3								! Indice e lunghezza vettore
REAL(dbl) :: r1											! r sul piano

! Inizio della subroutine vera e propria

! Prendo le lunghezze degli array
l1=SIZE(x)
l2=SIZE(y)
l3=SIZE(z)

size_if: IF (l1 /= l2 .OR. l2 /= l3 .OR. l3 /= l1) THEN !Check dimensioni
			error=1
			RETURN
		 ELSE
		 	error=0
		 END IF size_if


DO i=1,l1		 
r(i)=SQRT(x(i)**2+y(i)**2+z(i)**2)  ! Coordinata radiale

r_if: IF (r(i)==0.0D0) THEN				! Se raggio e' zero
						
		theta(i)=0.0D0
		phi(i)=0.0D0
		
ELSE
		theta(i)=ACOS(z(i)/r(i))
		r1=SQRT(x(i)**2+y(i)**2)
		
		r1_if: IF (r1==0.0D0) THEN 		! Se siamo sull'asse delle zeta
				
				phi(i)=0.0D0
				
			   ELSE
			   
			   	x_if: IF (x(i)<0.0D0) THEN
						phi(i)=Pi_D-ASIN(y(i)/r1)
					  ELSE
					  	phi(i)=ASIN(y(i)/r1)
				END IF x_if
				
		 END IF r1_if
		
END IF r_if

END DO
			
END SUBROUTINE cart_spher_r_dbl_v


!------------------------------------------------------------------------------

SUBROUTINE cart_spher_r_sgl (x,y,z,r,theta,phi,error)	! Scalari single precision	

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(sgl) ,INTENT(IN) :: x,y,z							! Posizione cartesiana	
REAL(sgl) ,INTENT(OUT) :: r,theta,phi					! Posizione Sferica													
INTEGER(lo), INTENT(OUT) :: error					! Errore

!Variabili interne alla procedura
REAL(sgl) :: r1											! r sul piano

! Inizio della subroutine vera e propria		 
error=0					! Lo assegno a priori

r=SQRT(x**2+y**2+z**2)  ! Coordinata radiale

r_if: IF (r==0.0D0) THEN				! Se raggio e' zero
						
		theta=0.0D0
		phi=0.0D0
		
ELSE
		theta=ACOS(z/r)
		r1=SQRT(x**2+y**2)
		
		r1_if: IF (r1==0.0D0) THEN 		! Se siamo sull'asse delle zeta
				
				phi=0.0D0
				
			   ELSE
			   
			   	x_if: IF (x<0.0D0) THEN
						phi=Pi-ASIN(y/r1)
					  ELSE
					  	phi=ASIN(y/r1)
				END IF x_if
				
		 END IF r1_if
		
END IF r_if
			
END SUBROUTINE cart_spher_r_sgl

!------------------------------------------------------------------------------

SUBROUTINE cart_spher_r_dbl (x,y,z,r,theta,phi,error)	! Scalari double precision	

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(dbl) ,INTENT(IN) :: x,y,z							! Posizione cartesiana	
REAL(dbl) ,INTENT(OUT) :: r,theta,phi					! Posizione Sferica													
INTEGER(lo), INTENT(OUT) :: error					! Errore

!Variabili interne alla procedura
REAL(dbl) :: r1											! r sul piano

! Inizio della subroutine vera e propria
error=0					! Lo assegno a priori

r=SQRT(x**2+y**2+z**2)  ! Coordinata radiale

r_if: IF (r==0.0D0) THEN				! Se raggio e' zero
						
		theta=0.0D0
		phi=0.0D0
		
ELSE
		theta=ACOS(z/r)
		r1=SQRT(x**2+y**2)
		
		r1_if: IF (r1==0.0D0) THEN 		! Se siamo sull'asse delle zeta
				
				phi=0.0D0
				
			   ELSE
			   
			   	x_if: IF (x<0.0D0) THEN
						phi=Pi_D-ASIN(y/r1)
					  ELSE
					  	phi=ASIN(y/r1)
				END IF x_if
				
		 END IF r1_if
		
END IF r_if
			
END SUBROUTINE cart_spher_r_dbl


!------------------------------------------------------------------------------

SUBROUTINE cart_spher_r_dbl1 (x,y,z,r,theta,phi,error)	! Scalari double precision	

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(dbl) ,INTENT(IN) :: x,y,z							! Posizione cartesiana	
REAL(dbl) ,INTENT(OUT) :: r,theta,phi					! Posizione Sferica													
INTEGER(lo), INTENT(OUT) :: error					! Errore

!Variabili interne alla procedura
REAL(dbl) :: r1											! r sul piano

! Inizio della subroutine vera e propria
error=0					! Lo assegno a priori

r=SQRT(x**2+y**2+z**2)  ! Coordinata radiale

r_if: IF (r==0.0D0) THEN				! Se raggio e' zero
						
		theta=0.0D0
		phi=0.0D0
		
ELSE
		theta=ACOS(z/r)
		r1=SQRT(x**2+y**2)
		
		r1_if: IF (r1==0.0D0) THEN 		! Se siamo sull'asse delle zeta
				
				phi=0.0D0
				
			   ELSE
			   
			   	phi=ATAN2(y,x)
				
		 END IF r1_if
		
END IF r_if
			
END SUBROUTINE cart_spher_r_dbl1



!------------------------------------------------------------------------------

SUBROUTINE spher_cart_r_dbl (r,theta,phi,x,y,z,error)	! Scalari double precision	

IMPLICIT NONE

! Dichiarazione dei dummy arguments
REAL(dbl) ,INTENT(OUT) :: x,y,z							! Posizione cartesiana	
REAL(dbl) ,INTENT(IN) :: r,theta,phi					! Posizione Sferica													
INTEGER(lo), INTENT(OUT) :: error						! Errore

! Inizio della subroutine vera e propria
error=0			

x=r*SIN(theta)*COS(phi)
y=r*SIN(theta)*SIN(phi)
z=r*COS(theta)
			
END SUBROUTINE spher_cart_r_dbl


!$$$$$$ !-------------------------------------------------------------------------------
!$$$$$$ !******************************************************************************
!$$$$$$ !******************************************************************************
!$$$$$$ !******************************************************************************
!$$$$$$ !7) FUNCTION LNF: logaritmo del fattoriale in doppia precisione
!$$$$$$ !******************************************************************************
!$$$$$$ !******************************************************************************
!$$$$$$ !******************************************************************************
!$$$$$$ FUNCTION lnf(z,ierr)
!$$$$$$ 
!$$$$$$ IMPLICIT NONE
!$$$$$$ 
!$$$$$$ ! Dichiarazione funzione
!$$$$$$ REAL(dbl) :: lnf
!$$$$$$ 
!$$$$$$ ! Argomento Funzione
!$$$$$$ REAL(dbl), INTENT(IN) :: z
!$$$$$$ INTEGER(lo) :: ierr
!$$$$$$ 
!$$$$$$ ! Variabili interne
!$$$$$$ INTEGER(lo) ::intz
!$$$$$$ 
!$$$$$$ 
!$$$$$$ 
!$$$$$$ !Funzione
!$$$$$$ intz=INT(Z)
!$$$$$$ 
!$$$$$$ IF (intz==0) THEN 
!$$$$$$     lnf=0.0D0
!$$$$$$ ELSE 
!$$$$$$     lnf=LOG(z)+dgamln(z,ierr)
!$$$$$$ END IF
!$$$$$$ !lnf=dgamln(1.0D0+z,ierr)
!$$$$$$ 
!$$$$$$ END FUNCTION lnf


!-------------------------------------------------------------------------------
!******************************************************************************
!******************************************************************************
!******************************************************************************
!7bis) FUNCTION LNF1: logaritmo del fattoriale in doppia precisione
!******************************************************************************
!******************************************************************************
!******************************************************************************
FUNCTION lnf(z,ierr)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: lnf

! Parametri e argomento funzione
REAL(dbl), PARAMETER :: v_c0(1:11) = (/  0.16427423239836267D5, -0.48589401600331902D5,  &   !Parametri dell'espansione in serie.
					      &        0.55557391003815523D5, -0.30964901015912058D5,  &
					      &        0.87287202992571788D4, -0.11714474574532352D4,  &
					      &        0.63103078123601037D2, -0.93060589791758878D0,  &
					      &        0.13919002438227877D-2,-0.45006835613027859D-8, & 
					      &        0.13069587914063262d-9   /)
REAL(dbl), INTENT(IN) :: z
INTEGER(lo), INTENT(INOUT) :: ierr

! Variabili interne
REAL(dbl) :: a,b,cp,z1 
INTEGER(lo) :: i

!Funzione
ierr=0
z1=z
a=1.0D0
cp=2.5066282746310005D0
b=z1+10.5D0
b=(z1+0.5D0)*LOG(b)-b

log_do: DO i=1,11

	z1=z1+1.0D0
    a=a+v_c0(i)/z1

END DO log_do

lnf=b+LOG(cp*a)

END FUNCTION lnf



!******************************************************************************
!******************************************************************************
!******************************************************************************
!7) FUNCTION LPOCH: logaritmo della funzione di pochhammer
!******************************************************************************
!******************************************************************************
!******************************************************************************
FUNCTION lpoch(x,n,ierr)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: lpoch

! Argomento Funzione
REAL(dbl), INTENT(IN) :: x,n
INTEGER(lo) :: ierr

! Variabili interne
REAL(dbl) :: somma

!Funzione
IF (ABS(n)<10**(-5)) THEN
	lpoch=1.0D0
ELSE
	somma=x+n
	lpoch=lnf(somma-1.0D0,ierr)-lnf(x-1.0D0,ierr)
END IF

END FUNCTION lpoch



!******************************************************************************
!******************************************************************************
!******************************************************************************
!8) FUNCTION FATT: FATTORIALE
!******************************************************************************
!******************************************************************************
!******************************************************************************
FUNCTION fatt(n)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: fatt

! Argomento Funzione
REAL(dbl), INTENT(IN) :: n

! Variabili interne
INTEGER(lo) :: ni,i
REAL(dbl) :: fact

!Funzione
ni=INT(n,lo)

f_if: IF (n<0.0D0) THEN !Flag di errore
	
	fatt=-1.0D0
	
ELSE IF (ni==0) THEN	!Caso n=0
	
	fatt=1.0D0
	
ELSE 					!Caso Generale

	fatt=1.0D0
	
	DO i=1,ni
		fatt=i*fatt
	END DO

END IF f_if

END FUNCTION fatt




!******************************************************************************
!******************************************************************************
!******************************************************************************
!9) FUNCTION FACT2: DOPPIO FATTORIALE
!******************************************************************************
!******************************************************************************
!******************************************************************************
FUNCTION fact2(n)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: fact2

! Argomento Funzione
REAL(dbl), INTENT(IN) :: n

! Variabili interne
INTEGER(lo) :: ni,i
REAL(dbl) :: fact

!Funzione
ni=INT(n,lo)

f_if: IF (n<-1.0D0) THEN
	
	fact2=-1.0D0
	
ELSE IF ((ni==-1) .OR. (ni==0)) THEN
	
	fact2=1.0D0
	
ELSE IF (MOD(ni,2)==0) THEN

	fact2=1.0D0
	
	DO i=2,ni,2
		fact2=i*fact2
	END DO

ELSE IF (MOD(ni,2)/=0) THEN

	fact2=1.0D0
	
	DO i=1,ni,2
		fact2=i*fact2
	END DO

END IF f_if

END FUNCTION fact2


!******************************************************************************
!******************************************************************************
!******************************************************************************
!9) Subroutine interp(x0,x1,x2,y1,y2)
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE interp(x0,x1,x2,y1,y2,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl), INTENT(IN) :: x0,x1,x2,y1	!Ordinata da calcolare,limiti ordinate,limite ascisse
REAL(dbl), INTENT(INOUT) :: y2			!Limite ordinate e output
INTEGER(lo), INTENT(OUT) :: error		!Flag di errore

! Dichiarazione variabili interne
REAL(dbl) :: mr,qr						!Pendenza ed intercetta

! Subroutine vera e propria
error=0

!Controllo gli argomenti
err_if: IF ((x1==x2) .OR. (x1>x2))THEN 
	error=1
	WRITE(*,10)
	10 FORMAT ("x1=x2 o x1>x2, la subroutine si ferma")
	RETURN
END IF err_if

!Calcolo i parametri e l'output
mr=(y2-y1)/(x2-x1)
qr=y1-mr*x1
y2=mr*x0+qr

END SUBROUTINE interp


!******************************************************************************
!******************************************************************************
!******************************************************************************
!10) Subroutine rot(alpha,beta,gamma,v_xyz,v_xyz1,error)
!******************************************************************************
!******************************************************************************
!******************************************************************************
!SUBROUTINE rot(alpha,beta,gamma,v_xyz,v_xyz1,error)
!
!IMPLICIT NONE
!
! Dichiarazione dummy arguments
!REAL(dbl), INTENT(IN) :: alpha,beta,gamma			!Angoli Eulero
!REAL(dbl), DIMENSION(:), INTENT(IN) :: v_xyz		!Coord in
!REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_xyz1		!Coord out
!INTEGER(lo), INTENT(OUT) :: error					!Flag di errore
!
! Dichiarazione variabili interne
!REAL(dbl), DIMENSION(3,3) :: m_rot					!Matrice Rotazione
!
! Subroutine vera e propria
!error=0
!
!m_rot(1,1)=COS(alpha)*COS(beta)*COS(gamma)-SIN(alpha)*SIN(gamma)
!m_rot(1,2)=COS(beta)*COS(gamma)*SIN(alpha)+COS(alpha)*SIN(gamma)
!m_rot(1,3)=-COS(gamma)*SIN(beta)
!m_rot(2,1)=-COS(gamma)*SIN(alpha)-COS(alpha)*COS(beta)*SIN(gamma)
!m_rot(2,2)=COS(alpha)*COS(gamma)-COS(beta)*SIN(alpha)*SIN(gamma)
!m_rot(2,3)=SIN(beta)*SIN(gamma)
!m_rot(3,1)=COS(alpha)*SIN(beta)
!m_rot(3,2)=SIN(alpha)*SIN(beta)
!m_rot(3,3)=COS(beta)
!
!CALL dgemv ('N',3,3,1.0D0,m_rot,3,v_xyz,1,0.0D0,v_xyz1,1)
!
!END SUBROUTINE rot

!******************************************************************************
!******************************************************************************
!******************************************************************************
!10) Subroutine weigths(nalpha,nbeta,v_walpha,v_xalpha,v_wbeta,v_xbeta,error)
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE weigths(nalpha,nbeta,v_walpha,v_xalpha,v_wbeta,v_xbeta,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nalpha,nbeta								!Numero punti integrazione
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_walpha,v_xalpha				!Vettori peso e coord per alpha
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_wbeta,v_xbeta					!Vettori peso e coord per beta
INTEGER(lo), INTENT(OUT) :: error										!Flag di errore

! Dichiarazione variabili interne
REAL(dbl), DIMENSION(20,20) :: m_walpha,m_xalpha,m_wbeta,m_xbeta		!Matrice Rotazione

! Subroutine vera e propria

!Controllo l'errore della subroutine
weigths_if : IF (( nalpha>20 ) .OR. ( nbeta>20 )) THEN
                WRITE(*,*)
                WRITE(*,10)
                10 FORMAT ("Steps dell'integrazione out of bounds")
                RETURN
END IF weigths_if



error=0
m_walpha=0.0D0
m_xalpha=0.0D0
m_wbeta=0.0D0
m_xbeta=0.0D0
!---------------------------------------------
!Matrice per alpha, dove la funzione peso e' 1.
!---------------------------------------------
m_xalpha(1,1)=3.1415926535897932D0
m_xalpha(2,1)=1.3277932893555754D0
m_xalpha(2,2)=4.9553920178240111D0
m_xalpha(3,1)=7.0812544800562610D-1
m_xalpha(3,2)=3.1415926535897932D0
m_xalpha(3,3)=5.5750598591739604D0
m_xalpha(4,1)=4.3625314334650658D-1
m_xalpha(4,2)=2.0735107047038175D0
m_xalpha(4,3)=4.2096746024757690D0
m_xalpha(4,4)=5.8469321638330799D0
m_xalpha(5,1)=2.9474470675775580D-1
m_xalpha(5,2)=1.4499414247782150D0
m_xalpha(5,3)=3.1415926535897932D0
m_xalpha(5,4)=4.8332438824013714D0
m_xalpha(5,5)=5.9884406004218307D0
m_xalpha(6,1)=2.1215327807272746D-1
m_xalpha(6,2)=1.0643421025827622D0
m_xalpha(6,3)=2.3919483715852460D0
m_xalpha(6,4)=3.8912369355943405D0
m_xalpha(6,5)=5.2188432045968243D0
m_xalpha(6,6)=6.0710320291068590D0
m_xalpha(7,1)=1.5988220870983761D-1
m_xalpha(7,2)=8.1200372850300619D-1
m_xalpha(7,3)=1.8665925075275247D0
m_xalpha(7,4)=3.1415926535897932D0
m_xalpha(7,5)=4.4165927996520617D0
m_xalpha(7,6)=5.4711815786765803D0
m_xalpha(7,7)=6.1233030984697489D0
m_xalpha(8,1)=1.2475309510033664D-1
m_xalpha(8,2)=6.3879110078588453D-1
m_xalpha(8,3)=1.4905838953733143D0
m_xalpha(8,4)=2.5653157283115897D0
m_xalpha(8,5)=3.7178695788679968D0
m_xalpha(8,6)=4.7926014118062722D0
m_xalpha(8,7)=5.6443942063937019D0
m_xalpha(8,8)=6.1584322120792498D0
m_xalpha(9,1)=1.0002755765490041D-1
m_xalpha(9,2)=5.1512346863989426D-1
m_xalpha(9,3)=1.2146294666957722D0
m_xalpha(9,4)=2.1229204807230464D0
m_xalpha(9,5)=3.1415926535897932D0
m_xalpha(9,6)=4.1602648264565401D0
m_xalpha(9,7)=5.0685558404838143D0
m_xalpha(9,8)=5.7680618385396922D0
m_xalpha(9,9)=6.1831577495246861D0
m_xalpha(10,1)=8.1975058317108093D-2
m_xalpha(10,2)=4.2391593591002604D-1
m_xalpha(10,3)=1.0071645450429653D0
m_xalpha(10,4)=1.7800408672936973D0
m_xalpha(10,5)=2.6738901239370640D0
m_xalpha(10,6)=3.6092951832425224D0
m_xalpha(10,7)=4.5031444398858892D0
m_xalpha(10,8)=5.2760207621366212D0
m_xalpha(10,9)=5.8592693712695604D0
m_xalpha(10,10)=6.2012102488624784D0
m_xalpha(11,1)=6.8396687627139341D-2
m_xalpha(11,2)=3.5480330688408202D-1
m_xalpha(11,3)=8.4775247687450612D-1
m_xalpha(11,4)=1.5108040675667751D0
m_xalpha(11,5)=2.2947978550244983D0
m_xalpha(11,6)=3.1415926535897932D0
m_xalpha(11,7)=3.9883874521550882D0
m_xalpha(11,8)=4.7723812396128113D0
m_xalpha(11,9)=5.4354328303050804D0
m_xalpha(11,10)=5.9283820002955045D0
m_xalpha(11,11)=6.2147886195524471D0
m_xalpha(12,1)=5.7928975987362026D-2
m_xalpha(12,2)=3.0122452299254974D-1
m_xalpha(12,3)=7.2287206836182955D-1
m_xalpha(12,4)=1.2964788830815699D0
m_xalpha(12,5)=1.9860159185781890D0
m_xalpha(12,6)=2.7481602974261530D0
m_xalpha(12,7)=3.5350250097534335D0
m_xalpha(12,8)=4.2971693886013975D0
m_xalpha(12,9)=4.9867064240980166D0
m_xalpha(12,10)=5.5603132388177569D0
m_xalpha(12,11)=5.9819607841870367D0
m_xalpha(12,12)=6.2252563311922245D0
m_xalpha(13,1)=4.9690399098315214D-2
m_xalpha(13,2)=2.5887226364513142D-1
m_xalpha(13,3)=6.2336081246349410D-1
m_xalpha(13,4)=1.1235926877657640D0
m_xalpha(13,5)=1.7326111217454157D0
m_xalpha(13,6)=2.4175865012264663D0
m_xalpha(13,7)=3.1415926535897932D0
m_xalpha(13,8)=3.8655988059531201D0
m_xalpha(13,9)=4.5505741854341708D0
m_xalpha(13,10)=5.1595926194138225D0
m_xalpha(13,11)=5.6598244947160924D0
m_xalpha(13,12)=6.0243130435344551D0
m_xalpha(13,13)=6.2334949080812713D0
m_xalpha(14,1)=4.3090685833326569D-2
m_xalpha(14,2)=2.2482844373581634D-1
m_xalpha(14,3)=5.4286307912680361D-1
m_xalpha(14,4)=9.8239831296901311D-1
m_xalpha(14,5)=1.5228913228348575D0
m_xalpha(14,6)=2.1390715796962990D0
m_xalpha(14,7)=2.8021280205467805D0
m_xalpha(14,8)=3.4810572866328060D0
m_xalpha(14,9)=4.1441137274832875D0
m_xalpha(14,10)=4.7602939843447290D0
m_xalpha(14,11)=5.3007869942105734D0
m_xalpha(14,12)=5.7403222280527829D0
m_xalpha(14,13)=6.0583568634437701D0
m_xalpha(14,14)=6.2400946213462599D0
m_xalpha(15,1)=3.7722617174954806D-2
m_xalpha(15,2)=1.9706144961855203D-1
m_xalpha(15,3)=4.7687308242109690D-1
m_xalpha(15,4)=8.6576723061849864D-1
m_xalpha(15,5)=1.3478306707186042D0
m_xalpha(15,6)=1.9033296772083993D0
m_xalpha(15,7)=2.5095227659417986D0
m_xalpha(15,8)=3.1415926535897932D0
m_xalpha(15,9)=3.7736625412377879D0
m_xalpha(15,10)=4.3798556299711871D0
m_xalpha(15,11)=4.9353546364609823D0
m_xalpha(15,12)=5.4174180765610878D0
m_xalpha(15,13)=5.8063122247584896D0
m_xalpha(15,14)=6.0861238575610344D0
m_xalpha(15,15)=6.2454626900046317D0
m_xalpha(16,1)=3.3297944765153212D-2
m_xalpha(16,2)=1.7412270033851633D-1
m_xalpha(16,3)=4.2213202745008155D-1
m_xalpha(16,4)=7.6841971381237148D-1
m_xalpha(16,5)=1.2004771833467960D0
m_xalpha(16,6)=1.7026905096809779D0
m_xalpha(16,7)=2.2569090072368731D0
m_xalpha(16,8)=2.8431020506847435D0
m_xalpha(16,9)=3.4400832564948430D0
m_xalpha(16,10)=4.0262762999427133D0
m_xalpha(16,11)=4.5804947974986085D0
m_xalpha(16,12)=5.0827081238327904D0
m_xalpha(16,13)=5.5147655933672150D0
m_xalpha(16,14)=5.8610532797295049D0
m_xalpha(16,15)=6.1090626068410701D0
m_xalpha(16,16)=6.2498873624144333D0
m_xalpha(17,1)=2.9608017515802154D-2
m_xalpha(17,2)=1.5495741845338888D-1
m_xalpha(17,3)=3.7623979483899766D-1
m_xalpha(17,4)=6.8639400027005689D-1
m_xalpha(17,5)=1.0754577713167543D0
m_xalpha(17,6)=1.5309278287139118D0
m_xalpha(17,7)=2.0381655258157074D0
m_xalpha(17,8)=2.5808680602204503D0
m_xalpha(17,9)=3.1415926535897932D0
m_xalpha(17,10)=3.7023172469591362D0
m_xalpha(17,11)=4.2450197813638791D0
m_xalpha(17,12)=4.7522574784656747D0
m_xalpha(17,13)=5.2077275358628322D0
m_xalpha(17,14)=5.5967913069095296D0
m_xalpha(17,15)=5.9069455123405888D0
m_xalpha(17,16)=6.1282278887261976D0
m_xalpha(17,17)=6.2535772896637843D0
m_xalpha(18,1)=2.6498804923070533D-2
m_xalpha(18,2)=1.3878315549110905D-1
m_xalpha(18,3)=3.3739930226494259D-1
m_xalpha(18,4)=6.1667905882802844D-1
m_xalpha(18,5)=9.6859372052814062D-1
m_xalpha(18,6)=1.3830207229940265D0
m_xalpha(18,7)=1.8480382296310620D0
m_xalpha(18,8)=2.3502687374168989D0
m_xalpha(18,9)=2.8752640954098987D0
m_xalpha(18,10)=3.4079212117696878D0
m_xalpha(18,11)=3.9329165697626876D0
m_xalpha(18,12)=4.4351470775485245D0
m_xalpha(18,13)=4.9001645841855600D0
m_xalpha(18,14)=5.3145915866514459D0
m_xalpha(18,15)=5.6665062483515580D0
m_xalpha(18,16)=5.9457860049146439D0
m_xalpha(18,17)=6.1444021516884774D0
m_xalpha(18,18)=6.2566865022565159D0
m_xalpha(19,1)=2.3854603598555350D-2
m_xalpha(19,2)=1.2500977692598067D-1
m_xalpha(19,3)=3.0424470174722993D-1
m_xalpha(19,4)=5.5695833261205538D-1
m_xalpha(19,5)=8.7661060738672053D-1
m_xalpha(19,6)=1.2549239363168119D0
m_xalpha(19,7)=1.6821006254103105D0
m_xalpha(19,8)=2.1470772027537888D0
m_xalpha(19,9)=2.6378111105068523D0
m_xalpha(19,10)=3.1415926535897932D0
m_xalpha(19,11)=3.6453741966727342D0
m_xalpha(19,12)=4.1361081044257977D0
m_xalpha(19,13)=4.6010846817692759D0
m_xalpha(19,14)=5.0282613708627746D0
m_xalpha(19,15)=5.4065746997928659D0
m_xalpha(19,16)=5.7262269745675311D0
m_xalpha(19,17)=5.9789406054323565D0
m_xalpha(19,18)=6.1581755302536058D0
m_xalpha(19,19)=6.2593307035810311D0
m_xalpha(20,1)=2.1587142319976703D-2
m_xalpha(20,2)=1.1318552858670486D-1
m_xalpha(20,3)=2.7572367544374244D-1
m_xalpha(20,4)=5.0542893941059704D-1
m_xalpha(20,5)=7.9692181911491876D-1
m_xalpha(20,6)=1.1433710829306258D0
m_xalpha(20,7)=1.5366566332996323D0
m_xalpha(20,8)=1.9675603506796556D0
m_xalpha(20,9)=2.4259822970515028D0
m_xalpha(20,10)=2.9011774969920140D0
m_xalpha(20,11)=3.3820078101875725D0
m_xalpha(20,12)=3.8572030101280836D0
m_xalpha(20,13)=4.3156249564999309D0
m_xalpha(20,14)=4.7465286738799542D0
m_xalpha(20,15)=5.1398142242489607D0
m_xalpha(20,16)=5.4862634880646677D0
m_xalpha(20,17)=5.7777563677689894D0
m_xalpha(20,18)=6.0074616317358440D0
m_xalpha(20,19)=6.1699997785928816D0
m_xalpha(20,20)=6.2615981648596098D0


m_walpha(1,1)=6.2831853071795865D0
m_walpha(2,1)=3.1415926535897932D0
m_walpha(2,2)=3.1415926535897932D0
m_walpha(3,1)=1.7453292519943296D0
m_walpha(3,2)=2.7925268031909273D0
m_walpha(3,3)=1.7453292519943296D0
m_walpha(4,1)=1.0928182259994402D0
m_walpha(4,2)=2.0487744275903530D0
m_walpha(4,3)=2.0487744275903530D0
m_walpha(4,4)=1.0928182259994402D0
m_walpha(5,1)=7.4432776153043700D-1
m_walpha(5,2)=1.5036563150382595D0
m_walpha(5,3)=1.7872171540421935D0
m_walpha(5,4)=1.5036563150382595D0
m_walpha(5,5)=7.4432776153043700D-1
m_walpha(6,1)=5.3823176663840207D-1
m_walpha(6,2)=1.1333659075855298D0
m_walpha(6,3)=1.4699949793658614D0
m_walpha(6,4)=1.4699949793658614D0
m_walpha(6,5)=1.1333659075855298D0
m_walpha(6,6)=5.3823176663840207D-1
m_walpha(7,1)=4.0678901846644394D-1
m_walpha(7,2)=8.7872040307216866D-1
m_walpha(7,3)=1.1995544815867014D0
m_walpha(7,4)=1.3130575009289585D0
m_walpha(7,5)=1.1995544815867014D0
m_walpha(7,6)=8.7872040307216866D-1
m_walpha(7,7)=4.0678901846644394D-1
m_walpha(8,1)=3.1801882594349384D-1
m_walpha(8,2)=6.9863062413641994D-1
m_walpha(8,3)=9.8553849407226550D-1
m_walpha(8,4)=1.1394047094376140D0
m_walpha(8,5)=1.1394047094376140D0
m_walpha(8,6)=9.8553849407226550D-1
m_walpha(8,7)=6.9863062413641994D-1
m_walpha(8,8)=3.1801882594349384D-1
m_walpha(9,1)=2.5533102140172596D-1
m_walpha(9,2)=5.6752293452347246D-1
m_walpha(9,3)=8.1873264926638200D-1
m_walpha(9,4)=9.8126728259911810D-1
m_walpha(9,5)=1.0374775315981894D0
m_walpha(9,6)=9.8126728259911810D-1
m_walpha(9,7)=8.1873264926638200D-1
m_walpha(9,8)=5.6752293452347246D-1
m_walpha(9,9)=2.5533102140172596D-1
m_walpha(10,1)=2.0945420548513033D-1
m_walpha(10,2)=4.6951526056054718D-1
m_walpha(10,3)=6.8828010698191944D-1
m_walpha(10,4)=8.4592634724050947D-1
m_walpha(10,5)=9.2841673332168683D-1
m_walpha(10,6)=9.2841673332168683D-1
m_walpha(10,7)=8.4592634724050947D-1
m_walpha(10,8)=6.8828010698191944D-1
m_walpha(10,9)=4.6951526056054718D-1
m_walpha(10,10)=2.0945420548513033D-1
m_walpha(11,1)=1.7488796148804153D-1
m_walpha(11,2)=3.9452236614603636D-1
m_walpha(11,3)=5.8524795808626295D-1
m_walpha(11,4)=7.3259981770514494D-1
m_walpha(11,5)=8.2562482636340274D-1
m_walpha(11,6)=8.5741944760180944D-1
m_walpha(11,7)=8.2562482636340274D-1
m_walpha(11,8)=7.3259981770514494D-1
m_walpha(11,9)=5.8524795808626295D-1
m_walpha(11,10)=3.9452236614603636D-1
m_walpha(11,11)=1.7488796148804153D-1
m_walpha(12,1)=1.4820569022249282D-1
m_walpha(12,2)=3.3595980092673639D-1
m_walpha(12,3)=5.0290090095070981D-1
m_walpha(12,4)=6.3826929524192654D-1
m_walpha(12,5)=7.3353843745694184D-1
m_walpha(12,6)=7.8271852879098584D-1
m_walpha(12,7)=7.8271852879098584D-1
m_walpha(12,8)=7.3353843745694184D-1
m_walpha(12,9)=6.3826929524192654D-1
m_walpha(12,10)=5.0290090095070981D-1
m_walpha(12,11)=3.3595980092673639D-1
m_walpha(12,12)=1.4820569022249282D-1
m_walpha(13,1)=1.2718425195861055D-1
m_walpha(13,2)=2.8940822712788102D-1
m_walpha(13,3)=4.3628399948471066D-1
m_walpha(13,4)=5.5966210442827737D-1
m_walpha(13,5)=6.5287336823995617D-1
m_walpha(13,6)=7.1088957674485286D-1
m_walpha(13,7)=7.3058225121100922D-1
m_walpha(13,8)=7.1088957674485286D-1
m_walpha(13,9)=6.5287336823995617D-1
m_walpha(13,10)=5.5966210442827737D-1
m_walpha(13,11)=4.3628399948471066D-1
m_walpha(13,12)=2.8940822712788102D-1
m_walpha(13,13)=1.2718425195861055D-1
m_walpha(14,1)=1.1033103857626982D-1
m_walpha(14,2)=2.5182405774691301D-1
m_walpha(14,3)=3.8176184894784863D-1
m_walpha(14,4)=4.9386831506522906D-1
m_walpha(14,5)=5.8288606647551246D-1
m_walpha(14,6)=6.4464998615473398D-1
m_walpha(14,7)=6.7627134062328629D-1
m_walpha(14,8)=6.7627134062328629D-1
m_walpha(14,9)=6.4464998615473398D-1
m_walpha(14,10)=5.8288606647551246D-1
m_walpha(14,11)=4.9386831506522906D-1
m_walpha(14,12)=3.8176184894784863D-1
m_walpha(14,13)=2.5182405774691301D-1
m_walpha(14,14)=1.1033103857626982D-1
m_walpha(15,1)=9.6614159129071119D-2
m_walpha(15,2)=2.2106145785079101D-1
m_walpha(15,3)=3.3665061978407636D-1
m_walpha(15,4)=4.3847421642935351D-1
m_walpha(15,5)=5.2235011551287746D-1
m_walpha(15,6)=5.8484203003381963D-1
m_walpha(15,7)=6.2339089654456458D-1
m_walpha(15,8)=6.3641831661047915D-1
m_walpha(15,9)=6.2339089654456458D-1
m_walpha(15,10)=5.8484203003381963D-1
m_walpha(15,11)=5.2235011551287746D-1
m_walpha(15,12)=4.3847421642935351D-1
m_walpha(15,13)=3.3665061978407636D-1
m_walpha(15,14)=2.2106145785079101D-1
m_walpha(15,15)=9.6614159129071119D-2
m_walpha(16,1)=8.5301967014861703D-2
m_walpha(16,2)=1.9557521346573255D-1
m_walpha(16,3)=2.9894928122825785D-1
m_walpha(16,4)=3.9153346052083872D-1
m_walpha(16,5)=4.6996965947265833D-1
m_walpha(16,6)=5.3142087863815935D-1
m_walpha(16,7)=5.7366554722553987D-1
m_walpha(16,8)=5.9517664602374486D-1
m_walpha(16,9)=5.9517664602374486D-1
m_walpha(16,10)=5.7366554722553987D-1
m_walpha(16,11)=5.3142087863815935D-1
m_walpha(16,12)=4.6996965947265833D-1
m_walpha(16,13)=3.9153346052083872D-1
m_walpha(16,14)=2.9894928122825785D-1
m_walpha(16,15)=1.9557521346573255D-1
m_walpha(16,16)=8.5301967014861703D-2
m_walpha(17,1)=7.5864130888491514D-2
m_walpha(17,2)=1.7423125005286554D-1
m_walpha(17,3)=2.6714893884282217D-1
m_walpha(17,4)=3.5149347239816092D-1
m_walpha(17,5)=4.2454342241352301D-1
m_walpha(17,6)=4.8394903131555572D-1
m_walpha(17,7)=5.2780045310765260D-1
m_walpha(17,8)=5.5468809807888326D-1
m_walpha(17,9)=5.6374771298367703D-1
m_walpha(17,10)=5.5468809807888326D-1
m_walpha(17,11)=5.2780045310765260D-1
m_walpha(17,12)=4.8394903131555572D-1
m_walpha(17,13)=4.2454342241352301D-1
m_walpha(17,14)=3.5149347239816092D-1
m_walpha(17,15)=2.6714893884282217D-1
m_walpha(17,16)=1.7423125005286554D-1
m_walpha(17,17)=7.5864130888491514D-2
m_walpha(18,1)=6.7908709294697567D-2
m_walpha(18,2)=1.5618286158496769D-1
m_walpha(18,3)=2.4009851271399466D-1
m_walpha(18,4)=3.1711878420264865D-1
m_walpha(18,5)=3.8501853706395925D-1
m_walpha(18,6)=4.4184274750877224D-1
m_walpha(18,7)=4.8595623899959872D-1
m_walpha(18,8)=5.1608979449347116D-1
m_walpha(18,9)=5.3137646772768331D-1
m_walpha(18,10)=5.3137646772768331D-1
m_walpha(18,11)=5.1608979449347116D-1
m_walpha(18,12)=4.8595623899959872D-1
m_walpha(18,13)=4.4184274750877224D-1
m_walpha(18,14)=3.8501853706395925D-1
m_walpha(18,15)=3.1711878420264865D-1
m_walpha(18,16)=2.4009851271399466D-1
m_walpha(18,17)=1.5618286158496769D-1
m_walpha(18,18)=6.7908709294697567D-2
m_walpha(19,1)=6.1141010928229008D-2
m_walpha(19,2)=1.4078804558342894D-1
m_walpha(19,3)=2.1690982823504019D-1
m_walpha(19,4)=2.8742437980586025D-1
m_walpha(19,5)=3.5049695403716089D-1
m_walpha(19,6)=4.0449250283415413D-1
m_walpha(19,7)=4.4801216790127013D-1
m_walpha(19,8)=4.7992867547209405D-1
m_walpha(19,9)=4.9941535055611331D-1
m_walpha(19,10)=5.0596747647288465D-1
m_walpha(19,11)=4.9941535055611331D-1
m_walpha(19,12)=4.7992867547209405D-1
m_walpha(19,13)=4.4801216790127013D-1
m_walpha(19,14)=4.0449250283415413D-1
m_walpha(19,15)=3.5049695403716089D-1
m_walpha(19,16)=2.8742437980586025D-1
m_walpha(19,17)=2.1690982823504019D-1
m_walpha(19,18)=1.4078804558342894D-1
m_walpha(19,19)=6.1141010928229008D-2
m_walpha(20,1)=5.5336035428638466D-2
m_walpha(20,2)=1.2755315358613732D-1
m_walpha(20,3)=1.9689004663186147D-1
m_walpha(20,4)=2.6162159955227133D-1
m_walpha(20,5)=3.2022291559736995D-1
m_walpha(20,6)=3.7131907330479027D-1
m_walpha(20,7)=4.1371205911317567D-1
m_walpha(20,8)=4.4640809313832121D-1
m_walpha(20,9)=4.6864075841638154D-1
m_walpha(20,10)=4.7988891882084600D-1
m_walpha(20,11)=4.7988891882084600D-1
m_walpha(20,12)=4.6864075841638154D-1
m_walpha(20,13)=4.4640809313832121D-1
m_walpha(20,14)=4.1371205911317567D-1
m_walpha(20,15)=3.7131907330479027D-1
m_walpha(20,16)=3.2022291559736995D-1
m_walpha(20,17)=2.6162159955227133D-1
m_walpha(20,18)=1.9689004663186147D-1
m_walpha(20,19)=1.2755315358613732D-1
m_walpha(20,20)=5.5336035428638466D-2

!---------------------------------------------
!Matrice per beta, dove la funzione peso e' Sin(beta)
!---------------------------------------------
m_xbeta(1,1)=1.5707963267948966D0
m_xbeta(2,1)=8.8712893670499358D-1
m_xbeta(2,2)=2.2544637168847997D0
m_xbeta(3,1)=5.5819508679658935D-1
m_xbeta(3,2)=1.5707963267948966D0
m_xbeta(3,3)=2.5833975667932039D0
m_xbeta(4,1)=3.8011976290004088D-1
m_xbeta(4,2)=1.1315088599348815D0
m_xbeta(4,3)=2.0100837936549118D0
m_xbeta(4,4)=2.7614728906897524D0
m_xbeta(5,1)=2.7435604676190318D-1
m_xbeta(5,2)=8.4480958885359663D-1
m_xbeta(5,3)=1.5707963267948966D0
m_xbeta(5,4)=2.2967830647361966D0
m_xbeta(5,5)=2.8672366068278901D0
m_xbeta(6,1)=2.0688502471741067D-1
m_xbeta(6,2)=6.5100566127776101D-1
m_xbeta(6,3)=1.2469442125820415D0
m_xbeta(6,4)=1.8946484410077517D0
m_xbeta(6,5)=2.4905869923120322D0
m_xbeta(6,6)=2.9347076288723826D0
m_xbeta(7,1)=1.6137945950861125D-1
m_xbeta(7,2)=5.1525636338817408D-1
m_xbeta(7,3)=1.0072943606181171D0
m_xbeta(7,4)=1.5707963267948966D0
m_xbeta(7,5)=2.1342982929716762D0
m_xbeta(7,6)=2.6263362902016192D0
m_xbeta(7,7)=2.9802131940811820D0
m_xbeta(8,1)=1.2930578661253905D-1
m_xbeta(8,2)=4.1707068133816903D-1
m_xbeta(8,3)=8.2732767891940417D-1
m_xbeta(8,4)=1.3142998193786654D0
m_xbeta(8,5)=1.8272928342111278D0
m_xbeta(8,6)=2.3142649746703891D0
m_xbeta(8,7)=2.7245219722516242D0
m_xbeta(8,8)=3.0122868669772542D0
m_xbeta(9,1)=1.0588124692379566D-1
m_xbeta(9,2)=3.4403557934869406D-1
m_xbeta(9,3)=6.8981465786634491D-1
m_xbeta(9,4)=1.1109864680874401D0
m_xbeta(9,5)=1.5707963267948966D0
m_xbeta(9,6)=2.0306061855023532D0
m_xbeta(9,7)=2.4517779957234483D0
m_xbeta(9,8)=2.7975570742410992D0
m_xbeta(9,9)=3.0357114066659976D0
m_xbeta(10,1)=8.8266870606002430D-2
m_xbeta(10,2)=2.8837448753348158D-1
m_xbeta(10,3)=5.8291104117038096D-1
m_xbeta(10,4)=9.4865486974379368D-1
m_xbeta(10,5)=1.3584534397258137D0
m_xbeta(10,6)=1.7831392138639795D0
m_xbeta(10,7)=2.1929377838459996D0
m_xbeta(10,8)=2.5586816124194123D0
m_xbeta(10,9)=2.8532181660563117D0
m_xbeta(10,10)=3.0533257829837908D0
m_xbeta(11,1)=7.4695232256887551D-2
m_xbeta(11,2)=2.4505418973274662D-1
m_xbeta(11,3)=4.9844561303943445D-1
m_xbeta(11,4)=8.1778888601126834D-1
m_xbeta(11,5)=1.1826856325847383D0
m_xbeta(11,6)=1.5707963267948966D0
m_xbeta(11,7)=1.9589070210050549D0
m_xbeta(11,8)=2.3238037675785249D0
m_xbeta(11,9)=2.6431470405503588D0
m_xbeta(11,10)=2.8965384638570466D0
m_xbeta(11,11)=3.0668974213329057D0
m_xbeta(12,1)=6.4021082203569426D-2
m_xbeta(12,2)=2.1071783641313667D-1
m_xbeta(12,3)=4.3071013629204451D-1
m_xbeta(12,4)=7.1119646407302626D-1
m_xbeta(12,5)=1.0366162020413269D0
m_xbeta(12,6)=1.3896358585027873D0
m_xbeta(12,7)=1.7519567950870059D0
m_xbeta(12,8)=2.1049764515484663D0
m_xbeta(12,9)=2.4303961895167670D0
m_xbeta(12,10)=2.7108825172977487D0
m_xbeta(12,11)=2.9308748171766566D0
m_xbeta(12,12)=3.0775715713862238D0
m_xbeta(13,1)=5.5476674127557389D-2
m_xbeta(13,2)=1.8306514205134980D-1
m_xbeta(13,3)=3.7565409616647830D-1
m_xbeta(13,4)=6.2348475941721876D-1
m_xbeta(13,5)=9.1451881007395918D-1
m_xbeta(13,6)=1.2351463238876992D0
m_xbeta(13,7)=1.5707963267948966D0
m_xbeta(13,8)=1.9064463297020941D0
m_xbeta(13,9)=2.2270738435158341D0
m_xbeta(13,10)=2.5181078941725745D0
m_xbeta(13,11)=2.7659385574233149D0
m_xbeta(13,12)=2.9585275115384434D0
m_xbeta(13,13)=3.0861159794622358D0
m_xbeta(14,1)=4.8532097368360352D-2
m_xbeta(14,2)=1.6048109778978809D-1
m_xbeta(14,3)=3.3035584489075966D-1
m_xbeta(14,4)=5.5060173287106814D-1
m_xbeta(14,5)=8.1178392414810371D-1
m_xbeta(14,6)=1.1030980990321637D0
m_xbeta(14,7)=1.4128319451121770D0
m_xbeta(14,8)=1.7287607084776163D0
m_xbeta(14,9)=2.0384945545576295D0
m_xbeta(14,10)=2.3298087294416895D0
m_xbeta(14,11)=2.5909909207187251D0
m_xbeta(14,12)=2.8112368086990336D0
m_xbeta(14,13)=2.9811115558000052D0
m_xbeta(14,14)=3.0930605562214329D0
m_xbeta(15,1)=4.2812235337337357D-2
m_xbeta(15,2)=1.4180692123709423D-1
m_xbeta(15,3)=2.9267389393658135D-1
m_xbeta(15,4)=4.8948266804606946D-1
m_xbeta(15,5)=7.2474988329063318D-1
m_xbeta(15,6)=9.8981031951733136D-1
m_xbeta(15,7)=1.2751663937224633D0
m_xbeta(15,8)=1.5707963267948966D0
m_xbeta(15,9)=1.8664262598673300D0
m_xbeta(15,10)=2.1517823340724619D0
m_xbeta(15,11)=2.4168427702991601D0
m_xbeta(15,12)=2.6521099855437238D0
m_xbeta(15,13)=2.8489187596532119D0
m_xbeta(15,14)=2.9997857323526990D0
m_xbeta(15,15)=3.0987804182524559D0
m_xbeta(16,1)=3.8045613558666141D-2
m_xbeta(16,2)=1.2619483942824516D-1
m_xbeta(16,3)=2.6101443634960295D-1
m_xbeta(16,4)=4.3778944819439625D-1
m_xbeta(16,5)=6.5051846472852565D-1
m_xbeta(16,6)=8.9218524581929117D-1
m_xbeta(16,7)=1.1550243500607072D0
m_xbeta(16,8)=1.4307618825479289D0
m_xbeta(16,9)=1.7108307710418644D0
m_xbeta(16,10)=1.9865683035290861D0
m_xbeta(16,11)=2.2494074077705021D0
m_xbeta(16,12)=2.4910741888612676D0
m_xbeta(16,13)=2.7038032053953970D0
m_xbeta(16,14)=2.8805782172401903D0
m_xbeta(16,15)=3.0153978141615481D0
m_xbeta(16,16)=3.1035470400311271D0
m_xbeta(17,1)=3.4031899584448591D-2
m_xbeta(17,2)=1.1301368197451478D-1
m_xbeta(17,3)=2.3417417851641769D-1
m_xbeta(17,4)=3.9372126805807085D-1
m_xbeta(17,5)=5.8679268868185191D-1
m_xbeta(17,6)=8.0765654042111924D-1
m_xbeta(17,7)=1.0499141317998898D0
m_xbeta(17,8)=1.3066886151777079D0
m_xbeta(17,9)=1.5707963267948966D0
m_xbeta(17,10)=1.8349040384120854D0
m_xbeta(17,11)=2.0916785217899034D0
m_xbeta(17,12)=2.3339361131686740D0
m_xbeta(17,13)=2.5547999649079413D0
m_xbeta(17,14)=2.7478713855317224D0
m_xbeta(17,15)=2.9074184750733755D0
m_xbeta(17,16)=3.0285789716152785D0
m_xbeta(17,17)=3.1075607540053446D0
m_xbeta(18,1)=3.0620681519636121D-2
m_xbeta(18,2)=1.0178608268265892D-1
m_xbeta(18,3)=2.1123248331639610D-1
m_xbeta(18,4)=3.5587764079887141D-1
m_xbeta(18,5)=5.3174490986019327D-1
m_xbeta(18,6)=7.3411247730572803D-1
m_xbeta(18,7)=9.5766912547259807D-1
m_xbeta(18,8)=1.1966625463457866D0
m_xbeta(18,9)=1.4450363174513588D0
m_xbeta(18,10)=1.6965563361384344D0
m_xbeta(18,11)=1.9449301072440066D0
m_xbeta(18,12)=2.1839235281171952D0
m_xbeta(18,13)=2.4074801762840652D0
m_xbeta(18,14)=2.6098477437296000D0
m_xbeta(18,15)=2.7857150127909218D0
m_xbeta(18,16)=2.9303601702733971D0
m_xbeta(18,17)=3.0398065709071343D0
m_xbeta(18,18)=3.1109719720701571D0
m_xbeta(19,1)=2.7697270683011191D-2
m_xbeta(19,2)=9.2145823656152250D-2
m_xbeta(19,3)=1.9147623255487096D-1
m_xbeta(19,4)=3.2315894870578621D-1
m_xbeta(19,5)=4.8391314684412630D-1
m_xbeta(19,6)=6.6981870377827689D-1
m_xbeta(19,7)=8.7643653076364654D-1
m_xbeta(19,8)=1.0989256462864816D0
m_xbeta(19,9)=1.3321531003996308D0
m_xbeta(19,10)=1.5707963267948966D0
m_xbeta(19,11)=1.8094395531901625D0
m_xbeta(19,12)=2.0426670073033116D0
m_xbeta(19,13)=2.2651561228261467D0
m_xbeta(19,14)=2.4717739498115163D0
m_xbeta(19,15)=2.6576795067456669D0
m_xbeta(19,16)=2.8184337048840070D0
m_xbeta(19,17)=2.9501164210349223D0
m_xbeta(19,18)=3.0494468299336410D0
m_xbeta(19,19)=3.1138953829067820D0
m_xbeta(20,1)=2.5172995233498733D-2
m_xbeta(20,2)=8.3808302693580758D-2
m_xbeta(20,3)=1.7434671211470294D-1
m_xbeta(20,4)=2.9469379653338350D-1
m_xbeta(20,5)=4.4212083888614678D-1
m_xbeta(20,6)=6.1334950870838540D-1
m_xbeta(20,7)=8.0464541030608077D-1
m_xbeta(20,8)=1.0119108980095357D0
m_xbeta(20,9)=1.2307736175843541D0
m_xbeta(20,10)=1.4566697916551874D0
m_xbeta(20,11)=1.6849228619346058D0
m_xbeta(20,12)=1.9108190360054392D0
m_xbeta(20,13)=2.1296817555802575D0
m_xbeta(20,14)=2.3369472432837125D0
m_xbeta(20,15)=2.5282431448814078D0
m_xbeta(20,16)=2.6994718147036465D0
m_xbeta(20,17)=2.8468988570564097D0
m_xbeta(20,18)=2.9672459414750903D0
m_xbeta(20,19)=3.0577843508962125D0
m_xbeta(20,20)=3.1164196583562945D0


m_wbeta(1,1)=2.0000000000000000D0
m_wbeta(2,1)=1.0000000000000000D0
m_wbeta(2,2)=1.0000000000000000D0
m_wbeta(3,1)=4.5584040803912226D-1
m_wbeta(3,2)=1.0883191839217555D0
m_wbeta(3,3)=4.5584040803912226D-1
m_wbeta(4,1)=2.2407061812762016D-1
m_wbeta(4,2)=7.7592938187237984D-1
m_wbeta(4,3)=7.7592938187237984D-1
m_wbeta(4,4)=2.2407061812762016D-1
m_wbeta(5,1)=1.2011199820616150D-1
m_wbeta(5,2)=5.0378251238821888D-1
m_wbeta(5,3)=7.5221097881123923D-1
m_wbeta(5,4)=5.0378251238821888D-1
m_wbeta(5,5)=1.2011199820616150D-1
m_wbeta(6,1)=6.9387748499877836D-2
m_wbeta(6,2)=3.2479855137722084D-1
m_wbeta(6,3)=6.0581370012290133D-1
m_wbeta(6,4)=6.0581370012290133D-1
m_wbeta(6,5)=3.2479855137722084D-1
m_wbeta(6,6)=6.9387748499877836D-2
m_wbeta(7,1)=4.2625750150901831D-2
m_wbeta(7,2)=2.1353015976164668D-1
m_wbeta(7,3)=4.5607388993189171D-1
m_wbeta(7,4)=5.7554040031111957D-1
m_wbeta(7,5)=4.5607388993189171D-1
m_wbeta(7,6)=2.1353015976164668D-1
m_wbeta(7,7)=4.2625750150901831D-2
m_wbeta(8,1)=2.7535633513767011D-2
m_wbeta(8,2)=1.4420409203022751D-1
m_wbeta(8,3)=3.3626447785280460D-1
m_wbeta(8,4)=4.9199579660320088D-1
m_wbeta(8,5)=4.9199579660320088D-1
m_wbeta(8,6)=3.3626447785280460D-1
m_wbeta(8,7)=1.4420409203022751D-1
m_wbeta(8,8)=2.7535633513767011D-2
m_wbeta(9,1)=1.8540690233345113D-2
m_wbeta(9,2)=1.0010472234198563D-1
m_wbeta(9,3)=2.4776782327152767D-1
m_wbeta(9,4)=4.0043355770226815D-1
m_wbeta(9,5)=4.6630641290174688D-1
m_wbeta(9,6)=4.0043355770226815D-1
m_wbeta(9,7)=2.4776782327152767D-1
m_wbeta(9,8)=1.0010472234198563D-1
m_wbeta(9,9)=1.8540690233345113D-2
m_wbeta(10,1)=1.2923574870973804D-2
m_wbeta(10,2)=7.1311631869472290D-2
m_wbeta(10,3)=1.8407934563768848D-1
m_wbeta(10,4)=3.1897307573622444D-1
m_wbeta(10,5)=4.1271237188564099D-1
m_wbeta(10,6)=4.1271237188564099D-1
m_wbeta(10,7)=3.1897307573622444D-1
m_wbeta(10,8)=1.8407934563768848D-1
m_wbeta(10,9)=7.1311631869472290D-2
m_wbeta(10,10)=1.2923574870973804D-2
m_wbeta(11,1)=9.2753069683002673D-3
m_wbeta(11,2)=5.2007339979733851D-2
m_wbeta(11,3)=1.3841860179438881D-1
m_wbeta(11,4)=2.5225098348326324D-1
m_wbeta(11,5)=3.5204407834712870D-1
m_wbeta(11,6)=3.9200737885437026D-1
m_wbeta(11,7)=3.5204407834712870D-1
m_wbeta(11,8)=2.5225098348326324D-1
m_wbeta(11,9)=1.3841860179438881D-1
m_wbeta(11,10)=5.2007339979733851D-2
m_wbeta(11,11)=9.2753069683002673D-3
m_wbeta(12,1)=6.8251304430123969D-3
m_wbeta(12,2)=3.8735670111072895D-2
m_wbeta(12,3)=1.0548056826517032D-1
m_wbeta(12,4)=1.9955864751696537D-1
m_wbeta(12,5)=2.9453802849521091D-1
m_wbeta(12,6)=3.5486195516856812D-1
m_wbeta(12,7)=3.5486195516856812D-1
m_wbeta(12,8)=2.9453802849521091D-1
m_wbeta(12,9)=1.9955864751696537D-1
m_wbeta(12,10)=1.0548056826517032D-1
m_wbeta(12,11)=3.8735670111072895D-2
m_wbeta(12,12)=6.8251304430123969D-3
m_wbeta(13,1)=5.1315088668852916D-3
m_wbeta(13,2)=2.9398056294239148D-2
m_wbeta(13,3)=8.1466283023251137D-2
m_wbeta(13,4)=1.5857380994293156D-1
m_wbeta(13,5)=2.4417546211714373D-1
m_wbeta(13,6)=3.1217069552828294D-1
m_wbeta(13,7)=3.3816836845453238D-1
m_wbeta(13,8)=3.1217069552828294D-1
m_wbeta(13,9)=2.4417546211714373D-1
m_wbeta(13,10)=1.5857380994293156D-1
m_wbeta(13,11)=8.1466283023251137D-2
m_wbeta(13,12)=2.9398056294239148D-2
m_wbeta(13,13)=5.1315088668852916D-3
m_wbeta(14,1)=3.9311912527090455D-3
m_wbeta(14,2)=2.2688833216019194D-2
m_wbeta(14,3)=6.3738471842993061D-2
m_wbeta(14,4)=1.2683268279318218D-1
m_wbeta(14,5)=2.0179416004657890D-1
m_wbeta(14,6)=2.7004165995687068D-1
m_wbeta(14,7)=3.1097300089164694D-1
m_wbeta(14,8)=3.1097300089164694D-1
m_wbeta(14,9)=2.7004165995687068D-1
m_wbeta(14,10)=2.0179416004657890D-1
m_wbeta(14,11)=1.2683268279318218D-1
m_wbeta(14,12)=6.3738471842993061D-2
m_wbeta(14,13)=2.2688833216019194D-2
m_wbeta(14,14)=3.9311912527090455D-3
m_wbeta(15,1)=3.0616594839236216D-3
m_wbeta(15,2)=1.7775620830836013D-2
m_wbeta(15,3)=5.0480612367418813D-2
m_wbeta(15,4)=1.0221293737820382D-1
m_wbeta(15,5)=1.6685635618251158D-1
m_wbeta(15,6)=2.3142033539202311D-1
m_wbeta(15,7)=2.7951722985166940D-1
m_wbeta(15,8)=2.9735049702682728D-1
m_wbeta(15,9)=2.7951722985166940D-1
m_wbeta(15,10)=2.3142033539202311D-1
m_wbeta(15,11)=1.6685635618251158D-1
m_wbeta(15,12)=1.0221293737820382D-1
m_wbeta(15,13)=5.0480612367418813D-2
m_wbeta(15,14)=1.7775620830836013D-2
m_wbeta(15,15)=3.0616594839236216D-3
m_wbeta(16,1)=2.4194677567615628D-3
m_wbeta(16,2)=1.4115268156854008D-2
m_wbeta(16,3)=4.0437893946503669D-2
m_wbeta(16,4)=8.3026647573217742D-2
m_wbeta(16,5)=1.3834195526951273D-1
m_wbeta(16,6)=1.9741148870253456D-1
m_wbeta(16,7)=2.4763632094635522D-1
m_wbeta(16,8)=2.7661095764826050D-1
m_wbeta(16,9)=2.7661095764826050D-1
m_wbeta(16,10)=2.4763632094635522D-1
m_wbeta(16,11)=1.9741148870253456D-1
m_wbeta(16,12)=1.3834195526951273D-1
m_wbeta(16,13)=8.3026647573217742D-2
m_wbeta(16,14)=4.0437893946503669D-2
m_wbeta(16,15)=1.4115268156854008D-2
m_wbeta(16,16)=2.4194677567615628D-3
m_wbeta(17,1)=1.9369674486454401D-3
m_wbeta(17,2)=1.1345551865279835D-2
m_wbeta(17,3)=3.2736482199317563D-2
m_wbeta(17,4)=6.7978435070800825D-2
m_wbeta(17,5)=1.1515990849948726D-1
m_wbeta(17,6)=1.6814565845293004D-1
m_wbeta(17,7)=2.1744637945414307D-1
m_wbeta(17,8)=2.5258319446354658D-1
m_wbeta(17,9)=2.6533484509169878D-1
m_wbeta(17,10)=2.5258319446354658D-1
m_wbeta(17,11)=2.1744637945414307D-1
m_wbeta(17,12)=1.6814565845293004D-1
m_wbeta(17,13)=1.1515990849948726D-1
m_wbeta(17,14)=6.7978435070800825D-2
m_wbeta(17,15)=3.2736482199317563D-2
m_wbeta(17,16)=1.1345551865279835D-2
m_wbeta(17,17)=1.9369674486454401D-3
m_wbeta(18,1)=1.5688436253924968D-3
m_wbeta(18,2)=9.2199875076378784D-3
m_wbeta(18,3)=2.6761374213526816D-2
m_wbeta(18,4)=5.6089752225160450D-2
m_wbeta(18,5)=9.6316257117859740D-2
m_wbeta(18,6)=1.4329022312108587D-1
m_wbeta(18,7)=1.8995726748033847D-1
m_wbeta(18,8)=2.2778381639235996D-1
m_wbeta(18,9)=2.4901247831663831D-1
m_wbeta(18,10)=2.4901247831663831D-1
m_wbeta(18,11)=2.2778381639235996D-1
m_wbeta(18,12)=1.8995726748033847D-1
m_wbeta(18,13)=1.4329022312108587D-1
m_wbeta(18,14)=9.6316257117859740D-2
m_wbeta(18,15)=5.6089752225160450D-2
m_wbeta(18,16)=2.6761374213526816D-2
m_wbeta(18,17)=9.2199875076378784D-3
m_wbeta(18,18)=1.5688436253924968D-3
m_wbeta(19,1)=1.2840821217594151D-3
m_wbeta(19,2)=7.5677017445384936D-3
m_wbeta(19,3)=2.2074656244751675D-2
m_wbeta(19,4)=4.6625489045122395D-2
m_wbeta(19,5)=8.0967083885420315D-2
m_wbeta(19,6)=1.2233034854344086D-1
m_wbeta(19,7)=1.6551794706902487D-1
m_wbeta(19,8)=2.0373683621316769D-1
m_wbeta(19,9)=2.3012138592799148D-1
m_wbeta(19,10)=2.3954893840956562D-1
m_wbeta(19,11)=2.3012138592799148D-1
m_wbeta(19,12)=2.0373683621316769D-1
m_wbeta(19,13)=1.6551794706902487D-1
m_wbeta(19,14)=1.2233034854344086D-1
m_wbeta(19,15)=8.0967083885420315D-2
m_wbeta(19,16)=4.6625489045122395D-2
m_wbeta(19,17)=2.2074656244751675D-2
m_wbeta(19,18)=7.5677017445384936D-3
m_wbeta(19,19)=1.2840821217594151D-3
m_wbeta(20,1)=1.0610419630623757D-3
m_wbeta(20,2)=6.2681947298298217D-3
m_wbeta(20,3)=1.8360759375151004D-2
m_wbeta(20,4)=3.9033364624202392D-2
m_wbeta(20,5)=6.8420974360331419D-2
m_wbeta(20,6)=1.0471375385635442D-1
m_wbeta(20,7)=1.4410932145222802D-1
m_wbeta(20,8)=1.8127933543040475D-1
m_wbeta(20,9)=2.1037592857364366D-1
m_wbeta(20,10)=2.2637732563479214D-1
m_wbeta(20,11)=2.2637732563479214D-1
m_wbeta(20,12)=2.1037592857364366D-1
m_wbeta(20,13)=1.8127933543040475D-1
m_wbeta(20,14)=1.4410932145222802D-1
m_wbeta(20,15)=1.0471375385635442D-1
m_wbeta(20,16)=6.8420974360331419D-2
m_wbeta(20,17)=3.9033364624202392D-2
m_wbeta(20,18)=1.8360759375151004D-2
m_wbeta(20,19)=6.2681947298298217D-3
m_wbeta(20,20)=1.0610419630623757D-3

!Assegno i pesi per l'output
v_xalpha=m_xalpha(nalpha,1:nalpha)
v_walpha=m_walpha(nalpha,1:nalpha)
v_xbeta=m_xbeta(nbeta,1:nbeta)
v_wbeta=m_wbeta(nbeta,1:nbeta)

END SUBROUTINE weigths



SUBROUTINE weigths1(nalpha,nbeta,v_walpha,v_xalpha,v_wbeta,v_xbeta,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nalpha,nbeta								!Numero punti integrazione
REAL(dbl), DIMENSION(:), INTENT(INOUT) :: v_walpha,v_xalpha				!Vettori peso e coord per alpha
REAL(dbl), DIMENSION(:), INTENT(INOUT) :: v_wbeta,v_xbeta					!Vettori peso e coord per beta
INTEGER(lo), INTENT(OUT) :: error										!Flag di errore

! Dichiarazione variabili interne

! Subroutine vera e propria

!Controllo l'errore della subroutine
weigths_if : IF (( nalpha>100 ) .OR. ( nbeta>100 )) THEN
                WRITE(*,*)
                WRITE(*,10)
                10 FORMAT ("Steps dell'integrazione out of bounds")
                RETURN
END IF weigths_if



error=0

!---------------------------------------------
!Matrice per alpha, dove la funzione peso e' 1.
!---------------------------------------------
	v_xalpha(1)=2.1832757774497611D-3
        v_walpha(1)=5.6023416145698184D-3
        v_xalpha(2)=1.1497862262221792D-2
        v_walpha(2)=1.3028289225575575D-2
        v_xalpha(3)=2.8232326149200534D-2
        v_walpha(3)=2.0434357370926012D-2
        v_xalpha(4)=5.2350702377344549D-2
        v_walpha(4)=2.7792915678578293D-2
        v_xalpha(5)=8.3796241863151145D-2
        v_walpha(5)=3.5085744882214234D-2
        v_xalpha(6)=1.2249446856775583D-1
        v_walpha(6)=4.2295412367258805D-2
        v_xalpha(7)=1.6835363106081901D-1
        v_walpha(7)=4.9404781813608457D-2
        v_xalpha(8)=2.2126497141654591D-1
        v_walpha(8)=5.6396978400262770D-2
        v_xalpha(9)=2.8110299681309014D-1
        v_walpha(9)=6.3255412500460179D-2
        v_xalpha(10)=3.4772578172761225D-1
        v_walpha(10)=6.9963814430584659D-2
        v_xalpha(11)=4.2097530633571288D-1
        v_walpha(11)=7.6506271484919367D-2
        v_xalpha(12)=5.0067783209032886D-1
        v_walpha(12)=8.2867265077102972D-2
        v_xalpha(13)=5.8664431418062425D-1
        v_walpha(13)=8.9031707290097443D-2
        v_xalpha(14)=6.7867085011593143D-1
        v_walpha(14)=9.4984976542466313D-2
        v_xalpha(15)=7.7653916346395678D-1
        v_walpha(15)=1.0071295220890957D-1
        v_xalpha(16)=8.8001712163803243D-1
        v_walpha(16)=1.0620204808097812D-1
        v_xalpha(17)=9.8885928652565835D-1
        v_walpha(17)=1.1143924457439217D-1
        v_xalpha(18)=1.1028074966627236D0
        v_walpha(18)=1.1641211959984361D-1
        v_xalpha(19)=1.2215914795781418D0
        v_walpha(19)=1.2110887802070008D-1
        v_xalpha(20)=1.3449294928596106D0
        v_walpha(20)=1.2551837962598257D-1
        v_xalpha(21)=1.4725289924217697D0
        v_walpha(21)=1.2963016555131031D-1
        v_xalpha(22)=1.6040873263927805D0
        v_walpha(22)=1.3343448308460509D-1
        v_xalpha(23)=1.7392924529741756D0
        v_walpha(23)=1.3692230879740383D-1
        v_xalpha(24)=1.8778236805717116D0
        v_walpha(24)=1.4008536994672600D-1
        v_xalpha(25)=2.0193524284419924D0
        v_walpha(25)=1.4291616409661711D-1
        v_xalpha(26)=2.1635430060508709D0
        v_walpha(26)=1.4540797691275535D-1
        v_xalpha(27)=2.3100534092951996D0
        v_walpha(27)=1.4755489808786776D-1
        v_xalpha(28)=2.4585361316994625D0
        v_walpha(28)=1.4935183536015226D-1
        v_xalpha(29)=2.6086389886632791D0
        v_walpha(29)=1.5079452659143512D-1
        v_xalpha(30)=2.7600059528047963D0
        v_walpha(30)=1.5187954987640225D-1
        v_xalpha(31)=2.9122779984186587D0
        v_walpha(31)=1.5260433165891736D-1
        v_xalpha(32)=3.0650939530456176D0
        v_walpha(32)=1.5296715283616980D-1
        v_xalpha(33)=3.2180913541339689D0
        v_walpha(33)=1.5296715283616980D-1
        v_xalpha(34)=3.3709073087609278D0
        v_walpha(34)=1.5260433165891736D-1
        v_xalpha(35)=3.5231793543747902D0
        v_walpha(35)=1.5187954987640225D-1
        v_xalpha(36)=3.6745463185163074D0
        v_walpha(36)=1.5079452659143512D-1
        v_xalpha(37)=3.8246491754801239D0
        v_walpha(37)=1.4935183536015226D-1
        v_xalpha(38)=3.9731318978843869D0
        v_walpha(38)=1.4755489808786776D-1
        v_xalpha(39)=4.1196423011287156D0
        v_walpha(39)=1.4540797691275535D-1
        v_xalpha(40)=4.2638328787375941D0
        v_walpha(40)=1.4291616409661711D-1
        v_xalpha(41)=4.4053616266078749D0
        v_walpha(41)=1.4008536994672600D-1
        v_xalpha(42)=4.5438928542054109D0
        v_walpha(42)=1.3692230879740383D-1
        v_xalpha(43)=4.6790979807868060D0
        v_walpha(43)=1.3343448308460509D-1
        v_xalpha(44)=4.8106563147578168D0
        v_walpha(44)=1.2963016555131031D-1
        v_xalpha(45)=4.9382558143199759D0
        v_walpha(45)=1.2551837962598257D-1
        v_xalpha(46)=5.0615938276014447D0
        v_walpha(46)=1.2110887802070008D-1
        v_xalpha(47)=5.1803778105168629D0
        v_walpha(47)=1.1641211959984361D-1
        v_xalpha(48)=5.2943260206539281D0
        v_walpha(48)=1.1143924457439217D-1
        v_xalpha(49)=5.4031681855415540D0
        v_walpha(49)=1.0620204808097812D-1
        v_xalpha(50)=5.5066461437156297D0
        v_walpha(50)=1.0071295220890957D-1
        v_xalpha(51)=5.6045144570636550D0
        v_walpha(51)=9.4984976542466313D-2
        v_xalpha(52)=5.6965409929989622D0
        v_walpha(52)=8.9031707290097443D-2
        v_xalpha(53)=5.7825074750892576D0
        v_walpha(53)=8.2867265077102972D-2
        v_xalpha(54)=5.8622100008438736D0
        v_walpha(54)=7.6506271484919367D-2
        v_xalpha(55)=5.9354595254519742D0
        v_walpha(55)=6.9963814430584659D-2
        v_xalpha(56)=6.0020823103664963D0
        v_walpha(56)=6.3255412500460179D-2
        v_xalpha(57)=6.0619203357630406D0
        v_walpha(57)=5.6396978400262770D-2
        v_xalpha(58)=6.1148316761187675D0
        v_walpha(58)=4.9404781813608457D-2
        v_xalpha(59)=6.1606908386118307D0
        v_walpha(59)=4.2295412367258805D-2
        v_xalpha(60)=6.1993890653164353D0
        v_walpha(60)=3.5085744882214234D-2
        v_xalpha(61)=6.2308346048022419D0
        v_walpha(61)=2.7792915678578293D-2
        v_xalpha(62)=6.2549529810303859D0
        v_walpha(62)=2.0434357370926012D-2
        v_xalpha(63)=6.2716874449173647D0
        v_walpha(63)=1.3028289225575575D-2
        v_xalpha(64)=6.2810020314021367D0
        v_walpha(64)=5.6023416145698184D-3

!---------------------------------------------
!Matrice per beta, dove la funzione peso e' Sin(beta)
!---------------------------------------------
	v_xbeta(1)=2.6972385680803485D-3
        v_wbeta(1)=1.2215020492883368D-5
        v_xbeta(2)=9.0356896742528891D-3
        v_wbeta(2)=7.3593756884554777D-5
        v_xbeta(3)=1.8980044582092850D-2
        v_wbeta(3)=2.2280825083203630D-4
        v_xbeta(4)=3.2505869741535165D-2
        v_wbeta(4)=4.9746303778877104D-4
        v_xbeta(5)=4.9580870527229776D-2
        v_wbeta(5)=9.3355976908733305D-4
        v_xbeta(6)=7.0164452102580136D-2
        v_wbeta(6)=1.5649451985242465D-3
        v_xbeta(7)=9.4207769719848453D-2
        v_wbeta(7)=2.4227211579765684D-3
        v_xbeta(8)=1.2165385904980514D-1
        v_wbeta(8)=3.5346162067368192D-3
        v_xbeta(9)=1.5243780016689957D-1
        v_wbeta(9)=4.9243241954268032D-3
        v_xbeta(10)=1.8648690544184348D-1
        v_wbeta(10)=6.6108206690576776D-3
        v_xbeta(11)=2.2372092778209790D-1
        v_wbeta(11)=8.6076738434502351D-3
        v_xbeta(12)=2.6405228707062770D-1
        v_wbeta(12)=1.0922372510290281D-2
        v_xbeta(13)=3.0738631306034678D-1
        v_wbeta(13)=1.3555698288758909D-2
        v_xbeta(14)=3.5362150314283185D-1
        v_wbeta(14)=1.6501173754010056D-2
        v_xbeta(15)=4.0264979351014949D-1
        v_wbeta(15)=1.9744620743810075D-2
        v_xbeta(16)=4.5435684231828653D-1
        v_wbeta(16)=2.3263864210844433D-2
        v_xbeta(17)=5.0862232355280492D-1
        v_wbeta(17)=2.7028616043254021D-2
        v_xbeta(18)=5.6532023039371640D-1
        v_wbeta(18)=3.1000570102339299D-2
        v_xbeta(19)=6.2431918697543585D-1
        v_wbeta(19)=3.5133734225510602D-2
        v_xbeta(20)=6.8548276753641007D-1
        v_wbeta(20)=3.9375017160733502D-2
        v_xbeta(21)=7.4866982204912256D-1
        v_wbeta(21)=4.3665078544981071D-2
        v_xbeta(22)=8.1373480751259168D-1
        v_wbeta(22)=4.7939438492362910D-2
        v_xbeta(23)=8.8052812417469800D-1
        v_wbeta(23)=5.2129830659580863D-2
        v_xbeta(24)=9.4889645602971519D-1
        v_wbeta(24)=5.6165769490572633D-2
        v_xbeta(25)=1.0186831150067474D0
        v_wbeta(25)=5.9976289496814327D-2
        v_xbeta(26)=1.0897283883272318D0
        v_wbeta(26)=6.3491802747443585D-2
        v_xbeta(27)=1.1618698885643601D0
        v_wbeta(27)=6.6646011060795933D-2
        v_xbeta(28)=1.2349429059845019D0
        v_wbeta(28)=6.9377802470636641D-2
        v_xbeta(29)=1.3087807627908954D0
        v_wbeta(29)=7.1633058012823838D-2
        v_xbeta(30)=1.3832151689234650D0
        v_wbeta(30)=7.3366295170601707D-2
        v_xbeta(31)=1.4580765790961145D0
        v_wbeta(31)=7.4542078614803438D-2
        v_xbeta(32)=1.5331945507746624D0
        v_wbeta(32)=7.5136137092773946D-2
        v_xbeta(33)=1.6083981028151308D0
        v_wbeta(33)=7.5136137092773946D-2
        v_xbeta(34)=1.6835160744936788D0
        v_wbeta(34)=7.4542078614803438D-2
        v_xbeta(35)=1.7583774846663282D0
        v_wbeta(35)=7.3366295170601707D-2
        v_xbeta(36)=1.8328118907988978D0
        v_wbeta(36)=7.1633058012823838D-2
        v_xbeta(37)=1.9066497476052913D0
        v_wbeta(37)=6.9377802470636641D-2
        v_xbeta(38)=1.9797227650254332D0
        v_wbeta(38)=6.6646011060795933D-2
        v_xbeta(39)=2.0518642652625614D0
        v_wbeta(39)=6.3491802747443585D-2
        v_xbeta(40)=2.1229095385830459D0
        v_wbeta(40)=5.9976289496814327D-2
        v_xbeta(41)=2.1926961975600780D0
        v_wbeta(41)=5.6165769490572633D-2
        v_xbeta(42)=2.2610645294150952D0
        v_wbeta(42)=5.2129830659580863D-2
        v_xbeta(43)=2.3278578460772016D0
        v_wbeta(43)=4.7939438492362910D-2
        v_xbeta(44)=2.3929228315406707D0
        v_wbeta(44)=4.3665078544981071D-2
        v_xbeta(45)=2.4561098860533832D0
        v_wbeta(45)=3.9375017160733502D-2
        v_xbeta(46)=2.5172734666143574D0
        v_wbeta(46)=3.5133734225510602D-2
        v_xbeta(47)=2.5762724231960768D0
        v_wbeta(47)=3.1000570102339299D-2
        v_xbeta(48)=2.6329703300369883D0
        v_wbeta(48)=2.7028616043254021D-2
        v_xbeta(49)=2.6872358112715067D0
        v_wbeta(49)=2.3263864210844433D-2
        v_xbeta(50)=2.7389428600796438D0
        v_wbeta(50)=1.9744620743810075D-2
        v_xbeta(51)=2.7879711504469614D0
        v_wbeta(51)=1.6501173754010056D-2
        v_xbeta(52)=2.8342063405294465D0
        v_wbeta(52)=1.3555698288758909D-2
        v_xbeta(53)=2.8775403665191655D0
        v_wbeta(53)=1.0922372510290281D-2
        v_xbeta(54)=2.9178717258076953D0
        v_wbeta(54)=8.6076738434502351D-3
        v_xbeta(55)=2.9551057481479498D0
        v_wbeta(55)=6.6108206690576776D-3
        v_xbeta(56)=2.9891548534228937D0
        v_wbeta(56)=4.9243241954268032D-3
        v_xbeta(57)=3.0199387945399881D0
        v_wbeta(57)=3.5346162067368192D-3
        v_xbeta(58)=3.0473848838699448D0
        v_wbeta(58)=2.4227211579765684D-3
        v_xbeta(59)=3.0714282014872131D0
        v_wbeta(59)=1.5649451985242465D-3
        v_xbeta(60)=3.0920117830625635D0
        v_wbeta(60)=9.3355976908733305D-4
        v_xbeta(61)=3.1090867838482581D0
        v_wbeta(61)=4.9746303778877104D-4
        v_xbeta(62)=3.1226126090077004D0
        v_wbeta(62)=2.2280825083203630D-4
        v_xbeta(63)=3.1325569639155403D0
        v_wbeta(63)=7.3593756884554777D-5
        v_xbeta(64)=3.1388954150217129D0
        v_wbeta(64)=1.2215020492883368D-5

END SUBROUTINE weigths1


END MODULE basicsubs
