MODULE gmmsubs

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! INDICE DELLE SUBROUTINE CONTENUTE NEL MODULO
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines
! 1) openinfo(filename, unit_in, h_eps, h_patt, nspheres, neps, neq, error)
!	COMMENT: opens the pricipal file end gets infos from it
!	VARS:
!		filename: name of the principal file in
!		unit_in: associated unit in
!		h_eps: linked list for dielectric constants names out
!		h_patt: linked list for the pattern out
!		nspheres: total number of spheres out
!		neps: total number of materials out
!		neq: number identical spheres groups out
!
!	ERROR: 0 ok and 1 error in one of the called subroutines
!
!2) intspheres(r,xyz,ns,error)
!	COMMENT: check if spheres intersection exists
!	VARS:
!		r: radii rank 1 array in
!		xyz: position rank 2 array
!		ns: total number of spheres
!
!	ERROR: 1 means intersecting spheres, 0 everyghing is ok
!3) fillarray(unit_in,ns,neps,h_patt,h_eps,v_r,m_xyz,v_unit_in,v_patt,v_eps,v_epsg,error)
!	COMMENT: fills arrays previously allocated
!	VARS:
!		unit_in: associated unit in
!		h_eps: linked list for dielectric constants names in
!		h_patt: linked list for the pattern in
!		ns: total number of spheres in
!		neps: total number of materials in
!		v_r: sphere raddi array out
!		m_xyz: cartesian sphere position array
!		v_unit_in: dielectric function opening unit array
!		v_patt:	equality pattern array
!		v_eps: dielectric functions name array
!		v_epsg: sphere dielectric function array
!
!	ERROR: 0 ok and 1 error in lldestructor
!
!4) checkepsilon(v_unit_in,neps,error)
!	COMMENT: checks if dielectric functions are homogenous and well tabulated
!	VARS:
!		v_unit_in: vector of input units in
!		neps: total number of materials in
!
!	ERROR: 0 ok, 1 different number of points in tabulation, 2 different tabulations length
!
!5) correps(lambda,e1,e2,r,par)
!	COMMENT: Corrects real and imaginary par of the dielectric function
!	VARS:
!		lambda: wavelength in question in
!		e1,e2: real and imaginary part of the dielectric function
!		r: cluster radius
!		par: vector containinge the 3 parameters necessary to correct the dielectric function
!
!	ERROR: no error flag
!
!6) lambdaepsilon
!	COMMENT: Gives an output matrix with corrected dielectric functions for the equivalent spheres
!	VARS:
!		neps,neq: number of materials and number of equal spheres in
! 		lambda: incident wavelength in
! 		v_eps,v_epseq:character arrays for materials and materials of equal spheres in
! 		v_req: array of radii for eqnals spheres in
! 		m_eps,m_par: matrices for dielectric functions ad correction parameters in
! 		m_epseq: matrix for corrected dielectric functions of equal spheres out
!
!
!	ERROR: no error flag
!
!
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

USE kinds
USE datatypes
USE operators
USE basicsubs

CONTAINS


!******************************************************************************
!******************************************************************************
!******************************************************************************
!1) SUBROUTINE inOPENINFO: per aprire il filone ed estrarne le informazioni
!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE openinfo(filename,unit_in,h_eps,h_patt,h_req,h_epseq,nspheres,neps,neq,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
INTEGER(lo), INTENT(IN) :: unit_in 						! Unita' apertura e errore
CHARACTER(len=*), INTENT(IN) :: filename					! Nome file principale
INTEGER(lo), INTENT(OUT) :: nspheres, neps, neq, error	! Numero sfere, numero materiali, num sfere identiche, errore
TYPE (NODE_CHAR), POINTER :: h_eps							! LL per le funzione dielettriche
TYPE (NODE_LONG_INT), POINTER :: h_patt						! Pattern sfere uguali
TYPE (NODE_DBL_REAL), POINTER :: h_req						! LL raggi per sfere identiche
TYPE (NODE_CHAR), POINTER :: h_epseq						! LL funzioni dielettriche per sfere identiche

! Dichiarazione variabili interne
REAL(dbl) :: x,y,z,r							! Posizione e raggio
INTEGER(lo) :: status_in						! Status apertura
CHARACTER(len=length) :: eps							! Nome file della funzione dielettrica
TYPE (RDBL_CHAR) :: input						! input per la routine patter
TYPE (NODE_RDBL_CHAR), POINTER :: h_aus=>NULL()	! LL ausiliaria per costruire il pattern

! Inizio della procedura vera e propria

! Apertura del file
OPEN(UNIT=unit_in,FILE=filename, ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ',IOSTAT=status_in)

! Controllo apertura file
openif: IF (status_in==0) THEN
			WRITE (*,*)
			WRITE (*,10) filename
			10 FORMAT ("Il file e' stato aperto con successo, il suo nome e' ", A)
			error=0
		ELSE
			WRITE (*,*)
            WRITE (*,20) status_in
			20 FORMAT ("Si e' verificato un errore di apertura in inputinfo, IOSTAT= ",I5 //,"Il programma termina ora...")
			error=1
            RETURN
        END IF openif

! Inizializzo i contatori
nspheres=0
neps=0
neq=0

! Ciclo per il reperimento delle informazioni
read_do: DO
			READ(unit_in,*,IOSTAT=status_in) x,y,z,r,eps
			IF (status_in<0) EXIT					! Esco a fine file

			nspheres=nspheres+1						! Aggiorno il numero delle sfere



			CALL c_addnode(h_eps,eps,error,neps)	! Costruisco la LL per le funzioni dielettriche

			! C_Addnode funziona?
			eps_if: IF (error==1) THEN
						WRITE (*,*)
						WRITE (*,30)
						30 FORMAT ("Si e' verificato un errore nella scansione delle funzioni dielettriche")
						RETURN
					ELSE
						error=0
			        END IF eps_if

			! Preparo l'input per la subroutine pattern
			input%value=r
			input%text=eps

			! Costruisco il patter di uni e dui ecc.
			CALL pattern(h_aus,input,h_patt,neq,h_req,h_epseq,error)

			! pattern funziona?
			patt_if: IF (error==1) THEN
						WRITE (*,*)
						WRITE (*,40)
						40 FORMAT ("Si e' verificato un errore nella routine per la costruzione del pattern")
						RETURN
					ELSE
						error=0
			        END IF patt_if

END DO read_do

! Disalloco la memoria della ll ausiliaria per la costruzione del pattern
CALL lldestructor(h_aus,error)

! pattern funziona?
dest_if: IF (error==1) THEN
			WRITE (*,*)
			WRITE (*,50)
			50 FORMAT ("Si e' verificato un errore nella routine lldesctructor")
			RETURN
		ELSE
			error=0
        END IF dest_if

! Torno all'inizio del file
REWIND(UNIT=unit_in)

END SUBROUTINE openinfo




!******************************************************************************
!******************************************************************************
!******************************************************************************
!1bis) SUBROUTINE OPENINFO_SHELL: per aprire il filone ed estrarne le informazioni
!generalizzato per farlo funzionare con le shell!
!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE openinfo_shell(filename,unit_in,h_eps,h_patt,h_req,h_epseq,h_req_shell,h_epseq_shell,nspheres,neps,neq,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
INTEGER(lo), INTENT(IN) :: unit_in 						! Unita' apertura e errore
CHARACTER(len=*), INTENT(IN) :: filename						! Nome file principale
INTEGER(lo), INTENT(OUT) :: nspheres, neps, neq, error		! Numero sfere, numero materiali, num sfere identiche, errore
TYPE (NODE_CHAR), POINTER :: h_eps							! LL per le funzione dielettriche
TYPE (NODE_LONG_INT), POINTER :: h_patt						! Pattern sfere uguali
TYPE (NODE_DBL_REAL), POINTER :: h_req, h_req_shell				! LL raggi per sfere identiche
TYPE (NODE_CHAR), POINTER :: h_epseq, h_epseq_shell				! LL funzioni dielettriche per sfere identiche

! Dichiarazione variabili interne
REAL(dbl) :: x,y,z,r,r_shell								! Posizione e raggio
INTEGER(lo) :: status_in								! Status apertura
CHARACTER(len=length) :: eps, eps_shell						! Nome file della funzione dielettrica
TYPE (RDBL_CHAR) :: input, input_shell						! input per la routine patter
TYPE (NODE_RDBL_CHAR), POINTER :: h_aus, h_aus_shell=>NULL()		! LL ausiliaria per costruire il pattern

! Inizio della procedura vera e propria

! Apertura del file
OPEN(UNIT=unit_in,FILE=filename, ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ',IOSTAT=status_in)

! Controllo apertura file
openif: IF (status_in==0) THEN
			WRITE (*,*)
			WRITE (*,10) filename
			10 FORMAT ("Il file e' stato aperto con successo, il suo nome e' ", A)
			error=0
		ELSE
			WRITE (*,*)
            WRITE (*,20) status_in
			20 FORMAT ("Si e' verificato un errore di apertura in inputinfo, IOSTAT= ",I5 //,"Il programma termina ora...")
			error=1
            RETURN
        END IF openif

! Inizializzo i contatori
nspheres=0
neps=0
neq=0

! Ciclo per il reperimento delle informazioni
read_do: DO
	READ(unit_in,*,IOSTAT=status_in) x,y,z,r,r_shell,eps,eps_shell
	IF (status_in<0) EXIT					! Esco a fine file

	nspheres=nspheres+1						! Aggiorno il numero delle sfere


	!Quit aggiungo le epsilon se e solo se sono diverse dalle precedenti
	CALL c_addnode(h_eps,eps,error,neps)		! Aggiungo eventualmente eps dal core
	CALL c_addnode(h_eps,eps_shell,error,neps)	! Aggiungo eventualmente eps dalla shell

	! C_Addnode funziona?
	eps_if: IF (error==1) THEN
		WRITE (*,*)
		WRITE (*,30)
		30 FORMAT ("Si e' verificato un errore nella scansione delle funzioni dielettriche")
		RETURN
	ELSE
		error=0
	END IF eps_if

	! Preparo l'input per la subroutine pattern
	input%value=r
	input%text=eps
	input_shell%value=r_shell
	input_shell%text=eps_shell

	! Costruisco il patter di uni e dui ecc.
	CALL pattern_shell(h_aus,h_aus_shell,input,input_shell,h_patt,neq,h_req,h_epseq,h_req_shell,h_epseq_shell,error)

	! pattern funziona?
	patt_if: IF (error==1) THEN
		WRITE (*,*)
		WRITE (*,40)
		40 FORMAT ("Si e' verificato un errore nella routine per la costruzione del pattern")
		RETURN
	ELSE
		error=0
	END IF patt_if

END DO read_do

! Disalloco la memoria della ll ausiliaria per la costruzione del pattern
CALL lldestructor(h_aus,error)
CALL lldestructor(h_aus_shell,error)

! pattern funziona?
dest_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,50)
	50 FORMAT ("Si e' verificato un errore nella routine lldesctructor")
	RETURN
ELSE
	error=0
END IF dest_if

! Torno all'inizio del file
REWIND(UNIT=unit_in)

END SUBROUTINE openinfo_shell





!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE INTSPHERES: per vedere se ci sono sfere intersecantisi
!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE intspheres(r,xyz,ns,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(:),INTENT(IN) :: r					! Array raggi delle sfere
REAL(dbl), DIMENSION(:,:),INTENT(IN) :: xyz 			! Array posizione sfere
INTEGER(lo), INTENT(IN) :: ns						! Numero di sfere
INTEGER(lo), INTENT(OUT) :: error					! Errore

! Dichiarazione variabili interne
INTEGER(lo), DIMENSION(ns,ns), TARGET ::	a			! Array di controllo
INTEGER(lo) :: i,j,somma,flag						! Indici somma e flag
INTEGER(lo), DIMENSION(:,:), POINTER :: p=>NULL()	! Puntatore su a

! Inizio della procedura vera e propria
a=0

! fill_do_out: DO i=1,ns
! 		fill_do_in: DO j=1,i
!
! 			IF (j==i) CYCLE
!
! 			IF (dist(xyz(i,:),xyz(j,:))<(r(i)+r(j))) THEN
! 				a(i,j)=1
! 			END IF
!
! 		END DO fill_do_in
! END DO fill_do_out

fill: FORALL (i=1:ns, j=1:ns,  (i/=j) .AND. (dist(xyz(i,:),xyz(j,:))<(r(i)+r(j))) .AND. (i<j) )
	a(i,j)=1
END FORALL fill



!La somma mi da quante sono le coppie di sfere che si intersecano
somma=SUM(a)

somma_if: IF (somma/=0) THEN

	error=1
	WRITE(*,10) somma
	10 FORMAT ("Ci sono ",I5, " coppie di sfere che si intersecano.")
	WRITE(*,20)
	20 FORMAT ("Se vuoi vedere le coppie (tanto tempo...) scrivi 1, se no 0: ", $)
	READ (*,*) flag

	flag_if:IF (flag==1) THEN

		p=>a

		out_do:	DO i=1,ns
			in_do:	DO j=i,ns
				IF (p(i,j)==1) THEN
					WRITE(*,30) i,j,dist(xyz(i,:),xyz(j,:)),r(i)+r(j)
					30 FORMAT("Coppia", 2I5,2ES15.6)
				END IF
			END DO in_do
		END DO out_do

		NULLIFY(p)

	END IF flag_if

ELSE

	error=0

END IF somma_if

END SUBROUTINE intspheres


!******************************************************************************
!******************************************************************************
!******************************************************************************
!2bis) SUBROUTINE INTSHELLS: per vedere se ci sono layers intersecantisi
!******************************************************************************
!******************************************************************************
!******************************************************************************

SUBROUTINE intshells(r,xyz,ns,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(:),INTENT(IN) :: r				! Array raggi della sfera e della shell
REAL(dbl), DIMENSION(:,:),INTENT(IN) :: xyz 			! Array posizione sfere
INTEGER(lo), INTENT(IN) :: ns					! Numero di sfere
INTEGER(lo), INTENT(OUT) :: error					! Errore

! Dichiarazione variabili interne
INTEGER(lo), DIMENSION(ns), TARGET ::	a			! Array di controllo
INTEGER(lo) :: i,j,somma,flag					! Indici somma e flag
INTEGER(lo), DIMENSION(:), POINTER :: p=>NULL()		! Puntatore su a

! Inizio della procedura vera e propria
a=0

fill_do: DO i=1,ns-1
     IF ( (dist(xyz(i,:),xyz(i+1,:))+r(i) ) >= r(i+1) ) THEN
          a(i)=1
     END IF
END DO fill_do


!La somma mi da quanti sono i layer non concentrici
somma=SUM(a)

somma_if:	IF (somma/=0) THEN

	error=1
	WRITE(*,10) somma
	10 FORMAT ("Ci sono ",I5, " layer non concentrici.")
	WRITE(*,20)
	20 FORMAT ("Se vuoi vedere quali layer sono scrivi 1, se no 0: ", $)
	READ (*,*) flag

	flag_if:	IF (flag==1) THEN

		p=>a

		scan_do:	DO i=1,ns
			IF (p(i)==1) THEN
				WRITE(*,30) i
				30 FORMAT("Layer", 2I5)
			END IF
		END DO scan_do

		NULLIFY(p)

	END IF flag_if

	ELSE

	error=0

END IF somma_if

END SUBROUTINE intshells


!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!3) SUBROUTINE FILL: farcisco un mucchio di array con le informazioni che mi servono
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE fillarray(unit_in,ns,neps,neq,h_patt,h_eps,h_req,h_epseq,v_r,v_req,m_xyz,v_unit_in,v_patt,v_eps,v_epsg,v_epseq,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
INTEGER(lo), INTENT(IN) :: unit_in, ns, neps,neq					! Unita' apertura file principale, num tot sfere, num materiali e num uguali
TYPE (NODE_long_INT), POINTER :: h_patt							! Linked List Pattern sfere uguali
TYPE (NODE_CHAR), POINTER :: h_eps,h_epseq						! LL per le funzione dielettriche e per sfere uguali
TYPE(NODE_DBL_REAL), POINTER :: h_req							! LL per i raggi delle sfere uguali

INTEGER(lo), INTENT(OUT) :: error								! errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_r,v_req					! Array raggi sfere e per sfere uguali
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_xyz						! Array coordinate cartesiane sfere
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_unit_in,v_patt			! Array unita' di apertura materiali e pattern uguaglianza
CHARACTER(len=*), DIMENSION(:), INTENT(OUT) :: v_eps,v_epsg,v_epseq		! Array nomi funzioni dielettriche, Array tutte le funzioni dielettriche

! Dichiarazione variabili interne
INTEGER(lo) :: status_in,i		 							! Variabile di status per leggere il file, indice, errore
TYPE (NODE_long_INT), POINTER :: curs_patt=>NULL()					! Linked List Pattern sfere uguali
TYPE (NODE_CHAR), POINTER :: curs_eps=>NULL()						! LL per le funzione dielettriche
TYPE (NODE_DBL_REAL), POINTER :: curs_req=>NULL()					! LL per raggi sfere uguali
TYPE (NODE_CHAR), POINTER :: curs_epseq=>NULL()						! LL per funzioni dielettriche sfere uguali
! Inizio della procedura vera e propria

! Riempio gli array
read_do: DO i=1,ns
			READ(unit_in,*,IOSTAT=status_in) m_xyz(i,1),m_xyz(i,2),m_xyz(i,3),v_r(i),v_epsg(i)
END DO read_do

! Riempio l'array uni_in
unit_do: DO i=1,neps
			v_unit_in(i)=unit_in+i
END DO unit_do


! Riavvolgo il file
REWIND(UNIT=unit_in)


! Riempio l'array v_patt da h_patt e distruggo h_patt
curs_patt=>h_patt					! Vettore di pattern uguaglianza

fill_do1: DO i=0,ns-1
			v_patt(ns-i)=curs_patt%value
			curs_patt=>curs_patt%next
		  END DO fill_do1

NULLIFY(curs_patt)				! Lo nullifichiamo per sicurezza

CALL lldestructor(h_patt,error)

! lldescructor funziona?
patt_if: IF (error==1) THEN
			WRITE (*,*)
			WRITE (*,10)
			10 FORMAT ("Si e' verificato un errore nella routine lldestructor")
			RETURN
		ELSE
			error=0
        END IF patt_if

! Riempio l'array v_eps da h_eps e distruggo h_eps
curs_eps=>h_eps						! Vettore nomi funzioni dielettriche

fill_do2: DO i=0,neps-1
			v_eps(neps-i)=curs_eps%value
			curs_eps=>curs_eps%next
		  END DO fill_do2

NULLIFY(curs_eps)					! Lo nullifichiamo per sicurezza

CALL lldestructor(h_eps,error)

! lldescructor funziona?
eps_if: IF (error==1) THEN
			WRITE (*,*)
			WRITE (*,20)
			20 FORMAT ("Si e' verificato un errore nella routine lldestructor")
			RETURN
		ELSE
			error=0
        END IF eps_if

! Riempio l'array v_req da h_req e distruggo h_req
curs_req=>h_req						! Vettore nomi funzioni dielettriche

fill_do3: DO i=0,neq-1
			v_req(neq-i)=curs_req%value
			curs_req=>curs_req%next
		  END DO fill_do3

NULLIFY(curs_req)					! Lo nullifichiamo per sicurezza

CALL lldestructor(h_req,error)

! lldescructor funziona?
req_if: IF (error==1) THEN
			WRITE (*,*)
			WRITE (*,30)
			30 FORMAT ("Si e' verificato un errore nella routine lldestructor")
			RETURN
		ELSE
			error=0
        END IF req_if

! Riempio l'array v_epseq da h_epseq e distruggo h_epseq
curs_epseq=>h_epseq						! Vettore nomi funzioni dielettriche

fill_do4: DO i=0,neq-1
			v_epseq(neq-i)=curs_epseq%value
			curs_epseq=>curs_epseq%next
		  END DO fill_do4

NULLIFY(curs_epseq)					! Lo nullifichiamo per sicurezza

CALL lldestructor(h_epseq,error)

! lldescructor funziona?
reps_if: IF (error==1) THEN
			WRITE (*,*)
			WRITE (*,40)
			40 FORMAT ("Si e' verificato un errore nella routine lldestructor")
			RETURN
		ELSE
			error=0
        END IF reps_if


END SUBROUTINE fillarray








!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!3) SUBROUTINE FILLARRAY_SHELL: farcisco un mucchio di array con le informazioni che mi servono, nel nuovo caso shell
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE fillarray_shell(unit_in,ns,neps,neq,h_patt,h_eps,h_req,h_epseq,h_req_shell,h_epseq_shell,v_r,v_req,v_r_shell,&
				&  v_req_shell,m_xyz,v_unit_in,v_patt,v_eps,v_epsg,v_epseq,v_epsg_shell,v_epseq_shell,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
INTEGER(lo), INTENT(IN) :: unit_in, ns, neps,neq					! Unita' apertura file principale, num tot sfere, num materiali e num uguali
TYPE (NODE_long_INT), POINTER :: h_patt							! Linked List Pattern sfere uguali
TYPE (NODE_CHAR), POINTER :: h_eps,h_epseq,h_epseq_shell				! LL per le funzione dielettriche e per sfere uguali
TYPE(NODE_DBL_REAL), POINTER :: h_req,h_req_shell					! LL per i raggi delle sfere uguali

INTEGER(lo), INTENT(OUT) :: error								! errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_r,v_req,v_r_shell,v_req_shell		! Array raggi sfere e per sfere uguali
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_xyz						! Array coordinate cartesiane sfere
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_unit_in,v_patt			! Array unita' di apertura materiali e pattern uguaglianza
CHARACTER(len=*), DIMENSION(:), INTENT(OUT) :: v_eps,v_epsg,v_epseq		! Array nomi funzioni dielettriche, Array tutte le funzioni dielettriche
CHARACTER(len=*), DIMENSION(:), INTENT(OUT) :: v_epsg_shell,v_epseq_shell	! Array nomi funzioni dielettriche, Array tutte le funzioni dielettriche

! Dichiarazione variabili interne
INTEGER(lo) :: status_in,i		 							! Variabile di status per leggere il file, indice, errore
TYPE (NODE_long_INT), POINTER :: curs_patt=>NULL()					! Linked List Pattern sfere uguali
TYPE (NODE_CHAR), POINTER :: curs_eps=>NULL()						! LL per le funzione dielettriche
TYPE (NODE_DBL_REAL), POINTER :: curs_req=>NULL(),curs_req_shell=>NULL()	! LL per raggi sfere uguali
TYPE (NODE_CHAR), POINTER :: curs_epseq=>NULL(),curs_epseq_shell=>NULL()! LL per funzioni dielettriche sfere uguali
! Inizio della procedura vera e propria

! Riempio gli array
read_do: DO i=1,ns
			READ(unit_in,*,IOSTAT=status_in) m_xyz(i,1),m_xyz(i,2),m_xyz(i,3),v_r(i),v_r_shell(i),v_epsg(i),v_epsg_shell(i)
END DO read_do

! Riempio l'array uni_in
unit_do: DO i=1,neps
			v_unit_in(i)=unit_in+i
END DO unit_do


! Riavvolgo il file
REWIND(UNIT=unit_in)




! Riempio l'array v_patt da h_patt e distruggo h_patt
curs_patt=>h_patt							! Vettore di pattern uguaglianza

fill_do1: DO i=0,ns-1
	v_patt(ns-i)=curs_patt%value
	curs_patt=>curs_patt%next
END DO fill_do1

NULLIFY(curs_patt)						! Lo nullifichiamo per sicurezza

CALL lldestructor(h_patt,error)

! lldescructor funziona?
patt_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,10)
	10 FORMAT ("Si e' verificato un errore nella routine lldestructor")
	RETURN
ELSE
	error=0
END IF patt_if





! Riempio l'array v_eps da h_eps e distruggo h_eps
curs_eps=>h_eps							! Vettore nomi funzioni dielettriche

fill_do2: DO i=0,neps-1
	v_eps(neps-i)=curs_eps%value
	curs_eps=>curs_eps%next
END DO fill_do2

NULLIFY(curs_eps)							! Lo nullifichiamo per sicurezza

CALL lldestructor(h_eps,error)

! lldescructor funziona?
eps_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,20)
	20 FORMAT ("Si e' verificato un errore nella routine lldestructor")
	RETURN
ELSE
	error=0
END IF eps_if




! Riempio l'array v_req da h_req e distruggo h_req
curs_req=>h_req							! Vettore nomi funzioni dielettriche

fill_do3: DO i=0,neq-1
	v_req(neq-i)=curs_req%value
	curs_req=>curs_req%next
END DO fill_do3

NULLIFY(curs_req)							! Lo nullifichiamo per sicurezza

CALL lldestructor(h_req,error)

! lldescructor funziona?
req_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,30)
	30 FORMAT ("Si e' verificato un errore nella routine lldestructor")
	RETURN
ELSE
	error=0
END IF req_if




! Riempio l'array v_epseq da h_epseq e distruggo h_epseq
curs_epseq=>h_epseq						! Vettore nomi funzioni dielettriche

fill_do4: DO i=0,neq-1
	v_epseq(neq-i)=curs_epseq%value
	curs_epseq=>curs_epseq%next
END DO fill_do4

NULLIFY(curs_epseq)						! Lo nullifichiamo per sicurezza

CALL lldestructor(h_epseq,error)

! lldescructor funziona?
reps_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,40)
	40 FORMAT ("Si e' verificato un errore nella routine lldestructor")
	RETURN
ELSE
	error=0
END IF reps_if




! Riempio l'array v_req_shell da h_req_shell e distruggo h_req_shell
curs_req_shell=>h_req_shell					! Vettore nomi funzioni dielettriche

fill_do5: DO i=0,neq-1
	v_req_shell(neq-i)=curs_req_shell%value
	curs_req_shell=>curs_req_shell%next
END DO fill_do5

NULLIFY(curs_req_shell)						! Lo nullifichiamo per sicurezza

CALL lldestructor(h_req_shell,error)

! lldescructor funziona?
req_shell_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,31)
	31 FORMAT ("Si e' verificato un errore nella routine lldestructor")
	RETURN
ELSE
	error=0
END IF req_shell_if




! Riempio l'array v_epseq_shell da h_epseq_shell e distruggo h_epseq_shell
curs_epseq_shell=>h_epseq_shell				! Vettore nomi funzioni dielettriche

fill_do6: DO i=0,neq-1
	v_epseq_shell(neq-i)=curs_epseq_shell%value
	curs_epseq_shell=>curs_epseq_shell%next
END DO fill_do6

NULLIFY(curs_epseq_shell)					! Lo nullifichiamo per sicurezza

CALL lldestructor(h_epseq_shell,error)

! lldescructor funziona?
reps_shell_if: IF (error==1) THEN
	WRITE (*,*)
	WRITE (*,41)
	41 FORMAT ("Si e' verificato un errore nella routine lldestructor")
	RETURN
ELSE
	error=0
END IF reps_shell_if


END SUBROUTINE fillarray_shell












!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!4) SUBROUTINE checkepsilon: controllo le tabulazioni delle funzioni dielettriche
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE checkepsilon(v_unit_in,neps,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
INTEGER(lo), INTENT(IN) :: neps								! Numero di funzioni dielettriche
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_unit_in			! Unita' di apertura
INTEGER(lo), INTENT(OUT) ::error								! Errore

! Dichiarazione variabili interne
REAL(dbl) :: eps1,eps2								! Var spazzatura per epsilon
REAL(dbl), DIMENSION(neps) :: v_lambda				! Vettore per tenere le lunghezze d'onda
INTEGER(lo) :: status_in,i						! Variabile di status lettura, indice
INTEGER(lo),DIMENSION(neps) :: ncount			! Array per vedere se la lunghezza dei file e' omogenea


! Inizio della procedura vera e propria
ncount=0

! Ciclo per vedere le lunghezze delle tabulazioni delle funzioni dielettriche
outer_do1: DO i=1,neps

			inner_do1: DO

 				READ (v_unit_in(i),*,IOSTAT=status_in) v_lambda(i),eps1,eps2
 				IF (status_in<0) CYCLE outer_do1
 				ncount(i)=ncount(i)+1

 			END DO inner_do1

END DO outer_do1

ncount=ncount-ncount(1)

!Le lunghezze sono uguali?
IF (DOT_PRODUCT(ncount,ncount)==0) THEN
	error=0
ELSE
	error=1
	RETURN
END IF

!Loop per il rewind dei file
rewind_do1: DO i=1,neps
				REWIND(v_unit_in(i))
END DO rewind_do1

! Mando avanti di uno perche' qui ho i parametri
forward_do: DO i=1,neps
				READ (v_unit_in(i),*,IOSTAT=status_in) v_lambda(i),eps1,eps2
END DO forward_do

! Ciclo per controllare che il passo delle lunghezze d'onda sia uguale
outer_do2: DO

				inner_do2: DO i=1,neps
								READ (v_unit_in(i),*,IOSTAT=status_in) v_lambda(i),eps1,eps2
				END DO inner_do2

				IF (status_in<0) EXIT
				v_lambda=v_lambda-v_lambda(1)

				lambda_if: IF (DOT_PRODUCT(v_lambda,v_lambda)<=1.0D-10) THEN
							error=0
						   ELSE
							error=2
				END IF lambda_if

				IF (error==2) EXIT

END DO outer_do2

!Loop per il rewind dei file
rewind_do2: DO i=1,neps
				REWIND(v_unit_in(i))
END DO rewind_do2

END SUBROUTINE checkepsilon


!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!5) SUBROUTINE correps: correzione della funzione dielettrica
!   Cambiamento 04/06/2009: da gamma=gamma_inf + par(3)/(r*1e-9) a gamma=gamma_inf + par(3)/(4.0D0*r*1e-9/3.0D0)
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE correps(lambda,e1,e2,r,par)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,r 				! Raggio della particella
REAL(dbl), DIMENSION(:), INTENT(IN) ::	par		! Array parametri
REAL(dbl), INTENT(INOUT) :: e1,e2				! Lunghezza d'onda, parte reale immaginaria funzione dielettrica

! Dichiarazione variabili interne
REAL(dbl) :: omega_drude,gamma_inf,gamma,omega,leff						! Parametri correzioni
REAL(dbl) :: corr1,corr2										! Correzioni funzione dielettrica

! Inizio della procedura vera e propria

! Parametri necessari alla correzione
leff=(4.0D0*r*1.0D-9)/3.0D0
omega_drude=SQRT((par(1)*(e**2))/(e0*me))
gamma_inf=1/par(2)
gamma=gamma_inf + par(3)/(leff)

! Correzione della funzione dielettrica complessa
omega=(2*pi_d*c)/(lambda*1e-9)
corr1=(omega_drude**2) * ( 1.0/(omega**2 + gamma_inf**2) - 1.0/(omega**2 + gamma**2) )
corr2=( (omega_drude**2) / omega ) * ( -gamma_inf/(omega**2 + gamma_inf**2) + gamma/(omega**2 + gamma**2) )

e1=e1+corr1
e2=e2+corr2

END SUBROUTINE correps






!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!5) SUBROUTINE correps: correzione della funzione dielettrica
!   Cambiamento 04/06/2009: da gamma=gamma_inf + par(3)/(r*1e-9) a gamma=gamma_inf + par(3)/(4.0D0*r*1e-9/3.0D0)
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE correps_shell(lambda,e1,e2,r,r_shell,par)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,r,r_shell			! Raggio della particella
REAL(dbl), DIMENSION(:), INTENT(IN) ::par			! Array parametri
REAL(dbl), INTENT(INOUT) :: e1,e2				! Lunghezza d'onda, parte reale immaginaria funzione dielettrica

! Dichiarazione variabili interne
REAL(dbl) :: omega_drude,gamma_inf,gamma,omega,leff,q						! Parametri correzioni
REAL(dbl) :: corr1,corr2										! Correzioni funzione dielettrica

! Inizio della procedura vera e propria

! Parametri necessari alla correzione
q=r/r_shell
leff=((4.0D0*r_shell*1.0D-9)/3.0D0)*(1.0D0-q**3)/(1.0D0+q**2)
omega_drude=SQRT((par(1)*(e**2))/(e0*me))
gamma_inf=1/par(2)

!Con questo if, risolvo il problema della correzione del caso di una particella solo core
IF (r==r_shell) THEN
	gamma=gamma_inf
ELSE
	gamma=gamma_inf + par(3)/(leff)
END IF

! Correzione della funzione dielettrica complessa
omega=(2*pi_d*c)/(lambda*1e-9)
corr1=(omega_drude**2) * ( 1.0/(omega**2 + gamma_inf**2) - 1.0/(omega**2 + gamma**2) )
corr2=( (omega_drude**2) / omega ) * ( -gamma_inf/(omega**2 + gamma_inf**2) + gamma/(omega**2 + gamma**2) )

e1=e1+corr1
e2=e2+corr2

END SUBROUTINE correps_shell



!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!6) SUBROUTINE lambdaepsilon: per ogni lunghezza d'onda costruisco due vettori per la parte reale e per la parte
!							  immaginaria della funzione dielettrica che mi serve, e cosi' poi posso fare i miei
!							  calcoli in relativa tranquillita'
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE lambdaepsilon(lambda,neps,neq,v_eps,v_epseq,v_req,m_eps,m_par,correction,m_epseq)

IMPLICIT NONE

! Dichiarazione dei dummy argument

INTEGER(lo), INTENT(IN) :: neps,neq							! Numero materiali e numero gruppi di sfere uguali
REAL(dbl), INTENT(IN) :: lambda									! Lunghezza d'onda
CHARACTER(len=*), INTENT(IN) :: correction						! Flag correzione
CHARACTER(len=*) ,DIMENSION(:), INTENT(IN) :: v_eps,v_epseq	! Array nomi funct diel
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req					! Array raggi sfere uguali
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_eps,m_par			! Matrici numeriche funzioni dielettriche e parametri

REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_epseq				! Matrice num funct diel corr sfere uguali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j											! Indici
REAL(dbl), DIMENSION(neq,3) :: m_pareq							! Matrice par eq

! Inizio della procedura vera e propria

out_do: DO i=1,neq
			in_do: DO j=1,neps

					IF (v_eps(j) /= v_epseq(i)) CYCLE in_do
					m_epseq(i,1:2)=m_eps(j,1:2)
					m_pareq(i,1:3)=m_par(j,1:3)

			END DO in_do
END DO out_do

!Non correggo la funzione dielettrica ed esco
IF (correction=='no') THEN
	RETURN
END IF

corr_do: DO i=1,neq
			CALL correps(lambda,m_epseq(i,1),m_epseq(i,2),v_req(i),m_pareq(i,:))
END DO corr_do

END SUBROUTINE lambdaepsilon



!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!6bis) SUBROUTINE lambdaepsilon_singleshell: per ogni lunghezza d'onda costruisco due vettori per la parte reale e per la parte
!							  immaginaria della funzione dielettrica che mi serve, e cosi' poi posso fare i miei
!							  calcoli in relativa tranquillita'
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE lambdaepsilon_singleshell(lambda,neps,neq,v_eps,v_epseq,v_req,m_eps,m_par,correction,m_epseq)

IMPLICIT NONE

! Dichiarazione dei dummy argument

INTEGER(lo), INTENT(IN) :: neps,neq							! Numero materiali e numero gruppi di sfere uguali
REAL(dbl), INTENT(IN) :: lambda									! Lunghezza d'onda
CHARACTER(len=*), INTENT(IN) :: correction						! Flag correzione
CHARACTER(len=*) ,DIMENSION(:), INTENT(IN) :: v_eps,v_epseq	! Array nomi funct diel
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req					! Array raggi sfere uguali
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_eps,m_par			! Matrici numeriche funzioni dielettriche e parametri

REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_epseq				! Matrice num funct diel corr sfere uguali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j											! Indici
REAL(dbl), DIMENSION(neq,3) :: m_pareq							! Matrice par eq

! Inizio della procedura vera e propria

out_do: DO i=1,neq
			in_do: DO j=1,neps

					IF (v_eps(j) /= v_epseq(i)) CYCLE in_do
					m_epseq(i,1:2)=m_eps(j,1:2)
					m_pareq(i,1:3)=m_par(j,1:3)

			END DO in_do
END DO out_do

!Non correggo la funzione dielettrica ed esco
IF (correction=='no') THEN
	RETURN
END IF

corr_do: DO i=1,neq

	IF (i==1) THEN
		CALL correps(lambda,m_epseq(i,1),m_epseq(i,2),v_req(i),m_pareq(i,:))
	ELSE IF (i==2) THEN
		CALL correps_shell(lambda,m_epseq(i,1),m_epseq(i,2),v_req(i-1),v_req(i),m_pareq(i,:))
	END IF

END DO corr_do

END SUBROUTINE lambdaepsilon_singleshell



!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************
!6tris) SUBROUTINE lambdaepsilon_shell: per ogni lunghezza d'onda costruisco due vettori per la parte reale e per la parte
!immaginaria della funzione dielettrica che mi serve, e cosi' poi posso fare i miei calcoli in relativa tranquillita',
!in piu correggo per le dimensioni, sia per il core sia per la shell, con mean free path diversi
!**********************************************************************************************************************
!**********************************************************************************************************************
!**********************************************************************************************************************

SUBROUTINE lambdaepsilon_shell(lambda,neps,neq,v_eps,v_epseq,v_epseq_shell,v_req,v_req_shell,m_eps,m_par,correction,&
					&correction_shell,m_epseq,m_epseq_shell)

IMPLICIT NONE

! Dichiarazione dei dummy argument

INTEGER(lo), INTENT(IN) :: neps,neq							! Numero materiali e numero gruppi di sfere uguali
REAL(dbl), INTENT(IN) :: lambda								! Lunghezza d'onda
CHARACTER(len=*), INTENT(IN) :: correction,correction_shell				! Flag correzione
CHARACTER(len=*) ,DIMENSION(:), INTENT(IN) :: v_eps,v_epseq,v_epseq_shell	! Array nomi funct diel
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req,v_req_shell				! Array raggi sfere uguali
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_eps,m_par					! Matrici numeriche funzioni dielettriche e parametri

REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_epseq,m_epseq_shell			! Matrice num funct diel corr sfere uguali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j											! Indici
REAL(dbl), DIMENSION(neq,3) :: m_pareq,m_pareq_shell						! Matrice par eq

! Inizio della procedura vera e propria


!Piazzo le funzioni dielettriche per il core
out_do: DO i=1,neq
			in_do: DO j=1,neps

					IF (v_eps(j) /= v_epseq(i)) CYCLE in_do
					m_epseq(i,1:2)=m_eps(j,1:2)
					m_pareq(i,1:3)=m_par(j,1:3)

			END DO in_do
END DO out_do


!Piazzo le funzioni dielettriche per la shell
out_do_shell: DO i=1,neq
			in_do_shell: DO j=1,neps

					IF (v_eps(j) /= v_epseq_shell(i)) CYCLE in_do_shell
					m_epseq_shell(i,1:2)=m_eps(j,1:2)
					m_pareq_shell(i,1:3)=m_par(j,1:3)

			END DO in_do_shell
END DO out_do_shell


!Secondo la flag correggo la funzione dielettrica oppure no
IF (correction=='yes') THEN
	corr_do: DO i=1,neq
		CALL correps(lambda,m_epseq(i,1),m_epseq(i,2),v_req(i),m_pareq(i,:))
	END DO corr_do
END IF

IF (correction_shell=='yes') THEN
	corr_do_shell: DO i=1,neq
		CALL correps_shell(lambda,m_epseq(i,1),m_epseq(i,2),v_req(i),v_req_shell(i),m_pareq_shell(i,:))
	END DO corr_do_shell
END IF

END SUBROUTINE lambdaepsilon_shell


END MODULE gmmsubs
