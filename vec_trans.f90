MODULE vec_trans

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! WIGNER 3JM SUBROUTINES INDEX
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines
! 1) j3minimum(j1,j2,m1,m2) f
!		COMMENT: finds the minimum value to start the w3jm recursion
!	
! 2) nw(j1,j2,m1,m2) f
!		COMMENT: number of calculated terms in the recursion
!
! 3) wdown0(j1,j2,m1,m2) f
!		COMMENT: calculates the first term of the backward recursion
!
! 4) wup0(j1,j2,m1,m2) f
!		COMMENT: calculates the first term of the upward recursion
!
! 5) cr(j1,j2,m1,m2) f
!		COMMENT: calculates the first coefficent for the recursion formula
!
! 6) dr(j1,j2,m1,m2) f
!		COMMENT: calculates the second coefficent for the recursion formula
!			
! 7) wigner3jm(j1,j2,m1,m2,j3min,j3max,v_w3jm) s
!		COMMENT: calculates the vector containing the values of 3jm simbols
!			
!		INPUT: REAL(dbl): j1,j2,m1,m2,
!		       INTEGER(lo): j3min,j3max
!
!		OUTPUT: REAL(DBL) 1-dim array: v_w3jm
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! GAUNT GENERAL FUNCTION AND SUBROUTINE INDEX
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines and Functions
! 1) idx(m,n,mu,nu,nstop) f
!		COMMENT: if the coefficient matrix is tought as a linear vector, this function calculates the the position
!				 the current value in the linear vector from indexes m,n,mu,nu and number of multipolar expansions 
!				 nstop
!	
! 2) c0sub(nstop,v_c0) s
!		COMMENT: given nstop, it calculates a 1-dim vector containing c0 coefficients for the calculations of vector
!				 translation coefficients
!
!		INPUT: INTEGER(lo): nstop
!
!		OUTPUT: REAL(dbl): 1-dim array: v_c0
!
!
! 2bis) c0xusub(nstop,v_c0xu) s
!		COMMENT: given nstop, it calculates a 1-dim vector containing c0 coefficients for the calculations of vector
!				 translation coefficients, with the last xu normalization
!
!		INPUT: INTEGER(lo): nstop
!
!		OUTPUT: REAL(dbl): 1-dim array: v_c0xu
!
!
!
! 3) qmaxsub(nstop,v_qmax) s
!		COMMENT: calculates qmax for gaunt and vector translation coefficients
!
!		INPUT: INTEGER(lo): nstop
!
!		OUTPUT: REAL(dbl): 1-dim array: v_qmax
!
! 4) qqmaxsub(nstop,v_qqmax) s
!		COMMENT: calculates Qmax for bq and vector translation coefficients
!
!		INPUT: INTEGER(lo): nstop
!
!		OUTPUT: REAL(dbl): 1-dim array: v_qqmax
!
! 5) idg(v_qmaxsum,nstop,m,n,mu,nu,q) f
!		COMMENT: calculates the index in order to recover a given gaunt coefficients from its storage vector
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! GAUNT COEFFICIENTS SUBROUTINES (CRUZAN) INDEX
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines and Functions
! 1) gaunt_cz(m,n,mu,nu,qmaxa,v_aq,error) s
!		COMMENT: calculates a series of gaunt coefficients for a quadruplet of given indexes
!	
!		INPUT: REAL(dbl): m,n,mu,nu
!			   INTEGER(lo): qmaxa
!
!		OUTPUT: REAL(dbl): 1-dim array: v_aq
!				INTEGER(lo): error
!
! 2) aqbq_cz(m,n,mu,nu,qmaxa,v_aq,qmaxb,v_bq,error) s
!		COMMENT: calculates a series of gaunt and b coefficients for a quadruplet of given indexes
!
!		INPUT: REAL(dbl): m,n,mu,nu
!			   INTEGER(lo): qmaxa,qmaxb
!
!		OUTPUT: REAL(dbl): 1-dim array: v_aq,v_bq
!				INTEGER(lo): error
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! GAUNT COEFFICIENTS SUBROUTINES (XU) INDEX
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Functions and Subroutines
!
! 1) f_alpha(n,nu,p) f
!		COMMENT: calculates the alpha parameters for Xu calculations
!
! 2) f_Ap(m,n,mu,nu,p) f
!		COMMENT: calculates the Ap parameters for Xu calculations
!
! 3) f_a0(m,n,mu,nu) f
!		COMMENT: starting value v_aq(q=0) for the backward recursion of gaunt coefficients
!
! 4) f_a1norm(m,n,mu,nu) f
!		COMMENT: normalized starting value v_aq(q=1) for the backward recursion of gaunt coefficients
!
! 5) f_a2norm(m,n,mu,nu) f
!		COMMENT: normalized starting value v_aq(q=2) for the backward recursion of gaunt coefficients
!
! 6) f_a2normr(m,n,mu,nu,a1norm) f
!		COMMENT: normalized starting value v_aq(q=2) for the backward recursion of gaunt coefficients
!				 calculated by recurrence (so a1norm is needed as input)
!
! 7) f_aqmax(m,n,mu,nu,qmax) f
!		COMMENT: starting value v_aq(q=qmax) for the forward recursion of gaunt coefficients
!				
! 8) f_aqmax_1(m,n,mu,nu,qmax) f
!		COMMENT: starting value v_aq(q=qmax) for the forward recursion of gaunt coefficients
!
! 9) gaunt_xu(m,n,mu,nu,qmax,v_aq,error) s
!		COMMENT: calculates a series of gaunt coefficients for a quadruplet of given indexes
!
!		INPUT: REAL(dbl): m,n,mu,nu
!			   INTEGER(lo): qmax
!
!		OUTPUT: REAL(dbl) 1-dim array: v_aq
!				INTEGER(lo): error
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! LEGENDRE SUBROUTINES INDEX
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines and Functions
! 1) legendre(nmin,nmax,m,mm,theta,v_leg,error) 
!		COMMENT: calculates a series of legendre functions in cos(theta) from nmin to nmax with fixed m
!
!		INPUT: INTEGER(lo): nmin,nmax,m,mm
!			   REAL(dbl): theta
!
!		OUTPUT: REAL(dbl): 1-dim array, v_leg
!				INTEGER(lo): error
!
! 2) pi_mn(nstop,theta,v_pimn,error) 
!		COMMENT: calculates a series of angular functions Pi_mn in cos(theta) n=1,nstop; m=-n,n
!
!		INPUT: INTEGER(lo): nstop
!			   REAL(dbl): theta
!
!		OUTPUT: REAL(dbl): 1-dim array, v_pimn
!				INTEGER(lo): error
!
!
! 3) tau_mn(nstop,theta,v_pimn,v_taumn,error) 
!		COMMENT: calculates a series of angular functions tau_mn in cos(theta) n=1,nstop; m=-n,n
!
!		INPUT: INTEGER(lo): nstop
!			   REAL(dbl): theta
!			   REAL(dbl): 1-dim array, v_pimn
!
!		OUTPUT: REAL(dbl): 1-dim array, v_taumn
!				INTEGER(lo): error
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! VECTOR TRANSLATION COEFFICIENTS SUBROUTINES INDEX
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines and Functions
! 1) vectrans_cz(m,n,mu,nu,kr,theta,phi,Avt,Bvt,error)
!		COMMENT: calculates vector translation coefficients Avt Bvt for given indexes and
!				 coordinates
!
!		INPUT: INTEGER(lo): n,m,mu,nu
!			   REAL(dbl): kr,theta,phi
!
!		OUTPUT: COMPLEX(dbl): Avt,Bvt
!				INTEGER(lo): error
!
! 2) vectrans_xu(m,n,mu,nu,kr,theta,phi,Avt,Bvt,error)
!		COMMENT: calculates vector translation coefficients Avt Bvt for given indexes and
!				 coordinates
!
!		INPUT: INTEGER(lo): n,m,mu,nu
!			   REAL(dbl): kr,theta,phi
!
!		OUTPUT: COMPLEX(dbl): Avt,Bvt
!				INTEGER(lo): error
!
! 3) fillblock_sub(nstop,v_qmax,v_qqmax,v_c0,v_aq_long,v_bq_long,kr,theta,phi,m_Aij,m_Bij,error)
!
!		INPUT: INTEGER(lo): nstop
!			   REAL(dbl): kr,theta,phi
!			   INTEGER(lo): 1-dim array: v_qmax,v_qqmax
!			   REAL(dbl): 1-dim array: v_c0,v_aq_long,v_bq_long
!
!		OUTPUT: COMPLEX(dbl): 2-dim array: m_Aij,m_Bij
!				INTEGER(lo): error
!
!
! 4) fill_AB_sub(nstop,k,v_patt,m_a,m_b,m_xyz,m_AB,error)
!
!		INPUT:	INTEGER(lo): nstop
!				REAL(dbl): k
!				INTEGER(lo): 1-dim array: v_patt
!				REAL(dbl): 2-dim array: m_xyz
!				COMPLEX(dbl): 2-dim: m_a,m_b
!
!		OUTPUT:	COMPLEX(dbl): 2-dim array: m_AB
!
! 5) field_exp_sub(nstop,ns,v_z,m_a,m_b,v_patt,v_pmn,v_qmn,v_p)
! 
! 		INPUT: 	INTEGER(lo): nstop,ns
!				INTEGER(lo): 1-dim array: v_patt
!				REAL(dbl): 1-dim array: v_z
!				COMPLEX(dbl): 2-dim array: m_a,m_b
!
!		OUTPUT:	COMPLEX(dbl): 1-dim array: v_pmn,v_qmn,v_p 
! 
!
! 6) cext_sub(nstop,ns,k,v_z,v_ab,v_cext)
! 
! 		INPUT:	INTEGER(lo): nstop,ns
!				REAL(dbl): k
!				REAL(dbl): 1-dim array: v_z
!				COMPLEX(dbl): 1-dim array: v_ab
!
!		OUTPUT:	REAL(dbl): 1-dim array: v_cext
! 
!
! 7) cext_sub_exp(nstop,ns,k,v_z,v_ab,m_cext)
! 
! 		INPUT:	INTEGER(lo): nstop,ns
!				REAL(dbl): k
!				REAL(dbl): 1-dim array: v_z
!				COMPLEX(dbl): 1-dim array: v_ab
!
!		OUTPUT:	REAL(dbl): 2-dim array: m_cext
! 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


USE kinds
USE basicsubs
USE sing_part
USE shared_data


CONTAINS
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! WIGNER 3JM SUBROUTINES
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!******************************************************************************
!1) FUNCTION j3min: per calcolare il minimo indice per i coeff di wigner
!******************************************************************************
FUNCTION j3minimum(j1,j2,m1,m2)

IMPLICIT NONE

! Dichiarazione funzione
INTEGER(lo) :: j3minimum

! Dichiarazione Argomenti
REAL(dbl), INTENT(IN) :: j1,j2,m1,m2	! Parametri Wigner

! Funzione vera e propria
j3minimum=INT(MAX(ABS(j1-j2),ABS(m1+m2)),lo)

END FUNCTION j3minimum



!******************************************************************************
!2) FUNCTION nw: numero di termini nella ricorsione di wigner
!******************************************************************************
FUNCTION nw(j1,j2,m1,m2)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: nw


! Dichiarazione Argomenti
REAL(dbl), INTENT(IN) :: j1,j2,m1,m2	! Parametri Wigner

! Funzione vera e propria
nw=j1+j2+1.0D0-MAX(ABS(j1-j2),ABS(m1+m2))

END FUNCTION nw



!******************************************************************************
!3) FUNCTION wdown0: per calcolare il primo termine della backward recursion
!******************************************************************************
FUNCTION wdown0(j1,j2,m1,m2)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: wdown0

! Dichiarazione Argomenti
REAL(dbl), INTENT(IN) :: j1,j2,m1,m2	! Parametri Wigner


! Dichiarazione variabili ausiliarie
REAL(dbl) :: logw
INTEGER(lo) :: r

! Funzione vera e propria

logw=0.5D0*(lnf(2.0D0*j1,r)+lnf(2.0D0*j2,r)+lnf(j1+j2-m1-m2,r)+lnf(j1+j2+m1+m2,r) - &
		&  	lnf(2.0D0*j1+2.0D0*j2+1.0D0,r)-lnf(j1-m1,r)-lnf(j1+m1,r)-lnf(j2-m2,r)-lnf(j2+m2,r))
			
wdown0=((-1.0D0)**INT(j1+j2+m1+m2,lo))*EXP(logw)

END FUNCTION wdown0


!******************************************************************************
!4) FUNCTION wup0: per calcolare il primo termine della upward recursion
!******************************************************************************
FUNCTION wup0(j1,j2,m1,m2)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: wup0

! Dichiarazione Argomenti
REAL(dbl), INTENT(IN) :: j1,j2,m1,m2 ! Parametri Wigner


! Dichiarazione variabili ausiliarie
REAL(dbl) :: logw
INTEGER(lo) :: r

! Funzione vera e propria
wup0_if : IF ( (INT(j1-j2)>=0) .AND. ( ABS(INT(j1-j2))>=ABS(INT(m1+m2))) ) THEN
			
			logw=0.5D0 * ( lnf(j1-m1,r)+lnf(j1+m1,r)+lnf(2.0D0*j1-2.0D0*j2,r)+lnf(2.0D0*j2,r) - &
			& lnf(j2-m2,r)-lnf(j2+m2,r)-lnf(j1-j2-m1-m2,r)-lnf(j1-j2+m1+m2,r)-lnf(2.0D0*j1+1.0D0,r) )
			
			wup0=((-1.0D0)**INT(j1+m1))*EXP(logw)
			
		  ELSE IF ( (INT(j1-j2)<0) .AND. ( ABS(INT(j1-j2))>=ABS(INT(m1+m2))) ) THEN

			logw=0.5D0 * ( lnf(j2-m2,r)+lnf(j2+m2,r)+lnf(2.0D0*j2-2.0D0*j1,r)+lnf(2.0D0*j1,r) - &
			& lnf(j1-m1,r)-lnf(j1+m1,r)-lnf(j2-j1-m1-m2,r)-lnf(j2-j1+m1+m2,r)-lnf(2.0D0*j2+1.0D0,r) )			
			
			wup0=((-1.0D0)**INT(j2+m2))*EXP(logw)

		  ELSE IF ( (INT(m1+m2)>0) .AND. (ABS(INT(j1-j2))<ABS(INT(m1+m2))) ) THEN

			logw=0.5D0 * ( lnf(j1+m1,r)+lnf(j2+m2,r)+lnf(j1+j2-m1-m2,r)+lnf(2.0D0*m1+2.0D0*m2,r) - &
			& lnf(j1-m1,r)-lnf(j2-m2,r)-lnf(j1-j2+m1+m2,r)-lnf(j2-j1+m1+m2,r)-lnf(j1+j2+m1+m2+1.0D0,r) )
			
			wup0=((-1.0D0)**INT(j2+m2))*EXP(logw)

		  ELSE IF ( (INT(m1+m2)<0) .AND. (ABS(INT(j1-j2))<ABS(INT(m1+m2))) ) THEN

			logw=0.5D0 * ( lnf(j1-m1,r)+lnf(j2-m2,r)+lnf(j1+j2+m1+m2,r)+lnf(-2.0D0*m1-2.0D0*m2,r) - &
			& lnf(j1+m1,r)-lnf(j2+m2,r)-lnf(j1-j2-m1-m2,r)-lnf(j2-j1-m1-m2,r)-lnf(j1+j2-m1-m2+1.0D0,r) )
			
			wup0=((-1.0D0)**INT(j1+m1))*EXP(logw)
			
END IF wup0_if

END FUNCTION wup0

!******************************************************************************
!5) FUNCTION cr: coefficente per il calcolo ricorsivo
!******************************************************************************

FUNCTION cr(j1,j2,j3,m1,m2,m3)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: cr

! Dichiarazione argomenti
REAL(dbl) :: j1,j2,j3,m1,m2,m3 ! Parametri Wigner

! Funzione vera e propria
cr=SQRT((j3**2-(j1-j2)**2)*( (j1+j2+1.0D0)**2-j3**2 )*( j3**2-m3**2 ))

END FUNCTION cr


!******************************************************************************
!6) FUNCTION dr: coefficente per il calcolo ricorsivo
!******************************************************************************

FUNCTION dr(j1,j2,j3,m1,m2,m3)

IMPLICIT NONE

! Dichiarazione funzione
REAL(dbl) :: dr

! Dichiarazione argomenti
REAL(dbl) :: j1,j2,j3,m1,m2,m3 ! Parametri Wigner

! Funzione vera e propria
dr=-(2.0D0*j3+1.0D0)*( j1*(j1+1.0D0)*m3-j2*(j2+1.0D0)*m3-j3*(j3+1.0D0)*(m2-m1) )

END FUNCTION dr

!******************************************************************************
!7) subroutine Wigner3jm: calcolo il vettore di simboli 3jm
!******************************************************************************

SUBROUTINE wigner3jm(j1,j2,m1,m2,j3min,j3max,v_w3jm)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl), INTENT(IN) :: j1,j2,m1,m2			! Parametri Wigner
INTEGER(lo), INTENT(IN) :: j3min,j3max		! Valore minimo e massimo di j3
REAL(dbl), DIMENSION(j3min:j3max), INTENT(OUT) :: v_w3jm	! Vettore contenente i coeff. di W.

! Dichiarazione variabili interne
REAL(dbl) :: cd1,cd2,cu1,cu2,cu3		! Coefficenti ricorsioni
REAL(dbl) :: j3,w3jm_temp			! j3 e storage temporanee
INTEGER(lo) :: j3int				! j3 tradotto intero

! Subroutine vera e propria

! Inizializzo gli indici per la downward recursion
j3=REAL(j3max,dbl)
j3int=j3max


! WRITE(*,*) "j3min,j3max:",j3min,j3max


! In questo if separo i casi in cui ho un vettore di lunghezza uno da quelli che 
! necessitano dell'uso della ricorsione
big_if: IF (j3min==j3max) THEN

	v_w3jm(j3max)=wdown0(j1,j2,m1,m2) ! Unico termine da calcolare
	
ELSE 
	
	! Si inizializza la ricorsione
	v_w3jm(j3max)=wdown0(j1,j2,m1,m2)
	v_w3jm(j3max-1)=-(dr(j1,j2,j1+j2,m1,m2,-m1-m2)/( (j1+j2+1)*cr(j1,j2,j1+j2,m1,m2,-m1-m2) ))*v_w3jm(j3max)
	

	! Ciclo per il calcolo ricorsivo
	down_do: DO 		
			
			!Condizione di uscita
			IF (j3int-2<j3min) EXIT
				
			!Primo coeff della ricorsione
			cd1=dr(j1,j2,j3-1.0D0,m1,m2,-m1-m2)/(j3*cr(j1,j2,j3-1.0D0,m1,m2,-m1-m2))
			cd2=((j3-1.0D0)*cr(j1,j2,j3,m1,m2,-m1-m2))/(j3*cr(j1,j2,j3-1.0D0,m1,m2,-m1-m2)) 						
			!Ricorsione
			v_w3jm(j3int-2)=-cd1*v_w3jm(j3int-1)-cd2*v_w3jm(j3int)
			
			!Aggiorno gli indici
			j3=j3-1.0D0
			j3int=INT(j3,lo)
			
	END DO down_do
	
	! Inizializzo gli indici per la upward recursion
	j3int=j3min
	j3=REAL(j3min,dbl)
	
	! Calcolo del primo termine di wigner dal basso
	v_w3jm(j3int)=wup0(j1,j2,m1,m2)
	
	! Calcolo del secondo termine di wigner dal basso
	! Pongo anche la condizione sul coefficienti nel caso ci sia signolarita'
	cu3_if: IF (j3min==0) THEN
			cu3=0.0D0
	ELSE
			cu3=dr(j1,j2,j3,m1,m2,-m1-m2)/(j3*cr(j1,j2,j3+1.0D0,m1,m2,-m1-m2))
	END IF cu3_if
	
	w3jm_temp=-cu3*v_w3jm(j3int)
	
	! If legato alla monotonia della successione
	up_if: IF (ABS(w3jm_temp)>ABS(v_w3jm(j3min))) THEN
	
			! Aggiorno gli indici e metto nell'array il secondo valore
			! in questo modo sono pronto per iniziale la upward recursion
			! a tre termini
			j3int=j3int+1
			v_w3jm(j3int)=w3jm_temp
					
			up_do: DO
					!Aggiorno gli indici
					j3int=j3int+1
					j3=REAL(j3int,dbl)
					
					IF (j3int-1==j3max) EXIT
					
! 					IF ((INT(-m1)==-1).AND.(INT(j1)==1).AND.(INT(m2)==1).AND.(INT(j2)==2)) THEN
! 					WRITE(*,*) "j3-1,cr1,cr2",j3-1,cr(j1,j2,j3,m1,m2,-m1-m2),cr(j1,j2,j3,m1,m2,-m1-m2)
! 					END IF
										
					!Primo e secondo coeff della ricorsione
					cu1=dr(j1,j2,j3-1.0D0,m1,m2,-m1-m2)/((j3-1.0D0)*cr(j1,j2,j3,m1,m2,-m1-m2))
					cu2=(j3*cr(j1,j2,j3-1.0D0,m1,m2,-m1-m2))/((j3-1.0D0)*cr(j1,j2,j3,m1,m2,-m1-m2)) 
				
					!Assegnazione temporanea della ricorsione
					w3jm_temp=-cu1*v_w3jm(j3int-1)-cu2*v_w3jm(j3int-2)

					IF ((ABS(w3jm_temp)<ABS(v_w3jm(j3int-1)))  .OR. ((j3int-1)==j3max) ) EXIT ! Cond. di uscita
					
					v_w3jm(j3int)=w3jm_temp	!Assegno perche' e' ok
					
			END DO up_do
			
	END IF up_if
			
END IF big_if

END SUBROUTINE wigner3jm












! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! GAUNT SUBROUTINES: COEFF. DI GAUNT SECONDO CRUZAN
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!******************************************************************************
!1) SUBROUTINE gaunt_cz: calcolo una serie di coefficienti di Gaunt per una data
!						 quadrupletta di indici
!******************************************************************************
SUBROUTINE gaunt_cz(m,n,mu,nu,qmaxa,v_aq,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl), INTENT(IN) :: n,m,mu,nu					! Quadrupletta di indici consueta
INTEGER(lo), INTENT(IN) :: qmaxa					! Upper bound del vettore di coeff di gaunt e di v_bq
REAL(dbl), DIMENSION(0:qmaxa), INTENT(OUT) :: v_aq	! Vettore coeff di gaunt
INTEGER(lo), INTENT(OUT) :: error					! Errore

! Dichiarazione variabili interne
INTEGER(lo) :: stat_a
INTEGER(lo) :: pmin,pmax,pmin0						! Up and low bounds dei vettori di wigner
INTEGER(lo) :: q,p,r									! Indici coeff di gaunt ed errore lnf
REAL(dbl) :: pr											! Corrispondente reale di p
REAL(dbl) :: logw,fac									! Per il fattoriale che serve per gaunt
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_w3jm,v_w3jm0	! Vettori coeff di wigner

!--------------------------
! Subroutine vera e propria
!--------------------------
error=0

! Check sull'errore
error_if: IF ((ABS(m)>n) .OR. (ABS(mu)>nu)) THEN
			WRITE(*,*)
			WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_cz:"
			WRITE(*,*) "|m|>n oppure |mu|>nu, la subroutine si ferma!!!"
			WRITE(*,*)
			error=1
			RETURN
END IF error_if

! Calcolo i bounds dei vettori di wigner
pmin=j3minimum(n,nu,m,mu)
pmax=INT(n+nu,lo)
pmin0=j3minimum(n,nu,0.0D0,0.0D0)

! Alloco i vettori di wigner e li calcolo
ALLOCATE(v_w3jm(pmin:pmax),v_w3jm0(pmin0:pmax),STAT=stat_a)
CALL wigner3jm(n,nu,m,mu,pmin,pmax,v_w3jm)
CALL wigner3jm(n,nu,0.0D0,0.0D0,pmin0,pmax,v_w3jm0)

! Entro nel ciclo per il calcolo dei coefficienti di gaunt
gaunt_do: DO q=0,qmaxa

			! Calcolo dell'indice p, sia reale che intero
			p=INT(n+nu,lo)-2*q
			pr=REAL(p,dbl)
			
			!Calcolo del fattoriale
			logw = 0.5D0* (lnf(n+m,r)+lnf(nu+mu,r)+lnf(pr-m-mu,r) - &
				 		   lnf(n-m,r)-lnf(nu-mu,r)-lnf(pr+m+mu,r))
			fac= EXP(logw)
 			
			! Calcolo del coefficiente di gaunt
			v_aq(q)=((-1.0D0)**INT(m+mu,lo))*(2.0D0*pr+1.0D0)*fac*v_w3jm(p)*v_w3jm0(p)
			
END DO gaunt_do

! Disalloco i vettori di wigner a lavoro finito
DEALLOCATE(v_w3jm,v_w3jm0)
			
END SUBROUTINE gaunt_cz


!******************************************************************************
!2) SUBROUTINE aqbq_cz: calcolo una serie di coefficienti di Gaunt e dei coeff bq
!per una data quadrupletta di indici
!******************************************************************************
SUBROUTINE aqbq_cz(m,n,mu,nu,qmaxa,v_aq,qmaxb,v_bq,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl), INTENT(IN) :: m,n,mu,nu					! Quadrupletta di indici consueta
INTEGER(lo), INTENT(IN) :: qmaxa,qmaxb			! Upper bound del vettore di coeff di gaunt e di v_bq
INTEGER(lo), INTENT(OUT) :: error					! Flag di errore
REAL(dbl), DIMENSION(0:qmaxa), INTENT(OUT) :: v_aq	! Vettore coeff di gaunt
REAL(dbl), DIMENSION(1:), INTENT(OUT) :: v_bq		! Vettore bq, non metto l'upper bound se non ho problemi d'allocazione
! Dichiarazione variabili interne
INTEGER(lo) :: stat_a
INTEGER(lo) :: pmin,pmax,pmin0						! Up and low bounds dei vettori di wigner
INTEGER(lo) :: q,p,r									! Indici coeff di gaunt ed errore lnf
REAL(dbl) :: pr											! Corrispondente reale di p
REAL(dbl) :: logw,fac									! Per il fattoriale che serve per gaunt
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_w3jm,v_w3jm0	! Vettori coeff di wigner

!--------------------------
! Subroutine vera e propria
!--------------------------
error=0

! Check sull'errore
error_if: IF ((ABS(m)>n) .OR. (ABS(mu)>nu)) THEN
			WRITE(*,*)
			WRITE(*,*) "Si e' verificato un errore nella subroutine aqbq_cz:"
			WRITE(*,*) "|m|>n oppure |mu|>nu, la subroutine si ferma!!!"
			WRITE(*,*)
			error=1
			RETURN
ELSE IF (SIZE(v_bq)==0) THEN
			WRITE(*,*)
			WRITE(*,*) "Si e' verificato un errore nella subroutine aqbq_cz:"
		  	WRITE(*,*) "v_bq ha lunghezza zero, la subroutine si ferma"
		  	WRITE(*,*)
		  	error=1
		  	RETURN
END IF error_if

! Calcolo i bounds dei vettori di wigner
pmin=j3minimum(n,nu,m,mu)
pmax=INT(n+nu,lo)
pmin0=j3minimum(n,nu,0.0D0,0.0D0)

! Alloco i vettori di wigner e li calcolo
ALLOCATE(v_w3jm(pmin:pmax),v_w3jm0(pmin0:pmax),STAT=stat_a)
CALL wigner3jm(n,nu,m,mu,pmin,pmax,v_w3jm)
CALL wigner3jm(n,nu,0.0D0,0.0D0,pmin0,pmax,v_w3jm0)

! Entro nel ciclo per il calcolo dei coefficienti di gaunt
gaunt_do: DO q=0,qmaxa

			! Calcolo dell'indice p, sia reale che intero
			p=INT(n+nu,lo)-2*q
			pr=REAL(p,dbl)
			
			!Calcolo del fattoriale
			logw = 0.5D0* (lnf(n+m,r)+lnf(nu+mu,r)+lnf(pr-m-mu,r) - &
				 		   lnf(n-m,r)-lnf(nu-mu,r)-lnf(pr+m+mu,r))
			fac= EXP(logw)
 			
			! Calcolo del coefficiente di gaunt
			v_aq(q)=((-1.0D0)**INT(m+mu,lo))*(2.0D0*pr+1.0D0)*fac*v_w3jm(p)*v_w3jm0(p)
			
END DO gaunt_do


!Entro nel ciclo per il calcolo dei coefficienti bq
qmaxb_if: IF (qmaxb==0) THEN
			v_bq(1)=0.0D0
ELSE

	bq_do: DO q=1,qmaxb

				! Calcolo dell'indice p, sia reale che intero
				p=INT(n+nu,lo)-2*q
				pr=REAL(p,dbl)
				
				!Calcolo del fattoriale
				logw = 0.5D0* (lnf(n+m,r)+lnf(nu+mu,r)+lnf(pr-m-mu+1.0D0,r) - &
					 		   lnf(n-m,r)-lnf(nu-mu,r)-lnf(pr+m+mu+1.0D0,r))
				fac= EXP(logw)
 			
				! Calcolo del coefficiente di gaunt
				v_bq(q)=((-1.0D0)**INT(m+mu,lo))*(2.0D0*pr+3.0D0)*fac*v_w3jm(p+1)*v_w3jm0(p)
			
	END DO bq_do
	
END IF qmaxb_if

! Disalloco i vettori di wigner a lavoro finito
DEALLOCATE(v_w3jm,v_w3jm0)
			
END SUBROUTINE aqbq_cz

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! GAUNT SUBROUTINES: COEFF. DI GAUNT SECONDO XU
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!******************************************************************************
!1) FUNCTION F_ALPHA:funzione per calcolare il parametro alpha(n,nu,p)
!******************************************************************************

FUNCTION f_alpha(n,nu,p)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_alpha

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: n,nu,p		! Indici consueti e indice costruito p 

! Funzione vera e propria
f_alpha= ((p**2 -(n-nu)**2)*(p**2-(n+nu+1)**2))/(4*(p**2)-1)

END FUNCTION f_alpha


!******************************************************************************
!2) FUNCTION f_Ap: funzione per il calcolo del coefficiente Ap(m,n,mu,nu,p)
!******************************************************************************

FUNCTION f_Ap(m,n,mu,nu,p)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_Ap

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu,p		! Indici consueti e indice costruito p

! Funzione vera e propria
f_Ap=p*(p-1)*(m-mu)-(m+mu)*(n-nu)*(n+nu+1)

END FUNCTION f_Ap


!***********************************************************************************
!3) FUNCTION f_a0: starting value per la backward (q crescente) recurrence di gaunt
!***********************************************************************************

FUNCTION f_a0(m,n,mu,nu)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_a0

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu		! Indici consueti

! Dichiarazione variabili interne
INTEGER(lo) :: ierr
REAL(dbl) :: logw,logp

! Funzione vera e propria
logw = lnf(n+nu-m-mu,ierr)-lnf(n-m,ierr)-lnf(nu-mu,ierr)
logp = lpoch(n+1.0D0,n,ierr)+lpoch(nu+1.0D0,nu,ierr)-lpoch(n+nu+1.0D0,n+nu,ierr)

f_a0=EXP(logw+logp)

END FUNCTION f_a0


!*********************************************************************************
!4) FUNCTION f_a1norm: coeff a(m,n,mu,nu,1) normalizzato per la backward recursion
!*********************************************************************************

FUNCTION f_a1norm(m,n,mu,nu)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_a1norm

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu			! Indici consueti

! Dichiarazione variabili interne
REAL(dbl) :: n4

! Funzione vera e propria
n4=n+nu-m-mu

f_a1norm=((2.0D0*n+2.0D0*nu-3.0D0)/2.0D0) * (1.0D0-((2.0D0*n+2.0D0*nu-1.0D0)/(n4*(n4-1.0D0))) * &
		 ((m-n)*(m-n+1.0D0)/(2.0D0*n-1.0D0)+(mu-nu)*(mu-nu+1.0D0)/(2.0D0*nu-1.0D0)) )

END FUNCTION f_a1norm


!*********************************************************************************
!5) FUNCTION f_a2norm: coeff a(m,n,mu,nu,2) normalizzato per la backward recursion
!*********************************************************************************

FUNCTION f_a2norm(m,n,mu,nu)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_a2norm

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu	! Indici Consueti

! Dichiarazione variabili interne
REAL(dbl) :: n4,n2nu,mn,munu

! Funzione vera e propria
n4=n+nu-m-mu
n2nu=2.0D0*n+2.0D0*nu
mn=m-n
munu=mu-nu

f_a2norm=((n2nu-1.0D0)*(n2nu-7.0D0)/4.0D0) * &
		 ( ((n2nu-3.0D0)/(n4*(n4-1.0D0))) * &
		 ( ((n2nu-5.0D0)/(2.0D0*(n4-2.0D0)*(n4-3.0D0))) * &
		 ( mn*(mn+1.0D0)*(mn+2.0D0)*(mn+3.0D0)/((2.0D0*n-1.0D0)*(2.0D0*n-3.0D0)) + &
		   2.0D0*mn*(mn+1.0D0)*munu*(munu+1.0D0)/((2.0D0*n-1.0D0)*(2.0D0*nu-1.0D0)) + &
		   munu*(munu+1.0D0)*(munu+2.0D0)*(munu+3.0D0)/((2.0D0*nu-1.0D0)*(2.0D0*nu-3.0D0)) &
		   ) - mn*(mn+1.0D0)/(2.0D0*n-1.0D0) - munu*(munu+1.0D0)/(2.0D0*nu-1.0D0) ) +0.5D0)
		   
		   
END FUNCTION f_a2norm


!*********************************************************************************
!6) FUNCTION f_a2normr: coeff a(m,n,mu,nu,2) normalizzato per la backward recursion
!						calcolato questa volta per ricorsione
!*********************************************************************************

FUNCTION f_a2normr(m,n,mu,nu,a1norm)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_a2normr

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu,a1norm			! Indici consueti e aq(q=1) normalizzato

! Dichiarazione variabili interne
REAL(dbl) :: c0,c1,c2
REAL(dbl) ::p,p1,p2
REAL(dbl) :: alphap1,alphap2,alphap3,alphap4,alphap5
REAL(dbl) :: Ap2,Ap3,Ap4,Ap5,Ap6

! Funzione vera e propria
p=n+nu-4
p1=p-m-mu
p2=p+m+mu

Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
Ap6=f_Ap(m,n,mu,nu,p+6.0D0)

alphap1=f_alpha(n,nu,p+1.0D0)
alphap2=f_alpha(n,nu,p+2.0D0)

Ap_if: IF ((INT(Ap4,lo)==0) .AND. (INT(Ap6,lo)==0)) THEN

		c0=(p+2.0D0)*(p1+1.0D0)*alphap1
		c1=(p+1.0D0)*(p2+2.0D0)*alphap2
		
		f_a2normr=(c1/c0)*a1norm
		
	   ELSE IF ((INT(Ap4,lo)==0) .AND. (INT(Ap6,lo)/=0)) THEN
	   
	   	Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
	   	Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
		Ap5=f_Ap(m,n,mu,nu,p+5.0D0)
		alphap5=f_alpha(n,nu,p+5.0D0)
		
		c0=(p+2.0D0)*(p+3.0D0)*(p+5.0D0)*(p1+1.0D0)*(p1+2.0D0)*(p1+4.0D0)*Ap6*alphap1
		c1=(p+5.0D0)*(p1+4.0D0)*Ap6*(Ap2*Ap3+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*alphap2)
		c2=(p+2.0D0)*(p2+3.0D0)*Ap2*(Ap5*Ap6+(p+4.0D0)*(p+6.0D0)*(p1+5.0D0)*(p2+5.0D0)*alphap5)
		
		f_a2normr=(c1/c0)*a1norm+(c2/c0)
		
	  ELSE
	  
	  	Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
	   	Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
	   	alphap3=f_alpha(n,nu,p+3.0D0)
	   	alphap4=f_alpha(n,nu,p+4.0D0)
	   	
	   	c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
	   	c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2 &
	   	 + (p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
	   	c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
	   	
	   	f_a2normr=(c1/c0)*a1norm+(c2/c0)
		   
END IF Ap_if   

END FUNCTION f_a2normr


!*********************************************************************************
!7) FUNCTION f_aqmax: coeff a(m,n,mu,nu,qmax) per la forward recursion
!*********************************************************************************

FUNCTION f_aqmax(m,n,mu,nu,qmax)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_aqmax

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu		! Indici consueti
INTEGER(lo), INTENT(IN) :: qmax		! qmax
! Dichiarazione variabili interne
REAL(dbl) :: pmin,logw,logp,qmaxr,Apmin
INTEGER(lo) :: pmin_i,ierr

! Funzione vera e propria
pmin=n+nu-2*REAL(qmax,dbl)
pmin_i=INT(pmin,lo)
qmaxr=REAL(qmax,dbl)


pmin_if: IF (pmin_i==INT(n-nu,lo)) THEN

	logw=lnf(n+m,ierr)+lnf(2.0D0*pmin+1.0D0,ierr) &
		-lnf(nu-mu,ierr)-lnf(n-nu,ierr)-lnf(pmin+m+mu,ierr)
	
	logp=lpoch(nu+1.0D0,nu,ierr)-lpoch(n+1.0D0,n+1.0D0,ierr)
	
! 	f_aqmax=((-1.0D0)**INT(mu,lo))*dpoch(nu+1.0D0,nu)*EXP(logw)/dpoch(n+1,n+1)	

	f_aqmax=((-1.0D0)**INT(mu,lo))*EXP(logw+logp)	

ELSE IF (pmin_i==INT(nu-n,lo)) THEN

	logw=lnf(nu+mu,ierr)+lnf(2.0D0*pmin+1.0D0,ierr) &
		-lnf(n-m,ierr)-lnf(nu-n,ierr)-lnf(pmin+m+mu,ierr)
		
	logp=lpoch(n+1.0D0,n,ierr)-lpoch(nu+1.0D0,nu+1.0D0,ierr)	
	
! 	f_aqmax=((-1.0D0)**INT(m,lo))*dpoch(n+1.0D0,n)*EXP(logw)/dpoch(nu+1,nu+1)

	f_aqmax=((-1.0D0)**INT(m,lo))*EXP(logw+logp)

ELSE IF (pmin_i==INT(m+mu,lo)) THEN

	logw=lpoch(qmaxr+1.0D0,qmaxr,ierr)+lnf(n+nu-qmaxr,ierr)+lnf(n+m,ierr)+lnf(nu+mu,ierr) &
		-lnf(n-qmaxr,ierr)-lnf(nu-qmaxr,ierr)-lnf(n-m,ierr)-lnf(nu-mu,ierr)-lnf(n+nu+pmin+1.0D0,ierr)

	f_aqmax=((-1.0D0)**INT(n+m-qmaxr,lo))*(2.0D0*pmin+1.0D0)*EXP(logw)

ELSE IF (pmin_I==INT(-m-mu,lo)) THEN

	logw=lpoch(qmaxr+1.0D0,qmaxr,ierr)+lnf(n+nu-qmaxr,ierr)+lnf(pmin-m-mu,ierr) &
		-lnf(n-qmaxr,ierr)-lnf(nu-qmaxr,ierr)-lnf(n+nu+pmin+1.0D0,ierr)
	
	f_aqmax=((-1.0D0)**INT(nu+mu-qmaxr,lo))*(2.0D0*pmin+1.0D0)*EXP(logw)

ELSE IF (pmin_i==INT(m+mu+1,lo)) THEN

	Apmin=f_Ap(m,n,mu,nu,pmin)

	logw=lpoch(qmaxr+1.0D0,qmaxr,ierr)+lnf(n+nu-qmaxr,ierr)+lnf(n+m,ierr)+lnf(nu+mu,ierr) &
		-lnf(n+nu+pmin+1.0D0,ierr)-lnf(n-qmaxr,ierr)-lnf(nu-qmax,ierr)-lnf(n-m,ierr)-lnf(nu-mu,ierr)
	
	f_aqmax=((-1.0D0)**INT(n+m-qmaxr,lo))*Apmin*(2.0D0*pmin+1.0D0)*EXP(logw) &
		   /(pmin-1.0D0)

ELSE IF (pmin_i==INT(-m-mu+1,lo)) THEN

	Apmin=f_Ap(m,n,mu,nu,pmin)

	logw=lpoch(qmaxr+1.0D0,qmaxr,ierr)+lnf(n+nu-qmaxr,ierr)+lnf(pmin-m-mu,ierr) &
		-lnf(n+nu+pmin+1.0D0,ierr)-lnf(n-qmaxr,ierr)-lnf(nu-qmaxr,ierr)
	
	f_aqmax=((-1.0D0)**INT(nu+mu-qmaxr,lo))*Apmin*(2.0D0*pmin+1.0D0)*EXP(logw) &
		   /(pmin-1.0D0)

END IF pmin_if

END FUNCTION f_aqmax


!*********************************************************************************
!8) FUNCTION f_aqmax_1: coeff a(m,n,mu,nu,qmax-1) per la forward recursion
!*********************************************************************************

FUNCTION f_aqmax_1(m,n,mu,nu,qmax)

IMPLICIT NONE

! Dichiarazione Funzione
REAL(dbl) :: f_aqmax_1

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu			! Indici Consueti
INTEGER(lo), INTENT(IN) :: qmax			! qmax
! Dichiarazione variabili interne
REAL(dbl) :: pmin,logw,qmaxr,Apmin2,f1,f2
INTEGER(lo) :: pmin_i,ierr

! Funzione vera e propria
pmin=n+nu-2*REAL(qmax,dbl)
pmin_i=INT(pmin,lo)
qmaxr=REAL(qmax,dbl)
Apmin2=f_Ap(m,n,mu,nu,pmin+2.0D0)

pmin_if: IF (pmin_i==INT(m+mu+1.0D0,lo)) THEN

	logw=lpoch(qmaxr+1.0D0,qmaxr,ierr)+lnf(n+nu-qmaxr,ierr)+lnf(n+m,ierr)+lnf(nu+mu,ierr) &
		-lnf(n+nu+pmin,ierr)-lnf(n-qmaxr+1.0D0,ierr)-lnf(nu-qmaxr+1.0D0,ierr)-lnf(n-m,ierr)-lnf(nu-mu,ierr)

	f1=((-1.0D0)**INT(m+n-qmaxr,lo))*Apmin2*(2.0D0*pmin+3.0D0)*(2.0D0*pmin+5.0D0) &
	  /(pmin*(n+nu+pmin+3.0D0))

	f2=(n-qmaxr)*(nu-qmaxr)*(2.0D0*qmaxr+1.0D0) &
	  /((pmin+m+mu+1.0D0)*(pmin+m+mu+2.0D0)*(2*qmaxr-1.0D0))
	  
	f_aqmax_1=f1*f2*EXP(logw)
	  	
ELSE IF (pmin_i==INT(-m-mu+1.0D0,lo)) THEN

	logw=lpoch(qmaxr+1.0D0,qmaxr+1.0D0,ierr)+lnf(n+nu-qmaxr,ierr)+lnf(pmin-m-mu,ierr) &
		-lnf(n+nu+pmin,ierr)-lnf(n-qmaxr+1.0D0,ierr)-lnf(nu-qmaxr+1.0D0,ierr)
		
	f1=((-1.0D0)**INT(nu+mu-qmaxr,lo))*Apmin2*(2.0D0*pmin+3.0D0)*(2.0D0*pmin+5.0D0) &
	  /(6.0D0*pmin*(n+nu+pmin+3.0D0))

	f2=(n-qmaxr)*(nu-qmaxr)/(2.0D0*qmaxr-1.0D0)

	f_aqmax_1=f1*f2*EXP(logw)

END IF pmin_if

END FUNCTION f_aqmax_1




!*********************************************************************************
!9) SUBROUTINE gaunt_xu: calcolo i coefficienti di gaunt tramite il formalismo di xu
!*********************************************************************************

SUBROUTINE gaunt_xu2(m,n,mu,nu,qmax,v_aq,error)

IMPLICIT NONE

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu						! Indici coefficienti di Gaunt Reali
INTEGER(lo), INTENT(IN) :: qmax						! Qmax, bound dell'indice
REAL(dbl), DIMENSION(0:qmax), INTENT(OUT) :: v_aq		! Vettore coefficienti di Gaunt
INTEGER(lo), INTENT(OUT) :: error						! Variabile di errore

! Dichiarazione variabili interne
INTEGER(lo), PARAMETER :: prec=14								! Precisione backward forward recursion
REAL(dbl) :: c0,c1,c2											! Coefficienti ricorsioni
REAL(dbl) ::p,p1,p2,pmin										! Coeff p per il calcolo coeff ricorsioni
REAL(dbl) :: alphap1,alphap2,alphap3,alphap4					! Coeff alpha per il calcolo coeff ricorsioni
REAL(dbl) :: Apmin,Ap2,Ap3,Ap4									! Coeff Ap per il calcolo coeff ricorsioni
REAL(dbl) :: aq_fwd												! Variabile swap per la forward recursion
REAL(dbl) :: res												! Differenza relativa tra i due
INTEGER(lo) :: q,qi,switch									! Indice e var di switch
LOGICAL :: test

!---------------------------------------------------------
! Subroutine vera e propria
!---------------------------------------------------------
error=0

! Check sull'errore
error_if: IF ((ABS(m)>n) .OR. (ABS(mu)>nu)) THEN
			WRITE(*,*)
			WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu:"
			WRITE(*,*) "|m|>n oppure |mu|>nu, la subroutine si ferma!!!"
			WRITE(*,*)
			error=1
			RETURN
END IF error_if

!::::::::::::::::::::::::::::::::::
!Struttura case, per i diversi qmax
!::::::::::::::::::::::::::::::::::
qmax_case: SELECT CASE (qmax)

CASE (0) !qmax=0

	v_aq(0)=f_a0(m,n,mu,nu)
	!WRITE(*,*) f_a0(m,n,mu,nu),v_aq(0)
	
CASE (1) !qmax=1

	v_aq(0)=f_a0(m,n,mu,nu)
	v_aq(1)=v_aq(0)*f_a1norm(m,n,mu,nu)
	
CASE (2) !qmax=2

	v_aq(0)=f_a0(m,n,mu,nu)
	v_aq(1)=v_aq(0)*f_a1norm(m,n,mu,nu)
	v_aq(2)=v_aq(0)*f_a2normr(m,n,mu,nu,v_aq(1)/v_aq(0))
	
CASE (3:) !qmax>2

	!:::::::::::::::::::::::::::::::::::::::::::::::
	!Struttura if, per i diversi casi di ricorsione
	!:::::::::::::::::::::::::::::::::::::::::::::::
	
	big_if: IF ((INT(m,lo)==0) .AND. (INT(mu,lo)==0)) THEN
	
		!'''''''''''''''''''''''''''''''''''
		!Caso mu=m=0 (1)
		!'''''''''''''''''''''''''''''''''''
		
		!BACKWARD
		
		!Primo valore per la backward recursion
		v_aq(0)=f_a0(m,n,mu,nu) 
	
		!Backward recursion
		uno_b_do: DO q=1,qmax
					
					p=n+nu-2.0D0*REAL(q,dbl)
					c0=f_alpha(n,nu,p+1.0D0)
					c1=f_alpha(n,nu,p+2.0D0)
					
					v_aq(q)=(c1/c0)*v_aq(q-1)
					
		END DO uno_b_do
		
	ELSE IF ((INT(mu,lo)==INT(m,lo)) .AND. (INT(nu,lo)==INT(n,lo))) THEN
	
		!'''''''''''''''''''''''''''''''''''
		!Caso mu=m e nu=n (2)
		!'''''''''''''''''''''''''''''''''''
		
		!BACKWARD
		
		!Primo valore per la backward recursion
		v_aq(0)=f_a0(m,n,mu,nu) 
	
		!Backward recursion
		due_b_do: DO q=1,qmax
					
					!Calcolo pre-coefficienti
					p=n+nu-2.0D0*REAL(q,dbl)
					p1=p-m-mu
					p2=p+m+mu
					
					!Calcolo coefficienti ricorsione
					c0=(p+2.0D0)*(p1+1.0D0)*f_alpha(n,nu,p+1.0D0)
					c1=(p+1.0D0)*(p2+2.0D0)*f_alpha(n,nu,p+2.0D0)
					
					!Ricorsione
					v_aq(q)=(c1/c0)*v_aq(q-1)
					
		END DO due_b_do
			
	ELSE IF (INT(mu,lo)==INT(-m,lo)) THEN		
	
		!'''''''''''''''''''''''''''''''''''
		!Caso mu=-m (3)
		!'''''''''''''''''''''''''''''''''''
		
		!-----------------------------------
		!BACKWARD
		!-----------------------------------
		
		!Primo valore per la backward recursion
		v_aq(0)=f_a0(m,n,mu,nu)
		v_aq(1)=f_a1norm(m,n,mu,nu)*v_aq(0) 
	
		!Backward recursion
		tre_b_do: DO q=2,qmax
					
					!Calcolo pre-coefficienti
					p=n+nu-2.0D0*REAL(q,dbl)
				
					!Calcolo coefficienti ricorsione
					c0=f_alpha(n,nu,p+1.0D0)
					c1=4.0D0*(m**2)+f_alpha(n,nu,p+2.0D0)+f_alpha(n,nu,p+3.0D0)
					c2=-f_alpha(n,nu,p+4.0D0)
					
					!Ricorsione
					v_aq(q)=(c1/c0)*v_aq(q-1)+(c2/c0)*v_aq(q-2)
					
		END DO tre_b_do
		
		!--------------------------------
		!FORWARD
		!--------------------------------
		
		!Primo valore per la forward recursion,errore relativo e suo swap
		aq_fwd=f_aqmax(m,n,mu,nu,qmax)
		res=ABS(aq_fwd-v_aq(qmax))/ABS(aq_fwd)
		
		!Se non ho precisione, sostituisco i valori
		tre_f_if: IF (res>(10.0D0**(-prec))) THEN 
		
			v_aq(qmax)=aq_fwd
			qi=1
			
			!Entro nel ciclo della sostituzione valori
			tre_f_do: DO q=qmax-1,0,-1
				
				tre_q_case:SELECT CASE (qmax-q)
				
				CASE(1) tre_q_case	!q=qmax-1
					
					!Calcolo v_aq(qmax-1)
					p=n+nu-2.0D0*REAL(q+2,dbl)
					c1=4.0D0*(m**2)+f_alpha(n,nu,p+2.0D0)+f_alpha(n,nu,p+3.0D0)
					c2=-f_alpha(n,nu,p+4.0D0)
					aq_fwd=-(c1/c2)*v_aq(qmax)
					
					!If a secondo che v_aq(qmax-1) sia zero o no
					zero1_if: IF ((aq_fwd/v_aq(qmax))<(10.0D0**(-8))) THEN
						v_aq(q)=aq_fwd	!Zero
						switch=1
						CYCLE tre_f_do
					ELSE
						res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)	!Diverso da zero
					END IF zero1_if
					
				CASE DEFAULT tre_q_case !Per tutti gli altri q
				
					!Calcolo v_aq(qmax-1)
					p=n+nu-2.0D0*REAL(q+2.0D0,dbl)
					c0=f_alpha(n,nu,p+1.0D0)
					c1=4.0D0*(m**2)+f_alpha(n,nu,p+2.0D0)+f_alpha(n,nu,p+3.0D0)
					c2=-f_alpha(n,nu,p+4.0D0)
					aq_fwd=-(c1/c2)*v_aq(q+1)+(c0/c2)*v_aq(q+2)
				
					!Case a seconda che il valore appena precedente sia zero o meno
					tre_switch_case: SELECT CASE (switch)
					
					CASE(1) tre_switch_case !Il valore precedente e' zero
					
						res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)
						
					CASE(0) tre_switch_case !Il valore precedente e' diverso da zero
					
						!If a secondo che v_aq(q) sia zero o no
						zero2_if: IF ((aq_fwd/v_aq(q+1))<(10.0D0**(-8))) THEN
							v_aq(q)=aq_fwd	!Zero
							switch=1
							CYCLE tre_f_do
						ELSE
							res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)	!Diverso da zero
						END IF zero2_if
						
					END SELECT tre_switch_case
				
				END SELECT tre_q_case

			!Adesso se la precisione e' raggiunta esco dal ciclo, se no sostituisco e rimango
			IF ((res<(10.0D0**(-prec))) .OR. (q==0)) EXIT
			
			!Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
			v_aq(q)=aq_fwd
			qi=q
			
			END DO tre_f_do

		! Check sul ciclo di sostituzione
		error_if1: IF (q==0) THEN
					WRITE(*,*)
					WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu:"
					WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
					WRITE(*,*) "e forward recursion non e' stata raggiunta"
					WRITE(*,*)
					error=1
					RETURN
		END IF error_if1
		
		
		END IF tre_f_if
			
	ELSE
	
		!'''''''''''''''''''''''''''''''''''
		!Caso generale (4)
		!'''''''''''''''''''''''''''''''''''
		
		!-----------------------------------
		!BACKWARD
		!-----------------------------------
		
		!Calcolo Ap4 per q=0, se e' uguale a zero chiamo cruzan ed esco dal ciclo
		Ap4=f_Ap(m,n,mu,nu,n+nu-REAL(qmax,dbl)+4.0D0)
	
		cz_if1:IF (Ap4==0) THEN
! 			WRITE(*,*) "cz Ap4"
			CALL gaunt_cz(m,n,mu,nu,qmax,v_aq,error)
			RETURN
		END IF cz_if1
	
		!Calcolo direttamente i primi tre valori della ricorsione
		v_aq(0)=f_a0(m,n,mu,nu)
		v_aq(1)=v_aq(0)*f_a1norm(m,n,mu,nu)
	
		gen_b_do: DO q=2,qmax
		
					!Calcolo pre-coeff. : questi li calcoli qui per poter uscire
					p=n+nu-2.0D0*REAL(q,dbl)
					Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
					
					IF (Ap2==0) EXIT	!Esco dal loop perche non posso fare la fwd recursion
					
					!Calcolo pre-coefficienti
					p1=p-m-mu
					p2=p+m+mu
					alphap1=f_alpha(n,nu,p+1.0D0)
					alphap2=f_alpha(n,nu,p+2.0D0)
					alphap3=f_alpha(n,nu,p+3.0D0)
					alphap4=f_alpha(n,nu,p+4.0D0)
					Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
					Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
					
					!Calcolo coefficienti ricorsione
					c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
					c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
					  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
					c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
				
					!Ricorsione
					v_aq(q)=(c1/c0)*v_aq(q-1)+(c2/c0)*v_aq(q-2)
			
		END DO gen_b_do
		
		cz_if2:IF (Ap2==0) THEN
! 			WRITE(*,*) "cz Ap2"
			CALL gaunt_cz(m,n,mu,nu,qmax,v_aq,error)
			RETURN
		END IF cz_if2
	
	
		!-----------------------------------
		!FORWARD
		!-----------------------------------
		
		!Calcolo pmin,Apmin e la mia variabile logica
		pmin=n+nu-2.0D0*REAL(qmax,dbl)
		Apmin=f_Ap(m,n,mu,nu,pmin)
		test=((INT(Apmin,lo)==0) .AND. &
		    &((INT(pmin,lo)==INT(m+mu+1.0D0,lo)).OR.(INT(pmin,lo)==INT(-m-mu+1.0D0,lo))))
		
		!........................................................
		!Se la mia variabile logica e' vera, Faccio il mio conto
		!........................................................
		Apmin_if: IF (test) THEN

			!Il valore per qmax allora e' zero		
			v_aq(qmax)=0.0D0
			
			!Calcolo il secondo valore, e se la precisione e' raggiunta esco
			aq_fwd=f_aqmax_1(m,n,mu,nu,qmax)
			res=ABS(aq_fwd-v_aq(qmax-1))/ABS(aq_fwd)
			IF (res<(10.0D0**(-prec))) THEN
				RETURN
			END IF
		
			!Assegno il secondo valore e faccio il ciclo
			v_aq(qmax-1)=aq_fwd
			qi=1
			
			Apmin_do: DO q=(qmax-2),0,-1
			
				!Calcolo pre-coefficienti
				p=n+nu-2.0D0*REAL(q+2,dbl)
				p1=p-m-mu
				p2=p+m+mu
				alphap1=f_alpha(n,nu,p+1.0D0)
				alphap2=f_alpha(n,nu,p+2.0D0)
				alphap3=f_alpha(n,nu,p+3.0D0)
				alphap4=f_alpha(n,nu,p+4.0D0)
				Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
				Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
				Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
				
				!Calcolo coefficienti ricorsione
				c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
				c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
				  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
				c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
						
				!Ricorsione e residuo
				aq_fwd=-(c1/c2)*v_aq(q+1)+(c0/c2)*v_aq(q+2)
				res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)
				
				IF (res<(10.0D0**(-prec))) EXIT
				
				v_aq(q)=aq_fwd
				qi=q
				
			END DO Apmin_do		
			
			! Check sul ciclo di sostituzione
			Apmin_error_if1: IF (qi==0) THEN
						WRITE(*,*)
						WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu, caso generale, Apmin=0:"
						WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
						WRITE(*,*) "e forward recursion non e' stata raggiunta"
						WRITE(*,*)
						error=1
						RETURN
			END IF Apmin_error_if1
			
			!Esco dalla subroutine gaunt_xu
			RETURN
			
		END IF Apmin_if
		
		!..............................................
		!Qui lavoro se la mia variabile logica e' falsa
		!Tutto e' identico al caso (3)
		!..............................................
		
		!Primo valore per la forward recursion,errore relativo e suo swap
		aq_fwd=f_aqmax(m,n,mu,nu,qmax)
		res=ABS(aq_fwd-v_aq(qmax))/ABS(aq_fwd)
		qi=1
		
		!Se non ho precisione, sostituisco i valori
		gen_f_if: IF (res>(10.0D0**(-prec))) THEN 
		
			v_aq(qmax)=aq_fwd
			
			!Entro nel ciclo della sostituzione valori
			gen_f_do: DO q=qmax-1,0,-1
				
				gen_q_case:SELECT CASE (qmax-q)
				
				CASE(1) gen_q_case	!q=qmax-1
					
					!Calcolo aq
					p=n+nu-2.0D0*REAL(q+2,dbl)
					p1=p-m-mu
					p2=p+m+mu
					alphap1=f_alpha(n,nu,p+1.0D0)
					alphap2=f_alpha(n,nu,p+2.0D0)
					alphap3=f_alpha(n,nu,p+3.0D0)
					alphap4=f_alpha(n,nu,p+4.0D0)
					Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
					Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
					Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
					c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
					  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
					c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
					aq_fwd=-(c1/c2)*v_aq(qmax) !E' qui che lo calcolo
						
					!If a secondo che v_aq(qmax-1) sia zero o no
					gen_zero1_if: IF ((aq_fwd/v_aq(qmax))<(10.0D0**(-8))) THEN
							v_aq(q)=aq_fwd	!Zero
							switch=1
							CYCLE gen_f_do
					ELSE
							res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)	!Diverso da zero
					END IF gen_zero1_if
					
				CASE DEFAULT gen_q_case !Per tutti gli altri q
				
					!Calcolo aq
					p=n+nu-2.0D0*REAL(q+2,dbl)
					p1=p-m-mu
					p2=p+m+mu
					alphap1=f_alpha(n,nu,p+1.0D0)
					alphap2=f_alpha(n,nu,p+2.0D0)
					alphap3=f_alpha(n,nu,p+3.0D0)
					alphap4=f_alpha(n,nu,p+4.0D0)
					Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
					Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
					Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
					c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
					c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
					  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
					c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
					aq_fwd=-(c1/c2)*v_aq(q+1)+(c0/c2)*v_aq(q+2)
				
					!Case a seconda che il valore appena precedente sia zero o meno
					gen_switch_case: SELECT CASE (switch)
					
					CASE(1) gen_switch_case !Il valore precedente e' zero
					
						res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)
						
					CASE(0) gen_switch_case !Il valore precedente e' diverso da zero
					
						!If a secondo che v_aq(q) sia zero o no
						gen_zero2_if: IF ((aq_fwd/v_aq(q+1))<(10.0D0**(-9))) THEN
							v_aq(q)=aq_fwd	!Zero
							switch=1
							CYCLE gen_f_do
						ELSE
							res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)	!Diverso da zero
						END IF gen_zero2_if
						
					END SELECT gen_switch_case
				
				END SELECT gen_q_case

			!Adesso se la precisione e' raggiunta esco dal ciclo, se no sostituisco e rimango
			IF ((res<(10.0D0**(-prec))) .OR. (q==0)) EXIT
			
			!Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
			v_aq(q)=aq_fwd
			qi=q
			
			END DO gen_f_do

			! Check sul ciclo di sostituzione
			gen_error_if1: IF (qi==0) THEN
						WRITE(*,*)
						WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu,caso generale:"
						WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
						WRITE(*,*) "e forward recursion non e' stata raggiunta"
						WRITE(*,*)
						error=1
						RETURN
			END IF gen_error_if1
		
		
		END IF gen_f_if

	END IF big_if
		
END SELECT qmax_case





END SUBROUTINE gaunt_xu2




!*********************************************************************************
!9) SUBROUTINE gaunt_xu2: calcolo i coefficienti di gaunt tramite il formalismo di xu
!*********************************************************************************

SUBROUTINE gaunt_xu(m,n,mu,nu,qmax,v_aq,error)

IMPLICIT NONE

! Dichiarazione argomenti
REAL(dbl), INTENT(IN) :: m,n,mu,nu						! Indici coefficienti di Gaunt Reali
INTEGER(lo), INTENT(IN) :: qmax						! Qmax, bound dell'indice
REAL(dbl), DIMENSION(0:qmax), INTENT(OUT) :: v_aq		! Vettore coefficienti di Gaunt
INTEGER(lo), INTENT(OUT) :: error						! Variabile di errore

! Dichiarazione variabili interno
REAL(dbl) ,DIMENSION(0:qmax) :: v_aq_cz							! Vettore ausiliario coeff cruzan
INTEGER(lo), DIMENSION(0:qmax) :: v_zero						! Vettore memorizzazione degli zeri
INTEGER(lo), PARAMETER :: prec=10,tol=8						! Precisione backward forward recursion
REAL(dbl) :: c0,c1,c2,c3										! Coefficienti ricorsioni
REAL(dbl) ::p,p1,p2,pmin										! Coeff p per il calcolo coeff ricorsioni
REAL(dbl) :: alphap1,alphap2,alphap3,alphap4,alphap5,alphap6	! Coeff alpha per il calcolo coeff ricorsioni
REAL(dbl) :: Apmin,Ap2,Ap3,Ap4,Ap5,Ap6							! Coeff Ap per il calcolo coeff ricorsioni
REAL(dbl) :: aq_fwd												! Variabile swap per la forward recursion
REAL(dbl) :: res												! Differenza relativa tra i due
INTEGER(lo) :: q,qi,switch,q4									! Indice e var di switch
LOGICAL :: test

!---------------------------------------------------------
! Subroutine vera e propria
!---------------------------------------------------------
error=0
v_zero=1

! Check sull'errore
error_if: IF ((ABS(m)>n) .OR. (ABS(mu)>nu)) THEN
			WRITE(*,*)
			WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu:"
			WRITE(*,*) "|m|>n oppure |mu|>nu, la subroutine si ferma!!!"
			WRITE(*,*)
			error=1
			RETURN
END IF error_if

!::::::::::::::::::::::::::::::::::
!Struttura case, per i diversi qmax
!::::::::::::::::::::::::::::::::::
qmax_case: SELECT CASE (qmax)

CASE (0) !qmax=0

	v_aq(0)=f_a0(m,n,mu,nu)
	!WRITE(*,*) f_a0(m,n,mu,nu),v_aq(0)
	
CASE (1) !qmax=1

	v_aq(0)=f_a0(m,n,mu,nu)
	v_aq(1)=v_aq(0)*f_a1norm(m,n,mu,nu)
	
	!Controllo gli zeri
	zg_if_1_0: IF (ABS(v_aq(1)/v_aq(0))<10.0D0**(-tol)) THEN
				v_aq(1)=0.0D0
				v_zero(1)=0
	END IF zg_if_1_0
	
CASE (2) !qmax=2

	v_aq(0)=f_a0(m,n,mu,nu)
	v_aq(1)=v_aq(0)*f_a1norm(m,n,mu,nu)
	
	!Controllo gli zeri
	zg_if_1_1: IF (ABS(v_aq(1)/v_aq(0))<10.0D0**(-tol)) THEN
				v_aq(1)=0.0D0
				v_zero(1)=0
	END IF zg_if_1_1
	
	v_aq(2)=v_aq(0)*f_a2normr(m,n,mu,nu,v_aq(1)/v_aq(0))
	
	!Controllo gli zeri
	zg_if_1_2: IF (ABS(v_aq(2)/v_aq(0))<10.0D0**(-tol)) THEN
				v_aq(2)=0.0D0
				v_zero(2)=0
	END IF zg_if_1_2
	
CASE (3:) !qmax>2

	!:::::::::::::::::::::::::::::::::::::::::::::::
	!Struttura if, per i diversi casi di ricorsione
	!:::::::::::::::::::::::::::::::::::::::::::::::
	
	big_if: IF ((INT(m,lo)==0) .AND. (INT(mu,lo)==0)) THEN
	
		!'''''''''''''''''''''''''''''''''''
		!Caso mu=m=0 (1)
		!'''''''''''''''''''''''''''''''''''
		
		!BACKWARD
		
		!Primo valore per la backward recursion
		v_aq(0)=f_a0(m,n,mu,nu) 
	
		!Backward recursion
		uno_b_do: DO q=1,qmax
					
					p=n+nu-2.0D0*REAL(q,dbl)
					c0=f_alpha(n,nu,p+1.0D0)
					c1=f_alpha(n,nu,p+2.0D0)
					
					v_aq(q)=(c1/c0)*v_aq(q-1)
					
					!Vedo se il q-esimo valore e' zero
					v_zero_if_1: IF (v_zero(q-1)==1) THEN
						zg_if_1: IF (ABS(v_aq(q)/v_aq(q-1))<10.0D0**(-tol)) THEN
							v_aq(q)=0.0D0
							v_zero(q)=0
						END IF zg_if_1 
					ELSE IF ((v_zero(q-1)==0).AND.(v_zero(q-2)/=0)) THEN
						zg_if1_1: IF (ABS(v_aq(q)/v_aq(q-2))<10.0D0**(-tol)) THEN
							v_aq(q)=0.0D0
							v_zero(q)=0
						END IF zg_if1_1
					END IF v_zero_if_1
					
		END DO uno_b_do
		
	ELSE IF ((INT(mu,lo)==INT(m,lo)) .AND. (INT(nu,lo)==INT(n,lo))) THEN
	
		!'''''''''''''''''''''''''''''''''''
		!Caso mu=m e nu=n (2)
		!'''''''''''''''''''''''''''''''''''
		
		!BACKWARD
		
		!Primo valore per la backward recursion
		v_aq(0)=f_a0(m,n,mu,nu) 
	
		!Backward recursion
		due_b_do: DO q=1,qmax
					
					!Calcolo pre-coefficienti
					p=n+nu-2.0D0*REAL(q,dbl)
					p1=p-m-mu
					p2=p+m+mu
					
					!Calcolo coefficienti ricorsione
					c0=(p+2.0D0)*(p1+1.0D0)*f_alpha(n,nu,p+1.0D0)
					c1=(p+1.0D0)*(p2+2.0D0)*f_alpha(n,nu,p+2.0D0)
					
					!Ricorsione
					v_aq(q)=(c1/c0)*v_aq(q-1)
					
					!Vedo se il q-esimo valore e' zero
					v_zero_if_2: IF (v_zero(q-1)==1) THEN
						zg_if_2: IF (ABS(v_aq(q)/v_aq(q-1))<10.0D0**(-tol)) THEN
							v_aq(q)=0.0D0
							v_zero(q)=0
						END IF zg_if_2 
					ELSE IF ((v_zero(q-1)==0).AND.(v_zero(q-2)/=0)) THEN
						zg_if1_2: IF (ABS(v_aq(q)/v_aq(q-2))<10.0D0**(-tol)) THEN
							v_aq(q)=0.0D0
							v_zero(q)=0
						END IF zg_if1_2
					END IF v_zero_if_2
					
		END DO due_b_do
			
	ELSE IF (INT(mu,lo)==INT(-m,lo)) THEN		
	
		!'''''''''''''''''''''''''''''''''''
		!Caso mu=-m (3)
		!'''''''''''''''''''''''''''''''''''
		
		!-----------------------------------
		!BACKWARD
		!-----------------------------------
		
		!Primo valore per la backward recursion
		v_aq(0)=f_a0(m,n,mu,nu)
		v_aq(1)=f_a1norm(m,n,mu,nu)*v_aq(0) 
	
		!Controllo gli zeri
		zg_if_3_0: IF (ABS(v_aq(1)/v_aq(0))<10.0D0**(-tol)) THEN
					v_aq(1)=0.0D0
					v_zero(1)=0
		END IF zg_if_3_0
	
		!Backward recursion
		tre_b_do: DO q=2,qmax
					
					!Calcolo pre-coefficienti
					p=n+nu-2.0D0*REAL(q,dbl)
				
					!Calcolo coefficienti ricorsione
					c0=f_alpha(n,nu,p+1.0D0)
					c1=4.0D0*(m**2)+f_alpha(n,nu,p+2.0D0)+f_alpha(n,nu,p+3.0D0)
					c2=-f_alpha(n,nu,p+4.0D0)
					
					!Ricorsione
					v_aq(q)=(c1/c0)*v_aq(q-1)+(c2/c0)*v_aq(q-2)
					
					!Vedo se il q-esimo valore e' zero
					v_zero_if_3: IF (v_zero(q-1)==1) THEN
						zg_if_3: IF (ABS(v_aq(q)/v_aq(q-1))<10.0D0**(-tol)) THEN
							v_aq(q)=0.0D0
							v_zero(q)=0
						END IF zg_if_3 
					ELSE IF ((v_zero(q-1)==0).AND.(v_zero(q-2)/=0)) THEN
						zg_if1_3: IF (ABS(v_aq(q)/v_aq(q-2))<10.0D0**(-tol)) THEN
							v_aq(q)=0.0D0
							v_zero(q)=0
						END IF zg_if1_3
					END IF v_zero_if_3
					
		END DO tre_b_do
		
		!--------------------------------
		!FORWARD
		!--------------------------------
		
		!Primo valore per la forward recursion,errore relativo e suo swap
		aq_fwd=f_aqmax(m,n,mu,nu,qmax)
		res=ABS(aq_fwd-v_aq(qmax))/ABS(aq_fwd)
		
		
		!Se non ho precisione, sostituisco i valori
		tre_f_if: IF (res>(10.0D0**(-prec))) THEN 
		
			v_aq(qmax)=aq_fwd
			qi=1
			
			!Entro nel ciclo della sostituzione valori
			tre_f_do: DO q=qmax-1,0,-1
				
				tre_q_case:SELECT CASE (qmax-q)
				
				CASE(1) tre_q_case	!q=qmax-1
					
					!Calcolo v_aq(qmax-1)
					p=n+nu-2.0D0*REAL(q+2,dbl)
					c1=4.0D0*(m**2)+f_alpha(n,nu,p+2.0D0)+f_alpha(n,nu,p+3.0D0)
					c2=-f_alpha(n,nu,p+4.0D0)
					aq_fwd=-(c1/c2)*v_aq(qmax)
					
					z_3_1_case: SELECT CASE (v_zero(q))
					CASE (0) z_3_1_case
						v_aq(q)=0.0D0
					CASE (1) z_3_1_case
						res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)
					END SELECT z_3_1_case
					
				CASE DEFAULT tre_q_case !Per tutti gli altri q
				
					!Calcolo v_aq(qmax-1)
					p=n+nu-2.0D0*REAL(q+2.0D0,dbl)
					c0=f_alpha(n,nu,p+1.0D0)
					c1=4.0D0*(m**2)+f_alpha(n,nu,p+2.0D0)+f_alpha(n,nu,p+3.0D0)
					c2=-f_alpha(n,nu,p+4.0D0)
					aq_fwd=-(c1/c2)*v_aq(q+1)+(c0/c2)*v_aq(q+2)
				
					z_3_2_case: SELECT CASE (v_zero(q))
					CASE (0) z_3_2_case
						v_aq(q)=0.0D0
					CASE (1) z_3_2_case
						res=ABS(aq_fwd-v_aq(q))/ABS(aq_fwd)
					END SELECT z_3_2_case
				
				END SELECT tre_q_case

			!Adesso se la precisione e' raggiunta esco dal ciclo, se no sostituisco e rimango
			IF ((res<(10.0D0**(-prec))) .OR. (q==0) .OR. (ABS(aq_fwd)<ABS(v_aq(q+1)))) EXIT
			
			!Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
			v_aq(q)=aq_fwd
			qi=q
			
			END DO tre_f_do

		! Check sul ciclo di sostituzione
		error_if1: IF (q==0) THEN
					WRITE(*,*)
					WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu caso mu=-m:"
					WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
					WRITE(*,*) "e forward recursion non e' stata raggiunta"
					WRITE(*,*)
					error=1
					RETURN
		END IF error_if1
		
		
		END IF tre_f_if
			
	ELSE
	
	
	
	
		!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		!Caso generale (4)
		!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		
		!---------------------------------------------------------------------------------------------------------
		!BACKWARD
		!---------------------------------------------------------------------------------------------------------
		
		!........................................................
		!Calcolo direttamente i primi due valori della ricorsione
		!........................................................
		v_aq(0)=f_a0(m,n,mu,nu)
		v_aq(1)=v_aq(0)*f_a1norm(m,n,mu,nu)
	
		!Vedo se il secondo valore e' zero
		zg1_if: IF (ABS(v_aq(1)/v_aq(0))<10.0D0**(-tol)) THEN
			v_aq(1)=0.0D0
			v_zero(1)=0
		END IF zg1_if
	
	
	
		!...........................................................
		!Calcolo il terzo valore della ricorsione in funzione di Ap4
		!...........................................................
		!Inizializzo i valori comuni per i coefficienti
		p=n+nu-2.0D0*REAL(2,dbl)
		p1=p-m-mu
		p2=p+m+mu
		alphap1=f_alpha(n,nu,p+1.0D0)
		alphap2=f_alpha(n,nu,p+2.0D0)
		Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
		Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
		Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
		
		!Con questo if decido se mi serve la ricorsione a 3 o 4 termini
		Ap4_2_if: IF (Ap4==0) THEN
		
			!Calcolo i restanti valori preliminari
			Ap5=f_Ap(m,n,mu,nu,p+5.0D0)
			Ap6=f_Ap(m,n,mu,nu,p+6.0D0)
			alphap5=f_alpha(n,nu,p+5.0D0)
			alphap6=f_alpha(n,nu,p+6.0D0)
			
			!Calcolo i coefficienti per la ricorsione ma non c3 perche' qui e solo qui non mi serve
			c0=(p+2.0D0)*(p+3.0D0)*(p+5.0D0)*(p1+1.0D0)*(p1+2.0D0)*(p1+4.0D0)*Ap6*alphap1
			c1=(p+5.0D0)*(p1+4.0D0)*Ap6*(Ap2*Ap3 + (p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*alphap2)
			c2=(p+2.0D0)*(p2+3.0D0)*Ap2*(Ap5*Ap6 + (p+4.0D0)*(p+6.0D0)*(p1+5.0D0)*(p2+5.0D0)*alphap5)

			!Calcolo il mio coefficiente
			v_aq(2)=(c1/c0)*v_aq(1)+(c2/c0)*v_aq(0)
		
			!Assegno l'indice segnaposto per Ap4=0
			q4=2
		
		ELSE

			!Calcolo i restanti valori preliminari		
			alphap3=f_alpha(n,nu,p+3.0D0)
			alphap4=f_alpha(n,nu,p+4.0D0)
		
			!Calcolo coefficienti ricorsione
			c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
			c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
			  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
			c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
		
			!Calcolo il mio coefficiente
			v_aq(2)=(c1/c0)*v_aq(1)+(c2/c0)*v_aq(0)
		
		END IF Ap4_2_if
	
		!Vedo se il terzo valore e' zero
		v_zero_if1: IF (v_zero(1)==1) THEN
			zg2_if: IF (ABS(v_aq(2)/v_aq(1))<10.0D0**(-tol)) THEN
				v_aq(2)=0.0D0
				v_zero(2)=0
			END IF zg2_if
		ELSE IF (v_zero(1)==0) THEN
			zg2_if1: IF (ABS(v_aq(2)/v_aq(0))<10.0D0**(-tol)) THEN
				v_aq(2)=0.0D0
				v_zero(2)=0
			END IF zg2_if1
		END IF v_zero_if1
	
	
		!...........................................................
		!Calcolo i restanti valori nel loop
		!...........................................................
		gen_bwd_do: DO q=3,qmax
		
			!Inizializzo i valori comuni per i coefficienti
			p=n+nu-2.0D0*REAL(q,dbl)
			p1=p-m-mu
			p2=p+m+mu
			alphap1=f_alpha(n,nu,p+1.0D0)
			alphap2=f_alpha(n,nu,p+2.0D0)
			Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
			Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
			Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
			
			!Con questo if decido se mi serve la ricorsione a 3 o 4 termini
			Ap4_bwd_if: IF (Ap4==0) THEN
			
				!Calcolo i restanti valori preliminari
				Ap5=f_Ap(m,n,mu,nu,p+5.0D0)
				Ap6=f_Ap(m,n,mu,nu,p+6.0D0)
				alphap5=f_alpha(n,nu,p+5.0D0)
				alphap6=f_alpha(n,nu,p+6.0D0)
				
				!Calcolo i coefficienti per la ricorsione ma non c3 perche' qui e solo qui non mi serve
				c0=(p+2.0D0)*(p+3.0D0)*(p+5.0D0)*(p1+1.0D0)*(p1+2.0D0)*(p1+4.0D0)*Ap6*alphap1
				c1=(p+5.0D0)*(p1+4.0D0)*Ap6*(Ap2*Ap3 + (p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*alphap2)
				c2=(p+2.0D0)*(p2+3.0D0)*Ap2*(Ap5*Ap6 + (p+4.0D0)*(p+6.0D0)*(p1+5.0D0)*(p2+5.0D0)*alphap5)
				c3=-(p+2.0D0)*(p+4.0D0)*(p+5.0D0)*(p2+3.0D0)*(p2+5.0D0)*(p2+6.0D0)*Ap2*alphap6
	
				!Calcolo il mio coefficiente
				v_aq(q)=(c1/c0)*v_aq(q-1)+(c2/c0)*v_aq(q-2)+(c3/c0)*v_aq(q-3)
			
				!Assegno l'indice segnaposto per Ap4=0
				q4=q
			
			ELSE
	
				!Calcolo i restanti valori preliminari		
				alphap3=f_alpha(n,nu,p+3.0D0)
				alphap4=f_alpha(n,nu,p+4.0D0)
			
				!Calcolo coefficienti ricorsione
				c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
				c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
				  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
				c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
			
				!Calcolo il mio coefficiente
				v_aq(q)=(c1/c0)*v_aq(q-1)+(c2/c0)*v_aq(q-2)
			
			END IF Ap4_bwd_if
			
			!Vedo se il q-esimo valore e' zero
			v_zero_ifq: IF (v_zero(q-1)==1) THEN
				zgq_if: IF (ABS(v_aq(q)/v_aq(q-1))<10.0D0**(-tol)) THEN
					v_aq(q)=0.0D0
					v_zero(q)=0
				END IF zgq_if
			ELSE IF ((v_zero(q-1)==0).AND.(v_zero(q-2)/=0)) THEN
				zgq_if1: IF (ABS(v_aq(q)/v_aq(q-2))<10.0D0**(-tol)) THEN
					v_aq(q)=0.0D0
					v_zero(q)=0
				END IF zgq_if1
			END IF v_zero_ifq
			
		END DO gen_bwd_do
		
		
	
	
	
	
		!---------------------------------------------------------------------------------
		!FORWARD
		!---------------------------------------------------------------------------------
		
		!Calcolo pmin,Apmin e la mia variabile logica
		pmin=n+nu-2.0D0*REAL(qmax,dbl)
		Apmin=f_Ap(m,n,mu,nu,pmin)
		test=((INT(Apmin,lo)==0) .AND. &
		    &((INT(pmin,lo)==INT(m+mu+1.0D0,lo)).OR.(INT(pmin,lo)==INT(-m-mu+1.0D0,lo))))
		
		!........................................................
		!Se la mia variabile logica e' vera, Faccio il mio conto
		!........................................................
		Apmin_if: IF (test) THEN

			!Il valore per qmax allora e' zero		
			v_aq(qmax)=0.0D0
			
			!Calcolo il secondo valore, e se la precisione e' raggiunta esco
			aq_fwd=f_aqmax_1(m,n,mu,nu,qmax)
			res=ABS(aq_fwd-v_aq(qmax-1))/ABS(aq_fwd)
			IF (res<(10.0D0**(-prec))) THEN
				RETURN
			END IF
		
			!Assegno il secondo valore e faccio il ciclo
			v_aq(qmax-1)=aq_fwd
			qi=1
			
			Apmin_do: DO q=qmax,2,-1
			
				!Calcolo pre-coefficienti
				p=n+nu-2.0D0*REAL(q,dbl)
				p1=p-m-mu
				p2=p+m+mu
				alphap1=f_alpha(n,nu,p+1.0D0)
				alphap2=f_alpha(n,nu,p+2.0D0)
				alphap3=f_alpha(n,nu,p+3.0D0)
				alphap4=f_alpha(n,nu,p+4.0D0)
				Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
				Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
				Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
				
				!Calcolo coefficienti ricorsione
				c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
				c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
				  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
				c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
						
				!Ricorsione e residuo
				aq_fwd=-(c1/c2)*v_aq(q-1)+(c0/c2)*v_aq(q)
				res=ABS(aq_fwd-v_aq(q-2))/ABS(aq_fwd)
				
				IF (res<(10.0D0**(-prec))) EXIT
				
				v_aq(q-2)=aq_fwd
				qi=q-2
				
			END DO Apmin_do		
			
			! Check sul ciclo di sostituzione
			Apmin_error_if1: IF (qi==0) THEN
						WRITE(*,*)
						WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu, caso generale, Apmin=0:"
						WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
						WRITE(*,*) "e forward recursion non e' stata raggiunta"
						WRITE(*,*)
						error=1
						RETURN
			END IF Apmin_error_if1
			
			!Esco dalla subroutine gaunt_xu
			RETURN
			
		END IF Apmin_if
		
		!..........................................................................
		!CASO GENERALE PER LA FORWARD RECURRENCE
		!..........................................................................
		
		!Primo valore per la forward recursion,errore relativo e suo swap
		aq_fwd=f_aqmax(m,n,mu,nu,qmax)
		res=ABS(aq_fwd-v_aq(qmax))/ABS(aq_fwd)
		qi=1
		
		gen_f_if: IF (res>(10.0D0**(-prec))) THEN 
		!Se non ho precisione, sostituisco i valori
		
			v_aq(qmax)=aq_fwd
			
			qi=qmax-1
			
			!Entro nel ciclo della sostituzione valori
			gen_f_do: DO 
				
				gen_q_case:SELECT CASE (qmax-qi)
				
									!$$$$$$$$$$$$$$$$
				CASE(1) gen_q_case	!q=qmax-1
									!$$$$$$$$$$$$$$$$
				
					
					!Calcolo Ap4 per qi+2 per vedere quale schema usare
					p=n+nu-2.0D0*REAL(qi+2,dbl)
					Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
					
					!Scelgo la ricorrenza a seconda del valore di Ap4
					Ap4_q1_if: IF (Ap4==0) THEN
					
						!Calcolo aq secondo la ricorrenza a 4 termini: uso qi+3 perche' il termine piu' alto e'
						!maggiore di 3 unita' rispetto a qi, pur essendo nullo e non comparendo nella ricorsione
						p=n+nu-2.0D0*REAL(qi+3,dbl)
						p1=p-m-mu
						p2=p+m+mu
						alphap5=f_alpha(n,nu,p+5.0D0)
						alphap6=f_alpha(n,nu,p+6.0D0)
						Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
						Ap5=f_Ap(m,n,mu,nu,p+5.0D0)
						Ap6=f_Ap(m,n,mu,nu,p+6.0D0)
						c2=(p+2.0D0)*(p2+3.0D0)*Ap2*(Ap5*Ap6 + (p+4.0D0)*(p+6.0D0)*(p1+5.0D0)*(p2+5.0D0)*alphap5)
						c3=-(p+2.0D0)*(p+4.0D0)*(p+5.0D0)*(p2+3.0D0)*(p2+5.0D0)*(p2+6.0D0)*Ap2*alphap6
						aq_fwd=-(c2/c3)*v_aq(qi+1) 
					
						!A seconda che il mio valore sia 0 o meno confronto i valori
						zAp41_case:SELECT CASE (v_zero(qi))
						CASE (0) zAp41_case
							v_aq(qi)=0.0D0
						CASE (1) zAp41_case
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
							IF (res<10.0D0**(-prec)) EXIT gen_f_do
						END SELECT zAp41_case
						
						!Qui calcolo il valore successivo dopo aver aggiornato qi:
						!Se v_aq(qi)=0 allora non chiamo cruzan, se no lo chamo e 
						!tengo un solo valore
						qi=qi-1
						
						zcz1_case:SELECT CASE (v_zero(qi))
						CASE (0) zcz1_case
							v_aq(qi)=0.0D0
							qi=qi-1
							CYCLE gen_f_do
						CASE (1) zcz1_case
							CALL gaunt_cz(m,n,mu,nu,qmax,v_aq_cz(qi),error)
							aq_fwd=v_aq_cz(qi)
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
						END SELECT zcz1_case
					
							!-----------------
					ELSE	!Qui Ap4/=0
							!-----------------
					
						!Calcolo aq
						p=n+nu-2.0D0*REAL(qi+2,dbl)
						p1=p-m-mu
						p2=p+m+mu
						alphap2=f_alpha(n,nu,p+2.0D0)
						alphap3=f_alpha(n,nu,p+3.0D0)
						alphap4=f_alpha(n,nu,p+4.0D0)
						Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
						Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
						Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
						c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
						  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
						c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
						aq_fwd=-(c1/c2)*v_aq(qi+1) !E' qui che lo calcolo
					
						!A seconda che il mio valore sia 0 o meno confronto i valori
						zAp4d1_case:SELECT CASE (v_zero(qi))
						CASE (0) zAp4d1_case
							v_aq(qi)=0.0D0
							qi=qi-1
							CYCLE gen_f_do
						CASE (1) zAp4d1_case
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
						END SELECT zAp4d1_case
					
					END IF Ap4_q1_if




					
					
									!$$$$$$$$$$$$$$$$
				CASE(2) gen_q_case	!q=qmax-2	
									!$$$$$$$$$$$$$$$$






					!Calcolo Ap4 per qi+2 per vedere quale schema usare
					p=n+nu-2.0D0*REAL(qi+2,dbl)
					Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
					
					!Scelgo la ricorrenza a seconda del valore di Ap4
					Ap4_q2_if: IF (Ap4==0) THEN
					
						!Calcolo aq secondo la ricorrenza a 4 termini: uso qi+3 perche' il termine piu' alto e'
						!maggiore di 3 unita' rispetto a qi, pur essendo nullo e non comparendo nella ricorsione
						p=n+nu-2.0D0*REAL(qi+3,dbl)
						p1=p-m-mu
						p2=p+m+mu
						alphap2=f_alpha(n,nu,p+2.0D0)
						alphap5=f_alpha(n,nu,p+5.0D0)
						alphap6=f_alpha(n,nu,p+6.0D0)
						Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
						Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
						Ap5=f_Ap(m,n,mu,nu,p+5.0D0)
						Ap6=f_Ap(m,n,mu,nu,p+6.0D0)
						c1=(p+5.0D0)*(p1+4.0D0)*Ap6*(Ap2*Ap3 + (p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*alphap2)
						c2=(p+2.0D0)*(p2+3.0D0)*Ap2*(Ap5*Ap6 + (p+4.0D0)*(p+6.0D0)*(p1+5.0D0)*(p2+5.0D0)*alphap5)
						c3=-(p+2.0D0)*(p+4.0D0)*(p+5.0D0)*(p2+3.0D0)*(p2+5.0D0)*(p2+6.0D0)*Ap2*alphap6
						aq_fwd=-(c1/c3)*v_aq(qi+2) -(c2/c3)*v_aq(qi+1) 
					
						!A seconda che il mio valore sia 0 o meno confronto i valori
						zAp42_case:SELECT CASE (v_zero(qi))
						CASE (0) zAp42_case
							v_aq(qi)=0.0D0
						CASE (1) zAp42_case
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
							IF (res<10.0D0**(-prec)) EXIT gen_f_do
						END SELECT zAp42_case
						
						!Qui calcolo il valore successivo dopo aver aggiornato qi:
						!Se v_aq(qi)=0 allora non chiamo cruzan, se no lo chamo e 
						!tengo un solo valore
						qi=qi-1
						
						zcz2_case:SELECT CASE (v_zero(qi))
						CASE (0) zcz2_case
							v_aq(qi)=0.0D0
							qi=qi-1
							CYCLE gen_f_do
						CASE (1) zcz2_case
							CALL gaunt_cz(m,n,mu,nu,qmax,v_aq_cz(qi),error)
							aq_fwd=v_aq_cz(qi)
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
						END SELECT zcz2_case
					
							!-----------------
					ELSE	!Qui Ap4/=0
							!-----------------
					
						!Calcolo aq
						p=n+nu-2.0D0*REAL(qi+2,dbl)
						p1=p-m-mu
						p2=p+m+mu
						alphap2=f_alpha(n,nu,p+2.0D0)
						alphap3=f_alpha(n,nu,p+3.0D0)
						alphap4=f_alpha(n,nu,p+4.0D0)
						Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
						Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
						Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
						c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
						c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
						  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
						c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
						aq_fwd=(c0/c2)*v_aq(qi+2)-(c1/c2)*v_aq(qi+1) !E' qui che lo calcolo
					
						!A seconda che il mio valore sia 0 o meno confronto i valori
						zAp4d2_case:SELECT CASE (v_zero(qi))
						CASE (0) zAp4d2_case
							v_aq(qi)=0.0D0
							qi=qi-1
							CYCLE gen_f_do
						CASE (1) zAp4d2_case
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
						END SELECT zAp4d2_case
					
					END IF Ap4_q2_if



					
										!$$$$$$$$$$$$$$$$$$$$$$
				CASE DEFAULT gen_q_case !Per tutti gli altri q
										!$$$$$$$$$$$$$$$$$$$$$$



				
				
					!Calcolo Ap4 per qi+2 per vedere quale schema usare
					p=n+nu-2.0D0*REAL(qi+2,dbl)
					Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
					
					!Scelgo la ricorrenza a seconda del valore di Ap4
					Ap4_qq_if: IF (Ap4==0) THEN
					
						!Calcolo aq secondo la ricorrenza a 4 termini: uso qi+3 perche' il termine piu' alto e'
						!maggiore di 3 unita' rispetto a qi, pur essendo nullo e non comparendo nella ricorsione
						p=n+nu-2.0D0*REAL(qi+3,dbl)
						p1=p-m-mu
						p2=p+m+mu
						alphap2=f_alpha(n,nu,p+2.0D0)
						alphap5=f_alpha(n,nu,p+5.0D0)
						alphap6=f_alpha(n,nu,p+6.0D0)
						Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
						Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
						Ap5=f_Ap(m,n,mu,nu,p+5.0D0)
						Ap6=f_Ap(m,n,mu,nu,p+6.0D0)
						c0=(p+2.0D0)*(p+3.0D0)*(p+5.0D0)*(p1+1.0D0)*(p1+2.0D0)*(p1+4.0D0)*Ap6*alphap1
						c1=(p+5.0D0)*(p1+4.0D0)*Ap6*(Ap2*Ap3 + (p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*alphap2)
						c2=(p+2.0D0)*(p2+3.0D0)*Ap2*(Ap5*Ap6 + (p+4.0D0)*(p+6.0D0)*(p1+5.0D0)*(p2+5.0D0)*alphap5)
						c3=-(p+2.0D0)*(p+4.0D0)*(p+5.0D0)*(p2+3.0D0)*(p2+5.0D0)*(p2+6.0D0)*Ap2*alphap6
						aq_fwd=(c0/c3)*v_aq(qi+3)-(c1/c3)*v_aq(qi+2) -(c2/c3)*v_aq(qi+1) 
					
						!A seconda che il mio valore sia 0 o meno confronto i valori
						zAp4q_case:SELECT CASE (v_zero(qi))
						CASE (0) zAp4q_case
							v_aq(qi)=0.0D0
						CASE (1) zAp4q_case
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
							IF (res<10.0D0**(-prec)) EXIT gen_f_do
						END SELECT zAp4q_case
						
						!Qui calcolo il valore successivo dopo aver aggiornato qi:
						!Se v_aq(qi)=0 allora non chiamo cruzan, se no lo chiamo e 
						!tengo un solo valore.L'if c'e' per non far sballare qi
						
						qi_if:IF (qi>0) THEN
						
							qi=qi-1
						
							zczq_case:SELECT CASE (v_zero(qi))
							CASE (0) zczq_case
								v_aq(qi)=0.0D0
								qi=qi-1
								CYCLE gen_f_do
							CASE (1) zczq_case
								CALL gaunt_cz(m,n,mu,nu,qmax,v_aq_cz(qi),error)
								aq_fwd=v_aq_cz(qi)
								res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
							END SELECT zczq_case
					
						END IF qi_if
					
							!-----------------
					ELSE	!Qui Ap4/=0
							!-----------------
					
						!Calcolo aq
						p=n+nu-2.0D0*REAL(qi+2,dbl)
						p1=p-m-mu
						p2=p+m+mu
						alphap2=f_alpha(n,nu,p+2.0D0)
						alphap3=f_alpha(n,nu,p+3.0D0)
						alphap4=f_alpha(n,nu,p+4.0D0)
						Ap2=f_Ap(m,n,mu,nu,p+2.0D0)
						Ap3=f_Ap(m,n,mu,nu,p+3.0D0)
						Ap4=f_Ap(m,n,mu,nu,p+4.0D0)
						c0=(p+2.0D0)*(p+3.0D0)*(p1+1.0D0)*(p1+2.0D0)*Ap4*alphap1
						c1=Ap2*Ap3*Ap4+(p+1.0D0)*(p+3.0D0)*(p1+2.0D0)*(p2+2.0D0)*Ap4*alphap2+ &
						  &(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*Ap2*alphap3
						c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*Ap2*alphap4
						aq_fwd=(c0/c2)*v_aq(qi+2)-(c1/c2)*v_aq(qi+1) !E' qui che lo calcolo
					
						!A seconda che il mio valore sia 0 o meno confronto i valori
						zAp4dq_case:SELECT CASE (v_zero(qi))
						CASE (0) zAp4dq_case
							v_aq(qi)=0.0D0
							qi=qi-1
							CYCLE gen_f_do
						CASE (1) zAp4dq_case
							res=ABS(aq_fwd-v_aq(qi))/ABS(aq_fwd)
						END SELECT zAp4dq_case
					
					END IF Ap4_qq_if
				
				END SELECT gen_q_case

			!Adesso se la precisione e' raggiunta esco dal ciclo, se no sostituisco e rimango
			IF ((res<(10.0D0**(-prec))) .OR. (qi==0) .OR. (ABS(aq_fwd)<ABS(v_aq(qi+1)))) EXIT
			
			!Sono nel ciclo, allora sostituisco eaggiorno indice e residuo
			v_aq(qi)=aq_fwd
			qi=qi-1
			
			END DO gen_f_do

			! Check sul ciclo di sostituzione
			gen_error_if1: IF (qi==0) THEN
						WRITE(*,*)
						WRITE(*,*) "Si e' verificato un errore nella subroutine gaunt_xu,caso generale:"
						WRITE(*,*) "la precisione richiesta per i coefficienti di Gaunt nella backward"
						WRITE(*,*) "e forward recursion non e' stata raggiunta"
						WRITE(*,*)
						error=1
						RETURN
			END IF gen_error_if1
		
		
		END IF gen_f_if

	END IF big_if
		
END SELECT qmax_case





END SUBROUTINE gaunt_xu















! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! LEGENDRE SUBROUTINES: COEFF. DI LEGENDRE SULLA FALSARIGA DI NUMERICAL RECIPES
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!******************************************************************************
!1) SUBROUTINE legendre: calcolo le mie serie di coefficienti di legendre
!******************************************************************************
SUBROUTINE legendre(nmin,nmax,m,mm,theta,v_leg,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nmin,nmax,m,mm				! Valori massimi e minimi di n, m e modulo di m
REAL(dbl), INTENT(IN) :: theta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(nmin:nmax), INTENT(OUT) :: v_leg	! Vettore di output valori funzioni di Legendre

!Dichiarazione variabili interne 
INTEGER(lo), PARAMETER :: prec=7						! Precisione su theta per 0 o pigreco
REAL(dbl) :: x											! x=Cos(thetain)
REAL(dbl) :: fact,somx2									! Semifattoriale e (1-x^2)			
REAL(dbl), DIMENSION(mm:nmax) :: v_legin			 	! Vettore interno legendre
REAL(dbl) :: mr,lr,logw,sig								! Indici reali,log dell'exp e segno
INTEGER(lo) :: i,l									! Indici


!Subroutine vera e propria
x=COS(theta)
error=0


!Controllo se nmin e nmax sono ok
n_check: IF (nmax<nmin) THEN
     WRITE(*,*) "nmax<nmin, la procedura Legendre si ferma"
     error=1
     RETURN
END IF n_check


!Controllo se ABS(m) e nmax sono ok
nm_check: IF (nmax<mm) THEN
     WRITE(*,*) "nmax<|m|, la procedura Legendre si ferma"
     error=1
     RETURN
END IF nm_check


!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((theta>PI_D) .OR. (theta<0)) THEN
     WRITE(*,*) "Theta e' al di fuori di [0,Pi], la procedura Legendre si ferma"
     error=1
     RETURN
END IF thetain_check


!Controllo se thetain e' zero, perche' cosi' i conti sono semplificati
!e non devo chiamare la routine
!m_if: IF (mm==0) THEN 

!     theta1_if: IF (ABS(theta)<10.0D0**(-prec)) THEN

!          !Theta=0
!          v_leg=1.0D0
!          RETURN

!          !Questo gruppo e' commentato perche' produce un errore, il caso e' da approfondire
!          !				
!          ! 			ELSE IF (ABS(Pi_D-theta)<10.0D0**(-prec)) THEN
!          ! 			
!          ! 				!Theta=+-Pigreco 				
!          ! 				DO i=nmin,nmax
!          ! 					v_leg(i)=-1.0D0**(i)
!          ! 				END DO
!          ! 				RETURN

!     END IF theta1_if

!ELSE 
!     theta2_if: IF (theta<10.0D0**(-prec)) THEN

!          !			La linea sottostante e' commentata per sicurezza, e uso quella
!          !			Immediatamente sottostante 			
!          !  			theta2_if: IF ((theta<10.0D0**(-prec)) .OR. (ABS(Pi_D-theta)<10.0D0**(-prec)) ) THEN

!          !Theta=0 or +- Pigreco
!          v_leg=0.0D0 		
!          RETURN

!     END IF theta2_if

!END IF m_if

!Alloco il vettore interno
m_check: IF (mm>nmin) THEN
     v_leg(nmin:(mm-1))=0.0D0
END IF m_check

!WRITE(*,*) theta,ABS(PIO2_D-theta)/PIO2_D


mr=REAL(mm,dbl)
!Questo if lo faccio perche' Pi/2 e' un valore particolare che mi cambia la ricorsione
!pi2_if: IF ((ABS(PIO2_D-theta)/PIO2_D)<10.0D0**(-7)) THEN

!     !	WRITE(*,*) "thetaqui"

!     !Inizializzo il vettore a zero,mi va sempre bene,indipendentemente dal numero di termini
!     v_legin=0.0D0

!     !Calcolo il primo termine della ricorsione: VA SEMPRE FATTO!!!
!     v_legin(mm)=fact2(2.0D0*mr-1.0D0)

!     !L'eventuale secondo termine e' gia ben inizializzato
!     rec_if: IF (mm+1<nmax) THEN

!          leg_do: DO l=mm+2,nmax,2
!               lr=REAL(l,dbl)
!               v_legin(l)=-(lr+mr-1.0D0)*v_legin(l-2)/(lr-mr)
!          END DO leg_do

!     END IF rec_if

!ELSE

!Calcolo il primo termine della ricorsione: VA SEMPRE FATTO!!!
v_legin(mm)=1.0D0
mzero_if:IF (mm>0) THEN

  somx2=SQRT((1.0D0-x)*(1.0D0+x))
  fact=1.0D0

  pmm_do: DO i=1,mm
	   v_legin(mm)=v_legin(mm)*fact*somx2
	   fact=fact+2.0D0
  END DO pmm_do

END IF mzero_if

!Procedo con la ricorsione
big_if: IF (mm+1==nmax) THEN

  v_legin(mm+1)=x*(2*mr+1.0D0)*v_legin(mm)

ELSE IF (mm+1<nmax) THEN

  v_legin(mm+1)=x*(2*mr+1.0D0)*v_legin(mm)

  leg_do1: DO l=mm+2,nmax
	   lr=REAL(l,dbl)
	   v_legin(l)=(x*(2*lr-1.0D0)*v_legin(l-1)-(lr+mr-1.0D0)*v_legin(l-2))/(lr-mr)
  END DO leg_do1

END IF big_if

!END IF pi2_if

! Assegno il valore a seconda del segno di m
mm_if: IF (mm<nmin) THEN

     sig_if1: IF (m<0) THEN

          sig=(-1.0D0)**mm

          assign_do1: DO l=nmin,nmax

               lr=REAL(l,dbl)
               logw=lnf(lr-mr,error)-lnf(lr+mr,error)
               v_leg(l)=sig*EXP(logw)*v_legin(l)
               !v_leg(l)=sig*m_fact(l-mm,l+mm)*v_legin(l)
               ! WRITE(*,*) "leg1",EXP(logw),m_fact(l-mm,l+mm)


          END DO assign_do1

     ELSE

          v_leg(nmin:nmax)=v_legin(nmin:nmax)

     END IF sig_if1

ELSE

     sig_if2: IF (m<0) THEN

          sig=(-1.0D0)**mm

          assign_do2: DO l=mm,nmax

               lr=REAL(l,dbl)
               logw=lnf(lr-mr,error)-lnf(lr+mr,error)
               v_leg(l)=sig*EXP(logw)*v_legin(l)
               !v_leg(l)=sig*m_fact(l-mm,l+mm)*v_legin(l)
               ! WRITE(*,*) "leg2",EXP(logw),m_fact(l-mm,l+mm)

          END DO assign_do2

     ELSE

          v_leg(mm:nmax)=v_legin(mm:nmax)

     END IF sig_if2

END IF mm_if

END SUBROUTINE legendre



!******************************************************************************
!2) SUBROUTINE pi_mn: calcolo le funzioni angolari Pi_mn
!******************************************************************************
SUBROUTINE pi_mn(nstop,theta,v_pimn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_pimn	! Vettore di output valori funzioni pi_mn

!Dichiarazione variabili interne 
INTEGER(lo), PARAMETER :: prec=7						! Precisione su theta per 0 o pigreco
INTEGER(lo) :: idx1,idx2,idx3,m,n					! Indice per riordinare i miei calcoli
REAL(dbl) :: fact,mr,nr									! Semifattoriale e m reale


!Subroutine vera e propria
error=0

!!$WRITE(*,*) 
!!$WRITE(*,*) "inside pi_mn"
!!$WRITE(*,*) "theta",theta
!!$WRITE(*,*) "nstop",nstop

!Controllo se nstop>0 
nm_check: IF (nstop<0) THEN
     WRITE(*,*) "nstop<0, la procedura pi_mn si ferma"
     error=1
     RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((theta>PI_D) .OR. (theta<0.0D0)) THEN
     WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
     error=1
     RETURN
END IF thetain_check

!Inizializzo a zero tutto il vettore
v_pimn=0.0D0

!Il calcolo dei valori cambia a seconda che theta sia 0 o pigreco
!theta_if :IF (theta<10.0D0**(-prec)) THEN

!     !Assegno i valori semplificati per theta=0
!     theta0_ndo: DO n=1,nstop
!          theta0_mdo: DO m=-n,n

!               !Calcolo per m=1 o m=-1
!               theta0_mcase: SELECT CASE (m)

!               CASE(-1)

!                    idx1=n*(n+1)+m
!                    v_pimn(idx1)=0.5D0

!               CASE(1)

!                    idx1=n*(n+1)+m
!                    nr=REAL(n,dbl)
!                    v_pimn(idx1)=0.5D0*nr*(nr+1.0D0)

!               END SELECT theta0_mcase

!          END DO theta0_mdo
!     END DO theta0_ndo

!     RETURN

!ELSE IF ((Pi_D-theta)<10.0D0**(-prec)) THEN

!     !Assegno i valori semplificati per theta=pigreco
!     thetaPi_ndo: DO n=1,nstop
!          thetaPi_mdo: DO m=-n,n

!               !Calcolo per m=1 o m=-1
!               thetaPi_mcase: SELECT CASE (m)

!               CASE(-1)

!                    idx1=n*(n+1)+m
!                    v_pimn(idx1)=((-1.0D0)**(n+1)) * 0.5D0

!               CASE(1)

!                    idx1=n*(n+1)+m
!                    nr=REAL(n,dbl)
!                    v_pimn(idx1)=((-1.0D0)**(n+1)) * 0.5D0 *nr*(nr+1.0D0)

!               END SELECT thetaPi_mcase

!          END DO thetaPi_mdo
!     END DO thetaPi_ndo

!     RETURN

!ELSE IF (ABS(PIO2_D-theta)<10.0D0**(-prec)) THEN

!     !Assegno i valori semplificati per theta=pigreco/2
!     pi2_mdo: DO m=1,nstop

!          idx1=m*(m+1)+m
!          mr=REAL(m,dbl)
!          v_pimn(idx1)=mr*fact2(2.0D0*mr-1.0D0)

!          pi2_ndo: DO n=m+2,nstop,2

!               nr=REAL(n,dbl)
!               idx1=(n-2)*(n-1)+m
!               idx2=n*(n+1)+m

!               v_pimn(idx2)=-(nr+mr-1.0D0)*v_pimn(idx1)/(nr-mr)

!          END DO pi2_ndo


!     END DO pi2_mdo

!ELSE	

!Qui c'e' il caso generale
m_do: DO m=1,nstop

  mr=REAL(m,dbl)

  !A seconda del valore di m,devo calcolare 1 o 2 o piu' valori
  m_case: SELECT CASE (nstop-m)

  CASE(0) m_case	

	   !******************************
	   !Caso in cui nstop=m
	   !******************************

	   one_case_ns: SELECT CASE(m)

	   CASE(1)			!m=1,m=nstop 

			idx1=3
			v_pimn(idx1)=1.0D0

	   CASE DEFAULT	!m/=1,m=nstop 

			idx1=nstop*(nstop+2)
			fact=fact2(2.0D0*mr-1.0D0)
			v_pimn(idx1)=mr*fact*(SIN(theta)**(m-1))

	   END SELECT one_case_ns

  CASE(1) m_case 

	   !******************************
	   !Caso in cui nstop=m+1
	   !******************************

	   one_case_ns1: SELECT CASE(m)

	   CASE(1)			!m=1,m=nstop-1 

			!Indici del vettore
			idx1=3
			idx2=7
			!Valori
			v_pimn(idx1)=1.0D0
			v_pimn(idx2)=3.0D0*COS(theta)*v_pimn(idx1)


	   CASE DEFAULT	!m/=1,m=nstop-1 

			!Indici del vettore
			idx1=(nstop**2)-1
			idx2=nstop*(nstop+1)+nstop-1
			!Calcolo i valori per (m,m) e (m,m+1)
			fact=fact2(2.0D0*mr-1.0D0)
			v_pimn(idx1)=mr*fact*(SIN(theta)**(m-1))
			v_pimn(idx2)=(2.0D0*mr+1.0D0)*COS(theta)*v_pimn(idx1)

	   END SELECT one_case_ns1

  CASE DEFAULT m_case

	   !******************************
	   !Caso generale in cui nstop>m+1
	   !******************************

	   !Calcolo gli starting values a seconda che sia m=1 o meno
	   start_case: SELECT CASE(m)

	   CASE(1)			!m=1 

			!Indici del vettore
			idx1=3
			idx2=7
			!Valori
			v_pimn(idx1)=1.0D0
			v_pimn(idx2)=3.0D0*COS(theta)*v_pimn(idx1)


	   CASE DEFAULT	!m/=1,m=nstop-1 

			!Indici del vettore
			idx1=m*(m+2)
			idx2=(m+1)*(m+2)+m
			!Calcolo i valori per (m,m) e (m,m+1)
			fact=fact2(2.0D0*mr-1.0D0)
			v_pimn(idx1)=mr*fact*(SIN(theta)**(m-1))
			v_pimn(idx2)=(2.0D0*mr+1.0D0)*COS(theta)*v_pimn(idx1)

	   END SELECT start_case

	   !Calcolo i restanti valori per ricorsione
	   n_do: DO n=m+2,nstop

			!Sistemo indici reali ed interi
			nr=REAL(n,dbl)
			idx1=(n-2)*(n-1)+m
			idx2=(n-1)*n+m
			idx3=n*(n+1)+m

			!Faccio il conto
			v_pimn(idx3)=((2.0D0*nr-1.0D0)/(nr-mr))*COS(theta)*v_pimn(idx2)- &
				 &((nr+mr-1.0D0)/(nr-mr))*v_pimn(idx1)

	   END DO n_do

  END SELECT m_case

END DO m_do

!END IF theta_if

!Calcolo dei valori per m negativo
nminus_do: DO n=1,nstop
     mminus_do: DO m=1,n

          !Sistemo gli indici
          nr=REAL(n,dbl)
          mr=REAL(m,dbl)
          idx1=n*(n+1)+m
          idx2=n*(n+1)-m

          !Faccio i conti
          fact=lnf(nr-mr,error)-lnf(nr+mr,error) 
          v_pimn(idx2)=((-1.0D0)**(m+1))*EXP(fact)*v_pimn(idx1)
          ! v_pimn(idx2)=((-1.0D0)**(m+1))*m_fact(n-m,n+m)*v_pimn(idx1)
          !WRITE(*,*) "pmn",EXP(fact),m_fact(n-m,n+m)


     END DO mminus_do
END DO nminus_do

END SUBROUTINE pi_mn


!******************************************************************************
!3) SUBROUTINE tau_mn: calcolo le funzioni angolari Tau_mn
!******************************************************************************
SUBROUTINE tau_mn(nstop,theta,v_pimn,v_taumn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_pimn			! Vettore di input valori funzioni pi_mn
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_taumn			! Vettore di output valori funzioni tau_mn

!Dichiarazione variabili interne 
INTEGER(lo), PARAMETER :: prec=7						! Precisione su theta per 0 o pigreco
INTEGER(lo) :: idx1,idx2,idx3,m,n						! Indice per riordinare i miei calcoli e non
REAL(dbl) :: fact,mr,nr,sig								! Semifattoriale e m,n reali,variabile segno


!Subroutine vera e propria
error=0

!Controllo se nstop>0 
nm_check: IF (nstop<0) THEN
     WRITE(*,*) "nstop<0, la procedura pi_mn si ferma"
     error=1
     RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((theta>PI_D) .OR. (theta<0.0D0)) THEN
     WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
     error=1
     RETURN
END IF thetain_check

!Inizializzo a zero tutto il vettore
v_taumn=0.0D0

!****************************************************************
!Calcolo dei valori per theta=0,Pi
!****************************************************************
!theta_if :IF ( (theta<10.0D0**(-prec)) .OR. ((PI_D-theta)<10.0D0**(-prec)) ) THEN

!     !Assegno il segno per i miei conti
!     sig_if: IF (theta<PIO2_D) THEN
!          sig=1.0D0
!     ELSE
!          sig=-1.0D0
!     END IF sig_if

!     !Assegno i valori semplificati per theta=0
!     theta0_ndo: DO n=1,nstop

!          nr=REAL(n,dbl)

!          theta0_mdo: DO m=-n,n

!               !Calcolo l'indice
!               idx1=n*(n+1)+m

!               !Calcolo per m=1 o m=-1
!               theta0_mcase: SELECT CASE (m)

!               CASE(-1)
!                    !m=-1
!                    v_taumn(idx1)=-0.5D0*(sig**(n))

!               CASE(1)
!                    !m=1
!                    v_taumn(idx1)=(sig**(n))*0.5D0*nr*(nr+1.0D0)

!               END SELECT theta0_mcase

!          END DO theta0_mdo
!     END DO theta0_ndo

!     !Chiudo la subroutine
!     RETURN

!END IF theta_if


!****************************************************************
!Calcolo dei valori nel caso generale
!****************************************************************
m_do: DO m=0,nstop

     mr=REAL(m,dbl)

     !A seconda del valore di m,devo calcolare 1 o 2 o piu' valori
     m_case: SELECT CASE (m)

     CASE(0) m_case	

          !******************************
          !Caso in cui m=0
          !******************************

!          m0_thetaif: IF (ABS(PIO2_D-theta)<10.0D0**(-prec)) THEN 

!               !******************************
!               !Caso theta=Pi/2
!               !******************************

!               !Starting value
!               idx1=2
!               v_taumn(idx1)=-1.0D0

!               !Ciclo per il calcolo
!               m0_pi2ndo: DO n=3,nstop,2

!                    !Indici
!                    nr=REAL(n,dbl)
!                    idx1=(n-2)*(n-1)
!                    idx2=n*(n+1)

!                    !Calcolo	
!                    v_taumn(idx2)=-(nr)*v_taumn(idx1)/(nr-1)

!               END DO m0_pi2ndo

!          ELSE

               !******************************
               !Caso generale theta/=Pi/2
               !******************************

               !Starting value
               idx1=2
               v_taumn(idx1)=-SIN(theta)

               !If a seconda del valore di nstop
               nstop_if: IF (nstop>1) THEN

                    !Secondo starting value
                    idx1=6
                    v_taumn(idx1)=-3.0D0*COS(theta)*SIN(theta)

                    !Ciclo per il calcolo
                    m0_ndo: DO n=3,nstop

                         !Assegnazione degli indici
                         nr=REAL(n,dbl)
                         idx1=(n-2)*(n-1)
                         idx2=(n-1)*n
                         idx3=n*(n+1)

                         !Calcolo
                         v_taumn(idx3)=COS(theta)*(2.0D0*nr-1.0D0)*v_taumn(idx2)/(nr-1.0D0) - &
                              & nr*v_taumn(idx1)/(nr-1.0D0)   

                    END DO m0_ndo

               END IF nstop_if

!          END IF m0_thetaif

     CASE DEFAULT 

          !******************************
          !Caso in cui m/=0
          !******************************

!          m_thetaif: IF (ABS(PIO2_D-theta)<10.0D0**(-prec)) THEN

!               !******************************
!               !Caso theta=Pi/2
!               !******************************

!               !Ciclo per l'assegnazione dei valori
!               pi2_ndo: DO n=m+1,nstop,2

!                    !Indici
!                    nr=REAL(n,dbl)
!                    idx1=(n-1)*n+m
!                    idx2=n*(n+1)+m

!                    !Calcolo
!                    v_taumn(idx2)=-(nr+mr)*v_pimn(idx1)/mr

!               END DO pi2_ndo

!          ELSE

               !******************************
               !Caso generale theta/=Pi/2
               !******************************

               !Starting value
               idx1=m*(m+1)+m
               v_taumn(idx1)=COS(theta)*v_pimn(idx1)

               !Ciclo per l'assegnazione dei valori
               ndo: DO n=m+1,nstop

                    !Indici
                    nr=REAL(n,dbl)
                    idx1=(n-1)*n+m
                    idx2=n*(n+1)+m

                    !Calcolo
                    v_taumn(idx2)=COS(theta)*v_pimn(idx2)*(nr/mr)-(nr+mr)*v_pimn(idx1)/mr

               END DO ndo

!          END IF m_thetaif

     END SELECT m_case

END DO m_do

!Calcolo dei valori per m negativo
nminus_do: DO n=1,nstop
     mminus_do: DO m=1,n

          !Sistemo gli indici
          nr=REAL(n,dbl)
          mr=REAL(m,dbl)
          idx1=n*(n+1)+m
          idx2=n*(n+1)-m

          !Faccio i conti
          fact=lnf(nr-mr,error)-lnf(nr+mr,error) 
          v_taumn(idx2)=((-1.0D0)**(m))*EXP(fact)*v_taumn(idx1)
          !v_taumn(idx2)=((-1.0D0)**(m))*m_fact(n-m,n+m)*v_taumn(idx1)
          ! WRITE(*,*) "tmn",EXP(fact),m_fact(n-m,n+m)

     END DO mminus_do
END DO nminus_do


END SUBROUTINE tau_mn












!******************************************************************************
!2) SUBROUTINE pi_mn: calcolo le funzioni angolari Pi_mn
!******************************************************************************
SUBROUTINE pi_mnsca(nstop,theta,v_pimn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_pimn	! Vettore di output valori funzioni pi_mn

!Dichiarazione variabili interne 
INTEGER(lo), PARAMETER :: prec=7						! Precisione su theta per 0 o pigreco
INTEGER(lo) :: idx1,idx2,idx3,m,n					! Indice per riordinare i miei calcoli
REAL(dbl) :: fact,mr,nr									! Semifattoriale e m reale


!Subroutine vera e propria
error=0

!Controllo se nstop>0 
nm_check: IF (nstop<0) THEN
			WRITE(*,*) "nstop<0, la procedura pi_mn si ferma"
            error=1
            RETURN
END IF nm_check

! !Controllo se thetain e' nei suoi limiti
! thetain_check: IF ((theta>PI_D) .OR. (theta<0.0D0)) THEN
! 				WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
!                 error=1
!                 RETURN
! END IF thetain_check

!Inizializzo a zero tutto il vettore
v_pimn=0.0D0

!Il calcolo dei valori cambia a seconda che theta sia 0 o pigreco
!theta_if :IF (theta<10.0D0**(-prec)) THEN

!	!Assegno i valori semplificati per theta=0
!	theta0_ndo: DO n=1,nstop
!		theta0_mdo: DO m=-n,n
		
!			!Calcolo per m=1 o m=-1
!			theta0_mcase: SELECT CASE (m)
			
!			CASE(-1)
			
!				idx1=n*(n+1)+m
!				v_pimn(idx1)=0.5D0
			
!			CASE(1)
			
!				idx1=n*(n+1)+m
!				nr=REAL(n,dbl)
!				v_pimn(idx1)=0.5D0*nr*(nr+1.0D0)
			
!			END SELECT theta0_mcase
	
!		END DO theta0_mdo
!	END DO theta0_ndo
	
!	RETURN

!ELSE IF ((Pi_D-theta)<10.0D0**(-prec)) THEN

!	!Assegno i valori semplificati per theta=pigreco
!	thetaPi_ndo: DO n=1,nstop
!		thetaPi_mdo: DO m=-n,n
		
!			!Calcolo per m=1 o m=-1
!			thetaPi_mcase: SELECT CASE (m)
			
!			CASE(-1)
			
!				idx1=n*(n+1)+m
!				v_pimn(idx1)=((-1.0D0)**(n+1)) * 0.5D0
			
!			CASE(1)
			
!				idx1=n*(n+1)+m
!				nr=REAL(n,dbl)
!				v_pimn(idx1)=((-1.0D0)**(n+1)) * 0.5D0 *nr*(nr+1.0D0)
			
!			END SELECT thetaPi_mcase
	
!		END DO thetaPi_mdo
!	END DO thetaPi_ndo

!	RETURN

!ELSE IF (ABS(PIO2_D-theta)<10.0D0**(-prec)) THEN

!	!Assegno i valori semplificati per theta=pigreco/2
!	pi2_mdo: DO m=1,nstop
	
!		idx1=m*(m+1)+m
!		mr=REAL(m,dbl)
!		v_pimn(idx1)=mr*fact2(2.0D0*mr-1.0D0)
	
!		pi2_ndo: DO n=m+2,nstop,2
		
!			nr=REAL(n,dbl)
!			idx1=(n-2)*(n-1)+m
!			idx2=n*(n+1)+m
		
!			v_pimn(idx2)=-(nr+mr-1.0D0)*v_pimn(idx1)/(nr-mr)
		
!		END DO pi2_ndo
	
	
!	END DO pi2_mdo

!ELSE	

!Qui c'e' il caso generale
m_do: DO m=1,nstop

mr=REAL(m,dbl)

!A seconda del valore di m,devo calcolare 1 o 2 o piu' valori
m_case: SELECT CASE (nstop-m)

CASE(0) m_case	

	!******************************
	!Caso in cui nstop=m
	!******************************

	one_case_ns: SELECT CASE(m)
				
	CASE(1)			!m=1,m=nstop 
	
		idx1=3
		v_pimn(idx1)=1.0D0

	CASE DEFAULT	!m/=1,m=nstop 
	
		idx1=nstop*(nstop+2)
		fact=fact2(2.0D0*mr-1.0D0)
		v_pimn(idx1)=mr*fact*(SIN(theta)**(m-1))
		
	END SELECT one_case_ns
	
CASE(1) m_case 

	!******************************
	!Caso in cui nstop=m+1
	!******************************

	one_case_ns1: SELECT CASE(m)
				
	CASE(1)			!m=1,m=nstop-1 
	
		!Indici del vettore
		idx1=3
		idx2=7
		!Valori
		v_pimn(idx1)=1.0D0
		v_pimn(idx2)=3.0D0*COS(theta)*v_pimn(idx1)
		
	
	CASE DEFAULT	!m/=1,m=nstop-1 
	
		!Indici del vettore
		idx1=(nstop**2)-1
		idx2=nstop*(nstop+1)+nstop-1
		!Calcolo i valori per (m,m) e (m,m+1)
		fact=fact2(2.0D0*mr-1.0D0)
		v_pimn(idx1)=mr*fact*(SIN(theta)**(m-1))
		v_pimn(idx2)=(2.0D0*mr+1.0D0)*COS(theta)*v_pimn(idx1)
			
	END SELECT one_case_ns1

CASE DEFAULT m_case
	
	!******************************
	!Caso generale in cui nstop>m+1
	!******************************

	!Calcolo gli starting values a seconda che sia m=1 o meno
	start_case: SELECT CASE(m)
				
	CASE(1)			!m=1 
	
		!Indici del vettore
		idx1=3
		idx2=7
		!Valori
		v_pimn(idx1)=1.0D0
		v_pimn(idx2)=3.0D0*COS(theta)*v_pimn(idx1)
		
	
	CASE DEFAULT	!m/=1,m=nstop-1 
	
		!Indici del vettore
		idx1=m*(m+2)
		idx2=(m+1)*(m+2)+m
		!Calcolo i valori per (m,m) e (m,m+1)
		fact=fact2(2.0D0*mr-1.0D0)
		v_pimn(idx1)=mr*fact*(SIN(theta)**(m-1))
		v_pimn(idx2)=(2.0D0*mr+1.0D0)*COS(theta)*v_pimn(idx1)
			
	END SELECT start_case

	!Calcolo i restanti valori per ricorsione
	n_do: DO n=m+2,nstop
		
		!Sistemo indici reali ed interi
		nr=REAL(n,dbl)
		idx1=(n-2)*(n-1)+m
		idx2=(n-1)*n+m
		idx3=n*(n+1)+m
		
		!Faccio il conto
		v_pimn(idx3)=((2.0D0*nr-1.0D0)/(nr-mr))*COS(theta)*v_pimn(idx2)- &
					&((nr+mr-1.0D0)/(nr-mr))*v_pimn(idx1)
	
	END DO n_do 

END SELECT m_case

END DO m_do
	
!END IF theta_if

!Calcolo dei valori per m negativo
nminus_do: DO n=1,nstop
	mminus_do: DO m=1,n

		!Sistemo gli indici
		nr=REAL(n,dbl)
		mr=REAL(m,dbl)
		idx1=n*(n+1)+m
		idx2=n*(n+1)-m
		
		!Faccio i conti
		fact=lnf(nr-mr,error)-lnf(nr+mr,error) 
		v_pimn(idx2)=((-1.0D0)**(m+1))*EXP(fact)*v_pimn(idx1)
		

	END DO mminus_do
END DO nminus_do

END SUBROUTINE pi_mnsca


!******************************************************************************
!3bis) SUBROUTINE tau_mn: calcolo le funzioni angolari Tau_mn
!******************************************************************************
SUBROUTINE tau_mnsca(nstop,theta,v_pimn,v_taumn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_pimn			! Vettore di input valori funzioni pi_mn
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_taumn			! Vettore di output valori funzioni tau_mn

!Dichiarazione variabili interne 
INTEGER(lo), PARAMETER :: prec=7						! Precisione su theta per 0 o pigreco
INTEGER(lo) :: idx1,idx2,idx3,m,n						! Indice per riordinare i miei calcoli e non
REAL(dbl) :: fact,mr,nr,sig								! Semifattoriale e m,n reali,variabile segno


!Subroutine vera e propria
error=0

!Controllo se nstop>0 
nm_check: IF (nstop<0) THEN
			WRITE(*,*) "nstop<0, la procedura pi_mn si ferma"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
!thetain_check: IF ((theta>PI_D) .OR. (theta<0.0D0)) THEN
!				WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
!                error=1
!                RETURN
!END IF thetain_check

!Inizializzo a zero tutto il vettore
v_taumn=0.0D0

!****************************************************************
!Calcolo dei valori per theta=0,Pi
!****************************************************************
!theta_if :IF ( (theta<10.0D0**(-prec)) .OR. ((PI_D-theta)<10.0D0**(-prec)) ) THEN

!	!Assegno il segno per i miei conti
!	sig_if: IF (theta<PIO2_D) THEN
!		sig=1.0D0
!	ELSE
!		sig=-1.0D0
!	END IF sig_if

!	!Assegno i valori semplificati per theta=0
!	theta0_ndo: DO n=1,nstop
	
!		nr=REAL(n,dbl)
	
!		theta0_mdo: DO m=-n,n
		
!			!Calcolo l'indice
!			idx1=n*(n+1)+m
			
!			!Calcolo per m=1 o m=-1
!			theta0_mcase: SELECT CASE (m)
			
!			CASE(-1)
!				!m=-1
!				v_taumn(idx1)=-0.5D0*(sig**(n))
				
!			CASE(1)
!				!m=1
!				v_taumn(idx1)=(sig**(n))*0.5D0*nr*(nr+1.0D0)
				
!			END SELECT theta0_mcase
	
!		END DO theta0_mdo
!	END DO theta0_ndo
	
!	!Chiudo la subroutine
!	RETURN

!END IF theta_if


!****************************************************************
!Calcolo dei valori nel caso generale
!****************************************************************
m_do: DO m=0,nstop

	mr=REAL(m,dbl)

	!A seconda del valore di m,devo calcolare 1 o 2 o piu' valori
	m_case: SELECT CASE (m)
	
	CASE(0) m_case	
	
		!******************************
		!Caso in cui m=0
		!******************************
	
!		m0_thetaif: IF (ABS(PIO2_D-theta)<10.0D0**(-prec)) THEN 
		
!			!******************************
!			!Caso theta=Pi/2
!			!******************************
		
!			!Starting value
!			idx1=2
!			v_taumn(idx1)=-1.0D0
			
!			!Ciclo per il calcolo
!			m0_pi2ndo: DO n=3,nstop,2
				
!				!Indici
!				nr=REAL(n,dbl)
!				idx1=(n-2)*(n-1)
!				idx2=n*(n+1)
				
!				!Calcolo	
!				v_taumn(idx2)=-(nr)*v_taumn(idx1)/(nr-1)
			
!			END DO m0_pi2ndo
		
!		ELSE
		
			!******************************
			!Caso generale theta/=Pi/2
			!******************************
		
			!Starting value
			idx1=2
			v_taumn(idx1)=-SIN(theta)
		
			!If a seconda del valore di nstop
			nstop_if: IF (nstop>1) THEN
			
				!Secondo starting value
				idx1=6
				v_taumn(idx1)=-3.0D0*COS(theta)*SIN(theta)
				
				!Ciclo per il calcolo
				m0_ndo: DO n=3,nstop
				
					!Assegnazione degli indici
					nr=REAL(n,dbl)
					idx1=(n-2)*(n-1)
					idx2=(n-1)*n
					idx3=n*(n+1)
					
					!Calcolo
					v_taumn(idx3)=COS(theta)*(2.0D0*nr-1.0D0)*v_taumn(idx2)/(nr-1.0D0) - &
								 & nr*v_taumn(idx1)/(nr-1.0D0)   
				
				END DO m0_ndo
			
			END IF nstop_if
		
!		END IF m0_thetaif
	
	CASE DEFAULT 
	
		!******************************
		!Caso in cui m/=0
		!******************************
	
!		m_thetaif: IF (ABS(PIO2_D-theta)<10.0D0**(-prec)) THEN
		
!			!******************************
!			!Caso theta=Pi/2
!			!******************************
		
!			!Ciclo per l'assegnazione dei valori
!			pi2_ndo: DO n=m+1,nstop,2
			
!				!Indici
!				nr=REAL(n,dbl)
!				idx1=(n-1)*n+m
!				idx2=n*(n+1)+m
				
!				!Calcolo
!				v_taumn(idx2)=-(nr+mr)*v_pimn(idx1)/mr
			
!			END DO pi2_ndo
			
!		ELSE
		
			!******************************
			!Caso generale theta/=Pi/2
			!******************************
		
			!Starting value
			idx1=m*(m+1)+m
			v_taumn(idx1)=COS(theta)*v_pimn(idx1)
			
			!Ciclo per l'assegnazione dei valori
			ndo: DO n=m+1,nstop
			
				!Indici
				nr=REAL(n,dbl)
				idx1=(n-1)*n+m
				idx2=n*(n+1)+m
				
				!Calcolo
				v_taumn(idx2)=COS(theta)*v_pimn(idx2)*(nr/mr)-(nr+mr)*v_pimn(idx1)/mr
			
			END DO ndo
					
!		END IF m_thetaif
		
	END SELECT m_case
	
END DO m_do

!Calcolo dei valori per m negativo
nminus_do: DO n=1,nstop
	mminus_do: DO m=1,n

		!Sistemo gli indici
		nr=REAL(n,dbl)
		mr=REAL(m,dbl)
		idx1=n*(n+1)+m
		idx2=n*(n+1)-m
		
		!Faccio i conti
		fact=lnf(nr-mr,error)-lnf(nr+mr,error) 
		v_taumn(idx2)=((-1.0D0)**(m))*EXP(fact)*v_taumn(idx1)
		
	END DO mminus_do
END DO nminus_do


END SUBROUTINE tau_mnsca


!******************************************************************************
!4) SUBROUTINE leg_mn: calcolo le funzioni angolari leg_mn
!******************************************************************************
SUBROUTINE leg_mn(nstop,theta,v_legmn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_legmn			! Vettore di output valori funzioni leg_mn

!Dichiarazione variabili interne 
INTEGER(lo), PARAMETER :: prec=7						! Precisione su theta per 0 o pigreco
INTEGER(lo) :: mm,m,n,beg,indx						! Indice per riordinare i miei calcoli e non
REAL(dbl) :: fact,mr,nr,sig								! Semifattoriale e m,n reali,variabile segno
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_leg			! Vettore temporaneo legendre

!Subroutine vera e propria
error=0

!Controllo se nstop>0 
nm_check: IF (nstop<0.0D0) THEN
			WRITE(*,*) "nstop<0, la procedura pi_mn si ferma"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((theta>PI_D) .OR. (theta<0.0D0)) THEN
				WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
                error=1
                RETURN
END IF thetain_check

!Inizializzo a zero tutto il vettore
v_legmn=0.0D0

!Riempio il vettore v_legtot
leg_mdo: DO m=-nstop,nstop

	mm=ABS(m)
	mr=REAL(m,dbl)

	ALLOCATE(v_leg(mm:nstop))
	CALL legendre(mm,nstop,m,mm,theta,v_leg,error)

	!If per l'errore su legendre
	errleg_if: IF (error/=0) THEN
		WRITE(*,*) "Errore in leg_mn, la procedura si arresta"
		RETURN
	END IF errleg_if

	!If per lo starting index
	begin_if: IF (mm==0) THEN
		beg=mm+1
	ELSE
		beg=mm
	END IF begin_if
	
	!Riempio il vettore finale
	leg_ndo: DO n=beg,nstop
	
		nr=REAL(n,dbl)		
		indx=n*(n+1)+m
		v_legmn(indx)=v_leg(n)
	
	END DO leg_ndo
		
	DEALLOCATE(v_leg)

END DO leg_mdo

END SUBROUTINE leg_mn


!******************************************************************************
!5) SUBROUTINE d_nmk: calcolo gli elementi ridotti di matrice secondo edmonds
!						 d(beta)^{nstop}_{m,k}	
!******************************************************************************
SUBROUTINE d_nmk(nstop,m,k,beta,v_d,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,m,k					! Numero di espansioni multipolari, d^{nstop}_{m,k}(beta)
REAL(dbl), INTENT(IN) :: beta							! Angolo Beta di Eulero
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_d		! Vettore di output valori funzioni pi_mn

!Dichiarazione variabili interne 
INTEGER(lo) :: ierr									! altra flag di errore
INTEGER(lo) :: nmin,modm,modk,MpK,MmK,emk,n			! Indice per i miei calcoli
REAL(dbl) :: logw,w,mr,kr,nr,x							! Fattoriale e indici reali, x=cos(beta)
REAL(dbl) :: c0,c1,c2									!Coefficienti ricorrenza


!Subroutine vera e propria
error=0

!Calcolo nmin
modm=ABS(m)
modk=ABS(k)
nmin=MAX(modm,modk)

!Controllo se nstop>nmin 
nm_check: IF (nstop<nmin) THEN
			WRITE(*,*) "nstop<nmin, la procedura d_nmk si ferma!"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((beta>PI_D) .OR. (beta<-PI_D)) THEN
				WRITE(*,*) "Beta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
                error=1
                RETURN
END IF thetain_check

!Inizializzo il vettore e alcune grandezze che fanno comodo
v_d=0.0D0		!Vettore

emk_if: IF (k>=m) THEN !Emk
	emk=1.0D0
ELSE
	emk=(-1.0D0)**(m-k)
END IF emk_if

MpK=ABS(m+k)	!moduli
MmK=ABS(m-k)

x=COS(beta)		!Coseno di beta

mr=REAL(m,dbl)	!Diventato reali
kr=REAL(k,dbl)

!If per vedere se m=k=0 come caso particolare
!In questo caso la ricorsione e' quella dei polinomi di legendre semplici
nmin_zero_if: IF (nmin==0) THEN

		!Primo valore
		v_d(nmin)=1

		!Se nstop>0 assegno altri valori
		nstop_zero_if: IF (nstop>0) THEN

			!Secondo valore
			v_d(nmin+1)=x

			!Ciclo su n
			nstop_zero_do: DO n=2,nstop
				nr=REAL(n,dbl)
				v_d(n)=((2.0D0*nr-1.0D0)*x*v_d(n-1)-(nr-1.0D0)*v_d(n-2))/nr
			END DO  nstop_zero_do

		END IF nstop_zero_if

ELSE

	!Calcolo il primo valore della ricorrenza
	logw=lnf(2.0D0*REAL(nmin,dbl),ierr) - lnf(REAL(MmK,dbl),ierr) -lnf(REAL(MpK,dbl),ierr)
	w=EXP(logw)
	v_d(nmin)=((-1.0D0)**(m+k))*emk*(2.0D0**(-nmin))*SQRT(w)*((SQRT(1-x))**MmK)*((SQRT(1+x))**MpK)

	!Ciclo per il calcolo
	nstop_do: DO n=nmin+1,nstop
	
				!Coefficienti della ricorrenza
				nr=REAL(n,dbl)
				c0=1.0D0/((nr-1.0D0)*SQRT(nr**2-mr**2)*SQRT(nr**2-kr**2))
				c1=(2.0D0*nr-1.0D0)*(nr*(nr-1.0D0)*x-mr*kr)
				c2=nr*SQRT((nr-1.0D0)**2-mr**2)*SQRT((nr-1.0D0)**2-kr**2)
				
				!Calcolo
				v_d(n)=(c1*v_d(n-1)-c2*v_d(n-2))*c0
				
	END DO nstop_do

END IF nmin_zero_if

END SUBROUTINE d_nmk



!******************************************************************************
!5) SUBROUTINE d_n_km: calcolo gli elementi ridotti di matrice secondo edmonds
!                        d(beta)^{nstop}_{k,m}  
!******************************************************************************
SUBROUTINE d_n_km(nstop,k,m,beta,v_d,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,k,m                  ! Numero di espansioni multipolari, d^{nstop}_{m,k}(beta)
REAL(dbl), INTENT(IN) :: beta                           ! Angolo Beta di Eulero
INTEGER(lo), INTENT(OUT) :: error                     ! Flag di errore
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_d       ! Vettore di output valori funzioni pi_mn

!Dichiarazione variabili interne 
INTEGER(lo) :: ierr,flag                              ! altra flag di errore
INTEGER(lo) :: nmin,modm,modk,n,mi,ki     ! Indice per i miei calcoli
REAL(dbl) :: logw,w,mr,kr,nr,x                          ! Fattoriale e indici reali, x=cos(beta)
REAL(dbl) :: c0,c1,c2,eps,eta                           !Coefficienti ricorrenza


!Subroutine vera e propria
error=0

!Calcolo nmin
modm=ABS(m)
modk=ABS(k)
nmin=MAX(modm,modk)

!Controllo se nstop>nmin 
nm_check: IF (nstop<nmin) THEN
            WRITE(*,*) "nstop<nmin, la procedura d_n_km si ferma!"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((beta>PI_D) .OR. (beta<-PI_D)) THEN
                WRITE(*,*) "Beta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
                error=1
                RETURN
END IF thetain_check

!Inizializzo il vettore e alcune grandezze che fanno comodo
v_d=0.0D0       !Vettore


x=COS(beta)     !Coseno di beta
eps=COS(beta/2.0D0)
eta=SIN(beta/2.0D0)

mr=REAL(m,dbl)  !Diventato reali
kr=REAL(k,dbl)

!If per vedere se m=k=0 come caso particolare
nmin_zero_if: IF (nmin==0) THEN

        !Primo valore
        v_d(nmin)=1.0D0
        
        !Se nstop>0 assegno altri valori
        nstop_zero_if: IF (nstop>0) THEN
            
            !Secondo valore
            v_d(nmin+1)=x
            
            !Ciclo su n
            nstop_zero_do: DO n=2,nstop
                nr=REAL(n,dbl)
                v_d(n)=((2.0D0*nr-1.0D0)*x*v_d(n-1)-(nr-1.0D0)*v_d(n-2))/nr
            END DO  nstop_zero_do

        END IF nstop_zero_if
RETURN

END IF nmin_zero_if


mag_if: IF (ABS(m)>=ABS(k)) THEN

    ki=k
    mi=m

    sig_if1: IF (mi>=0) THEN

        flag=1

    ELSE

        flag=2
        ki=-ki
        mi=-mi

    END IF sig_if1

ELSE

    ki=m
    mi=k

    sig_if2: IF (mi>=0) THEN

        flag=3

    ELSE

        flag=4
        ki=-ki
        mi=-mi

    END IF sig_if2

END IF mag_if


!Calcolo il primo valore della ricorrenza
logw=lnf(2.0D0*REAL(mi,dbl),ierr) - lnf(REAL(mi+ki,dbl),ierr) -lnf(REAL(mi-ki,dbl),ierr)
w=EXP(logw)
v_d(nmin)=((-1.0D0)**(mi-ki))*SQRT(w)*(eps**(mi+ki))*(eta**(mi-ki))

!Ciclo per il calcolo
nstop_do: DO n=nmin+1,nstop
    
            !Coefficienti della ricorrenza
            nr=REAL(n,dbl)
            c0=1.0D0/((nr-1.0D0)*SQRT(nr**2-kr**2)*SQRT(nr**2-mr**2))
            c1=(2.0D0*nr-1.0D0)*(nr*(nr-1.0D0)*x-kr*mr)
            c2=nr*SQRT((nr-1.0D0)**2-kr**2)*SQRT((nr-1.0D0)**2-mr**2)
            
            !Calcolo
            v_d(n)=(c1*v_d(n-1)-c2*v_d(n-2))*c0
                
END DO nstop_do

!Moltiplico per il segno
final_if: IF ((flag==2).OR.(flag==3)) THEN
    v_d=((-1.0D0)**(k-m))*v_d
END IF final_if

END SUBROUTINE d_n_km




!******************************************************************************
!6) SUBROUTINE d_nm12: calcolo gli elementi ridotti di matrice per k=1 e k=-1, in
! maniera da poter calcolare la sezione di estinzione
!******************************************************************************
SUBROUTINE d_nm12(nstop,beta,v_dnm1,v_dnm2,error)


IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: beta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_dnm1,v_dnm2	! Vettore di output valori funzioni d_mn1

!Dichiarazione variabili interne 
INTEGER(lo) :: mm,m,n,beg,indx						! Indice per riordinare i miei calcoli e non
REAL(dbl) :: mr,nr										! Semifattoriale e m,n reali,variabile segno
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_d1,v_d2		! Vettore temporaneo red rot mat el

!Subroutine vera e propria
error=0

!Controllo se nstop>0 
nm_check: IF (nstop<0) THEN
			WRITE(*,*) "nstop<0, la procedura d_nm2 si ferma"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
betain_check: IF ((beta>PI_D) .OR. (beta<-PI_D)) THEN
				WRITE(*,*) "Beta e' al di fuori dei suoi bounds, la procedura d_nm2 si ferma"
                error=1
                RETURN
END IF betain_check

!Inizializzo a zero tutto il vettore
v_dnm1=0.0D0
v_dnm2=0.0D0

!Riempio il vettore v_legtot
d_mdo: DO m=-nstop,nstop

	mm=ABS(m)
	mr=REAL(m,dbl)

	ALLOCATE(v_d1(0:nstop),v_d2(0:nstop))
	CALL d_nmk(nstop,m,1,beta,v_d1,error)
	CALL d_nmk(nstop,m,-1,beta,v_d2,error)

	!If per l'errore su legendre
	errd_if: IF (error/=0) THEN
		WRITE(*,*) "Errore in d_nm12, la procedura si arresta"
		RETURN
	END IF errd_if

	!If per lo starting index
	begin_if: IF (mm==0) THEN
		beg=mm+1
	ELSE
		beg=mm
	END IF begin_if
	
	!Riempio il vettore finale
	d_ndo: DO n=beg,nstop
	
		nr=REAL(n,dbl)		
		indx=n*(n+1)+m
		v_dnm1(indx)=v_d1(n)
		v_dnm2(indx)=v_d2(n)
	
	END DO d_ndo
		
	DEALLOCATE(v_d1,v_d2)

END DO d_mdo

END SUBROUTINE d_nm12






!******************************************************************************
!5tris) SUBROUTINE field_expRandom_sub: calcolo i coefficienti del campo incidente
! per una direzione arbitraria rispetto al sistema di riferimento del cluster che e'
! ruotato degli angoli alpha bete gamma.
!******************************************************************************
SUBROUTINE field_expRandom_sub(nstop,ns,k,betap,m_xyz,alpha,beta,gamma,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,betap							!Vettore d'onda e polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Matrice coordinate
REAL(dbl), INTENT(IN) :: alpha,beta,gamma					!Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p				!Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj,n,m,dimens								!Indici e dimensioni
REAL(dbl) :: nr,mr,x,y,z,d,theta,phi							!Indici reali e coordinate
REAL(dbl) :: kd													!Phase 
COMPLEX(dbl) :: pshift,expg1,expg2								!Phase shift				
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2			!Elementi ridotti di matrice di rotazione


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))

!Funzione angolare Pi_mn
CALL d_nm12(nstop,beta,v_dnm1,v_dnm2,error)

dnm_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da field_expRandom_sub: il programma termina ora...") 
       			STOP
END IF dnm_if

!Calcolo gli esponenziali
expg1=EXP(gamma*(0.0D0,1.0D0))
expg2=EXP(-gamma*(0.0D0,1.0D0))

!Inizializzo j e i vettori output
j=0
jj=0
v_p=(0.0D0,0.0D0)

i_do: DO i=1,ns

	jj=0

	!Assegno x y e z
	x=m_xyz(i,1)
	y=m_xyz(i,2)
	z=m_xyz(i,3)
	
	!Calcolo phase e phase shift
	kd=k*(SIN(beta)*(x*COS(alpha)+y*SIN(alpha)) + z*COS(beta))
	pshift=EXP(CMPLX(0.0D0,kd))

	n_do: DO n=1,nstop
	
		!n reale
		nr=REAL(n,dbl)
	
		m_do: DO m=-n,n
		
			!m reale
			mr=REAL(m,dbl)	
		
			!Incremento j
			j=j+1
			jj=jj+1
		
			v_p(2*j-1)=0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
			v_p(2*j)=  0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )
				
		END DO m_do
	END DO n_do
END DO i_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE field_expRandom_sub



!******************************************************************************
!5tris) SUBROUTINE field_expRandom_sub: calcolo i coefficienti del campo incidente
! per una direzione arbitraria rispetto al sistema di riferimento del cluster che e'
! ruotato degli angoli alpha bete gamma.
!******************************************************************************
SUBROUTINE field_expRandom_shell_sub(nstop,ns,k,betap,v_xyz,alpha,beta,gamma,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns					!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,betap						!Vettore d'onda e polarizzazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_xyz				!Matrice coordinate
REAL(dbl), INTENT(IN) :: alpha,beta,gamma					!Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p				!Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj,n,m,dimens						!Indici e dimensioni
REAL(dbl) :: nr,mr,x,y,z,d,theta,phi					!Indici reali e coordinate
REAL(dbl) :: kd									!Phase 
COMPLEX(dbl) :: pshift,expg1,expg2						!Phase shift				
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2			!Elementi ridotti di matrice di rotazione


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))

!Funzione angolare Pi_mn
CALL d_nm12(nstop,beta,v_dnm1,v_dnm2,error)

dnm_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da field_expRandom_sub: il programma termina ora...") 
       			STOP
END IF dnm_if

!Calcolo gli esponenziali
expg1=EXP(gamma*(0.0D0,1.0D0))
expg2=EXP(-gamma*(0.0D0,1.0D0))

!Inizializzo j e i vettori output
j=0
jj=0
v_p=(0.0D0,0.0D0)



jj=0

!Assegno x y e z
x=v_xyz(1)
y=v_xyz(2)
z=v_xyz(3)

!WRITE(*,*) "alpha,beta,gamma",alpha,beta,gamma

!Calcolo phase e phase shift
kd=k*(SIN(beta)*(x*COS(alpha)+y*SIN(alpha)) + z*COS(beta))
pshift=EXP(CMPLX(0.0D0,kd))

n_do: DO n=1,nstop

	!n reale
	nr=REAL(n,dbl)

	m_do: DO m=-n,n
	
		!m reale
		mr=REAL(m,dbl)
	
		!Incremento j
		j=j+1
		jj=jj+1
	
		v_p(2*j-1)=0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
				& ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
		v_p(2*j)=  0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
				& ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )
				
!		WRITE(*,*) "v_p_int",v_p(2*j-1)
!		WRITE(*,*) "v_q_int",v_p(2*j)
			
	END DO m_do
END DO n_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE field_expRandom_shell_sub



!******************************************************************************
!5tris-bis) SUBROUTINE field_expRandom_dip_sub: calcolo i coefficienti del campo incidente,ma li randomizzo per la 
! fase e per l'orientamento,e questo mi serve per un tipo di calcoli.
!******************************************************************************
SUBROUTINE field_expRandom_dip_sub(nstop,ns,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns				!N espansioni multipolari e N sfere
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p			!Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error				!Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj,n,m,dimens					!Indici e dimensioni
REAL(dbl) :: alpha,beta,gamma						!Angoli di eulero 
REAL(dbl) :: nr,mr,x,y,z,d,theta,phi				!Indici reali e coordinate
REAL(dbl) :: kd								!Phase 
COMPLEX(dbl) :: pshift,expg1,expg2					!Phase shift
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2		!Elementi ridotti di matrice di rotazione


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))


!Inizializzo j e i vettori output
j=0
jj=0
v_p=(0.0D0,0.0D0)

i_do: DO i=1,ns

	!Calcolo gli angoli randomizzati
	CALL RANDOM_NUMBER(alpha)
	alpha=alpha*2*Pi_D

	CALL RANDOM_NUMBER(beta)
	beta=beta*Pi_D

	CALL RANDOM_NUMBER(gamma)
	gamma=gamma*2*Pi_D

	!Funzione angolare Pi_mn
	CALL d_nm12(nstop,beta,v_dnm1,v_dnm2,error)

	dnm_if: IF (error/=0) THEN																	
	    			WRITE(*,10) 
	    			10 FORMAT ("Errore in d_nm12 chiamata da field_expRandom_dip_sub: il programma termina ora...") 
	       			STOP
	END IF dnm_if

	!Calcolo gli esponenziali
	expg1=EXP(gamma*(0.0D0,1.0D0))
	expg2=EXP(-gamma*(0.0D0,1.0D0))

	jj=0
	
	!Calcolo phase e phase shift
	CALL RANDOM_NUMBER(kd)
	kd=kd*2*Pi_D
	pshift=EXP(CMPLX(0.0D0,kd))

	n_do: DO n=1,nstop
	
		!n reale
		nr=REAL(n,dbl)
	
		m_do: DO m=-n,n
		
			!m reale
			mr=REAL(m,dbl)	
		
			!Incremento j
			j=j+1
			jj=jj+1
		
			v_p(2*j-1)=0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
			v_p(2*j)=  0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )
				
		END DO m_do
	END DO n_do
END DO i_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE field_expRandom_dip_sub



!******************************************************************************
!5tris-tris) SUBROUTINE field_expPhased_dip_sub: calcolo i coefficienti del campo incidente
! per una direzione arbitraria rispetto al sistema di riferimento del cluster che e'
! ruotato degli angoli alpha bete gamma, qui però suppongo che tutti i dipoli siano in fase,
! quindi metto a zero fittiziamente tutte le coordinate.
!******************************************************************************
SUBROUTINE field_expPhased_dip_sub(nstop,ns,k,betap,alpha,beta,gamma,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns				!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,betap					!Vettore d'onda e polarizzazione
REAL(dbl), INTENT(IN) :: alpha,beta,gamma				!Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p			!Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error				!Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj,n,m,dimens					!Indici e dimensioni
REAL(dbl) :: nr,mr,x,y,z,d,theta,phi				!Indici reali e coordinate
REAL(dbl) :: kd								!Phase 
COMPLEX(dbl) :: pshift,expg1,expg2					!Phase shift				
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2		!Elementi ridotti di matrice di rotazione


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))

!Funzione angolare Pi_mn
CALL d_nm12(nstop,beta,v_dnm1,v_dnm2,error)

dnm_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da field_expRandom_sub: il programma termina ora...") 
       			STOP
END IF dnm_if

!Calcolo gli esponenziali
expg1=EXP(gamma*(0.0D0,1.0D0))
expg2=EXP(-gamma*(0.0D0,1.0D0))

!Inizializzo j e i vettori output
j=0
jj=0
v_p=(0.0D0,0.0D0)

i_do: DO i=1,ns

	jj=0

	!Assegno x y e z
	x=0.0D0
	y=0.0D0
	z=0.0D0
	
	!Calcolo phase e phase shift
	kd=k*(SIN(beta)*(x*COS(alpha)+y*SIN(alpha)) + z*COS(beta))
	pshift=EXP(CMPLX(0.0D0,kd))

	n_do: DO n=1,nstop
	
		!n reale
		nr=REAL(n,dbl)
	
		m_do: DO m=-n,n
		
			!m reale
			mr=REAL(m,dbl)	
		
			!Incremento j
			j=j+1
			jj=jj+1
		
			v_p(2*j-1)=0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
			v_p(2*j)=  0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )
				
		END DO m_do
	END DO n_do
END DO i_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE field_expPhased_dip_sub




!******************************************************************************
!5quater) SUBROUTINE fnmp_sub: calcolo parte dei coefficienti di espansione del
!campo incidente, la parte indipendente dalla particella e dalla lunghezza d'onda
!per risparmiare tempo visto che la faccenda comincia a diventare complicata
!******************************************************************************
SUBROUTINE fmnp_sub(nstop,alpha,beta,v_fnmp1,v_fnmp2,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop                           !N espansioni multipolari
REAL(dbl), INTENT(IN) :: alpha,beta                         !Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_fnmp1,v_fnmp2  !Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error                        !Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: j,n,m,dimens                              !Indici e dimensioni
REAL(dbl) :: nr,mr                                         !Indici reali e coordinate
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2      !Elementi ridotti di matrice di rotazione


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))

!Funzione angolare Pi_mn
CALL d_nm12(nstop,beta,v_dnm1,v_dnm2,error)

dnm_if: IF (error/=0) THEN                                                                  
                WRITE(*,10) 
                10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da field_expRandom_sub: il programma termina ora...") 
                STOP
END IF dnm_if

!Inizializzo j e i vettori output
j=0
v_fnmp1=(0.0D0,0.0D0)
v_fnmp2=(0.0D0,0.0D0)

    n_do: DO n=1,nstop
    
        !n reale
        nr=REAL(n,dbl)
    
        m_do: DO m=-n,n
        
            !m reale
            mr=REAL(m,dbl)  
        
            !Incremento j
            j=j+1
        
            v_fnmp1(2*j-1)=-0.5D0*SQRT(2.0D0*nr+1.0D0)*EXP(-(0.0D0,1.0D0)*mr*alpha)*v_dnm2(j)
            v_fnmp1(2*j)=  -v_fnmp1(2*j-1)
                
            v_fnmp2(2*j-1)=0.5D0*SQRT(2.0D0*nr+1.0D0)*EXP(-(0.0D0,1.0D0)*mr*alpha)*v_dnm1(j)
            v_fnmp2(2*j)=  v_fnmp2(2*j-1)

        END DO m_do
    END DO n_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE fmnp_sub



!******************************************************************************
!5quinties) SUBROUTINE field_expfnmp_sub: calcolo i coefficienti del campo incidente
! non polarizzati come moltiplicazione del ritardo di fase per fnmp
!******************************************************************************
SUBROUTINE field_expfnmp_sub(nstop,ns,k,m_xyz,alpha,beta,v_fnmp,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns                        !N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k                                  !Vettore d'onda e polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz              !Matrice coordinate
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_fnmp            !Vettore espansione campo parziale
REAL(dbl), INTENT(IN) :: alpha,beta                         !Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p              !Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error                        !Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,low_b,up_b,dimens2               !Indici, dimensioni e bounds
REAL(dbl) :: x,y,z                                  !Indici reali e coordinate
REAL(dbl) :: kd                                     !Phase 
COMPLEX(dbl) :: pshift                              !Phase shift                


!Inizio subroutine vera e propria
error=0
dimens2=2*nstop*(nstop+2)

!Inizializzo j e i vettori output
v_p=(0.0D0,0.0D0)

i_do: DO i=1,ns

    !Assegno x y e z
    x=m_xyz(i,1)
    y=m_xyz(i,2)
    z=m_xyz(i,3)
    
    !Calcolo phase e phase shift
    kd=k*(SIN(beta)*(x*COS(alpha)+y*SIN(alpha)) + z*COS(beta))
    pshift=EXP(CMPLX(0.0D0,kd))

    !Calcolo bounds
    low_b=1+(i-1)*dimens2
    up_b=i*dimens2

    !Riempio v_p
    v_p(low_b:up_b)=pshift*v_fnmp    

END DO i_do

END SUBROUTINE field_expfnmp_sub



!******************************************************************************
!6) SUBROUTINE cext_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cext_sub(nstop,ns,k,beta,v_z,v_ab,norm,v_cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_z					!Vettore coordinata z
REAL(dbl), INTENT(IN) :: k,beta								!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab				!Vettori coefficienti espansioni campo, non corretti e corretti  
CHARACTER(len=3), INTENT(IN) :: norm						!Flag normalizzazione

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cext				!Vettore sezioni di estinzione
! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,l_ib,up_ib,l_ia,up_ia,prec,elem
REAL(dbl) :: nr
COMPLEX(dbl) :: somma


!Funzione vera e propria
elem=2*nstop*(nstop+2)

norm_if: IF (norm=='old') THEN

		i_old_do: DO i=1,ns
		
			!Inizializzo la variabile somma
			somma=(0.0D0,0.0D0)
			
			!Mi porto sulla giusta posizione nel vettore v_ab, per la giusta sfera
			prec=(i-1)*elem
			
			n_old_do: DO n=1,nstop
			
				!Indice reale n
				nr=REAL(n,dbl)
				
				!Costruisco gli indici per recuperare i coefficienti
				l_ib=prec+2*(n*(n+1)-1)
				up_ib=prec+2*(n*(n+1)+1)
				l_ia=l_ib-1
				up_ia=up_ib-1
			
				!Costruisco la somma
				somma=somma+(2.0D0*nr+1.0D0)*((v_ab(up_ia)+v_ab(up_ib))*EXP((0.0D0,1.0D0)*beta) - &
											& nr*(nr+1.0D0)*(v_ab(l_ia)-v_ab(l_ib))*EXP(-(0.0D0,1.0D0)*beta))
			
			END DO n_old_do
			
			!Calcolo la funzione
			v_cext(i)=(2*Pi_d/(k**2))*REAL(EXP(-(0.0D0,1.0D0)*k*v_z(i))*somma,dbl)
			
		END DO i_old_do
		
ELSE
		
		i_new_do: DO i=1,ns
		
			!Inizializzo la variabile somma
			somma=(0.0D0,0.0D0)
			
			!Mi porto sulla giusta posizione nel vettore v_ab, per la giusta sfera
			prec=(i-1)*elem
			
			n_new_do: DO n=1,nstop
			
				!Indice reale n
				nr=REAL(n,dbl)
				
				!Costruisco gli indici per recuperare i coefficienti
				l_ib=prec+2*(n*(n+1)-1)
				up_ib=prec+2*(n*(n+1)+1)
				l_ia=l_ib-1
				up_ia=up_ib-1
			
				!Costruisco la somma
				somma=somma+SQRT(2.0D0*nr+1.0D0)*( (v_ab(up_ia)+v_ab(up_ib))*EXP((0.0D0,1.0D0)*beta) - &
												&  (v_ab(l_ia)-v_ab(l_ib))*EXP(-(0.0D0,1.0D0)*beta))
			
			END DO n_new_do
			
			!Calcolo la funzione
			v_cext(i)=(2*Pi_d/(k**2))*REAL(EXP(-(0.0D0,1.0D0)*k*v_z(i))*somma,dbl)
			
		END DO i_new_do

END IF norm_if

END SUBROUTINE cext_sub



!******************************************************************************
!6bis) SUBROUTINE cext_random_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cext_random_sub(nstop,ns,k,v_p,v_ab,v_cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k									!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p,v_ab			!Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cext				!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j
COMPLEX(dbl) :: somma

!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma+CONJG(v_p(2*j-1))*v_ab(2*j-1)+CONJG(v_p(2*j))*v_ab(2*j)
		
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_cext(i)=(4*Pi_d/(k**2))*REAL(somma,dbl)
	
END DO sphere_do

END SUBROUTINE cext_random_sub



!******************************************************************************
!6tris) SUBROUTINE csca_random_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE csca_random_sub(nstop,ns,k,v_ab_sca,v_ab,v_csca)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k									!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab_sca,v_ab		!Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_csca				!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j
COMPLEX(dbl) :: somma

!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma+CONJG(v_ab(2*j-1))*v_ab_sca(2*j-1)+CONJG(v_ab(2*j))*v_ab_sca(2*j)
		
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_csca(i)=(4*Pi_d/(k**2))*REAL(somma,dbl)
	
END DO sphere_do

END SUBROUTINE csca_random_sub



!******************************************************************************
!6quinties) SUBROUTINE cabs_random_sub: calcolo la sezione di 
!assorbimento per ciascuna sfera.
!******************************************************************************
SUBROUTINE cabs_random_sub(ref_index,nstop,ns,k,v_dc,v_req,m_epseq,v_cabs)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns				!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,ref_index				!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc 			!Vettori coefficienti espansioni campo, non corretti e corretti  
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq		!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			!Raggi equivalenti

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs			!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j,error
COMPLEX(dbl) :: somma,z,mc
REAL(dbl) :: radius,xpar
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1


!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Calcolo gli argomenti per le funzioni di riccati bessel
    radius=v_req(v_patt(i))
    mc=SQRT(CMPLX(m_epseq(v_patt(i),1),m_epseq(v_patt(i),2),dbl))
    mc=mc/ref_index
    xpar=k*radius
    z=mc*xpar

!	WRITE(*,*) CMPLX(m_epseq(v_patt(i),1),m_epseq(v_patt(i),2),dbl)
!	WRITE(*,*) 'k',k
!	WRITE(*,*) 'radius',radius
!	WRITE(*,*) 'ref_index',ref_index
!	WRITE(*,*) 'mc',mc
!	WRITE(*,*) 'xpar',xpar
!	WRITE(*,*) 'z',z

	!Calcolo la funzione di Riccati Bessel
    CALL psi_z_sub(nstop,z,v_psi,error)

	!Calcolo la derivata della funzione di riccati bessel
	der_do: DO n=1,nstop 
		v_psi1(n)=v_psi(n-1)-REAL(n,dbl)*v_psi(n)/z
	END DO der_do

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma + v_psi(n)*CONJG(v_psi1(n))*( ( ABS(v_dc(2*j) ) )**2 ) - &
            &             v_psi1(n)*CONJG(v_psi(n))*( ( ABS(v_dc(2*j-1) ) )**2 )
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_cabs(i)=-(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma/mc,dbl)
	
END DO sphere_do

END SUBROUTINE cabs_random_sub




!******************************************************************************
!6bis-bis-bis) SUBROUTINE cabs_random_sub: calcolo la sezione di 
!assorbimento per ciascuna sfera.
!******************************************************************************
SUBROUTINE cabs_random_sub_dip(ref_index,nstop,ns_ant,k,v_dc,v_req,m_epseq,v_cabs)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns_ant					!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,ref_index					!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc 	!Vettori coefficienti espansioni campo, non corretti e corretti  
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq		!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			!Raggi equivalenti

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs			!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j,error
COMPLEX(dbl) :: somma,z,mc
REAL(dbl) :: radius,xpar
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1


!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=2,ns_ant+1

	!Calcolo gli argomenti per le funzioni di riccati bessel
	radius=v_req(v_patt(i))
	mc=SQRT(CMPLX(m_epseq(v_patt(i),1),m_epseq(v_patt(i),2)))
	mc=mc/ref_index
	xpar=k*radius
	z=mc*xpar

	!Calcolo la funzione di Riccati Bessel
	CALL psi_z_sub(nstop,z,v_psi,error)

	!Calcolo la derivata della funzione di riccati bessel
	der_do: DO n=1,nstop 
		v_psi1(n)=v_psi(n-1)-REAL(n,dbl)*v_psi(n)/z
	END DO der_do

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma + v_psi(n)*CONJG(v_psi1(n))*( ( ABS(v_dc(2*j) ) )**2 ) - &
            &             v_psi1(n)*CONJG(v_psi(n))*( ( ABS(v_dc(2*j-1) ) )**2 )
		
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_cabs(i-1)=-(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma/mc,dbl)
	
END DO sphere_do

END SUBROUTINE cabs_random_sub_dip





!******************************************************************************
!6bis-bis-bis) SUBROUTINE cabs_random_sub: calcolo la sezione di 
!assorbimento per ciascuna sfera.
!******************************************************************************
SUBROUTINE cabs_random_sub1(ref_index,nstop,ns,k,v_ab,v_req,m_epseq,v_cabs)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns					!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,ref_index					!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab 			!Vettori coefficienti espansioni campo, non corretti e corretti  
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq		!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			!Raggi equivalenti

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs			!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j,error
COMPLEX(dbl) :: somma,z,mc
REAL(dbl) :: radius,xpar
REAL(dbl), DIMENSION(0:nstop) :: v_psi
REAL(dbl), DIMENSION(1:nstop) :: v_psi1,nomD,nomC
COMPLEX(dbl), DIMENSION(0:nstop) :: v_cpsi,v_D,v_C
COMPLEX(dbl), DIMENSION(1:nstop) :: v_cpsi1,denD,denC


!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Calcolo gli argomenti per le funzioni di riccati bessel
    radius=v_req(v_patt(i))
    mc=SQRT(CMPLX(m_epseq(v_patt(i),1),m_epseq(v_patt(i),2)))
    mc=mc/ref_index
    xpar=k*radius
    z=mc*CMPLX(xpar,0.0D0)
    
	!Calcolo la funzione di Riccati Bessel complessa
    CALL psi_z_sub(nstop,z,v_cpsi,error)

	!Calcolo la derivata della funzione di riccati bessel complessa
    cder_do: DO n=1,nstop 
		v_cpsi1(n)=v_cpsi(n-1)-REAL(n,dbl)*v_cpsi(n)/z
	END DO cder_do
    
	!Calcolo la funzione di Riccati Bessel reale
    CALL psi_d_sub(nstop,xpar,v_psi,error)

	!Calcolo la derivata della funzione di riccati bessel reale
    der_do: DO n=1,nstop 
		v_psi1(n)=v_psi(n-1)-REAL(n,dbl)*v_psi(n)/xpar
	END DO der_do
    
	
	!Calcolo i numeratori	
	nomden_DO: DO n=1,nstop
    
		nomD(n)=REAL( (0.0D0,1.0D0)*mc*v_cpsi(n)*CONJG(v_cpsi1(n)) ,dbl)
        nomC(n)=REAL( (0.0D0,1.0D0)*CONJG(mc)*v_cpsi(n)*CONJG(v_cpsi1(n)) ,dbl)
        denD(n)=ABS( mc*v_cpsi(n)*v_psi1(n) - v_psi(n)*v_cpsi1(n) )*ABS( mc*v_cpsi(n)*v_psi1(n) - v_psi(n)*v_cpsi1(n) )
        denC(n)=ABS( v_cpsi(n)*v_psi1(n) - mc*v_psi(n)*v_cpsi1(n) )*ABS( v_cpsi(n)*v_psi1(n) - mc*v_psi(n)*v_cpsi1(n) )

	END DO nomden_DO

	v_D=nomD/denD
    v_C=nomC/denC
	

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma + v_D(n)*(ABS(v_ab(2*j-1))**2) + v_C(n)*(ABS(v_ab(2*j))**2)
		
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_cabs(i)=(4*Pi_d/(k**2))*REAL(somma,dbl)
	
END DO sphere_do

END SUBROUTINE cabs_random_sub1





!******************************************************************************
!6tris) SUBROUTINE cext_rm_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cext_rm_sub(nstop,ns,dimvp,k,m_p,m_ab,v_cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns,dimvp                  !N espansioni multipolari,N sfere e N termini noti
REAL(dbl), INTENT(IN) :: k                                  !Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_p,m_ab          !Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cext              !Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j,jj,dimvp1
COMPLEX(dbl) :: somma

!Comincia la subroutine vera e propria
dimvp1=dimvp/2

direction_do: DO j=1,dimvp1

    !Inizializzo la variabile somma
    somma=(0.0D0,0.0D0)
    jj=0
    
    !Ciclo per il calcolo della sezione di estinzione per una direzione fissa.
    sphere_do: DO i=1,ns
    
            n_do: DO n=1,nstop
    
            m_do: DO m=-n,n
        
                !Aggiorno l'indice
                jj=jj+1
                
                !Costruisco la somma
                somma=somma+CONJG(m_p(2*jj-1,2*j))*m_ab(2*jj-1,2*j)+CONJG(m_p(2*jj,2*j))*m_ab(2*jj,2*j) &
                     &     +CONJG(m_p(2*jj-1,2*j-1))*m_ab(2*jj-1,2*j-1)+CONJG(m_p(2*jj,2*j-1))*m_ab(2*jj,2*j-1) 
        
            END DO m_do
    
        END DO n_do
    
    END DO sphere_do



    !Calcolo la sezione d'estinzione cumulativa, una direzione, luce non polarizzata
    v_cext(j)=(4*Pi_d/(k**2))*REAL(somma,dbl)
    
END DO direction_do

END SUBROUTINE cext_rm_sub


!******************************************************************************
!6quater) SUBROUTINE cext_tot_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cext_tot_sub(nalpha,nbeta,v_walpha,v_wbeta,v_cext,cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nalpha,nbeta
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_cext,v_walpha,v_wbeta     !Vettore cext e pesi per gauss
REAL(dbl), INTENT(OUT) :: cext                                     !Sezione di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj

!Inizio subroutine vera e propria

!Inizializzo cext e indici
cext=0.0D0
jj=0

!Ciclo per il calcolo di cext, integrando secondo gauss
alpha_do: DO i=1,nalpha
    beta_do: DO j=1,nbeta

        jj=jj+1
        cext=cext+v_walpha(i)*v_wbeta(j)*v_cext(jj)

    END DO beta_do
END DO alpha_do

cext=cext/(4.0D0*(Pi_d))

END SUBROUTINE cext_tot_sub



!******************************************************************************
!6tris) SUBROUTINE cext_rm_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cext_nopol_sub(nstop,ns,k,v_p1,v_p2,v_ab1,v_ab2,v_cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns                                !N espansioni multipolari,N sfere e N termini noti
REAL(dbl), INTENT(IN) :: k                                          !Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p1,v_p2,v_ab1,v_ab2     !Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cext              !Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j,jj,dimvp1
COMPLEX(dbl) :: somma

!Comincia la subroutine vera e propria



!Ciclo per il calcolo della sezione di estinzione per una direzione fissa.
j=0
sphere_do: DO i=1,ns

        !Inizializzo la variabile somma
        somma=(0.0D0,0.0D0)

        n_do: DO n=1,nstop
    
        m_do: DO m=-n,n
      
            !Aggiorno l'indice
            j=j+1
                
            !Costruisco la somma
            somma=somma+CONJG(v_p1(2*j-1))*v_ab1(2*j-1)+CONJG(v_p1(2*j))*v_ab1(2*j) &
                 &     +CONJG(v_p2(2*j-1))*v_ab2(2*j-1)+CONJG(v_p2(2*j))*v_ab2(2*j) 
        
        END DO m_do
    
    END DO n_do
    
    !Calcolo la sezione per una sfere, luce non polarizzata
     v_cext(i)=(4*Pi_d/(k**2))*REAL(somma,dbl)

END DO sphere_do

END SUBROUTINE cext_nopol_sub




!******************************************************************************
!7) SUBROUTINE cext_sub_exp: calcolo le sezioni di estinzione per ciascuna sfera
!conservando la sezione calcolata per ciascun ordine multipolare
!******************************************************************************
SUBROUTINE cext_sub_exp(nstop,ns,k,beta,v_z,v_ab,norm,m_cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_z					!Vettore coordinata z
REAL(dbl), INTENT(IN) :: k,beta								!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab				!Vettori coefficienti espansioni campo, non corretti e corretti  
CHARACTER(len=3), INTENT(IN) :: norm						!Flag normalizzazione

REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_cext			!Vettore sezioni di estinzione
! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,l_ib,up_ib,l_ia,up_ia,prec,elem
REAL(dbl) :: nr
COMPLEX(dbl) :: somma


!Funzione vera e propria
elem=2*nstop*(nstop+2)

norm_if: IF (norm=='old') THEN

		i_old_do: DO i=1,ns
		
			!Inizializzo la variabile somma
			somma=(0.0D0,0.0D0)
			
			!Mi porto sulla giusta posizione nel vettore v_ab, per la giusta sfera
			prec=(i-1)*elem
			
			n_old_do: DO n=1,nstop
			
				!Indice reale n
				nr=REAL(n,dbl)
				
				!Costruisco gli indici per recuperare i coefficienti
				l_ib=prec+2*(n*(n+1)-1)
				up_ib=prec+2*(n*(n+1)+1)
				l_ia=l_ib-1
				up_ia=up_ib-1
			
				!Riempio la matrice delle sezioni di estinzione
				somma=somma+(2.0D0*nr+1.0D0)*((v_ab(up_ia)+v_ab(up_ib))*EXP((0.0D0,1.0D0)*beta) - &
											& nr*(nr+1.0D0)*(v_ab(l_ia)-v_ab(l_ib))*EXP(-(0.0D0,1.0D0)*beta))
				m_cext(n,i)=(2*Pi_d/(k**2))*REAL(EXP(-(0.0D0,1.0D0)*k*v_z(i))*somma,dbl)
			
			END DO n_old_do
			
		END DO i_old_do
		
ELSE
		
		i_new_do: DO i=1,ns
		
			!Inizializzo la variabile somma
			somma=(0.0D0,0.0D0)
			
			!Mi porto sulla giusta posizione nel vettore v_ab, per la giusta sfera
			prec=(i-1)*elem
			
			n_new_do: DO n=1,nstop
			
				!Indice reale n
				nr=REAL(n,dbl)
				
				!Costruisco gli indici per recuperare i coefficienti
				l_ib=prec+2*(n*(n+1)-1)
				up_ib=prec+2*(n*(n+1)+1)
				l_ia=l_ib-1
				up_ia=up_ib-1
			
				!Riempio la matrice delle sezioni di estinzione
				somma=somma+SQRT(2.0D0*nr+1.0D0)*( (v_ab(up_ia)+v_ab(up_ib))*EXP((0.0D0,1.0D0)*beta) - &
												&  (v_ab(l_ia)-v_ab(l_ib))*EXP(-(0.0D0,1.0D0)*beta))
				m_cext(n,i)=(2*Pi_d/(k**2))*REAL(EXP(-(0.0D0,1.0D0)*k*v_z(i))*somma,dbl)
			
			END DO n_new_do
			
		END DO i_new_do

END IF norm_if

END SUBROUTINE cext_sub_exp



!******************************************************************************
!5) SUBROUTINE dmncmn_sub: calcolo i coefficienti di espansione del campo 
!interno per ciascuna sfera
!******************************************************************************
SUBROUTINE dmncmn_sub(nstop,ns,v_patt,m_da,m_cb,v_amnbmn,v_dmncmn)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore coeff espansioni campo esterno
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_da,m_cb		!Matrici rapporto coefficienti sfera singola d/a,c/b
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt			!Vettore pattern sfere identiche
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_dmncmn			!Vettori coefficienti espansioni campo interno

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,n,m


!Funzione vera e propria

!Inizializzo j e i vettori output
j=0
v_dmncmn=(0.0D0,0.0D0)

i_do: DO i=1,ns
	n_do: DO n=1,nstop
		m_do: DO m=-n,n
		
			!Incremento j
			j=j+1	
		
			!Calcolo dei coefficienti
			v_dmncmn(2*j-1)=m_da(n,v_patt(i))*v_amnbmn(2*j-1)
			v_dmncmn(2*j)=m_cb(n,v_patt(i))*v_amnbmn(2*j)
					
		END DO m_do
	END DO n_do
END DO i_do

END SUBROUTINE dmncmn_sub



!******************************************************************************
!5bis) SUBROUTINE dmncmn_sub_dip: calcolo i coefficienti di espansione del campo 
!interno per ciascuna sfera
!******************************************************************************
SUBROUTINE dmncmn_sub_dip(nstop,ns,v_patt,m_da,m_cb,v_amnbmn,v_dmncmn)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore coeff espansioni campo esterno
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_da,m_cb		!Matrici rapporto coefficienti sfera singola d/a,c/b
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt			!Vettore pattern sfere identiche
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_dmncmn			!Vettori coefficienti espansioni campo interno

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,n,m


!Funzione vera e propria

!Inizializzo j e i vettori output
j=0
v_dmncmn=(0.0D0,0.0D0)

i_do: DO i=1,ns
	n_do: DO n=1,nstop
		m_do: DO m=-n,n
		
			!Incremento j
			j=j+1	
		
			!Calcolo dei coefficienti
			v_dmncmn(2*j-1)=m_da(n,v_patt(i)-1)*v_amnbmn(2*j-1)
			v_dmncmn(2*j)=m_cb(n,v_patt(i)-1)*v_amnbmn(2*j)
					
		END DO m_do
	END DO n_do
END DO i_do

END SUBROUTINE dmncmn_sub_dip



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! SPARSE SUBROUTINES: funzioni e subroutines generali
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!******************************************************************************
!1) SUBROUTINE qsum_sparse: calcolo la somma di tutti i qmax nel caso sparse
!******************************************************************************
SUBROUTINE qsum_sparse(nstop,sumq,sumqq)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop              ! N espansioni multipolari
INTEGER(lo), INTENT(OUT) :: sumq,sumqq        ! vettore coeff qqmax

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,low_nu                  ! Indici

! Funzione vera e propria
sumq=0
sumqq=0

n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
   
            sumq=sumq+MIN(n,nu) + 1                         
            sumqq=sumqq+MIN(n,nu)                         

        END DO nu_do        
    END DO m_do
END DO n_do

END SUBROUTINE qsum_sparse



!******************************************************************************
!1bis) SUBROUTINE cgsum_sparse: calcolo la lunghezza del vettore che contiene in 
!fila tutti i coefficienti di clebsch-gordan
!******************************************************************************
SUBROUTINE cgsum_sparse(nstop,sumcg)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop		! N espansioni multipolari
INTEGER(lo), INTENT(OUT) :: sumcg		! vettore coeff qqmax

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,low_nu		! Indici

! Funzione vera e propria
sumcg=0

n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop

            sumcg=sumcg+n+nu-ABS(n-nu)+1

        END DO nu_do        
    END DO m_do
END DO n_do

END SUBROUTINE cgsum_sparse



!****************************************************************************************
!2) SUBROUTINE fill_aqbq_sparse:calcolo i vettori, tutti dritti, di coefficienti aq bq
!                 nel caso m=mu
!****************************************************************************************
SUBROUTINE fill_aqbq_sparse(nstop,v_aq_long,v_bq_long)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                              ! N espansioni multipolari e bounds vettori
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_aq_long          ! Vettore per tutti i coeff di gaunt
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_bq_long         ! Vettore per tutti i coeff di gaunt

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,error,qmax
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq,low_nu
REAL(dbl) :: mr,nr,mur,nur
REAL(dbl), DIMENSION(:), ALLOCATABLE :: v_aq,v_bq
REAL(dbl) :: p,p1,p2                                        !I soliti p,p1,p2
REAL(dbl) :: alphap1,alphap2,alphap3,alphap4                !I soliti coefficienti alpha
REAL(dbl) :: Ap2,Ap3,Ap4                                    !I soliti coefficienti Ap
REAL(dbl) :: c0,c1,c2                                       !I coefficienti ricorsivi
INTEGER(lo) :: q,Ap2i                                     !Come sopra ma integer

! Funzione vera e propria

!Inizializzo indici e bounds

!aq
i=0
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

n_do: DO n=1,nstop
    m_do: DO m=-n,n
        
        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
            
                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo qmax
                qmax=MIN(n,nu)

                !---------------------------------
                !Parte di storage per aq
                !---------------------------------
                
                !Aggiorno i
                i=i+1
                
                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_aq=1
                ELSE
                    low_aq=up_aq+1  
                END IF ia_if
                
                up_aq=low_aq+qmax

!                WRITE(*,*) "Bounds aq",low_aq,up_aq

                !Alloco la memoria per il mio v_aq
                ALLOCATE(v_aq(0:qmax)) 

                !Archivio i dati nel vettore lungo
                CALL gaunt_xu(-mr,nr,mur,nur,qmax,v_aq,error)
                v_aq_long(low_aq:up_aq)=v_aq
                
                
                !--------------------------------------------
                !Parte di storage per bq: se sono qui Qmax/=0
                !--------------------------------------------
                                
                !Aggiorno i bounds
                ib_if: IF (i==1) THEN
                    low_bq=1
                ELSE
                    low_bq=up_bq+1  
                END IF ib_if
                
                up_bq=low_bq+qmax-1

!                WRITE(*,*) "Bounds bq",low_bq,up_bq
                          
                !Alloco la memoria per il mio v_aq
                ALLOCATE(v_bq(1:qmax))
                
                
                !***************************************************************************
                !Calcolo bq che mi serve
                !***************************************************************************
                mmu_if: IF (m==0) THEN
                
                    !----------------
                    !Cosi' avro' B==0
                    !----------------
                    v_bq=0.0D0
                    
                ELSE    
                
                    !----------------
                    !Caso generale
                    !----------------
                    bq_case: SELECT CASE (qmax)
                
                    CASE(1) bq_case
                                        
                        !Calcolo coefficienti parziali
                        p=nr+nur-2.0D0
                        p1=p+mr-mur
                        Ap3=f_Ap(-mr,nr,mur,nur,p+3.0D0)
                        
                        !Calcolo Bq
                        v_bq(1)=v_aq(0)*(2.0D0*p+3.0D0)*Ap3/((p+3.0D0)*(p1+2.0D0))  
                        
                    CASE DEFAULT bq_case
                    
                        !Calcolo il primo valore di bq
                        p=nr+nur-2.0D0
                        p1=p+mr-mur
                        Ap3=f_Ap(-mr,nr,mur,nur,p+3.0D0)
                        v_bq(1)=v_aq(0)*(2.0D0*p+3.0D0)*Ap3/((p+3.0D0)*(p1+2.0D0))
                    
                        !Comincio il ciclo do per il calcolo di tutti i bq
                        bq_do: DO q=2,qmax
                        
                            !Calcolo preliminarmente p ed Ap2
                            p=nr+nur-2.0D0*REAL(q,dbl)
                            Ap2=f_Ap(-mr,nr,mur,nur,p+2.0D0)
                            Ap2i=INT(Ap2,lo)
                        
                            ap_case: SELECT CASE (Ap2i)
                            
                            CASE(0) ap_case
                            
                                !Calcolo coefficienti parziali
                                p1=p+mr-mur
                                p2=p-mr+mur
                                Ap3=f_Ap(-mr,nr,mur,nur,p+3.0D0)
                                Ap4=f_Ap(-mr,nr,mur,nur,p+4.0D0)
                                alphap3=f_alpha(nr,nur,p+3.0D0)
                                alphap4=f_alpha(nr,nur,p+4.0D0)
                                
                                !Calcolo coefficienti ricorsione
                                c0=(2.0D0*p+3.0D0)/((p+3.0D0)*(p1+2.0D0)*Ap4)
                                c1=Ap3*Ap4+(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*alphap3
                                c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*alphap4
                                
                                !Calcolo bq
                                v_bq(q)=c0*(c1*v_aq(q-1)+c2*v_aq(q-2))
                    
                                IF (v_bq(1)==0.0D0) THEN
                                    WRITE(*,*) m,n,mu,nu,c0,c1,c2
                                    WRITE(*,*) Ap3,Ap4,f_Ap(-mr,nr,mur,nur,p+5.0D0)
                                    WRITE(*,*)
                                END IF
                                                
                            CASE DEFAULT ap_case
                            
                                !Calcolo coefficienti parziali
                                p1=p+mr-mur
                                p2=p-mr+mur
                                alphap1=f_alpha(nr,nur,p+1.0D0)
                                alphap2=f_alpha(nr,nur,p+2.0D0)
                                
                                !Calcolo coefficienti ricorsione
                                c0=(2.0D0*p+3.0D0)/Ap2
                                c1=(p+2.0D0)*(p1+1.0D0)*alphap1
                                c2=-(p+1.0D0)*(p2+2.0D0)*alphap2
                                
                                
                                v_bq(q)=c0*(c1*v_aq(q)+c2*v_aq(q-1))
                               
                            
                            END SELECT ap_case
                            
                        END DO bq_do
                    
                    END SELECT bq_case
            
                END IF mmu_if

                v_bq_long(low_bq:up_bq)=v_bq
            
                !***************************************************************************
                !Fine calcolo bq
                !***************************************************************************
                
                !Disalloco v_aq e v_bq
                DEALLOCATE(v_aq,v_bq)
                                    
                !Esco ad un certo punto
!                IF ( i==( nstop*(1+2*nstop*(3+nstop)) )/3 ) EXIT n_do

        END DO nu_do        
    END DO m_do
END DO n_do

END SUBROUTINE fill_aqbq_sparse




!******************************************************************************
!2bis) SUBROUTINE fill_index: filling indexes and derived stuff once for all
!so i don't lose time when I fill up the matrix, or hopefully i make things
!easier for an openmp implementation
!******************************************************************************
SUBROUTINE fill_index(nstop,m_index,m_Apmn,v_jABij_template,v_iABij_template,v_mask,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(OUT) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_mask					! Masking vector for AB
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax				!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr						!Indici reali b
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq				!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes						!Estremi per i vett e n bes

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
error=0
v_iABij_template(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

n_do: DO n=1,nstop
	m_do: DO m=-n,n

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		nu_do: DO nu=low_nu,nstop

			!Rendo reali gli indici
			mr=REAL(m,dbl)
			nr=REAL(n,dbl)
			mur=REAL(m,dbl)
			nur=REAL(nu,dbl)

			!Calcolo qmax
			qmax=MIN(n,nu)

			!Aggiorno i
			i=i+1


			!Aggiorno i bounds
			ia_if: IF (i==1) THEN
				low_aq=1
			ELSE
				low_aq=up_aq+1  
			END IF ia_if

			ib_if: IF (i==1) THEN
				low_bq=1
			ELSE
				low_bq=up_bq+1  
			END IF ib_if

			up_aq=low_aq+qmax

			up_bq=low_bq+qmax-1

			!Calcolo Pmin e Pmax
			pmin=n+nu -2*qmax
			pmax=n+nu
			nbes=pmax-pmin+2

			!-----------
			!Calcolo Avt
			!-----------
			!Calcolo del fattore numerico sotto sommatoria per Avt
			DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				m_Apmn(i,p+1)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
			END DO

			!--------------------------------------------------------------
			!Infilo tutto quello che devo infilare nella mia matrice indice
			!--------------------------------------------------------------
			m_index(i,1)=n
			m_index(i,2)=m
			m_index(i,3)=nu
			m_index(i,4)=qmax
			m_index(i,5)=pmin
			m_index(i,6)=pmax
			m_index(i,7)=nbes
			m_index(i,8)=low_aq
			m_index(i,9)=up_aq
			m_index(i,10)=low_bq
			m_index(i,11)=up_bq

			!-----------------------------------------
			!Updatin indexing and filling mask vector
			!-----------------------------------------
			v_jABij_template(i)=nu*(nu+1)+m
			v_mask(i)=v_oner(n+nu)

		END DO nu_do

		!-----------------------------------
		!Dove comincia la prossima riga?Qui!
		!-----------------------------------
		v_iABij_template(n*(n+1)+m+1)=i+1

	END DO m_do
END DO n_do

END SUBROUTINE fill_index





!******************************************************************************************
!2tris) SUBROUTINE fill_TU_sparse_pardiso: riempio un blocco delle matrici Aij e Bij, TUij ma
!seguendo la data structure nuova o meglio, ora che posso, scrivo TU come matrice unica in
!fashion compressed sparse row, e poi scrivo Aij e Bij, che mi servono comunque, come al
!solito
!*****************************************************************************************
SUBROUTINE fill_index_TU(nstop,m_index,m_index_TU,m_Apmn,v_jABij,v_iABij,v_jTUij,v_iTUij,error)



IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_index,m_index_TU			! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(OUT) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij,v_jTUij,v_iTUij		! Vettori sparse colonne e righe
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: m,m1,n,nu,i,i1,i2,q,p,low_nu,qmax				!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr						!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq				!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes						!Estremi per i vett e n bes

! Funzione vera e propria

!Inizializzo indici e bounds
error=0

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1


!---------------------------------------------------------------------------------------
!Loop for the openMP implementation of the loops for the calculation of A and B
!---------------------------------------------------------------------------------------
AB_n_do: DO n=1,nstop
	AB_m_do: DO m=-n,n

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		AB_nu_do: DO nu=low_nu,nstop


			!Rendo reali gli indici
			mr=REAL(m,dbl)
			nr=REAL(n,dbl)
			mur=REAL(m,dbl)
			nur=REAL(nu,dbl)

			!Calcolo qmax
			qmax=MIN(n,nu)

			!Aggiorno i
			i=i+1

			!Aggiorno i bounds
			ia_if: IF (i==1) THEN
				low_aq=1
			ELSE
				low_aq=up_aq+1  
			END IF ia_if

			ib_if: IF (i==1) THEN
				low_bq=1
			ELSE
				low_bq=up_bq+1  
			END IF ib_if

			up_aq=low_aq+qmax

			up_bq=low_bq+qmax-1

			!Calcolo Pmin e Pmax
			pmin=n+nu -2*qmax
			pmax=n+nu
			nbes=pmax-pmin+2

			!-----------
			!Calcolo Avt
			!-----------
			!Calcolo del fattore numerico sotto sommatoria per Avt
			DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				m_Apmn(i,p+1)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
			END DO


			!--------------------------------------------------------------
			!Infilo tutto quello che devo infilare nella mia matrice indice
			!--------------------------------------------------------------
			m_index(i,1)=n
			m_index(i,2)=m
			m_index(i,3)=nu
			m_index(i,4)=qmax
			m_index(i,5)=pmin
			m_index(i,6)=pmax
			m_index(i,7)=nbes
			m_index(i,8)=low_aq
			m_index(i,9)=up_aq
			m_index(i,10)=low_bq
			m_index(i,11)=up_bq

			!-----------------------------------------
			!Aggiorno gli indici per lo storage sparse
			!-----------------------------------------
			v_jABij(i)=nu*(nu+1)+m

		END DO AB_nu_do

		!-----------------------------------
		!Dove comincia la prossima riga?Qui!
		!-----------------------------------
		v_iABij(n*(n+1)+m+1)=i+1

	END DO AB_m_do
END DO AB_n_do





!------------------------------------------------------------------------
! Loop for the indexes of the TU matrix, it is trickier because I have m1
!------------------------------------------------------------------------
!Indici
i=0			!Indice per TU
i1=0			!indice per A e B
i2=0			!indice per A e B
v_iTUij(1)=1

n_do: DO n=1,nstop
	m_do: DO m=-n,n

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		!Questo ciclo e' sembre sulle righe. Per TU le righe si sdoppiano, anche le colonne ma risolvo con i+2
		m1_do: DO m1=1,2

			nu_do: DO nu=low_nu,nstop


				!Aggiorno gli indici
				i=i+2						!Salto di due in due per quello di TU
!				WRITE(*,*) "n,m,nu,i",n,m,nu,i


				!Assegno i valori per la matrice TU, con un if su m1
				m1_if: IF (m1==1) THEN

					i1=i1+1

					!Index filling: 3 rows (n,m,m1) one column (nu) and one special progressive (i1)
					m_index_TU(i-1,1)=n
					m_index_TU(i,1)=n
					m_index_TU(i-1,2)=m
					m_index_TU(i,2)=m
					m_index_TU(i-1,3)=m1
					m_index_TU(i,3)=m1
					m_index_TU(i-1,4)=nu
					m_index_TU(i,4)=nu
					!Special treatment here, since i1 and i2 alernate in filling the 5th space
					m_index_TU(i-1,5)=i1
					m_index_TU(i,5)=i1

				ELSE

					i2=i2+1

					!Index filling: 3 rows (n,m,m1) one column (nu) and one special progressive (i2)
					m_index_TU(i-1,1)=n
					m_index_TU(i,1)=n
					m_index_TU(i-1,2)=m
					m_index_TU(i,2)=m
					m_index_TU(i-1,3)=m1
					m_index_TU(i,3)=m1
					m_index_TU(i-1,4)=nu
					m_index_TU(i,4)=nu
					!Special treatment here, since i1 and i2 alernate in filling the 5th space
					m_index_TU(i-1,5)=i2
					m_index_TU(i,5)=i2

				END IF m1_if

				!-----------------------------------------
				!Sparse storage indexes
				!-----------------------------------------
				v_jTUij(i-1)=2*(nu*(nu+1)+m-1)+1
				v_jTUij(i)=2*(nu*(nu+1)+m-1)+2

			END DO nu_do

			!-----------------------------------
			!Starting index of the next row
			!-----------------------------------
			v_iTUij(2*(n*(n+1)+m-1)+m1+1)=i+1

		END DO m1_do
	END DO m_do
END DO n_do


END SUBROUTINE fill_index_TU





!******************************************************************************
!3) SUBROUTINE c0_sparse(nstop,v_c0_sparse): coefficiente di normalizzazione per
!              i vector translation coefficients, calcolato nel caso sparse.
!******************************************************************************
SUBROUTINE c0_sparse(nstop,v_c0_sparse)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                   ! N espansioni multipolari
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_c0_sparse  ! vettore coeff c0

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,r,low_nu
REAL(dbl) :: mr,nr,mur,nur,logw,sqrtw

! Funzione vera e propria
i=0

n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
                
                i=i+1

                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)
                                
                logw =lnf(nr+mr,r)+lnf(nur-mur,r)-lnf(nr-mr,r)-lnf(nur+mur,r)
                
                sqrtw=EXP(logw)*((2*nr+1.0D0)*(2*nur+1.0D0))/(nr*(nr+1.0D0)*nur*(nur+1.0D0))
            
                v_c0_sparse(i)=0.5D0*((-1.0D0)**m)*SQRT(sqrtw)
            
        END DO nu_do        
    END DO m_do
END DO n_do

END SUBROUTINE c0_sparse


!******************************************************************************
!3) SUBROUTINE fillblock_sparse: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse.
!******************************************************************************
SUBROUTINE fillblock_sparse(nstop,v_c0,v_aq,v_bq,kr,v_Aij,v_Bij,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                                  ! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr                                         ! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               ! Vettore norm,gaunt e bq
INTEGER(lo), INTENT(OUT) :: error                                 ! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij            ! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij          ! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax                   !Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr                               !Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq                  !Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes                             !Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB                               !Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h              !Vett comp fun hankel
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn              !Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_h(0:(2*nstop+1)))

!Calcolo infine la mia funzione di Hankel
CALL hankel1_d_sub((2*nstop+1),kr,v_h,error)


n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
            

                
                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo qmax
                qmax=MIN(n,nu)
                

                !Aggiorno i
                i=i+1

                
                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_aq=1
                ELSE
                    low_aq=up_aq+1  
                END IF ia_if

                ib_if: IF (i==1) THEN
                    low_bq=1
                ELSE
                    low_bq=up_bq+1  
                END IF ib_if
                
                up_aq=low_aq+qmax

                up_bq=low_bq+qmax-1
                    
                !Calcolo Pmin e Pmax
                pmin=n+nu -2*qmax
                pmax=n+nu
                nbes=pmax-pmin+2
                
                     
                !Adesso posso allocare tutti ma proprio tutti i miei vettori
!                ALLOCATE(v_h(0:pmax+1),v_Apmn(pmin:pmax+1))
                ALLOCATE(v_Apmn(pmin:pmax+1))
                    
                !Inizializzo v_Apmn
                v_Apmn=0.0D0
                    
                    
                !-------
!                !Calcolo infine la mia funzione di Hankel
!                CALL hankel1_d_sub(pmax+1,kr,v_h,error)
                    
                    
                !-----------
                !Calcolo Avt
                !-----------
                !Calcolo del fattore numerico sotto sommatoria per Avt
                DO q=0,qmax
                    pr=nr+nur-2.0D0*REAL(q,dbl)
                    p=INT(pr,lo)
                    v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
                END DO
                    
              
                !Calcolo della sommatoria e di Avt
                sommaA=(0.0D0,0.0D0)
                sommaA_do:DO q=0,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_h(p)
                END DO sommaA_do
                    
                v_Aij(i)=v_c0(i)*sommaA
                    
                !-----------
                !Calcolo Bvt
                !-----------    
                !Calcolo della sommatoria e di Bvt
                sommaB=(0.0D0,0.0D0)
                sommaB_do:DO q=1,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
                END DO sommaB_do
                    
                v_Bij(i)=v_c0(i)*sommaB
                    
                DEALLOCATE(v_Apmn)    
                
                !-----------
                !Aggiorno gli indici per lo storage sparse
                !-----------   
                v_jABij(i)=nu*(nu+1)+m
      
        END DO nu_do

       !-----------
       !Dove comincia la prossima riga?Qui!
       !-----------   
       v_iABij(n*(n+1)+m+1)=i+1 
        
    END DO m_do
END DO n_do

END SUBROUTINE fillblock_sparse



!******************************************************************************
!3bis) SUBROUTINE fillblock_sparse_dip: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fillblock_sparse_dip(nstop,v_c0,v_aq,v_bq,kr,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                                  			! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr                                         			! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
INTEGER(lo), INTENT(OUT) :: error                                 			! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca      ! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij         			! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax                   !Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr                               !Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq                  !Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes                             !Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca           !Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h              !Vett comp fun hankel
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn              !Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1



n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
            

                
                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo qmax
                qmax=MIN(n,nu)
                

                !Aggiorno i
                i=i+1

                
                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_aq=1
                ELSE
                    low_aq=up_aq+1  
                END IF ia_if

                ib_if: IF (i==1) THEN
                    low_bq=1
                ELSE
                    low_bq=up_bq+1  
                END IF ib_if
                
                up_aq=low_aq+qmax

                up_bq=low_bq+qmax-1
                    
                !Calcolo Pmin e Pmax
                pmin=n+nu -2*qmax
                pmax=n+nu
                nbes=pmax-pmin+2
                
                     
                !Adesso posso allocare tutti ma proprio tutti i miei vettori
                ALLOCATE(v_h(0:pmax+1),v_Apmn(pmin:pmax+1))
                    
                !Inizializzo v_Apmn
                v_Apmn=0.0D0
                    
                    
                !-------
                !Calcolo infine la mia funzione di Hankel
                CALL hankel1_d_sub(pmax+1,kr,v_h,error)
                    
                    
                !-----------
                !Calcolo Avt
                !-----------
                !Calcolo del fattore numerico sotto sommatoria per Avtfill_aqbq_sparse
                DO q=0,qmax
                    pr=nr+nur-2.0D0*REAL(q,dbl)
                    p=INT(pr,lo)
                    v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
                END DO
                    
              
                !Calcolo della sommatoria e di Avt
                sommaA=(0.0D0,0.0D0)
                sommaAsca=(0.0D0,0.0D0)
                sommaA_do:DO q=0,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_h(p)
                            sommaAsca=sommaAsca+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
                END DO sommaA_do
                    
                v_Aij(i)=v_c0(i)*sommaA
                v_Aij_sca(i)=v_c0(i)*sommaAsca
                    
                !-----------
                !Calcolo Bvt
                !-----------    
                !Calcolo della sommatoria e di Bvt
                sommaB=(0.0D0,0.0D0)
                sommaBsca=(0.0D0,0.0D0)
                sommaB_do:DO q=1,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
                            sommaBsca=sommaBsca+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
                END DO sommaB_do
                    
                v_Bij(i)=v_c0(i)*sommaB
                v_Bij_sca(i)=v_c0(i)*sommaBsca
                    
                DEALLOCATE(v_h,v_Apmn)    
                
                !-----------
                !Aggiorno gli indici per lo storage sparse
                !-----------   
                v_jABij(i)=nu*(nu+1)+m
      
        END DO nu_do

       !-----------
       !Dove comincia la prossima riga?Qui!
       !-----------   
       v_iABij(n*(n+1)+m+1)=i+1 
        
    END DO m_do
END DO n_do

END SUBROUTINE fillblock_sparse_dip






!******************************************************************************
!3bisbis) SUBROUTINE fillblock_sparse_dip_nosca: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fillblock_sparse_dip_nosca(nstop,v_c0,v_aq,v_bq,kr,v_Aij,v_Bij,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                                  			! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr                                         			! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
INTEGER(lo), INTENT(OUT) :: error                                 			! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij					        ! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij         			! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax                   !Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr                               !Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq                  !Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes                             !Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB					            !Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h              !Vett comp fun hankel
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn              !Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1



n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
            

                
                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo qmax
                qmax=MIN(n,nu)
                

                !Aggiorno i
                i=i+1

                
                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_aq=1
                ELSE
                    low_aq=up_aq+1  
                END IF ia_if

                ib_if: IF (i==1) THEN
                    low_bq=1
                ELSE
                    low_bq=up_bq+1  
                END IF ib_if
                
                up_aq=low_aq+qmax

                up_bq=low_bq+qmax-1
                    
                !Calcolo Pmin e Pmax
                pmin=n+nu -2*qmax
                pmax=n+nu
                nbes=pmax-pmin+2
                
                     
                !Adesso posso allocare tutti ma proprio tutti i miei vettori
                ALLOCATE(v_h(0:pmax+1),v_Apmn(pmin:pmax+1))
                    
                !Inizializzo v_Apmn
                v_Apmn=0.0D0
                    
                    
                !-------
                !Calcolo infine la mia funzione di Hankel
                CALL hankel1_d_sub(pmax+1,kr,v_h,error)
                    
                    
                !-----------
                !Calcolo Avt
                !-----------
                !Calcolo del fattore numerico sotto sommatoria per Avt
                DO q=0,qmax
                    pr=nr+nur-2.0D0*REAL(q,dbl)
                    p=INT(pr,lo)
                    v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
                END DO
                    
              
                !Calcolo della sommatoria e di Avt
                sommaA=(0.0D0,0.0D0)
                sommaA_do:DO q=0,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_h(p)
                END DO sommaA_do
                    
                v_Aij(i)=v_c0(i)*sommaA
                    
                !-----------
                !Calcolo Bvt
                !-----------    
                !Calcolo della sommatoria e di Bvt
                sommaB=(0.0D0,0.0D0)
                sommaB_do:DO q=1,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
                END DO sommaB_do
                    
                v_Bij(i)=v_c0(i)*sommaB
                    
                DEALLOCATE(v_h,v_Apmn)    
                
                !-----------
                !Aggiorno gli indici per lo storage sparse
                !-----------   
                v_jABij(i)=nu*(nu+1)+m
      
        END DO nu_do

       !-----------
       !Dove comincia la prossima riga?Qui!
       !-----------   
       v_iABij(n*(n+1)+m+1)=i+1 
        
    END DO m_do
END DO n_do

END SUBROUTINE fillblock_sparse_dip_nosca





!******************************************************************************
!3tris) SUBROUTINE fillblock_sparse_sca: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,kr,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                                  			! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr                                         			! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
INTEGER(lo), INTENT(OUT) :: error                                 			! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca      ! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij         			! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax                   !Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr                               !Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq                  !Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes                             !Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca           !Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h              !Vett comp fun hankel
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn              !Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_h(0:(2*nstop+1)))

!Calcolo infine la mia funzione di Hankel
CALL hankel1_d_sub((2*nstop+1),kr,v_h,error)

n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop


                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo qmax
                qmax=MIN(n,nu)


                !Aggiorno i
                i=i+1


                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_aq=1
                ELSE
                    low_aq=up_aq+1  
                END IF ia_if

                ib_if: IF (i==1) THEN
                    low_bq=1
                ELSE
                    low_bq=up_bq+1  
                END IF ib_if

                up_aq=low_aq+qmax

                up_bq=low_bq+qmax-1

                !Calcolo Pmin e Pmax
                pmin=n+nu -2*qmax
                pmax=n+nu
                nbes=pmax-pmin+2


!                !Adesso posso allocare tutti ma proprio tutti i miei vettori
!                ALLOCATE(v_h(0:pmax+1),v_Apmn(pmin:pmax+1))
                ALLOCATE(v_Apmn(pmin:pmax+1))

                !Inizializzo v_Apmn
                v_Apmn=0.0D0


!                !-------
!                !Calcolo infine la mia funzione di Hankel
!                CALL hankel1_d_sub(pmax+1,kr,v_h,error)


                !-----------
                !Calcolo Avt
                !-----------
                !Calcolo del fattore numerico sotto sommatoria per Avt
                DO q=0,qmax
                    pr=nr+nur-2.0D0*REAL(q,dbl)
                    p=INT(pr,lo)
                    v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
                END DO


                !Calcolo della sommatoria e di Avt
                sommaA=(0.0D0,0.0D0)
                sommaAsca=(0.0D0,0.0D0)
                sommaA_do:DO q=0,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_h(p)
                            sommaAsca=sommaAsca+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
                END DO sommaA_do

                v_Aij(i)=v_c0(i)*sommaA
                v_Aij_sca(i)=v_c0(i)*sommaAsca

                !-----------
                !Calcolo Bvt
                !-----------    
                !Calcolo della sommatoria e di Bvt
                sommaB=(0.0D0,0.0D0)
                sommaBsca=(0.0D0,0.0D0)
                sommaB_do:DO q=1,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
                            sommaBsca=sommaBsca+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
                END DO sommaB_do

                v_Bij(i)=v_c0(i)*sommaB
                v_Bij_sca(i)=v_c0(i)*sommaBsca

                DEALLOCATE(v_Apmn)
!                DEALLOCATE(v_h,v_Apmn)

                !-----------
                !Aggiorno gli indici per lo storage sparse
                !-----------   
                v_jABij(i)=nu*(nu+1)+m

        END DO nu_do

       !-----------
       !Dove comincia la prossima riga?Qui!
       !-----------   
       v_iABij(n*(n+1)+m+1)=i+1

    END DO m_do
END DO n_do

END SUBROUTINE fillblock_sparse_sca




!******************************************************************************
!3quater) SUBROUTINE fillblock_sparse_sca: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fillblock_sparse_sca_fast(nstop,v_c0,v_aq,v_bq,kr,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				&    v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr								! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca		! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij				! Vettori sparse colonne e righe
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax			!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr					!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq			!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes,NNZab				!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca		!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h			!Vett comp fun hankel
INTEGER(lo) :: chunk=1					!Thread id, number of total threads and chunk size

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
NNZab=nstop*(1+2*nstop*(3+nstop))	!Computing NNZab,so i can have the cycle on the whole number of filled matrix values
NNZab=NNZab/3

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of hankel functions
ALLOCATE(v_h(0:(2*nstop+1)))
CALL hankel1_d_sub((2*nstop+1),kr,v_h,error)

h_if: IF (error/=0) THEN                                                                  
	WRITE (*,50) error
	50 FORMAT ("Si e' verificato un errore in fillblock_sparse_sca_fast, error= ",I5 //,"Il programma termina ora...") 
	STOP
END IF h_if




!This is a flattened do for the block filling, as in the standard routine, it is meant to be openMP friendly

!$OMP PARALLEL PRIVATE(n,m,nu,qmax,pmin,pmax,low_aq,up_aq,low_bq,up_bq,q,p,pr,sommaA,sommaAsca,sommaB,sommaBsca)

!Choosing the chunk size
chunk=NNZab/(200*2)
chunk_if: IF (chunk==0) THEN
	chunk=1
END IF chunk_if
!WRITE (*,*) 'Chunk length =',chunk

!$OMP DO SCHEDULE(guided,chunk)
flat_do: DO i=1,NNZab

	!----------------------------------------------------------------------------------------
	!Assigning all the needed indexes, which will be private, I assume, in the openMP version
	!----------------------------------------------------------------------------------------
	n=m_index(i,1)
	m=m_index(i,2)
	nu=m_index(i,3)
	qmax=m_index(i,4)
	pmin=m_index(i,5)
	pmax=m_index(i,6)
	low_aq=m_index(i,8)
	up_aq=m_index(i,9)
	low_bq=m_index(i,10)
	up_bq=m_index(i,11)

	!-----------
	!Calcolo Avt
	!-----------
	!Calcolo della sommatoria e di Avt
	sommaA=(0.0D0,0.0D0)
	sommaAsca=(0.0D0,0.0D0)
	sommaA_do:DO q=0,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
!		sommaA=sommaA+((0.0D0,1.0D0)**p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_h(p)
!		sommaAsca=sommaAsca+((0.0D0,1.0D0)**p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
		sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_h(p)
		sommaAsca=sommaAsca+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
	END DO sommaA_do

	v_Aij(i)=v_c0(i)*sommaA
	v_Aij_sca(i)=v_c0(i)*sommaAsca

	!-----------
	!Calcolo Bvt
	!----------- 
	!Calcolo della sommatoria e di Bvt
	sommaB=(0.0D0,0.0D0)
	sommaBsca=(0.0D0,0.0D0)
	sommaB_do:DO q=1,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
!		sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
!		sommaBsca=sommaBsca+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
		sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_h(p+1)
		sommaBsca=sommaBsca+v_ic(p+1)*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
	END DO sommaB_do

	v_Bij(i)=v_c0(i)*sommaB
	v_Bij_sca(i)=v_c0(i)*sommaBsca

END DO flat_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!Making a copy of the block filling pattern

v_jABij=v_jABij_template
v_iABij=v_iABij_template

END SUBROUTINE fillblock_sparse_sca_fast



!******************************************************************************
!3.5) SUBROUTINE fillblock_sparse_sca: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fillblock_sparse_fast(nstop,v_c0,v_aq,v_bq,kr,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				&    v_Aij,v_Bij,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr								! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij					! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij				! Vettori sparse colonne e righe
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax			!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr					!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq			!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes,NNZab				!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB					!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h			!Vett comp fun hankel
INTEGER(lo) :: chunk=1					!Thread id, number of total threads and chunk size

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
NNZab=nstop*(1+2*nstop*(3+nstop))	!Computing NNZab,so i can have the cycle on the whole number of filled matrix values
NNZab=NNZab/3

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of hankel functions
ALLOCATE(v_h(0:(2*nstop+1)))
CALL hankel1_d_sub((2*nstop+1),kr,v_h,error)

!Making a copy of the block filling pattern
v_jABij=v_jABij_template
v_iABij=v_iABij_template

!This is a flattened do for the block filling, as in the standard routine, it is meant to be openMP friendly

!$OMP PARALLEL PRIVATE(n,m,nu,qmax,pmin,pmax,low_aq,up_aq,low_bq,up_bq,q,p,pr,sommaA,sommaB)

!Choosing the chunk size
chunk=NNZab/(200*2)
chunk_if: IF (chunk==0) THEN
	chunk=1
END IF chunk_if
!WRITE (*,*) 'Chunk length =',chunk

!$OMP DO SCHEDULE(guided,chunk)
flat_do: DO i=1,NNZab

	!----------------------------------------------------------------------------------------
	!Assigning all the needed indexes, which will be private, I assume, in the openMP version
	!----------------------------------------------------------------------------------------
	n=m_index(i,1)
	m=m_index(i,2)
	nu=m_index(i,3)
	qmax=m_index(i,4)
	pmin=m_index(i,5)
	pmax=m_index(i,6)
	low_aq=m_index(i,8)
	up_aq=m_index(i,9)
	low_bq=m_index(i,10)
	up_bq=m_index(i,11)

	!-----------
	!Calcolo Avt
	!-----------
	!Calcolo della sommatoria e di Avt
	sommaA=(0.0D0,0.0D0)
	sommaA_do:DO q=0,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
!		sommaA=sommaA+((0.0D0,1.0D0)**p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_h(p)
!		sommaAsca=sommaAsca+((0.0D0,1.0D0)**p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
		sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_h(p)
	END DO sommaA_do

	v_Aij(i)=v_c0(i)*sommaA

	!-----------
	!Calcolo Bvt
	!----------- 
	!Calcolo della sommatoria e di Bvt
	sommaB=(0.0D0,0.0D0)
	sommaB_do:DO q=1,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
!		sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
!		sommaBsca=sommaBsca+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
		sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_h(p+1)
	END DO sommaB_do

	v_Bij(i)=v_c0(i)*sommaB

END DO flat_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE fillblock_sparse_fast


!******************************************************************************
!3.6) SUBROUTINE fillblock_sparse_onlysca: riempio un blocco solo della matrice di
!scattering
!******************************************************************************
SUBROUTINE fillblock_sparse_onlysca(nstop,v_c0,v_aq,v_bq,kr,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                                  			! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr                                         			! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
INTEGER(lo), INTENT(OUT) :: error                                 			! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij_sca,v_Bij_sca      			! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij         			! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax                   !Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr                               !Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq                  !Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes                             !Estremi per i vett e n bes
COMPLEX(dbl) :: sommaAsca,sommaBsca           !Somme parziali per vec trans
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_j                 !Vett comp fun sferiche bessel
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn              !Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_j(0:(2*nstop+1)))

!Calcolo infine la mia funzione di Hankel
CALL besselj_d_sub((2*nstop+1),kr,v_j,error)

n_do: DO n=1,nstop
    m_do: DO m=-n,n

        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop
            

                
                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo qmax
                qmax=MIN(n,nu)
                

                !Aggiorno i
                i=i+1

                
                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_aq=1
                ELSE
                    low_aq=up_aq+1  
                END IF ia_if

                ib_if: IF (i==1) THEN
                    low_bq=1
                ELSE
                    low_bq=up_bq+1  
                END IF ib_if
                
                up_aq=low_aq+qmax

                up_bq=low_bq+qmax-1
                    
                !Calcolo Pmin e Pmax
                pmin=n+nu -2*qmax
                pmax=n+nu
                nbes=pmax-pmin+2
                
                     
                !Adesso posso allocare tutti ma proprio tutti i miei vettori
                ALLOCATE(v_Apmn(pmin:pmax+1))
                    
                !Inizializzo v_Apmn
                v_Apmn=0.0D0
                    
                    
                !-------
                !Calcolo infine la mia funzione di Hankel
                CALL besselj_d_sub(pmax+1,kr,v_j,error)
                    
                    
                !-----------
                !Calcolo Avt
                !-----------
                !Calcolo del fattore numerico sotto sommatoria per Avt
                DO q=0,qmax
                    pr=nr+nur-2.0D0*REAL(q,dbl)
                    p=INT(pr,lo)
                    v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
                END DO
                    
              
                !Calcolo della sommatoria e di Avt
                sommaAsca=(0.0D0,0.0D0)
                sommaA_do:DO q=0,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaAsca=sommaAsca+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_j(p)
                END DO sommaA_do
                    
                v_Aij_sca(i)=v_c0(i)*sommaAsca
                    
                !-----------
                !Calcolo Bvt
                !-----------    
                !Calcolo della sommatoria e di Bvt
                sommaBsca=(0.0D0,0.0D0)
                sommaB_do:DO q=1,qmax
                            pr=nr+nur-2.0D0*REAL(q,dbl)
                            p=INT(pr,lo)
                            sommaBsca=sommaBsca+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_j(p+1)
                END DO sommaB_do
                    
                v_Bij_sca(i)=v_c0(i)*sommaBsca
                    
                DEALLOCATE(v_Apmn)    
                
                !-----------
                !Aggiorno gli indici per lo storage sparse
                !-----------   
                v_jABij(i)=nu*(nu+1)+m
      
        END DO nu_do

       !-----------
       !Dove comincia la prossima riga?Qui!
       !-----------   
       v_iABij(n*(n+1)+m+1)=i+1 
        
    END DO m_do
END DO n_do

END SUBROUTINE fillblock_sparse_onlysca


!******************************************************************************
!3.7) SUBROUTINE fillblock_sparse_sca: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fillblock_sparse_onlysca_fast(nstop,v_c0,v_aq,v_bq,kr,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				&    v_Aij,v_Bij,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr								! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij					! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij				! Vettori sparse colonne e righe
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax			!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr					!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq			!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes,NNZab				!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB					!Somme parziali per vec trans
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_j			!Vett comp fun bessel
INTEGER(lo) :: chunk=1					!Thread id, number of total threads and chunk size

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
NNZab=nstop*(1+2*nstop*(3+nstop))	!Computing NNZab,so i can have the cycle on the whole number of filled matrix values
NNZab=NNZab/3

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of hankel functions
ALLOCATE(v_j(0:(2*nstop+1)))
CALL besselj_d_sub((2*nstop+1),kr,v_j,error)

!This is a flattened do for the block filling, as in the standard routine, it is meant to be openMP friendly

!$OMP PARALLEL PRIVATE(n,m,nu,qmax,pmin,pmax,low_aq,up_aq,low_bq,up_bq,q,p,pr,sommaA,sommaB)

!Choosing the chunk size
chunk=NNZab/(200*2)
chunk_if: IF (chunk==0) THEN
	chunk=1
END IF chunk_if
!WRITE (*,*) 'Chunk length =',chunk

!$OMP DO SCHEDULE(guided,chunk)
flat_do: DO i=1,NNZab

	!----------------------------------------------------------------------------------------
	!Assigning all the needed indexes, which will be private, I assume, in the openMP version
	!----------------------------------------------------------------------------------------
	n=m_index(i,1)
	m=m_index(i,2)
	nu=m_index(i,3)
	qmax=m_index(i,4)
	pmin=m_index(i,5)
	pmax=m_index(i,6)
	low_aq=m_index(i,8)
	up_aq=m_index(i,9)
	low_bq=m_index(i,10)
	up_bq=m_index(i,11)

	!-----------
	!Calcolo Avt
	!-----------
	!Calcolo della sommatoria e di Avt
	sommaA=(0.0D0,0.0D0)
	sommaA_do:DO q=0,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
		sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_j(p)
	END DO sommaA_do

	v_Aij(i)=v_c0(i)*sommaA

	!-----------
	!Calcolo Bvt
	!----------- 
	!Calcolo della sommatoria e di Bvt
	sommaB=(0.0D0,0.0D0)
	sommaB_do:DO q=1,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
		sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_j(p+1)
	END DO sommaB_do

	v_Bij(i)=v_c0(i)*sommaB

END DO flat_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!Making a copy of the block filling pattern
v_jABij=v_jABij_template
v_iABij=v_iABij_template

END SUBROUTINE fillblock_sparse_onlysca_fast


!******************************************************************************
!4) FUNCTION id: per calcolare l'indice nei miei vettori di storage sparse
!******************************************************************************
FUNCTION id(n,m,k)

IMPLICIT NONE

! Dichiarazione funzione
INTEGER(lo) :: id

! Dichiarazione argomenti
INTEGER(lo), INTENT(IN) :: n,m,k            ! Indici coeff di gaunt e numero esp multipolari

! Dichiarazione variabili interne
INTEGER(lo) :: somma

! Funzione vera e propria
somma=4*(n**3)+2*n
somma=somma/3

id=somma+(2*n+1)*(n+m)+k

END FUNCTION id




!******************************************************************************
!5) SUBROUTINE fillblock_Dkmn:riempio un blocco ij dei coefficienti Dkmn  
!******************************************************************************
SUBROUTINE fillblock_Dkmn(nstop,theta,v_Dkmn,v_jDkmn,v_iDkmn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                               ! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta                                   ! Angolo theta coordinate relative
INTEGER(lo), INTENT(OUT) :: error                              ! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_Dkmn                   ! Vettore sparse di output valori Dkmn
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jDkmn,v_iDkmn      ! Vettori indici colonne e righe

!Dichiarazione variabili interne 
INTEGER(lo) :: i,n,m,k,lowb,upb,nrow,ncol,modm,maxmk      !Indici
INTEGER(lo) :: imk,ikm,i_k_m,i_m_k,imm,i_m_m,im_m,i_mm    !Indici blocco
REAL(dbl), DIMENSION(0:nstop) :: v_d,v_d1                   ! Vettori fattori Dkmn


!Subroutine vera e propria
error=0
v_d=0.0D0
v_d1=0.0D0

!Check di errore su theta
thetain_check: IF ((theta>PI_D) .OR. (theta<-PI_D)) THEN
                WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura fillblock_Dkmn si ferma"
                error=1
                RETURN
END IF thetain_check

!Inizializzo degli indici
i=0
v_iDkmn(1)=1
!v_Dkmn=0.0D0

!Assegno gli indici di righe e colonne prima di tutto
block_do: DO n=1,nstop

            !Bounds
            lowb=n**2
            upb=n*(n+2)

            row_do: DO nrow=lowb,upb
                col_do: DO ncol=lowb,upb

                    i=i+1
                    v_jDkmn(i)=ncol

!                    WRITE(*,*) lowb,upb,i

                END DO col_do

                v_iDkmn(nrow+1)=i+1

            END DO row_do

END DO block_do


!Assegno i valori lungo la diagonale per m=k=0
CALL d_nmk(nstop,0,0,theta,v_d,error)

diag_zero_do: DO n=1,nstop

    imm=id(n,0,0)

!    WRITE(*,*) "Diagonale zero:", imm

    v_Dkmn(imm)=v_d(n)
!    WRITE(*,*) v_d
!    WRITE(*,*)

END DO diag_zero_do


!Calcolo i valori che stanno sulle due diagonali principali:simmetria di 2
diag_do: DO m=1,nstop

    !Chiamata per ciascuna delle due diagonali
    CALL d_n_km(nstop,m,m,theta,v_d,error)
    CALL d_n_km(nstop,m,-m,theta,v_d1,error)
!            WRITE(*,*) m,k,v_d
!            WRITE(*,*)
!            WRITE(*,*) -m,k,v_d1
    !Faccio le assegnazioni saltando di blocco in blocco
    n_diag_do: DO n=m,nstop

            !Calcolo gli indici
            imm=id(n,m,m)
            i_m_m=id(n,-m,-m)
            im_m=id(n,m,-m)
            i_mm=id(n,-m,m)

!            WRITE(*,*) "Diagonale: ",imm,i_m_m,im_m,i_mm

            !Assegno i valori
            v_Dkmn(imm)=v_d(n)
            v_Dkmn(i_m_m)=v_d(n)
            v_Dkmn(im_m)=v_d1(n)
            v_Dkmn(i_mm)=v_d1(n)

    END DO n_diag_do

END DO diag_do


!Calcolo il caso generale
m_do: DO m=-(nstop-1),(nstop-1)
    k_do: DO k=ABS(m)+1,nstop

            CALL d_n_km(nstop,m,k,theta,v_d,error)
!            WRITE(*,*) m,k,v_d
!            WRITE(*,*)
            modm=ABS(m)
            maxmk=MAX(ABS(m),k)


            !Faccio le assegnazioni saltando di blocco in blocco
            n_do: DO n=maxmk,nstop

                !Sempre indici
                imk=id(n,m,k)
                ikm=id(n,k,m)
                i_k_m=id(n,-k,-m)
                i_m_k=id(n,-m,-k)
                
!                WRITE(*,*) "Generale: ", imk,ikm,i_k_m,i_m_k

                !Sempre valori
                v_Dkmn(ikm)=v_d(n)
                v_Dkmn(i_k_m)=((-1.0D0)**(m+k))*v_d(n)
                v_Dkmn(i_m_k)=v_d(n)
                v_Dkmn(imk)=v_Dkmn(i_k_m)

            END DO n_do
    END DO k_do
END DO m_do

END SUBROUTINE fillblock_Dkmn





!******************************************************************************
!6) SUBROUTINE nnb_sparse:calcolo il numero di blocchi diversi da zero. I blocchi
!              diagonali qui non sono inclusi
!******************************************************************************
SUBROUTINE nnb_sparse(ns,m_xyz,v_r,fint,nnb)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                      ! N sfere
REAL(dbl), INTENT(IN) :: fint                        ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r           ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz       ! Matrix posizione
INTEGER(lo), INTENT(OUT) :: nnb                    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f

! Funzione vera e propria
nnb=0
row_do: DO i=1,ns
    col_do: DO j=1,ns

        IF (i==j) CYCLE col_do

        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        xj=m_xyz(j,1)
        yj=m_xyz(j,2)
        zj=m_xyz(j,3)

        dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

        !Calcolo il coeff d'interazione
        f=(v_r(i)+v_r(j))/dij

        !Se la distanza e' piccola allora il blocco non e' zero
        nnb_if: IF (f>=fint) THEN
            nnb=nnb+1
        END IF nnb_if

    END DO col_do
END DO row_do

END SUBROUTINE nnb_sparse



!***********************************************************************************************************************************
!6bis) SUBROUTINE nnb_sparse_ant:calcolo il numero di blocchi diversi da zero. I blocchi diagonali non sono inclusi e scarto sempre 
! le ultime ndip righe.
!***********************************************************************************************************************************
SUBROUTINE nnb_sparse_ant(ns,ns_ant,m_xyz,v_r,fint,nnb_ant)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns				! N sfere
INTEGER(lo), INTENT(IN) :: ns_ant			! N dipoli
REAL(dbl), INTENT(IN) :: fint				! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r		! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz		! Matrix posizione
INTEGER(lo), INTENT(OUT) :: nnb_ant			! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f

! Subroutine

!First loop for the interaction among spheres
nnb_ant=0
row_ant_do: DO i=1,ns_ant
	col_ant_do: DO j=1,ns_ant

		IF (i==j) CYCLE col_ant_do

		!Coordinate e distanze
		xi=m_xyz(i,1)
		yi=m_xyz(i,2)
		zi=m_xyz(i,3)

		xj=m_xyz(j,1)
		yj=m_xyz(j,2)
		zj=m_xyz(j,3)

		dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

		!Calcolo il coeff d'interazione
		f=(v_r(i)+v_r(j))/dij

		!Se la distanza e' piccola allora il blocco non e' zero
		nnb_ant_if: IF (f>=fint) THEN
			nnb_ant=nnb_ant+1
		END IF nnb_ant_if

	END DO col_ant_do
END DO row_ant_do
nnb_ant=nnb_ant/2

!Sec ond loop for the interaction among spheres and dipoles
row_do: DO i=1,ns_ant
	col_do: DO j=ns_ant+1,ns

		IF (i==j) CYCLE col_do

		!Coordinate e distanze
		xi=m_xyz(i,1)
		yi=m_xyz(i,2)
		zi=m_xyz(i,3)

		xj=m_xyz(j,1)
		yj=m_xyz(j,2)
		zj=m_xyz(j,3)

		dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

		!Calcolo il coeff d'interazione
		f=(v_r(i)+v_r(j))/dij

		!Se la distanza e' piccola allora il blocco non e' zero
		nnb_if: IF (f>=fint) THEN
			nnb_ant=nnb_ant+1
		END IF nnb_if

	END DO col_do
END DO row_do

END SUBROUTINE nnb_sparse_ant




!***********************************************************************************************************************************
!6bis) SUBROUTINE nnb_sparse_dip_per:calcolo il numero di blocchi diversi da zero. I blocchi diagonali non sono inclusi e scarto sempre 
! le ultime ndip righe.
!***********************************************************************************************************************************
SUBROUTINE nnb_sparse_dip_per(ns,ndip,m_xyz,v_r,fint,nnb)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                      ! N sfere
INTEGER(lo), INTENT(IN) :: ndip                    ! N dipoli
REAL(dbl), INTENT(IN) :: fint                        ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r           ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz       ! Matrix posizione
INTEGER(lo), INTENT(OUT) :: nnb                    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,icell,jcell
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f,sigx,sigy

! Funzione vera e propria
nnb=0
row_do: DO i=1,ns
    col_do: DO j=1,ns

        IF ((i==j) .OR. (i>(ns-ndip))) CYCLE col_do

        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        !Coordinate per la seconda sfere, comprese complicate traslazioni
        cell_row_do_nanb: DO icell=0,nacell,nacell
             cell_col_do_nanb: DO jcell=0,nbcell,nbcell
             
                    sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
                    sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2)))      

                    xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
                    yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
                    zj=m_xyz(j,3)

                    dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                    !Calcolo il coeff d'interazione
                    f=(v_r(i)+v_r(j))/dij

                    !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                    !WRITE(*,*) i,j,f,dij
                    IF (f>=fint) EXIT cell_row_do_nanb

              END DO cell_col_do_nanb
        END DO cell_row_do_nanb

        !Se la distanza e' piccola allora il blocco non e' zero
        nnb_if: IF (f>=fint) THEN
            nnb=nnb+1
        END IF nnb_if

    END DO col_do
END DO row_do

END SUBROUTINE nnb_sparse_dip_per




!******************************************************************************
!6bis) SUBROUTINE nnb_sparse_dip:calcolo il numero di blocchi diversi da zero. I blocchi
!              diagonali qui non sono inclusi,e scarto sempre la prima riga!!!
!******************************************************************************
SUBROUTINE nnb_sparse_dip_real(ns,ndip,nnb)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,ndip                 ! N sfere
INTEGER(lo), INTENT(OUT) :: nnb                    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j

! Funzione vera e propria
nnb=0
row_do: DO i=1,ns

	IF (i<ndip) CYCLE row_do

	col_do: DO j=1,ns

		IF (j>ndip) CYCLE col_do
		nnb=nnb+1

	END DO col_do
END DO row_do

END SUBROUTINE nnb_sparse_dip_real






!******************************************************************************
!7) SUBROUTINE fill_jBlock_iBlock_sparse:riempio i vettori indici riga e colonna
! per la disposizione dei blocchi di matrice.
!******************************************************************************
SUBROUTINE fill_jBlock_iBlock_sparse(ns,m_xyz,v_r,fint,v_jBlock,v_iBlock)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                                 ! N sfere
REAL(dbl), INTENT(IN) :: fint                                   ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                      ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                  ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(OUT) :: v_jBlock,v_iBlock    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,next
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f

! Funzione vera e propria
next=1
!Per il primo blocco e' sempre uno
v_iBlock(1)=1

row_do: DO i=1,ns
	col_do: DO j=1,ns

		IF (i>=j) CYCLE col_do

		!Coordinate e distanze
		xi=m_xyz(i,1)
		yi=m_xyz(i,2)
		zi=m_xyz(i,3)

		xj=m_xyz(j,1)
		yj=m_xyz(j,2)
		zj=m_xyz(j,3)

		dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

		!Calcolo il coeff d'interazione
		f=(v_r(i)+v_r(j))/dij

		!Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
		IF (f>=fint) THEN

			!Assegno i valori di colonna
			v_jBlock(next)=j
			next=next+1

		END IF

	END DO col_do

	v_iBlock(i+1)=next

END DO row_do

END SUBROUTINE fill_jBlock_iBlock_sparse




!******************************************************************************
!7bis) SUBROUTINE fill_jBlock_iBlock_sparse_dip:riempio i vettori indici riga e colonna
! per la disposizione dei blocchi di matrice, escludo la prima riga.
!******************************************************************************
SUBROUTINE fill_jBlock_iBlock_sparse_dip(ns,ndip,m_xyz,v_r,fint,v_jBlock,v_iBlock)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                                 ! N sfere
INTEGER(lo), INTENT(IN) :: ndip                               ! N dipoli
REAL(dbl), INTENT(IN) :: fint                                   ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                      ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                  ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(OUT) :: v_jBlock,v_iBlock    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,next
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f

! Funzione vera e propria
next=1
!Per il primo blocco e' sempre uno
v_iBlock(1)=1

row_do: DO i=1,ns
    col_do: DO j=1,ns

        IF ((i==j) .OR. (i>(ns-ndip))) CYCLE col_do

        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        xj=m_xyz(j,1)
        yj=m_xyz(j,2)
        zj=m_xyz(j,3)

        dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

        !Calcolo il coeff d'interazione
        f=(v_r(i)+v_r(j))/dij

        !Se la distanza e' grande o sono sulla prima riga allora il blocco e' zero e non tocco gli indici
        IF (f<fint) CYCLE col_do

        !Assegno i valori di colonna
        v_jBlock(next)=j
        next=next+1

    END DO col_do

    v_iBlock(i+1)=next

END DO row_do

END SUBROUTINE fill_jBlock_iBlock_sparse_dip





!******************************************************************************
!7bis-bis) SUBROUTINE fill_jBlock_iBlock_sparse_dip_per:riempio i vettori indici riga e colonna
! per la disposizione dei blocchi di matrice, escludo le ultime ndip righe per strutture
! periodiche
!******************************************************************************
SUBROUTINE fill_jBlock_iBlock_sparse_dip_per(ns,ndip,m_xyz,v_r,fint,v_jBlock,v_iBlock)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                                 ! N sfere
INTEGER(lo), INTENT(IN) :: ndip                               ! N dipoli
REAL(dbl), INTENT(IN) :: fint                                   ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                      ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                  ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(OUT) :: v_jBlock,v_iBlock    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,next,icell,jcell
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f,sigx,sigy

! Funzione vera e propria
next=1
!Per il primo blocco e' sempre uno
v_iBlock(1)=1

row_do: DO i=1,ns
    col_do: DO j=1,ns

        IF ((i==j) .OR. (i>(ns-ndip))) CYCLE col_do

        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        !Coordinate per la seconda sfere, comprese complicate traslazioni
        cell_row_do_nanb: DO icell=0,nacell,nacell
             cell_col_do_nanb: DO jcell=0,nbcell,nbcell
             
                    sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
                    sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2)))      

                    xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
                    yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
                    zj=m_xyz(j,3)

                    dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                    !Calcolo il coeff d'interazione
                    f=(v_r(i)+v_r(j))/dij

                    !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                    !WRITE(*,*) i,j,f,dij
                    IF (f>=fint) EXIT cell_row_do_nanb

              END DO cell_col_do_nanb
        END DO cell_row_do_nanb

!        !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
!        IF (f<fint) CYCLE col_do

        !Assegno i valori di colonna
        fint_if: IF (f>=fint) THEN
        	v_jBlock(next)=j
        	next=next+1
        END IF fint_if

    END DO col_do

    v_iBlock(i+1)=next

END DO row_do

END SUBROUTINE fill_jBlock_iBlock_sparse_dip_per






!******************************************************************************
!7bis) SUBROUTINE fill_jBlock_iBlock_sparse_dip_real:riempio i vettori indici riga e colonna
! per la disposizione dei blocchi di matrice, escludo la prima riga.
!******************************************************************************
SUBROUTINE fill_jBlock_iBlock_sparse_dip_real(ns,ndip,v_jBlock,v_iBlock)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,ndip                            ! N sfere
INTEGER(lo), DIMENSION(:),INTENT(OUT) :: v_jBlock,v_iBlock    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,next

! Funzione vera e propria
next=1
!Per il primo blocco e' sempre uno
v_iBlock(1)=1

row_do: DO i=1,ns

    col_do: DO j=1,ns

        IF (j>ndip) CYCLE col_do
	IF (i<=ndip) EXIT col_do

        !Assegno i valori di colonna
        v_jBlock(next)=j
        next=next+1

    END DO col_do

    v_iBlock(i+1)=next

END DO row_do

END SUBROUTINE fill_jBlock_iBlock_sparse_dip_real







!******************************************************************************
!8) SUBROUTINE fill_D_PHI_sparse:riempio tutta la struttura per i blocchi di 
! rotazione Dij e per Exp[phi_ij]
!******************************************************************************
SUBROUTINE fill_D_PHI_sparse(ns,nstop,nnb,m_xyz,v_jBlock,v_iBlock,m_Dnkm,m_jDnkm,m_iDnkm,m_exphi)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop,nnb                           ! N sfere, N multipolar exp e N nnz blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         ! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Dnkm                    ! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jDnkm,m_iDnkm       ! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop,1:nnb), INTENT(OUT) :: m_exphi !Matrice Esponenziali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!I cycle col_do if I am in the lower triangular part of the matrix
		IF (i>j) CYCLE col_do

		!-----Diagnostics over the cycle------------------
		!WRITE(*,*) "D,i,j",i,j
		!-------------------------------------------------

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire Dij
		CALL fillblock_Dkmn(nstop,thetaij,m_Dnkm(:,next),m_jDnkm(:,next),m_iDnkm(:,next),error)

		!Riempio anche la colonna degli esponenziali
		phi_do: DO m=-nstop,nstop

			m_exphi(m,next)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

		END DO phi_do

	END DO col_do

END DO row_do

END SUBROUTINE fill_D_PHI_sparse




!******************************************************************************
!8) SUBROUTINE update_D_PHI_sparse: updating the rotation matrix to account for
! the changing position of the dipole
!******************************************************************************
SUBROUTINE update_D_PHI_sparse(ns,nstop,nnb,m_xyz,v_jBlock,v_iBlock,m_Dnkm,m_jDnkm,m_iDnkm,m_exphi)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop,nnb                           ! N sfere, N multipolar exp e N nnz blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         ! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Dnkm                    ! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jDnkm,m_iDnkm       ! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop,1:nnb), INTENT(OUT) :: m_exphi !Matrice Esponenziali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!I cycle col_do if I am in the lower triangular part of the matrix
		IF ((i>j) .OR. (j<=ns_ant)) CYCLE col_do

		!-----Diagnostics over the cycle------------------
		!WRITE(*,*) "D,i,j",i,j
		!-------------------------------------------------

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire Dij
		CALL fillblock_Dkmn(nstop,thetaij,m_Dnkm(:,next),m_jDnkm(:,next),m_iDnkm(:,next),error)

		!Riempio anche la colonna degli esponenziali
		phi_do: DO m=-nstop,nstop

			m_exphi(m,next)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

		END DO phi_do

	END DO col_do

END DO row_do

END SUBROUTINE update_D_PHI_sparse




!******************************************************************************
!8bis) SUBROUTINE fill_D_PHI_sparse_singleblock:riempio tutta la struttura per  
! un singolo blocco di rotazione Dij e per Exp[phi_ij]
!******************************************************************************
SUBROUTINE fill_D_PHI_sparse_singleblock(nstop,theta,phi,v_Dnkm,v_jDnkm,v_iDnkm,v_exphi,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! N sfere, N multipolar exp e N nnz blocks
REAL(dbl), INTENT(IN) :: theta,phi						! Angoli per la singola matrice di rotazione
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_Dnkm				! Vettore per la matrice Dij
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jDnkm,v_iDnkm		! Vettori per gli indici della matrice Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop), INTENT(OUT) :: v_exphi	! Matrice Esponenziali
INTEGER(lo),  INTENT(OUT) :: error					! Vettori per gli indici della matrice Dij

! Dichiarazione variabili interne
INTEGER(lo) :: m

!Comincio la subroutine vera e propria
error=0

	WRITE(*,*) "qui1a"
!Chiamo la subroutine per riempire Dij
CALL fillblock_Dkmn(nstop,theta,v_Dnkm,v_jDnkm,v_iDnkm,error)
	WRITE(*,*) "qui2a"
!Check di errore su fillblock_Dkmn
thetain_check: IF (error==1) THEN
                WRITE(*,*) "errore in fillblock_Dkmn chiamata da fill_D_PHI_sparse_singleblock: il programma si ferma"
                error=1
                STOP
END IF thetain_check
	WRITE(*,*) "qui3a"
!Riempio anche la colonna degli esponenziali
phi_do: DO m=-nstop,nstop
	v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phi)
END DO phi_do
	WRITE(*,*) "qui4a"
	
END SUBROUTINE fill_D_PHI_sparse_singleblock







!******************************************************************************
!8tris) SUBROUTINE fill_D_PHI_sparse:riempio tutta la struttura per i blocchi di 
! rotazione Dij e per Exp[phi_ij]
!******************************************************************************
SUBROUTINE fill_D_PHI_sparse_rhs(ns,nstop,nnb,m_xyz,v_xyz_dip,v_jBlock,v_iBlock,m_Dnkm,m_jDnkm,m_iDnkm,m_exphi)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop,nnb                           ! N sfere, N multipolar exp e N nnz blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_xyz_dip                    ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         ! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Dnkm                    ! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jDnkm,m_iDnkm       ! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop,1:nnb), INTENT(OUT) :: m_exphi !Matrice Esponenziali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

        !Calcolo le coordinate sferiche relative
        xij=m_xyz(i,1)-v_xyz_dip(1)
        yij=m_xyz(i,2)-v_xyz_dip(2)
        zij=m_xyz(i,3)-v_xyz_dip(3)
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        !Chiamo la subroutine per riempire Dij
        CALL fillblock_Dkmn(nstop,thetaij,m_Dnkm(:,next),m_jDnkm(:,next),m_iDnkm(:,next),error)

        !Riempio anche la colonna degli esponenziali
        phi_do: DO m=-nstop,nstop

            m_exphi(m,next)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

        END DO phi_do

    END DO col_do

END DO row_do

END SUBROUTINE fill_D_PHI_sparse_rhs


!******************************************************************************
!9) SUBROUTINE fill_AB_sparse:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse(ns,nstop,k,m_xyz,v_c0,v_aq,v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop                               			! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k                                          			! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      			! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         			! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij					  	! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij       			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

        !Calcolo le coordinate sferiche relative
        xij=m_xyz(i,1)-m_xyz(j,1)
        yij=m_xyz(i,2)-m_xyz(j,2)
        zij =m_xyz(i,3)-m_xyz(j,3)
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        !Chiamo la subroutine per riempire Dij
        CALL fillblock_sparse(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &
                            & m_jABij(:,next),m_iABij(:,next),error)

    END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse



!******************************************************************************
!9) SUBROUTINE fill_AB_sparse:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_rhs(ns,nstop,k,m_xyz,v_xyz_dip,v_c0,v_aq,v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop                               			! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k                                          			! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      			! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_xyz_dip                    			! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         			! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij					  	! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij       			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

        !Calcolo le coordinate sferiche relative
        xij=m_xyz(i,1)-v_xyz_dip(1)
        yij=m_xyz(i,2)-v_xyz_dip(2)
        zij=m_xyz(i,3)-v_xyz_dip(3)
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        !Chiamo la subroutine per riempire Dij
        CALL fillblock_sparse(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &
                            & m_jABij(:,next),m_iABij(:,next),error)

    END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_rhs



!******************************************************************************
!9.3) SUBROUTINE fill_AB_sparse_dip:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_dip&
&(ns,ndip,nstop,k,m_xyz,v_c0,v_aq,v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,ndip,nstop                             			! N sfere, N multipoli,N multipolar exp 
REAL(dbl) ,INTENT(IN) :: k                                          			! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      			! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         			! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca  		! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij       			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

        !Calcolo le coordinate sferiche relative
        xij=m_xyz(i,1)-m_xyz(j,1)
        yij=m_xyz(i,2)-m_xyz(j,2)
        zij =m_xyz(i,3)-m_xyz(j,3)
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        dip_if: IF (i<=(ns-ndip)) THEN
        
        	!Chiamo la subroutine per riempire Aij e Aijsca
        	CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &  
                            & m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)
        ELSE
        
        	!Chiamo la subroutine per riempire solamente Aijsca, dal momento che sono nella zona dipoli
        	CALL fillblock_sparse_onlysca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij_sca(:,next),m_Bij_sca(:,next),&  
                            & m_jABij(:,next),m_iABij(:,next),error)
                            
        END IF dip_if

    END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_dip


!******************************************************************************
!9.3) SUBROUTINE fill_AB_sparse_dip_fast:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_dip_fast(ns,ns_ant,nstop,k,m_xyz,v_c0,v_aq,v_bq,&
				   m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				   v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,ns_ant,nstop						! N sfere, N multipoli,N multipolar exp 
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock				! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca  		! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij				! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!I cycle col_do in I am in the lower triangular part of the matrix
		IF (i>j) CYCLE col_do

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		ant_if: IF (i<=ns_ant) THEN

			!Filling both Aij and Aijsca
			CALL fillblock_sparse_sca_fast(nstop,v_c0,v_aq,v_bq,k*rij,&
			m_index,m_Apmn,v_jABij_template,v_iABij_template, &
			m_Aij(:,next),m_Bij(:,next),m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

		ELSE

			!Filling Aijsca for the scattering coefficient computation
			CALL fillblock_sparse_onlysca_fast(nstop,v_c0,v_aq,v_bq,k*rij,&
			m_index,m_Apmn,v_jABij_template,v_iABij_template, &
			m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

		END IF ant_if

	END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_dip_fast



!******************************************************************************
!9.4) SUBROUTINE update_AB_sparse_dip_fast: updating the sparse block matrix
!structure to acconunt for the displacement of the dipoles
!******************************************************************************
SUBROUTINE update_AB_sparse_dip_fast(ns,ns_ant,nstop,k,m_xyz,v_c0,v_aq,v_bq,&
				   m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				   v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,ns_ant,nstop						! N sfere, N multipoli,N multipolar exp 
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock				! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca		! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij				! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!I cycle col_do in I am in the lower triangular part of the matrix
		IF ((i>j) .OR. j<=ns_ant) CYCLE col_do

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		ant_if: IF (i<=ns_ant) THEN

			!Filling both Aij and Aijsca
			CALL fillblock_sparse_sca_fast(nstop,v_c0,v_aq,v_bq,k*rij,&
			m_index,m_Apmn,v_jABij_template,v_iABij_template, &
			m_Aij(:,next),m_Bij(:,next),m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

		ELSE

			!Filling Aijsca for the scattering coefficient computation
			CALL fillblock_sparse_onlysca_fast(nstop,v_c0,v_aq,v_bq,k*rij,&
			m_index,m_Apmn,v_jABij_template,v_iABij_template, &
			m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

		END IF ant_if

	END DO col_do

END DO row_do

END SUBROUTINE update_AB_sparse_dip_fast



!******************************************************************************
!8bis) SUBROUTINE fill_AB_sparse_dip_per:riempio tutta la struttura per i blocchi di 
! traslazione Aij e per Bij in un contesto periodico
!******************************************************************************
SUBROUTINE fill_AB_sparse_dip_per &
&(ns,ndip,nstop,k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,ndip,nstop                             			! N sfere, N multipoli,N multipolar exp 
REAL(dbl) ,INTENT(IN) :: k                                          			! vettore d'onda 
REAL(dbl), INTENT(IN) :: fint										! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r							! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      			! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         			! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca  	! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij       			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m,icell,jcell
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij,sigx,sigy,f,dij,xi,yi,zi,xj,yj,zj

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

		!Coordinate e distanze
		xi=m_xyz(i,1)
		yi=m_xyz(i,2)
		zi=m_xyz(i,3)

		!Coordinate per la seconda sfere, comprese complicate traslazioni
		cell_row_do_nanb: DO icell=0,nacell,nacell
			cell_col_do_nanb: DO jcell=0,nbcell,nbcell

				sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
				sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2))) 

				xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
				yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
				zj=m_xyz(j,3)

				dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

				!Calcolo il coeff d'interazione
				f=(v_r(i)+v_r(j))/dij

				!Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
				IF (f>=fint) EXIT cell_row_do_nanb

			END DO cell_col_do_nanb
		END DO cell_row_do_nanb

		!Calcolo le coordinate sferiche relative
		xij=xi-xj
		yij=yi-yj
		zij =zi-zj

		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        dip_if: IF (i<=(ns-ndip)) THEN
        
        	!Chiamo la subroutine per riempire Aij e Aijsca
        	CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &  
                            & m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)
        ELSE
        
        	!Chiamo la subroutine per riempire solamente Aijsca, dal momento che sono nella zona dipoli
        	CALL fillblock_sparse_onlysca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij_sca(:,next),m_Bij_sca(:,next),&  
                            & m_jABij(:,next),m_iABij(:,next),error)
                            
        END IF dip_if

    END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_dip_per





!******************************************************************************
!8bis) SUBROUTINE fill_AB_sparse_dip:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_dip_nosca(ns,nstop,k,m_xyz,v_c0,v_aq,v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop                               			! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k                                          			! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      			! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         			! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij					  	! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij       			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

        !Calcolo le coordinate sferiche relative
        xij=m_xyz(i,1)-m_xyz(j,1)
        yij=m_xyz(i,2)-m_xyz(j,2)
        zij =m_xyz(i,3)-m_xyz(j,3)
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        !Chiamo la subroutine per riempire Dij
        CALL fillblock_sparse_dip_nosca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &  
                            & m_jABij(:,next),m_iABij(:,next),error)
                            

    END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_dip_nosca



!******************************************************************************
!8tris) SUBROUTINE fill_AB_sparse_sca:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_sca(ns,nstop,k,m_xyz,v_c0,v_aq,v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop                               			! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k                                          			! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      			! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         			! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               			! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca  	! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij       			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire Dij
		CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &  
				    & m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)


	END DO col_do
END DO row_do

END SUBROUTINE fill_AB_sparse_sca


!******************************************************************************
!8quater) SUBROUTINE fill_AB_sparse_sca:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_sca_fast(ns,nstop,k,m_xyz,v_c0,v_aq,v_bq,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				&  v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop							! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock				! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca		! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij				! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!I cycle col_do in I am in the lower triangular part of the matrix
		IF (i>j) CYCLE col_do

		!-----Diagnostics over the cycle------------------
		!WRITE(*,*) "AB,i,j",i,j
		!-------------------------------------------------

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire ABij
		CALL fillblock_sparse_sca_fast(nstop,v_c0,v_aq,v_bq,k*rij,m_index,m_Apmn,v_jABij_template,v_iABij_template, &
		m_Aij(:,next),m_Bij(:,next),m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

	END DO col_do
END DO row_do

END SUBROUTINE fill_AB_sparse_sca_fast



!******************************************************************************
!8quinties) SUBROUTINE fill_AB_sparse_sca:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_fast(ns,nstop,k,m_xyz,v_c0,v_aq,v_bq,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				&  v_jBlock,v_iBlock,m_Aij,m_Bij,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop							! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock				! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij				! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij				! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire ABij
		CALL fillblock_sparse_fast(nstop,v_c0,v_aq,v_bq,k*rij,m_index,m_Apmn,v_jABij_template,v_iABij_template, &
		m_Aij(:,next),m_Bij(:,next),m_jABij(:,next),m_iABij(:,next),error)

	END DO col_do
END DO row_do

END SUBROUTINE fill_AB_sparse_fast


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! PERIODIC PROBLEM SUBROUTINES
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!******************************************************************************
!1) SUBROUTINE nnb_sparse_per:calcolo il numero di blocchi diversi da zero. I blocchi
!              diagonali qui non sono inclusi
!******************************************************************************
SUBROUTINE nnb_sparse_per(ns,m_xyz,v_r,fint,nnb)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                      ! N sfere
REAL(dbl), INTENT(IN) :: fint                        ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r           ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz       ! Matrix posizione
INTEGER(lo), INTENT(OUT) :: nnb                    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,icell,jcell
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f,sigx,sigy

! Funzione vera e propria
nnb=0
row_do: DO i=1,ns
    col_do: DO j=1,ns

        IF (i==j) CYCLE col_do

        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        !Coordinate per la seconda sfere, comprese complicate traslazioni
        cell_row_do_nanb: DO icell=0,nacell,nacell
             cell_col_do_nanb: DO jcell=0,nbcell,nbcell
             
                    sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
                    sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2)))      

                    xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
                    yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
                    zj=m_xyz(j,3)

                    dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                    !Calcolo il coeff d'interazione
                    f=(v_r(i)+v_r(j))/dij

                    !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                    !WRITE(*,*) i,j,f,dij
                    IF (f>=fint) EXIT cell_row_do_nanb

              END DO cell_col_do_nanb
        END DO cell_row_do_nanb

        !Se la distanza e' piccola allora il blocco non e' zero
        nnb_if: IF (f>=fint) THEN
            nnb=nnb+1
        END IF nnb_if

    END DO col_do
END DO row_do

END SUBROUTINE nnb_sparse_per


!******************************************************************************
!2) SUBROUTINE fill_jBlock_iBlock_sparse_per:riempio i vettori indici riga e colonna
! per la disposizione dei blocchi di matrice.
!******************************************************************************
SUBROUTINE fill_jBlock_iBlock_sparse_per(ns,m_xyz,v_r,fint,v_jBlock,v_iBlock)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns                                 ! N sfere
REAL(dbl), INTENT(IN) :: fint                                   ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                      ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                  ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(OUT) :: v_jBlock,v_iBlock    ! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,next,icell,jcell
REAL(dbl) :: xi,yi,zi,xj,yj,zj,dij,f,sigx,sigy

! Funzione vera e propria
next=1
!Per il primo blocco e' sempre uno
v_iBlock(1)=1

row_do: DO i=1,ns
    col_do: DO j=1,ns

        IF (i==j) CYCLE col_do

        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        !Coordinate per la seconda sfere, comprese complicate traslazioni
        cell_row_do_nanb: DO icell=0,nacell,nacell
             cell_col_do_nanb: DO jcell=0,nbcell,nbcell
             
                    sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
                    sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2)))      

                    xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
                    yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
                    zj=m_xyz(j,3)

                    dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                    !Calcolo il coeff d'interazione
                    f=(v_r(i)+v_r(j))/dij

                    !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                    !WRITE(*,*) i,j,f,dij
                    IF (f>=fint) EXIT cell_row_do_nanb

              END DO cell_col_do_nanb
        END DO cell_row_do_nanb

!        !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
!        IF (f<fint) CYCLE col_do

        !Assegno i valori di colonna
        fint_if: IF (f>=fint) THEN
        	v_jBlock(next)=j
        	next=next+1
        END IF fint_if

    END DO col_do

    v_iBlock(i+1)=next

END DO row_do

END SUBROUTINE fill_jBlock_iBlock_sparse_per


!******************************************************************************
!3) SUBROUTINE fill_D_PHI_sparse_per:riempio tutta la struttura per i blocchi di 
! rotazione Dij e per Exp[phi_ij]
!******************************************************************************
SUBROUTINE fill_D_PHI_sparse_per(ns,nstop,nnb,m_xyz,fint,v_r,v_jBlock,v_iBlock,m_Dnkm,m_jDnkm,m_iDnkm,m_exphi)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop,nnb                           ! N sfere, N multipolar exp e N nnz blocks
REAL(dbl), INTENT(IN) :: fint                                       ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                          ! Vettore raggi
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         ! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Dnkm                    ! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jDnkm,m_iDnkm       ! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop,1:nnb), INTENT(OUT) :: m_exphi !Matrice Esponenziali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m,icell,jcell
REAL(dbl) :: xij,yij,zij,xi,yi,zi,dij,xj,yj,zj,rij,thetaij,phiij,sigx,sigy,f

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)
        
        !Coordinate e distanze
        xi=m_xyz(i,1)
        yi=m_xyz(i,2)
        zi=m_xyz(i,3)

        !Coordinate per la seconda sfere, comprese complicate traslazioni
        cell_row_do_nanb: DO icell=0,nacell,nacell
             cell_col_do_nanb: DO jcell=0,nbcell,nbcell
             
                    sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
                    sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2))) 
                    
                    xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
                    yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
                    zj=m_xyz(j,3)

                    dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                    !Calcolo il coeff d'interazione
                    f=(v_r(i)+v_r(j))/dij

                    !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                    IF (f>=fint) EXIT cell_row_do_nanb

              END DO cell_col_do_nanb
        END DO cell_row_do_nanb

        !Calcolo le coordinate sferiche relative
        xij=xi-xj
        yij=yi-yj
        zij =zi-zj
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        !Chiamo la subroutine per riempire Dij
        CALL fillblock_Dkmn(nstop,thetaij,m_Dnkm(:,next),m_jDnkm(:,next),m_iDnkm(:,next),error)

        !Riempio anche la colonna degli esponenziali
        phi_do: DO m=-nstop,nstop

            m_exphi(m,next)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

        END DO phi_do

    END DO col_do

END DO row_do

END SUBROUTINE fill_D_PHI_sparse_per



!******************************************************************************
!4) SUBROUTINE fill_AB_sparse_per_sca:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_per_sca(ns,nstop,k,fint,v_r,m_xyz,v_c0,v_aq,&
					  & v_bq,v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop								! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k										! vettore d'onda
REAL(dbl), INTENT(IN) :: fint										! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r							! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz							! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock					! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq						! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca		! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij				! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m,icell,jcell
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,dij,f,thetaij,phiij,sigx,sigy

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!Coordinate e distanze
		xi=m_xyz(i,1)
		yi=m_xyz(i,2)
		zi=m_xyz(i,3)

		!Coordinate per la seconda sfere, comprese complicate traslazioni
		cell_row_do_nanb: DO icell=0,nacell,nacell
			cell_col_do_nanb: DO jcell=0,nbcell,nbcell

				sigx=SIGN(1.0D0,(m_xyz(i,1)-m_xyz(j,1)))
				sigy=SIGN(1.0D0,(m_xyz(i,2)-m_xyz(j,2))) 

				xj=m_xyz(j,1)+sigx*REAL(icell,dbl)*ax+sigy*REAL(jcell,dbl)*bx
				yj=m_xyz(j,2)+sigy*REAL(jcell,dbl)*by
				zj=m_xyz(j,3)

				dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

				!Calcolo il coeff d'interazione
				f=(v_r(i)+v_r(j))/dij

				!Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
				IF (f>=fint) EXIT cell_row_do_nanb

			END DO cell_col_do_nanb
		END DO cell_row_do_nanb

		!Calcolo le coordinate sferiche relative
		xij=xi-xj
		yij=yi-yj
		zij =zi-zj

		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire Dij
		CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &  
						& m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

	END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_per_sca







! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! PERIODIC PROBLEM SUBROUTINES: SECONDA FORMALIZZAZIONE
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!******************************************************************************
!1) SUBROUTINE fillblock_expPhi:riempio un blocco ij dei coefficienti Dkmn  
!******************************************************************************
!SUBROUTINE fillblock_expPhi(nstop,phi,v_expPhi,v_expPhim,v_jexpPhi,v_iexpPhi,error)

!IMPLICIT NONE

!Dichiarazione dummy arguments
!INTEGER(lo), INTENT(IN) :: nstop                               ! Numero di espansioni multipolari
!REAL(dbl), INTENT(IN) :: phi                                     ! Angolo theta coordinate relative
!INTEGER(lo), INTENT(OUT) :: error                              ! Flag di errore
!REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_expPhi,v_expPhim       ! Vettore sparse di output valori Dkmn
!INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jexpPhi,v_iexpPhi        ! Vettori indici colonne e righe

!Dichiarazione variabili interne 
!INTEGER(lo) :: n,m      !Indici

!Subroutine vera e propria
!error=0

!Check di errore su theta
!phiin_check: IF ((phi>2.0D0*PI_D) .OR. (phi<-2.0D0*PI_D)) THEN
!                WRITE(*,*) "Phi e' al di fuori dei suoi bounds, la procedura fillblock_Dkmn si ferma"
!                error=1
!                RETURN
!END IF phiin_check

!Riempio le matrici sparse,siccome sono diagonali il processo e' oltremodo semplice
!n_do: DO i=1,nstop
!    m_do: DO m=-n,n

        !Riempio la matrice,niente fronzoli perche' e' molto semplice

!        v_expPhi(n*(n+1)+m) =  EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phi )
!        v_expPhim(n*(n+1)-m) =  (-1.0D0**(m))*v_expPhi(n*(n+1)+m)
!        v_jexpPhi(n*(n+1)+m) = n*(n+1)+m
!       v_iexpPhi(n*(n+1)+m) = n*(n+1)+m

        
!    END DO m_do
!END DO n_do

!Qui setto l'ultimo valore per v_iexpPhi, quello che dice dove comincerebbe la linea successiva
!v_iexpPhi(nstop*(nstop+2)+1) = nstop*(nstop+2)+1

!END SUBROUTINE fillblock_expPhi



!******************************************************************************
!2) SUBROUTINE fill_AB_per:subroutine per riempire  la matrice densa periodica,
!              quella super compressa con la bella idea!!!
!******************************************************************************
SUBROUTINE fill_AB_per(k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,m_a,m_b,m_AB)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k                                          ! vettore d'onda
REAL(dbl), INTENT(IN) :: fint                                       ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                          ! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               ! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b                 ! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB                   ! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr   ! Indici
INTEGER(lo) :: n,m,mu,nu,NNZab,NNZd,na,nb,blockside,matrixside    ! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2  ! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,xj0,yj0,dij,f,thetaij,phiij,sigx,sigy
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_Aij,v_Bij              ! Matrice blocchi Aij  Bij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij         ! Matrici indici Aij Bij
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Dnkm                      ! Matrice blocchi Dij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jDnkm,v_iDnkm         ! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop) :: v_exphi                     ! Matrice blocchi PHIij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_rhs                  ! Matrice Right Hand Side

! Funzione vera e propria
   
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!Lato di un blocco
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside


!-------------------------------------------------------------------------------------------------------------
!PRIMA PARTE DELLA SUBROUTINE,RIEMPIO, CON SOMMA TUTTI I POSTI A MENO DI UNA SIMMETRIA!!!
!-------------------------------------------------------------------------------------------------------------

!Alloco tutti i vettori necessari ai miei calcoli
ALLOCATE(v_Aij(1:NNZab),v_Bij(1:NNZab),v_jABij(1:NNZab),v_iABij(1:blockside+1),STAT=status_arr)
ALLOCATE(v_Dnkm(1:NNZd),v_jDnkm(1:NNZd),v_iDnkm(1:blockside+1),STAT=status_arr)
ALLOCATE(m_rhs(1:blockside,1:8),STAT=status_arr)

fill_row_block_do: DO icell=1,ns
    fill_col_block_do: DO jcell=icell,ns

        !Offset indici per il riempimento della matrice
        low_row=2*blockside*(icell-1)
        low_col=2*blockside*(jcell-1)

        fill_row_cell_do: DO na=0,nacell-1
            fill_col_cell_do: DO nb=0,nbcell-1

                !Se sono su un blocco diagonale, non sommo subito i fattori diagonali che rompono le simmetrie
                IF ( (icell==jcell) .AND. (na==0) .AND. (nb==0) ) CYCLE fill_col_cell_do

                !Coordinate e distanze
                xi=m_xyz(icell,1)
                yi=m_xyz(icell,2)
                zi=m_xyz(icell,3)

                xj0=m_xyz(jcell,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
                yj0=m_xyz(jcell,2) + REAL(nb,dbl)*by
                zj=m_xyz(jcell,3)

                !Coordinate per la seconda sfere, comprese complicate traslazioni
                cell_row_do_nanb: DO icell1=0,nacell,nacell
                    cell_col_do_nanb: DO jcell1=0,nbcell,nbcell
             
                        sigx=SIGN(1.0D0,(xi-xj0))
                        sigy=SIGN(1.0D0,(yi-yj0)) 

                        xj=xj0 + sigx*REAL(icell1,dbl)*ax + sigy*REAL(jcell1,dbl)*bx
                        yj=yj0 + sigy*REAL(jcell1,dbl)*by

                        dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                       !Calcolo il coeff d'interazione
                       f=(v_r(icell)+v_r(jcell))/dij

                       !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                       IF (f>=fint) EXIT cell_row_do_nanb

                   END DO cell_col_do_nanb
               END DO cell_row_do_nanb

			   IF (f<fint) CYCLE fill_col_cell_do

               !Calcolo le coordinate sferiche relative
                xij=xi-xj
                yij=yi-yj
                zij=zi-zj
        
               !Conversione coordinate da cartesiane a sferiche
               CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
               !Chiamo la subroutine per riempire Aij e Bij, Dij e PHIij
               CALL fillblock_sparse(nstop,v_c0,v_aq,v_bq,k*rij,v_Aij,v_Bij,v_jABij,v_iABij,error)
               CALL fillblock_Dkmn(nstop,thetaij,v_Dnkm,v_jDnkm,v_iDnkm,error)
               
               !Riempio anche la colonna degli esponenziali
                phi_do: DO m=-nstop,nstop

                    v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

                END DO phi_do

                !Qui alloco la matrice che uso per i vettori Right Hand Side
                m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)

               !Adesso comincio il ciclo su mu e nu, con tutte le moltiplicazioni del caso
               fill_nu_do: DO nu=1,nstop
                fill_mu_do: DO mu= -nu,nu

		    		m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)

                    !Qui estraggo una colonna (riga e' uguale), della mat diag e poi la densifico
                    m_rhs(nu*(nu+1)+mu,1) = v_exphi(mu)

                    !Faccio tutte le mie moltiplicazioni vettore matrice del caso, Qui rotazione
                    CALL amudz(blockside,m_rhs(:,1),m_rhs(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
                    
                    !Due traslazioni
                    CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,3),v_Aij,v_jABij,v_iABij)   !Per Aij
                    CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,4),v_Bij,v_jABij,v_iABij)   !Per Bij

                    !Aggiustamento della parita'
                    par_ndo: DO n=1,nstop
                        par_mdo: DO m=-n,n

                            m_rhs(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,3)
                            m_rhs(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,4)

                        END DO par_mdo 
                    END DO par_ndo


                    !Rotazioni per ciascuna delle due branche
                    CALL amudz(blockside,m_rhs(:,3),m_rhs(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte A
                    CALL amudz(blockside,m_rhs(:,4),m_rhs(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte B

                    !Estraggo le sottomatrici e poi moltiplico
                    !Aggiustamento della parita'
                    i=0
                    phi_ndo: DO n=nu,nstop
                        
                        diag_if: IF (n==nu) THEN

                            phi_mdo1: DO m=-n,-mu
                                i=n*(n+1)+m
                                m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
                                m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
                            END DO phi_mdo1

                        ELSE

                            phi_mdo2: DO m=-n,n
                                i=n*(n+1)+m
                                m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
                                m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
                            END DO phi_mdo2

                        END IF diag_if 

                    END DO phi_ndo

                    !Adesso riempio dai due vettori i pezzi di matrice che mi servono!!!
                    fill_ndo: DO n=nu,nstop
                        
                        fill_diag_if: IF (n==nu) THEN

                            !Qui riempio il blocco diagonale    
                            fill_mdo1: DO m=-n,-mu

                                iab=2*( n*(n+1)+m )
                                jab=2*( nu*(nu+1)+mu )
                                iab2=n*(n+1)+m
                                jab2=nu*(nu+1)+mu

                                !Indici per A e B
                                ia1=low_row+iab-1
                                ja1=low_col+jab-1
                                ib1=low_row+iab
                                jb1=low_col+jab-1

                                !Assegno i valori passando dai vettori che tengono a e b alla matrice
                                !Metto la somma per via della nuova teoria
                                m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
                                m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)

                            END DO fill_mdo1

                        ELSE

                            !Qui riempio il blocco fuori dalla diagonale diagonale
                            fill_mdo2: DO m=-n,n

                                iab=2*( n*(n+1)+m )
                                jab=2*( nu*(nu+1)+mu )
                                iab2=n*(n+1)+m
                                jab2=nu*(nu+1)+mu
                                
                                !Indici per A e B
                                ia1=low_row+iab-1
                                ja1=low_col+jab-1
                                ib1=low_row+iab
                                jb1=low_col+jab-1

                                !Assegno i valori passando dai vettori che tengono a e b alla matrice
                                !Metto la somma per via della nuova teoria
                                m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
                                m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)

                            END DO fill_mdo2



                        END IF fill_diag_if

                    END DO fill_ndo


                END DO fill_mu_do
               END DO fill_nu_do


            END DO fill_col_cell_do
        END DO fill_row_cell_do
    END DO fill_col_block_do
END DO fill_row_block_do

!-------------------------------------------------------------------------------------------------------------
!SECONDA PARTE DELLA SUBROUTINE: FATTE LE SOMME SFRUTTO LE SIMMETRIE
!-------------------------------------------------------------------------------------------------------------

sim_row_cell_do: DO icell=1,ns
    sim_col_cell_do: DO jcell=icell,ns
    
		!Offset indici per il riempimento della matrice
        low_row=2*blockside*(icell-1)
        low_col=2*blockside*(jcell-1)	

        sim_nu_do: DO nu=1,nstop
            sim_mu_do: DO mu=-nu,nu
                sim_n_do: DO n=nu,nstop

                    sim_if: IF (nu==n) THEN

                        sim_m_do1: DO m=-n,-mu
                            
                            !Indici per A e B originali
                            ia1=low_row+2*( n*(n+1)+m )-1
                            ja1=low_col+2*( nu*(nu+1)+mu )-1
                            ib1=ia1+1
                            jb1=ja1

                            !Indici per i primi doppioni
                            ia2=ia1+1
                            ja2=ja1+1
                            ib2=ib1-1
                            jb2=jb1+1
                            m_AB(ia2,ja2)=m_AB(ia1,ja1)
                            m_AB(ib2,jb2)=m_AB(ib1,jb1)

                            !Indici per prima simmetria da scambio di indici n=nu, m=-mu
                            ia3=low_row+2*( nu*(nu+1)-mu )-1
                            ja3=low_col+2*( n*(n+1)-m )-1
                            ib3=ia3+1
                            jb3=ja3
                            m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
                            m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)

                            !Indici per i doppioni della simmetria di indici
                            ia4=ia3+1
                            ja4=ja3+1
                            ib4=ib3-1
                            jb4=jb3+1
                            m_AB(ia4,ja4)=m_AB(ia3,ja3)
                            m_AB(ib4,jb4)=m_AB(ib3,jb3)

                            sim_diag_if1: IF (icell/=jcell) THEN

                                !Traslo tutti gli indici
                                ia5=ia1-low_row+low_col
                                ja5=ja1-low_col+low_row
                                ib5=ib1-low_row+low_col
                                jb5=jb1-low_col+low_row
                                m_AB(ia5,ja5)=(-1.0D0**(n+nu))*m_AB(ia1,ja1)
                                m_AB(ib5,jb5)=(-1.0D0**(n+nu+1))*m_AB(ib1,jb1)

                                ia6=ia2-low_row+low_col
                                ja6=ja2-low_col+low_row
                                ib6=ib2-low_row+low_col
                                jb6=jb2-low_col+low_row
                                m_AB(ia6,ja6)=(-1.0D0**(n+nu))*m_AB(ia2,ja2)
                                m_AB(ib6,jb6)=(-1.0D0**(n+nu+1))*m_AB(ib2,jb2)
                                
                                ia7=ia3-low_row+low_col
                                ja7=ja3-low_col+low_row
                                ib7=ib3-low_row+low_col
                                jb7=jb3-low_col+low_row
                                m_AB(ia7,ja7)=(-1.0D0**(n+nu))*m_AB(ia3,ja3)
                                m_AB(ib7,jb7)=(-1.0D0**(n+nu+1))*m_AB(ib3,jb3)

                                ia8=ia4-low_row+low_col
                                ja8=ja4-low_col+low_row
                                ib8=ib4-low_row+low_col
                                jb8=jb4-low_col+low_row
                                m_AB(ia8,ja8)=(-1.0D0**(n+nu))*m_AB(ia4,ja4)
                                m_AB(ib8,jb8)=(-1.0D0**(n+nu+1))*m_AB(ib4,jb4)

                            END IF sim_diag_if1

                        END DO sim_m_do1

                    ELSE

                        sim_m_do2: DO m=-n,n

                            !Indici per A e B originali
                            ia1=low_row+2*( n*(n+1)+m )-1
                            ja1=low_col+2*( nu*(nu+1)+mu )-1
                            ib1=ia1+1
                            jb1=ja1

                            !Indici per i primi doppioni
                            ia2=ia1+1
                            ja2=ja1+1
                            ib2=ib1-1
                            jb2=jb1+1
                            m_AB(ia2,ja2)=m_AB(ia1,ja1)
                            m_AB(ib2,jb2)=m_AB(ib1,jb1)
                                
!                                 WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
! 				101 FORMAT (8I4)
                            !Indici per prima simmetria da scambio di indici n=nu, m=-mu
                            ia3=low_row+2*( nu*(nu+1)-mu )-1
                            ja3=low_col+2*( n*(n+1)-m )-1
                            ib3=ia3+1
                            jb3=ja3
                            m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
                            m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
!                                 
!                                 WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
! 				102 FORMAT (8I4)
                            !Indici per i doppioni della simmetria di indici
                            ia4=ia3+1
                            ja4=ja3+1
                            ib4=ib3-1
                            jb4=jb3+1
                            m_AB(ia4,ja4)=m_AB(ia3,ja3)
                            m_AB(ib4,jb4)=m_AB(ib3,jb3)
!                                 
!                                 WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
! 				103 FORMAT (8I4)
                            sim_diag_if2: IF (icell/=jcell) THEN

                                !Traslo tutti gli indici
                                ia5=ia1-low_row+low_col
                                ja5=ja1-low_col+low_row
                                ib5=ib1-low_row+low_col
                                jb5=jb1-low_col+low_row
                                m_AB(ia5,ja5)=(-1.0D0**(n+nu))*m_AB(ia1,ja1)
                                m_AB(ib5,jb5)=(-1.0D0**(n+nu+1))*m_AB(ib1,jb1)

                                ia6=ia2-low_row+low_col
                                ja6=ja2-low_col+low_row
                                ib6=ib2-low_row+low_col
                                jb6=jb2-low_col+low_row
                                m_AB(ia6,ja6)=(-1.0D0**(n+nu))*m_AB(ia2,ja2)
                                m_AB(ib6,jb6)=(-1.0D0**(n+nu+1))*m_AB(ib2,jb2)
                                
                                ia7=ia3-low_row+low_col
                                ja7=ja3-low_col+low_row
                                ib7=ib3-low_row+low_col
                                jb7=jb3-low_col+low_row
                                m_AB(ia7,ja7)=(-1.0D0**(n+nu))*m_AB(ia3,ja3)
                                m_AB(ib7,jb7)=(-1.0D0**(n+nu+1))*m_AB(ib3,jb3)

                                ia8=ia4-low_row+low_col
                                ja8=ja4-low_col+low_row
                                ib8=ib4-low_row+low_col
                                jb8=jb4-low_col+low_row
                                m_AB(ia8,ja8)=(-1.0D0**(n+nu))*m_AB(ia4,ja4)
                                m_AB(ib8,jb8)=(-1.0D0**(n+nu+1))*m_AB(ib4,jb4)

                            END IF sim_diag_if2

                        END DO sim_m_do2

                    END IF sim_if

                END DO sim_n_do
            END DO sim_mu_do
        END DO sim_nu_do
    END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!ULTIMA PARTE DELLA SUBROUTINE, AGGIUNGO I TERMINI DIAGONALI!!!
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
    diagn_do: DO n=1,nstop
        diagm_do: DO m=-n,n
        
            j=j+1
            m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
            m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))


        END DO diagm_do
    END DO diagn_do
END DO diagi_do

                                
!                                 WRITE(*,*) "QUId"


END SUBROUTINE fill_AB_per




!******************************************************************************
!2.1) SUBROUTINE fill_AB_per_sca:subroutine per riempire  la matrice densa periodica,
!              quella super compressa con la bella idea, e pure lo scattering!!!
!******************************************************************************
SUBROUTINE fill_AB_per_sca_old(k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,m_a,m_b,m_AB,m_ABJ)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k                                          ! vettore d'onda
REAL(dbl), INTENT(IN) :: fint                                       ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                          ! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               ! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b                 ! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_ABJ             ! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr   ! Indici
INTEGER(lo) :: n,m,mu,nu,NNZab,NNZd,na,nb,blockside,matrixside    ! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2  ! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,xj0,yj0,dij,f,thetaij,phiij,sigx,sigy
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca  	! Matrice blocchi Aij  Bij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij         			! Matrici indici Aij Bij
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Dnkm                      			! Matrice blocchi Dij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jDnkm,v_iDnkm         			! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop) :: v_exphi                     			! Matrice blocchi PHIij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_rhs,m_rhs_sca           			! Matrice Right Hand Side

! Funzione vera e propria
   
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!Lato di un blocco
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside


!-------------------------------------------------------------------------------------------------------------
!PRIMA PARTE DELLA SUBROUTINE,RIEMPIO, CON SOMMA TUTTI I POSTI A MENO DI UNA SIMMETRIA!!!
!-------------------------------------------------------------------------------------------------------------

!Alloco tutti i vettori necessari ai miei calcoli
ALLOCATE(v_Aij(1:NNZab),v_Bij(1:NNZab),v_Aij_sca(1:NNZab),v_Bij_sca(1:NNZab), &
	   & v_jABij(1:NNZab),v_iABij(1:blockside+1),STAT=status_arr)
ALLOCATE(v_Dnkm(1:NNZd),v_jDnkm(1:NNZd),v_iDnkm(1:blockside+1),STAT=status_arr)
ALLOCATE(m_rhs(1:blockside,1:8),m_rhs_sca(1:blockside,1:8),STAT=status_arr)

fill_row_block_do: DO icell=1,ns
    fill_col_block_do: DO jcell=icell,ns

        !Offset indici per il riempimento della matrice
        low_row=2*blockside*(icell-1)
        low_col=2*blockside*(jcell-1)

        fill_row_cell_do: DO na=0,nacell-1
            fill_col_cell_do: DO nb=0,nbcell-1

                !Se sono su un blocco diagonale, non sommo subito i fattori diagonali che rompono le simmetrie
                IF ( (icell==jcell) .AND. (na==0) .AND. (nb==0) ) CYCLE fill_col_cell_do

                !Coordinate e distanze
                xi=m_xyz(icell,1)
                yi=m_xyz(icell,2)
                zi=m_xyz(icell,3)

                xj0=m_xyz(jcell,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
                yj0=m_xyz(jcell,2) + REAL(nb,dbl)*by
                zj=m_xyz(jcell,3)

                !Coordinate per la seconda sfere, comprese complicate traslazioni
                cell_row_do_nanb: DO icell1=0,nacell,nacell
                    cell_col_do_nanb: DO jcell1=0,nbcell,nbcell
             
                        sigx=SIGN(1.0D0,(xi-xj0))
                        sigy=SIGN(1.0D0,(yi-yj0)) 

                        xj=xj0 + sigx*REAL(icell1,dbl)*ax + sigy*REAL(jcell1,dbl)*bx
                        yj=yj0 + sigy*REAL(jcell1,dbl)*by

                        dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                       !Calcolo il coeff d'interazione
                       f=(v_r(icell)+v_r(jcell))/dij

                       !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                       IF (f>=fint) EXIT cell_row_do_nanb

                   END DO cell_col_do_nanb
               END DO cell_row_do_nanb

			   IF (f<fint) CYCLE fill_col_cell_do

               !Calcolo le coordinate sferiche relative
                xij=xi-xj
                yij=yi-yj
                zij=zi-zj
        
               !Conversione coordinate da cartesiane a sferiche
               CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
               !Chiamo la subroutine per riempire Aij e Bij, Dij e PHIij
               CALL fillblock_sparse_dip(nstop,v_c0,v_aq,v_bq,k*rij,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)
               CALL fillblock_Dkmn(nstop,thetaij,v_Dnkm,v_jDnkm,v_iDnkm,error)
               
               !Riempio anche la colonna degli esponenziali
                phi_do: DO m=-nstop,nstop

                    v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

                END DO phi_do

                !Qui alloco la matrice che uso per i vettori Right Hand Side
                m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
                m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

               !Adesso comincio il ciclo su mu e nu, con tutte le moltiplicazioni del caso
               fill_nu_do: DO nu=1,nstop
                fill_mu_do: DO mu= -nu,nu

				    m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
                    m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

                    !Qui estraggo una colonna (riga e' uguale), della mat diag e poi la densifico
                    m_rhs(nu*(nu+1)+mu,1) = v_exphi(mu)
                    m_rhs_sca(nu*(nu+1)+mu,1) = v_exphi(mu)

                    !Faccio tutte le mie moltiplicazioni vettore matrice del caso, Qui rotazione
                    CALL amudz(blockside,m_rhs(:,1),m_rhs(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
                    CALL amudz(blockside,m_rhs_sca(:,1),m_rhs_sca(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
                    
                    !Due traslazioni
                    CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,3),v_Aij,v_jABij,v_iABij)   !Per Aij
                    CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,4),v_Bij,v_jABij,v_iABij)   !Per Bij
                    CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,3),v_Aij_sca,v_jABij,v_iABij)   !Per AJij
                    CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,4),v_Bij_sca,v_jABij,v_iABij)   !Per BJij                    

                    !Aggiustamento della parita'
                    par_ndo: DO n=1,nstop
                        par_mdo: DO m=-n,n

                            m_rhs(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,3)
                            m_rhs(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,4)
                            m_rhs_sca(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs_sca(n*(n+1)+m,3)
                            m_rhs_sca(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs_sca(n*(n+1)+m,4)                            

                        END DO par_mdo 
                    END DO par_ndo


                    !Rotazioni per ciascuna delle due branche
                    CALL amudz(blockside,m_rhs(:,3),m_rhs(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte A
                    CALL amudz(blockside,m_rhs(:,4),m_rhs(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte B
                    CALL amudz(blockside,m_rhs_sca(:,3),m_rhs_sca(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte AJ
                    CALL amudz(blockside,m_rhs_sca(:,4),m_rhs_sca(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte BJ

                    !Estraggo le sottomatrici e poi moltiplico
                    !Aggiustamento della parita'
                    i=0
                    phi_ndo: DO n=nu,nstop
                        
                        diag_if: IF (n==nu) THEN

                            phi_mdo1: DO m=-n,-mu
                                i=n*(n+1)+m
                                m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
                                m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
                                m_rhs_sca(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,5)
                                m_rhs_sca(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,6)
                            END DO phi_mdo1

                        ELSE

                            phi_mdo2: DO m=-n,n
                                i=n*(n+1)+m
                                m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
                                m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
                                m_rhs_sca(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,5)
                                m_rhs_sca(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,6)
                            END DO phi_mdo2

                        END IF diag_if 

                    END DO phi_ndo

                    !Adesso riempio dai due vettori i pezzi di matrice che mi servono!!!
                    fill_ndo: DO n=nu,nstop
                        
                        fill_diag_if: IF (n==nu) THEN

                            !Qui riempio il blocco diagonale    
                            fill_mdo1: DO m=-n,-mu

                                iab=2*( n*(n+1)+m )
                                jab=2*( nu*(nu+1)+mu )
                                iab2=n*(n+1)+m
                                jab2=nu*(nu+1)+mu

                                !Indici per A e B
                                ia1=low_row+iab-1
                                ja1=low_col+jab-1
                                ib1=low_row+iab
                                jb1=low_col+jab-1

                                !Assegno i valori passando dai vettori che tengono a e b alla matrice
                                !Metto la somma per via della nuova teoria
                                m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
                                m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)
                                m_ABJ(ia1,ja1)=m_ABJ(ia1,ja1)+m_rhs_sca(iab2,7)
                                m_ABJ(ib1,jb1)=m_ABJ(ib1,jb1)+m_rhs_sca(iab2,8)
                                
                            END DO fill_mdo1

                        ELSE

                            !Qui riempio il blocco fuori dalla diagonale diagonale
                            fill_mdo2: DO m=-n,n

                                iab=2*( n*(n+1)+m )
                                jab=2*( nu*(nu+1)+mu )
                                iab2=n*(n+1)+m
                                jab2=nu*(nu+1)+mu
                                
                                !Indici per A e B
                                ia1=low_row+iab-1
                                ja1=low_col+jab-1
                                ib1=low_row+iab
                                jb1=low_col+jab-1

                                !Assegno i valori passando dai vettori che tengono a e b alla matrice
                                !Metto la somma per via della nuova teoria
                                m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
                                m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)
                                m_ABJ(ia1,ja1)=m_ABJ(ia1,ja1)+m_rhs_sca(iab2,7)
                                m_ABJ(ib1,jb1)=m_ABJ(ib1,jb1)+m_rhs_sca(iab2,8)

                            END DO fill_mdo2



                        END IF fill_diag_if

                    END DO fill_ndo


                END DO fill_mu_do
               END DO fill_nu_do


            END DO fill_col_cell_do
        END DO fill_row_cell_do
    END DO fill_col_block_do
END DO fill_row_block_do

!-------------------------------------------------------------------------------------------------------------
!SECONDA PARTE DELLA SUBROUTINE: FATTE LE SOMME SFRUTTO LE SIMMETRIE
!-------------------------------------------------------------------------------------------------------------

sim_row_cell_do: DO icell=1,ns
    sim_col_cell_do: DO jcell=icell,ns
    
		!Offset indici per il riempimento della matrice
        low_row=2*blockside*(icell-1)
        low_col=2*blockside*(jcell-1)

        sim_nu_do: DO nu=1,nstop
            sim_mu_do: DO mu=-nu,nu
                sim_n_do: DO n=nu,nstop

                    sim_if: IF (nu==n) THEN

                        sim_m_do1: DO m=-n,-mu
                            
                            !Indici per A e B originali
                            ia1=low_row+2*( n*(n+1)+m )-1
                            ja1=low_col+2*( nu*(nu+1)+mu )-1
                            ib1=ia1+1
                            jb1=ja1

                            !Indici per i primi doppioni
                            ia2=ia1+1
                            ja2=ja1+1
                            ib2=ib1-1
                            jb2=jb1+1
                            m_AB(ia2,ja2)=m_AB(ia1,ja1)
                            m_AB(ib2,jb2)=m_AB(ib1,jb1)
                            m_ABJ(ia2,ja2)=m_ABJ(ia1,ja1)
                            m_ABJ(ib2,jb2)=m_ABJ(ib1,jb1)

                            !Indici per prima simmetria da scambio di indici n=nu, m=-mu
                            ia3=low_row+2*( nu*(nu+1)-mu )-1
                            ja3=low_col+2*( n*(n+1)-m )-1
                            ib3=ia3+1
                            jb3=ja3
                            m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
                            m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
                            m_ABJ(ia3,ja3)=((-1.0D0)**(m+mu))*m_ABJ(ia1,ja1)
                            m_ABJ(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_ABJ(ib1,jb1)

                            !Indici per i doppioni della simmetria di indici
                            ia4=ia3+1
                            ja4=ja3+1
                            ib4=ib3-1
                            jb4=jb3+1
                            m_AB(ia4,ja4)=m_AB(ia3,ja3)
                            m_AB(ib4,jb4)=m_AB(ib3,jb3)
                            m_ABJ(ia4,ja4)=m_ABJ(ia3,ja3)
                            m_ABJ(ib4,jb4)=m_ABJ(ib3,jb3)

                            sim_diag_if1: IF (icell/=jcell) THEN

                                !Traslo tutti gli indici
                                ia5=ia1-low_row+low_col
                                ja5=ja1-low_col+low_row
                                ib5=ib1-low_row+low_col
                                jb5=jb1-low_col+low_row
                                m_AB(ia5,ja5)=(-1.0D0**(n+nu))*m_AB(ia1,ja1)
                                m_AB(ib5,jb5)=(-1.0D0**(n+nu+1))*m_AB(ib1,jb1)
                                m_ABJ(ia5,ja5)=(-1.0D0**(n+nu))*m_ABJ(ia1,ja1)
                                m_ABJ(ib5,jb5)=(-1.0D0**(n+nu+1))*m_ABJ(ib1,jb1)

                                ia6=ia2-low_row+low_col
                                ja6=ja2-low_col+low_row
                                ib6=ib2-low_row+low_col
                                jb6=jb2-low_col+low_row
                                m_AB(ia6,ja6)=(-1.0D0**(n+nu))*m_AB(ia2,ja2)
                                m_AB(ib6,jb6)=(-1.0D0**(n+nu+1))*m_AB(ib2,jb2)
                                m_ABJ(ia6,ja6)=(-1.0D0**(n+nu))*m_ABJ(ia2,ja2)
                                m_ABJ(ib6,jb6)=(-1.0D0**(n+nu+1))*m_ABJ(ib2,jb2)
                                
                                ia7=ia3-low_row+low_col
                                ja7=ja3-low_col+low_row
                                ib7=ib3-low_row+low_col
                                jb7=jb3-low_col+low_row
                                m_AB(ia7,ja7)=(-1.0D0**(n+nu))*m_AB(ia3,ja3)
                                m_AB(ib7,jb7)=(-1.0D0**(n+nu+1))*m_AB(ib3,jb3)
                                m_ABJ(ia7,ja7)=(-1.0D0**(n+nu))*m_ABJ(ia3,ja3)
                                m_ABJ(ib7,jb7)=(-1.0D0**(n+nu+1))*m_ABJ(ib3,jb3)

                                ia8=ia4-low_row+low_col
                                ja8=ja4-low_col+low_row
                                ib8=ib4-low_row+low_col
                                jb8=jb4-low_col+low_row
                                m_AB(ia8,ja8)=(-1.0D0**(n+nu))*m_AB(ia4,ja4)
                                m_AB(ib8,jb8)=(-1.0D0**(n+nu+1))*m_AB(ib4,jb4)
                                m_ABJ(ia8,ja8)=(-1.0D0**(n+nu))*m_ABJ(ia4,ja4)
                                m_ABJ(ib8,jb8)=(-1.0D0**(n+nu+1))*m_ABJ(ib4,jb4)

                            END IF sim_diag_if1

                        END DO sim_m_do1

                    ELSE

                        sim_m_do2: DO m=-n,n

                            !Indici per A e B originali
                            ia1=low_row+2*( n*(n+1)+m )-1
                            ja1=low_col+2*( nu*(nu+1)+mu )-1
                            ib1=ia1+1
                            jb1=ja1

                            !Indici per i primi doppioni
                            ia2=ia1+1
                            ja2=ja1+1
                            ib2=ib1-1
                            jb2=jb1+1
                            m_AB(ia2,ja2)=m_AB(ia1,ja1)
                            m_AB(ib2,jb2)=m_AB(ib1,jb1)
                            m_ABJ(ia2,ja2)=m_ABJ(ia1,ja1)
                            m_ABJ(ib2,jb2)=m_ABJ(ib1,jb1)
                                
!                                 WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
! 				101 FORMAT (8I4)
                            !Indici per prima simmetria da scambio di indici n=nu, m=-mu
                            ia3=low_row+2*( nu*(nu+1)-mu )-1
                            ja3=low_col+2*( n*(n+1)-m )-1
                            ib3=ia3+1
                            jb3=ja3
                            m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
                            m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
                            m_ABJ(ia3,ja3)=((-1.0D0)**(m+mu))*m_ABJ(ia1,ja1)
                            m_ABJ(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_ABJ(ib1,jb1)
!                                 
!                                 WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
! 				102 FORMAT (8I4)
                            !Indici per i doppioni della simmetria di indici
                            ia4=ia3+1
                            ja4=ja3+1
                            ib4=ib3-1
                            jb4=jb3+1
                            m_AB(ia4,ja4)=m_AB(ia3,ja3)
                            m_AB(ib4,jb4)=m_AB(ib3,jb3)
                            m_ABJ(ia4,ja4)=m_ABJ(ia3,ja3)
                            m_ABJ(ib4,jb4)=m_ABJ(ib3,jb3)
!                                 
!                                 WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
! 				103 FORMAT (8I4)
                            sim_diag_if2: IF (icell/=jcell) THEN

                                !Traslo tutti gli indici
                                ia5=ia1-low_row+low_col
                                ja5=ja1-low_col+low_row
                                ib5=ib1-low_row+low_col
                                jb5=jb1-low_col+low_row
                                m_AB(ia5,ja5)=(-1.0D0**(n+nu))*m_AB(ia1,ja1)
                                m_AB(ib5,jb5)=(-1.0D0**(n+nu+1))*m_AB(ib1,jb1)
                                m_ABJ(ia5,ja5)=(-1.0D0**(n+nu))*m_ABJ(ia1,ja1)
                                m_ABJ(ib5,jb5)=(-1.0D0**(n+nu+1))*m_ABJ(ib1,jb1)

                                ia6=ia2-low_row+low_col
                                ja6=ja2-low_col+low_row
                                ib6=ib2-low_row+low_col
                                jb6=jb2-low_col+low_row
                                m_AB(ia6,ja6)=(-1.0D0**(n+nu))*m_AB(ia2,ja2)
                                m_AB(ib6,jb6)=(-1.0D0**(n+nu+1))*m_AB(ib2,jb2)
                                m_ABJ(ia6,ja6)=(-1.0D0**(n+nu))*m_ABJ(ia2,ja2)
                                m_ABJ(ib6,jb6)=(-1.0D0**(n+nu+1))*m_ABJ(ib2,jb2)
                                
                                ia7=ia3-low_row+low_col
                                ja7=ja3-low_col+low_row
                                ib7=ib3-low_row+low_col
                                jb7=jb3-low_col+low_row
                                m_AB(ia7,ja7)=(-1.0D0**(n+nu))*m_AB(ia3,ja3)
                                m_AB(ib7,jb7)=(-1.0D0**(n+nu+1))*m_AB(ib3,jb3)
                                m_ABJ(ia7,ja7)=(-1.0D0**(n+nu))*m_ABJ(ia3,ja3)
                                m_ABJ(ib7,jb7)=(-1.0D0**(n+nu+1))*m_ABJ(ib3,jb3)

                                ia8=ia4-low_row+low_col
                                ja8=ja4-low_col+low_row
                                ib8=ib4-low_row+low_col
                                jb8=jb4-low_col+low_row
                                m_AB(ia8,ja8)=(-1.0D0**(n+nu))*m_AB(ia4,ja4)
                                m_AB(ib8,jb8)=(-1.0D0**(n+nu+1))*m_AB(ib4,jb4)
                                m_ABJ(ia8,ja8)=(-1.0D0**(n+nu))*m_ABJ(ia4,ja4)
                                m_ABJ(ib8,jb8)=(-1.0D0**(n+nu+1))*m_ABJ(ib4,jb4)

                            END IF sim_diag_if2

                        END DO sim_m_do2

                    END IF sim_if

                END DO sim_n_do
            END DO sim_mu_do
        END DO sim_nu_do
    END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!ULTIMA PARTE DELLA SUBROUTINE, AGGIUNGO I TERMINI DIAGONALI!!!
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
    diagn_do: DO n=1,nstop
        diagm_do: DO m=-n,n
        
            j=j+1
            m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
            m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))
            m_ABJ(2*j-1,2*j-1)=m_ABJ(2*j-1,2*j-1)+1.0D0
            m_ABJ(2*j,2*j)=m_ABJ(2*j,2*j)+1.0D0


        END DO diagm_do
    END DO diagn_do
END DO diagi_do

                                
!                                 WRITE(*,*) "QUId"


END SUBROUTINE fill_AB_per_sca_old



!******************************************************************************
!2.2) SUBROUTINE fillblock_sparse_per_sca:filling sparse block for periodic calc.
!******************************************************************************
SUBROUTINE fillblock_sparse_per_sca(nstop,v_c0,v_aq,v_bq,kr,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
				&    v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr								! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca		! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij				! Vettori sparse colonne e righe
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,q,p,low_nu,qmax			!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr					!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq			!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes,NNZab				!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca		!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h			!Vett comp fun hankel
INTEGER(lo) :: chunk=1					!Thread id, number of total threads and chunk size

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
NNZab=nstop*(1+2*nstop*(3+nstop))	!Computing NNZab,so i can have the cycle on the whole number of filled matrix values
NNZab=NNZab/3

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of hankel functions
ALLOCATE(v_h(0:(2*nstop+1)))
CALL hankel1_d_sub((2*nstop+1),kr,v_h,error)

!Checking the call for hankel
error_h_if: IF (error==1) THEN
		WRITE(*,10) 
		10 FORMAT ("Error in hankel1_d_sub called by fillblock_sparse_per_sc, STOPPING") 
		STOP
END IF error_h_if

!Making a copy of the block filling pattern
v_jABij=v_jABij_template
v_iABij=v_iABij_template

!This is a flattened do for the block filling, as in the standard routine, it is meant to be openMP friendly
flat_do: DO i=1,NNZab

	!----------------------------------------------------------------------------------------
	!Assigning all the needed indexes, which will be private, I assume, in the openMP version
	!----------------------------------------------------------------------------------------
	n=m_index(i,1)
	m=m_index(i,2)
	nu=m_index(i,3)
	qmax=m_index(i,4)
	pmin=m_index(i,5)
	pmax=m_index(i,6)
	low_aq=m_index(i,8)
	up_aq=m_index(i,9)
	low_bq=m_index(i,10)
	up_bq=m_index(i,11)

	!-----------
	!Calcolo Avt
	!-----------
	!Calcolo della sommatoria e di Avt
	sommaA=(0.0D0,0.0D0)
	sommaAsca=(0.0D0,0.0D0)
	sommaA_do:DO q=0,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
!		sommaA=sommaA+((0.0D0,1.0D0)**p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_h(p)
!		sommaAsca=sommaAsca+((0.0D0,1.0D0)**p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
		sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_h(p)
		sommaAsca=sommaAsca+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*REAL(v_h(p),dbl)
	END DO sommaA_do

	v_Aij(i)=v_c0(i)*sommaA
	v_Aij_sca(i)=v_c0(i)*sommaAsca

	!-----------
	!Calcolo Bvt
	!----------- 
	!Calcolo della sommatoria e di Bvt
	sommaB=(0.0D0,0.0D0)
	sommaBsca=(0.0D0,0.0D0)
	sommaB_do:DO q=1,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
!		sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
!		sommaBsca=sommaBsca+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
		sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_h(p+1)
		sommaBsca=sommaBsca+v_ic(p+1)*v_bq(low_bq+(q-1))*REAL(v_h(p+1),dbl)
	END DO sommaB_do

	v_Bij(i)=v_c0(i)*sommaB
	v_Bij_sca(i)=v_c0(i)*sommaBsca

END DO flat_do

END SUBROUTINE fillblock_sparse_per_sca



!******************************************************************************
!2.3) SUBROUTINE fillblock_dense_per_sca:filling a dense block for periodic structures
!******************************************************************************
SUBROUTINE fillblock_dense_per_sca(k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,&
				   m_index,m_Apmn,v_jABij_template,v_iABij_template,na,nb,nacell,nbcell,icell,jcell,&
				   f,m_block,m_block_sca)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda
REAL(dbl), INTENT(IN) :: fint								! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r						! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
INTEGER(lo), INTENT(IN) :: na,nb,nacell,nbcell,icell,jcell				! Vettore norm,gaunt e bq calling the present subroutine
REAL(dbl), INTENT(OUT) :: f								! Interaction factor
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_block,m_block_sca			! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell1,jcell1,status_arr					! Indici
INTEGER(lo) :: n,m,mu,nu,NNZab,NNZd,blockside,matrixside				! Indici
INTEGER(lo) :: ia1,ia2,iab,jab,iab2,jab2  ! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2
INTEGER(lo) :: ib1,ib2
INTEGER(lo) :: jb1,jb2
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,xj0,yj0,dij,thetaij,phiij,sigx,sigy
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca		! Matrice blocchi Aij  Bij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij				! Matrici indici Aij Bij
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Dnkm						! Matrice blocchi Dij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jDnkm,v_iDnkm				! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop) :: v_exphi					! Matrice blocchi PHIij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_rhs,m_rhs_sca				! Matrice Right Hand Side

! Subroutine
   
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!block and matrix side
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside

!Initializing outputs
m_block=(0.0D0,0.0D0)
m_block_sca=(0.0D0,0.0D0)


!Allocation of necessary vectors
ALLOCATE(v_Aij(1:NNZab),v_Bij(1:NNZab),v_Aij_sca(1:NNZab),v_Bij_sca(1:NNZab), &
	   & v_jABij(1:NNZab),v_iABij(1:blockside+1),STAT=status_arr)
ALLOCATE(v_Dnkm(1:NNZd),v_jDnkm(1:NNZd),v_iDnkm(1:blockside+1),STAT=status_arr)
ALLOCATE(m_rhs(1:blockside,1:8),m_rhs_sca(1:blockside,1:8),STAT=status_arr)

!Getting coordinates of the interacting spheres
xi=m_xyz(icell,1)
yi=m_xyz(icell,2)
zi=m_xyz(icell,3)

xj0=m_xyz(jcell,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
yj0=m_xyz(jcell,2) + REAL(nb,dbl)*by
zj=m_xyz(jcell,3)

!Cheking if I have to work directly with the sphere, or with its image
cell_row_do_nanb: DO icell1=0,nacell,nacell
	cell_col_do_nanb: DO jcell1=0,nbcell,nbcell

		sigx=SIGN(1.0D0,(xi-xj0))
		sigy=SIGN(1.0D0,(yi-yj0)) 

		xj=xj0 + sigx*REAL(icell1,dbl)*ax + sigy*REAL(jcell1,dbl)*bx
		yj=yj0 + sigy*REAL(jcell1,dbl)*by

		dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

		!Interaction coefficient
		f=(v_r(icell)+v_r(jcell))/dij

		!If spheres are interacting, I jump off the loops
		IF (f>=fint) EXIT cell_row_do_nanb

	END DO cell_col_do_nanb
END DO cell_row_do_nanb

!If spheres are not interacting, what should I do? I Stop the subroutine, and i cycle the loop
! that is calling this subroutine
exit_if: IF (f<fint) THEN

	DEALLOCATE(v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij)
	DEALLOCATE(v_Dnkm,v_jDnkm,v_iDnkm)
	DEALLOCATE(m_rhs,m_rhs_sca)

	RETURN

END IF exit_if

!Relative spherical coordinates
!WRITE(*,*)"xi,xj0,xj: ",xi,xj0,xj
!WRITE(*,*)"yi,yj: ",yi,yj
!WRITE(*,*)"zi,zj: ",zi,zj
xij=xi-xj
yij=yi-yj
zij=zi-zj

!Cartesian to spherical
CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

!Filling the relevant sparse blocks
!CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

CALL fillblock_sparse_per_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
			 v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

!--------------------------------------
!Diagnostics
!--------------------------------------
!IF ((na==1).AND.(nb==0)) THEN
!	WRITE(*,*) 'k,rij',k,rij
!	WRITE(*,*) 'v_jABij',v_jABij
!	WRITE(*,*)
!	WRITE(*,*) 'v_iABij',v_iABij
!	WRITE(*,*)
!	WRITE(*,*) 'v_Aij',v_Aij
!	WRITE(*,*)
!	WRITE(*,*) 'v_Bij',v_Bij
!END IF
!--------------------------------------

CALL fillblock_Dkmn(nstop,thetaij,v_Dnkm,v_jDnkm,v_iDnkm,error)
!--------------------------------------
!Diagnostics
!--------------------------------------
!IF ((na==1).AND.(nb==0)) THEN
!	WRITE(*,*)  'theta,phi',thetaij,phiij
!	WRITE(*,*) 'v_jDkmn',v_jDnkm
!	WRITE(*,*)
!	WRITE(*,*) 'v_iDkmn',v_iDnkm
!	WRITE(*,*)
!	WRITE(*,*) 'v_Dij',v_Dnkm
!END IF
!--------------------------------------
!Filling the phase matrix (diagonal, so using a vector)
phi_do: DO m=-nstop,nstop
	v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )
END DO phi_do

!Allocating column vectors that will be going to build the dense blocks of the matrix
!Each column of the matrix, will contain the updated version, multiplication after multiplication
m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Here we are going to build 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Beginning the loops over mu and mu
fill_nu_do: DO nu=1,nstop
	fill_mu_do: DO mu= -nu,nu

		m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
		m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

		!Qui estraggo una colonna (riga e' uguale), della mat diag e poi la densifico
		m_rhs(nu*(nu+1)+mu,1) = v_exphi(mu)
		m_rhs_sca(nu*(nu+1)+mu,1) = v_exphi(mu)

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
		!IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1)) THEN
		!	WRITE(*,*) "Multiplication 1",blockside
		!	do i=1,blockside
		!		WRITE(*,*) m_rhs_sca(i,1)
		!	end do
		!END IF
		!WRITE(*,*)
		!--------------------------------------

		!Faccio tutte le mie moltiplicazioni vettore matrice del caso, Qui rotazione
		CALL amudz(blockside,m_rhs(:,1),m_rhs(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
		CALL amudz(blockside,m_rhs_sca(:,1),m_rhs_sca(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
		!IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1)) THEN
		!	WRITE(*,*) "Multiplication 2",blockside
		!	do i=1,blockside
		!		WRITE(*,*) m_rhs_sca(i,2)
		!	end do
		!END IF
		!WRITE(*,*)
		!---------------------------------------

		!Due traslazioni
		CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,3),v_Aij,v_jABij,v_iABij)   !Per Aij
		CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,4),v_Bij,v_jABij,v_iABij)   !Per Bij
		CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,3),v_Aij_sca,v_jABij,v_iABij)   !Per AJij
		CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,4),v_Bij_sca,v_jABij,v_iABij)   !Per BJij

		!Aggiustamento della parita'
		par_ndo: DO n=1,nstop
			par_mdo: DO m=-n,n

				m_rhs(n*(n+1)+m,3)=v_oner(m)*m_rhs(n*(n+1)+m,3)
				m_rhs(n*(n+1)+m,4)=v_oner(m)*m_rhs(n*(n+1)+m,4)
				m_rhs_sca(n*(n+1)+m,3)=v_oner(m)*m_rhs_sca(n*(n+1)+m,3)
				m_rhs_sca(n*(n+1)+m,4)=v_oner(m)*m_rhs_sca(n*(n+1)+m,4)

			END DO par_mdo 
		END DO par_ndo

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
		!IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1)) THEN
		!	WRITE(*,*) "Multiplication 3,4",blockside
		!	do i=1,blockside
		!		WRITE(*,*) m_rhs_sca(i,3),m_rhs_sca(i,4)
		!	end do
		!END IF
		!WRITE(*,*)
		!--------------------------------------

		!Rotazioni per ciascuna delle due branche
		CALL amudz(blockside,m_rhs(:,3),m_rhs(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte A
		CALL amudz(blockside,m_rhs(:,4),m_rhs(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte B
		CALL amudz(blockside,m_rhs_sca(:,3),m_rhs_sca(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte AJ
		CALL amudz(blockside,m_rhs_sca(:,4),m_rhs_sca(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte BJ

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
		!IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1)) THEN
		!	WRITE(*,*) "Multiplication 5,6",blockside
		!	do i=1,blockside
		!		WRITE(*,*) m_rhs_sca(i,5),m_rhs_sca(i,6)
		!	end do
		!END IF
		!WRITE(*,*)
		!--------------------------------------

		!----------------------------------------------------------------------------------------------
		!Beginning nested loops on (n,m), i.e. I am cycling on the rows, i.e. i am filling a column
		!of the final matrix. I am not sure this is usefol or smart anyway
		!----------------------------------------------------------------------------------------------

		!Estraggo le sottomatrici e poi moltiplico
		!Aggiustamento della parita'
		i=0
		phi_ndo: DO n=nu,nstop

			diag_if: IF (n==nu) THEN

				phi_mdo1: DO m=-n,-mu
					i=n*(n+1)+m
					m_rhs(i,7)=v_exphi(-m)*v_oner(m)*m_rhs(i,5)
					m_rhs(i,8)=v_exphi(-m)*v_oner(m)*m_rhs(i,6)
					m_rhs_sca(i,7)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,5)
					m_rhs_sca(i,8)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,6)
				END DO phi_mdo1

			ELSE

				phi_mdo2: DO m=-n,n
					i=n*(n+1)+m
					m_rhs(i,7)=v_exphi(-m)*v_oner(m)*m_rhs(i,5)
					m_rhs(i,8)=v_exphi(-m)*v_oner(m)*m_rhs(i,6)
					m_rhs_sca(i,7)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,5)
					m_rhs_sca(i,8)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,6)
				END DO phi_mdo2

			END IF diag_if 

			!--------------------------------------
			!Diagnostics
			!--------------------------------------
			!IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1)) THEN
			!	WRITE(*,*) "Multiplication 7,8",blockside
			!	do i=1,blockside
			!		WRITE(*,*) m_rhs_sca(i,7),m_rhs_sca(i,8)
			!	end do
			!END IF
			!WRITE(*,*)
			!---------------------------------------

		END DO phi_ndo

		!Filling the dense block (then i complete with the simmetry)
		fill_ndo: DO n=nu,nstop

			fill_diag_if: IF (n==nu) THEN

				!Diagonal block n==nu
				fill_mdo1: DO m=-n,-mu

					iab=2*( n*(n+1)+m )
					jab=2*( nu*(nu+1)+mu )
					iab2=n*(n+1)+m
					jab2=nu*(nu+1)+mu

					!AB indexes: just the top row of the small block of four, then simmetry
					!AB indexes not translated for the previous block, given the input
					ia1=iab-1
					ja1=jab-1
					ib1=iab
					jb1=jab-1

!					WRITE(*,*) 'm,n,nu,ia1,ib1',m,n,nu,ia1,ib1

					!No sum just one block
					m_block(ia1,ja1)=m_rhs(iab2,7)
					m_block(ib1,jb1)=m_rhs(iab2,8)
					m_block_sca(ia1,ja1)=m_rhs_sca(iab2,7)
					m_block_sca(ib1,jb1)=m_rhs_sca(iab2,8)

				END DO fill_mdo1

			ELSE

				!Off diagonal block: n not nu
				fill_mdo2: DO m=-n,n

					iab=2*( n*(n+1)+m )
					jab=2*( nu*(nu+1)+mu )
					iab2=n*(n+1)+m
					jab2=nu*(nu+1)+mu

					!AB indexes: just the top row of the small block of four, then simmetry
					!AB indexes not translated for the previous block, given the input
					ia1=iab-1
					ja1=jab-1
					ib1=iab
					jb1=jab-1

!					WRITE(*,*) 'm,n,nu,ia1,ib1',m,n,nu,ia1,ib1

					!No sum just one block
					m_block(ia1,ja1)=m_rhs(iab2,7)
					m_block(ib1,jb1)=m_rhs(iab2,8)
					m_block_sca(ia1,ja1)=m_rhs_sca(iab2,7)
					m_block_sca(ib1,jb1)=m_rhs_sca(iab2,8)

				END DO fill_mdo2

			END IF fill_diag_if

		END DO fill_ndo

	END DO fill_mu_do
END DO fill_nu_do

!--------------------------------------
!Diagnostics
!--------------------------------------
!IF ((na==1).AND.(nb==0)) THEN
!	WRITE(*,*) blockside
!	do i=1,blockside
!		WRITE(*,1111) m_block_sca(i,:)
!		1111 FORMAT (12ES14.6)
!	end do

!	WRITE(*,*)
!	do i=1,matrixside
!		WRITE(*,1112) m_block(i,:)
!		1112 FORMAT (12ES14.6)
!	end do
!END IF
!-------------------------------------

END SUBROUTINE fillblock_dense_per_sca










!******************************************************************************
!2.3.1) SUBROUTINE fillblock_dense_per_sca_asym:filling a dense block for periodic structures
!******************************************************************************
SUBROUTINE fillblock_dense_per_sca_asym(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
				   m_index,m_Apmn,v_jABij_template,v_iABij_template,na,nb,icell,jcell,&
				   m_block,m_block_sca)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r						! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
INTEGER(lo), INTENT(IN) :: na,nb,icell,jcell						! Vettore norm,gaunt e bq calling the present subroutine
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_block,m_block_sca			! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell1,jcell1,status_arr					! Indici
INTEGER(lo) :: n,m,mu,nu,NNZab,NNZd,blockside,matrixside				! Indici
INTEGER(lo) :: ia1,ia2,iab,jab,iab2,jab2  ! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2
INTEGER(lo) :: ib1,ib2
INTEGER(lo) :: jb1,jb2
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,xj0,yj0,dij,thetaij,phiij,sigx,sigy
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca		! Matrice blocchi Aij  Bij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij				! Matrici indici Aij Bij
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Dnkm						! Matrice blocchi Dij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jDnkm,v_iDnkm				! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop) :: v_exphi					! Matrice blocchi PHIij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_rhs,m_rhs_sca				! Matrice Right Hand Side

! Subroutine
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!block and matrix side
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside

!Initializing outputs
m_block=(0.0D0,0.0D0)
m_block_sca=(0.0D0,0.0D0)


!Allocation of necessary vectors
ALLOCATE(v_Aij(1:NNZab),v_Bij(1:NNZab),v_Aij_sca(1:NNZab),v_Bij_sca(1:NNZab), &
	   & v_jABij(1:NNZab),v_iABij(1:blockside+1),STAT=status_arr)
ALLOCATE(v_Dnkm(1:NNZd),v_jDnkm(1:NNZd),v_iDnkm(1:blockside+1),STAT=status_arr)
ALLOCATE(m_rhs(1:blockside,1:8),m_rhs_sca(1:blockside,1:8),STAT=status_arr)

!Getting coordinates of the interacting spheres
xi=m_xyz(icell,1)
yi=m_xyz(icell,2)
zi=m_xyz(icell,3)

xj=m_xyz(jcell,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
yj=m_xyz(jcell,2) + REAL(nb,dbl)*by
zj=m_xyz(jcell,3)


!Relative spherical coordinates
!WRITE(*,*)"xi,xj0,xj: ",xi,xj0,xj
!WRITE(*,*)"yi,yj: ",yi,yj
!WRITE(*,*)"zi,zj: ",zi,zj
xij=xi-xj
yij=yi-yj
zij=zi-zj

!Cartesian to spherical
CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

!Filling the relevant sparse blocks
!CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

CALL fillblock_sparse_per_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_index,m_Apmn,v_jABij_template,v_iABij_template,&
			 v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

!--------------------------------------
!Diagnostics
!--------------------------------------
!IF ((icell==1).AND.(jcell==2)) THEN
!	WRITE(*,*) 'icell,jcell,na,rij',icell,jcell,na,rij
!	WRITE(*,*) 'v_jABij',v_jABij
!	WRITE(*,*)
!	WRITE(*,*) 'v_iABij',v_iABij
!	WRITE(*,*)
!	WRITE(*,*) 'v_Aij',v_Aij
!	WRITE(*,*)
!	WRITE(*,*) 'v_Bij',v_Bij
!END IF
!--------------------------------------

CALL fillblock_Dkmn(nstop,thetaij,v_Dnkm,v_jDnkm,v_iDnkm,error)
!--------------------------------------
!Diagnostics
!--------------------------------------
!IF ((icell==1).AND.(jcell==2)) THEN
!	WRITE(*,*)  'theta,phi',thetaij,phiij
!	WRITE(*,*) 'v_jDkmn',v_jDnkm
!	WRITE(*,*)
!	WRITE(*,*) 'v_iDkmn',v_iDnkm
!	WRITE(*,*)
!	WRITE(*,*) 'v_Dij',v_Dnkm
!END IF
!--------------------------------------
!Filling the phase matrix (diagonal, so using a vector)
phi_do: DO m=-nstop,nstop
	v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )
END DO phi_do

!Allocating column vectors that will be going to build the dense blocks of the matrix
!Each column of the matrix, will contain the updated version, multiplication after multiplication
m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Here we are going to build 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Beginning the loops over mu and mu
fill_nu_do: DO nu=1,nstop
	fill_mu_do: DO mu= -nu,nu

		m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
		m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

		!Qui estraggo una colonna (riga e' uguale), della mat diag e poi la densifico
		m_rhs(nu*(nu+1)+mu,1) = v_exphi(mu)
		m_rhs_sca(nu*(nu+1)+mu,1) = v_exphi(mu)

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
!		IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1).AND.(icell==1).AND.(jcell==2)) THEN
!			WRITE(*,*) "Multiplication 1",blockside
!			do i=1,blockside
!				WRITE(*,*) m_rhs_sca(i,1)
!			end do
!			do i=1,blockside
!				WRITE(*,*) m_rhs(i,1)
!			end do
!			WRITE(*,*)
!		END IF
		!--------------------------------------

		!Faccio tutte le mie moltiplicazioni vettore matrice del caso, Qui rotazione
		CALL amudz(blockside,m_rhs(:,1),m_rhs(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
		CALL amudz(blockside,m_rhs_sca(:,1),m_rhs_sca(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
!		IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1).AND.(icell==1).AND.(jcell==2)) THEN
!			WRITE(*,*) "Multiplication 2",blockside
!			do i=1,blockside
!				WRITE(*,*) m_rhs_sca(i,2)
!			end do
!			do i=1,blockside
!				WRITE(*,*) m_rhs(i,2)
!			end do
!			WRITE(*,*)
!		END IF
		!---------------------------------------

		!Due traslazioni
		CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,3),v_Aij,v_jABij,v_iABij)   !Per Aij
		CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,4),v_Bij,v_jABij,v_iABij)   !Per Bij
		CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,3),v_Aij_sca,v_jABij,v_iABij)   !Per AJij
		CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,4),v_Bij_sca,v_jABij,v_iABij)   !Per BJij

		!Aggiustamento della parita'
		par_ndo: DO n=1,nstop
			par_mdo: DO m=-n,n

				m_rhs(n*(n+1)+m,3)=v_oner(m)*m_rhs(n*(n+1)+m,3)
				m_rhs(n*(n+1)+m,4)=v_oner(m)*m_rhs(n*(n+1)+m,4)
				m_rhs_sca(n*(n+1)+m,3)=v_oner(m)*m_rhs_sca(n*(n+1)+m,3)
				m_rhs_sca(n*(n+1)+m,4)=v_oner(m)*m_rhs_sca(n*(n+1)+m,4)

			END DO par_mdo 
		END DO par_ndo

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
!		IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1).AND.(icell==1).AND.(jcell==2)) THEN
!			WRITE(*,*) "Multiplication 3,4",blockside
!			do i=1,blockside
!				WRITE(*,*) m_rhs_sca(i,3),m_rhs_sca(i,4)
!			end do
!			do i=1,blockside
!				WRITE(*,*) m_rhs(i,3),m_rhs(i,4)
!			end do
!			WRITE(*,*)
!		END IF
		!--------------------------------------

		!Rotazioni per ciascuna delle due branche
		CALL amudz(blockside,m_rhs(:,3),m_rhs(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte A
		CALL amudz(blockside,m_rhs(:,4),m_rhs(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte B
		CALL amudz(blockside,m_rhs_sca(:,3),m_rhs_sca(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte AJ
		CALL amudz(blockside,m_rhs_sca(:,4),m_rhs_sca(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte BJ

		!--------------------------------------
		!Diagnostics
		!--------------------------------------
!		IF ((na==1).AND.(nb==0).AND.(nu==1).AND.(mu==-1).AND.(icell==1).AND.(jcell==2)) THEN
!			WRITE(*,*) "Multiplication 5,6",blockside
!			do i=1,blockside
!				WRITE(*,*) m_rhs_sca(i,5),m_rhs_sca(i,6)
!			end do
!			do i=1,blockside
!				WRITE(*,*) m_rhs(i,5),m_rhs(i,6)
!			end do
!			WRITE(*,*)
!		END IF
		!--------------------------------------

		!----------------------------------------------------------------------------------------------
		!Beginning nested loops on (n,m), i.e. I am cycling on the rows, i.e. i am filling a column
		!of the final matrix. I am not sure this is usefol or smart anyway
		!----------------------------------------------------------------------------------------------

		!Estraggo le sottomatrici e poi moltiplico
		!Aggiustamento della parita'
		i=0
		phi_ndo: DO n=nu,nstop

			diag_if: IF (n==nu) THEN

				phi_mdo1: DO m=-n,-mu
					i=n*(n+1)+m
					m_rhs(i,7)=v_exphi(-m)*v_oner(m)*m_rhs(i,5)
					m_rhs(i,8)=v_exphi(-m)*v_oner(m)*m_rhs(i,6)
					m_rhs_sca(i,7)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,5)
					m_rhs_sca(i,8)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,6)
				END DO phi_mdo1

			ELSE

				phi_mdo2: DO m=-n,n
					i=n*(n+1)+m
					m_rhs(i,7)=v_exphi(-m)*v_oner(m)*m_rhs(i,5)
					m_rhs(i,8)=v_exphi(-m)*v_oner(m)*m_rhs(i,6)
					m_rhs_sca(i,7)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,5)
					m_rhs_sca(i,8)=v_exphi(-m)*v_oner(m)*m_rhs_sca(i,6)
				END DO phi_mdo2

			END IF diag_if 

			!--------------------------------------
			!Diagnostics
			!--------------------------------------
!			IF ((nu==1).AND.(mu==-1).AND.(icell==1).AND.(jcell==2)) THEN
!				WRITE(*,*) "Multiplication 7,8",blockside
!				do i=1,blockside
!					WRITE(*,*) m_rhs_sca(i,7),m_rhs_sca(i,8)
!				end do
!				do i=1,blockside
!					WRITE(*,*) m_rhs(i,7),m_rhs(i,8)
!				end do
!				WRITE(*,*)
!			END IF
			!---------------------------------------

		END DO phi_ndo

		!Filling the dense block (then i complete with the simmetry)
		fill_ndo: DO n=nu,nstop

			fill_diag_if: IF (n==nu) THEN

				!Diagonal block n==nu
				fill_mdo1: DO m=-n,-mu

					iab=2*( n*(n+1)+m )
					jab=2*( nu*(nu+1)+mu )
					iab2=n*(n+1)+m
					jab2=nu*(nu+1)+mu

					!AB indexes: just the top row of the small block of four, then simmetry
					!AB indexes not translated for the previous block, given the input
					ia1=iab-1
					ja1=jab-1
					ib1=iab
					jb1=jab-1

!					WRITE(*,*) 'm,n,nu,ia1,ib1',m,n,nu,ia1,ib1

					!No sum just one block
					m_block(ia1,ja1)=m_rhs(iab2,7)
					m_block(ib1,jb1)=m_rhs(iab2,8)
					m_block_sca(ia1,ja1)=m_rhs_sca(iab2,7)
					m_block_sca(ib1,jb1)=m_rhs_sca(iab2,8)

				END DO fill_mdo1

			ELSE

				!Off diagonal block: n not nu
				fill_mdo2: DO m=-n,n

					iab=2*( n*(n+1)+m )
					jab=2*( nu*(nu+1)+mu )
					iab2=n*(n+1)+m
					jab2=nu*(nu+1)+mu

					!AB indexes: just the top row of the small block of four, then simmetry
					!AB indexes not translated for the previous block, given the input
					ia1=iab-1
					ja1=jab-1
					ib1=iab
					jb1=jab-1

!					WRITE(*,*) 'm,n,nu,ia1,ib1',m,n,nu,ia1,ib1

					!No sum just one block
					m_block(ia1,ja1)=m_rhs(iab2,7)
					m_block(ib1,jb1)=m_rhs(iab2,8)
					m_block_sca(ia1,ja1)=m_rhs_sca(iab2,7)
					m_block_sca(ib1,jb1)=m_rhs_sca(iab2,8)

				END DO fill_mdo2

			END IF fill_diag_if

		END DO fill_ndo

	END DO fill_mu_do
END DO fill_nu_do

!--------------------------------------
!Diagnostics
!--------------------------------------
!IF ((na==1).AND.(nb==0).AND.(icell==1).AND.(jcell==2)) THEN
!	WRITE(*,*)"na,nb", na,nb
!	WRITE(*,*)"icell,jcell", icell,jcell
!	WRITE(*,*)"Blockside", blockside
!	do i=1,2*blockside
!!		WRITE(*,*)'row',i
!		WRITE(*,1111) m_block_sca(i,:)
!		1111 FORMAT (12ES14.6)
!	end do

!	WRITE(*,*)
!	do i=1,matrixside
!		WRITE(*,1112) m_block(i,:)
!		1112 FORMAT (12ES14.6)
!	end do
!END IF
!-------------------------------------

END SUBROUTINE fillblock_dense_per_sca_asym



!******************************************************************************
!2.4) SUBROUTINE fill_AB_per_sca_test:subroutine per riempire  la matrice densa periodica,
!              quella super compressa con la bella idea, e pure lo scattering!!!
!******************************************************************************
SUBROUTINE fill_AB_per_sca(k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,m_a,m_b,m_AB,m_ABJ)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k                                          ! vettore d'onda
REAL(dbl), INTENT(IN) :: fint                                       ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                          ! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               ! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b                 ! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_ABJ             ! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr   ! Indici
INTEGER(lo) :: n,m,mu,nu,NNZab,NNZd,na,nb,blockside,matrixside    ! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2  ! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,xj0,yj0,dij,f,thetaij,phiij,sigx,sigy
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca  	! Matrice blocchi Aij  Bij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij         			! Matrici indici Aij Bij
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Dnkm                      			! Matrice blocchi Dij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jDnkm,v_iDnkm         			! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop) :: v_exphi                     			! Matrice blocchi PHIij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_rhs,m_rhs_sca           			! Matrice Right Hand Side

! Funzione vera e propria
   
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!Lato di un blocco
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside


!-------------------------------------------------------------------------------------------------------------
!PRIMA PARTE DELLA SUBROUTINE,RIEMPIO, CON SOMMA TUTTI I POSTI A MENO DI UNA SIMMETRIA!!!
!-------------------------------------------------------------------------------------------------------------

!Alloco tutti i vettori necessari ai miei calcoli
ALLOCATE(v_Aij(1:NNZab),v_Bij(1:NNZab),v_Aij_sca(1:NNZab),v_Bij_sca(1:NNZab), &
	   & v_jABij(1:NNZab),v_iABij(1:blockside+1),STAT=status_arr)
ALLOCATE(v_Dnkm(1:NNZd),v_jDnkm(1:NNZd),v_iDnkm(1:blockside+1),STAT=status_arr)
ALLOCATE(m_rhs(1:blockside,1:8),m_rhs_sca(1:blockside,1:8),STAT=status_arr)

fill_row_block_do: DO icell=1,ns
	fill_col_block_do: DO jcell=1,ns

		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)
		low_col=2*blockside*(jcell-1)

		fill_row_cell_do: DO na=0,nacell-1
			fill_col_cell_do: DO nb=0,nbcell-1

				!Se sono su un blocco diagonale, non sommo subito i fattori diagonali che rompono le simmetrie
				IF ( (icell==jcell) .AND. (na==0) .AND. (nb==0) ) CYCLE fill_col_cell_do

				!Coordinate e distanze
				xi=m_xyz(icell,1)
				yi=m_xyz(icell,2)
				zi=m_xyz(icell,3)

				xj0=m_xyz(jcell,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
				yj0=m_xyz(jcell,2) + REAL(nb,dbl)*by
				zj=m_xyz(jcell,3)

				!Coordinate per la seconda sfere, comprese complicate traslazioni
				cell_row_do_nanb: DO icell1=0,nacell,nacell
					cell_col_do_nanb: DO jcell1=0,nbcell,nbcell

						sigx=SIGN(1.0D0,(xi-xj0))
						sigy=SIGN(1.0D0,(yi-yj0)) 

						xj=xj0 + sigx*REAL(icell1,dbl)*ax + sigy*REAL(jcell1,dbl)*bx
						yj=yj0 + sigy*REAL(jcell1,dbl)*by

						dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

						!Calcolo il coeff d'interazione
						f=(v_r(icell)+v_r(jcell))/dij

						!If spheres are interacting, I jump off the loops
						IF (f>=fint) EXIT cell_row_do_nanb

					END DO cell_col_do_nanb
				END DO cell_row_do_nanb

				!If spheres are not interacting, i cycle the loop
				IF (f<fint) CYCLE fill_col_cell_do

				!Calcolo le coordinate sferiche relative
				!WRITE(*,*)"xi,xj0,xj: ",xi,xj0,xj
				!WRITE(*,*)"yi,yj: ",yi,yj
				!WRITE(*,*)"zi,zj: ",zi,zj
				xij=xi-xj
				yij=yi-yj
				zij=zi-zj

				!Conversione coordinate da cartesiane a sferiche
				CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

				!Chiamo la subroutine per riempire Aij e Bij, Dij e PHIij
				CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,&
				v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)
				CALL fillblock_Dkmn(nstop,thetaij,v_Dnkm,v_jDnkm,v_iDnkm,error)

				!Riempio anche la colonna degli esponenziali
				phi_do: DO m=-nstop,nstop

					v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

				END DO phi_do

				!Qui alloco la matrice che uso per i vettori Right Hand Side
				m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
				m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

				!Adesso comincio il ciclo su mu e nu, con tutte le moltiplicazioni del caso
				fill_nu_do: DO nu=1,nstop
					fill_mu_do: DO mu= -nu,nu

						m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
						m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

						!Qui estraggo una colonna (riga e' uguale), della mat diag e poi la densifico
						m_rhs(nu*(nu+1)+mu,1) = v_exphi(mu)
						m_rhs_sca(nu*(nu+1)+mu,1) = v_exphi(mu)

						!Faccio tutte le mie moltiplicazioni vettore matrice del caso, Qui rotazione
						CALL amudz(blockside,m_rhs(:,1),m_rhs(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
						CALL amudz(blockside,m_rhs_sca(:,1),m_rhs_sca(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm

						!Due traslazioni
						CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,3),v_Aij,v_jABij,v_iABij)   !Per Aij
						CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,4),v_Bij,v_jABij,v_iABij)   !Per Bij
						CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,3),v_Aij_sca,v_jABij,v_iABij)   !Per AJij
						CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,4),v_Bij_sca,v_jABij,v_iABij)   !Per BJij

						!Aggiustamento della parita'
						par_ndo: DO n=1,nstop
							par_mdo: DO m=-n,n

								m_rhs(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,3)
								m_rhs(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,4)
								m_rhs_sca(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs_sca(n*(n+1)+m,3)
								m_rhs_sca(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs_sca(n*(n+1)+m,4)

							END DO par_mdo 
						END DO par_ndo


						!Rotazioni per ciascuna delle due branche
						CALL amudz(blockside,m_rhs(:,3),m_rhs(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte A
						CALL amudz(blockside,m_rhs(:,4),m_rhs(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte B
						CALL amudz(blockside,m_rhs_sca(:,3),m_rhs_sca(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte AJ
						CALL amudz(blockside,m_rhs_sca(:,4),m_rhs_sca(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte BJ

						!Estraggo le sottomatrici e poi moltiplico
						!Aggiustamento della parita'
						i=0
						phi_ndo: DO n=nu,nstop

						diag_if: IF (n==nu) THEN

								phi_mdo1: DO m=-n,-mu
									i=n*(n+1)+m
									m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
									m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
									m_rhs_sca(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,5)
									m_rhs_sca(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,6)
								END DO phi_mdo1

							ELSE

								phi_mdo2: DO m=-n,n
									i=n*(n+1)+m
									m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
									m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
									m_rhs_sca(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,5)
									m_rhs_sca(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,6)
								END DO phi_mdo2

							END IF diag_if 

						END DO phi_ndo

						!Adesso riempio dai due vettori i pezzi di matrice che mi servono!!!
						fill_ndo: DO n=nu,nstop

							fill_diag_if: IF (n==nu) THEN

								!Qui riempio il blocco diagonale    
								fill_mdo1: DO m=-n,-mu

									iab=2*( n*(n+1)+m )
									jab=2*( nu*(nu+1)+mu )
									iab2=n*(n+1)+m
									jab2=nu*(nu+1)+mu

									!Indici per A e B
									ia1=low_row+iab-1
									ja1=low_col+jab-1
									ib1=low_row+iab
									jb1=low_col+jab-1

									!Assegno i valori passando dai vettori che tengono a e b alla matrice
									!Metto la somma per via della nuova teoria
									m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
									m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)
									m_ABJ(ia1,ja1)=m_ABJ(ia1,ja1)+m_rhs_sca(iab2,7)
									m_ABJ(ib1,jb1)=m_ABJ(ib1,jb1)+m_rhs_sca(iab2,8)
		
								END DO fill_mdo1

							ELSE

								!Qui riempio il blocco fuori dalla diagonale diagonale
								fill_mdo2: DO m=-n,n

									iab=2*( n*(n+1)+m )
									jab=2*( nu*(nu+1)+mu )
									iab2=n*(n+1)+m
									jab2=nu*(nu+1)+mu
		
									!Indici per A e B
									ia1=low_row+iab-1
									ja1=low_col+jab-1
									ib1=low_row+iab
									jb1=low_col+jab-1

									!Assegno i valori passando dai vettori che tengono a e b alla matrice
									!Metto la somma per via della nuova teoria
									m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
									m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)
									m_ABJ(ia1,ja1)=m_ABJ(ia1,ja1)+m_rhs_sca(iab2,7)
									m_ABJ(ib1,jb1)=m_ABJ(ib1,jb1)+m_rhs_sca(iab2,8)

								END DO fill_mdo2

							END IF fill_diag_if

						END DO fill_ndo

					END DO fill_mu_do
				END DO fill_nu_do

			END DO fill_col_cell_do
		END DO fill_row_cell_do

	END DO fill_col_block_do
END DO fill_row_block_do

!-------------------------------------------------------------------------------------------------------------
!SECONDA PARTE DELLA SUBROUTINE: FATTE LE SOMME SFRUTTO LE SIMMETRIE
!-------------------------------------------------------------------------------------------------------------

sim_row_cell_do: DO icell=1,ns
	sim_col_cell_do: DO jcell=1,ns
 
		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)
		low_col=2*blockside*(jcell-1)

		sim_nu_do: DO nu=1,nstop
			sim_mu_do: DO mu=-nu,nu
				sim_n_do: DO n=nu,nstop

					sim_if: IF (nu==n) THEN

						sim_m_do1: DO m=-n,-mu

							!Indici per A e B originali
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Indici per i primi doppioni
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_ABJ(ia2,ja2)=m_ABJ(ia1,ja1)
							m_ABJ(ib2,jb2)=m_ABJ(ib1,jb1)

							!Indici per prima simmetria da scambio di indici n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_ABJ(ia3,ja3)=((-1.0D0)**(m+mu))*m_ABJ(ia1,ja1)
							m_ABJ(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_ABJ(ib1,jb1)

							!Indici per i doppioni della simmetria di indici
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_ABJ(ia4,ja4)=m_ABJ(ia3,ja3)
							m_ABJ(ib4,jb4)=m_ABJ(ib3,jb3)

						END DO sim_m_do1

					ELSE

						sim_m_do2: DO m=-n,n

							!Indici per A e B originali
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Indici per i primi doppioni
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_ABJ(ia2,ja2)=m_ABJ(ia1,ja1)
							m_ABJ(ib2,jb2)=m_ABJ(ib1,jb1)

							!WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
							!101 FORMAT (8I4)
							!Indici per prima simmetria da scambio di indici n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_ABJ(ia3,ja3)=((-1.0D0)**(m+mu))*m_ABJ(ia1,ja1)
							m_ABJ(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_ABJ(ib1,jb1)

							!WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
							!102 FORMAT (8I4)
							!Indici per i doppioni della simmetria di indici
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_ABJ(ia4,ja4)=m_ABJ(ia3,ja3)
							m_ABJ(ib4,jb4)=m_ABJ(ib3,jb3)

							!WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
							!103 FORMAT (8I4)


						END DO sim_m_do2

					END IF sim_if

				END DO sim_n_do
			END DO sim_mu_do
		END DO sim_nu_do
	END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!ULTIMA PARTE DELLA SUBROUTINE, AGGIUNGO I TERMINI DIAGONALI!!!
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
	diagn_do: DO n=1,nstop
		diagm_do: DO m=-n,n

			j=j+1
			m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
			m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))
			m_ABJ(2*j-1,2*j-1)=m_ABJ(2*j-1,2*j-1)+1.0D0
			m_ABJ(2*j,2*j)=m_ABJ(2*j,2*j)+1.0D0

		END DO diagm_do
	END DO diagn_do
END DO diagi_do

                                
!                                 WRITE(*,*) "QUId"


END SUBROUTINE fill_AB_per_sca









!******************************************************************************
!2.5) SUBROUTINE fill_AB_per_sca_fast:filling AB with fast routines openMP friendly
!******************************************************************************
SUBROUTINE fill_AB_per_sca_fast(k,m_fint,v_r,m_xyz,v_c0,v_aq,v_bq,&
				m_index,m_Apmn,v_jABij_template,v_iABij_template,m_a,m_b,m_AB,m_AB_sca)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_fint						! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r						! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b					! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_AB_sca				! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr			! Indici
INTEGER(lo) :: n,m,mu,nu,na,nb,blockside,ncore,omp_get_num_procs,chunk	! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2		! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col,high_row, high_col,low_i,high_i
INTEGER(lo), DIMENSION(1:nacell*nbcell,1:2) :: m_nanb
REAL(dbl) :: fint
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_f						!Interaction factor vector
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:,:) :: m_block,m_block_sca			! Temporary block matrix storage

! Subroutine

!Self commenting :-)
blockside=nstop*(nstop+2)
ncore=1!omp_get_num_procs()

!Filling the linearized index vector
i=0
na_do: DO na=0,nacell-1
	nb_do: DO nb=0,nbcell-1
		i=i+1
		m_nanb(i,1)=na
		m_nanb(i,2)=nb
	END DO nb_do
END DO na_do

!Allocating temporary vectors according to the number of cores
ALLOCATE(m_block(1:2*blockside,1:2*blockside,1:ncore),m_block_sca(1:2*blockside,1:2*blockside,1:ncore),STAT=status_arr)
m_block=0.0D0
m_block_sca=0.0D0
ALLOCATE(v_f(1:ncore),STAT=status_arr)

!-------------------------------------------------------------------------------------------------------------
!First Part: filling direct values, what remains is filled by simmetry
!-------------------------------------------------------------------------------------------------------------
fill_row_block_do: DO icell=1,ns
	fill_col_block_do: DO jcell=1,ns

		!Getting the appropriate interaction coefficient for truncation
		fint=m_fint(icell,jcell)

		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)+1
		low_col=2*blockside*(jcell-1)+1
		high_row=2*blockside*(icell)
		high_col=2*blockside*(jcell)

		!Bounds for the unrolled loop
		low_i=1
		high_i=ncore

		!Here I unroll the loop so as to preserve the openMP implementation
		unroll_do: DO

			!$OMP PARALLEL PRIVATE(i)
			chunk=1		!Choosing the chunk size
			!$OMP DO SCHEDULE(guided,chunk)
			linear_do: DO i=low_i,high_i

!				WRITE(*,*) 'i:',i

				!Diagonal factors at the end to preserve simmetries
				cell_if: IF ( (icell==jcell) .AND. (i==1)) THEN

!					WRITE(*,*) 'qui1,i,low_i,high_i,i-low_i+1: ',i,low_i,high_i,i-low_i+1

					m_block(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)
					m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)

!					WRITE(*,*) 'qui2'

				ELSE

!					WRITE(*,*) 'qui3'

				!Subroutine to fill a whole dense block
				CALL fillblock_dense_per_sca(k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,&
							     m_index,m_Apmn,v_jABij_template,v_iABij_template,&
							     m_nanb(i,1),m_nanb(i,2),nacell,nbcell,icell,jcell,v_f(i-low_i+1),&
							     m_block(:,:,i-low_i+1),m_block_sca(:,:,i-low_i+1))

!					WRITE(*,*) 'qui4'

				END IF cell_if

			END DO linear_do
			!$OMP END DO NOWAIT
			!$OMP END PARALLEL

			!Here i do the sum to fill the compact AB and AB_sca matrix
			fill_do: DO i=low_i,high_i

				!Not updating for irrelevant blocks
!				IF (((i==1) .AND. icell==jcell) .OR. (v_f(i-low_i+1)<fint)) CYCLE fill_do

				!Filling the hankel block
				m_AB(low_row:high_row,low_col:high_col)=&
				m_AB(low_row:high_row,low_col:high_col)+m_block(1:2*blockside,1:2*blockside,i-low_i+1)

				!Filling the bessel block
				m_AB_sca(low_row:high_row,low_col:high_col)=&
				m_AB_sca(low_row:high_row,low_col:high_col)+m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)

			END DO fill_do

			!If the job is finished, I leave the loop
			IF (high_i==nacell*nbcell) EXIT unroll_do

			!Refresh the indexes for the openMP loop
			low_i=high_i+1
			high_i=high_i+ncore
			high_if: IF (high_i>=(nacell*nbcell)) THEN
				high_i=(nacell*nbcell)
			END IF high_if

		END DO unroll_do

	END DO fill_col_block_do
END DO fill_row_block_do

!Freeing some memory
DEALLOCATE(m_block,m_block_sca)




!-------------------------------------------------------------------------------------------------------------
! Filling the remainder of the matrix exploiting symmetries
!-------------------------------------------------------------------------------------------------------------
sim_row_cell_do: DO icell=1,ns
	sim_col_cell_do: DO jcell=1,ns
 
		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)
		low_col=2*blockside*(jcell-1)

		!-------------------------------------------------------------------------------------
		!Filling the (icell,jcell) block exploiting the symmetries
		!-------------------------------------------------------------------------------------
		sim_nu_do: DO nu=1,nstop
			sim_mu_do: DO mu=-nu,nu
				sim_n_do: DO n=nu,nstop


					!-------------------------------------------------------------
					!Filling the diagonal block inside the block
					!-------------------------------------------------------------
					sim_if: IF (nu==n) THEN

						sim_m_do1: DO m=-n,-mu

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

						END DO sim_m_do1

					!-------------------------------------------------------------
					!Filling the off diagonal block inside the block
					!-------------------------------------------------------------
					ELSE

						sim_m_do2: DO m=-n,n

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
							!101 FORMAT (8I4)
							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
							!102 FORMAT (8I4)
							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

							!WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
							!103 FORMAT (8I4)


						END DO sim_m_do2

					END IF sim_if

				END DO sim_n_do
			END DO sim_mu_do
		END DO sim_nu_do

		!--------------------------------------------------------
		! Symmetry cannot be exploited if icell=jcell, of course
		!--------------------------------------------------------
		IF (icell==jcell) CYCLE sim_col_cell_do

		!------------------------------------------------------------------------------------------------
		!Filling the (jcell,icell) block exploiting the symmetry Aij=(-1)^(n+nu)Aji,Bij=(-1)^(n+nu+1)Bji
		!------------------------------------------------------------------------------------------------
		ij_n_do: DO n=1,nstop
			ij_m_do: DO m=-n,n
				ij_nu_do: DO nu=1,nstop
					ij_mu_do: DO mu=-nu,nu

						!Original index (ia1,ja1) (ib1,jb1), for the already filled part
						ia1=low_row+2*( n*(n+1)+m )-1
						ja1=low_col+2*( nu*(nu+1)+mu )-1
						ib1=ia1+1
						jb1=ja1

						!Index (ia2,ja2) (ib2,jb2) for the the mapping (i,j)->(j,i), low_col and low_row are exchanged
						ia2=low_col+2*( n*(n+1)+m )-1
						ja2=low_row+2*( nu*(nu+1)+mu )-1
						ib2=ia2+1
						jb2=ja2

						!Index (ia3,ja3) (ib3,jb3) for the second simmetry, i.e. the small block of four for (j,i)
						ia3=ia2+1
						ja3=ja2+1
						ib3=ib2-1
						jb3=jb2+1

						!Mapping the values from (i,j) to (j,i)
						m_AB(ia2,ja2)=v_oner(n+nu)*m_AB(ia1,ja1)
						m_AB(ib2,jb2)=v_oner(n+nu+1)*m_AB(ib1,jb1)
						m_AB_sca(ia2,ja2)=v_oner(n+nu)*m_AB_sca(ia1,ja1)
						m_AB_sca(ib2,jb2)=v_oner(n+nu+1)*m_AB_sca(ib1,jb1)

						!Filling with the simmetry on the block of four
						m_AB(ia3,ja3)=m_AB(ia2,ja2)
						m_AB(ib3,jb3)=m_AB(ib2,jb2)
						m_AB_sca(ia3,ja3)=m_AB_sca(ia2,ja2)
						m_AB_sca(ib3,jb3)=m_AB_sca(ib2,jb2)

					END DO ij_mu_do
				END DO ij_nu_do
			END DO ij_m_do
		END DO ij_n_do


	END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!Adding the diagonal terms
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
	diagn_do: DO n=1,nstop
		diagm_do: DO m=-n,n

			j=j+1
			m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
			m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))
			m_AB_sca(2*j-1,2*j-1)=m_AB_sca(2*j-1,2*j-1)+1.0D0
			m_AB_sca(2*j,2*j)=m_AB_sca(2*j,2*j)+1.0D0

		END DO diagm_do
	END DO diagn_do
END DO diagi_do

END SUBROUTINE fill_AB_per_sca_fast




!******************************************************************************
!2.6) SUBROUTINE fill_AB_per_sca_fast:filling AB with fast routines openMP friendly
!******************************************************************************
SUBROUTINE fill_AB_per_sca_fast_asym(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
				m_index,m_Apmn,v_jABij_template,v_iABij_template,m_a,m_b,m_AB,m_AB_sca)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r						! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b					! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_AB_sca				! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr			! Indici
INTEGER(lo) :: n,m,mu,nu,na,nb,blockside,ncore,omp_get_num_procs,chunk	! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2		! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col,high_row, high_col,low_i,high_i
INTEGER(lo), DIMENSION(1:(namax-namin+1)*(nbmax-nbmin+1),1:2) :: m_nanb
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:,:) :: m_block,m_block_sca			! Temporary block matrix storage

! Subroutine

!Self commenting :-)
blockside=nstop*(nstop+2)
ncore=1!omp_get_num_procs()

!Filling the linearized index vector
i=0
na_do: DO na=namin,namax
	nb_do: DO nb=nbmin,nbmax
		i=i+1
		m_nanb(i,1)=na
		m_nanb(i,2)=nb
	END DO nb_do
END DO na_do

!Allocating temporary vectors according to the number of cores
ALLOCATE(m_block(1:2*blockside,1:2*blockside,1:ncore),m_block_sca(1:2*blockside,1:2*blockside,1:ncore),STAT=status_arr)
m_block=0.0D0
m_block_sca=0.0D0

!-------------------------------------------------------------------------------------------------------------
!First Part: filling direct values, what remains is filled by simmetry
!-------------------------------------------------------------------------------------------------------------
fill_row_block_do: DO icell=1,ns
	fill_col_block_do: DO jcell=icell,ns

		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)+1
		low_col=2*blockside*(jcell-1)+1
		high_row=2*blockside*(icell)
		high_col=2*blockside*(jcell)

		!Bounds for the unrolled loop
		low_i=1
		high_i=ncore

		!Here I unroll the loop so as to preserve the openMP implementation
		unroll_do: DO

			!$OMP PARALLEL PRIVATE(i)
			chunk=1		!Choosing the chunk size
			!$OMP DO SCHEDULE(guided,chunk)
			linear_do: DO i=low_i,high_i

!				WRITE(*,*) 'i,na,nb:',i,m_nanb(i,1),m_nanb(i,2)

				!Diagonal factors at the end to preserve simmetries
				cell_if: IF ( (icell==jcell) .AND. (m_nanb(i,1)==0) .AND. (m_nanb(i,2)==0) ) THEN

!					WRITE(*,*) 'qui1,i,low_i,high_i,i-low_i+1: ',i,low_i,high_i,i-low_i+1

					m_block(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)
					m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)

!					WRITE(*,*) 'qui2'

				ELSE

!					WRITE(*,*) 'qui3'

				!Subroutine to fill a whole dense block
				CALL fillblock_dense_per_sca_asym(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
							     m_index,m_Apmn,v_jABij_template,v_iABij_template,&
							     m_nanb(i,1),m_nanb(i,2),icell,jcell,&
							     m_block(:,:,i-low_i+1),m_block_sca(:,:,i-low_i+1))

!					WRITE(*,*) 'qui4'

				END IF cell_if

			END DO linear_do
			!$OMP END DO NOWAIT
			!$OMP END PARALLEL

			!Here i do the sum to fill the compact AB and AB_sca matrix
			fill_do: DO i=low_i,high_i

				!Filling the hankel block
				m_AB(low_row:high_row,low_col:high_col)=&
				m_AB(low_row:high_row,low_col:high_col)+m_block(1:2*blockside,1:2*blockside,i-low_i+1)

				!Filling the bessel block
				m_AB_sca(low_row:high_row,low_col:high_col)=&
				m_AB_sca(low_row:high_row,low_col:high_col)+m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)

			END DO fill_do

			!If the job is finished, I leave the loop
			IF (high_i==(namax-namin+1)*(nbmax-nbmin+1)) EXIT unroll_do

			!Refresh the indexes for the openMP loop
			low_i=high_i+1
			high_i=high_i+ncore
			high_if: IF (high_i>=((namax-namin+1)*(nbmax-nbmin+1))) THEN
				high_i=((namax-namin+1)*(nbmax-nbmin+1))
			END IF high_if

		END DO unroll_do

	END DO fill_col_block_do
END DO fill_row_block_do

!Freeing some memory
DEALLOCATE(m_block,m_block_sca)




!-------------------------------------------------------------------------------------------------------------
! Filling the remainder of the matrix exploiting symmetries
!-------------------------------------------------------------------------------------------------------------
sim_row_cell_do: DO icell=1,ns
	sim_col_cell_do: DO jcell=icell,ns
 
		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)
		low_col=2*blockside*(jcell-1)

		!-------------------------------------------------------------------------------------
		!Filling the (icell,jcell) block exploiting the symmetries
		!-------------------------------------------------------------------------------------
		sim_nu_do: DO nu=1,nstop
			sim_mu_do: DO mu=-nu,nu
				sim_n_do: DO n=nu,nstop


					!-------------------------------------------------------------
					!Filling the diagonal block inside the block
					!-------------------------------------------------------------
					sim_if: IF (nu==n) THEN

						sim_m_do1: DO m=-n,-mu

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

						END DO sim_m_do1

					!-------------------------------------------------------------
					!Filling the off diagonal block inside the block
					!-------------------------------------------------------------
					ELSE

						sim_m_do2: DO m=-n,n

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
							!101 FORMAT (8I4)
							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
							!102 FORMAT (8I4)
							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

							!WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
							!103 FORMAT (8I4)


						END DO sim_m_do2

					END IF sim_if

				END DO sim_n_do
			END DO sim_mu_do
		END DO sim_nu_do

		!--------------------------------------------------------
		! Symmetry cannot be exploited if icell=jcell, of course
		!--------------------------------------------------------
		IF (icell==jcell) CYCLE sim_col_cell_do

		!------------------------------------------------------------------------------------------------
		!Filling the (jcell,icell) block exploiting the symmetry Aij=(-1)^(n+nu)Aji,Bij=(-1)^(n+nu+1)Bji
		!------------------------------------------------------------------------------------------------
		ij_n_do: DO n=1,nstop
			ij_m_do: DO m=-n,n
				ij_nu_do: DO nu=1,nstop
					ij_mu_do: DO mu=-nu,nu

						!Original index (ia1,ja1) (ib1,jb1), for the already filled part
						ia1=low_row+2*( n*(n+1)+m )-1
						ja1=low_col+2*( nu*(nu+1)+mu )-1
						ib1=ia1+1
						jb1=ja1

						!Index (ia2,ja2) (ib2,jb2) for the the mapping (i,j)->(j,i), low_col and low_row are exchanged
						ia2=low_col+2*( n*(n+1)+m )-1
						ja2=low_row+2*( nu*(nu+1)+mu )-1
						ib2=ia2+1
						jb2=ja2

						!Index (ia3,ja3) (ib3,jb3) for the second simmetry, i.e. the small block of four for (j,i)
						ia3=ia2+1
						ja3=ja2+1
						ib3=ib2-1
						jb3=jb2+1

						!Mapping the values from (i,j) to (j,i)
						m_AB(ia2,ja2)=v_oner(n+nu)*m_AB(ia1,ja1)
						m_AB(ib2,jb2)=v_oner(n+nu+1)*m_AB(ib1,jb1)
						m_AB_sca(ia2,ja2)=v_oner(n+nu)*m_AB_sca(ia1,ja1)
						m_AB_sca(ib2,jb2)=v_oner(n+nu+1)*m_AB_sca(ib1,jb1)

						!Filling with the simmetry on the block of four
						m_AB(ia3,ja3)=m_AB(ia2,ja2)
						m_AB(ib3,jb3)=m_AB(ib2,jb2)
						m_AB_sca(ia3,ja3)=m_AB_sca(ia2,ja2)
						m_AB_sca(ib3,jb3)=m_AB_sca(ib2,jb2)

					END DO ij_mu_do
				END DO ij_nu_do
			END DO ij_m_do
		END DO ij_n_do


	END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!Adding the diagonal terms
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
	diagn_do: DO n=1,nstop
		diagm_do: DO m=-n,n

			j=j+1
			m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
			m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))
			m_AB_sca(2*j-1,2*j-1)=m_AB_sca(2*j-1,2*j-1)+1.0D0
			m_AB_sca(2*j,2*j)=m_AB_sca(2*j,2*j)+1.0D0

		END DO diagm_do
	END DO diagn_do
END DO diagi_do

END SUBROUTINE fill_AB_per_sca_fast_asym




!******************************************************************************
!2.7) SUBROUTINE fill_AB_per_sca_fast:filling AB with fast routines openMP friendly
!******************************************************************************
SUBROUTINE fill_AB_per_sca_fast_asym_dip(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
				m_index,m_Apmn,v_jABij_template,v_iABij_template,m_a,m_b,m_AB,m_AB_sca)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r						! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b					! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_AB_sca				! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr			! Indici
INTEGER(lo) :: n,m,mu,nu,na,nb,blockside,ncore,omp_get_num_procs,chunk	! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2		! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col,high_row, high_col,low_i,high_i
INTEGER(lo), DIMENSION(1:(namax-namin+1)*(nbmax-nbmin+1),1:2) :: m_nanb
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:,:) :: m_block,m_block_sca			! Temporary block matrix storage

! Subroutine

!Self commenting :-)
blockside=nstop*(nstop+2)
ncore=1!omp_get_num_procs()

!Filling the linearized index vector
i=0
na_do: DO na=namin,namax
	nb_do: DO nb=nbmin,nbmax
		i=i+1
		m_nanb(i,1)=na
		m_nanb(i,2)=nb
	END DO nb_do
END DO na_do

!Allocating temporary vectors according to the number of cores
ALLOCATE(m_block(1:2*blockside,1:2*blockside,1:ncore),m_block_sca(1:2*blockside,1:2*blockside,1:ncore),STAT=status_arr)
m_block=0.0D0
m_block_sca=0.0D0

!-------------------------------------------------------------------------------------------------------------
!First Part: filling direct values, what remains is filled by simmetry
!-------------------------------------------------------------------------------------------------------------
fill_row_block_do: DO icell=1,ns
	fill_col_block_do: DO jcell=icell,ns

		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)+1
		low_col=2*blockside*(jcell-1)+1
		high_row=2*blockside*(icell)
		high_col=2*blockside*(jcell)

		!Bounds for the unrolled loop
		low_i=1
		high_i=ncore

		!Here I unroll the loop so as to preserve the openMP implementation
		unroll_do: DO

			!$OMP PARALLEL PRIVATE(i)
			chunk=1		!Choosing the chunk size
			!$OMP DO SCHEDULE(guided,chunk)
			linear_do: DO i=low_i,high_i

!				WRITE(*,*) 'i,na,nb:',i,m_nanb(i,1),m_nanb(i,2)

				!Diagonal factors at the end to preserve simmetries
				cell_if: IF ( (icell==jcell) .AND. (m_nanb(i,1)==0) .AND. (m_nanb(i,2)==0) ) THEN

!					WRITE(*,*) 'qui1,i,low_i,high_i,i-low_i+1: ',i,low_i,high_i,i-low_i+1

					m_block(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)
					m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)

!					WRITE(*,*) 'qui2'

				ELSE

!					WRITE(*,*) 'qui3'

				!Subroutine to fill a whole dense block
				CALL fillblock_dense_per_sca_asym(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
							     m_index,m_Apmn,v_jABij_template,v_iABij_template,&
							     m_nanb(i,1),m_nanb(i,2),icell,jcell,&
							     m_block(:,:,i-low_i+1),m_block_sca(:,:,i-low_i+1))

!					WRITE(*,*) 'qui4'

				END IF cell_if

			END DO linear_do
			!$OMP END DO NOWAIT
			!$OMP END PARALLEL

			!Here i do the sum to fill the compact AB and AB_sca matrix
			fill_do: DO i=low_i,high_i

				!Filling the hankel block
				m_AB(low_row:high_row,low_col:high_col)=&
				m_AB(low_row:high_row,low_col:high_col)+m_block(1:2*blockside,1:2*blockside,i-low_i+1)

				!Filling the bessel block
				m_AB_sca(low_row:high_row,low_col:high_col)=&
				m_AB_sca(low_row:high_row,low_col:high_col)+m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)

			END DO fill_do

			!If the job is finished, I leave the loop
			IF (high_i==(namax-namin+1)*(nbmax-nbmin+1)) EXIT unroll_do

			!Refresh the indexes for the openMP loop
			low_i=high_i+1
			high_i=high_i+ncore
			high_if: IF (high_i>=((namax-namin+1)*(nbmax-nbmin+1))) THEN
				high_i=((namax-namin+1)*(nbmax-nbmin+1))
			END IF high_if

		END DO unroll_do

	END DO fill_col_block_do
END DO fill_row_block_do

!Freeing some memory
DEALLOCATE(m_block,m_block_sca)




!-------------------------------------------------------------------------------------------------------------
! Filling the remainder of the matrix exploiting symmetries
!-------------------------------------------------------------------------------------------------------------
sim_row_cell_do: DO icell=1,ns
	sim_col_cell_do: DO jcell=icell,ns
 
		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)
		low_col=2*blockside*(jcell-1)

		!-------------------------------------------------------------------------------------
		!Filling the (icell,jcell) block exploiting the symmetries
		!-------------------------------------------------------------------------------------
		sim_nu_do: DO nu=1,nstop
			sim_mu_do: DO mu=-nu,nu
				sim_n_do: DO n=nu,nstop


					!-------------------------------------------------------------
					!Filling the diagonal block inside the block
					!-------------------------------------------------------------
					sim_if: IF (nu==n) THEN

						sim_m_do1: DO m=-n,-mu

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

						END DO sim_m_do1

					!-------------------------------------------------------------
					!Filling the off diagonal block inside the block
					!-------------------------------------------------------------
					ELSE

						sim_m_do2: DO m=-n,n

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
							!101 FORMAT (8I4)
							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
							!102 FORMAT (8I4)
							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

							!WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
							!103 FORMAT (8I4)


						END DO sim_m_do2

					END IF sim_if

				END DO sim_n_do
			END DO sim_mu_do
		END DO sim_nu_do

		!--------------------------------------------------------
		! Symmetry cannot be exploited if icell=jcell, of course
		!--------------------------------------------------------
		IF (icell==jcell) CYCLE sim_col_cell_do

		!------------------------------------------------------------------------------------------------
		!Filling the (jcell,icell) block exploiting the symmetry Aij=(-1)^(n+nu)Aji,Bij=(-1)^(n+nu+1)Bji
		!------------------------------------------------------------------------------------------------
		ij_n_do: DO n=1,nstop
			ij_m_do: DO m=-n,n
				ij_nu_do: DO nu=1,nstop
					ij_mu_do: DO mu=-nu,nu

						!Original index (ia1,ja1) (ib1,jb1), for the already filled part
						ia1=low_row+2*( n*(n+1)+m )-1
						ja1=low_col+2*( nu*(nu+1)+mu )-1
						ib1=ia1+1
						jb1=ja1

						!Index (ia2,ja2) (ib2,jb2) for the the mapping (i,j)->(j,i), low_col and low_row are exchanged
						ia2=low_col+2*( n*(n+1)+m )-1
						ja2=low_row+2*( nu*(nu+1)+mu )-1
						ib2=ia2+1
						jb2=ja2

						!Index (ia3,ja3) (ib3,jb3) for the second simmetry, i.e. the small block of four for (j,i)
						ia3=ia2+1
						ja3=ja2+1
						ib3=ib2-1
						jb3=jb2+1

						!Mapping the values from (i,j) to (j,i)
						m_AB(ia2,ja2)=v_oner(n+nu)*m_AB(ia1,ja1)
						m_AB(ib2,jb2)=v_oner(n+nu+1)*m_AB(ib1,jb1)
						m_AB_sca(ia2,ja2)=v_oner(n+nu)*m_AB_sca(ia1,ja1)
						m_AB_sca(ib2,jb2)=v_oner(n+nu+1)*m_AB_sca(ib1,jb1)

						!Filling with the simmetry on the block of four
						m_AB(ia3,ja3)=m_AB(ia2,ja2)
						m_AB(ib3,jb3)=m_AB(ib2,jb2)
						m_AB_sca(ia3,ja3)=m_AB_sca(ia2,ja2)
						m_AB_sca(ib3,jb3)=m_AB_sca(ib2,jb2)

					END DO ij_mu_do
				END DO ij_nu_do
			END DO ij_m_do
		END DO ij_n_do


	END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!Adding the diagonal terms
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
	diagn_do: DO n=1,nstop
		diagm_do: DO m=-n,n


			j=j+1
			ant_if: IF (i<=ns_ant) THEN
				m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
				m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))
			END IF ant_if

			m_AB_sca(2*j-1,2*j-1)=m_AB_sca(2*j-1,2*j-1)+1.0D0
			m_AB_sca(2*j,2*j)=m_AB_sca(2*j,2*j)+1.0D0

		END DO diagm_do
	END DO diagn_do
END DO diagi_do

END SUBROUTINE fill_AB_per_sca_fast_asym_dip



!*********************************************************************************************************************************
!2.8) SUBROUTINE update_AB_per_sca_fast:updating AB with fast routines openMP friendly. Here I update only when i do the mapping
!displacing the dipole
!*********************************************************************************************************************************
SUBROUTINE update_AB_per_sca_fast_asym_dip(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
				m_index,m_Apmn,v_jABij_template,v_iABij_template,m_a,m_b,m_AB,m_AB_sca)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k								! vettore d'onda
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r						! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index					! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_jABij_template,v_iABij_template		! Vettori sparse colonne e righe template, da copiare a piedi pari
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b					! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_AB_sca				! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr			! Indici
INTEGER(lo) :: n,m,mu,nu,na,nb,blockside,ncore,omp_get_num_procs,chunk	! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2		! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col,high_row, high_col,low_i,high_i
INTEGER(lo), DIMENSION(1:(namax-namin+1)*(nbmax-nbmin+1),1:2) :: m_nanb
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:,:) :: m_block,m_block_sca			! Temporary block matrix storage

! Subroutine

!Self commenting :-)
blockside=nstop*(nstop+2)
ncore=1!omp_get_num_procs()

!Filling the linearized index vector
i=0
na_do: DO na=namin,namax
	nb_do: DO nb=nbmin,nbmax
		i=i+1
		m_nanb(i,1)=na
		m_nanb(i,2)=nb
	END DO nb_do
END DO na_do

!Allocating temporary vectors according to the number of cores
ALLOCATE(m_block(1:2*blockside,1:2*blockside,1:ncore),m_block_sca(1:2*blockside,1:2*blockside,1:ncore),STAT=status_arr)
m_block=0.0D0
m_block_sca=0.0D0

!-------------------------------------------------------------------------------------------------------------
!First Part: updating direct values, what remains is filled by simmetry
!-------------------------------------------------------------------------------------------------------------
fill_row_block_do: DO icell=1,ns_ant
	fill_col_block_do: DO jcell=ns,ns

		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)+1
		low_col=2*blockside*(jcell-1)+1
		high_row=2*blockside*(icell)
		high_col=2*blockside*(jcell)

		!Bounds for the unrolled loop
		low_i=1
		high_i=ncore

		!Here I unroll the loop so as to preserve the openMP implementation
		unroll_do: DO

			!$OMP PARALLEL PRIVATE(i)
			chunk=1		!Choosing the chunk size
			!$OMP DO SCHEDULE(guided,chunk)
			linear_do: DO i=low_i,high_i

!				WRITE(*,*) 'i,na,nb:',i,m_nanb(i,1),m_nanb(i,2)

				!Diagonal factors at the end to preserve simmetries
				cell_if: IF ( (icell==jcell) .AND. (m_nanb(i,1)==0) .AND. (m_nanb(i,2)==0) ) THEN

!					WRITE(*,*) 'qui1,i,low_i,high_i,i-low_i+1: ',i,low_i,high_i,i-low_i+1

					m_block(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)
					m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)=(0.0D0,0.0D0)

!					WRITE(*,*) 'qui2'

				ELSE

!					WRITE(*,*) 'qui3'

				!Subroutine to fill a whole dense block
				CALL fillblock_dense_per_sca_asym(k,v_r,m_xyz,v_c0,v_aq,v_bq,&
							     m_index,m_Apmn,v_jABij_template,v_iABij_template,&
							     m_nanb(i,1),m_nanb(i,2),icell,jcell,&
							     m_block(:,:,i-low_i+1),m_block_sca(:,:,i-low_i+1))

!					WRITE(*,*) 'qui4'

				END IF cell_if

			END DO linear_do
			!$OMP END DO NOWAIT
			!$OMP END PARALLEL

			!Here i do the sum to fill the compact AB and AB_sca matrix
			fill_do: DO i=low_i,high_i

				!Filling the hankel block
				m_AB(low_row:high_row,low_col:high_col)=&
				m_AB(low_row:high_row,low_col:high_col)+m_block(1:2*blockside,1:2*blockside,i-low_i+1)

				!Filling the bessel block
				m_AB_sca(low_row:high_row,low_col:high_col)=&
				m_AB_sca(low_row:high_row,low_col:high_col)+m_block_sca(1:2*blockside,1:2*blockside,i-low_i+1)

			END DO fill_do

			!If the job is finished, I leave the loop
			IF (high_i==(namax-namin+1)*(nbmax-nbmin+1)) EXIT unroll_do

			!Refresh the indexes for the openMP loop
			low_i=high_i+1
			high_i=high_i+ncore
			high_if: IF (high_i>=((namax-namin+1)*(nbmax-nbmin+1))) THEN
				high_i=((namax-namin+1)*(nbmax-nbmin+1))
			END IF high_if

		END DO unroll_do

	END DO fill_col_block_do
END DO fill_row_block_do

!Freeing some memory
DEALLOCATE(m_block,m_block_sca)




!-------------------------------------------------------------------------------------------------------------
! Filling the remainder of the matrix exploiting symmetries
!-------------------------------------------------------------------------------------------------------------
sim_row_cell_do: DO icell=1,ns_ant
	sim_col_cell_do: DO jcell=ns,ns
 
		!Offset indici per il riempimento della matrice
		low_row=2*blockside*(icell-1)
		low_col=2*blockside*(jcell-1)

		!-------------------------------------------------------------------------------------
		!Filling the (icell,jcell) block exploiting the symmetries
		!-------------------------------------------------------------------------------------
		sim_nu_do: DO nu=1,nstop
			sim_mu_do: DO mu=-nu,nu
				sim_n_do: DO n=nu,nstop


					!-------------------------------------------------------------
					!Filling the diagonal block inside the block
					!-------------------------------------------------------------
					sim_if: IF (nu==n) THEN

						sim_m_do1: DO m=-n,-mu

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

						END DO sim_m_do1

					!-------------------------------------------------------------
					!Filling the off diagonal block inside the block
					!-------------------------------------------------------------
					ELSE

						sim_m_do2: DO m=-n,n

							!Original index (ia1,ja1) (ib1,jb1), for the already filled part
							ia1=low_row+2*( n*(n+1)+m )-1
							ja1=low_col+2*( nu*(nu+1)+mu )-1
							ib1=ia1+1
							jb1=ja1

							!Index (ia2,ja2) (ib2,jb2) for the first simmetry, i.e the small block of four
							ia2=ia1+1
							ja2=ja1+1
							ib2=ib1-1
							jb2=jb1+1
							m_AB(ia2,ja2)=m_AB(ia1,ja1)
							m_AB(ib2,jb2)=m_AB(ib1,jb1)
							m_AB_sca(ia2,ja2)=m_AB_sca(ia1,ja1)
							m_AB_sca(ib2,jb2)=m_AB_sca(ib1,jb1)

							!WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
							!101 FORMAT (8I4)
							!Index (ia3,ja3) (ib3,jb3) for the third simmetry, i.e. n=nu, m=-mu
							ia3=low_row+2*( nu*(nu+1)-mu )-1
							ja3=low_col+2*( n*(n+1)-m )-1
							ib3=ia3+1
							jb3=ja3
							m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
							m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
							m_AB_sca(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB_sca(ia1,ja1)
							m_AB_sca(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB_sca(ib1,jb1)

							!WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
							!102 FORMAT (8I4)
							!Index (ia4,ja4) (ib4,jb4) for the third simmetry, i.e. the small block of four for n=nu, m=-mu
							ia4=ia3+1
							ja4=ja3+1
							ib4=ib3-1
							jb4=jb3+1
							m_AB(ia4,ja4)=m_AB(ia3,ja3)
							m_AB(ib4,jb4)=m_AB(ib3,jb3)
							m_AB_sca(ia4,ja4)=m_AB_sca(ia3,ja3)
							m_AB_sca(ib4,jb4)=m_AB_sca(ib3,jb3)

							!WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
							!103 FORMAT (8I4)


						END DO sim_m_do2

					END IF sim_if

				END DO sim_n_do
			END DO sim_mu_do
		END DO sim_nu_do

		!--------------------------------------------------------
		! Symmetry cannot be exploited if icell=jcell, of course
		!--------------------------------------------------------
		IF (icell==jcell) CYCLE sim_col_cell_do

		!------------------------------------------------------------------------------------------------
		!Filling the (jcell,icell) block exploiting the symmetry Aij=(-1)^(n+nu)Aji,Bij=(-1)^(n+nu+1)Bji
		!------------------------------------------------------------------------------------------------
		ij_n_do: DO n=1,nstop
			ij_m_do: DO m=-n,n
				ij_nu_do: DO nu=1,nstop
					ij_mu_do: DO mu=-nu,nu

						!Original index (ia1,ja1) (ib1,jb1), for the already filled part
						ia1=low_row+2*( n*(n+1)+m )-1
						ja1=low_col+2*( nu*(nu+1)+mu )-1
						ib1=ia1+1
						jb1=ja1

						!Index (ia2,ja2) (ib2,jb2) for the the mapping (i,j)->(j,i), low_col and low_row are exchanged
						ia2=low_col+2*( n*(n+1)+m )-1
						ja2=low_row+2*( nu*(nu+1)+mu )-1
						ib2=ia2+1
						jb2=ja2

						!Index (ia3,ja3) (ib3,jb3) for the second simmetry, i.e. the small block of four for (j,i)
						ia3=ia2+1
						ja3=ja2+1
						ib3=ib2-1
						jb3=jb2+1

						!Mapping the values from (i,j) to (j,i)
						m_AB(ia2,ja2)=v_oner(n+nu)*m_AB(ia1,ja1)
						m_AB(ib2,jb2)=v_oner(n+nu+1)*m_AB(ib1,jb1)
						m_AB_sca(ia2,ja2)=v_oner(n+nu)*m_AB_sca(ia1,ja1)
						m_AB_sca(ib2,jb2)=v_oner(n+nu+1)*m_AB_sca(ib1,jb1)

						!Filling with the simmetry on the block of four
						m_AB(ia3,ja3)=m_AB(ia2,ja2)
						m_AB(ib3,jb3)=m_AB(ib2,jb2)
						m_AB_sca(ia3,ja3)=m_AB_sca(ia2,ja2)
						m_AB_sca(ib3,jb3)=m_AB_sca(ib2,jb2)

					END DO ij_mu_do
				END DO ij_nu_do
			END DO ij_m_do
		END DO ij_n_do


	END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!Adding the diagonal terms, but only for ns and the scattering coefficient matrix
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
	diagn_do: DO n=1,nstop
		diagm_do: DO m=-n,n

			j=j+1
			ant_if: IF (i==ns) THEN
				m_AB_sca(2*j-1,2*j-1)=m_AB_sca(2*j-1,2*j-1)+1.0D0
				m_AB_sca(2*j,2*j)=m_AB_sca(2*j,2*j)+1.0D0
			END IF ant_if

		END DO diagm_do
	END DO diagn_do
END DO diagi_do

END SUBROUTINE update_AB_per_sca_fast_asym_dip



!******************************************************************************
!2bis) SUBROUTINE fill_AB_per_dip:subroutine per riempire  la matrice densa periodica,
!              quella super compressa con la bella idea, e pure lo scattering!!!
!******************************************************************************
SUBROUTINE fill_AB_per_dip(k,fint,v_r,m_xyz,v_c0,v_aq,v_bq,m_a,m_b,m_AB,m_ABJ)

IMPLICIT NONE

! Dichiarazione dummy arguments
REAL(dbl) ,INTENT(IN) :: k                                          ! vettore d'onda
REAL(dbl), INTENT(IN) :: fint                                       ! Coefficiente di interazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_r                          ! Vettore raggi 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq               ! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_a,m_b                 ! Matrici contenenti i coeff di sfera singola
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_AB,m_ABJ             ! Matrice di output



! Dichiarazione variabili interne
INTEGER(lo) :: i,j,error,icell,jcell,icell1,jcell1,status_arr   ! Indici
INTEGER(lo) :: n,m,mu,nu,NNZab,NNZd,na,nb,blockside,matrixside    ! Indici
INTEGER(lo) :: ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8,iab,jab,iab2,jab2  ! Indici per lo sfruttamento delle simm.
INTEGER(lo) :: ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
INTEGER(lo) :: ib1,ib2,ib3,ib4,ib5,ib6,ib7,ib8
INTEGER(lo) :: jb1,jb2,jb3,jb4,jb5,jb6,jb7,jb8
INTEGER(lo) :: low_row, low_col
REAL(dbl) :: xij,yij,zij,rij,xi,yi,zi,xj,yj,zj,xj0,yj0,dij,f,thetaij,phiij,sigx,sigy
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca  	! Matrice blocchi Aij  Bij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij         			! Matrici indici Aij Bij
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Dnkm                      			! Matrice blocchi Dij
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jDnkm,v_iDnkm         			! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop) :: v_exphi                     			! Matrice blocchi PHIij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_rhs,m_rhs_sca           			! Matrice Right Hand Side

! Funzione vera e propria
   
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!Lato di un blocco
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside


!-------------------------------------------------------------------------------------------------------------
!PRIMA PARTE DELLA SUBROUTINE,RIEMPIO, CON SOMMA TUTTI I POSTI A MENO DI UNA SIMMETRIA!!!
!-------------------------------------------------------------------------------------------------------------

!Alloco tutti i vettori necessari ai miei calcoli
ALLOCATE(v_Aij(1:NNZab),v_Bij(1:NNZab),v_Aij_sca(1:NNZab),v_Bij_sca(1:NNZab), &
	   & v_jABij(1:NNZab),v_iABij(1:blockside+1),STAT=status_arr)
ALLOCATE(v_Dnkm(1:NNZd),v_jDnkm(1:NNZd),v_iDnkm(1:blockside+1),STAT=status_arr)
ALLOCATE(m_rhs(1:blockside,1:8),m_rhs_sca(1:blockside,1:8),STAT=status_arr)

fill_row_block_do: DO icell=1,ns
    fill_col_block_do: DO jcell=1,ns

		IF (icell==1) CYCLE fill_col_block_do

        !Offset indici per il riempimento della matrice
        low_row=2*blockside*(icell-1)
        low_col=2*blockside*(jcell-1)

        fill_row_cell_do: DO na=0,nacell-1
            fill_col_cell_do: DO nb=0,nbcell-1

                !Se sono su un blocco diagonale, non sommo subito i fattori diagonali che rompono le simmetrie
                IF ( (icell==jcell) .AND. (na==0) .AND. (nb==0) ) CYCLE fill_col_cell_do

                !Coordinate e distanze
                xi=m_xyz(icell,1)
                yi=m_xyz(icell,2)
                zi=m_xyz(icell,3)

                xj0=m_xyz(jcell,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
                yj0=m_xyz(jcell,2) + REAL(nb,dbl)*by
                zj=m_xyz(jcell,3)

                !Coordinate per la seconda sfere, comprese complicate traslazioni
                cell_row_do_nanb: DO icell1=0,nacell,nacell
                    cell_col_do_nanb: DO jcell1=0,nbcell,nbcell
             
                        sigx=SIGN(1.0D0,(xi-xj0))
                        sigy=SIGN(1.0D0,(yi-yj0)) 

                        xj=xj0 + sigx*REAL(icell1,dbl)*ax + sigy*REAL(jcell1,dbl)*bx
                        yj=yj0 + sigy*REAL(jcell1,dbl)*by

                        dij=SQRT((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

                       !Calcolo il coeff d'interazione
                       f=(v_r(icell)+v_r(jcell))/dij

                       !Se la distanza e' grande allora il blocco e' zero e non tocco gli indici
                       IF (f>=fint) EXIT cell_row_do_nanb

                   END DO cell_col_do_nanb
               END DO cell_row_do_nanb

			   IF (f<fint) CYCLE fill_col_cell_do

               !Calcolo le coordinate sferiche relative
                xij=xi-xj
                yij=yi-yj
                zij=zi-zj
        
               !Conversione coordinate da cartesiane a sferiche
               CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
               !Chiamo la subroutine per riempire Aij e Bij, Dij e PHIij
               CALL fillblock_sparse_dip(nstop,v_c0,v_aq,v_bq,k*rij,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)
               CALL fillblock_Dkmn(nstop,thetaij,v_Dnkm,v_jDnkm,v_iDnkm,error)
               
               !Riempio anche la colonna degli esponenziali
                phi_do: DO m=-nstop,nstop

                    v_exphi(m)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

                END DO phi_do

                !Qui alloco la matrice che uso per i vettori Right Hand Side
                m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
                m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

               !Adesso comincio il ciclo su mu e nu, con tutte le moltiplicazioni del caso
               fill_nu_do: DO nu=1,nstop
                fill_mu_do: DO mu= -nu,nu

				    m_rhs(1:blockside,1:8)=(0.0D0,0.0D0)
                    m_rhs_sca(1:blockside,1:8)=(0.0D0,0.0D0)

                    !Qui estraggo una colonna (riga e' uguale), della mat diag e poi la densifico
                    m_rhs(nu*(nu+1)+mu,1) = v_exphi(mu)
                    m_rhs_sca(nu*(nu+1)+mu,1) = v_exphi(mu)

                    !Faccio tutte le mie moltiplicazioni vettore matrice del caso, Qui rotazione
                    CALL amudz(blockside,m_rhs(:,1),m_rhs(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
                    CALL amudz(blockside,m_rhs_sca(:,1),m_rhs_sca(:,2),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm
                    
                    !Due traslazioni
                    CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,3),v_Aij,v_jABij,v_iABij)   !Per Aij
                    CALL amuzz(blockside,m_rhs(:,2),m_rhs(:,4),v_Bij,v_jABij,v_iABij)   !Per Bij
                    CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,3),v_Aij_sca,v_jABij,v_iABij)   !Per AJij
                    CALL amuzz(blockside,m_rhs_sca(:,2),m_rhs_sca(:,4),v_Bij_sca,v_jABij,v_iABij)   !Per BJij                    

                    !Aggiustamento della parita'
                    par_ndo: DO n=1,nstop
                        par_mdo: DO m=-n,n

                            m_rhs(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,3)
                            m_rhs(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs(n*(n+1)+m,4)
                            m_rhs_sca(n*(n+1)+m,3)=((-1.0D0)**m)*m_rhs_sca(n*(n+1)+m,3)
                            m_rhs_sca(n*(n+1)+m,4)=((-1.0D0)**m)*m_rhs_sca(n*(n+1)+m,4)                            

                        END DO par_mdo 
                    END DO par_ndo


                    !Rotazioni per ciascuna delle due branche
                    CALL amudz(blockside,m_rhs(:,3),m_rhs(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte A
                    CALL amudz(blockside,m_rhs(:,4),m_rhs(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte B
                    CALL amudz(blockside,m_rhs_sca(:,3),m_rhs_sca(:,5),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte AJ
                    CALL amudz(blockside,m_rhs_sca(:,4),m_rhs_sca(:,6),v_Dnkm,v_jDnkm,v_iDnkm)  !Per Dnkm, parte BJ

                    !Estraggo le sottomatrici e poi moltiplico
                    !Aggiustamento della parita'
                    i=0
                    phi_ndo: DO n=nu,nstop
                        
                        diag_if: IF (n==nu) THEN

                            phi_mdo1: DO m=-n,-mu
                                i=n*(n+1)+m
                                m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
                                m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
                                m_rhs_sca(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,5)
                                m_rhs_sca(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,6)
                            END DO phi_mdo1

                        ELSE

                            phi_mdo2: DO m=-n,n
                                i=n*(n+1)+m
                                m_rhs(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,5)
                                m_rhs(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs(i,6)
                                m_rhs_sca(i,7)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,5)
                                m_rhs_sca(i,8)=v_exphi(-m)*((-1.0D0)**m)*m_rhs_sca(i,6)
                            END DO phi_mdo2

                        END IF diag_if 

                    END DO phi_ndo

                    !Adesso riempio dai due vettori i pezzi di matrice che mi servono!!!
                    fill_ndo: DO n=nu,nstop
                        
                        fill_diag_if: IF (n==nu) THEN

                            !Qui riempio il blocco diagonale    
                            fill_mdo1: DO m=-n,-mu

                                iab=2*( n*(n+1)+m )
                                jab=2*( nu*(nu+1)+mu )
                                iab2=n*(n+1)+m
                                jab2=nu*(nu+1)+mu

                                !Indici per A e B
                                ia1=low_row+iab-1
                                ja1=low_col+jab-1
                                ib1=low_row+iab
                                jb1=low_col+jab-1

                                !Assegno i valori passando dai vettori che tengono a e b alla matrice
                                !Metto la somma per via della nuova teoria
                                m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
                                m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)
                                m_ABJ(ia1,ja1)=m_ABJ(ia1,ja1)+m_rhs_sca(iab2,7)
                                m_ABJ(ib1,jb1)=m_ABJ(ib1,jb1)+m_rhs_sca(iab2,8)
                                
                            END DO fill_mdo1

                        ELSE

                            !Qui riempio il blocco fuori dalla diagonale diagonale
                            fill_mdo2: DO m=-n,n

                                iab=2*( n*(n+1)+m )
                                jab=2*( nu*(nu+1)+mu )
                                iab2=n*(n+1)+m
                                jab2=nu*(nu+1)+mu
                                
                                !Indici per A e B
                                ia1=low_row+iab-1
                                ja1=low_col+jab-1
                                ib1=low_row+iab
                                jb1=low_col+jab-1

                                !Assegno i valori passando dai vettori che tengono a e b alla matrice
                                !Metto la somma per via della nuova teoria
                                m_AB(ia1,ja1)=m_AB(ia1,ja1)+m_rhs(iab2,7)
                                m_AB(ib1,jb1)=m_AB(ib1,jb1)+m_rhs(iab2,8)
                                m_ABJ(ia1,ja1)=m_ABJ(ia1,ja1)+m_rhs_sca(iab2,7)
                                m_ABJ(ib1,jb1)=m_ABJ(ib1,jb1)+m_rhs_sca(iab2,8)

                            END DO fill_mdo2



                        END IF fill_diag_if

                    END DO fill_ndo


                END DO fill_mu_do
               END DO fill_nu_do


            END DO fill_col_cell_do
        END DO fill_row_cell_do
    END DO fill_col_block_do
END DO fill_row_block_do

!-------------------------------------------------------------------------------------------------------------
!SECONDA PARTE DELLA SUBROUTINE: FATTE LE SOMME SFRUTTO LE SIMMETRIE
!-------------------------------------------------------------------------------------------------------------

sim_row_cell_do: DO icell=1,ns
    sim_col_cell_do: DO jcell=1,ns
    
		IF (icell==1) CYCLE sim_col_cell_do

		!Offset indici per il riempimento della matrice
        low_row=2*blockside*(icell-1)
        low_col=2*blockside*(jcell-1)

        sim_nu_do: DO nu=1,nstop
            sim_mu_do: DO mu=-nu,nu
                sim_n_do: DO n=nu,nstop

                    sim_if: IF (nu==n) THEN

                        sim_m_do1: DO m=-n,-mu
                            
                            !Indici per A e B originali
                            ia1=low_row+2*( n*(n+1)+m )-1
                            ja1=low_col+2*( nu*(nu+1)+mu )-1
                            ib1=ia1+1
                            jb1=ja1

                            !Indici per i primi doppioni
                            ia2=ia1+1
                            ja2=ja1+1
                            ib2=ib1-1
                            jb2=jb1+1
                            m_AB(ia2,ja2)=m_AB(ia1,ja1)
                            m_AB(ib2,jb2)=m_AB(ib1,jb1)
                            m_ABJ(ia2,ja2)=m_ABJ(ia1,ja1)
                            m_ABJ(ib2,jb2)=m_ABJ(ib1,jb1)

                            !Indici per prima simmetria da scambio di indici n=nu, m=-mu
                            ia3=low_row+2*( nu*(nu+1)-mu )-1
                            ja3=low_col+2*( n*(n+1)-m )-1
                            ib3=ia3+1
                            jb3=ja3
                            m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
                            m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
                            m_ABJ(ia3,ja3)=((-1.0D0)**(m+mu))*m_ABJ(ia1,ja1)
                            m_ABJ(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_ABJ(ib1,jb1)

                            !Indici per i doppioni della simmetria di indici
                            ia4=ia3+1
                            ja4=ja3+1
                            ib4=ib3-1
                            jb4=jb3+1
                            m_AB(ia4,ja4)=m_AB(ia3,ja3)
                            m_AB(ib4,jb4)=m_AB(ib3,jb3)
                            m_ABJ(ia4,ja4)=m_ABJ(ia3,ja3)
                            m_ABJ(ib4,jb4)=m_ABJ(ib3,jb3)

                        END DO sim_m_do1

                    ELSE

                        sim_m_do2: DO m=-n,n

                            !Indici per A e B originali
                            ia1=low_row+2*( n*(n+1)+m )-1
                            ja1=low_col+2*( nu*(nu+1)+mu )-1
                            ib1=ia1+1
                            jb1=ja1

                            !Indici per i primi doppioni
                            ia2=ia1+1
                            ja2=ja1+1
                            ib2=ib1-1
                            jb2=jb1+1
                            m_AB(ia2,ja2)=m_AB(ia1,ja1)
                            m_AB(ib2,jb2)=m_AB(ib1,jb1)
                            m_ABJ(ia2,ja2)=m_ABJ(ia1,ja1)
                            m_ABJ(ib2,jb2)=m_ABJ(ib1,jb1)
                                
!                                 WRITE(*,101) ia1,ja1,ib1,jb1,ia2,ja2,ib2,jb2
! 				101 FORMAT (8I4)
                            !Indici per prima simmetria da scambio di indici n=nu, m=-mu
                            ia3=low_row+2*( nu*(nu+1)-mu )-1
                            ja3=low_col+2*( n*(n+1)-m )-1
                            ib3=ia3+1
                            jb3=ja3
                            m_AB(ia3,ja3)=((-1.0D0)**(m+mu))*m_AB(ia1,ja1)
                            m_AB(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_AB(ib1,jb1)
                            m_ABJ(ia3,ja3)=((-1.0D0)**(m+mu))*m_ABJ(ia1,ja1)
                            m_ABJ(ib3,jb3)=((-1.0D0)**(m+mu+1))*m_ABJ(ib1,jb1)
!                                 
!                                 WRITE(*,102) ia1,ja1,ib1,jb1,ia3,ja3,ib3,jb3
! 				102 FORMAT (8I4)
                            !Indici per i doppioni della simmetria di indici
                            ia4=ia3+1
                            ja4=ja3+1
                            ib4=ib3-1
                            jb4=jb3+1
                            m_AB(ia4,ja4)=m_AB(ia3,ja3)
                            m_AB(ib4,jb4)=m_AB(ib3,jb3)
                            m_ABJ(ia4,ja4)=m_ABJ(ia3,ja3)
                            m_ABJ(ib4,jb4)=m_ABJ(ib3,jb3)
!                                 
!                                 WRITE(*,103) ia1,ja1,ib1,jb1,ia4,ja4,ib4,jb4
! 				103 FORMAT (8I4)


                        END DO sim_m_do2

                    END IF sim_if

                END DO sim_n_do
            END DO sim_mu_do
        END DO sim_nu_do
    END DO sim_col_cell_do
END DO sim_row_cell_do

!-------------------------------------------------------------------------------------------------------------
!ULTIMA PARTE DELLA SUBROUTINE, AGGIUNGO I TERMINI DIAGONALI!!!
!-------------------------------------------------------------------------------------------------------------
j=0
diagi_do: DO i=1,ns
    diagn_do: DO n=1,nstop
        diagm_do: DO m=-n,n
        
            j=j+1
            m_AB(2*j-1,2*j-1)=m_AB(2*j-1,2*j-1)+1.0D0/m_a(n,v_patt(i))
            m_AB(2*j,2*j)=m_AB(2*j,2*j)+1.0D0/m_b(n,v_patt(i))
            m_ABJ(2*j-1,2*j-1)=m_ABJ(2*j-1,2*j-1)+1.0D0
            m_ABJ(2*j,2*j)=m_ABJ(2*j,2*j)+1.0D0


        END DO diagm_do
    END DO diagn_do
END DO diagi_do

                                
!                                 WRITE(*,*) "QUId"


END SUBROUTINE fill_AB_per_dip









! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! SPARSE BLAS
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!******************************************************************************
!1)  FUNCTION Getelm:estraggo un elemento sparse
!******************************************************************************


COMPLEX(dbl) function getelm (i,j,a,ja,ia,iadd,sorted) 
!c-----------------------------------------------------------------------
!c     purpose:
!c     -------- 
!c     this function returns the element a(i,j) of a matrix a, 
!c     for any pair (i,j).  the matrix is assumed to be stored 
!c     in compressed sparse row (csr) format. getelm performs a
!c     binary search in the case where it is known that the elements 
!c     are sorted so that the column indices are in increasing order. 
!c     also returns (in iadd) the address of the element a(i,j) in 
!c     arrays a and ja when the search is successsful (zero if not).
!c----- 
!c     first contributed by noel nachtigal (mit). 
!c     recoded jan. 20, 1991, by y. saad [in particular
!c     added handling of the non-sorted case + the iadd output] 
!c-----------------------------------------------------------------------
!c     parameters:
!c     ----------- 
!c on entry: 
!c---------- 
!c     i      = the row index of the element sought (input).
!c     j      = the column index of the element sought (input).
!c     a      = the matrix a in compressed sparse row format (input).
!c     ja     = the array of column indices (input).
!c     ia     = the array of pointers to the rows' data (input).
!c     sorted = logical indicating whether the matrix is knonw to 
!c              have its column indices sorted in increasing order 
!c              (sorted=.true.) or not (sorted=.false.).
!c              (input). 
!c on return:
!c----------- 
!c     getelm = value of a(i,j). 
!c     iadd   = address of element a(i,j) in arrays a, ja if found,
!c              zero if not found. (output) 
!c
!c     note: the inputs i and j are not checked for validity. 
!c-----------------------------------------------------------------------
!c     noel m. nachtigal october 28, 1990 -- youcef saad jan 20, 1991.
!c----------------------------------------------------------------------- 
      integer(lo) i, ia(:), iadd, j, ja(:)
      COMPLEX(dbl) a(:)
      logical sorted 
!c
!c     local variables.
!c
      integer(lo) ibeg, iend, imid, k
!c
!c     initialization 
!c
      iadd = 0 
      getelm = 0.0
      ibeg = ia(i)
      iend = ia(i+1)-1
!c
!c     case where matrix is not necessarily sorted
!c     
      if (.not. sorted) then 
!c
!c scan the row - exit as soon as a(i,j) is found
!c
         do 5  k=ibeg, iend
            if (ja(k) .eq.  j) then
               iadd = k 
               goto 20 
            endif
 5       continue
!c     
!c     end unsorted case. begin sorted case
!c     
      else
!c     
!c     begin binary search.   compute the middle index.
!c     
 10      imid = ( ibeg + iend ) / 2
!c     
!c     test if  found
!c     
         if (ja(imid).eq.j) then
            iadd = imid 
            goto 20
         endif
         if (ibeg .ge. iend) goto 20
!c     
!c     else     update the interval bounds. 
!c     
         if (ja(imid).gt.j) then
            iend = imid -1
         else 
            ibeg = imid +1
         endif
         goto 10  
!c     
!c     end both cases
!c     
      endif
!c     
 20   if (iadd .ne. 0) getelm = a(iadd) 
!c
      return
!c--------end-of-getelm--------------------------------------------------
!c-----------------------------------------------------------------------
end function getelm


!******************************************************************************
!2)  FUNCTION amudd: real sparse matrix - real vector multiplication:
!			a) Saad version: universal but slower
!			b) MKL: faster but requires MKL libraries
!			Choose one and comment out the other
!******************************************************************************

SUBROUTINE amudd (n, x, y, a,ja,ia) 
REAL(dbl) :: x(:), y(:), a(:) 
INTEGER(lo) :: n, ja(:), ia(:)


!#######################################################################
!#######################################################################
!a)Saad Implementation
!#######################################################################
!#######################################################################
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      REAL(dbl) :: t
      INTEGER(lo) i, k
!-----------------------------------------------------------------------
      outer_do: do i = 1,n
!
!     compute the inner product of row i with vector x
! 
         t = 0.0d0
         inner_do: do k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
         END DO inner_do
!
!     store result in y(i) 
!
         y(i) = t
      END DO outer_do

!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------




!#######################################################################
!#######################################################################
!b)MKL Implementation
!#######################################################################
!#######################################################################
!call mkl_dcsrgemv('N',n,a,ia,ja,x,y)

end SUBROUTINE amudd

!******************************************************************************
!3)  FUNCTION amuzz: complex sparse matrix - complex vector multiplication:
!			a) Saad version: universal but slower
!			b) MKL: faster but requires MKL libraries
!			Choose one and comment out the other
!******************************************************************************

SUBROUTINE amuzz (n, x, y, a,ja,ia) 
COMPLEX(dbl) :: x(:), y(:), a(:) 
INTEGER(lo) :: n, ja(:), ia(:)




!#######################################################################
!#######################################################################
!a)Saad Implementation
!#######################################################################
!#######################################################################
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      COMPLEX(dbl) :: t
      INTEGER(lo) i, k
!-----------------------------------------------------------------------
      outer_do: do i = 1,n
!
!     compute the inner product of row i with vector x
! 
         t = 0.0d0
         inner_do: do k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
         END DO inner_do
!
!     store result in y(i) 
!
         y(i) = t
      END DO outer_do

!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------




!#######################################################################
!#######################################################################
!b)MKL Implementation
!#######################################################################
!#######################################################################
!call mkl_zcsrgemv('N',n,a,ia,ja,x,y)

end SUBROUTINE amuzz


!******************************************************************************
!4)  FUNCTION amudz: real sparse matrix - complex vector multiplication:
!			a) Saad version: universal but slower
!			b) MKL: faster but requires MKL libraries
!			Choose one and comment out the other
!******************************************************************************

SUBROUTINE amudz (n, x, y, a,ja,ia) 
COMPLEX(dbl) :: x(:), y(:)
REAL(dbl) :: a(:)
INTEGER(lo) :: n, ja(:), ia(:)


!!#######################################################################
!!#######################################################################
!!a)Saad Implementation
!!#######################################################################
!!#######################################################################
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      COMPLEX(dbl) :: t
      INTEGER(lo) i, k
!-----------------------------------------------------------------------
      outer_do: do i = 1,n
!
!     compute the inner product of row i with vector x
! 
         t = 0.0d0
         inner_do: do k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
         END DO inner_do
!
!     store result in y(i) 
!
         y(i) = t
      END DO outer_do

!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------




!#######################################################################
!#######################################################################
!b)MKL Implementation
!#######################################################################
!#######################################################################
!Local MKL variables
!REAL(dbl), DIMENSION(1:n) :: xr,xim,yr,yim

!!Separating real and imaginary part of the x vector
!xr=REAL(x,dbl)
!xim=AIMAG(x)

!!Separate matrix vector multiplication for the real and imaginary part
!call mkl_dcsrgemv('N',n,a,ia,ja,xr,yr)
!call mkl_dcsrgemv('N',n,a,ia,ja,xim,yim)

!!Putting together real and imaginary part for the final result
!y=CMPLX(yr,yim,dbl)


end SUBROUTINE amudz




!-------------------------------------
!zdot_bcg(n,zx,zy)
!-------------------------------------

function zdot_bcg(n,zx,zy)

! complex inner product function

implicit none

COMPLEX(dbl) :: zdot_bcg
integer,       intent(in):: n
COMPLEX(dbl),intent(in):: zx(n),zy(n)

zdot_bcg=DOT_PRODUCT(zx,zy)
!zdot_bcg = sum( conjg(zx) * zy )

end function zdot_bcg

!-------------------------------------
!dnorm2_bcg(n,zx)
!-------------------------------------
function dnorm2_bcg(n,zx)

! l2 norm function

implicit none
REAL(dbl) :: dnorm2_bcg
integer,       intent(in):: n
COMPLEX(dbl),intent(in):: zx(n)

dnorm2_bcg = sqrt( abs( zdot_bcg(n, zx, zx) ) )

end function dnorm2_bcg


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! extracts the submatrix A(i1:i2,j1:j2) and puts the result in 
! matrix ao,iao,jao
!---- In place: ao,jao,iao may be the same as a,ja,ia.
!-------------- 
! on input
!---------
! n = row dimension of the matrix 
! i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
!          extracted. 
! j1,j2 = two integers with j2 .ge. j1 indicating the range of columns 
!         to be extracted.
!         * There is no checking whether the input values for i1, i2, j1,
!           j2 are between 1 and n. 
! a,
! ja,
! ia    = matrix in compressed sparse row format. 
!
! job   = job indicator: if job .ne. 1 then the real values in a are NOT
!         extracted, only the column indices (i.e. data structure) are.
!         otherwise values as well as column indices are extracted...
!         
! on output
!-------------- 
! nr    = number of rows of submatrix 
! nc    = number of columns of submatrix 
!     * if either of nr or nc is nonpositive the code will quit.
!
! ao,
! jao,iao = extracted matrix in general sparse format with jao containing
!   the column indices,and iao being the pointer to the beginning 
!   of the row,in arrays a,ja.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      subroutine submatc (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      integer(lo) n,job,i1,i2,j1,j2,nr,nc,ia(:),ja(:),jao(:),iao(:)
      COMPLEX(dbl) a(:),ao(:) 

      nr = i2-i1+1
      nc = j2-j1+1
!     
      if ( nr .le. 0 .or. nc .le. 0) return
!     
      klen = 0
!     
!     simple procedure. proceeds row-wise...
!     
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1
!-----------------------------------------------------------------------
         do 60 k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) = j - j1+1
            endif
 60      continue
 100  continue
      iao(nr+1) = klen+1
      return
      end subroutine submatc
!------------end-of submat---------------------------------------------- 
!----------------------------------------------------------------------- 

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! OLD SUBROUTINES
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!******************************************************************************
!2bis) SUBROUTINE c0xusub(nstop,v_c0xu): coefficente per i vector translation coefficent
!******************************************************************************
SUBROUTINE c0xusub(nstop,v_c0xu)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop				! N espansioni multipolari
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_c0xu	! vettore coeff c0

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,r
REAL(dbl) :: mr,nr,mur,nur,logw,sqrtw

! Funzione vera e propria
i=0

n_do: DO n=1,nstop
	m_do: DO m=-n,n
		nu_do: DO nu=1,nstop
			mu_do: DO mu=-nu,nu
			
				mr=REAL(m,dbl)
				nr=REAL(n,dbl)
				mur=REAL(mu,dbl)
				nur=REAL(nu,dbl)
								
				logw =lnf(nr+mr,r)+lnf(nur-mur,r)-lnf(nr-mr,r)-lnf(nur+mur,r)
				
				sqrtw=EXP(logw)*((2*nr+1.0D0)*(2*nur+1.0D0))/(nr*(nr+1.0D0)*nur*(nur+1.0D0))
			
				i=i+1
			
				v_c0xu(i)=0.5D0*((-1.0D0)**m)*SQRT(sqrtw)
			
			END DO mu_do
		END DO nu_do		
	END DO m_do
END DO n_do

END SUBROUTINE c0xusub




!******************************************************************************
!3) SUBROUTINE qmaxsub: indice superiore della sommatoria per i vector translation coefficent
!******************************************************************************
SUBROUTINE qmaxsub(nstop,v_qmax)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop				! N espansioni multipolari
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_qmax	! vettore coeff qmax

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,thirdi
REAL(dbl) :: mr,nr,mur,nur,thirdr

! Funzione vera e propria
i=0

n_do: DO n=1,nstop
	m_do: DO m=-n,n
		nu_do: DO nu=1,nstop
			mu_do: DO mu=-nu,nu
			
				mr=REAL(m,dbl)
				nr=REAL(n,dbl)
				mur=REAL(mu,dbl)
				nur=REAL(nu,dbl)
				
				thirdr=(nr+nur-ABS(mr-mur))/2.0D0
				thirdi=INT(thirdr)
				
				i=i+1
				
				v_qmax(i)=MIN(n,nu,thirdi)
								
			END DO mu_do
		END DO nu_do		
	END DO m_do
END DO n_do

END SUBROUTINE qmaxsub


!******************************************************************************
!4) SUBROUTINE qqmax: indice superiore della sommatoria per i vector translation coefficent
!******************************************************************************
SUBROUTINE qqmaxsub(nstop,v_qqmax)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop				! N espansioni multipolari
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_qqmax	! vettore coeff qqmax

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,thirdi
REAL(dbl) :: mr,nr,mur,nur,thirdr

! Funzione vera e propria
i=0

n_do: DO n=1,nstop
	m_do: DO m=-n,n
		nu_do: DO nu=1,nstop
			mu_do: DO mu=-nu,nu
			
				mr=REAL(m,dbl)
				nr=REAL(n,dbl)
				mur=REAL(mu,dbl)
				nur=REAL(nu,dbl)
				
				thirdr=(1.0D0+nr+nur-ABS(mr-mur))/2.0D0
				thirdi=INT(thirdr)
				
				i=i+1
				
				v_qqmax(i)=MIN(n,nu,thirdi)
				
			END DO mu_do
		END DO nu_do		
	END DO m_do
END DO n_do

END SUBROUTINE qqmaxsub


!******************************************************************************
!6bis) SUBROUTINE fill_aqbq_sub2: infilo tutti i coefficienti di gaunt in un solo vettore lineare
!******************************************************************************
SUBROUTINE fill_aqbq_sub2(nstop,sumq,sumqq,elem2,v_qmax,v_qqmax,v_aq_long,v_bq_long)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,sumq,sumqq,elem2				! N espansioni multipolari e bounds vettori
INTEGER(lo), DIMENSION(1:elem2), INTENT(IN) :: v_qmax,v_qqmax	! vettore coeff qqmax
REAL(dbl), DIMENSION(1:sumq), INTENT(OUT) :: v_aq_long			! Vettore per tutti i coeff di gaunt
REAL(dbl), DIMENSION(1:sumqq), INTENT(OUT) :: v_bq_long			! Vettore per tutti i coeff di gaunt

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,error
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq
REAL(dbl) :: mr,nr,mur,nur
REAL(dbl), DIMENSION(:), ALLOCATABLE :: v_aq,v_bq
REAL(dbl) :: p,p1,p2										!I soliti p,p1,p2
REAL(dbl) :: alphap1,alphap2,alphap3,alphap4				!I soliti coefficienti alpha
REAL(dbl) :: Ap2,Ap3,Ap4									!I soliti coefficienti Ap
REAL(dbl) :: c0,c1,c2										!I coefficienti ricorsivi
INTEGER(lo) :: q,Ap2i										!Come sopra ma integer

! Funzione vera e propria

!Inizializzo indici e bounds

!aq
i=0
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

n_do: DO n=1,nstop
	m_do: DO m=-n,n
		nu_do: DO nu=1,nstop
			mu_do: DO mu=-nu,nu
			
				mr=REAL(m,dbl)
				nr=REAL(n,dbl)
				mur=REAL(mu,dbl)
				nur=REAL(nu,dbl)

				!---------------------------------
				!Parte di storage per aq
				!---------------------------------
				
				!Aggiorno i
				i=i+1
				
				!Aggiorno i bounds
				ia_if: IF (i==1) THEN
					low_aq=1
				ELSE
					low_aq=up_aq+1	
				END IF ia_if
				
				up_aq=low_aq+v_qmax(i)

!				WRITE(*,*) "Bounds aq 1",low_aq,up_aq,v_qmax(i)
						
				!Alloco la memoria per il mio v_aq
				ALLOCATE(v_aq(0:v_qmax(i)))	

				!Archivio i dati nel vettore lungo
				CALL gaunt_xu(-mr,nr,mur,nur,v_qmax(i),v_aq,error)
				v_aq_long(low_aq:up_aq)=v_aq
				
				!Se Qmax=0 allora salto l'assegnazione per bq e ricomincio
				!il ciclo successivo
				IF (v_qqmax(i)==0) THEN
					DEALLOCATE(v_aq) 
					CYCLE mu_do
				END IF
				
				!--------------------------------------------
				!Parte di storage per bq: se sono qui Qmax/=0
				!--------------------------------------------
				
				
				!Aggiorno i bounds
				ib_if: IF (i==1) THEN
					low_bq=1
				ELSE
					low_bq=up_bq+1	
				END IF ib_if
				
				up_bq=low_bq+v_qqmax(i)-1
				
! 				WRITE(*,*) "Bounds bq 1",low_bq,up_bq,v_qqmax(i)
						
				!Alloco la memoria per il mio v_aq
				ALLOCATE(v_bq(1:v_qqmax(i)))
				
				
				!***************************************************************************
				!Calcolo bq che mi serve
				!***************************************************************************
				mmu_if: IF (((m==0).AND.(mu==0)).OR.((mu==-m).AND.(n==nu))) THEN
				
					!----------------
					!Cosi' avro' B==0
					!----------------
					v_bq=0.0D0
					
				ELSE	
				
					!----------------
					!Caso generale
					!----------------
					bq_case: SELECT CASE (v_qqmax(i))
				
					CASE(1) bq_case
					
! 						WRITE(*,*) "qui"
					
						!Calcolo coefficienti parziali
						p=nr+nur-2.0D0
						p1=p+mr-mur
						Ap3=f_Ap(-mr,nr,mur,nur,p+3.0D0)
						
						!Calcolo Bq
						v_bq(1)=v_aq(0)*(2.0D0*p+3.0D0)*Ap3/((p+3.0D0)*(p1+2.0D0))	
						
					CASE DEFAULT bq_case
					
						!Calcolo il primo valore di bq
						p=nr+nur-2.0D0
						p1=p+mr-mur
						Ap3=f_Ap(-mr,nr,mur,nur,p+3.0D0)
						v_bq(1)=v_aq(0)*(2.0D0*p+3.0D0)*Ap3/((p+3.0D0)*(p1+2.0D0))
					
						!Comincio il ciclo do per il calcolo di tutti i bq
						bq_do: DO q=2,v_qqmax(i)
						
							!Calcolo preliminarmente p ed Ap2
							p=nr+nur-2.0D0*REAL(q,dbl)
							Ap2=f_Ap(-mr,nr,mur,nur,p+2.0D0)
							Ap2i=INT(Ap2,lo)
						
							ap_case: SELECT CASE (Ap2i)
							
							CASE(0) ap_case
							
								!Calcolo coefficienti parziali
								p1=p+mr-mur
								p2=p-mr+mur
								Ap3=f_Ap(-mr,nr,mur,nur,p+3.0D0)
								Ap4=f_Ap(-mr,nr,mur,nur,p+4.0D0)
								alphap3=f_alpha(nr,nur,p+3.0D0)
								alphap4=f_alpha(nr,nur,p+4.0D0)
								
								!Calcolo coefficienti ricorsione
								c0=(2.0D0*p+3.0D0)/((p+3.0D0)*(p1+2.0D0)*Ap4)
								c1=Ap3*Ap4+(p+2.0D0)*(p+4.0D0)*(p1+3.0D0)*(p2+3.0D0)*alphap3
								c2=-(p+2.0D0)*(p+3.0D0)*(p2+3.0D0)*(p2+4.0D0)*alphap4
								
								!Calcolo bq
								v_bq(q)=c0*(c1*v_aq(q-1)+c2*v_aq(q-2))
					
								IF (v_bq(1)==0.0D0) THEN
									WRITE(*,*) m,n,mu,nu,c0,c1,c2
									WRITE(*,*) Ap3,Ap4,f_Ap(-mr,nr,mur,nur,p+5.0D0)
									WRITE(*,*)
								END IF
												
							CASE DEFAULT ap_case
							
								!Calcolo coefficienti parziali
								p1=p+mr-mur
								p2=p-mr+mur
								alphap1=f_alpha(nr,nur,p+1.0D0)
								alphap2=f_alpha(nr,nur,p+2.0D0)
								
								!Calcolo coefficienti ricorsione
								c0=(2.0D0*p+3.0D0)/Ap2
								c1=(p+2.0D0)*(p1+1.0D0)*alphap1
								c2=-(p+1.0D0)*(p2+2.0D0)*alphap2
								
								!Calcolo bq
								IF ((q==v_qqmax(i)) .AND. (v_qqmax(i)>v_qmax(i))) THEN
									v_bq(q)=c0*c2*v_aq(q-1)
								ELSE
									v_bq(q)=c0*(c1*v_aq(q)+c2*v_aq(q-1))
								END IF
							
							END SELECT ap_case
							
						END DO bq_do
					
					END SELECT bq_case
			
				END IF mmu_if
			
				!***************************************************************************
				!Fine calcolo bq
				!***************************************************************************
								
! 				WRITE(*,*)	"long,short:", SIZE(v_bq_long(low_bq:up_bq)),SIZE(v_bq)			
! 				WRITE(*,*)	"low,up:", low_bq,up_bq
				!Salvo v_bq nel vettore
				v_bq_long(low_bq:up_bq)=v_bq
				
! 				WRITE(*,*)	v_bq
! 				WRITE(*,*)	v_bq_long(low_bq:up_bq)
! 				WRITE(*,*)
				
				!Disalloco v_aq e v_bq
				DEALLOCATE(v_aq,v_bq)
									
				!Esco ad un certo punto
				IF (i==(nstop*(nstop+2))**2) EXIT n_do
					
				
				
			END DO mu_do
		END DO nu_do		
	END DO m_do
END DO n_do

END SUBROUTINE fill_aqbq_sub2



!******************************************************************************
!3) SUBROUTINE fillblock_sub2: riempio un blocco delle matrici Aij e Bij
!******************************************************************************
SUBROUTINE fillblock_sub2(nstop,v_qmax,v_qqmax,v_c0,v_aq_long,v_bq_long,kr,theta,phi,m_Aij,m_Bij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop									! N espansioni multipolari
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_qmax,v_qqmax			! vettore coeff qqmax
REAL(dbl), INTENT(IN) :: kr,theta,phi								! Parametri coordinate relative
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq_long,v_bq_long		! Vettore per tutti i coeff di gaunt
INTEGER(lo), INTENT(OUT) :: error 								! Flag di errore
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij			! Blocco ij delle matrici A e B

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,j,l,q,pi						!Indici per il ciclo, reali o immaginari
REAL(dbl) :: mr,nr,mur,nur,p								
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq					!Bounds per i vettori aq_long,bq_long
INTEGER(lo):: elem,elem2									!Limiti di uscita per il ciclo
INTEGER(lo) :: pmin,pmax,nbes								!Estremi per i vett e n bes
COMPLEX(dbl) :: coeff2,sommaA,sommaB						!Coefficienti exp e somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h				!Vett comp fun hankel
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_leg				!Vettori Legendre
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn				!Vettori Bessel jn ed yn
LOGICAL :: test

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=1
j=0
l=0
elem=nstop*(nstop+2)
elem2=(nstop*(nstop+2))**2

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

! IF ((kr==1.0D0).AND.(theta==0.0D0).AND.(phi==0.0D0)) THEN
! WRITE(*,*) "vqq2" ,v_qqmax
! END IF
! 
! IF ((kr==1.0D0).AND.(theta==0.0D0).AND.(phi==0.0D0)) THEN
! WRITE(*,*) "Low_bq,Up_bq iniziali, 2", low_bq,up_bq
! END IF

n_do: DO n=1,nstop
	m_do: DO m=-n,n
		nu_do: DO nu=1,nstop
			mu_do: DO mu=-nu,nu
			
				test=(m==0 .AND. mu==0) .OR. (mu==-m .AND. nu==n) .OR. &
				&(m==n .AND. mu==-nu) .OR. (m==-n .AND. mu==nu)
		

				!Aggiorno i
				l=l+1
				j=j+1
				
				!Aggiorno i bounds
				ia_if: IF (l==1) THEN
					low_aq=1
				ELSE
					low_aq=up_aq+1	
				END IF ia_if
				
				up_aq=low_aq+v_qmax(l)				
		
! 				WRITE(*,*) "Bounds aq 2",low_aq,up_aq,v_qmax(l)
			
				mr=REAL(m,dbl)
				nr=REAL(n,dbl)
				mur=REAL(mu,dbl)
				nur=REAL(nu,dbl)

				Bvt_if: IF ((v_qqmax(l)==0)) THEN 
				
					!-----------------------------------------------------------------------
					!Calcolo Bvt
					!-----------------------------------------------------------------------
					!Bvt e' sempre zero in questi casi
					m_Bij(i,j)=(0.0D0,0.0D0)
				
				
					!-----------------------------------------------------------------------
					!Calcolo Avt
					!-----------------------------------------------------------------------
					! Calcolo Pmin e Pmax
					pmin=n+nu -2*v_qmax(l)
					pmax=n+nu
					nbes=pmax-pmin+1
					
					!-------
					!Calcolo coeff1 e coeff2
					coeff2=EXP(CMPLX(0.0D0,(mur-mr)*phi))
					
					!-------
					!Adesso posso allocare tutti ma proprio tutti i miei vettori
					ALLOCATE(v_h(0:pmax),v_leg(pmin:pmax),&
							&v_Apmn(pmin:pmax))
				
					!Inizializzo v_Apmn
					v_Apmn=0.0D0
				
					!-------
					!Calcolo la mia funzione di legendre
					CALL legendre(pmin,pmax,mu-m,ABS(mu-m),theta,v_leg,error)
					
					! Check sull'errore in legendre
					error_leg_if: IF (error==1) THEN
								WRITE(*,*)
								WRITE(*,*) "Si e' verificato un errore nella subroutine fillblock_sub:"
								WRITE(*,*) "nel calcolo delle funzioni di Legendre, la subroutine si ferma!!!"
								WRITE(*,*)
								RETURN
					END IF error_leg_if
									
					!-------
					!Calcolo infine la mia funzione di Hankel
					CALL hankel1_d_sub(pmax,kr,v_h,error)
										
					!------
					!Calcolo del fattore numerico sotto sommatoria per Avt
					DO q=0,v_qmax(l)
						p=nr+nur-2.0D0*REAL(q,dbl)
						pi=INT(p,lo)
						v_Apmn(pi)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-p*(p+1.0D0)
					END DO
					
					!------
					!Calcolo della sommatoria e di Avt
					sommaA=(0.0D0,0.0D0)
					sommaAB0_do:DO q=0,v_qmax(l)
								p=nr+nur-2.0D0*REAL(q,dbl)
								pi=INT(p,lo)
								sommaA=sommaA+((0.0D0,1.0D0)**pi)*v_Apmn(pi)*v_aq_long(low_aq+q)*v_h(pi)*v_leg(pi)
				! 				WRITE(*,*) (0.0D0,1.0D0)**pi,v_leg(pi),v_h(pi),v_aq(q)
					END DO sommaAB0_do
				! 	WRITE(*,*)
					
					!Calcolo il mio coefficiente
					m_Aij(i,j)=v_c0(l)*coeff2*sommaA
					
					!Disalloco i miei vettore
					DEALLOCATE(v_h,v_leg,v_Apmn)					
					
					!Aggiorno l'indice e vedo eventualmente se uscire

					IF (l==elem2) EXIT n_do
					
					
! 					IF ((kr==1.0D0).AND.(theta==0.0D0).AND.(phi==0.0D0)) THEN
! 					WRITE(*,*) low_bq,up_bq
! 					END IF
					
					
					!Cambio della matrice se finisco quella precedente
					!Rientro ad inizio riga
					shift_if: IF (j==elem) THEN
								i=i+1
								j=0
					END IF shift_if
					
					
					
						!######################################################
				ELSE	!QUI COMINCIA IL CASO GENERALE Avt Bvt
						!######################################################
				
					!-------------------------
					!Calcolo componenti comuni
					!-------------------------

					
					!Aggiorno i bounds
					ib_if: IF (l==1) THEN
						low_bq=1
					ELSE
						low_bq=up_bq+1	
					END IF ib_if
					
					up_bq=low_bq+v_qqmax(l)-1

! 					WRITE(*,*) "Bounds bq 2",low_bq,up_bq,v_qqmax(l)
					
					! Calcolo Pmin e Pmax
					pmin=n+nu -2*v_qqmax(l)
					pmax=n+nu
					nbes=pmax-pmin+2
				
					!Calcolo coeff2
					coeff2=EXP(CMPLX(0.0D0,(mur-mr)*phi))
					
					!Adesso posso allocare tutti ma proprio tutti i miei vettori
					ALLOCATE(v_h(0:pmax+1),v_leg(pmin:pmax+1),&
						   & v_Apmn(pmin:pmax+1))
					
					!Inizializzo v_Apmn
					v_Apmn=0.0D0
					
					!-------
					!Calcolo la mia funzione di legendre
					CALL legendre(pmin,pmax+1,mu-m,ABS(mu-m),theta,v_leg,error)
					
					! Check sull'errore in legendre
					error_leg_if1: IF (error==1) THEN
								WRITE(*,*)
								WRITE(*,*) "Si e' verificato un errore nella subroutine fillblock_sub,caso generale:"
								WRITE(*,*) "nel calcolo delle funzioni di Legendre, la subroutine si ferma!!!"
								WRITE(*,*)
								RETURN
					END IF error_leg_if1
					
					!-------
					!Calcolo infine la mia funzione di Hankel
					CALL hankel1_d_sub(pmax+1,kr,v_h,error)
					
					
					!-----------
					!Calcolo Avt
					!-----------
					!Calcolo del fattore numerico sotto sommatoria per Avt
					DO q=0,v_qmax(l)
						p=nr+nur-2.0D0*REAL(q,dbl)
						pi=INT(p,lo)
						v_Apmn(pi)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-p*(p+1.0D0)
					END DO
					
					
					!Calcolo della sommatoria e di Avt
					sommaA=0
					sommaA_do:DO q=0,v_qmax(l)
								p=nr+nur-2.0D0*REAL(q,dbl)
								pi=INT(p,lo)
								sommaA=sommaA+((0.0D0,1.0D0)**pi)*v_Apmn(pi)*v_aq_long(low_aq+q)*v_h(pi)*v_leg(pi)
					END DO sommaA_do
					
					m_Aij(i,j)=v_c0(l)*coeff2*sommaA
					
					
					!-----------
					!Calcolo Bvt
					!-----------	
					
! 					IF ((kr==1.0D-2).AND.(theta==0.0D0).AND.(phi==0.0D0).AND.(m==3).AND.(n==4).AND.(mu==3).AND.(nu==3)) THEN
! 					 WRITE(*,*) v_h(pmin:pmax+1)
! 					 WRITE(*,*)
! 					END IF
					
					!Calcolo della sommatoria e di Bvt
					sommaB=0
					sommaB_do:DO q=1,v_qqmax(l)
								p=nr+nur-2.0D0*REAL(q,dbl)
								pi=INT(p,lo)
								sommaB=sommaB+((0.0D0,1.0D0)**(pi+1))*v_bq_long(low_bq+(q-1))*v_h(pi+1)*v_leg(pi+1)
! 						IF ((kr==1.0D0).AND.(theta==0.0D0).AND.(phi==0.0D0).AND.(m==2).AND.(n==2).AND.(mu==2).AND.(nu==2)) THEN
! 						WRITE(*,*) "Somma B2",sommaB,low_bq+(q-1)
! 						END IF
					END DO sommaB_do
					
					m_Bij(i,j)=v_c0(l)*coeff2*sommaB
					
					DEALLOCATE(v_h,v_leg,v_Apmn)	
					
					!Aggiorno l'indice e vedo eventualmente se uscire se ho finito il giro
					IF (l==elem2) EXIT n_do
					
					!Aggiorno i bounds
					
! 					IF ((kr==1.0D0).AND.(theta==0.0D0).AND.(phi==0.0D0)) THEN
! 					WRITE(*,*) "Bounds 2",low_bq,up_bq
! 					END IF
					
					!Shifto la colonna della matrice se finisco quella precedente
					shift_if1: IF (j==elem) THEN
								i=i+1
								j=0
					END IF shift_if1
					
				END IF Bvt_if
				
			END DO mu_do
		END DO nu_do		
	END DO m_do
END DO n_do

END SUBROUTINE fillblock_sub2






! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MACKOWSKI SUBROUTINES: metto qui tutte le subroutines di mackowski
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!*********************************************************************************
!1) SUBROUTINE clebsch_mk_sub: calcolo i coefficienti di gaunt tramite il formalismo di xu
!*********************************************************************************

SUBROUTINE clebsch_mk_sub(m,n,k,l,nmin,nmax,v_cg,error)

IMPLICIT NONE

! Dichiarazione argomenti
INTEGER(lo), INTENT(IN) :: m,n,k,l					! Indici coefficienti di di Clebsch-Gordan
INTEGER(lo), INTENT(IN) :: nmin,nmax					! Bound dell'indice
REAL(dbl), DIMENSION(nmin:nmax), INTENT(OUT) :: v_cg			! Vettore coefficienti di Gaunt
INTEGER(lo), INTENT(OUT) :: error						! Variabile di errore

! Dichiarazione variabili interno
INTEGER(lo) :: nhalf,w							! Indice switch ricorsioni e indice
INTEGER(lo) :: ierr								! Errore dummy per il logaritmo
REAL(dbl) :: Ac,Bc,Cc,Dc,Ec,Fc						! Coefficienti ricorsioni
REAL(dbl) :: p,p1,p2,pmin							! Coeff p per il calcolo coeff ricorsioni
REAL(dbl) :: num,den								! numeratore,denominatore
REAL(dbl) :: mr,nr,kr,lr,wr,nminr,nmaxr					! indici reali

!---------------------------------------------------------
! Subroutine vera e propria
!---------------------------------------------------------
error=0
ierr=0
nhalf=(nmin+nmax)/2


nmax_error_if: IF (nmax /= (n+l)) THEN
	WRITE(*,*)
	WRITE(*,*) "Si e' verificato un errore nella subroutine clebsch_mk_sub"
	WRITE(*,*) "nmax /= n+l"
	WRITE(*,*)
	error=1
	RETURN
END IF nmax_error_if

nmin_error_if: IF (nmin /= MAX(ABS(n-l),ABS(k+m))) THEN
	WRITE(*,*)
	WRITE(*,*) "Si e' verificato un errore nella subroutine clebsch_mk_sub"
	WRITE(*,*) "nmin /= MAX(|n-l|,|k+m|)"
	WRITE(*,*)
	error=1
	RETURN
END IF nmin_error_if


km_error_if: IF ((ABS(k)>l) .OR. (ABS(m)>n)) THEN
	WRITE(*,*)
	WRITE(*,*) "Si e' verificato un errore nella subroutine clebsch_mk_sub"
	WRITE(*,*) "k o m fuori dai bounds"
	WRITE(*,*)
	error=1
	RETURN
END IF km_error_if

mr=REAL(m,dbl)
nr=REAL(n,dbl)
kr=REAL(k,dbl)
lr=REAL(l,dbl)
nminr=REAL(nmin,dbl)
nmaxr=REAL(nmax,dbl)

!************
!PRIMO VALORE
!************
!Calcolo il valore per C(Nmin), i valori li calcolo man mano ed esco a seconda di Nmax-Nmin
nmin_if: IF (ABS(n-l)>=ABS(k+m)) THEN

	nl_if: IF (n>=l) THEN

		num=lnf(nr+mr,ierr)+lnf(nr-mr,ierr)+lnf(2.0D0*lr,ierr)+lnf(2.0D0*(nr-lr)+1.0D0,ierr)
		den=lnf(2.0D0*nr+1.0D0,ierr)+lnf(lr+kr,ierr)+lnf(lr-kr,ierr)+lnf(nr-lr+kr+mr,ierr)+lnf(nr-lr-kr-mr,ierr)
		v_cg(nmin)=((-1.0D0)**(k+l))*SQRT(EXP(num-den))

	ELSE

		num=lnf(lr+kr,ierr)+lnf(lr-kr,ierr)+lnf(2.0D0*nr,ierr)+lnf(2.0D0*(lr-nr)+1.0D0,ierr)
		den=lnf(2.0D0*lr+1.0D0,ierr)+lnf(nr+mr,ierr)+lnf(nr-mr,ierr)+lnf(lr-nr+mr+kr,ierr)+lnf(lr-nr-mr-kr,ierr)
		v_cg(nmin)=((-1.0D0)**(m+n))*SQRT(EXP(num-den))

	END IF nl_if

ELSE

	mk_if: IF ((m+k)>=0) THEN

		num=lnf(2.0D0*(kr+mr)+1.0D0,ierr)+lnf(nr+lr-kr-mr,ierr)+lnf(nr+mr,ierr)+lnf(lr+kr,ierr)
		den=lnf(nr+lr+kr+mr+1.0D0,ierr)+lnf(nr-lr+kr+mr,ierr)+lnf(lr-nr+kr+mr,ierr)+lnf(nr-mr,ierr)+lnf(lr-kr,ierr)
		v_cg(nmin)=((-1.0D0)**(m+n))*SQRT(EXP(num-den))

	ELSE

		num=lnf(-2.0D0*(kr+mr)+1.0D0,ierr)+lnf(nr+lr+kr+mr,ierr)+lnf(nr-mr,ierr)+lnf(lr-kr,ierr)
		den=lnf(nr+lr-kr-mr+1.0D0,ierr)+lnf(nr-lr-kr-mr,ierr)+lnf(lr-nr-kr-mr,ierr)+lnf(nr+mr,ierr)+lnf(lr+kr,ierr)
		v_cg(nmin)=((-1.0D0)**(k+l))*SQRT(EXP(num-den))

	END IF mk_if

END IF nmin_if

!Se nmax=nmin calcolo il primo sopra e ora esco
zero_if: IF ((nmax-nmin)==0) THEN
	RETURN
END IF zero_if


!**************
!SECONDO VALORE
!**************
!Calcolo il secondo valore, in direzione upward e fuori dal ciclo
wr=nminr+1.0D0
Ac=(4.0D0*(wr**2)*(2.0D0*wr+1.0D0)*(2.0D0*wr-1.0D0))/&					!Coeff A
  &((wr+kr+mr)*(wr-kr-mr)*(lr-nr+wr)*(nr-lr+wr)*(nr+lr-wr+1.0D0)*(nr+lr+wr+1.0D0))
Ac=SQRT(Ac)
!WRITE(*,*) m,n,k,l
!WRITE(*,*) "A", Ac

b_if: IF (n==l) THEN
	Bc=(mr-kr)/2.0D0
ELSE
	Bc=((mr-kr)*wr*(wr-1.0D0)-(kr+mr)*nr*(nr+1.0D0)+(kr+mr)*lr*(lr+1.0D0))/&		!Coeff B
	  &(2.0D0*wr*(wr-1.0D0))
END IF b_if
!WRITE(*,*) "B", Bc

v_cg(nmin+1)=Ac*Bc*v_cg(nmin)

!Se nmax=nmin+1 calcolo anche il secondo ed esco
one_if: IF ((nmax-nmin)==1) THEN
	RETURN
END IF one_if


!************
!TERZO VALORE
!************
!Calcolo il terzo valore, direzione backward e fuori dal ciclo
num=lnf(2.0D0*nr,ierr)+lnf(2.0D0*lr,ierr)+lnf(nr+lr+kr+mr,ierr)+lnf(nr+lr-kr-mr,ierr)
den=lnf(2.0D0*(nr+lr),ierr)+lnf(nr+mr,ierr)+lnf(nr-mr,ierr)+lnf(lr+kr,ierr)+lnf(lr-kr,ierr)
v_cg(nmax)=SQRT(EXP(num-den))

!Se nmax=nmin+2 calcolo anche il secondo ed esco
two_if: IF ((nmax-nmin)==2) THEN
	RETURN
END IF two_if


!*************
!QUARTO VALORE
!*************
wr=nmaxr-1.0D0
Dc=(4.0D0*((wr+1.0D0)**2)*(2.0D0*wr+1.0D0)*(2.0D0*wr+3.0D0))/&
  &((wr+kr+mr+1.0D0)*(wr-kr-mr+1.0D0)*(lr-nr+wr+1.0D0)*(nr-lr+wr+1.0D0)*(nr+lr-wr)*(nr+lr+wr+2.0D0))
Dc=SQRT(Dc)
!WRITE(*,*) m,n,k,l
!WRITE(*,*) "D", Dc

Ec=((mr-kr)*(wr+2.0D0)*(wr+1.0D0)-(kr+mr)*(nr+1.0D0)*nr+(kr+mr)*(lr+1.0D0)*lr)/&
  &(2.0D0*(wr+2.0D0)*(wr+1.0D0))
!WRITE(*,*) "E", Ec
v_cg(nmax-1)=Dc*Ec*v_cg(nmax)

!Se nmax=nmin+3 calcolo anche il secondo ed esco
three_if: IF ((nmax-nmin)==3) THEN
	RETURN
END IF three_if

!********************************
!CICLI UPWARD E DOWNWARD GENERALI
!********************************

!ciclo upward
up_do: DO w=nmin+2,nhalf

	wr=REAL(w,dbl)
	Ac=(4.0D0*(wr**2)*(2.0D0*wr+1.0D0)*(2.0D0*wr-1.0D0))/&					!Coeff A
	  &((wr+kr+mr)*(wr-kr-mr)*(lr-nr+wr)*(nr-lr+wr)*(nr+lr-wr+1.0D0)*(nr+lr+wr+1.0D0))
	Ac=SQRT(Ac)
!	WRITE(*,*) m,n,k,l,Ac

	b_if_loop: IF (n==l) THEN
		Bc=(mr-kr)/2.0D0
	ELSE
		Bc=((mr-kr)*wr*(wr-1.0D0)-(kr+mr)*nr*(nr+1.0D0)+(kr+mr)*lr*(lr+1.0D0))/&		!Coeff B
		  &(2.0D0*wr*(wr-1.0D0))
	END IF b_if_loop
!	WRITE(*,*) Bc

	Cc=((wr-kr-mr-1.0D0)*(wr+kr+mr-1.0D0)*(lr-nr+wr-1.0D0)*(nr-lr+wr-1.0D0)*(nr+lr-wr+2.0D0)*(nr+lr+wr))/&
	  &(4.0D0*((wr-1.0D0)**2)*(2.0D0*wr-3.0D0)*(2.0D0*wr-1.0D0))
	Cc=SQRT(Cc)
!	WRITE(*,*) Cc

	v_cg(w)=Ac*(Bc*v_cg(w-1)-Cc*v_cg(w-2))

END DO up_do


!ciclo downward
down_do: DO w=nmax-2,nhalf+1,-1

	wr=REAL(w,dbl)
	Dc=(4.0D0*((wr+1.0D0)**2)*(2.0D0*wr+1.0D0)*(2.0D0*wr+3.0D0))/&
	  &((wr+kr+mr+1.0D0)*(wr-kr-mr+1.0D0)*(lr-nr+wr+1.0D0)*(nr-lr+wr+1.0D0)*(nr+lr-wr)*(nr+lr+wr+2.0D0))
	Dc=SQRT(Dc)
!	WRITE(*,*) m,n,k,l,Dc

	Ec=((mr-kr)*(wr+2.0D0)*(wr+1.0D0)-(kr+mr)*(nr+1.0D0)*nr+(kr+mr)*(lr+1.0D0)*lr)/&
	  &(2.0D0*(wr+2.0D0)*(wr+1.0D0))
!	WRITE(*,*) Ec

	Fc=((wr-kr-mr+2.0D0)*(wr+kr+mr+2.0D0)*(lr-nr+wr+2.0D0)*(nr-lr+wr+2.0D0)*(nr+lr-wr-1.0D0)*(nr+lr+wr+3.0D0))/&
	  &(4.0D0*((wr+2.0D0)**2)*(2.0D0*wr+5.0D0)*(2.0D0*wr+3.0D0))
	Fc=SQRT(Fc)
!	WRITE(*,*) Fc


	v_cg(w)=Dc*(Ec*v_cg(w+1)-Fc*v_cg(w+2))

END DO down_do

END SUBROUTINE clebsch_mk_sub


!******************************************************************************
!2) SUBROUTINE fillblock_sparse_mk: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse di mackowski, anche per lo scattering
!******************************************************************************
SUBROUTINE fillblock_sparse_mk(nstop,v_cg,kr,v_Aij,v_Bij,v_Aij_sca,v_Bij_sca,v_jABij,v_iABij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop								! N espansioni multipolari
REAL(dbl), INTENT(IN) :: kr									! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_cg						! Vettore clebsch-gordan
INTEGER(lo), INTENT(OUT) :: error								! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Aij,v_Bij,v_Aij_sca,v_Bij_sca	! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij				! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: i,m,n,nu,w,low_nu,nmin,nmax			!Indici per il ciclo
REAL(dbl) :: nr,nur							!Indici reali 
INTEGER(lo) :: low_cg,up_cg						!Bounds per i vettori aq,bq
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca,c0		!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h			!Vett comp fun hankel

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!cg
low_cg=1
up_cg=1

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_h(0:(2*nstop+1)))

!Calcolo infine la mia funzione di Hankel
CALL hankel1_d_sub((2*nstop+1),kr,v_h,error)

n_do: DO n=1,nstop
	m_do: DO m=-n,n

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		nu_do: DO nu=low_nu,nstop


			!Rendo reali gli indici
			nr=REAL(n,dbl)
			nur=REAL(nu,dbl)

			!Calcolo qmax
			nmin=ABS(n-nu)
			nmax=n+nu

			!Aggiorno i
			i=i+1


			!Aggiorno i bounds
			ia_if: IF (i==1) THEN
				low_cg=1
			ELSE
				low_cg=up_cg+1
			END IF ia_if

			up_cg=low_cg+nmax-nmin


			!Calcolo della sommatoria e di Avt
			sommaA=(0.0D0,0.0D0)
			sommaAsca=(0.0D0,0.0D0)
			sommaA_do:DO w=nmin,nmax,2
				sommaA=sommaA+((0.0D0,1.0D0)**w)*v_cg(low_cg-nmin+w)*v_h(w)
				sommaAsca=sommaAsca+((0.0D0,1.0D0)**w)*v_cg(low_cg-nmin+w)*REAL(v_h(w),dbl)
			END DO sommaA_do

			c0=-((-1.0D0)**(m))*((0.0D0,1.0D0)**(n-nu))*SQRT((2.0D0*nr+1.0D0)*(2.0D0*nur+1.0D0))
!			c0=-((-1.0D0)**(m))*SQRT((2.0D0*nr+1.0D0)*(2.0D0*nur+1.0D0))
			v_Aij(i)=c0*sommaA
			v_Aij_sca(i)=c0*sommaAsca

			!-----------
			!Calcolo Bvt
			!-----------    
			!Calcolo della sommatoria e di Bvt
			sommaB=(0.0D0,0.0D0)
			sommaBsca=(0.0D0,0.0D0)
			sommaB_do:DO w=nmin+1,nmax-1,2
				sommaB=sommaB+((0.0D0,1.0D0)**w)*v_cg(low_cg-nmin+w)*v_h(w)
				sommaBsca=sommaBsca+((0.0D0,1.0D0)**w)*v_cg(low_cg-nmin+w)*REAL(v_h(w),dbl)
			END DO sommaB_do

			v_Bij(i)=c0*sommaB
			v_Bij_sca(i)=c0*sommaBsca

			!-----------
			!Aggiorno gli indici per lo storage sparse
			!-----------   
			v_jABij(i)=nu*(nu+1)+m

		END DO nu_do

		!-----------
		!Dove comincia la prossima riga?Qui!
		!-----------   
		v_iABij(n*(n+1)+m+1)=i+1

	END DO m_do
END DO n_do

END SUBROUTINE fillblock_sparse_mk


!******************************************************************************
!3) SUBROUTINE fill_AB_sparse_sca_mk:riempio tutta la struttura per i blocchi di 
! traslazione Aij e per Bij, formalismo di mackowski
!******************************************************************************
SUBROUTINE fill_AB_sparse_sca_mk(ns,nstop,k,m_xyz,v_cg,v_jBlock,v_iBlock,m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop								! N sfere, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k										! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz							! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock					! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_cg							! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca		! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij				! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

	lowb=v_iBlock(i)
	upb=v_iBlock(i+1)-1

	col_do: DO next=lowb,upb

		!Adesso so sia riga che colonna
		j=v_jBlock(next)

		!Calcolo le coordinate sferiche relative
		xij=m_xyz(i,1)-m_xyz(j,1)
		yij=m_xyz(i,2)-m_xyz(j,2)
		zij =m_xyz(i,3)-m_xyz(j,3)
		CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

		!Chiamo la subroutine per riempire Dij
		CALL fillblock_sparse_mk(nstop,v_cg,k*rij,m_Aij(:,next),m_Bij(:,next), &  
					& m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)


	END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_sca_mk


!******************************************************************************
!4) SUBROUTINE D_k_mn: calcolo gli elementi ridotti di matrice secondo mackowski
!   D(beta)^{k}_{m,n}: la ricorrenza e' in n.In sostanza e' come per Edmonds,
!   ma cambio il segno con (-1)^{m+k}  
!******************************************************************************
SUBROUTINE D_k_mn(nstop,m,k,beta,v_d,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,m,k                  ! Numero di espansioni multipolari, d^{nstop}_{m,k}(beta)
REAL(dbl), INTENT(IN) :: beta                           ! Angolo Beta di Eulero
INTEGER(lo), INTENT(OUT) :: error                     ! Flag di errore
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_d       ! Vettore di output valori funzioni pi_mn

!Dichiarazione variabili interne 
INTEGER(lo) :: ierr                                   ! altra flag di errore
INTEGER(lo) :: nmin,modm,modk,MpK,MmK,emk,n           ! Indice per i miei calcoli
REAL(dbl) :: logw,w,mr,kr,nr,x                          ! Fattoriale e indici reali, x=cos(beta)
REAL(dbl) :: c0,c1,c2                                   !Coefficienti ricorrenza


!Subroutine vera e propria
error=0

!Calcolo nmin
modm=ABS(m)
modk=ABS(k)
nmin=MAX(modm,modk)

!Controllo se nstop>nmin 
nm_check: IF (nstop<nmin) THEN
            WRITE(*,*) "nstop<nmin, la procedura d_nmk si ferma!"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
thetain_check: IF ((beta>PI_D) .OR. (beta<-PI_D)) THEN
                WRITE(*,*) "Beta e' al di fuori dei suoi bounds, la procedura pi_mn si ferma"
                error=1
                RETURN
END IF thetain_check

!Inizializzo il vettore e alcune grandezze che fanno comodo
v_d=0.0D0       !Vettore

emk_if: IF (k>=m) THEN !Emk
    emk=1.0D0
ELSE
    emk=(-1.0D0)**(m-k)
END IF emk_if

MpK=ABS(m+k)    !moduli
MmK=ABS(m-k)

x=COS(beta)     !Coseno di beta

mr=REAL(m,dbl)  !Diventato reali
kr=REAL(k,dbl)

!If per vedere se m=k=0 come caso particolare
nmin_zero_if: IF (nmin==0) THEN

        !Primo valore
        v_d(nmin)=1
        
        !Se nstop>0 assegno altri valori
        nstop_zero_if: IF (nstop>0) THEN
            
            !Secondo valore
            v_d(nmin+1)=x
            
            !Ciclo su n
            nstop_zero_do: DO n=2,nstop
                nr=REAL(n,dbl)
                v_d(n)=((2.0D0*nr-1.0D0)*x*v_d(n-1)-(nr-1.0D0)*v_d(n-2))/nr
            END DO  nstop_zero_do

        END IF nstop_zero_if

ELSE

    !Calcolo il primo valore della ricorrenza
    logw=lnf(2.0D0*REAL(nmin,dbl),ierr) - lnf(REAL(MmK,dbl),ierr) -lnf(REAL(MpK,dbl),ierr)
    w=EXP(logw)
    v_d(nmin)=emk*(2.0D0**(-nmin))*SQRT(w)*((SQRT(1-x))**MmK)*((SQRT(1+x))**MpK)

    !Ciclo per il calcolo
    nstop_do: DO n=nmin+1,nstop
    
                !Coefficienti della ricorrenza
                nr=REAL(n,dbl)
                c0=1.0D0/((nr-1.0D0)*SQRT(nr**2-mr**2)*SQRT(nr**2-kr**2))
                c1=(2.0D0*nr-1.0D0)*(nr*(nr-1.0D0)*x-mr*kr)
                c2=nr*SQRT((nr-1.0D0)**2-mr**2)*SQRT((nr-1.0D0)**2-kr**2)
                
                !Calcolo
                v_d(n)=(c1*v_d(n-1)-c2*v_d(n-2))*c0
                
    END DO nstop_do

END IF nmin_zero_if

END SUBROUTINE D_k_mn




!******************************************************************************
!5) SUBROUTINE fillblock_Dkmn:riempio un blocco ij dei coefficienti Dkmn
!secondo il formalismo di mackowski,che poi e' tutto quasi uguale
!******************************************************************************
SUBROUTINE fillblock_Dkmn_mk(nstop,theta,v_Dkmn,v_jDkmn,v_iDkmn,error)

IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop                               ! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: theta                                   ! Angolo theta coordinate relative
INTEGER(lo), INTENT(OUT) :: error                              ! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_Dkmn                   ! Vettore sparse di output valori Dkmn
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jDkmn,v_iDkmn      ! Vettori indici colonne e righe

!Dichiarazione variabili interne 
INTEGER(lo) :: i,n,m,k,lowb,upb,nrow,ncol,modm,maxmk      !Indici
INTEGER(lo) :: imk,ikm,i_k_m,i_m_k,imm,i_m_m,im_m,i_mm    !Indici blocco
REAL(dbl), DIMENSION(0:nstop) :: v_d,v_d1                   ! Vettori fattori Dkmn


!Subroutine vera e propria
error=0
v_d=0.0D0
v_d1=0.0D0

!Check di errore su theta
thetain_check: IF ((theta>PI_D) .OR. (theta<-PI_D)) THEN
	WRITE(*,*) "Theta e' al di fuori dei suoi bounds, la procedura fillblock_Dkmn_mk si ferma"
	error=1
	RETURN
END IF thetain_check

!Inizializzo degli indici
i=0
v_iDkmn(1)=1
!v_Dkmn=0.0D0

!Assegno gli indici di righe e colonne prima di tutto
!Lo faccio cosi perche' so gia di avere i blocchi quadrati diagonali
block_do: DO n=1,nstop

	!Bounds
	lowb=n**2
	upb=n*(n+2)

	row_do: DO nrow=lowb,upb
		col_do: DO ncol=lowb,upb

			i=i+1
			v_jDkmn(i)=ncol

		END DO col_do

		v_iDkmn(nrow+1)=i+1

	END DO row_do
END DO block_do


!Assegno i valori lungo la diagonale per m=k=0
CALL D_k_mn(nstop,0,0,theta,v_d,error)

diag_zero_do: DO n=1,nstop
	imm=id(n,0,0)
	v_Dkmn(imm)=v_d(n)
END DO diag_zero_do


!Calcolo i valori che stanno sulle due diagonali principali:simmetria di 2
diag_do: DO m=1,nstop

	!Chiamata per ciascuna delle due diagonali
	CALL D_k_mn(nstop,m,m,theta,v_d,error)
	CALL D_k_mn(nstop,m,-m,theta,v_d1,error)

	!Faccio le assegnazioni saltando di blocco in blocco
	n_diag_do: DO n=m,nstop

		!Calcolo gli indici
		imm=id(n,m,m)
		i_m_m=id(n,-m,-m)
		im_m=id(n,m,-m)
		i_mm=id(n,-m,m)

		!Assegno i valori
		v_Dkmn(imm)=v_d(n)
		v_Dkmn(i_m_m)=v_d(n)
		v_Dkmn(im_m)=v_d1(n)
		v_Dkmn(i_mm)=v_d1(n)

	END DO n_diag_do
END DO diag_do


!Calcolo il caso generale
m_do: DO m=-(nstop-1),(nstop-1)
	k_do: DO k=ABS(m)+1,nstop

		CALL D_k_mn(nstop,m,k,theta,v_d,error)
		modm=ABS(m)
		maxmk=MAX(ABS(m),k)

		!Faccio le assegnazioni saltando di blocco in blocco
		n_do: DO n=maxmk,nstop

			!Sempre indici
			imk=id(n,m,k)
			ikm=id(n,k,m)
			i_k_m=id(n,-k,-m)
			i_m_k=id(n,-m,-k)

			!Sempre valori
			v_Dkmn(ikm)=v_d(n)
			v_Dkmn(i_k_m)=((-1.0D0)**(m-k))*v_d(n)
			v_Dkmn(i_m_k)=v_d(n)
			v_Dkmn(imk)=v_Dkmn(i_k_m)

		END DO n_do
	END DO k_do
END DO m_do

END SUBROUTINE fillblock_Dkmn_mk



!******************************************************************************
!6) SUBROUTINE fill_D_PHI_sparse_mk:riempio tutta la struttura per i blocchi di 
! rotazione Dij e per Exp[phi_ij]
!******************************************************************************
SUBROUTINE fill_D_PHI_sparse_mk(ns,nstop,nnb,m_xyz,v_jBlock,v_iBlock,m_Dnkm,m_jDnkm,m_iDnkm,m_exphi)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ns,nstop,nnb                           ! N sfere, N multipolar exp e N nnz blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz                      ! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock,v_iBlock         ! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Dnkm                    ! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jDnkm,m_iDnkm       ! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop,1:nnb), INTENT(OUT) :: m_exphi !Matrice Esponenziali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,rij,thetaij,phiij

! Funzione vera e propria
row_do: DO i=1,ns

    lowb=v_iBlock(i)
    upb=v_iBlock(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock(next)

        !Calcolo le coordinate sferiche relative
        xij=m_xyz(i,1)-m_xyz(j,1)
        yij=m_xyz(i,2)-m_xyz(j,2)
        zij =m_xyz(i,3)-m_xyz(j,3)
        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)
        
        !Chiamo la subroutine per riempire Dij
        CALL fillblock_Dkmn_mk(nstop,thetaij,m_Dnkm(:,next),m_jDnkm(:,next),m_iDnkm(:,next),error)

        !Riempio anche la colonna degli esponenziali
        phi_do: DO m=-nstop,nstop

            m_exphi(m,next)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )

        END DO phi_do

    END DO col_do

END DO row_do

END SUBROUTINE fill_D_PHI_sparse_mk




!******************************************************************************
!7) SUBROUTINE D_nm12_mk: calcolo gli elementi ridotti di matrice per k=1 e k=-1, in
! maniera da poter calcolare la sezione di estinzione
!******************************************************************************
SUBROUTINE D_nm12_mk(nstop,beta,v_dnm1,v_dnm2,error)


IMPLICIT NONE

!Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: beta							! Theta di cos(theta)
INTEGER(lo), INTENT(OUT) :: error						! Flag di errore
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_dnm1,v_dnm2	! Vettore di output valori funzioni d_mn1

!Dichiarazione variabili interne 
INTEGER(lo) :: mm,m,n,beg,indx						! Indice per riordinare i miei calcoli e non
REAL(dbl) :: mr,nr										! Semifattoriale e m,n reali,variabile segno
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_d1,v_d2		! Vettore temporaneo red rot mat el

!Subroutine vera e propria
error=0

!Controllo se nstop>0 
nm_check: IF (nstop<0) THEN
			WRITE(*,*) "nstop<0, la procedura d_nm2 si ferma"
            error=1
            RETURN
END IF nm_check

!Controllo se thetain e' nei suoi limiti
betain_check: IF ((beta>PI_D*1.001) .OR. (beta<-PI_D*1.001)) THEN
				WRITE(*,*) "Beta e' al di fuori dei suoi bounds, la procedura d_nm2 si ferma"
                error=1
                RETURN
END IF betain_check

!Inizializzo a zero tutto il vettore
v_dnm1=0.0D0
v_dnm2=0.0D0

!Riempio il vettore v_legtot
d_mdo: DO m=-nstop,nstop

	mm=ABS(m)
	mr=REAL(m,dbl)

	ALLOCATE(v_d1(0:nstop),v_d2(0:nstop))
	CALL D_k_mn(nstop,m,1,beta,v_d1,error)
	CALL D_k_mn(nstop,m,-1,beta,v_d2,error)

	!If per l'errore su legendre
	errd_if: IF (error/=0) THEN
		WRITE(*,*) "Errore in d_nm12_mk, la procedura si arresta"
		RETURN
	END IF errd_if

	!If per lo starting index
	begin_if: IF (mm==0) THEN
		beg=mm+1
	ELSE
		beg=mm
	END IF begin_if
	
	!Riempio il vettore finale
	d_ndo: DO n=beg,nstop
	
		nr=REAL(n,dbl)		
		indx=n*(n+1)+m
		v_dnm1(indx)=v_d1(n)
		v_dnm2(indx)=v_d2(n)
	
	END DO d_ndo
		
	DEALLOCATE(v_d1,v_d2)

END DO d_mdo

END SUBROUTINE D_nm12_mk






!******************************************************************************
!8) SUBROUTINE field_expRandom_sub: calcolo i coefficienti del campo incidente
! per una direzione arbitraria rispetto al sistema di riferimento del cluster che e'
! ruotato degli angoli alpha bete gamma.
!******************************************************************************
SUBROUTINE field_expRandom_sub_mk(nstop,ns,k,betap,m_xyz,alpha,beta,gamma,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k,betap							!Vettore d'onda e polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Matrice coordinate
REAL(dbl), INTENT(IN) :: alpha,beta,gamma					!Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p				!Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj,n,m,dimens								!Indici e dimensioni
REAL(dbl) :: nr,mr,x,y,z,d,theta,phi							!Indici reali e coordinate
REAL(dbl) :: kd													!Phase 
COMPLEX(dbl) :: pshift,expg1,expg2								!Phase shift				
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2			!Elementi ridotti di matrice di rotazione


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))

!Funzione angolare Pi_mn
CALL d_nm12_mk(nstop,beta,v_dnm1,v_dnm2,error)

dnm_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da field_expRandom_sub: il programma termina ora...") 
       			STOP
END IF dnm_if

!Calcolo gli esponenziali
expg1=EXP(gamma*(0.0D0,1.0D0))
expg2=EXP(-gamma*(0.0D0,1.0D0))

!Inizializzo j e i vettori output
j=0
jj=0
v_p=(0.0D0,0.0D0)

i_do: DO i=1,ns

	jj=0

	!Assegno x y e z
	x=m_xyz(i,1)
	y=m_xyz(i,2)
	z=m_xyz(i,3)
	
	!Calcolo phase e phase shift
	kd=k*(SIN(beta)*(x*COS(alpha)+y*SIN(alpha)) + z*COS(beta))
	pshift=EXP(CMPLX(0.0D0,kd))

	n_do: DO n=1,nstop
	
		!n reale
		nr=REAL(n,dbl)
	
		m_do: DO m=-n,n
		
			!m reale
			mr=REAL(m,dbl)	
		
			!Incremento j
			j=j+1
			jj=jj+1
		
			v_p(2*j-1)=-((0.0D0,1.0D0)**(n+1))*((-1.0D0)**m)*0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
			v_p(2*j)=-((0.0D0,1.0D0)**(n+1))*((-1.0D0)**m)*0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
					& ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )

!			v_p(2*j-1)=((-1.0D0)**m)*0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
!					& ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
!			v_p(2*j)=((-1.0D0)**m)*0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
!					& ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )

		END DO m_do
	END DO n_do
END DO i_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE field_expRandom_sub_mk

!****************************************************************************************
!2bis) SUBROUTINE fill_cg_sparse:calcolo i vettori, tutti dritti, di coefficienti di
! clebsch-gordan nel caso si abbia mu=-m
!****************************************************************************************
SUBROUTINE fill_cg_sparse(nstop,v_cg_long)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop					! N espansioni multipolari e bounds vettori
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cg_long		! Vettore per tutti i coeff di gaunt

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,error,nmin,nmax
INTEGER(lo) :: low_cg,up_cg,low_nu
REAL(dbl) :: mr,nr,mur,nur
REAL(dbl), DIMENSION(:), ALLOCATABLE :: v_cg,v_cg_m1

! Funzione vera e propria

!Inizializzo indici e bounds

!cg
error=0
i=0
low_cg=1
up_cg=1

n_do: DO n=1,nstop
    m_do: DO m=-n,n
        
        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop

                !Rendo reali gli indici
                mr=REAL(m,dbl)
                nr=REAL(n,dbl)
                mur=REAL(m,dbl)
                nur=REAL(nu,dbl)

                !Calcolo nmin,nmax
                nmin=ABS(n-nu)
                nmax=n+nu

                !---------------------------------
                !Parte di storage per aq
                !---------------------------------

                !Aggiorno i
                i=i+1

                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_cg=1
                ELSE
                    low_cg=up_cg+1  
                END IF ia_if
                
                up_cg=low_cg+nmax-nmin

!                WRITE(*,*) "Bounds aq",low_aq,up_aq

                !Alloco la memoria per il mio v_aq
                ALLOCATE(v_cg(nmin:nmax)) 

                !Archivio i dati nel vettore lungo
                CALL clebsch_mk_sub(-m,n,m,nu,nmin,nmax,v_cg,error)

                error_if: IF (error/=0) THEN
                  	WRITE(*,*)
                  	WRITE(*,*) "Si e' verificato un errore nella subroutine fill_cg_sparse"
                  	WRITE(*,*)
                  	STOP
                END IF error_if

                v_cg_long(low_cg:up_cg)=v_cg

                !Disalloco v_cg
                DEALLOCATE(v_cg)

        END DO nu_do
    END DO m_do
END DO n_do

END SUBROUTINE fill_cg_sparse




!****************************************************************************************
!2bis) SUBROUTINE fill_cg_sparse_norm:calcolo i vettori, tutti dritti, di coefficienti di
! clebsch-gordan nel caso si abbia mu=-m, e poi tutto e' moltiplicato per i coeff di cg
! per m=-1,mu=1, in maniera tale da essere piu' svelti a riempire la matrice. qui lo faccio
! malissimo con un numero assurdo di chiamate,ma che palle,va bene cosi!
!****************************************************************************************
SUBROUTINE fill_cg_sparse_norm(nstop,v_cg_long)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop					! N espansioni multipolari e bounds vettori
REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cg_long		! Vettore per tutti i coeff di gaunt

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,mu,nu,i,error,nmin,nmax
INTEGER(lo) :: low_cg,up_cg,low_nu
REAL(dbl) :: mr,nr,mur,nur
REAL(dbl), DIMENSION(:), ALLOCATABLE :: v_cg,v_cg_norm

! Funzione vera e propria

!Inizializzo indici e bounds

!cg
error=0
i=0
low_cg=1
up_cg=1

n_do: DO n=1,nstop
    m_do: DO m=-n,n
        
        !Calcolo il lower bounds per il ciclo successivo
        IF (m==0) THEN
            low_nu=1
        ELSE
            low_nu=ABS(m)
        END IF

        nu_do: DO nu=low_nu,nstop

                !Calcolo nmin,nmax
                nmin=ABS(n-nu)
                nmax=n+nu

                !---------------------------------
                !Parte di storage per aq
                !---------------------------------

                !Aggiorno i
                i=i+1

                !Aggiorno i bounds
                ia_if: IF (i==1) THEN
                    low_cg=1
                ELSE
                    low_cg=up_cg+1  
                END IF ia_if
                
                up_cg=low_cg+nmax-nmin

!                WRITE(*,*) "Bounds aq",low_aq,up_aq

                !Alloco la memoria per il mio v_aq
                ALLOCATE(v_cg(nmin:nmax),v_cg_norm(nmin:nmax)) 

                !Archivio i dati nel vettore lungo
                CALL clebsch_mk_sub(-1,n,1,nu,nmin,nmax,v_cg_norm,error)
                CALL clebsch_mk_sub(-m,n,m,nu,nmin,nmax,v_cg,error)

                error_if: IF (error/=0) THEN
                  	WRITE(*,*)
                  	WRITE(*,*) "Si e' verificato un errore nella subroutine fill_cg_sparse"
                  	WRITE(*,*)
                  	STOP
                END IF error_if

                v_cg_long(low_cg:up_cg)=v_cg*v_cg_norm

                !Disalloco v_cg
                DEALLOCATE(v_cg,v_cg_norm)

        END DO nu_do
    END DO m_do
END DO n_do

END SUBROUTINE fill_cg_sparse_norm

!******************************************************************************
!6bis) SUBROUTINE cext_random_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cext_random_sub_mk(nstop,ns,k,v_p,v_ab,v_cext)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k									!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p,v_ab			!Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cext				!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j
COMPLEX(dbl) :: somma

!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma+CONJG(v_p(2*j-1))*v_ab(2*j-1)+CONJG(v_p(2*j))*v_ab(2*j)
		
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_cext(i)=-(4*Pi_d/(k**2))*REAL(somma,dbl)
	
END DO sphere_do

END SUBROUTINE cext_random_sub_mk


!******************************************************************************
!6bis) SUBROUTINE cext_random_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE cabs_random_sub_mk(nstop,ns,k,v_ab,v_cabs)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns				!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k						!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab			!Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs			!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j
COMPLEX(dbl) :: somma

!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)
		
	n_do: DO n=1,nstop
	
		m_do: DO m=-n,n
		
			!Aggiorno l'indice
			j=j+1
				
			!Costruisco la somma
			somma=somma+(1.0D0+1.0D0/(CONJG(m_a(n,v_patt(i)))))*CONJG(v_ab(2*j-1))*v_ab(2*j-1)+&
				 &    (1.0D0+1.0D0/(CONJG(m_b(n,v_patt(i)))))*CONJG(v_ab(2*j))*v_ab(2*j)
		
		END DO m_do
	
	END DO n_do
	
	!Calcolo la funzione
	v_cabs(i)=-(4*Pi_d/(k**2))*REAL(somma,dbl)
	
END DO sphere_do

END SUBROUTINE cabs_random_sub_mk








!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!SUBROUTINES FOR THE SINGLE SHELL CODE
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================




!******************************************************************************
!1) SUBROUTINE fill_TU_sparse: riempio un blocco delle matrici Aij e Bij, TUij ma
!seguendo la data structure usata di solito, quindi dividendo in blocchi
!******************************************************************************
SUBROUTINE fill_TU_sparse(nstop,v_c0,v_aq,v_bq,m_t1,m_u1,m_t2,m_u2,kr,&
				& v_precond_shell,v_sAij,v_sBij,v_sT1ij,v_sT2ij,v_sU1ij,v_sU2ij,v_jTUij,v_iTUij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop									! N espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: kr									! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq						! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_t1,m_u1,m_t2,m_u2				! Matrici ausiliarie per sistema lineare

INTEGER(lo), INTENT(OUT) :: error									! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_precond_shell					! Vettore preconditioning
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sT1ij,v_sT2ij,v_sU1ij,v_sU2ij		! Vettori blocco sparse ij delle matrici A e B
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sAij,v_sBij					! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jTUij,v_iTUij					! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,j,q,p,low_nu,qmax							!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr										!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq							!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes									!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca						!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_j							!Vett comp fun besselj per argomento complesso
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn							!Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iTUij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_j(0:(2*nstop+1)))

!Calcolo infine la mia funzione di Hankel
CALL besselj_z_sub((2*nstop+1),kr,v_j,error)

j=0
n_do: DO n=1,nstop
	m_do: DO m=-n,n

		!Aggiorno l'indice per il vettore di preconditioning
		j=j+1

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		nu_do: DO nu=low_nu,nstop


			!Rendo reali gli indici
			mr=REAL(m,dbl)
			nr=REAL(n,dbl)
			mur=REAL(m,dbl)
			nur=REAL(nu,dbl)

			!Calcolo qmax
			qmax=MIN(n,nu)

			!Aggiorno i
			i=i+1


			!Aggiorno i bounds
			ia_if: IF (i==1) THEN
				low_aq=1
			ELSE
				low_aq=up_aq+1  
			END IF ia_if

			ib_if: IF (i==1) THEN
				low_bq=1
			ELSE
				low_bq=up_bq+1  
			END IF ib_if

			up_aq=low_aq+qmax

			up_bq=low_bq+qmax-1

			!Calcolo Pmin e Pmax
			pmin=n+nu -2*qmax
			pmax=n+nu
			nbes=pmax-pmin+2

			!Adesso posso allocare tutti ma proprio tutti i miei vettori
			ALLOCATE(v_Apmn(pmin:pmax+1))

			!Inizializzo v_Apmn
			v_Apmn=0.0D0


			!-----------
			!Calcolo Avt
			!-----------
			!Calcolo del fattore numerico sotto sommatoria per Avt
			DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
			END DO


			!Calcolo della sommatoria e di Avt
			sommaA=(0.0D0,0.0D0)
			sommaAsca=(0.0D0,0.0D0)
			sommaA_do:DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_j(p)
			END DO sommaA_do

			v_sAij(i)=v_c0(i)*sommaA
			v_sT1ij(i)=v_c0(i)*sommaA*m_t1(n,nu)
			v_sU2ij(i)=v_c0(i)*sommaA*m_u2(n,nu)

			!-----------
			!Calcolo Bvt
			!-----------
			!Calcolo della sommatoria e di Bvt
			sommaB=(0.0D0,0.0D0)
			sommaBsca=(0.0D0,0.0D0)
			sommaB_do:DO q=1,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_j(p+1)
			END DO sommaB_do

			v_sBij(i)=v_c0(i)*sommaB
			v_sT2ij(i)=v_c0(i)*sommaB*m_t2(n,nu)
			v_sU1ij(i)=v_c0(i)*sommaB*m_u1(n,nu)

			DEALLOCATE(v_Apmn)

			!-------------------------------------
			!Riempio il vettore di preconditioning
			!-------------------------------------
			precond_if: IF (n==nu) THEN
				v_precond_shell(2*j-1)=(1.0D0,0.0D0)/v_sU2ij(i)
				v_precond_shell(2*j)=(1.0D0,0.0D0)/v_sT1ij(i)
			END IF precond_if
				

			!-----------------------------------------
			!Aggiorno gli indici per lo storage sparse
			!-----------------------------------------
			v_jTUij(i)=nu*(nu+1)+m

		END DO nu_do

		!-----------------------------------
		!Dove comincia la prossima riga?Qui!
		!-----------------------------------
		v_iTUij(n*(n+1)+m+1)=i+1

	END DO m_do
END DO n_do

END SUBROUTINE fill_TU_sparse



!******************************************************************************************
!2) SUBROUTINE fill_TU_sparse_pardiso: riempio un blocco delle matrici Aij e Bij, TUij ma
!seguendo la data structure nuova o meglio, ora che posso, scrivo TU come matrice unica in
!fashion compressed sparse row, e poi scrivo Aij e Bij, che mi servono comunque, come al
!solito
!*****************************************************************************************
SUBROUTINE fill_TU_sparse_pardiso(nstop,v_c0,v_aq,v_bq,m_t1,m_u1,m_t2,m_u2,kr,&
				&v_sAij,v_sBij,v_jABij,v_iABij,v_sTUij,v_jTUij,v_iTUij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop									! N espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: kr									! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq						! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_t1,m_u1,m_t2,m_u2				! Matrici ausiliarie per sistema lineare

INTEGER(lo), INTENT(OUT) :: error									! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sAij,v_sBij					! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij					! Vettori sparse colonne e righe
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sTUij    					! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jTUij,v_iTUij					! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,i1,i2,j,q,p,low_nu,qmax,m1,m2					!Indici per il ciclo
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq							!Bounds per i vettori aq,bq
REAL(dbl) :: mr,nr,mur,nur,pr										!Indici reali 
INTEGER(lo) :: pmin,pmax,nbes									!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAsca,sommaBsca						!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_j							!Vett comp fun besselj per argomento complesso
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn							!Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_j(0:(2*nstop+1)))

!Calcolo infine la mia funzione di Hankel
CALL besselj_z_sub((2*nstop+1),kr,v_j,error)


!------------------------------------------------------------------------
!Questo e' il ciclo per riempire le matrici A e B che mi servono comunque
!------------------------------------------------------------------------
j=0
AB_n_do: DO n=1,nstop
	AB_m_do: DO m=-n,n

		!Aggiorno l'indice per il vettore di preconditioning
		j=j+1

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		AB_nu_do: DO nu=low_nu,nstop


			!Rendo reali gli indici
			mr=REAL(m,dbl)
			nr=REAL(n,dbl)
			mur=REAL(m,dbl)
			nur=REAL(nu,dbl)

			!Calcolo qmax
			qmax=MIN(n,nu)

			!Aggiorno i
			i=i+1

			!Aggiorno i bounds
			ia_if: IF (i==1) THEN
				low_aq=1
			ELSE
				low_aq=up_aq+1  
			END IF ia_if

			ib_if: IF (i==1) THEN
				low_bq=1
			ELSE
				low_bq=up_bq+1  
			END IF ib_if

			up_aq=low_aq+qmax

			up_bq=low_bq+qmax-1

			!Calcolo Pmin e Pmax
			pmin=n+nu -2*qmax
			pmax=n+nu
			nbes=pmax-pmin+2

			!Adesso posso allocare tutti ma proprio tutti i miei vettori
			ALLOCATE(v_Apmn(pmin:pmax+1))

			!Inizializzo v_Apmn
			v_Apmn=0.0D0


			!-----------
			!Calcolo Avt
			!-----------
			!Calcolo del fattore numerico sotto sommatoria per Avt
			DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
			END DO


			!Calcolo della sommatoria e di Avt
			sommaA=(0.0D0,0.0D0)
			sommaAsca=(0.0D0,0.0D0)
			sommaA_do:DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_j(p)
			END DO sommaA_do

			v_sAij(i)=v_c0(i)*sommaA

			!-----------
			!Calcolo Bvt
			!-----------
			!Calcolo della sommatoria e di Bvt
			sommaB=(0.0D0,0.0D0)
			sommaBsca=(0.0D0,0.0D0)
			sommaB_do:DO q=1,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_j(p+1)
			END DO sommaB_do

			v_sBij(i)=v_c0(i)*sommaB

			DEALLOCATE(v_Apmn)

			!-----------------------------------------
			!Aggiorno gli indici per lo storage sparse
			!-----------------------------------------
			v_jABij(i)=nu*(nu+1)+m

		END DO AB_nu_do

		!-----------------------------------
		!Dove comincia la prossima riga?Qui!
		!-----------------------------------
		v_iABij(n*(n+1)+m+1)=i+1

	END DO AB_m_do
END DO AB_n_do





!------------------------------------------------------------------------
!Questo e' il ciclo per riempire la matrici TU, quindi ho anche i cicli
!su m1 e m2
!------------------------------------------------------------------------
!Indici
i=0			!Indice per TU
i1=0			!indice per A e B
i2=0			!indice per A e B
v_iTUij(1)=1

n_do: DO n=1,nstop
	m_do: DO m=-n,n

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		m1_do: DO m1=1,2

			nu_do: DO nu=low_nu,nstop


				!Aggiorno gli indici
				i=i+2						!Salto di due in due per quello di TU


				!Assegno i valori per la matrice TU, con un if su m1
				m1_if: IF (m1==1) THEN
					i1=i1+1
					v_sTUij(i-1)=v_sAij(i1)*m_u2(n,nu)
					v_sTUij(i)=v_sBij(i1)*m_t2(n,nu)
					! WRITE(*,1004) v_sAij(i1),v_sBij(i1)
					! 1004 FORMAT (4ES14.5)
				ELSE
					i2=i2+1
					v_sTUij(i-1)=v_sBij(i2)*m_u1(n,nu)
					v_sTUij(i)=v_sAij(i2)*m_t1(n,nu)
					! WRITE(*,1005) v_sBij(i2),v_sAij(i2)
					! 1005 FORMAT (4ES14.5)
				END IF m1_if


				!-----------------------------------------
				!Aggiorno gli indici per lo storage sparse
				!-----------------------------------------
				v_jTUij(i-1)=2*(nu*(nu+1)+m-1)+1
				v_jTUij(i)=2*(nu*(nu+1)+m-1)+2

			END DO nu_do

			!-----------------------------------
			!Dove comincia la prossima riga?Qui!
			!-----------------------------------
			v_iTUij(2*(n*(n+1)+m-1)+m1+1)=i+1

		END DO m1_do
	END DO m_do
END DO n_do


END SUBROUTINE fill_TU_sparse_pardiso


!******************************************************************************
!3quater) SUBROUTINE fillblock_sparse_sca: riempio un blocco delle matrici Aij e Bij,ma
!utilizzando il formalismo sparse, nel caso dei conti con dipoli.
!******************************************************************************
SUBROUTINE fill_TU_sparse_pardiso_fast(nstop,v_c0,v_aq,v_bq,m_t1,m_u1,m_t2,m_u2,&
				     & kr,m_index,m_index_TU,m_Apmn,v_sAij,v_sBij,v_sTUij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: kr								! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Vettore norm,gaunt e bq
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index,m_index_TU				! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sAij,v_sBij,v_sTUij			! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_t1,m_u1,m_t2,m_u2				! Matrici ausiliarie per sistema lineare

! Dichiarazione variabili interne
INTEGER(lo) :: m,m1,n,nu,i,iab,q,p,low_nu,qmax			!Indici per il ciclo
REAL(dbl) :: mr,nr,mur,nur,pr					!Indici reali 
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq			!Bounds per i vettori aq,bq
INTEGER(lo) :: pmin,pmax,nbes,NNZab				!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB					!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_j			!Vett comp fun Bessel
INTEGER(lo) :: chunk=1					!Thread id, number of total threads and chunk size

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
NNZab=nstop*(1+2*nstop*(3+nstop))	!Computing NNZab,so i can have the cycle on the whole number of filled matrix values
NNZab=NNZab/3

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of Bessel functions
ALLOCATE(v_j(0:(2*nstop+1)))
CALL besselj_z_sub((2*nstop+1),kr,v_j,error)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!This is a flattened do for the block filling, as in the standard routine, it is meant to be openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!$OMP PARALLEL PRIVATE(n,m,nu,qmax,pmin,pmax,low_aq,up_aq,low_bq,up_bq,q,p,pr,sommaA,sommaB)

!Choosing the chunk size
chunk=NNZab/(200*2)
chunk_if: IF (chunk==0) THEN
	chunk=1
END IF chunk_if
!WRITE (*,*) 'Chunk length =',chunk

!$OMP DO SCHEDULE(guided,chunk)
flat_do: DO i=1,NNZab

	!----------------------------------------------------------------------------------------
	!Assigning all the needed indexes, which will be private, I assume, in the openMP version
	!----------------------------------------------------------------------------------------
	n=m_index(i,1)
	m=m_index(i,2)
	nu=m_index(i,3)
	qmax=m_index(i,4)
	pmin=m_index(i,5)
	pmax=m_index(i,6)
	low_aq=m_index(i,8)
	up_aq=m_index(i,9)
	low_bq=m_index(i,10)
	up_bq=m_index(i,11)

	!-----------
	!Calcolo Avt
	!-----------
	!Calcolo della sommatoria e di Avt
	sommaA=(0.0D0,0.0D0)
	sommaA_do:DO q=0,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
		sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_j(p)
	END DO sommaA_do

	v_sAij(i)=v_c0(i)*sommaA

	!-----------
	!Calcolo Bvt
	!----------- 
	!Calcolo della sommatoria e di Bvt
	sommaB=(0.0D0,0.0D0)
	sommaB_do:DO q=1,qmax
		p=n+nu-2*q
		pr=REAL(p,dbl)
		sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_j(p+1)
	END DO sommaB_do

	v_sBij(i)=v_c0(i)*sommaB

END DO flat_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL


!$OMP PARALLEL PRIVATE(n,m,m1,nu,iab)

!Choosing the chunk size
chunk=4*NNZab/(200*2)
chunk_if_TU: IF (chunk==0) THEN
	chunk=1
END IF chunk_if_TU


!$OMP DO SCHEDULE(guided,chunk)
flat_TU_do: DO i=1,4*NNZab,2

	!----------------------------------------------------------------------------------------
	!Assigning all the needed indexes, which will be private, I assume, in the openMP version
	!----------------------------------------------------------------------------------------
	n=m_index_TU(i,1)
	n=m_index_TU(i+1,1)
	m=m_index_TU(i,2)
	m=m_index_TU(i+1,2)
	m1=m_index_TU(i,3)
	m1=m_index_TU(i+1,3)
	nu=m_index_TU(i,4)
	nu=m_index_TU(i+1,4)
	iab=m_index_TU(i,5)
	iab=m_index_TU(i+1,5)


	!assigning values to TU matrix elements
	m1_if: IF (m1==1) THEN
		v_sTUij(i)=v_sAij(iab)*m_u2(n,nu)
		v_sTUij(i+1)=v_sBij(iab)*m_t2(n,nu)
	ELSE
		v_sTUij(i)=v_sBij(iab)*m_u1(n,nu)
		v_sTUij(i+1)=v_sAij(iab)*m_t1(n,nu)
	END IF m1_if

END DO flat_TU_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE fill_TU_sparse_pardiso_fast


!******************************************************************************
!3) SUBROUTINE fill_jBlock_iBlock_sparse_ss_dip_rhs:riempio i vettori indici riga e colonna
! per la disposizione dei blocchi di matrice per il rhs della dip single shell
!******************************************************************************
SUBROUTINE fill_jBlock_iBlock_sparse_ss_dip_rhs(ndip,v_jBlock_rhs,v_iBlock_rhs)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ndip						! N dipoli
INTEGER(lo), DIMENSION(:),INTENT(OUT) :: v_jBlock_rhs,v_iBlock_rhs			! Numero di blocchi diversi da zero


! Dichiarazione variabili interne
INTEGER(lo) :: i,j,next

! Funzione vera e propria
next=1
!Per il primo blocco e' sempre uno
v_iBlock_rhs(1)=1

!Adesso il ciclo per riempire i vettori del pattern della matrice a blocchi
row_do: DO i=1,ndip
	col_do: DO j=1,ndip

		IF (i/=j) CYCLE col_do

		!Assegno i valori di colonna
		v_jBlock_rhs(next)=j
		next=next+1

	END DO col_do

	v_iBlock_rhs(i+1)=next

END DO row_do

END SUBROUTINE fill_jBlock_iBlock_sparse_ss_dip_rhs



!**********************************************************************************************************
!4) SUBROUTINE fill_D_PHI_sparse:riempio tutta la struttura per i blocchi di 
! rotazione Dij e per Exp[phi_ij]
!**********************************************************************************************************
SUBROUTINE fill_D_PHI_sparse_ss_dip_rhs(ndip,nstop,nnb_rhs,m_xyz,v_jBlock_rhs,v_iBlock_rhs,dip_flag,m_Dnkm,m_jDnkm,m_iDnkm,m_exphi)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ndip,nstop,nnb_rhs			        ! N dipoli, N multipolar exp
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock_rhs,v_iBlock_rhs	! Numero di blocchi diversi da zero
INTEGER(lo), INTENT(IN) :: dip_flag                                   ! dipole position flag

REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Dnkm			! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jDnkm,m_iDnkm		! Matrici indici Dij
COMPLEX(dbl), DIMENSION(-nstop:nstop,1:nnb_rhs), INTENT(OUT) :: m_exphi	!Matrice Esponenziali

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,xc,yc,zc,rij,thetaij,phiij


!Funzione vera e propria

!Assigning either core or shell coordinates depending on dip_flag
pos_if: IF (dip_flag==0) THEN
     xc=m_xyz(1,1)
     yc=m_xyz(1,2)
     zc=m_xyz(1,3)
ELSE
     xc=m_xyz(2,1)
     yc=m_xyz(2,2)
     zc=m_xyz(2,3)
END IF pos_if

row_do: DO i=1,ndip

     lowb=v_iBlock_rhs(i)
     upb=v_iBlock_rhs(i+1)-1

     col_do: DO next=lowb,upb

          !Adesso so sia riga che colonna
          j=v_jBlock_rhs(next)

          !Calcolo le coordinate sferiche relative
          xij=xc-m_xyz(i+2,1)
          yij=yc-m_xyz(i+2,2)
          zij=zc-m_xyz(i+2,3)

          CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

          !Chiamo la subroutine per riempire Dij
          CALL fillblock_Dkmn(nstop,thetaij,m_Dnkm(:,next),m_jDnkm(:,next),m_iDnkm(:,next),error)

          !Riempio anche la colonna degli esponenziali
          phi_do: DO m=-nstop,nstop
               m_exphi(m,next)=EXP( (0.0D0,1.0D0)*REAL(m,dbl)*phiij )
          END DO phi_do

     END DO col_do
END DO row_do

END SUBROUTINE fill_D_PHI_sparse_ss_dip_rhs



!******************************************************************************
!5) SUBROUTINE fill_AB_sparse:riempio tutta la struttura per i blocchi di 
! rotazione Aij e per Bij
!******************************************************************************
SUBROUTINE fill_AB_sparse_ss_dip_rhs(ndip,nstop,k,m_xyz,v_c0,v_aq,v_bq,v_jBlock_rhs,v_iBlock_rhs,dip_flag,&
                                   & m_Aij,m_Bij,m_Aij_sca,m_Bij_sca,m_jABij,m_iABij)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: ndip,nstop								! N dipoles, N multipolar exp e N nnz blocks
REAL(dbl) ,INTENT(IN) :: k									! vettore d'onda 
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz							! Matrix posizione
INTEGER(lo), DIMENSION(:),INTENT(IN) :: v_jBlock_rhs,v_iBlock_rhs				! Numero di blocchi diversi da zero
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq						! Vettore norm,gaunt e bq
INTEGER(lo), INTENT(IN) :: dip_flag                                                           ! dipole position flag

COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_Aij,m_Bij,m_Aij_sca,m_Bij_sca			! Matrice blocchi Dij
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_jABij,m_iABij	             			! Matrici indici Dij

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,lowb,upb,next,error,m
REAL(dbl) :: xij,yij,zij,xc,yc,zc,rij,thetaij,phiij


!Funzione vera e propria

!Assigning either core or shell coordinates depending on dip_flag
pos_if: IF (dip_flag==0) THEN
     xc=m_xyz(1,1)
     yc=m_xyz(1,2)
     zc=m_xyz(1,3)
ELSE
     xc=m_xyz(2,1)
     yc=m_xyz(2,2)
     zc=m_xyz(2,3)
END IF pos_if

row_do: DO i=1,ndip

    lowb=v_iBlock_rhs(i)
    upb=v_iBlock_rhs(i+1)-1

    col_do: DO next=lowb,upb

        !Adesso so sia riga che colonna
        j=v_jBlock_rhs(next)

        !Calcolo le coordinate sferiche relative
        xij=xc-m_xyz(i+2,1)
        yij=yc-m_xyz(i+2,2)
        zij=zc-m_xyz(i+2,3)

        CALL cart_spher_r_dbl1(xij,yij,zij,rij,thetaij,phiij,error)

        !Chiamo la subroutine per riempire Dij
        ! CALL fillblock_sparse_onlysca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next), &
        !                     & m_jABij(:,next),m_iABij(:,next),error)

        CALL fillblock_sparse_sca(nstop,v_c0,v_aq,v_bq,k*rij,m_Aij(:,next),m_Bij(:,next),&
             & m_Aij_sca(:,next),m_Bij_sca(:,next),m_jABij(:,next),m_iABij(:,next),error)

    END DO col_do

END DO row_do

END SUBROUTINE fill_AB_sparse_ss_dip_rhs

!******************************************************************************************
!6) SUBROUTINE fill_TU_sparse_pardiso_dip: riempio un blocco delle matrici Aij e Bij, TUij ma
!seguendo la data structure nuova o meglio, ora che posso, scrivo TU come matrice unica in
!fashion compressed sparse row, e poi scrivo Aij e Bij, che mi servono comunque, come al
!solito
!*****************************************************************************************
SUBROUTINE fill_TU_sparse_pardiso_dip(nstop,v_c0,v_aq,v_bq,v_qa,v_qb,m_t1,m_u1,m_t2,m_u2,kr,tflag,&
				&v_sAij,v_sBij,v_sAijh,v_sBijh,v_jABij,v_iABij,v_sTUij,v_jTUij,v_iTUij,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! N espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: kr							! Distanza radiale kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq				! Vettore norm,gaunt e bq
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_qa,v_qb				! Vettore coeff ausiliari espansioni
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_t1,m_u1,m_t2,m_u2		! Matrici ausiliarie per sistema lineare
INTEGER(lo), INTENT(IN) :: tflag							! Numero sfere eq e lato matrice

INTEGER(lo), INTENT(OUT) :: error							! Flag di errore
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sAij,v_sBij			! Vettori blocco sparse ij delle matrici A e B
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sAijh,v_sBijh			!Vettori blocco sparse ij delle matrici A e B con h_n
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jABij,v_iABij			! Vettori sparse colonne e righe
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sTUij    			! Vettori blocco sparse ij delle matrici A e B
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jTUij,v_iTUij			! Vettori sparse colonne e righe

! Dichiarazione variabili interne
INTEGER(lo) :: m,n,nu,i,i1,i2,j,q,p,low_nu,qmax,m1,m2			!Indici per il ciclo
INTEGER(lo) :: nnab									!Non zero elements ABij
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq					!Bounds per i vettori aq,bq
REAL(dbl) :: mr,nr,mur,nur,pr								!Indici reali 
INTEGER(lo) :: pmin,pmax,nbes							!Estremi per i vett e n bes
COMPLEX(dbl) :: sommaA,sommaB,sommaAh,sommaBh					!Somme parziali per vec trans
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_j,v_y,v_h			!Vett comp j_n,y_n,h_n per argomento complesso
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_Apmn					!Vettori fattori Apmn

! Funzione vera e propria

!Inizializzo indici e bounds

!Indici
i=0
v_iABij(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Calcolo l'upper bound per i miei vettori AB ausiliari
nnab=SIZE(v_sAij,1)

!Adesso posso allocare tutti ma proprio tutti i miei vettori
ALLOCATE(v_j(0:(2*nstop+1)),v_y(0:(2*nstop+1)),v_h(0:(2*nstop+1)))

!Calcolo infine la mie funzioni di bessel e di hankel
CALL besselj_z_sub((2*nstop+1),kr,v_j,error)
CALL bessely_z_sub((2*nstop+1),kr,v_y,error)
v_h=v_j+(0.0D0,1.0D0)*v_y

!Check delle funzioni
!WRITE(*,*) "kr: ", kr
!WRITE(*,*)

!WRITE(*,*) "v_j"
!DO i=1,2*nstop+1
!WRITE(*,*) v_j(i)
!END DO
!WRITE(*,*)

!WRITE(*,*) "v_y"
!DO i=1,2*nstop+1
!WRITE(*,*) v_y(i)
!END DO
!WRITE(*,*)

!WRITE(*,*) "v_h"
!DO i=1,2*nstop+1
!WRITE(*,*) v_h(i)
!!40 FORMAT (2ES21.12)
!END DO
!WRITE(*,*)


!------------------------------------------------------------------------
!Questo e' il ciclo per riempire le matrici A e B che mi servono comunque
!------------------------------------------------------------------------
j=0
AB_n_do: DO n=1,nstop
	AB_m_do: DO m=-n,n

		!Aggiorno l'indice per il vettore di preconditioning
		j=j+1

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		AB_nu_do: DO nu=low_nu,nstop


			!Rendo reali gli indici
			mr=REAL(m,dbl)
			nr=REAL(n,dbl)
			mur=REAL(m,dbl)
			nur=REAL(nu,dbl)

			!Calcolo qmax
			qmax=MIN(n,nu)

			!Aggiorno i
			i=i+1

			!Aggiorno i bounds
			ia_if: IF (i==1) THEN
				low_aq=1
			ELSE
				low_aq=up_aq+1  
			END IF ia_if

			ib_if: IF (i==1) THEN
				low_bq=1
			ELSE
				low_bq=up_bq+1  
			END IF ib_if

			up_aq=low_aq+qmax

			up_bq=low_bq+qmax-1

			!Calcolo Pmin e Pmax
			pmin=n+nu -2*qmax
			pmax=n+nu
			nbes=pmax-pmin+2

			!Adesso posso allocare tutti ma proprio tutti i miei vettori
			ALLOCATE(v_Apmn(pmin:pmax+1))

			!Inizializzo v_Apmn
			v_Apmn=0.0D0


			!-----------
			!Calcolo Avt
			!-----------
			!Calcolo del fattore numerico sotto sommatoria per Avt
			DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				v_Apmn(p)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
			END DO


			!Calcolo della sommatoria e di Avt
			sommaA=(0.0D0,0.0D0)
			sommaAh=(0.0D0,0.0D0)
			sommaA_do:DO q=0,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				sommaA=sommaA+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_j(p)
				sommaAh=sommaAh+((0.0D0,1.0D0)**p)*v_Apmn(p)*v_aq(low_aq+q)*v_h(p)
			END DO sommaA_do

			v_sAij(i)=v_c0(i)*sommaA
			v_sAijh(i)=v_c0(i)*sommaAh

			!-----------
			!Calcolo Bvt
			!-----------
			!Calcolo della sommatoria e di Bvt
			sommaB=(0.0D0,0.0D0)
			sommaBh=(0.0D0,0.0D0)
			sommaB_do:DO q=1,qmax
				pr=nr+nur-2.0D0*REAL(q,dbl)
				p=INT(pr,lo)
				sommaB=sommaB+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_j(p+1)
				sommaBh=sommaBh+((0.0D0,1.0D0)**(p+1))*v_bq(low_bq+(q-1))*v_h(p+1)
			END DO sommaB_do

			v_sBij(i)=v_c0(i)*sommaB
			v_sBijh(i)=v_c0(i)*sommaBh

			DEALLOCATE(v_Apmn)

			!-----------------------------------------
			!Aggiorno gli indici per lo storage sparse
			!-----------------------------------------
			v_jABij(i)=nu*(nu+1)+m

		END DO AB_nu_do

		!-----------------------------------
		!Dove comincia la prossima riga?Qui!
		!-----------------------------------
		v_iABij(n*(n+1)+m+1)=i+1

	END DO AB_m_do
END DO AB_n_do





!------------------------------------------------------------------------
!Questo e' il ciclo per riempire la matrici TU, quindi ho anche i cicli
!su m1 e m2
!------------------------------------------------------------------------
!Indici
i=0			!Indice per TU
i1=0			!indice per A e B
i2=0			!indice per A e B
v_iTUij(1)=1

n_do: DO n=1,nstop
	m_do: DO m=-n,n

		!Calcolo il lower bounds per il ciclo successivo
		IF (m==0) THEN
			low_nu=1
		ELSE
			low_nu=ABS(m)
		END IF

		m1_do: DO m1=1,2

			nu_do: DO nu=low_nu,nstop


				!Aggiorno gli indici
				i=i+2						!Salto di due in due per quello di TU

				tflag_if: IF (tflag==0) THEN

					!Assegno i valori per la matrice TU, con un if su m1
					m1_if: IF (m1==1) THEN
						i1=i1+1
						v_sTUij(i-1)=v_sAij(i1)*m_u2(n,nu)
						v_sTUij(i)=v_sBij(i1)*m_t2(n,nu)
					ELSE
						i2=i2+1
						v_sTUij(i-1)=v_sBij(i2)*m_u1(n,nu)
						v_sTUij(i)=v_sAij(i2)*m_t1(n,nu)
					END IF m1_if

				ELSE 

					!Assegno i valori per la matrice TU, con un if su m1
					m1_if1: IF (m1==1) THEN
						i1=i1+1
!						v_sTUij(i-1)=(v_sAij(i1)+v_qa(nu)*v_sAijh(i1))*m_u2(n,nu)
!						v_sTUij(i)=(v_sBij(i1)+v_qb(nu)*v_sBijh(i1))*m_t2(n,nu)
						v_sTUij(i-1)=(v_sAij(i1)-v_qa(nu)*v_sAijh(i1))*m_u2(n,n) !Good ones
						v_sTUij(i)=(v_sBij(i1)-v_qb(nu)*v_sBijh(i1))*m_t2(n,n)
					ELSE
						i2=i2+1
!						v_sTUij(i-1)=(v_sBij(i1)+v_qb(nu)*v_sBijh(i1))*m_u1(n,nu)
!						v_sTUij(i)=(v_sAij(i1)+v_qa(nu)*v_sAijh(i1))*m_t1(n,nu)
						v_sTUij(i-1)=(v_sBij(i2)-v_qa(nu)*v_sBijh(i2))*m_u1(n,n) !Good ones
						v_sTUij(i)=(v_sAij(i2)-v_qb(nu)*v_sAijh(i2))*m_t1(n,n)
					END IF m1_if1

				END IF tflag_if

				!-----------------------------------------
				!Aggiorno gli indici per lo storage sparse
				!-----------------------------------------
				v_jTUij(i-1)=2*(nu*(nu+1)+m-1)+1
				v_jTUij(i)=2*(nu*(nu+1)+m-1)+2

			END DO nu_do

			!-----------------------------------
			!Dove comincia la prossima riga?Qui!
			!-----------------------------------
			v_iTUij(2*(n*(n+1)+m-1)+m1+1)=i+1

		END DO m1_do
	END DO m_do
END DO n_do


END SUBROUTINE fill_TU_sparse_pardiso_dip



!******************************************************************************
!7) SUBROUTINE cext_random_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE csca_shell_sub(nstop,ns,k,v_ab_host,v_csca)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k							!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab_host				!Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_csca				!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j
REAL(dbl) :: somma

!Comincia la subroutine vera e propria
j=0
		
sphere_do: DO i=1,ns

	!Inizializzo la variabile somma
	somma=(0.0D0,0.0D0)

	n_do: DO n=1,nstop

		m_do: DO m=-n,n

			!Aggiorno l'indice
			j=j+1

			!Costruisco la somma
			somma=somma+(ABS(v_ab_host(2*j-1))**2)+(ABS(v_ab_host(2*j))**2)

		END DO m_do

	END DO n_do

	!Calcolo la funzione
	v_csca(i)=(4*Pi_d/(k**2))*somma

END DO sphere_do

END SUBROUTINE csca_shell_sub

!******************************************************************************
!8) SUBROUTINE cabs_shell_sub: sezione assorbimento per sistema core shell
!******************************************************************************
SUBROUTINE cabs_shell_sub(lambda,k,ref_index,v_req,m_epseq,nstop,neq,blockside,matrixside,&
     & v_ab_shell1,v_dc_shell1,v_dc_core,v_cabs,error,v_cabs_n) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index,k					! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req					! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq					! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop							! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside,blockside				! Numero sfere eq e lato matrice
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab_shell1,v_dc_shell1	! Espansioni shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc_core				! espansioni core

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs				! Matrici coeff. singola sfera
REAL(dbl), DIMENSION(:), INTENT(OUT), OPTIONAL :: v_cabs_n              ! Absorption multipole contribution
INTEGER(lo), INTENT(OUT) :: error							! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x										! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc									! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq) :: v_k									! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq) :: v_m									! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx								! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z22						! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z11,v_csi_z22						! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11,v_psi1_z22						! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z11,v_csi1_z22						! Vettori Csi' complessi
COMPLEX(dbl), DIMENSION(1:blockside) :: v_f1,v_f2,v_f3,v_f4,v_f5,v_f6,v_f7,v_f8		! Vettori prodotto di passaggio
INTEGER(lo) :: m,n,nu,i,j											! Indici
REAL(dbl) :: nr													! Indice complessificato
COMPLEX(dbl) :: somma_core, somma_coreshell								! Somme parziali cabs
COMPLEX(dbl) :: n_sum_c=0.0D0                                                                !Partial Sum for the n-th multipole


! Inizio della procedura vera e propria

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
     !Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
     !da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
     v_epsc(j)=CMPLX(m_epseq(neq+1-j,1),m_epseq(neq+1-j,2))
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Calcolo il vettore per i wavevectors
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq)=(2*pi_d*SQRT(v_epsc))/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq)=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)

!Calcolo la matrice mx
mx_out_do: DO i=0,neq
     mx_in_do: DO j=1,neq

          m_mx(i,j)=v_m(i)*CMPLX(v_x(j),KIND=dbl)

     END DO mx_in_do
END DO mx_out_do

!Loop per calcolare le matrici coefficiente a e b 

!Calcolo riccati bessel di prima specie Complessa
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,10)
10   FORMAT ("Errore: errore in psi_d_sub chiamata da cabs_shell_sub")
     RETURN
END IF psi_z_if

!Calcolo riccati bessel di terza specie complesso
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(2,2),v_csi_z22,error)

csi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,20)
20   FORMAT ("Errore: errore in csi_z_sub chiamata da cabs_shell_sub")
     RETURN
END IF csi_z_if

!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
der_do: DO n=1,nstop

     nr=REAL(n,dbl)
     !Derivate di Psi
     v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
     v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

     !Derivate di csi
     v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
     v_csi1_z22(n)=v_csi_z22(n-1)-nr*v_csi_z22(n)/m_mx(2,2)

END DO der_do

!output di controllo per i vettori di riccatibessel
!WRITE(*,*) "Argomento funzione 11: ", m_mx(1,1)
!WRITE(*,*) "Psi11: ", (v_psi_z11(n), n=1,nstop)
!WRITE(*,*) "Psi-11: ", (v_psi1_z11(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi11: ", (v_csi_z11(n), n=1,nstop)
!WRITE(*,*) "Csi1-11: ", (v_csi1_z11(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Argomento funzione 22: ", m_mx(2,2)
!WRITE(*,*) "Psi22: ", (v_psi_z22(n), n=1,nstop)
!WRITE(*,*) "Psi1-22: ", (v_psi1_z22(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi12: ", (v_csi_z12(n), n=1,nstop)
!WRITE(*,*) "Csi1-12: ", (v_csi1_z12(n), n=1,nstop)
!WRITE(*,*)

!Calcolo i valori per i vettori iniziali,che sono tutto quello che mi serve per una shell concentrica
j=0
factor_out_do: DO n=1,nstop
     factor_in_do: DO m=-n,n

          !Aggiorno l'indice
          j=j+1

          !Riempio i vettori di passaggio per calcolare la sezione di assorbimento core shell
          v_f1(j)=(ABS(v_ab_shell1(2*j))**2)*v_csi_z11(n)*CONJG(v_csi1_z11(n))
          v_f2(j)=(ABS(v_ab_shell1(2*j-1))**2)*v_csi1_z11(n)*CONJG(v_csi_z11(n))
          v_f3(j)=(ABS(v_dc_shell1(2*j))**2)*v_psi_z11(n)*CONJG(v_psi1_z11(n))
          v_f4(j)=(ABS(v_dc_shell1(2*j-1))**2)*v_psi1_z11(n)*CONJG(v_psi_z11(n))


          v_f5(j)=v_ab_shell1(2*j-1)*CONJG(v_dc_shell1(2*j-1))*v_csi1_z11(n)*CONJG(v_psi_z11(n))
          v_f6(j)=v_dc_shell1(2*j-1)*CONJG(v_ab_shell1(2*j-1))*v_psi1_z11(n)*CONJG(v_csi_z11(n))
          v_f7(j)=v_ab_shell1(2*j)*CONJG(v_dc_shell1(2*j))*v_csi_z11(n)*CONJG(v_psi1_z11(n))
          v_f8(j)=v_dc_shell1(2*j)*CONJG(v_ab_shell1(2*j))*v_psi_z11(n)*CONJG(v_csi1_z11(n))

     END DO factor_in_do
END DO factor_out_do



!Inizializzo la variabile somma
somma_core=(0.0D0,0.0D0)
somma_coreshell=(0.0D0,0.0D0)
j=0
n_do: DO n=1,nstop

     n_sum_c=(0.0D0,0.0D0)

     m_do: DO m=-n,n

          !Aggiorno l'indice
          j=j+1

          !Costruisco la somma per il core
          somma_core=somma_core + v_psi_z22(n)*CONJG(v_psi1_z22(n))*( (ABS(v_dc_core(2*j)))**2 ) - &
               & v_psi1_z22(n)*CONJG(v_psi_z22(n))*( (ABS(v_dc_core(2*j-1)))**2 )

          !Costruisco la somma per il core shell
          somma_coreshell=somma_coreshell+v_f1(j)-v_f2(j)+v_f3(j)-v_f4(j)+v_f5(j)+v_f6(j)-v_f7(j)-v_f8(j)

          !Partial sum

          n_sum_c=n_sum_c+v_f1(j)-v_f2(j)+v_f3(j)-v_f4(j)+v_f5(j)+v_f6(j)-v_f7(j)-v_f8(j)

     END DO m_do

     multi_if: IF (PRESENT(v_cabs_n)) THEN
          v_cabs_n(n)=-(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*n_sum_c/v_m(1),dbl)
     END IF multi_if

END DO n_do

!Calcolo la funzione per il core
v_cabs(1)=-(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma_core/v_m(2),dbl)

!Calcolo la funzione per la sola shell, per sottrazione
v_cabs(2)=-(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma_coreshell/v_m(1),dbl)
!WRITE(*,*)
!WRITE(*,*) "k", k
!WRITE(*,*) "m", v_m(1)
!WRITE(*,*) "somma_coreshell",somma_coreshell
!WRITE(*,*) "presa parte reale", REAL((0.0D0,1.0D0)*somma_coreshell/v_m(1),dbl)
!WRITE(*,*) "Sezione", v_cabs(2)

END SUBROUTINE cabs_shell_sub



!******************************************************************************
!7bis) SUBROUTINE cext_random_sub: calcolo le sezioni di estinzione per ciascuna sfera
!******************************************************************************
SUBROUTINE rad_shell_sub_dip(nstop,ns,k,v_ab_host,v_csca)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns						!N espansioni multipolari e N sfere
REAL(dbl), INTENT(IN) :: k							!Vettore d'onda e polarizzazione
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab_host				!Vettori coefficienti espansioni campo, non corretti e corretti  


REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_csca				!Vettore sezioni di estinzione

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,j
REAL(dbl) :: somma

!Comincia la subroutine vera e propria
j=0

sphere_do: DO i=1,ns

     !Inizializzo la variabile somma
     somma=(0.0D0,0.0D0)

     n_do: DO n=1,nstop

          m_do: DO m=-n,n

               !Aggiorno l'indice
               j=j+1

               !Costruisco la somma
               somma=somma+(ABS(v_ab_host(2*j-1))**2)+(ABS(v_ab_host(2*j))**2)

          END DO m_do

     END DO n_do

     !Calcolo la funzione
     v_csca(i)=(4*Pi_d/(k**2))*somma

END DO sphere_do

END SUBROUTINE rad_shell_sub_dip




!******************************************************************************
!9) SUBROUTINE cabs_shell_sub: sezione assorbimento per sistema core shell
!******************************************************************************
SUBROUTINE cabs_shell_sub_dip(lambda,k,ref_index,v_req,m_epseq,nstop,neq,blockside,matrixside,&
     & v_ab_shell1,v_dc_shell1,v_ab_shell2,v_dc_shell2,v_cabs,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index,k					! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req					! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq					! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop							! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside,blockside				! Numero sfere eq e lato matrice
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab_shell1,v_dc_shell1		! Espansioni shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_ab_shell2,v_dc_shell2		! espansioni core

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs				! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error							! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq-1) :: v_x										! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq-1) :: v_epsc									! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_k									! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_m									! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq-1,1:neq-1) :: m_mx								! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z12						! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z11,v_csi_z12						! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11,v_psi1_z12						! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z11,v_csi1_z12						! Vettori Csi' complessi
COMPLEX(dbl), DIMENSION(1:blockside) :: v_f1,v_f2,v_f3,v_f4,v_f5,v_f6,v_f7,v_f8		! Vettori prodotto di passaggio
COMPLEX(dbl), DIMENSION(1:blockside) :: v_g1,v_g2,v_g3,v_g4,v_g5,v_g6,v_g7,v_g8		! Vettori prodotto di passaggio
INTEGER(lo) :: m,n,nu,i,j											! Indici
REAL(dbl) :: nr													! Indice complessificato
COMPLEX(dbl) :: somma_shellint, somma_shellext								! Somme parziali cabs



! Inizio della procedura vera e propria

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq-1
     !Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
     !da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
     v_epsc(j)=CMPLX(m_epseq((neq-1)+1-j,1),m_epseq((neq-1)+1-j,2))
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Calcolo il vettore per i wavevectors
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq-1)=(2*pi_d*SQRT(v_epsc(1:neq-1)))/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq-1)=SQRT(v_epsc(1:neq-1))/CMPLX(ref_index,KIND=dbl)

!Calcolo la matrice mx
mx_out_do: DO i=0,neq-1
     mx_in_do: DO j=1,neq-1

          m_mx(i,j)=v_m(i)*CMPLX(v_x(j),KIND=dbl)

     END DO mx_in_do
END DO mx_out_do

!Loop per calcolare le matrici coefficiente a e b 

!Calcolo riccati bessel di prima specie Complessa
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(1,2),v_psi_z12,error)

psi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,10)
10   FORMAT ("Errore: errore in psi_d_sub chiamata da cabs_shell_sub")
     RETURN
END IF psi_z_if

!Calcolo riccati bessel di terza specie complesso
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(1,2),v_csi_z12,error)

csi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,20)
20   FORMAT ("Errore: errore in csi_z_sub chiamata da cabs_shell_sub")
     RETURN
END IF csi_z_if

!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
der_do: DO n=1,nstop

     nr=REAL(n,dbl)
     !Derivate di Psi
     v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
     v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)

     !Derivate di csi
     v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
     v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)

END DO der_do


!Calcolo i valori per i vettori iniziali,che sono tutto quello che mi serve per una shell concentrica
j=0
factor_out_do: DO n=1,nstop
     factor_in_do: DO m=-n,n

          !Aggiorno l'indice
          j=j+1

          !Riempio i vettori di passaggio per calcolare la sezione di assorbimento core shell, confine esterno
          v_f1(j)=(ABS(v_ab_shell1(2*j))**2)*v_csi_z11(n)*CONJG(v_csi1_z11(n))
          v_f2(j)=(ABS(v_ab_shell1(2*j-1))**2)*v_csi1_z11(n)*CONJG(v_csi_z11(n))
          v_f3(j)=(ABS(v_dc_shell1(2*j))**2)*v_psi_z11(n)*CONJG(v_psi1_z11(n))
          v_f4(j)=(ABS(v_dc_shell1(2*j-1))**2)*v_psi1_z11(n)*CONJG(v_psi_z11(n))


          v_f5(j)=v_ab_shell1(2*j-1)*CONJG(v_dc_shell1(2*j-1))*v_csi1_z11(n)*CONJG(v_psi_z11(n))
          v_f6(j)=v_dc_shell1(2*j-1)*CONJG(v_ab_shell1(2*j-1))*v_psi1_z11(n)*CONJG(v_csi_z11(n))
          v_f7(j)=v_ab_shell1(2*j)*CONJG(v_dc_shell1(2*j))*v_csi_z11(n)*CONJG(v_psi1_z11(n))
          v_f8(j)=v_dc_shell1(2*j)*CONJG(v_ab_shell1(2*j))*v_psi_z11(n)*CONJG(v_csi1_z11(n))


          !Riempio i vettori di passaggio per calcolare la sezione di assorbimento core shell, confine interno
          v_g1(j)=(ABS(v_ab_shell2(2*j))**2)*v_csi_z12(n)*CONJG(v_csi1_z12(n))
          v_g2(j)=(ABS(v_ab_shell2(2*j-1))**2)*v_csi1_z12(n)*CONJG(v_csi_z12(n))
          v_g3(j)=(ABS(v_dc_shell2(2*j))**2)*v_psi_z12(n)*CONJG(v_psi1_z12(n))
          v_g4(j)=(ABS(v_dc_shell2(2*j-1))**2)*v_psi1_z12(n)*CONJG(v_psi_z12(n))


          v_g5(j)=v_ab_shell2(2*j-1)*CONJG(v_dc_shell2(2*j-1))*v_csi1_z12(n)*CONJG(v_psi_z12(n))
          v_g6(j)=v_dc_shell2(2*j-1)*CONJG(v_ab_shell2(2*j-1))*v_psi1_z12(n)*CONJG(v_csi_z12(n))
          v_g7(j)=v_ab_shell2(2*j)*CONJG(v_dc_shell2(2*j))*v_csi_z12(n)*CONJG(v_psi1_z12(n))
          v_g8(j)=v_dc_shell2(2*j)*CONJG(v_ab_shell2(2*j))*v_psi_z12(n)*CONJG(v_csi1_z12(n))

     END DO factor_in_do
END DO factor_out_do



!Inizializzo la variabile somma
somma_shellint=(0.0D0,0.0D0)
somma_shellext=(0.0D0,0.0D0)
j=0
n_do: DO n=1,nstop

     m_do: DO m=-n,n

          !Aggiorno l'indice
          j=j+1

          !Costruisco la somma per il core shell interno
          somma_shellint=somma_shellint+v_g1(j)-v_g2(j)+v_g3(j)-v_g4(j)+v_g5(j)+v_g6(j)-v_g7(j)-v_g8(j)

          !Costruisco la somma per il core shell esterno
          somma_shellext=somma_shellext+v_f1(j)-v_f2(j)+v_f3(j)-v_f4(j)+v_f5(j)+v_f6(j)-v_f7(j)-v_f8(j)

     END DO m_do

END DO n_do

!Calcolo la funzione per il core. Di solito davanti c'e' il meno, ma qui lo ho tolto perche' cambia l'orientamento della superficie
!attraverso la quale calcolo il flusso del vettore di poynting
v_cabs(1)=(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma_shellint/v_m(1),dbl)

!Calcolo la funzione per la sola shell, per sottrazione
v_cabs(2)=(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma_shellext/v_m(1),dbl)

END SUBROUTINE cabs_shell_sub_dip


!******************************************************************************
!8) SUBROUTINE cabs_shell_sub: sezione assorbimento per sistema core shell
!******************************************************************************
SUBROUTINE cabs_shell_sub_coredip(lambda,k,ref_index,v_req,m_epseq,nstop,neq,blockside,matrixside,&
     & v_p,v_dc_core,v_cabs,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index,k					! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req					! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq					! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop							! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside,blockside				! Numero sfere eq e lato matrice
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p	! Espansioni shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc_core				! espansioni core

REAL(dbl), DIMENSION(:), INTENT(OUT) :: v_cabs				! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error							! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x										! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc									! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq) :: v_k									! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq) :: v_m									! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx								! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z22								! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z22								! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z22								! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z22								! Vettori Csi' complessi
COMPLEX(dbl), DIMENSION(1:blockside) :: v_f1,v_f2,v_f3,v_f4,v_f5,v_f6,v_f7,v_f8		! Vettori prodotto di passaggio
INTEGER(lo) :: m,n,nu,i,j											! Indici
REAL(dbl) :: nr													! Indice complessificato
COMPLEX(dbl) :: somma_core, somma_coreshell								! Somme parziali cabs



! Inizio della procedura vera e propria

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq-1
     !Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
     !da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
     v_epsc(j)=CMPLX(m_epseq((neq-1)+1-j,1),m_epseq((neq-1)+1-j,2))
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Calcolo il vettore per i wavevectors
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq-1)=(2*pi_d*SQRT(v_epsc(1:neq-1)))/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq-1)=SQRT(v_epsc(1:neq-1))/CMPLX(ref_index,KIND=dbl)

!Calcolo la matrice mx
mx_out_do: DO i=0,neq-1
     mx_in_do: DO j=1,neq-1

          m_mx(i,j)=v_m(i)*CMPLX(v_x(j),KIND=dbl)

     END DO mx_in_do
END DO mx_out_do

!Loop per calcolare le matrici coefficiente a e b 

!Calcolo riccati bessel di prima specie Complessa
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,10)
10   FORMAT ("Errore: errore in psi_d_sub chiamata da cabs_shell_sub")
     RETURN
END IF psi_z_if

!Calcolo riccati bessel di terza specie complesso
CALL csi_z_sub(nstop,m_mx(2,2),v_csi_z22,error)

csi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,20)
20   FORMAT ("Errore: errore in csi_z_sub chiamata da cabs_shell_sub")
     RETURN
END IF csi_z_if

!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
der_do: DO n=1,nstop

     nr=REAL(n,dbl)
     !Derivate di Psi
     v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

     !Derivate di csi
     v_csi1_z22(n)=v_csi_z22(n-1)-nr*v_csi_z22(n)/m_mx(2,2)

END DO der_do

!output di controllo per i vettori di riccatibessel
!WRITE(*,*) "Argomento funzione 11: ", m_mx(1,1)
!WRITE(*,*) "Psi11: ", (v_psi_z11(n), n=1,nstop)
!WRITE(*,*) "Psi-11: ", (v_psi1_z11(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi11: ", (v_csi_z11(n), n=1,nstop)
!WRITE(*,*) "Csi1-11: ", (v_csi1_z11(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Argomento funzione 22: ", m_mx(2,2)
!WRITE(*,*) "Psi22: ", (v_psi_z22(n), n=1,nstop)
!WRITE(*,*) "Psi1-22: ", (v_psi1_z22(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi12: ", (v_csi_z12(n), n=1,nstop)
!WRITE(*,*) "Csi1-12: ", (v_csi1_z12(n), n=1,nstop)
!WRITE(*,*)

!Calcolo i valori per i vettori iniziali,che sono tutto quello che mi serve per una shell concentrica
j=0
factor_out_do: DO n=1,nstop
     factor_in_do: DO m=-n,n

          !Aggiorno l'indice
          j=j+1

          !Riempio i vettori di passaggio per calcolare la sezione di assorbimento core shell
          v_f1(j)=(ABS(v_p(2*j))**2)*v_csi_z22(n)*CONJG(v_csi1_z22(n))
          v_f2(j)=(ABS(v_p(2*j-1))**2)*v_csi1_z22(n)*CONJG(v_csi_z22(n))
          v_f3(j)=(ABS(v_dc_core(2*j))**2)*v_psi_z22(n)*CONJG(v_psi1_z22(n))
          v_f4(j)=(ABS(v_dc_core(2*j-1))**2)*v_psi1_z22(n)*CONJG(v_psi_z22(n))


          v_f5(j)=v_p(2*j-1)*CONJG(v_dc_core(2*j-1))*v_csi1_z22(n)*CONJG(v_psi_z22(n))
          v_f6(j)=v_dc_core(2*j-1)*CONJG(v_p(2*j-1))*v_psi1_z22(n)*CONJG(v_csi_z22(n))
          v_f7(j)=v_p(2*j)*CONJG(v_dc_core(2*j))*v_csi_z22(n)*CONJG(v_psi1_z22(n))
          v_f8(j)=v_dc_core(2*j)*CONJG(v_p(2*j))*v_psi_z22(n)*CONJG(v_csi1_z22(n))

     END DO factor_in_do
END DO factor_out_do



!Inizializzo la variabile somma
somma_core=(0.0D0,0.0D0)
j=0
n_do: DO n=1,nstop

     m_do: DO m=-n,n

          !Aggiorno l'indice
          j=j+1

          !Costruisco la somma per il core
          somma_core=somma_core+v_f1(j)-v_f2(j)+v_f3(j)-v_f4(j)+v_f5(j)+v_f6(j)-v_f7(j)-v_f8(j)

     END DO m_do

END DO n_do

!Calcolo la funzione per il core
v_cabs(1)=(4*Pi_d/(k**2))*REAL((0.0D0,1.0D0)*somma_core/v_m(2),dbl)

!WRITE(*,*)
!WRITE(*,*) "k", k
!WRITE(*,*) "m", v_m(1)
!WRITE(*,*) "somma_coreshell",somma_coreshell
!WRITE(*,*) "presa parte reale", REAL((0.0D0,1.0D0)*somma_coreshell/v_m(1),dbl)
!WRITE(*,*) "Sezione", v_cabs(2)

END SUBROUTINE cabs_shell_sub_coredip


!******************************************************************************
!5tris-bis) SUBROUTINE field_expRandom_dip_ss_sub: calcolo i coefficienti del campo incidente,ma li randomizzo per la 
! fase e per l'orientamento,e questo mi serve per un tipo di calcoli. qui in piu' va bene per gli spettri, cosi' non mi
! dipende da lambda e lo uso per la single shell
!******************************************************************************
SUBROUTINE field_expRandom_dip_ss_sub(nstop,ns,v_alpha,v_beta,v_gamma,v_kd,v_p,error)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN):: nstop,ns				!N espansioni multipolari e N sfere
REAL(dbl), DIMENSION(:) :: v_alpha,v_beta,v_gamma,v_kd	!Angoli di eulero 
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_p			!Vettori coefficienti espansioni campo, non corretti e corretti  
INTEGER(lo) , INTENT(OUT) :: error				!Flag di errore

! Dichiarazione variabili interne
INTEGER(lo) :: i,j,jj,n,m,dimens					!Indici e dimensioni
REAL(dbl) :: nr,mr,x,y,z,d,theta,phi				!Indici reali e coordinate
REAL(dbl) :: kd								!Phase 
COMPLEX(dbl) :: pshift,expg1,expg2					!Phase shift
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_dnm1,v_dnm2		!Elementi ridotti di matrice di rotazione
REAL(dbl) :: alpha,beta,gamma		!Angoli di eulero 

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Alloco i vettori delle funzioni angolari
ALLOCATE(v_dnm1(1:dimens),v_dnm2(1:dimens))


!Inizializzo j e i vettori output
j=0
jj=0
v_p=(0.0D0,0.0D0)

i_do: DO i=1,ns

     !Calcolo gli angoli randomizzati
     alpha=v_alpha(i)*Pi_D

     beta=v_beta(i)*Pi_D

     gamma=v_gamma(i)*2*Pi_D

     !Funzione angolare Pi_mn
     CALL d_nm12(nstop,beta,v_dnm1,v_dnm2,error)

     dnm_if: IF (error/=0) THEN																	
          WRITE(*,10) 
10        FORMAT ("Errore in d_nm12 chiamata da field_expRandom_dip_ss_sub: il programma termina ora...") 
          STOP
     END IF dnm_if

     !Calcolo gli esponenziali
     expg1=EXP(gamma*(0.0D0,1.0D0))
     expg2=EXP(-gamma*(0.0D0,1.0D0))

     jj=0

     !Calcolo phase e phase shift
     kd=v_kd(i)
     pshift=EXP(CMPLX(0.0D0,kd))

     n_do: DO n=1,nstop

          !n reale
          nr=REAL(n,dbl)

          m_do: DO m=-n,n

               !m reale
               mr=REAL(m,dbl)	

               !Incremento j
               j=j+1
               jj=jj+1

               v_p(2*j-1)=0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
                    & ( v_dnm1(jj)*expg2-v_dnm2(jj)*expg1 )
               v_p(2*j)=  0.5D0*SQRT(2.0D0*nr+1.0D0)*pshift*EXP(-(0.0D0,1.0D0)*mr*alpha) * &
                    & ( v_dnm1(jj)*expg2+v_dnm2(jj)*expg1 )

          END DO m_do
     END DO n_do
END DO i_do

DEALLOCATE(v_dnm1,v_dnm2)


END SUBROUTINE field_expRandom_dip_ss_sub





!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!SUBROUTINES FOR THE SINGLE SHELL CODE FOLLOWING BORGHESE FORMALISM
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================



!******************************************************************************************
!1) SUBROUTINE fill_index_AB: filling al the relevant indexes for the sparse structure of 
!submatrixes composing the overall coefficient matrix. So matrix ad submatrix can be readily
!filled even in the oponMP implementation
!*****************************************************************************************
SUBROUTINE fill_index_AB(nstop,m_index,m_index_AABB,m_Apmn,v_jAB,v_iAB,v_jAABB,v_iAABB,error)



IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! Multipolar Expansions
INTEGER(lo), DIMENSION(:,:), INTENT(OUT) :: m_index,m_index_AABB			! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(OUT) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_jAB,v_iAB,v_jAABB,v_iAABB		! Row and columns CSR matrix Storage
INTEGER(lo), INTENT(OUT) :: error							! Error flag

! Dichiarazione variabili interne
INTEGER(lo) :: m,m1,n,nu,i,i1,i2,q,p,low_nu,qmax			!Loop indexes
REAL(dbl) :: mr,nr,mur,nur,pr						!Real loop indexes
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq				!aq bq vector Bounds
INTEGER(lo) :: pmin,pmax,nbes						!Bessel Hankel vector bounds

! Subroutine

!Inizializzo indici e bounds
error=0

!Index initialization
i=0
v_iAB(1)=1

!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1


!---------------------------------------------------------------------------------------
!Loop for the openMP implementation of the loops for the calculation of A and B
!---------------------------------------------------------------------------------------
AB_n_do: DO n=1,nstop
     AB_m_do: DO m=-n,n

          !Lower bound for the next loop
          IF (m==0) THEN
               low_nu=1
          ELSE
               low_nu=ABS(m)
          END IF

          AB_nu_do: DO nu=low_nu,nstop


               !Real indexes
               mr=REAL(m,dbl)
               nr=REAL(n,dbl)
               mur=REAL(m,dbl)
               nur=REAL(nu,dbl)

               !qmax
               qmax=MIN(n,nu)

               !Aggiorno i
               i=i+1

               !Updating the bounds
               ia_if: IF (i==1) THEN
                    low_aq=1
               ELSE
                    low_aq=up_aq+1  
               END IF ia_if

               ib_if: IF (i==1) THEN
                    low_bq=1
               ELSE
                    low_bq=up_bq+1  
               END IF ib_if

               up_aq=low_aq+qmax

               up_bq=low_bq+qmax-1

               !Pmin and Pmax
               pmin=n+nu -2*qmax
               pmax=n+nu
               nbes=pmax-pmin+2

               !-----------
               !Avt
               !-----------
               !Avt factor under the sum symbol for VTC computation
               DO q=0,qmax
                    pr=nr+nur-2.0D0*REAL(q,dbl)
                    p=INT(pr,lo)
                    m_Apmn(i,p+1)=nr*(nr+1.0D0)+nur*(nur+1.0D0)-pr*(pr+1.0D0)
               END DO


               !--------------------------------------------------------------
               !Filling the matrix containing all the indexes
               !--------------------------------------------------------------
               m_index(i,1)=n
               m_index(i,2)=m
               m_index(i,3)=nu
               m_index(i,4)=qmax
               m_index(i,5)=pmin
               m_index(i,6)=pmax
               m_index(i,7)=nbes
               m_index(i,8)=low_aq
               m_index(i,9)=up_aq
               m_index(i,10)=low_bq
               m_index(i,11)=up_bq

               !-----------------------------------------
               !Updating indexes for sparse storage
               !-----------------------------------------
               v_jAB(i)=nu*(nu+1)+m

          END DO AB_nu_do

          !-----------------------------------
          !Here begins the very next line
          !-----------------------------------
          v_iAB(n*(n+1)+m+1)=i+1

     END DO AB_m_do
END DO AB_n_do





!------------------------------------------------------------------------
! Loop for the indexes of the full AABB matrix, it is trickier because I have m1
!------------------------------------------------------------------------
!Indexes
i=0			!Running index for AABB
i1=0			!Running index for A e B
i2=0			!Running index for A e B
v_iAABB(1)=1

n_do: DO n=1,nstop
     m_do: DO m=-n,n

          !Lower bound for the next loop
          IF (m==0) THEN
               low_nu=1
          ELSE
               low_nu=ABS(m)
          END IF

          !Looping on every row of AABB. The columns are also doubled, but I handle that with i+2
          m1_do: DO m1=1,2

               nu_do: DO nu=low_nu,nstop


                    !Updating indexes
                    i=i+2						!Junping by 2 for AABB


                    !Filling AABB indexes, I know they are right, but it's i bit obscure
                    m1_if: IF (m1==1) THEN

                         i1=i1+1

                         !Index filling: 3 rows (n,m,m1) one column (nu) and one special progressive (i1)
                         m_index_AABB(i-1,1)=n
                         m_index_AABB(i,1)=n
                         m_index_AABB(i-1,2)=m
                         m_index_AABB(i,2)=m
                         m_index_AABB(i-1,3)=m1
                         m_index_AABB(i,3)=m1
                         m_index_AABB(i-1,4)=nu
                         m_index_AABB(i,4)=nu
                         !Special treatment here, since i1 and i2 alernate in filling the 5th space
                         m_index_AABB(i-1,5)=i1
                         m_index_AABB(i,5)=i1

                    ELSE

                         i2=i2+1

                         !Index filling: 3 rows (n,m,m1) one column (nu) and one special progressive (i2)
                         m_index_AABB(i-1,1)=n
                         m_index_AABB(i,1)=n
                         m_index_AABB(i-1,2)=m
                         m_index_AABB(i,2)=m
                         m_index_AABB(i-1,3)=m1
                         m_index_AABB(i,3)=m1
                         m_index_AABB(i-1,4)=nu
                         m_index_AABB(i,4)=nu
                         !Special treatment here, since i1 and i2 alernate in filling the 5th space
                         m_index_AABB(i-1,5)=i2
                         m_index_AABB(i,5)=i2

                    END IF m1_if

                    !-----------------------------------------
                    !Sparse storage indexes
                    !-----------------------------------------
                    v_jAABB(i-1)=2*(nu*(nu+1)+m-1)+1
                    v_jAABB(i)=2*(nu*(nu+1)+m-1)+2

               END DO nu_do

               !-----------------------------------
               !Starting index of the next row
               !-----------------------------------
               v_iAABB(2*(n*(n+1)+m-1)+m1+1)=i+1

          END DO m1_do
     END DO m_do
END DO n_do


END SUBROUTINE fill_index_AB



!******************************************************************************
!2) SUBROUTINE fill_M: Filling the coefficient matrix for the linear system
! to solve all the expansion coefficients 
!******************************************************************************
SUBROUTINE fill_M(nstop,kr,m_xyz,v_c0,v_aq,v_bq,v_Rn,v_Vn,v_Zn,m_index,m_index_AABB,m_Apmn,v_iAB,v_jAB,v_iAABB,v_jAABB,&
     & v_sA,v_sB,v_sAABB,v_iM,v_jM,v_sM,error)

IMPLICIT NONE

! Dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! Multipolar expansions
COMPLEX(dbl), INTENT(IN) :: kr								! Radial distance kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Norm,gaunt e bq vectors
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Core and host Coordinates
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_Rn,v_Vn,v_Zn				! M auxiliary vectors
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index,m_index_AABB			! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iAB,v_jAB,v_iAABB,v_jAABB			! CSR row column vectors for A,B,AB
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_iM,v_jM					! CSR row column vectors for M
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sA,v_sB,v_sAABB,v_sM			! CSR value matrixes
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Internal variables
INTEGER(lo) :: m,m1,n,nu,i,j,jjt,jjb,iab,q,p,low_nu,qmax,expb				!Loop indexes
REAL(dbl) :: mr,nr,mur,nur,pr,dz							!Loop indexes real type, and delta z coordinate
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq						!Bounds for aq,bq vectors
INTEGER(lo) :: pmin,pmax,nbes								!Bounds for Bessel vectors
INTEGER(lo) :: nnzAB,nnzAABB,nnzM							!Elements of (A,B), AB, M matrix
INTEGER(lo) :: blockside,AABBside,Mside						!Sides of (A,B), AB, M matrix
COMPLEX(dbl) :: sommaA,sommaB								!Partial Sums for VTC
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_j						!J Bessel
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_n,v_ABsign,v_index_i,v_index_jj		!n vector and and signs for A and B
INTEGER(lo) :: chunk=1								!Thread id, number of total threads and chunk size

! Subroutine

!Index and bounds initialization

!nnzAB: non-zero matrix elements of A or B
!blockside: side of A,B blocks
i=0
nnzAB=nstop*(1+2*nstop*(3+nstop))
nnzAB=nnzAB/3
blockside=nstop*(nstop+2)

!nnzABAB: non-zero matrix elements of AB
!AABBside: side of the off diagonal AB block
nnzAABB=4*nnzAB
AABBside=2*blockside

!nnzM: non-zero matrix elements of M
!Mside: side of the coefficient matrix
nnzM=2*nnzAABB+2*AABBside	!2 off diagonal blocks plus long diagonal
Mside=2*AABBside


!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of Bessel functions
ALLOCATE(v_j(0:(2*nstop+1)))
CALL besselj_z_sub((2*nstop+1),kr,v_j,error)

!Allocation of the matrix for the signs of A and B, and for the running index i
ALLOCATE(v_ABsign(1:nnzAABB),v_index_i(1:(AABBSide+nnzAABB)),v_index_jj(1:(AABBSide+nnzAABB)))


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Loop for A and B: openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!CALL OMP_SET_NUM_THREADS(1)

!$OMP PARALLEL PRIVATE(n,m,nu,qmax,pmin,pmax,low_aq,up_aq,low_bq,up_bq,q,p,pr,sommaA,sommaB)

!Choosing the chunk size
chunk=nnzAB/(200*2)
chunk_if_AB: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_AB
!WRITE (*,*) 'Chunk length =',chunk

!$OMP DO SCHEDULE(guided,chunk)
AB_do: DO i=1,nnzAB

     !----------------------------------------------------------------------------------------
     !Assigning all the needed indexes, which will be private, I assume, in the openMP version
     !----------------------------------------------------------------------------------------
     n=m_index(i,1)
     m=m_index(i,2)
     nu=m_index(i,3)
     qmax=m_index(i,4)
     pmin=m_index(i,5)
     pmax=m_index(i,6)
     low_aq=m_index(i,8)
     up_aq=m_index(i,9)
     low_bq=m_index(i,10)
     up_bq=m_index(i,11)

     !-----------
     !Calcolo Avt
     !-----------
     !Calcolo della sommatoria e di Avt
     sommaA=(0.0D0,0.0D0)
     sommaA_do:DO q=0,qmax
          p=n+nu-2*q
          pr=REAL(p,dbl)
          sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_j(p)
     END DO sommaA_do

     v_sA(i)=v_c0(i)*sommaA

     !-----------
     !Calcolo Bvt
     !----------- 
     !Calcolo della sommatoria e di Bvt
     sommaB=(0.0D0,0.0D0)
     sommaB_do:DO q=1,qmax
          p=n+nu-2*q
          pr=REAL(p,dbl)
          sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_j(p+1)
     END DO sommaB_do

     v_sB(i)=v_c0(i)*sommaB

END DO AB_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL


!CALL OMP_SET_NUM_THREADS(1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Loop for AB: openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!$OMP PARALLEL PRIVATE(n,m,m1,nu,iab)

!Choosing the chunk size
chunk=nnzAABB/(200*2)
chunk_if_AABB: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_AABB


!$OMP DO SCHEDULE(guided,chunk)
AABB_do: DO i=1,nnzAABB,2

     !----------------------------------------------------------------------------------------
     !Assigning all the needed indexes, which will be private, I assume, in the openMP version
     !----------------------------------------------------------------------------------------
     n=m_index_AABB(i,1)
     n=m_index_AABB(i+1,1)
     m=m_index_AABB(i,2)
     m=m_index_AABB(i+1,2)
     m1=m_index_AABB(i,3)
     m1=m_index_AABB(i+1,3)
     nu=m_index_AABB(i,4)
     nu=m_index_AABB(i+1,4)
     iab=m_index_AABB(i,5)
     iab=m_index_AABB(i+1,5)


     !assigning values to TU matrix elements
     m1_if: IF (m1==1) THEN
          v_sAABB(i)=v_sA(iab)
          v_sAABB(i+1)=v_sB(iab)
          v_ABsign(i)=0
          v_ABsign(i+1)=1
          !		WRITE(*,1004) v_sAABB(i),v_sAABB(i+1)
          !		1004 FORMAT (4ES14.5)
     ELSE
          v_sAABB(i)=v_sB(iab)
          v_sAABB(i+1)=v_sA(iab)
          v_ABsign(i)=1
          v_ABsign(i+1)=0
          !		WRITE(*,1005) v_sAABB(i),v_sAABB(i+1)
          !		1005 FORMAT (4ES14.5)
     END IF m1_if

END DO AABB_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!WRITE(*,*)
!WRITE(*,*) "A,B"
!WRITE(*,1006) v_sA(1),v_sB(1)
!1006 FORMAT (4ES14.5)
!WRITE(*,1007) v_sA(2),v_sB(2)
!1007 FORMAT (4ES14.5)
!WRITE(*,1008) v_sA(3),v_sB(3)
!1008 FORMAT (4ES14.5)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Loops for M: openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!Loop for M row vector filling
v_iM(1)=1				!Diagonal is filled,therefore first line begins at 1

row_M_do: DO i=2,Mside
     half_row_if: IF (i<=(AABBside+1)) THEN
          v_iM(i)=v_iM(i-1)+(v_iAABB(i)-v_iAABB(i-1)+1)
     ELSE
          v_iM(i)=v_iM(i-1)+(v_iAABB(i-AABBside)-v_iAABB(i-AABBside-1)+1)
     END IF half_row_if
END DO row_M_do
v_iM(Mside+1)=nnzM+1		!+1 Last row begins at non zero elements+1

!WRITE(*,*) v_iM

!WRITE(*,*) 'qui2'


!Loop for M column vector filling
col_M_do: DO i=1,Mside

     half_col_if: IF (i<=AABBside) THEN		!I am in the first half

          !		WRITE(*,*) 'qui1a',v_iAABB(i),v_iAABB(i+1)-1,v_iM(i)+1,v_iM(i+1)-1
          !Diagonal term, the first of the row
          v_jM(v_iM(i))=i
          !Remainder of the line from the of diagonal block AB + a shift of half the matrix side
          v_jM(v_iM(i)+1:v_iM(i+1)-1)=AABBside+v_jAABB(v_iAABB(i):v_iAABB(i+1)-1)
          !		WRITE(*,*) 'qui1b'

     ELSE					!I am in the second half

          !		WRITE(*,*) 'qui2a',v_iAABB(i-AABBside),v_iAABB(i+1-AABBside)-1,v_iM(i),v_iM(i+1)-2
          !Diagonal term, now it is the last of the row...
          v_jM(v_iM(i+1)-1)=i
          !Beginning of the line from the of diagonal block AB
          v_jM(v_iM(i):v_iM(i+1)-2)=v_jAABB(v_iAABB(i-AABBside):(v_iAABB(i+1-AABBside)-1))
          !		WRITE(*,*) 'qui2b',v_iM(i),v_iM(i+1)-2

     END IF half_col_if

END DO col_M_do
!WRITE(*,*) 'qui3'
!Loop for M values filling
dz=m_xyz(1,3)-m_xyz(2,3)		!dz coordinate

!I build a cursor to get the diagonal correctly
!Not very elegant but...
ALLOCATE(v_n(1:2*nstop*(nstop+2)))
i=0
cursor_n_do: DO n=1,nstop
     cursor_m_do: DO m=-n,n

          i=i+1
          v_n(2*i-1)=2*n-1
          v_n(2*i)=2*n

     END DO cursor_m_do
END DO cursor_n_do

!WRITE(*,*) 'v_iM',v_iM
!WRITE(*,*)
!WRITE(*,*) 'v_jM',v_jM

! CALL OMP_SET_NUM_THREADS(1)

!---------------------------------------------------------
!Finally Filling the top half of the matrix
!---------------------------------------------------------

!Initializing indexing
i=0
jjt=1

!Top index initialization
val_i_top_do: DO j=1,AABBside+nnzAABB

     diag_i_top_if: IF (i==v_jM(j)-1) THEN
          i=i+1
          v_index_i(j)=i
          v_index_jj(j)=jjt
          ! WRITE(*,*)"A,i,j,jjt",i,j,jjt
     ELSE
          v_index_i(j)=i
          v_index_jj(j)=jjt
          ! WRITE(*,*)"B,i,j,jjt",i,j,jjt
          jjt=jjt+1
     END IF diag_i_top_if

END DO val_i_top_do

!$OMP PARALLEL PRIVATE(i,j,jjt,n,nu)

!Choosing the chunk size
chunk=nnzM/(200*2)
chunk_if_top_M: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_top_M


!$OMP DO SCHEDULE(guided,chunk)
val_m_M_top_do: DO j=1,AABBside+nnzAABB

     !Setting Indexes
     i=v_index_i(j)
     jjt=v_index_jj(j)
     !	WRITE(*,*)"Top,i,j,n,jjt",i,j,v_n(i),jjt

     diag_m_top_if: IF (i==v_jM(j)) THEN			!I am on the diagonal
          v_sM(j)=v_Rn(v_n(i))
     ELSE							!I am off diagonal, the remainder of the row
          v_sM(j)=v_sAABB(jjt)

     END IF diag_m_top_if

END DO val_m_M_top_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!---------------------------------------------------------
!Finally Filling the bottom half of the matrix
!---------------------------------------------------------

i=AABBside+1
jjb=0

!Bottom index initialization
val_i_bottom_do: DO j=AABBside+nnzAABB+1,nnzM

     diag_i_bottom_if: IF (i==v_jM(j)) THEN
          v_index_i(j-(AABBside+nnzAABB))=i
          v_index_jj(j-(AABBside+nnzAABB))=jjb
          ! WRITE(*,*)"C,i,j,jjb",i,j,jjb
          i=i+1
     ELSE
          jjb=jjb+1
          v_index_i(j-(AABBside+nnzAABB))=i
          v_index_jj(j-(AABBside+nnzAABB))=jjb
          ! WRITE(*,*)"D,i,j,jjb",i,j,jjb
     END IF diag_i_bottom_if

END DO val_i_bottom_do

!$OMP PARALLEL PRIVATE(i,j,jjb,n,nu)

!Choosing the chunk size
chunk=nnzM/(200*2)
chunk_if_bottom_M: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_bottom_M

!$OMP DO SCHEDULE(guided,chunk)
val_m_M_bottom_do: DO j=AABBside+nnzAABB+1,nnzM

     !Initializing indexes
     i=v_index_i(j-(AABBside+nnzAABB))
     jjb=v_index_jj(j-(AABBside+nnzAABB))
     n=m_index_AABB(jjb,1)
     nu=m_index_AABB(jjb,4)
     ! WRITE(*,*)"Bottom,i,j,n,jjb",i,j,v_n(i-AABBSide),jjb

     diag_m_bottom_if: IF (i==v_jM(j)) THEN	!I am on the diagonal
          v_sM(j)=v_Vn(v_n(i-AABBside))
     ELSE					!I am off diagonal, the remainder of the row
          v_sM(j)=v_Zn(v_n(i-AABBside))*(((-1.0D0)**(v_ABsign(jjb)))*((-1.0D0)**(n+nu)))*v_sAABB(jjb)

     END IF diag_m_bottom_if

END DO val_m_M_bottom_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

DEALLOCATE(v_n,v_ABsign,v_index_i,v_index_jj)

END SUBROUTINE fill_M







!******************************************************************************
!2.1) SUBROUTINE fill_M_dip: Filling the coefficient matrix for the linear system
! to solve all the expansion coefficients 
!******************************************************************************
SUBROUTINE fill_M_dip(nstop,kr,m_xyz,v_c0,v_aq,v_bq,v_Rn,v_Vn,v_Zn,m_index,m_index_AABB,m_Apmn,v_iAB,v_jAB,v_iAABB,v_jAABB,&
     & v_sA,v_sB,v_sAABB,v_iM,v_jM,v_sM,error)

IMPLICIT NONE

! Dummy arguments
INTEGER(lo), INTENT(IN) :: nstop							! Multipolar expansions
COMPLEX(dbl), INTENT(IN) :: kr								! Radial distance kr_ij
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_c0,v_aq,v_bq					! Norm,gaunt e bq vectors
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						! Core and host Coordinates
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_Rn,v_Vn,v_Zn				! M auxiliary vectors
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_index,m_index_AABB			! Index matrix
REAL(lo), DIMENSION(:,:), INTENT(IN) :: m_Apmn					! Apmn coefficient matrix
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iAB,v_jAB,v_iAABB,v_jAABB			! CSR row column vectors for A,B,AB
INTEGER(lo), DIMENSION(:), INTENT(OUT) :: v_iM,v_jM					! CSR row column vectors for M
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_sA,v_sB,v_sAABB,v_sM			! CSR value matrixes
INTEGER(lo), INTENT(OUT) :: error							! Flag di errore

! Internal variables
INTEGER(lo) :: m,m1,n,nu,i,j,jjt,jjb,iab,q,p,low_nu,qmax,expb				!Loop indexes
REAL(dbl) :: mr,nr,mur,nur,pr,dz							!Loop indexes real type, and delta z coordinate
INTEGER(lo) :: low_aq,up_aq,low_bq,up_bq						!Bounds for aq,bq vectors
INTEGER(lo) :: pmin,pmax,nbes								!Bounds for Bessel vectors
INTEGER(lo) :: nnzAB,nnzAABB,nnzM							!Elements of (A,B), AB, M matrix
INTEGER(lo) :: blockside,AABBside,Mside						!Sides of (A,B), AB, M matrix
COMPLEX(dbl) :: sommaA,sommaB								!Partial Sums for VTC
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_j						!J Bessel
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_n,v_ABsign,v_index_i,v_index_jj		!n vector and and signs for A and B
INTEGER(lo) :: chunk=1								!Thread id, number of total threads and chunk size

! Subroutine

!Index and bounds initialization

!nnzAB: non-zero matrix elements of A or B
!blockside: side of A,B blocks
i=0
nnzAB=nstop*(1+2*nstop*(3+nstop))
nnzAB=nnzAB/3
blockside=nstop*(nstop+2)

!nnzABAB: non-zero matrix elements of AB
!AABBside: side of the off diagonal AB block
nnzAABB=4*nnzAB
AABBside=2*blockside

!nnzM: non-zero matrix elements of M
!Mside: side of the coefficient matrix
nnzM=2*nnzAABB+2*AABBside	!2 off diagonal blocks plus long diagonal
Mside=2*AABBside


!aq
low_aq=1
up_aq=1

!bq
low_bq=1
up_bq=1

!Allocation and computation of Bessel functions
ALLOCATE(v_j(0:(2*nstop+1)))
CALL besselj_z_sub((2*nstop+1),kr,v_j,error)

!Allocation of the matrix for the signs of A and B, and for the running index i
ALLOCATE(v_ABsign(1:nnzAABB),v_index_i(1:(AABBSide+nnzAABB)),v_index_jj(1:(AABBSide+nnzAABB)))


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Loop for A and B: openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!CALL OMP_SET_NUM_THREADS(1)

!$OMP PARALLEL PRIVATE(n,m,nu,qmax,pmin,pmax,low_aq,up_aq,low_bq,up_bq,q,p,pr,sommaA,sommaB)

!Choosing the chunk size
chunk=nnzAB/(200*2)
chunk_if_AB: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_AB
!WRITE (*,*) 'Chunk length =',chunk

!$OMP DO SCHEDULE(guided,chunk)
AB_do: DO i=1,nnzAB

     !----------------------------------------------------------------------------------------
     !Assigning all the needed indexes, which will be private, I assume, in the openMP version
     !----------------------------------------------------------------------------------------
     n=m_index(i,1)
     m=m_index(i,2)
     nu=m_index(i,3)
     qmax=m_index(i,4)
     pmin=m_index(i,5)
     pmax=m_index(i,6)
     low_aq=m_index(i,8)
     up_aq=m_index(i,9)
     low_bq=m_index(i,10)
     up_bq=m_index(i,11)

     !-----------
     !Calcolo Avt
     !-----------
     !Calcolo della sommatoria e di Avt
     sommaA=(0.0D0,0.0D0)
     sommaA_do:DO q=0,qmax
          p=n+nu-2*q
          pr=REAL(p,dbl)
          sommaA=sommaA+v_ic(p)*m_Apmn(i,p+1)*v_aq(low_aq+q)*v_j(p)
     END DO sommaA_do

     v_sA(i)=v_c0(i)*sommaA

     !-----------
     !Calcolo Bvt
     !----------- 
     !Calcolo della sommatoria e di Bvt
     sommaB=(0.0D0,0.0D0)
     sommaB_do:DO q=1,qmax
          p=n+nu-2*q
          pr=REAL(p,dbl)
          sommaB=sommaB+v_ic(p+1)*v_bq(low_bq+(q-1))*v_j(p+1)
     END DO sommaB_do

     v_sB(i)=v_c0(i)*sommaB

END DO AB_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL


!CALL OMP_SET_NUM_THREADS(1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Loop for AB: openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!$OMP PARALLEL PRIVATE(n,m,m1,nu,iab)

!Choosing the chunk size
chunk=nnzAABB/(200*2)
chunk_if_AABB: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_AABB


!$OMP DO SCHEDULE(guided,chunk)
AABB_do: DO i=1,nnzAABB,2

     !----------------------------------------------------------------------------------------
     !Assigning all the needed indexes, which will be private, I assume, in the openMP version
     !----------------------------------------------------------------------------------------
     n=m_index_AABB(i,1)
     n=m_index_AABB(i+1,1)
     m=m_index_AABB(i,2)
     m=m_index_AABB(i+1,2)
     m1=m_index_AABB(i,3)
     m1=m_index_AABB(i+1,3)
     nu=m_index_AABB(i,4)
     nu=m_index_AABB(i+1,4)
     iab=m_index_AABB(i,5)
     iab=m_index_AABB(i+1,5)


     !assigning values to TU matrix elements
     m1_if: IF (m1==1) THEN
          v_sAABB(i)=v_sA(iab)
          v_sAABB(i+1)=v_sB(iab)
          v_ABsign(i)=0
          v_ABsign(i+1)=1
          !		WRITE(*,1004) v_sAABB(i),v_sAABB(i+1)
          !		1004 FORMAT (4ES14.5)
     ELSE
          v_sAABB(i)=v_sB(iab)
          v_sAABB(i+1)=v_sA(iab)
          v_ABsign(i)=1
          v_ABsign(i+1)=0
          !		WRITE(*,1005) v_sAABB(i),v_sAABB(i+1)
          !		1005 FORMAT (4ES14.5)
     END IF m1_if

END DO AABB_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!WRITE(*,*)
!WRITE(*,*) "A,B"
!WRITE(*,1006) v_sA(1),v_sB(1)
!1006 FORMAT (4ES14.5)
!WRITE(*,1007) v_sA(2),v_sB(2)
!1007 FORMAT (4ES14.5)
!WRITE(*,1008) v_sA(3),v_sB(3)
!1008 FORMAT (4ES14.5)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Loops for M: openMP friendly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!Loop for M row vector filling
v_iM(1)=1				!Diagonal is filled,therefore first line begins at 1

row_M_do: DO i=2,Mside
     half_row_if: IF (i<=(AABBside+1)) THEN
          v_iM(i)=v_iM(i-1)+(v_iAABB(i)-v_iAABB(i-1)+1)
     ELSE
          v_iM(i)=v_iM(i-1)+(v_iAABB(i-AABBside)-v_iAABB(i-AABBside-1)+1)
     END IF half_row_if
END DO row_M_do
v_iM(Mside+1)=nnzM+1		!+1 Last row begins at non zero elements+1

!WRITE(*,*) v_iM

!WRITE(*,*) 'qui2'


!Loop for M column vector filling
col_M_do: DO i=1,Mside

     half_col_if: IF (i<=AABBside) THEN		!I am in the first half

          !		WRITE(*,*) 'qui1a',v_iAABB(i),v_iAABB(i+1)-1,v_iM(i)+1,v_iM(i+1)-1
          !Diagonal term, the first of the row
          v_jM(v_iM(i))=i
          !Remainder of the line from the of diagonal block AB + a shift of half the matrix side
          v_jM(v_iM(i)+1:v_iM(i+1)-1)=AABBside+v_jAABB(v_iAABB(i):v_iAABB(i+1)-1)
          !		WRITE(*,*) 'qui1b'

     ELSE					!I am in the second half

          !		WRITE(*,*) 'qui2a',v_iAABB(i-AABBside),v_iAABB(i+1-AABBside)-1,v_iM(i),v_iM(i+1)-2
          !Diagonal term, now it is the last of the row...
          v_jM(v_iM(i+1)-1)=i
          !Beginning of the line from the of diagonal block AB
          v_jM(v_iM(i):v_iM(i+1)-2)=v_jAABB(v_iAABB(i-AABBside):(v_iAABB(i+1-AABBside)-1))
          !		WRITE(*,*) 'qui2b',v_iM(i),v_iM(i+1)-2

     END IF half_col_if

END DO col_M_do
!WRITE(*,*) 'qui3'
!Loop for M values filling
dz=m_xyz(1,3)-m_xyz(2,3)		!dz coordinate

!I build a cursor to get the diagonal correctly
!Not very elegant but...
ALLOCATE(v_n(1:2*nstop*(nstop+2)))
i=0
cursor_n_do: DO n=1,nstop
     cursor_m_do: DO m=-n,n

          i=i+1
          v_n(2*i-1)=2*n-1
          v_n(2*i)=2*n

     END DO cursor_m_do
END DO cursor_n_do

!WRITE(*,*) 'v_iM',v_iM
!WRITE(*,*)
!WRITE(*,*) 'v_jM',v_jM

! CALL OMP_SET_NUM_THREADS(1)

!---------------------------------------------------------
!Finally Filling the top half of the matrix
!---------------------------------------------------------

!Initializing indexing
i=0
jjt=1

!Top index initialization
val_i_top_do: DO j=1,AABBside+nnzAABB

     diag_i_top_if: IF (i==v_jM(j)-1) THEN
          i=i+1
          v_index_i(j)=i
          v_index_jj(j)=jjt
          ! WRITE(*,*)"A,i,j,jjt",i,j,jjt
     ELSE
          v_index_i(j)=i
          v_index_jj(j)=jjt
          ! WRITE(*,*)"B,i,j,jjt",i,j,jjt
          jjt=jjt+1
     END IF diag_i_top_if

END DO val_i_top_do

!$OMP PARALLEL PRIVATE(i,j,jjt,n,nu)

!Choosing the chunk size
chunk=nnzM/(200*2)
chunk_if_top_M: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_top_M


!$OMP DO SCHEDULE(guided,chunk)
val_m_M_top_do: DO j=1,AABBside+nnzAABB

     !Setting Indexes
     i=v_index_i(j)
     jjt=v_index_jj(j)
     !	WRITE(*,*)"Top,i,j,n,jjt",i,j,v_n(i),jjt

     diag_m_top_if: IF (i==v_jM(j)) THEN			!I am on the diagonal
          v_sM(j)=v_Zn(v_n(i))
     ELSE							!I am off diagonal, the remainder of the row
          v_sM(j)=v_Vn(v_n(i))*v_sAABB(jjt)

     END IF diag_m_top_if

END DO val_m_M_top_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!---------------------------------------------------------
!Finally Filling the bottom half of the matrix
!---------------------------------------------------------

i=AABBside+1
jjb=0

!Bottom index initialization
val_i_bottom_do: DO j=AABBside+nnzAABB+1,nnzM

     diag_i_bottom_if: IF (i==v_jM(j)) THEN
          v_index_i(j-(AABBside+nnzAABB))=i
          v_index_jj(j-(AABBside+nnzAABB))=jjb
          ! WRITE(*,*)"C,i,j,jjb",i,j,jjb
          i=i+1
     ELSE
          jjb=jjb+1
          v_index_i(j-(AABBside+nnzAABB))=i
          v_index_jj(j-(AABBside+nnzAABB))=jjb
          ! WRITE(*,*)"D,i,j,jjb",i,j,jjb
     END IF diag_i_bottom_if

END DO val_i_bottom_do

!$OMP PARALLEL PRIVATE(i,j,jjb,n,nu)

!Choosing the chunk size
chunk=nnzM/(200*2)
chunk_if_bottom_M: IF (chunk==0) THEN
     chunk=1
END IF chunk_if_bottom_M

!$OMP DO SCHEDULE(guided,chunk)
val_m_M_bottom_do: DO j=AABBside+nnzAABB+1,nnzM

     !Initializing indexes
     i=v_index_i(j-(AABBside+nnzAABB))
     jjb=v_index_jj(j-(AABBside+nnzAABB))
     n=m_index_AABB(jjb,1)
     nu=m_index_AABB(jjb,4)
     ! WRITE(*,*)"Bottom,i,j,n,jjb",i,j,v_n(i-AABBSide),jjb

     diag_m_bottom_if: IF (i==v_jM(j)) THEN	!I am on the diagonal
          v_sM(j)=1.0D0
     ELSE					!I am off diagonal, the remainder of the row
          v_sM(j)=v_Rn(v_n(i-AABBside))*(((-1.0D0)**(v_ABsign(jjb)))*((-1.0D0)**(n+nu)))*v_sAABB(jjb)

     END IF diag_m_bottom_if

END DO val_m_M_bottom_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

DEALLOCATE(v_n,v_ABsign,v_index_i,v_index_jj)

END SUBROUTINE fill_M_dip

















!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE field_exp_shell_borghese: subroutine to get all the expansions
! coefficients for the shell, in every reference frame
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE field_exp_shell_borghese(nstop,neq,matrixside,lambda,ref_index,r_ih,v_req,m_epseq,&
     & v_p,v_ab_shell2,v_dc_shell1,v_Rn,v_Vn,v_Zn,&
     & v_ab_shell1,v_dc_shell2,v_dc_core,v_ab_host,error)

IMPLICIT NONE

! Dummy argument
INTEGER(lo), INTENT(IN) :: nstop,neq,matrixside			! Multipole expansion, equal spheres and matrixsize
REAL(dbl), INTENT(IN) :: lambda,ref_index,r_ih				! Wavelength, ref index and core-host distance
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req				! Equivalent radii vector
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq				! Equivalent dielectric functions matrix
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p				! Incident fiyeld expansion coefficients
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc_shell1,v_ab_shell2	! Vectors: d_mn1,c_mn1,a_mn2,b_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_Rn,v_Vn,v_Zn		! Auxiliary coefficients
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_shell1,v_dc_shell2	! Vectors: a_mn1,b_mn1,d_mn2,c_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_host,v_dc_core		! Vectors: a_mn,b_mn,d_mn,c_mn
INTEGER(lo), INTENT(OUT) :: error					! Error Flag


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x						! Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc						! Complex dielectric function
COMPLEX(dbl), DIMENSION(0:neq) :: v_k						! Comples wave vector
COMPLEX(dbl), DIMENSION(0:neq) :: v_m						! Normalized complex ref index
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx					! Argument Psi and Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z01,v_psi_z11,v_psi_z12,v_psi_z22	! Complex psi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12		! Complex csi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z01,v_psi1_z11,v_psi1_z12,v_psi1_z22	! Complex psi derivatives
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z11,v_csi1_z12			! Complex pcsi derivatives
COMPLEX(dbl), DIMENSION(1:nstop) :: v_qintd,v_qintc				! Auxiliary coefficients
COMPLEX(dbl), DIMENSION(1:2*nstop) :: v_Vscatn,v_Zscatn				! Auxiliary coefficients
INTEGER(lo) :: m,n,nu,i,j							! Indexes
REAL(dbl) :: nr									! Real indexes



! Subroutine

!Complex dielectric function
v_epsc_do: DO j=1,neq
     !Using neq-j,because i use out to in numbering in Ngo paper (0=host,1=shell,2=core), but in m_epsqe
     !i do the opposite (1=core,2=shell). so i use Ngo numbering here
     v_epsc(j)=CMPLX(m_epseq(neq+1-j,1),m_epseq(neq+1-j,2),KIND=dbl)
END DO v_epsc_do

!Size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Wavevectors
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq)=(2*pi_d*SQRT(v_epsc))/lambda

!Normalized refractive indexes
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq)=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)

!Arguments
mx_out_do: DO i=0,neq
     mx_in_do: DO j=1,neq

          m_mx(i,j)=v_m(i)*CMPLX(v_x(j),KIND=dbl)

     END DO mx_in_do
END DO mx_out_do


!Riccati bessel psi
CALL psi_z_sub(nstop,m_mx(0,1),v_psi_z01,error)
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(1,2),v_psi_z12,error)
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,10)
10   FORMAT ("Error: error in psi_z_sub called by field_exp_shell_borghese")
     RETURN
END IF psi_z_if

!Riccati bessel csi
CALL csi_z_sub(nstop,m_mx(0,1),v_csi_z01,error)
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(1,2),v_csi_z12,error)

csi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,20)
20   FORMAT ("Error: error in csi_z_sub called by field_exp_shell_borghese")
     RETURN
END IF csi_z_if

!Riccati bessel derivatives
der_do: DO n=1,nstop

     nr=REAL(n,dbl)
     !D_Psi
     v_psi1_z01(n)=v_psi_z01(n-1)-nr*v_psi_z01(n)/m_mx(0,1)
     v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
     v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)
     v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

     !Derivate di csi
     v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
     v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)

END DO der_do


!Auxiliary vectors to get the expansion coefficients for the fields of the internal core
qint_do:DO n=1,nstop

     v_qintd(n)= ( v_k(2)*v_psi1_z12(n)*v_csi_z12(n) - v_k(2)*v_psi_z12(n)*v_csi1_z12(n) ) / &
          & ( v_k(1)*v_psi1_z22(n)*v_csi_z12(n) - v_k(2)*v_psi_z22(n)*v_csi1_z12(n) )

     v_qintc(n)= ( v_k(2)*v_psi1_z12(n)*v_csi_z12(n) - v_k(2)*v_psi_z12(n)*v_csi1_z12(n) ) / &
          & ( v_k(2)*v_psi1_z22(n)*v_csi_z12(n) - v_k(1)*v_psi_z22(n)*v_csi1_z12(n) )

END DO qint_do

!-----------------------------
!Vectors: a_mn2,b_mn2
!-----------------------------
j=0
dmn2_n_do: DO n=1,nstop
     dmn2_m_do: DO m=-n,n

          !amn2,bmn2
          j=j+1
          v_dc_shell2(2*j-1)=-v_ab_shell2(2*j-1)*v_Rn(2*n-1)
          v_dc_shell2(2*j)=-v_ab_shell2(2*j)*v_Rn(2*n)

     END DO dmn2_m_do
END DO dmn2_n_do


!-----------------------------
!Vectors: a_mn1,b_mn1
!-----------------------------
j=0
amn1_n_do: DO n=1,nstop
     amn1_m_do: DO m=-n,n

          !amn1,bmn1
          j=j+1
          v_ab_shell1(2*j-1)=(v_p(2*j-1)-v_dc_shell1(2*j-1)*v_Vn(2*n-1))/v_Zn(2*n-1)
          v_ab_shell1(2*j)=(v_p(2*j)-v_dc_shell1(2*j)*v_Vn(2*n))/v_Zn(2*n)

     END DO amn1_m_do
END DO amn1_n_do

!---------------------------------------------------
!Expansion coefficients for the core
!---------------------------------------------------
j=0
dmn_n_do: DO n=1,nstop
     dmn_m_do: DO m=-n,n

          !dmn,cmn
          j=j+1
          v_dc_core(2*j-1)=v_qintd(n)*v_dc_shell2(2*j-1)
          v_dc_core(2*j)=v_qintc(n)*v_dc_shell2(2*j)

     END DO dmn_m_do
END DO dmn_n_do



!---------------------------------------------------
!Expansion coefficients for the scattered field
!---------------------------------------------------
!Coefficients for the calculation of scatterd field
!Tn
v_Vscatn(1:2*nstop:2)=(0.0D0,1.0D0)*( (v_k(0)/v_k(1))*v_psi_z01(1:nstop)*v_psi1_z11 - v_psi1_z01*v_psi_z11(1:nstop))
v_Vscatn(2:2*nstop:2)=(0.0D0,1.0D0)*( v_psi_z01(1:nstop)*v_psi1_z11 - (v_k(0)/v_k(1))*v_psi1_z01*v_psi_z11(1:nstop))

!Un
v_Zscatn(1:2*nstop:2)=(0.0D0,1.0D0)*( v_psi1_z01*v_csi_z11(1:nstop) - (v_k(0)/v_k(1))*v_psi_z01(1:nstop)*v_csi1_z11)
v_Zscatn(2:2*nstop:2)=(0.0D0,1.0D0)*( (v_k(0)/v_k(1))*v_psi1_z01*v_csi_z11(1:nstop) - v_psi_z01(1:nstop)*v_csi1_z11 )

j=0
amn_n_do: DO n=1,nstop
     amn_m_do: DO m=-n,n

          !amn,bmn
          j=j+1
          v_ab_host(2*j-1)=v_ab_shell1(2*j-1)*v_Zscatn(2*n-1) + v_dc_shell1(2*j-1)*v_Vscatn(2*n-1)

          v_ab_host(2*j)=v_ab_shell1(2*j)*v_Zscatn(2*n) + v_dc_shell1(2*j)*v_Vscatn(2*n)

     END DO amn_m_do
END DO amn_n_do

END SUBROUTINE field_exp_shell_borghese



!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE coeff_shell: si calcolano i vettori dei coefficienti di per passare
! tra i diversi coefficienti della shell
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE field_exp_shell_dip_borghese(nstop,neq,matrixside,lambda,ref_index,r_ih,v_req,m_epseq,&
     & v_p,v_ab_shell2,v_dc_shell1,v_Rn,v_Vn,v_Zn,&
     & v_ab_shell1,v_dc_shell2,v_dc_core,v_ab_host,error)

IMPLICIT NONE

! Dummy argument
INTEGER(lo), INTENT(IN) :: nstop,neq,matrixside			! Multipole expansion, equal spheres and matrixsize
REAL(dbl), INTENT(IN) :: lambda,ref_index,r_ih				! Wavelength, ref index and core-host distance
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req				! Equivalent radii vector
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq				! Equivalent dielectric functions matrix
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p				! Incident fiyeld expansion coefficients
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc_shell1,v_ab_shell2	! Vectors: d_mn1,c_mn1,a_mn2,b_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_Rn,v_Vn,v_Zn		! Auxiliary coefficients

COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_shell1,v_dc_shell2	! Vectors: a_mn1,b_mn1,d_mn2,c_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_host,v_dc_core		! Vectors: a_mn,b_mn,d_mn,c_mn
INTEGER(lo), INTENT(OUT) :: error					! Error Flag


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq-1) :: v_x							! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq-1) :: v_epsc						! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_k							! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_m							! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq-1,1:neq-1) :: m_mx						! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z12,v_psi_z22			! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12,v_csi_z22	        ! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11						! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z01,v_csi1_z11				! Vettori Csi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_qintd,v_qintc					! Vettore coefficienti ausiliari
INTEGER(lo) :: m,n,nu,i,j								! Indici
REAL(dbl) :: nr										! Indice complessificato



! Inizio della procedura vera e propria

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq-1
     !Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
     !da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
     v_epsc(j)=CMPLX(m_epseq((neq-1)+1-j,1),m_epseq((neq-1)+1-j,2),KIND=dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
!Calcolo il vettore per i size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Calcolo il vettore per i wavevectors
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq-1)=(2*pi_d*SQRT(v_epsc))/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq-1)=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)

!Calcolo la matrice mx
mx_out_do: DO i=0,neq-1
     mx_in_do: DO j=1,neq-1

          m_mx(i,j)=v_m(i)*CMPLX(v_x(j),KIND=dbl)

     END DO mx_in_do
END DO mx_out_do

!Loop per calcolare le matrici coefficiente a e b 

!Calcolo riccati bessel di prima specie Complessa
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(1,2),v_psi_z12,error)
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,10)
10   FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_shell")
     RETURN
END IF psi_z_if

!Calcolo riccati bessel di terza specie complesso
CALL csi_z_sub(nstop,m_mx(0,1),v_csi_z01,error)
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(1,2),v_csi_z12,error)
CALL csi_z_sub(nstop,m_mx(2,2),v_csi_z22,error)

csi_z_if: IF (error/=0) THEN
     error=1
     WRITE(*,20)
20   FORMAT ("Errore: errore in csi_z_sub chiamata da coeff_shell")
     RETURN
END IF csi_z_if

!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
der_do: DO n=1,nstop

     nr=REAL(n,dbl)
     !Derivate di Psi
     v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)

     !Derivate di csi
     v_csi1_z01(n)=v_csi_z01(n-1)-nr*v_csi_z01(n)/m_mx(0,1)
     v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)

END DO der_do


!Cofficients for the computation of overall scattering coefficients
qint_do:DO n=1,nstop

     v_qintd(n)= (0.0D0,1.0D0) / &
          & ( (v_k(1)/v_k(0))*v_csi1_z01(n)*v_csi_z11(n) - v_csi_z01(n)*v_csi1_z11(n) )

     v_qintc(n)= (0.0D0,1.0D0) / &
          & ( v_csi1_z01(n)*v_csi_z11(n) - (v_k(1)/v_k(0))*v_csi_z01(n)*v_csi1_z11(n) )

END DO qint_do

!-----------------------------
!Vectors: a_mn1,b_mn1
!-----------------------------
j=0
dmn1_n_do: DO n=1,nstop
     dmn1_m_do: DO m=-n,n

          !amn1,bmn1
          j=j+1
          v_ab_shell1(2*j-1)=-v_dc_shell1(2*j-1)/v_Rn(2*n-1)
          v_ab_shell1(2*j)=-v_dc_shell1(2*j)/v_Rn(2*n)

     END DO dmn1_m_do
END DO dmn1_n_do

!-----------------------------
!Vectors: d_mn2,c_mn2
!-----------------------------
j=0
dmn2_n_do: DO n=1,nstop
     dmn2_m_do: DO m=-n,n

          !dmn2,cmn2
          j=j+1
          v_dc_shell2(2*j-1)=(v_p(2*j-1)-v_ab_shell2(2*j-1)*v_Zn(2*n-1))/v_Vn(2*n-1)
          v_dc_shell2(2*j)=(v_p(2*j)-v_ab_shell2(2*j)*v_Zn(2*n))/v_Vn(2*n)

     END DO dmn2_m_do
END DO dmn2_n_do



!---------------------------------------------------
!Vectors a_mn e b_mn at the host interface
!---------------------------------------------------
j=0
amn_n_do: DO n=1,nstop
     amn_m_do: DO m=-n,n

          !amn bmn host
          j=j+1
          v_ab_host(2*j-1)=v_qintd(n)*v_dc_shell1(2*j-1)
          v_ab_host(2*j)=v_qintc(n)*v_dc_shell1(2*j)

     END DO amn_m_do
END DO amn_n_do



!---------------------------------------------------
!Vectors  d_mn e c_mn at the core interface
!---------------------------------------------------
j=0
dmn_n_do: DO n=1,nstop
     dmn_m_do: DO m=-n,n

          !dmn2 and cmn2
          j=j+1
          v_dc_core(2*j-1)=-( v_ab_shell2(2*j-1)*v_csi_z12(n) - v_dc_shell2(2*j-1)*v_psi_z12(n) - v_p(2*j-1)*v_csi_z22(n))&
               & / (v_psi_z22(n))
          v_dc_core(2*j)=-( v_ab_shell2(2*j)*v_k(2)*v_csi_z12(n) - v_dc_shell2(2*j)*v_k(2)*v_psi_z12(n) -  & 
               v_p(2*j)*v_k(1)*v_csi_z22(n)) / (v_k(1)*v_psi_z22(n))

     END DO dmn_m_do
END DO dmn_n_do



END SUBROUTINE field_exp_shell_dip_borghese




END MODULE vec_trans










