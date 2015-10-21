MODULE sing_part
	

USE kinds
USE basicsubs

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! INDICE DELLE SUBROUTINE CONTENUTE NEL MODULO
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Indice Subroutines
!1) nstop_sp(lambda,ref_index,v_req,m_epseq,nstop,nmx)
!	COMMENT: calculates number of multipolar expansions for single particles
!
!	VARS:
! 		lambda,ref_index: actual wavelength and refraction index of the embedding medium INPUT
!		v_req: vector containing radii of equivalent spheres INPUT
! 		m_epseq: matrix containing dielectric function for equivalent spheres INPUT
! 		nstop,nmx: numbers for multipolar expansion of single particles OUTPUT
!
!	ERROR: no error flags
!
!2) coeff_sp(lambda,ref_index,v_req,m_epseq,nstop,nmx,neq,m_a,m_b,error)
!	COMMENT: calculates coefficents a and b of multipolar expansions for single particles
!
!	VARS:
! 		lambda,ref_index: actual wavelength and refraction index of the embedding medium INPUT
!		v_req: vector containing radii of equivalent spheres INPUT
! 		m_epseq: matrix containing dielectric function for equivalent spheres INPUT
! 		nstop,nmx: numbers for multipolar expansion of single particles INPUT
!		neq: number of equivalent sphers in order to allocate arrays INPUT
!		m_a,m_b: matrice of multipolar expansion coefficents OUTPUT
!
!	ERROR: 
!		error=0 OK
!		error=1 Underflow in computing bessel functions
!
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


CONTAINS



!******************************************************************************
!******************************************************************************
!******************************************************************************
!1) SUBROUTINE besselj_d_sub(nstop,x,v_besj,error): 
!			   Si calcolano le funzioni sferiche di bessel del primo ordine
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE besselj_d_sub(nstop,x,v_besj,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: x								!Argomento reale
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_besj	!Vettore funzioni di bessel
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n									!Numero di coeff pn e indice
REAL(dbl) :: pn,nr										!Coefficiente rapporto e indice reale
COMPLEX(dbl) :: z										!x complessificato

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
v_besj=0.0D0
z=CMPLX(x,KIND=dbl)

!Controllo che l'argomento sia maggiore di zero
zero1_if: IF (x<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento in besselj_d_sub e' minore di zero")
				RETURN
END IF zero1_if

!Se x=0 allora metto solo un valore
zero2_if: IF (x==0.0D0) THEN
				v_besj(0)=1.0D0
				RETURN
END IF zero2_if

nstop_case: SELECT CASE (nstop)

CASE (0)

	!Calcolo il mio valore iniziale
	v_besj(0)=SIN(x)/x

CASE DEFAULT

	!Calcolo i valori iniziali a secondo dell'argomento
	v_besj(0)=SIN(x)/x

	!Calcolo i rapporti	
	CALL rn_d_sub(1,nstop,x,v_besj(1:nstop),error)

	!If su errore di rn
	rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore in rn_d_sub chiamata da besselj_d_sub")
				RETURN
	END IF rn_if
	
	!Calcolo bessel
	rn_do: DO n=1,nstop
		v_besj(n)=v_besj(n-1)/v_besj(n)
	END DO rn_do

END SELECT nstop_case

END SUBROUTINE besselj_d_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE besselj_d_sub(nstop,x,v_besj,error): 
!			   Si calcolano le funzioni sferiche di bessel del primo ordine
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE besselj_z_sub(nstop,z,v_besj,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop					!Numero di espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: z						!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_besj	!Vettore funzioni di bessel
INTEGER(lo), INTENT(OUT) :: error					!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n							!Numero di coeff pn e indice
REAL(dbl) :: pn,nr							!Coefficiente rapporto e indice reale

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
v_besj=(0.0D0,0.0D0)

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(z,dbl)<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento z besselj_z_sub e' errato")
				RETURN
END IF zero_if

!Controllo che la parte immaginaria dell'argomento non sia troppo grande
max_if: IF (ABS(AIMAG(z))>700.0D0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: l'argomento z besselj_z_sub e' > 700,overflow")
				RETURN
END IF max_if

!Se x=0 allora metto solo un valore
zero2_if: IF (ABS(z)==0.0D0) THEN
				v_besj(0)=1.0D0
				RETURN
END IF zero2_if

nstop_case: SELECT CASE (nstop)

CASE (0)

	!Calcolo il mio valore iniziale
	v_besj(0)=SIN(z)/z

CASE DEFAULT

	!Calcolo i valori iniziali a secondo dell'argomento
	v_besj(0)=SIN(z)/z

	!Calcolo i rapporti	
	CALL rn_z_sub(1,nstop,z,v_besj(1:nstop),error)

	!If su errore di rn
	rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore in rn_z_sub chiamata da besselj_z_sub")
				RETURN
	END IF rn_if
	
	!Calcolo bessel
	rn_do: DO n=1,nstop
		v_besj(n)=v_besj(n-1)/v_besj(n)
	END DO rn_do

END SELECT nstop_case

END SUBROUTINE besselj_z_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE bessely_d_sub(nstop,x,v_besy,error): 
!			   Si calcolano le funzioni sferiche di bessel del secondo ordine
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE bessely_d_sub(nstop,x,v_besy,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: x								!Argomento reale
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_besy	!Vettore funzioni di bessel
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n										!Indice
REAL(dbl) :: nr											!Espansione per sinx/x e indice reale

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
v_besy=0.0D0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (x<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento in bessely_d_sub e' minore di zero")
				RETURN
END IF zero_if

!Calcolo le mie funzioni a seconda del valore di nstop
nstop_case: SELECT CASE (nstop)

CASE (0) nstop_case

	v_besy(0)=-COS(x)/x
	
CASE (1) nstop_case

	v_besy(0)=-COS(x)/x
	v_besy(1)=-COS(x)/(x**2)-SIN(x)/x

CASE DEFAULT

	!Calcolo i valori iniziali a secondo dell'argomento
	v_besy(0)=-COS(x)/x
	v_besy(1)=-COS(x)/(x**2)-SIN(x)/x
	
	!Calcolo le funzioni sferiche di Bessel del primo ordine
	bes_do: DO n=2,nstop
	
		nr=REAL(n,dbl)
		v_besy(n)=(2.0D0*nr-1.0D0)*v_besy(n-1)/x - v_besy(n-2)
	
	END DO bes_do

END SELECT nstop_case

END SUBROUTINE bessely_d_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE bessely_z_sub(nstop,x,v_besy,error): 
!Si calcolano le funzioni sferiche di bessel del secondo ordine su arg complex
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE bessely_z_sub(nstop,z,v_besy,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Numero di espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: z							!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_besy		!Vettore funzioni di bessel
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n								!Indice
COMPLEX(dbl) :: nc								!Espansione per sinx/x e indice reale

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
v_besy=(0.0D0,0.0D0)

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(z,dbl)<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento in bessely_d_sub e' minore di zero")
				RETURN
END IF zero_if

!Calcolo le mie funzioni a seconda del valore di nstop
nstop_case: SELECT CASE (nstop)

CASE (0) nstop_case

	v_besy(0)=-COS(z)/z
	
CASE (1) nstop_case

	v_besy(0)=-COS(z)/z
	v_besy(1)=-COS(z)/(z**2)-SIN(z)/z

CASE DEFAULT

	!Calcolo i valori iniziali a secondo dell'argomento
	v_besy(0)=-COS(z)/z
	v_besy(1)=-COS(z)/(z**2)-SIN(z)/z
	
	!Calcolo le funzioni sferiche di Bessel del primo ordine
	bes_do: DO n=2,nstop
	
		nc=CMPLX(n,KIND=dbl)
		v_besy(n)=(2.0D0*nc-1.0D0)*v_besy(n-1)/z - v_besy(n-2)
	
	END DO bes_do

END SELECT nstop_case

END SUBROUTINE bessely_z_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!5) SUBROUTINE hankel1_d_sub(nstop,x,v_besy,error): 
!			   Si calcolano le funzioni sferiche di hankel del primo ordine
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE hankel1_d_sub(nstop,x,v_hank1,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop							!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: x									!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_hank1	!Vettore funzioni di hankel
INTEGER(lo), INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
REAL(dbl), DIMENSION(0:nstop) :: v_besj,v_besy 				!Vettori Besselj e Bessely

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
v_besy=0.0D0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (x<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento in hankel1_sub e' minore di zero")
				RETURN
END IF zero_if

!Calcolo della funzione sferiche di bessel prima e secondo genere
CALL besselj_d_sub(nstop,x,v_besj,error)

besj_if: IF (error/=0) THEN !Controllo besselj non sbagli
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore in besselj_d_sub chiamata da hankel1_d_sub")
				RETURN
END IF besj_if

CALL bessely_d_sub(nstop,x,v_besy,error)

besy_if: IF (error/=0) THEN !Controllo besselj non sbagli
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore in bessely_d_sub chiamata da hankel1_d_sub")
				RETURN
END IF besy_if

!Calcolo le funzioni sferiche di hankel1
v_hank1=CMPLX(v_besj,v_besy,KIND=dbl)

END SUBROUTINE hankel1_d_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!6) SUBROUTINE hankel1_z_sub(nstop,z,v_besy,error): 
! Si calcolano le funzioni sferiche di hankel del primo ordine per arg complex
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE hankel1_z_sub(nstop,z,v_hank1,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop					!Numero di espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: z						!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_hank1	!Vettore funzioni di hankel
INTEGER(lo), INTENT(OUT) :: error					!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl), DIMENSION(0:nstop) :: v_besj,v_besy		!Vettori Besselj e Bessely

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
v_besj=(0.0D0,0.0D0)
v_besy=(0.0D0,0.0D0)

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(z,dbl)<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento in hankel1_z_sub e' minore di zero")
				RETURN
END IF zero_if

!Calcolo della funzione sferiche di bessel prima e secondo genere
CALL besselj_z_sub(nstop,z,v_besj,error)

besj_if: IF (error/=0) THEN !Controllo besselj non sbagli
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore in besselj_z_sub chiamata da hankel1_z_sub")
				RETURN
END IF besj_if

CALL bessely_z_sub(nstop,z,v_besy,error)

besy_if: IF (error/=0) THEN !Controllo besselj non sbagli
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore in bessely_z_sub chiamata da hankel1_z_sub")
				RETURN
END IF besy_if

!Calcolo le funzioni sferiche di hankel1
v_hank1=v_besj+(0.0D0,1.0D0)*v_besy

END SUBROUTINE hankel1_z_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!7) SUBROUTINE eta_d_sub(nstop,r,v_eta,error): 
!			   Si calcolano le funzioni eta(n,r)=j(n,r)/r
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE eta_d_sub(nstop,r,v_eta,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: r								!Argomento reale
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_eta		!Vettore funzioni eta
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n									!Numero di coeff pn e indice
REAL(dbl) :: pn,nr										!Coefficiente rapporto e indice reale

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
!Metto v_eta(0)=0 anche se non e' vero,comunque non lo devo usare,e' per
!uniformita' di inizializzazione
v_eta=0.0D0

!Controllo che l'argomento sia maggiore di zero
zero1_if: IF (r<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: r in eta_d_sub e' minore di zero")
				RETURN
END IF zero1_if

!Se z=0 allora metto solo un valore
zero2_if: IF (r==0.0D0) THEN
				v_eta(1)=1.0D0/3.0D0
				RETURN
END IF zero2_if


nstop_case: SELECT CASE (nstop)

CASE (1)

	!Calcolo il mio valore iniziale
	in1_if: IF (r>0.1) THEN
		v_eta(1)=SIN(r)/(r**3)-COS(r)/(r**2)
	ELSE
		
		!Ciclo do per l'espansione in serie
		v_eta(1)=0.0D0
		in1_do: DO n=4,0,-1
			
			nr=REAL(n,dbl)
			v_eta(1)=v_eta(1)+((-1.0D0)**n)*2.0D0*(nr+1.0D0)*(r**(2*n))/fatt(2.0D0*nr+3.0D0)
		
		END DO in1_do

	END IF in1_if

CASE DEFAULT

	!Calcolo il mio valore iniziale
	in2_if: IF (r>0.1) THEN
	
		v_eta(1)=SIN(r)/(r**3)-COS(r)/(r**2)
		
	ELSE
		
		!Ciclo do per l'espansione in serie
		v_eta(1)=(0.0D0,0.0D0)
		in2_do: DO n=0,4
			
			nr=REAL(n,dbl)
			v_eta(1)=v_eta(1)+((-1.0D0)**n)*2.0D0*(nr+1.0D0)*(r**(2*n))/fatt(2.0D0*nr+3.0D0)
		
		END DO in2_do

	END IF in2_if
	
	!Calcolo i rapporti	
	CALL rn_d_sub(2,nstop,r,v_eta(2:nstop),error)

	!If su errore di rn
	rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore in rn_d_sub chiamata da eta_d_sub")
				RETURN
	END IF rn_if
	
	!Calcolo bessel
	rn_do: DO n=2,nstop
		v_eta(n)=v_eta(n-1)/v_eta(n)
	END DO rn_do

END SELECT nstop_case

END SUBROUTINE eta_d_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!8) SUBROUTINE eta_z_sub(nstop,z,v_eta,error): 
!			   Si calcolano le funzioni eta(n,r)=j(n,r)/r
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE eta_z_sub(nstop,z,v_eta,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Numero di espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: z							!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_eta		!Vettore funzioni eta
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n									!Numero di coeff pn e indice
REAL(dbl) :: pn,nr										!Coefficiente rapporto e indice reale

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0
!Metto v_eta(0)=0 anche se non e' vero,comunque non lo devo usare,e' per
!uniformita' di inizializzazione
v_eta=(0.0D0,0.0D0)

!Controllo che l'argomento sia maggiore di zero
zero1_if: IF (REAL(z,dbl)<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: Re(z) in eta_z_sub e' minore di zero")
				RETURN
END IF zero1_if

!Se z=0 allora metto solo un valore
zero2_if: IF (z==(0.0D0,0.0D0)) THEN
				v_eta(1)=1.0D0/3.0D0
				RETURN
END IF zero2_if


nstop_case: SELECT CASE (nstop)

CASE (1)

	!Calcolo il mio valore iniziale
	in1_if: IF (ABS(z)>0.1) THEN
		v_eta(1)=SIN(z)/(z**3)-COS(z)/(z**2)
	ELSE
		
		!Ciclo do per l'espansione in serie
		v_eta(1)=(0.0D0,0.0D0)
		in1_do: DO n=4,0,-1
			
			nr=REAL(n,dbl)
			v_eta(1)=v_eta(1)+((-1.0D0)**n)*2.0D0*(nr+1.0D0)*(z**(2*n))/fatt(2.0D0*nr+3.0D0)
		
		END DO in1_do

	END IF in1_if

CASE DEFAULT

	!Calcolo il mio valore iniziale
	in2_if: IF (ABS(z)>0.1) THEN
	
		v_eta(1)=SIN(z)/(z**3)-COS(z)/(z**2)
		
	ELSE
		
		!Ciclo do per l'espansione in serie
		v_eta(1)=(0.0D0,0.0D0)
		in2_do: DO n=0,4
			
			nr=REAL(n,dbl)
			v_eta(1)=v_eta(1)+((-1.0D0)**n)*2.0D0*(nr+1.0D0)*(z**(2*n))/fatt(2.0D0*nr+3.0D0)
		
		END DO in2_do

	END IF in2_if
	
	!Calcolo i rapporti	
	CALL rn_z_sub(2,nstop,z,v_eta(2:nstop),error)

	!If su errore di rn
	rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore in rn_z_sub chiamata da eta_z_sub")
				RETURN
	END IF rn_if
	
	!Calcolo bessel
	rn_do: DO n=2,nstop
		v_eta(n)=v_eta(n-1)/v_eta(n)
	END DO rn_do

END SELECT nstop_case

END SUBROUTINE eta_z_sub




!********************************************************************************
!********************************************************************************
!********************************************************************************
!9) SUBROUTINE rn_d_sub(low_b,up_b,x,v_rn,error): Si calcolano i coeff 
!rn=psi(n-1,x)/psi(n,x), tali coeff servono e per il calcolo delle funzioni di
!riccati-bessel di prima specie, e per il calcolo di a(n) e b(n). Qui il calcolo
!e' per argomento reale
!********************************************************************************
!********************************************************************************
!********************************************************************************
SUBROUTINE rn_d_sub(low_b,up_b,x,v_rn,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: low_b,up_b						!Lower and upper bound array
REAL(dbl), INTENT(IN) :: x								!Argomento reale
REAL(dbl), DIMENSION(low_b:up_b), INTENT(OUT) :: v_rn	!Vettore rn
INTEGER(lo), INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n,inc,digit=15			!Indice di partenza, indice,incremento,cifre significative
REAL(dbl) :: acc								!Accuracy	
REAL(dbl) :: rn,nr								!Pn e indice complesso
COMPLEX(dbl) :: z								!x complessificato


!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (x<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento x rn_d_sub e' errato")
				RETURN
END IF zero_if

!Controllo sui bounds
bounds_if: IF (up_b<low_b) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: in rn_d_sub up_b<low_b")
				RETURN
END IF bounds_if

!Controllo sui up_b
low_if: IF (low_b<1) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: in rn_d_sub il lower bound dev'essere > di 1")
				RETURN
END IF low_if


!Calcolo dell'indice di partenza per il coeff. r(n)
inc=NINT(MAX(4.0D0*(ABS(x)**(1.0D0/3.0D0)),5.0D0),lo)
nmax=up_b+inc
z=CMPLX(x,KIND=dbl)

acc_do: DO 
	acc=l_sub(nmax,z)-l_sub(up_b,z)
	IF (acc>digit) EXIT
	nmax=nmax+inc
END DO acc_do

!Calcolo primo valore dell'iterazione r(nmax+1)
rn=(2.0D0*REAL(nmax,dbl)+2.0D0)/x

!Ciclo per il calcolo di pn e per il suo storage nel vettore v_besj
!fino alla sua successiva sostituzione
rn_do: DO n=nmax,low_b,-1

	!Calcolo pn(n)
	nr=REAL(n,dbl)
	rn=(2.0D0*nr+1.0D0)/x - 1.0D0/rn

	!Se n>nmax allora torno all'inizio
	IF (n>up_b) CYCLE rn_do
	
	v_rn(n)=rn
	
END DO rn_do



END SUBROUTINE rn_d_sub



!********************************************************************************
!********************************************************************************
!********************************************************************************
!10) SUBROUTINE rn_z_sub(low_b,up_b,x,v_rn,error): Si calcolano i coeff 
!rn=psi(n-1,x)/psi(n,x), tali coeff servono e per il calcolo delle funzioni di
!riccati-bessel di prima specie, e per il calcolo di a(n) e b(n). Qui il calcolo
!e' per argomento complesso
!********************************************************************************
!********************************************************************************
!********************************************************************************
SUBROUTINE rn_z_sub(low_b,up_b,x,v_rn,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: low_b,up_b						!Lower and upper bound array
COMPLEX(dbl), INTENT(IN) :: x								!Argomento reale
COMPLEX(dbl), DIMENSION(low_b:up_b), INTENT(OUT) :: v_rn	!Vettore rn
INTEGER(lo), INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n,inc,digit=15			!Indice di partenza, indice,incremento,cifre significative
REAL(dbl) :: acc								!Accuracy	
COMPLEX(dbl) :: rn,nc							!Pn e indice complesso

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(x,dbl)<0.0D0) THEN
!zero_if: IF ((REAL(x,dbl)<0.0D0) .OR. (-REAL((0.0D0,1.0D0)*x,dbl)<0.0D0)) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento x rn_z_sub e' errato")
				RETURN
END IF zero_if

!Controllo sui bounds
bounds_if: IF (up_b<low_b) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: in rn_z_sub up_b<low_b")
				RETURN
END IF bounds_if

!Controllo sui up_b
low_if: IF (low_b<1) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: in rn_z_sub il lower bound dev'essere > di 1")
				RETURN
END IF low_if


!Calcolo dell'indice di partenza per il coeff. r(n)
inc=NINT(MAX(4.0D0*(ABS(x)**(1.0D0/3.0D0)),5.0D0),lo)
nmax=up_b+inc

acc_do: DO 
	acc=l_sub(nmax,x)-l_sub(up_b,x)
	IF (acc>digit) EXIT
	nmax=nmax+inc
END DO acc_do

!Calcolo primo valore dell'iterazione r(nmax+1)
rn=((2.0D0,0.0D0)*CMPLX(nmax,KIND=dbl)+(2.0D0,0.0D0))/x

!Ciclo per il calcolo di pn e per il suo storage nel vettore v_besj
!fino alla sua successiva sostituzione
rn_do: DO n=nmax,low_b,-1

	!Calcolo pn(n)
	nc=CMPLX(n,KIND=dbl)
	rn=((2.0D0,0.0D0)*nc+(1.0D0,0.0D0))/x - (1.0D0,0.0D0)/rn

	!Se n>nmax allora torno all'inizio
	IF (n>up_b) CYCLE 
	
	v_rn(n)=rn
	
END DO rn_do

END SUBROUTINE rn_z_sub



!********************************************************************************
!********************************************************************************
!********************************************************************************
!9) SUBROUTINE qn_d_sub(low_b,up_b,r,v_rn,error): Si calcolano i coeff 
!qn(r)=r*rn(r), tali coeff servono al calcolo del campo interno alla sfera senza
!avere divergenze in zero
!********************************************************************************
!********************************************************************************
!********************************************************************************
SUBROUTINE qn_d_sub(low_b,up_b,r,v_qn,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: low_b,up_b						!Lower and upper bound array
REAL(dbl), INTENT(IN) :: r									!Argomento reale
REAL(dbl), DIMENSION(low_b:up_b), INTENT(OUT) :: v_qn		!Vettore rn
INTEGER(lo), INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n,inc,digit=15			!Indice di partenza, indice,incremento,cifre significative
REAL(dbl) :: acc								!Accuracy	
REAL(dbl) :: qn,nr								!Indice reale
COMPLEX(dbl) :: z								!Argomento complessificato

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (r<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento r qn_d_sub e' errato")
				RETURN
END IF zero_if

!Controllo sui bounds
bounds_if: IF (up_b<low_b) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: in qn_d_sub up_b<low_b")
				RETURN
END IF bounds_if

!Controllo sui up_b
low_if: IF (low_b<1) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: in qn_d_sub il lower bound dev'essere > di 1")
				RETURN
END IF low_if

!Caso particolare per z=0.0
zero1_if: IF (r==0.0D0) THEN

	!Ciclo per il vettore per z=0.0
	zero_do: DO n=low_b,up_b
	
		nr=REAL(n,dbl)
		v_qn(n)=2.0D0*nr+1.0D0
		
	END DO zero_do
	RETURN

END IF zero1_if


!Calcolo dell'indice di partenza per il coeff. r(n)
inc=NINT(MAX(4.0D0*(ABS(r)**(1.0D0/3.0D0)),5.0D0),lo)
nmax=up_b+inc
z=CMPLX(r,KIND=dbl)

acc_do: DO 
	acc=l_sub(nmax,z)-l_sub(up_b,z)
	IF (acc>digit) EXIT
	nmax=nmax+inc
END DO acc_do

!Calcolo primo valore dell'iterazione r(nmax+1)
qn=(2.0D0*REAL(nmax,dbl)+2.0D0)

!Ciclo per il calcolo di pn e per il suo storage nel vettore v_besj
!fino alla sua successiva sostituzione
qn_do: DO n=nmax,low_b,-1

	!Calcolo pn(n)
	nr=REAL(n,dbl)
	qn=(2.0D0*nr+1.0D0) - (r**2)/qn

	!Se n>nmax allora torno all'inizio
	IF (n>up_b) CYCLE qn_do
	
	v_qn(n)=qn
	
END DO qn_do

END SUBROUTINE qn_d_sub






!********************************************************************************
!********************************************************************************
!********************************************************************************
!9bis) SUBROUTINE qn_z_sub(low_b,up_b,z,v_rn,error): Si calcolano i coeff 
!qn(z)=z*rn(z), tali coeff servono al calcolo del campo interno alla sfera senza
!avere divergenze in zero
!********************************************************************************
!********************************************************************************
!********************************************************************************
SUBROUTINE qn_z_sub(low_b,up_b,z,v_qn,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: low_b,up_b						!Lower and upper bound array
COMPLEX(dbl), INTENT(IN) :: z								!Argomento complesso
COMPLEX(dbl), DIMENSION(low_b:up_b), INTENT(OUT) :: v_qn	!Vettore rn
INTEGER(lo), INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: nmax,n,inc,digit=15			!Indice di partenza, indice,incremento,cifre significative
REAL(dbl) :: acc								!Accuracy	
COMPLEX(dbl) :: qn,nr							!Indice complesso

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(z,dbl)<0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento z qn_z_sub e' errato")
				RETURN
END IF zero_if

!Controllo sui bounds
bounds_if: IF (up_b<low_b) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: in qn_z_sub up_b<low_b")
				RETURN
END IF bounds_if

!Controllo sui up_b
low_if: IF (low_b<1) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: in qn_z_sub il lower bound dev'essere > di 1")
				RETURN
END IF low_if

!Caso particolare per z=0.0
zero1_if: IF (z==(0.0D0,0.0D0)) THEN

	!Ciclo per il vettore per z=0.0
	zero_do: DO n=low_b,up_b
	
		nr=REAL(n,dbl)
		v_qn(n)=2.0D0*nr+1.0D0
		
	END DO zero_do
	RETURN

END IF zero1_if


!Calcolo dell'indice di partenza per il coeff. r(n)
inc=NINT(MAX(4.0D0*(ABS(z)**(1.0D0/3.0D0)),5.0D0),lo)
nmax=up_b+inc

acc_do: DO 
	acc=l_sub(nmax,z)-l_sub(up_b,z)
	IF (acc>digit) EXIT
	nmax=nmax+inc
END DO acc_do

!Calcolo primo valore dell'iterazione r(nmax+1)
qn=(2.0D0*CMPLX(nmax,KIND=dbl)+2.0D0)

!Ciclo per il calcolo di pn e per il suo storage nel vettore v_besj
!fino alla sua successiva sostituzione
qn_do: DO n=nmax,low_b,-1

	!Calcolo pn(n)
	nr=REAL(n,dbl)
	qn=(2.0D0*nr+1.0D0) - (z**2)/qn

	!Se n>nmax allora torno all'inizio
	IF (n>up_b) CYCLE qn_do
	
	v_qn(n)=qn
	
END DO qn_do

END SUBROUTINE qn_z_sub







!******************************************************************************
!******************************************************************************
!******************************************************************************
!10) SUBROUTINE psi_d_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di secondo ordine: chi(0)...chi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE psi_d_sub(nstop,x,v_psi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
REAL(dbl), INTENT(IN) :: x								!Argomento reale
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_psi		!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n										!Indice												!Pn e indice reale
REAL(dbl) :: nr											!Pn e indice complesso

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(x,dbl)<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento x psi_d_sub e' errato")
				RETURN
END IF zero_if

!Calcolo i valori di partenza
v_psi(0)=SIN(x)

!Chiamo la mia subroutine r(n) reale
CALL rn_d_sub(1,nstop,x,v_psi(1:nstop),error)

rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in rn_d_sub chiamata da psi_d_sub")
				RETURN
END IF rn_if

!Ricorsione
psi_do: DO n=1,nstop
	v_psi(n)=v_psi(n-1)/v_psi(n)
END DO psi_do

END SUBROUTINE psi_d_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!11) SUBROUTINE psi_z_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di primo tipo: psi(0)...psi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE psi_z_sub(nstop,z,v_psi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
COMPLEX(dbl), INTENT(IN) :: z								!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_psi		!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n										!Indice												!Pn e indice reale
COMPLEX(dbl) :: nr											!Pn e indice complesso

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(z,dbl)<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento x psi_d_sub e' errato")
				RETURN
END IF zero_if

!Controllo che la parte immaginaria dell'argomento non sia troppo grande
max_if: IF (ABS(AIMAG(z))>700.0D0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: l'argomento z psi_d_sub e' > 700,overflow")
				RETURN
END IF max_if

!Calcolo i valori di partenza
v_psi(0)=SIN(z)

!Chiamo la mia subroutine r(n) reale
CALL rn_z_sub(1,nstop,z,v_psi(1:nstop),error)

rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in rn_z_sub chiamata da psi_z_sub")
				RETURN
END IF rn_if

!Ricorsione
psi_do: DO n=1,nstop
	v_psi(n)=v_psi(n-1)/v_psi(n)
END DO psi_do

END SUBROUTINE psi_z_sub





!******************************************************************************
!******************************************************************************
!******************************************************************************
!12) SUBROUTINE SUBROUTINE chi_d_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di secondo ordine: chi(0)...chi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE chi_d_sub(nstop,x,v_chi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
REAL(dbl), INTENT(IN) :: x								!Argomento reale
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_chi		!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n										!Indice												!Pn e indice reale
REAL(dbl) :: nr											!Pn e indice complesso

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(x,dbl)<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento x chi_d_sub e' errato")
				RETURN
END IF zero_if

!Calcolo i valori di partenza
v_chi(0)=COS(x)
v_chi(1)=COS(x)/x+SIN(x)

!Ricorsione
chi_do: DO n=2,nstop

	nr=REAL(n,dbl)
	v_chi(n)=(2.0D0*nr-1.0D0)*v_chi(n-1)/x-v_chi(n-2)

END DO chi_do

END SUBROUTINE chi_d_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!12bis) SUBROUTINE SUBROUTINE chi_z_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di secondo ordine: chi(0)...chi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE chi_z_sub(nstop,z,v_chi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
COMPLEX(dbl), INTENT(IN) :: z								!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_chi		!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n										!Indice												!Pn e indice reale
REAL(dbl) :: nr											!Pn e indice complesso

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Calcolo i valori di partenza
v_chi(0)=COS(z)
v_chi(1)=COS(z)/z+SIN(z)

!Ricorsione
chi_do: DO n=2,nstop

	nr=REAL(n,dbl)
	v_chi(n)=(2.0D0*nr-1.0D0)*v_chi(n-1)/z-v_chi(n-2)

END DO chi_do

END SUBROUTINE chi_z_sub


!******************************************************************************
!******************************************************************************
!******************************************************************************
!13) SUBROUTINE csi_d_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di terzo ordine: csi(0)...csi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE csi_d_sub(nstop,x,v_csi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
REAL(dbl), INTENT(IN) :: x								!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_csi	!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n										!Indice
REAL(dbl) :: nr											!Pn e indice complesso
REAL(dbl) ,DIMENSION(0:nstop) :: v_psi,v_chi			!Vettori psi(0,x)...psi(nstop,x) e chi(0,x)...chi(nstop,x)
!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Controllo che l'argomento sia maggiore di zero
zero_if: IF (REAL(x,dbl)<=0.0D0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: l'argomento x csi_d_sub e' errato")
				RETURN
END IF zero_if

!Chiamo la subroutine per psi
CALL psi_d_sub(nstop,x,v_psi,error)

psi_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in psi_d_sub chiamata da csi_d_sub")
				RETURN
END IF psi_if

!Chiamo la subroutine per chii
CALL chi_d_sub(nstop,x,v_chi,error)

chi_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in chi_d_sub chiamata da csi_d_sub")
				RETURN
END IF chi_if

!Calcolo csi
v_csi=CMPLX(v_psi,-v_chi,KIND=dbl)


END SUBROUTINE csi_d_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!13bis) SUBROUTINE csi_z_sub(nstop,z,v_csi,error): Calcolo le funzioni 
!di riccati bessel di terzo ordine: csi(0)...csi(nstop). L'argomento delle
!funzioni e' complesso.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE csi_z_sub(nstop,z,v_csi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop					!Ordine funzione calcolate
COMPLEX(dbl), INTENT(IN) :: z						!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_csi	!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error					!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n							!Indice
REAL(dbl) :: nr								!Pn e indice complesso
COMPLEX(dbl) ,DIMENSION(0:nstop) :: v_psi,v_chi			!Vettori psi(0,x)...psi(nstop,x) e chi(0,x)...chi(nstop,x)
!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Chiamo la subroutine per psi
CALL psi_z_sub(nstop,z,v_psi,error)

psi_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in psi_z_sub chiamata da csi_z_sub")
				RETURN
END IF psi_if

!Chiamo la subroutine per chii
CALL chi_z_sub(nstop,z,v_chi,error)

chi_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in chi_z_sub chiamata da csi_z_sub")
				RETURN
END IF chi_if

!Calcolo csi
v_csi=v_psi-CMPLX(0.0D0,1.0D0,KIND=dbl)*v_chi

END SUBROUTINE csi_z_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!18) SUBROUTINE chi_z_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di secondo ordine: chi(0)...chi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE dlog_d_sub(nstop,x,v_dlog,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
REAL(dbl), INTENT(IN) :: x							!Argomento reale
REAL(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_dlog		!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n							!Indice
REAL(dbl) :: nr,dlog_aus						!Indice reale e calcolo ausiliario

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Calcolo i valori di partenza
dlog_aus=0.0D0

!Ricorsione ausiliaria
aus_do: DO n=nstop+15,nstop+2

	nr=REAL(n,dbl)
	dlog_aus=nr/x-1/(dlog_aus+nr/x)		!D[n-1]=n/x-1/(D[n]+n/x)

END DO aus_do

!Calcolo a parte il primo valore vero
nr=REAL(nstop+1,dbl)
v_dlog(nstop)=nr/x-1/(dlog_aus+nr/x)

!Ricorsione
dlog_do: DO n=nstop,1

	nr=REAL(n,dbl)
	v_dlog(n-1)=nr/x-1/(v_dlog(n)-nr/x)		!D[n-1]=n/x-1/(D[n]+n/x)

END DO dlog_do

END SUBROUTINE dlog_d_sub




!******************************************************************************
!******************************************************************************
!******************************************************************************
!18bis) SUBROUTINE chi_z_sub(nstop,x,v_chi,error): Calcolo le funzioni 
!di riccati bessel di secondo ordine: chi(0)...chi(nstop). L'argomento delle
!funzioni e' reale.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE dlog_z_sub(nstop,z,v_dlog,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: nstop						!Ordine funzione calcolate
COMPLEX(dbl), INTENT(IN) :: z							!Argomento reale
COMPLEX(dbl), DIMENSION(0:nstop), INTENT(OUT) :: v_dlog		!Vettore riccati-bessel2
INTEGER(lo), INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: n							!Indice
COMPLEX(dbl) :: nr,dlog_aus						!Indice reale e calcolo ausiliario

!Inizio subroutine vera e propria

!Inizializzo l'array e l'errore
error=0

!Calcolo i valori di partenza
dlog_aus=(0.0D0,0.0D0)

!Ricorsione ausiliaria
aus_do: DO n=nstop+15,nstop+2

	nr=REAL(n,dbl)
	dlog_aus=nr/z-1/(dlog_aus+nr/z)		!D[n-1]=n/x-1/(D[n]+n/x)

END DO aus_do

!Calcolo a parte il primo valore vero
nr=REAL(nstop+1,dbl)
v_dlog(nstop)=nr/z-1/(dlog_aus+nr/z)

!Ricorsione
dlog_do: DO n=nstop,1

	nr=REAL(n,dbl)
	v_dlog(n-1)=nr/z-1/(v_dlog(n)-nr/z)		!D[n-1]=n/x-1/(D[n]+n/x)

END DO dlog_do

END SUBROUTINE dlog_z_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!14) l_sub(n,x): calcolo il numero di cifre significative perse in una forward
!recursion per le funzioni sferiche di bessel ad argomento complesso
!******************************************************************************
!******************************************************************************
!******************************************************************************
FUNCTION l_sub(n,x)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo) ,INTENT(IN) :: n							!Ordine funzione calcolate
COMPLEX(dbl), INTENT(IN) :: x							!Argomento reale
REAL(dbl) :: l_sub									!Funzione

!Variabili interne
COMPLEX(dbl) :: radc
COMPLEX(dbl) :: nc
REAL(dbl) :: nr

!Inizio funzione vera e propria
nr=REAL(n,dbl)
nc=CMPLX(nr,KIND=dbl)
radc=SQRT((1.0D0,0.0D0)-(x/nr)**2)

!Calcolo la funzione
IF ((REAL(x,dbl)<0.0D0) .OR. (-REAL((0.0D0,1.0D0)*x,dbl)<0.0D0)) THEN
	l_sub=-1.0D0
END IF

l_sub=(ABS(REAL((0.0D0,1.0D0)*x,dbl))-LOG(2.0D0)-nr*REAL(LOG(x/nc)+radc-LOG((1.0D0,0.0D0)+radc),dbl)) &
&	   /LOG(10.0D0)

END FUNCTION l_sub

!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
!15) SUBROUTINE coeff_sp2: si calcolano le matrici dei coefficienti di espansione
!di singola particella m_a e m_b, che verranno utilizzati per risolvere il sistema lineare.
!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
SUBROUTINE coeff_sp2(lambda,ref_index,v_req,m_epseq,nstop,neq,m_a,m_b,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index				! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop  					! Numero Espansioni multipolari	
INTEGER(lo), INTENT(IN) :: neq						! Numero sfere eq

COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_a,m_b	! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error		       			! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x						! Vettore Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc					! Vettore complesso funct diel
COMPLEX(dbl), DIMENSION(neq) :: v_m,v_mx				! Vettore ref index normalizzato e vettore mx
REAL(dbl), DIMENSION(0:nstop) :: v_psi,v_chi			! Matrici funzioni bessel
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi				! Matrici funzioni bessel
COMPLEX(dbl), DIMENSION(1:nstop) :: v_rn				! Matrici funzioni bessel
INTEGER(lo) :: n,j									! indici
COMPLEX(dbl) :: nc										! Indice complessificato



! Inizio della procedura vera e propria
		
!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
			v_epsc(j)=CMPLX(m_epseq(j,1),m_epseq(j,2),dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x=(2*pi_d*ref_index*v_req)/lambda
!WRITE(*,*) "x: ", v_x

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)
v_mx=v_m*CMPLX(v_x,KIND=dbl)
!WRITE(*,*) "mx: ", v_mx

m_a=0.0D0
m_b=0.0D0
!Loop per calcolare le matrici coefficiente a e b 
m_ab_loop_out: DO j=1,neq

		!Calcolo riccati bessel di prima specie
		CALL psi_d_sub(nstop,v_x(j),v_psi,error)
		
		psi_if: IF (error/=0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_sp2")
				RETURN
		END IF psi_if	
		
		CALL chi_d_sub(nstop,v_x(j),v_chi,error)
		
		chi_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in chi_d_sub chiamata da coeff_sp2")
				RETURN
		END IF chi_if
		
		
		v_csi=CONJG(CMPLX(v_psi,v_chi,KIND=dbl))
		
		
		CALL rn_z_sub(1,nstop,v_mx(j),v_rn,error)

		rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in rn_z_sub chiamata da coeff_sp2")
				RETURN
		END IF rn_if


	m_ab_loop_in: DO n=1,nstop
	
		nc=CMPLX(n,KIND=dbl)
	
		m_a(n,j)=( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_psi(n) - v_psi(n-1) ) / &
			   & ( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_csi(n) - v_csi(n-1) )
		
		m_b(n,j)=( v_rn(n)*v_m(j)*v_psi(n)-v_psi(n-1) ) / &
			   & ( v_rn(n)*v_m(j)*v_csi(n)-v_csi(n-1) )
					
	END DO m_ab_loop_in 

END DO m_ab_loop_out
				
END SUBROUTINE coeff_sp2


!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
!15bis) SUBROUTINE coeff_sp2_dip: si calcolano le matrici dei coefficienti di espansione
!di singola particella m_a e m_b, che verranno utilizzati per risolvere il sistema lineare,
! ma qui i conti li faccio per il dipolo
!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
SUBROUTINE coeff_sp2_dip(lambda,ref_index,v_req,m_epseq,nstop,neq,m_a,m_b,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index				! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop  					! Numero Espansioni multipolari	
INTEGER(lo), INTENT(IN) :: neq						! Numero sfere eq

COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_a,m_b	! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error		       			! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x						! Vettore Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc					! Vettore complesso funct diel
COMPLEX(dbl), DIMENSION(neq) :: v_m,v_mx				! Vettore ref index normalizzato e vettore mx
REAL(dbl), DIMENSION(0:nstop) :: v_psi,v_chi			! Matrici funzioni bessel
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi				! Matrici funzioni bessel
COMPLEX(dbl), DIMENSION(1:nstop) :: v_rn				! Matrici funzioni bessel
INTEGER(lo) :: n,j									! indici
COMPLEX(dbl) :: nc										! Indice complessificato



! Inizio della procedura vera e propria
		
!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
			v_epsc(j)=CMPLX(m_epseq(j,1),m_epseq(j,2),KIND=dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x=(2*pi_d*ref_index*v_req)/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)
v_mx=v_m*CMPLX(v_x,KIND=dbl)

m_a=0.0D0
m_b=0.0D0
!Loop per calcolare le matrici coefficiente a e b 
m_ab_loop_out: DO j=1,neq

		!Calcolo riccati bessel di prima specie
		CALL psi_d_sub(nstop,v_x(j),v_psi,error)
		
		psi_if: IF (error/=0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_sp2")
				RETURN
		END IF psi_if	
		
		CALL chi_d_sub(nstop,v_x(j),v_chi,error)
		
		chi_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in chi_d_sub chiamata da coeff_sp2")
				RETURN
		END IF chi_if
		
		
		v_csi=CONJG(CMPLX(v_psi,v_chi,KIND=dbl))
		
		
		CALL rn_z_sub(1,nstop,v_mx(j),v_rn,error)

		rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in rn_z_sub chiamata da coeff_sp2")
				RETURN
		END IF rn_if

	dip_if: IF (j==neq) THEN

		m_ab_loop_in_dip: DO n=1,1
			nc=CMPLX(n,KIND=dbl)

			m_a(n,j)=( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_psi(n) - v_psi(n-1) ) / &
				   & ( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_csi(n) - v_csi(n-1) )

			m_b(n,j)=( v_rn(n)*v_m(j)*v_psi(n)-v_psi(n-1) ) / &
				   & ( v_rn(n)*v_m(j)*v_csi(n)-v_csi(n-1) )
		END DO m_ab_loop_in_dip


	ELSE

		m_ab_loop_in: DO n=1,nstop
			nc=CMPLX(n,KIND=dbl)

			m_a(n,j)=( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_psi(n) - v_psi(n-1) ) / &
				   & ( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_csi(n) - v_csi(n-1) )

			m_b(n,j)=( v_rn(n)*v_m(j)*v_psi(n)-v_psi(n-1) ) / &
				   & ( v_rn(n)*v_m(j)*v_csi(n)-v_csi(n-1) )
		END DO m_ab_loop_in

	END IF dip_if

END DO m_ab_loop_out
				
END SUBROUTINE coeff_sp2_dip



!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
!15) SUBROUTINE coeff_sp3: si calcolano le matrici dei coefficienti di espansione
!di singola particella m_a, m_b,m_c e m_d che verranno utilizzati per risolvere il sistema lineare.
!******************************************************************************************
!******************************************************************************************
!******************************************************************************************
SUBROUTINE coeff_sp3(lambda,ref_index,v_req,m_epseq,nstop,neq,m_a,m_b,m_c,m_d,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index					! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req				! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq				! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop  					! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq						! Numero sfere eq

COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_a,m_b,m_c,m_d	! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error		       			! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x					! Vettore Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc				! Vettore complesso funct diel
COMPLEX(dbl), DIMENSION(neq) :: v_m,v_mx				! Vettore ref index normalizzato e vettore mx
REAL(dbl), DIMENSION(0:nstop) :: v_psi,v_chi			! Matrici funzioni bessel
REAL(dbl), DIMENSION(1:nstop) :: v_psi1				! Matrici funzioni bessel
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi,v_psiz		! Matrici funzioni bessel
COMPLEX(dbl), DIMENSION(1:nstop) :: v_rn,v_psi1z,v_csi1	! Matrici funzioni bessel
INTEGER(lo) :: n,j							! indici
REAL(dbl) :: nr								! Indice complessificato
COMPLEX(dbl) :: nc							! Indice complessificato



! Inizio della procedura vera e propria
		
!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
			v_epsc(j)=CMPLX(m_epseq(j,1),m_epseq(j,2),KIND=dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x=(2*pi_d*ref_index*v_req)/lambda


!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)
v_mx=v_m*CMPLX(v_x,KIND=dbl)
!WRITE(*,*) "x: ", v_x
!WRITE(*,*)
!WRITE(*,*) "m: ", v_m
!WRITE(*,*)
!WRITE(*,*) "mx: ", v_mx

m_a=0.0D0
m_b=0.0D0
m_c=0.0D0
m_d=0.0D0
!Loop per calcolare le matrici coefficiente a e b 
m_abcd_loop_out: DO j=1,neq

		!Calcolo riccati bessel di prima specie
		CALL psi_d_sub(nstop,v_x(j),v_psi,error)
		
		psi_if: IF (error/=0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_sp3")
				RETURN
		END IF psi_if	


		!Calcolo riccati bessel di prima specie Complessa
		CALL psi_z_sub(nstop,v_mx(j),v_psiz,error)

		psiz_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in psi_z_sub chiamata da coeff_sp3")
				RETURN
		END IF psiz_if

		CALL chi_d_sub(nstop,v_x(j),v_chi,error)

		chi_if: IF (error/=0) THEN
				error=1
				WRITE(*,25)
				25 FORMAT ("Errore: errore in chi_d_sub chiamata da coeff_sp3")
				RETURN
		END IF chi_if
		
		
		v_csi=CONJG(CMPLX(v_psi,v_chi,KIND=dbl))
		
		
		CALL rn_z_sub(1,nstop,v_mx(j),v_rn,error)

		rn_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in rn_z_sub chiamata da coeff_sp2")
				RETURN
		END IF rn_if


		!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
		der_do: DO n=1,nstop

			nr=REAL(n,dbl)
			nc=CMPLX(n,KIND=dbl)
			!Derivate di psi(x)
			v_psi1(n)=v_psi(n-1)-nr*v_psi(n)/v_x(j)

			!Derivate di psi(x)
			v_psi1z(n)=v_psiz(n-1)-nc*v_psiz(n)/v_mx(j)

			!Derivate di csi(x)
			v_csi1(n)=v_csi(n-1)-nc*v_csi(n)/v_x(j)

		END DO der_do


	m_abcd_loop_in: DO n=1,nstop
	
		nc=CMPLX(n,KIND=dbl)

		m_a(n,j)=( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_psi(n) - v_psi(n-1) ) / &
			   & ( (v_rn(n)/v_m(j) + nc*((1.0D0,0.0D0)-(1.0D0,0.0D0)/(v_m(j)**2))/v_x(j))*v_csi(n) - v_csi(n-1) )

		m_b(n,j)=( v_rn(n)*v_m(j)*v_psi(n)-v_psi(n-1) ) / &
			   & ( v_rn(n)*v_m(j)*v_csi(n)-v_csi(n-1) )


		m_d(n,j)=( v_m(j)*CMPLX(0.0D0,1.0D0,KIND=dbl) ) / ( v_m(j)*v_psiz(n)*v_csi1(n) - v_csi(n)*v_psi1z(n) )

		m_c(n,j)=( v_m(j)*CMPLX(0.0D0,1.0D0,KIND=dbl) ) / ( v_psiz(n)*v_csi1(n) - v_m(j)*v_csi(n)*v_psi1z(n) )

	END DO m_abcd_loop_in 

END DO m_abcd_loop_out
				
END SUBROUTINE coeff_sp3




!******************************************************************************
!******************************************************************************
!******************************************************************************
!16) SUBROUTINE coeff_sp2: si calcolano le matrici dei coefficienti di espansione
!						di singola particella m_a e m_b, che verranno utilizza-
!						-ti per risolvere il sistema lineare.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE coeff_ad_ca(lambda,ref_index,v_req,m_epseq,nstop,neq,m_da,m_cb,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index				! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop						! Numero Espansioni multipolari	
INTEGER(lo), INTENT(IN) :: neq						! Numero sfere eq

COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_da,m_cb	! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error					! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x						! Vettore Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc					! Vettore complesso funct diel
COMPLEX(dbl), DIMENSION(neq) :: v_m,v_mx				! Vettore ref index normalizzato e vettore mx
REAL(dbl), DIMENSION(0:nstop) :: v_psi_d				! Vettore Psi reale
REAL(dbl), DIMENSION(1:nstop) :: v_rn_d				! Vettore rn reale
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z				! Vettore Psi complesso
COMPLEX(dbl), DIMENSION(1:nstop) :: v_rn_z				! Vettore rn complesso
INTEGER(lo) :: n,j									! Indici
REAL(dbl) :: nr											! Indice complessificato



! Inizio della procedura vera e propria		

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
			v_epsc(j)=CMPLX(m_epseq(j,1),m_epseq(j,2),KIND=dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x=(2*pi_d*ref_index*v_req)/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)
v_mx=v_m*CMPLX(v_x,KIND=dbl)


!Loop per calcolare le matrici coefficiente a e b 
m_ab_loop_out: DO j=1,neq

		!Calcolo riccati bessel di prima specie reale
		CALL psi_d_sub(nstop,v_x(j),v_psi_d,error)

		psi_d_if: IF (error/=0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_da_cb")
				RETURN
		END IF psi_d_if	

		!Calcolo riccati bessel di prima specie complesso
		CALL psi_z_sub(nstop,v_mx(j),v_psi_z,error)
		
		psi_z_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in psi_z_sub chiamata da coeff_da_cb")
				RETURN
		END IF psi_z_if
		
		!Calcolo derivata logaritmica reale
		CALL rn_d_sub(1,nstop,v_x(j),v_rn_d,error)

		rn_d_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in rn_d_sub chiamata da coeff_da_cb")
				RETURN
		END IF rn_d_if

		!Ciclo per il calcolo effettivo della derivata reale
		dd_do: DO n=1,nstop
			nr=REAL(n,dbl)
			v_rn_d(n)=v_rn_d(n)-nr/v_x(j)
		END DO dd_do

		!Calcolo derivata logaritmica complessa
		CALL rn_z_sub(1,nstop,v_mx(j),v_rn_z,error)

		rn_z_if: IF (error/=0) THEN
				error=1
				WRITE(*,40)
				40 FORMAT ("Errore: errore in rn_z_sub chiamata da coeff_da_cb")
				RETURN
		END IF rn_z_if

		!Ciclo per il calcolo effettivo della derivata complessa
		dz_do: DO n=1,nstop
			nr=REAL(n,dbl)
			v_rn_z(n)=v_rn_z(n)-nr/v_mx(j)
		END DO dz_do

		DO n=1,nstop

		!Calcolo di m_da e m_cb
		m_da(n,j)=(0.0D0,1.0D0)*v_m(j) / &
			    & (v_psi_d(n)*v_psi_z(n)*(v_m(j)*v_rn_d(n)-v_rn_z(n)))

		m_cb(n,j)=(0.0D0,1.0D0)*v_m(j) / &
			    & (v_psi_d(n)*v_psi_z(n)*(v_rn_d(n)-v_m(j)*v_rn_z(n)))

!		WRITE(*,*) "da",n,j,m_da(n,j)
!		WRITE(*,*) "cb",n,j,m_cb(n,j)

		END DO

END DO m_ab_loop_out

END SUBROUTINE coeff_ad_ca


!******************************************************************************
!******************************************************************************
!******************************************************************************
!16) SUBROUTINE coeff_ad_ca_dip: calcolo i coefficienti per poi calcolare 
! i coefficienti di espansione per il campo interno.
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE coeff_ad_ca_dip(lambda,ref_index,v_req,m_epseq,nstop,neq,m_da,m_cb,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index				! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
INTEGER(lo), INTENT(IN) :: nstop						! Numero Espansioni multipolari	
INTEGER(lo), INTENT(IN) :: neq						! Numero sfere eq

COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_da,m_cb	! Matrici coeff. singola sfera
INTEGER(lo), INTENT(OUT) :: error					! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x						! Vettore Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc					! Vettore complesso funct diel
COMPLEX(dbl), DIMENSION(neq) :: v_m,v_mx				! Vettore ref index normalizzato e vettore mx
REAL(dbl), DIMENSION(0:nstop) :: v_psi_d				! Vettore Psi reale
REAL(dbl), DIMENSION(1:nstop) :: v_rn_d				! Vettore rn reale
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z				! Vettore Psi complesso
COMPLEX(dbl), DIMENSION(1:nstop) :: v_rn_z				! Vettore rn complesso
INTEGER(lo) :: n,j									! Indici
REAL(dbl) :: nr											! Indice complessificato



! Inizio della procedura vera e propria		

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
			v_epsc(j)=CMPLX(m_epseq(j,1),m_epseq(j,2),KIND=dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
v_x=(2*pi_d*ref_index*v_req)/lambda

!Calcolo l'indice di rifrazione normalizzato e il vettore mx
v_m=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)
v_mx=v_m*CMPLX(v_x,KIND=dbl)

m_da=0.0D0
m_cb=0.0D0
!Loop per calcolare le matrici coefficiente a e b 
m_ab_loop_out: DO j=1,neq

		!Calcolo riccati bessel di prima specie reale
		CALL psi_d_sub(nstop,v_x(j),v_psi_d,error)

		psi_d_if: IF (error/=0) THEN
				error=1
				WRITE(*,10)
				10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_da_cb")
				RETURN
		END IF psi_d_if

		!Calcolo riccati bessel di prima specie complesso
		CALL psi_z_sub(nstop,v_mx(j),v_psi_z,error)

		psi_z_if: IF (error/=0) THEN
				error=1
				WRITE(*,20)
				20 FORMAT ("Errore: errore in psi_z_sub chiamata da coeff_da_cb")
				RETURN
		END IF psi_z_if

		!Calcolo derivata logaritmica reale
		CALL rn_d_sub(1,nstop,v_x(j),v_rn_d,error)

		rn_d_if: IF (error/=0) THEN
				error=1
				WRITE(*,30)
				30 FORMAT ("Errore: errore in rn_d_sub chiamata da coeff_da_cb")
				RETURN
		END IF rn_d_if

		!Ciclo per il calcolo effettivo della derivata reale
		dd_do: DO n=1,nstop
			nr=REAL(n,dbl)
			v_rn_d(n)=v_rn_d(n)-nr/v_x(j)
		END DO dd_do

		!Calcolo derivata logaritmica complessa
		CALL rn_z_sub(1,nstop,v_mx(j),v_rn_z,error)

		rn_z_if: IF (error/=0) THEN
				error=1
				WRITE(*,40)
				40 FORMAT ("Errore: errore in rn_z_sub chiamata da coeff_da_cb")
				RETURN
		END IF rn_z_if

		!Ciclo per il calcolo effettivo della derivata complessa
		dz_do: DO n=1,nstop
			nr=REAL(n,dbl)
			v_rn_z(n)=v_rn_z(n)-nr/v_mx(j)
		END DO dz_do

		dip_if: IF (j==neq) THEN

			DO n=1,1
				!Calcolo di m_da e m_cb
				m_da(n,j)=(0.0D0,1.0D0)*v_m(j) / &
					    & (v_psi_d(n)*v_psi_z(n)*(v_m(j)*v_rn_d(n)-v_rn_z(n)))

				m_cb(n,j)=(0.0D0,1.0D0)*v_m(j) / &
					    & (v_psi_d(n)*v_psi_z(n)*(v_rn_d(n)-v_m(j)*v_rn_z(n)))
			END DO

		ELSE

			DO n=1,nstop
				!Calcolo di m_da e m_cb
				m_da(n,j)=(0.0D0,1.0D0)*v_m(j) / &
					    & (v_psi_d(n)*v_psi_z(n)*(v_m(j)*v_rn_d(n)-v_rn_z(n)))

				m_cb(n,j)=(0.0D0,1.0D0)*v_m(j) / &
					    & (v_psi_d(n)*v_psi_z(n)*(v_rn_d(n)-v_m(j)*v_rn_z(n)))
			END DO

		END IF dip_if



END DO m_ab_loop_out

END SUBROUTINE coeff_ad_ca_dip




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
!******************************************************************************
!******************************************************************************
!1) SUBROUTINE coeff_shell: si calcolano i vettori dei coefficienti di per passare
! tra i diversi coefficienti della shell
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE coeff_shell(lambda,ref_index,v_req,m_epseq,v_p,nstop,neq,matrixside,v_qa,v_qb,v_gamma, &
			   & m_t1,m_u1,m_t2,m_u2,v_dc0,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index				! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p			! Vettore espensioni campo incidente
INTEGER(lo), INTENT(IN) :: nstop					! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside			! Numero sfere eq e lato matrice

COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_qa,v_qb				! Matrici coeff. singola sfera
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_gamma				! Vettore ausiliario RHS
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_t1,m_u1,m_t2,m_u2		! Matrici ausiliarie per sistema lineare
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_dc0				! Starting Values d_mn,c_mn

INTEGER(lo), INTENT(OUT) :: error					! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x										! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc									! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq) :: v_k									! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq) :: v_m									! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx								! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z01,v_psi_z11,v_psi_z12,v_psi_z22			! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12					! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z01,v_psi1_z11,v_psi1_z12,v_psi1_z22		! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z01,v_csi1_z11,v_csi1_z12				! Vettori Csi' complessi
INTEGER(lo) :: m,n,nu,i,j											! Indici
REAL(dbl) :: nr													! Indice complessificato



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
v_m(1:neq)=SQRT(v_epsc)/CMPLX(ref_index)

!Calcolo la matrice mx
mx_out_do: DO i=0,neq
	mx_in_do: DO j=1,neq

		m_mx(i,j)=v_m(i)*CMPLX(v_x(j))

	END DO mx_in_do
END DO mx_out_do

!Loop per calcolare le matrici coefficiente a e b 

!Calcolo riccati bessel di prima specie Complessa
CALL psi_z_sub(nstop,m_mx(0,1),v_psi_z01,error)
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(1,2),v_psi_z12,error)
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
		error=1
		WRITE(*,10)
		10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_shell")
		RETURN
END IF psi_z_if	

!Calcolo riccati bessel di terza specie complesso
CALL csi_z_sub(nstop,m_mx(0,1),v_csi_z01,error)
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(1,2),v_csi_z12,error)





csi_z_if: IF (error/=0) THEN
		error=1
		WRITE(*,20)
		20 FORMAT ("Errore: errore in csi_z_sub chiamata da coeff_shell")
		RETURN
END IF csi_z_if

!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
der_do: DO n=1,nstop

	nr=REAL(n,dbl)
	!Derivate di Psi
	v_psi1_z01(n)=v_psi_z01(n-1)-nr*v_psi_z01(n)/m_mx(0,1)
	v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
	v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)
	v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

	!Derivate di csi
	v_csi1_z01(n)=v_csi_z01(n-1)-nr*v_csi_z01(n)/m_mx(0,1)
	v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
	v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)

END DO der_do


!output di controllo per i vettori di riccatibessel
!WRITE(*,*) "Argomento funzione 01: ", m_mx(0,1)
!WRITE(*,*) "Psi01: ", (v_psi_z01(n), n=1,nstop)
!WRITE(*,*) "Psi1-01: ", (v_psi1_z01(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi01: ", (v_csi_z01(n), n=1,nstop)
!WRITE(*,*) "Csi1-01: ", (v_csi1_z01(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Argomento funzione 11: ", m_mx(1,1)
!WRITE(*,*) "Psi11: ", (v_psi_z11(n), n=1,nstop)
!WRITE(*,*) "Psi-11: ", (v_psi1_z11(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi11: ", (v_csi_z11(n), n=1,nstop)
!WRITE(*,*) "Csi1-11: ", (v_csi1_z11(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Argomento funzione 12: ", m_mx(1,2)
!WRITE(*,*) "Psi12: ", (v_psi_z12(n), n=1,nstop)
!WRITE(*,*) "Psi1-12: ", (v_psi1_z12(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Csi12: ", (v_csi_z12(n), n=1,nstop)
!WRITE(*,*) "Csi1-12: ", (v_csi1_z12(n), n=1,nstop)
!WRITE(*,*)
!WRITE(*,*) "Argomento funzione 22: ", m_mx(2,2)
!WRITE(*,*) "Psi22: ", (v_psi_z22(n), n=1,nstop)
!WRITE(*,*) "Psi1-22: ", (v_psi1_z22(n), n=1,nstop)
!WRITE(*,*)


!Calcolo i fattori Qa e Qb
eq_if: IF (v_epsc(1)==v_epsc(2)) THEN

	v_qa=(0.0D0,0.0D0)
	v_qb=(0.0D0,0.0D0)

ELSE

	q_do: DO n=1,nstop

		!Calcolo di qa e qb
		v_qb(n)=( v_k(1)*v_psi1_z12(n)*v_psi_z22(n) - v_k(2)*v_psi_z12(n)*v_psi1_z22(n) ) / &
			& ( v_k(2)*v_csi_z12(n)*v_psi1_z22(n) - v_k(1)*v_csi1_z12(n)*v_psi_z22(n) )

		v_qa(n)=( v_k(2)*v_psi1_z12(n)*v_psi_z22(n) - v_k(1)*v_psi_z12(n)*v_psi1_z22(n) ) / &
			& ( v_k(1)*v_csi_z12(n)*v_psi1_z22(n) - v_k(2)*v_csi1_z12(n)*v_psi_z22(n) )

	END DO q_do

END IF eq_if


!Calcolo il vettore ausiliaro RHS gamma
v_gamma=v_k(1)*( v_csi1_z01*v_psi_z01(1:nstop) - v_psi1_z01*v_csi_z01(1:nstop) )

!Calcolo le matrici ausiliarie che servono per la matrice dei coefficienti del sistema lineare
n_do: DO n=1,nstop
	nu_do: DO nu=1,nstop

	m_t1(n,nu)=v_k(0)*v_csi1_z01(n)*(v_psi_z11(n)+v_qb(nu)*v_csi_z11(n)) - &
		   & v_k(1)*v_csi_z01(n)*(v_psi1_z11(n)+v_qb(nu)*v_csi1_z11(n))

	m_t2(n,nu)=v_k(1)*v_csi1_z01(n)*(v_psi_z11(n)+v_qb(nu)*v_csi_z11(n)) - &
		   & v_k(0)*v_csi_z01(n)*(v_psi1_z11(n)+v_qb(nu)*v_csi1_z11(n))

	m_u1(n,nu)=v_k(0)*v_csi1_z01(n)*(v_psi_z11(n)+v_qa(nu)*v_csi_z11(n)) - &
		   & v_k(1)*v_csi_z01(n)*(v_psi1_z11(n)+v_qa(nu)*v_csi1_z11(n))

	m_u2(n,nu)=v_k(1)*v_csi1_z01(n)*(v_psi_z11(n)+v_qa(nu)*v_csi_z11(n)) - &
		   & v_k(0)*v_csi_z01(n)*(v_psi1_z11(n)+v_qa(nu)*v_csi1_z11(n))

	END DO nu_do
END DO n_do

!Calcolo i valori per i vettori iniziali,che sono tutto quello che mi serve per una shell concentrica
j=0
dc_out_do: DO n=1,nstop
	dc_in_do: DO m=-n,n

	!Aggiorno l'indice
	j=j+1

	!Riempio il vettore dei valori iniziali
	v_dc0(2*j-1)=v_gamma(n)*v_p(2*j-1)/m_u2(n,n)		!d_mn
	v_dc0(2*j)=v_gamma(n)*v_p(2*j)/m_t1(n,n)		!c_mn

!	WRITE(*,*) "v_d"
!	WRITE(*,*) v_gamma(n)/m_u2(n,n)
!	WRITE(*,*)
!	WRITE(*,*) "Quello che dovrebbe essere i"
!	WRITE(*,*) v_gamma/v_k(1)
!	WRITE(*,*)

	END DO dc_in_do
END DO dc_out_do

END SUBROUTINE coeff_shell





!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE coeff_shell_borghese: Calculation of auxiliary coefficients
! to calculate the coefficent matrix to be solved
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE coeff_shell_borghese(lambda,ref_index,v_req,m_epseq,nstop,neq,v_Rn,v_Vn,v_Zn,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index			! Wavelength and ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vector of equals r
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Dielectri Function Matrix
INTEGER(lo), INTENT(IN) :: nstop				! Multipolar expansions
INTEGER(lo), INTENT(IN) :: neq				! Number of equal spheres

COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Rn,v_Vn,v_Zn	! Auxiliary coefficients
INTEGER(lo), INTENT(OUT) :: error				! Error flag


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x							! Size Parameters Vector
COMPLEX(dbl), DIMENSION(neq) :: v_epsc							! Dielectric Function
COMPLEX(dbl), DIMENSION(0:neq) :: v_k							! Wave vectors
COMPLEX(dbl), DIMENSION(0:neq) :: v_m							! Normalized refractive index
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx						! Riccati bessel functions arguments
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z12,v_psi_z22			! Riccati bessel
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12			! 
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11,v_psi1_z12,v_psi1_z22			! 
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z01,v_csi1_z11,v_csi1_z12			! 
INTEGER(lo) :: m,n,nu,i,j								! Indexes
REAL(dbl) :: nr										! 



!Subroutine begins

!Remapping of the dielectric function
v_epsc_do: DO j=1,neq
			!using neq-j,to match Ngo formalism, (0=host,1=shell,2=core), instead of mine (1=core,2=shell)
			v_epsc(j)=CMPLX(m_epseq(neq+1-j,1),m_epseq(neq+1-j,2))
END DO v_epsc_do


!Size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Wavevectors
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq)=(2*pi_d*SQRT(v_epsc))/lambda

!Normalized ref index
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq)=SQRT(v_epsc)/CMPLX(ref_index)

!Function argument matrix
mx_out_do: DO i=0,neq
	mx_in_do: DO j=1,neq

		m_mx(i,j)=v_m(i)*CMPLX(v_x(j))

	END DO mx_in_do
END DO mx_out_do


!Riccati Bessel Calculation psi
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(1,2),v_psi_z12,error)
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
		error=1
		WRITE(*,10)
		10 FORMAT ("Error: psi_z_sub called by coeff_shell_borghese")
		RETURN
END IF psi_z_if	

!Riccati Bessel Calculation csi
CALL csi_z_sub(nstop,m_mx(0,1),v_csi_z01,error)
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(1,2),v_csi_z12,error)

csi_z_if: IF (error/=0) THEN
		error=1
		WRITE(*,20)
		20 FORMAT ("Error: csi_z_sub called by coeff_shell_borghese")
		RETURN
END IF csi_z_if


!Riccati Bessel derivatives
der_do: DO n=1,nstop

	nr=REAL(n,dbl)
	!Psi
	v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
	v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)
	v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

	!Csi
	v_csi1_z01(n)=v_csi_z01(n-1)-nr*v_csi_z01(n)/m_mx(0,1)
	v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
	v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)

END DO der_do


!Rn
v_Rn(1:2*nstop:2)=( v_k(1)*v_csi_z12(1:nstop)*v_psi1_z22 - v_k(2)*v_csi1_z12*v_psi_z22(1:nstop) ) / &
                      &( v_k(2)*v_psi1_z12*v_psi_z22(1:nstop) - v_k(1)*v_psi_z12(1:nstop)*v_psi1_z22 )
v_Rn(2:2*nstop:2)=( v_k(2)*v_csi_z12(1:nstop)*v_psi1_z22 - v_k(1)*v_csi1_z12*v_psi_z22(1:nstop) ) / &
                      &( v_k(1)*v_psi1_z12*v_psi_z22(1:nstop) - v_k(2)*v_psi_z12(1:nstop)*v_psi1_z22 )

!Tn
v_Vn(1:2*nstop:2)=(0.0D0,1.0D0)*( (v_k(0)/v_k(1))*v_csi_z01(1:nstop)*v_psi1_z11 - v_csi1_z01*v_psi_z11(1:nstop))
v_Vn(2:2*nstop:2)=(0.0D0,1.0D0)*( v_csi_z01(1:nstop)*v_psi1_z11 - (v_k(0)/v_k(1))*v_csi1_z01*v_psi_z11(1:nstop))

!Un
v_Zn(1:2*nstop:2)=(0.0D0,1.0D0)*( v_csi1_z01*v_csi_z11(1:nstop) - (v_k(0)/v_k(1))*v_csi_z01(1:nstop)*v_csi1_z11)
v_Zn(2:2*nstop:2)=(0.0D0,1.0D0)*( (v_k(0)/v_k(1))*v_csi1_z01*v_csi_z11(1:nstop) - v_csi_z01(1:nstop)*v_csi1_z11 )

END SUBROUTINE coeff_shell_borghese






!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE coeff_shell_dip: si calcolano i vettori dei coefficienti di per passare
! tra i diversi coefficienti della shell
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE coeff_shell_dip(lambda,ref_index,v_req,m_epseq,v_p,nstop,neq,matrixside,tflag,v_qa,v_qb,v_gamma, &
			   & m_t1,m_u1,m_t2,m_u2,v_dc0,error) 

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index				! Lunghezza d'onda in questione e ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p			! Vettore espensioni campo incidente
INTEGER(lo), INTENT(IN) :: nstop					! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside,tflag		! Numero sfere eq e lato matrice

COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_qa,v_qb				! Matrici coeff. singola sfera
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_gamma				! Vettore ausiliario RHS
COMPLEX(dbl), DIMENSION(:,:), INTENT(OUT) :: m_t1,m_u1,m_t2,m_u2		! Matrici ausiliarie per sistema lineare
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_dc0				! Starting Values d_mn,c_mn

INTEGER(lo), INTENT(OUT) :: error					! Controllo errore


! Dichiarazione variabili interne
!Metto dappertutto neq-1 per tagliare via la parte col dipolo
REAL(dbl), DIMENSION(neq-1) :: v_x										! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq-1) :: v_epsc									! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_k									! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_m									! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq-1,1:neq-1) :: m_mx								! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z12,v_psi_z22					! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12,v_csi_z22			! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11,v_psi1_z12,v_psi1_z22				! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z01,v_csi1_z11,v_csi1_z12,v_csi1_z22		! Vettori Csi' complessi
INTEGER(lo) :: m,n,nu,i,j											! Indici
REAL(dbl) :: nr													! Indice complessificato



! Inizio della procedura vera e propria

!Costruisco il vettore complesso, neq -> (neq-1)
v_epsc_do: DO j=1,neq-1
			!Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
			!da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
			v_epsc(j)=CMPLX(m_epseq(neq-j,1),m_epseq(neq-j,2),KIND=dbl)
END DO v_epsc_do


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
		10 FORMAT ("Errore: errore in psi_d_sub chiamata da coeff_shell")
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
		20 FORMAT ("Errore: errore in csi_z_sub chiamata da coeff_shell")
		RETURN
END IF csi_z_if

!Calcolo le derivate corrispondenti,nessu problema numerico in vista credo
der_do: DO n=1,nstop

	nr=REAL(n,dbl)
	!Derivate di Psi
	v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
	v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)
	v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

	!Derivate di csi
	v_csi1_z01(n)=v_csi_z01(n-1)-nr*v_csi_z01(n)/m_mx(0,1)
	v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
	v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)
	v_csi1_z22(n)=v_csi_z22(n-1)-nr*v_csi_z22(n)/m_mx(2,2)

END DO der_do



!Calcolo i fattori Qa e Qb
eq_if: IF (v_epsc(1)==v_epsc(2)) THEN

	v_qa=(0.0D0,0.0D0)
	v_qb=(0.0D0,0.0D0)

ELSE

	q_do: DO n=1,nstop

		!Calcolo di qa e qb
		v_qb(n)=( v_k(1)*v_psi1_z11(n)*v_csi_z01(n) - v_k(0)*v_psi_z11(n)*v_csi1_z01(n) ) / &
			& ( v_k(0)*v_csi_z11(n)*v_csi1_z01(n) - v_k(1)*v_csi1_z11(n)*v_csi_z01(n) )

		v_qa(n)=( v_k(0)*v_psi1_z11(n)*v_csi_z01(n) - v_k(1)*v_psi_z11(n)*v_csi1_z01(n) ) / &
			& ( v_k(1)*v_csi_z11(n)*v_csi1_z01(n) - v_k(0)*v_csi1_z11(n)*v_csi_z01(n) )

	END DO q_do

END IF eq_if


!Calcolo il vettore ausiliaro RHS gamma: in questo caso ho il meno davanti, perche' per via delle sostituzioni del dipolo
!funziona cosi'
v_gamma=-v_k(1)*( v_psi1_z22*v_csi_z22(1:nstop) - v_csi1_z22*v_psi_z22(1:nstop) )


!Calcolo le matrici ausiliarie che servono per la matrice dei coefficienti del sistema lineare
!Avro' due casi: tflag=0, caso consueto. tflag=1, caso in cui la traslazione e' piu' grande del raggio.
tflag_if: IF (tflag==0) THEN

	n_do: DO n=1,nstop
		nu_do: DO nu=1,nstop

		m_t1(n,nu)=v_k(2)*v_psi1_z22(n)*(v_psi_z12(n)+v_qb(nu)*v_csi_z12(n)) - &
			   & v_k(1)*v_psi_z22(n)*(v_psi1_z12(n)+v_qb(nu)*v_csi1_z12(n))

		m_t2(n,nu)=v_k(1)*v_psi1_z22(n)*(v_psi_z12(n)+v_qb(nu)*v_csi_z12(n)) - &
			   & v_k(2)*v_psi_z22(n)*(v_psi1_z12(n)+v_qb(nu)*v_csi1_z12(n))

		m_u1(n,nu)=v_k(2)*v_psi1_z22(n)*(v_psi_z12(n)+v_qa(nu)*v_csi_z12(n)) - &
			   & v_k(1)*v_psi_z22(n)*(v_psi1_z12(n)+v_qa(nu)*v_csi1_z12(n))

		m_u2(n,nu)=v_k(1)*v_psi1_z22(n)*(v_psi_z12(n)+v_qa(nu)*v_csi_z12(n)) - &
			   & v_k(2)*v_psi_z22(n)*(v_psi1_z12(n)+v_qa(nu)*v_csi1_z12(n))

		END DO nu_do
	END DO n_do

ELSE IF (tflag==1) THEN

	n_do1: DO n=1,nstop
		nu_do1: DO nu=1,nstop

		m_t1(n,nu)=v_k(2)*v_psi1_z22(n)*v_psi_z12(n)-v_k(1)*v_psi_z22(n)*v_psi1_z12(n)

		m_t2(n,nu)=v_k(1)*v_psi1_z22(n)*v_psi_z12(n)-v_k(2)*v_psi_z22(n)*v_psi1_z12(n)

		m_u1(n,nu)=v_k(2)*v_psi1_z22(n)*v_psi_z12(n)-v_k(1)*v_psi_z22(n)*v_psi1_z12(n)

		m_u2(n,nu)=v_k(1)*v_psi1_z22(n)*v_psi_z12(n)-v_k(2)*v_psi_z22(n)*v_psi1_z12(n)

		END DO nu_do1
	END DO n_do1

END IF tflag_if


!Calcolo i valori per i vettori iniziali,che sono tutto quello che mi serve per una shell concentrica
!Alla fine, formalmente, il tutto e' analogo al caso classico e quello che cambia e' solo il RHS, quindi e' tutto incluso
j=0
dc_out_do: DO n=1,nstop
	dc_in_do: DO m=-n,n

	!Aggiorno l'indice
	j=j+1

	!Riempio il vettore dei valori iniziali
	v_dc0(2*j-1)=v_gamma(n)*v_p(2*j-1)/m_u2(n,n)		!d_mn
	v_dc0(2*j)=v_gamma(n)*v_p(2*j)/m_t1(n,n)		!c_mn

	END DO dc_in_do
END DO dc_out_do

END SUBROUTINE coeff_shell_dip



!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE coeff_shell_dip: Calculation of auxiliary coefficients
! to calculate the coefficent matrix to be solved
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE coeff_shell_dip_borghese(lambda,ref_index,v_req,m_epseq,nstop,neq,v_Rn,v_Vn,v_Zn,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index			! Wavelength and ref index
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vector of equals r
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Dielectri Function Matrix
INTEGER(lo), INTENT(IN) :: nstop				! Multipolar expansions
INTEGER(lo), INTENT(IN) :: neq			        ! Number of equal spheres

COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_Rn,v_Vn,v_Zn	! Auxiliary coefficients
INTEGER(lo), INTENT(OUT) :: error				! Error flag


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x							! Size Parameters Vector
COMPLEX(dbl), DIMENSION(neq-1) :: v_epsc							! Dielectric Function
COMPLEX(dbl), DIMENSION(0:neq) :: v_k							! Wave vectors
COMPLEX(dbl), DIMENSION(0:neq) :: v_m							! Normalized refractive index
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx						! Riccati bessel functions arguments
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z12,v_psi_z22			! Riccati bessel
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12			! 
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11,v_psi1_z12,v_psi1_z22			! 
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z01,v_csi1_z11,v_csi1_z12			! 
INTEGER(lo) :: m,n,nu,i,j								! Indexes
REAL(dbl) :: nr										! 

!Subroutine begins

!Remapping of the dielectric function, neq -> (neq-1), since neq includes dipoles
v_epsc_do: DO j=1,neq-1
			!Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
			!da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
			v_epsc(j)=CMPLX(m_epseq(neq-j,1),m_epseq(neq-j,2),KIND=dbl)
END DO v_epsc_do

!Size parameters
v_x(1)=(2*pi_d*ref_index*v_req(2))/lambda
v_x(2)=(2*pi_d*ref_index*v_req(1))/lambda

!Wavevectors, again neq->neq-1
v_k(0)=(2*pi_d*ref_index)/lambda
v_k(1:neq-1)=(2*pi_d*SQRT(v_epsc))/lambda

!Normalized ref index, again neq->neq-1
v_m(0)=(1.0D0,0.0D0)
v_m(1:neq-1)=SQRT(v_epsc)/CMPLX(ref_index,KIND=dbl)

!Function argument matrix, again neq->neq-1
mx_out_do: DO i=0,neq-1
	mx_in_do: DO j=1,neq-1

		m_mx(i,j)=v_m(i)*CMPLX(v_x(j),KIND=dbl)

	END DO mx_in_do
END DO mx_out_do



!Riccati Bessel Calculation psi
CALL psi_z_sub(nstop,m_mx(1,1),v_psi_z11,error)
CALL psi_z_sub(nstop,m_mx(1,2),v_psi_z12,error)
CALL psi_z_sub(nstop,m_mx(2,2),v_psi_z22,error)

psi_z_if: IF (error/=0) THEN
		error=1
		WRITE(*,10)
		10 FORMAT ("Error: psi_z_sub called by coeff_shell_borghese")
		RETURN
END IF psi_z_if	

!Riccati Bessel Calculation csi
CALL csi_z_sub(nstop,m_mx(0,1),v_csi_z01,error)
CALL csi_z_sub(nstop,m_mx(1,1),v_csi_z11,error)
CALL csi_z_sub(nstop,m_mx(1,2),v_csi_z12,error)

csi_z_if: IF (error/=0) THEN
		error=1
		WRITE(*,20)
		20 FORMAT ("Error: csi_z_sub called by coeff_shell_borghese")
		RETURN
END IF csi_z_if


!Riccati Bessel derivatives
der_do: DO n=1,nstop

	nr=REAL(n,dbl)
	!Psi
	v_psi1_z11(n)=v_psi_z11(n-1)-nr*v_psi_z11(n)/m_mx(1,1)
	v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)
	v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

	!Csi
	v_csi1_z01(n)=v_csi_z01(n-1)-nr*v_csi_z01(n)/m_mx(0,1)
	v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)
	v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)

END DO der_do



!Rn
v_Rn(1:2*nstop:2)=(v_k(1)*v_csi_z11(1:nstop)*v_csi1_z01 - v_k(0)*v_csi1_z11*v_csi_z01(1:nstop)) / &
                      &(v_k(0)*v_psi1_z11*v_csi_z01(1:nstop) - v_k(1)*v_psi_z11(1:nstop)*v_csi1_z01)
v_Rn(2:2*nstop:2)=(v_k(0)*v_csi_z11(1:nstop)*v_csi1_z01 - v_k(1)*v_csi1_z11*v_csi_z01(1:nstop)) / &
                      &(v_k(1)*v_psi1_z11*v_csi_z01(1:nstop) - v_k(0)*v_psi_z11(1:nstop)*v_csi1_z01)

!Vn
v_Vn(1:2*nstop:2)=(0.0D0,1.0D0)*( (v_k(2)/v_k(1))*v_psi_z22(1:nstop)*v_psi1_z12 - v_psi1_z22*v_psi_z12(1:nstop))
v_Vn(2:2*nstop:2)=(0.0D0,1.0D0)*( v_psi_z22(1:nstop)*v_psi1_z12 - (v_k(2)/v_k(1))*v_psi1_z22*v_psi_z12(1:nstop))

!Zn
v_Zn(1:2*nstop:2)=(0.0D0,1.0D0)*( v_psi1_z22*v_csi_z12(1:nstop) - (v_k(2)/v_k(1))*v_psi_z22(1:nstop)*v_csi1_z12)
v_Zn(2:2*nstop:2)=(0.0D0,1.0D0)*( (v_k(2)/v_k(1))*v_psi1_z22*v_csi_z12(1:nstop) - v_psi_z22(1:nstop)*v_csi1_z12)




END SUBROUTINE coeff_shell_dip_borghese



END MODULE sing_part















