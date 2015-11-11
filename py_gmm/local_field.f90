MODULE local_field

USE kinds
USE datatypes
USE basicsubs
USE sing_part
USE shared_data
USE vec_trans

CONTAINS
!******************************************************************************
!******************************************************************************
!******************************************************************************
!1) SUBROUTINE emn_sub: calcola i fattori di normalizzazione:li calcolo una volta
! per tutte visto che servono praticamente solo una volta
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE emn_sub(nstop,v_emn,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
INTEGER(lo), INTENT(IN) :: nstop							!Numero di espansioni multipolari
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_emn			!Vettore  normalizzazione
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
INTEGER(lo) :: m,n,i											!Indici e dimensione
REAL(dbl) :: nr,mr,logw											!Indici reali

!Inizio subroutine vera e propria
error=0
i=0
n_do: DO n=1,nstop

	nr=REAL(n,dbl)

	m_do: DO m=-n,n
	
		i=i+1
		mr=REAL(m,dbl)
		
		logw=lnf(nr-mr,error)-lnf(nr+mr,error)
		v_emn(i)=((0.0D0,1.0D0)**n)*SQRT(EXP(logw)*(2*nr+1.0D0)/(nr*(nr+1.0D0)))
	
	END DO m_do
END DO n_do

END SUBROUTINE emn_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE int_field: calcolo il campo interno alla sfera in coordinate
! sferiche locali
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE int_field(v_dmncmn,v_emn,nstop,rho,theta,phi,Er,Etheta,Ephi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn,v_emn	!Vettore campo interno e normalizzazione
INTEGER(lo), INTENT(IN) :: nstop							!Numero di espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: rho								!k*r con r coordinata sferica
REAL(dbl), INTENT(IN) :: theta,phi 							!Coordinate sferiche locali
COMPLEX(dbl), INTENT(OUT) :: Er,Etheta,Ephi					!Campo calcolato in coordinate sferiche
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_legmn	!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_eta,v_qn			!Vettori funzioni radiali
INTEGER(lo) :: n,m,i,dimens									!Indici e dimensione
REAL(dbl) :: nr,mr												!Indici reali

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_legmn(1:dimens),v_eta(0:nstop),v_qn(1:nstop))

!Funzione angolare Pi_mn
CALL pi_mn(nstop,theta,v_pimn,error)

pimn_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF pimn_if

!Funzione angolare Tau_mn
CALL tau_mn(nstop,theta,v_pimn,v_taumn,error)

taumn_if: IF (error/=0) THEN																	
    			WRITE(*,20) 
    			20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF taumn_if

!Funzione radiale eta_n
CALL eta_z_sub(nstop,rho,v_eta,error)

eta_if: IF (error/=0) THEN																	
    			WRITE(*,30) 
    			30 FORMAT ("Si e' verificato un errore in eta_z_sub chiamata da int_field: il programma termina ora...") 
       			STOP
END IF eta_if

!Funzione radiale q_n
CALL qn_z_sub(1,nstop,rho,v_qn,error)

qn_if: IF (error/=0) THEN																	
    			WRITE(*,40) 
    			40 FORMAT ("Si e' verificato un errore in qn_z_sub chiamata da int_field: il programma termina ora...") 
       			STOP
END IF qn_if

!Funzione angolare legendre
CALL leg_mn(nstop,theta,v_legmn,error)

legmn_if: IF (error/=0) THEN																	
    			WRITE(*,50) 
    			50 FORMAT ("Si e' verificato un errore in leg_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF legmn_if


!Calcolo del campo radiale
i=0
Er=(0.0D0,0.0D0)
Ern_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Erm_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Er=Er+v_emn(i)*nr*(nr+1.0D0)*v_dmncmn(2*i-1)*v_legmn(i)*v_eta(n)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Erm_do
END DO Ern_do

Er=(0.0D0,-1.0D0)*Er

!Calcolo del campo theta
i=0
Etheta=(0.0D0,0.0D0)
Ethetan_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ethetam_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Etheta=Etheta+v_emn(i)*( (0.0D0,-1.0D0)*v_dmncmn(2*i-1)*(v_qn(n)-nr)*v_taumn(i) + &
			 & rho*v_dmncmn(2*i)*v_pimn(i)) * v_eta(n)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ethetam_do
END DO Ethetan_do

!Calcolo del campo phi
i=0
Ephi=(0.0D0,0.0D0)
Ephin_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ephim_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Ephi=Ephi+v_emn(i)*( (0.0D0,1.0D0)*rho*v_dmncmn(2*i)*v_taumn(i) + &
			 & v_dmncmn(2*i-1)*(v_qn(n)-nr)*v_pimn(i)) * v_eta(n)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ephim_do
END DO Ephin_do


DEALLOCATE(v_pimn,v_taumn,v_legmn,v_eta,v_qn)

END SUBROUTINE int_field



!******************************************************************************
!******************************************************************************
!******************************************************************************
!2c) SUBROUTINE inc_field: calcolo il campo incidente in coord sferiche
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE inc_field(v_pmnqmn,v_emn,nstop,rho,theta,phi,Er,Etheta,Ephi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_pmnqmn,v_emn	!Vettore campo interno e normalizzazione
INTEGER(lo), INTENT(IN) :: nstop							!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: rho								!k*r con r coordinata sferica
REAL(dbl), INTENT(IN) :: theta,phi 							!Coordinate sferiche locali
COMPLEX(dbl), INTENT(OUT) :: Er,Etheta,Ephi					!Campo calcolato in coordinate sferiche
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_legmn	!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_eta,v_qn			!Vettori funzioni radiali
INTEGER(lo) :: n,m,i,dimens									!Indici e dimensione
REAL(dbl) :: nr,mr												!Indici reali

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_legmn(1:dimens),v_eta(0:nstop),v_qn(1:nstop))

!Funzione angolare Pi_mn
CALL pi_mn(nstop,theta,v_pimn,error)

pimn_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF pimn_if

!Funzione angolare Tau_mn
CALL tau_mn(nstop,theta,v_pimn,v_taumn,error)

taumn_if: IF (error/=0) THEN																	
    			WRITE(*,20) 
    			20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF taumn_if

!Funzione radiale eta_n
CALL eta_d_sub(nstop,rho,v_eta,error)

eta_if: IF (error/=0) THEN																	
    			WRITE(*,30) 
    			30 FORMAT ("Si e' verificato un errore in eta_z_sub chiamata da int_field: il programma termina ora...") 
       			STOP
END IF eta_if

!Funzione radiale q_n
CALL qn_d_sub(1,nstop,rho,v_qn,error)

qn_if: IF (error/=0) THEN																	
    			WRITE(*,40) 
    			40 FORMAT ("Si e' verificato un errore in qn_z_sub chiamata da int_field: il programma termina ora...") 
       			STOP
END IF qn_if

!Funzione angolare legendre
CALL leg_mn(nstop,theta,v_legmn,error)

legmn_if: IF (error/=0) THEN																	
    			WRITE(*,50) 
    			50 FORMAT ("Si e' verificato un errore in leg_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF legmn_if


!Calcolo del campo radiale
i=0
Er=(0.0D0,0.0D0)
Ern_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Erm_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Er=Er+v_emn(i)*nr*(nr+1.0D0)*v_pmnqmn(2*i-1)*v_legmn(i)*v_eta(n)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Erm_do
END DO Ern_do

Er=(0.0D0,-1.0D0)*Er

!Calcolo del campo theta
i=0
Etheta=(0.0D0,0.0D0)
Ethetan_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ethetam_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Etheta=Etheta+v_emn(i)*( (0.0D0,-1.0D0)*v_pmnqmn(2*i-1)*(v_qn(n)-nr)*v_taumn(i) + &
			 & rho*v_pmnqmn(2*i)*v_pimn(i)) * v_eta(n)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ethetam_do
END DO Ethetan_do

!Calcolo del campo phi
i=0
Ephi=(0.0D0,0.0D0)
Ephin_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ephim_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Ephi=Ephi+v_emn(i)*( (0.0D0,1.0D0)*rho*v_pmnqmn(2*i)*v_taumn(i) + &
			 & v_pmnqmn(2*i-1)*(v_qn(n)-nr)*v_pimn(i)) * v_eta(n)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ephim_do
END DO Ephin_do


DEALLOCATE(v_pimn,v_taumn,v_legmn,v_eta,v_qn)

END SUBROUTINE inc_field






!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE ext_field: calcolo il campo esterno alla sfera in coordinate
! sferiche locali
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE ext_field(v_amnbmn,v_emn,nstop,rho,theta,phi,Er,Etheta,Ephi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn,v_emn	!Vettore campo interno e normalizzazione
INTEGER(lo), INTENT(IN) :: nstop							!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: rho								!k*r con r coordinata sferica
REAL(dbl), INTENT(IN) :: theta,phi 							!Coordinate sferiche locali
COMPLEX(dbl), INTENT(OUT) :: Er,Etheta,Ephi					!Campo calcolato in coordinate sferiche
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_legmn	!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h,v_dh				!Vettori funzioni radiali
INTEGER(lo) :: n,m,i,dimens									!Indici e dimensione
REAL(dbl) :: nr,mr												!Indici reali

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

!Input monitoring
!!$WRITE(*,*) 
!!$WRITE(*,*) "inside ext_field"
!!$WRITE(*,*) "v_amnbmn",v_amnbmn
!!$WRITE(*,*) "v_emn",v_emn
!!$WRITE(*,*) "nstop",nstop
!!$WRITE(*,*) "rho",rho
!!$WRITE(*,*) "theta",theta
!!$WRITE(*,*) "phi",phi


ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_legmn(1:dimens),v_h(0:nstop),v_dh(1:nstop))

!Funzione angolare Pi_mn
CALL pi_mn(nstop,theta,v_pimn,error)

pimn_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da ext_field: il programma termina ora...") 
       			STOP
END IF pimn_if

!!$WRITE(*,*)
!!$DO i=1,dimens
!!$   WRITE(*,*) "i,pi_mn",i,v_pimn(i)
!!$END DO


!Funzione angolare Tau_mn
CALL tau_mn(nstop,theta,v_pimn,v_taumn,error)

taumn_if: IF (error/=0) THEN																	
    			WRITE(*,20) 
    			20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da ext_field: il programma termina ora...") 
       			STOP
END IF taumn_if

!!$WRITE(*,*)
!!$DO i=1,dimens
!!$   WRITE(*,*) "i,tau_mn",i,v_taumn(i)
!!$END DO

!Funzione angolare legendre
CALL leg_mn(nstop,theta,v_legmn,error)

legmn_if: IF (error/=0) THEN																	
    			WRITE(*,40) 
    			40 FORMAT ("Si e' verificato un errore in leg_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF legmn_if

!!$WRITE(*,*)
!!$DO i=1,dimens
!!$   WRITE(*,*) "i,leg_mn",i,v_legmn(i)
!!$END DO

!Funzione radiale hankel1_n
CALL hankel1_d_sub(nstop,rho,v_h,error)

hankel_if: IF (error/=0) THEN																	
    			WRITE(*,30) 
    			30 FORMAT ("Si e' verificato un errore in hankel1_d_sub chiamata da ext_field: il programma termina ora...") 
       			STOP
END IF hankel_if

!Funzione radiale (1/rho)*D[rho*h_n[rho],rho]
dh_do: DO n=1,nstop
	nr=REAL(n,dbl)
	v_dh(n)=v_h(n-1)-nr*v_h(n)/rho
END DO dh_do

!!$WRITE(*,*)
!!$DO i=1,nstop
!!$   WRITE(*,*) "i,hn",i,v_h(i)
!!$END DO
!!$WRITE(*,*)
!!$DO i=0,nstop
!!$   WRITE(*,*) "i,hn",i,v_dh(i)
!!$END DO

!Calcolo del campo radiale
i=0
Er=(0.0D0,0.0D0)
Ern_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Erm_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Er=Er+v_emn(i)*nr*(nr+1.0D0)*v_amnbmn(2*i-1)*v_legmn(i)*(v_h(n)/rho)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Erm_do
END DO Ern_do

Er=(0.0D0,1.0D0)*Er

!Calcolo del campo theta
i=0
Etheta=(0.0D0,0.0D0)
Ethetan_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ethetam_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Etheta=Etheta+v_emn(i)*( (0.0D0,1.0D0)*v_amnbmn(2*i-1)*v_taumn(i)*v_dh(n)  -&
			 & v_amnbmn(2*i)*v_pimn(i)*v_h(n))*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ethetam_do
END DO Ethetan_do

!Calcolo del campo phi
i=0
Ephi=(0.0D0,0.0D0)
Ephin_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ephim_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Ephi=Ephi-v_emn(i)*( v_amnbmn(2*i-1)*v_pimn(i)*v_dh(n) + &
			 & v_amnbmn(2*i)*v_taumn(i)*v_h(n))*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ephim_do
END DO Ephin_do


DEALLOCATE(v_pimn,v_taumn,v_legmn,v_h,v_dh)

END SUBROUTINE ext_field



!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE Exyz_sub: dato un punto (x,y,z) nello spazio, calcolo il campo
! elettrico totale in coordinate cartesiane,sia che sia dentro sia che sia fuori
! una sfera
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Exyz_sub(v_kinc,v_Einc,flaginc,nstop,ratio,lambda,x,y,z,v_amnbmn,&
			 &v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,Ex,Ey,Ez,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc		!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc						!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz			!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn		!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn		!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn			!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq		!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index						!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt		!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio				!Numero di espansioni multipolari e raggi non calcolo
REAL(dbl), INTENT(IN) :: x,y,z							!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda							!Lunghezza d'onda incidente
COMPLEX(dbl), INTENT(OUT) :: Ex,Ey,Ez					!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Exinc,Eyinc,Ezinc,pshift					!Onda incidente cartesiana
COMPLEX(dbl) :: Er,Etheta,Ephi						!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: rhoc,mc								!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase								!Rho reale
REAL(dbl) :: r,theta,phi 							!Coordinate sferiche locali
REAL(dbl) :: xl,yl,zl	 							!Coordinate cartesiane locali
REAL(dbl) :: d,radius								!Distanza punto-sfere,raggio sfera
INTEGER(lo) :: n,m,i,dimens,ns,low_b,up_b,exit_flag			!Indici e dimensione
REAL(dbl) :: nr,mr									!Indici reali

!Inizio subroutine vera e propria

!Calcolo il numero di sfere e dimens
dimens=nstop*(nstop+2)
ns=SIZE(m_xyz(:,1))
exit_flag=0

!--------------------------------------------
!Se sono qui allora non sono in nessuna sfera
!--------------------------------------------
Ex=(0.0D0,0.0D0)
Ey=(0.0D0,0.0D0)
Ez=(0.0D0,0.0D0)

!Input monitoring
!!$IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$
!!$   WRITE(*,*) 'v_kinc',v_kinc
!!$   WRITE(*,*) 'v_Einc',v_Einc
!!$   WRITE(*,*) 'flaginc',flaginc
!!$   WRITE(*,*) 'nstop',nstop
!!$   WRITE(*,*) 'ratio',ratio
!!$   WRITE(*,*) 'lambda',lambda
!!$   WRITE(*,*) 'x,y,z',x,y,z
!!$   WRITE(*,*) 'v_amnbmn',v_amnbmn
!!$   WRITE(*,*) 'v_dmncmn',v_dmncmn
!!$   WRITE(*,*) 'v_emn',v_emn
!!$   WRITE(*,*) 'm_xyz',m_xyz
!!$   WRITE(*,*) 'm_eps',m_epseq
!!$   WRITE(*,*) 'v_req',v_req
!!$   WRITE(*,*) 'ref_index',ref_index
!!$   
!!$END IF


!Faccio un ciclo per vedere se interseco sulle sfere
int_do: DO i=1,ns

			!Calcolo raggio e distanza
			radius=v_req(v_patt(i))
			d=SQRT((x-m_xyz(i,1))**2 + (y-m_xyz(i,2))**2 + (z-m_xyz(i,3))**2)
			
			!Se non sono in una sfera, controllo se sono in un'altra
			IF (d>radius) CYCLE int_do
			
			!-------------------------------------------
			!Se sono qui allora sono nella i-esima sfera
			!-------------------------------------------
						
			!Calcolo le coordinate sferiche locali
			xl=x-m_xyz(i,1)
			yl=y-m_xyz(i,2)
			zl=z-m_xyz(i,3)
			CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)
			
			c_s_if: IF (error/=0) THEN
				error=1
				WRITE(*,10) 
    			10 FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")                           
       			RETURN
			END IF c_s_if
			
			!Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
			mc=SQRT(CMPLX(m_epseq(v_patt(i),1),m_epseq(v_patt(i),2),KIND=dbl))
			rhoc=2*Pi_D*mc*r/lambda
			
			!Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
			low_b=2*dimens*(i-1)+1
			up_b=2*dimens*i
			CALL int_field(v_dmncmn(low_b:up_b),v_emn,nstop,rhoc,theta,phi,Er,Etheta,Ephi,error)
			
			intfield_if: IF (error/=0) THEN
				error=1
				WRITE(*,20) 
    			20 FORMAT ("Errore in int_field chiamata da Exyz_sub: il programma termina ora...")                           
       			RETURN
			END IF intfield_if
			
			!Conversione in coordinate cartesiane e uscita dalla subroutine
			Ex=Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
			Ey=Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
			Ez=Er*COS(theta) - Etheta*SIN(theta)
			
			exit_flag=1

END DO int_do



ext_do: DO i=1,ns

        IF (exit_flag==1) EXIT ext_do 

	!Calcolo le coordinate sferiche locali
	xl=x-m_xyz(i,1)
	yl=y-m_xyz(i,2)
	zl=z-m_xyz(i,3)
	CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

	c_s_if1: IF (error/=0) THEN
		WRITE(*,30) 
		30 FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")                           
		RETURN
	END IF c_s_if1
    
!!$        !Input monitoring
!!$        IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "xl,yl,zl",xl,yl,zl
!!$          WRITE(*,*) "r,theta,phi",r,theta,phi
!!$        END IF

    IF ( ( r/( v_req(v_patt(i)) ) ) > ratio) CYCLE ext_do
	
!!$        !Input monitoring
!!$        IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "computing external field, ns:",i
!!$        END IF

	!Calcolo  rho per chiamare la subroutine per il calcolo del campo
	rhor=2*Pi_D*ref_index*r/lambda
	!Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
	low_b=2*dimens*(i-1)+1
	up_b=2*dimens*i

!!$       !Input monitoring
!!$       IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "v_amnbmn(low_b:up_b)",v_amnbmn(low_b:up_b)
!!$          WRITE(*,*) "v_emn",v_emn
!!$          WRITE(*,*) "nstop",nstop
!!$          WRITE(*,*) "rhor",rhor
!!$          WRITE(*,*) "theta",theta
!!$          WRITE(*,*) "phi",phi
!!$       END IF

	CALL ext_field(v_amnbmn(low_b:up_b),v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)
	
	extfield_if: IF (error/=0) THEN
		WRITE(*,40) 
		40 FORMAT ("Errore in ext_field chiamata da Exyz_sub: il programma termina ora...")                           
		RETURN
	END IF extfield_if
	
!!$      !Input monitoring
!!$       IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "Er,Etheta,Ephi",Er,Etheta,Ephi
!!$       END IF

	!Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
	Ex=Ex+Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
	Ey=Ey+Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
	Ez=Ez+Er*COS(theta) - Etheta*SIN(theta)

END DO ext_do

!Adesso, se la flag e' quella attesa, aggiungo anche il campo incidente, tanto il calcolo e' banale			
inc_if: IF ((flaginc=="yes") .AND. (exit_flag==0)) THEN

	phase=x*v_kinc(1)+y*v_kinc(2)+z*v_kinc(3)
	pshift=EXP((0.0D0,1.0D0)*phase)
	Exinc=CMPLX(v_Einc(1),0.0D0)*pshift
	Eyinc=CMPLX(v_Einc(2),0.0D0)*pshift
	Ezinc=CMPLX(v_Einc(3),0.0D0)*pshift
	Ex=Ex+Exinc
	Ey=Ey+Eyinc
	Ez=Ez+Ezinc

END IF inc_if

END SUBROUTINE Exyz_sub





!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE Exyz_per_sub: dato un punto (x,y,z) nello spazio, calcolo il campo
! elettrico totale in coordinate cartesiane,sia che sia dentro sia che sia fuori
! una sfera.Qui lo calcolo per la soluzione periodica compatta
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Exyz_per_sub(v_kinc,v_Einc,flaginc,nstop,ratio,lambda,x,y,z,v_amnbmn,&
			 &v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,Ex,Ey,Ez,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc	!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc				!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz			!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn		!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn		!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn			!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq		!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index					!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt		!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio				!Numero di espansioni multipolari e raggi non calcolo
REAL(dbl), INTENT(IN) :: x,y,z					!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda					!Lunghezza d'onda incidente
COMPLEX(dbl), INTENT(OUT) :: Ex,Ey,Ez				!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error				!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Exinc,Eyinc,Ezinc,pshift				!Onda incidente cartesiana
COMPLEX(dbl) :: Er,Etheta,Ephi					!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: rhoc,mc							!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase							!Rho reale
REAL(dbl) :: r,theta,phi 						!Coordinate sferiche locali
REAL(dbl) :: xg,yg,zg,xl,yl,zl					!Coordinate cartesiane generali e locali
REAL(dbl) :: d,radius							!Distanza punto-sfere,raggio sfera
INTEGER(lo) :: n,m,i,dimens,ns,low_b,up_b			!Indici e dimensione
INTEGER(lo) :: na,nb							!Indici e dimensione
REAL(dbl) :: nr,mr							!Indici reali

!Inizio subroutine vera e propria

!Calcolo il numero di sfere e dimens
dimens=nstop*(nstop+2)
ns=SIZE(m_xyz(:,1))

!Faccio un ciclo per vedere se interseco sulle sfere
int_do: DO i=1,ns

	row_do: DO na=namin,namax
		col_do: DO nb=nbmin,nbmax

			!Calcolo raggio e distanza
			radius=v_req(v_patt(i))

			xg=m_xyz(i,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
			yg=m_xyz(i,2) + REAL(nb,dbl)*by
			zg=m_xyz(i,3)

			d=SQRT((x-xg)**2 + (y-yg)**2 + (z-zg)**2)

			!Se non sono in una sfera, controllo se sono in un'altra
			IF (d>radius) CYCLE col_do

			!-------------------------------------------
			!Se sono qui allora sono nella i-esima sfera
			!-------------------------------------------

			!Calcolo le coordinate sferiche locali
			xl=x-xg
			yl=y-yg
			zl=z-zg
			CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

			c_s_if: IF (error/=0) THEN
				error=1
				WRITE(*,10) 
			10 FORMAT ("Errore in cart_spher chiamata da Exyz_per_sub: il programma termina ora...")                           
				RETURN
			END IF c_s_if

			!Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
			mc=SQRT(CMPLX(m_epseq(v_patt(i),1),m_epseq(v_patt(i),2),KIND=dbl))
			rhoc=2*Pi_D*mc*r/lambda

			!Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
			low_b=2*dimens*(i-1)+1
			up_b=2*dimens*i
			CALL int_field(v_dmncmn(low_b:up_b),v_emn,nstop,rhoc,theta,phi,Er,Etheta,Ephi,error)

			intfield_if: IF (error/=0) THEN
				error=1
				WRITE(*,20) 
			20 FORMAT ("Errore in int_field chiamata da Exyz_sub: il programma termina ora...")                           
				RETURN
			END IF intfield_if

			!Conversione in coordinate cartesiane e uscita dalla subroutine
			Ex=Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
			Ey=Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
			Ez=Er*COS(theta) - Etheta*SIN(theta)

			RETURN

		END DO col_do
	END DO row_do
END DO int_do

!--------------------------------------------
!Se sono qui allora non sono in nessuna sfera
!--------------------------------------------
Ex=(0.0D0,0.0D0)
Ey=(0.0D0,0.0D0)
Ez=(0.0D0,0.0D0)

ext_do: DO i=1,ns

!	row_do_ext: DO na=0,nacell-1
	row_do_ext: DO na=namin,namax
!		col_do_ext: DO nb=0,nbcell-1
		col_do_ext: DO nb=nbmin,nbmax

!			WRITE(*,*) 'nacell', nacell
!			WRITE(*,*) 'nbcell', nbcell

			xg=m_xyz(i,1) + REAL(na,dbl)*ax + REAL(nb,dbl)*bx
			yg=m_xyz(i,2) + REAL(nb,dbl)*by
			zg=m_xyz(i,3)

			!Calcolo le coordinate sferiche locali
			xl=x-xg
			yl=y-yg
			zl=z-zg
			CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

			c_s_if1: IF (error/=0) THEN
				WRITE(*,30) 
				30 FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")                           
				RETURN
			END IF c_s_if1

			IF ( v_req(v_patt(i))>1 ) THEN
				IF ( ( r/( v_req(v_patt(i)) ) ) > ratio) CYCLE col_do_ext
			ELSE
				IF ( ( r/( v_req(v_patt(i)) ) ) > ratio*1000) CYCLE col_do_ext
			END IF

			!Calcolo  rho per chiamare la subroutine per il calcolo del campo
			rhor=2*Pi_D*ref_index*r/lambda
			!Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
			low_b=2*dimens*(i-1)+1
			up_b=2*dimens*i
			CALL ext_field(v_amnbmn(low_b:up_b),v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)

			extfield_if: IF (error/=0) THEN
				WRITE(*,40) 
				40 FORMAT ("Errore in ext_field chiamata da Exyz_sub: il programma termina ora...")                           
				RETURN
			END IF extfield_if

			!Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
			Ex=Ex+Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
			Ey=Ey+Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
			Ez=Ez+Er*COS(theta) - Etheta*SIN(theta)

		END DO col_do_ext
	END DO row_do_ext
END DO ext_do

!Adesso, se la flag e' quella attesa, aggiungo anche il campo incidente, tanto il calcolo e' banale			
inc_if: IF (flaginc=="yes") THEN

	phase=x*v_kinc(1)+y*v_kinc(2)+z*v_kinc(3)
	pshift=EXP((0.0D0,1.0D0)*phase)
	Exinc=CMPLX(v_Einc(1),0.0D0)*pshift
	Eyinc=CMPLX(v_Einc(2),0.0D0)*pshift
	Ezinc=CMPLX(v_Einc(3),0.0D0)*pshift
	Ex=Ex+Exinc
	Ey=Ey+Eyinc
	Ez=Ez+Ezinc

END IF inc_if

END SUBROUTINE Exyz_per_sub






!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE Efield_sub: calcolo il campo su tutto un piano in coordinate  
! cartesiane
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efield_sub(v_kinc,v_Einc,flaginc,xmin,xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep,&
     & nstop,ratio,lambda,v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,m_E,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc		!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc						!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn			!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn				!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq			!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req				!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index							!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt			!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio    				!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: xmin,xmax,ymin,ymax,zmin,zmax		!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda								!Lunghezza d'onda incidente
INTEGER(lo), INTENT(IN) :: xstep,ystep,zstep				!Numero di step
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_E			!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,pshift					!Campo calcolato in coordinate sferiche e phase shift c inc
COMPLEX(dbl) :: rhoc,mc								!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase								!Rho reale e fase campo inc
REAL(dbl) :: r,theta,phi 							!Coordinate sferiche locali
REAL(dbl) :: x,y,z,dx,dy,dz							!Coordinate cartesiane locali,delta x,y,z
REAL(dbl) :: d,radius								!Distanza punto-sfere,raggio sfera
INTEGER(lo) :: i,j,chunk      						!Indici e dimensione
REAL(dbl) :: nr,mr								!Indici reali
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pos_a, v_pos_b			!First and second dimension position vectors

!Inizio subroutine vera e propria

!Controllo di avere un piano ben definito
plane_if: IF ((xmin/=xmax) .AND. (ymin/=ymax) .AND. (zmin/=zmax)) THEN
     error=1
     WRITE(*,10) 
10   FORMAT ("Errore in Efield_sub: nessun piano definito, il programma termina ora...")                           
     RETURN
END IF plane_if

!Controllo di avere i bounds ben definiti
bounds_if: IF ((xmin>xmax) .OR. (ymin>ymax) .OR. (zmin>zmax)) THEN
     error=1
     WRITE(*,20) 
20   FORMAT ("Errore in Efield_sub: bounds mal definiti, il programma termina ora...")                           
     RETURN
END IF bounds_if

!Controllo di avere gli steps ben definiti
steps_if: IF ((xstep<2) .OR. (ystep<2) .OR. (zstep<2)) THEN
     error=1
     WRITE(*,30) 
30   FORMAT ("Errore in Efield_sub: steps mal definiti, il programma termina ora...")                           
     RETURN
END IF steps_if

!Calcolo i delta a seconda degli step
dx=(xmax-xmin)/REAL(xstep-1,dbl)
dy=(ymax-ymin)/REAL(ystep-1,dbl)
dz=(zmax-zmin)/REAL(zstep-1,dbl)

!Comincio i conti del campo e un grande if 
x=xmin
y=ymin
z=zmin


!Filling the spatial grid for the calculation
grid_if: IF (xmin==xmax) THEN

     !Allocating and initializing the space vector
     ALLOCATE(v_pos_a(1:ystep))				!a=y
     ALLOCATE(v_pos_b(1:zstep))				!b=z
     v_pos_a=0.0D0
     v_pos_b=0.0D0

     !Filling first dimension
     y=ymin
     y_yz_grid_do: DO i=1,ystep
          v_pos_a(i)=y
          y=y+dy
     END DO y_yz_grid_do

     !Filling second dimension
     z=zmin
     z_yz_grid_do: DO i=1,zstep
          v_pos_b(i)=z
          z=z+dz
     END DO z_yz_grid_do

ELSE IF (ymin==ymax) THEN

     !Allocating and initializing the space vector
     ALLOCATE(v_pos_a(1:xstep))				!a=x
     ALLOCATE(v_pos_b(1:zstep))				!b=z
     v_pos_a=0.0D0
     v_pos_b=0.0D0

     !Filling first dimension
     x=xmin
     x_xz_grid_do: DO i=1,xstep
          v_pos_a(i)=x
          x=x+dx
     END DO x_xz_grid_do

     !Filling second dimension
     z=zmin
     z_xz_grid_do: DO i=1,zstep
          v_pos_b(i)=z
          z=z+dz
     END DO z_xz_grid_do

ELSE 

     !Allocating and initializing the space vector
     ALLOCATE(v_pos_a(1:xstep))				!a=x
     ALLOCATE(v_pos_b(1:ystep))				!b=y
     v_pos_a=0.0D0
     v_pos_b=0.0D0

     !Filling first dimension
     x=xmin
     x_xy_grid_do: DO i=1,xstep
          v_pos_a(i)=x
          x=x+dx
     END DO x_xy_grid_do

     !Filling second dimension
     y=ymin
     y_xy_grid_do: DO i=1,ystep
          v_pos_b(i)=y
          y=y+dy
     END DO y_xy_grid_do

END IF grid_if



!Loop for the calculation of the field in each grid point
xyz_if: IF (xmin==xmax) THEN

     !Comincio i conti del campo e un grande if 

     !$OMP PARALLEL PRIVATE(i,j)

     !Choosing the chunk size
     chunk=32
     chunk_if: IF (chunk==0) THEN
          chunk=1
     END IF chunk_if

     !$OMP DO SCHEDULE(guided,chunk)
     y_yz_do: DO i=1,ystep
          z_yz_do: DO j=1,zstep

               m_E(i,j,1)=v_pos_a(i)		!y
               m_E(i,j,2)=v_pos_b(j)		!z
               CALL Exyz_sub(v_kinc,v_Einc,flaginc,nstop,&
                    &ratio,lambda,x,v_pos_a(i),v_pos_b(j),v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

               Exyz1_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,40) 
40                  FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
               END IF Exyz1_if

          END DO z_yz_do
     END DO y_yz_do
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL


ELSE IF (ymin==ymax) THEN

     !$OMP PARALLEL PRIVATE(i,j)

     !Choosing the chunk size
     chunk=32
     chunk_if1: IF (chunk==0) THEN
          chunk=1
     END IF chunk_if1

     !$OMP DO SCHEDULE(guided,chunk)
     x_xz_do: DO i=1,xstep
          z_xz_do: DO j=1,zstep

               m_E(i,j,1)=v_pos_a(i)		!x
               m_E(i,j,2)=v_pos_b(j)		!z
               CALL Exyz_sub(v_kinc,v_Einc,flaginc,nstop,&
                    &ratio,lambda,v_pos_a(i),y,v_pos_b(j),v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

               Exyz2_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,50) 
50                  FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
               END IF Exyz2_if

          END DO z_xz_do
     END DO x_xz_do
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL
     

ELSE 

     !$OMP PARALLEL PRIVATE(i,j)

     !Choosing the chunk size
     chunk=32
     chunk_if2: IF (chunk==0) THEN
          chunk=1
     END IF chunk_if2

     !$OMP DO SCHEDULE(guided,chunk)
     x_xy_do: DO i=1,xstep
          y_xy_do: DO j=1,ystep

               m_E(i,j,1)=v_pos_a(i)		!x
               m_E(i,j,2)=v_pos_b(j)		!y
               CALL Exyz_sub(v_kinc,v_Einc,flaginc,nstop,&
                    &ratio,lambda,v_pos_a(i),v_pos_b(j),z,v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

               Exyz3_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,60) 
60                  FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
               END IF Exyz3_if

          END DO y_xy_do
     END DO x_xy_do
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

END IF xyz_if


final_if: IF (error/=0) THEN
     error=1
     WRITE(*,70) 
70   FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
     RETURN
END IF final_if

END SUBROUTINE Efield_sub





!******************************************************************************
!******************************************************************************
!******************************************************************************
!4bis) SUBROUTINE Efield_per_sub: calcolo il campo su tutto un piano in coordinate  
! cartesiane, nel caso periodico compatto
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efield_per_sub(v_kinc,v_Einc,flaginc,xmin,xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep,&
& nstop,ratio,lambda,v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,m_E,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc		!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc					!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn			!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn				!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq			!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req				!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index							!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt			!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio    				!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: xmin,xmax,ymin,ymax,zmin,zmax		!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda								!Lunghezza d'onda incidente
INTEGER(lo), INTENT(IN) :: xstep,ystep,zstep				!Numero di step
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_E			!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,pshift					!Campo calcolato in coordinate sferiche e phase shift c inc
COMPLEX(dbl) :: rhoc,mc								!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase								!Rho reale e fase campo inc
REAL(dbl) :: r,theta,phi 							!Coordinate sferiche locali
REAL(dbl) :: x,y,z,dx,dy,dz							!Coordinate cartesiane locali,delta x,y,z
REAL(dbl) :: d,radius								!Distanza punto-sfere,raggio sfera
INTEGER(lo) :: i,j								!Indici e dimensione
REAL(dbl) :: nr,mr								!Indici reali
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pos_a, v_pos_b			!First and second dimension position vectors

!Inizio subroutine vera e propria

!Controllo di avere un piano ben definito
plane_if: IF ((xmin/=xmax) .AND. (ymin/=ymax) .AND. (zmin/=zmax)) THEN
	error=1
	WRITE(*,10) 
	10 FORMAT ("Errore in Efield_sub: nessun piano definito, il programma termina ora...")                           
	RETURN
END IF plane_if

!Controllo di avere i bounds ben definiti
bounds_if: IF ((xmin>xmax) .OR. (ymin>ymax) .OR. (zmin>zmax)) THEN
	error=1
	WRITE(*,20) 
	20 FORMAT ("Errore in Efield_sub: bounds mal definiti, il programma termina ora...")                           
	RETURN
END IF bounds_if

!Controllo di avere gli steps ben definiti
steps_if: IF ((xstep<2) .OR. (ystep<2) .OR. (zstep<2)) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efield_sub: steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if

!Calcolo i delta a seconda degli step
dx=(xmax-xmin)/REAL(xstep-1,dbl)
dy=(ymax-ymin)/REAL(ystep-1,dbl)
dz=(zmax-zmin)/REAL(zstep-1,dbl)

!Comincio i conti del campo e un grande if 
x=xmin
y=ymin
z=zmin

!Filling the spatial grid for the calculation
grid_if: IF (xmin==xmax) THEN

	!Allocating and initializing the space vector
	ALLOCATE(v_pos_a(1:ystep))				!a=y
	ALLOCATE(v_pos_b(1:zstep))				!b=z
	v_pos_a=0.0D0
	v_pos_b=0.0D0

	!Filling first dimension
	y=ymin
	y_yz_grid_do: DO i=1,ystep
		v_pos_a(i)=y
		y=y+dy
	END DO y_yz_grid_do

	!Filling second dimension
	z=zmin
	z_yz_grid_do: DO i=1,zstep
		v_pos_b(i)=z
		z=z+dz
	END DO z_yz_grid_do

ELSE IF (ymin==ymax) THEN

	!Allocating and initializing the space vector
	ALLOCATE(v_pos_a(1:xstep))				!a=x
	ALLOCATE(v_pos_b(1:zstep))				!b=z
	v_pos_a=0.0D0
	v_pos_b=0.0D0

	!Filling first dimension
	x=xmin
	x_xz_grid_do: DO i=1,xstep
		v_pos_a(i)=x
		x=x+dx
	END DO x_xz_grid_do

	!Filling second dimension
	z=zmin
	z_xz_grid_do: DO i=1,zstep
		v_pos_b(i)=z
		z=z+dz
	END DO z_xz_grid_do

ELSE 

	!Allocating and initializing the space vector
	ALLOCATE(v_pos_a(1:xstep))				!a=x
	ALLOCATE(v_pos_b(1:ystep))				!b=y
	v_pos_a=0.0D0
	v_pos_b=0.0D0

	!Filling first dimension
	x=xmin
	x_xy_grid_do: DO i=1,xstep
		v_pos_a(i)=x
		x=x+dx
	END DO x_xy_grid_do

	!Filling second dimension
	y=ymin
	y_xy_grid_do: DO i=1,ystep
		v_pos_b(i)=y
		y=y+dy
	END DO y_xy_grid_do

END IF grid_if



!Loop for the calculation of the field in each grid point
xyz_if: IF (xmin==xmax) THEN

	!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static,1)
	y_yz_do: DO i=1,ystep
		z_yz_do: DO j=1,zstep

			m_E(i,j,1)=v_pos_a(i)		!y
			m_E(i,j,2)=v_pos_b(j)		!z
			CALL Exyz_per_sub(v_kinc,v_Einc,flaginc,nstop,&
			&ratio,lambda,x,v_pos_a(i),v_pos_b(j),v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
			& m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

			Exyz1_if: IF (error/=0) THEN
				error=1
				WRITE(*,40) 
			40 FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
			END IF Exyz1_if

		END DO z_yz_do
	END DO y_yz_do
	!OMP END PARALLEL DO

ELSE IF (ymin==ymax) THEN

	!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static,1)
	x_xz_do: DO i=1,xstep
		z_xz_do: DO j=1,zstep

			m_E(i,j,1)=v_pos_a(i)		!x
			m_E(i,j,2)=v_pos_b(j)		!z
			CALL Exyz_per_sub(v_kinc,v_Einc,flaginc,nstop,&
			&ratio,lambda,v_pos_a(i),y,v_pos_b(j),v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
			& m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

			Exyz2_if: IF (error/=0) THEN
				error=1
				WRITE(*,50) 
			50 FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
			END IF Exyz2_if

		END DO z_xz_do
	END DO x_xz_do
	!OMP END PARALLEL DO

ELSE 

	!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static,1)
	x_xy_do: DO i=1,xstep
		y_xy_do: DO j=1,ystep

			m_E(i,j,1)=v_pos_a(i)		!x
			m_E(i,j,2)=v_pos_b(j)		!y
			CALL Exyz_per_sub(v_kinc,v_Einc,flaginc,nstop,&
			&ratio,lambda,v_pos_a(i),v_pos_b(j),z,v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
			& m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

			Exyz3_if: IF (error/=0) THEN
				error=1
				WRITE(*,60) 
			60 FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
			END IF Exyz3_if

		END DO y_xy_do
	END DO x_xy_do
	!OMP END PARALLEL DO

END IF xyz_if


final_if: IF (error/=0) THEN
	error=1
	WRITE(*,70) 
	70 FORMAT ("Errore in Exyz_sub chiamata da Efield_sub: il programma termina ora...")
	RETURN
END IF final_if

END SUBROUTINE Efield_per_sub





!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE Exyz_sub: dato un punto (x,y,z) nello spazio, calcolo il campo
! elettrico totale in coordinate cartesiane,sia che sia dentro sia che sia fuori
! una sfera
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Exyz_ext_sub(nstop,ratio,lambda,x,y,z,v_amnbmn,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,Ex,Ey,Ez,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz			!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn		!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn		!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn			!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq		!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index						!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt		!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio				!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: x,y,z							!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda							!Lunghezza d'onda incidente
COMPLEX(dbl), INTENT(OUT) :: Ex,Ey,Ez					!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Er,Etheta,Ephi						!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: rhoc,mc								!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor									!Rho reale
REAL(dbl) :: r,theta,phi 							!Coordinate sferiche locali
REAL(dbl) :: xl,yl,zl	 							!Coordinate cartesiane locali
REAL(dbl) :: d,radius								!Distanza punto-sfere,raggio sfera
INTEGER(lo) :: n,m,i,dimens,ns,low_b,up_b			!Indici e dimensione
REAL(dbl) :: nr,mr									!Indici reali

!Inizio subroutine vera e propria

!Calcolo il numero di sfere e dimens
dimens=nstop*(nstop+2)
ns=SIZE(m_xyz(:,1))

!--------------------------------------------
!Se sono qui allora non sono in nessuna sfera
!--------------------------------------------
Ex=(0.0D0,0.0D0)
Ey=(0.0D0,0.0D0)
Ez=(0.0D0,0.0D0)

ext_do: DO i=1,ns

	!Calcolo le coordinate sferiche locali
	xl=x-m_xyz(i,1)
	yl=y-m_xyz(i,2)
	zl=z-m_xyz(i,3)
	CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

	c_s_if1: IF (error/=0) THEN
		WRITE(*,30) 
		30 FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")                           
		RETURN
	END IF c_s_if1
	
	!Calcolo  rho per chiamare la subroutine per il calcolo del campo
	rhor=2*Pi_D*ref_index*r/lambda
	!Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
	low_b=2*dimens*(i-1)+1
	up_b=2*dimens*i
	CALL ext_field(v_amnbmn(low_b:up_b),v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)
	
	extfield_if: IF (error/=0) THEN
		WRITE(*,40) 
		40 FORMAT ("Errore in ext_field chiamata da Exyz_sub: il programma termina ora...")                           
		RETURN
	END IF extfield_if
	
	!Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
	Ex=Ex+Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
	Ey=Ey+Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
	Ez=Ez+Er*COS(theta) - Etheta*SIN(theta)

END DO ext_do

END SUBROUTINE Exyz_ext_sub


!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE Efield_sub: calcolo il campo su tutto un piano in coordinate  
! cartesiane
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_sub(k,phimin,phimax,phistep,thetastep,betap,v_amnbmn,m_xyz,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k						!Wavevector
REAL(dbl), INTENT(IN) :: phimin,phimax					!Angoli per il far field
REAL(dbl), INTENT(INOUT) :: betap						!Polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				!Intervalli di calcolo
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC			!Quantita' scattering
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: scatot					!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: s1,s2,s3,s4,sum1,sum2,sum3,sum4					!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,phase,dtheta,dphi					!Piani e angoli di scattering
REAL(dbl) :: nr,mr								!Indici reali
REAL(dbl) :: logw,fr								!logaritmo per cmn e frazione
REAL(dbl) :: s11,s12,dcpar,dcnorm,dc,cs,ss,dc1					!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,chunk			!Indici e dimensione
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_psimn, v_phimn, v_thetamn,	v_ximn	!Vettori funzioni angolari
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta					!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_cmn			!Vettori funzioni angolari
REAL(dbl), DIMENSION(thetastep+1,phistep+1,2) :: m_thetaphi			!Angle storage matrix

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_cmn(1:dimens))
ALLOCATE(v_psimn(1:side),v_phimn(1:side),v_thetamn(1:side),v_ximn(1:side))

!Calcolo i delta a seconda degli step
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=180.0D0/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0


!Filling the angular matrix
theta=0.0D0
fill_theta_do: DO itheta=1,thetastep+1

	phi=phimin

	fill_phi_do: DO iphi=1,phistep+1

		m_thetaphi(itheta,iphi,1)=theta
		m_thetaphi(itheta,iphi,2)=phi

		phi=phi+dphi

	END DO fill_phi_do

	theta=theta+dtheta

END DO fill_theta_do





!coefficienti di normalizzazione per le funzioni angolari
error=0
i=0
n_do_cmn: DO n=1,nstop

	nr=REAL(n,dbl)

	m_do_cmn: DO m=-n,n
	
		i=i+1
		mr=REAL(m,dbl)
		
		logw=lnf(nr-mr,error)-lnf(nr+mr,error)
		v_cmn(i)=SQRT(EXP(logw)*(2*nr+1.0D0)/(nr*(nr+1.0D0)))
	
	END DO m_do_cmn
END DO n_do_cmn


!Controllo di avere un piano ben definito
phi_if: IF (phimax<phimin) THEN
	error=1
	WRITE(*,10) 
	10 FORMAT ("Errore in Efar_sub: phi bounds non corretti, il programma termina ora...")
	RETURN
END IF phi_if


!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if


!Initializing the FIRSPRIVATE allocatable arrays for safety
v_pimn=0.0D0
v_taumn=0.0D0
v_psimn=0.0D0
v_phimn=0.0D0
v_thetamn=0.0D0
v_ximn=0.0D0

!Comincio i conti del campo e un grande if 
!$OMP PARALLEL PRIVATE(v_delta,theta,phi,betap,error,j,jm,j0,l,m,n,itheta,iphi,phase,s1,s2,s3,s4,sum1,sum2,sum3,sum4,mr,fr,s11,&
!$OMP s12,dcpar,dcnorm,cs,ss,v_pimn,v_taumn,v_psimn,v_phimn,v_thetamn,v_ximn)

!Choosing the chunk size
chunk=thetastep/10
chunk_if: IF (chunk==0) THEN
	chunk=1
END IF chunk_if

!$OMP DO SCHEDULE(guided,chunk)
theta_do: DO itheta=1,thetastep+1

	!Assigning theta
	theta=m_thetaphi(itheta,1,1)

	!Funzione angolare Pi_mn
	CALL pi_mnsca(nstop,theta,v_pimn,error)

	pimn_if: IF (error/=0) THEN
				WRITE(*,11) 
				11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF pimn_if

	!Funzione angolare Tau_mn
	CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)

	taumn_if: IF (error/=0) THEN
				WRITE(*,20) 
				20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF taumn_if

	!Normalizzo le funzioni angolari
	v_taumn=v_cmn*v_taumn
	v_pimn=v_cmn*v_pimn
		
	!Inizializzo il ciclo per i vettori di funzioni che mi servono nel calcolo di S
	j=0
	jm=0
	l_do: DO l=1,ns

		j0=0

		n_do: DO n=1,nstop
			m_do: DO m=-n,n 

				j0=n*(n+1)+m
				j=(l-1)*dimens+n*(n+1)+m	
				jm=(l-1)*dimens+n*(n+1)-m

				v_psimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) + &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))

				v_phimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) - &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))

				v_thetamn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) - &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))

				v_ximn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) + &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))

			END DO m_do
		END DO n_do
	END DO l_do


	phi_do: DO iphi=1,phistep+1

		!Assigning phi
		phi=m_thetaphi(itheta,iphi,2)

		!Riempio il vettore degli esponenziali
		delta_do: DO l=1,ns

			phase=m_xyz(l,1)*SIN(theta)*COS(phi) + m_xyz(l,2)*SIN(theta)*SIN(phi) + m_xyz(l,3)*COS(theta)
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)
			
			!If ((itheta==1) .AND. (iphi)==1) THEN
			!	WRITE(*,*) v_delta(l)
			!END IF
		
		END DO delta_do

		!Coordinate angolari
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi
		
		!Inizializzo il ciclo per le somme di S
		j=0
		jm=0
		s1=(0.0D0,0.0D0)
		s2=(0.0D0,0.0D0)
		s3=(0.0D0,0.0D0)
		s4=(0.0D0,0.0D0)
		l_do_phi: DO l=1,ns

			j0=0
			sum1=(0.0D0,0.0D0)
			sum2=(0.0D0,0.0D0)
			sum3=(0.0D0,0.0D0)
			sum4=(0.0D0,0.0D0)
			n_do_phi: DO n=1,nstop
				m_do_phi: DO m=0,n 

					j0=n*(n+1)+m
					j=(l-1)*dimens+n*(n+1)+m	
					jm=(l-1)*dimens+n*(n+1)-m
					mr=REAL(m,dbl)

					fr_if: IF (m==0) THEN 
						fr=0.5D0
					ELSE
						fr=1.0D0
					END IF fr_if

					cs=COS((mr-1.0D0)*phi+betap)
					ss=SIN((mr-1.0D0)*phi+betap)

					sum1=sum1+fr*(v_ximn(j)*cs + (0.0D0,1.0D0)*v_thetamn(j)*ss)
					sum2=sum2+fr*(v_psimn(j)*cs + (0.0D0,1.0D0)*v_phimn(j)*ss)
					sum3=sum3+fr*((0.0D0,1.0D0)*v_phimn(j)*cs - v_psimn(j)*ss)
					sum4=sum4+fr*(-(0.0D0,1.0D0)*v_thetamn(j)*cs + v_ximn(j)*ss)

				END DO m_do_phi
			END DO n_do_phi

			!Elementi dell'amplitude scattering matrix
			s1=s1+v_delta(l)*sum1
			s2=s2+v_delta(l)*sum2
			s3=s3+v_delta(l)*sum3
			s4=s4+v_delta(l)*sum4
			
		END DO l_do_phi

		!Elementi della matrice di muller
		s11=0.5D0*((ABS(s1)**2)+(ABS(s2)**2)+(ABS(s3)**2)+(ABS(s4)**2))
		s12=0.5D0*((ABS(s2)**2)-(ABS(s1)**2)+(ABS(s4)**2)-(ABS(s3)**2))

		!Sezione angolare differenziale parallela,normale e complessiva
		dcpar= (s11+s12)/(k**2)
		dcnorm=(s11-s12)/(k**2)
		m_SC(itheta,iphi,7)=((COS(phi-betap))**2)*dcpar + ((SIN(phi-betap))**2)*dcnorm


		!Coordinate angolari
		m_SC(itheta,iphi,3)=s1
		m_SC(itheta,iphi,4)=s2
		m_SC(itheta,iphi,5)=s3
		m_SC(itheta,iphi,6)=s4

	END DO phi_do
END DO theta_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!Integrating the scattering cross section
scatot=0.0D0
sca_theta_do: DO itheta=1,thetastep+1
	sca_phi_do: DO iphi=1,phistep+1

		scatot=scatot+m_SC(itheta,iphi,7)*SIN(m_thetaphi(itheta,iphi,1))*dtheta*dphi

	END DO sca_phi_do
END DO sca_theta_do

END SUBROUTINE Efar_sub








!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE Efield_sub_w: calcolo lo scattering e poi integro con gauss
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_sub_w(k,phistep,thetastep,betap,v_amnbmn,m_xyz,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k						!Wavevector
REAL(dbl), INTENT(INOUT) :: betap						!Polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				!Intervalli di calcolo
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC			!Quantita' scattering
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: scatot					!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: s1,s2,s3,s4,sum1,sum2,sum3,sum4					!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,phase,dtheta,dphi					!Piani e angoli di scattering
REAL(dbl) :: nr,mr								!Indici reali
REAL(dbl) :: logw,fr								!logaritmo per cmn e frazione
REAL(dbl) :: s11,s12,dcpar,dcnorm,dc,cs,ss,dc1					!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta			!Indici e dimensione
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_psimn, v_phimn, v_thetamn,	v_ximn	!Vettori funzioni angolari
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta					!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_cmn			!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_wphi,v_xphi,v_wtheta,v_xtheta	!Vettori peso

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_cmn(1:dimens))
ALLOCATE(v_psimn(1:side),v_phimn(1:side),v_thetamn(1:side),v_ximn(1:side))
ALLOCATE(v_wphi(1:phistep),v_xphi(1:phistep),v_wtheta(1:thetastep),v_xtheta(1:thetastep))

!Alloco e calcolo i vettori per contenere gli angoli phi e theta e per i successivi pesi
CALL weigths1(phistep,thetastep,v_wphi,v_xphi,v_wtheta,v_xtheta,error)

!coefficienti di normalizzazione per le funzioni angolari
error=0
i=0
n_do_cmn: DO n=1,nstop

	nr=REAL(n,dbl)

	m_do_cmn: DO m=-n,n
	
		i=i+1
		mr=REAL(m,dbl)
		
		logw=lnf(nr-mr,error)-lnf(nr+mr,error)
		v_cmn(i)=SQRT(EXP(logw)*(2*nr+1.0D0)/(nr*(nr+1.0D0)))
	
	END DO m_do_cmn
END DO n_do_cmn


!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if

scatot=0.0D0

theta_do: DO itheta=1,thetastep

	!Funzione angolare Pi_mn
	CALL pi_mnsca(nstop,v_xtheta(itheta),v_pimn,error)
	
	pimn_if: IF (error/=0) THEN																	
				WRITE(*,11) 
				11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF pimn_if
	
	!Funzione angolare Tau_mn
	CALL tau_mnsca(nstop,v_xtheta(itheta),v_pimn,v_taumn,error)
	
	taumn_if: IF (error/=0) THEN																	
				WRITE(*,20) 
				20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF taumn_if

	!Normalizzo le funzioni angolari
	v_taumn=v_cmn*v_taumn
	v_pimn=v_cmn*v_pimn
		
	!Inizializzo il ciclo per i vettori di funzioni che mi servono nel calcolo di S
	j=0
	jm=0
	l_do: DO l=1,ns
	
		j0=0
	
		n_do: DO n=1,nstop
			m_do: DO m=-n,n 

				j0=n*(n+1)+m
				j=(l-1)*dimens+n*(n+1)+m	
				jm=(l-1)*dimens+n*(n+1)-m
				
				v_psimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) + &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))
				
				v_phimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) - &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))
			
				v_thetamn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) - &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))
			
				v_ximn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) + &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))
			
			END DO m_do
		END DO n_do
	END DO l_do

	

	phi_do: DO iphi=1,phistep
	
		!Riempio il vettore degli esponenziali
		delta_do: DO l=1,ns
		
			phase=m_xyz(l,1)*SIN(v_xtheta(itheta))*COS(v_xphi(iphi)) + m_xyz(l,2)*SIN(v_xtheta(itheta))*SIN(v_xphi(iphi)) + &
			& m_xyz(l,3)*COS(v_xtheta(itheta))
			
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)
		
		END DO delta_do
	
		!Coordinate angolari
		m_SC(itheta,iphi,1)=v_xtheta(itheta)
		m_SC(itheta,iphi,2)=v_xphi(iphi)
		
		!Inizializzo il ciclo per le somme di S
		j=0
		jm=0
		s1=(0.0D0,0.0D0)
		s2=(0.0D0,0.0D0)
		s3=(0.0D0,0.0D0)
		s4=(0.0D0,0.0D0)
		l_do_phi: DO l=1,ns
		
			j0=0
			sum1=(0.0D0,0.0D0)
			sum2=(0.0D0,0.0D0)
			sum3=(0.0D0,0.0D0)
			sum4=(0.0D0,0.0D0)
			n_do_phi: DO n=1,nstop
				m_do_phi: DO m=0,n 
	
					j0=n*(n+1)+m
					j=(l-1)*dimens+n*(n+1)+m	
					jm=(l-1)*dimens+n*(n+1)-m
					mr=REAL(m,dbl)
					
					fr_if: IF (m==0) THEN 
						fr=0.5D0
					ELSE
						fr=1.0D0
					END IF fr_if
					
					cs=COS((mr-1.0D0)*v_xphi(iphi)+betap)
					ss=SIN((mr-1.0D0)*v_xphi(iphi)+betap)
					
					sum1=sum1+fr*(v_ximn(j)*cs + (0.0D0,1.0D0)*v_thetamn(j)*ss)
					sum2=sum2+fr*(v_psimn(j)*cs + (0.0D0,1.0D0)*v_phimn(j)*ss)
					sum3=sum3+fr*((0.0D0,1.0D0)*v_phimn(j)*cs - v_psimn(j)*ss)
					sum4=sum4+fr*(-(0.0D0,1.0D0)*v_thetamn(j)*cs + v_ximn(j)*ss)
				
				END DO m_do_phi
			END DO n_do_phi
			
			!Elementi dell'amplitude scattering matrix
			s1=s1+v_delta(l)*sum1
			s2=s2+v_delta(l)*sum2
			s3=s3+v_delta(l)*sum3
			s4=s4+v_delta(l)*sum4
			
		END DO l_do_phi
		
		!Elementi della matrice di muller
		s11=0.5D0*((ABS(s1)**2)+(ABS(s2)**2)+(ABS(s3)**2)+(ABS(s4)**2))
		s12=0.5D0*((ABS(s2)**2)-(ABS(s1)**2)+(ABS(s4)**2)-(ABS(s3)**2))
		
		!Sezione angolare differenziale parallela,normale e complessiva
		dcpar= (s11+s12)/(k**2)
		dcnorm=(s11-s12)/(k**2)
		dc=((COS(v_xphi(iphi)-betap))**2)*dcpar + ((SIN(v_xphi(iphi)-betap))**2)*dcnorm
		
		! scatot=scatot+dc*SIN(theta)*dtheta*dphi
		scatot=scatot+dc*v_wphi(iphi)*v_wtheta(itheta)
				
		!Coordinate angolari
		m_SC(itheta,iphi,3)=s1
		m_SC(itheta,iphi,4)=s2
		m_SC(itheta,iphi,5)=s3
		m_SC(itheta,iphi,6)=s4
		m_SC(itheta,iphi,7)=dc

	END DO phi_do
	
END DO theta_do

END SUBROUTINE Efar_sub_w
















!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE Efar_sub_poynting: calcolo lo scattering e poi integro con gauss
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_sub_poynting(k,phimin,phimax,phistep,thetastep,betap,v_amnbmn,m_xyz,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k							!Wavevector
REAL(dbl), INTENT(IN) :: phimin,phimax					!Angoli per il far field
REAL(dbl), INTENT(INOUT) :: betap						!Polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				!Intervalli di calcolo
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC			!Quantita' scattering
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: scatot						!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,Htheta,Hphi,sum1,sum2,sum3,sum4,dccomplex					!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,phase,dtheta,dphi								!Piani e angoli di scattering
REAL(dbl) :: dc												!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,chunk						!Indici e dimensione
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta								!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn							!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_amn_taumn,v_amn_pimn,v_bmn_taumn,v_bmn_pimn,v_emn,v_ab_tot	!Vettori conti campo
REAL(dbl), DIMENSION(thetastep+1,phistep+1,2) :: m_thetaphi						!Angle storage matrix

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_emn(1:dimens),v_ab_tot(1:2*dimens))
ALLOCATE(v_amn_pimn(1:dimens),v_amn_taumn(1:dimens),v_bmn_pimn(1:dimens),v_bmn_taumn(1:dimens))

!Calcolo i delta a seconda degli step
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=180.0D0/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0


!Filling the angular matrix
theta=0.0D0
fill_theta_do: DO itheta=1,thetastep+1

	phi=phimin

	fill_phi_do: DO iphi=1,phistep+1

		m_thetaphi(itheta,iphi,1)=theta
		m_thetaphi(itheta,iphi,2)=phi

		phi=phi+dphi

	END DO fill_phi_do

	theta=theta+dtheta

END DO fill_theta_do


!coefficienti di normalizzazione per le funzioni angolari
error=0
CALL emn_sub(nstop,v_emn,error)

!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if

scatot=0.0D0

!Comincio i conti del campo e un grande if 

!Initializing the FIRSPRIVATE allocatable arrays for safety
v_pimn=0.0D0
v_taumn=0.0D0
v_ab_tot=0.0D0
v_amn_pimn=0.0D0
v_amn_taumn=0.0D0
v_bmn_pimn=0.0D0
v_bmn_taumn=0.0D0

!Comincio i conti del campo e un grande if 
!$OMP PARALLEL PRIVATE(phase,v_delta,theta,phi,error,l,m,n,j0,j,itheta,iphi,Etheta,Ephi,Htheta,Hphi,sum1,sum2,sum3,sum4,dccomplex,&
!$OMP dc,v_pimn,v_taumn,v_ab_tot,v_amn_pimn,v_amn_taumn,v_bmn_pimn,v_bmn_taumn)

!Choosing the chunk size
chunk=thetastep/10
chunk_if: IF (chunk==0) THEN
	chunk=1
END IF chunk_if

!$OMP DO SCHEDULE(guided,chunk)
theta_do: DO itheta=1,thetastep+1

	theta=m_thetaphi(itheta,1,1)

	!Funzione angolare Pi_mn
	CALL pi_mnsca(nstop,theta,v_pimn,error)
	
	pimn_if: IF (error/=0) THEN
		WRITE(*,11) 
		11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
		STOP
	END IF pimn_if
	
	!Funzione angolare Tau_mn
	CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)
	
	taumn_if: IF (error/=0) THEN
		WRITE(*,20) 
		20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
		STOP
	END IF taumn_if


	phi_do: DO iphi=1,phistep+1

		phi=m_thetaphi(itheta,iphi,2)

		!Riempio il vettore degli esponenziali
		delta_do: DO l=1,ns

			phase=m_xyz(l,1)*SIN(theta)*COS(phi) + m_xyz(l,2)*SIN(theta)*SIN(phi) + m_xyz(l,3)*COS(theta)
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)

		END DO delta_do

		!Coordinate angolari
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi

		!Qui calcolo i coefficienti totali del campo,questo approccio mi dovrebbe dare il campo totale nel far field
		!e quindi dovrei calcolare tutto unitariamente e senza problemi!

		!inizializzo i coefficienti totali
		v_ab_tot=0.0D0

		tot_sca_ndo: DO n=1,nstop
			tot_sca_mdo: DO m=-n,n
				tot_sca_ldo: DO l=1,ns

					j0=n*(n+1)+m
					j=(l-1)*dimens+n*(n+1)+m
					v_ab_tot(2*j0-1)=v_ab_tot(2*j0-1)+v_delta(l)*v_amnbmn(2*j-1)
					v_ab_tot(2*j0)=v_ab_tot(2*j0)+v_delta(l)*v_amnbmn(2*j)

				END DO tot_sca_ldo
			END DO tot_sca_mdo
		END DO tot_sca_ndo

	
		!Qui calcolo i vettori misti taumn pimn amn bmn per il calcolo dei campi

		v_amn_taumn=v_ab_tot(1:2*dimens:2)*v_taumn
		v_bmn_taumn=v_ab_tot(2:2*dimens:2)*v_taumn
		v_amn_pimn=v_ab_tot(1:2*dimens:2)*v_pimn
		v_bmn_pimn=v_ab_tot(2:2*dimens:2)*v_pimn

		!Inizializzo il ciclo per le somme di S
		Etheta=(0.0D0,0.0D0)
		Ephi=(0.0D0,0.0D0)
		Htheta=(0.0D0,0.0D0)
		Hphi=(0.0D0,0.0D0)

		j0=0
		n_do_phi: DO n=1,nstop

			sum1=(0.0D0,0.0D0)
			sum2=(0.0D0,0.0D0)
			sum3=(0.0D0,0.0D0)
			sum4=(0.0D0,0.0D0)

			m_do_phi: DO m=-n,n 

				j0=n*(n+1)+m

				sum1=sum1+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum2=sum2+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum3=sum3+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum4=sum4+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)

			END DO m_do_phi

			Etheta=Etheta+((0.0D0,-1.0D0)**(n+1))*sum1
			Ephi=Ephi+((0.0D0,-1.0D0)**n)*sum2
			Htheta=Htheta+((0.0D0,-1.0D0)**n)*sum3
			Hphi=Hphi+((0.0D0,-1.0D0)**(n+1))*sum4

		END DO n_do_phi

		!Correggo dove occore per il segno
		Etheta=-Etheta
		Ephi=-Ephi
		Hphi=-Hphi

		!Sezione angolare differenziale complessa e reale
		dccomplex=(Etheta*CONJG(Hphi)-Ephi*CONJG(Htheta))/(k**2)
		dc=REAL(dccomplex,dbl)

		!Coordinate angolari
		m_SC(itheta,iphi,3)=REAL(Etheta,dbl)
		m_SC(itheta,iphi,4)=REAL(Ephi,dbl)
		m_SC(itheta,iphi,5)=REAL(Htheta,dbl)
		m_SC(itheta,iphi,6)=REAL(Hphi,dbl)
		m_SC(itheta,iphi,7)=dc

	END DO phi_do
END DO theta_do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!Integrating the scattering cross section
scatot=0.0D0
sca_theta_do: DO itheta=1,thetastep+1
	sca_phi_do: DO iphi=1,phistep+1

		scatot=scatot+m_SC(itheta,iphi,7)*SIN(m_thetaphi(itheta,iphi,1))*dtheta*dphi

	END DO sca_phi_do
END DO sca_theta_do

END SUBROUTINE Efar_sub_poynting




































!******************************************************************************
!******************************************************************************
!******************************************************************************
!5) SUBROUTINE Efar_sub_poynting: calcolo lo scattering e poi integro con gauss
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_sub_poynting_per(k,phimin,phimax,phistep,thetamin,thetamax,thetastep,betap,v_amnbmn,m_xyz,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k							!Wavevector
REAL(dbl), INTENT(IN) :: phimin,phimax					!Angoli per il far field
REAL(dbl), INTENT(IN) :: thetamin,thetamax				!Angoli per il far field
REAL(dbl), INTENT(INOUT) :: betap						!Polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				!Intervalli di calcolo
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC			!Quantita' scattering
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: scatot						!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,Htheta,Hphi,sum1,sum2,sum3,sum4,dccomplex		!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,phase,dtheta,dphi							!Piani e angoli di scattering
COMPLEX(dbl) :: qa,qb,suma,sumb								!Basi delle somme esponenziali
REAL(dbl) :: Sa,Sb										!Somme parziali
REAL(dbl) :: dc,x,y,z										!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,na,nb			!Indici e dimensione
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta							!Vettori funzioni angolari
COMPLEX(dbl), DIMENSION(0:(nacell-1),0:(nbcell-1),1:ns) :: v_delta_nanb		!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn				!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_amn_taumn,v_amn_pimn,v_bmn_taumn,v_bmn_pimn,v_emn,v_ab_tot	!Vettori conti campo

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_emn(1:dimens),v_ab_tot(1:2*dimens))
ALLOCATE(v_amn_pimn(1:dimens),v_amn_taumn(1:dimens),v_bmn_pimn(1:dimens),v_bmn_taumn(1:dimens))

!Calcolo i delta a seconda degli step
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=(thetamax-thetamin)/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0

!coefficienti di normalizzazione per le funzioni angolari
error=0
CALL emn_sub(nstop,v_emn,error)

!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")
	RETURN
END IF steps_if

scatot=0.0D0

!Comincio i conti del campo e un grande if 
theta=0.0D0
phi=phimin

theta_do: DO itheta=1,thetastep+1

	!Funzione angolare Pi_mn
	CALL pi_mnsca(nstop,theta,v_pimn,error)
	
	pimn_if: IF (error/=0) THEN	
				WRITE(*,11) 
				11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF pimn_if
	
	!Funzione angolare Tau_mn
	CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)
	
	taumn_if: IF (error/=0) THEN																	
				WRITE(*,20) 
				20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF taumn_if


	phi=phimin


	phi_do: DO iphi=1,phistep+1

		!Calcolo le somme compatte per i fattori di fase,aiutano a fare velocemente i calcoli
		Sa=ax*SIN(theta)*COS(phi)
		Sb=bx*SIN(theta)*COS(phi)+by*SIN(theta)*SIN(phi)
		qa=EXP(CMPLX(0.0D0,-k*Sa))
		qb=EXP(CMPLX(0.0D0,-k*Sb))

		suma_if: IF (Sa==0.0D0) THEN
			suma=CMPLX(REAL(nacell,dbl),0.0D0,KIND=dbl)
		ELSE
			suma=((qa**nacell)-(1.0D0,0.0D0))/(qa-(1.0D0,0.0D0))
		END IF suma_if

		sumb_if: IF (Sb==0.0D0) THEN
			sumb=CMPLX(REAL(nbcell,dbl),0.0D0,KIND=dbl)
		ELSE
			sumb=((qb**nbcell)-(1.0D0,0.0D0))/(qb-(1.0D0,0.0D0))
		END IF sumb_if

!		WRITE(*,*) "theta,phi",theta,phi
!		WRITE(*,*) "qa: ",qa
!		WRITE(*,*) "qb: ",qb
!		WRITE(*,*) "sa: ",sa
!		WRITE(*,*) "sb: ",sb
!		WRITE(*,*) "suma: ",suma
!		WRITE(*,*) "sumb: ",sumb
!		WRITE(*,*)

		!Riempio il vettore degli esponenziali
		delta_do: DO l=1,ns
		
			phase=m_xyz(l,1)*SIN(theta)*COS(phi) + m_xyz(l,2)*SIN(theta)*SIN(phi) + &
			& m_xyz(l,3)*COS(theta)
			
!			phase=m_xyz(l,3)*COS(theta)
			
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)
		
		END DO delta_do


!		!Riempio il vettore degli esponenziali
!		delta_do: DO l=1,ns
!			delta_na_do: DO na=0,nacell-1
!				delta_nb_do: DO nb=0,nbcell-1

!					x=m_xyz(l,1)+na*ax+nb*bx
!					y=m_xyz(l,2)+nb*by
!					z=m_xyz(l,3)

!					phase=x*SIN(theta)*COS(phi) + y*SIN(theta)*SIN(phi) + z*COS(theta)

!					v_delta_nanb(na,nb,l)=EXP(-(0.0D0,1.0D0)*k*phase)

!				END DO delta_nb_do
!			END DO delta_na_do
!		END DO delta_do


		!Coordinate angolari
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi
		
		!Qui calcolo i coefficienti totali del campo,questo approccio mi dovrebbe dare il campo totale nel far field
		!e quindi dovrei calcolare tutto unitariamente e senza problemi!
		
		!inizializzo i coefficienti totali
		v_ab_tot=0.0D0

		tot_sca_ndo: DO n=1,nstop
			tot_sca_mdo: DO m=-n,n
				tot_sca_ldo: DO l=1,ns

					j0=n*(n+1)+m
					j=(l-1)*dimens+n*(n+1)+m
					v_ab_tot(2*j0-1)=v_ab_tot(2*j0-1)+v_delta(l)*v_amnbmn(2*j-1)
					v_ab_tot(2*j0)=v_ab_tot(2*j0)+v_delta(l)*v_amnbmn(2*j)

				END DO tot_sca_ldo
			END DO tot_sca_mdo
		END DO tot_sca_ndo


!		tot_sca_ndo: DO n=1,nstop
!			tot_sca_mdo: DO m=-n,n
!				tot_sca_ldo: DO l=1,ns
!					tot_sca_nado: DO na=0,nacell-1
!						tot_sca_nbdo: DO nb=0,nbcell-1

!							j0=n*(n+1)+m
!							j=(l-1)*dimens+n*(n+1)+m
!							v_ab_tot(2*j0-1)=v_ab_tot(2*j0-1)+v_delta_nanb(na,nb,l)*v_amnbmn(2*j-1)
!							v_ab_tot(2*j0)=v_ab_tot(2*j0)+v_delta_nanb(na,nb,l)*v_amnbmn(2*j)

!						END DO tot_sca_nbdo
!					END DO tot_sca_nado
!				END DO tot_sca_ldo
!			END DO tot_sca_mdo
!		END DO tot_sca_ndo


		v_ab_tot=v_ab_tot*suma*sumb

		!Qui calcolo i vettori misti taumn pimn amn bmn per il calcolo dei campi
		
		v_amn_taumn=v_ab_tot(1:2*dimens:2)*v_taumn
		v_bmn_taumn=v_ab_tot(2:2*dimens:2)*v_taumn
		v_amn_pimn=v_ab_tot(1:2*dimens:2)*v_pimn
		v_bmn_pimn=v_ab_tot(2:2*dimens:2)*v_pimn
		
		!Inizializzo il ciclo per le somme di S
		Etheta=(0.0D0,0.0D0)
		Ephi=(0.0D0,0.0D0)
		Htheta=(0.0D0,0.0D0)
		Hphi=(0.0D0,0.0D0)
			
		j0=0
		n_do_phi: DO n=1,nstop
		
			sum1=(0.0D0,0.0D0)
			sum2=(0.0D0,0.0D0)
			sum3=(0.0D0,0.0D0)
			sum4=(0.0D0,0.0D0)
			
			m_do_phi: DO m=-n,n 
				
				j0=n*(n+1)+m
				
				sum1=sum1+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum2=sum2+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum3=sum3+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum4=sum4+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
			
			END DO m_do_phi
			
			Etheta=Etheta+((0.0D0,-1.0D0)**(n+1))*sum1
			Ephi=Ephi+((0.0D0,-1.0D0)**n)*sum2
			Htheta=Htheta+((0.0D0,-1.0D0)**n)*sum3
			Hphi=Hphi+((0.0D0,-1.0D0)**(n+1))*sum4
			
		END DO n_do_phi
		
		!Correggo dove occore per il segno
		Etheta=-Etheta
		Ephi=-Ephi
		Hphi=-Hphi
				
		!Sezione angolare differenziale complessa e reale
		dccomplex=(Etheta*CONJG(Hphi)-Ephi*CONJG(Htheta))/(k**2)
		dc=REAL(dccomplex,dbl)
		scatot=scatot+dc*SIN(theta)*dtheta*dphi
!		scatot=scatot+dc*v_wphi(iphi)*v_wtheta(itheta)
				
		!Coordinate angolari
		m_SC(itheta,iphi,3)=REAL(Etheta,dbl)
		m_SC(itheta,iphi,4)=REAL(Ephi,dbl)
		m_SC(itheta,iphi,5)=REAL(Htheta,dbl)
		m_SC(itheta,iphi,6)=REAL(Hphi,dbl)
		m_SC(itheta,iphi,7)=dc
		
		phi=phi+dphi

	END DO phi_do
	
	theta=theta+dtheta
	
END DO theta_do

END SUBROUTINE Efar_sub_poynting_per













!******************************************************************************
!******************************************************************************
!******************************************************************************
!5) SUBROUTINE Efar_sub_poynting: calcolo lo scattering e poi integro con gauss
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Rad_sub_poynting_per(k,phiin,thetain,v_amnbmn,m_xyz,rad,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k							!Wavevector
REAL(dbl), INTENT(IN) :: phiin							!Angoli per il far field
REAL(dbl), INTENT(IN) :: thetain						!Angoli per il far field
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: rad							!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,Htheta,Hphi,sum1,sum2,sum3,sum4,dccomplex		!Somme parziali e finali per gli elementi S
REAL(dbl) :: phi,theta,phase,dtheta,dphi							!Piani e angoli di scattering
COMPLEX(dbl) :: qa,qb,suma,sumb								!Basi delle somme esponenziali
REAL(dbl) :: Sa,Sb										!Somme parziali
REAL(dbl) :: dc,x,y,z										!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,na,nb			!Indici e dimensione
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta							!Vettori funzioni angolari
COMPLEX(dbl), DIMENSION(0:(nacell-1),0:(nbcell-1),1:ns) :: v_delta_nanb		!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn				!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_amn_taumn,v_amn_pimn,v_bmn_taumn,v_bmn_pimn,v_emn,v_ab_tot	!Vettori conti campo

!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_emn(1:dimens),v_ab_tot(1:2*dimens))
ALLOCATE(v_amn_pimn(1:dimens),v_amn_taumn(1:dimens),v_bmn_pimn(1:dimens),v_bmn_taumn(1:dimens))


!coefficienti di normalizzazione per le funzioni angolari
error=0
CALL emn_sub(nstop,v_emn,error)

theta=thetain*Pi_d/180.0D0
phi=phiin*Pi_d/180.0D0

rad=0.0D0

!Funzione angolare Pi_mn
CALL pi_mnsca(nstop,theta,v_pimn,error)

pimn_if: IF (error/=0) THEN	
			WRITE(*,11) 
			11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
			STOP
END IF pimn_if

!Funzione angolare Tau_mn
CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)

taumn_if: IF (error/=0) THEN																	
			WRITE(*,20) 
			20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
			STOP
END IF taumn_if

!Calcolo le somme compatte per i fattori di fase,aiutano a fare velocemente i calcoli
Sa=ax*SIN(theta)*COS(phi)
Sb=bx*SIN(theta)*COS(phi)+by*SIN(theta)*SIN(phi)
qa=EXP(CMPLX(0.0D0,-k*Sa))
qb=EXP(CMPLX(0.0D0,-k*Sb))

suma_if: IF (Sa==0.0D0) THEN
	suma=CMPLX(REAL(nacell,dbl),0.0D0,KIND=dbl)
ELSE
	suma=((qa**nacell)-(1.0D0,0.0D0))/(qa-(1.0D0,0.0D0))
END IF suma_if

sumb_if: IF (Sb==0.0D0) THEN
	sumb=CMPLX(REAL(nbcell,dbl),0.0D0,KIND=dbl)
ELSE
	sumb=((qb**nbcell)-(1.0D0,0.0D0))/(qb-(1.0D0,0.0D0))
END IF sumb_if



!Riempio il vettore degli esponenziali
delta_do: DO l=1,ns

	phase=m_xyz(l,1)*SIN(theta)*COS(phi) + m_xyz(l,2)*SIN(theta)*SIN(phi) + &
	& m_xyz(l,3)*COS(theta)

	v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)

END DO delta_do


!Qui calcolo i coefficienti totali del campo,questo approccio mi dovrebbe dare il campo totale nel far field
!e quindi dovrei calcolare tutto unitariamente e senza problemi!

!inizializzo i coefficienti totali
v_ab_tot=0.0D0

tot_sca_ndo: DO n=1,nstop
	tot_sca_mdo: DO m=-n,n
		tot_sca_ldo: DO l=1,ns

			j0=n*(n+1)+m
			j=(l-1)*dimens+n*(n+1)+m
			v_ab_tot(2*j0-1)=v_ab_tot(2*j0-1)+v_delta(l)*v_amnbmn(2*j-1)
			v_ab_tot(2*j0)=v_ab_tot(2*j0)+v_delta(l)*v_amnbmn(2*j)

		END DO tot_sca_ldo
	END DO tot_sca_mdo
END DO tot_sca_ndo


v_ab_tot=v_ab_tot*suma*sumb

!Qui calcolo i vettori misti taumn pimn amn bmn per il calcolo dei campi

v_amn_taumn=v_ab_tot(1:2*dimens:2)*v_taumn
v_bmn_taumn=v_ab_tot(2:2*dimens:2)*v_taumn
v_amn_pimn=v_ab_tot(1:2*dimens:2)*v_pimn
v_bmn_pimn=v_ab_tot(2:2*dimens:2)*v_pimn

!Inizializzo il ciclo per le somme di S
Etheta=(0.0D0,0.0D0)
Ephi=(0.0D0,0.0D0)
Htheta=(0.0D0,0.0D0)
Hphi=(0.0D0,0.0D0)
	
j0=0
n_do_phi: DO n=1,nstop

	sum1=(0.0D0,0.0D0)
	sum2=(0.0D0,0.0D0)
	sum3=(0.0D0,0.0D0)
	sum4=(0.0D0,0.0D0)
	
	m_do_phi: DO m=-n,n 

		j0=n*(n+1)+m
		
		sum1=sum1+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
		sum2=sum2+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
		sum3=sum3+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
		sum4=sum4+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
	
	END DO m_do_phi
	
	Etheta=Etheta+((0.0D0,-1.0D0)**(n+1))*sum1
	Ephi=Ephi+((0.0D0,-1.0D0)**n)*sum2
	Htheta=Htheta+((0.0D0,-1.0D0)**n)*sum3
	Hphi=Hphi+((0.0D0,-1.0D0)**(n+1))*sum4

END DO n_do_phi

!Correggo dove occore per il segno
Etheta=-Etheta
Ephi=-Ephi
Hphi=-Hphi

!Sezione angolare differenziale complessa e reale
dccomplex=(Etheta*CONJG(Hphi)-Ephi*CONJG(Htheta))/(k**2)
rad=REAL(dccomplex,dbl)

END SUBROUTINE Rad_sub_poynting_per







!******************************************************************************
!******************************************************************************
!******************************************************************************
!6) SUBROUTINE Efield_sub_per: calcolo il campo su tutto un piano in coordinate  
! cartesiane
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_sub_per(k,phimin,phimax,phistep,thetastep,betap,v_amnbmn,m_xyz,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k								!Wavevector
REAL(dbl), INTENT(IN) :: phimin,phimax					!Angoli per il far field
REAL(dbl), INTENT(INOUT) :: betap							!Polarizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz			!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn	!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep			!Intervalli di calcolo
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC	!Quantita' scattering
REAL(dbl), INTENT(OUT) :: scatot						!Sezione di scattering totale
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: s1,s2,s3,s4,sum1,sum2,sum3,sum4								!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,phase,dtheta,dphi											!Piani e angoli di scattering
REAL(dbl) :: nr,mr															!Indici reali
REAL(dbl) :: logw,fr															!logaritmo per cmn e frazione
REAL(dbl) :: x,y,z																!Coordinate traslate sfere
REAL(dbl) :: s11,s12,dcpar,dcnorm,dc,cs,ss										!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,na,nb						!Indici e dimensione
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_psimn, v_phimn, v_thetamn,	v_ximn	!Vettori funzioni angolari
COMPLEX(dbl), DIMENSION(1:(ns*nacell*nbcell)) :: v_delta							!Vettori funzioni angolari															!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_cmn					!Vettori funzioni angolari


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_cmn(1:dimens))
ALLOCATE(v_psimn(1:side),v_phimn(1:side),v_thetamn(1:side),v_ximn(1:side))

!Calcolo i delta a seconda degli step
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=180.0D0/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0

!coefficienti di normalizzazione per le funzioni angolari
error=0
i=0
n_do_cmn: DO n=1,nstop

	nr=REAL(n,dbl)

	m_do_cmn: DO m=-n,n
	
		i=i+1
		mr=REAL(m,dbl)
		
		logw=lnf(nr-mr,error)-lnf(nr+mr,error)
		v_cmn(i)=SQRT(EXP(logw)*(2*nr+1.0D0)/(nr*(nr+1.0D0)))
	
	END DO m_do_cmn
END DO n_do_cmn


!Controllo di avere un piano ben definito
phi_if: IF (phimax<phimin) THEN
	error=1
	WRITE(*,10) 
	10 FORMAT ("Errore in Efar_sub: phi bounds non corretti, il programma termina ora...")
	RETURN
END IF phi_if


!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if


!Comincio i conti del campo e un grande if 
theta=0.0D0
phi=phimin

scatot=0.0D0

theta_do: DO itheta=1,thetastep+1

	!Funzione angolare Pi_mn
	CALL pi_mnsca(nstop,theta,v_pimn,error)
	
	pimn_if: IF (error/=0) THEN																	
				WRITE(*,11) 
				11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF pimn_if
	
	!Funzione angolare Tau_mn
	CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)
	
	taumn_if: IF (error/=0) THEN																	
				WRITE(*,20) 
				20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
				STOP
	END IF taumn_if

	!Normalizzo le funzioni angolari
	v_taumn=v_cmn*v_taumn
	v_pimn=v_cmn*v_pimn
		
	!Inizializzo il ciclo per i vettori di funzioni che mi servono nel calcolo di S
	j=0
	jm=0
	l_do: DO l=1,ns
	
		j0=0
	
		n_do: DO n=1,nstop
			m_do: DO m=-n,n 

				j0=n*(n+1)+m
				j=(l-1)*dimens+n*(n+1)+m	
				jm=(l-1)*dimens+n*(n+1)-m
				
				v_psimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) + &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))
				
				v_phimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) - &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))
			
				v_thetamn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) - &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))
			
				v_ximn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) + &
				& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))
			
			END DO m_do
		END DO n_do
	END DO l_do

	phi=phimin

	phi_do: DO iphi=1,phistep+1
	
		!Riempio il vettore degli esponenziali
		j=0
		delta_do_l: DO l=1,ns
			delta_do_na: DO na=0,nacell-1
				delta_do_nb: DO nb=0,nbcell-1
				
					j=j+1
				
					x=m_xyz(l,1)+ax*REAL(na,dbl)+bx*REAL(nb,dbl)
					y=m_xyz(l,2)+by*REAL(nb,dbl)
					z=m_xyz(l,3)
				
					phase=x*SIN(theta)*COS(phi) + y*SIN(theta)*SIN(phi) + z*COS(theta)
					v_delta(j)=EXP(-(0.0D0,1.0D0)*k*phase)
					
					!If ((itheta==1) .AND. (iphi)==1) THEN
					!	WRITE(*,*) v_delta(j)
					!END IF
				
				END DO delta_do_nb
			END DO delta_do_na
		END DO delta_do_l
	
		!Coordinate angolari
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi
		
		!Inizializzo il ciclo per le somme di S
		j0=0
		jm=0
		s1=(0.0D0,0.0D0)
		s2=(0.0D0,0.0D0)
		s3=(0.0D0,0.0D0)
		s4=(0.0D0,0.0D0)
		l_do_phi: DO l=1,ns
			na_do_phi: DO na=0,nacell-1
				nb_do_phi: DO nb=0,nbcell-1
					
					j0=j0+1
					
					sum1=(0.0D0,0.0D0)
					sum2=(0.0D0,0.0D0)
					sum3=(0.0D0,0.0D0)
					sum4=(0.0D0,0.0D0)
					n_do_phi: DO n=1,nstop
						m_do_phi: DO m=0,n 
							
							j=(l-1)*dimens+n*(n+1)+m	
							mr=REAL(m,dbl)
							
							fr_if: IF (m==0) THEN 
								fr=0.5D0
							ELSE
								fr=1.0D0
							END IF fr_if
							
							cs=COS((mr-1.0D0)*phi+betap)
							ss=SIN((mr-1.0D0)*phi+betap)
							
							sum1=sum1+fr*(v_ximn(j)*cs + (0.0D0,1.0D0)*v_thetamn(j)*ss)
							sum2=sum2+fr*(v_psimn(j)*cs + (0.0D0,1.0D0)*v_phimn(j)*ss)
							sum3=sum3+fr*((0.0D0,1.0D0)*v_phimn(j)*cs - v_psimn(j)*ss)
							sum4=sum4+fr*(-(0.0D0,1.0D0)*v_thetamn(j)*cs + v_ximn(j)*ss)
							
						END DO m_do_phi
					END DO n_do_phi
					
					!Elementi dell'amplitude scattering matrix
					s1=s1+v_delta(j0)*sum1
					s2=s2+v_delta(j0)*sum2
					s3=s3+v_delta(j0)*sum3
					s4=s4+v_delta(j0)*sum4
					
				END DO nb_do_phi
			END DO na_do_phi
		END DO l_do_phi
		
		!Elementi della matrice di muller
		s11=0.5D0*((ABS(s1)**2)+(ABS(s2)**2)+(ABS(s3)**2)+(ABS(s4)**2))
		s12=0.5D0*((ABS(s2)**2)-(ABS(s1)**2)+(ABS(s4)**2)-(ABS(s3)**2))
		
		!Sezione angolare differenziale parallela,normale e complessiva
		dcpar= (s11+s12)/(k**2)
		dcnorm=(s11-s12)/(k**2)
		dc=((COS(phi-betap))**2)*dcpar + ((SIN(phi-betap))**2)*dcnorm
		scatot=scatot+dc*SIN(theta)*dtheta*dphi
				
		!Coordinate angolari
		m_SC(itheta,iphi,3)=s1
		m_SC(itheta,iphi,4)=s2
		m_SC(itheta,iphi,5)=s3
		m_SC(itheta,iphi,6)=s4
		m_SC(itheta,iphi,7)=dc
		
		phi=phi+dphi

	END DO phi_do
	
	theta=theta+dtheta
	
END DO theta_do

scatot=scatot/nacell

END SUBROUTINE Efar_sub_per




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
!1) SUBROUTINE ext_field_shell: calcolo il campo esterno alla shell in coordinate
! sferiche locali: contributo outgoing, calcolato su argomenti complessi
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE ext_field_shell(v_amnbmn,v_emn,nstop,rho,theta,phi,Er,Etheta,Ephi,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn,v_emn	!Vettore campo interno e normalizzazione
INTEGER(lo), INTENT(IN) :: nstop					!Numero di espansioni multipolari
COMPLEX(dbl), INTENT(IN) :: rho					!k*r con r coordinata sferica
REAL(dbl), INTENT(IN) :: theta,phi 					!Coordinate sferiche locali
COMPLEX(dbl), INTENT(OUT) :: Er,Etheta,Ephi			!Campo calcolato in coordinate sferiche
INTEGER(lo) , INTENT(OUT) :: error				!Flag di errore

!Dichiarazione variabili interne
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_legmn	!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_h,v_dh			!Vettori funzioni radiali
INTEGER(lo) :: n,m,i,dimens							!Indici e dimensione
REAL(dbl) :: nr,mr								!Indici reali
REAL(dbl) :: nc,mc								!Indici complessi
!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)

ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_legmn(1:dimens),v_h(0:nstop),v_dh(1:nstop))

!Funzione angolare Pi_mn
CALL pi_mn(nstop,theta,v_pimn,error)

pimn_if: IF (error/=0) THEN																	
    			WRITE(*,10) 
    			10 FORMAT ("Si e' verificato un errore in pi_mn chiamata da ext_field: il programma termina ora...") 
       			STOP
END IF pimn_if

!Funzione angolare Tau_mn
CALL tau_mn(nstop,theta,v_pimn,v_taumn,error)

taumn_if: IF (error/=0) THEN																	
    			WRITE(*,20) 
    			20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da ext_field: il programma termina ora...") 
       			STOP
END IF taumn_if

!Funzione angolare legendre
CALL leg_mn(nstop,theta,v_legmn,error)

legmn_if: IF (error/=0) THEN																	
    			WRITE(*,40) 
    			40 FORMAT ("Si e' verificato un errore in leg_mn chiamata da int_field: il programma termina ora...") 
       			STOP
END IF legmn_if

!Funzione radiale hankel1_n
CALL hankel1_z_sub(nstop,rho,v_h,error)

hankel_if: IF (error/=0) THEN																	
    			WRITE(*,30) 
    			30 FORMAT ("Si e' verificato un errore in hankel1_z_sub chiamata da ext_field: il programma termina ora...") 
       			STOP
END IF hankel_if

!Funzione radiale (1/rho)*D[rho*h_n[rho],rho]
dh_do: DO n=1,nstop
	nc=CMPLX(n,KIND=dbl)
	v_dh(n)=v_h(n-1)-nc*v_h(n)/rho
END DO dh_do




!Calcolo del campo radiale
i=0
Er=(0.0D0,0.0D0)
Ern_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Erm_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Er=Er+v_emn(i)*nr*(nr+1.0D0)*v_amnbmn(2*i-1)*v_legmn(i)*(v_h(n)/rho)*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Erm_do
END DO Ern_do

Er=(0.0D0,1.0D0)*Er

!Calcolo del campo theta
i=0
Etheta=(0.0D0,0.0D0)
Ethetan_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ethetam_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Etheta=Etheta+v_emn(i)*( (0.0D0,1.0D0)*v_amnbmn(2*i-1)*v_taumn(i)*v_dh(n)  -&
			 & v_amnbmn(2*i)*v_pimn(i)*v_h(n))*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ethetam_do
END DO Ethetan_do

!Calcolo del campo phi
i=0
Ephi=(0.0D0,0.0D0)
Ephin_do: DO n=1,nstop

	nr=REAL(n,dbl)

	Ephim_do: DO m=-n,n
		
		i=i+1
		mr=REAL(m,dbl)
		
		Ephi=Ephi-v_emn(i)*( v_amnbmn(2*i-1)*v_pimn(i)*v_dh(n) + &
			 & v_amnbmn(2*i)*v_taumn(i)*v_h(n))*EXP((0.0D0,1.0D0)*mr*phi)

	END DO Ephim_do
END DO Ephin_do


DEALLOCATE(v_pimn,v_taumn,v_legmn,v_h,v_dh)

END SUBROUTINE ext_field_shell



!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE Exyz_ss_sub: dato un punto (x,y,z) nello spazio, calcolo il campo
! elettrico totale in coordinate cartesiane,sia che sia dentro sia che sia fuori
! una sfera, nel caso di una shell isolata di qualunque simmetria
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Exyz_ss_sub(v_kinc,v_Einc,flaginc,nstop,ratio,lambda,x,y,z,v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,&
			 &v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,Ex,Ey,Ez,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc				!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc							!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn					!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn_shell,v_dmncmn_shell		!Vettori campo shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn					!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn						!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq					!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req						!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index								!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt					!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio							!Numero di espansioni multipolari e raggi non calcolo
REAL(dbl), INTENT(IN) :: x,y,z								!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda								!Lunghezza d'onda incidente
COMPLEX(dbl), INTENT(OUT) :: Ex,Ey,Ez							!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Exinc,Eyinc,Ezinc,pshift							!Onda incidente cartesiana
COMPLEX(dbl) :: Er,Etheta,Ephi								!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: Er_in,Etheta_in,Ephi_in,Er_out,Etheta_out,Ephi_out		!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: rhoc,mc										!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase										!Rho reale
REAL(dbl) :: r,theta,phi 									!Coordinate sferiche locali
REAL(dbl) :: xl,yl,zl	 									!Coordinate cartesiane locali
REAL(dbl) :: dcore,dshell,rcore,rshell							!Distanza punto-centrocore,punto-centroshell,r core e shell
INTEGER(lo) :: n,m,i,dimens,ns								!Indici e dimensione
REAL(dbl) :: nr,mr										!Indici reali

!Inizio subroutine vera e propria

!Calcolo il numero di sfere e dimens
dimens=nstop*(nstop+2)
ns=SIZE(m_xyz(:,1))

WRITE(*,*) "kinc",v_kinc
WRITE(*,*) "einc",v_Einc
WRITE(*,*) "flaginc",flaginc
WRITE(*,*) "nstop",nstop
WRITE(*,*) "ratio",ratio
WRITE(*,*) "lambda",lambda
WRITE(*,*) "x,y,z",x,y,z
WRITE(*,*)
WRITE(*,*) "m_xyz",m_xyz
WRITE(*,*)
WRITE(*,*) "m_epseq",m_epseq
WRITE(*,*)
WRITE(*,*) "v_req",v_req
WRITE(*,*)
WRITE(*,*) "ref_index",ref_index
WRITE(*,*)
WRITE(*,*) "ns",ns
WRITE(*,*)

!Calcolo raggio e distanza dal centro del core e dal centro della shell
rcore=v_req(1)
rshell=v_req(2)
dcore=SQRT((x-m_xyz(1,1))**2 + (y-m_xyz(1,2))**2 + (z-m_xyz(1,3))**2)
dshell=SQRT((x-m_xyz(2,1))**2 + (y-m_xyz(2,2))**2 + (z-m_xyz(2,3))**2)

WRITE(*,*) "rcore,rshell",rcore,rshell
WRITE(*,*) "dcore,dshell",dcore,dshell
WRITE(*,*)

!Inizializzo i campi
Ex=(0.0D0,0.0D0)
Ey=(0.0D0,0.0D0)
Ez=(0.0D0,0.0D0)

!Qui faccio un'if a seconda di dove mi trovo, dentro al core, dentro alla shell, o fuori.
position_if: IF (dcore<rcore) THEN

	!-------------------------------------------
	!Se sono qui sono dentro al core
	!-------------------------------------------
				
	!Calcolo le coordinate sferiche locali,nel sistema di riferimento del core
	xl=x-m_xyz(1,1)
	yl=y-m_xyz(1,2)
	zl=z-m_xyz(1,3)
	
	WRITE(*,*) "xl,yl,zl",xl,yl,zl
	WRITE(*,*)
	
	CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

	WRITE(*,*) "r,theta,phi",r,theta,phi
	WRITE(*,*)

	c_s_if: IF (error/=0) THEN
		error=1
		WRITE(*,100) 
		100 FORMAT ("Errore in cart_spher chiamata da Exyz_ss_sub: il programma termina ora...")
			RETURN
	END IF c_s_if

	!Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
	mc=SQRT(CMPLX(m_epseq(1,1),m_epseq(1,2),KIND=dbl))
	rhoc=2*Pi_D*mc*r/lambda
	
	WRITE(*,*) "mc,rhoc",mc,rhoc
	WRITE(*,*)

	!Calcolo il campo in coordinate sferiche locali
	CALL int_field(v_dmncmn,v_emn,nstop,rhoc,theta,phi,Er,Etheta,Ephi,error)

	WRITE(*,*) "Er,Etheta,Ephi",Er,Etheta,Ephi
	WRITE(*,*)

	intfield_if: IF (error/=0) THEN
		error=1
		WRITE(*,20) 
		20 FORMAT ("Errore in int_field chiamata da Exyz_ss_sub: il programma termina ora...")
			RETURN
	END IF intfield_if

	!Conversione in coordinate cartesiane e uscita dalla subroutine
	Ex=Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
	Ey=Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
	Ez=Er*COS(theta) - Etheta*SIN(theta)

	WRITE(*,*) "Ez,Ex,Ey",Ez,Ex,Ey
	WRITE(*,*)

	RETURN

ELSE IF (dshell<rshell) THEN

	!-------------------------------------------
	!se sono qui sono dentro alla shell
	!-------------------------------------------

	!Calcolo le coordinate sferiche locali, sempre nel sistema di riferimento del core
	!quindi uso i coefficienti di espansione del core, quelli con l'indice 2.
	xl=x-m_xyz(1,1)
	yl=y-m_xyz(1,2)
	zl=z-m_xyz(1,3)
	CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

	c_s_if0: IF (error/=0) THEN
		error=1
		WRITE(*,101) 
		101 FORMAT ("Errore in cart_spher chiamata da Exyz_ss_sub: il programma termina ora...")
			RETURN
	END IF c_s_if0

	!Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
	mc=SQRT(CMPLX(m_epseq(2,1),m_epseq(2,2),KIND=dbl))
	rhoc=2*Pi_D*mc*r/lambda


	!----------------------------------------------------------------------------------------------
	!Qui calcolo il campo della shell, quindi sommo campo interno e campo esterno
	!----------------------------------------------------------------------------------------------

	!Calcolo il campo nella shell in coordinate sferiche locali, il contributo ingoing
	CALL int_field(v_dmncmn_shell,v_emn,nstop,rhoc,theta,phi,Er_in,Etheta_in,Ephi_in,error)

	intfield_shell_if: IF (error/=0) THEN
		error=1
		WRITE(*,11) 
		11 FORMAT ("Errore in int_field chiamata da Exyz_ss_sub: il programma termina ora...")
			RETURN
	END IF intfield_shell_if

	!Calcolo il campo nella shell in coordinate sferiche locali, il contributo outgoing
	CALL ext_field_shell(v_amnbmn_shell,v_emn,nstop,rhoc,theta,phi,Er_out,Etheta_out,Ephi_out,error)

	extfield_shell_if: IF (error/=0) THEN
		WRITE(*,21) 
		21 FORMAT ("Errore in ext_field chiamata da Exyz_ss_sub: il programma termina ora...")
		RETURN
	END IF extfield_shell_if

!	Er_out=(zerodbl,zerodbl)
!	Etheta_out=(zerodbl,zerodbl)
!	Ephi_out=(zerodbl,zerodbl)

	!Conversione in coordinate cartesiane e uscita dalla subroutine
	Ex=(Er_in+Er_out)*SIN(theta)*COS(phi) + (Etheta_in+Etheta_out)*COS(theta)*COS(phi) - (Ephi_in+Ephi_out)*SIN(phi)
	Ey=(Er_in+Er_out)*SIN(theta)*SIN(phi) + (Etheta_in+Etheta_out)*COS(theta)*SIN(phi) + (Ephi_in+Ephi_out)*COS(phi)
	Ez=(Er_in+Er_out)*COS(theta) - (Etheta_in+Etheta_out)*SIN(theta)
!	Ex=(zerodbl,zerodbl)
!	Ey=(zerodbl,zerodbl)
!	Ez=(zerodbl,zerodbl)
	RETURN

ELSE

	!--------------------------------------------
	!Se sono qui allora sono fuori dalla sfera
	!--------------------------------------------

	!Calcolo le coordinate sferiche locali
	xl=x-m_xyz(2,1)
	yl=y-m_xyz(2,2)
	zl=z-m_xyz(2,3)
	CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

	c_s_if1: IF (error/=0) THEN
		WRITE(*,102) 
		102 FORMAT ("Errore in cart_spher chiamata da Exyz_ss_sub: il programma termina ora...")
		RETURN
	END IF c_s_if1

	!Calcolo  rho per chiamare la subroutine per il calcolo del campo
	rhor=2*Pi_D*ref_index*r/lambda

	!Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
	CALL ext_field(v_amnbmn,v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)

	extfield_if: IF (error/=0) THEN
		WRITE(*,22) 
		22 FORMAT ("Errore in ext_field_shell chiamata da Exyz_ss_sub: il programma termina ora...")
		RETURN
	END IF extfield_if

	!Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
	Ex=Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
	Ey=Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
	Ez=Er*COS(theta) - Etheta*SIN(theta)

END IF position_if 

!Adesso, se la flag e' quella attesa, aggiungo anche il campo incidente, tanto il calcolo e' banale
inc_if: IF (flaginc=="yes") THEN

!    WRITE(*,*) 'v_kinc',v_kinc

	phase=x*v_kinc(1)+y*v_kinc(2)+z*v_kinc(3)
	pshift=EXP((0.0D0,1.0D0)*phase)
	Exinc=CMPLX(v_Einc(1),0.0D0)*pshift
	Eyinc=CMPLX(v_Einc(2),0.0D0)*pshift
	Ezinc=CMPLX(v_Einc(3),0.0D0)*pshift
	
	WRITE(*,*) 'v_Exinc',Exinc
	WRITE(*,*) 'v_Eyinc',Eyinc
	WRITE(*,*) 'v_Eyinc',Eyinc
	
	Ex=Ex+Exinc
	Ey=Ey+Eyinc
	Ez=Ez+Ezinc

END IF inc_if

END SUBROUTINE Exyz_ss_sub



!******************************************************************************
!******************************************************************************
!******************************************************************************
!3) SUBROUTINE Efield_ss sub: calcolo il campo su tutto un piano in coordinate  
! cartesiane, per il caso della shell
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efield_ss_sub(v_kinc,v_Einc,flaginc,m_euler,xmin,xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep,&
     & nstop,ratio,lambda,v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,m_E,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc				!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc							!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_euler								!Matrice di eulero
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn					!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn					!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn_shell,v_dmncmn_shell		!Vettori campo shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn						!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq					!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req						!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index								!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt					!Vettore pattern sfere identiche
INTEGER(lo), INTENT(IN) :: nstop,ratio    						!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: xmin,xmax,ymin,ymax,zmin,zmax				!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda								!Lunghezza d'onda incidente
INTEGER(lo), INTENT(IN) :: xstep,ystep,zstep						!Numero di step
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_E					!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,pshift					!Campo calcolato in coordinate sferiche e phase shift c inc
COMPLEX(dbl) :: rhoc,mc							!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase							!Rho reale e fase campo inc
REAL(dbl) :: r,theta,phi 						!Coordinate sferiche locali
REAL(dbl) :: x,y,z,dx,dy,dz,dtheta,dphi				!Coordinate cartesiane locali,delta x,y,z
REAL(dbl) :: d,radius							!Distanza punto-sfere,raggio sfera
REAL(dbl), DIMENSION(1:3) :: v_xyz					!Vettore coordinate da ruotare
INTEGER(lo) :: i,j							!Indici e dimensione
REAL(dbl) :: nr,mr							!Indici reali

!Inizio subroutine vera e propria

!Controllo di avere un piano ben definito
plane_if: IF ((xmin/=xmax) .AND. (ymin/=ymax) .AND. (zmin/=zmax)) THEN
     error=1
     WRITE(*,10) 
10   FORMAT ("Errore in Efield_sub: nessun piano definito, il programma termina ora...")                           
     RETURN
END IF plane_if

!Controllo di avere i bounds ben definiti
bounds_if: IF ((xmin>xmax) .OR. (ymin>ymax) .OR. (zmin>zmax)) THEN
     error=1
     WRITE(*,20) 
20   FORMAT ("Errore in Efield_sub: bounds mal definiti, il programma termina ora...")                           
     RETURN
END IF bounds_if

!Controllo di avere gli steps ben definiti
steps_if: IF ((xstep<2) .OR. (ystep<2) .OR. (zstep<2)) THEN
     error=1
     WRITE(*,30) 
30   FORMAT ("Errore in Efield_sub: steps mal definiti, il programma termina ora...")                           
     RETURN
END IF steps_if

!Calcolo i delta a seconda degli step
dx=(xmax-xmin)/REAL(xstep-1,dbl)
dy=(ymax-ymin)/REAL(ystep-1,dbl)
dz=(zmax-zmin)/REAL(zstep-1,dbl)
dtheta=Pi_D/REAL(180-1,dbl)
dphi=2.0D0*Pi_D/REAL(180-1,dbl)

!Comincio i conti del campo e un grande if 
x=xmin
y=ymin
z=zmin
theta=0.0D0
phi=0.0D0

xyz_if: IF ((xmin==xmax).AND.(ymin==ymax).AND.(zmin==zmax)) THEN !Se e' tutto uguale,e' una flag per calcolare sulla sfera



     theta_do: DO i=1,180

          phi=0.0D0

          phi_do: DO j=1,360

!Inizializzo i vettori da ruotare
               v_xyz(1)=(v_req(1)*1.05)*SIN(theta)*COS(phi)
               v_xyz(2)=(v_req(1)*1.05)*SIN(theta)*SIN(phi)
               v_xyz(3)=(v_req(1)*1.05)*COS(theta)

               v_xyz=MATMUL(m_euler,v_xyz)

               m_E(i,j,1)=theta
               m_E(i,j,2)=phi
               CALL Exyz_ss_sub(v_kinc,v_Einc,flaginc,nstop,&
                    &ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,&
                    &v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

               Exyz0_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,35) 
35                  FORMAT ("Errore in Exyz_ss_sub chiamata da Efield_sub: il programma termina ora...")                           
                    RETURN
               END IF Exyz0_if


               phi=phi+dphi

          END DO phi_do

          theta=theta+dtheta

     END DO theta_do



ELSE IF (xmin==xmax) THEN


     y_yz_do: DO i=1,ystep

          z=zmin

          z_yz_do: DO j=1,zstep

!Inizializzo i vettori da ruotare
               v_xyz(1)=x
               v_xyz(2)=y
               v_xyz(3)=z

               v_xyz=MATMUL(m_euler,v_xyz)

               m_E(i,j,1)=y
               m_E(i,j,2)=z
               CALL Exyz_ss_sub(v_kinc,v_Einc,flaginc,nstop,&
                    &ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,&
                    &v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

               Exyz1_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,40) 
40                  FORMAT ("Errore in Exyz_ss_sub chiamata da Efield_sub: il programma termina ora...")                           
                    RETURN
               END IF Exyz1_if


!			!Calcolo raggio e distanza dal centro del core e dal centro della shell
!			rcore=v_req(1)
!			rshell=v_req(2)
!			dcore=SQRT((v_xyz(1)-m_xyz(1,1))**2 + (v_xyz(2)-m_xyz(1,2))**2 + (v_xyz(3)-m_xyz(1,3))**2)
!			dshell=SQRT((v_xyz(1)-m_xyz(2,1))**2 + (v_xyz(2)-m_xyz(2,2))**2 + (v_xyz(3)-m_xyz(2,3))**2)

!			!Adesso, se la flag e' quella attesa, aggiungo anche il campo incidente, tanto il calcolo e' banale
!			inc_if: IF ( (flaginc=="yes").AND.(dshell>rshell) ) THEN

!				phase=x*v_kinc(1)+y*v_kinc(2)+z*v_kinc(3)
!				pshift=EXP((0.0D0,1.0D0)*phase)
!				Exinc=CMPLX(v_Einc(1),0.0D0)*pshift
!				Eyinc=CMPLX(v_Einc(2),0.0D0)*pshift
!				Ezinc=CMPLX(v_Einc(3),0.0D0)*pshift
!				Ex=Ex+Exinc
!				Ey=Ex+Eyinc
!				Ez=Ex+Ezinc

!			END IF inc_if





               z=z+dz

          END DO z_yz_do

          y=y+dy

     END DO y_yz_do

ELSE IF (ymin==ymax) THEN

     x_xz_do: DO i=1,xstep

          z=zmin

          z_xz_do: DO j=1,zstep

!Inizializzo i vettori da ruotare
               v_xyz(1)=x
               v_xyz(2)=y
               v_xyz(3)=z

               v_xyz=MATMUL(m_euler,v_xyz)

               m_E(i,j,1)=x
               m_E(i,j,2)=z
               CALL Exyz_ss_sub(v_kinc,v_Einc,flaginc,nstop,&
                    & ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,&
                    & v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

               Exyz2_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,50) 
50                  FORMAT ("Errore in Exyz_ss_sub chiamata da Efield_sub: il programma termina ora...")                           
                    RETURN
               END IF Exyz2_if

               z=z+dz

          END DO z_xz_do

          x=x+dx

     END DO x_xz_do

ELSE 

     x_xy_do: DO i=1,xstep

          y=ymin

          y_xy_do: DO j=1,ystep

!Inizializzo i vettori da ruotare
               v_xyz(1)=x
               v_xyz(2)=y
               v_xyz(3)=z

               v_xyz=MATMUL(m_euler,v_xyz)

               m_E(i,j,1)=x
               m_E(i,j,2)=y
               CALL Exyz_ss_sub(v_kinc,v_Einc,flaginc,nstop,&
                    & ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,&
                    & v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
                    & m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)


               Exyz3_if: IF (error/=0) THEN
                    error=1
                    WRITE(*,60) 
60                  FORMAT ("Errore in Exyz_ss_sub chiamata da Efield_sub: il programma termina ora...")                           
                    RETURN
               END IF Exyz3_if

               y=y+dy

          END DO y_xy_do

          x=x+dx

     END DO x_xy_do

END IF xyz_if

END SUBROUTINE Efield_ss_sub


!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE Efar_ss_sub: calcolo il campo su tutto un piano in coordinate  
! cartesiane
!
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_ss_sub(k,phimin,phimax,phistep,thetastep,betap,v_amnbmn,v_xyz,m_euler,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k							!Wavevector
REAL(dbl), INTENT(IN) :: phimin,phimax					!Angoli per il far field
REAL(dbl), INTENT(INOUT) :: betap						!Polarizzazione
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
REAL(dbl), DIMENSION(:,:) ::m_euler						!Matrice di eulero
INTEGER(lo), INTENT(IN) :: phistep,thetastep				!Intervalli di calcolo
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC			!Quantita' scattering
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: scatot						!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: s1,s2,s3,s4,sum1,sum2,sum3,sum4							!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,theta1,r1,phi1,phase,dtheta,dphi					!Piani e angoli di scattering
REAL(dbl) :: nr,mr											!Indici reali
REAL(dbl) :: logw,fr											!logaritmo per cmn e frazione
REAL(dbl) :: s11,s12,dcpar,dcnorm,dc,cs,ss,dc1							!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta					!Indici e dimensione
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_psimn, v_phimn, v_thetamn,	v_ximn	!Vettori funzioni angolari
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta								!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_cmn				!Vettori funzioni angolari
REAL(dbl), DIMENSION(1:3) :: v_xyz_rot									!Vettori funzioni angolari


!Inizio subroutine vera e propria
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_cmn(1:dimens))
ALLOCATE(v_psimn(1:side),v_phimn(1:side),v_thetamn(1:side),v_ximn(1:side))

!Calcolo i delta a seconda degli step
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=180.0D0/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0

!coefficienti di normalizzazione per le funzioni angolari
error=0
i=0
n_do_cmn: DO n=1,nstop

	nr=REAL(n,dbl)

	m_do_cmn: DO m=-n,n
	
		i=i+1
		mr=REAL(m,dbl)
		
		logw=lnf(nr-mr,error)-lnf(nr+mr,error)
		v_cmn(i)=SQRT(EXP(logw)*(2*nr+1.0D0)/(nr*(nr+1.0D0)))
	
	END DO m_do_cmn
END DO n_do_cmn


!Controllo di avere un piano ben definito
phi_if: IF (phimax<phimin) THEN
	error=1
	WRITE(*,10) 
	10 FORMAT ("Errore in Efar_ss_sub: phi bounds non corretti, il programma termina ora...")
	RETURN
END IF phi_if


!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if


!Inizializzo gli angoli
scatot=0.0D0

!Inizializzo gli angoli
theta=0.0D0
phi=phimin

theta_do: DO itheta=1,thetastep+1

	phi=phimin

	phi_do: DO iphi=1,phistep+1

		v_xyz_rot(1)=SIN(theta)*COS(phi)
		v_xyz_rot(2)=SIN(theta)*SIN(phi)
		v_xyz_rot(3)=COS(theta)

		v_xyz_rot=MATMUL(m_euler,v_xyz_rot)

		CALL cart_spher_r_dbl1 (v_xyz_rot(1),v_xyz_rot(2),v_xyz_rot(3),r1,theta1,phi1,error)

		!Funzione angolare Pi_mn
		CALL pi_mnsca(nstop,theta1,v_pimn,error)
		
		pimn_if: IF (error/=0) THEN
					WRITE(*,11) 
					11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
					STOP
		END IF pimn_if
		
		!Funzione angolare Tau_mn
		CALL tau_mnsca(nstop,theta1,v_pimn,v_taumn,error)
		
		taumn_if: IF (error/=0) THEN
					WRITE(*,20) 
					20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
					STOP
		END IF taumn_if

		!Normalizzo le funzioni angolari
		v_taumn=v_cmn*v_taumn
		v_pimn=v_cmn*v_pimn
			
		!Inizializzo il ciclo per i vettori di funzioni che mi servono nel calcolo di S
		j=0
		jm=0
		l_do: DO l=1,ns
		
			j0=0
		
			n_do: DO n=1,nstop
				m_do: DO m=-n,n 

					j0=n*(n+1)+m
					j=(l-1)*dimens+n*(n+1)+m	
					jm=(l-1)*dimens+n*(n+1)-m
					
					v_psimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) + &
					& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))
					
					v_phimn(j)=v_amnbmn(2*j-1)*v_taumn(j0) + v_amnbmn(2*j)*v_pimn(j0) - &
					& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_taumn(j0) - v_amnbmn(2*jm)*v_pimn(j0))
				
					v_thetamn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) - &
					& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))
				
					v_ximn(j)=v_amnbmn(2*j-1)*v_pimn(j0) + v_amnbmn(2*j)*v_taumn(j0) + &
					& ((-1.0D0)**m) * (v_amnbmn(2*jm-1)*v_pimn(j0) - v_amnbmn(2*jm)*v_taumn(j0))
				
				END DO m_do
			END DO n_do
		END DO l_do

		!Riempio il vettore degli esponenziali
		delta_do: DO l=1,ns
		
			phase=v_xyz(1)*SIN(theta1)*COS(phi1) + v_xyz(2)*SIN(theta1)*SIN(phi1) + v_xyz(3)*COS(theta1)
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)
			
			!If ((itheta==1) .AND. (iphi)==1) THEN
			!	WRITE(*,*) v_delta(l)
			!END IF
		
		END DO delta_do

		!Coordinate angolari
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi
		
		!Inizializzo il ciclo per le somme di S
		j=0
		jm=0
		s1=(0.0D0,0.0D0)
		s2=(0.0D0,0.0D0)
		s3=(0.0D0,0.0D0)
		s4=(0.0D0,0.0D0)
		l_do_phi: DO l=1,ns
		
			j0=0
			sum1=(0.0D0,0.0D0)
			sum2=(0.0D0,0.0D0)
			sum3=(0.0D0,0.0D0)
			sum4=(0.0D0,0.0D0)
			n_do_phi: DO n=1,nstop
				m_do_phi: DO m=0,n 
	
					j0=n*(n+1)+m
					j=(l-1)*dimens+n*(n+1)+m	
					jm=(l-1)*dimens+n*(n+1)-m
					mr=REAL(m,dbl)
					
					fr_if: IF (m==0) THEN 
						fr=0.5D0
					ELSE
						fr=1.0D0
					END IF fr_if
					
					cs=COS((mr-1.0D0)*phi1+betap)
					ss=SIN((mr-1.0D0)*phi1+betap)
					
					sum1=sum1+fr*(v_ximn(j)*cs + (0.0D0,1.0D0)*v_thetamn(j)*ss)
					sum2=sum2+fr*(v_psimn(j)*cs + (0.0D0,1.0D0)*v_phimn(j)*ss)
					sum3=sum3+fr*((0.0D0,1.0D0)*v_phimn(j)*cs - v_psimn(j)*ss)
					sum4=sum4+fr*(-(0.0D0,1.0D0)*v_thetamn(j)*cs + v_ximn(j)*ss)
				
				END DO m_do_phi
			END DO n_do_phi
			
			!Elementi dell'amplitude scattering matrix
			s1=s1+v_delta(l)*sum1
			s2=s2+v_delta(l)*sum2
			s3=s3+v_delta(l)*sum3
			s4=s4+v_delta(l)*sum4
			
		END DO l_do_phi
		
		!Elementi della matrice di muller
		s11=0.5D0*((ABS(s1)**2)+(ABS(s2)**2)+(ABS(s3)**2)+(ABS(s4)**2))
		s12=0.5D0*((ABS(s2)**2)-(ABS(s1)**2)+(ABS(s4)**2)-(ABS(s3)**2))
		
		!Sezione angolare differenziale parallela,normale e complessiva
		dcpar= (s11+s12)/(k**2)
		dcnorm=(s11-s12)/(k**2)
		dc=((COS(phi1-betap))**2)*dcpar + ((SIN(phi1-betap))**2)*dcnorm
		scatot=scatot+dc*SIN(theta1)*dtheta*dphi
				
		!Coordinate angolari
		m_SC(itheta,iphi,3)=s1
		m_SC(itheta,iphi,4)=s2
		m_SC(itheta,iphi,5)=s3
		m_SC(itheta,iphi,6)=s4
		m_SC(itheta,iphi,7)=dc
		
		phi=phi+dphi

	END DO phi_do
	
	theta=theta+dtheta
	
END DO theta_do

END SUBROUTINE Efar_ss_sub







!******************************************************************************
!******************************************************************************
!******************************************************************************
!4) SUBROUTINE Efar_ss_sub_poynting: calcolo lo scattering per la shell,pure
!integro, ma faccio con il vettore di poynting cosi' non ho scazzi coi
!sistemi di riferimento e le polarizzazioni
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efar_ss_sub_poynting(k,phimin,phimax,phistep,thetastep,v_amnbmn,v_xyz,m_euler,m_SC,scatot,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k							!Wavevector
REAL(dbl), INTENT(IN) :: phimin,phimax					!Angoli per il far field
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_xyz				!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn			!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				!Intervalli di calcolo
REAL(dbl), DIMENSION(:,:) ::m_euler						!Matrice di eulero
COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_SC			!Quantita' scattering
INTEGER(lo) , INTENT(OUT) :: error					!Flag di errore
REAL(dbl), INTENT(OUT) :: scatot						!Sezione di scattering totale

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,Htheta,Hphi,sum1,sum2,sum3,sum4,dccomplex		!Somme parziali e finali per gli elementi S
REAL(dbl) :: theta,phi,phase,r1,theta1,phi1,dtheta,dphi				!Piani e angoli di scattering
REAL(dbl) :: dc											!sezione diff
INTEGER(lo) :: i,l,n,m,j0,j,jm,blockside,side,iphi,itheta				!Indici e dimensione
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta							!Vettori funzioni angolari
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn				!Vettori funzioni angolari
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_amn_taumn,v_amn_pimn,v_bmn_taumn,v_bmn_pimn,v_emn,v_ab_tot	!Vettori conti campo
REAL(dbl), DIMENSION(1:3) :: v_xyz_rot							!Vettori funzioni angolari

!Inizio subroutine vera e propria
blockside=nstop*(nstop+2)
side=ns*blockside
ALLOCATE(v_pimn(1:blockside),v_taumn(1:blockside),v_emn(1:blockside),v_ab_tot(1:2*blockside))
ALLOCATE(v_amn_pimn(1:blockside),v_amn_taumn(1:blockside),v_bmn_pimn(1:blockside),v_bmn_taumn(1:blockside))

!Calcolo i delta a seconda degli step
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=180.0D0/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0

!coefficienti di normalizzazione per le funzioni angolari
error=0
CALL emn_sub(nstop,v_emn,error)

!Controllo di avere gli steps ben definiti
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efar_sub: theta steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if

scatot=0.0D0

!Comincio i conti del campo e un grande if 
theta=0.0D0
phi=phimin

theta_do: DO itheta=1,thetastep+1

	phi=phimin

	phi_do: DO iphi=1,phistep+1


		v_xyz_rot(1)=SIN(theta)*COS(phi)
		v_xyz_rot(2)=SIN(theta)*SIN(phi)
		v_xyz_rot(3)=COS(theta)

		v_xyz_rot=MATMUL(m_euler,v_xyz_rot)

		CALL cart_spher_r_dbl1 (v_xyz_rot(1),v_xyz_rot(2),v_xyz_rot(3),r1,theta1,phi1,error)

		!Funzione angolare Pi_mn
		CALL pi_mnsca(nstop,theta1,v_pimn,error)

		pimn_if: IF (error/=0) THEN
					WRITE(*,11) 
					11 FORMAT ("Si e' verificato un errore in pi_mn chiamata da Efar_sub: il programma termina ora...") 
					STOP
		END IF pimn_if

		!Funzione angolare Tau_mn
		CALL tau_mnsca(nstop,theta1,v_pimn,v_taumn,error)

		taumn_if: IF (error/=0) THEN	
					WRITE(*,20) 
					20 FORMAT ("Si e' verificato un errore in tau_mn chiamata da Efar_sub: il programma termina ora...") 
					STOP
		END IF taumn_if

		!Riempio il vettore degli esponenziali
		delta_do: DO l=1,ns

			phase=v_xyz(1)*SIN(theta1)*COS(phi1) + v_xyz(2)*SIN(theta1)*SIN(phi1) + &
			& v_xyz(3)*COS(theta1)
			
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)

		END DO delta_do
	
		!Coordinate angolari
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi
		
		!Qui calcolo i coefficienti totali del campo,questo approccio mi dovrebbe dare il campo totale nel far field
		!e quindi dovrei calcolare tutto unitariamente e senza problemi!
		
		!inizializzo i coefficienti totali
		v_ab_tot=0.0D0
		
		tot_sca_ndo: DO n=1,nstop
			tot_sca_mdo: DO m=-n,n
				tot_sca_ldo: DO l=1,ns

					j0=n*(n+1)+m
					j=(l-1)*blockside+n*(n+1)+m
					v_ab_tot(2*j0-1)=v_ab_tot(2*j0-1)+v_delta(l)*v_amnbmn(2*j-1)
					v_ab_tot(2*j0)=v_ab_tot(2*j0)+v_delta(l)*v_amnbmn(2*j)

				END DO tot_sca_ldo
			END DO tot_sca_mdo
		END DO tot_sca_ndo
		
		
		!Qui calcolo i vettori misti taumn pimn amn bmn per il calcolo dei campi
		v_amn_taumn=v_ab_tot(1:2*blockside:2)*v_taumn
		v_bmn_taumn=v_ab_tot(2:2*blockside:2)*v_taumn
		v_amn_pimn=v_ab_tot(1:2*blockside:2)*v_pimn
		v_bmn_pimn=v_ab_tot(2:2*blockside:2)*v_pimn
		
		!Inizializzo il ciclo per le somme di S
		Etheta=(0.0D0,0.0D0)
		Ephi=(0.0D0,0.0D0)
		Htheta=(0.0D0,0.0D0)
		Hphi=(0.0D0,0.0D0)
			
		j0=0
		n_do_phi: DO n=1,nstop
		
			sum1=(0.0D0,0.0D0)
			sum2=(0.0D0,0.0D0)
			sum3=(0.0D0,0.0D0)
			sum4=(0.0D0,0.0D0)
			
			m_do_phi: DO m=-n,n 
				
				j0=n*(n+1)+m
				
				sum1=sum1+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum2=sum2+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum3=sum3+(v_emn(j0)*(v_bmn_taumn(j0)+v_amn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
				sum4=sum4+(v_emn(j0)*(v_amn_taumn(j0)+v_bmn_pimn(j0)))*EXP((0.0D0,1.0D0)*REAL(m,dbl)* phi)
			
			END DO m_do_phi
			
			Etheta=Etheta+((0.0D0,-1.0D0)**(n+1))*sum1
			Ephi=Ephi+((0.0D0,-1.0D0)**n)*sum2
			Htheta=Htheta+((0.0D0,-1.0D0)**n)*sum3
			Hphi=Hphi+((0.0D0,-1.0D0)**(n+1))*sum4
			
		END DO n_do_phi
		
		!Correggo dove occore per il segno
		Etheta=-Etheta
		Ephi=-Ephi
		Hphi=-Hphi
				
		!Sezione angolare differenziale complessa e reale
		dccomplex=(Etheta*CONJG(Hphi)-Ephi*CONJG(Htheta))/(k**2)
		dc=REAL(dccomplex,dbl)
		scatot=scatot+dc*SIN(theta)*dtheta*dphi
!		scatot=scatot+dc*v_wphi(iphi)*v_wtheta(itheta)
				
		!Coordinate angolari
		m_SC(itheta,iphi,3)=REAL(Etheta,dbl)
		m_SC(itheta,iphi,4)=REAL(Ephi,dbl)
		m_SC(itheta,iphi,5)=REAL(Htheta,dbl)
		m_SC(itheta,iphi,6)=REAL(Hphi,dbl)
		m_SC(itheta,iphi,7)=dc
		
		phi=phi+dphi

	END DO phi_do
	
	theta=theta+dtheta
	
END DO theta_do

END SUBROUTINE Efar_ss_sub_poynting









!******************************************************************************
!******************************************************************************
!******************************************************************************
!5) SUBROUTINE Exyz_ss_sub_dip: dato un punto (x,y,z) nello spazio, calcolo il campo
! elettrico totale in coordinate cartesiane,sia che sia dentro sia che sia fuori
! una sfera, nel caso di una shell isolata di qualunque simmetria, nel caso i
! dipoli siano dentro la sfera
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Exyz_ss_sub_dip(v_kinc,v_Einc,flaginc,dip_flag,nstop,ratio,lambda,x,y,z,v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,&
     &v_dmncmn,v_amnbmn_rhs,v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,Ex,Ey,Ez,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc				        !Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc							!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz						!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn					!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn_shell,v_dmncmn_shell		        !Vettori campo shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn,v_amnbmn_rhs		        	!Vettore campo interno e dipoli
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn						!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq					!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req						!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index							!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt					!Vettore pattern sfere identiche	
INTEGER(lo), INTENT(IN) :: nstop,ratio,dip_flag					!Numero di espansioni multipolari e raggi non calcolo

REAL(dbl), INTENT(IN) :: x,y,z								!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda								!Lunghezza d'onda incidente
COMPLEX(dbl), INTENT(OUT) :: Ex,Ey,Ez							!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error							!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Exinc,Eyinc,Ezinc,pshift							!Onda incidente cartesiana
COMPLEX(dbl) :: Er,Etheta,Ephi								!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: Er_in,Etheta_in,Ephi_in,Er_out,Etheta_out,Ephi_out		!Campo calcolato in coordinate sferiche
COMPLEX(dbl) :: rhoc,mc										!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase										!Rho reale
REAL(dbl) :: r,theta,phi 									!Coordinate sferiche locali
REAL(dbl) :: xl,yl,zl	 									!Coordinate cartesiane locali
REAL(dbl) :: dcore,dshell,rcore,rshell							!Distanza punto-centrocore,punto-centroshell,r core e shell
INTEGER(lo) :: n,m,i,dimens,ns,low_b,up_b						!Indici e dimensione
REAL(dbl) :: nr,mr										!Indici reali
REAL(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_xyz_int						!Vettore posizione senza shell

!Inizio subroutine vera e propria

!Calcolo il numero di sfere e dimens
dimens=nstop*(nstop+2)
ns=SIZE(m_xyz(:,1))

!Adesso ricolloco m_xyz saltando la shell
ALLOCATE(m_xyz_int(1:(ns-1),1:3))

m_xyz_int(1,:)=m_xyz(1,:)

int_do: DO i=3,(ns-1)+1
     m_xyz_int(i-1,:)=m_xyz(i,:)
END DO int_do

!Calcolo raggio e distanza dal centro del core e dal centro della shell
rcore=v_req(1)
rshell=v_req(2)
dcore=SQRT((x-m_xyz(1,1))**2 + (y-m_xyz(1,2))**2 + (z-m_xyz(1,3))**2)
dshell=SQRT((x-m_xyz(2,1))**2 + (y-m_xyz(2,2))**2 + (z-m_xyz(2,3))**2)

!Inizializzo i campi
Ex=(0.0D0,0.0D0)
Ey=(0.0D0,0.0D0)
Ez=(0.0D0,0.0D0)

!Qui faccio un'if a seconda di dove mi trovo, dentro al core, dentro alla shell, o fuori.
position_if: IF (dcore<rcore) THEN

     !-------------------------------------------
     !Se sono qui sono dentro al core
     !-------------------------------------------

     !Calcolo le coordinate sferiche locali,nel sistema di riferimento del core
     xl=x-m_xyz(1,1)
     yl=y-m_xyz(1,2)
     zl=z-m_xyz(1,3)
     CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

     c_s_if: IF (error/=0) THEN
          error=1
          WRITE(*,100) 
100       FORMAT ("Errore in cart_spher chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF c_s_if

     !Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
     mc=SQRT(CMPLX(m_epseq(1,1),m_epseq(1,2),KIND=dbl))
     rhoc=2*Pi_D*mc*r/lambda

     !Calcolo il campo in coordinate sferiche locali
     CALL int_field(v_dmncmn,v_emn,nstop,rhoc,theta,phi,Er,Etheta,Ephi,error)

     intfield_if: IF (error/=0) THEN
          error=1
          WRITE(*,20) 
20        FORMAT ("Errore in int_field chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF intfield_if

     !Conversione in coordinate cartesiane e uscita dalla subroutine
     Ex=Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
     Ey=Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
     Ez=Er*COS(theta) - Etheta*SIN(theta)

     !Voglio includere il campo dei dipoli?
     inc_if: IF ((flaginc=="yes") .AND. (dip_flag==0)) THEN

          !Qui faccio il ciclo sui dipoli invece, i margini del ciclo sono cosi perche' parto dalla seconda componente di amnbmn_rhs,
          !e quindi salto la shell, fino a ns-1, perche' il vettore e' lungo (ns-1) blocchi
          ext_do: DO i=2,ns-1

               !Calcolo le coordinate sferiche locali
               xl=x-m_xyz_int(i,1)
               yl=y-m_xyz_int(i,2)
               zl=z-m_xyz_int(i,3)
               CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

               c_s_if1: IF (error/=0) THEN
                    WRITE(*,30) 
30                  FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")
                    RETURN
               END IF c_s_if1

               !Calcolo  rho per chiamare la subroutine per il calcolo del campo, ovviamente l'indice di rifrazione e' quello interno
               !del core
               rhor=2*Pi_D*SQRT(m_epseq(1,1))*r/lambda
               !Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
               low_b=2*dimens*((i-1)-1)+1
               up_b=2*dimens*((i)-1)
               CALL ext_field(v_amnbmn_rhs(low_b:up_b),v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)

               extfield_if: IF (error/=0) THEN
                    WRITE(*,40) 
40                  FORMAT ("Errore in ext_field chiamata da Exyz_ss_sub_dip: il programma termina ora...")
                    RETURN
               END IF extfield_if

               !Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
               Ex=Ex+Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
               Ey=Ey+Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
               Ez=Ez+Er*COS(theta) - Etheta*SIN(theta)

          END DO ext_do

     END IF inc_if

     RETURN

ELSE IF (dshell<rshell) THEN

     !-------------------------------------------
     !se sono qui sono dentro alla shell
     !-------------------------------------------

     !Calcolo le coordinate sferiche locali, sempre nel sistema di riferimento del core
     !quindi uso i coefficienti di espansione del core, quelli con l'indice 2.
     xl=x-m_xyz(1,1)
     yl=y-m_xyz(1,2)
     zl=z-m_xyz(1,3)
     CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

     c_s_if0: IF (error/=0) THEN
          error=1
          WRITE(*,101) 
101       FORMAT ("Errore in cart_spher chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF c_s_if0

     !Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
     mc=SQRT(CMPLX(m_epseq(2,1),m_epseq(2,2),KIND=dbl))
     rhoc=2*Pi_D*mc*r/lambda


     !----------------------------------------------------------------------------------------------
     !Qui calcolo il campo della shell, quindi sommo campo interno e campo esterno
     !----------------------------------------------------------------------------------------------

     !Calcolo il campo nella shell in coordinate sferiche locali, il contributo ingoing
     CALL int_field(v_dmncmn_shell,v_emn,nstop,rhoc,theta,phi,Er_in,Etheta_in,Ephi_in,error)

     intfield_shell_if: IF (error/=0) THEN
          error=1
          WRITE(*,11) 
11        FORMAT ("Errore in int_field chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF intfield_shell_if

     !Calcolo il campo nella shell in coordinate sferiche locali, il contributo outgoing
     CALL ext_field_shell(v_amnbmn_shell,v_emn,nstop,rhoc,theta,phi,Er_out,Etheta_out,Ephi_out,error)

     extfield_shell_if: IF (error/=0) THEN
          WRITE(*,21) 
21        FORMAT ("Errore in ext_field chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF extfield_shell_if

     !	Er_out=(zerodbl,zerodbl)
     !	Etheta_out=(zerodbl,zerodbl)
     !	Ephi_out=(zerodbl,zerodbl)

     !Conversione in coordinate cartesiane e uscita dalla subroutine
     Ex=(Er_in+Er_out)*SIN(theta)*COS(phi) + (Etheta_in+Etheta_out)*COS(theta)*COS(phi) - (Ephi_in+Ephi_out)*SIN(phi)
     Ey=(Er_in+Er_out)*SIN(theta)*SIN(phi) + (Etheta_in+Etheta_out)*COS(theta)*SIN(phi) + (Ephi_in+Ephi_out)*COS(phi)
     Ez=(Er_in+Er_out)*COS(theta) - (Etheta_in+Etheta_out)*SIN(theta)
     !	Ex=(zerodbl,zerodbl)
     !	Ey=(zerodbl,zerodbl)
     !	Ez=(zerodbl,zerodbl)
     RETURN

ELSE

     !--------------------------------------------
     !Se sono qui allora sono fuori dalla sfera
     !--------------------------------------------

     !Calcolo le coordinate sferiche locali
     xl=x-m_xyz(2,1)
     yl=y-m_xyz(2,2)
     zl=z-m_xyz(2,3)
     CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

     c_s_if2: IF (error/=0) THEN
          WRITE(*,102) 
102       FORMAT ("Errore in cart_spher chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF c_s_if2

     !Calcolo  rho per chiamare la subroutine per il calcolo del campo
     rhor=2*Pi_D*ref_index*r/lambda

     !Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
     CALL ext_field(v_amnbmn,v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)

     extfield_if1: IF (error/=0) THEN
          WRITE(*,22) 
22        FORMAT ("Errore in ext_field_shell chiamata da Exyz_ss_sub: il programma termina ora...")
          RETURN
     END IF extfield_if1

     !Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
     Ex=Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
     Ey=Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
     Ez=Er*COS(theta) - Etheta*SIN(theta)

     !Voglio includere il campo dei dipoli?
     inc_out_if: IF ((flaginc=="yes") .AND. (dip_flag==1)) THEN

          !Qui faccio il ciclo sui dipoli invece, i margini del ciclo sono cosi perche' parto dalla seconda componente di amnbmn_rhs,
          !e quindi salto la shell, fino a ns-1, perche' il vettore e' lungo (ns-1) blocchi
          ext_out_do: DO i=2,ns-1

               !Calcolo le coordinate sferiche locali
               xl=x-m_xyz_int(i,1)
               yl=y-m_xyz_int(i,2)
               zl=z-m_xyz_int(i,3)
               CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

               c_s_out_if1: IF (error/=0) THEN
                    WRITE(*,35) 
35                  FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")
                    RETURN
               END IF c_s_out_if1

               !Calcolo  rho per chiamare la subroutine per il calcolo del campo, ovviamente l'indice di rifrazione e' quello interno
               !del core
               rhor=2*Pi_D*SQRT(m_epseq(1,1))*r/lambda
               !Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
               low_b=2*dimens*((i-1)-1)+1
               up_b=2*dimens*((i)-1)
               CALL ext_field(v_amnbmn_rhs(low_b:up_b),v_emn,nstop,rhor,theta,phi,Er,Etheta,Ephi,error)

               extfield_out_if: IF (error/=0) THEN
                    WRITE(*,45) 
45                  FORMAT ("Errore in ext_field chiamata da Exyz_ss_sub_dip: il programma termina ora...")
                    RETURN
               END IF extfield_out_if

               !Conversione in coordinate cartesiane e somma su i contributi di tutte le sfere
               Ex=Ex+Er*SIN(theta)*COS(phi) + Etheta*COS(theta)*COS(phi) - Ephi*SIN(phi)
               Ey=Ey+Er*SIN(theta)*SIN(phi) + Etheta*COS(theta)*SIN(phi) + Ephi*COS(phi)
               Ez=Ez+Er*COS(theta) - Etheta*SIN(theta)

          END DO ext_out_do

     END IF inc_out_if

END IF position_if

END SUBROUTINE Exyz_ss_sub_dip


!******************************************************************************
!******************************************************************************
!******************************************************************************
!6) SUBROUTINE Efield_ss sub_dip: calcolo il campo su tutto un piano in coordinate  
! cartesiane, per il caso della shell con dentro i dipoli
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE Efield_ss_sub_dip(v_kinc,v_Einc,flaginc,dip_flag,m_euler,xmin,xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep,&
& nstop,ratio,lambda,v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,v_amnbmn_rhs,&
& v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt,m_E,error)

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), DIMENSION(1:3), INTENT(IN) :: v_kinc,v_Einc				!Vettore campo e vettore d'onda
CHARACTER(len=3), INTENT(IN) :: flaginc						!Flag per aggiungere il campo incidente
REAL(dbl), DIMENSION(:,:) ::m_euler						!Matrice di eulero
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_xyz					!Posizione sfere
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn				!Vettore campo esterno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dmncmn,v_amnbmn_rhs			!Vettore campo interno
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_amnbmn_shell,v_dmncmn_shell		!Vettori campo shell
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_emn					!Vettore normalizzazione
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_epseq				!Funzioni dielettriche corrette
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req					!Raggi equivalenti
REAL(dbl), INTENT(IN) :: ref_index						!Indice di rifrazione matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_patt				!Vettore pattern sfere identiche
INTEGER(lo), INTENT(IN) :: nstop,ratio,dip_flag				!Numero di espansioni multipolari
REAL(dbl), INTENT(IN) :: xmin,xmax,ymin,ymax,zmin,zmax				!Coordinate punto nello spazio
REAL(dbl), INTENT(IN) :: lambda							!Lunghezza d'onda incidente
INTEGER(lo), INTENT(IN) :: xstep,ystep,zstep					!Numero di step

COMPLEX(dbl), DIMENSION(:,:,:), INTENT(OUT) :: m_E				!Campo cartesiano
INTEGER(lo) , INTENT(OUT) :: error						!Flag di errore

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,pshift					!Campo calcolato in coordinate sferiche e phase shift c inc
COMPLEX(dbl) :: rhoc,mc							!Rho complesso,indice di rifrazione complesso
REAL(dbl) :: rhor,phase							!Rho reale e fase campo inc
REAL(dbl) :: r,theta,phi 						!Coordinate sferiche locali
REAL(dbl) :: x,y,z,dx,dy,dz						!Coordinate cartesiane locali,delta x,y,z
REAL(dbl) :: d,radius							!Distanza punto-sfere,raggio sfera
REAL(dbl), DIMENSION(1:3) :: v_xyz					!Vettore coordinate da ruotare
INTEGER(lo) :: i,j							!Indici e dimensione
REAL(dbl) :: nr,mr							!Indici reali

!Inizio subroutine vera e propria

!Controllo di avere un piano ben definito
plane_if: IF ((xmin/=xmax) .AND. (ymin/=ymax) .AND. (zmin/=zmax)) THEN
	error=1
	WRITE(*,10) 
	10 FORMAT ("Errore in Efield_sub: nessun piano definito, il programma termina ora...")                           
	RETURN
END IF plane_if

!Controllo di avere i bounds ben definiti
bounds_if: IF ((xmin>xmax) .OR. (ymin>ymax) .OR. (zmin>zmax)) THEN
	error=1
	WRITE(*,20) 
	20 FORMAT ("Errore in Efield_sub: bounds mal definiti, il programma termina ora...")                           
	RETURN
END IF bounds_if

!Controllo di avere gli steps ben definiti
steps_if: IF ((xstep<2) .OR. (ystep<2) .OR. (zstep<2)) THEN
	error=1
	WRITE(*,30) 
	30 FORMAT ("Errore in Efield_sub: steps mal definiti, il programma termina ora...")                           
	RETURN
END IF steps_if

!Calcolo i delta a seconda degli step
dx=(xmax-xmin)/REAL(xstep-1,dbl)
dy=(ymax-ymin)/REAL(ystep-1,dbl)
dz=(zmax-zmin)/REAL(zstep-1,dbl)

!Comincio i conti del campo e un grande if 
x=xmin
y=ymin
z=zmin

xyz_if: IF (xmin==xmax) THEN

	
	y_yz_do: DO i=1,ystep
	
		z=zmin
	
		z_yz_do: DO j=1,zstep

			!Inizializzo i vettori da ruotare
			v_xyz(1)=x
			v_xyz(2)=y
			v_xyz(3)=z

			v_xyz=MATMUL(m_euler,v_xyz)

			m_E(i,j,1)=y
			m_E(i,j,2)=z
			CALL Exyz_ss_sub_dip(v_kinc,v_Einc,flaginc,dip_flag,nstop,&
				&ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,v_amnbmn_rhs,&
				&v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
				& m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

			Exyz1_if: IF (error/=0) THEN
				error=1
				WRITE(*,40) 
			40 FORMAT ("Errore in Exyz_ss_sub_dip chiamata da Efield_sub_dip: il programma termina ora...")
				RETURN
			END IF Exyz1_if

			z=z+dz

		END DO z_yz_do
		
		y=y+dy
		
	END DO y_yz_do

ELSE IF (ymin==ymax) THEN

	x_xz_do: DO i=1,xstep
	
		z=zmin
	
		z_xz_do: DO j=1,zstep

			!Inizializzo i vettori da ruotare
			v_xyz(1)=x
			v_xyz(2)=y
			v_xyz(3)=z

			v_xyz=MATMUL(m_euler,v_xyz)
		
			m_E(i,j,1)=x
			m_E(i,j,2)=z
			CALL Exyz_ss_sub_dip(v_kinc,v_Einc,flaginc,dip_flag,nstop,&
				& ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,v_amnbmn_rhs,&
				& v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
				& m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)

			Exyz2_if: IF (error/=0) THEN
				error=1
				WRITE(*,50) 
			50 FORMAT ("Errore in Exyz_ss_sub_dip chiamata da Efield_sub_dip: il programma termina ora...")
				RETURN
			END IF Exyz2_if

			z=z+dz

		END DO z_xz_do
		
		x=x+dx
		
	END DO x_xz_do

ELSE

	x_xy_do: DO i=1,xstep
	
		y=ymin
	
		y_xy_do: DO j=1,ystep

			!Inizializzo i vettori da ruotare
			v_xyz(1)=x
			v_xyz(2)=y
			v_xyz(3)=z

			v_xyz=MATMUL(m_euler,v_xyz)

			m_E(i,j,1)=x
			m_E(i,j,2)=y
			CALL Exyz_ss_sub_dip(v_kinc,v_Einc,flaginc,dip_flag,nstop,&
				& ratio,lambda,v_xyz(1),v_xyz(2),v_xyz(3),v_amnbmn,v_amnbmn_shell,v_dmncmn_shell,v_dmncmn,v_amnbmn_rhs,&
				& v_emn,m_xyz,m_epseq,v_req,ref_index,v_patt, &
				& m_E(i,j,3),m_E(i,j,4),m_E(i,j,5),error)


			Exyz3_if: IF (error/=0) THEN
				error=1
				WRITE(*,60) 
			60 FORMAT ("Errore in Exyz_ss_sub_dip chiamata da Efield_sub_dip: il programma termina ora...")
				RETURN
			END IF Exyz3_if

			y=y+dy

		END DO y_xy_do
		
		x=x+dx
		
	END DO x_xy_do

END IF xyz_if

END SUBROUTINE Efield_ss_sub_dip

END MODULE local_field
