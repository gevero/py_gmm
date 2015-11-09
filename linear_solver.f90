MODULE linear_solver

USE KINDS
USE SHARED_DATA
USE VEC_TRANS

CONTAINS

!******************************************************************************
!1) SUBROUTINE rot1_sub: faccio la prima rotazione
!******************************************************************************
SUBROUTINE rot1_sub(nstop,next,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                         ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_x1

! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_x1(1:blockside))

!Moltiplico per l'esponenziale
i=0
n_do: DO n=1,nstop
    m_do: DO m=-n,n

        i=i+1
        v_x1(i)=m_exPhi(m,next)*v_x(i)

    END DO m_do
END DO n_do

!Routine di moltiplicazione
CALL amudz(blockside,v_x1,v_y,m_Dij(:,next),m_jDij(:,next),m_iDij(:,next))

DEALLOCATE(v_x1)

END SUBROUTINE rot1_sub



!******************************************************************************
!1) SUBROUTINE rot1_sub: faccio la prima rotazione per il right hand side
!******************************************************************************
SUBROUTINE rot1_sub_rhs(nstop,next,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                         ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_x1

! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_x1(1:blockside))

!Moltiplico per l'esponenziale
i=0
n_do: DO n=1,nstop
    m_do: DO m=-n,n

        i=i+1
        v_x1(i)=m_exPhi_rhs(m,next)*v_x(i)

    END DO m_do
END DO n_do

!Routine di moltiplicazione
CALL amudz(blockside,v_x1,v_y,m_Dij_rhs(:,next),m_jDij_rhs(:,next),m_iDij_rhs(:,next))

DEALLOCATE(v_x1)

END SUBROUTINE rot1_sub_rhs




!******************************************************************************
!2) SUBROUTINE rot2_sub: faccio la seconda rotazione
!******************************************************************************
SUBROUTINE rot2_sub(nstop,next,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x                      ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                         ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,blockside

! Funzione vera e propria
blockside=nstop*(nstop+2)

!Cambio di segno
i=0
n_do: DO n=1,nstop
    m_do: DO m=-n,n

        i=i+1
        v_x(i)=v_oner(m)*v_x(i)

    END DO m_do
END DO n_do

!Routine di moltiplicazione
CALL amudz(blockside,v_x,v_y,m_Dij(:,next),m_jDij(:,next),m_iDij(:,next))

!Cambio di segno
i=0
n_do1: DO n=1,nstop
    m_do1: DO m=-n,n

        i=i+1
        v_y(i)=m_exPhi(-m,next)*v_oner(m)*v_y(i)

    END DO m_do1
END DO n_do1

END SUBROUTINE rot2_sub





!******************************************************************************
!2bis) SUBROUTINE rot2_sub_rhs: faccio la seconda rotazione per il right hand side
!******************************************************************************
SUBROUTINE rot2_sub_rhs(nstop,next,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x                      ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                         ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,blockside

! Funzione vera e propria
blockside=nstop*(nstop+2)

!Cambio di segno
i=0
n_do: DO n=1,nstop
    m_do: DO m=-n,n

        i=i+1
        v_x(i)=((-1.0D0)**m)*v_x(i)

    END DO m_do
END DO n_do

!Routine di moltiplicazione
CALL amudz(blockside,v_x,v_y,m_Dij_rhs(:,next),m_jDij_rhs(:,next),m_iDij_rhs(:,next))

!Cambio di segno
i=0
n_do1: DO n=1,nstop
    m_do1: DO m=-n,n

        i=i+1
        v_y(i)=m_exPhi_rhs(-m,next)*((-1.0D0)**m)*v_y(i)

    END DO m_do1
END DO n_do1

END SUBROUTINE rot2_sub_rhs





!******************************************************************************
!3) SUBROUTINE trasl_sub: faccio la traslazione
!******************************************************************************
SUBROUTINE trasl_sub(nstop,next,v_a,v_b,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_a,v_b                     ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
CALL amuzz(blockside,v_a,v_y1,m_sAij(:,next),m_jABij(:,next),m_iABij(:,next))
CALL amuzz(blockside,v_b,v_y2,m_sBij(:,next),m_jABij(:,next),m_iABij(:,next))

!Calcolo il vettore
v_y=v_y1+v_y2

DEALLOCATE(v_y1,v_y2)

END SUBROUTINE trasl_sub



!******************************************************************************
!3bis) SUBROUTINE trasl_sub_sca: faccio la traslazione
!******************************************************************************
SUBROUTINE trasl_sub_sca(nstop,next,v_a,v_b,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_a,v_b                     ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
CALL amuzz(blockside,v_a,v_y1,m_sAij_sca(:,next),m_jABij(:,next),m_iABij(:,next))
CALL amuzz(blockside,v_b,v_y2,m_sBij_sca(:,next),m_jABij(:,next),m_iABij(:,next))

!Calcolo il vettore
v_y=v_y1+v_y2

DEALLOCATE(v_y1,v_y2)

END SUBROUTINE trasl_sub_sca





!******************************************************************************
!3tris) SUBROUTINE trasl_sub_rhs: faccio la traslazione
!******************************************************************************
SUBROUTINE trasl_sub_rhs(nstop,next,v_a,v_b,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop,next                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_a,v_b                     ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
CALL amuzz(blockside,v_a,v_y1,m_sAij_rhs(:,next),m_jABij_rhs(:,next),m_iABij_rhs(:,next))
CALL amuzz(blockside,v_b,v_y2,m_sBij_rhs(:,next),m_jABij_rhs(:,next),m_iABij_rhs(:,next))

!Calcolo il vettore
v_y=v_y1+v_y2

DEALLOCATE(v_y1,v_y2)

END SUBROUTINE trasl_sub_rhs



!******************************************************************************
!3quater) SUBROUTINE trasl_sub_shell1: faccio la traslazione
!******************************************************************************
SUBROUTINE trasl_sub_shell1(nstop,v_d,v_c,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop				! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_d,v_c	! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y		! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
CALL amuzz(blockside,v_d,v_y1,v_sU1ij,v_jTUij,v_iTUij)
CALL amuzz(blockside,v_c,v_y2,v_sT1ij,v_jTUij,v_iTUij)

!Calcolo il vettore
v_y=v_y1+v_y2

DEALLOCATE(v_y1,v_y2)

END SUBROUTINE trasl_sub_shell1


!******************************************************************************
!3quinties) SUBROUTINE trasl_sub_shell2: faccio la traslazione
!******************************************************************************
SUBROUTINE trasl_sub_shell2(nstop,v_d,v_c,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop				! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_d,v_c	! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y		! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
CALL amuzz(blockside,v_d,v_y1,v_sU2ij,v_jTUij,v_iTUij)
CALL amuzz(blockside,v_c,v_y2,v_sT2ij,v_jTUij,v_iTUij)

!Calcolo il vettore
v_y=v_y1+v_y2

DEALLOCATE(v_y1,v_y2)

END SUBROUTINE trasl_sub_shell2



!******************************************************************************
!3sexsties) SUBROUTINE trasl_sub_shell2: faccio la traslazione
!******************************************************************************
SUBROUTINE trasl_sub_shell(nstop,v_a,v_b,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop				! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_a,v_b	! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y		! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
CALL amuzz(blockside,v_a,v_y1,v_sAij,v_jABij,v_iABij)
CALL amuzz(blockside,v_b,v_y2,v_sBij,v_jABij,v_iABij)

!Calcolo il vettore
v_y=v_y1+v_y2

DEALLOCATE(v_y1,v_y2)

END SUBROUTINE trasl_sub_shell



!******************************************************************************
!4) SUBROUTINE matvec(matrixside,v_x,v_y)
!******************************************************************************
SUBROUTINE matvec(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "matrixside:",matrixside
!WRITE(*,*) "v_x",v_x
!WRITE(*,*) "v_y",v_y

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "v_bin", v_bin
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ns

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Tratto unicamente la porzione j-esima del vettore incognita
    v_a1=v_ain(lowbi:upbi)
    v_b1=v_bin(lowbi:upbi)

!	WRITE(*,*)"Inp: ",dnorm2_bcg(upbi-lowbi,v_a1),dnorm2_bcg(upbi-lowbi,v_b1)

    !Qui sono sulla diagonale, quindi mi basta dividere
    !per i coeff di particella singola
    cursor=0
    diag_ndo: DO n=1,nstop
        diag_mdo: DO m=-n,n

            cursor=cursor+1
            v_a4(cursor)=v_a1(cursor)/m_a(n,v_patt(i))
            v_b4(cursor)=v_b1(cursor)/m_b(n,v_patt(i))

        END DO diag_mdo
    END DO diag_ndo

!	WRITE(*,*)"Diag: ",dnorm2_bcg(upbi-lowbi,v_a4),dnorm2_bcg(upbi-lowbi,v_b4)

    !Faccio la somma per l'update dei vettori
    v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
    v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

!	WRITE(*,*) "v_a4 diag", v_a4
!	WRITE(*,*) "v_b4 diag", v_b4

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock(i),v_iBlock(i+1)-1

        !Indice colonna
        j=v_jBlock(next)

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)


!        WRITE(*,*) "v_a1", v_a1
!	WRITE(*,*) "v_b1", v_b1

        !Qui ho la moltiplicazione fattorizzata

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

!        WRITE(*,*) "v_a2", v_a2
!	WRITE(*,*) "v_b2", v_b2
!		WRITE(*,*)"rot1: ",dnorm2_bcg(upb-lowb,v_a2),dnorm2_bcg(upbi-lowbi,v_b2)

        !Traslazione
        CALL trasl_sub(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub(nstop,next,v_b2,v_a2,v_b3)

!        WRITE(*,*) "v_a3", v_a3
!	WRITE(*,*) "v_b3", v_b3
!		WRITE(*,*)"trasl: ",dnorm2_bcg(upb-lowb,v_a3),dnorm2_bcg(upbi-lowbi,v_b3)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

!        WRITE(*,*) "v_a4", v_a4
!	WRITE(*,*) "v_b4", v_b4
!		WRITE(*,*)"rot2: ",dnorm2_bcg(upb-lowb,v_a4),dnorm2_bcg(upbi-lowbi,v_b4)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec





!******************************************************************************
!4bis) SUBROUTINE matvec_sca(matrixside,v_x,v_y)
!******************************************************************************
SUBROUTINE matvec_sca(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ns

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Tratto unicamente la porzione j-esima del vettore incognita
    v_a1=v_ain(lowbi:upbi)
    v_b1=v_bin(lowbi:upbi)


    !Qui sono sulla diagonale, quindi mi basta dividere
    !per i coeff di particella singola
    cursor=0
    diag_ndo: DO n=1,nstop
        diag_mdo: DO m=-n,n

            cursor=cursor+1
            v_a4(cursor)=v_a1(cursor)
            v_b4(cursor)=v_b1(cursor)

        END DO diag_mdo
    END DO diag_ndo

    !Faccio la somma per l'update dei vettori
    v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
    v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock(i),v_iBlock(i+1)-1

        !Indice colonna
        j=v_jBlock(next)

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

        !Traslazione
        CALL trasl_sub_sca(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub_sca(nstop,next,v_b2,v_a2,v_b3)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_sca






!******************************************************************************
!4bis) SUBROUTINE matvec_sca(matrixside,v_x,v_y)
!******************************************************************************
SUBROUTINE matvec_sca_dip(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ns

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Tratto unicamente la porzione j-esima del vettore incognita
    v_a1=v_ain(lowbi:upbi)
    v_b1=v_bin(lowbi:upbi)


    !Qui sono sulla diagonale, quindi mi basta dividere
    !per i coeff di particella singola
    cursor=0
    diag_ndo: DO n=1,nstop
        diag_mdo: DO m=-n,n

            cursor=cursor+1
            v_a4(cursor)=v_a1(cursor)
            v_b4(cursor)=v_b1(cursor)

        END DO diag_mdo
    END DO diag_ndo

    !Faccio la somma per l'update dei vettori
    v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
    v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock_sca(i),v_iBlock_sca(i+1)-1

        !Indice colonna
        j=v_jBlock_sca(next)

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

        !Traslazione
        CALL trasl_sub_sca(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub_sca(nstop,next,v_b2,v_a2,v_b3)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_sca_dip






!******************************************************************************
!4tris) SUBROUTINE matvec_rhs(matrixside,v_x,v_y): per il right hand side
!******************************************************************************
SUBROUTINE matvec_rhs(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ns

!    WRITE(*,*) "ns:",i
!    WRITE(*,*)

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock_rhs(i),v_iBlock_rhs(i+1)-1

        !Indice colonna
        j=v_jBlock_rhs(next)

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

        WRITE(*,*) "v_a1"
        WRITE(*,*) (v_a1(n),n=lowb,upb)
        WRITE(*,*)

        !Prima rotazione
        CALL rot1_sub_rhs(nstop,next,v_a1,v_a2)
        CALL rot1_sub_rhs(nstop,next,v_b1,v_b2)

        WRITE(*,*) "v_a2"
        WRITE(*,*) (v_a2(n),n=lowb,upb)
        WRITE(*,*)

        !Traslazione
        CALL trasl_sub_rhs(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub_rhs(nstop,next,v_b2,v_a2,v_b3)

        WRITE(*,*) "v_a3"
        WRITE(*,*) (v_a3(n),n=lowb,upb)
        WRITE(*,*)

        !Rotazione all'indietro
        CALL rot2_sub_rhs(nstop,next,v_a3,v_a4)
        CALL rot2_sub_rhs(nstop,next,v_b3,v_b4)

        WRITE(*,*) "v_a4"
        WRITE(*,*) (v_a4(n),n=lowb,upb)
        WRITE(*,*)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

WRITE(*,*)"v_aout"
WRITE(*,*) (v_aout(n),n=1,matrixside/2)
WRITE(*,*)

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_rhs





!******************************************************************************
!4quadris) SUBROUTINE matvec_rhs(matrixside,v_x,v_y): per il right hand side
!******************************************************************************
SUBROUTINE matvec_rhs_dip(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,(ns-ndip)

!    WRITE(*,*) "ns:",i
!    WRITE(*,*)

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock(i),v_iBlock(i+1)-1

        !Indice colonna
        j=v_jBlock(next)

        IF (j<=(ns-ndip)) CYCLE next_do

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

!        WRITE(*,*) "v_a1"
!        WRITE(*,*) (v_a1(n),n=lowb,upb)
!        WRITE(*,*)

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

!        WRITE(*,*) "v_a2"
!        WRITE(*,*) (v_a2(n),n=lowb,upb)
!        WRITE(*,*)

        !Traslazione
        CALL trasl_sub(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub(nstop,next,v_b2,v_a2,v_b3)

!        WRITE(*,*) "v_a3"
!        WRITE(*,*) (v_a3(n),n=lowb,upb)
!        WRITE(*,*)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

!        WRITE(*,*) "v_a4"
!        WRITE(*,*) (v_a4(n),n=lowb,upb)
!        WRITE(*,*)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!WRITE(*,*)"v_aout"
!WRITE(*,*) (v_aout(n),n=1,matrixside/2)
!WRITE(*,*)

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_rhs_dip



!******************************************************************************
!4quinties) SUBROUTINE matvec_sca_dipant(matrixside,v_x,v_y), per sistemare
!la sezione di scattering solo per l'antenna
!******************************************************************************
SUBROUTINE matvec_sca_dipant(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ns-ndip

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Tratto unicamente la porzione j-esima del vettore incognita
    v_a1=v_ain(lowbi:upbi)
    v_b1=v_bin(lowbi:upbi)


    !Qui sono sulla diagonale, quindi mi basta dividere
    !per i coeff di particella singola
    cursor=0
    diag_ndo: DO n=1,nstop
        diag_mdo: DO m=-n,n

            cursor=cursor+1
            v_a4(cursor)=v_a1(cursor)
            v_b4(cursor)=v_b1(cursor)

        END DO diag_mdo
    END DO diag_ndo

    !Faccio la somma per l'update dei vettori
    v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
    v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock_sca(i),v_iBlock_sca(i+1)-1

        !Indice colonna
        j=v_jBlock_sca(next)

        IF (j>(ns-ndip)) CYCLE next_do

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

        !Traslazione
        CALL trasl_sub_sca(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub_sca(nstop,next,v_b2,v_a2,v_b3)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_sca_dipant



!******************************************************************************
!4sexties) SUBROUTINE matvec_dip(matrixside,v_x,v_y)
!******************************************************************************
SUBROUTINE matvec_dip(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ns-ndip

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Tratto unicamente la porzione j-esima del vettore incognita
    v_a1=v_ain(lowbi:upbi)
    v_b1=v_bin(lowbi:upbi)

!	WRITE(*,*)"Inp: ",dnorm2_bcg(upbi-lowbi,v_a1),dnorm2_bcg(upbi-lowbi,v_b1)

    !Qui sono sulla diagonale, quindi mi basta dividere
    !per i coeff di particella singola
    cursor=0
    diag_ndo: DO n=1,nstop
        diag_mdo: DO m=-n,n

            cursor=cursor+1
            v_a4(cursor)=v_a1(cursor)/m_a(n,v_patt(i))
            v_b4(cursor)=v_b1(cursor)/m_b(n,v_patt(i))

        END DO diag_mdo
    END DO diag_ndo


!	WRITE(*,*)"Diag: ",dnorm2_bcg(upbi-lowbi,v_a4),dnorm2_bcg(upbi-lowbi,v_b4)

    !Faccio la somma per l'update dei vettori
    v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
    v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock(i),v_iBlock(i+1)-1

        !Indice colonna
        j=v_jBlock(next)

        IF ((j>(ns-ndip)) .OR. ((ns-ndip)==1)) CYCLE next_do

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

!		WRITE(*,*)"rot1: ",dnorm2_bcg(upb-lowb,v_a2),dnorm2_bcg(upbi-lowbi,v_b2)

        !Traslazione
        CALL trasl_sub(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub(nstop,next,v_b2,v_a2,v_b3)

!		WRITE(*,*)"trasl: ",dnorm2_bcg(upb-lowb,v_a3),dnorm2_bcg(upbi-lowbi,v_b3)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

!		WRITE(*,*)"rot2: ",dnorm2_bcg(upb-lowb,v_a4),dnorm2_bcg(upbi-lowbi,v_b4)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_dip



!******************************************************************************
!4septies) SUBROUTINE matvec_shell(matrixside,v_x,v_y)
!******************************************************************************
SUBROUTINE matvec_shell(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_din,v_cin,v_dout,v_cout

! Funzione vera e propria
blockside=nstop*(nstop+2)

! Alloco i vettori
ALLOCATE(v_din(1:blockside),v_cin(1:blockside),v_dout(1:blockside),v_cout(1:blockside))

!Scrivo i miei due vettori in entrata
v_din=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_cin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_dout=(0.0D0,0.0D0)
v_cout=(0.0D0,0.0D0)

!Comincio la moltiplicazione

!Traslazione
CALL trasl_sub_shell2(nstop,v_din,v_cin,v_dout)
CALL trasl_sub_shell1(nstop,v_din,v_cin,v_cout)


!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_dout
v_y(2:matrixside:2)=v_cout

END SUBROUTINE matvec_shell



!******************************************************************************
!4septies) SUBROUTINE matvec_trasl_shell(matrixside,v_x,v_y), per traslare i
!coefficienti das un sistema di riferimento all'altro all'interno della shell
!******************************************************************************
SUBROUTINE matvec_trasl_shell(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout

! Funzione vera e propria
blockside=nstop*(nstop+2)

! Alloco i vettori
ALLOCATE(v_ain(1:blockside),v_bin(1:blockside),v_aout(1:blockside),v_bout(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)

!Comincio la moltiplicazione

!WRITE(*,*) "lower bounds ain",LBOUND(v_ain, 1)
!WRITE(*,*) "upper bounds ain",UBOUND(v_ain, 1)
!WRITE(*,*) "lower bounds bin",LBOUND(v_bin, 1)
!WRITE(*,*) "upper bounds bin",UBOUND(v_bin, 1)

!Traslazione
CALL trasl_sub_shell(nstop,v_ain,v_bin,v_aout)
CALL trasl_sub_shell(nstop,v_bin,v_ain,v_bout)

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_trasl_shell



!******************************************************************************
!5) SUBROUTINE precond(matrixside,m_x)
!******************************************************************************
SUBROUTINE precond(matrixside,v_x)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x                         ! Vettore input

! Dichiarazione variabili interne
INTEGER(lo) :: n,m,j,i

! Funzione vera e propria
j=0
ns_do: DO i=1,ns
    n_do: DO n=1,nstop
        m_do: DO m=-n,n

        j=j+1

        v_x(2*j-1)=v_x(2*j-1)*m_a(n,v_patt(i))
        v_x(2*j)=v_x(2*j)*m_b(n,v_patt(i))

        END DO m_do
    END DO n_do
END DO ns_do

END SUBROUTINE precond




!******************************************************************************
!5bis) SUBROUTINE precond(matrixside,m_x)
!******************************************************************************
SUBROUTINE precond_dip(matrixside,v_x)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x                         ! Vettore input

! Dichiarazione variabili interne
INTEGER(lo) :: n,m,j,i

! Funzione vera e propria
j=0
ns_do: DO i=1,(ns-ndip)
    n_do: DO n=1,nstop
        m_do: DO m=-n,n

        j=j+1

        v_x(2*j-1)=v_x(2*j-1)*m_a(n,v_patt(i))
        v_x(2*j)=v_x(2*j)*m_b(n,v_patt(i))

        END DO m_do
    END DO n_do
END DO ns_do

END SUBROUTINE precond_dip


!******************************************************************************
!5bis) SUBROUTINE precond(matrixside,m_x)
!******************************************************************************
SUBROUTINE precond_shell(matrixside,v_x)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x                         ! Vettore input

! Dichiarazione variabili interne
INTEGER(lo) :: n,m,j,i

! Funzione vera e propria
j=0
n_do: DO n=1,nstop
	m_do: DO m=-n,n

		j=j+1

		v_x(2*j-1)=v_precond_shell(2*j-1)*v_x(2*j-1)
		v_x(2*j)=v_precond_shell(2*j)*v_x(2*j)

	END DO m_do
END DO n_do


END SUBROUTINE precond_shell



!******************************************************************************
!6) zbcg2(print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                  mxmatvec,work,info)
!******************************************************************************
subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,toler, &
                  mxmatvec,work,info)

! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)
!
! Improved "vanilla" BiCGStab(2) iterative method
!
! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on original work of D.R.Fokkema:
!
! subroutine zbistbl v1.1 1998
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! }}
!
! Your bug reports, comments, etc. are welcome:
! m.a.botchev@math.utwente.nl
!
! ------------------------------
! Description of the parameters:
! ------------------------------
!
! print_resid (input) LOGICAL. If print_resid=.true. the number of
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve
!
! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not.
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix
!            vector multiplications allowed to be done.  On output:
!            if info>=0, mxmatvec is set to the actual number of matrix
!            vector multiplications done
!
! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
implicit none

! Parameters:
LOGICAL, INTENT(in)   :: print_resid,nonzero_x
INTEGER(lo), INTENT(in)   :: l, n
INTEGER(lo), INTENT(inout):: mxmatvec
INTEGER(lo), INTENT(out)  :: info
COMPLEX(dbl), DIMENSION(:) ,INTENT(inout):: x(:)
COMPLEX(dbl), DIMENSION(:) ,INTENT(in)   :: rhs(:)
REAL(dbl), INTENT(inout) :: toler
COMPLEX(dbl), DIMENSION(:,:), INTENT(out)  :: work(:,:)

! Local variables:
COMPLEX(dbl), DIMENSION(1:l+1) :: y0,yl,zy0,zyl
COMPLEX(dbl), DIMENSION(1:l+1,1:l+1)  :: matrix_z
LOGICAL    :: rcmp, xpdt
INTEGER(lo)    :: i, j, k, nmatvec,ii,jj
COMPLEX(dbl) :: alpha, beta, omega, rho0, rho1, sigma
COMPLEX(dbl) :: varrho, hatgamma
REAL(dbl)    :: rnrm0, rnrm
REAL(dbl)    :: mxnrmx, mxnrmr
COMPLEX(dbl) :: kappa0, kappal

! Aliases for the parts of the work array:
INTEGER(lo)          :: rr, r, u, xp, bp

! Constants:
REAL(dbl),    parameter :: delta = 1d-2
COMPLEX(dbl), parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

    ! Functions:
    !REAL(sgl)     :: dnorm2_bcg
    !COMPLEX(dbl) :: zdot_bcg


info = 0

if (l<1 .or. l>2) info = -2
if (n<1) info = -3
if (toler<=0d0) info = -9
if (mxmatvec<0) info = -10

rr = 1
r = rr+1
u = r+(l+1)
xp = u+(l+1)
bp = xp+1

if (info/=0) return

!!~~~~~~~~~~~~~~~~~~~~
!WRITE(*,*) 'work,n1'
!Do,jj=1,9
!Do ii=1,n
!WRITE(*,300) work(ii,jj)
!300 FORMAT (2ES21.12)
!end do
!end do
!WRITE(*,*)
!!~~~~~~~~~~~~~~~~~~~
! Initialize first residual
nmatvec=0
!WRITE(*,*) 'x',x
if (nonzero_x) then
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n,x)
   call matvec (n, x, work(1:n,r))

!   !~~~~~~~~~~~~~~~~~~~~
!   WRITE(*,*) 'work,n2'
!   Do,jj=1,9
!   Do ii=1,n
!   WRITE(*,400) work(ii,jj)
!   400 FORMAT (2ES21.12)
!   end do
!   end do
!   WRITE(*,*)
!   if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
!   !~~~~~~~~~~~~~~~~~~~~

   work(1:n,r) = rhs - work(1:n,r)
   nmatvec = 1

!   !~~~~~~~~~~~~~~~~~~~~
!   WRITE(*,*) 'rhs'
!   Do ii=1,n
!   WRITE(*,450) rhs(ii)
!   450 FORMAT (2ES21.12)
!   end do
!   WRITE(*,*)
!   !~~~~~~~~~~~~~~~~~~~~

!   !~~~~~~~~~~~~~~~~~~~~
!   WRITE(*,*) 'work,n3'
!   Do,jj=1,9
!   Do ii=1,n
!   WRITE(*,500) work(ii,jj)
!   500 FORMAT (2ES21.12)
!   end do
!   end do
!   WRITE(*,*)
!   if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
!   !~~~~~~~~~~~~~~~~~~~~

else
   work(1:n,r) = rhs
   nmatvec = 0
end if
call precond (n, work(1:n,r))

! Initialize iteration loop

work(1:n,rr) = work(1:n,r)
work(1:n,bp) = work(1:n,r)
work(1:n,xp) = x
x = zzero

rnrm0 = dnorm2_bcg (n, work(1:n,r))
rnrm = rnrm0

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zzero
omega = zone
sigma = zone
rho0  = zone

! Iterate



do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)

! =====================
! The BiCG part ---
! =====================

   rho0 = -omega*rho0
   do k=1,l
!      WRITE(*,*) 'work(1:n,rr),k',work(1:n,rr),k
!      WRITE(*,*) 'work(1:n,r+k-1),k',work(1:n,r+k-1),k
      rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
!      WRITE(*,*) 'rho1,k',rho1,k
      if (rho0.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
      enddo
      call matvec (n, work(1:n,u+k-1), work(1:n,u+k))
      call precond (n, work(1:n,u+k))
      nmatvec = nmatvec+1
!     if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))

      sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
      if (sigma.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      alpha = rho1/sigma
      x(1:n) = alpha*work(1:n,u) + x(1:n)
      do j=0,k-1
         work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
      enddo
      call matvec (n, work(1:n,r+k-1), work(1:n,r+k))
      call precond (n, work(1:n,r+k))
      nmatvec = nmatvec+1
!     if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
      rnrm = dnorm2_bcg (n, work(1:n,r))
      mxnrmx = max (mxnrmx, rnrm)
      mxnrmr = max (mxnrmr, rnrm)
   enddo

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

   do i=1,l+1
      do j=1,i
         matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
      end do
   end do

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
   do j=2,l+1
      matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
   end do

!  small vectors y0 and yl

y0(1) = -zone
y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
y0(l+1) = zzero

yl(1) = zzero
yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
yl(l+1) = -zone

!  --- Convex combination

! compute Z*y0 and Z*yl
zy0 = zzero
zyl = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
   zyl = zyl + matrix_z(:,j)*yl(j)
end do

kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)

hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 )

y0 = y0 - (hatgamma*kappa0/kappal)*yl


!  --- Update

omega = y0(l+1)

do j=1,l
   work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
   x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
   work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
enddo

! y0 has changed; compute Z*y0 once more
zy0 = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
end do

rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )

! ================================
! The reliable update part ---
! ================================

mxnrmx = max (mxnrmx, rnrm)
mxnrmr = max (mxnrmr, rnrm)
xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
if (rcmp) then
   call matvec (n, x, work(1:n,r))
   call precond (n, work(1:n,r))
   nmatvec = nmatvec + 1
   work(1:n,r) =  work(1:n,bp) - work(1:n,r)
   mxnrmr = rnrm
!  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
   if (xpdt) then

      work(1:n,xp) = x(1:n) + work(1:n,xp)
      x = zzero
      work(1:n,bp) = work(1:n,r)

      mxnrmx = rnrm
   endif
endif

if (print_resid) print *,nmatvec,' ',rnrm

enddo

! =========================
! End of iterations ---
! =========================

x(1:n) = x(1:n) + work(1:n,xp)

if (rnrm>toler*rnrm0) info = 1

! compute the true residual:

! --------------------- One matvec can be saved by commenting out this:
!$$$$$$ call matvec (n, x, work(1:n,r) )
!$$$$$$ work(1:n,r) = rhs(1:n) - work(1:n,r)
!$$$$$$ call precond (n, work(1:n,r))
!$$$$$$ rnrm = dnorm2_bcg(n,work(1:n,r))
!$$$$$$ nmatvec = nmatvec+1
! --------------------- One matvec can be saved by commenting out this^

!$$$$$$ toler = rnrm/rnrm0
!$$$$$$ mxmatvec = nmatvec

end subroutine zbcg2






!******************************************************************************
!6bis) zbcg2(print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                  mxmatvec,work,info)
!******************************************************************************
subroutine zbcg2_dip (print_resid,l,n,x,nonzero_x,rhs,toler, &
                  mxmatvec,work,info)

! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)
!
! Improved "vanilla" BiCGStab(2) iterative method
!
! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on original work of D.R.Fokkema:
!
! subroutine zbistbl v1.1 1998
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! }}
!
! Your bug reports, comments, etc. are welcome:
! m.a.botchev@math.utwente.nl
!
! ------------------------------
! Description of the parameters:
! ------------------------------
!
! print_resid (input) LOGICAL. If print_resid=.true. the number of
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve
!
! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not.
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix
!            vector multiplications allowed to be done.  On output:
!            if info>=0, mxmatvec is set to the actual number of matrix
!            vector multiplications done
!
! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
implicit none

! Parameters:
LOGICAL, INTENT(in)   :: print_resid,nonzero_x
INTEGER(lo), INTENT(in)   :: l, n
INTEGER(lo), INTENT(inout):: mxmatvec
INTEGER(lo), INTENT(out)  :: info
COMPLEX(dbl), DIMENSION(:) ,INTENT(inout):: x(:)
COMPLEX(dbl), DIMENSION(:) ,INTENT(in)   :: rhs(:)
REAL(dbl), INTENT(inout) :: toler
COMPLEX(dbl), DIMENSION(:,:), INTENT(out)  :: work(:,:)

! Local variables:
COMPLEX(dbl), DIMENSION(1:l+1) :: y0,yl,zy0,zyl
COMPLEX(dbl), DIMENSION(1:l+1,1:l+1)  :: matrix_z
LOGICAL    :: rcmp, xpdt
INTEGER(lo)    :: i, j, k, nmatvec
COMPLEX(dbl) :: alpha, beta, omega, rho0, rho1, sigma
COMPLEX(dbl) :: varrho, hatgamma
REAL(dbl)    :: rnrm0, rnrm
REAL(dbl)    :: mxnrmx, mxnrmr
COMPLEX(dbl) :: kappa0, kappal

! Aliases for the parts of the work array:
INTEGER(lo)          :: rr, r, u, xp, bp

! Constants:
REAL(dbl),    parameter :: delta = 1d-2
COMPLEX(dbl), parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

    ! Functions:
    !REAL(sgl)     :: dnorm2_bcg
    !COMPLEX(dbl) :: zdot_bcg


info = 0

if (l<1 .or. l>2) info = -2
if (n<1) info = -3
if (toler<=0d0) info = -9
if (mxmatvec<0) info = -10

rr = 1
r = rr+1
u = r+(l+1)
xp = u+(l+1)
bp = xp+1

if (info/=0) return

! Initialize first residual
nmatvec=0


if (nonzero_x) then
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n,x)
   call matvec_dip (n, x, work(1:n,r))
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
   work(1:n,r) = rhs - work(1:n,r)
   nmatvec = 1
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
else
   work(1:n,r) = rhs
   nmatvec = 0
end if
call precond_dip (n, work(1:n,r))

! Initialize iteration loop

work(1:n,rr) = work(1:n,r)
work(1:n,bp) = work(1:n,r)
work(1:n,xp) = x
x = zzero

rnrm0 = dnorm2_bcg (n, work(1:n,r))
rnrm = rnrm0

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zzero
omega = zone
sigma = zone
rho0  = zone

! Iterate



do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)

! =====================
! The BiCG part ---
! =====================

   rho0 = -omega*rho0
   do k=1,l
      rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
      if (rho0.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
      enddo
      call matvec_dip (n, work(1:n,u+k-1), work(1:n,u+k))
      call precond_dip (n, work(1:n,u+k))
      nmatvec = nmatvec+1
     if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))

      sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
      if (sigma.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      alpha = rho1/sigma
      x(1:n) = alpha*work(1:n,u) + x(1:n)
      do j=0,k-1
         work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
      enddo
      call matvec_dip (n, work(1:n,r+k-1), work(1:n,r+k))
      call precond_dip (n, work(1:n,r+k))
      nmatvec = nmatvec+1
     if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
      rnrm = dnorm2_bcg (n, work(1:n,r))
      mxnrmx = max (mxnrmx, rnrm)
      mxnrmr = max (mxnrmr, rnrm)
   enddo

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

   do i=1,l+1
      do j=1,i
         matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
      end do
   end do

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
   do j=2,l+1
      matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
   end do

!  small vectors y0 and yl

y0(1) = -zone
y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
y0(l+1) = zzero

yl(1) = zzero
yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
yl(l+1) = -zone

!  --- Convex combination

! compute Z*y0 and Z*yl
zy0 = zzero
zyl = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
   zyl = zyl + matrix_z(:,j)*yl(j)
end do

kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)

hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 )

y0 = y0 - (hatgamma*kappa0/kappal)*yl


!  --- Update

omega = y0(l+1)

do j=1,l
   work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
   x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
   work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
enddo

! y0 has changed; compute Z*y0 once more
zy0 = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
end do

rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )

! ================================
! The reliable update part ---
! ================================

mxnrmx = max (mxnrmx, rnrm)
mxnrmr = max (mxnrmr, rnrm)
xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
if (rcmp) then
   call matvec_dip (n, x, work(1:n,r))
   call precond_dip (n, work(1:n,r))
   nmatvec = nmatvec + 1
   work(1:n,r) =  work(1:n,bp) - work(1:n,r)
   mxnrmr = rnrm
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
   if (xpdt) then

      work(1:n,xp) = x(1:n) + work(1:n,xp)
      x = zzero
      work(1:n,bp) = work(1:n,r)

      mxnrmx = rnrm
   endif
endif

if (print_resid) print *,nmatvec,' ',rnrm

enddo

! =========================
! End of iterations ---
! =========================

x(1:n) = x(1:n) + work(1:n,xp)

if (rnrm>toler*rnrm0) info = 1

! compute the true residual:

! --------------------- One matvec can be saved by commenting out this:
!$$$$$$ call matvec_dip (n, x, work(1:n,r) )
!$$$$$$ work(1:n,r) = rhs(1:n) - work(1:n,r)
!$$$$$$ call precond (n, work(1:n,r))
!$$$$$$ rnrm = dnorm2_bcg(n,work(1:n,r))
!$$$$$$ nmatvec = nmatvec+1
! --------------------- One matvec can be saved by commenting out this^

!$$$$$$ toler = rnrm/rnrm0
!$$$$$$ mxmatvec = nmatvec

end subroutine zbcg2_dip






!******************************************************************************
!6bis) zbcg2(print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                  mxmatvec,work,info)
!******************************************************************************
subroutine zbcg2_shell (print_resid,l,n,x,nonzero_x,rhs,toler, &
                  mxmatvec,work,info)

! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)
!
! Improved "vanilla" BiCGStab(2) iterative method
!
! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on original work of D.R.Fokkema:
!
! subroutine zbistbl v1.1 1998
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! }}
!
! Your bug reports, comments, etc. are welcome:
! m.a.botchev@math.utwente.nl
!
! ------------------------------
! Description of the parameters:
! ------------------------------
!
! print_resid (input) LOGICAL. If print_resid=.true. the number of
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve
!
! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not.
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix
!            vector multiplications allowed to be done.  On output:
!            if info>=0, mxmatvec is set to the actual number of matrix
!            vector multiplications done
!
! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
implicit none

! Parameters:
LOGICAL, INTENT(in)   :: print_resid,nonzero_x
INTEGER(lo), INTENT(in)   :: l, n
INTEGER(lo), INTENT(inout):: mxmatvec
INTEGER(lo), INTENT(out)  :: info
COMPLEX(dbl), DIMENSION(:) ,INTENT(inout):: x(:)
COMPLEX(dbl), DIMENSION(:) ,INTENT(in)   :: rhs(:)
REAL(dbl), INTENT(inout) :: toler
COMPLEX(dbl), DIMENSION(:,:), INTENT(out)  :: work(:,:)

! Local variables:
COMPLEX(dbl), DIMENSION(1:l+1) :: y0,yl,zy0,zyl
COMPLEX(dbl), DIMENSION(1:l+1,1:l+1)  :: matrix_z
LOGICAL    :: rcmp, xpdt
INTEGER(lo)    :: i, j, k, nmatvec
COMPLEX(dbl) :: alpha, beta, omega, rho0, rho1, sigma
COMPLEX(dbl) :: varrho, hatgamma
REAL(dbl)    :: rnrm0, rnrm
REAL(dbl)    :: mxnrmx, mxnrmr
COMPLEX(dbl) :: kappa0, kappal

! Aliases for the parts of the work array:
INTEGER(lo)          :: rr, r, u, xp, bp

! Constants:
REAL(dbl),    parameter :: delta = 1d-2
COMPLEX(dbl), parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

    ! Functions:
    !REAL(sgl)     :: dnorm2_bcg
    !COMPLEX(dbl) :: zdot_bcg


info = 0

if (l<1 .or. l>2) info = -2
if (n<1) info = -3
if (toler<=0d0) info = -9
if (mxmatvec<0) info = -10

rr = 1
r = rr+1
u = r+(l+1)
xp = u+(l+1)
bp = xp+1

if (info/=0) return

! Initialize first residual
nmatvec=0


if (nonzero_x) then
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n,x)
   call matvec_shell (n, x, work(1:n,r))
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
   work(1:n,r) = rhs - work(1:n,r)
   nmatvec = 1
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
else
   work(1:n,r) = rhs
   nmatvec = 0
end if
call precond_shell (n, work(1:n,r))

! Initialize iteration loop

work(1:n,rr) = work(1:n,r)
work(1:n,bp) = work(1:n,r)
work(1:n,xp) = x
x = zzero

rnrm0 = dnorm2_bcg (n, work(1:n,r))
rnrm = rnrm0

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zzero
omega = zone
sigma = zone
rho0  = zone

! Iterate



do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)

! =====================
! The BiCG part ---
! =====================

   rho0 = -omega*rho0
   do k=1,l
      rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
      if (rho0.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
      enddo
      call matvec_shell (n, work(1:n,u+k-1), work(1:n,u+k))
      call precond_shell (n, work(1:n,u+k))
      nmatvec = nmatvec+1
     if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))

      sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
      if (sigma.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      alpha = rho1/sigma
      x(1:n) = alpha*work(1:n,u) + x(1:n)
      do j=0,k-1
         work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
      enddo
      call matvec_shell (n, work(1:n,r+k-1), work(1:n,r+k))
      call precond_shell (n, work(1:n,r+k))
      nmatvec = nmatvec+1
     if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
      rnrm = dnorm2_bcg (n, work(1:n,r))
      mxnrmx = max (mxnrmx, rnrm)
      mxnrmr = max (mxnrmr, rnrm)
   enddo

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

   do i=1,l+1
      do j=1,i
         matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
      end do
   end do

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
   do j=2,l+1
      matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
   end do

!  small vectors y0 and yl

y0(1) = -zone
y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
y0(l+1) = zzero

yl(1) = zzero
yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
yl(l+1) = -zone

!  --- Convex combination

! compute Z*y0 and Z*yl
zy0 = zzero
zyl = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
   zyl = zyl + matrix_z(:,j)*yl(j)
end do

kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)

hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 )

y0 = y0 - (hatgamma*kappa0/kappal)*yl


!  --- Update

omega = y0(l+1)

do j=1,l
   work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
   x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
   work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
enddo

! y0 has changed; compute Z*y0 once more
zy0 = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
end do

rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )

! ================================
! The reliable update part ---
! ================================

mxnrmx = max (mxnrmx, rnrm)
mxnrmr = max (mxnrmr, rnrm)
xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
if (rcmp) then
   call matvec_shell (n, x, work(1:n,r))
   call precond_shell (n, work(1:n,r))
   nmatvec = nmatvec + 1
   work(1:n,r) =  work(1:n,bp) - work(1:n,r)
   mxnrmr = rnrm
  if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
   if (xpdt) then

      work(1:n,xp) = x(1:n) + work(1:n,xp)
      x = zzero
      work(1:n,bp) = work(1:n,r)

      mxnrmx = rnrm
   endif
endif

if (print_resid) print *,nmatvec,' ',rnrm

enddo

! =========================
! End of iterations ---
! =========================

x(1:n) = x(1:n) + work(1:n,xp)

if (rnrm>toler*rnrm0) info = 1

! compute the true residual:

! --------------------- One matvec can be saved by commenting out this:
!$$$$$$ call matvec_dip (n, x, work(1:n,r) )
!$$$$$$ work(1:n,r) = rhs(1:n) - work(1:n,r)
!$$$$$$ call precond (n, work(1:n,r))
!$$$$$$ rnrm = dnorm2_bcg(n,work(1:n,r))
!$$$$$$ nmatvec = nmatvec+1
! --------------------- One matvec can be saved by commenting out this^

!$$$$$$ toler = rnrm/rnrm0
!$$$$$$ mxmatvec = nmatvec

end subroutine zbcg2_shell






!$$$$$$ !*************************************************************
!$$$$$$ !SUBROUTINE zbistbl (l, n, x, b, mv, solve, tol, &
!$$$$$$ !                 & mxmv, work, ldw, rwork, ldrw, iwork, info)
!$$$$$$ !*************************************************************
!$$$$$$
!$$$$$$ SUBROUTINE zbistbl (l,n,x,b,tol,mxmv,work,ldw,rwork,ldrw,iwork,info)
!$$$$$$
!$$$$$$ ! subroutine zbistbl v1.1 1998
!$$$$$$ !
!$$$$$$ ! Copyright (c) 1995-1998 by D.R. Fokkema.
!$$$$$$ ! Permission to copy all or part of this work is granted,
!$$$$$$ ! provided that the copies are not made or distributed
!$$$$$$ ! for resale, and that the copyright notice and this
!$$$$$$ ! notice are retained.
!$$$$$$ !
!$$$$$$ ! THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE AUTHOR
!$$$$$$ ! PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESSED OR IMPLIED,
!$$$$$$ ! REGARDING THE WORK, INCLUDING WARRANTIES WITH RESPECT TO ITS
!$$$$$$ ! MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE.
!$$$$$$ !
!$$$$$$ ! Parameter list (by M.Botchev, October, 1997):
!$$$$$$ !
!$$$$$$ ! l      == (input) INTEGER BiCGstab's dimension
!$$$$$$ ! n      == (input) INTEGER size of the system to solve
!$$$$$$ ! x      == (input/output) DOUBLE PRECISION array dimension n
!$$$$$$ !           initial guess on input, solution on output
!$$$$$$ ! b      == (input) DOUBLE PRECISION array dimension n
!$$$$$$ !           right-hand side (rhs) vector
!$$$$$$ ! mv     == (input) EXTERNAL name of matrix vector subroutine
!$$$$$$ !           to deliver y:=A*x by CALL mv(n,x,y)
!$$$$$$ ! solve  == (input) EXTERNAL name of subroutine to perform
!$$$$$$ !           a preconditioner solve by CALL solve(n,x). This call
!$$$$$$ !           delivers solution to M*x=x, M is the left preconditioning
!$$$$$$ !           matrix, x is rhs on input and solution on output.
!$$$$$$ ! tol    == (input) DOUBLE PRECISION tolerance for the relative
!$$$$$$ !           residual norm stopping criterion: ||r_k|| <= tol*||r_0||
!$$$$$$ ! mxmv   == (input/output) INTEGER.  On input: maximum number of matrix
!$$$$$$ !           vector multiplications allowed to be done.  On input:
!$$$$$$ !           actual number of matrix vector multiplications done
!$$$$$$ ! work   == (workspace) DOUBLE PRECISION array dimension (n,3+2*(l+1))
!$$$$$$ ! ldw    == (input) INTEGER leading dimension of work (??)
!$$$$$$ ! rwork  == (workspace) DOUBLE PRECISION array dimension (l+1,3+2*(l+1))
!$$$$$$ ! ldrw   == (input) INTEGER leading dimension of rwork (??)
!$$$$$$ ! iwork  == (workspace) INTEGER array dimension (l+1)
!$$$$$$ ! info   == (output) INTEGER.  Info = 0 in case of normal computations
!$$$$$$ !           and
!$$$$$$ !           info = m < 0 - means paramater number m has an illegal value
!$$$$$$ !           info = 1 - means no convergence achieved (stopping criterion
!$$$$$$ !           is not fulfilled)
!$$$$$$ !           info = 2 - means breakdown of the algorithm (try to enlarge
!$$$$$$ !           parameter l to get rid of this)
!$$$$$$
!$$$$$$       implicit none
!$$$$$$ !
!$$$$$$ !     .. Parameters ..
!$$$$$$ !
!$$$$$$ INTEGER(lo), INTENT(INOUT) :: l, n, mxmv, ldw, ldrw
!$$$$$$ INTEGER(lo), DIMENSION(:), INTENT(INOUT) :: iwork
!$$$$$$ INTEGER(lo), INTENT(OUT) :: info
!$$$$$$ COMPLEX(dbl), INTENT(IN) :: b(:)
!$$$$$$ COMPLEX(dbl), INTENT(INOUT) :: x(:)
!$$$$$$ REAL(dbl), INTENT(INOUT) :: tol
!$$$$$$ COMPLEX(dbl), DIMENSION(:,:), INTENT(INOUT) :: work(:,:), rwork(:,:)
!$$$$$$ !
!$$$$$$ !     .. Matrix ..
!$$$$$$ !
!$$$$$$ !
!$$$$$$ !     .. Local ..
!$$$$$$ !
!$$$$$$       logical rcmp, xpdt
!$$$$$$       integer i, j, k, nmv
!$$$$$$       complex*16 alpha, beta, omega, rho0, rho1, sigma
!$$$$$$       complex*16 varrho, hatgamma
!$$$$$$       real*8 rnrm0, rnrm
!$$$$$$       real*8 mxnrmx, mxnrmr
!$$$$$$       complex*16 kappa0, kappal
!$$$$$$ !
!$$$$$$ !     .. Work Aliases ..
!$$$$$$ !
!$$$$$$       integer z, zz, y0, yl, y
!$$$$$$       integer rr, r, u, xp, bp
!$$$$$$ !
!$$$$$$ !     .. Constants ..
!$$$$$$ !
!$$$$$$       real*8 delta
!$$$$$$       parameter (delta = 1d-2)
!$$$$$$       complex*16 zzero, zone
!$$$$$$       parameter (zzero = (0d0,0d0), zone = (1d0,0d0))
!$$$$$$ !
!$$$$$$ !     .. BLAS and LAPACK ..
!$$$$$$ !
!$$$$$$ !     subroutine zaxpy
!$$$$$$ !     subroutine zcopy
!$$$$$$ !     subroutine zgemv
!$$$$$$ !     subroutine zgetrf
!$$$$$$ !     subroutine zgetrs
!$$$$$$ !     subroutine zlacpy
!$$$$$$ !     subroutine zlaset
!$$$$$$ !     subroutine zsymv
!$$$$$$ !     subroutine zhemv
!$$$$$$ !     function zdotc
!$$$$$$ !     function dznrm2
!$$$$$$ !
!$$$$$$       real*8 dznrm2
!$$$$$$       complex*16 zdotc
!$$$$$$ !
!$$$$$$ !     .. Intrinsic ..
!$$$$$$ !
!$$$$$$       intrinsic abs, max, sqrt
!$$$$$$ !
!$$$$$$ !     ===========================
!$$$$$$ !     .. Executable Statements ..
!$$$$$$ !     ===========================
!$$$$$$ !
!$$$$$$       info = 0
!$$$$$$
!$$$$$$       if (l.lt.1) info = -1
!$$$$$$       if (n.lt.1) info = -2
!$$$$$$       if (tol.le.0d0) info = -7
!$$$$$$       if (mxmv.lt.0) info = -8
!$$$$$$
!$$$$$$       rr = 1
!$$$$$$       r = rr+1
!$$$$$$       u = r+(l+1)
!$$$$$$       xp = u+(l+1)
!$$$$$$       bp = xp+1
!$$$$$$       if (bp*n.gt.ldw) info = -10
!$$$$$$
!$$$$$$       z = 1
!$$$$$$       zz = z+(l+1)
!$$$$$$       y0 = zz+(l+1)
!$$$$$$       yl = y0+1
!$$$$$$       y = yl+1
!$$$$$$       if (y*(l+1).gt.ldrw) info = -12
!$$$$$$
!$$$$$$       if (info.ne.0) return
!$$$$$$ !
!$$$$$$ !     --- Initialize first residual
!$$$$$$ !
!$$$$$$       call matvec (n, x, work(1:n,r))
!$$$$$$       do i=1,n
!$$$$$$          work(i,r) = b(i) - work(i,r)
!$$$$$$       enddo
!$$$$$$       call precond (n, work(1:n,r))
!$$$$$$ !
!$$$$$$ !     --- Initialize iteration loop
!$$$$$$ !
!$$$$$$       nmv = 0
!$$$$$$
!$$$$$$       call zcopy (n, work(1:n,r), 1, work(1:n,rr), 1)
!$$$$$$       call zcopy (n, work(1:n,r), 1, work(1:n,bp), 1)
!$$$$$$       call zcopy (n, x, 1, work(1:n,xp), 1)
!$$$$$$       call zlaset ('n', n, 1, zzero, zzero, x, 1)
!$$$$$$       rnrm0 = dznrm2 (n, work(1:n,r), 1)
!$$$$$$       rnrm = rnrm0
!$$$$$$
!$$$$$$       mxnrmx = rnrm0
!$$$$$$       mxnrmr = rnrm0
!$$$$$$       rcmp = .false.
!$$$$$$       xpdt = .false.
!$$$$$$
!$$$$$$       alpha = zzero
!$$$$$$       omega = zone
!$$$$$$       sigma = zone
!$$$$$$       rho0 = zone
!$$$$$$ !
!$$$$$$ !     --- Iterate
!$$$$$$ !
!$$$$$$       do while (rnrm.gt.tol*rnrm0 .and. nmv.lt.mxmv)
!$$$$$$ !
!$$$$$$ !     =====================
!$$$$$$ !     --- The BiCG part ---
!$$$$$$ !     =====================
!$$$$$$ !
!$$$$$$          rho0 = -omega*rho0
!$$$$$$          do k=1,l
!$$$$$$             rho1 = zdotc (n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
!$$$$$$             if (rho0.eq.zzero) then
!$$$$$$                info = 2
!$$$$$$                return
!$$$$$$             endif
!$$$$$$             beta = alpha*(rho1/rho0)
!$$$$$$             rho0 = rho1
!$$$$$$             do j=0,k-1
!$$$$$$                do i=1,n
!$$$$$$                   work(i,u+j) = work(i,r+j) - beta*work(i,u+j)
!$$$$$$                enddo
!$$$$$$             enddo
!$$$$$$             call matvec (n, work(1:n,u+k-1), work(1:n,u+k))
!$$$$$$             call precond (n, work(1:n,u+k))
!$$$$$$             nmv = nmv+1
!$$$$$$             sigma = zdotc (n, work(1:n,rr), 1, work(1:n,u+k), 1)
!$$$$$$             if (sigma.eq.zzero) then
!$$$$$$                info = 2
!$$$$$$                return
!$$$$$$             endif
!$$$$$$             alpha = rho1/sigma
!$$$$$$             call zaxpy (n, alpha, work(1:n,u), 1, x, 1)
!$$$$$$             do j=0,k-1
!$$$$$$                call zaxpy (n, (-alpha), work(1:n,u+j+1), 1,&
!$$$$$$      &              work(1,r+j), 1)
!$$$$$$             enddo
!$$$$$$             call matvec (n, work(1:n,r+k-1), work(1:n,r+k))
!$$$$$$             call precond (n, work(1:n,r+k))
!$$$$$$             nmv = nmv+1
!$$$$$$             rnrm = dznrm2 (n, work(1:n,r), 1)
!$$$$$$             mxnrmx = max (mxnrmx, rnrm)
!$$$$$$             mxnrmr = max (mxnrmr, rnrm)
!$$$$$$          enddo
!$$$$$$ !
!$$$$$$ !     ==================================
!$$$$$$ !     --- The convex polynomial part ---
!$$$$$$ !     ==================================
!$$$$$$ !
!$$$$$$ !        --- Z = R'R
!$$$$$$ !
!$$$$$$          do i=1,l+1
!$$$$$$             call zgemv ('!', n, l+1-(i-1), zone, work(1:n,r+i-1),&
!$$$$$$      &           n, work(1:n,r+i-1), 1, zzero, rwork(i,z+i-1), 1)
!$$$$$$             if (l-(i-1).ne.0) then
!$$$$$$                call zcopy (l-(i-1), rwork(i+1,z+i-1), 1,&
!$$$$$$      &              rwork(i,z+i), l+1)
!$$$$$$                call zlacgv (l-(i-1), rwork(i,z+i), l+1)
!$$$$$$             endif
!$$$$$$          enddo
!$$$$$$          call zlacpy ('a', l+1, l+1, rwork(1:n,z), l+1,&
!$$$$$$      &        rwork(1,zz), l+1)
!$$$$$$          call zgetrf (l-1, l-1, rwork(2,zz+1), l+1,&
!$$$$$$      &        iwork, info)
!$$$$$$ !
!$$$$$$ !        --- tilde r0 and tilde rl (small vectors)
!$$$$$$ !
!$$$$$$          rwork(1,y0) = -zone
!$$$$$$          call zcopy (l-1, rwork(2,z), 1, rwork(2,y0), 1)
!$$$$$$          call zgetrs ('n', l-1, 1, rwork(2,zz+1), l+1, iwork,&
!$$$$$$      &        rwork(2,y0), l+1, info)
!$$$$$$          rwork(l+1,y0) = zzero
!$$$$$$
!$$$$$$          rwork(1,yl) = zzero
!$$$$$$          call zcopy (l-1, rwork(2,z+l), 1, rwork(2,yl), 1)
!$$$$$$          call zgetrs ('n', l-1, 1, rwork(2,zz+1), l+1, iwork,&
!$$$$$$      &        rwork(2,yl), l+1, info)
!$$$$$$          rwork(l+1,yl) = -zone
!$$$$$$ !
!$$$$$$ !        --- Convex combination
!$$$$$$ !
!$$$$$$          call zhemv ('u', l+1, zone, rwork(1:n,z), l+1,&
!$$$$$$      &        rwork(1,y0), 1, zzero, rwork(1:n,y), 1)
!$$$$$$          kappa0 = sqrt(zdotc (l+1, rwork(1:n,y0), 1,&
!$$$$$$      &        rwork(1,y), 1))
!$$$$$$
!$$$$$$          call zhemv ('u', l+1, zone, rwork(1:n,z), l+1,&
!$$$$$$      &        rwork(1,yl), 1, zzero, rwork(1:n,y), 1)
!$$$$$$          kappal = sqrt(zdotc (l+1, rwork(1:n,yl), 1,&
!$$$$$$      &        rwork(1,y), 1))
!$$$$$$
!$$$$$$          call zhemv ('u', l+1, zone, rwork(1:n,z), l+1,&
!$$$$$$      &        rwork(1,y0), 1, zzero, rwork(1:n,y), 1)
!$$$$$$          varrho =zdotc (l+1, rwork(1:n,yl), 1, rwork(1:n,y), 1)&
!$$$$$$      &        / (kappa0*kappal)
!$$$$$$
!$$$$$$          hatgamma = varrho/abs(varrho)*max(abs(varrho),7d-1)&
!$$$$$$      &        * (kappa0/kappal)
!$$$$$$
!$$$$$$          call zaxpy (l+1, (-hatgamma), rwork(1:n,yl), 1,&
!$$$$$$      &        rwork(1:n,y0), 1)
!$$$$$$ !
!$$$$$$ !        --- Update
!$$$$$$ !
!$$$$$$          omega = rwork(l+1,y0)
!$$$$$$
!$$$$$$          call zgemv ('n', n, l, (-zone), work(1:n,u+1), n,&
!$$$$$$      &        rwork(2,y0), 1, zone, work(1:n,u), 1)
!$$$$$$          call zgemv ('n', n, l, zone, work(1:n,r), n,&
!$$$$$$      &        rwork(2,y0), 1, zone, x, 1)
!$$$$$$          call zgemv ('n', n, l, (-zone), work(1:n,r+1), n,&
!$$$$$$      &        rwork(2,y0), 1, zone, work(1:n,r), 1)
!$$$$$$
!$$$$$$          call zhemv ('u', l+1, zone, rwork(1:n,z), l+1,&
!$$$$$$      &        rwork(1,y0), 1, zzero, rwork(1:n,y), 1)
!$$$$$$          rnrm = sqrt (zdotc (l+1, rwork(1:n,y0), 1,&
!$$$$$$      &        rwork(1:n,y), 1))
!$$$$$$ !
!$$$$$$ !     ================================
!$$$$$$ !     --- The reliable update part ---
!$$$$$$ !     ================================
!$$$$$$ !
!$$$$$$          mxnrmx = max (mxnrmx, rnrm)
!$$$$$$          mxnrmr = max (mxnrmr, rnrm)
!$$$$$$          xpdt = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.mxnrmx)
!$$$$$$          rcmp = ((rnrm.lt.delta*mxnrmr.and.rnrm0.lt.mxnrmr)&
!$$$$$$      &        .or.xpdt)
!$$$$$$          if (rcmp) then
!$$$$$$             call matvec (n, x, work(1:n,r))
!$$$$$$             call precond (n, work(1:n,r))
!$$$$$$             do i=1,n
!$$$$$$                work(i,r) =  work(i,bp) - work(i,r)
!$$$$$$             enddo
!$$$$$$             mxnrmr = rnrm
!$$$$$$             if (xpdt) then
!$$$$$$                call zaxpy (n, zone, x, 1, work(1:n,xp), 1)
!$$$$$$                call zlaset ('n', n, 1, zzero, zzero, x, 1)
!$$$$$$                call zcopy (n, work(1:n,r), 1, work(1:n,bp), 1)
!$$$$$$                mxnrmx = rnrm
!$$$$$$             endif
!$$$$$$          endif
!$$$$$$       enddo
!$$$$$$ !
!$$$$$$ !     =========================
!$$$$$$ !     --- End of iterations ---
!$$$$$$ !     =========================
!$$$$$$ !
!$$$$$$       call zaxpy (n, zone, work(1:n,xp), 1, x, 1)
!$$$$$$ !
!$$$$$$ !     --- Check stopping criterion
!$$$$$$ !
!$$$$$$       call matvec (n, x, work(1:n,r))
!$$$$$$       do i=1,n
!$$$$$$          work(i,r) = b(i) - work(i,r)
!$$$$$$       enddo
!$$$$$$       call precond (n, work(1:n,r))
!$$$$$$       rnrm = dznrm2 (n, work(1:n,r), 1)
!$$$$$$       if (rnrm.gt.tol*rnrm0) info = 1
!$$$$$$ !
!$$$$$$ !     --- Return
!$$$$$$ !
!$$$$$$       tol = rnrm/rnrm0
!$$$$$$       mxmv = nmv
!$$$$$$
!$$$$$$       return
!$$$$$$
!$$$$$$       end SUBROUTINE zbistbl



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
SUBROUTINE field_exp_shell(lambda,ref_index,v_req,m_epseq,v_p,v_dc_shell2,nstop,neq,matrixside,r_ih,v_qa,v_qb,&
				 & v_ab_shell2,v_dc_shell1,v_ab_shell1,v_dc_core,v_ab_host,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index,r_ih			! Lunghezza d'onda in questione, ref index,dist centri host core
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p			! Vettore espensioni campo incidente
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc_shell2		! Vettore d_mn2,c_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_qa,v_qb		! Matrici coeff. singola sfera
INTEGER(lo), INTENT(IN) :: nstop					! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside			! Numero sfere eq e lato matrice


COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_shell2			! Vettore a_mn2,b_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_shell1,v_dc_shell1	! Vettore a_mn1,b_mn1,d_mn1,c_mn1
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_host,v_dc_core		! Vettore a_mn,b_mn,d_mn,c_mn
INTEGER(lo), INTENT(OUT) :: error							! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq) :: v_x								! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq) :: v_epsc							! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq) :: v_k							! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq) :: v_m							! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq,1:neq) :: m_mx						! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z01,v_psi_z11,v_psi_z12,v_psi_z22	! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12			! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z12,v_psi1_z22				! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z12						! Vettori Csi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_qintd,v_qintc					! Vettore coefficienti ausiliari
INTEGER(lo) :: m,n,nu,i,j									! Indici
REAL(dbl) :: nr											! Indice complessificato



! Inizio della procedura vera e propria

!Costruisco il vettore complesso
v_epsc_do: DO j=1,neq
			!Metto neq-j,perche mentre nel formalismo numero da fuori a dentro, (0=host,1=shell,2=core), in m_epsqe numero
			!da dentro a fuori (1=core,2=shell). Cosi mi riporto al formalismo dell'articolo
			v_epsc(j)=CMPLX(m_epseq(neq+1-j,1),m_epseq(neq+1-j,2),KIND=dbl)
END DO v_epsc_do

!Calcolo il vettore per i size parameters
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
	v_psi1_z12(n)=v_psi_z12(n-1)-nr*v_psi_z12(n)/m_mx(1,2)
	v_psi1_z22(n)=v_psi_z22(n-1)-nr*v_psi_z22(n)/m_mx(2,2)

	!Derivate di csi
	v_csi1_z12(n)=v_csi_z12(n-1)-nr*v_csi_z12(n)/m_mx(1,2)

END DO der_do


!Calcolo il vettore ausiliario v_qint, per calcolare i coefficienti all'interfaccia interna
qint_do:DO n=1,nstop

	v_qintd(n)= ( v_k(2)*v_psi1_z12(n)*v_csi_z12(n) - v_k(2)*v_psi_z12(n)*v_csi1_z12(n) ) / &
		   & ( v_k(1)*v_psi1_z22(n)*v_csi_z12(n) - v_k(2)*v_psi_z22(n)*v_csi1_z12(n) )

	v_qintc(n)= ( v_k(2)*v_psi1_z12(n)*v_csi_z12(n) - v_k(2)*v_psi_z12(n)*v_csi1_z12(n) ) / &
		   & ( v_k(2)*v_psi1_z22(n)*v_csi_z12(n) - v_k(1)*v_psi_z22(n)*v_csi1_z12(n) )

END DO qint_do

!WRITE(*,*)
!WRITE(*,*) "qint",v_qint
!WRITE(*,*)

!-----------------------------
!Calcolo i Vettore a_mn2,b_mn2
!-----------------------------
j=0
amn2_n_do: DO n=1,nstop
	amn2_m_do: DO m=-n,n

		!Calcolo di amn2 e bmn2
		j=j+1
		v_ab_shell2(2*j-1)=-v_qa(n)*v_dc_shell2(2*j-1)
		v_ab_shell2(2*j)=-v_qb(n)*v_dc_shell2(2*j)

	END DO amn2_m_do
END DO amn2_n_do



!---------------------------------------------------
!Calcolo i vettori a_mn1 e b_mn1,d_mn1 e c_mn1
!ma uso l'if,perche se non traslo, allora no problem
!---------------------------------------------------
!WRITE(*,*) "lower bounds dc2",LBOUND(v_dc_shell2, 1)
!WRITE(*,*) "upper bounds dc2",UBOUND(v_dc_shell2, 1)
!WRITE(*,*) "lower bounds ab2",LBOUND(v_ab_shell2, 1)
!WRITE(*,*) "upper bounds ab2",UBOUND(v_ab_shell2, 1)
!WRITE(*,*) "lower bounds dc1",LBOUND(v_dc_shell1, 1)
!WRITE(*,*) "upper bounds dc1",UBOUND(v_dc_shell1, 1)
!WRITE(*,*) "lower bounds ab1",LBOUND(v_ab_shell1, 1)
!WRITE(*,*) "upper bounds ab1",UBOUND(v_ab_shell1, 1)

r_ih_if: IF (r_ih /= 0.0D0) THEN

	!Calcolo i vettori a_mn1 e b_mn1 tramite i vectore translation coefficients
	CALL matvec_trasl_shell(matrixside,v_ab_shell2,v_ab_shell1)

	!Calcolo i vettori d_mn1 e c_mn1 tramite i vectore translation coefficients
	CALL matvec_trasl_shell(matrixside,v_dc_shell2,v_dc_shell1)

ELSE


	v_ab_shell1=v_ab_shell2
	v_dc_shell1=v_dc_shell2

END IF r_ih_if




!---------------------------------------------------
!Calcolo i vettori d_mn e c_mn all'interfaccia interna
!---------------------------------------------------
j=0
dmn_n_do: DO n=1,nstop
	dmn_m_do: DO m=-n,n

		!Calcolo di amn2 e bmn2
		j=j+1
		v_dc_core(2*j-1)=v_qintd(n)*v_dc_shell2(2*j-1)
		v_dc_core(2*j)=v_qintc(n)*v_dc_shell2(2*j)

	END DO dmn_m_do
END DO dmn_n_do



!---------------------------------------------------
!Calcolo i vettori a_mn e b_mn all'interfaccia esterna
!---------------------------------------------------
j=0
amn_n_do: DO n=1,nstop
	amn_m_do: DO m=-n,n

		!Calcolo di amn2 e bmn2
		j=j+1
		v_ab_host(2*j-1)=( v_ab_shell1(2*j-1)*v_csi_z11(n) - v_dc_shell1(2*j-1)*v_psi_z11(n) + v_p(2*j-1)*v_psi_z01(n))&
				 & / (v_csi_z01(n))
		v_ab_host(2*j)=( v_ab_shell1(2*j)*v_k(0)*v_csi_z11(n) - v_dc_shell1(2*j)*v_k(0)*v_psi_z11(n) +  &
				     v_p(2*j)*v_k(1)*v_psi_z01(n)) / (v_k(1)*v_csi_z01(n))

	END DO amn_m_do
END DO amn_n_do



END SUBROUTINE field_exp_shell



!******************************************************************************
!******************************************************************************
!******************************************************************************
!2) SUBROUTINE coeff_shell: si calcolano i vettori dei coefficienti di per passare
! tra i diversi coefficienti della shell
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE field_exp_shell_dip(lambda,ref_index,v_req,m_epseq,v_p,v_dc_shell1,nstop,neq,matrixside,r_ih,v_qa,v_qb,tflag,&
				 & v_ab_shell1,v_dc_shell2,v_ab_shell2,v_dc_core,v_ab_host,error)

IMPLICIT NONE

! Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: lambda,ref_index,r_ih			! Lunghezza d'onda in questione, ref index,dist centri host core
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_req			! Vettore raggi equivalenti
REAL(dbl), DIMENSION(:,:), INTENT(IN) ::m_epseq			! Matrice funzioni dielettriche
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_p			! Vettore espensioni campo incidente
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_dc_shell1		! Vettore d_mn1,c_mn1
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_qa,v_qb		! Matrici coeff. singola sfera
INTEGER(lo), INTENT(IN) :: nstop					! Numero Espansioni multipolari
INTEGER(lo), INTENT(IN) :: neq,matrixside			! Numero sfere eq e lato matrice
INTEGER(lo), INTENT(IN) :: tflag					! Flag tipo di traslazione VTT


COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_shell1			! Vettore a_mn2,b_mn2
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_shell2,v_dc_shell2	! Vettore a_mn1,b_mn1,d_mn1,c_mn1
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_ab_host,v_dc_core		! Vettore a_mn,b_mn,d_mn,c_mn
INTEGER(lo), INTENT(OUT) :: error							! Controllo errore


! Dichiarazione variabili interne
REAL(dbl), DIMENSION(neq-1) :: v_x								! Vettori Size Parameters
COMPLEX(dbl), DIMENSION(neq-1) :: v_epsc							! Vettori complesso funct diel
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_k							! Vettori complesso wave vector
COMPLEX(dbl), DIMENSION(0:neq-1) :: v_m							! Vettore ref index normalizzato
COMPLEX(dbl),DIMENSION(0:neq-1,1:neq-1) :: m_mx						! Argomenti Psi e Csi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_psi_z11,v_psi_z12,v_psi_z22			! Vettori Psicomplessi
COMPLEX(dbl), DIMENSION(0:nstop) :: v_csi_z01,v_csi_z11,v_csi_z12,v_csi_z22	! Vettori Csi complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_psi1_z11						! Vettori Psi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_csi1_z01,v_csi1_z11				! Vettori Csi' complessi
COMPLEX(dbl), DIMENSION(1:nstop) :: v_qintd,v_qintc					! Vettore coefficienti ausiliari
INTEGER(lo) :: m,n,nu,i,j									! Indici
REAL(dbl) :: nr											! Indice complessificato



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

	!Derivate di csi
	v_csi1_z01(n)=v_csi_z01(n-1)-nr*v_csi_z01(n)/m_mx(0,1)
	v_csi1_z11(n)=v_csi_z11(n-1)-nr*v_csi_z11(n)/m_mx(1,1)

END DO der_do


!Calcolo il vettore ausiliario v_qint, per calcolare i coefficienti all'interfaccia interna
qint_do:DO n=1,nstop

	v_qintd(n)= ( v_k(0)*v_psi1_z11(n)*v_csi_z11(n) - v_k(0)*v_psi_z11(n)*v_csi1_z11(n) ) / &
		   & ( v_k(1)*v_csi1_z01(n)*v_csi_z11(n) - v_k(0)*v_csi_z01(n)*v_csi1_z11(n) )

	v_qintc(n)= ( v_k(0)*v_psi1_z11(n)*v_csi_z11(n) - v_k(0)*v_psi_z11(n)*v_csi1_z11(n) ) / &
		   & ( v_k(0)*v_csi1_z01(n)*v_csi_z11(n) - v_k(1)*v_csi_z01(n)*v_csi1_z11(n) )

END DO qint_do

!-----------------------------
!Calcolo i Vettore a_mn2,b_mn2
!-----------------------------
j=0
amn2_n_do: DO n=1,nstop
	amn2_m_do: DO m=-n,n

		!Calcolo di amn2 e bmn2
		j=j+1
		v_ab_shell1(2*j-1)=-v_qa(n)*v_dc_shell1(2*j-1)
		v_ab_shell1(2*j)=-v_qb(n)*v_dc_shell1(2*j)

	END DO amn2_m_do
END DO amn2_n_do

r_ih_if: IF (r_ih /= 0.0D0) THEN

	flag_trasl_if: IF (tflag==0) THEN	!Caso in cui uso gli stessi coefficienti secondo i Vector Translation Theorems

		!Calcolo i vettori a_mn1 e b_mn1 tramite i vectore translation coefficients
		CALL matvec_singleblock_axial(matrixside,v_sAij,v_sBij,v_iABij,v_jABij,v_ab_shell1,v_ab_shell2)

		!Calcolo i vettori d_mn1 e c_mn1 tramite i vectore translation coefficients
		CALL matvec_singleblock_axial(matrixside,v_sAij,v_sBij,v_iABij,v_jABij,v_dc_shell1,v_dc_shell2)

	ELSE						!Caso in cui uso i coeff diversi, secondo i Vector Translation Theorems

		!Calcolo i vettori a_mn1 e b_mn1 tramite i vectore translation coefficients
		CALL matvec_singleblock_axial(matrixside,v_sAijh,v_sBijh,v_iABij,v_jABij,v_ab_shell1,v_ab_shell2)

		!Calcolo i vettori d_mn1 e c_mn1 tramite i vectore translation coefficients
		CALL matvec_singleblock_axial(matrixside,v_sAij,v_sBij,v_iABij,v_jABij,v_dc_shell1,v_dc_shell2)

	END IF flag_trasl_if

ELSE


	v_ab_shell2=v_ab_shell1
	v_dc_shell2=v_dc_shell1

END IF r_ih_if




!---------------------------------------------------
!Calcolo i vettori d_mn e c_mn all'interfaccia interna
!---------------------------------------------------
j=0
dmn_n_do: DO n=1,nstop
	dmn_m_do: DO m=-n,n

		!Calcolo di amn2 e bmn2
		j=j+1
		v_ab_host(2*j-1)=-v_qintd(n)*v_dc_shell1(2*j-1)
		v_ab_host(2*j)=-v_qintc(n)*v_dc_shell1(2*j)

	END DO dmn_m_do
END DO dmn_n_do



!---------------------------------------------------
!Calcolo i vettori a_mn e b_mn all'interfaccia esterna
!---------------------------------------------------
j=0
amn_n_do: DO n=1,nstop
	amn_m_do: DO m=-n,n

		!Calcolo di amn2 e bmn2
		j=j+1
		v_dc_core(2*j-1)=-( v_ab_shell2(2*j-1)*v_csi_z12(n) - v_dc_shell2(2*j-1)*v_psi_z12(n) - v_p(2*j-1)*v_csi_z22(n))&
				 & / (v_psi_z22(n))
		v_dc_core(2*j)=-( v_ab_shell2(2*j)*v_k(2)*v_csi_z12(n) - v_dc_shell2(2*j)*v_k(2)*v_psi_z12(n) -  &
				     v_p(2*j)*v_k(1)*v_csi_z22(n)) / (v_k(1)*v_psi_z22(n))

	END DO amn_m_do
END DO amn_n_do



END SUBROUTINE field_exp_shell_dip


!******************************************************************************
!4quadris) SUBROUTINE matvec_rhs(matrixside,v_x,v_y): per il right hand side
!******************************************************************************
SUBROUTINE matvec_rhs_ss_dip(matrixside,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Funzione vera e propria

!WRITE(*,*) "QUI2",matrixside
!WRITE(*,*) "QUI3",v_y
!WRITE(*,*) "QUI4"

blockside=nstop*(nstop+2)
veclength=matrixside/2

! Alloco i vettori
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:veclength),v_bout(1:veclength))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
&        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Comincio la moltiplicazione
ns_do: DO i=1,ndip

!    WRITE(*,*) "ns:",i
!    WRITE(*,*)

    !Upper e lover bounds
    lowbi=1+blockside*(i-1)
    upbi=blockside*i

    !Faccio i fuori diagonale
    next_do: DO next=v_iBlock_rhs(i),v_iBlock_rhs(i+1)-1

        !Indice colonna
        j=v_jBlock_rhs(next)

!        IF (j<=(ns-ndip)) CYCLE next_do

        !Upper e lover bounds
        lowb=1+blockside*(j-1)
        upb=blockside*j

        !Tratto unicamente la porzione j-esima del vettore incognita
        v_a1=v_ain(lowb:upb)
        v_b1=v_bin(lowb:upb)

        !Qui ho la moltiplicazione fattorizzata

!        WRITE(*,*) "v_a1"
!        WRITE(*,*) (v_a1(n),n=lowb,upb)
!        WRITE(*,*)

        !Prima rotazione
        CALL rot1_sub(nstop,next,v_a1,v_a2)
        CALL rot1_sub(nstop,next,v_b1,v_b2)

!        WRITE(*,*) "v_a2"
!        WRITE(*,*) (v_a2(n),n=lowb,upb)
!        WRITE(*,*)

        !Traslazione
        CALL trasl_sub(nstop,next,v_a2,v_b2,v_a3)
        CALL trasl_sub(nstop,next,v_b2,v_a2,v_b3)

!        WRITE(*,*) "v_a3"
!        WRITE(*,*) (v_a3(n),n=lowb,upb)
!        WRITE(*,*)

        !Rotazione all'indietro
        CALL rot2_sub(nstop,next,v_a3,v_a4)
        CALL rot2_sub(nstop,next,v_b3,v_b4)

!        WRITE(*,*) "v_a4"
!        WRITE(*,*) (v_a4(n),n=lowb,upb)
!        WRITE(*,*)

        !Faccio la somma per l'update dei vettori
        v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
        v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

    END DO next_do
END DO ns_do

!WRITE(*,*)"v_aout"
!WRITE(*,*) (v_aout(n),n=1,matrixside/2)
!WRITE(*,*)

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_rhs_ss_dip









!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!SUBROUTINES DI MOLTIPLICAZIONE MATRICE VETTORE SPARSE DI CARATTERE GENERALE CON LE MATRICI CHE FIGURANO COME INPUT
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================
!===================================================================================================================================


!******************************************************************************
!1) SUBROUTINE matvec_singleblock_axial, per traslare i
!coefficienti das un sistema di riferimento all'altro nel caso di una
!traslazione puramente assiale
!******************************************************************************
SUBROUTINE matvec_singleblock_axial(matrixside,v_sAij,v_sBij,v_iABij,v_jABij,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: matrixside                               ! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_sAij,v_sBij               ! Contenuti matrici traslazioni
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iABij,v_jABij             ! Indici matrici traslazioni
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x                         ! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y                        ! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout

! Funzione vera e propria
blockside=nstop*(nstop+2)

! Alloco i vettori
ALLOCATE(v_ain(1:blockside),v_bin(1:blockside),v_aout(1:blockside),v_bout(1:blockside))

!Scrivo i miei due vettori in entrata
v_ain=v_x(1:matrixside:2)
!WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside:2)
!WRITE(*,*) "QUI3"
!Inizializzo i vettori
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)

!Comincio la moltiplicazione

!WRITE(*,*) "lower bounds ain",LBOUND(v_ain, 1)
!WRITE(*,*) "upper bounds ain",UBOUND(v_ain, 1)
!WRITE(*,*) "lower bounds bin",LBOUND(v_bin, 1)
!WRITE(*,*) "upper bounds bin",UBOUND(v_bin, 1)

!Traslazione
CALL trasl_singleblock(nstop,v_sAij,v_sBij,v_iABij,v_jABij,v_ain,v_bin,v_aout)
CALL trasl_singleblock(nstop,v_sAij,v_sBij,v_iABij,v_jABij,v_bin,v_ain,v_bout)

!Assegno finalmente i vettori
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_singleblock_axial


!******************************************************************************
!2) SUBROUTINE trasl_singleblock: faccio la traslazione generale,
!cioe' con le matrici che figurano come inputs, per un singolo blocco
!******************************************************************************
SUBROUTINE trasl_singleblock(nstop,v_sAij,v_sBij,v_iABij,v_jABij,v_a,v_b,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop					! N multipolar exp e posizione sparse
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_sAij,v_sBij			! Contenuti matrici traslazioni
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iABij,v_jABij		! Indici matrici traslazioni
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_a,v_b			! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y				! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_y1,v_y2
! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_y1(1:blockside),v_y2(1:blockside))


!Routine di moltiplicazione
!!$WRITE(*,*) 'qui1'
CALL amuzz(blockside,v_a,v_y1,v_sAij,v_jABij,v_iABij)
!!$WRITE(*,*) 'qui2'
CALL amuzz(blockside,v_b,v_y2,v_sBij,v_jABij,v_iABij)
!!$WRITE(*,*) 'qui3,y1,y2',v_y1,v_y2

!Calcolo il vettore
v_y=v_y1+v_y2
!!$WRITE(*,*) 'qui4'
DEALLOCATE(v_y1,v_y2)
!!$WRITE(*,*) 'qui5'
END SUBROUTINE trasl_singleblock


!******************************************************************************
!3) SUBROUTINE rot1_singleblock: faccio la prima rotazione generale,
!cioe' con le matrici che figurano come inputs, per un singolo blocco
!******************************************************************************
SUBROUTINE rot1_singleblock(nstop,v_sDij,v_iDij,v_jDij,v_exPhi,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! N multipolar exp
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_sDij				! Vettore coefficienti matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iDij,v_jDij		! Vettori indici matrice
COMPLEX(dbl), DIMENSION(-nstop:nstop), INTENT(IN) :: v_exPhi			! Vettore angoli phi
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x				! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y				! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_x1

! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_x1(1:blockside))


!D Matrix Parameters------------------------------------
!!$WRITE(*,*) 'inside the rotation'
!!$WRITE(*,*) "nstop",nstop
!!$WRITE(*,*) "v_exphi",v_exPhi
!!$WRITE(*,*) "v_iD",v_iDij
!!$WRITE(*,*) "v_jD",v_jDij
!!$WRITE(*,*) "v_D",v_sDij
!!$WRITE(*,*) "first time v_x",v_x



!Moltiplico per l'esponenziale
i=0
n_do: DO n=1,nstop
    m_do: DO m=-n,n

        i=i+1
        v_x1(i)=v_exPhi(m)*v_x(i)

    END DO m_do
END DO n_do

!!$WRITE(*,*) "second time v_x1",v_x1

!Routine di moltiplicazione
CALL amudz(blockside,v_x1,v_y,v_sDij,v_jDij,v_iDij)

!!$WRITE(*,*) "third time v_y",v_y

DEALLOCATE(v_x1)

END SUBROUTINE rot1_singleblock


!******************************************************************************
!2) SUBROUTINE rot2_singleblock: faccio la prima rotazione generale,
!cioe' con le matrici che figurano come inputs, per un singolo blocco
!******************************************************************************
SUBROUTINE rot2_singleblock(nstop,v_sDij,v_iDij,v_jDij,v_exPhi,v_x,v_y)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: nstop						! N multipolar exp
REAL(dbl), DIMENSION(:), INTENT(IN) :: v_sDij				! Vettore coefficienti matrice
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iDij,v_jDij		! Vettori indici matrice
COMPLEX(dbl), DIMENSION(-nstop:nstop), INTENT(IN) :: v_exPhi			! Vettore angoli phi
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x				! Vettore input
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y				! Vettore output

! Dichiarazione variabili interne
INTEGER(lo) :: i,n,m,blockside
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_x1

! Funzione vera e propria
blockside=nstop*(nstop+2)
ALLOCATE(v_x1(1:blockside))

!Cambio di segno
i=0
n_do: DO n=1,nstop
    m_do: DO m=-n,n

        i=i+1
        v_x(i)=v_oner(m)*v_x(i)

    END DO m_do
END DO n_do

!Routine di moltiplicazione
CALL amudz(blockside,v_x,v_y,v_sDij,v_jDij,v_iDij)

!Cambio di segno
i=0
n_do1: DO n=1,nstop
    m_do1: DO m=-n,n

        i=i+1
        v_y(i)=v_exPhi(-m)*v_oner(m)*v_y(i)

    END DO m_do1
END DO n_do1

END SUBROUTINE rot2_singleblock



!******************************************************************************
!4) SUBROUTINE matvec_gen(matrixside,v_x,v_y)
!******************************************************************************
SUBROUTINE matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
     v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
     v_x,v_y,error)

IMPLICIT NONE

! Dummy arguments
INTEGER(lo), INTENT(IN) :: nstop					! nstop
INTEGER(lo), INTENT(IN) :: blockside					! Size of one block
INTEGER(lo), INTENT(IN) :: imin,imax,jmin,jmax			! Matrix Boundaries
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iBlock,v_jBlock		! Row and column CSR Blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_Dij				! Dmn and phase CSR Blocks
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_exphi			! Dmn and phase CSR Blocks
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_iDij,m_jDij		! Dmn row and column CSR Blocks
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_Aij,m_Bij			! A,B CSR Blocks
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_iABij,m_jABij		! A,B row and column CSR Blocks
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_diag			! Coefficient matrix diagonal
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_x				! Input vector
COMPLEX(dbl), DIMENSION(:), INTENT(OUT) :: v_y				! Output vector
INTEGER(lo), INTENT(OUT) :: error					! Error flag

! Internal variables
INTEGER(lo) :: matrixside,matrixside_col,matrixside_vec,veclength,i,j,n,m,next,cursor,upb,lowb,upbi,lowbi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ain,v_bin,v_aout,v_bout
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_a1,v_b1,v_a2,v_b2,v_a3,v_b3,v_a4,v_b4

! Subroutine
error=0

!Checking matrix side and vector length
matrixside=2*(imax-imin+1)*blockside
matrixside_col=2*(jmax-jmin+1)*blockside
veclength=SIZE(v_x,1)

!!$WRITE(*,*) "matrixside,matrixside_col,veclength: ",matrixside,matrixside_col,veclength

length_if: IF (matrixside_col /= veclength) THEN
     WRITE(*,*)
     WRITE(*,10)
10   FORMAT ("Mismatch in matvec_gen between matrix and vector")
     error=1
     STOP
END IF length_if

veclength=veclength/2



! Allocating vectors
ALLOCATE(v_ain(1:veclength),v_bin(1:veclength),v_aout(1:matrixside/2),v_bout(1:matrixside/2))
ALLOCATE(v_a1(1:blockside),v_b1(1:blockside),v_a2(1:blockside),v_b2(1:blockside),&
     &        v_a3(1:blockside),v_b3(1:blockside),v_a4(1:blockside),v_b4(1:blockside))

!!$WRITE(*,*) 'matrixside,matrixside_col,veclength n1',matrixside,matrixside_col,veclength

!Assigning electric and magnetic input vectors
v_ain=v_x(1:matrixside_col:2)
!!$WRITE(*,*) "v_ain", v_ain
v_bin=v_x(2:matrixside_col:2)
!!$WRITE(*,*) "v_bin", v_bin

!!$WRITE(*,*) 'matrixside,matrixside_col,veclength n2',matrixside,matrixside_col,veclength

!Initializing intermediate vectors
v_aout=(0.0D0,0.0D0)
v_bout=(0.0D0,0.0D0)
v_a1=(0.0D0,0.0D0)
v_a2=(0.0D0,0.0D0)
v_a3=(0.0D0,0.0D0)
v_a4=(0.0D0,0.0D0)
v_b1=(0.0D0,0.0D0)
v_b2=(0.0D0,0.0D0)
v_b3=(0.0D0,0.0D0)
v_b4=(0.0D0,0.0D0)

!Beginning multiplications of the relevant subsection of the matrix
ns_do: DO i=imin,imax

     !Diagonal multiplication by identity matrix,to be done if and only if imin<=i<=imax and jmin<=i<=jmax
     diag_if: IF ( (imin<=i) .AND. (i<=imax) .AND. (jmin<=i) .AND. (i<=jmax) ) THEN

          !Upper e lower bounds, keeping in mind that vectors have an offset
          lowbi=1+blockside*(i-imin)
          upbi=blockside*(i-imin+1)

          !Taking the right portion of the vectors
          v_a1=v_ain(lowbi:upbi)
          v_b1=v_bin(lowbi:upbi)

          cursor=0
          diag_ndo: DO n=1,nstop
               diag_mdo: DO m=-n,n
                    cursor=cursor+1
                    v_a4(cursor)=v_a1(cursor)*v_diag(2*blockside*(i-1)+2*cursor-1)
                    v_b4(cursor)=v_b1(cursor)*v_diag(2*blockside*(i-1)+2*cursor)
               END DO diag_mdo
          END DO diag_ndo


          !Updating output
          v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
          v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

     END IF diag_if

     !Off diagonal blocks
     next_do: DO next=v_iBlock(i),v_iBlock(i+1)-1

          !Columns
          j=v_jBlock(next)

          !I cycle next_do if I am in the lower triangular part of the matrix
          IF (i>j) CYCLE next_do

          !-------------------------------------------------------------------------------------------
          !(i,j) Block multiplication, i.e. direct multiplication
          !-------------------------------------------------------------------------------------------

          !Skipping or exiting loop if I am before or after the selected columns
          IF (j<jmin) CYCLE next_do
          IF (j>jmax) EXIT next_do

          !-----Diagnostics over the cycle------------------
!!$          WRITE(*,*) "matvec dir,i,j",i,j
          !-------------------------------------------------

          !Upper e lower bounds, keeping in mind that vectors have an offset
          lowbi=1+blockside*(i-imin)
          upbi=blockside*(i-imin+1)
          lowb=1+blockside*(j-jmin)
          upb=blockside*(j-jmin+1)

          !getting the jth part of the
          v_a1=v_ain(lowb:upb)
          v_b1=v_bin(lowb:upb)

          !Matrix Vector Multiplication along the Mackowski Rotation Translation Rotation scheme

          !First Rotation
!!$          WRITE(*,*) "nstop",nstop
          CALL rot1_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_a1,v_a2)
          CALL rot1_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_b1,v_b2)

!!$          WRITE(*,*) "first rotation done"
!!$
!!$          WRITE(*,*) 'inside matvec'
!!$          WRITE(*,*) "nstop",nstop
!!$          WRITE(*,*) "m_exphi",m_exphi
!!$          WRITE(*,*) "m_iD",m_iDij
!!$          WRITE(*,*) "m_jD",m_jDij
!!$          WRITE(*,*) "m_D",m_Dij
!!$
!!$          WRITE(*,*) "v_a1", v_a1
!!$          WRITE(*,*) "v_a2", v_a2
!!$          WRITE(*,*) "v_b1", v_b1
!!$          WRITE(*,*) "v_b2", v_b2

          !Translation
          CALL trasl_singleblock(nstop,m_Aij(:,next),m_Bij(:,next),m_iABij(:,next),m_jABij(:,next),v_a2,v_b2,v_a3)
!!$          WRITE(*,*) "first translation done"
          CALL trasl_singleblock(nstop,m_Aij(:,next),m_Bij(:,next),m_iABij(:,next),m_jABij(:,next),v_b2,v_a2,v_b3)
!!$          WRITE(*,*) "second translation done"

          !Second Rotation
          CALL rot2_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_a3,v_a4)
          CALL rot2_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_b3,v_b4)

!!$          WRITE(*,*) "second rotation done"

          !Updating the vector subsection
          v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
          v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

          !-------------------------------------------------------------------------------------------
          !(j,i) Block multiplication, i.e. multiplication by simmetry rules
          !-------------------------------------------------------------------------------------------

          !Skipping or exiting loop if I am before or after the selected columns (now "i" addresses the column)
          IF (i<jmin) CYCLE next_do
          IF (i>jmax) EXIT next_do

          !-----Diagnostics over the cycle------------------
!!$          WRITE(*,*) "matvec sim,j,i",j,i
          !-------------------------------------------------

          !Upper e lower bounds, keeping in mind that vectors have an offset (now "i" and "j" switched)
          lowbi=1+blockside*(j-imin)
          upbi=blockside*(j-imin+1)
          lowb=1+blockside*(i-jmin)
          upb=blockside*(i-jmin+1)

          !getting the jth part of the
          v_a1=v_ain(lowb:upb)
          v_b1=v_bin(lowb:upb)

          !Matrix Vector Multiplication along the Mackowski Rotation Translation Rotation scheme

          !First Rotation
          CALL rot1_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_a1,v_a2)
          CALL rot1_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_b1,v_b2)

          !Translation
          CALL trasl_singleblock(nstop,v_mask*m_Aij(:,next),-v_mask*m_Bij(:,next),&
               m_iABij(:,next),m_jABij(:,next),v_a2,v_b2,v_a3)
          CALL trasl_singleblock(nstop,v_mask*m_Aij(:,next),-v_mask*m_Bij(:,next),&
               m_iABij(:,next),m_jABij(:,next),v_b2,v_a2,v_b3)

          !Second Rotation
          CALL rot2_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_a3,v_a4)
          CALL rot2_singleblock(nstop,m_Dij(:,next),m_iDij(:,next),m_jDij(:,next),m_exphi(:,next),v_b3,v_b4)

          !Updating the vector subsection
          v_aout(lowbi:upbi)=v_aout(lowbi:upbi)+v_a4
          v_bout(lowbi:upbi)=v_bout(lowbi:upbi)+v_b4

     END DO next_do
END DO ns_do

!Final update for the output vector
v_y(1:matrixside:2)=v_aout
v_y(2:matrixside:2)=v_bout

END SUBROUTINE matvec_gen



!******************************************************************************
!5) SUBROUTINE zbistbl(l,n,blockside,imin,imax,jmin,jmax,&
!		    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,&
!		    x,b,tol,mxmv,work,ldw,rwork,ldrw,iwork,info)
!******************************************************************************

SUBROUTINE zbistbl (l,n,blockside,imin,imax,jmin,jmax,&
		    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
		    x,b,tol,mxmv,work,ldw,rwork,ldrw,iwork,info)

! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.

! THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE AUTHOR
! PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESSED OR IMPLIED,
! REGARDING THE WORK, INCLUDING WARRANTIES WITH RESPECT TO ITS
! MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE.

IMPLICIT NONE

!Dummy arguments
INTEGER(lo), INTENT(INOUT) :: l,n,mxmv,ldw,ldrw,info
INTEGER(lo), INTENT(IN) :: blockside					! Size of one block
INTEGER(lo), INTENT(IN) :: imin,imax,jmin,jmax			! Matrix Boundaries
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iBlock,v_jBlock		! Row and column CSR Blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_Dij				! Dmn and phase CSR Blocks
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_exphi			! Dmn and phase CSR Blocks
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_iDij,m_jDij		! Dmn row and column CSR Blocks
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_Aij,m_Bij			! A,B CSR Blocks
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_iABij,m_jABij		! A,B row and column CSR Blocks
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_diag			! Coefficient matrix diagonal
INTEGER(lo), DIMENSION(:) ,INTENT(INOUT) :: iwork(:) 			! iwork(l+1)
COMPLEX(dbl), DIMENSION(:) ,INTENT(INOUT) :: x(:)
COMPLEX(dbl), DIMENSION(:) ,INTENT(IN):: b(:)
REAL(dbl), INTENT(INOUT) :: tol
COMPLEX(dbl) ,DIMENSION(:,:), INTENT(OUT):: work(:,:), rwork(:,:) 	! work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1))

!Internal variables
LOGICAL :: rcmp, xpdt
INTEGER(lo) :: i, j, k, nmv
COMPLEX(dbl) :: alpha, beta, omega, rho0, rho1, sigma
COMPLEX(dbl) :: varrho, hatgamma
REAL(dbl) :: rnrm0, rnrm
REAL(dbl) :: mxnrmx, mxnrmr
COMPLEX(dbl) :: kappa0, kappal

!Work Aliases
INTEGER(lo) :: z, zz, y0, yl, y,error
INTEGER(lo) :: rr, r, u, xp, bp,ii,jj

!Constants ..
REAL(dbl), PARAMETER :: delta=1.0D-2
COMPLEX(dbl), PARAMETER :: zzero=(0.0D0,0.0D0),zone=(1.0D0,0.0D0)

!BLAS and LAPACK
!subroutine zaxpy
!subroutine zcopy
!subroutine zgemv
!subroutine zgetrf
!subroutine zgetrs
!subroutine zlacpy
!subroutine zlaset
!subroutine zsymv
!subroutine zhemv
!function zdotc
!function dznrm2
REAL(dbl) :: dznrm2
COMPLEX(dbl) :: zdotc

!===========================
!.. Executable Statements ..
!===========================

!Checking errors
info = 0

IF (l<1) info = -1
IF (info/=0) RETURN

IF (n<1) info = -2
IF (info/=0) RETURN

IF (tol<=0.0D0) info = -7
IF (info/=0) RETURN

IF (mxmv<0) info = -8
IF (info/=0) RETURN

rr = 1
r = rr+1
u = r+(l+1)
xp = u+(l+1)
bp = xp+1
IF (bp*n>ldw) info = -10
IF (info/=0) RETURN

z = 1
zz = z+(l+1)
y0 = zz+(l+1)
yl = y0+1
y = yl+1
IF (y*(l+1)>ldrw) info = -12
IF (info/=0) RETURN


!Initialize first residual
CALL matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
		v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
		x,work(1:n,r),error) !mv(n, x,work(1:n,r))
DO i=1,n
	work(i,r) = b(i) - work(i,r)
END DO
CALL precond(n, work(1:n,r))

!   !~~~~~~~~~~~~~~~~~~~~
   WRITE(*,*) 'work,n2'
   Do,jj=1,9
   Do ii=1,n
   WRITE(*,400) work(ii,jj)
   400 FORMAT (2ES21.12)
   end do
   end do
   WRITE(*,*)
   STOP
   !~~~~~~~~~~~~~~~~~~~~

!Initialize iteration loop
nmv = 0

CALL zcopy (n, work(1:n,r), 1, work(1:n,rr), 1)
CALL zcopy (n, work(1:n,r), 1, work(1:n,bp), 1)
CALL zcopy (n, x, 1, work(1:n,xp), 1)
CALL zlaset ('n', n, 1, zzero, zzero, x, 1)
rnrm0 = dznrm2 (n, work(1:n,r), 1)
rnrm = rnrm0

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zzero
omega = zone
sigma = zone
rho0 = zone

!Iterate
do while ((rnrm.gt.tol*rnrm0) .and. (nmv.lt.mxmv))

	!=====================
	!--- The BiCG part ---
	!=====================
	rho0 = -omega*rho0
	do k=1,l

		rho1 = zdotc (n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
		if (rho0.eq.zzero) then
			info = 2
			return
		endif

		beta = alpha*(rho1/rho0)
		rho0 = rho1
		do j=0,k-1
			do i=1,n
				work(i,u+j) = work(i,r+j) - beta*work(i,u+j)
			enddo
		enddo

		call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
			    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
			    work(1:n,u+k-1),work(1:n,u+k),error)!mv(n,work(1:n,u+k-1),work(1:n,u+k))
		call precond (n, work(1:n,u+k))

		nmv = nmv+1
		sigma = zdotc (n, work(1:n,rr), 1, work(1:n,u+k), 1)
		if (sigma.eq.zzero) then
			info = 2
			return
		endif

		alpha = rho1/sigma
		call zaxpy (n, alpha, work(1:n,u), 1, x, 1)
		do j=0,k-1
			call zaxpy (n, (-alpha), work(1:n,u+j+1), 1, work(1:n,r+j), 1)
		enddo

		call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
			    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
			    work(1:n,r+k-1),work(1:n,r+k),error)!mv(n,work(1:n,r+k-1),work(1:n,r+k))
		call precond (n, work(1:n,r+k))
		nmv = nmv+1
		rnrm = dznrm2 (n, work(1:n,r), 1)
		WRITE(*,*) "n1,mv,res",nmv,rnrm
		mxnrmx = max (mxnrmx, rnrm)
		mxnrmr = max (mxnrmr, rnrm)

	enddo


	!==================================
	!--- The convex polynomial part ---
	!==================================

	!--- Z = R'R
	do i=1,l+1
		call zgemv ('c', n, l+1-(i-1), zone, work(1:n,r+i-1), n, work(1:n,r+i-1), 1, zzero, rwork(i,z+i-1), 1)
		if (l-(i-1).ne.0) then
			call zcopy (l-(i-1), rwork(i+1,z+i-1), 1, rwork(i,z+i), l+1)
			call zlacgv (l-(i-1), rwork(i,z+i), l+1)
		endif
	enddo
	call zlacpy ('a', l+1, l+1, rwork(1,z), l+1, rwork(1,zz), l+1)
	call zgetrf (l-1, l-1, rwork(2,zz+1), l+1, iwork, info)


	!        --- tilde r0 and tilde rl (small vectors)
	rwork(1,y0) = -zone
	call zcopy (l-1, rwork(2,z), 1, rwork(2,y0), 1)
	call zgetrs ('n', l-1, 1, rwork(2,zz+1), l+1, iwork, rwork(2,y0), l+1, info)
	rwork(l+1,y0) = zzero

	rwork(1,yl) = zzero
	call zcopy (l-1, rwork(2,z+l), 1, rwork(2,yl), 1)
	call zgetrs ('n', l-1, 1, rwork(2,zz+1), l+1, iwork, rwork(2,yl), l+1, info)
	rwork(l+1,yl) = -zone


	!        --- Convex combination
	call zhemv ('u', l+1, zone, rwork(1,z), l+1, rwork(1,y0), 1, zzero, rwork(1,y), 1)
	kappa0 = sqrt(zdotc (l+1, rwork(1,y0), 1, rwork(1,y), 1))

	call zhemv ('u', l+1, zone, rwork(1,z), l+1, rwork(1,yl), 1, zzero, rwork(1,y), 1)
	kappal = sqrt(zdotc (l+1, rwork(1,yl), 1, rwork(1,y), 1))

	call zhemv ('u', l+1, zone, rwork(1,z), l+1, rwork(1,y0), 1, zzero, rwork(1,y), 1)
	varrho = zdotc (l+1, rwork(1,yl), 1, rwork(1,y), 1) / (kappa0*kappal)

	hatgamma = varrho/abs(varrho)*max(abs(varrho),7d-1) * (kappa0/kappal)

	call zaxpy (l+1, (-hatgamma), rwork(1,yl), 1, rwork(1,y0), 1)


	!        --- Update
	omega = rwork(l+1,y0)

	call zgemv ('n', n, l, (-zone), work(1:n,u+1), n, rwork(2,y0), 1, zone, work(1:n,u), 1)
	call zgemv ('n', n, l, zone, work(1:n,r), n, rwork(2,y0), 1, zone, x, 1)
	call zgemv ('n', n, l, (-zone), work(1:n,r+1), n, rwork(2,y0), 1, zone, work(1:n,r), 1)

	call zhemv ('u', l+1, zone, rwork(1,z), l+1, rwork(1,y0), 1, zzero, rwork(1,y), 1)
	rnrm = sqrt (zdotc (l+1, rwork(1,y0), 1, rwork(1,y), 1))
	WRITE(*,*) "n2,mv,res",nmv,rnrm


	!================================
	!--- The reliable update part ---
	!================================
	mxnrmx = max (mxnrmx, rnrm)
	mxnrmr = max (mxnrmr, rnrm)
	xpdt = (rnrm.lt.delta*rnrm0.and.rnrm0.lt.mxnrmx)
	rcmp = ((rnrm.lt.delta*mxnrmr.and.rnrm0.lt.mxnrmr).or.xpdt)
	if (rcmp) then
		call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
			    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
			    x,work(1:n,r),error)!mv(n,x,work(1:n,r))
		call precond (n, work(1:n,r))
		do i=1,n
			work(i,r) =  work(i,bp) - work(i,r)
		enddo
		mxnrmr = rnrm
		if (xpdt) then
			call zaxpy (n, zone, x, 1, work(1:n,xp), 1)
			call zlaset ('n', n, 1, zzero, zzero, x, 1)
			call zcopy (n, work(1:n,r), 1, work(1:n,bp), 1)
			mxnrmx = rnrm
		endif
	endif

	WRITE(*,*) "n3,mv,res",nmv,rnrm

enddo

!=========================
!--- End of iterations ---
!=========================
call zaxpy (n, zone, work(1:n,xp), 1, x, 1)

if (rnrm.gt.tol*rnrm0) info = 1

!-------------------------------------------------------------------------------------------
!Computing the true residual: this is commented out to save one matvec
!-------------------------------------------------------------------------------------------
!call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
!		v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,&
!		x,work(1:n,r),error)!mv(n,x,work(1:n,r))
!do i=1,n
!	work(i,r) = b(i) - work(i,r)
!enddo
!call precond (n, work(1:n,r))
!rnrm = dznrm2 (n, work(1:n,r), 1)
!-------------------------------------------------------------------------------------------

!Return
tol = rnrm/rnrm0
mxmv = nmv
RETURN

END SUBROUTINE zbistbl



!******************************************************************************
!6) zbcg2(print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                  mxmatvec,work,info)
!******************************************************************************
subroutine zbcg2_gen(nstop,print_resid,nonzero_x,&
		     l,n,blockside,imin,imax,jmin,jmax,&
		     v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
		     x,rhs,toler,mxmatvec,work,info)
! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)
!
! Improved "vanilla" BiCGStab(2) iterative method
!
! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
!
! {{ This code based on original work of D.R.Fokkema:
!
! subroutine zbistbl v1.1 1998
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.
!
! }}
!
! Your bug reports, comments, etc. are welcome:
! m.a.botchev@math.utwente.nl
!
! ------------------------------
! Description of the parameters:
! ------------------------------
!
! print_resid (input) LOGICAL. If print_resid=.true. the number of
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration
!
! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices
!
! n          (input) INTEGER size of the linear system to solve
!
! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output
!
! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector
!
! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
!
! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not.
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided
!
! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of
!            toler is set to the actually achieved residual reduction
!
! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix
!            vector multiplications allowed to be done.  On output:
!            if info>=0, mxmatvec is set to the actual number of matrix
!            vector multiplications done
!
! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
!
! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)
!
! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
implicit none

! Parameters:
LOGICAL, INTENT(in)   :: print_resid,nonzero_x
INTEGER(lo), INTENT(in)   :: l, n,nstop
INTEGER(lo), INTENT(inout):: mxmatvec
INTEGER(lo), INTENT(out)  :: info
COMPLEX(dbl), DIMENSION(:) ,INTENT(inout):: x(:)
COMPLEX(dbl), DIMENSION(:) ,INTENT(in)   :: rhs(:)
REAL(dbl), INTENT(inout) :: toler
COMPLEX(dbl), DIMENSION(:,:), INTENT(out)  :: work(:,:)
INTEGER(lo), INTENT(IN) :: blockside					! Size of one block
INTEGER(lo), INTENT(IN) :: imin,imax,jmin,jmax			! Matrix Boundaries
INTEGER(lo), DIMENSION(:), INTENT(IN) :: v_iBlock,v_jBlock		! Row and column CSR Blocks
REAL(dbl), DIMENSION(:,:), INTENT(IN) :: m_Dij				! Dmn and phase CSR Blocks
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_exphi			! Dmn and phase CSR Blocks
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_iDij,m_jDij		! Dmn row and column CSR Blocks
COMPLEX(dbl), DIMENSION(:,:), INTENT(IN) :: m_Aij,m_Bij			! A,B CSR Blocks
INTEGER(lo), DIMENSION(:,:), INTENT(IN) :: m_iABij,m_jABij		! A,B row and column CSR Blocks
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_diag			! Coefficient matrix diagonal

! Local variables:
COMPLEX(dbl), DIMENSION(1:l+1) :: y0,yl,zy0,zyl
COMPLEX(dbl), DIMENSION(1:l+1,1:l+1)  :: matrix_z
LOGICAL    :: rcmp, xpdt
INTEGER(lo)    :: i, j, k, nmatvec,ii,jj
COMPLEX(dbl) :: alpha, beta, omega, rho0, rho1, sigma
COMPLEX(dbl) :: varrho, hatgamma
REAL(dbl)    :: rnrm0, rnrm
REAL(dbl)    :: mxnrmx, mxnrmr
COMPLEX(dbl) :: kappa0, kappal

! Aliases for the parts of the work array:
INTEGER(lo)          :: rr, r, u, xp, bp,error

! Constants:
REAL(dbl),    parameter :: delta = 1d-2
COMPLEX(dbl), parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

! Functions:
!REAL(dbl)     :: dnorm2_bcg
!COMPLEX(dbl) :: zdot_bcg

!!$WRITE(*,*) "nstop just inside solver",nstop

info = 0

if (l<1 .or. l>2) info = -2
if (n<1) info = -3
if (toler<=0d0) info = -9
if (mxmatvec<0) info = -10

rr = 1
r = rr+1
u = r+(l+1)
xp = u+(l+1)
bp = xp+1

if (info/=0) return

!~~~~~~~~~~~~~~~~~~~~
!!$WRITE(*,*) 'work,n1'
!!$Do,jj=1,9
!!$Do ii=1,n
!!$WRITE(*,300) work(ii,jj)
!!$300 FORMAT (2ES21.12)
!!$end do
!!$end do
!!$WRITE(*,*)
!~~~~~~~~~~~~~~~~~~~

! Initialize first residual
nmatvec=0
!!$WRITE(*,*) 'x',x
if (nonzero_x) then
!  if (print_resid) WRITE(*,*) "n1,nmv,res",nmatvec,dnorm2_bcg (n,x)
   call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
		   v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
		   x,work(1:n,r),error)

   !~~~~~~~~~~~~~~~~~~~~
!!$   WRITE(*,*) 'work,n2'
!!$   Do,jj=1,9
!!$   Do ii=1,n
!!$   WRITE(*,400) work(ii,jj)
!!$   400 FORMAT (2ES21.12)
!!$   end do
!!$   end do
!!$   WRITE(*,*)
!!$   if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
   !~~~~~~~~~~~~~~~~~~~~

   work(1:n,r) = rhs - work(1:n,r)
   nmatvec = 1

   !~~~~~~~~~~~~~~~~~~~~
!!$   WRITE(*,*) 'rhs'
!!$   Do ii=1,n
!!$   WRITE(*,450) rhs(ii)
!!$   450 FORMAT (2ES21.12)
!!$   end do
!!$   WRITE(*,*)
   !~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~
!!$   WRITE(*,*) 'work,n3'
!!$   Do,jj=1,9
!!$   Do ii=1,n
!!$   WRITE(*,500) work(ii,jj)
!!$   500 FORMAT (2ES21.12)
!!$   end do
!!$   end do
!!$   WRITE(*,*)
!!$   if (print_resid) print *,nmatvec,' ',dnorm2_bcg (n, work(1:n,r))
!   !~~~~~~~~~~~~~~~~~~~~

else
   work(1:n,r) = rhs
   nmatvec = 0
end if
CALL precond_gen(blockside,imin,imax,v_diag,work(1:n,r))

! Initialize iteration loop

work(1:n,rr) = work(1:n,r)
work(1:n,bp) = work(1:n,r)
work(1:n,xp) = x
x = zzero

rnrm0 = dnorm2_bcg (n, work(1:n,r))
rnrm = rnrm0

mxnrmx = rnrm0
mxnrmr = rnrm0
rcmp = .false.
xpdt = .false.

alpha = zzero
omega = zone
sigma = zone
rho0  = zone

! Iterate

!!$if (print_resid) WRITE(*,*) "tol",toler*rnrm0

do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)

! =====================
! The BiCG part ---
! =====================

   rho0 = -omega*rho0
   do k=1,l
      rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
      if (rho0.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      beta = alpha*(rho1/rho0)
      rho0 = rho1
      do j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
      enddo
      call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
			    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
			    work(1:n,u+k-1),work(1:n,u+k),error)
       CALL precond_gen(blockside,imin,imax,v_diag,work(1:n,u+k))
      nmatvec = nmatvec+1
!!$     if (print_resid) WRITE(*,*) "n2,nmv,res",nmatvec,dnorm2_bcg (n, work(1:n,r))

      sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
      if (sigma.eq.zzero) then
         info = 2
         toler = rnrm/rnrm0
         mxmatvec = nmatvec
         return
      endif
      alpha = rho1/sigma
      x(1:n) = alpha*work(1:n,u) + x(1:n)
      do j=0,k-1
         work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
      enddo
      call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
			    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
			    work(1:n,r+k-1),work(1:n,r+k),error)
       CALL precond_gen(blockside,imin,imax,v_diag,work(1:n,r+k))
      nmatvec = nmatvec+1
      rnrm = dnorm2_bcg (n, work(1:n,r))
!!$      if (print_resid) WRITE(*,*) "n3,nmv,res",nmatvec,rnrm
      mxnrmx = max (mxnrmx, rnrm)
      mxnrmr = max (mxnrmr, rnrm)
   enddo

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

   do i=1,l+1
      do j=1,i
         matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
      end do
   end do

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
   do j=2,l+1
      matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
   end do

!  small vectors y0 and yl

y0(1) = -zone
y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
y0(l+1) = zzero

yl(1) = zzero
yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
yl(l+1) = -zone


!  --- Convex combination

! compute Z*y0 and Z*yl
zy0 = zzero
zyl = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
   zyl = zyl + matrix_z(:,j)*yl(j)
end do

kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)

hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 )

y0 = y0 - (hatgamma*kappa0/kappal)*yl

!  --- Update

omega = y0(l+1)

do j=1,l
   work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
   x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
   work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
enddo

! y0 has changed; compute Z*y0 once more
zy0 = zzero
do j=1,l+1
   zy0 = zy0 + matrix_z(:,j)*y0(j)
end do

rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
!if (print_resid) WRITE(*,*) "n3.1,nmv,res",nmatvec,rnrm

! ================================
! The reliable update part ---
! ================================

mxnrmx = max (mxnrmx, rnrm)
mxnrmr = max (mxnrmr, rnrm)
xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
if (rcmp) then
   call matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
			    v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_Aij,m_Bij,m_iABij,m_jABij,v_diag,&
			    x,work(1:n,r),error)
    CALL precond_gen(blockside,imin,imax,v_diag,work(1:n,r))
   nmatvec = nmatvec + 1
   work(1:n,r) =  work(1:n,bp) - work(1:n,r)
   mxnrmr = rnrm
!  if (print_resid) WRITE(*,*) "n3,nmv,res",nmatvec,dnorm2_bcg (n, work(1:n,r))
   if (xpdt) then

      work(1:n,xp) = x(1:n) + work(1:n,xp)
      x = zzero
      work(1:n,bp) = work(1:n,r)

      mxnrmx = rnrm
   endif
endif

!!$if (print_resid) WRITE(*,*) "n4,nmv,res",nmatvec,rnrm

enddo

! =========================
! End of iterations ---
! =========================

x(1:n) = x(1:n) + work(1:n,xp)

if (rnrm>toler*rnrm0) info = 1

end subroutine zbcg2_gen




!******************************************************************************
!7) SUBROUTINE precond_gen(blockside,imin,imax,v_diag,v_x)
!******************************************************************************
SUBROUTINE precond_gen(blockside,imin,imax,v_diag,v_x)

IMPLICIT NONE

! Dichiarazione dummy arguments
INTEGER(lo), INTENT(IN) :: blockside,imin,imax		! Matrix blockside and block position
COMPLEX(dbl), DIMENSION(:), INTENT(IN) :: v_diag		! Matrix diagonal
COMPLEX(dbl), DIMENSION(:), INTENT(INOUT) :: v_x		! Unknowns vector



! Subroutine
v_x=v_x/v_diag(1+2*blockside*(imin-1):2*blockside*imax)

END SUBROUTINE precond_gen


END MODULE linear_solver
