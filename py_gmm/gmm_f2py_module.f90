!--------------------------------------------------------------
!Precisions ect...
!--------------------------------------------------------------
MODULE local_kinds

IMPLICIT NONE

SAVE

INTEGER, PARAMETER :: sgl = SELECTED_REAL_KIND(p=6,r=37)
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=15,r=307)
!INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(32)
INTEGER, PARAMETER :: vshort = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: short = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: long = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: lo = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: length=30

REAL(sgl), PARAMETER :: PI=3.1415926535897932384626433832795028841971
REAL(sgl), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858
REAL(sgl), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394
REAL(sgl), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967
REAL(sgl), PARAMETER :: EULER=0.5772156649015328606065120900824024310422

REAL(dbl), PARAMETER :: PI_D=3.141592653589793238462643383279502884197
REAL(dbl), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858
REAL(dbl), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394
REAL(dbl), PARAMETER :: SQRT2_D=1.41421356237309504880168872420969807856967
REAL(dbl), PARAMETER :: EULER_D=0.5772156649015328606065120900824024310422
REAL(dbl), PARAMETER :: e=1.6021765314e-19										! electron charge
REAL(dbl), PARAMETER :: me=9.109382616e-31										! electron mass
REAL(dbl), PARAMETER :: e0=8.854187817e-12									  ! vacuum dielectri constant
REAL(dbl), PARAMETER :: c=2.99792458e8											  ! light speed in vacuum
REAL(dbl), PARAMETER :: zerodbl=0.0D0   											! double precision zero

END MODULE

!--------------------------------------------------------------
!Actual subroutines
!--------------------------------------------------------------
MODULE gmm_f2py_module

USE datatypes			! derived data types
USE operators			! operator overloading
USE basicsubs			! "low level" subroutines
USE gmmsubs				! "high level" subroutines
USE sing_part			! single particle subroutines
USE vec_trans			! vector translation coefficients subroutines
USE local_field			! local field computation subroutines
USE shared_data			! shared data subroutines
USE linear_solver		! linear solver subroutines

CONTAINS


!--------------------------------------------------------------------------------
!1) SUBROUTINE emn_sub: field normalization factors
!--------------------------------------------------------------------------------
SUBROUTINE emn(nstop,v_emn,error)
USE local_kinds				! precisions

IMPLICIT NONE

! dummy argument
INTEGER(lo), INTENT(IN) :: nstop																! multipolar expansions
COMPLEX(dbl), DIMENSION(nstop*(nstop+2)), INTENT(OUT) :: v_emn	! normalization coefficients
INTEGER(lo) , INTENT(OUT) :: error														  ! error flags

! internal variables
INTEGER(lo) :: m,n,i											! indexes and dimensions
REAL(dbl) :: nr,mr,logw										! indexes real kind

! subroutine begins here
error=0
i=0
n_do: DO n=1,nstop

	nr=REAL(n,dbl)

	m_do: DO m=-n,n

		i=i+1
		mr=REAL(m,dbl)

		logw=lnf(nr-mr,error)-lnf(nr+mr,error)
		v_emn(i)=((0.0D0,1.0D0)**n)*SQRT(EXP(logw)*(2*nr+1.0D0)/(nr*(nr+1.0D0)))
!		WRITE(*,*) m,n,v_emn(i)

	END DO m_do
END DO n_do

END SUBROUTINE emn

!--------------------------------------------------------------------------------
!2) SUBROUTINE shell_coefficients: shell multipolar coefficients
!--------------------------------------------------------------------------------
SUBROUTINE expansion_coefficients(ns,m_xyz_inp,v_r_inp,m_eps_inp,fint,&
                                 &ref_index,lambda_inp,alpha0,beta0,gamma0,nstop,&
								 							   &quasi_static_flag,&
                                 &m_coeff,v_cext,v_csca,v_cabs)

  USE local_kinds				! Kinds numerici, come single e double precision


  IMPLICIT NONE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Dummy arguments declaration
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL(dbl), DIMENSION(ns,3), INTENT(IN) :: m_xyz_inp			! Sphere coordinates
  REAL(dbl), DIMENSION(ns), INTENT(IN) :: v_r_inp				! Sphere radii
  REAL(dbl), DIMENSION(ns,2), INTENT(IN) :: m_eps_inp			! Sphere e1,e2
  REAL(dbl), INTENT(IN) :: fint                                 ! Interaction parameter
  REAL(dbl), INTENT(IN) :: ref_index							! Matrix Refractive Index
  REAL(dbl), INTENT(INOUT) :: lambda_inp							! Vacuum Wavelength
  REAL(dbl), INTENT(INOUT) :: alpha0,beta0,gamma0				! Euler angles for incident light
  INTEGER(lo), INTENT(IN) :: ns,nstop							! Spheres and Multipolar expansions
  CHARACTER(len=3), INTENT(IN) :: quasi_static_flag				!Flag to force quasi-static incident field
  COMPLEX(dbl), DIMENSION(2*nstop*(nstop+2)*ns,2), INTENT(OUT):: m_coeff! Expansion coefficients
  REAL(dbl), DIMENSION(ns), INTENT(OUT) :: v_cext,v_csca,v_cabs         ! Cross Sections

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Internal Variables
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! File opening
  INTEGER(lo) :: status_in, status_arr		! Status apertura file principale e allocazione array
  INTEGER(lo) :: neq,matrixside,i,j,m		! Numero sfere totali, numero materiali e numero gruppi di sfere uguali
  INTEGER(lo) :: nnb,NNZab,NNZd,blockside	! Non zero blocks, elements Aij Bij, elements Dij, block side
  INTEGER(lo) :: error			        ! Flag di errore

  REAL(dbl) :: t1,t2				! Tempi
  REAL(dbl) :: lambda0,lambda1,lambda2,k,m1,m2	! Lunghezza d'onda
  REAL(dbl) :: lambdamin,lambdamax,lambdastep	! Lunghezza d'onda

  ! Single particle variables
  COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_da,m_cb! Matrici coefficienti calcolo campo interno

  ! General internal variables
  INTEGER(lo) :: sumqsp,sumqqsp				!Somme di qmax e Qmax nel caso sparse
  INTEGER(lo), ALLOCATABLE, DIMENSION(:,:) :: m_index	! Index matrix
  REAL(lo), ALLOCATABLE, DIMENSION(:,:) :: m_Apmn	! Apmn coefficient matrix
  REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_c0,v_aq,v_bq!Vettori coefficienti c0,gaunt e bq per Aij Bij

  COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_p,v_p1,v_ab,v_ab1,v_ab_sca,v_dc,v_diag! Vettori espansioni campo e vettore noto
  INTEGER(lo) :: info,n								! Flag errore

  CHARACTER(len=3) :: norm='new'						! Flag normalizz fill_AB
  CHARACTER(len=3) :: norm1='div'						! Flag normalizz matrice

  ! Linear solver variables
  INTEGER(lo) :: lvanilla=2,mxmv,mxmv0,cont		! Bicgstab dimension and max iterations
  REAL(dbl) :: tol,tol0					! Bicgstab tolerance
  COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: work	! Auxiliary array for bicgstab
  LOGICAL :: logic1,logic2

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Dirty hacks to interface fortran and python properly
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  neq=ns
  ALLOCATE (v_patt(1:ns),STAT=status_arr)
  DO i=1,ns
     v_patt(i)=i
  END DO



!!$  WRITE(*,*) 'ns',ns
!!$  WRITE(*,*) 'm_xyz_inp',m_xyz_inp
!!$  WRITE(*,*) 'v_r_inp',v_r_inp
!!$  WRITE(*,*) 'm_eps_inp',m_eps_inp
!!$  WRITE(*,*) 'fint',fint
!!$  WRITE(*,*) 'ref_index',ref_index
!!$  WRITE(*,*) 'lambda_inp',lambda_inp
!!$  WRITE(*,*) 'alpha0,beta0,gamma0',alpha0,beta0,gamma0
!!$  WRITE(*,*) 'nstop',nstop
!!$  WRITE(*,*) 'v_patt',v_patt

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Checking sphere intersection
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL intspheres(v_r_inp,m_xyz_inp,ns,error)

  err_intspheres: IF (error/=0) THEN
     WRITE(*,*)
     WRITE(*,12)
12   FORMAT ("Intersection spheres, the program stops.")
     STOP
  END IF err_intspheres


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Non zero blocks in coeff matrix
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL nnb_sparse(ns,m_xyz_inp,v_r_inp,fint,nnb) !Divided by two because i store only the upper triangular part
  nnb=nnb/2
!!$  WRITE(*,*) "Non zero blocks:", nnb

  !Initializing v_ic to eliminate the load on the intrinsic library pow
  ic_do: DO i=0,400
     v_ic(i)=(0.0D0,1.0D0)**i
  END DO ic_do

  !Initializing v_ir to eliminate the load on the intrinsic library pow
  ir_do: DO i=-200,200
     v_oner(i)=(-1.0D0)**i
  END DO ir_do


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! If I have blocks, I start filling them now
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  zero_Blockif1: IF (nnb/=0) THEN

     !nnb multiplied by 2,because i store the whole skeleton of the matrix for flexibility
     ALLOCATE (v_jBlock(1:nnb),v_iBlock(1:ns+1),STAT=status_arr)
     CALL fill_jBlock_iBlock_sparse(ns,m_xyz_inp,v_r_inp,fint,v_jBlock,v_iBlock)

     !--------------Filling diagnostics-----------------------------
!!$     WRITE(*,*) 'v_iBlock',v_iBlock
!!$     WRITE(*,*) 'v_jBlock',v_jBlock
     !--------------------------------------------------------------

     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !Elements for Aij e Dij, matrix sizes and so on
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !Non zero elements Aij Bij
     NNZab=nstop*(1+2*nstop*(3+nstop))
     NNZab=NNZab/3
     !Non zero elements Dij
     NNZd=nstop*(11+12*nstop+4*(nstop**2))
     NNZd=NNZd/3
     !Lato di un blocco
     blockside=nstop*(nstop+2)
     matrixside=2*ns*blockside

     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! Vector allocation
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ALLOCATE (m_sAij(1:NNZab,1:nnb),m_sBij(1:NNZab,1:nnb),m_sAij_sca(1:NNZab,1:nnb),&
          m_sBij_sca(1:NNZab,1:nnb),m_jABij(1:NNZab,1:nnb),m_iABij(1:blockside+1,1:nnb),STAT=status_arr)
     ALLOCATE (m_Dij(1:NNZd,1:nnb),m_jDij(1:NNZd,1:nnb),m_iDij(1:blockside+1,1:nnb),&
          &m_exphi(-nstop:nstop,1:nnb),STAT=status_arr)
     ALLOCATE(work(1:matrixside,1:3+2*(lvanilla+1)),STAT=status_arr)
     !Filling the rotational part
     CALL fill_D_PHI_sparse(ns,nstop,nnb,m_xyz_inp,v_jBlock,v_iBlock,m_Dij,m_jDij,m_iDij,m_exphi)

!!$     WRITE(*,*) 'm_Dij',m_Dij
!!$     WRITE(*,*) 'm_jDij',m_jDij
!!$     WRITE(*,*) 'm_iDij',m_iDij
!!$     WRITE(*,*) 'm_exphi',m_exphi

     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !Aij Bij wavelength independent coefficients
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !Coefficienti c0
     ALLOCATE(v_c0(1:NNZab))
     CALL c0_sparse(nstop,v_c0)

     !Coefficienti aq e bq
     CALL qsum_sparse(nstop,sumqsp,sumqqsp)
     ALLOCATE(v_aq(1:sumqsp))
     ALLOCATE(v_bq(1:sumqqsp))
     CALL fill_aqbq_sparse(nstop,v_aq,v_bq)

     !Calculation support coefficients for ABij once for all
     ALLOCATE(m_index(1:NNZab,1:11),STAT=status_arr)
     ALLOCATE(m_Apmn(1:NNZab,1:(2*nstop)+1),STAT=status_arr)
     ALLOCATE(v_jABij_template(1:NNZab),v_iABij_template(1:blockside+1),v_mask(1:NNZab),STAT=status_arr)
     CALL fill_index(nstop,m_index,m_Apmn,v_jABij_template,v_iABij_template,v_mask,error)
     !		Write(*,600) m_index(:,4)
600  FORMAT (I5)
  END IF zero_Blockif1

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Allocating what is left, i.e. expansion
  !coefficients, cross sections, etc...
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Lato di un blocco
  blockside=nstop*(nstop+2)
  matrixside=2*ns*blockside
  imin=1
  imax=ns
  jmin=1
  jmax=ns

  ALLOCATE (v_p(matrixside),v_ab(matrixside),v_diag(matrixside),&
       v_dc(matrixside),v_ab_sca(matrixside),v_ab1(matrixside),STAT=status_arr)

  !quasi static approximation
  qs_if: IF (quasi_static_flag=='yes') THEN
	lambda_inp=MAXVAL(v_r_inp)*1.0D2
  END IF qs_if

  !Wavevector
  k=2*Pi_d*ref_index/(lambda_inp)



  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Single particle coefficients
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ALLOCATE(m_a(1:nstop,1:neq),m_b(1:nstop,1:neq),m_da(1:nstop,1:neq),m_cb(1:nstop,1:neq),STAT=status_arr)

  CALL coeff_sp2(lambda_inp,ref_index,v_r_inp,m_eps_inp,nstop,neq,m_a,m_b,error)

  coeff_sp_if: IF (error/=0) THEN
     WRITE(*,70)
70   FORMAT ("Error in coeff_sp2, stopping the program...")
     STOP
  END IF coeff_sp_if

  CALL coeff_ad_ca(lambda_inp,ref_index,v_r_inp,m_eps_inp,nstop,neq,m_da,m_cb,error)

  coeff_dacb_if: IF (error/=0) THEN
     WRITE(*,71)
71   FORMAT ("Error in coeff_ad_ca, stopping the program a...")
     STOP
  END IF coeff_dacb_if


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Filling the coefficient matrix
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  zero_Blockif2: IF (nnb/=0) THEN
     CALL fill_AB_sparse_sca_fast(ns,nstop,k,m_xyz_inp,v_c0,v_aq,v_bq,m_index,m_Apmn,&
          &v_jABij_template,v_iABij_template,v_jBlock,v_iBlock,m_sAij,m_sBij,&
          &m_sAij_sca,m_sBij_sca,m_jABij,m_iABij)
  END IF zero_Blockif2

  !Filling the diagonal of the matrix
  j=0
  DO i=1,ns
     DO n=1,nstop
        DO m=-n,n
           j=j+1
           v_diag(2*j-1)=1.0D0/m_a(n,v_patt(i))
           v_diag(2*j)=1.0D0/m_b(n,v_patt(i))
        END DO
     END DO
  END DO

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Filling the RHS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Incident field expansion coefficients
  CALL field_expRandom_sub(nstop,ns,k,zerodbl,m_xyz_inp,alpha0,beta0,gamma0,v_p,error)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Solving the linear system
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  zero_Blockif3: IF (nnb/=0) THEN

     !Initial guess
     v_ab=v_p/v_diag

     !BicgSTAB linear solver: initializing parameters
     mxmv0=5000
     tol0=1.0D-15
     work=(1.0d0,1.0d0)
     logic1=.TRUE.
     logic2=.TRUE.
     mxmv=mxmv0
     tol=tol0
     info=0
!!$     WRITE(*,*) "l,lambda,mxmv,tol,matrixside,info: " ,lvanilla,lambda_inp,mxmv,tol,matrixside,info

!!$     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$     ! INPUT DIAGNOSTICS
!!$     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$     !		!Logic Variables----------------------------------------
!!$     WRITE(*,*) "logic1",logic1
!!$     WRITE(*,*) "logic2",logic2
!!$     !Integer Parameters-------------------------------------
!!$     WRITE(*,*) "nstop",nstop
!!$     WRITE(*,*) "lvanilla",lvanilla
!!$     WRITE(*,*) "matrixside",matrixside
!!$     WRITE(*,*) "blockside",blockside
!!$     WRITE(*,*) "imin,imax,jmin,jmax,mxmv,info"
!!$     WRITE(*,200) imin,imax,jmin,jmax,mxmv,info
!!$200  FORMAT (6I10)
!!$     !Block Matrix Parameters--------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_iBlock"
!!$     DO i=1,ns+1
!!$        WRITE(*,205) v_iBlock(i)
!!$205     FORMAT (I)
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_jBlock"
!!$     DO i=1,nnb
!!$        WRITE(*,210) v_jBlock(i)
!!$210     FORMAT (I)
!!$     END DO
!!$
!!$     !D Matrix Parameters------------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_iD"
!!$     DO i=1,nnb
!!$        DO j=1,ns+1
!!$           WRITE(*,215) m_iDij(j,i)
!!$215        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_jD"
!!$     DO i=1,nnb
!!$        DO j=1,NNZd
!!$           WRITE(*,220) m_jDij(j,i)
!!$220        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_D"
!!$     DO i=1,nnb
!!$        DO j=1,NNZd
!!$           WRITE(*,225) m_Dij(j,i)
!!$225        FORMAT (1ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     !AB Matrix Parameters-----------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_iAB"
!!$     DO i=1,nnb
!!$        DO j=1,ns+1
!!$           WRITE(*,230) m_iABij(j,i)
!!$230        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_jAB"
!!$     DO i=1,nnb
!!$        DO j=1,NNZab
!!$           WRITE(*,235) m_jABij(j,i)
!!$235        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_A"
!!$     DO i=1,nnb
!!$        DO j=1,NNZab
!!$           WRITE(*,240) m_sAij(j,i)
!!$240        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_B"
!!$     DO i=1,nnb
!!$        DO j=1,NNZab
!!$           WRITE(*,245) m_sBij(j,i)
!!$245        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     !Phi Matrix Parameters----------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_phi"
!!$     DO i=1,nnb
!!$        DO j=-nstop,nstop
!!$           WRITE(*,250) m_exphi(j,i)
!!$250        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     !Matrix Diagonal---------------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_diag"
!!$     DO i=1,matrixside
!!$        WRITE(*,255) v_diag(i)
!!$255     FORMAT (2ES16.7)
!!$     END DO
!!$
!!$     !Unknown vector-----------------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_ab"
!!$     DO i=1,matrixside
!!$        WRITE(*,260) v_ab(i)
!!$260     FORMAT (2ES16.7)
!!$     END DO
!!$
!!$     !RHS----------------------------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_p"
!!$     DO i=1,matrixside
!!$        WRITE(*,265) v_p(i)
!!$265     FORMAT (2ES16.7)
!!$     END DO
!!$
!!$     !workspace----------------------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "work"
!!$     DO j=1,3+2*(lvanilla+1)
!!$        DO i=1,matrixside
!!$           WRITE(*,270) work(i,j)
!!$270        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     WRITE(*,*)
!!$     WRITE(*,*) "tolerance"
!!$     WRITE(*,275) tol
!!$275  FORMAT (1ES16.7)
!!$
!!$     !BicgSTAB linear solver: calling the solver
!!$     WRITE(*,*) "nstop just outside the solver",nstop
     CALL zbcg2_gen(nstop,logic1,logic2,&
          lvanilla,matrixside,blockside,imin,imax,jmin,jmax,&
          v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_sAij,m_sBij,m_iABij,m_jABij,v_diag,&
          v_ab,v_p,tol,mxmv,work,info)


!!$     WRITE(*,*) 'qui1'


     info_if: IF (info /= 0) THEN
        WRITE(*,*) "Problem in bicgstab", info
        STOP
     END IF info_if


  ELSE	!No interaction among the nanoparticles
     v_ab=v_p/v_diag
  END IF zero_Blockif3

!!$  WRITE(*,*) 'qui2'

  !Computing v_ab_sca
  v_diag=1.0D0
  CALL matvec_gen(nstop,blockside,imin,imax,jmin,jmax,&
       v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_sAij_sca,m_sBij_sca,m_iABij,m_jABij,v_diag,&
       v_ab,v_ab_sca,error)

!!$  WRITE(*,*) 'qui3'

  !Internal field coefficients
  CALL dmncmn_sub(nstop,ns,v_patt,m_da,m_cb,v_ab,v_dc)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Filling m_coeff with all the coefficients
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !ab
  m_coeff(:,1)=v_ab
  !dc
  m_coeff(:,2)=v_dc

!!$  WRITE(*,*) 'qui4'

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Cross sections
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL cext_random_sub(nstop,ns,k,v_p,v_ab,v_cext)
  CALL csca_random_sub(nstop,ns,k,v_ab_sca,v_ab,v_csca)
  CALL cabs_random_sub(ref_index,nstop,ns,k,v_dc,v_r_inp,m_eps_inp,v_cabs)

!!$  WRITE(*,*) 'qui5'

  DEALLOCATE(v_patt)
  DEALLOCATE (v_jBlock,v_iBlock)
  DEALLOCATE (m_sAij,m_sBij,m_sAij_sca,m_sBij_sca,m_jABij,m_iABij)
  DEALLOCATE (m_Dij,m_jDij,m_iDij,m_exphi,STAT=status_arr)
  DEALLOCATE(work)
  DEALLOCATE(v_c0)
  DEALLOCATE(v_aq)
  DEALLOCATE(v_bq)
  DEALLOCATE(m_index)
  DEALLOCATE(m_Apmn)
  DEALLOCATE(v_jABij_template,v_iABij_template,v_mask)
  DEALLOCATE (v_p,v_ab,v_diag,v_dc,v_ab_sca,v_ab1)
  DEALLOCATE(m_a,m_b,m_da,m_cb)

  RETURN

  END SUBROUTINE expansion_coefficients




  !-------------------------------------------------------------------------------
  !3) Exyz_sub: calculating field in the (x,y,z) coordinates
  !-------------------------------------------------------------------------------
  SUBROUTINE Exyz(flaginc,ns,nstop,ratio,lambda_inp,alpha0,beta0,gamma0,x,y,z,v_amnbmn,&
       &v_dmncmn,v_emn,m_xyz_inp,m_eps_inp,v_r_inp,ref_index,&
	   &quasi_static_flag,&
	   &Ex,Ey,Ez,Eabs,error)

    USE local_kinds

    IMPLICIT NONE

    !Dichiarazione dei dummy argument
    CHARACTER(len=3), INTENT(IN) :: flaginc			!Flag per aggiungere il campo incidente
    REAL(dbl), DIMENSION(ns,3), INTENT(IN) :: m_xyz_inp		!Posizione sfere
    COMPLEX(dbl), DIMENSION(2*nstop*(nstop+2)*ns), INTENT(IN) :: v_amnbmn		!Vettore campo esterno
    COMPLEX(dbl), DIMENSION(2*nstop*(nstop+2)*ns), INTENT(IN) :: v_dmncmn		!Vettore campo interno
    COMPLEX(dbl), DIMENSION(nstop*(nstop+2)), INTENT(IN) :: v_emn		!Vettore normalizzazione
    REAL(dbl), DIMENSION(ns,2), INTENT(IN) :: m_eps_inp		!Funzioni dielettriche corrette
    REAL(dbl), DIMENSION(ns), INTENT(IN) :: v_r_inp		!Raggi equivalenti
    REAL(dbl), INTENT(IN) :: ref_index				!Indice di rifrazione matrice
    INTEGER(lo), INTENT(IN) :: nstop,ratio,ns			!Numero di espansioni multipolari e raggi non calcolo
    REAL(dbl), INTENT(IN) :: x,y,z				!Coordinate punto nello spazio
    REAL(dbl), INTENT(INOUT) :: lambda_inp	        		! Vacuum Wavelength
    REAL(dbl), INTENT(IN) :: alpha0,beta0,gamma0		! Euler angles for incident light
	CHARACTER(len=3), INTENT(IN) :: quasi_static_flag				!Flag to force quasi-static incident field
    COMPLEX(dbl), INTENT(OUT) :: Ex,Ey,Ez			!Campo cartesiano
    REAL(dbl), INTENT(OUT) :: Eabs      			!Field module
    INTEGER(lo) , INTENT(OUT) :: error				!Flag di errore

    !Dichiarazione variabili interne
    COMPLEX(dbl) :: Exinc,Eyinc,Ezinc,pshift			!Onda incidente cartesiana
    REAL(dbl), DIMENSION(1:3) :: v_Einc=0.0D0,v_kinc=0.0D0	!Vettore campo e vettore d'onda per l'onda piana incidente
    COMPLEX(dbl) :: Er,Etheta,Ephi				!Campo calcolato in coordinate sferiche
    COMPLEX(dbl) :: rhoc,mc					!Rho complesso,indice di rifrazione complesso
    REAL(dbl) :: rhor,phase,k					!Rho reale
    REAL(dbl) :: r,theta,phi 					!Coordinate sferiche locali
    REAL(dbl) :: xl,yl,zl	 				!Coordinate cartesiane locali
    REAL(dbl) :: d,radius					!Distanza punto-sfere,raggio sfera
    INTEGER(lo) :: n,m,i,dimens,low_b,up_b,exit_flag		!Indici e dimensione
    REAL(dbl) :: nr,mr						!Indici reali

    !Inizio subroutine vera e propria

    !Calcolo il numero di sfere e dimens
    dimens=nstop*(nstop+2)
    exit_flag=0

    !--------------------------------------------
    !Se sono qui allora non sono in nessuna sfera
    !--------------------------------------------
    Ex=(0.0D0,0.0D0)
    Ey=(0.0D0,0.0D0)
    Ez=(0.0D0,0.0D0)

    ! Calcolo i miei vettori campo e vettore d'onda
    !Faccio k

	!quasi static approximation
	qs_if: IF (quasi_static_flag=='yes') THEN
		lambda_inp=MAXVAL(v_r_inp)*1.0D2
	END IF qs_if

    k=2*Pi_d*ref_index/lambda_inp
    v_kinc(1)=k*COS(alpha0)*SIN(beta0)
    v_kinc(2)=k*SIN(alpha0)*SIN(beta0)
    v_kinc(3)=k*COS(beta0)

    !Calcolo E
    v_Einc(1)=COS(alpha0)*COS(beta0)*COS(gamma0)-SIN(alpha0)*SIN(gamma0)
    v_Einc(2)=SIN(alpha0)*COS(beta0)*COS(gamma0)-COS(alpha0)*SIN(gamma0)
    v_Einc(3)=-SIN(beta0)*COS(gamma0)

    !Input monitoring
!!$    IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$
!!$       WRITE(*,*) 'v_kinc',v_kinc
!!$       WRITE(*,*) 'v_Einc',v_Einc
!!$       WRITE(*,*) 'flaginc',flaginc
!!$       WRITE(*,*) 'nstop',nstop
!!$       WRITE(*,*) 'ratio',ratio
!!$       WRITE(*,*) 'lambda_inp',lambda_inp
!!$       WRITE(*,*) 'x,y,z',x,y,z
!!$       WRITE(*,*) 'v_amnbmn',v_amnbmn
!!$       WRITE(*,*) 'v_dmncmn',v_dmncmn
!!$       WRITE(*,*) 'v_emn',v_emn
!!$       WRITE(*,*) 'm_xyz_inp',m_xyz_inp
!!$       WRITE(*,*) 'm_eps',m_eps_inp
!!$       WRITE(*,*) 'v_req',v_r_inp
!!$       WRITE(*,*) 'ref_index',ref_index
!!$
!!$    END IF

    !Faccio un ciclo per vedere se interseco sulle sfere
    int_do: DO i=1,ns

       !Calcolo raggio e distanza
       radius=v_r_inp(i)
       d=SQRT((x-m_xyz_inp(i,1))**2 + (y-m_xyz_inp(i,2))**2 + (z-m_xyz_inp(i,3))**2)

              !Se non sono in una sfera, controllo se sono in un'altra
       IF (d>radius) CYCLE int_do

       !-------------------------------------------
       !Se sono qui allora sono nella i-esima sfera
       !-------------------------------------------

!!$       WRITE(*,*) "sono nella sfera"

       !Calcolo le coordinate sferiche locali
       xl=x-m_xyz_inp(i,1)
       yl=y-m_xyz_inp(i,2)
       zl=z-m_xyz_inp(i,3)
       CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

       c_s_if: IF (error/=0) THEN
          error=1
          WRITE(*,10)
10        FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")
          RETURN
       END IF c_s_if

       !Calcolo indice di rifrazione complesso e poi rho per chiamare la subroutine per il calcolo del campo
       mc=SQRT(CMPLX(m_eps_inp(i,1),m_eps_inp(i,2),KIND=dbl))
       rhoc=2*Pi_D*mc*r/lambda_inp

       !Calcolo i bound per i coeff interni, e il campo in coordinate sferiche locali
       low_b=2*dimens*(i-1)+1
       up_b=2*dimens*i
       CALL int_field(v_dmncmn(low_b:up_b),v_emn,nstop,rhoc,theta,phi,Er,Etheta,Ephi,error)

       intfield_if: IF (error/=0) THEN
          error=1
          WRITE(*,20)
20        FORMAT ("Errore in int_field chiamata da Exyz_sub: il programma termina ora...")
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
       xl=x-m_xyz_inp(i,1)
       yl=y-m_xyz_inp(i,2)
       zl=z-m_xyz_inp(i,3)
       CALL cart_spher_r_dbl1(xl,yl,zl,r,theta,phi,error)

       c_s_if1: IF (error/=0) THEN
          WRITE(*,30)
30        FORMAT ("Errore in cart_spher chiamata da Exyz_sub: il programma termina ora...")
          RETURN
       END IF c_s_if1

!!$       !Input monitoring
!!$       IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "xl,yl,zl",xl,yl,zl
!!$          WRITE(*,*) "r,theta,phi",r,theta,phi
!!$       END IF

       IF ( ( r/( v_r_inp(i) ) ) > ratio) CYCLE ext_do

!!$       !Input monitoring
!!$       IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "computing external field, ns:",i
!!$       END IF

       !Calcolo  rho per chiamare la subroutine per il calcolo del campo
       rhor=2*Pi_D*ref_index*r/lambda_inp
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

!!$       !Input monitoring
!!$       IF ((ABS(x)<1.0D-10) .AND. ((ABS(y)<1.0D-10))) THEN
!!$          WRITE(*,*) "Er,Etheta,Ephi",Er,Etheta,Ephi
!!$       END IF

       extfield_if: IF (error/=0) THEN
          WRITE(*,40)
40        FORMAT ("Errore in ext_field chiamata da Exyz_sub: il programma termina ora...")
          RETURN
       END IF extfield_if

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

    !Field module
    Eabs=SQRT((ABS(Ex))**2+(ABS(Ey))**2+(ABS(Ez))**2)

  END SUBROUTINE Exyz


!-------------------------------------------------------------------------------
!4) Efar_sub: calculating scattering matrix and far field radiation pattern
!-------------------------------------------------------------------------------
SUBROUTINE Efar(k,ns,nstop,phimin,phimax,phistep,thetastep,betap,v_amnbmn,m_xyz_inp,m_SC,scatot,error)

USE local_kinds				! Kinds numerici, come single e double precision

IMPLICIT NONE

!Dichiarazione dei dummy argument
REAL(dbl), INTENT(IN) :: k						!Wavevector
INTEGER(lo), INTENT(IN) :: nstop,ns			!Numero di espansioni multipolari e raggi non calcolo
REAL(dbl), INTENT(IN) :: phimin,phimax					! Azimuth far field angles
REAL(dbl), INTENT(INOUT) :: betap						! Incident polarization angle
REAL(dbl), DIMENSION(ns,3), INTENT(IN) :: m_xyz_inp		!Posizione sfere
COMPLEX(dbl), DIMENSION(2*nstop*(nstop+2)*ns), INTENT(IN) :: v_amnbmn		!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				! Azimuth and polar steps
COMPLEX(dbl), DIMENSION(thetastep+1,phistep+1,7), INTENT(OUT) :: m_SC			! Scattering quantities
INTEGER(lo) , INTENT(OUT) :: error					! Error flag
REAL(dbl), INTENT(OUT) :: scatot					! Scattering cross section

!Dichiarazione variabili interne
COMPLEX(dbl) :: s1,s2,s3,s4,sum1,sum2,sum3,sum4					! Final and partial scatt. matrix elements
REAL(dbl) :: theta,phi,phase,dtheta,dphi					!
REAL(dbl) :: nr,mr								! Real indexes
REAL(dbl) :: logw,fr								! Auxiliary quantities
REAL(dbl) :: s11,s12,dcpar,dcnorm,dc,cs,ss,dc1					! Auxiliary quantities
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,chunk			! Auxiliary quantities
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_psimn, v_phimn, v_thetamn,	v_ximn	! Angular functions
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta					! Angular functions
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn,v_cmn			!Angular functions
REAL(dbl), DIMENSION(thetastep+1,phistep+1,2) :: m_thetaphi			! Angle storage matrix

! Subroutine
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_cmn(1:dimens))
ALLOCATE(v_psimn(1:side),v_phimn(1:side),v_thetamn(1:side),v_ximn(1:side))

! Angular increments
dphi=(phimax-phimin)/REAL(phistep,dbl)
dtheta=180.0D0/REAL(thetastep,dbl)
dtheta=dtheta*Pi_d/180.0D0
dphi=dphi*Pi_d/180.0D0


! Filling the angular matrix
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





	! Angular function normalization coefficients
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


	! Checking for well defined angles
	phi_if: IF (phimax<phimin) THEN
		error=1
		WRITE(*,10)
		10 FORMAT ("Error in Efar_sub: wrong phi bounds...")
		RETURN
	END IF phi_if


	! Checking for well defined angles
	steps_if: IF (thetastep<2) THEN
		error=1
		WRITE(*,30)
		30 FORMAT ("Error in Efar_sub: wrong theta steps...")
		RETURN
	END IF steps_if


	!Initializing the FIRSPRIVATE allocatable arrays for safety
	v_pimn=0.0D0
	v_taumn=0.0D0
	v_psimn=0.0D0
	v_phimn=0.0D0
	v_thetamn=0.0D0
	v_ximn=0.0D0

	! loop over theta
	theta_do: DO itheta=1,thetastep+1

		!Assigning theta
		theta=m_thetaphi(itheta,1,1)

		! Pi_mn angular function
		CALL pi_mnsca(nstop,theta,v_pimn,error)

		pimn_if: IF (error/=0) THEN
					WRITE(*,11)
					11 FORMAT ("Error in pi_mn called by Efar_sub: terminating...")
					STOP
		END IF pimn_if

		! Tau_mn angular function
		CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)

		taumn_if: IF (error/=0) THEN
					WRITE(*,20)
					20 FORMAT ("Error in tau_mn called by Efar_sub: terminating...")
					STOP
		END IF taumn_if

		! Normalizing angular functions
		v_taumn=v_cmn*v_taumn
		v_pimn=v_cmn*v_pimn

		! Initializing auxiliary vectors for scattering matrix
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

			! phi
			phi=m_thetaphi(itheta,iphi,2)

			! filling in phase terms
			delta_do: DO l=1,ns

				phase=m_xyz_inp(l,1)*SIN(theta)*COS(phi) + m_xyz_inp(l,2)*SIN(theta)*SIN(phi) + m_xyz_inp(l,3)*COS(theta)
				v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)

				!If ((itheta==1) .AND. (iphi)==1) THEN
				!	WRITE(*,*) v_delta(l)
				!END IF

			END DO delta_do

			! angular coordinates
			m_SC(itheta,iphi,1)=theta
			m_SC(itheta,iphi,2)=phi

			! loop to get S matrix terms
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

				! Amplitude scattering matrix elements
				s1=s1+v_delta(l)*sum1
				s2=s2+v_delta(l)*sum2
				s3=s3+v_delta(l)*sum3
				s4=s4+v_delta(l)*sum4

			END DO l_do_phi

			! muller matrix elements
			s11=0.5D0*((ABS(s1)**2)+(ABS(s2)**2)+(ABS(s3)**2)+(ABS(s4)**2))
			s12=0.5D0*((ABS(s2)**2)-(ABS(s1)**2)+(ABS(s4)**2)-(ABS(s3)**2))

			! parallel, normal and averange differenzial cross section
			dcpar= (s11+s12)/(k**2)
			dcnorm=(s11-s12)/(k**2)
			m_SC(itheta,iphi,7)=((COS(phi-betap))**2)*dcpar + ((SIN(phi-betap))**2)*dcnorm


			! Scattering matrix elements
			m_SC(itheta,iphi,3)=s1
			m_SC(itheta,iphi,4)=s2
			m_SC(itheta,iphi,5)=s3
			m_SC(itheta,iphi,6)=s4

		END DO phi_do
	END DO theta_do

	! Integrating the scattering cross section
	scatot=0.0D0
	sca_theta_do: DO itheta=1,thetastep+1
		sca_phi_do: DO iphi=1,phistep+1

			scatot=scatot+m_SC(itheta,iphi,7)*SIN(m_thetaphi(itheta,iphi,1))*dtheta*dphi

		END DO sca_phi_do
	END DO sca_theta_do

END SUBROUTINE Efar


!-------------------------------------------------------------------------------
!5) Efar_sub_poynting: calculating scattering matrix and far field radiation pattern
! with explicit poynting vector
!-------------------------------------------------------------------------------
SUBROUTINE Efar_poynting(k,ns,nstop,phimin,phimax,phistep,thetastep,betap,v_amnbmn,m_xyz_inp,m_SC,scatot,error)

USE local_kinds				! Kinds numerici, come single e double precision

IMPLICIT NONE

! dummy arguments
REAL(dbl), INTENT(IN) :: k							!Wavevector
INTEGER(lo), INTENT(IN) :: nstop,ns			!Numero di espansioni multipolari e raggi non calcolo
REAL(dbl), INTENT(IN) :: phimin,phimax					! azimuth angles
REAL(dbl), INTENT(INOUT) :: betap						! polarization angle
REAL(dbl), DIMENSION(ns,3), INTENT(IN) :: m_xyz_inp		!Posizione sfere
COMPLEX(dbl), DIMENSION(2*nstop*(nstop+2)*ns), INTENT(IN) :: v_amnbmn		!Vettore campo esterno
INTEGER(lo), INTENT(IN) :: phistep,thetastep				! azimuth and polar steps
COMPLEX(dbl), DIMENSION(thetastep+1,phistep+1,7), INTENT(OUT) :: m_SC			! scattering outputs
INTEGER(lo) , INTENT(OUT) :: error					! error flag
REAL(dbl), INTENT(OUT) :: scatot						! total scattering cross section

!Dichiarazione variabili interne
COMPLEX(dbl) :: Etheta,Ephi,Htheta,Hphi,sum1,sum2,sum3,sum4,dccomplex					! fields and scattering partial sums
REAL(dbl) :: theta,phi,phase,dtheta,dphi								! angle and angular increments
REAL(dbl) :: dc												! differential cross section
INTEGER(lo) :: i,l,n,m,j0,j,jm,dimens,side,iphi,itheta,chunk						! auxiliary variables
COMPLEX(dbl), DIMENSION(1:ns) :: v_delta								! phase terms
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_pimn,v_taumn							! angular functions
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_amn_taumn,v_amn_pimn,v_bmn_taumn,v_bmn_pimn,v_emn,v_ab_tot	! angular functions
REAL(dbl), DIMENSION(thetastep+1,phistep+1,2) :: m_thetaphi						! angle storage

! subroutine
dimens=nstop*(nstop+2)
side=ns*dimens
ALLOCATE(v_pimn(1:dimens),v_taumn(1:dimens),v_emn(1:dimens),v_ab_tot(1:2*dimens))
ALLOCATE(v_amn_pimn(1:dimens),v_amn_taumn(1:dimens),v_bmn_pimn(1:dimens),v_bmn_taumn(1:dimens))

! angular increments
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


! field normalization coefficients
error=0
CALL emn_sub(nstop,v_emn,error)

! well defined theta steps?
steps_if: IF (thetastep<2) THEN
	error=1
	WRITE(*,30)
	30 FORMAT ("Error in Efar_sub: ill defined theta steps., Terminating...")
	RETURN
END IF steps_if

scatot=0.0D0


!Initializing allocatable arrays for safety
v_pimn=0.0D0
v_taumn=0.0D0
v_ab_tot=0.0D0
v_amn_pimn=0.0D0
v_amn_taumn=0.0D0
v_bmn_pimn=0.0D0
v_bmn_taumn=0.0D0

! angular functions
theta_do: DO itheta=1,thetastep+1

	theta=m_thetaphi(itheta,1,1)

	! Pi_mn angular function
	CALL pi_mnsca(nstop,theta,v_pimn,error)

	pimn_if: IF (error/=0) THEN
		WRITE(*,11)
		11 FORMAT ("Error in pi_mn called by Efar_sub: stopping...")
		STOP
	END IF pimn_if

	! Tau_mn angular function
	CALL tau_mnsca(nstop,theta,v_pimn,v_taumn,error)

	taumn_if: IF (error/=0) THEN
		WRITE(*,20)
		20 FORMAT ("Error in tau_mn called by Efar_sub: stopping...")
		STOP
	END IF taumn_if


	phi_do: DO iphi=1,phistep+1

		phi=m_thetaphi(itheta,iphi,2)

		! phase terms
		delta_do: DO l=1,ns

			phase=m_xyz_inp(l,1)*SIN(theta)*COS(phi) + m_xyz_inp(l,2)*SIN(theta)*SIN(phi) + m_xyz_inp(l,3)*COS(theta)
			v_delta(l)=EXP(-(0.0D0,1.0D0)*k*phase)

		END DO delta_do

		! angular coordinates
		m_SC(itheta,iphi,1)=theta
		m_SC(itheta,iphi,2)=phi

		! total scattering expansion coefficients in the far field
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


		! mixed vectors for far field calculations

		v_amn_taumn=v_ab_tot(1:2*dimens:2)*v_taumn
		v_bmn_taumn=v_ab_tot(2:2*dimens:2)*v_taumn
		v_amn_pimn=v_ab_tot(1:2*dimens:2)*v_pimn
		v_bmn_pimn=v_ab_tot(2:2*dimens:2)*v_pimn

		! computing the fields
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

		! flipping the fields where necessary
		Etheta=-Etheta
		Ephi=-Ephi
		Hphi=-Hphi

		! complex and real differential cross section
		dccomplex=(Etheta*CONJG(Hphi)-Ephi*CONJG(Htheta))/(k**2)
		dc=REAL(dccomplex,dbl)

		! fields and cross section
		m_SC(itheta,iphi,3)=REAL(Etheta,dbl)
		m_SC(itheta,iphi,4)=REAL(Ephi,dbl)
		m_SC(itheta,iphi,5)=REAL(Htheta,dbl)
		m_SC(itheta,iphi,6)=REAL(Hphi,dbl)
		m_SC(itheta,iphi,7)=dc

	END DO phi_do
END DO theta_do

! integrating the scattering cross section
scatot=0.0D0
sca_theta_do: DO itheta=1,thetastep+1
	sca_phi_do: DO iphi=1,phistep+1

		scatot=scatot+m_SC(itheta,iphi,7)*SIN(m_thetaphi(itheta,iphi,1))*dtheta*dphi

	END DO sca_phi_do
END DO sca_theta_do

END SUBROUTINE Efar_poynting

!--------------------------------------------------------------------------------
!6) SUBROUTINE dip_coefficients: dipole multipolar coefficients
!--------------------------------------------------------------------------------
SUBROUTINE dip_coefficients(ns,ndip,m_xyz_inp,v_r_inp,m_eps_inp,fint,q0,&
                               &ref_index,lambda_inp,alpha0,beta0,gamma0,nstop,&
							 &quasi_static_flag,&
                               &m_coeff,v_W)

USE local_kinds				! Kinds numerici, come single e double precision


IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Dummy arguments declaration
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dbl), DIMENSION(ns,3), INTENT(IN) :: m_xyz_inp			! Sphere coordinates
REAL(dbl), DIMENSION(ns), INTENT(IN) :: v_r_inp			! Sphere radii
REAL(dbl), DIMENSION(ns,2), INTENT(IN) :: m_eps_inp			! Sphere e1,e2
REAL(dbl), INTENT(IN) :: fint,q0                                      ! Interaction parameter, q0
REAL(dbl), INTENT(IN) :: ref_index					! Matrix Refractive Index
REAL(dbl), INTENT(INOUT) :: lambda_inp					! Vacuum Wavelength
REAL(dbl), INTENT(INOUT) :: alpha0,beta0,gamma0			! Euler angles for incident light
INTEGER(lo), INTENT(IN) :: ns,ndip,nstop				! Spheres, dipoles and Multipolar expansions
CHARACTER(len=3), INTENT(IN) :: quasi_static_flag				!Flag to force quasi-static incident field
COMPLEX(dbl), DIMENSION(2*nstop*(nstop+2)*ns,2), INTENT(OUT):: m_coeff! Expansion coefficients
REAL(dbl), DIMENSION(6), INTENT(OUT) :: v_W         ! Radiated power and efficiencies

!~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Internal Variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~

! File opening
INTEGER(lo) :: status_in, status_arr		! Status apertura file principale e allocazione array
INTEGER(lo) :: neq,ns_ant,matrixside,i,j,m	! Numero sfere totali, numero materiali e numero gruppi di sfere uguali
INTEGER(lo) :: matrixside_ant         	! Numero sfere totali, numero materiali e numero gruppi di sfere uguali
INTEGER(lo) :: nnb,nnb_ant,NNZab,NNZd,blockside	! Non zero blocks, elements Aij Bij, elements Dij, block side
INTEGER(lo) :: error			        ! Flag di errore

REAL(dbl) :: t1,t2				! Tempi
REAL(dbl) :: lambda0,lambda1,lambda2,k,m1,m2	! Lunghezza d'onda
REAL(dbl) :: lambdamin,lambdamax,lambdastep	! Lunghezza d'onda

! Single particle variables
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_da,m_cb! Matrici coefficienti calcolo campo interno

! General internal variables
INTEGER(lo) :: sumqsp,sumqqsp				!Somme di qmax e Qmax nel caso sparse
INTEGER(lo), ALLOCATABLE, DIMENSION(:,:) :: m_index	! Index matrix
REAL(lo), ALLOCATABLE, DIMENSION(:,:) :: m_Apmn	! Apmn coefficient matrix
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_c0,v_aq,v_bq!Vettori coefficienti c0,gaunt e bq per Aij Bij

COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_p,v_p1,v_ab,v_ab1,v_ab_sca,v_dc,v_diag! Vettori espansioni campo e vettore noto
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_ab_sca_ant					! Vettori espansioni campo
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_csca,v_csca_ant,v_csca_dip,v_cext,v_cabs		! Sezioni
INTEGER(lo) :: info,n								! Flag errore
REAL(dbl) :: dipsca

CHARACTER(len=3) :: norm='new'						! Flag normalizz fill_AB
CHARACTER(len=3) :: norm1='div'						! Flag normalizz matrice

! Linear solver variables
INTEGER(lo) :: lvanilla=2,mxmv,mxmv0,cont		! Bicgstab dimension and max iterations
REAL(dbl) :: tol,tol0					! Bicgstab tolerance
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: work	! Auxiliary array for bicgstab
LOGICAL :: logic1,logic2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Dirty hacks to interface fortran and python properly
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
neq=ns
ALLOCATE (v_patt(1:ns),STAT=status_arr)
DO i=1,ns
   v_patt(i)=i
END DO


!!$
!!$  WRITE(*,*) 'ns',ns
!!$  WRITE(*,*) 'm_xyz_inp',m_xyz_inp
!!$  WRITE(*,*) 'v_r_inp',v_r_inp
!!$  WRITE(*,*) 'm_eps_inp',m_eps_inp
!!$  WRITE(*,*) 'fint',fint
!!$  WRITE(*,*) 'ref_index',ref_index
!!$  WRITE(*,*) 'lambda_inp',lambda_inp
!!$  WRITE(*,*) 'alpha0,beta0,gamma0',alpha0,beta0,gamma0
!!$  WRITE(*,*) 'nstop',nstop
!!$  WRITE(*,*) 'v_patt',v_patt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Checking sphere intersection
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALL intspheres(v_r_inp,m_xyz_inp,ns,error)

err_intspheres: IF (error/=0) THEN
   WRITE(*,*)
   WRITE(*,12)
12   FORMAT ("Intersection spheres, the program stops.")
   STOP
END IF err_intspheres


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Non zero blocks in coeff matrix
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Assigning new values to spheres and dipole numbers
ns_ant=ns-ndip		!Number of real spheres
CALL nnb_sparse(ns,m_xyz_inp,v_r_inp,fint,nnb) !Divided by two because i store only the upper triangular part
nnb=nnb/2
!!$  WRITE(*,*) "Non zero blocks:", nnb

!Non zero blocks fincluding incident dipoles, computed already correctly for the new storage scheme
CALL nnb_sparse_ant(ns,ns_ant,m_xyz_inp,v_r_inp,fint,nnb_ant)

!Initializing v_ic to eliminate the load on the intrinsic library pow
ic_do: DO i=0,400
   v_ic(i)=(0.0D0,1.0D0)**i
END DO ic_do

!Initializing v_ir to eliminate the load on the intrinsic library pow
ir_do: DO i=-200,200
   v_oner(i)=(-1.0D0)**i
END DO ir_do



!nnb multiplied by 2,because i store the whole skeleton of the matrix for flexibility
ALLOCATE (v_jBlock(1:nnb),v_iBlock(1:ns+1),STAT=status_arr)
CALL fill_jBlock_iBlock_sparse(ns,m_xyz_inp,v_r_inp,fint,v_jBlock,v_iBlock)

!--------------Filling diagnostics-----------------------------
!!$  WRITE(*,*) 'v_iBlock',v_iBlock
!!$  WRITE(*,*) 'v_jBlock',v_jBlock
!--------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Elements for Aij e Dij, matrix sizes and so on
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Non zero elements Aij Bij
NNZab=nstop*(1+2*nstop*(3+nstop))
NNZab=NNZab/3
!Non zero elements Dij
NNZd=nstop*(11+12*nstop+4*(nstop**2))
NNZd=NNZd/3
!Lato di un blocco
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside
matrixside_ant=2*ns_ant*blockside

!---------------------------------
!Allocating matrixes and vectors
!---------------------------------
!Coefficient Matrix
ALLOCATE (m_Dij(1:NNZd,1:nnb),m_jDij(1:NNZd,1:nnb),m_iDij(1:blockside+1,1:nnb),&
         &m_exphi(-nstop:nstop,1:nnb),STAT=status_arr)
ALLOCATE(m_sAij(1:NNZab,1:nnb_ant),m_sBij(1:NNZab,1:nnb_ant),&
        &m_sAij_sca(1:NNZab,1:nnb),m_sBij_sca(1:NNZab,1:nnb),&
        &m_jABij(1:NNZab,1:nnb),m_iABij(1:blockside+1,1:nnb),STAT=status_arr)

!Expansion coefficient vectors
ALLOCATE (v_p(matrixside),v_p1(matrixside),v_ab(matrixside),v_ab1(matrixside),v_ab_sca(matrixside),&
          v_ab_sca_ant(matrixside),v_dc(matrixside),v_diag(matrixside),STAT=status_arr)

!Cross sections
ALLOCATE(v_cext(1:ns),v_cabs(1:ns),v_csca(1:ns),v_csca_ant(1:ns),v_csca_dip(1:ns),STAT=status_arr)

!Workspace for the linear solver
ALLOCATE(work(1:matrixside_ant,1:3+2*(lvanilla+1)),STAT=status_arr)


!--------------------------------------
!Filling rotational coefficient matrix
!--------------------------------------
CALL fill_D_PHI_sparse(ns,nstop,nnb,m_xyz_inp,v_jBlock,v_iBlock,m_Dij,m_jDij,m_iDij,m_exphi)

!!$  WRITE(*,*) 'm_Dij',m_Dij
!!$  WRITE(*,*) 'm_jDij',m_jDij
!!$  WRITE(*,*) 'm_iDij',m_iDij
!!$  WRITE(*,*) 'm_exphi',m_exphi

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Aij Bij wavelength independent coefficients
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Coefficienti c0
ALLOCATE(v_c0(1:NNZab))
CALL c0_sparse(nstop,v_c0)

!Coefficienti aq e bq
CALL qsum_sparse(nstop,sumqsp,sumqqsp)
ALLOCATE(v_aq(1:sumqsp))
ALLOCATE(v_bq(1:sumqqsp))
CALL fill_aqbq_sparse(nstop,v_aq,v_bq)

!Calculation support coefficients for ABij once for all
ALLOCATE(m_index(1:NNZab,1:11),STAT=status_arr)
ALLOCATE(m_Apmn(1:NNZab,1:(2*nstop)+1),STAT=status_arr)
ALLOCATE(v_jABij_template(1:NNZab),v_iABij_template(1:blockside+1),v_mask(1:NNZab),STAT=status_arr)
CALL fill_index(nstop,m_index,m_Apmn,v_jABij_template,v_iABij_template,v_mask,error)
!		Write(*,600) m_index(:,4)
600 FORMAT (I5)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Allocating what is left, i.e. expansion
!coefficients, cross sections, etc...
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Lato di un blocco
blockside=nstop*(nstop+2)
matrixside=2*ns*blockside
imin=1
imax=ns
jmin=1
jmax=ns

!quasi static approximation
qs_if: IF (quasi_static_flag=='yes') THEN
lambda_inp=MAXVAL(v_r_inp)*1.0D2
END IF qs_if

!Wavevector
k=2*Pi_d*ref_index/(lambda_inp)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Single particle coefficients
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ALLOCATE(m_a(1:nstop,1:neq-ndip),m_b(1:nstop,1:neq-ndip),&
        &m_da(1:nstop,1:neq-ndip),m_cb(1:nstop,1:neq-ndip),STAT=status_arr)

CALL coeff_sp2(lambda_inp,ref_index,v_r_inp(1:ns-ndip),m_eps_inp(1:ns-ndip,:),nstop,neq-ndip,m_a,m_b,error)

coeff_sp_if: IF (error/=0) THEN
   WRITE(*,70)
70   FORMAT ("Error in coeff_sp2, stopping the program...")
   STOP
END IF coeff_sp_if

CALL coeff_ad_ca(lambda_inp,ref_index,v_r_inp(1:ns-ndip),m_eps_inp(1:ns-ndip,:),nstop,neq-ndip,m_da,m_cb,error)

coeff_dacb_if: IF (error/=0) THEN
   WRITE(*,71)
71   FORMAT ("Error in coeff_ad_ca, stopping the program a...")
   STOP
END IF coeff_dacb_if


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Filling the coefficient matrix
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALL fill_AB_sparse_dip_fast(ns,ns_ant,nstop,k,m_xyz_inp,v_c0,v_aq,v_bq,&
     m_index,m_Apmn,v_jABij_template,v_iABij_template,&
     v_jBlock,v_iBlock,m_sAij,m_sBij,m_sAij_sca,m_sBij_sca,m_jABij,m_iABij)

!Filling the diagonal of the matrix
j=0
v_diag=1.0D0
DO i=1,ns_ant
   DO n=1,nstop
      DO m=-n,n
         j=j+1
         v_diag(2*j-1)=1.0D0/m_a(n,v_patt(i))
         v_diag(2*j)=1.0D0/m_b(n,v_patt(i))
      END DO
   END DO
END DO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Filling the RHS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Incident field expansion coefficients
CALL field_expRandom_sub(nstop,ns,k,zerodbl,m_xyz_inp,alpha0,beta0,gamma0,v_p,error)


!Saving a copy of the fictious incident field
v_p1=v_p

!Computing the expansion coefficients of the dipoles in their own original reference frames
j=0
v_ab=0.0D0
ns_rhs_do: DO i=1,ns
   n_rhs_do: DO n=1,nstop
      m_rhs_do: DO m=-n,n
         j=j+1
         IF ((i>ns_ant) .AND. (n==1)) THEN
            v_ab(2*j-1)=a1mie*v_p(2*j-1)
         END IF
      END DO m_rhs_do
   END DO n_rhs_do
END DO ns_rhs_do



!Computing the RHS by traslating the dipoles in the antenna sphere reference frame
CALL matvec_gen(nstop,blockside,1,ns_ant,ns_ant+1,ns,&
     v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_sAij,m_sBij,m_iABij,m_jABij,v_diag,&
     v_ab(matrixside_ant+1:matrixside),v_p(1:matrixside_ant),error)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Solving the linear system
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!Initial guess
j=0
DO i=1,ns_ant
   DO n=1,nstop
      DO m=-n,n
         j=j+1
         v_ab(2*j-1)=-m_a(n,v_patt(i))*v_p(2*j-1)
         v_ab(2*j)=-m_b(n,v_patt(i))*v_p(2*j)
      END DO
   END DO
END DO

solve_if: IF (ns_ant/=1) THEN

   !BicgSTAB inputs
   work=(1.0d0,1.0d0)
   logic1=.TRUE.
   logic2=.TRUE.
   mxmv=20000
   tol=1.0D-15
!!$     WRITE(*,*) "l,lambda,mxmv,tol: " ,lvanilla,lambda_inp,mxmv,tol
!!$     WRITE(*,*)
!!$
!!$     !~~~~~~~~~~~~~~~~~~~
!!$     ! INPUT DIAGNOSTICS
!!$     !~~~~~~~~~~~~~~~~~~~
!!$     !------Logic Variables--------
!!$     WRITE(*,*) "logic1",logic1
!!$     WRITE(*,*) "logic2",logic2
!!$     !------Integer Parameters---------
!!$     WRITE(*,*) "nstop",nstop
!!$     WRITE(*,*) "lvanilla",lvanilla
!!$     WRITE(*,*) "matrixside",matrixside
!!$     WRITE(*,*) "blockside",blockside
!!$     !------Block Matrix Parameters-----------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_iBlock"
!!$     DO i=1,ns+1
!!$        WRITE(*,205) v_iBlock(i)
!!$205     FORMAT (I)
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_jBlock"
!!$     DO i=1,nnb
!!$        WRITE(*,210) v_jBlock(i)
!!$210     FORMAT (I)
!!$     END DO
!!$
!!$     !------D Matrix Parameters------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_iD"
!!$     DO i=1,nnb
!!$        DO j=1,ns+1
!!$           WRITE(*,215) m_iDij(j,i)
!!$215        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_jD"
!!$     DO i=1,nnb
!!$        DO j=1,NNZd
!!$           WRITE(*,220) m_jDij(j,i)
!!$220        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_D"
!!$     DO i=1,nnb
!!$        DO j=1,NNZd
!!$           WRITE(*,225) m_Dij(j,i)
!!$225        FORMAT (1ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     !------AB Matrix Parameters---------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_iAB"
!!$     DO i=1,nnb
!!$        DO j=1,ns+1
!!$           WRITE(*,230) m_iABij(j,i)
!!$230        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_jAB"
!!$     DO i=1,nnb
!!$        DO j=1,NNZab
!!$           WRITE(*,235) m_jABij(j,i)
!!$235        FORMAT (I)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_A"
!!$     DO i=1,nnb
!!$        DO j=1,NNZab
!!$           WRITE(*,240) m_sAij(j,i)
!!$240        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_B"
!!$     DO i=1,nnb
!!$        DO j=1,NNZab
!!$           WRITE(*,245) m_sBij(j,i)
!!$245        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     !------Phi Matrix Parameters-----------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "m_phi"
!!$     DO i=1,nnb
!!$        DO j=-nstop,nstop
!!$           WRITE(*,250) m_exphi(j,i)
!!$250        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     !------Matrix Diagonal-----------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_diag"
!!$     DO i=1,matrixside
!!$        WRITE(*,255) v_diag(i)
!!$255     FORMAT (2ES16.7)
!!$     END DO
!!$
!!$     !------Unknown vector------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_ab"
!!$     DO i=1,matrixside
!!$        WRITE(*,260) v_ab(i)
!!$260     FORMAT (2ES16.7)
!!$     END DO
!!$
!!$     !RHS----------------------------------------------------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "v_p"
!!$     DO i=1,matrixside
!!$        WRITE(*,265) v_p(i)
!!$265     FORMAT (2ES16.7)
!!$     END DO
!!$
!!$     !------workspace-----------
!!$     WRITE(*,*)
!!$     WRITE(*,*) "work"
!!$     DO j=1,3+2*(lvanilla+1)
!!$        DO i=1,matrixside_ant
!!$           WRITE(*,270) work(i,j)
!!$270        FORMAT (2ES16.7)
!!$        END DO
!!$     END DO
!!$
!!$     WRITE(*,*)
!!$     WRITE(*,*) "tolerance"
!!$     WRITE(*,275) tol
!!$275  FORMAT (1ES16.7)


   !BicgSTAB linear solver: calling the solver
   CALL zbcg2_gen(nstop,logic1,logic2,&
        lvanilla,matrixside_ant,blockside,1,ns_ant,1,ns_ant,&
        v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_sAij,m_sBij,m_iABij,m_jABij,v_diag,&
        v_ab(1:matrixside_ant),-v_p(1:matrixside_ant),tol,mxmv,work,info)

!!$     WRITE(*,*) 'just outside the solver'

   info_if: IF (info /= 0) THEN
      WRITE(*,*) "Problem in bicgstab", info
      STOP
   END IF info_if

END IF solve_if

!----------------------------------------------------------
! Cross sections, or better radiated and dissipated powers
!----------------------------------------------------------

!!$  WRITE(*,*) 'before first matvec'

!Expansion coefficients for Total radiated power
v_diag=1.0D0
CALL matvec_gen(nstop,blockside,1,ns,1,ns,&
     v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_sAij_sca,m_sBij_sca,m_iABij,m_jABij,v_diag,&
     v_ab,v_ab_sca,error)

!!$  WRITE(*,*) 'before second matvec'

!Expansion coefficients for antenna radiated power
CALL matvec_gen(nstop,blockside,1,ns_ant,1,ns_ant,&
     v_iBlock,v_jBlock,m_Dij,m_iDij,m_jDij,m_exphi,m_sAij_sca,m_sBij_sca,m_iABij,m_jABij,v_diag,&
     v_ab(1:matrixside_ant),v_ab_sca_ant(1:matrixside_ant),error)

!!$  WRITE(*,*) 'after second matvec'

!Isolated dipoles radiated power (subroutines mutuated by the shell calculations)
v_csca_dip=0.0D0
CALL rad_shell_sub_dip(nstop,ns-ns_ant,k,v_ab(matrixside_ant+1:matrixside),v_csca_dip(ns_ant+1:ns))
! dipsca=(4*Pi_d/(k**2))*SUM(v_csca_dip)
dipsca=SUM(v_csca_dip)

!!$  WRITE(*,*) 'after isolated dipole'

!Antenna extinction
v_cext=0.0D0
CALL cext_random_sub(nstop,ns_ant,k,-v_p(1:matrixside_ant),v_ab(1:matrixside_ant),v_cext(1:ns_ant))

!!$  WRITE(*,*) 'after ext'

!Antenna radiated power
v_csca_ant=0.0D0
CALL csca_random_sub(nstop,ns_ant,k,v_ab_sca_ant(1:matrixside_ant),v_ab(1:matrixside_ant),v_csca_ant(1:ns_ant))

!!$  WRITE(*,*) 'after scaant'

!Antenna dissipated power
v_cabs=0.0D0
CALL dmncmn_sub(nstop,ns_ant,v_patt,m_da,m_cb,v_ab(1:matrixside_ant),v_dc(1:matrixside_ant))
CALL cabs_random_sub(ref_index,nstop,ns_ant,k,v_dc(1:matrixside_ant),&
     v_r_inp(1:neq-1),m_eps_inp(1:neq-1,:),v_cabs(1:ns_ant))

!!$  WRITE(*,*) 'after abs'

!Total radiated power
CALL csca_random_sub(nstop,ns,k,v_ab_sca,v_ab,v_csca)

!!$  WRITE(*,*) 'after sca'

!-------------
!Final outputs
!-------------
v_W(1)=SUM(v_csca)/dipsca !W_rad_tot: Power radiated by the dipole antenna system
v_W(2)=SUM(v_csca_ant)/dipsca !W_rad_ant: power radiated by the antenna alone
v_W(3)=(SUM(v_cext)-SUM(v_csca_ant))/dipsca !W_abs_1: Power dissipated: extinction method calculation
v_W(4)=SUM(v_cabs)/dipsca !W_abs_2: Power dissipated, absorption calculation
v_W(5)=(SUM(v_cabs)+dipsca*(1-q0)*100.0D0)/dipsca !W_nr non radiative recomb rate
v_W(6)=(SUM(v_csca)/dipsca)/(SUM(v_cabs)/dipsca+SUM(v_csca)/dipsca+(1-q0)/q0) !q, final quantum efficiency

!!$  WRITE(*,*) 'after outputs W'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Filling m_coeff with all the coefficients
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!ab
m_coeff(:,1)=v_ab
!dc
m_coeff(:,2)=v_dc

!!$  WRITE(*,*) 'after outputs coeff'

DEALLOCATE(v_patt)
!!$  WRITE(*,*) "1"
DEALLOCATE (v_jBlock,v_iBlock)
!!$  WRITE(*,*) "2"
DEALLOCATE (m_sAij,m_sBij,m_sAij_sca,m_sBij_sca,m_jABij,m_iABij)
!!$  WRITE(*,*) "3"
DEALLOCATE (m_Dij,m_jDij,m_iDij,m_exphi,STAT=status_arr)
!!$  WRITE(*,*) "4"
DEALLOCATE(work)
!!$  WRITE(*,*) "5"
DEALLOCATE(v_c0)
!!$  WRITE(*,*) "6"
DEALLOCATE(v_aq)
!!$  WRITE(*,*) "7"
DEALLOCATE(v_bq)
!!$  WRITE(*,*) "8"
DEALLOCATE(m_index)
!!$  WRITE(*,*) "9"
DEALLOCATE(m_Apmn)
!!$  WRITE(*,*) "10"
DEALLOCATE(v_jABij_template,v_iABij_template,v_mask)
!!$  WRITE(*,*) "11"
DEALLOCATE (v_p,v_p1,v_ab,v_diag,v_dc,v_ab_sca,v_ab_sca_ant,v_ab1)
!!$  WRITE(*,*) "12"
DEALLOCATE(v_csca,v_cabs,v_csca_ant,v_cext,v_csca_dip)
!!$  WRITE(*,*) "13"
DEALLOCATE(m_a,m_b,m_da,m_cb)
!!$  WRITE(*,*) "14"

!!$    WRITE(*,*) 'after deallocation'

END SUBROUTINE dip_coefficients





END MODULE
