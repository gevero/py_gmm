MODULE shared_data

USE KINDS

IMPLICIT NONE
SAVE
!********************************************************************************************************************
!********************************************************************************************************************
!********************************************************************************************************************
! Shared Variables
!********************************************************************************************************************
!********************************************************************************************************************
!********************************************************************************************************************
INTEGER(lo) :: nstop,ns,nsc,ncell,nacell,nbcell,ndip,ns_ant,imin,imax,jmin,jmax,namin,namax,nbmin,nbmax
REAL(dbl) :: ax,ay,az,bx,by,bz
COMPLEX(dbl) :: a1mie=-2.0D0*(0.0D0,1.0D0)/3.0D0
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jBlock,v_iBlock,v_jBlock_rhs,v_iBlock_rhs,v_jBlock_sca,v_iBlock_sca
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_sAij,m_sBij,m_sAij_sca,m_sBij_sca,m_sAij_rhs,m_sBij_rhs
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_sT1ij,v_sT2ij,v_sU1ij,v_sU2ij,v_sTUij,v_sAij,v_sBij,v_sAijh,v_sBijh
INTEGER(lo), ALLOCATABLE, DIMENSION(:,:) :: m_jABij,m_iABij,m_jABij_rhs,m_iABij_rhs
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_jABij,v_iABij,v_jABij_template,v_iABij_template,v_jTUij,v_iTUij,v_jDij,v_iDij
REAL(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_Dij,m_Dij_rhs,m_Dij_sca
REAL(dbl), ALLOCATABLE, DIMENSION(:) :: v_sDij
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_mask
INTEGER(lo), ALLOCATABLE, DIMENSION(:,:) :: m_jDij,m_iDij,m_jDij_rhs,m_iDij_rhs,m_jDij_sca,m_iDij_sca
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_precond_shell					! Vettore preconditioning codice shell
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_exphi, m_exphi_rhs
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_exphi
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_a,m_b,m_c,m_d           !Matrici coefficienti particella singola
INTEGER(lo), ALLOCATABLE ,DIMENSION(:) :: v_patt                      !Array per il pattern di uguaglianza
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:,:) :: m_AB,m_ABJ,m_AB_swap
REAL(dbl), DIMENSION(-200:200) :: v_oner				!Real vector for the 1.0D0 exponential
REAL(dbl), DIMENSION(0:120,0:120) :: m_fact			!Real vector for the 1.0D0 exponential
COMPLEX(dbl), DIMENSION(0:400) :: v_ic					!Complex vector for the (0.0D0,1.0D0) exponential

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Core shell codes shared data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
INTEGER(lo), ALLOCATABLE, DIMENSION(:) :: v_iAB,v_jAB,v_iAABB,v_jAABB,v_iM,v_jM,v_iM_dip,v_jM_dip
COMPLEX(dbl), ALLOCATABLE, DIMENSION(:) :: v_sA,v_sB,v_sAABB,v_sM,v_sM_dip
END MODULE shared_data
