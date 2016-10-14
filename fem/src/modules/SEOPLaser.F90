!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Geoffry Schrank
! *  Email:   geoffry.schrank@duke.edu
! *  Web:     NA
! *  Address: 311 Research Drive
! *           Durham, NC
! *           27710 U.S.A.
! *
! *  Original Date: 29 Sept. 2016
! *
!------------------------------------------------------------------------------
! * This solver should implement a module to model the absorption of laser light
! * though an SEOP cell. The output function will be the pump rate of the laser
! * scalar field over the model. The modeled equation is taken from Fink et al.
! * PRA 72 (2005); equation A6.

! * The template for this module is based on the ModelPDE.F90. I will add notes
! * as I figure out what I'm doing.
!-------------------------------------------------------------------------------
SUBROUTINE SEOPLaser( Model,Solver,dt,TransientSimulation )
    !------------------------------------------------------------------------------
    !Use standard Elmer tools
    USE DefUtils

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver !Solver and model data structs, I don't know what they do
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt !Probably time delta definition
    LOGICAL :: TransientSimulation !Yes/No Transient Simulation?
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element !pointer to element struct, elements of mesh?
    REAL(KIND=dp) :: Norm
    INTEGER :: n, nb, nd, t, active !element counters?
    INTEGER :: iter, maxiter !loop counters?
    LOGICAL :: Found !Yes/No the guy is in the SIF
    !------------------------------------------------------------------------------
    !ListGetInterger looks to be a function that retrieves information from the
    !SIF, presumably only for integers
    !GetSolverParams() looks to be a function that function that
    maxiter = ListGetInteger( GetSolverParams(),&
        'Nonlinear System Max Iterations',Found,minv=1)

    IF(.NOT. Found ) maxiter = 1

    ! Nonlinear iteration loop:
    !--------------------------
    DO iter=1,maxiter

        ! System assembly:
        !----------------
        CALL DefaultInitialize()
        Active = GetNOFActive()
        DO t=1,Active
            Element => GetActiveElement(t)
            n  = GetElementNOFNodes()
            nd = GetElementNOFDOFs()
            nb = GetElementNOFBDOFs()
            CALL LocalMatrix(  Element, n, nd+nb )
        END DO

        CALL DefaultFinishBulkAssembly()

        Active = GetNOFBoundaryElements()
        DO t=1,Active
            Element => GetBoundaryElement(t)
            IF(ActiveBoundaryElement()) THEN
                n  = GetElementNOFNodes()
                nd = GetElementNOFDOFs()
                nb = GetElementNOFBDOFs()
                CALL LocalMatrixBC(  Element, n, nd+nb )
            END IF
        END DO

        CALL DefaultFinishBoundaryAssembly()
        CALL DefaultFinishAssembly()
        CALL DefaultDirichletBCs()

        ! And finally, solve:
        !--------------------
        Norm = DefaultSolve()

        IF( Solver % Variable % NonlinConverged == 1 ) EXIT

    END DO

CONTAINS

    ! Assembly of the matrix entries arising from the bulk elements
    !------------------------------------------------------------------------------
    SUBROUTINE LocalMatrix( Element, n, nd )
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: nRb(n), Beta(n), &
            K, direction(3,n),N(3), Weight
        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
        REAL(KIND=dp) :: STIFF(nd,nd), MASS(nd,nd), FORCE(nd,nd), LOAD(nd,nd)

        LOGICAL :: Stat,Found
        INTEGER :: i,t,p,q,dim
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(ValueList_t), POINTER ::  Material
        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------

        dim = CoordinateSystemDimension()

        CALL GetElementNodes( Nodes )

        STIFF = 0._dp

        !Material Properties Def
        !Eventually some of this will all be passed from another solver
        Material => GetMaterial()
        nRb(1:n)=GetReal(Material,'rubidium number density',Found)

        DO i =1,dim
            direction(i,1:n)=GetReal(Material,&
                'propogation direction'//TRIM(I2S(i)), Found)
        END DO

        Beta(n) = BetaCalc(Model,n) !Calculate Beta from (A6), I have a feeling this isn't expressed correctly

        ! Numerical integration:
        !-----------------------
        IP = GaussPoints( Element )
        DO t=1,IP %
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

        ! The source term at the integration point:
        !------------------------------------------



        N = MATMUL(direction(:,1:n),Basis(1:n))
        K = SUM(Basis(1:n)*Beta(1:n)) !Need to write the full non-linear version of the code.

        Weight = IP % s(t) * DetJ



        DO p=1,nd
            DO q=1,nd
                ! Absorption Term (N*grad(u),grad(v)):
                ! -----------------------------------

                STIFF (p,q) = STIFF(p,q) + Weight * &
                    * SUM(N(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

                ! Polarization Term (K*u,v)
                ! -----------------------------------
                STIFF(p,q) = STIFF(p,q) + Weight * K*Basis(q) * Basis(p)

            END DO
        END DO

        FORCE(1:nd) = 0._dp
        MASS(1:nd) = 0._dp
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
SUBROUTINE LocalMatrixBC( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BC,'field flux', Found )
    Coeff(1:n) = GetReal( BC,'robin coefficient', Found )
    Ext_t(1:n) = GetReal( BC,'external field', Found )

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ

        ! Evaluate terms at the integration point:
        !------------------------------------------

        ! Given flux:
        ! -----------
        F = SUM(Basis(1:n)*flux(1:n))

        ! Robin condition (C*(u-u_0)):
        ! ---------------------------
        C = SUM(Basis(1:n)*coeff(1:n))
        Ext = SUM(Basis(1:n)*ext_t(1:n))

        DO p=1,nd
            DO q=1,nd
                STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
        END DO

        FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
    END DO
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

! Perform static condensation in case bubble dofs are present
!------------------------------------------------------------------------------
SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
        Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SEOPLaser
!------------------------------------------------------------------------------
