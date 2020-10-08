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
!/*****************************************************************************/
! *
! * A prototype solver for advection-diffusion-reaction equation,
! * This equation is generic and intended for education purposes
! * but may also serve as a starting point for more complex solvers.
! *
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE SESolver( Model,Solver,dt,TransientSimulation )
    !------------------------------------------------------------------------------
    USE DefUtils
    USE Adaptive

    IMPLICIT NONE

    INTERFACE
        FUNCTION InsideResidual(Model, Element, Mesh,&
            Quant, Perm, Fnorm) RESULT( Indicator)
            !------------------------------------------------------------------------------

            USE Types
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: Perm(:)
            REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
            TYPE(Mesh_t), POINTER    :: Mesh
            TYPE(Element_t), POINTER :: Element
        END FUNCTION InsideResidual

        FUNCTION EdgeResidual( Model, Edge, Mesh, Quant, Perm ) RESULT( Indicator )
            !------------------------------------------------------------------------------
            USE Types

            IMPLICIT NONE

            TYPE(Model_t) :: Model
            INTEGER :: Perm(:)
            REAL(KIND=dp) :: Quant(:), Indicator(2)
            TYPE( Mesh_t ), POINTER    :: Mesh
            TYPE( Element_t ), POINTER :: Edge
        END FUNCTION EdgeResidual

        FUNCTION BoundaryResidual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
            !------------------------------------------------------------------------------
            USE Types

            IMPLICIT NONE
            !------------------------------------------------------------------------------
            TYPE(Model_t) :: Model
            INTEGER :: Perm(:)
            REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
            TYPE( Mesh_t ), POINTER    :: Mesh
            TYPE( Element_t ), POINTER :: Edge
        END FUNCTION BoundaryResidual

    END INTERFACE

    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    REAL(KIND=dp), POINTER :: XePol(:)
    INTEGER, POINTER :: Perms(:)
    REAL(KIND=dp) :: Norm
    INTEGER :: n, nb, nd, t, active
    INTEGER :: iter, maxiter
    LOGICAL :: Found
    !------------------------------------------------------------------------------

    CALL DefaultStart()
    !Check all arguments related to SESolver in the SIF
    CALL ArguementCheck()

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
                CALL LocalMatrixBC(  Element, n, nd )
            END IF
        END DO

        CALL DefaultFinishBoundaryAssembly()
        CALL DefaultFinishAssembly()
        CALL DefaultDirichletBCs()

        ! And finally, solve:
        !--------------------
        Norm = DefaultSolve()

        XePol=> Solver % Variable % Values
        Perms=> Solver % Variable % Perm

        IF (ListGetLogical(Solver % Values, 'Adaptive Mesh Refinement', Found)) &
            CALL RefineMesh(Model,Solver,Xepol,Perms, &
            InsideResidual, EdgeResidual, BoundaryResidual)

        IF( DefaultConverged() ) EXIT

    END DO

    CALL DefaultFinish()
  
CONTAINS

    ! Assembly of the matrix entries arising from the bulk elements
    !------------------------------------------------------------------------------
    SUBROUTINE LocalMatrix( Element, n, nd )
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: diff_coeff(n), conv_coeff(n),react_coeff(n), &
            time_coeff(n), D,C,R, rho,Velo(3,n),a(3), Weight
        REAL(KIND=dp) :: spinexchangerate(n), spinrelaxationrate(n),&
            alkalipolarization(n)
        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
        REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
        LOGICAL :: Stat,found
        INTEGER :: ind, i,t,p,q,dim
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(ValueList_t), POINTER :: BodyForce, Material
        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------

        dim = CoordinateSystemDimension()

        CALL GetElementNodes( Nodes )
        MASS  = 0._dp
        STIFF = 0._dp
        FORCE = 0._dp
        LOAD = 0._dp

        !Diffusion coefficient
        Material => GetMaterial()
        diff_coeff(1:n)=GetReal(Material,'xe diffusivity',Found)

        DO ind = 1, n
            IF (diff_coeff(ind) .LT. 0) THEN
                CALL FATAL('SESolver',&
                    'diffusion coefficient is less than 0, this is not physically possible')
            END IF
        END DO

        alkalipolarization(1:n)=GetReal(Material, 'alkali polarization', found)
        CALL FoundCheck(found, 'alkali polarization', 'fatal')

        DO ind = 1, n
            IF (alkalipolarization(ind) .LT. 0 .or. alkalipolarization(ind) .GT. 1) THEN
                CALL FATAL('SESolver',&
                    'alkali polarization is not bound by 0 and 1')
            END IF
        END DO

        spinrelaxationrate(1:n)=GetReal(Material, 'spin relaxation rate', found)
        CALL FoundCheck(found, 'spin relaxation rate', 'fatal')

        DO ind = 1, n
            IF (spinrelaxationrate(ind) .LT. 0) THEN
                CALL FATAL('SESolver',&
                    'spin relaxation rate is less than 0, this is not physically possible')
            END IF
        END DO

        spinexchangerate(1:n)=GetReal(Material, 'spin exchange rate', found)
        CALL FoundCheck(found, 'spin exchange rate', 'fatal')

        DO ind = 1, n
            IF (spinexchangerate(ind) .LT. 0) THEN
                CALL FATAL('SESolver',&
                    'spin exchange rate is less than 0, this is not physically possible')
            END IF
        END DO

        Load(1:n) = spinexchangerate*alkalipolarization


        !The reaction coefficient is the spin-exchange rate +
        !The spin-destruction rate.

        react_coeff = spinexchangerate + spinrelaxationrate

        time_coeff(1:n) = 1
        Velo = 0._dp
        DO i=1,dim
            Velo(i,1:n)=GetReal(Material,&
                'convection velocity '//TRIM(I2S(i)),found)
        END DO

        ! Numerical integration:
        !-----------------------
        IP = GaussPoints( Element )
        DO t=1,IP % n
            ! Basis function values & derivatives at the integration point:
            !--------------------------------------------------------------
            stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis, dBasisdx )

            ! The source term at the integration point:
            !------------------------------------------
            LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

            rho = SUM(Basis(1:n)*time_coeff(1:n))
            a = MATMUL(Velo(:,1:n),Basis(1:n))
            D = SUM(Basis(1:n)*diff_coeff(1:n))
            C = SUM(Basis(1:n))
            R = SUM(Basis(1:n)*react_coeff(1:n))

            Weight = IP % s(t) * DetJ

            ! diffusion term (D*grad(u),grad(v)):
            ! -----------------------------------
            STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
                D * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

            DO p=1,nd
                DO q=1,nd
                    ! advection term (C*grad(u),v)
                    ! -----------------------------------
                    STIFF (p,q) = STIFF(p,q) + Weight * &
                        C * SUM(a(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

                    ! reaction term (R*u,v)
                    ! -----------------------------------
                    STIFF(p,q) = STIFF(p,q) + Weight * R*Basis(q) * Basis(p)

                    ! time derivative (rho*du/dt,v):
                    ! ------------------------------
                    MASS(p,q) = MASS(p,q) + Weight * rho * Basis(q) * Basis(p)
                END DO
            END DO

            FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
        END DO

        IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
        CALL CondensateP( nd-nb, nb, STIFF, FORCE )
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
        INTEGER :: ind, i,t,p,q,dim
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

        Flux = 0.0d0
        Coeff(1:n) = 0.0d0
        Ext_t(1:n) = 0.0d0

        Coeff(1:n)  = GetReal( BC,'T1 Coefficient', Found )
        !CALL FoundCheck(Found, 'T1 Coefficient', 'warn')

        IF (.NOT. Found) THEN
            CALL INFO('SESolver', 'T1 Coefficient is not found', level = 15)
        END IF

        DO ind = 1, n
            IF (Coeff(ind) .LT. 0) THEN
                CALL FATAL('SESolver',&
                    'T1 coefficient is less than 0, this is not physically possible')
            END IF
        END DO

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
            F = SUM(Basis(1:n)*Flux(1:n))

                  ! Robin condition (C*(u-u_0)):
                  ! ---------------------------
            C = SUM(Basis(1:n)*Coeff(1:n))
            Ext = SUM(Basis(1:n)*Ext_t(1:n))

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

    !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE SESolver
!------------------------------------------------------------------------------
SUBROUTINE ArguementCheck()
    USE defutils
    IMPLICIT NONE

    REAL(KIND=dp) :: cellradius, T1, T1pressure, T1temperature, binaryexrate,&
        xemolserate, n2molserate, hemolserate, binaryrelrate, xevdWrelrate, herelax,&
        n2relax
    TYPE(ValueList_t), POINTER :: Constants
    LOGICAL :: found
    !---------------------------------------------------------------

    Constants => GetConstants()

    cellradius = GetConstReal(Constants, 'cell radius', found)
    IF (found) THEN
        IF (cellradius .lt. 0) THEN
            CALL Fatal('SESolver', 'Cell radius is less than 0')
        END IF
    END IF

    T1 = GetConstReal(Constants, 'T1', found)
    IF (found) THEN
        IF (T1 .lt. 0) THEN
            CALL Fatal('SESolver', 'T1 is less than 0')
        END IF
    END IF

    T1pressure = GetConstReal(Constants, 'T1 pressure', found)
    IF (found) THEN
        IF (T1pressure .lt. 0) THEN
            CALL Fatal('SESolver', 'T1 pressure is less than 0')
        END IF
    END IF

    T1temperature = GetConstReal(Constants, 'T1 temperature', found)
    IF (found) THEN
        IF (T1temperature .lt. 0) THEN
            CALL Fatal('SESolver', 'T1 temperature is less than 0')
        END IF
    END IF

END SUBROUTINE ArguementCheck


SUBROUTINE FoundCheck(found,name,warn_fatal_flag)
    !------------------------------------------------------------------------------
    USE DefUtils

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: found
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: warn_fatal_flag
    CHARACTER(len=len(name)+28) :: outputstring

    !Putting together the text to be printed with the warning or fatal warning.

    outputstring=TRIM('The parameter '//name//' was not found')


    IF (.NOT. found) THEN
        IF (warn_fatal_flag .EQ. 'warn') THEN
            CALL Warn('SESolver', outputstring)
        ELSE
            CALL Fatal('SESolver', outputstring)
        END IF
    END IF

!-------------------------------------------------------------------------
END SUBROUTINE FoundCheck

FUNCTION InsideResidual(Model, Element, Mesh,&
    Quant, Perm, Fnorm) RESULT( Indicator )
    !------------------------------------------------------------------------------

    USE Defutils

    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
    TYPE(Mesh_t), POINTER    :: Mesh
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: f, hK, detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:),&
        spinexchangerate(:), alkalipolarization(:), solution(:),&
        diff_coeff(:), spinrelaxationrate(:), Velo(:,:)
    LOGICAL :: stat, found
    INTEGER :: n, dim, j, numnodes

    numnodes = Element % TYPE % NumberOfNodes

    dim=CoordinateSystemDimension()

    ALLOCATE(Basis(numnodes), dBasisdx(numnodes,dim),ddBasisddx(numnodes,dim,dim),&
        spinexchangerate(numnodes), alkalipolarization(numnodes),&
        solution(numnodes), diff_coeff(numnodes), spinrelaxationrate(numnodes),&
        Velo(dim,numnodes))

    CALL GetElementNodes( Nodes )

    Indicator = 0.0d0
    Fnorm = 0.0d0
    hK = Element % hK

            !Diffusion coefficient
    Material => GetMaterial()
    diff_coeff(1:numnodes)=GetReal(Material,'xe diffusivity',found)

    alkalipolarization(1:numnodes)=GetReal(Material, 'alkali polarization', found)
    CALL FoundCheck(found, 'alkali polarization', 'fatal')

    spinrelaxationrate(1:numnodes)=GetReal(Material, 'spin relaxation rate', found)
    CALL FoundCheck(found, 'spin relaxation rate', 'fatal')

    spinexchangerate(1:numnodes)=GetReal(Material, 'spin exchange rate', found)
    CALL FoundCheck(found, 'spin exchange rate', 'fatal')

    !The reaction coefficient is the spin-exchange rate +
    !The spin-destruction rate.

    Velo = 0._dp
    DO j=1,dim
        Velo(j,1:numnodes)=GetReal(Material,&
            'convection velocity '//TRIM(I2S(j)),found)
    END DO

    CALL GetScalarLocalSolution(solution, 'XePol', Element)

    IP = GaussPoints(Element)

    DO n = 1, IP%n

        stat = ElementInfo(Element, Nodes, IP % u(n), IP % v(n),&
            IP % w(n), detJ, Basis, dBasisdx, ddBasisddx, .TRUE.)



        DO j=1,dim
            !+dDXe.dsolution

            Indicator= SUM(solution(1:numnodes)*dBasisdx(1:numnodes,j))*&
                SUM(diff_coeff(1:numnodes)*dBasisdx(1:numnodes,j))
            !+u.dsolution
            Indicator=Indicator+&
                SUM( Velo(j,1:numnodes) * Basis(1:numnodes) ) * &
                SUM( solution(1:numnodes) * dBasisdx(1:numnodes,j) )
        END DO

        !+(spinex+spinrelax)*solution

        Indicator=Indicator+&
            (SUM(spinexchangerate(1:numnodes)*Basis(1:numnodes))+&
            SUM(spinrelaxationrate(1:numnodes)*Basis(1:numnodes)))*&
            SUM(solution(1:numnodes)*Basis(1:numnodes))

        !-spinex*alkalipol
        Indicator=Indicator-&
            SUM(spinexchangerate(1:numnodes)*Basis(numnodes))*&
            SUM(alkalipolarization(1:numnodes)*Basis(numnodes))

        Indicator = Indicator**2.0D0*detJ*IP % s(n)

    END DO

    Fnorm = SQRT(Fnorm)
    Indicator = hK*SQRT(Indicator)

    DEALLOCATE(Basis, dBasisdx,ddBasisddx,&
        spinexchangerate, alkalipolarization,&
        solution, diff_coeff, spinrelaxationrate,&
        Velo)

!----------------------------------------------------------------------------------
END FUNCTION InsideResidual
!----------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
FUNCTION EdgeResidual( Model, Edge, Mesh, Quant, Perm ) RESULT( Indicator )
    !------------------------------------------------------------------------------
    USE Defutils

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2)
    TYPE( Mesh_t ), POINTER    :: Mesh
    TYPE( Element_t ), POINTER :: Edge

    Indicator=0.0d0

!---------------------------------------------------------------------------------
END FUNCTION EdgeResidual
    !---------------------------------------------------------------------------------

FUNCTION BoundaryResidual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
    !------------------------------------------------------------------------------
    USE DefUtils

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
    TYPE( Mesh_t ), POINTER    :: Mesh
    TYPE( Element_t ), POINTER :: Edge

    Indicator=0.0d0
!---------------------------------------------------------------------------------
END FUNCTION BoundaryResidual
    !--------------------------------------------------------------------------------


