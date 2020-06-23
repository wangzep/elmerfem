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
! * Essentially what we are at here is to solve the equation for optical pumping rate
! * presented by Fink et al. The equation is
! *
! * diff(opt,z)-beta*N_Rb*opt*(1-opt/(opt+sigma_SD))
! *
! * where opt is the optical pumping rate, beta is a complicated function that we will define in
! * the helper function file OPUtil, N_Rb is the rubidium number denisty, which we will get from
! * another solver or which will be specified by the user
! * and sigma_SD is the rubidium spin destruction rate which is either specified by the user or
! * calculated as a function of the temperature solver.
! *
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE OPSolver( Model,Solver,dt,TransientSimulation )
    !------------------------------------------------------------------------------
    USE DefUtils
    !USE OPUtil

    ! * Boiler-plate begining code for Elmer Solvers. See Elmer Solver manaul for details.
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    REAL(KIND=dp) :: Norm
    INTEGER :: n, nb, nd, t, active
    INTEGER :: iter, maxiter, newton_interations
    LOGICAL :: found, newton = .FALSE.
    !------------------------------------------------------------------------------

    ! * Again, this is boiler plate code for Elmer Solvers

    CALL DefaultStart()

    maxiter = ListGetInteger( GetSolverParams(),&
        'Nonlinear System Max Iterations',found,minv=2)
    IF(.NOT. found ) maxiter = 2

    ! Nonlinear iteration loop:
    !--------------------------
    DO iter=1,maxiter

        ! System assembly:
        !----------------
        CALL DefaultInitialize()
        ! * Number of active cells
        Active = GetNOFActive()

        ! * Loop over the number of active elements to make the local matrix.
        DO t=1,Active
            Element => GetActiveElement(t)
            n  = GetElementNOFNodes()
            nd = GetElementNOFDOFs()
            nb = GetElementNOFBDOFs()
            IF (newton) THEN
                CALL Warn('OPSolver',&
                    "Newton's Method has not yet been implemented in this solver. Continuing with Picard Interations.")

                !Assemble the matrix
                CALL LocalMatrix(Element, n, nd)
            ELSE

                ! Assemble the matrix
                CALL LocalMatrix(  Element, n, nd+nb )
            END IF

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
        REAL(KIND=dp) ::  laser_direction(3,n),directionterm(3), Weight, nonlinearterm
        REAL(KIND=dp) :: alkali_density(n), spin_destruction_rate(n), alkali_line_width(n), &
            beta(n),previous_solution(n), temp(n)

        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ!,LoadAtIP
        REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)

        LOGICAL :: Stat,found

        INTEGER :: i,t,p,q,dim
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(ValueList_t), POINTER :: BodyForce, Material, Constants
        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------

        dim = CoordinateSystemDimension()

        CALL GetElementNodes( Nodes )
        MASS  = 0._dp
        STIFF = 0._dp
        FORCE = 0._dp
        LOAD = 0._dp

        !BodyForce => GetBodyForce()
        !IF ( ASSOCIATED(BodyForce) ) &
        !    Load(1:n) = GetReal( BodyForce,'field source', found )

        !Get the previous solution
        CALL GetScalarLocalSolution(previous_solution, 'optrate', Element)

        !Pointers to the SIF Constants and Material information
        Constants => GetConstants()
        Material => GetMaterial()

        !Let's get the relevent information
        alkali_density(1:n)=GetReal(Material,'alkali density',found)
        spin_destruction_rate(1:n)=GetReal(Material,'spin destrution rate',found)
        alkali_line_width(1:n)=GetReal(Constants,'alkali line width',found)
        beta(1:n)=GetReal(Material,'beta',found)

        laser_direction = 0._dp
        DO i=1,dim
            laser_direction(i,1:n)=GetReal(Material,&
                'laser direction '//TRIM(I2S(i)),found)
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

            directionterm = MATMUL(laser_direction(:,1:n),Basis(1:n))


            !Put together the nonlinear term
            temp=spin_destruction_rate+previous_solution
            temp=previous_solution/temp
            temp=1-temp
            temp=alkali_line_width*temp
            temp=alkali_density*temp
            temp=beta*temp

            nonlinearterm = SUM(Basis(1:n)*temp(1:n))


            Weight = IP % s(t) * DetJ

            ! diffusion term (D*grad(u),grad(v)):
            ! -----------------------------------
            !STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
            !    D * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

            DO p=1,nd
                DO q=1,nd
                    ! advection term (C*grad(u),v)
                    ! -----------------------------------
                    STIFF (p,q) = STIFF(p,q) + Weight * &
                        SUM(directionterm(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

                    ! reaction term (R*u,v)
                    ! -----------------------------------
                    STIFF(p,q) = STIFF(p,q) + Weight * nonlinearterm * Basis(q) * Basis(p)

                  ! Mass matrix is identically zero in order to solve the equation instaniously at each time point, that is
                  !the equation has no knowledge of its history.
                  ! ------------------------------
                  !MASS(p,q) = MASS(p,q) + Weight * rho * Basis(q) * Basis(p)
                END DO
            END DO

            !FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
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
        LOGICAL :: Stat,found
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

        Flux(1:n)  = GetReal( BC,'field flux', found )
        Coeff(1:n) = GetReal( BC,'robin coefficient', found )
        Ext_t(1:n) = GetReal( BC,'external field', found )

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

!------------------------------------------------------------------------------
END SUBROUTINE OPSolver
!------------------------------------------------------------------------------
