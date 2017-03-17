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
! *****************************************************************************
!
!******************************************************************************
! *
! *  Authors: Mikko, Lyly, Juha Ruokolainen, Thomas Zwinger
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 04 Apr 2004
! *  Updated by: Geoffry Schrank
! *  Update Date: 06 Feb 2017
! *
! ****************************************************************************

!------------------------------------------------------------------------------
!> A prototype solver that heavily relies on the Advection-Reaction Solver code
!> to solver for laser absorption in SEOP processes.
!>  Advection-reaction equation solver for scalar fields with discontinous Galerkin method.
!> \ingroup Solvers. Stable Verision.
!------------------------------------------------------------------------------
SUBROUTINE SEOPSolver( Model,Solver,dt,Transient )
    !------------------------------------------------------------------------------
    USE DefUtils

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
    TYPE(Solver_t), TARGET :: Solver  !< Linear & nonlinear equation solver options
    REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
    LOGICAL :: Transient    !< Steady state or transient simulation
    !------------------------------------------------------------------------------
    !    Local variables
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC, BodyForce, Material,&
        Equation, Constants, SolverParams

    TYPE(Element_t), POINTER :: Element, Face, &
        ParentElement, LeftParent, RightParent

    LOGICAL :: AllocationsDone = .FALSE., Found, Stat, &
        LimitSolution, FoundLowerLimit, FoundUpperLimit
    INTEGER :: Active, DIM,NonLinearIterMin,NonlinearIterMax,iter,&
        CorrectedLowerLimit,CorrectedUpperLimit
    INTEGER :: n1,n2, k, n, t, istat, i, j, dummyInt, NumberOfFAces, Indexes(128)
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: Norm,RelativeChange,at,at0,totat,st,totst,&
          OriginalValue
#else
        REAL(KIND=dp) :: Norm,RelativeChange,at,at0,totat,st,totst,CPUTime,RealTime,&
        OriginalValue
#endif
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), &
        FORCE(:), Velo(:,:), MeshVelo(:,:), Gamma(:), Ref(:), &
        UpperLimit(:), LowerLimit(:), nRb(:), spin_destruction(:), RbPol_Term(:), &
        Absorption_Term(:), Flux(:)

    !Laser Configuration Variables Declaration------------------------------------------
    REAL(KIND=dp) ::rubidium_wavelength,rubidium_freq_width,laser_wavelength,&
        laser_linewidth,oscillator_strength,Beta
    !---------------------------------------

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER ::  Var
    CHARACTER(LEN=MAX_NAME_LEN) ::  VariableName, SolverName, ExpVariableName

    SAVE MASS, STIFF, LOAD, FORCE, Velo, MeshVelo, Gamma, &
        AllocationsDone, DIM, VariableName, SolverName, &
        UpperLimit, LowerLimit, nRb, spin_destruction, RbPol_Term, Absorption_Term, Flux
    !*******************************************************************************

    TYPE( Element_t ), POINTER :: Faces(:)
    TYPE(Nodes_t) :: ElementNodes

    INTEGER :: DOFs
    !*******************************************************************************
    SolverName = 'SEOPSolver ('// TRIM(Solver % Variable % Name) // ')'
    VariableName = TRIM(Solver % Variable % Name)
    WRITE(Message,'(A,A)')&
        'AdvectionReactionSolver for variable ', VariableName
    CALL INFO(SolverName,Message,Level=1)

    Mesh => GetMesh()

    IF ( CoordinateSystemDimension() == 2 ) THEN
        Faces => Mesh % Edges
        NumberOfFaces = Mesh % NumberOfEdges
    ELSE
        Faces => Mesh % Faces
        NumberOfFaces = Mesh % NumberOfFaces
    END IF

    ! Initialize & allocate some permanent storage, this is done first time only:
    !----------------------------------------------------------------------------
    IF ( .NOT. AllocationsDone ) THEN

        N = 2 * MAX(Mesh % MaxElementDOFs, Mesh % MaxElementNodes )
        ALLOCATE( FORCE(N), MASS(n,n), STIFF(N,N), LOAD(N),  &
            Velo(3,N), MeshVelo( 3,N ), Gamma(n), &
            UpperLimit(n), LowerLimit(n), nRb(n), spin_destruction(n), RbPol_Term(n), &
            Absorption_Term(n), Flux(n), STAT = istat )

        IF ( istat /= 0 ) THEN
            CALL FATAL(SolverName,'Memory allocation error.' )
        ELSE
            CALL INFO(SolverName,'Memory allocation done',Level=1 )
        END IF
        DIM = CoordinateSystemDimension()
        AllocationsDone = .TRUE.
    END IF

    !------------------------------------------------------------------------------
    !    Read physical and numerical constants and initialize
    !------------------------------------------------------------------------------
    Constants => GetConstants()
    SolverParams => GetSolverParams()

    !------------------Get Laser Parameters----------------------------------------
    rubidium_wavelength = GetConstReal(Model % Constants,'rubidium wavelength',Found)
    laser_wavelength = GetConstReal(Model % Constants,'laser wavelength',Found)
    laser_linewidth = GetConstReal(Model % Constants,'laser line width',Found)
    rubidium_freq_width = GetConstReal(Model % Constants,'rubidium frequency width',Found)
    oscillator_strength = GetConstReal(Model % Constants,'oscillator strength', Found)
    !------------------------------------------------------------------------------
    NonlinearIterMax = GetInteger(   SolverParams, &
        'Nonlinear System Max Iterations', Found )
    IF ( .NOT.Found ) THEN
        CALL WARN(SolverName,'No > Nonlinear System Max Iterations < found. Setting 1')
        NonlinearIterMax = 1
    END IF

    NonlinearIterMin = GetInteger(   SolverParams, &
        'Nonlinear System Min Iterations', Found )
    IF ( .NOT.Found ) THEN
        CALL WARN(SolverName,'No >Nonlinear System Min Iterations< found. Setting 1')
        NonlinearIterMin = 1
    ELSE IF (NonlinearIterMin > NonlinearIterMax) THEN
        CALL WARN(SolverName,&
            '>Nonlinear System Min Iterations< is exceeding >Nonlinear System Max Iterations<.')
        CALL WARN(SolverName,&
            'First is being reset to the latter.')
        NonlinearIterMin = NonlinearIterMax
    END IF

    ExpVariableName = GetString(SolverParams , 'Exported Variable 1', Found )
    IF (.NOT.Found) &
        CALL FATAL(SolverName,'No value > Exported Variable 1 < found in Solver')

    LimitSolution = GetLogical( SolverParams, &
        'Limit Solution', Found )
    IF ( .NOT.Found ) &
        LimitSolution = .FALSE.

    IF (LimitSolution) THEN
        CALL INFO(SolverName, 'Keyword > Limit Solution < found. Solution will be limited',Level=1)
    ELSE
        CALL INFO(SolverName, 'No keyword > Limit Solution < found. Solution will not be limited',Level=1)
    END IF
    !Calculate Beta for the this laser configuration----------------

    Beta = BetaCalc(rubidium_wavelength,rubidium_freq_width,laser_wavelength,&
        laser_linewidth,oscillator_strength)

    !----------------------------------------------------------------

    !------------------------------------------------------------------------------
    !       non-linear system iteration loop
    !------------------------------------------------------------------------------
    totat = 0; totst = 0
    DO iter=1,NonlinearIterMax
        ! Assembly of the bulk elements:
        !-------------------------------
        !------------------------------------------------------------------------------
        ! print out some information
        !------------------------------------------------------------------------------
        at  = CPUTime()
        at0 = RealTime()

        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, '-------------------------------------',Level=4 )
        WRITE( Message,'(A,I3,A,I3)') &
            'Nonlinear iteration no.', iter,' of max',NonlinearIterMax
        CALL Info( SolverName, Message, Level=4 )
        CALL Info( SolverName, '-------------------------------------',Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, 'Starting Assembly...', Level=4 )


        CALL DefaultInitialize()
        Active = GetNOFActive()
        DO t = 1, Active
            !------------------------------------------------------------------------------
            ! write some info on status of assembly
            !------------------------------------------------------------------------------
            IF ( RealTime() - at0 > 1.0 ) THEN
                WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                    (Active-t) / &
                    (1.0*Active)), ' % done'

                CALL Info( SolverName, Message, Level=5 )

                at0 = RealTime()
            END IF
            !------------------------------------------------------------------------------
            ! assign pointers and get number of nodes in element
            !------------------------------------------------------------------------------
            Element => GetActiveElement( t )
            n = GetElementNOfNodes( Element )
            Material => GetMaterial()
            BodyForce => GetBodyForce( Element )
            Equation => GetEquation()
            IF (LimitSolution) THEN
                dummyInt = GetElementDOFs( Indexes )
                UpperLimit(1:n) = GetReal(Material,TRIM(VariableName) // ' Upper Limit', FoundUpperLimit)
                LowerLimit(1:n) = GetReal(Material,TRIM(VariableName) // ' Lower Limit', FoundLowerLimit)
                DO i=1,n
                    IF (FoundUpperLimit) THEN
                        OriginalValue = Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
                        IF (OriginalValue > UpperLimit(i)) THEN
                            Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) ) = UpperLimit(i)
                            CorrectedUpperLimit = CorrectedUpperLimit + 1
                        END IF
                        IF (OriginalValue < LowerLimit(i)) THEN
                            Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) ) = LowerLimit(i)
                            CorrectedLowerLimit = CorrectedLowerLimit + 1
                        END IF
                    END IF
                END DO
            END IF
            !------------------------------------------------------------------------------
            ! the body force (r.h.s) = source
            !------------------------------------------------------------------------------
            LOAD(1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Source', Found )
            IF (.NOT.Found) THEN
                !WRITE(Message,'(A,A,A)') 'Body Force >',TRIM(VariableName) // ' Source','< not found'
                !CALL INFO(SolverName, Message, Level=42)
                LOAD(1:n)  = 0.0d0
            END IF
            !------------------------------------------------------------------------------
            ! Get convection and mesh velocity
            CALL GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
                Equation,Solver,Model,Element)

            !------------------------------------------------------------------------------
            !Get the Material Properties---------------------------------------------------

            nRb(1:n)=GetReal(Material,'rubidium number density',Found)
            spin_destruction(1:n) = GetReal(Material,'spin destruction rate',Found)

            Absorption_Term(1:n) = Beta*nRb(1:n)

            !-----------------------Get the Local Solution for the Non-Linear Term----------
            !CALL GetScalarLocalSolution(Flux)
            !Use Advect-React Solver's method to get the local solution, as was done to limit
            !Variable Values.

            DO i = 1,n
                Flux(i) = Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
            END DO

            !---------------------Non-linear Term---------------------------------------------------------
            DO i = 1,n
                RbPol_Term(i) = (1.)-(Flux(i)/(Flux(i)+spin_destruction(i)))
                Absorption_Term(i) = Absorption_Term(i)*RbPol_Term(i)
            END DO

            !-----------------------
            ! get reaction constant
            !-----------------------
            Gamma(1:n) = Absorption_Term(1:n)

            !Gamma(1:n)  = GetReal( Material, TRIM(VariableName) // ' Gamma', Found )
            !IF (.NOT.Found) THEN
            !    WRITE(Message,'(A,A,A)') 'Material Property >',TRIM(VariableName) // ' Gamma','< not found'
            !    CALL INFO(SolverName, Message, Level=42)
            !    Gamma(1:n)  = 0.0d0
            !END IF

            CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, Velo, MeshVelo, Gamma, Element, n )
            IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
            CALL DefaultUpdateEquations( STIFF, FORCE )
        END DO


        ! Assembly of the face terms:
        !----------------------------
        FORCE = 0.0d0
        DO t=1,NumberOfFaces
            Face => Faces(t)
            IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

            LeftParent  => Face % BoundaryInfo % Left
            RightParent => Face % BoundaryInfo % Right
            IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
                n  = GetElementNOFNodes( Face )
                n1 = GetElementNOFNodes( LeftParent )
                n2 = GetElementNOFNodes( RightParent )

                !------------------------------------------------------------------------------
                ! Get convection and mesh velocity
                !------------------------------------------------------------------------------
                Material => GetMaterial( LeftParent )
                Equation => GetEquation( LeftParent )
                CALL GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
                    Equation,Solver,Model,Face)

                CALL LocalJumps( STIFF,Face,n,LeftParent,n1,RightParent,n2,Velo )
                CALL DefaultUpdateEquations( STIFF, FORCE, Face )
            END IF
        END DO

        CALL DefaultFinishBulkAssembly()


        ! Loop over the boundary elements:
        !---------------------------------
        DO t=1,Mesh % NumberOfBoundaryElements
            Element => GetBoundaryElement(t)
            IF( .NOT. ActiveBoundaryElement() )  CYCLE

            !------------------------------------------------------------------------------
            ! Check that the dimension of element is suitable for BCs
            !------------------------------------------------------------------------------
            IF( .NOT. PossibleFluxElement(Element) ) CYCLE

            ParentElement => Element % BoundaryInfo % Left
            IF ( .NOT. ASSOCIATED( ParentElement ) ) &
                ParentElement => Element % BoundaryInfo % Right

            n  = GetElementNOFNodes( Element )
            n1 = GetElementNOFnodes( ParentElement )

            !------------------------------------------------------------------------------
            ! Get convection and mesh velocity
            !------------------------------------------------------------------------------
            Material => GetMaterial( ParentElement )
            Equation => GetEquation( ParentElement )
            CALL GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
                Equation,Solver,Model,Element)

            BC => GetBC()
            LOAD = 0.0d0
            Found = .FALSE.
            IF ( ASSOCIATED(BC) ) THEN
                LOAD(1:n) = GetReal( BC, Solver % Variable % Name, Found )
            END IF

            MASS = 0.0d0
            CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD, &
                Element, n, ParentElement, n1, Velo, Found )

            IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
            CALL DefaultUpdateEquations( STIFF, FORCE )
        END DO

        CALL DefaultFinishAssembly()
        CALL Info( SolverName, 'Assembly done', Level=4 )

        !------------------------------------------------------------------------------
        !     Solve the system and check for convergence
        !------------------------------------------------------------------------------
        at = CPUTime() - at
        st = CPUTime()

        ! Solve the system:
        !------------------
        Norm = DefaultSolve()

        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
        CALL Info( SolverName, Message, Level=4 )
        WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
        CALL Info( SolverName, Message, Level=4 )

        RelativeChange = Solver % Variable % NonlinChange

        IF ( Solver % Variable % NonlinConverged == 1 )  THEN
            WRITE(Message,'(A,I6,A,I6,A)') &
                'Nonlinear iteration converged after ', iter, &
                ' out of max ',NonlinearIterMax,' iterations'
            CALL INFO(SolverName,Message)
            EXIT
        ELSE IF ((iter .EQ. NonlinearIterMax) .AND. (NonlinearIterMax > 1)) THEN
            CALL WARN(SolverName,'Maximum nonlinear iterations reached, but system not converged')
        END IF

    !-----------------------------------------------
    END DO ! End of nonlinear iteration loop
    !----------------------------------------------

    !------------------------------------------------------------------------------
    ! limit solution
    !------------------------------------------------------------------------------
    CorrectedUpperLimit = 0
    CorrectedLowerLimit = 0
    IF (LimitSolution) THEN
        Active = GetNOFActive()
        DO t = 1, Active
            Element => GetActiveElement( t )
            n = GetElementNOfNodes( Element )
            Material => GetMaterial()
            dummyInt = GetElementDOFs( Indexes )
            UpperLimit(1:n) = GetReal(Material,TRIM(VariableName) // ' Upper Limit', FoundUpperLimit)
            LowerLimit(1:n) = GetReal(Material,TRIM(VariableName) // ' Lower Limit', FoundLowerLimit)
            DO i=1,n
                IF (FoundUpperLimit) THEN
                    OriginalValue = Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
                    IF (OriginalValue > UpperLimit(i)) THEN
                        Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) ) = &
                            MIN(UpperLimit(i), OriginalValue)
                        CorrectedUpperLimit = CorrectedUpperLimit + 1
                    END IF
                    IF (OriginalValue < LowerLimit(i)) THEN
                        Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) ) = &
                            MAX(LowerLimit(i), OriginalValue )
                        CorrectedLowerLimit = CorrectedLowerLimit + 1
                    END IF
                END IF
            END DO
        END DO
        WRITE(Message,'(a,i10)') 'Limited values for upper limit: ', CorrectedUpperLimit
        CALL Info( SolverName, Message, Level=3 )
        WRITE(Message,'(a,i10)') 'Limited values for lower limit: ', CorrectedLowerLimit
        CALL Info( SolverName, Message, Level=3 )
    END IF


    ! Average the elemental results to nodal values:
    !-----------------------------------------------
    Var => VariableGet( Mesh % Variables,TRIM(ExpVariableName))
    IF ( ASSOCIATED( Var ) ) THEN
        n1 = Mesh % NumberOfNodes
        ALLOCATE( Ref(n1) )
        Ref = 0
        Var % TYPE = Variable_on_nodes

        IF ( ASSOCIATED( Var % Perm, Solver % Variable % Perm ) ) THEN
            ALLOCATE( Var % Perm(SIZE(Solver % Variable % Perm))  )
            Var % Perm = 0
            DO i = 1,n1
                Var % Perm(i) = i
            END DO
        END IF

        Var % Values = 0.0d0
        DO t=1,Active
            Element => GetActiveElement(t)
            n = GetElementDOFs( Indexes )
            n = GetElementNOFNodes()
            DO i=1,n
                k = Element % NodeIndexes(i)
                Var % Values(k) = Var % Values(k) + &
                    Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
                Ref(k) = Ref(k) + 1
            END DO
        END DO

        WHERE( Ref > 0 )
            Var % Values(1:n1) = Var % Values(1:n1) / Ref
        END WHERE
        DEALLOCATE( Ref )
    ELSE
        WRITE(Message,'(A,A,A)') 'Exported Variable >',TRIM(VariableName) //   ' Nodal Result','< not found'
        CALL FATAL(SolverName,Message)
    END IF

CONTAINS

    !------------------------------------------------------------------------------
    SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Velo, Mvelo, Gamma, Element, n)
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), &
            LOAD(:), Velo(:,:), Gamma(:), MVelo(:,:)
        INTEGER :: n
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
        REAL(KIND=dp) :: detJ,U,V,W,S,A,L,cu(3),g,divMVelo
        LOGICAL :: Stat
        INTEGER :: i,p,q,t,dim
        TYPE(GaussIntegrationPoints_t) :: IntegStuff
        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------
        dim = CoordinateSystemDimension()
        FORCE = 0.0d0
        STIFF = 0.0d0
        MASS  = 0.0d0
        CALL GetElementNodes( Nodes, Element )
        !------------------------------------------------------------------------------
        !      Numerical integration
        !------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Element )

        DO t=1,IntegStuff % n
            U = IntegStuff % u(t)
            V = IntegStuff % v(t)
            W = IntegStuff % w(t)
            S = IntegStuff % s(t)
            !------------------------------------------------------------------------------
            !        Basis function values & derivatives at the integration point
            !------------------------------------------------------------------------------
            stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                Basis, dBasisdx )

            S = S * detJ
            L = SUM( LOAD(1:n) *  Basis(1:n) )
            g = SUM( Basis(1:n) * Gamma(1:n) )

            ! This term was missing for the ALE formulation
            divMVelo = 0.0_dp
            DO i=1,dim
                divMVelo = divMVelo + SUM( MVelo(i,1:n) * dBasisdx(1:n,i) )
            END DO

            cu   = 0.0d0
            DO i=1,dim
                cu(i) = SUM( Basis(1:n) * Velo(i,1:n) )
            END DO
            !------------------------------------------------------------------------------
            !        The advection-reaction equation: dc/dt + grad(u . c) + gamma c = s
            !------------------------------------------------------------------------------
            DO p=1,n
                DO q=1,n
                    MASS(p,q)  = MASS(p,q)  + s * Basis(q) * Basis(p)
                    STIFF(p,q) = STIFF(p,q) + s * (g + divMVelo) * Basis(q) * Basis(p)
                    DO i=1,dim
                        STIFF(p,q) = STIFF(p,q) - s * cu(i) * Basis(q) * dBasisdx(p,i)
                    END DO
                END DO
            END DO
            FORCE(1:n) = FORCE(1:n) + s*L*Basis(1:n)
        !------------------------------------------------------------------------------
        END DO
    !------------------------------------------------------------------------------
    END SUBROUTINE LocalMatrix
    !------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW(Face, nFace, Parent, nParent, U, V, W, Basis)
        !------------------------------------------------------------------------------
        IMPLICIT NONE
        TYPE(Element_t) :: Face, Parent
        INTEGER :: nFace, nParent
        REAL( KIND=dp ) :: U, V, W, Basis(:)
        !------------------------------------------------------------------------------
        INTEGER :: i,j
        REAL(KIND=dp) :: ParentU(nFace), ParentV(nFace), ParentW(nFace)
        !------------------------------------------------------------------------------
        DO i = 1,nFace
            DO j = 1,nParent
                IF ( Face % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
                    ParentU(i) = Parent % TYPE % NodeU(j)
                    ParentV(i) = Parent % TYPE % NodeV(j)
                    ParentW(i) = Parent % TYPE % NodeW(j)
                    EXIT
                END IF
            END DO
        END DO
        U = SUM( Basis(1:nFace) * ParentU(1:nFace) )
        V = SUM( Basis(1:nFace) * ParentV(1:nFace) )
        W = SUM( Basis(1:nFace) * ParentW(1:nFace) )
    !------------------------------------------------------------------------------
    END SUBROUTINE FindParentUVW
    !------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF,Face,n,LeftParent,n1,RightParent,n2,Velo )
        !------------------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=dp) :: STIFF(:,:), Velo(:,:)
        INTEGER :: n,n1,n2
        TYPE(Element_t), POINTER :: Face, LeftParent, RightParent
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: FaceBasis(n), FacedBasisdx(n,3)
        REAL(KIND=dp) :: LeftBasis(n1), LeftdBasisdx(n1,3)
        REAL(KIND=dp) :: RightBasis(n2), RightdBasisdx(n2,3)
        REAL(KIND=dp) :: Jump(n1+n2), Average(n1+n2)
        REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
        LOGICAL :: Stat
        INTEGER :: i, j, p, q, dim, t, nFace, nParent
        TYPE(GaussIntegrationPoints_t) :: IntegStuff
        REAL(KIND=dp) :: hE, Normal(3), cu(3), LeftOut(3)

        TYPE(Nodes_t) :: FaceNodes, LeftParentNodes, RightParentNodes
        SAVE FaceNodes, LeftParentNodes, RightParentNodes
        !------------------------------------------------------------------------------
        dim = CoordinateSystemDimension()
        STIFF = 0.0d0

        CALL GetElementNodes( FaceNodes, Face )
        CALL GetElementNodes( LeftParentNodes,  LeftParent )
        CALL GetElementNodes( RightParentNodes, RightParent )
        !------------------------------------------------------------------------------
        !     Numerical integration over the edge
        !------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Face )

        LeftOut(1) = SUM(LeftParentNodes % x(1:n1)) / n1
        LeftOut(2) = SUM(LeftParentNodes % y(1:n1)) / n1
        LeftOut(3) = SUM(LeftParentNodes % z(1:n1)) / n1
        LeftOut(1) = SUM(FaceNodes % x(1:n)) / n - LeftOut(1)
        LeftOut(2) = SUM(FaceNodes % y(1:n)) / n - LeftOut(2)
        LeftOut(3) = SUM(FaceNodes % z(1:n)) / n - LeftOut(3)

        DO t=1,IntegStuff % n
            U = IntegStuff % u(t)
            V = IntegStuff % v(t)
            W = IntegStuff % w(t)
            S = IntegStuff % s(t)

            ! Basis function values & derivatives at the integration point:
            !--------------------------------------------------------------
            stat = ElementInfo( Face, FaceNodes, U, V, W, detJ, &
                FaceBasis, FacedBasisdx )

            S = S * detJ

            Normal = NormalVector( Face, FaceNodes, U, V, .FALSE. )
            IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

            ! Find basis functions for the parent elements:
            ! ---------------------------------------------
            CALL FindParentUVW( Face,n,LeftParent,n1,U,V,W,FaceBasis )
            stat = ElementInfo( LeftParent, LeftParentNodes, U, V, W, detJ, &
                LeftBasis, LeftdBasisdx )

            CALL FindParentUVW( Face,n,RightParent,n2,U,V,W,FaceBasis )
            stat = ElementInfo( RightParent, RightParentNodes, U, V, W, detJ, &
                RightBasis, RightdBasisdx )

            ! Integrate jump terms:
            ! ---------------------
            Jump(1:n1) = LeftBasis(1:n1)
            Jump(n1+1:n1+n2) = -RightBasis(1:n2)

            Average(1:n1) = LeftBasis(1:n1) / 2
            Average(n1+1:n1+n2) = RightBasis(1:n2) / 2

            cu = 0.0d0
            DO i=1,dim
                cu(i) = SUM( Velo(i,1:n) * FaceBasis(1:n) )
            END DO
            Udotn = SUM( Normal * cu )

            DO p=1,n1+n2
                DO q=1,n1+n2
                    STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
                    STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
                END DO
            END DO
        END DO
    !------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
    !------------------------------------------------------------------------------



    !------------------------------------------------------------------------------
    SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, ParentElement, np, Velo, InFlowBC )
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:)
        INTEGER :: n, np
        LOGICAL :: InFlowBC
        TYPE(Element_t), POINTER :: Element, ParentElement
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
        REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3)
        INTEGER :: i,j,p,q,t,dim

        REAL(KIND=dp) :: Normal(3), g, L, Udotn, cu(3), cu1(3), detJ,U,V,W,S
        LOGICAL :: Stat, Inflow
        TYPE(GaussIntegrationPoints_t) :: IntegStuff

        TYPE(Nodes_t) :: Nodes, ParentNodes
        SAVE Nodes, ParentNodes
        !------------------------------------------------------------------------------
        dim = CoordinateSystemDimension()
        FORCE = 0.0d0
        STIFF = 0.0d0

        CALL GetElementNodes( Nodes, Element )
        CALL GetElementNodes( ParentNodes, ParentElement )

        Normal = NormalVector( Element, Nodes, 0.0d0, 0.0d0, .TRUE. )
        DO i=1,3
            cu(i) = SUM( Velo(i,1:n) ) / n
        END DO
        Inflow = InFlowBC .AND. SUM( Normal * cu ) < 0.0d0

        ! Numerical integration:
        !-----------------------
        IntegStuff = GaussPoints( Element )

        DO t=1,IntegStuff % n
            U = IntegStuff % u(t)
            V = IntegStuff % v(t)
            W = IntegStuff % w(t)
            S = IntegStuff % s(t)

            Normal = NormalVector( Element, Nodes, U, V, .TRUE. )

            ! Basis function values & derivatives at the integration point:
            ! -------------------------------------------------------------
            stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                Basis, dBasisdx )
            S = S * detJ

            CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
            stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
                detJ, ParentBasis, ParentdBasisdx )

            L = SUM( LOAD(1:n) * Basis(1:n) )
            cu = 0.0d0
            DO i=1,dim
                cu(i)  = SUM( Velo(i,1:n) * Basis(1:n) )
            END DO
            Udotn = SUM( Normal * cu )

            DO p = 1,np
                IF ( Inflow ) THEN
                    FORCE(p) = FORCE(p) - s * Udotn*L*ParentBasis(p)
                ELSE
                    DO q=1,np
                        STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
                    END DO
                END IF
            END DO
        END DO
    !------------------------------------------------------------------------------
    END SUBROUTINE LocalMatrixBoundary
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE GetLocalALEVelocity(Velo,MeshVelo,SolverName,Material,&
        Equation,Solver,Model,Element)
        IMPLICIT NONE
        !------------------------------------------------------------------------------
        REAL (KIND=dp) :: Velo(:,:), MeshVelo(:,:)
        CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
        TYPE(ValueList_t), POINTER :: Material,Equation
        TYPE(Solver_t), TARGET :: Solver
        TYPE(Model_t), TARGET :: Model
        TYPE(Element_t),POINTER :: Element
        !------------------------------------------------------------------------------
        CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, FlowSolName
        INTEGER :: i,j,k,N,FlowDOFs
        INTEGER, POINTER :: FlowPerm(:)
        REAL(KIND=dp), POINTER :: FlowSolution(:)
        TYPE(Variable_t), POINTER :: FlowSol
        LOGICAL :: Found
        !------------------------------------------------------------------------------
        n  = GetElementNOFNodes( Element )
        ConvectionFlag = GetString( Equation, 'Convection', Found )
        !IF (.NOT. Found) &
        !    CALL FATAL(SolverName, 'No string for keyword >Convection< found in Equation')
        !Velo = 0.0d00
        !Get the laser direction. I've commented out all the computed velocity stuff because
        !I don't need it.
        !---------------------------------------------------
        !IF ( ConvectionFlag == 'constant' ) THEN
        Velo(1,1:N) = GetReal( Material, 'Laser Direction 1', Found, Element )
        IF (.NOT.Found) Velo(1,1:N) = 0.0d0
        Velo(2,1:N) = GetReal( Material, 'Laser Direction 2', Found, Element )
        IF (.NOT.Found) Velo(2,1:N) = 0.0d0
        Velo(3,1:N) = GetReal( Material, 'Laser Direction 3', Found, Element )
        IF (.NOT.Found) Velo(3,1:N) = 0.0d0
           ! computed velocity
           !------------------Equation => GetEquation()
        !ELSE IF (ConvectionFlag == 'computed' ) THEN
        !    FlowSolName =  GetString( Equation,'Flow Solution Name', Found)
        !    IF(.NOT.Found) THEN
        !        CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Equation<')
        !        CALL WARN(SolverName,'Taking default value >Flow Solution<')
        !        WRITE(FlowSolName,'(A)') 'Flow Solution'
        !    END IF
        !    FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
        !    IF ( ASSOCIATED( FlowSol ) ) THEN
        !        FlowPerm     => FlowSol % Perm
        !        FlowDOFs     =  FlowSol % DOFs
        !        FlowSolution => FlowSol % Values
        !    ELSE
        !        WRITE(Message,'(A,A,A)') &
        !            'Convection flag set to >computed<, but no variable >',FlowSolName,'< found'
        !        CALL FATAL(SolverName,Message)
        !    END IF


        !    DO i=1,n
        !        k = FlowPerm(Element % NodeIndexes(i))
        !        IF ( k > 0 ) THEN
                    ! Pressure = FlowSolution(FlowDOFs*k)

        !            SELECT CASE( FlowDOFs )
         !               CASE(2)
          !                  Velo(1,i) = FlowSolution( FlowDOFs*k-1 )
           !                 Velo(2,i) = 0.0d0
            !                Velo(3,i) = 0.0d0
             !           CASE(3)
              !              Velo(1,i) = FlowSolution( FlowDOFs*k-2 )
               !             Velo(2,i) = FlowSolution( FlowDOFs*k-1 )
                !            Velo(3,i) = 0.0d0
                 !       CASE(4)
                  !          Velo(1,i) = FlowSolution( FlowDOFs*k-3 )
                   !         Velo(2,i) = FlowSolution( FlowDOFs*k-2 )
                    !        Velo(3,i) = FlowSolution( FlowDOFs*k-1 )
                    !END SELECT
                !END IF
            !END DO
        !ELSE IF (ConvectionFlag == 'none' ) THEN
        !    Velo = 0.0d0
        !ELSE
        !    WRITE(Message,'(A,A,A)') 'Convection flag >', ConvectionFlag ,'< not recognised'
        !    CALL FATAL(SolverName,Message)
        !END IF
        CALL INFO(SolverName,Message,Level=10)
        !-------------------------------------------------
        ! Add mesh velocity (if any) for ALE formulation
        !-------------------------------------------------
        MeshVelo = 0.0d0
        CALL GetVectorLocalSolution( MeshVelo, 'Mesh Velocity', Element)
        IF (ANY( MeshVelo /= 0.0d0 )) THEN
            DO i=1,FlowDOFs
                Velo(i,1:N) = Velo(i,1:N) - MeshVelo(i,1:N)
            END DO
        END IF
    END SUBROUTINE GetLocalALEVelocity

    FUNCTION BetaCalc(rubidium_wavelength,rubidium_freq_width,laser_wavelength,&
        laser_linewidth,oscillator_strength) RESULT (Beta)
        USE DefUtils
        IMPLICIT NONE

        !------------------------------------------------------------------------!
        REAL(KIND=dp) :: TWO, C, LOG2
        REAL(KIND=dp) :: electron_radius, oscillator_strength, laser_wavelength, &
            laser_linewidth, frequency_shift,&
            absorb_laser_ratio,laser_frequency,&
            rubidium_frequency, laser_freq_width,&
            rubidium_freq_width, rubidium_wavelength,w_prime, w_dprime,&
            w_input_real,w_input_imaginary, Beta
        LOGICAL :: FLAG, Found
        !------------------------------------------------------------------------!
        !Declare constants-------------------------------------------------------
        TWO = 2.0D0
        C = 299792458.0D0
        LOG2 = LOG(TWO)
        electron_radius = 2.8179403267e-15
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------

        !Define spectral overlap function (eq. A7, A10 of Fink's paper)
        rubidium_frequency = C/rubidium_wavelength

        !The second equation here uses the dispresion relationship
        laser_frequency = C/laser_wavelength

        laser_freq_width = (C*laser_linewidth)/(laser_wavelength)**2


        absorb_laser_ratio = rubidium_freq_width/laser_freq_width

        frequency_shift = 2*(laser_frequency - rubidium_frequency)/laser_freq_width

        !So... A8 is the Faddeeva function written in a funny way. Note, that in order to recover
        !The Faddeeva function normally, we need to divide the arguement of erfc in A8 by -i.
        !When you do that you get that the real part is -sqrt(ln(2))*s and the imaginary part
        !is sqrt(ln(2))*r. This is very confusing, and I'm going to make a note in my lab note-
        !book regarding the maninpulation.

        !So I don't think the below assignments are wrong!
        w_input_real = -SQRT(LOG2)*(frequency_shift)
        w_input_imaginary = SQRT(LOG2)*(absorb_laser_ratio)


        CALL WOFZ(w_input_real,w_input_imaginary,w_prime,w_dprime,FLAG)

        Beta = 2*DSQRT(PI*LOG2)*(electron_radius*oscillator_strength*&
            laser_wavelength**2*w_prime)/laser_linewidth

    END FUNCTION

    !We need the FADDEEVA function from equation (A8). Thankfully, someone in the
    !interweb has coded a nice function to calculate it for me. The ACM has collected
    !aligorithms to do this very thing. The following code is copied line for line
    !from ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
    !THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !VOL. 16, NO. 1, PP. 47. - http://www.netlib.org/toms/68

    !      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
    !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !      VOL. 16, NO. 1, PP. 47.
    SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
        !
        !  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
        !  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
        !  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
        !  MEANS SQRT(-1).
        !  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
        !  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
        !  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
        !  OF THE FUNCTION.
        !  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
        !
        !
        !  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
        !     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
        !                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
        !                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
        !                FLOATING-POINT ARITHMETIC
        !     RMAXEXP  = LN(RMAX) - LN(2)
        !     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
        !                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
        !  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
        !  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
        !
        !
        !  PARAMETER LIST
        !     XI     = REAL      PART OF Z
        !     YI     = IMAGINARY PART OF Z
        !     U      = REAL      PART OF W(Z)
        !     V      = IMAGINARY PART OF W(Z)
        !     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
        !              OCCUR OR NOT; TYPE LOGICAL;
        !              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
        !              MEANING :
        !              FLAG=.FALSE. : NO ERROR CONDITION
        !              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
        !                             BECOMES INACTIVE
        !  XI, YI      ARE THE INPUT-PARAMETERS
        !  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
        !
        !  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
        !
        !  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
        !  PUT TO 0 UPON UNDERFLOW;
        !
        !  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
        !  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
        !


        USE DefUtils

        IMPLICIT REAL(KIND=dp) (A-H, O-Z)

        REAL(KIND=dp) ::  KAPN, NP1
        INTEGER :: NU, I, J

        !DOUBLE PRECISION FACTOR,RMAMXREAL,RMAXEXP,RMAXGONI

        LOGICAL A, B, FLAG
        PARAMETER (FACTOR   = 1.12837916709551257388D0,&
            RMAXREAL = 0.5D+154,&
            RMAXEXP  = 708.50306146160D0,&
            RMAXGONI = 3.53711887601422D+15)

        FLAG = .FALSE.

        XABS = DABS(XI)
        YABS = DABS(YI)
        X    = XABS/6.3
        Y    = YABS/4.4

        !
        !     THE FOLLOWING IF-STATEMENT PROTECTS
        !     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
        !
        IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100

        QRHO = X**2 + Y**2

        XABSQ = XABS**2
        XQUAD = XABSQ - YABS**2
        YQUAD = 2*XABS*YABS

        A     = QRHO.LT.0.085264D0

        IF (A) THEN
            !
            !  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
            !  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
            !  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
            !  ACCURACY
            !
            QRHO  = (1-0.85*Y)*DSQRT(QRHO)
            N     = IDNINT(6 + 72*QRHO)
            J     = 2*N+1
            XSUM  = 1.0/J
            YSUM  = 0.0D0
            DO 10 I=N, 1, -1
                J    = J - 2
                XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
                YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
                XSUM = XAUX + 1.0/J
10          CONTINUE
            U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
            V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
            DAUX =  DEXP(-XQUAD)
            U2   =  DAUX*DCOS(YQUAD)
            V2   = -DAUX*DSIN(YQUAD)

            U    = U1*U2 - V1*V2
            V    = U1*V2 + V1*U2

        ELSE
            !
            !  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
            !  CONTINUED FRACTION
            !  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
            !  ACCURACY
            !
            !  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
            !  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
            !  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
            !  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
            !  TO OBTAIN THE REQUIRED ACCURACY
            !  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
            !  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
            !

            IF (QRHO.GT.1.0) THEN
                H    = 0.0D0
                KAPN = 0
                QRHO = DSQRT(QRHO)
                NU   = IDINT(3 + (1442/(26*QRHO+77)))
            ELSE
                QRHO = (1-Y)*DSQRT(1-QRHO)
                H    = 1.88*QRHO
                H2   = 2*H
                KAPN = IDNINT(7  + 34*QRHO)
                NU   = IDNINT(16 + 26*QRHO)
            ENDIF

            B = (H.GT.0.0)

            IF (B) QLAMBDA = H2**KAPN

            RX = 0.0
            RY = 0.0
            SX = 0.0
            SY = 0.0

            DO 11 N=NU, 0, -1
                NP1 = N + 1
                TX  = YABS + H + NP1*RX
                TY  = XABS - NP1*RY
                C   = 0.5/(TX**2 + TY**2)
                RX  = C*TX
                RY  = C*TY
                IF ((B).AND.(N.LE.KAPN)) THEN
                    TX = QLAMBDA + SX
                    SX = RX*TX - RY*SY
                    SY = RY*TX + RX*SY
                    QLAMBDA = QLAMBDA/H2
                ENDIF
11          CONTINUE

            IF (H.EQ.0.0) THEN
                U = FACTOR*RX
                V = FACTOR*RY
            ELSE
                U = FACTOR*SX
                V = FACTOR*SY
            END IF

            IF (YABS.EQ.0.0) U = DEXP(-XABS**2)

        END IF


        !
        !  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
        !

        IF (YI.LT.0.0) THEN

            IF (A) THEN
                U2    = 2*U2
                V2    = 2*V2
            ELSE
                XQUAD =  -XQUAD

                !
                !         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
                !         AGAINST OVERFLOW
                !
                IF ((YQUAD.GT.RMAXGONI).OR.&
                    (XQUAD.GT.RMAXEXP)) GOTO 100

                W1 =  2*DEXP(XQUAD)
                U2  =  W1*DCOS(YQUAD)
                V2  = -W1*DSIN(YQUAD)
            END IF

            U = U2 - U
            V = V2 - V
            IF (XI.GT.0.0) V = -V
        ELSE
            IF (XI.LT.0.0) V = -V
        END IF

        RETURN

100     FLAG = .TRUE.
        RETURN

    END


!------------------------------------------------------------------------------
END SUBROUTINE SEOPSolver
!------------------------------------------------------------------------------
