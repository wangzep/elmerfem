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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found, PiolaVersion, SecondOrder, LorentzConductivity
  TYPE(ValueList_t), POINTER :: SolverParams
 
  LorentzConductivity = ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .or. &
    ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity")
  IF(LorentzConductivity) CALL INFO("WhitneyAVSolver_Init0", "Material is moving: has Lorentz force present")

  SolverParams => GetSolverParams()
  IF ( .NOT.ListCheckPresent(SolverParams, "Element") ) THEN
    PiolaVersion = GetLogical(SolverParams, &
        'Use Piola Transform', Found )   
    SecondOrder = GetLogical(SolverParams, 'Quadratic Approximation', Found)
    IF (PiolaVersion) THEN
      IF ( Transient .OR. LorentzConductivity ) THEN
        IF (SecondOrder) THEN
          CALL ListAddString( SolverParams, &
              "Element", "n:1 e:2 -brick b:6 -prism b:2 -quad_face b:4 -tri_face b:2" )  
        ELSE
          CALL ListAddString( SolverParams, "Element", "n:1 e:1 -brick b:3 -quad_face b:2" )
        END IF
      ELSE
        IF (SecondOrder) THEN
          CALL ListAddString( SolverParams, "Element", &
              "n:0 e:2 -brick b:6 -prism b:2 -quad_face b:4 -tri_face b:2" )  
        ELSE
          CALL ListAddString( SolverParams, "Element", "n:0 e:1 -brick b:3 -quad_face b:2" )
        END IF
      END IF
    ELSE
      IF ( Transient .OR. LorentzConductivity ) THEN
        CALL ListAddString( SolverParams, "Element", "n:1 e:1" )
      ELSE
        CALL ListAddString( SolverParams, "Element", "n:0 e:1" )
      END IF

      IF(GetString(SolverParams,'Linear System Solver')=='block') THEN
        CALL ListAddString( SolverParams, "Element", "n:1 e:1" )
        CALL ListAddLogical( SolverParams, "Optimize Bandwidth", .FALSE.)
      END IF
    END IF
  END IF
  CALL ListAddLogical( SolverParams,'Use Global Mass Matrix',.TRUE.) 
!------------------------------------------------------------------------------
END SUBROUTINE WhitneyAVSolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve vector potential A
! 
!> sigma @A/@t + rot (1/mu) rot A = J^s + curl(M^s) - sigma grad(V^s)
!>   -div(sigma*(@A/@t+grad(V))=0 .
!
!>  using edge elements (Nedelec/Whitney basis of lowest degree)+nodal basis for V.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  USE CircuitUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element, Edge

  REAL(KIND=dp) :: Norm, PrevDT=-1, RelChange
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, BodyParams

  INTEGER :: n,nb,nd,t,istat,i,j,k,l,nNodes,Active,FluxCount=0, &
          NoIterationsMax,NoIterationsMin, nsize

  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), POINTER :: VecPot(:)
  REAL (KIND=DP), POINTER :: Cwrk(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:), Tcoef(:,:,:), &
                                GapLength(:), AirGapMu(:), LamThick(:), &
                                LamCond(:), Wbase(:), RotM(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:), PrevSol(:)

  CHARACTER(LEN=MAX_NAME_LEN):: LaminateStackModel, CoilType
  LOGICAL :: LaminateStack, CoilBody

  INTEGER, POINTER :: Perm(:)
  INTEGER, ALLOCATABLE :: FluxMap(:)
  LOGICAL, ALLOCATABLE :: TreeEdges(:)
  LOGICAL :: Stat, EigenAnalysis, TG, DoneAssembly=.FALSE., &
         SkipAssembly, ConstantSystem, ConstantBulk, FixJ, FoundRelax, &
         PiolaVersion, SecondOrder, LFact, LFactFound, EdgeBasis

  REAL(KIND=dp) :: Relax

  TYPE(Variable_t), POINTER :: Var, FixJVar, CoordVar
  TYPE(Matrix_t), POINTER :: A
  TYPE(ListMatrix_t), POINTER :: BasicCycles(:)
  TYPE(ValueList_t), POINTER :: CompParams
  
  INTEGER :: n_n, n_e
  INTEGER, POINTER :: Vperm(:), Aperm(:)
  REAL(KIND=dp), POINTER :: Avals(:), Vvals(:)

  SAVE STIFF, LOAD, MASS, FORCE, Tcoef, GapLength, AirGapMu, &
       Acoef, Cwrk, LamThick, LamCond, Wbase, RotM, AllocationsDone
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN	

  CALL Info('WhitneyAVSolver','-------------------------------------------',Level=8 )
  CALL Info('WhitneyAVSolver','Solving the AV equations with edge elements',Level=5 )

  PiolaVersion = GetLogical( GetSolverParams(), 'Use Piola Transform', Found )
  SecondOrder = GetLogical( GetSolverParams(), 'Quadratic Approximation', Found )
  IF (PiolaVersion) THEN
    CALL Info('WhitneyAVSolver', &
        'Using Piola Transformed element basis functions',Level=4)
    CALL Info('WhitneyAVSolver', &
        'Currently only homogeneous Dirichlet and Neumann BCs for the vector potential supported',Level=4)
    CALL Info('WhitneyAVSolver', &
        'The option > Use Tree Gauge < is not available',Level=4)
    IF (SecondOrder) &
         CALL Info('WhitneyAVSolver', &
        'Using quadratic approximation, pyramidical elements are not yet available',Level=4)
  END IF

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()
  nNodes = Mesh % NumberOfNodes
  Perm => Solver % Variable % Perm
  Vecpot => Solver % Variable % Values

  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough
     ALLOCATE( FORCE(N), LOAD(7,N), STIFF(N,N), &
          MASS(N,N), Tcoef(3,3,N), GapLength(N), &
          AirGapMu(N), Acoef(N), LamThick(N), &
          LamCond(N), Wbase(N), RotM(3,3,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'WhitneyAVSolver', 'Memory allocation error.' )
     END IF

     IF(GetString(GetSolverParams(),'Linear System Solver')=='block') THEN
       n = Mesh % NumberOfNodes
       n_n = COUNT(Perm(1:n)>0)
       n_e = COUNT(Perm(n+1:)>0)
       ALLOCATE( Avals(n_e), Vvals(n_n) )
       Vvals = Vecpot(1:n)
       Avals = Vecpot(n+1:)

       ALLOCATE(Aperm(SIZE(Perm)),Vperm(SIZE(Perm)))
       Aperm = 0; Vperm = 0
       DO i=1,SIZE(AVals)
         Aperm(i+n) = i
       END DO
       DO i=1,SIZE(Vvals)
         Vperm(i)=i
       END DO

       CALL VariableAdd(Mesh % Variables,Mesh,Solver, & 
                 GetVarName(Solver % Variable)//' 1',1,Vvals,Vperm)

       CALL VariableAdd(Mesh % Variables,Mesh,Solver, & 
                 GetVarName(Solver % Variable)//' 2',1,Avals,Aperm)
     END IF

     NULLIFY( Cwrk )

     AllocationsDone = .TRUE.
  END IF

  ConstantSystem = GetLogical( GetSolverParams(), &
       'Constant System', Found )

  ConstantBulk = GetLogical( GetSolverParams(), &
       'Constant Bulk System', Found )

  SkipAssembly = DoneAssembly.AND.(ConstantBulk.OR.ConstantSystem)

  FixJ = GetLogical(GetSolverParams(),'Fix input Current Density', Found)
  IF (.NOT. Found .AND. .NOT. Transient ) THEN
    ! Only fix the current density if there is one
    FixJ = ListCheckPrefixAnyBodyForce(Model, 'Current Density')
  END IF
  IF (FixJ) THEN
    CALL JfixPotentialSolver(Model,Solver,dt,Transient)
    FixJVar => VariableGet(Mesh % Variables, 'Jfix')
  END IF

  ! 
  ! Use vec.pot. dofs only for convergence:
  ! ----------------------------------------
  CALL ListAddInteger(GetSolverParams(),'Norm Permutation',nNodes+1)

  ! 
  ! Resolve internal non.linearities, if requeted:
  ! ----------------------------------------------
  NoIterationsMax = GetInteger( GetSolverParams(), &
	      'Nonlinear System Max Iterations',Found)
  IF(.NOT. Found) NoIterationsMax = 1

  NoIterationsMin = GetInteger( GetSolverParams(), &
	      'Nonlinear System Min Iterations',Found)
  IF(.NOT. Found) NoIterationsMin = 1

  LFact = GetLogical( GetSolverParams(),'Linear System Refactorize', LFactFound )
  IF ( dt /= PrevDT .AND. LFactFound .AND. .NOT. LFact ) THEN
    CALL ListAddLogical( getSolverParams(), 'Linear System Refactorize', .TRUE. )
  END IF
  EdgeBasis = .NOT.LFactFound .AND. GetLogical( GetSolverParams(), 'Edge Basis', Found )

  DO i=1,NoIterationsMax
    IF( NoIterationsMax > 1 ) THEN
      CALL Info('WhitneyAVSolver','Nonlinear iteration: '//TRIM(I2S(i)),Level=8 )
    END IF
    IF( DoSolve(i) ) THEN
      IF(i>=NoIterationsMin) EXIT
    END IF
    IF( EdgeBasis ) CALL ListAddLogical(GetSolverParams(),'Linear System Refactorize',.FALSE.)
  END DO
  IF ( EdgeBasis ) CALL ListRemove( GetSolverParams(), 'Linear System Refactorize' )

  IF ( dt /= PrevDT .AND. LFactFound .AND. .NOT. LFact ) THEN
    CALL ListAddLogical( GetSolverParams(), 'Linear System Refactorize', .FALSE. )
  END IF
  PrevDT = dt

  CALL CalculateLumped(Model % NumberOfBodyForces)

  CoordVar => VariableGet(Mesh % Variables,'Coordinates')
  IF(ASSOCIATED(CoordVar)) THEN
    DO i=1,Mesh % NumberOfNodes
      j = 3*(CoordVar % Perm(i)-1)
      CoordVar % Values(j+1) = Mesh % Nodes % x(i)
      CoordVar % Values(j+2) = Mesh % Nodes % y(i)
      CoordVar % Values(j+3) = Mesh % Nodes % z(i)
    END DO
  END IF

  CALL Info('WhitneyAVSolver','All done',Level=8 )
  CALL Info('WhitneyAVSolver','-------------------------------------------',Level=8 )


CONTAINS

!------------------------------------------------------------------------------
  LOGICAL FUNCTION DoSolve(IterNo) RESULT(Converged)
!------------------------------------------------------------------------------
   CHARACTER(LEN=MAX_NAME_LEN) :: potname
   INTEGER :: i,j,k,t,n,nd,nb,IterNo

   REAL(KIND=dp)::TOL,Norm,PrevNorm, NonLinError, LinTol, RelTol, BaseTol
   LOGICAL :: Found, FoundMagnetization, CalculateNonlinearResidual, LFactFound
   LOGICAL :: AdaptiveTols

   TYPE(Matrix_t), POINTER :: MMatrix
   REAL(KIND=dp), POINTER :: Mx(:), Mb(:), Mr(:)
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpRVec, TmpRHSVec   
   CHARACTER(LEN=MAX_NAME_LEN) :: ConvergenceType

   SAVE TmpRHSVec, TmpRVec
  !-----------------
  !System assembly:
  !-----------------

  A => GetMatrix()

300 CONTINUE

  IF ( SkipAssembly) THEN
    A % RHS = A % BulkRHS
    A % Values = A % BulkValues
    IF ( ConstantBulk    ) GOTO 100
    IF ( ConstantSystem  ) GOTO 200
  END IF

  CALL DefaultInitialize()
  Active = GetNOFActive()
  DO t=1,active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes() ! kulmat
     nd = GetElementNOFDOFs()  ! vapausasteet
     nb = GetElementNOFBDOFs()  ! sisäiset vapausasteet

     LOAD = 0.0d0
     BodyForce => GetBodyForce()
     FoundMagnetization = .FALSE.
     IF ( ASSOCIATED(BodyForce) ) THEN

       CALL GetRealVector( BodyForce, Load(1:3,1:n), 'Current Density', Found )
       CALL GetRealVector( BodyForce, Load(4:6,1:n), &
                'Magnetization', FoundMagnetization )
	Load(7,1:n) = GetReal( BodyForce, 'Electric Potential', Found )
     END IF

     Material => GetMaterial( Element )

     IF(ASSOCIATED(Material).AND..NOT.FoundMagnetization) THEN
       CALL GetRealVector( Material, Load(4:6,1:n), &
                'Magnetization', FoundMagnetization )
     END IF
     
     CoilBody = .FALSE.
     CompParams => GetComponentParams( Element )
     CoilType = ''
     RotM = 0._dp
     IF (ASSOCIATED(CompParams)) THEN
       CoilType = GetString(CompParams, 'Coil Type', Found)
       IF (Found) THEN
         SELECT CASE (CoilType)
         CASE ('stranded')
            CoilBody = .TRUE.
         CASE ('massive')
            CoilBody = .TRUE.
         CASE ('foil winding')
            CoilBody = .TRUE.
            CALL GetElementRotM(Element, RotM, n)
         CASE DEFAULT
            CALL Fatal ('WhitneyAVSolver', 'Non existent Coil Type Chosen!')
         END SELECT
       END IF
     END IF

     Acoef = 0.0d0
     Tcoef = 0.0d0
     Material => GetMaterial( Element )
     IF ( ASSOCIATED(Material) ) THEN
       CALL GetReluctivity(Material,Acoef,n)

!------------------------------------------------------------------------------
!      Read conductivity values (might be a tensor)
!------------------------------------------------------------------------------

       Tcoef = GetElectricConductivityTensor(Element,n,'re',CoilBody,CoilType)

       LaminateStackModel = GetString( Material, 'Laminate Stack Model', LaminateStack )
       IF (.NOT. LaminateStack) LaminateStackModel = ''
     END IF


     LamThick=0d0
     LamCond=0d0
     IF (LaminateStack) THEN
       SELECT CASE(LaminateStackModel)
       CASE('low-frequency model')
         LamThick(1:n) = GetReal( Material, 'Laminate Thickness', Found )
         IF (.NOT. Found) CALL Fatal('WhitneyAVSolver', 'Laminate Thickness not found!')

         LamCond(1:n) = GetReal( Material, 'Laminate Stack Conductivity', Found )
         IF (.NOT. Found) CALL Fatal('WhitneyAVSolver', 'Laminate Stack Conductivity not found!')

       CASE DEFAULT
         CALL WARN('WhitneyAVSolver', 'Nonexistent Laminate Stack Model chosen!')
       END SELECT
     END IF


     !Get element local matrix and rhs vector:
     !----------------------------------------
     CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, &
          Tcoef, Acoef, LaminateStack, LaminateStackModel, &
          LamThick, LamCond, CoilBody, CoilType, RotM, &
          Element, n, nd+nb, PiolaVersion, SecondOrder )    

     !Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     IF (Transient) CALL DefaultUpdateMass(MASS)
     CALL DefaultUpdateEquations(STIFF,FORCE)
  END DO

  CALL DefaultFinishBulkAssembly(BulkUpdate=ConstantBulk)

100 CONTINUE

  !
  ! Robin type of BC in terms of H:
  !--------------------------------
  Active = GetNOFBoundaryElements()
  DO t=1,Active
     Element => GetBoundaryElement(t)
     IF (.NOT. ActiveBoundaryElement()) CYCLE
     BC=>GetBC()
     IF (.NOT. ASSOCIATED(BC) ) CYCLE
     
     SELECT CASE(GetElementFamily())
     CASE(1)
       CYCLE
     CASE(2)
       k = GetBoundaryEdgeIndex(Element,1); Element => Mesh % Edges(k)
     CASE(3,4)
       k = GetBoundaryFaceIndex(Element)  ; Element => Mesh % Faces(k)
     END SELECT

     Model % CurrentElement => Element
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     CALL GetRealVector( BC, Load(1:3,1:n), 'Magnetic Field Strength', Found )

     Acoef(1:n) = GetReal( BC, 'Magnetic Transfer Coefficient', Found ) !???

     Load(4,1:n) = GetReal( BC, 'Electric Current Density', Found )
     IF (.NOT. Found) &
          Load(4,1:n) = GetReal( BC, 'Electric Flux', Found )       
     Load(5,1:n) = GetReal( BC, 'Electric Transfer Coefficient', Found )

     !If air gap length keyword is detected, use air gap boundary condition
     GapLength=GetConstReal( BC, 'Air Gap Length', Found)
     IF (Found) THEN
       AirGapMu=GetConstReal( BC, 'Air Gap Relative Permeability', Found)
       IF (.NOT. Found) AirGapMu=1d0 ! if not found default to "air" property
       CALL LocalMatrixAirGapBC(STIFF,FORCE,LOAD,GapLength,AirGapMu,Element,n,nd )
     ELSE
       CALL LocalMatrixBC(STIFF,FORCE,LOAD,Acoef,Element,n,nd )
     END IF

     CALL DefaultUpdateEquations(STIFF,FORCE,Element)
  END DO

  CALL DefaultFinishBoundaryAssembly(BulkUpdate=ConstantSystem)

  DoneAssembly = .TRUE.

200 CONTINUE

  ! This is now automatically invoked as the time integration is set gloabl in the Solver_init
  ! IF ( Transient ) CALL Default1stOrderTimeGlobal()
  CALL DefaultFinishAssembly()

  ! The following is now commented out intentionally, as ensuring the uniqueness of the scalar
  ! potential via using a Dirichlet constraint should not be a prerequisite for obtaining 
  ! a solution (although not unique) if the data satisfies a compatibility condition.
  ! ----------------------------------------------------------------------------------------
  IF (.FALSE.) THEN
     ! Check that the auxliary potential fixed somehow,
     ! if not disable the variable (set to zero)
     ! ------------------------------------------------
     IF  (FixJ) THEN
        Found = .FALSE.
        potname = Solver % Variable % Name(1:Solver % Variable % NameLen)

        Found = ListCheckPresentAnyBC(Model,potname)

        IF (.NOT. Found ) THEN
           DO i=1,Mesh % NumberOfNodes
              J = Solver % Variable % Perm(i)
              IF(J>0) THEN
                 CALL ZeroRow(A,j)
                 CALL SetMatrixElement(A,j,j,1._dp)
                 A % RHS(j)=0._dp
              END IF
           END DO
        END IF
     END IF
  END IF

  !
  ! Dirichlet BCs in terms of vector potential A:
  ! ---------------------------------------------
  CALL DefaultDirichletBCs()

  ! Dirichlet BCs in terms of magnetic flux density B:
  ! --------------------------------------------------
  CALL DirichletAfromB()
  CALL ConstrainUnused(A)

  !
  ! Gauge tree, if requested or using direct solver:
  ! ------------------------------------------------
  TG=GetLogical(GetSolverParams(), 'Use tree gauge', Found)
  IF (.NOT. Found) TG=GetString(GetSolverParams(), &
     'Linear System Solver',Found)=='direct'

  IF (TG) THEN
    CALL GaugeTree()
    WRITE(Message,*) 'Volume tree edges: ', &
           TRIM(i2s(COUNT(TreeEdges))),     &
             ' of total: ',Mesh % NumberOfEdges
    CALL Info('WhitneyAVSolver: ', Message, Level=5)
  END IF

  IF (DefaultLineSearch(Converged)) RETURN
  IF ( Converged ) GOTO 10

  ! The following gives the user an option to adapt the linear system convergence tolerance
  ! adaptively:
  !-------------------------------------------------------------------------------------------
  AdaptiveTols = ListGetLogical(GetSolverParams(), 'Linear System Adaptive Tolerance', Found)
  IF (AdaptiveTols) THEN
     IF (iterno == 1) THEN
        BaseTol = ListGetConstReal(GetSolverParams(), 'Linear System Base Tolerance', Found)
        IF (Found) CALL ListAddConstReal(GetSolverParams(), &
             'Linear System Convergence Tolerance', BaseTol)
     ELSE
        IF (.NOT. ListCheckPresent(GetSolverParams(), 'Linear System Relative Tolerance')) THEN
           RelTol = 1.0d-2
        ELSE
           RelTol = ListGetConstReal(GetSolverParams(), 'Linear System Relative Tolerance')
        END IF
        IF (.NOT. ListCheckPresent(GetSolverParams(), 'Linear System Base Tolerance')) THEN
           BaseTol = 1.0d-3
        ELSE
           BaseTol = ListGetConstReal(GetSolverParams(), 'Linear System Base Tolerance')
        END IF
        LinTol = MIN(BaseTol, RelTol * Solver % Variable % NonlinChange)
        CALL ListAddConstReal(GetSolverParams(),'Linear System Convergence Tolerance', &
             LinTol)
     END IF
  END IF

  norm = DefaultSolve()
  Converged = Solver % Variable % NonlinConverged==1

10 CONTINUE

  IF ( ALLOCATED(FluxMap) ) DEALLOCATE(FluxMap)
  IF ( ALLOCATED(TreeEdges) ) DEALLOCATE(TreeEdges)

! CALL WriteResults  ! debugging helper


!------------------------------------------------------------------------------
 END FUNCTION DoSolve
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
 

   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF
   END IF

   RotM = 0._dp
   DO j = 1, n
     DO k=1,RotMvar % DOFs
       RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
             RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE ConstrainUnused(A)
!------------------------------------------------------------------------------
  TYPE(Matrix_t) :: A
!------------------------------------------------------------------------------
  INTEGER :: i,j,n

  REAL(kind=dp), ALLOCATABLE :: dDiag(:)
!------------------------------------------------------------------------------

  n = A % NumberOfRows
  ALLOCATE(dDiag(n)); dDiag = 0._dp

  DO i=1,n
    j = A % Diag(i)
    IF (j>0) dDiag(i) = A % Values(j)
  END DO
  IF (ParEnv % PEs>1) CALL ParallelSumVector(A, dDiag)

  n = SIZE(Solver % Variable % Perm)
  DO i=1,n
    J = Solver % Variable % Perm(i)
    IF(j>0) THEN
      IF ( dDiag(j)==0._dp ) THEN
        CALL ZeroRow(A,j)
        CALL SetMatrixElement(A,j,j,1._dp)
        A % RHS(j)=0._dp
        IF(ALLOCATED(A % ConstrainedDOF)) A % ConstrainedDOF(j)=.TRUE.
      END IF
    END IF
  END DO
!------------------------------------------------------------------------------
 END SUBROUTINE ConstrainUnused
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumped(nbf)
!------------------------------------------------------------------------------
   INTEGER::nbf
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: torq,a(nbf),u(nbf),IMoment,IA,zforce,zzforce
   INTEGER :: i,bfid,n,nd,EdgeBasisDegree
   LOGICAL :: Found, CalcTorque,CalcPotential,CalcInertia
   TYPE(ValueList_t),POINTER::Params
   TYPE(ELement_t), POINTER :: Element, Parent
!------------------------------------------------------------------------------

   CalcTorque = ListCheckPresentAnyBody(Model,'r inner')
   CalcPotential = ListGetLogicalAnyBodyForce( Model,'Calculate Potential')
   CalcInertia = ListGetLogicalAnyBody( Model,'Calculate Inertial Moment')

   IF(.NOT. (CalcTorque .OR. CalcPotential .OR. CalcInertia ) ) RETURN

   EdgeBasisDegree = 1
   IF (SecondOrder) EdgeBasisDegree = 2

   U=0._dp; a=0._dp; torq=0._dp; IMoment=0._dp;IA=0; zforce=0
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     IF( CalcTorque ) THEN
       CALL Torque(Torq,Element,n,nd,EdgeBasisDegree)
       CALL AxialForce(zforce,Element,n,nd,EdgeBasisDegree)
     END IF

     IF( CalcPotential ) THEN
       Params=>GetBodyForce(Element)
       IF(ASSOCIATED(Params)) THEN
         bfid=GetBodyForceId(Element)
         IF(GetLogical(Params,'Calculate Potential',Found)) &
             CALL Potential(u(bfid),a(bfid),Element,n,nd,EdgeBasisDegree)
       END IF
     END IF

     IF( CalcInertia ) THEN
       Params=>GetBodyParams(Element)
       IF(ASSOCIATED(Params)) THEN
         IF(GetLogical(Params,'Calculate Inertial Moment',Found)) &
             CALL InertialMoment(IMoment,IA,Element,n,nd)
       END IF
     END IF
   END DO

   zzforce = 0
   IF(ListGetLogicalAnyBC(Model,'Calculate Axial Force')) THEN
     DO i=1,Mesh % NumberOFBoundaryElements
       Element => GetBoundaryElement(i)
       IF (.NOT.GetLogical(GetBC(), 'Calculate Axial Force', Found ) ) CYCLE

       Parent => Element % BoundaryInfo % Left
       n  = GetELementNofNodes(Parent)
       nd = GetELementNofDOFs(Parent)
       CALL AxialForceSurf(zzforce,Element,n,nd,EdgeBasisDegree)
     END DO
   END IF

   IF( CalcPotential ) THEN
     DO i=1,nbf
       a(i) = ParallelReduction(a(i))
       u(i) = ParallelReduction(u(i))
     END DO
     
     DO i=1,nbf
       IF(a(i)>0) THEN
         CALL ListAddConstReal(Model % Simulation,'res: Potential / bodyforce ' &
             //TRIM(i2s(i)),u(i)/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
             //TRIM(i2s(i)),a(i))
       END IF
     END DO
   END IF

   IF( CalcTorque ) THEN
     Torq = ParallelReduction(Torq)
     CALL ListAddConstReal(Model % Simulation,'res: Air Gap Torque', Torq)
     zforce = ParallelReduction(zforce)
     CALL ListAddConstReal(Model % Simulation,'res: Axial Force(vol)', zforce)

     IF(ListGetLogicalAnyBC(Model,'Calculate Axial Force')) THEN
       zzforce = ParallelReduction(zzforce)
       CALL ListAddConstReal(Model % Simulation,'res: Axial force(surf)', zzforce )
     END IF
   END IF

   IF( CalcInertia ) THEN
     IMoment = ParallelReduction(IMoment)
     IA = ParallelReduction(IA)
     CALL ListAddConstReal(Model % Simulation,'res: Inertial Volume', IA)
     CALL ListAddConstReal(Model % Simulation,'res: Inertial Moment', IMoment)
   END IF

!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumped
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE InertialMoment(U,A,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    REAL(KIND=dp)::U,a
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), DetJ,x,y,r,Density(n)
    INTEGER :: t
    LOGICAL :: stat,Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
    IF(.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )

      x = SUM(Nodes % x(1:n)*Basis(1:n))
      y = SUM(Nodes % y(1:n)*Basis(1:n))
      r = SQRT(x**2+y**2)
      A = A + IP % s(t)*detJ
      U = U + IP % s(t)*detJ*R*SUM(Density(1:n)*Basis(1:n))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE InertialMoment
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Torque(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    INTEGER :: n,nd,EdgeBasisDegree
    REAL(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(nd),x,y,r,r0,r1,Br,Bp,Bx,By,B(3,nd),Wbasis(nd,3),RotWBasis(nd,3)
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    r0 = GetCReal(GetBodyParams(),'r inner',Found)
    r1 = GetCReal(GetBodyParams(),'r outer',Found)
    IF (.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n
    r = SQRT(x**2+y**2)
    IF (r<r0.OR.r>r1) RETURN

    CALL GetLocalSolution(POT, UElement=Element)
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      IF (PiolaVersion) THEN
        stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
             DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
             BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
      END IF

      x = SUM(Nodes % x(1:n)*Basis(1:n))
      y = SUM(Nodes % y(1:n)*Basis(1:n))
      r = SQRT(x**2+y**2)

      Bx =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,1))
      By =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,2))
      Br =  x/r*Bx + y/r*By
      Bp = -y/r*Bx + x/r*By
      U = U + IP % s(t)*detJ*r*Br*Bp/(PI*4.0d-7*(r1-r0))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Torque
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AxialForce(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    INTEGER :: n,nd,EdgeBasisDegree
    REAL(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(nd),x,y,r,r0,r1,Bx,By,Bz,B(3,nd),Wbasis(nd,3),RotWBasis(nd,3)
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    r0 = GetCReal(GetBodyParams(),'r inner',Found)
    r1 = GetCReal(GetBodyParams(),'r outer',Found)
    IF (.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n
    r = SQRT(x**2+y**2)
    IF (r<r0.OR.r>r1) RETURN

    CALL GetLocalSolution(POT, UElement=Element)
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      IF (PiolaVersion) THEN
        stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
             DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
             BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
      END IF

      x = SUM(Nodes % x(1:n)*Basis(1:n))
      y = SUM(Nodes % y(1:n)*Basis(1:n))
      r = SQRT(x**2+y**2)
      x = x/r; y = y/r

      Bx =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,1))
      By =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,2))
      Bz =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,3))
      U = U + IP % s(t)*detJ*1*(Bx*Bz*x+By*Bz*y)/(PI*4.0d-7*(r1-r0))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AxialForce
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AxialForceSurf(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    INTEGER :: n,nd,EdgeBasisDegree
    REAL(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, Pdetj, uu,v,w, &
             POT(nd),x,y,r,r0,r1,Wbasis(nd,3),RotWBasis(nd,3)
    REAL(KIND=dp) :: B(3,nd), Bx, By, Bz
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes, PNodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    CALL GetElementNodes( Nodes, Element )
    Parent => Element % BoundaryInfo % Left
    CALL GetElementNodes( PNodes, Parent )

    CALL GetLocalSolution(POT, UElement=Parent )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      IF (PiolaVersion) THEN
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )

        CALL GetParentUVW(Element,GetElementNOFNodes(Element),Parent,n,uu,v,w,Basis)
 
        stat = EdgeElementInfo( Parent, PNodes, uu, v, w, &
              DetF = PDetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
              BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)

      ELSE

        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )

        CALL GetParentUVW(Element,GetElementNOFNodes(Element),Parent,n,uu,v,w,Basis)
        stat = ElementInfo( Parent, PNodes, uu,v,w, pdetJ, Basis, dBasisdx )

        CALL GetEdgeBasis(Parent,WBasis,RotWBasis,Basis,dBasisdx)
      END IF

      x = SUM(Basis(1:n) * PNodes % x(1:n))
      y = SUM(Basis(1:n) * PNodes % y(1:n))
      r = SQRT(x**2 + y**2)
      x=x/r; y=y/r

      Bx =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,1)) 
      By =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,2))
      Bz =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,3))
      U = U + IP % s(t) * detJ * (Bx*Bz + By*Bz) /(PI*4.0d-7) !/ 2
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AxialForceSurf
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Potential( U, A, Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: U,A
    INTEGER :: n, nd, EdgeBasisDegree
    TYPE(Element_t) :: Element

    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3),DetJ,POT(nd),pPOT(nd), &
          dPOT(nd),wBasis(nd,3),rotWBasis(nd,3),Wpot(nd),w(3)
    INTEGER :: t
    LOGICAL :: stat, WbaseFound
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    CALL GetElementNodes( Nodes )

    CALL GetLocalSolution(POT,UElement=Element)
    CALL GetLocalSolution(pPOT,tstep=-1,UElement=Element)
    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      dPot = (POT - pPOT)/dt
    ELSE
      dPot = 1.5_dp*POT - 2*pPOT
      CALL GetLocalSolution(pPOT,tstep=-2,UElement=Element)
      dPot = (dPOT + 0.5_dp*pPOT)/dt
    END IF

    CALL GetLocalSolution(Wpot,'W',UElement=Element)
    W = [0._dp, 0._dp, 1._dp]
    WbaseFound = ANY(Wpot(1:n)/=0._dp)

    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      IF (PiolaVersion) THEN
        stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
             DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, dBasisdx = dBasisdx, &
             BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)         
      ELSE
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
      END IF

      IF(WBaseFound) W = MATMUL(Wpot(1:n),dBasisdx(1:n,:))

      A = A + IP % s(t) * detJ
      U = U + IP % s(t) * detJ * SUM(dPot(n+1:nd)*MATMUL(WBasis(1:nd-n,:),w))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Potential
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE GaugeTree()
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry
    TYPE(ListMatrix_t), POINTER :: Alist(:)
    INTEGER :: i,j,k,l,n,Start
    LOGICAL, ALLOCATABLE :: Done(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Edge, Boundary
!------------------------------------------------------------------------------

    IF ( .NOT. ALLOCATED(TreeEdges) ) THEN
      ALLOCATE(TreeEdges(Mesh % NumberOfEdges)); TreeEdges = .FALSE.
    END IF

    n = Mesh % NumberOfNodes
    ALLOCATE(Done(n)); Done=.FALSE.

    ! 
    ! Skip Dirichlet BCs in terms of A:
    ! ---------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(i)
      IF ( .NOT.ActiveBoundaryElement()) CYCLE
      IF ( GetElementFamily()==1 ) CYCLE
      BC => GetBC()
      IF (.NOT.ASSOCIATED(BC)) CYCLE
      IF (.NOT.ListCheckPresent( BC, &
           TRIM(Solver % Variable % Name)//' {e}')) CYCLE
 
      j=1; k=GetBoundaryEdgeIndex(Boundary,j)
      DO WHILE(k>0)
        Edge => Mesh % Edges(k)
        Done(Edge % NodeIndexes) = .TRUE.
        j=j+1; k=GetBoundaryEdgeIndex(Boundary,j)
      END DO
    END DO

    ! 
    ! Skip Dirichlet BCs in terms of B:
    ! ---------------------------------
    DO i=1,FluxCount
      j = FluxMap(i)
      IF ( Perm(j+n)<=0 ) CYCLE
      Edge => Mesh % Edges(j)
      Done(Edge % NodeIndexes)=.TRUE.
    END DO

    !
    ! node -> edge list
    ! -----------------
    Alist => NULL()
    n = Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
      Edge => Mesh % Edges(i)
      IF ( Perm(i+n)<=0 ) CYCLE
      DO j=1,Edge % TYPE % NumberOfNodes
        k=Edge % NodeIndexes(j)
        Aentry=>List_GetMatrixIndex(Alist,k,i)
      END DO
    END DO

    !
    ! generate the tree for all (perhaps disconnected) parts:
    ! -------------------------------------------------------
    DO WHILE(.NOT.ALL(Done))
      DO Start=1,n
        IF (.NOT. Done(Start)) EXIT
      END DO
      CALL DepthFirstSearch(Alist,Done,Start)
    END DO
    DEALLOCATE(Done)
    CALL List_FreeMatrix(SIZE(Alist),Alist)
!------------------------------------------------------------------------------
  END SUBROUTINE GaugeTree
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GaugeTreeFluxBC()
!------------------------------------------------------------------------------
!   TYPE(Mesh_t) :: Mesh
!   LOGICAL, ALLOCATABLE :: TreeEdges(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry, Ltmp
    TYPE(ListMatrix_t), POINTER :: Alist(:)
    INTEGER :: i,j,k,l,n,Start,nCount,fixedge
    LOGICAL, ALLOCATABLE :: Done(:)
    INTEGER, ALLOCATABLE :: NodeList(:)
    TYPE(Element_t), POINTER :: Edge, Boundary
!------------------------------------------------------------------------------

    IF ( .NOT. ALLOCATED(TreeEdges) ) THEN
      ALLOCATE(TreeEdges(Mesh % NumberOfEdges)); TreeEdges = .FALSE.
    END IF

    n = Mesh % NumberOfNodes
    ALLOCATE(Done(n)); Done=.FALSE.

    !
    ! list the candidate nodes:
    ! -------------------------
    DO i=1,FluxCount
      j = FluxMap(i)
      Edge => Mesh % Edges(j)
      Done(Edge % NodeIndexes)=.TRUE.
    END DO

    ALLOCATE(NodeList(COUNT(Done)))
    nCount = 0
    DO i=1,n
      IF ( Done(i) ) THEN
        nCount = nCount+1
        NodeList(nCount)=i
      END IF
    END DO

    Done=.FALSE.
    DO i=1,FluxCount
      IF ( TreeEdges(FluxMap(i)) ) THEN
        Edge => Mesh % Edges(FluxMap(i))
        Done(Edge % NodeIndexes)=.TRUE.
      END IF
    END DO

    ! 
    ! Skip Dirichlet BCs in terms of A:
    ! ---------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(i)
      IF ( .NOT.ActiveBoundaryElement()) CYCLE
      IF ( GetElementFamily()==1 ) CYCLE
      BC => GetBC()
      IF (.NOT.ASSOCIATED(BC)) CYCLE
      IF (.NOT.ListCheckPresent( BC, &
           TRIM(Solver % Variable % Name)//' {e}')) CYCLE
 
      j=1; k=GetBoundaryEdgeIndex(Boundary,j)
      DO WHILE(k>0)
        Edge => Mesh % Edges(k)
        TreeEdges(k) = .TRUE.
        Done(Edge % NodeIndexes) = .TRUE.
        j=j+1; k=GetBoundaryEdgeIndex(Boundary,j)
      END DO
    END DO

    ! node -> edge list
    ! -----------------
    Alist => NULL()
    DO i=1,FluxCount
      j = FluxMap(i)
      IF ( Perm(j+n)<=0 ) CYCLE

      Edge => Mesh % Edges(j)
      DO k=1,Edge % TYPE % NumberOfNodes
        l=Edge % NodeIndexes(k)
        Aentry=>List_GetMatrixIndex(Alist,l,j)
      END DO
    END DO
 
    ! generate the tree for all (perhaps disconnected) parts:
    ! -------------------------------------------------------
    DO WHILE(.NOT.ALL(Done(NodeList)))
      DO i=1,nCount
        Start = NodeList(i)
        IF ( .NOT. Done(Start) ) EXIT
      END DO
      CALL BreadthFirstSearch(Alist,Done,start,nCount,NodeList)
    END DO
    DEALLOCATE(Done,NodeList)
    CALL List_FreeMatrix(SIZE(Alist),Alist)
!------------------------------------------------------------------------------
  END SUBROUTINE GaugeTreeFluxBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BreadthFirstSearch(Alist,done,start,nCount,NodeList)
!------------------------------------------------------------------------------
    INTEGER :: start,nCount,NodeList(:)
    LOGICAL :: Done(:)
    TYPE(ListMatrix_t) :: Alist(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry, Ltmp, Btmp
    INTEGER :: i,j,k,l,n,m,ll,IF,bcycle
    TYPE(Element_t), POINTER :: Edge,Edge1,Boundary
    LOGICAL, ALLOCATABLE :: DoneL(:)
    INTEGER, ALLOCATABLE :: Fifo(:), Previous(:), FiFo1(:)
!------------------------------------------------------------------------------
   ALLOCATE(DoneL(Mesh % NumberOfEdges)); DoneL=.FALSE.
   ALLOCATE(Fifo(FluxCount),FiFo1(FluxCount))
   ALLOCATE(Previous(Mesh % NumberOfNodes)); Previous=0;

   IF = 0; m=0
   DO i=1,nCount
     j = NodeList(i)
     IF ( Done(j) ) THEN
       m=m+1; fifo1(m)=j
       IF=IF+1; fifo(IF)=j
     END IF
   END DO

   IF ( IF>0 ) THEN
     DO WHILE(m>0)
       j = Fifo1(m); m=m-1

       Aentry => Alist(j) % Head
       DO WHILE(ASSOCIATED(Aentry))
         k = Aentry % Index
         Aentry => Aentry % Next

         Edge => Mesh % Edges(k)
         IF (.NOT. TreeEdges(k) .OR. DoneL(k) ) CYCLE
         DoneL(k)=.TRUE.

         l = Edge % NodeIndexes(1)
         IF (l==j) l=Edge % NodeIndexes(2)

         IF=IF+1; Fifo(IF)=l
         m=m+1; Fifo1(m)=l
         Previous(l)=j
       END DO
     END DO
     Start = l
   END IF
   
   IF ( IF==0 ) THEN
     Done(Start)=.TRUE.
     IF=1; fifo(IF)=start;
   END IF

   Bcycle=0;
   ALLOCATE(BasicCycles(FluxCount))
   BasicCycles(:) % Degree = 0
   DO i=1,FluxCount
     BasicCycles(i) % Head => NULL()
   END DO

   DO WHILE(IF>0)
     j = Fifo(IF); IF=IF-1

     Aentry => Alist(j) % Head
     DO WHILE(ASSOCIATED(Aentry))
       k = Aentry % Index
       Aentry => Aentry % Next

       Edge => Mesh % Edges(k)
       IF ( DoneL(k) ) CYCLE
       DoneL(k)=.TRUE.

       l = Edge % NodeIndexes(1)
       IF (l==j) l=Edge % NodeIndexes(2)

       IF ( Done(l) ) THEN
         ! Generate fundamental cycle
         bcycle = bcycle+1
         CALL AddToCycle(bcycle,k)

         m = j
         DO WHILE(m/=Previous(l))
           Ltmp => Alist(m) % Head
           DO WHILE(ASSOCIATED(Ltmp))
             Edge1 => Mesh % Edges(Ltmp % Index)
             IF ( ANY(Edge1 % NodeIndexes(1:2)==Previous(m)) ) THEN
               CALL AddToCycle(bcycle,Ltmp % Index); EXIT
             END IF
             Ltmp=>Ltmp % Next
           END DO
           IF ( ANY(Edge1 % NodeIndexes(1:2) == l) ) EXIT
           m = Previous(m)
         END DO

         IF ( ALL(Edge1 % NodeIndexes(1:2) /= l) ) THEN
           ltmp => Alist(l) % Head
           DO WHILE(ASSOCIATED(ltmp))
             edge1 => Mesh % Edges(Ltmp % Index)
             IF ( ANY(Edge1 % NodeIndexes(1:2)==Previous(l)) ) THEN
               CALL AddToCycle(bcycle,Ltmp % Index); EXIT
             END IF
             ltmp=>ltmp % Next
           END DO
         END IF
       ELSE
         IF (.NOT.TreeEdges(k)) CALL SetDOFToValue(Solver,k,0._dp)
         IF=IF+1; Fifo(IF)=l
         Previous(l)=j
         Done(l)=.TRUE.
         TreeEdges(k) = .TRUE.
       END IF
     END DO
   END DO
   DEALLOCATE(Fifo, Fifo1, DoneL)
!------------------------------------------------------------------------------
  END SUBROUTINE BreadthFirstSearch
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddToCycle(bcycle,index)
    INTEGER :: bcycle,index
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Btmp

    ALLOCATE(Btmp); Btmp % Next => BasicCycles(bcycle) % Head;
    Btmp % Index = index; BasicCycles(bcycle) % Head => Btmp
    BasicCycles(bcycle) % Degree=BasicCycles(bcycle) % Degree+1
!------------------------------------------------------------------------------
  END SUBROUTINE AddToCycle
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DepthFirstSearch(Alist,done,i)
!------------------------------------------------------------------------------
    TYPE(ListMatrix_t) :: Alist(:)
    INTEGER :: i
    LOGICAL :: Done(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry
    INTEGER :: j,k,l,n
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

    ! To give better matrix conditioning some directional heuristics
    ! could be added,e.g. select the order of going through the nodes
    ! edge list here:

    Done(i) = .TRUE.

    Aentry => Alist(i) % Head
    DO WHILE(ASSOCIATED(Aentry))
      k = Aentry % Index
      Aentry => Aentry % Next

      Edge => Mesh % Edges(k)
      IF (ALL(Done(Edge % NodeIndexes))) CYCLE

      IF ( .NOT. TreeEdges(k)) CALL SetDOFToValue(Solver,k,0._dp)
      TreeEdges(k)=.TRUE.
      DO l=1,2
        n = Edge % NodeIndexes(l)
        IF (.NOT. Done(n)) CALL DepthFirstSearch(Alist,done,n)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE DepthFirstSearch
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( MASS, STIFF, FORCE, LOAD, &
            Tcoef, Acoef, LaminateStack, LaminateStackModel, &
            LamThick, LamCond, CoilBody, CoilType, RotM, &
            Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), MASS(:,:)
    REAL(KIND=dp) :: LOAD(:,:), Tcoef(:,:,:), Acoef(:), &
                     LamThick(:), LamCond(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Aloc(nd), JAC(nd,nd), mu, muder, B_ip(3), Babs
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), A, Acoefder(n), C(3,3), &
                     RotMLoc(3,3), RotM(3,3,n), velo(3), omega_velo(3,n), lorentz_velo(3,n)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ, L(3), G(3), M(3), FixJPot(nd)
    REAL(KIND=dp) :: LocalLamThick, LocalLamCond

    CHARACTER(LEN=MAX_NAME_LEN):: LaminateStackModel, CoilType

    LOGICAL :: Stat, Found, Newton, Cubic, HBCurve, LaminateStack, CoilBody, &
        HasVelocity, HasLorenzVelocity, HasAngularVelocity
    INTEGER :: t, i, j, p, q, np, siz, EdgeBasisDegree
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:),  &
           CubicCoeff(:)=>NULL(),HB(:,:)=>NULL()
    TYPE(ValueListEntry_t), POINTER :: Lst
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS  = 0.0d0

    JAC = 0._dp
    Newton = .FALSE.

    FixJpot = 0._dp
    IF (FixJ) THEN
      FixJPot(1:n) = FixJVar % Values(FixJVar % Perm(Element % NodeIndexes))
    END IF

    HasVelocity = .FALSE.
    IF(ASSOCIATED(BodyForce)) THEN
      CALL GetRealVector( BodyForce, omega_velo, 'Angular velocity', HasAngularVelocity)
      CALL GetRealVector( BodyForce, lorentz_velo, 'Lorentz velocity', HasLorenzVelocity)
      HasVelocity = HasAngularVelocity .OR. HasLorenzVelocity
    END IF

    CALL GetConstRealArray( Material, HB, 'H-B curve', HBCurve )
    siz = 0
    Cval => NULL()
    IF ( HBCurve ) THEN
      siz = SIZE(HB,1)
      IF(siz>1) THEN
        Bval=>HB(:,1)
        Hval=>HB(:,2)
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        IF (Cubic.AND..NOT.ASSOCIATED(CubicCoeff)) THEN
          ALLOCATE(CubicCoeff(siz))
          CALL CubicSpline(siz,Bval,Hval,CubicCoeff)
        END IF
        Cval=>CubicCoeff
        HBCurve = .TRUE.
      END IF
    END IF

    IF(siz<=1) THEN
      Lst => ListFind(Material,'H-B Curve',HBcurve)
      IF(HBcurve) THEN
        Cval => Lst % CubicCoeff
        Bval => Lst % TValues
        Hval => Lst % FValues(1,1,:)
      END IF
    END IF

    IF(HBCurve) THEN
      CALL GetScalarLocalSolution(Aloc)
      Newton = GetLogical( GetSolverParams(),'Newton-Raphson iteration',Found)
    END IF

    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree )

    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
    DO t=1,IP % n
       IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
               RotBasis = RotWBasis, dBasisdx = dBasisdx, &
               BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
       ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )

          CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
       END IF

       A = SUM( Basis(1:n) * Acoef(1:n) )
       mu = A

       ! Compute convection type term coming from rotation
       ! -------------------------------------------------
       IF(HasVelocity) THEN
         velo = 0.0_dp
         IF( HasAngularVelocity ) THEN
           DO i=1,n
             velo(1:3) = velo(1:3) + CrossProduct(omega_velo(1:3,i), [ &
                 basis(i) * Nodes % x(i), &
                 basis(i) * Nodes % y(i), &
                 basis(i) * Nodes % z(i)])
           END DO
         END IF
         IF( HasLorenzVelocity ) THEN
           velo(1:3) = velo(1:3) + [ &
               SUM(basis(1:n)*lorentz_velo(1,1:n)), &
               SUM(basis(1:n)*lorentz_velo(2,1:n)), &
               SUM(basis(1:n)*lorentz_velo(3,1:n))]
         END IF
       END IF

       ! Compute the conductivity tensor
       ! -------------------------------
       DO i=1,3
         DO j=1,3
           C(i,j) = SUM( Tcoef(i,j,1:n) * Basis(1:n) )
           IF(CoilBody .AND. CoilType /= 'massive') RotMLoc(i,j) = SUM( RotM(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       ! Transform the conductivity tensor (in case of a foil winding):
       ! --------------------------------------------------------------
       IF (CoilBody .AND. CoilType /= 'massive') C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))

       M = MATMUL( LOAD(4:6,1:n), Basis(1:n) )
       L = MATMUL( LOAD(1:3,1:n), Basis(1:n) )
       L = L-MATMUL(FixJPot(1:n), dBasisdx(1:n,:))

       LocalLamThick = SUM( Basis(1:n) * LamThick(1:n) )
       LocalLamCond = SUM( Basis(1:n) * LamCond(1:n) )

       ! Add -C * grad(V^s), where C is a tensor
       ! -----------------------------------------
       L = L-MATMUL(C, MATMUL(LOAD(7,1:n), dBasisdx(1:n,:)))

       IF ( HBCurve ) THEN
         B_ip = MATMUL( Aloc(np+1:nd), RotWBasis(1:nd-np,:) )
         babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )
         mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
         IF ( Newton ) THEN
           muder=(DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/babs
         END IF
       END IF

       ! ------------------------------------------------------------------
       ! Compute element stiffness matrix and force vector.
       ! If we calculate a coil, the nodal degrees of freedom are not used.
       ! ------------------------------------------------------------------
       IF (.NOT. CoilBody) THEN
         ! --------------------------------------------------------
         ! The constraint equation involving the scalar potential:
         !     -div(C*(dA/dt+grad(V)))=0
         ! --------------------------------------------------------
         IF ( Transient ) THEN

           DO i=1,np
             p = i
             DO j=1,np
               q = j

	       ! Compute the conductivity term <C grad V,grad v> for stiffness 
               ! matrix (anisotropy taken into account)
               ! -------------------------------------------

               IF ( SUM(C) /= 0._dp ) THEN
                 STIFF(p,q) = STIFF(p,q) + SUM(MATMUL(C, dBasisdx(q,:)) * dBasisdx(p,:))*detJ*IP % s(t)
               END IF
             END DO
             DO j=1,nd-np
               q = j+np

               ! Compute the conductivity term <C A,grad v> for 
               ! mass matrix (anisotropy taken into account)
               ! -------------------------------------------
               MASS(p,q) = MASS(p,q) + SUM(MATMUL(C, Wbasis(j,:))*dBasisdx(i,:))*detJ*IP % s(t)

	       ! Compute the conductivity term <C grad V, eta> for 
               ! stiffness matrix (anisotropy taken into account)
               ! ------------------------------------------------
               STIFF(q,p) = STIFF(q,p) + SUM(MATMUL(C, dBasisdx(i,:))*WBasis(j,:))*detJ*IP % s(t)
             END DO
           END DO
         ELSE
           ! ---------------------------------------------------------------
           ! This is the steady state branch. Adding the scalar pontential 
           ! solver as done in the following is reasonable only in the case 
           ! where the electrical conductivity is nonzero.
           ! ------------------------------------------------------------------
           IF (.NOT. LaminateStack ) THEN
             DO i=1,np
               p = i
               DO j=1,np
                 q = j

                 ! Compute the conductivity term <C grad V,grad v> for stiffness 
                 ! matrix (anisotropy taken into account)
                 ! -------------------------------------------
                 STIFF(p,q) = STIFF(p,q) + SUM(MATMUL(C, dBasisdx(j,:)) * dBasisdx(i,:))*detJ*IP % s(t)
               END DO

               DO j=1,nd-np
                 q = j+np
                 ! The equation for the vector potential:
                 ! Compute the conductivity term <C grad V, eta> for 
                 ! stiffness matrix (anisotropy taken into account)
                 ! ------------------------------------------------
                 STIFF(q,p) = STIFF(q,p) + SUM(MATMUL(C, dBasisdx(i,:))*WBasis(j,:))*detJ*IP % s(t)
               END DO
               IF ( HasVelocity ) THEN
                 DO j=1,nd-np
                   q = j+np
                   STIFF(p,q) = STIFF(p,q) - SUM(MATMUL(C,CrossProduct(velo, RotWBasis(j,:)))*dBasisdx(i,:))*detJ*IP % s(t)
                 END DO
               END IF
             END DO
           END IF
         END IF

       END IF ! (.NOT. CoilBody)

       !-----------------------------------------------------------------
       ! The equations for the H(curl)-conforming part, i.e. the equation 
       ! for the vector potential
       !    C*dA/dt + curl(nu*curl(A)) + C*grad(V) =  
       !    J^s + curl(M^s) - C*grad(V^s),
       ! with the term C*grad(V) already handled above.
       ! -----------------------------------------------------------------
       DO i = 1,nd-np
         p = i+np
         FORCE(p) = FORCE(p) + (SUM(L*WBasis(i,:)) + &
            SUM(M*RotWBasis(i,:)))*detJ*IP%s(t) 
         DO j = 1,nd-np
           q = j+np
           STIFF(p,q) = STIFF(p,q) + mu * SUM(RotWBasis(i,:)*RotWBasis(j,:))*detJ*IP%s(t) 
           IF( HasVelocity ) THEN
             STIFF(p,q) = STIFF(p,q) &
                 - SUM(WBasis(i,:)*MATMUL(C,CrossProduct(velo, RotWBasis(j,:))))*detJ*IP%s(t)
           END IF
           IF ( Newton ) THEN
             JAC(p,q) = JAC(p,q) + muder * SUM(B_ip(:)*RotWBasis(j,:)) * &
                 SUM(B_ip(:)*RotWBasis(i,:))*detJ*IP % s(t)/Babs
           END IF

           ! Compute the conductivity term <C A,eta> for 
           ! mass matrix (anisotropy taken into account)
           ! This is not used in the case of stranded coil:
           ! ----------------------------------------------
           IF (CoilType /= 'stranded') MASS(p,q) = MASS(p,q) + & 
               SUM(MATMUL(C, WBasis(j,:))*WBasis(i,:) )*detJ*IP % s(t)

           ! Compute the low frequency eddy term for laminate stack model.
           ! Note that the conductivity term <C A, eta> above can be used to 
           ! introduce the anisotropic effect in the laminate stack. However, 
           ! in classical approach of the Low-Frequency model it is set 
           ! to zero (this is left to the user to deside).
           ! -------------------------------------------------------------------
           IF (LaminateStackModel=='low-frequency model') THEN
               MASS(p,q) = MASS(p,q) + LocalLamCond * LocalLamthick**2/12d0 * & 
                           SUM(RotWBasis(i,:)*RotWBasis(j,:)) * detJ*IP % s(t)
           END IF

         END DO
       END DO
    END DO

    IF ( Newton ) THEN
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + JAC
      FORCE(1:nd) = FORCE(1:nd) + MATMUL(JAC,Aloc)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, LOAD, Bcoef, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: LOAD(:,:), Bcoef(:)
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,L(3),Normal(3),w0(3),w1(3)
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), B, F, TC
    LOGICAL :: Stat
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, k, ii,jj, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes, Element )

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n
       IF ( PiolaVersion ) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
               BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
       ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
       END IF

       B  = SUM(Basis(1:n) * Bcoef(1:n))
       L  = MATMUL(LOAD(1:3,1:n), Basis(1:n))

       F  = SUM(LOAD(4,1:n)*Basis(1:n))
       TC = SUM(LOAD(5,1:n)*Basis(1:n))

       Normal = NormalVector(Element,Nodes,IP % U(t), IP % V(t),.TRUE.)

       ! Compute element stiffness matrix and force vector:
       !---------------------------------------------------
       DO p=1,np
         FORCE(p) = FORCE(p) + F*Basis(p)*detJ*IP % s(t)
         DO q=1,np
           STIFF(p,q) = STIFF(p,q) + TC * &
                  Basis(p)*Basis(q)*detJ*IP % s(T)
         END DO
       END DO

       DO i = 1,nd-np
         w0 = -CrossProduct(Wbasis(i,:),Normal)
         p = i+np
         FORCE(p) = FORCE(p) + SUM(L*w0)*detJ*IP % s(t)
         DO j = 1,nd-np
           w1 = CrossProduct(Wbasis(j,:),Normal)
           q = j+np
           STIFF(p,q) = STIFF(p,q) - B * &
              SUM(w1*w0)*detJ*IP % s(t)
         END DO
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  FUNCTION LocalFluxBC( LOAD, Element, n, nd ) RESULT(Bn)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: LOAD(:,:), Bn
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Edge, Parent
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,L(3),Ln
    REAL(KIND=dp) :: Normal(3)
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,  Element )
    !
    ! Integrate (B,n) over boundary face:
    ! -----------------------------------
    IP = GaussPoints(Element)
    Bn = 0._dp
    DO t=1,IP % n
      stat = ElementInfo( Element,Nodes,IP % U(t),IP % V(t), &
                 IP % W(t),detJ,Basis,dBasisdx )

      Normal=NormalVector(Element,Nodes,IP % u(t),ip % v(t),.TRUE.)
      Ln = SUM(LOAD(4,1:n)*Basis(1:n))
      L  = MATMUL(LOAD(1:3,1:n), Basis(1:n))
      Bn = Bn + Detj * IP % S(t) * (Ln+SUM(L*Normal))
    END DO
!------------------------------------------------------------------------------
  END FUNCTION LocalFluxBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(  STIFF, FORCE, LOAD, GapLength, AirGapMu, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: LOAD(:,:), GapLength(:), AirGapMu(:)
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Normal(3)
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), localGapLength, muAir, muVacuum
    LOGICAL :: Stat
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )

    EdgeBasisDegree = 1
    IF (SecondOrder) EdgeBasisDegree = 2

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp

    muVacuum = 4 * PI * 1d-7

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n
       IF ( PiolaVersion ) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
               DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
               BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
       ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )

          CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
       END IF

       localGapLength  = SUM(Basis(1:n) * GapLength(1:n))
       muAir  = SUM(Basis(1:n) * AirGapMu(1:n))
 
       DO i = 1,nd-np
         p = i+np
         DO j = 1,nd-np
           q = j+np
           STIFF(p,q) = STIFF(p,q) + localGapLength / (muAir*muVacuum) * &
              SUM(RotWBasis(i,:)*RotWBasis(j,:))*detJ*IP%s(t)
         END DO
       END DO  
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixAirGapBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE DirichletAfromB()
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: s,p(3),q(3),cx(3),r,xmin,ymin,zmin,xmax,ymax,zmax
    TYPE(ListMatrixEntry_t), POINTER :: Ltmp
    TYPE(Matrix_t), POINTER :: Smat
    TYPE(Nodes_t),SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: BC

    LOGICAL :: Found, Found1,Found2,Found3,L1,L2,L3
    INTEGER :: i,j,k,l,m,t,ii,Faces,n,nd,Active,je1,je2,pe1,pe2

    TYPE(Element_t), POINTER :: Element, Edge, Edge1
    REAL(KIND=dp), ALLOCATABLE :: Bn(:)
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER, ALLOCATABLE :: dMap(:),FaceMap(:)
    LOGICAL, ALLOCATABLE :: FluxBoundaryEdge(:), CycleEdges(:), UsedFaces(:)
!------------------------------------------------------------------------------
    ALLOCATE(FluxBoundaryEdge(Mesh % NumberOFEdges)); FluxBoundaryEdge=.FALSE.

    Active = GetNOFBoundaryElements()
    DO t=1,Active
       Element => GetBoundaryElement(t)

       IF (.NOT. ActiveBoundaryElement()) CYCLE
       IF ( GetElementFamily()==1 ) CYCLE
       BC=>GetBC()
       IF (.NOT. ASSOCIATED(BC) ) CYCLE

       Found = ListCheckPrefix(BC,'Magnetic Flux Density')

       IF ( Found ) THEN
         SELECT CASE(GetElementFamily())
         CASE(2)
           CYCLE !what would it mean in 2D,at least with only B_z solved?
         CASE(3,4)
           k = GetBoundaryFaceIndex(Element); Element => Mesh % Faces(k)
         END SELECT
         FluxBoundaryEdge(Element % EdgeIndexes)=.TRUE.
       END IF
    END DO

    FluxCount = COUNT(FluxBoundaryEdge)
    IF ( FluxCount==0 ) THEN
      DEALLOCATE(FluxBoundaryEdge); RETURN
    END IF

    IF (.NOT.ALLOCATED(FluxMap) ) ALLOCATE(FluxMap(FluxCount))
    FluxCount = 0
    FluxMap   = 0
    DO i=1,Mesh % NumberOfEdges
      IF ( FluxBoundaryEdge(i) ) THEN
        FluxCount = FluxCount+1
        FluxMap(FluxCount) = i
      END IF
    END DO
    DEALLOCATE(FluxBoundaryEdge)

    DO i=1,FluxCount
      Edge => Mesh % Edges(FluxMap(i))
      Edge % BoundaryInfo % Left => NULL()
      Edge % BoundaryInfo % Right => NULL()
    END DO

    ALLOCATE(FaceMap(Mesh % NumberOfFaces)); FaceMap=0
    Faces = 0
    DO t=1,Active
      Element => GetBoundaryElement(t)

      IF (.NOT. ActiveBoundaryElement()) CYCLE
      IF ( GetElementFamily()==1 ) CYCLE
      BC=>GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      Found = ListCheckPrefix(BC,'Magnetic Flux Density')
      IF ( .NOT. Found ) CYCLE

      k = GetBoundaryFaceIndex(Element); Element=>Mesh % Faces(k)
      Faces = Faces+1
      FaceMap(k) = Faces

      DO i=1,Element % TYPE % NumberOfNodes
        Edge => Mesh % Edges(Element % EdgeIndexes(i))
        IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Left)) THEN
           Edge % BoundaryInfo % Left => Element
        ELSE IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Right)) THEN
           Edge % BoundaryInfo % Right => Element
        END IF
      END DO
    END DO


    ! Make gauge tree for the boundary:
    ! ---------------------------------
    CALL GaugeTreeFluxBC()

    WRITE(Message,*) 'Boundary tree edges: ', &
      TRIM(i2s(COUNT(TreeEdges(FluxMap)))),   &
             ' of total: ',TRIM(i2s(FluxCount))
    CALL Info('WhitneyAVSolver: ', Message, Level=5)

    ! Get (B,n) for BC faces:
    ! -----------------------
    ALLOCATE(Bn(Faces))
    DO t=1,Active
      Element => GetBoundaryElement(t)

      IF (.NOT.ActiveBoundaryElement()) CYCLE
      IF ( GetElementFamily()==1 ) CYCLE
      BC=>GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      n  = GetElementNOFNodes(Element)
      CALL GetRealVector(BC,Load(1:3,1:n),'Magnetic Flux Density',Found1)
      LOAD(4,1:n) = GetReal(BC,'Magnetic Flux Density {n}',Found)

      IF (Found.OR.Found1) THEN
        k = GetBoundaryFaceIndex(Element)
        Element => Mesh % Faces(k)
        nd = GetElementNOFDOFs(Element)
        Bn(FaceMap(k))=LocalFluxBC(LOAD,Element,n,nd)
      END IF
    END DO

    !
    ! Calculate value for free edges using the fundamental loop basis
    ! generated by GaugeTreeFluxBC():
    ! ---------------------------------------------------------------
    ALLOCATE(CycleEdges(Mesh % NumberOFEdges), UsedFaces(Faces))
    CycleEdges = .FALSE.
    ALLOCATE(dMap(MAXVAL(BasicCycles(:) % Degree)))

    Smat => GetMatrix()
    DO i=1,SIZE(BasicCycles)
      IF (BasicCycles(i) % Degree<=0 ) CYCLE

      ! 
      ! Extract loop edge indices: 
      ! --------------------------
      j = 0
      Ltmp => BasicCycles(i) % Head
      DO WHILE(ASSOCIATED(Ltmp))
        j = j + 1
        dMap(j) = Ltmp % Index; Ltmp => Ltmp % Next
      END DO
      IF ( j<= 0 ) CYCLE

      !
      ! Orient edges to form a polygonal path:
      ! --------------------------------------
      Edge  => Mesh % Edges(dMap(j))
      Edge1 => Mesh % Edges(dMap(j-1))
      IF ( ANY(Edge % NodeIndexes(1)==Edge1 % NodeIndexes) ) THEN
        l = Edge % NodeIndexes(1)
        Edge % NodeIndexes(1) = Edge % NodeIndexes(2)
        Edge % NodeIndexes(2) = l
      END IF
 
      DO k=j-1,1,-1
        Edge1 => Mesh % Edges(dMap(k))
        IF (Edge % NodeIndexes(2)==Edge1 % NodeIndexes(2)) THEN
          l = Edge1 % NodeIndexes(1)
          Edge1 % NodeIndexes(1) = Edge1 % NodeIndexes(2)
          Edge1 % NodeIndexes(2) = l
        END IF
        Edge => Edge1
      END DO

      !
      ! Try to find which way is inside...
      ! ----------------------------------
      Edge => Mesh % Edges(dMap(j))
      Element => Edge % BoundaryInfo % Left
      IF ( j==3 ) THEN
        m = 0
        DO k=1,3
          DO l=1,3
            IF (dMap(l)==Element % EdgeIndexes(k)) m=m+1
          END DO
        END DO
        L1 = m==3
        IF ( .NOT. L1 ) Element=>Edge % BoundaryInfo % Right
        S = Bn(FaceMap(Element % ElementIndex))
      ELSE
        ! If not a triangle, try a (planar) polygonal test. This
        ! will fail for general 3D paths. We'll spot the failure
        ! later by trial and error...Might be preferable to skip
        ! this altogether? Dunno....
        ! ------------------------------------------------------
        xmin=HUGE(xmin); xmax=-HUGE(xmax);
        ymin=HUGE(ymin); ymax=-HUGE(ymax);
        zmin=HUGE(zmin); zmax=-HUGE(zmax);
        DO k=1,j
          Edge1 => Mesh % Edges(dMap(k))
          DO l=1,2
            m = Edge1 % NodeIndexes(l)
            xmin = MIN(xmin,Mesh % Nodes % x(m))
            ymin = MIN(ymin,Mesh % Nodes % y(m))
            zmin = MIN(zmin,Mesh % Nodes % z(m))

            xmax = MAX(xmax,Mesh % Nodes % x(m))
            ymax = MAX(ymax,Mesh % Nodes % y(m))
            zmax = MAX(zmax,Mesh % Nodes % z(m))
          END DO
        END DO
        L1 = xmax-xmin > ymax-ymin
        L2 = xmax-xmin > zmax-zmin
        L3 = ymax-ymin > zmax-zmin
        IF ( l1 ) THEN
          l=1
          IF ( l3 ) THEN
            m=2; n=3
          ELSE
            m=3; n=2
          END IF
        ELSE
          IF ( l2 ) THEN
            l=1; m=2; n=3
          ELSE
            l=3; m=1; n=2
          END IF
        END IF
        cx(l) = SUM(Mesh % Nodes % x(Element % NodeIndexes))/3._dp
        cx(m) = SUM(Mesh % Nodes % y(Element % NodeIndexes))/3._dp
        cx(n) = SUM(Mesh % Nodes % z(Element % NodeIndexes))/3._dp

        L1 = .FALSE.
        DO k=j,1,-1
          Edge1 => Mesh % Edges(dMap(k))
          je1 = Edge1 % NodeIndexes(1)
          je2 = Edge1 % NodeIndexes(2)
          p(l) = Mesh % Nodes % x(je1)
          p(m) = Mesh % Nodes % y(je1)
          p(n) = Mesh % Nodes % z(je1)

          q(l) = Mesh % Nodes % x(je2)
          q(m) = Mesh % Nodes % y(je2)
          q(n) = Mesh % Nodes % z(je2)

          IF ((q(2)>cx(2)).NEQV.(p(2)>cx(2))) THEN
            IF (cx(1)<(p(1)-q(1))*(cx(2)-q(2))/(p(2)-q(2))+q(1)) L1=.NOT.L1
          END IF
        END DO
        IF (.NOT.L1) THEN
          IF (ASSOCIATED(Edge % BoundaryInfo % Right)) &
            Element=>Edge % BoundaryInfo % Right
        END IF

        ! Compute integral of (B,n) inside the cycle path
        ! -----------------------------------------------
        CycleEdges(dMap(1:j))=.TRUE.
        DO m=1,2
          S=0; UsedFaces=.FALSE.;
          IF( FloodFill(Element,CycleEdges, &
                       FaceMap,UsedFaces,Bn,S) )EXIT

          ! the in/out guess was wrong, try the other way:
          ! ----------------------------------------------
          IF (ASSOCIATED(Edge % BoundaryInfo % Right,Element)) THEN
            Element => Edge % BoundaryInfo % Left
          ELSE
            Element => Edge % BoundaryInfo % Right
          END IF
        END DO
        CycleEdges(dMap(1:j))=.FALSE.
      END IF

      !
      ! Orient edge to parent triangle...
      ! ---------------------------------
      je1 = Edge % NodeIndexes(1)
      je2 = Edge % NodeIndexes(2)
      EdgeMap => GetEdgeMap(GetElementFamily(Element))
      DO t=1,Element % TYPE % NumberOfEdges
        pe1 = Element % NodeIndexes(EdgeMap(t,1))
        pe2 = Element % NodeIndexes(EdgeMap(t,2))
        IF (pe1==je1.AND.pe2==je2 .OR. pe1==je2.AND.pe2==je1) EXIT
      END DO
      IF ( pe1/=je1 ) S=-S

      !
      ! ...because we now know how to orient against outward normal:
      ! ------------------------------------------------------------
      CALL GetElementNodes(Nodes,Element)
      p = NormalVector(Element,Nodes,0._dp,0._dp)
      q = NormalVector(Element,Nodes,0._dp,0._dp,.TRUE.)
      IF ( SUM(p*q)<0 ) S=-S

      !
      ! Check whether some edges in the path have nonzero values,
      ! if so, substract from integral:
      ! ---------------------------------------------------------
      DO k=j-1,1,-1
        l = Perm(dMap(k)+nNodes)
        R = Smat % RHS(l)/Smat % Values(Smat % Diag(l))
        IF ( R==0 ) CYCLE

        Edge1 => Mesh % Edges(dMap(k))
        pe1=Edge1 % NodeIndexes(1)
        pe2=Edge1 % NodeIndexes(2)
        IF ( pe2<pe1 ) R=-R; S=S-R
      END DO

      !
      ! ...and finally we should have the edge value:
      ! ---------------------------------------------
      IF ( je2<je1 ) S=-S
      CALL SetDOFtoValue(Solver,dMap(j),S)
    END DO
    DEALLOCATE(dMap, CycleEdges, FaceMap, UsedFaces, Bn)
    CALL List_FreeMatrix(SIZE(BasicCycles), BasicCycles)
!------------------------------------------------------------------------------
  END SUBROUTINE DirichletAfromB 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE FUNCTION FloodFill(Element,CycleEdges, &
          FaceMap,UsedFaces,Bn,CycleSum) RESULT(Found)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: e, Element
    REAL(KIND=dp) :: CycleSum, Bn(:)
    INTEGER :: i,j,n, FaceMap(:)
    LOGICAL :: CycleEdges(:), UsedFaces(:), Found

    Found=.FALSE.
    IF (.NOT.ASSOCIATED(Element)) RETURN

    n=FaceMap(Element % ElementIndex)
    IF (UsedFaces(n)) THEN
      Found=.TRUE.; RETURN
    END IF
    UsedFaces(n)=.TRUE.
    CycleSum = CycleSum+Bn(n)

    DO i=1,Element % TYPE % NumberOfEdges
      j = Element % EdgeIndexes(i)
      IF ( CycleEdges(j) ) CYCLE

      e => Mesh % Edges(j) % BoundaryInfo % Right
      IF(.NOT.FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum)) RETURN

      e => Mesh % Edges(j) % BoundaryInfo % Left
      IF(.NOT.FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum)) RETURN
    END DO
    Found=.TRUE.; RETURN
!------------------------------------------------------------------------------
  END FUNCTION FloodFill
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE WhitneyAVSolver
!------------------------------------------------------------------------------
