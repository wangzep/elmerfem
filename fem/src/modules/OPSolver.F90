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
    !TYPE(ValueList_t), POINTER :: Constants, Equation
    REAL(KIND=dp) :: Norm
    REAL(KIND=dp) :: beta
    INTEGER :: n, nb, nd, t, active
    INTEGER :: iter, maxiter, newton_interations
    LOGICAL :: found, newton = .FALSE., beta_solve=.FALSE.
    !------------------------------------------------------------------------------

    ! * Again, this is boiler plate code for Elmer Solvers

    CALL DefaultStart()

    !Get the beta parameter. Done outside the loop to save on processing.
    !We just need it the first time.

    beta_solve=GetLogical(Model%Constants, 'Compute Beta', found)
    CALL FoundCheck(found, 'Compute Beta', 'warn')

    IF (beta_solve) THEN
        beta=BetaCalc()
    ELSE
        beta=GetConstReal(Model%Constants,'beta',found)
        CALL FoundCheck(found, 'beta', 'fatal')
    END IF

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
                CALL LocalMatrix(Element, n, nd, beta)
            ELSE

                ! Assemble the matrix
                CALL LocalMatrix(  Element, n, nd+nb, beta )
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
    SUBROUTINE LocalMatrix( Element, n, nd, betain )
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        REAL(KIND=dp), INTENT(IN) :: betain
        !------------------------------------------------------------------------------
        REAL(KIND=dp) ::  laser_direction(3,n),directionterm(3), Weight, nonlinearterm
        REAL(KIND=dp) :: alkali_density(n), spin_destruction_rate(n),  &
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
        CALL FoundCheck(found, 'alkali density', 'fatal')

        spin_destruction_rate(1:n)=GetReal(Material,'spin destruction rate',found)
        CALL FoundCheck(found, 'spin destruction', 'fatal')

        beta(1:n)=betain

        laser_direction = 0._dp
        DO i=1,dim
            laser_direction(i,1:n)=GetReal(Material,&
                'laser direction '//TRIM(I2S(i)),found)
                CALL FoundCheck(found, 'laser direction', 'fatal')
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
                CALL Warn('OPSolver', outputstring)
            ELSE
                CALL Fatal('OPSolver', outputstring)
            END IF
        END IF

    !-------------------------------------------------------------------------
    END SUBROUTINE FoundCheck
    !-------------------------------------------------------------------------

    FUNCTION BetaCalc(Model, n, x) RESULT(beta)
        USE DefUtils
        IMPLICIT NONE
        TYPE(Model_t), OPTIONAL :: Model
        REAL(KIND=dp), INTENT(IN), OPTIONAL ::n, x(3)
        REAL(KIND=dp) :: beta

        !------------------------------------------------------------------------!
        REAL(KIND=dp) :: RPi

        REAL(KIND=dp) :: TWO, LOG2
        REAL(KIND=dp) :: speed_of_light, plank_constant, oscillator_strength,&
            electron_radius

        REAL(KIND=dp) :: alkali_wavelength,laser_wavelength, laser_linewidth, alkali_freq_width

        REAL(KIND=dp) :: laser_frequency, alkali_frequency, laser_freq_width, absorb_laser_ratio,&
            frequency_shift

        REAL(KIND=dp) :: w_prime, w_dprime, w_input_real, w_input_imaginary

        TYPE(ValueList_t), POINTER :: Constants

        LOGICAL :: FLAG, found, found1, found2, found3, found4

        Constants => GetConstants()

        !------------------------------------------------------------------------!
        !Declare constants-------------------------------------------------------
        TWO = 2.0D0
        LOG2 = LOG(TWO)
        RPi=4.D0*DATAN(1.D0)
        !-------------------------------------------------------------------------

        !Get the information about the rubidium and the laser from the SIF
        !file---------------------------------------------------------------------

        alkali_wavelength = GetConstReal(Constants,'alkali wavelength',found)
        CALL FoundCheck(found, 'alkali wavelength', 'fatal')

        laser_wavelength = GetConstReal(Constants,'laser wavelength',found)
        CALL FoundCheck(found, 'laser wavelength', 'fatal')

        laser_linewidth = GetConstReal(Constants,'laser line width',found)
        CALL FoundCheck(found, 'laser line width', 'fatal')

        alkali_freq_width = GetConstReal(Constants,'alkali frequency width',found)
        CALL FoundCheck(found, 'alkali frequency width', 'fatal')

        oscillator_strength = GetConstReal(Constants,'oscillator strength', found1)
        electron_radius = GetConstReal(Constants, 'electron radius', found2)
        plank_constant = GetConstReal(Constants, 'Planks Constant',found3)
        speed_of_light = GetConstReal(Constants,'speed of light',found4)

        IF (.NOT. (found1 .AND. found2 .AND. found3 .AND. found4)) THEN
            oscillator_strength=1.0D0/3.0D0
            plank_constant = 6.62607004D-34
            speed_of_light = 299792458.0D0
            electron_radius = 2.8179403267e-15

            CALL Warn('BetaCalc',&
                'One or more of the constants are not listed in the SIF. Using default values SI units.')
        END IF

        !For testing---------------------------------!!!!!!!!!!!!!!!!!!!!!!!
        !rubidium_wavelength = 800e-9
        !rubidium_freq_width = 126.65e9

        !laser_wavelength = 800e-9
        !laser_linewidth = 2e-9

        !oscillator_strength = 1.0/3.0


        !-------------------------------------------------------------------------

        !Define spectral overlap function (eq. A7, A10 of Fink's paper)
        alkali_frequency = speed_of_light/alkali_wavelength

        !The second equation here uses the dispresion relationship
        laser_frequency = speed_of_light/laser_wavelength

        laser_freq_width = (speed_of_light*laser_linewidth)/(laser_wavelength)**2

        absorb_laser_ratio = alkali_freq_width/laser_freq_width

        frequency_shift = 2*(laser_frequency - alkali_frequency)/laser_freq_width

        !So... A8 is the Faddeeva function written in a funny way. Note, that in order to recover
        !The Faddeeva function normally, we need to divide the arguement of erfc in A8 by -i.
        !When you do that you get that the real part is -sqrt(ln(2))*s and the imaginary part
        !is sqrt(ln(2))*r. This is very confusing, and I'm going to make a note in my lab note-
        !book regarding the maninpulation.

        !So I don't think the below assignments are wrong!
        w_input_real = -SQRT(LOG2)*(frequency_shift)
        w_input_imaginary = SQRT(LOG2)*(absorb_laser_ratio)


        CALL WOFZ(w_input_real,w_input_imaginary,w_prime,w_dprime,FLAG)

        beta = 2*DSQRT(RPi*LOG2)*(electron_radius*oscillator_strength*&
            laser_wavelength**2*w_prime)/laser_linewidth
            !I added a cos term to assure that the laser optical pumping rate goes to zero at the boundaries
            !of the laser beam area. I was getting some computational artifacts that appeared to be due to the
            !discontinity at the border.
            !The factor of 3.14/2 is to account for the difference in area under the curve of a cos and a Heaviside
            !function.
            !initopt = 3.14/2*beta*power/(area*plank_constant*laser_frequency)*COS((3.14**(3/2)*SQRT(x(1)**2 + x(2)**2))/(SQRT(area)))


            RETURN

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
10              CONTINUE
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
11              CONTINUE

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

100         FLAG = .TRUE.
            RETURN

        END SUBROUTINE WOFZ


    !------------------------------------------------------------------------------
    END SUBROUTINE OPSolver
!------------------------------------------------------------------------------
