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
    !Use the element description model to try to construct a limit for non-linear term
    USE elementdescription

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
    LOGICAL :: found, newton = .FALSE., beta_solve=.FALSE.,&
        laser_power_use=.FALSE.
    TYPE(ValueList_t), POINTER :: SolverParams
    !------------------------------------------------------------------------------
    !Create a Flag structure that is somewhat easier to pass to LocalMatrix.
    !This will have all the flags for turning stuff on or off that can't be
    !contained in the matrix loop for various reasons.

    TYPE flags
        LOGICAL limitrate
        LOGICAL limitratewarning
        LOGICAL limitalkalidensitywarning
    END TYPE Flags

    !Define initial value in Flags, should all be false:
    TYPE(flags) :: modflags

    modflags % limitrate = .FALSE.
    modflags % limitratewarning = .FALSE.
    modflags % limitalkalidensitywarning = .FALSE.


    !------------------------------------------------------------------------------
    ! Factor to convert alkali density from kg/m^3 to part/m^3, I think this is only valid for rubidium
    !REAL, PARAMETER :: rb_density_conversion_factor = 7.043279D24

    ! * Again, this is boiler plate code for Elmer Solvers

    CALL DefaultStart()

    !Check to make sure that the typical OP variables in the SIF are within bounds
    CALL ArgumentCheck()

    !Get information from the solver part of the sif and set the rate limiter flag

    modflags % limitrate = GetLogical(Model% Solver % Values, 'Limit OP Rate', found)

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
                CALL LocalMatrix(Element, n, nd, beta, modflags)
            ELSE

                ! Assemble the matrix
                CALL LocalMatrix(  Element, n, nd+nb, beta, modflags)
            END IF

        END DO

        CALL DefaultFinishBulkAssembly()

        laser_power_use=GetLogical(Model%Constants, 'Compute Laser Power', found)
        CALL FoundCheck(found, 'Compute Laser Power', 'warn')


        !Calculate the inital optical pumping rate using laser power.

        IF (laser_power_use) THEN
            Active = GetNOFBoundaryElements()
            DO t=1,Active
                Element => GetBoundaryElement(t)
                IF(ActiveBoundaryElement()) THEN
                    n  = GetElementNOFNodes()
                    nd = GetElementNOFDOFs()
                    CALL CalcLaserPowerBC( Model, Element, n, nd,beta )
                END IF
            END DO
        END IF

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
    SUBROUTINE LocalMatrix( Element, n, nd, betain, modflags)
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        REAL(KIND=dp), INTENT(IN) :: betain
        !------------------------------------------------------------------------------
        REAL(KIND=dp) ::  laser_direction(3,n),directionterm(3), Weight, nonlinearterm
        REAL(KIND=dp) :: alkali_density(n), spin_destruction_rate(n),  &
            beta(n),previous_solution(n), temp(n), theta(n), elementdiam

        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ!,LoadAtIP
        REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)

        LOGICAL :: Stat,found

        LOGICAL :: convert_density=.FALSE., skewlight=.FALSE.

        INTEGER :: i,t,p,q,dim, ind
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(ValueList_t), POINTER :: BodyForce, Material, Constants
        TYPE(Nodes_t) :: Nodes
        TYPE(flags) :: modflags
        SAVE Nodes
        !------------------------------------------------------------------------------
        dim = CoordinateSystemDimension()

        CALL GetElementNodes( Nodes )
        MASS  = 0._dp
        STIFF = 0._dp
        FORCE = 0._dp
        LOAD = 0._dp
        theta(1:n) = 0._dp

        !BodyForce => GetBodyForce()
        !IF ( ASSOCIATED(BodyForce) ) &
        !    Load(1:n) = GetReal( BodyForce,'field source', found )

        !Get the previous solution
        CALL GetScalarLocalSolution(previous_solution, 'optrate', Element)

        DO ind = 1, n
            IF (previous_solution(ind) .LT. 0) THEN
                previous_solution(ind) = 0

                !If the rate is not limited then call fatal if the optrate drops
                !below zero.

                IF (.NOT. modflags % limitrate) THEN
                    CALL FATAL('OPSolver',&
                        'optrate is less than 0, this is not physically possible')
                END IF
            END IF
        END DO

        !Pointers to the SIF Constants and Material information
        Constants => GetConstants()
        Material => GetMaterial()

        !Let's get the relevent information

        !Alkali Density-------------------------------------------------

        alkali_density(1:n)=GetReal(Material,'alkali density',found)
        CALL FoundCheck(found, 'alkali density', 'fatal')

        !Check and correct the alaklidensity if needed
        DO ind = 1, n
            IF (alkali_density(ind) .LT. 0) THEN

                alkali_density(ind) = 0
                !Show a warning, but only the first time it is thrown
                IF (.NOT. modflags % limitalkalidensitywarning) THEN

                    CALL WARN('OPSolver',&
                        'alkali density is less than 0, this is not physically possible')
                    !Change the value of the flag so that the warning is not thrown
                    !the next time through.
                    modflags % limitalkalidensitywarning = .TRUE.

                END IF

            END IF

        END DO

        !        convert_density=GetLogical(Material, 'Convert Density')
        !        CALL FoundCheck(found, 'Convert Density', 'warn')

        !        IF (convert_density) THEN
        !            alkali_density = rb_density_conversion_factor*alkali_density
        !        END IF



        ! Spin destruction Rate--------------------------------------

        spin_destruction_rate(1:n)=GetReal(Material,'spin destruction rate',found)
        CALL FoundCheck(found, 'spin destruction', 'fatal')

        DO ind = 1, n
            IF (spin_destruction_rate(ind) .LT. 0) THEN
                CALL FATAL('OPSolver',&
                    'spin destruction rate is less than 0, this is not physically possible')
            END IF
        END DO

        !beta
        beta(1:n)=betain

        DO ind = 1, n
            IF (beta(ind) .LT. 0) THEN
                CALL FATAL('OPSolver',&
                    'beta is less than 0, this is not physically possible')
            END IF
        END DO

        !Direction of the laser beam
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


            !Put together the nonlinear term----------------------------

            !Do we use skew light?
            skewlight=GetLogical(Material, 'Use Skew Light', found)
            CALL FoundCheck(found, 'Use Skew Light', 'warn')


            IF (.NOT. skewlight) THEN
                temp=spin_destruction_rate+previous_solution
                temp=previous_solution/temp


            ELSE
                            !Skew light term
                theta(1:n)=GetReal(Material, 'theta', found)
                CALL FoundCheck(found, 'theta', 'fatal')

                DO ind = 1, n
                    IF (theta(ind) .LT. 0) THEN
                        CALL FATAL('OPSolver',&
                            'theta is less than 0, this is not physically possible')
                    END IF
                    !Check if theta is greater than pi/2
                    IF (theta(ind) .GT. 1.57079632679) THEN
                        CALL FATAL('OPSolver',&
                            'theta is greater than pi/2, this is not physically possible')
                    END IF
                END DO

                temp=spin_destruction_rate+previous_solution
                temp=(previous_solution*(1-SIN(theta)**2))/temp

            END IF


            temp=1-temp

            !Assure that the temp term is positive and if not
            !Make it a small positive number

            DO ind = 1, n
                IF (temp(ind) .LT. 1D-8) THEN
                    temp(ind)=1D-8
                    CALL WARN('OPSolver',&
                        'Derivative of oprate is greater than 0, reset to 0')
                END IF
            END DO

            temp=alkali_density*temp
            temp=beta*temp

            !If the rate is limited, check to see if temp*elementdiam < 1. This would be an indicator that the solution will be
            !than zero and that temp is too big for the geometry. Essentially, this is a hack
            !to account for rapidly changing optical pumping rates that change more quickly
            !than our mesh size can account for, thus delivering really unrealistic solutions
            !In and ideal world, we would apply mesh refinements, but Elmer doesn't support that
            !really well.


            IF (modflags % limitrate) THEN
                elementdiam = ElementDiameter(Element,Nodes)
                DO ind = 1, n
                    IF (temp(ind)*elementdiam .GT. 1) THEN
                        temp = 1/elementdiam

                        !Show warning, but only the first time it is thrown
                        IF (.NOT. modflags % limitratewarning) THEN
                            CALL WARN('OPSolver',&
                                'The nonlinear term exceeds estmiated possible value. Resetting value')
                            !Change flag so that the warning is not shown the next time through
                            modflags % limitratewarning = .TRUE.
                        END IF

                    END IF
                END DO
            END IF


            !Construct the nonlinear term for the local matrix.

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
    SUBROUTINE CalcLaserPowerBC( Model, Element, n, nd , beta)
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        !        REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
        !        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
        !        REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
        LOGICAL :: Stat,found, found1, found2
        INTEGER :: dim
        !        INTEGER :: i,t,p,q,dim
        !        TYPE(GaussIntegrationPoints_t) :: IP

        REAL(KIND=dp) :: area, laserpower, beta, oprate
        TYPE(ValueList_t), POINTER :: BC, Constants
        TYPE(Model_t) :: Model

        REAL(KIND=dp) :: plank_constant, speed_of_light, laser_frequency,&
            laser_wavelength

        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------
        BC => GetBC()
        IF (.NOT.ASSOCIATED(BC) ) RETURN

        !Constants=GetConstants()

        laser_wavelength = GetConstReal(Model%Constants,'laser wavelength',found)
        CALL FoundCheck(found, 'laser wavelength', 'fatal')
        IF (laser_wavelength .LT. 0) THEN
            CALL Fatal('CalcLasePowerBC',&
                'laser wavelength is less than 0, that is not physically possible')
        END IF



        dim = CoordinateSystemDimension()

        plank_constant = GetConstReal(Model%Constants, 'planks constant',found1)
        IF (plank_constant .lt. 0) THEN
            CALL Fatal('CalcLaserPowerBC',&
                'planks constant is less than 0, that is not physically possible')
        END IF

        speed_of_light = GetConstReal(Model%Constants,'speed of light',found2)
        IF (speed_of_light .lt. 0) THEN
            CALL Fatal('CalcLaserPowerBC',&
                'speed of light is less than 0, that is not physically possible')
        END IF


        IF (.NOT. (found1 .AND. found2)) THEN

            plank_constant = 6.62607004D-34
            speed_of_light = 299792458.0D0

            CALL Warn('BetaCalc',&
                'One or more of the constants are not listed in the SIF. Using default values SI units.')
        END IF

        area=GetConstReal(Model%Constants, 'Laser Spot Size', found)
        CALL FoundCheck(found, 'Laser Spot Size', 'fatal')
        IF (area .lt. 0) THEN
            CALL Fatal('CalcLaserPowerBC',&
                'laser spot size is less than 0, that is not physically possible')
        END IF

        laserpower=GetConstReal(BC, 'Laser Power', found)
        IF (.NOT. found) RETURN !If the BC isn't there for this element, we shouldn't set it
        IF (laserpower .lt. 0) THEN
            CALL Fatal('CalcLaserPowerBC',&
                'Laser Power is less than 0, that is not physically possible')
        END IF


        laser_frequency = speed_of_light/laser_wavelength

        !Add the optrate entry

        oprate=beta*laserpower/(area*plank_constant*laser_frequency)

        CALL ListAddConstReal(BC, 'optrate', oprate)

    !------------------------------------------------------------------------------
    END SUBROUTINE CalcLaserPowerBC
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    FUNCTION BetaCalc(Model, n, x) RESULT(beta)
        !------------------------------------------------------------------------------
        !Default function (no change of speed of light, plank's constant or electron
        !radius) gives answer in m^2. Default oscilator strength is for Rb = 1/3
        !Output value for:
        !alkali wavelength = Real 800.0D-9,
        !alkali frequency width = Real 126.65D9
        !alkali temperature = Real 433.15
        !laser wavelength = Real 800.0D-9
        !laser line width = Real 2.0D-9
        ! is 7.8480143564070032e-19 m^2.
        !------------------------------------------------------------------------------
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

        IF (alkali_wavelength .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'alakali wavelength is less that zero, this is not physically possible')
        END IF

        laser_wavelength = GetConstReal(Constants,'laser wavelength',found)
        CALL FoundCheck(found, 'laser wavelength', 'fatal')

        IF (laser_wavelength .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'laser wavelength is less that zero, this is not physically possible')
        END IF

        laser_linewidth = GetConstReal(Constants,'laser line width',found)
        CALL FoundCheck(found, 'laser line width', 'fatal')

        IF (laser_linewidth .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'laser linewidth is less that zero, this is not physically possible')
        END IF

        alkali_freq_width = GetConstReal(Constants,'alkali frequency width',found)
        CALL FoundCheck(found, 'alkali frequency width', 'fatal')

        IF (alkali_freq_width .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'alkali frequency width is less that zero, this is not physically possible')
        END IF

        oscillator_strength = GetConstReal(Constants,'oscillator strength', found1)

        IF (oscillator_strength .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'oscillator strength is less that zero, this is not physically possible')
        END IF

        electron_radius = GetConstReal(Constants, 'electron radius', found2)

        IF (electron_radius .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'electron radius is less that zero, this is not physically possible')
        END IF

        plank_constant = GetConstReal(Constants, 'planks constant',found3)

        IF (plank_constant .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'planks constant is less that zero, this is not physically possible')
        END IF

        speed_of_light = GetConstReal(Constants,'speed of light',found4)

        IF (speed_of_light .LT. 0) THEN
            CALL Fatal('BetaCalc',&
                'speed_of_light is less that zero, this is not physically possible')
        END IF

        IF (.NOT. (found1 .AND. found2 .AND. found3 .AND. found4)) THEN
            oscillator_strength=1.0D0/3.0D0
            plank_constant = 6.62607004D-34
            speed_of_light = 299792458.0D0
            electron_radius = 2.8179403267D-15

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

        laser_freq_width = (speed_of_light*laser_linewidth)/(laser_wavelength)**2.0D0

        absorb_laser_ratio = alkali_freq_width/laser_freq_width

        frequency_shift = 2.0D0*(laser_frequency - alkali_frequency)/laser_freq_width

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

        RETURN
    !------------------------------------------------------------------------------
    END FUNCTION
    !------------------------------------------------------------------------------

    !We need the FADDEEVA function from equation (A8). Thankfully, someone in the
        !interweb has coded a nice function to calculate it for me. The ACM has collected
        !aligorithms to do this very thing. The following code is copied line for line
        !from ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
        !THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
        !VOL. 16, NO. 1, PP. 47. - http://www.netlib.org/toms/68

        !      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
        !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
        !      VOL. 16, NO. 1, PP. 47.
    !--------------------------------------------------------------------------
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

    !--------------------------------------------------------------------------
    END SUBROUTINE WOFZ
        !--------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    !    FUNCTION GetAlkaliConcentration(Element, n) RESULT(alkali_density)
    !        USE DefUtils
    !        !--------------------------------------------------------------------------
    !        INTEGER, INTENT(IN) :: n
    !        TYPE(Element_t), POINTER, INTENT(IN) :: Element
    !        !--------------------------------------------------------------------------
    !        REAL(KIND=dp) :: alkali_density(n)
    !        TYPE(ValueList_t), POINTER :: Materials
    !        LOGICAL :: convert_density = .FALSE.
    !        INTEGER :: ind = 0
    !
    !        alkali_density = 0.0D0
    !
    !        Materials => GetMaterial()
    !
    !        CALL GetScalarLocalSolution(alkali_density, 'Concentration', Element)
    !
    !        !Check to make sure we actually found the solution
    !        DO ind = 1, n
    !
    !
    !            IF (alkali_density(ind) .EQ. 0) Call Fatal('GetAlkaliConcentration',&
    !                'Concentration variable not found. Check name of AdvectDiff variable.')
    !
    !        END DO
    !
    !        IF (convert_density) THEN
    !            alkali_density(1:n) = 7.043279D24*alkali_density
    !        END IF
    !
    !
    !
    !    END FUNCTION GetAlkaliConcentration



!------------------------------------------------------------------------------
END SUBROUTINE OPSolver
!------------------------------------------------------------------------------

!FUNCTION calculaterbmumdensitym(Model,n,Temp) RESULT(RbNumDensity_m)
!    !-------------------------------------------------------------------------
!    !Calculates Rb number density in m^-3 using Killian equation as presented
!    !in Fink et al. 2005.
!    !n_Rb=10^(9.55-4132/T)/kT
!    !where T is the temperature in Kelvin and k is Boltzman's constant.
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: model
!    INTEGER :: n
!    REAL(KIND=dp) :: Temp !Dummy Variable. Not actually used
!    REAL(KIND=dp) :: RbNumDensity_m, Temperature
!    LOGICAL :: found=.FALSE.
!    !------------------------------------------------------------------------
!    TYPE(ValueList_t), POINTER :: Materials
!    !-----------------------------------------------------------------------
!
!    Materials => GetConstants()
!
!    Temperature=GetConstReal(Materials, 'alkali temperature', found)
!    IF (.NOT. found) CALL Fatal('RbNumDensity',&
!        'Temperature not found')
!
!
!    RbNumDensity_m=(10**(9.55D0-4132.0D0/Temperature))/(1.380648521D-23*Temperature)
!
!
!END FUNCTION calculaterbmumdensitym

!FUNCTION CalculateSpinDestructionRate(Model,n,argument)&
!    RESULT(SpinDestrucionRate)
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: model
!    INTEGER :: n
!    REAL(KIND=dp) :: argument(3)
!    REAL(KIND=dp) :: Concentration, Temperature, Pressure, SpinDestrucionRate
!    !----------------------------------------------------------------------
!
!    REAL(KIND=dp) :: alkali_alkali_spin_destruction_rate, xe_spin_destruction_rate,&
!        he_spin_destruction_rate, n2_spin_destruction_rate, G1
!    REAL(KIND=dp) :: xe_fraction, n2_fraction, he_fraction, xe_numberdensity,&
!        he_numberdensity, n2_numberdensity,tot_numberdensity, ref_pressure
!    REAL(KIND=dp) :: loschmidt
!    !-----------------------------------------------------------------
!    INTEGER :: ind
!    !-----------------------------------------------------------------
!    LOGICAL :: convert_density=.FALSE., he_term_included=.FALSE.
!    LOGICAL :: n2_term_included=.FALSE., vanderWall_term_included=.FALSE.
!    LOGICAL :: xe_term_included=.FALSE., local_temperature_used = .FALSE.
!    LOGICAL :: local_pressure_used = .FALSE., found=.FALSE.
!    !--------------------------------------------------------
!    TYPE(ValueList_t), POINTER :: Materials, Constants
!    !--------------------------------------------------------------
!
!
!
!    REAL, PARAMETER :: rb_density_conversion_factor = 7.043279D24
!
!    SpinDestrucionRate = 0.0D0
!
!    Materials => GetMaterial()
!    Constants => GetConstants()
!
!    !!!!!!!!!!!!!!!!!!!Alkali-Alkali Spin Destruction Term!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    alkali_alkali_spin_destruction_rate=&
!        GetConstReal(Materials, 'alkali-alkali spin destruction rate', found)
!    CALL FoundCheck(found, 'alkali-alkali spin destruction rate','fatal')
!
!    ! Assign the Concentration
!    Concentration=argument(1)
!
!    !Check if we should convert the density (we probably do want to
!    convert_density=GetLogical(Materials, 'convert density', found)
!    CALL FoundCheck(found, 'convert density', 'warn')
!
!    !Convert the alkali density to particles/m^3
!    IF (convert_density) THEN
!        Concentration = rb_density_conversion_factor*Concentration
!    END IF
!
!    SpinDestrucionRate=alkali_alkali_spin_destruction_rate*Concentration
!
!    !Inclusion of other terms in the spin-destruction rate: He-Rb, N2-Rb,
!    !and Van der Walls terms
!
!    he_term_included=GetLogical(Materials,&
!        'Helium Spin Destruction Term Included',found)
!
!    n2_term_included = GetLogical(Materials,&
!        'Nitrogen Spin Destruction Term Included',found)
!
!    xe_term_included = GetLogical(Materials,&
!        'Xenon Spin Destruction Term Included',found)
!
!    vanderWall_term_included = GetLogical(Materials,&
!        'vander Wall Spin Destruction Term Included',found)
!
!    IF (he_term_included .OR. n2_term_included .OR. &
!        xe_term_included .OR. vanderWall_term_included) THEN
!
!        !Get the temperature
!
!        Temperature = argument(2)
!
!        IF (Temperature .EQ. 0) Call Fatal('GetSpinDestructionRate',&
!            'Temperature variable not found.')
!
!
!        !Get the pressure
!
!        Pressure = argument(3)
!
!        !Check to make sure we actually found the solution
!
!
!        !        IF (Pressure .EQ. 0) Call Fatal('GetSpinDestructionRate',&
!        !            'Pressure variable not found. Check name of N-S variable.')
!
!
!        !Get the gas fractions
!        xe_fraction=GetConstReal(Materials, 'xe fraction', found)
!        CALL FoundCheck(found, 'xe fraction', 'fatal')
!
!        n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
!        CALL FoundCheck(found, 'n2 fraction' , 'fatal')
!
!        he_fraction = GetConstReal(Materials, 'he fraction', found)
!        CALL FoundCheck(found, 'he fraction', 'fatal')
!
!
!        !Call fatal if the gas fractions don't add to 1
!
!        IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
!            CALL Fatal('GetSpinDestructionRate', &
!                'Gas fractions do not add to 1')
!        END IF
!
!        !Calculate pressure in amagats
!
!        !Get Loschmidt's number if defined in constants
!
!        loschmidt=GetConstReal(Constants, 'loschmidts constant', found)
!        !CALL FoundCheck(found, 'loschmidts constant' , 'warn')
!        IF (.NOT. found) THEN
!            loschmidt= 2.6867811D25
!        END IF
!
!        tot_numberdensity = ((Pressure)/101325.0D0)*(273.15D0/Temperature)*loschmidt
!
!        xe_numberdensity=tot_numberdensity*xe_fraction
!
!        n2_numberdensity=tot_numberdensity*n2_fraction
!
!        he_numberdensity=tot_numberdensity*he_fraction
!
!        !Implement the xenon contribution to spin destruction rate
!        IF (xe_term_included) THEN
!
!            xe_spin_destruction_rate = 0.0D0
!            xe_spin_destruction_rate = GetConstReal(Materials,&
!                'xe spin destruction rate',found)
!            CALL FoundCheck(found, 'xe spin destruction rate', 'fatal')
!
!            SpinDestrucionRate= SpinDestrucionRate+&
!                xe_spin_destruction_rate*xe_numberdensity
!
!        END IF
!
!        !Implement helium contribution to spin destruction rate
!        IF (he_term_included) THEN
!
!            he_spin_destruction_rate = 0.0D0
!            he_spin_destruction_rate = GetConstReal(Materials,&
!                'he spin destruction rate',found)
!            CALL FoundCheck(found, 'he spin destruction rate', 'fatal')
!
!            SpinDestrucionRate = SpinDestrucionRate+&
!                he_spin_destruction_rate*he_numberdensity
!        END IF
!
!        !Implement N2 contribution to spin destruction rate
!        IF (n2_term_included) THEN
!
!            n2_spin_destruction_rate = 0.0D0
!            n2_spin_destruction_rate = GetConstReal(Materials,&
!                'n2 spin destruction rate',found)
!            CALL FoundCheck(found, 'n2 spin destruction rate', 'fatal')
!
!            SpinDestrucionRate=SpinDestrucionRate+&
!                n2_spin_destruction_rate*n2_numberdensity
!        END IF
!
!        !Implement van der Walls contribution to spin destruction rate
!        !See Nelson's 2001 thesis for details.
!        IF (vanderWall_term_included) THEN
!
!
!            G1 = 0.0D0
!            G1 = GetConstReal(Materials,'short-very short transition density',found)
!            CALL FoundCheck(found, 'short-very short transition density', 'fatal')
!
!            SpinDestrucionRate=SpinDestrucionRate+&
!                (0.385D0+0.642D0*1.0D0/(1.0D0+G1/tot_numberdensity))&
!                *6469.0D0/(xe_fraction+1.1D0*n2_fraction+3.2D0*he_fraction)
!        END IF
!    END IF
!
!END FUNCTION CalculateSpinDestructionRate
!
!FUNCTION CalculateRbPol(Model,n,argument) RESULT(RbPol)
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: argument(4), sdargument(3)
!    REAL(KIND=dp) :: spindestrucionrate, optrate
!    REAL(KIND=dp) :: RbPol
!    REAL(KIND=dp) :: CalculateSpinDestructionRate
!    !--------------------------------------------------------------------------------
!
!    optrate=argument(1)
!
!    sdargument= (/argument(2),argument(3),argument(4)/)
!
!    spindestrucionrate=CalculateSpinDestructionRate(Model, n, sdargument)
!
!    RbPol=optrate/(optrate+spindestrucionrate)
!
!END FUNCTION CalculateRbPol

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

SUBROUTINE ArgumentCheck()

    !Function that just checks that the SIF are positive
    !Tried to just look at terms that are related to OPSolver
    USE defutils

    IMPLICIT NONE
    REAL(KIND=dp) :: alkaliwavelength, alkalifreqwidth, laserwavelength,&
        laserlinewidth, laserspotsize
    TYPE(ValueList_t), POINTER :: Constants
    LOGICAL :: found
    !---------------------------------------------------------------

    Constants => GetConstants()

    alkaliwavelength = GetConstReal(Constants, 'alkali wavelength', found)
    IF (found) THEN
        IF (alkaliwavelength .lt. 0) THEN
            CALL Fatal('OPUtil', 'Alakli wavelength is less than 0')
        END IF
    END IF

    alkalifreqwidth = GetConstReal(Constants, 'alkali frequencey width', found)
    IF (found) THEN
        IF (alkalifreqwidth .lt. 0) THEN
            CALL Fatal('OPUtil', 'Alakli wavelength is less than 0')
        END IF
    END IF

    laserwavelength = GetConstReal(Constants, 'laser wavelength', found)
    IF (found) THEN
        IF (laserwavelength .lt. 0) THEN
            CALL Fatal('OPUtil', 'Laser wave length is less than 0')
        END IF
    END IF

    laserspotsize = GetConstReal(Constants, 'Laser Spot Size', found)
    IF (found) THEN
        IF (laserspotsize .lt. 0) THEN
            CALL Fatal('OPUtil', 'Laser spot size is less than 0')
        END IF
    END IF


END SUBROUTINE ArgumentCheck

!--------------------------------------------------------------
!From here down will need to moved to the eventual new XePol Solver
!FUNCTION CalculateSpinExchangeRate(Model,n,Argument)&
!    RESULT(SpinExchangeRate)
!
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: Argument
!    REAL(KIND=dp) :: Concentration
!    REAL(KIND=dp) :: SpinExchangeRate
!    REAL(KIND=dp) :: binaryexchangerate
!    TYPE(ValueList_t), POINTER :: Materials
!    LOGICAL :: found
!    !-----------------------------------------------------------
!
!    !Eventaully I might import more terms, so this is anticpating that eventuallity
!    Concentration=Argument
!
!    Materials=GetMaterial()
!
!    binaryexchangerate=GetConstReal(Materials, 'binary exchange rate', found)
!    CALL FoundCheck(found, 'binaryexchangerate', 'fatal')
!
!    SpinExchangeRate=binaryexchangerate*Concentration
!
!END FUNCTION CalculateSpinExchangeRate
!
!FUNCTION CalculateSpinRelaxationRate(Model,n,Argument)&
!    RESULT(SpinRelaxationRate)
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: Argument(2)
!    REAL(KIND=dp) :: Pressure,Temperature
!    REAL(KIND=dp) :: he_fraction=0, xe_fraction=0, n2_fraction=0
!    REAL(KIND=dp) :: SpinRelaxationRate
!    REAL(kind=dp) :: binary_term=0, vdWterm=0
!    TYPE(ValueList_t), POINTER :: Materials
!    LOGICAL :: found
!    !------------------------------------------------------------------------------
!
!    Materials=>GetMaterial()
!
!    Pressure=Argument(1)
!    Temperature=Argument(2)
!
!    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
!    CALL FoundCheck(found, 'xe fraction', 'fatal')
!
!    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
!    CALL FoundCheck(found, 'n2 fraction' , 'fatal')
!
!    he_fraction = GetConstReal(Materials, 'he fraction', found)
!    CALL FoundCheck(found, 'he fraction', 'fatal')
!
!    binary_term=(5D-6)*xe_fraction*((Pressure)/101325)*(273.15/Temperature)
!
!    vdWterm=6.72D-5*(1/(1+0.25*he_fraction/xe_fraction+1.05*n2_fraction/xe_fraction))
!
!    SpinRelaxationRate=binary_term+vdWterm
!
!END FUNCTION CalculateSpinRelaxationRate
!
!REAL(KIND=selected_real_kind(12)) FUNCTION CalculateXenonDiffusion(Model,n,Argument) RESULT(D_Xe)
!    !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
!    !I need to go back and document this better (probably when I revamp the XePol
!    !solver.
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: Argument(2)
!    !REAL(KIND=dp) :: D_Xe
!    !-------------------------------------------------------
!    REAL(KIND=dp) :: Pressure,Temperature
!    REAL(KIND=dp) :: pressure_atm=0
!
!    REAL(KIND=dp) :: mixKesp=0, tprime=0, omega=0, sigma=0, Mtot=0
!    !------------------------------------------------------------
!    !I believe these are all gotten from Lightfoot, sans the masses
!    REAL(Kind=dp), PARAMETER :: A=1.858D-7, sigmaHe=2.576, sigmaXe=4.009,&
!        KespHe=10.02, KespXe=234.7, massXe=131.29, massHe=4.003
!
!    !------------------------------------------------------------
!
!    !Getting assignments
!    Pressure=Argument(1)
!    Temperature=Argument(2)
!
!    !Convert to atm
!    pressure_atm=(Pressure)/101325
!
!    !Actually doing the calcuation.
!    mixKesp=sqrt(KespHe*KespXe)
!    tprime = Temperature/mixKesp
!
!    omega = (1.06036/tprime**(0.15610))+(0.19300/exp(0.47635*tprime))+&
!        (1.03587/exp(1.52296*tprime))+(1.76474/exp(3.89411*tprime))
!
!    sigma = 0.5*(sigmaHe+sigmaXe)
!
!    Mtot = sqrt(1/massHe+1/massXe)
!
!    D_Xe = (A*Temperature**(3/2)*Mtot)/(pressure_atm*omega*sigma**2)
!
!
!
!END FUNCTION CalculateXenonDiffusion
!
!FUNCTION CalculateDecayRate(Model,n,argument) RESULT(decayrate)
!
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: argument(2)
!    REAL(KIND=dp) :: pressure, temperature
!    REAL(KIND=dp) :: decayrate
!    !--------------------------------------------
!    REAL(KIND=dp) :: cell_radius=0, cellT1=0,&
!        pressureT1=0, temperatureT1=0, T1vals(2)
!    REAL(KIND=SELECTED_REAL_KIND(12)) :: initialD_Xe=0, CalculateXenonDiffusion
!
!    TYPE(ValueList_t), POINTER :: Constants
!
!    LOGICAL :: found, check
!    !------------------------------------------------------------
!
!    Constants=>GetConstants()
!
!    cell_radius=GetConstReal(Constants, 'cell radius', found)
!    CALL FoundCheck(found, 'cell radius', 'fatal')
!
!    cellT1=GetConstReal(Constants, 'T1', found)
!    CALL FoundCheck(found, 'T1', 'fatal')
!
!    pressureT1=GetConstReal(Constants, 'T1 Pressure',found)
!    CALL FoundCheck(found, 'T1 Pressure', 'fatal')
!
!    temperatureT1=GetConstReal(Constants, 'T1 Temperature', found)
!    CALL FoundCheck(found, 'T1 Temperature', 'fatal')
!
!    T1vals = (/pressureT1,temperatureT1/)
!
!    initialD_Xe=CalculateXenonDiffusion(Model, 1, T1vals)
!
!    check = (cellT1>cell_radius**2/(2*initialD_Xe))
!
!    IF (check) THEN
!
!        decayrate=(initialD_Xe/cell_radius)*(1+(cell_radius/sqrt(initialD_Xe*cellT1*60))/&
!            tan(cell_radius/sqrt(initialD_Xe*cellT1*60)))
!
!    ELSE
!
!        CALL Fatal('CalculateDecayRate',&
!            'The cell T1 is shorter that what is possible given the cell radius')
!
!    END IF
!
!END FUNCTION CalculateDecayRate
!
!
!FUNCTION CalculateXeGamma(Model,n,argument)&
!    RESULT(XeGamma)
!
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: argument(3)
!    REAL(KIND=dp) :: concentration,pressure,temperature
!    REAL(KIND=dp) :: spinexchangerate
!    REAL(KIND=dp) :: binaryexchangerate
!    !-------------------------------------------------------
!    REAL(KIND=dp) :: he_fraction=0, xe_fraction=0, n2_fraction=0
!    REAL(KIND=dp) :: spinrelaxationrate
!    REAL(kind=dp) :: binary_term=0, vdWterm=0
!    !----------------------------------------------------------
!    REAL(KIND=dp) :: XeGamma
!    !-------------------------------------------------------
!    TYPE(ValueList_t), POINTER :: Materials
!    LOGICAL :: found
!    !-----------------------------------------------------------
!
!    !----Spin Exchange Rate-------------------------------------------------------
!
!    !Eventaully I might import more terms, so this is anticpating that eventuallity
!    concentration=argument(1)
!
!    Materials=GetMaterial()
!
!    binaryexchangerate=GetConstReal(Materials, 'binary exchange rate', found)
!    CALL FoundCheck(found, 'binaryexchangerate', 'fatal')
!
!    spinexchangerate=binaryexchangerate*Concentration
!
!    !------------------------------------------------------------------------------
!
!    !-------------------Spin Relaxation Rate--------------------------------------
!    pressure=argument(2)
!    temperature=argument(3)
!
!    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
!    CALL FoundCheck(found, 'xe fraction', 'fatal')
!
!    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
!    CALL FoundCheck(found, 'n2 fraction' , 'fatal')
!
!    he_fraction = GetConstReal(Materials, 'he fraction', found)
!    CALL FoundCheck(found, 'he fraction', 'fatal')
!
!    binary_term=(5D-6)*xe_fraction*((Pressure)/101325)*(273.15/Temperature)
!
!    vdWterm=6.72D-5*(1/(1+0.25*he_fraction/xe_fraction+1.05*n2_fraction/xe_fraction))
!
!    spinrelaxationrate=binary_term+vdWterm
!
!
!    XeGamma=spinexchangerate+spinrelaxationrate
!
!END FUNCTION CalculateXeGamma
!
!FUNCTION CalculateXeField(Model,n,argument)&
!    RESULT(XeField)
!
!    USE DefUtils
!    IMPLICIT None
!    TYPE(Model_t) :: Model
!    INTEGER :: n
!    REAL(KIND=dp) :: argument(4)
!    !-------------------------------------------------
!    REAL(KIND=dp) :: concentration,pressure,temperature, optrate
!    REAL(KIND=dp) :: alkali_alkali_spin_destruction_rate, xe_spin_destruction_rate
!    REAL(KIND=dp) :: xe_fraction, n2_fraction, he_fraction, xe_pressure,&
!        he_pressure, n2_pressure,tot_pressure_amg, ref_pressure
!    REAL(KIND=dp) :: SpinDestrucionRate
!    !-----------------------------------------------------------------
!    INTEGER :: ind
!    !-----------------------------------------------------------------
!    LOGICAL :: convert_density=.FALSE., he_term_included=.FALSE.
!    LOGICAL :: n2_term_included=.FALSE., vanderWall_term_included=.FALSE.
!    LOGICAL :: xe_term_included=.FALSE., local_temperature_used = .FALSE.
!    LOGICAL :: local_pressure_used = .FALSE., found=.FALSE.
!    !-----------------------------------------------------------
!    REAL(KIND=dp) :: SpinExchangeRate
!    REAL(KIND=dp) :: binaryexchangerate
!    !----------------------------------------------------------
!    REAL(KIND=dp) :: XeField
!    !--------------------------------------------------------
!    TYPE(ValueList_t), POINTER :: Materials
!    !--------------------------------------------------------------
!
!    !Spin Destruction Rate Calculation
!
!    REAL, PARAMETER :: rb_density_conversion_factor = 7.043279D24
!
!    SpinDestrucionRate = 0.0D0
!
!    Materials => GetMaterial()
!
!    !!!!!!!!!!!!!!!!!!!Alkali-Alkali Spin Destruction Term!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    alkali_alkali_spin_destruction_rate=&
!        GetConstReal(Materials, 'alkali-alkali spin destruction rate', found)
!    CALL FoundCheck(found, 'alkali-alkali spin destruction rate','fatal')
!
!    ! Assign the Concentration
!    Concentration=argument(1)
!
!    !Check if we should convert the density (we probably do want to
!    convert_density=GetLogical(Materials, 'convert density', found)
!    CALL FoundCheck(found, 'convert density', 'warn')
!
!    !Convert the alkali density to particles/m^3
!    IF (convert_density) THEN
!        Concentration = rb_density_conversion_factor*Concentration
!    END IF
!
!    SpinDestrucionRate=alkali_alkali_spin_destruction_rate*Concentration
!
!    !Inclusion of other terms in the spin-destruction rate: He-Rb, N2-Rb,
!    !and Van der Walls terms
!
!    he_term_included=GetLogical(Materials,&
!        'Helium Spin Destruction Term Included',found)
!
!    n2_term_included = GetLogical(Materials,&
!        'Nitrogen Spin Destruction Term Included',found)
!
!    xe_term_included = GetLogical(Materials,&
!        'Xenon Spin Destruction Term Included',found)
!
!    vanderWall_term_included = GetLogical(Materials,&
!        'vander Wall Spin Destruction Term Included',found)
!
!    IF (he_term_included .OR. n2_term_included .OR. &
!        xe_term_included .OR. vanderWall_term_included) THEN
!
!        !Get the temperature
!
!        Temperature = argument(2)
!
!        IF (Temperature .EQ. 0) Call Fatal('GetSpinDestructionRate',&
!            'Temperature variable not found.')
!
!
!        !Get the pressure
!
!        Pressure = argument(3)
!
!        !Check to make sure we actually found the solution
!
!
!        !        IF (Pressure .EQ. 0) Call Fatal('GetSpinDestructionRate',&
!        !            'Pressure variable not found. Check name of N-S variable.')
!
!
!        !Get the gas fractions
!        xe_fraction=GetConstReal(Materials, 'xe fraction', found)
!        CALL FoundCheck(found, 'xe fraction', 'fatal')
!
!        n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
!        CALL FoundCheck(found, 'n2 fraction' , 'fatal')
!
!        he_fraction = GetConstReal(Materials, 'he fraction', found)
!        CALL FoundCheck(found, 'he fraction', 'fatal')
!
!
!        !Call fatal if the gas fractions don't add to 1
!
!        IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
!            CALL Fatal('GetSpinDestructionRate', &
!                'Gas fractions do not add to 1')
!        END IF
!
!        !Calculate pressure in amagats
!
!        tot_pressure_amg = ((Pressure)/101325)*(273.15/Temperature)
!
!        xe_pressure=tot_pressure_amg*xe_fraction
!
!        n2_pressure=tot_pressure_amg*n2_fraction
!
!        he_pressure=tot_pressure_amg*he_fraction
!
!        !Implement the xenon contribution to spin destruction rate
!        IF (xe_term_included) THEN
!
!            xe_spin_destruction_rate = 0.0D0
!            xe_spin_destruction_rate = GetConstReal(Materials,&
!                'xe spin destruction rate',found)
!            CALL FoundCheck(found, 'xe spin destruction rate', 'fatal')
!
!            SpinDestrucionRate= SpinDestrucionRate+&
!                xe_spin_destruction_rate*xe_pressure
!
!        END IF
!
!        !Implement helium contribution to spin destruction rate
!        IF (he_term_included) THEN
!            SpinDestrucionRate = SpinDestrucionRate+&
!                24.6*(1+(Temperature+273.15-90)/94.6)*he_pressure
!        END IF
!
!        !Implement N2 contribution to spin destruction rate
!        IF (n2_term_included) THEN
!            SpinDestrucionRate=SpinDestrucionRate+&
!                170*(1+(Temperature+273.15-90)/194.36)*n2_pressure
!        END IF
!
!        !Implement van der Walls contribution to spin destruction rate
!        IF (vanderWall_term_included) THEN
!            SpinDestrucionRate=SpinDestrucionRate+&
!                6469/(xe_fraction+1.1*n2_fraction+3.2*he_fraction)
!        END IF
!    END IF
!
!    !-------------Spin Exchange Rate---------------------------------
!
!
!    !Eventaully I might import more terms, so this is anticpating that eventuallity
!
!    binaryexchangerate=GetConstReal(Materials, 'binary exchange rate', found)
!    CALL FoundCheck(found, 'binaryexchangerate', 'fatal')
!
!    SpinExchangeRate=binaryexchangerate*Concentration
!
!    optrate=argument(4)
!
!    XeField=SpinExchangeRate*(optrate/(SpinDestrucionRate+optrate))
!
!END FUNCTION CalculateXeField
