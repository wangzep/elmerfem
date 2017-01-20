!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation,
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
!------------------------------------------------------------------------------
SUBROUTINE SEOPSolver( Model,Solver,dt,TransientSimulation )
    !------------------------------------------------------------------------------
    USE DefUtils

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
    INTEGER :: iter, maxiter
    LOGICAL :: Found

    !Laser Configuration Variables------------------------------------------
    REAL(KIND=dp) ::rubidium_wavelength,rubidium_freq_width,laser_wavelength,&
        laser_linewidth,oscillator_strength,Beta


    rubidium_wavelength = GetConstReal(Model % Constants,'rubidium wavelength',Found)
    laser_wavelength = GetConstReal(Model % Constants,'laser wavelength',Found)
    laser_linewidth = GetConstReal(Model % Constants,'laser line width',Found)
    rubidium_freq_width = GetConstReal(Model % Constants,'rubidium frequency width',Found)
    oscillator_strength = GetConstReal(Model % Constants,'oscillator strength', Found)
    !------------------------------------------------------------------------------

    maxiter = ListGetInteger( GetSolverParams(),&
        'Nonlinear System Max Iterations',Found,minv=1)
    IF(.NOT. Found ) maxiter = 1

    !Calculate Beta for the this laser configuration----------------



    Beta = BetaCalc(rubidium_wavelength,rubidium_freq_width,laser_wavelength,&
        laser_linewidth,oscillator_strength)
    !----------------------------------------------------------------------

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

            CALL LocalMatrix(  Element, n, nd + nb, Beta )

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
    SUBROUTINE LocalMatrix( Element, n, nd , Beta)
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: Beta, Absorption_Term(n), nRb(n), spin_destruction(n), &
            D,C,R, direction(3,n),a(3), Weight, Flux(n), RbPol_Term(n),one
        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
        !REAL(KIND=dp) :: rubidium_wavelength
        REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd)
        LOGICAL :: Stat,Found
        INTEGER :: i,t,p,q,dim
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(ValueList_t), POINTER :: BodyForce, Material
        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------

        dim = CoordinateSystemDimension()

        CALL GetScalarLocalSolution(Flux)
        CALL GetElementNodes( Nodes )

        MASS  = 0._dp
        STIFF = 0._dp
        FORCE = 0._dp
        RbPol_Term = 0._dp
        Absorption_Term = 0._dp


        Material => GetMaterial()
        nRb(1:n)=GetReal(Material,'rubidium number density',Found)
        spin_destruction(1:n) = GetReal(Material,'spin destruction rate',Found)

        Absorption_Term = Beta*nRb

        !Absorption_Term = nRb !This is just to make testing easier. Switch back to the above term when I'm done.

        DO i=1,n
            IF(Flux(i)<0) Flux(i) = 0 !This prevents non-phyical values of photon flux being used to make this calculation
            RbPol_Term(i) = (1.)-(Flux(i)/(Flux(i)+spin_destruction(i)))
            Absorption_Term(i) = Absorption_Term(i)*RbPol_Term(i)
        END DO

        direction = 0._dp

        DO i=1,dim
            direction(i,1:n)=GetReal(Material,&
                'laser direction '//TRIM(I2S(i)),Found)
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


            a = MATMUL(direction(:,1:n),Basis(1:n))
            R = SUM(Basis(1:n)*Absorption_Term(1:n))

            Weight = IP % s(t) * DetJ

            DO p=1,nd
                DO q=1,nd
                    ! Spacial derivative term
                    ! -----------------------------------
                    STIFF (p,q) = STIFF(p,q) + Weight * &
                        SUM(a(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

                    ! Absorption Term
                    ! -----------------------------------
                    STIFF(p,q) = STIFF(p,q) + Weight * R * Basis(q) * Basis(p)


                END DO
            END DO
        END DO

        IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
        CALL LCondensate( nd-nb, nb, STIFF, FORCE )
        CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
    END SUBROUTINE LocalMatrix
    !------------------------------------------------------------------------------


    ! Assembly of the matrix entries arising from the Neumann and Robin conditions
    !There are no possible Neumann or Robin conditions for this equation, set all to
    !zero.
    !------------------------------------------------------------------------------
    SUBROUTINE LocalMatrixBC( Element, n, nd )
        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: Weight
        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
        REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
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

        IP = GaussPoints( Element )

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
