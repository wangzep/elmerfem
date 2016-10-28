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
!******************************************************************************
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
            CALL LocalMatrix(  Element, n, nd+nb, Model, Solver )
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
    SUBROUTINE LocalMatrix( Element, n, nd, Model, Solver )

        !------------------------------------------------------------------------------
        INTEGER :: n, nd
        TYPE(Element_t), POINTER :: Element
        !------------------------------------------------------------------------------
        REAL(KIND=dp) :: nRb(n), Beta(n), &
            K, Weight, spin_destruction(n), SOL(n)
        REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
        REAL(KIND=dp) :: STIFF(nd,nd), MASS(nd,nd), FORCE(nd), LOAD(nd,nd)

        LOGICAL :: Found, stat
        INTEGER :: i,t,p,q,dim
        REAL(KIND=dp), POINTER :: direction(:,:)
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(ValueList_t), POINTER ::  Material
        TYPE(Solver_t) :: Solver
        TYPE(Model_t) :: Model
        TYPE(Nodes_t) :: Nodes
        SAVE Nodes
        !------------------------------------------------------------------------------

        dim = CoordinateSystemDimension()

        CALL GetElementNodes( Nodes )
        CALL GetScalarLocalSolution(SOL,UElement=Element)

        STIFF = 0._dp

        !Material Properties Def
        !Eventually some of this will all be passed from another solver
        Material => GetMaterial()
        nRb(1:n)=GetReal(Material,'rubidium number density',Found)
        spin_destruction(1:n) = GetReal(Material,'spin destruction rate',Found)

        !DO i =1,dim
        !    direction(i,1:n)=GetReal(Material,&
        !        'propogation direction'//TRIM(I2S(i)), Found)
        !END DO
        CALL GetConstRealArray(Solver % Values, direction, 'laser direction',Found)
        IF(.NOT.Found) CALL Fatal('SEOPLaser','Unable to find laser direction')

        Beta(1:n) = BetaCalc(Model,n) !Calculate Beta from (A6), I have a feeling this isn't expressed correctly

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



            !N = MATMUL(direction(:,1:n),Basis(1:n)) I'm not sure this line is needed
            K = SUM(Basis(1:n)*Beta(1:n)*(1-(SOL(1:n)/(SOL(1:n)+spin_destruction(1:n))))) !I think this properly implements the Picard manuver
            Weight = IP % s(t) * DetJ



            DO p=1,nd
                DO q=1,nd
                    ! Absorption Term (N*grad(u),grad(v)):
                    ! -----------------------------------

                    STIFF (p,q) = STIFF(p,q) + Weight * &
                         SUM(direction(1:dim,1)*dBasisdx(q,1:dim)) * Basis(p)

                    ! Polarization Term (K*u,v)
                    ! -----------------------------------
                    STIFF(p,q) = STIFF(p,q) + Weight * K*Basis(q) * Basis(p)

                END DO
            END DO

            FORCE(1:nd) = 0._dp
            MASS(1:nd,1:nd) = 0._dp
        END DO

        !IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
        IF(TransientSimulation)CALL Fatal('SEOPLaser','SEOPLaser Solver cannot be transient')
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


    !User defined function used in conjuction with SEOPLaser.F90. This function
    !calculates the Beta parameter of the PDE presented in Fink et al. PRA 72 (2005)
    !This function will get the Rubidium parameters from the SIF file and return
    !the Beta parameter. It encodes equation (A7)-(A10) from Fink et al.

    FUNCTION BetaCalc(Model, n) RESULT(Beta)
        USE DefUtils
        IMPLICIT NONE
        TYPE(Model_t) :: Model
        INTEGER :: n

        !------------------------------------------------------------------------!
        REAL(KIND=dp) :: TWO, C, LOG2
        REAL(KIND=dp) :: electron_radius, oscillator_strength, laser_wavelength, &
            laser_linewidth, frequency_shift,&
            absorb_laser_ratio,laser_frequency,&
            rubidium_frequency, laser_freq_width,&
            rubidium_freq_width, rubidium_wavelength,w_prime, w_dprime,&
            w_input_real,w_input_imaginary, Beta(n)
        LOGICAL :: FLAG, Found
        !------------------------------------------------------------------------!
        !Declare constants-------------------------------------------------------
        TWO = 2.0D0
        C = 299792458.0D0
        LOG2 = LOG(TWO)
        electron_radius = 2.8179403267e-15
        !-------------------------------------------------------------------------

        !Get the information about the rubidium and the laser from the SIF
        !file---------------------------------------------------------------------


        rubidium_wavelength = GetConstReal(Model % Constants,'rubidium wavelength',Found)
        laser_wavelength = GetConstReal(Model % Constants,'laser wavelength',Found)
        laser_linewidth = GetConstReal(Model % Constants,'laser line width',Found)
        rubidium_freq_width = GetConstReal(Model % Constants,'rubidium frequency width',Found)
        oscillator_strength = GetConstReal(Model % Constants,'oscillator strength', Found)

            !For testing---------------------------------!!!!!!!!!!!!!!!!!!!!!!!
        !rubidium_wavelength = 800e-9
        !rubidium_freq_width = 126.65e9

        !laser_wavelength = 800e-9
        !laser_linewidth = 2e-9

        !oscillator_strength = 1.0/3.0


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

        !print *,'wprime',w_prime

        Beta(1:n) = 2*DSQRT(PI*LOG2)*(electron_radius*oscillator_strength*&
            laser_wavelength**2*w_prime)/laser_linewidth

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
END SUBROUTINE SEOPLaser
!------------------------------------------------------------------------------
