MODULE OPUtil


    IMPLICIT NONE


    INTERFACE

        FUNCTION calculaterbextconc(Model,n,Temperature) RESULT(RbNumDensity_m)

            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: model
            INTEGER :: n
            REAL(KIND=dp) :: Temperature
            REAL(KIND=dp) :: RbNumDensity_m
        END

        FUNCTION CalculateSpinDestructionRate(Model,n,argument)&
            RESULT(SpinDestrucionRate)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: model
            INTEGER :: n
            REAL(KIND=dp) :: argument(3)
            REAL(KIND=dp) :: SpinDestrucionRate
        END

        FUNCTION CalculateRbPol(Model,n,argument) RESULT(RbPol)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: Model
            REAL(KIND=dp) :: argument(4)
            INTEGER :: n
            REAL(KIND=dp) :: RbPol
        END

        FUNCTION CalculateSpinExchangeRate(Model,n,Argument)&
            RESULT(SpinExchangeRate)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: Argument(3)
            REAL(KIND=dp) :: SpinExchangeRate
        END

        FUNCTION CalculateSpinRelaxationRate(Model,n,Argument)&
            RESULT(SpinRelaxationRate)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: Argument(3)
            REAL(KIND=dp) :: SpinRelaxationRate
        END

        FUNCTION CalculateDecayRate(Model,n,argument) RESULT(decayrate)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: argument(2)
            REAL(KIND=dp) :: decayrate
        END

        FUNCTION CalculateXenonDiffusion(Model,n,Argument) RESULT(D_Xe)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: Argument(3)
            REAL(KIND=dp) :: D_Xe
        END

        FUNCTION CalculateRubidiumDiffusion(Model,n,Argument) RESULT(D_Rb)
            !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: Argument(3)
            REAL(KIND=dp) :: D_Rb
        END

        FUNCTION calculaterbmumdensitym(Model,n,Temp) RESULT(RbNumDensity_m)
            USE DefUtils
            IMPLICIT None
            TYPE(Model_t) :: model
            INTEGER :: n
            REAL(KIND=dp) :: Temp
            REAL(KIND=dp) :: RbNumDensity_m
        END

        FUNCTION calculaterbfrequency(Model,n,pressure) RESULT(freqwidth)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: pressure
            REAL(KIND=dp) :: freqwidth
        END

        FUNCTION calculatelaserheating(Model,n,arguments)RESULT(heating)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: heating
        END

        FUNCTION calculateevaprate(Model,n,arguments)Result(evaprate)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: evaprate
        END

        FUNCTION calculateheattransfercoef(Model,n,arguments)RESULT(heatranscoef)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: heatranscoef
        END

        FUNCTION calculateviscosity(Model,n,arguments)RESULT(viscositytot)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: viscositytot
        END

        FUNCTION calculateheatcapratio(Model,n,arguments)RESULT(heatcapratio)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: heatcapratio
        END

        FUNCTION calculatecp(Model,n,arguments)RESULT(cp)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: cp
        END

        FUNCTION calculatethermalconductivity(Model,n,arguments)RESULT(ktot)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments(3)
            REAL(KIND=dp) :: ktot
        END

        FUNCTION calculatedensity(Model,n,arguments)RESULT(density)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments(3)
            REAL(KIND=dp) :: density
        END

        FUNCTION CalculateDiffusion(Concentration, Pressure, Temperature,&
            mass1, mass2, sigma1,sigma2, Kesp1, Kesp2) RESULT(diffusioncoef)
            !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
            USE DefUtils
            IMPLICIT None
            REAL(KIND=dp) :: diffusioncoef
            REAL(KIND=dp) :: Concentration,Pressure,Temperature
            REAL(KIND=dp) :: mass1, mass2, sigma1, sigma2, Kesp1,Kesp2
        END

        SUBROUTINE ArgumentCheck(Concentration, Pressure, Temperature, Caller)
            USE defutils

            IMPLICIT NONE
            REAL(KIND=dp), INTENT(IN) :: Pressure, Temperature
            CHARACTER(len=*), INTENT(IN) :: Caller
            REAL(KIND=dp) :: Concentration
        END

        SUBROUTINE GasFracCheck(xe_fraction, he_fraction, n2_fraction)
            USE defutils

            IMPLICIT NONE
            REAL(KIND=dp), INTENT(IN) :: xe_fraction, he_fraction, n2_fraction
        END

        SUBROUTINE FoundCheck(found,name,warn_fatal_flag)
            !------------------------------------------------------------------------------
            USE DefUtils

            IMPLICIT NONE

            LOGICAL, INTENT(IN) :: found
            CHARACTER(len=*), INTENT(IN) :: name
            CHARACTER(len=*), INTENT(IN) :: warn_fatal_flag
        END



    END INTERFACE
END MODULE

FUNCTION calculaterbextconc(Model,n,Temperature) RESULT(RbNumDensity_m)
    !-------------------------------------------------------------------------
    !Calculates Rb number density in m^-3 using Killian equation as presented
    !in Fink et al. 2005.
    !n_Rb=10^(9.55-4132/T)/kT
    !where T is the temperature in Kelvin and k is Boltzman's constant.
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: Temperature
    REAL(KIND=dp) :: RbNumDensity_m
    LOGICAL :: found=.FALSE.
    !------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    !-----------------------------------------------------------------------

    IF (Temperature .lt. 0) THEN
        CALL Fatal('calculaterbextconc',&
            'Temperature is less than 0, this is not physically possible')
    END IF

    RbNumDensity_m=(10**(9.55D0-4132.0D0/Temperature))/(1.380648521D-23*Temperature)

    IF (RbNumDensity_m .lt. 0) THEN
        CALL Fatal('calculaterbextconc',&
            'Calculated Rb Number Density is less than 0, this is not physically possible')
    END IF

END FUNCTION calculaterbextconc

FUNCTION CalculateSpinDestructionRate(Model,n,argument)&
    RESULT(SpinDestrucionRate)
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: argument(3)
    REAL(KIND=dp) :: Concentration, Temperature, Pressure, SpinDestrucionRate
    !----------------------------------------------------------------------

    REAL(KIND=dp) :: alkali_alkali_spin_destruction_rate, xe_spin_destruction_rate,&
        he_spin_destruction_rate, n2_spin_destruction_rate, G1
    REAL(KIND=dp) :: xe_fraction, n2_fraction, he_fraction, xe_numberdensity,&
        he_numberdensity, n2_numberdensity,tot_numberdensity, ref_pressure
    REAL(KIND=dp) :: loschmidt
    !-----------------------------------------------------------------
    INTEGER :: ind
    !-----------------------------------------------------------------
    LOGICAL :: convert_density=.FALSE., he_term_included=.FALSE.
    LOGICAL :: n2_term_included=.FALSE., vanderWall_term_included=.FALSE.
    LOGICAL :: xe_term_included=.FALSE., local_temperature_used = .FALSE.
    LOGICAL :: local_pressure_used = .FALSE., found=.FALSE.
    !--------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials, Constants
    !--------------------------------------------------------------

    SpinDestrucionRate = 0.0D0

    Materials => GetMaterial()
    Constants => GetConstants()

    ! Assign the Concentration
    Concentration=argument(1)


    !Get the pressure

    Pressure = argument(2)

    !Get Reference Pressure

    ref_pressure=GetConstReal(Materials, 'Reference Pressure', found)

    IF (found) THEN
        Pressure=ref_pressure+Pressure
    END IF

    !Get the temperature

    Temperature = argument(3)


    !Check the arguments are reasonable

    CALL ArgumentCheck(Concentration, Pressure, Temperature,&
        'CalculateSpinDestructionRate')
    !!!!!!!!!!!!!!!!!!!Alkali-Alkali Spin Destruction Term!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    alkali_alkali_spin_destruction_rate=&
        GetConstReal(Materials, 'alkali-alkali spin destruction rate', found)
    CALL FoundCheck(found, 'alkali-alkali spin destruction rate','fatal')


    IF (alkali_alkali_spin_destruction_rate .lt. 0) THEN
        CALL Fatal('OPUtil', 'Alkali-alkali spin destruction rate is less than 0')
    END IF

    !    !Check if we should convert the density (we probably do want to
    !    convert_density=GetLogical(Materials, 'convert density', found)
    !    CALL FoundCheck(found, 'convert density', 'warn')
    !
    !    !Convert the alkali density to particles/m^3
    !    IF (convert_density) THEN
    !        Concentration = rb_density_conversion_factor*Concentration
    !    END IF

    SpinDestrucionRate=alkali_alkali_spin_destruction_rate*Concentration

    !Inclusion of other terms in the spin-destruction rate: He-Rb, N2-Rb,
    !and Van der Walls terms

    he_term_included=GetLogical(Materials,&
        'Helium Spin Destruction Term Included',found)

    n2_term_included = GetLogical(Materials,&
        'Nitrogen Spin Destruction Term Included',found)

    xe_term_included = GetLogical(Materials,&
        'Xenon Spin Destruction Term Included',found)

    vanderWall_term_included = GetLogical(Materials,&
        'vander Wall Spin Destruction Term Included',found)

    IF (he_term_included .OR. n2_term_included .OR. &
        xe_term_included .OR. vanderWall_term_included) THEN

        !Check to make sure we actually found the solution


        !        IF (Pressure .EQ. 0) Call Fatal('GetSpinDestructionRate',&
        !            'Pressure variable not found. Check name of N-S variable.')


        !Get the gas fractions
        xe_fraction=GetConstReal(Materials, 'xe fraction', found)
        CALL FoundCheck(found, 'xe fraction', 'fatal')

        n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
        CALL FoundCheck(found, 'n2 fraction' , 'fatal')

        he_fraction = GetConstReal(Materials, 'he fraction', found)
        CALL FoundCheck(found, 'he fraction', 'fatal')

        CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)


        !Calculate pressure in amagats

        !Get Loschmidt's number if defined in constants

        loschmidt=GetConstReal(Constants, 'loschmidts constant', found)
        CALL FoundCheck(found, 'loschmidts constant' , 'warn')
        IF (.NOT. found) THEN
            loschmidt= 2.6867811D25
        END IF

        tot_numberdensity = ((Pressure)/101325.0D0)*(273.15D0/Temperature)*loschmidt

        xe_numberdensity=tot_numberdensity*xe_fraction

        n2_numberdensity=tot_numberdensity*n2_fraction

        he_numberdensity=tot_numberdensity*he_fraction

        !Implement the xenon contribution to spin destruction rate
        IF (xe_term_included) THEN

            xe_spin_destruction_rate = 0.0D0
            xe_spin_destruction_rate = GetConstReal(Materials,&
                'xe spin destruction rate',found)
            CALL FoundCheck(found, 'xe spin destruction rate', 'fatal')

            IF (xe_spin_destruction_rate .lt. 0) THEN
                CALL Fatal('OPUtil', 'Xe spin destruction rate is less than 0')
            END IF

            SpinDestrucionRate= SpinDestrucionRate+&
                xe_spin_destruction_rate*xe_numberdensity

        END IF

        !Implement helium contribution to spin destruction rate
        IF (he_term_included) THEN

            he_spin_destruction_rate = 0.0D0
            he_spin_destruction_rate = GetConstReal(Materials,&
                'he spin destruction rate',found)
            CALL FoundCheck(found, 'he spin destruction rate', 'fatal')

            IF (he_spin_destruction_rate .lt. 0) THEN
                CALL Fatal('OPUtil', 'He spin destruction rate is less than 0')
            END IF

            SpinDestrucionRate = SpinDestrucionRate+&
                he_spin_destruction_rate*he_numberdensity
        END IF

        !Implement N2 contribution to spin destruction rate
        IF (n2_term_included) THEN

            n2_spin_destruction_rate = 0.0D0
            n2_spin_destruction_rate = GetConstReal(Materials,&
                'n2 spin destruction rate',found)
            CALL FoundCheck(found, 'n2 spin destruction rate', 'fatal')

            IF (n2_spin_destruction_rate .lt. 0) THEN
                CALL Fatal('OPUtil', 'N2 spin destruction rate is less than 0')
            END IF

            SpinDestrucionRate=SpinDestrucionRate+&
                n2_spin_destruction_rate*n2_numberdensity
        END IF

        !Implement van der Walls contribution to spin destruction rate
        !See Nelson's 2001 thesis for details.
        IF (vanderWall_term_included) THEN


            G1 = 0.0D0
            G1 = GetConstReal(Materials,'short-very short transition density',found)
            CALL FoundCheck(found, 'short-very short transition density', 'fatal')

            SpinDestrucionRate=SpinDestrucionRate+&
                (0.385D0+0.642D0*1.0D0/(1.0D0+G1/tot_numberdensity))&
                *6469.0D0/(xe_fraction+1.1D0*n2_fraction+3.2D0*he_fraction)
        END IF
    END IF

    IF (SpinDestrucionRate .lt. 0) THEN
        CALL Fatal('CalculateSpinDestructionRate',&
            'Calculated spin destruction rate is less than 0')
    END IF

END FUNCTION CalculateSpinDestructionRate

FUNCTION CalculateRbPol(Model,n,argument) RESULT(RbPol)
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: argument(4), sdargument(3)
    REAL(KIND=dp) :: spindestrucionrate, optrate
    REAL(KIND=dp) :: RbPol
    REAL(KIND=dp) :: CalculateSpinDestructionRate
    !--------------------------------------------------------------------------------

    optrate=argument(1)

    IF (optrate .lt. 0) THEN
        optrate = 0
        !CALL Fatal('CalculateRbPol',&
         !   'optrate is less than 0. This is not physical.')
    END IF

    sdargument= (/argument(2),argument(3),argument(4)/)

    spindestrucionrate=CalculateSpinDestructionRate(Model, n, sdargument)

    RbPol=optrate/(optrate+spindestrucionrate)

    IF (RbPol .gt. 1 .or. RbPol .lt. 0) THEN
        CALL Fatal('CalculateRbPol',&
            'Calculated RbPol is not bound between 0 and 1')
    END IF

END FUNCTION CalculateRbPol

FUNCTION CalculateSpinExchangeRate(Model,n,Argument)&
    RESULT(SpinExchangeRate)

    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: Argument(3)
    REAL(KIND=dp) :: Concentration, Temperature, Pressure
    REAL(KIND=dp) :: SpinExchangeRate
    REAL(KIND=dp) :: binaryexchangerate, xemolcrate, n2molcrate, hemolcrate
    REAL(KIND=dp) :: xe_fraction,n2_fraction, he_fraction, loschmidt,&
        tot_numberdensity, xe_numberdensity, n2_numberdensity, he_numberdensity,&
        ref_pressure
    TYPE(ValueList_t), POINTER :: Materials, Constants
    LOGICAL :: found, is131 = .FALSE.
        !-----------------------------------------------------------


    Constants=>GetConstants()
    Materials=>GetMaterial()

    Concentration=Argument(1)

    Temperature=Argument(3)

    IF (Temperature .EQ. 0.0D0) Call Fatal('CalculateSpinExchangeRate',&
        'Temperature variable not found.')

    Pressure=Argument(2)

        !Get Reference Pressure

    ref_pressure=GetConstReal(Materials, 'Reference Pressure', found)

    IF (found) THEN
        Pressure=ref_pressure+Pressure
    END IF

    CALL ArgumentCheck(Concentration, Pressure, Temperature,&
        'CalculateSpinExchangeRate')

    !Binary component
    binaryexchangerate=GetConstReal(Materials, 'binary exchange rate', found)
    CALL FoundCheck(found, 'binaryexchangerate', 'fatal')


    IF (binaryexchangerate .lt. 0) THEN
        CALL Fatal('SESolver', 'Binary exchange rate is less than 0')
    END IF

    SpinExchangeRate=binaryexchangerate*Concentration


    !Molecular component

    !Get the gas fractions
    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    !Calculate pressure in amagats

    !Get Loschmidt's number if defined in constants

    loschmidt=GetConstReal(Model % Constants, 'loschmidts constant', found)
    CALL FoundCheck(found, 'loschmidts constant' , 'warn')
    IF (.NOT. found) THEN
        loschmidt= 2.6867811D25
    END IF

    tot_numberdensity = ((Pressure)/101325.0D0)*(273.15D0/Temperature)*loschmidt

    xe_numberdensity=tot_numberdensity*xe_fraction

    n2_numberdensity=tot_numberdensity*n2_fraction

    he_numberdensity=tot_numberdensity*he_fraction

    !Get molecular exchange rates

    xemolcrate=GetConstReal(Materials, 'xe molecular se rate', found)
    CALL FoundCheck(found, 'xe molecular se rate', 'fatal')

    IF (xemolcrate .lt. 0) THEN
        CALL Fatal('OPUtil', 'Xe molecular se rate is less than 0')
    END IF

    hemolcrate=GetConstReal(Materials, 'he molecular se rate', found)
    CALL FoundCheck(found, 'xe molecular se rate', 'fatal')

    IF (hemolcrate .lt. 0) THEN
        CALL Fatal('OPUtil', 'He molecular se rate is less than 0')
    END IF

    n2molcrate=GetConstReal(Materials, 'n2 molecular se rate', found)
    CALL FoundCheck(found, 'xe molecular se rate', 'fatal')

    IF (n2molcrate .lt. 0) THEN
        CALL Fatal('OPUtil', 'N2 molecular se rate is less than 0')
    END IF

    SpinExchangeRate=SpinExchangeRate+(1.0D0/(xe_numberdensity/xemolcrate&
        +he_numberdensity/hemolcrate+n2_numberdensity/n2molcrate))*Concentration

        !Addition for 131Xe rate, just divide by 2.28
    is131 = GetLogical(Materials, '131Xe', found)
    IF (is131) THEN
        SpinExchangeRate=SpinExchangeRate/2.28
    END IF

    IF (SpinExchangeRate .lt. 0) THEN
        CALL Fatal('CalculateSpinExchangeRate',&
            'Calculated spin exchange rate is less than 0')
    END IF

END FUNCTION CalculateSpinExchangeRate

FUNCTION CalculateSpinRelaxationRate(Model,n,Argument)&
    RESULT(SpinRelaxationRate)
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: Argument(3)
    REAL(KIND=dp) :: Concentration,Pressure,Temperature
    REAL(KIND=dp) :: he_fraction=0, xe_fraction=0, n2_fraction=0, ref_pressure
    REAL(KIND=dp) :: he_ratio_term, n2_ratio_term, xe_vdW_term, xe_binary,&
        loschmidt
    REAL(KIND=dp) :: SpinRelaxationRate
    REAL(kind=dp) :: binary_term=0, vdWterm=0
    TYPE(ValueList_t), POINTER :: Materials
    LOGICAL :: found
    !------------------------------------------------------------------------------

    Materials=>GetMaterial()

    Concentration=Argument(1)

    Pressure=Argument(2)

        !Get Reference Pressure

    ref_pressure=GetConstReal(Materials, 'Reference Pressure', found)

    IF (found) THEN
        Pressure=ref_pressure+Pressure
    END IF

    Temperature=Argument(3)

    CALL ArgumentCheck(Concentration, Pressure, Temperature,&
        'CalculateSpinRelaxationRate')

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    loschmidt=GetConstReal(Model % Constants, 'loschmidts constant', found)
    CALL FoundCheck(found, 'loschmidts constant' , 'warn')
    IF (.NOT. found) THEN
        loschmidt= 2.6867811D25
    END IF

    xe_binary=GetConstReal(Materials, 'binary spin relaxation', found)
    CALL FoundCheck(found, 'binary spin relaxation', 'fatal')

    IF (xe_binary .lt. 0) THEN
        CALL Fatal('OPUtil', 'Binary spin relaxation rate is less than 0')
    END IF

    binary_term=xe_binary*xe_fraction*&
        ((Pressure)/101325.0D0)*(273.15D0/Temperature)*loschmidt

    xe_vdW_term=GetConstReal(Materials, 'xe van der Waal spin relaxation', found)
    CALL FoundCheck(found, 'xe van der Waal spin relaxation', 'fatal')

    IF (xe_vdW_term .lt. 0) THEN
        CALL Fatal('SESolver', 'vdW spin relaxation is less than 0')
    END IF


    he_ratio_term=GetConstReal(Materials, 'spin relaxation he r', found)
    CALL FoundCheck(found, 'spin relaxation he r', 'fatal')
    IF (he_ratio_term .lt. 0) THEN
        CALL Fatal('SESolver', 'spin relaxation he r is less than 0')
    END IF

    n2_ratio_term=GetConstReal(Materials, 'spin relaxation n2 r', found)
    CALL FoundCheck(found, 'spin relaxation n2 r', 'fatal')

    IF (n2_ratio_term .lt. 0) THEN
        CALL Fatal('SESolver', 'spin relaxation n2 r is less than 0')
    END IF

    vdWterm=xe_vdW_term*(1.0D0/(1.0D0+he_ratio_term*he_fraction/xe_fraction+&
        n2_ratio_term*n2_fraction/xe_fraction))

    SpinRelaxationRate=binary_term+vdWterm

    IF (SpinRelaxationRate .lt. 0) THEN
        CALL Fatal('CalculateSpinRelaxationRate',&
            'Calculated spin relaxation rate is less than 0')
    END IF

END FUNCTION CalculateSpinRelaxationRate

FUNCTION CalculateXenonDiffusion(Model,n,Argument) RESULT(D_Xe)
    !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
    !I need to go back and document this better (probably when I revamp the XePol
    !solver.
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: Argument(3)
    REAL(KIND=dp) :: D_Xe
    !-------------------------------------------------------
    REAL(KIND=dp) :: Concentration,Pressure,Temperature
    REAL(KIND=dp) :: ref_pressure
    !------------------------------------------------------------
    !From Lightfoot Table E.1 page 864
    REAL(KIND=dp), PARAMETER ::sigmaHe=2.576D0, sigmaXe=4.009D0,&
        KespHe=10.02D0, KespXe=234.7D0
    REAL(KIND=dp), PARAMETER ::massXe=131.29D0, massHe=4.003D0
    REAL(KIND=dp) :: CalculateDiffusion
    TYPE(ValueList_t), POINTER :: Materials
    LOGICAL :: found
    !------------------------------------------------------------

    !Getting assignments
    Concentration=Argument(1)

    Pressure=Argument(2)

    !Get the refeerence pressure for Perfect Gas compressibility
    !model runs.
    Materials=>GetMaterial()
    ref_pressure=GetConstReal(Materials, 'Reference Pressure', found)

    IF (ref_pressure .lt. 0) THEN
        CALL Fatal('CalculateXenonDiffusion',&
            'Reference Pressure is less than 0')
    END IF

    IF (found) THEN
        Pressure=ref_pressure+Pressure
    END IF

    !Get temperature assignment
    Temperature=Argument(3)

    CALL ArgumentCheck(Concentration, Pressure, Temperature,&
        'CalculateXenonDiffusion')

    D_Xe=CalculateDiffusion(Concentration, Pressure, Temperature,&
        massHe, massXe, sigmaHe, sigmaXe, KespHe, KespXe)

    IF (D_Xe .lt. 0) THEN
        CALL Fatal('CalculateXenonDiffusion',&
            'Calculated Xenon Diffusion is less than 0')
    END IF

END FUNCTION CalculateXenonDiffusion

FUNCTION CalculateRubidiumDiffusion(Model,n,Argument) RESULT(D_Rb)
    !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: Argument(3)
    REAL(KIND=dp) :: D_Rb
    !-------------------------------------------------------
    REAL(KIND=dp) :: Concentration,Pressure,Temperature
    REAL(KIND=dp) :: ref_pressure
    !------------------------------------------------------------
    !From Lightfoot Table E.1 page 864 and
    !JOURNAL OF RESEARCH of the National Bureau of Stondards Vol. 84, No.6,
    !November- December 1979, Mountain and Haan for the Lenard-Jones parameters of Rb
    REAL(KIND=dp), PARAMETER ::sigmaHe=2.576D0, sigmaRb=4.48D0,&
        KespHe=10.02D0, KespRb=383D0
    REAL(KIND=dp), PARAMETER ::massRb=85.4678D0, massHe=4.003D0
    REAL(KIND=dp) :: CalculateDiffusion
    TYPE(ValueList_t), POINTER :: Materials
    LOGICAL :: found
    !------------------------------------------------------------

    !Getting assignments
    Concentration=Argument(1)

    Pressure=Argument(2)

    !Get the refeerence pressure for Perfect Gas compressibility
    !model runs.
    Materials=>GetMaterial()
    ref_pressure=GetConstReal(Materials, 'Reference Pressure', found)

    IF (ref_pressure .lt. 0) THEN
        CALL Fatal('CalculateXenonDiffusion',&
            'Reference Pressure is less than 0')
    END IF

    IF (found) THEN
        Pressure=ref_pressure+Pressure
    END IF

    !Get temperature assignment
    Temperature=Argument(3)

    CALL ArgumentCheck(Concentration, Pressure, Temperature,&
        'CalculateRubidiumDiffusion')

    D_Rb=CalculateDiffusion(Concentration, Pressure, Temperature,&
        massHe, massRb, sigmaHe, sigmaRb, KespHe, KespRb)
END FUNCTION CalculateRubidiumDiffusion

FUNCTION CalculateDecayRate(Model,n,argument) RESULT(decayrate)

    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: argument(3)!Dummy for Elmer, doesn't actually do anything.
    REAL(KIND=dp) :: decayrate
    !--------------------------------------------
    REAL(KIND=dp) :: cell_radius=0, cellT1=0,&
        pressureT1=0, temperatureT1=0, T1vals(3), minT1, rpi,&
        dumb=0
    !--------------------------------------------
            !From Lightfoot Table E.1 page 864
    REAL(KIND=dp), PARAMETER ::sigmaHe=2.576D0, sigmaXe=4.009D0,&
        KespHe=10.02D0, KespXe=234.7D0
    REAL(KIND=dp), PARAMETER ::massXe=131.29D0, massHe=4.003D0
    REAL(KIND=dp) :: initialD_Xe=0, CalculateDiffusion
    TYPE(ValueList_t), POINTER :: Constants

    LOGICAL :: found, check
    !------------------------------------------------------------

    rpi=4.D0*DATAN(1.D0)

    Constants=>GetConstants()

    cell_radius=GetConstReal(Constants, 'cell radius', found)
    CALL FoundCheck(found, 'cell radius', 'fatal')

    cellT1=GetConstReal(Constants, 'T1', found)
    CALL FoundCheck(found, 'T1', 'fatal')

    pressureT1=GetConstReal(Constants, 'T1 Pressure',found)
    CALL FoundCheck(found, 'T1 Pressure', 'fatal')

    temperatureT1=GetConstReal(Constants, 'T1 Temperature', found)
    CALL FoundCheck(found, 'T1 Temperature', 'fatal')

    initialD_Xe=CalculateDiffusion(dumb, pressureT1, temperatureT1,&
        massHe, massXe, sigmaHe, sigmaXe, KespHe, KespXe)

    !Check to make sure that we are longer than the minumum T1
    !See https://arxiv.org/abs/1911.01574 for this equation.

    minT1 = cell_radius**2.0D0/((rpi**2.0D0)*initialD_Xe)

    check = (cellT1>minT1)

    IF (check) THEN

        !See https://arxiv.org/abs/1911.01574 for this equation.
        decayrate=(initialD_Xe/cell_radius)*(1-(cell_radius/sqrt(initialD_Xe*cellT1))/&
            TAN(cell_radius/sqrt(initialD_Xe*cellT1)))

    ELSE
        PRINT *, 'The minimum T1 for a cell of this radius is:', minT1
        CALL Fatal('CalculateDecayRate',&
            'The cell T1 is shorter that what is possible given the cell radius')

    END IF

END FUNCTION CalculateDecayRate

FUNCTION calculaterbmumdensitym(Model,n,Temp) RESULT(RbNumDensity_m)
    !-------------------------------------------------------------------------
    !Calculates Rb number density in m^-3 using Killian equation as presented
    !in Fink et al. 2005.
    !n_Rb=10^(9.55-4132/T)/kT
    !where T is the temperature in Kelvin and k is Boltzman's constant.
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: Temp(3) !Dummy Variable. Not actually used
    REAL(KIND=dp) :: RbNumDensity_m, Temperature
    LOGICAL :: found=.FALSE.
    !------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    !-----------------------------------------------------------------------

    Materials => GetConstants()

    Temperature=GetConstReal(Materials, 'alkali temperature', found)
    IF (.NOT. found) CALL Fatal('RbNumDensity',&
        'Temperature not found')


    RbNumDensity_m=(10**(9.55D0-4132.0D0/Temperature))/(1.380648521D-23*Temperature)


END FUNCTION calculaterbmumdensitym

!-------------------------------------------------------------------
FUNCTION calculaterbfrequency(Model,n,pressure) RESULT(freqwidth)
    !--------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: pressure
    REAL(KIND=dp) :: freqwidth
    !-----------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials, Constants
    REAL(KIND=dp) :: xe_fraction, n2_fraction, he_fraction
    LOGICAL :: found

    Materials=>GetMaterial()
    Constants=>GetConstants()

    !Get the gas fractions
    xe_fraction = GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    pressure = GetConstReal(Materials, 'frequency pressure', found)
    CALL FoundCheck(found, 'frequency pressure', 'fatal')

    freqwidth=pressure/(18.9D0*xe_fraction+17.8D0*n2_fraction+18.0D0*he_fraction)

END FUNCTION calculaterbfrequency

!--------------------------------------------------------------------------------
FUNCTION calculatelaserheating(Model,n,arguments)RESULT(heating)
    !---------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(4), sdargument(3)
    REAL(KIND=dp) :: heating
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials, Constants
    REAL(KIND=dp) :: laser_wavelength, speed_of_light, plank_constant,&
        concentration, spin_destruction_rate, alkali_polarization,&
        laser_frequency, heatlimit
    REAL(KIND=dp) :: CalculateSpinDestructionRate, CalculateRbPol
    LOGICAL :: found, found1, found2, limitheat = .FALSE.

    sdargument= (/arguments(2),arguments(3),arguments(4)/)

    Materials=>GetMaterial()
    Constants=>GetConstants()

    laser_wavelength = GetConstReal(Constants,'laser wavelength',found)
    CALL FoundCheck(found, 'laser wavelength', 'fatal')

    spin_destruction_rate = GetConstReal(Materials, 'spin destruction rate',found)
    CALL FoundCheck(found, 'spin destruction rate', 'fatal')

    IF (spin_destruction_rate .EQ. 0) THEN
        spin_destruction_rate=CalculateSpinDestructionRate(Model, n, sdargument)
    END IF

    alkali_polarization = GetConstReal(Materials, 'alkali polarization',found)
    CALL FoundCheck(found, 'alkali polarization', 'fatal')

    IF (alkali_polarization .EQ. 0) THEN
        alkali_polarization=CalculateRbPol(Model, n, arguments)
    END IF

    plank_constant = GetConstReal(Constants, 'planks constant',found1)
    speed_of_light = GetConstReal(Constants,'speed of light',found2)

    IF (.NOT. (found1 .AND. found2)) THEN

        plank_constant = 6.62607004D-34
        speed_of_light = 299792458.0D0

        CALL Warn('calculatelaserheating',&
            'One or more of the constants are not listed in the SIF. Using default values SI units.')
    END IF

    concentration=arguments(2)

    IF (concentration .lt. 0) THEN
        CALL INFO('calculatelaserheating',&
            'Concentration is less than 0', level = 6)
        concentration = 0
    END IF

    laser_frequency = speed_of_light/laser_wavelength

    heating=plank_constant*laser_frequency*concentration*spin_destruction_rate*&
        alkali_polarization

    IF (heating .lt. 0) THEN
        CALL Fatal('calculatelaserheating',&
            'Calcualted laser heating is less than 0')
    END IF

    !Code to limit the amount of heat. This might be necessary if the heating is greater
    !than what can be supported by the mesh size.

    limitheat = GetLogical(Materials, 'Limit Heat', found)
    CALL FoundCheck(found, 'Limit Heat', 'warn')

    heatlimit = GetConstReal(Materials, 'Heat Limit', found)

    IF (limitheat) THEN
        IF (heating .gt. heatlimit) THEN
            heating = heatlimit
        END IF
    END IF

!----------------------------------------------------------------------------------
END FUNCTION calculatelaserheating
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
FUNCTION calculateevaprate(Model,n,arguments)Result(evaprate)
    !---------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: evaprate
    !----------------------------------------------------------------------------
    REAL(KIND=dp) :: gasconstants, avagradosnumber, boltzmansconstant,&
        rpi, alpha, mass, temperature, loschmidt, prefactor
    TYPE(ValueList_t), POINTER :: Constants, BC
    LOGICAL :: found, found1, found2, found3
    !----------------------------------------------------------------------------

    rpi=4.D0*DATAN(1.D0)

    BC=> GetBC()
    Constants=> GetConstants()

    temperature=arguments(3)

    IF (temperature .lt. 0) THEN
        CALL Fatal('calculate evaprate',&
            'Temperature is less than 0')
    END IF

    !Adjust to keep below boiling point of Rb.

    IF (temperature .GT. 960) THEN
        temperature = 960
    END IF

    !Gas constant and avagrdos number are not used in this formulation
    !gasconstants=GetConstReal(Constants, 'gas constant', found1)
    !avagradosnumber=GetConstReal(Constants, 'avagrados number', found2)
    boltzmansconstant=GetConstReal(Constants, 'boltzmans constant', found)

    IF (.NOT. (found)) THEN
        !gasconstants=8.31446261815324D0
        !avagradosnumber=6.02214076D23
        boltzmansconstant=1.38064852D-23

        CALL Warn('calculateevaprate',&
            'One or more of the constants are not listed in the SIF. Using default values SI units.')
    END IF

    alpha=GetConstReal(BC, 'alpha', found)
    CALL FoundCheck(found, 'alpha', 'fatal')

    mass=GetConstReal(BC, 'evap atomic mass', found)
    CALL FoundCheck(found, 'evap atomic mass', 'fatal')


    !Using the Herz-Knudsen equation h=alpha*(psat-p)*Na/sqrt(2*pi*M_Rb*k*T)
    !This calculates the prefactor for the (psat-p) term, and then is is used
    ! in the mass transfer coeficient. Note, we are not calculating mass, but
    !particle flux

    evaprate=(alpha)/&
        (sqrt(2*rpi*mass*boltzmansconstant*temperature))

    !Put in a prefactor to convert number density to pressure in pascal
    !To do this use pressure=(numdensity/loscmidt)*(Temperature/273.15D0)*101375.0D0

    loschmidt=GetConstReal(Model % Constants, 'loschmidts constant', found)
    CALL FoundCheck(found, 'loschmidts constant' , 'warn')
    IF (.NOT. found) THEN
        loschmidt= 2.6867811D25
    END IF

    prefactor=(temperature/273.15D0)*101375.0D0/loschmidt

    evaprate=prefactor*evaprate

    IF (evaprate .lt. 0) THEN
        CALL Fatal('calculateevaprate',&
            'Calculated evaporation rate is less than 0')
    END IF

!---------------------------------------------------------------------------------
END FUNCTION calculateevaprate
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
FUNCTION calculateheattransfercoef(Model,n,arguments)RESULT(heatranscoef)
    !---------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: heatranscoef
    !----------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp) :: airheattrans, glassthermalconst, glassthickness
    LOGICAL :: found
    !----------------------------------------------------------------------------

    BC=> GetBC()

    airheattrans=GetConstReal(BC, 'air heat transfer coefficient', found)
    CALL FoundCheck(found, 'air heat transfer coefficient', 'fatal')

    glassthermalconst=GetConstReal(BC, 'glass thermal conductivity', found)
    CALL FoundCheck(found, 'glass thermal conductivity', 'fatal')

    glassthickness=GetConstReal(BC, 'glass thickness', found)
    CALL FoundCheck(found, 'glass thickness', 'fatal')

    heatranscoef = (1/airheattrans + glassthickness/glassthermalconst)**(-1)

    IF (heatranscoef .lt. 0) THEN
        CALL Fatal('calcualteheatransfercoef',&
            'Calculated heat transfer coefficient is less than 0')
    END IF

END FUNCTION calculateheattransfercoef
!--------------------------------------------------------------------------------

FUNCTION calculateviscosity(Model,n,arguments)RESULT(viscositytot)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: viscositytot
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp) :: n2_fraction, xe_fraction, he_fraction
    !From Lightfoot page 26, eq. 1.4-14
    REAL(KIND=dp), PARAMETER :: A=2.6693D-6
    !From Lightfoot Table E.1 page 864
    REAL(KIND=dp), PARAMETER :: sigmaHe=2.576D0, sigmaXe=4.009D0, sigmaN2=3.667D0,&
        KespHe=10.02D0, KespXe=234.7D0, KespN2=99.8D0
    REAL(KIND=dp), PARAMETER :: massXe=131.29D0, massHe=4.003D0, massN2=28.0134D0
    REAL(KIND=dp) :: omega(3), tprime(3), mass(3), Kesp(3), sigma(3), frac(3)
    REAL(KIND=dp) :: phid(3,3)
    REAL(KIND=dp) :: viscosity(3), denom
    INTEGER :: i,j
    LOGICAL :: found
    !------------------------------------------------------------------------------

    temperature = arguments(3)

    Materials=>GetMaterial()

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    !Make vector assignments in order: (1) He, (2) N2, (3) Xe
    frac = (/he_fraction,n2_fraction,xe_fraction/)
    mass = (/massHe,massN2,massXe/)
    sigma = (/sigmaHe,sigmaN2,sigmaXe/)
    Kesp = (/KespHe,KespN2, KespXe/)


    !Calculate the tprimes and collsion integral, Lightfoot Table E.2 footnote
    !page 865
    !Calculate indivual viscosities; Lightfoot eq. 1.4-14
    DO i = 1,3
        tprime(i) = temperature/Kesp(i)
        omega(i)  = 1.16145/tprime(i)**(0.14874D0)+0.52487/EXP(0.77370D0*tprime(i))+&
            2.16178/EXP(2.43787*tprime(i))

        viscosity(i) = A * SQRT(mass(i)*temperature)/((sigma(i)**2)*omega(i))
    END DO

    !Dimensionaless overlap functions Lightfoot eq. 1.4-16
    DO i = 1, 3
        DO j = 1, 3
            phid(i,j) = (8.0D0)**(-0.5D0)*(1+mass(i)/mass(j))**(-0.5D0)*&
                (1+(viscosity(i)/viscosity(j))**(0.5D0)*(mass(j)/mass(i))**(0.25D0))**2.0D0
        END DO
    END DO

    viscositytot = 0.0D0

    !Lightfoot eq. 1.4-15
    DO i = 1, 3
        !instantiate denom as 0
        denom = 0.0D0

        DO j = 1, 3
            denom = denom + phid(i,j)*frac(j)
        END DO

        viscositytot = viscositytot + frac(i)*viscosity(i)/denom

    END DO

    IF (viscositytot .lt. 0) THEN
        CALL Fatal('calculateviscosity',&
            'Calculated viscosity is less than 0')
    END IF

END FUNCTION calculateviscosity

FUNCTION calculateheatcapratio(Model,n,arguments)RESULT(heatcapratio)
    !Defines the Heat Capacity ratio as a function of gas fraction, Note* Xe and He (gamma = 1.666) are assumed to be perfect
    !monotonic gasses and N2 is assumed to be a perfect diatomic gas (gamma = 1.4).
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: heatcapratio
    !-----------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: n2_fraction, xe_fraction, he_fraction
    REAL(KIND=dp) :: avgmolarmass, xemassfrac, n2massfrac, hemassfrac
    LOGICAL :: found
    !------------------------------------------------------------------------------

    Materials=>GetMaterial()

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    avgmolarmass = xe_fraction*131.293D0+n2_fraction*28.0134D0+he_fraction*4.002602D0

    xemassfrac = xe_fraction*131.293D0/avgmolarmass
    n2massfrac = n2_fraction*28.0134D0/avgmolarmass
    hemassfrac = he_fraction*4.0026202D0/avgmolarmass

    heatcapratio=(5.0D0/2.0D0*(xemassfrac+hemassfrac)+7.0D0/2.0D0*n2massfrac)/&
        (3.0D0/2.0D0*(xemassfrac+hemassfrac)+5.0D0/2.0D0*n2massfrac)

    IF (heatcapratio .lt. 0) THEN
        CALL Fatal('calculateheatcapratio', &
            'Calculated heat capacity ration is less than 0')
    END IF

END FUNCTION calculateheatcapratio


FUNCTION calculatecp(Model,n,arguments)RESULT(cp)
    !Define heat capacity at constant pressure a function of gas fraction,
    !He = 5196.118 J/kg*K N2 = 1039.67 J/kg*K Xe = 158.31 J/kg*K, from NIST Chemistry Webbook
    !Need to use mass fraction instead of mole or volume fraction.
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: cp
        !---------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: n2_fraction, xe_fraction, he_fraction
    REAL(KIND=dp) :: avgmolarmass, xemassfrac, n2massfrac, hemassfrac
    LOGICAL :: found
    !------------------------------------------------------------------------------

    Materials=>GetMaterial()

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    avgmolarmass = xe_fraction*131.293D0+n2_fraction*28.0134D0+he_fraction*4.002602D0

    xemassfrac = xe_fraction*131.293D0/avgmolarmass
    n2massfrac = n2_fraction*28.0134D0/avgmolarmass
    hemassfrac = he_fraction*4.0026202D0/avgmolarmass

    cp = xemassfrac*158.31D0+n2massfrac*1039.67D0+5196.118D0*hemassfrac

    IF (cp .lt. 0) THEN
        CALL Fatal('calculatecp',&
            'Calculated cp is less than 0')
    END IF

END FUNCTION calculatecp

FUNCTION calculatethermalconductivity(Model,n,arguments)RESULT(ktot)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: ktot
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp) :: n2_fraction, xe_fraction, he_fraction
    !From Lightfoot page 26, eq. 9.3-13, modified to go from cal/(cm*s*K) to W/(m*K)
    !1 cal/(s*cm*K)=418.4 W/(m*K). So 1.9891e-4*418.4 = 0.083223944
    REAL(KIND=dp), PARAMETER :: A=0.083223944
    !From Lightfoot Table E.1 page 864
    REAL(KIND=dp), PARAMETER :: sigmaHe=2.576D0, sigmaXe=4.009D0, sigmaN2=3.667D0,&
        KespHe=10.02D0, KespXe=234.7D0, KespN2=99.8D0
    REAL(KIND=dp), PARAMETER :: massXe=131.29D0, massHe=4.003D0, massN2=28.0134D0
    REAL(KIND=dp) :: omega(3), tprime(3), mass(3), Kesp(3), sigma(3), frac(3)
    REAL(KIND=dp) :: phid(3,3)
    REAL(KIND=dp) :: k(3), denom
    INTEGER :: i,j
    LOGICAL :: found
    !------------------------------------------------------------------------------

    temperature = arguments(3)

    IF (temperature .lt. 0) THEN
        CALL Warn('calculatethermalconductivity',&
            'Temperature is less than 0')

        temperature = 0
    END IF

    Materials=>GetMaterial()

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    CALL GasFracCheck(xe_fraction, he_fraction, n2_fraction)

    !Make vector assignments in order: (1) He, (2) N2, (3) Xe
    frac = (/he_fraction,n2_fraction,xe_fraction/)
    mass = (/massHe,massN2,massXe/)
    sigma = (/sigmaHe,sigmaN2,sigmaXe/)
    Kesp = (/KespHe,KespN2, KespXe/)


    !Calculate the tprimes and collsion integral, Lightfoot Table E.2 footnote
    !page 865
    !Calculate indivual viscosities; Lightfoot eq. 9.3-13
    DO i = 1,3
        tprime(i) = temperature/Kesp(i)
        omega(i)  = 1.16145/tprime(i)**(0.14874D0)+0.52487/EXP(0.77370D0*tprime(i))+&
            2.16178/EXP(2.43787*tprime(i))

        k(i) = A * SQRT(temperature/mass(i))/((sigma(i)**2)*omega(i))
    END DO

    !Dimensionaless overlap functions Lightfoot eq. 9.3-17
    DO i = 1, 3
        DO j = 1, 3
            phid(i,j) = (8.0D0)**(-0.5D0)*(1+mass(i)/mass(j))**(-0.5D0)*&
                (1+(k(i)/k(j))**(0.5D0)*(mass(j)/mass(i))**(0.25D0))**2.0D0
        END DO
    END DO

    ktot = 0.0D0

    !Lightfoot eq. 1.4-15
    DO i = 1, 3
        !instantiate denom as 0
        denom = 0.0D0

        DO j = 1, 3
            denom = denom + phid(i,j)*frac(j)
        END DO

        ktot = ktot + frac(i)*k(i)/denom

    END DO

    IF (ktot .lt. 0) THEN
        CALL Fatal('calculatethermalconductivity', &
            'Calculated thermal conductivity is less than 0')
    END IF

END FUNCTION calculatethermalconductivity

FUNCTION calculatedensity(Model,n,arguments)RESULT(density)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments(3)
    REAL(KIND=dp) :: density
    !----------------------------------------------------------
    REAL(KIND=dp) :: RConst, Concentration, Pressure, Temperature
    REAL(KIND=dp) :: specheatratio, cp, refpressure
    REAL(KIND=dp) :: calculateheatcapratio, calculatecp
    TYPE(ValueList_t), POINTER :: Materials
    LOGICAL :: found
    !-----------------------------------------------------------

    Materials=> GetMaterial()

    Concentration = arguments(1)
    Pressure = arguments(2)
    Temperature = arguments(3)

    refpressure = GetConstReal(Materials, 'Reference Pressure', found)

    IF (found) THEN
        Pressure=Pressure+refpressure
    END IF

    CALL ArgumentCheck(Concentration, Pressure, Temperature, 'calculatedensity')

    specheatratio = calculateheatcapratio(Model, n, arguments)
    cp = calculatecp(Model, n, arguments)

    RConst = ((specheatratio-1)/specheatratio)*cp
    density = Pressure/(RConst*Temperature)

END FUNCTION calculatedensity


FUNCTION CalculateDiffusion(Concentration, Pressure, Temperature,&
    mass1, mass2, sigma1,sigma2, Kesp1, Kesp2) RESULT(diffusioncoef)
    !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
    USE DefUtils
    IMPLICIT None
    REAL(KIND=dp) :: diffusioncoef
    REAL(KIND=dp) :: Concentration,Pressure,Temperature
    REAL(KIND=dp) :: mass1, mass2, sigma1, sigma2, Kesp1,Kesp2
    !------------------------------------------------------------------------
    REAL(KIND=dp) :: pressure_atm=0

    REAL(KIND=dp) :: mixKesp=0, tprime=0, omega=0, sigma=0, Mtot=0
    !------------------------------------------------------------
    !From Lightfoot pg. 526
    REAL(KIND=dp), PARAMETER :: A=1.8583D-7
    !------------------------------------------------------------
    !Convert to atm
    pressure_atm=(Pressure)/101325.0D0

    !Actually doing the calcuation.
    !Lightfoot eq. 17.3-15
    mixKesp=sqrt(Kesp1*Kesp2)
    !Lightfoot Table E.2 footnote page 865
    tprime = Temperature/mixKesp

    !Approximation to the collision integral from Lightfoot Table E.2 footnote, page 865.
    omega = (1.06036D0/tprime**(0.15610D0))+(0.19300D0/exp(0.47635D0*tprime))+&
        (1.03587D0/exp(1.52296D0*tprime))+(1.76474D0/exp(3.89411D0*tprime))

    !Lightfoot eq.17.3-14
    sigma = 0.5D0*(sigma1+sigma2)

    !Lightfoot eq. 17.3-12
    Mtot = sqrt(1/mass1+1/mass2)

    !Light eq. 17.3-12
    diffusioncoef = (A*Temperature**(3.0D0/2.0D0)*Mtot)/(pressure_atm*omega*sigma**2D0)

END FUNCTION CalculateDiffusion

SUBROUTINE ArgumentCheck(Concentration, Pressure, Temperature, Caller)
    USE defutils

    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Pressure, Temperature
    CHARACTER(len=*), INTENT(IN) :: Caller
    REAL(KIND=dp) :: Concentration
    !-----------------------------------------------------------------
    IF (Concentration .lt. 0) THEN
        !CALL INFO(Caller, 'Concentration is less than 0', level = 6)
        Concentration = 0
    END IF

    IF (Pressure .lt. 0 .or. Pressure .eq. 0) THEN
        CALL Fatal(Caller, 'Pressure is less than 0')
    END IF

    IF (Temperature .lt. 0 .or. Temperature .eq. 0) THEN
        CALL Fatal(Caller, 'Temperature is less than 0')
    END IF

END SUBROUTINE ArgumentCheck

SUBROUTINE GasFracCheck(xe_fraction, he_fraction, n2_fraction)
    USE defutils

    IMPLICIT NONE
    REAL(Kind=dp), INTENT(IN) :: xe_fraction, he_fraction, n2_fraction
    !---------------------------------------------------------------
    IF (xe_fraction .lt. 0 .or. xe_fraction .gt. 1) THEN
        CALL Fatal('OPUtil', 'Xenon fraction is not bound by 0 and 1')
    END IF

    IF (n2_fraction .lt. 0 .or. n2_fraction .gt. 1) THEN
        CALL Fatal('OPUtil', 'N2 fraction is not bound by 0 and 1')
    END IF

    IF (he_fraction .lt. 0 .or. he_fraction .gt. 1) THEN
        CALL Fatal('OPUtil', 'He fraction is not bound by 0 and 1')
    END IF

    !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

END SUBROUTINE GasFracCheck


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
            CALL Warn('OPUtil', outputstring)
        ELSE
            CALL Fatal('OPUtil', outputstring)
        END IF
    END IF

!-------------------------------------------------------------------------
END SUBROUTINE FoundCheck
