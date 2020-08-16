MODULE OPUtil


    IMPLICIT NONE


    INTERFACE

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
            REAL(KIND=dp) :: Argument(2)
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
            REAL(KIND=dp) :: Argument(2)
            REAL(KIND=dp) :: D_Xe
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

        FUNCTION calculateviscosity(Model,n,arguments)RESULT(viscosity)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: arguments
            REAL(KIND=dp) :: viscosity
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
        tot_numberdensity, xe_numberdensity, n2_numberdensity, he_numberdensity
    TYPE(ValueList_t), POINTER :: Materials, Constants
    LOGICAL :: found
    !-----------------------------------------------------------

    !Eventaully I might import more terms, so this is anticpating that eventuallity
    Concentration=Argument(1)

    Temperature=Argument(2)


    IF (Temperature .EQ. 0) Call Fatal('GetSpinDestructionRate',&
        'Temperature variable not found.')

    Pressure=Argument(3)

    Constants=>GetConstants()
    Materials=>GetMaterial()


    !Binary component
    binaryexchangerate=GetConstReal(Materials, 'binary exchange rate', found)
    CALL FoundCheck(found, 'binaryexchangerate', 'fatal')

    SpinExchangeRate=binaryexchangerate*Concentration


    !Molecular component

    !Get the gas fractions
    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')


    !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

    !Calculate pressure in amagats

    !Get Loschmidt's number if defined in constants

    loschmidt=GetConstReal(Model % Constants, 'loschmidts constant', found)
    CALL FoundCheck(found, 'loschmidts constant' , 'warn')
    IF (.NOT. found) THEN
        loschmidt= 2.6867811D25
    END IF



    tot_numberdensity = ((Pressure)/101325)*(273.15/Temperature)*loschmidt

    xe_numberdensity=tot_numberdensity*xe_fraction

    n2_numberdensity=tot_numberdensity*n2_fraction

    he_numberdensity=tot_numberdensity*he_fraction

    !Get molecular exchange rates

    xemolcrate=GetConstReal(Materials, 'xe molecular se rate', found)
    CALL FoundCheck(found, 'xe molecular se rate', 'fatal')

    hemolcrate=GetConstReal(Materials, 'he molecular se rate', found)
    CALL FoundCheck(found, 'xe molecular se rate', 'fatal')

    n2molcrate=GetConstReal(Materials, 'n2 molecular se rate', found)
    CALL FoundCheck(found, 'xe molecular se rate', 'fatal')


    SpinExchangeRate=SpinExchangeRate+(1/(xe_numberdensity/xemolcrate&
        +he_numberdensity/hemolcrate+n2_numberdensity/n2molcrate))*Concentration

END FUNCTION CalculateSpinExchangeRate

FUNCTION CalculateSpinRelaxationRate(Model,n,Argument)&
    RESULT(SpinRelaxationRate)
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: Argument(2)
    REAL(KIND=dp) :: Pressure,Temperature
    REAL(KIND=dp) :: he_fraction=0, xe_fraction=0, n2_fraction=0
    REAL(KIND=dp) :: he_ratio_term, n2_ratio_term, xe_vdW_term, xe_binary,&
        loschmidt
    REAL(KIND=dp) :: SpinRelaxationRate
    REAL(kind=dp) :: binary_term=0, vdWterm=0
    TYPE(ValueList_t), POINTER :: Materials
    LOGICAL :: found
    !------------------------------------------------------------------------------

    Materials=>GetMaterial()

    Pressure=Argument(1)
    Temperature=Argument(2)

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

        !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

    loschmidt=GetConstReal(Model % Constants, 'loschmidts constant', found)
    CALL FoundCheck(found, 'loschmidts constant' , 'warn')
    IF (.NOT. found) THEN
        loschmidt= 2.6867811D25
    END IF

    xe_binary=GetConstReal(Materials, 'binary spin relaxation', found)
    CALL FoundCheck(found, 'binary spin relaxation', 'fatal')

    binary_term=xe_binary*xe_fraction*&
        ((Pressure)/101325)*(273.15/Temperature)*loschmidt

    xe_vdW_term=GetConstReal(Materials, 'xe van der Waal spin relaxation', found)
    CALL FoundCheck(found, 'xe van der Waal spin relaxation', 'fatal')
    he_ratio_term=GetConstReal(Materials, 'spin relaxation he r', found)
    CALL FoundCheck(found, 'spin relaxation he r', 'fatal')
    n2_ratio_term=GetConstReal(Materials, 'spin relaxation n2 r', found)
    CALL FoundCheck(found, 'spin relaxation n2 r', 'fatal')

    vdWterm=xe_vdW_term*(1/(1+he_ratio_term*he_fraction/xe_fraction+&
        n2_ratio_term*n2_fraction/xe_fraction))

    SpinRelaxationRate=binary_term+vdWterm

END FUNCTION CalculateSpinRelaxationRate

FUNCTION CalculateXenonDiffusion(Model,n,Argument) RESULT(D_Xe)
    !Implements terms from Bird, Stewart, and Lightfoot. Diffusion in m^2/s.
    !I need to go back and document this better (probably when I revamp the XePol
    !solver.
    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: Argument(2)
    REAL(KIND=dp) :: D_Xe
    !-------------------------------------------------------
    REAL(KIND=dp) :: Pressure,Temperature
    REAL(KIND=dp) :: pressure_atm=0

    REAL(KIND=dp) :: mixKesp=0, tprime=0, omega=0, sigma=0, Mtot=0
    !------------------------------------------------------------
    !I believe these are all gotten from Lightfoot, sans the masses
    REAL(Kind=dp), PARAMETER :: A=1.858D-7, sigmaHe=2.576, sigmaXe=4.009,&
        KespHe=10.02, KespXe=234.7, massXe=131.29, massHe=4.003

    !------------------------------------------------------------

    !Getting assignments
    Pressure=Argument(1)
    Temperature=Argument(2)

    !Convert to atm
    pressure_atm=(Pressure)/101325

    !Actually doing the calcuation.
    mixKesp=sqrt(KespHe*KespXe)
    tprime = Temperature/mixKesp

    omega = (1.06036/tprime**(0.15610))+(0.19300/exp(0.47635*tprime))+&
        (1.03587/exp(1.52296*tprime))+(1.76474/exp(3.89411*tprime))

    sigma = 0.5*(sigmaHe+sigmaXe)

    Mtot = sqrt(1/massHe+1/massXe)

    D_Xe = (A*Temperature**(3/2)*Mtot)/(pressure_atm*omega*sigma**2)



END FUNCTION CalculateXenonDiffusion

FUNCTION CalculateDecayRate(Model,n,argument) RESULT(decayrate)

    USE DefUtils
    IMPLICIT None
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: argument(2)
    REAL(KIND=dp) :: pressure, temperature
    REAL(KIND=dp) :: decayrate
    !--------------------------------------------
    REAL(KIND=dp) :: cell_radius=0, cellT1=0,&
        pressureT1=0, temperatureT1=0, T1vals(2), minT1, rpi
    REAL(KIND=SELECTED_REAL_KIND(12)) :: initialD_Xe=0, CalculateXenonDiffusion

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

    T1vals = (/pressureT1,temperatureT1/)

    initialD_Xe=CalculateXenonDiffusion(Model, 1, T1vals)

    !Check to make sure that we are longer than the minumum T1

    minT1 = cell_radius**2/((rpi**2)*initialD_Xe)

    check = (cellT1>minT1)

    IF (check) THEN

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
    REAL(KIND=dp) :: Temp !Dummy Variable. Not actually used
    REAL(KIND=dp) :: RbNumDensity_m, Temperature
    LOGICAL :: found=.FALSE.
    !------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    !-----------------------------------------------------------------------

    Materials => GetConstants()

    Temperature=GetConstReal(Materials, 'alkali temperature', found)
    IF (.NOT. found) CALL Fatal('RbNumDensity',&
        'Temperature not found')


    RbNumDensity_m=(10**(9.55-4132/Temperature))/(1.380648521D-23*Temperature)


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


        !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

    pressure = GetConstReal(Materials, 'frequency pressure', found)
    CALL FoundCheck(found, 'frequency pressure', 'fatal')

    freqwidth=pressure/(18.9*xe_fraction+17.8*n2_fraction+18*he_fraction)

END FUNCTION calculaterbfrequency

!--------------------------------------------------------------------------------
FUNCTION calculatelaserheating(Model,n,arguments)RESULT(heating)
    !---------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments
    REAL(KIND=dp) :: heating
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials, Constants
    REAL(KIND=dp) :: laser_wavelength, speed_of_light, plank_constant,&
        concentration, spin_destruction_rate, alkali_polarization,&
        laser_frequency
    LOGICAL :: found, found1, found2

    Materials=>GetMaterial()
    Constants=>GetConstants()

    laser_wavelength = GetConstReal(Constants,'laser wavelength',found)
    CALL FoundCheck(found, 'laser wavelength', 'fatal')

    spin_destruction_rate = GetConstReal(Materials, 'spin destruction rate',found)
    CALL FoundCheck(found, 'spin destruction rate', 'fatal')

    alkali_polarization = GetConstReal(Materials, 'alkali polarization',found)
    CALL FoundCheck(found, 'alkali polarization', 'fatal')

    plank_constant = GetConstReal(Constants, 'planks constant',found1)
    speed_of_light = GetConstReal(Constants,'speed of light',found2)

    IF (.NOT. (found1 .AND. found2)) THEN

        plank_constant = 6.62607004D-34
        speed_of_light = 299792458.0D0

        CALL Warn('calculatelaserheating',&
            'One or more of the constants are not listed in the SIF. Using default values SI units.')
    END IF

    concentration=arguments

    laser_frequency = speed_of_light/laser_wavelength

    heating=plank_constant*laser_frequency*concentration*spin_destruction_rate*&
        alkali_polarization

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
    REAL(KIND=dp) :: arguments
    REAL(KIND=dp) :: evaprate
    !----------------------------------------------------------------------------
    REAL(KIND=dp) :: gasconstants, avagradosnumber, boltzmansconstant,&
        rpi, alpha, mass, temperature
    TYPE(ValueList_t), POINTER :: Materials, Constants
    LOGICAL :: found, found1, found2, found3
    !----------------------------------------------------------------------------

    rpi=4.D0*DATAN(1.D0)

    Materials=> GetMaterial()
    Constants=> GetConstants()

    temperature=arguments

    gasconstants=GetConstReal(Constants, 'gas constant', found1)
    avagradosnumber=GetConstReal(Constants, 'avagrados number', found2)
    boltzmansconstant=GetConstReal(Constants, 'boltzmans constant', found3)

    IF (.NOT. (found1 .OR. found2 .OR. found3)) THEN
        gasconstants=8.31446261815324
        avagradosnumber=6.02214076e23
        boltzmansconstant=1.38064852e-23

        CALL Warn('calculateevaprate',&
            'One or more of the constants are not listed in the SIF. Using default values SI units.')
    END IF

    alpha=GetConstReal(Materials, 'alpha', found)
    CALL FoundCheck(found, 'alpha', 'fatal')

    mass=GetConstReal(Materials, 'evap atomic mass', found)

    evaprate=(alpha*gasconstants*temperature)/(avagradosnumber*&
        sqrt(2*rpi*mass/avagradosnumber*boltzmansconstant*temperature))

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
    REAL(KIND=dp) :: arguments
    REAL(KIND=dp) :: heatranscoef
    !----------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: airheattrans, glassthermalconst, glassthickness
    LOGICAL :: found
    !----------------------------------------------------------------------------

    Materials=>GetMaterial()

    airheattrans=GetConstReal(Materials, 'air heat transfer coefficient', found)
    CALL FoundCheck(found, 'air heat transfer coefficient', 'fatal')

    glassthermalconst=GetConstReal(Materials, 'glass thermal conductivity', found)
    CALL FoundCheck(found, 'glass thermal conductivity', 'fatal')

    glassthickness=GetConstReal(Materials, 'glass thickness', found)
    CALL FoundCheck(found, 'glass thickness', 'fatal')

    heatranscoef = (1/airheattrans + glassthickness/glassthermalconst)**(-1)

END FUNCTION calculateheattransfercoef
!--------------------------------------------------------------------------------

FUNCTION calculateviscosity(Model,n,arguments)RESULT(viscosity)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments
    REAL(KIND=dp) :: viscosity
    !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp) :: n2_fraction, xe_fraction, he_fraction
    REAL(KIND=dp) :: viscosity_xe, viscosity_n2, viscosity_he
    LOGICAL :: found
    !------------------------------------------------------------------------------

    temperature = arguments

    Materials=>GetMaterial()

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

    !Calculate the individual viscosities
    viscosity_n2 = 1.66*(273.15+111)/(temperature+111)*(temperature/273.15)**(3/2)

    viscosity_he = 1.87*(273.15+79.4)/(temperature+79.4)*(temperature/273.15)**(3/2)

    viscosity_xe = 2.12*(273.15+252)/(temperature+252)*(temperature/273.15)**(3/2)

    !Add viscosities
    viscosity=(he_fraction*viscosity_he**(1/3)+n2_fraction*viscosity_n2**(1/3)&
        +xe_fraction*viscosity_xe**(1/3))**3

END FUNCTION calculateviscosity

FUNCTION calculateheatcapratio(Model,n,arguments)RESULT(heatcapratio)
    !Defines the Heat Capacity ratio as a function of gas fraction, Note* Xe and He (gamma = 1.666) are assumed to be perfect
    !monotonic gasses and N2 is assumed to be a perfect diatomic gas (gamma = 1.4).
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments
    REAL(KIND=dp) :: heatcapratio
    !-----------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Materials
    REAL(KIND=dp) :: n2_fraction, xe_fraction, he_fraction
    LOGICAL :: found
    !------------------------------------------------------------------------------

    Materials=>GetMaterial()

    xe_fraction=GetConstReal(Materials, 'xe fraction', found)
    CALL FoundCheck(found, 'xe fraction', 'fatal')

    n2_fraction = GetConstReal(Materials, 'n2 fraction', found)
    CALL FoundCheck(found, 'n2 fraction' , 'fatal')

    he_fraction = GetConstReal(Materials, 'he fraction', found)
    CALL FoundCheck(found, 'he fraction', 'fatal')

    !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

    heatcapratio=5/3*(xe_fraction+he_fraction)+7/5*n2_fraction

END FUNCTION calculateheatcapratio


FUNCTION calculatecp(Model,n,arguments)RESULT(cp)
    !Define heat capacity at constant pressure a function of gas fraction, He = 5196.118 N2 = 1039.67 Xe = 158.31 J/kg*K, from NIST Chemistry Webbook
    !Need to use mass fraction instead of mole or volume fraction.
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: arguments
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

    !Call fatal if the gas fractions don't add to 1

    IF (ABS(1-he_fraction-n2_fraction-xe_fraction)>1e-5) THEN
        CALL Fatal('GetSpinDestructionRate', &
            'Gas fractions do not add to 1')
    END IF

    avgmolarmass = xe_fraction*131.293+n2_fraction*28.0134+he_fraction*4.002602

    xemassfrac = xe_fraction*131.293/avgmolarmass
    n2massfrac = n2_fraction*28.0134/avgmolarmass
    hemassfrac = he_fraction*4.0026202/avgmolarmass

    cp = xemassfrac*158.31+n2massfrac*1039.67+5196.118*hemassfrac
END FUNCTION calculatecp

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
