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

    binary_term=(5D-6)*xe_fraction*((Pressure)/101325)*(273.15/Temperature)

    vdWterm=6.72D-5*(1/(1+0.25*he_fraction/xe_fraction+1.05*n2_fraction/xe_fraction))

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
        pressureT1=0, temperatureT1=0, T1vals(2)
    REAL(KIND=SELECTED_REAL_KIND(12)) :: initialD_Xe=0, CalculateXenonDiffusion

    TYPE(ValueList_t), POINTER :: Constants

    LOGICAL :: found, check
    !------------------------------------------------------------

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

    check = (cellT1>cell_radius**2/(2*initialD_Xe))

    IF (check) THEN

        decayrate=(1/cell_radius)*(1+(cell_radius/sqrt(initialD_Xe*cellT1))/&
            tan(cell_radius/sqrt(initialD_Xe*cellT1)))

    ELSE

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
