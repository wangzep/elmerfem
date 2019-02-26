!----------------------------------------------------------------------------
! File: C-ECalc
! Written by: Geoff Schrank
! Date : 29 Jan 2019
!----------------------------------------------------------------------------

!User defined function used to calculate viscosity using Chapmans-Enkogs theory
!as described in Bird, Steward, and Lightfoot. The function returns the viscosity
!in Pa-s.

FUNCTION Viscosity(Model, n, temperature) RESULT(mumix)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n

    !------------------------------------------------------------------------!
    REAL(KIND=dp) :: viscositycoef, omega
    REAL(KIND=dp), POINTER, DIMENSION(:) :: mass,sigma, K_esp, frac
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp), ALLOCATABLE,DIMENSION(:) :: mu
    REAL(KIND=dp) :: mumix

    INTEGER :: ind,numofparams
 
    !------------------------------------------------------------------------!
    !Declare constants-------------------------------------------------------
    !This coefficients and parameters are taken from Bird, Steward,
    !and Lightfoot.
    viscositycoef=2.6693d-5

    !-------------------------------------------------------------------------
    CALL GETINFO(Model,frac,mass,sigma,K_esp,numofparams)
    

    ALLOCATE(mu(numofparams))
    !Calculate viscosity for each of the three gas components
    !mu = 2.6693e-5*sqrt(m*T)/(simga^2*omega_mu)
    DO ind=1,numofparams
      CALL OMEGACALC(temperature,k_esp(ind),omega)	
      mu(ind)=viscositycoef*(mass(ind)*temperature)**(0.5)/(sigma(ind)**(2.0)*omega)
    END DO

    !Calculate the viscosity of the mixture   
    CALL MIXCALC(mu,mass,frac,numofparams,mumix)
    !Don't make a memory leak!
    DEALLOCATE(mu)
    !To get to the right units. mu is in g*cm/s and needs to be in kg*m/s
    !Conversion factor is 1 g*cm/s=.1 kg*m/s. 
    mumix=0.1*mumix

    RETURN

END FUNCTION

!User defined function used to calculate viscosity using Chapmans-Enkogs theory
!as described in Bird, Steward, and Lightfoot. The function returns the thermal
!conductivity in W/(m K).

FUNCTION ThermalCond(Model, n, temperature) RESULT(kmix)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    INTEGER :: numparams

    !------------------------------------------------------------------------!
    REAL(KIND=dp) :: thermcoef, omega
    REAL(KIND=dp), POINTER, DIMENSION(:) :: mass,sigma,K_esp,frac
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: k
    REAL(KIND=dp) :: kmix

    INTEGER :: numofparams,ind
    !------------------------------------------------------------------------!
    !Declare constants-------------------------------------------------------
    !These coefficients and parameters are taken from Bird, Steward,
    !and Lightfoot.
    thermcoef=1.9891d-4

    !-------------------------------------------------------------------------

    !Get the information about the gasses from the SIF
    !-------------------------------------------------------------------------
    CALL GETINFO(Model,frac,mass,sigma,K_esp,numofparams)
    ALLOCATE(k(numofparams))
    !-------------------------------------------------------------------------

    !Calculate viscosity for each of the three gas components

    DO ind=1,numofparams
       CALL OMEGACALC(temperature,K_esp(ind),omega)	
       k(ind)=thermcoef*(temperature/mass(ind))**(0.5)/(sigma(ind)**(2.0)*omega)
    END DO

    !Calculate the viscosity of the mixture   
    CALL MIXCALC(k,mass,frac,numofparams,kmix)
    
    DEALLOCATE(k)
    !k is in cal/(cm*s*C) and it needs to be in W/(m*K).
    !The converstion factor is 1 cal/(cm*s*C)=418.4 W/(m*K)
    kmix=418.4*kmix

    RETURN

END FUNCTION

SUBROUTINE GETINFO(Model,xfrac,mass,sigma,K_esp,numofparams)
    
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(ValueList_t), POINTER :: material

    REAL(KIND=dp),POINTER,DIMENSION(:) :: xfrac,mass,sigma,K_esp
    INTEGER :: numofparams, numcheck

    LOGICAL :: Found

    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
       CALL Fatal('getGasInfo', 'No material found')
    END IF
    !Get the information about the N2frac and Xefrac from the SIF
    !-------------------------------------------------------------------------
    xfrac = GetConstReal(material,'Gas Frac',Found)
    IF(.NOT. Found) THEN
        CALL Fatal('getGasInfo', 'Gas Frac not found')
    END IF
    numofparams=SIZE(xfrac)

    mass = GetConstReal(material,'Gas Mass',Found)
    IF(.NOT. Found) THEN
        CALL Fatal('getGasInfo', 'Gas Mass not found')
    END IF
    numcheck=SIZE(mass)

    IF(numcheck.NE.numofparams) THEN
        CALL Fatal('getGasInfo','mass size not equal to xfrac')
    END IF

    sigma = GetConstReal(material,'Gas Sigma',Found)
    IF(.NOT. Found) THEN
        CALL Fatal('getGasInfo', 'Gas Sigma not found')
    END IF
    numcheck=SIZE(sigma)

    IF(numcheck.NE.numofparams) THEN
        CALL Fatal('getGasInfo','sigma size not equal to xfrac')
    END IF
    K_esp = GetConstReal(material,'Gas K_Esp',Found)
    IF(.NOT. Found) THEN
        CALL Fatal('getGasInfo', 'Gas K_ESP not found')
    END IF
    numcheck=SIZE(mass)

    IF(numcheck.NE.numofparams) THEN
        CALL Fatal('getGasInfo','Gas K_Esp size not equal to xfrac')
    END IF
    !-------------------------------------------------------------------------
END SUBROUTINE



!Calculates a Lenard-Jones parameter given by 
!omega=1.16145/Tp^0.14874+0.52487/exp(0.77320*Tp)+2.16178/exp(2.43787*Tp)
!where Tp=kT/eps.
SUBROUTINE OMEGACALC(temperature, K_esp, omega)
    USE DefUtils
    IMPLICIT NONE 

    REAL(KIND=dp) :: temperature, K_esp, omega, tempprime
	
    tempprime=temperature/K_esp
    omega=1.16145/tempprime**(0.14874)+0.52487/exp(0.77320*tempprime)+2.16178/exp(2.43787*tempprime)
END SUBROUTINE




SUBROUTINE MIXCALC(tc,mass,xfrac,numparams,tcmix)
!Calculate the gas mixture viscosity using Wilke's equation, again
    !from Bird, Steward, and Lightfoot.

    USE DefUtils
    IMPLICIT NONE

    REAL(KIND=dp) :: tc(numparams),mass(numparams),xfrac(numparams),&
                     psi(numparams,numparams)
    REAL(KIND=dp) :: frontbit, middlebit,lastbit,denominator
    REAL(KIND=dp) :: tcmix

    INTEGER :: ind, jnd,numparams 
    
    !Initialize some of our variables and constants
    tcmix=0d0
    frontbit=8.0**(-0.5)

    !Calculate psi(i,j) for all components. I'm doing it in pieces because
    !it is easier to keep track of the parantheses.
    DO ind=1,numparams
    
       DO jnd=1,numparams

         middlebit=(1.0+mass(ind)/mass(jnd))**(-0.5)

         lastbit=(1.0+(tc(ind)/tc(jnd))**(0.5)*(mass(jnd)/mass(ind))**(0.25))**(2)

         psi(ind,jnd)=frontbit*middlebit*lastbit
	
       END DO
    
    END DO

    !Now do the actual sum
    DO ind=1,numparams
    
       denominator=0d0
    
       DO jnd=1,numparams
    
          denominator=denominator+xfrac(jnd)*psi(ind,jnd)
    
       END DO

       tcmix=tcmix+((xfrac(ind)*tc(ind))/denominator)

    END DO
END SUBROUTINE


