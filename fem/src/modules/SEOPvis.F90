!----------------------------------------------------------------------------
! File: SEOPvis
! Written by: Geoff Schrank
! Date : 29 Jan 2019
!----------------------------------------------------------------------------

!User defined function used to calculate viscosity using Chapmans-Enkogs theory
!as described in Bird, Steward, and Lightfoot. The function returns the viscosity
!in Pa-s.

FUNCTION SEOPvis(Model, n, temperature) RESULT(mu)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n

    !------------------------------------------------------------------------!
    REAL(KIND=dp) :: N2Frac, XeFrac
    REAL(KIND=dp) :: viscositycoef, massHe, massXe, massN2, sigmaHe, sigmaN2, &
                     sigmaXe, k_espHe,k_espN2,k_espXe
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp) :: muHe, muN2, muXe
    REAL(KIND=dp) :: omega,tempprime
    REAL(KIND=dp) :: mucomp(3),mass(3),psi(3,3),xfrac(3),denominator, frontcoef,&
                     middlebit, lastbit
    REAL(KIND=dp) :: mu
    INTEGER :: ind,jnd
    LOGICAL :: FLAG, Found
    !------------------------------------------------------------------------!
    !Declare constants-------------------------------------------------------
    !These coefficients and parameters are taken from Bird, Steward,
    !and Lightfoot.
    viscositycoef=2.6693d-5
    massHe=4.003
    massN2=28.013
    massXe=131.29
    sigmaHe=2.576
    sigmaN2=3.667
    sigmaXe=4.009
    k_espHe=10.2
    k_espN2=99.8
    k_espXe=234.7
    frontcoef=8.0**(-0.5)

    !Initialize mu so that it isn't some random memory location.    
    mu=0d0

    !-------------------------------------------------------------------------

    !Get the information about the N2frac and Xefrac from the SIF
    !-------------------------------------------------------------------------
    N2Frac = GetConstReal(Model % Constants,'N2Frac',Found)
    XeFrac = GetConstReal(Model % Constants,'XeFrac',Found)
    !-------------------------------------------------------------------------

    !Calculate viscosity for each of the three gas components
    !mu = 2.6693e-5*sqrt(m*T)/(simga^2*omega_mu)
    CALL OMEGACALC(temperature,k_espHe,omega)	
    muHe=viscositycoef*(massHe*temperature)**(0.5)/(sigmaHe**(2.0)*omega)

    CALL OMEGACALC(temperature,k_espXe,omega)
    muXe=viscositycoef*(massXe*temperature)**(0.5)/(sigmaXe**(2.0)*omega)

    CALL OMEGACALC(temperature,k_espN2,omega)
    muN2=viscositycoef*(massN2*temperature)**(0.5)/(sigmaN2**(2.0)*omega)

    !Calculate the gas mixture viscosity using Wilke's equation, again
    !from Bird, Steward, and Lightfoot.
    mucomp(1)=muHe
    mucomp(2)=muXe
    mucomp(3)=muN2

    mass(1)=massHe
    mass(2)=massXe
    mass(3)=massN2

    xfrac(1)=1-N2Frac-XeFrac
    xfrac(2)=N2Frac
    xfrac(3)=XeFrac

    !Calculate psi(i,j) for all components. I'm doing it in pieces because
    !it is easier to keep track of the parantheses.
    DO ind=1,3
    
    DO jnd=1,3

    middlebit=(1.0+mass(ind)/mass(jnd))**(-0.5)

    lastbit=(1.0+(mucomp(ind)/mucomp(jnd))**(0.5)*(mass(jnd)/mass(ind))**(0.25))**(2)

    psi(ind,jnd)=frontcoef*middlebit*lastbit
	
    END DO
    
    END DO

    !Now do the actual sum
    DO ind=1,3
    
    denominator=0d0
    
    DO jnd=1,3
    
    denominator=denominator+xfrac(jnd)*psi(ind,jnd)
    
    END DO

    mu=mu+xfrac(ind)*mucomp(ind)/denominator

    END DO
    
    !To get to the right units. mu is in g*cm/s and needs to be in kg*m/s
    !Conversion factor is 1 g*cm/s=.1 kg*m/s.
    mu=0.1*mu

    RETURN

END FUNCTION



!Calculates a Lenard-Jones parameter given by 
!omega=1.16145/Tp^0.14874+0.52487/exp(0.77320*Tp)+2.16178/exp(2.43787*Tp)
!where Tp=kT/eps.
SUBROUTINE OMEGACALC(temperature, k_esp, omega)
    USE DefUtils
    IMPLICIT NONE 

    REAL(KIND=dp) :: temperature, k_esp, omega, tempprime
	
    tempprime=temperature/k_esp
    omega=1.16145/tempprime**(0.14874)+0.52487/exp(0.77320*tempprime)+2.16178/exp(2.43787*tempprime)
END SUBROUTINE

