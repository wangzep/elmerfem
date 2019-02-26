!----------------------------------------------------------------------------
! File: SEOPvis
! Written by: Geoff Schrank
! Date : 29 Jan 2019
!----------------------------------------------------------------------------

!User defined function used to calculate viscosity using Chapmans-Enkogs theory
!as described in Bird, Steward, and Lightfoot. The function returns the viscosity
!in g/(cm*s).

FUNCTION SEOPtherm(Model, n, temperature) RESULT(k)
    USE DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n

    !------------------------------------------------------------------------!
    REAL(KIND=dp) :: N2Frac, XeFrac
    REAL(KIND=dp) :: thermcoef, massHe, massXe, massN2, sigmaHe, sigmaN2, &
                     sigmaXe, k_espHe,k_espN2,k_espXe
    REAL(KIND=dp) :: temperature
    REAL(KIND=dp) :: kHe, kN2, kXe
    REAL(KIND=dp) :: omega,tempprime
    REAL(KIND=dp) :: kcomp(3),mass(3),psi(3,3),xfrac(3),denominator, frontcoef,&
                     middlebit, lastbit
    REAL(KIND=dp) :: k
    INTEGER :: ind,jnd
    LOGICAL :: FLAG, Found
    !------------------------------------------------------------------------!
    !Declare constants-------------------------------------------------------
    !These coefficients and parameters are taken from Bird, Steward,
    !and Lightfoot.
    thermcoef=1.9891d-4
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

    !Initialize k so that it isn't some random memory location.    
    k=0d0

    !-------------------------------------------------------------------------

    !Get the information about the N2frac and Xefrac from the SIF
    !-------------------------------------------------------------------------
    N2Frac = GetConstReal(Model % Constants,'N2Frac',Found)
    XeFrac = GetConstReal(Model % Constants,'XeFrac',Found)
    !-------------------------------------------------------------------------

    !Calculate viscosity for each of the three gas components
    CALL OMEGACALC(temperature,k_espHe,omega)	
    kHe=thermcoef*(temperature/massHe)**(0.5)/(sigmaHe**(2.0)*omega)

    CALL OMEGACALC(temperature,k_espXe,omega)
    kXe=thermcoef*(temperature/massXe)**(0.5)/(sigmaXe**(2.0)*omega)

    CALL OMEGACALC(temperature,k_espN2,omega)
    kN2=thermcoef*(temperature/massN2)**(0.5)/(sigmaN2**(2.0)*omega)

    !Calculate the gas mixture viscosity using Wilke's equation, again
    !from Bird, Steward, and Lightfoot.
    kcomp(1)=kHe
    kcomp(2)=kXe
    kcomp(3)=kN2

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

    lastbit=(1.0+(kcomp(ind)/kcomp(jnd))**(0.5)*(mass(jnd)/mass(ind))**(0.25))**(2)

    psi(ind,jnd)=frontcoef*middlebit*lastbit
	
    END DO
    
    END DO

    !Now do the actual sum
    DO ind=1,3
    
    denominator=0d0
    
    DO jnd=1,3
    
    denominator=denominator+xfrac(jnd)*psi(ind,jnd)
    
    END DO

    k=k+xfrac(ind)*kcomp(ind)/denominator

    END DO
    !k is in cal/(cm*s*C) and it needs to be in W/(m*K).
    !The converstion factor is 1 cal/(cm*s*C)=418.4 W/(m*K)
    k=418.4*k



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
    omega=1.16145/tempprime**(0.14874) + 0.52487/exp(0.77320*tempprime) +&
    2.16178/exp(2.43787*tempprime)
END SUBROUTINE

