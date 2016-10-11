!----------------------------------------------------------------------------
! File: BetaCalc
! Written by: Geoff Schrank
! Date : 30 Sept, 2016
!----------------------------------------------------------------------------

!User defined function used in conjuction with SEOPLaser.F90. This function
!calculates the Beta parameter of the PDE presented in Fink et al. PRA 72 (2005)
!This function will get the Rubidium parameters from the SIF file and return
!the Beta parameter. It encodes equation (A7)-(A10) from Fink et al.

FUNCTION BetaCalc(Model, n) RESULT(Beta(n))
    IMPLICIT NONE
    USE DefUtils
    TYPE(Model_t) :: Model
    INTEGER :: n

    !------------------------------------------------------------------------!
    REAL(KIND=dp) :: C,PI
    REAL(KIND=dp) :: electron_radius, oscillator_strength, laser_wavelength, &
        spectral_overlap, laser_linewidth, frequency_shift,&
        absorb_laser_ratio,rubidium_linewidth,laser_frequency,&
        rubidium_frequency, laser_freq_width,&
        rubidium_freq_width, w, w_prime
    REAL(KIND=dp) :: Beta(n)
    !------------------------------------------------------------------------!
    !Declare constants-------------------------------------------------------
    C = 2.998e8
    PI=4.D0*DATAN(1.D0)
    !-------------------------------------------------------------------------

    !Get the material information about the rubidium and the laser from the SIF
    !file---------------------------------------------------------------------
    !Material => GetMaterial()

    !rubidium_wavelength = GetReal(Material,'rubidium wavelength',Found)
    !laser_wavelength = GetReal(Material,'laser wavelength',Found)
    !laser_linewidth = GetReal(Material,'laser line width',Found)
    !rubidium_linewidth = GetReal(Material,'rubidium linewidth',Found)
    !For testing
rubidium_wavelength = 795e-9
laser_wavelength = 795e-9
laser_linewidth = 0.3e-9
rubidium_linewidth = 0.01e-9

    !-------------------------------------------------------------------------

    !Define spectral overlap function (eq. A7, A10 of Fink's paper)
    rubidium_frequency = C/rubidium_wavelength
    rubidium_freq_width = C*(1/(rubidium_wavelength-rubidium_linewidth/2) -&
        1/(rubidium_wavelength+rubidium_linewidth/2)

    laser_frequency = C/laser_wavelength
    laser_freq_width = C*(1/(laser_wavelength-laser_linewidth/2) -&
        1/(laser_wavelength+laser_linewidth/2))


    absorb_laser_ratio = rubidium_freq_width/laser_freq_width
    frequency_shift = 2*(laser_frequency - rubidium_frequency)/laser_freq_width

    w = EXP(LOG(2)*(COMPLEX(absorb_laser_ratio,frequency_shift))**2&
        )*ERFC(SQRT(LOG(2))*COMPLEX(absorb_laser_ratio,frequency_shift))

    w_prime = REAL(w)

    Beta(1:n) = 2*SQRT(PI*LOG(2))*(electron_radius*oscillator_strength*&
        laser_wavelength**2*w_prime)/laser_linewidth


END FUNCTION
