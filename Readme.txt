The code is devoted to the solution of the linearized gravity-wave equations
using a multilayer approach.

I. PROGRAM STRUCTURE AND INPUT/OUTPUT

1. List of files
==================
The following files are included in the code:
- Waves.f: main file
- WavesAux.f: file containing auxiliary routines, such as cubic spline
              interpolation and sorting routines
- WavesPrint.f: file for printing output data
- IRI.f: file for reading and handling IRI data
- lapack.f, lapack1.f, lapack2.f, lapack3.f: files containing LAPACK routines
- fft.f: file containing Fast Fourier Transform (FFT) routines (in its present form,
         the hard-coded variable controlling the spatial Fourier transform over
         altitude, DoSpatialFourierTransform, is set to false)
- nrlmsise00_modified.f: file for the neutral atmosphere empirical model included
                         in sami2py
- hwm93.f, hwm07e_modified.f90, hwm14.f90: files for the horizontal wind models
                                           included in sami2py

The code is compiled with gfortran using the make command, which creates the
executable file waves.x.


2. Input data
==============
The input parameters and their descriptions are provided in the file
WaveData.namelist. In addition, the output of the International Reference
Ionosphere (IRI) code available at
                  https://ccmc.gsfc.nasa.gov/models/IRI~2016/
is used. From the data produced by this code, we read
- the date (year, month, and day) and the time,
- the geographic latitude and longitude,
- the magnetic dip angle,
- the solar radio flux f10.7 and its 81-day average,
- the number density of O+ ions as a function of altitude.

The IRI output data are stored in the directory 'IRIData/', which contains the
following files:

    File name       Location                            Lat  Long   Height[km]
---------------------------------------------------------------------------
- iri_outputA.txt  Jicamarca(Peru)	                    −12   283     300
- iri_outputB.txt  Arecibo(Puerto Rico)	                +18   293     300
- iri_outputC.txt  Millstone Hill(USA)	                +42   288     300
- iri_outputD.txt  Saint-Santin(France)	                +44     2     300
- iri_outputE.txt  EISCAT TromsøEISCAT Tromsø (Auroral) +70    19     300
- iri_outputF.txt  Svalbard archipelago(Norway)	        +80    15     300

The altitude grid starts at zmin = 80 km and ends at zmax = 500 km, with a step
size of dz = 1.0 km. The user may generate a custom iri_output.txt file by running
the IRI code and specifying the corresponding file name in WaveData.namelist.

The IRI data are subsequently used in a manner analogous to that in the SAMI2 model
of the Naval Research Laboratory. In SAMI2, the neutral atmospheric parameters,
namely the neutral number density, total mass density, and temperature, are
specified using the Mass Spectrometer Incoherent Scatter (MSIS) model, whereas
the meridional and zonal winds are specified using the Horizontal Wind Model (HWM)

The derivatives of the background parameters are computed using central finite
differences. Prior to applying the finite-difference calculations, the background
parameters are smoothed by means of cubic spline interpolation with regularization.


3. Output files
================
The output data are written to the following directories
- AtmosPars: atmospheric parameters,
- EigenValWaveNumSF: eigenvalues and wavenumbers for a single-frequency wave,
- WaveResultsSF: wave parameters for a single-frequency wave,
- WaveResultsWP: wave parameters for a time-dependent wave packet,

1. Atmospheric parameters versus altitude are written to the directory AtmosPars
   in the following files:
   - the temperature in 'Temperature.dat'
   - the mass density in 'MassDensity.dat'
   - the pressure in 'Pressure.dat'·
   - the horizontal velocity in the southward direction in 'UVelocity.dat'
   - the atmospheric scale height in 'ScaleHeightAtmos.dat'
   - the density scale height in 'ScaleHeightDensity.dat'
   - the specific heat constant in 'HeatConstantVolume.dat'
   - the ratio of specific heats in 'RatioSpecificHeats.dat'
   - the sound speed in 'SoundSpeed.dat'

   - the temperature derivative in 'TemperatureDerivative.dat'
   - the horizontal-velocity derivative in 'UVelocityDerivative.dat'
   - the mass-density derivative in 'MassDensityDerivative.dat'
   - the pressure derivative in 'PressureDerivative.dat'
   - the density-scale height derivative in 'HScaleDensDerivative.dat'

   - the number density of O+ ions in 'AtmosPars/IonNumberDensity.dat'
   - the neutral-ion frequency in 'AtmosPars/FrequencyNeutralIon.dat'
   - the ion-neutral frequency in 'AtmosPars/FrequencyIonNeutral.dat'
   - the diffusion velocity in 'AtmosPars/DiffusionVelocity.dat'
   - the ion number density derivative in 'AtmosPars/IonNumberDensityDerivative.dat'

2. Eigenvalues and wavenumbers for a single-frequency wave versus altitude are
   written to the directory 'EigenValWaveNumSF' in the following files:
   - the real part of the eigenvalues for all modes in
    '            'EigenValRealLinModelX.dat',
     where X denotes TypeLinModel, i.e., the linearization model.
   - the real part of the eigenvalues for ascending and descending gravity waves in
                 'EigenValRealAscDescGWLinModelX.dat'
   - the imaginary part of the eigenvalues for all modes in
  '              'EigenValImagLinModelX.dat'

   - the imaginary part of the wavenumbers for all modes in
                 'WavenumberImagLinModelX.dat'
   - the imaginary part of the wavenumbers for ascending and descending gravity in
                 'WavenumberImagAscDescGWLinModelX.dat'
   - the real part of the wavenumbers for all modes in
                 'WavenumberRealLinModelX.dat'

3. Wave parameters computed by the Global-Matrix Method for Amplitudes (GMMA),
   for a single-frequency wave versus altitude are written to the directory
   'WaveResultsSF' in the following files:
   - the horizontal velocity in
                 'UVelocitySolMetGMMALinModelX.dat'
   - the sum of ascending and descending horizontal velocities in
                 'UVelocityAscPlusDescSolMetGMMALinModelX.dat'
   - the ascending horizontal velocity in
                 'UVelocityAscSolMetGMMALinModelX.dat'
   - the desscending horizontal velocity in
                 'UVelocityDescSolMetGMMALinModelX.dat'

   - the vertical velocity in
                 'WVelocitySolMetGMMALinModelX.dat'
   - the sum of ascending and descending vertical velocities in
                 'WVelocityAscPlusDescSolMetGMMALinModelX.dat'
   - the ascending vertical velocity in
                 'WVelocityAscSolMetGMMALinModelX.dat'
   - the desscending vertical velocity in
                 'WVelocityDescSolMetGMMALinModelX.dat'

   - the temperature in
                 'TemperatureSolMetGMMALinModelX.dat'
   - the sum of ascending and descending wave temperatures in
                 'TemperatureAscPlusDescSolMetGMMALinModelX.dat'
   - the ascending wave temperature in
                 'TemperatureAscSolMetGMMALinModelX.dat'
   - the desscending wave temperature in
                 'TemperatureDescSolMetGMMALinModelX.dat'

   - the pressure in 'WaveResultsSFSolMetGMMA/PressureLinModelX.dat'
   - the mass density in 'WaveResultsSFSolMetGMMA/MassDensityLinModelX.dat'

4. Wave parameters computed by the Global-Matrix Method for Nodal Values (GMMN)
   for a single-frequency wave versus altitude are written to the directory
   'WaveResultsSF' in the following files:
   - the horizontal velocity in
                 'UVelocitySolMetGMMNLinModelX.dat'
   - the vertical velocity in
                 'WVelocitySolMetGMMNLinModelX.dat'
   - the temperature in
                 'TemperatureSolMetGMMNLinModelX.dat'
   - the pressure in
                 'PressureSolMetGMMNLinModelX.dat'
   - the mass density in
                 'MassDensitySolMetGMMNLinModelX.dat'

5. Wave parameters computed by the Scattering-Matrix Method for Amplitudes (SMMA)
   for a single-frequency wave versus altitude are written to the directory
   'WaveResultsSF' in the following files:
   - the horizontal velocity in
                 'UVelocitySolMetSMMALinModelX.dat'
   - the vertical velocity in
                 'WVelocitySolMetSMMALinModelX.dat'
   - the temperature in
                 'TemperatureSolMetSMMALinModelX.dat'
   - the pressure in
                 'PressureSolMetSMMALinModelX.dat'
   - the mass density in
                 'MassDensitySolMetSMMALinModelX.dat'

6. Wave parameters for a time-dependent wave packet are written to the directory
  'WaveResultsWP' in the following files:
   - the horizontal velocity versus time and altitude in the Gnuplot-format file in
                 'UVelocityWP.dat'
   - the vertical velocity versus time and altitude in the Gnuplot-format file in
                 'WVelocityWP.dat'
   - the temperature versus time and altitude in the Gnuplot-format file in
                 'TemperatureWP.dat'
   - the source function versus time and altitude in the Gnuplot-format file in
                 'SrcFctWP.dat'
   - the plots generated with Gnuplot for velocities, temperature, and source
     function in 'UVelocityWP.png', 'WVelocityWP.png', 'TemperatureWP.png', and
     'SrcFctWP.png'

    - the maximum horizontal velocity over altitude versus time in
                 'UVelocityWPMax.dat'
    - the maximum vertical velocity over altitude versus time in
                 'WVelocityWPMax.dat'
    - the maximum temperature over altitude versus time in
                 'TemperatureWPMax.dat'
    - the maximum of the source function versus time in




II. THEORETICAL BACKGROUND AND NUMERICAL MODEL

1. Linearized hydrodynamic equations
=====================================
The model solves the linearized hydrodynamic equations for the neutral
atmosphere, including the effects of viscosity, thermal conduction and ion-drag.
The equations are derived under the following main assumptions:
- The geographic and geomagnetic coordinate systems are assumed to be identical.
- The gravity wave propagates in the meridional direction; the x-coordinate
  is positive southward, while the z-coordinate is positive upward.
- All background (unperturbed) quantities vary only with altitude z, whereas
  all perturbations vary harmonically in time and in the x-direction.

The ion-drag terms are treated in an approximate manner in order to decouple 
the hydrodynamic and ion equations. To achieve this, we adopt the following 
assumptions:
- In the ion continuity equation, the perturbed production and loss terms are
  neglected.
- In the ion momentum equation, ion inertia and ion–ion collisions are neglected,
  and only transport parallel to the magnetic field lines is retained. Under these
  assumptions, the ion momentum equation reduces to the ambipolar diffusion
  equation.
- To decouple the ion continuity equation from the diffusion equation, fast
  field-aligned diffusion is assumed, meaning that the field-aligned diffusion is
  sufficiently strong for the relative ion perturbation and the perturbed diffusion
  velocity to remain nearly constant along a magnetic field line.

The linearized equations are transformed into a linear system of ordinary
differential equations with variable coefficients that depend on the background
atmospheric parameters and their vertical derivatives.


2. Multilayer method
=====================
To integrate the linearized equations, a multilayer method is employed. In
this approach, the atmosphere is divided into a sequence of thin layers. Within
each layer, a linear system of ordinary differential equations with constant
coefficients is solved analytically. The wave solutions in neighboring layers are
then matched by imposing continuity conditions on the variables at the layer
interfaces.

In each layer,
- the dynamic (molecular) viscosity and the thermal conductivity are assumed to
  depend on temperature, and
- the vertical derivatives of the background temperature, wind velocity, and
  background dynamic viscosity are included in the formulation.
The second assumption implies the use of a numerical multilayer method, in
which the coefficients of the differential equations are approximated by their
values at the midpoint of each layer. As a consequence, the height derivatives of
the atmospheric parameters, evaluated at the layer midpoints, are explicitly
included in the resulting system of equations.


3. Linearization models
========================
We consider two linearization models:
- General model, that accounts on the altitude derivatives of the background
  velocity u0, temperature T0, density scale height HDens, and dynamic
  viscosity mu0,

- Simplified model for an isothermal (T0 = constant), homogeneous (kinematic
  viscosity = constant) and and windless atmosphere (u0 = constant) without
  ion drag.

The simplified model is used solely for testing purposes. In the code, the selected
linearization model is specified by the integer variable TypeLinModel, which takes
the following values:
- TypeLinModel = 1: General model
- TypeLinModel = 2: Simplified model

Note that the general model is consistent in form with the linearized model of
- Vadas and Nicolls (Vadas, S. L., and M. J. Nicolls (2012), The phases and
  amplitudes of gravity waves propagating and dissipating in the thermosphere:
  Theory, Journal of Geophysical Research: Space Physics, 117, A05322, 
  doi:10.1029/2011JA017426), and
- Knight et al. (Knight, H., Broutman, D., and Eckermann, S. (2024),
  Compressible and anelastic governing-equation solution methods for thermospheric
  gravity waves with realistic background parameters, Theoretical and Computational
  Fluid Dynamics, 38, 479–509, doi:10.1007/s00162-024-00709-x).

4. Solution methods
====================
The solution methods are based on the matrix exponential formalism and
encompass two main approaches:
- global matrix methods, and
- scattering matrix methods.
Both approaches aim to determine either
- the amplitudes of the characteristic solutions, or
- the nodal (grid-point) values of the state vector.
Ascending and descending gravity-wave modes are distinguished according to the
criterion that the real parts of the eigenvalues of the characteristic
equation for ascending modes are smaller than those for descending modes.

In global matrix methods, ascending and descending modes may be defined
- only at the upper and lower boundaries of the domain, or
- independently within each layer.
In contrast, scattering matrix methods require the explicit identification
of the mode type within each layer.

In the code, the selected solution method is specified by the integer variable
TypeSolMet, which takes the following values:
- TypeSolMet = 1: Global-Matrix Method for Amplitudes (GMMA)
- TypeSolMet = 2: Global-Matrix Method for Nodal (grid-point) Values (GMMN)
                  of the state vector
- TypeSolMet = 3: Scattering-Matrix Method for Amplitudes (SMMA)


5. Boundary values
===================
The state vector consists of the perturbed horizontal velocity, the perturbed
vertical velocity, the perturbed temperature, and their vertical derivatives.

At the lower boundary, we assume that only ascending wave modes transport
energy upward. Accordingly, we impose that a selected component of the state
vector is finite and that its first and second derivatives with respect to
height are zero. The component used for this boundary condition is specified by
the integer variable qNEV, which takes the following values:
- qNEV = 1: perturbed horizontal velocity
- qNEV = 2: perturbed vertical velocity
- qNEV = 3: perturbed temperature

The amplitudes of the perturbations at the lower boundary are prescribed as
follows:
- The amplitude of the perturbed horizontal velocity is set to
  epsU * UTMax, where epsU is an input parameter and UTMax is the maximum
  southward neutral horizontal velocity over the altitude range.

- The amplitude of the perturbed vertical velocity is specified by the input
  parameter w0BC.

- The amplitude of the perturbed temperature is set to
  epsT * TempMax, where epsT is an input parameter and TempMax is the
  maximum neutral temperature over the altitude range.

At the upper boundary, we assume that there is no downward energy propagation (the
amplitudes of all descending wave modes are set to zero).


6. Wave period
===============
To define practical lower and upper bounds for the wave period, the code evaluates
a simplified dispersion relation that neglects viscous and dissipative effects.
This inviscid formulation provides conservative estimates of the characteristic
wave frequencies and is used only to determine an admissible range of periods for
the simulations.

The dispersion relation is evaluated for two prescribed vertical wavelengths. 
The minimum vertical wavelength is set to LambdaZMin = 125 km, and the maximum vertical
wavelength is set to LambdaZMax = 250 km. These values are hard coded and chosen to
represent typical vertical scales of interest for the modeled wave dynamics.

From the corresponding solutions, the smallest and largest wave periods are obtained
by converting the maximum and minimum frequencies, respectively.

Because dissipative processes are neglected in this simplified formulation, the
resulting periods are expected to be underestimated. In realistic conditions,
viscosity and related effects reduce oscillation frequencies and therefore lead to
longer wave periods. To account for this, the computed minimum and maximum periods
are rounded upward to the nearest multiple of ten minutes. These rounded values are
then used as the final bounds for the wave period range in the simulations.


7. Source function
===================
The model handle in a first step a single-frequency wave, and in a second step
a time-dependent wave packet. In the latter case, the source function depends
explicitly on time, and the perturbed quantities are not monochromatic waves
with a single specified angular frequency. For time-dependent wave packets, the
governing equations are treated in the frequency domain by applying a Fourier
transform in time.

8. Frequency and time discretization for the Fourier transform
================================================================
In the code, we do not use a Fast Fourier Transform (FFT). Instead, we
perform a direct (discrete) Fourier transform by explicitly discretizing
the Fourier integral. The discretization parameters are chosen to resolve a
source that is localized both in frequency and in time.

8.1 Frequency band centered on Omega0
--------------------------------------
The transform is performed over a frequency band centered on the reference
frequency Omega0. We introduce a frequency standard deviation SigmaOmega
defined by
                   SigmaOmega = Omega0 / RatioSigmaOmega
The effective frequency band of interest is chosen to cover approximately
plus/minus three standard deviations of the source spectrum:
                   OmegaMin = Omega0 - 3.d0 * SigmaOmega
                   OmegaMax = Omega0 + 3.d0 * SigmaOmega
The total frequency interval and the frequency step are, respectively,
                   DeltaOmega = OmegaMax - OmegaMin
                   DOmega     = DeltaOmega / ( NFFT - 1 )
The discrete frequency grid OmegaK(i) is then constructed as
                   OmegaK(i) = OmegaMin + (i - 1) * DOmega,   i = 1,...,NFFT

8.2 Time interval and time step
--------------------------------
The time grid is chosen to be consistent with the assumed temporal localization
of the source and with the chosen frequency bandwidth. First, we define the
time standard deviation SigmaTime as
                   SigmaTime = 1.d0 / SigmaOmega
The time interval DeltaTime is chosen to contain NPeriod periods of the
reference frequency. The period is given by
                   PeriodTime = 2.d0 * Pi / Omega0
and therefore
                   DeltaTime = NPeriod * PeriodTime
The time step DTime is then defined as
                   DTime = DeltaTime / ( NFFT - 1 )
The time grid TimeK(i) is defined on the interval [TimeMin, TimeMax], chosen as
                   TimeMin = 0.d0
                   TimeMax = DeltaTime
by
                   TimeK(i) = TimeMin + (i - 1) * DTime,   i = 1,...,NFFT
We also define a time shift by
                   TimeShift = TimeK( NFFT / 2 )
which places the reference time close to the center of the time window.

8.3 Consistency condition
--------------------------
In the code, RatioSigmaOmega and NPeriod are input parameters. The discretization
is considered adequate if the following conditions are satisfied:
- The time window is long enough for the chosen frequency bandwidth. This requires
  that the number of time standard deviations contained in the interval satisfies
                    KSigma = DeltaTime / SigmaTime > 6,
  which ensures that the main part of the source is well contained within the time
  window and that truncation effects are negligible.
- The Nyquist condition
                   DTime < Pi / OmegaMax,
  and the dual Nyquist condition
                   DOmega < 2*Pi/DeltaTime
  are satisfied. These conditions ensure that the temporal signal is properly
  sampled in time and frequency, prevent aliasing effects, and guarantee a
  consistent discrete Fourier transform between the time and frequency domains.
In summary, the frequency grid OmegaK(i) is symmetric around Omega0, while
the time grid TimeK(i) covers the main part of the source in time, with TimeShift
located near the center of the window. This choice of parameters is appropriate
for a direct discrete Fourier transform and does not rely on FFT-specific
grid constraints.

8.4 Comment
------------
Since the Fourier transform is computed by direct discretization of the
Fourier integral, no FFT-specific constraint is imposed between the frequency
step DOmega and the time step DTime. In particular, the grid relation
                     DOmega = 2.d0 * Pi / ( NFFT * DTime )
-which is characteristic of discrete Fourier transforms based on periodic
sampling and FFT algorithms-is not required here. Instead, the frequency
and time grids are chosen independently, based on the desired frequency band
and time window needed to resolve a source that is localized in both domains.


9. Imaginary frequency shift
=============================
The imaginary frequency shift is determined using a heuristic approach that 
combines two criteria: application of the Layerwise Causality (LC) condition 
at selected altitude levels, and a Source-Function Reconstruction (SFR) test.

9.1 Layerwise causality condition at a selected set of altitude levels
-----------------------------------------------------------------------
First, a subset of altitude levels is chosen from the full altitude grid. Only
these selected altitudes are used in this step of the test. The total number of
selected altitude levels is recorded.

For each selected altitude level, the algorithm starts with an initial value of the
imaginary frequency shift. This value is then increased gradually in fixed steps,
but never beyond a prescribed upper limit. After each increase, the algorithm checks
whether the layerwise causality condition is satisfied at that altitude level for
all test frequencies within the specified frequency range. The check is performed
using a predefined tolerance.

If, for any selected altitude level, it is not possible to satisfy the layerwise
causality condition for all test frequencies by increasing the imaginary frequency
shift within the allowed interval, the subset-based layerwise causality test is 
declared to have failed.

If the causality condition can be satisfied at every selected altitude level, the
procedure collects the smallest admissible frequency-shift value found at each
altitude. The final frequency shift returned by the procedure is then chosen as the
largest of these values. This ensures that the causality condition is satisfied 
simultaneously at all selected altitude levels.

9.2. Source-function reconstruction test
-----------------------------------------
For a selected value of the frequency shift, the algorithm reconstructs the source
function in the time domain by applying an inverse Fourier transform to the
frequency-domain representation that has been shifted into the complex plane. This
reconstructed signal is then compared with the original source function.

To ensure that the comparison is independent of absolute amplitude, both the
original and the reconstructed source functions are first normalized by their
respective maximum real values. The algorithm then computes a relative
root-mean-square difference between these two normalized time signals over a
discrete set of time samples.

If the relative error is below a prescribed tolerance, the frequency shift is
considered valid, meaning that the reconstructed source function is sufficiently
close to the original one.

If the error exceeds the tolerance, the algorithm relaxes the frequency shift by
reducing its value in discrete steps and repeats the reconstruction and comparison.
This relaxation continues until either the error falls below the tolerance or a
predefined lower limit for the frequency shift is reached.

If the error remains above the tolerance even at this lower limit, the
source-function reconstruction test is declared to have failed, indicating that no
admissible frequency shift exists within the allowed interval.

9.3. Algorithm
---------------
The algorithm proceeds as follows.

- Step 1: Initialization and source-function reconstruction test at the minimum
  ------------------------------------------------------------------------------
  frequency shift.
  ----------------
  We first define a minimum and a maximum admissible value for the imaginary
  frequency shift, denoted by AbsOmegaImMin and AbsOmegaImMax, respectively.
  The minimum value is chosen equal to the initial discrete step DOmegaIm, and the 
  maximum value is obtained  by multiplying this step by a fixed integer factor.
  The discrete step DOmegaIm represents  the smallest increment by which the
  imaginary frequency shift is modified during the algorithm.

  As a first consistency check, we apply the source-function reconstruction test to
  the minimum frequency shift AbsOmegaImMin. If the source-function reconstruction
  test fails at this minimum value, the algorithm is stopped immediately, since no
  smaller admissible frequency shift is allowed. If the test is successful, the
  minimum frequency shift is accepted as admissible. Because the reconstruction
  error behaves monotonically with respect to the frequency shift, we assume that
  all smaller values would also satisfy the reconstruction test.

- Step 2: Source-function reconstruction test at the maximum frequency shift and
  -------------------------------------------------------------------------------
  refinement of the admissible interval
  --------------------------------------
  Next, we apply the source-function reconstruction test to the initially chosen
  maximum frequency shift AbsOmegaImMax. In this case, the test is allowed to
  decrease the frequency shift down to the previously accepted minimum value,
  using the same discrete step.

  Since the reconstruction test has already succeeded at the minimum frequency
  shift AbsOmegaImMin, this second test is guaranteed to succeed as well: if
  necessary, the procedure can always relax the frequency shift down to the minimum
  admissible value. The reconstruction routine then returns the largest value of
  the frequency shift AbsOmegaImMax for which the test remains satisfied.

  We replace the initial maximum frequency shift by this returned value. By
  construction, every frequency shift between the minimum and the updated maximum
  now satisfies the source-function reconstruction test.

- Step 3: Construction and adjustment of an internal frequency-shift grid
  ------------------------------------------------------------------------
  The next stage of the algorithm requires an internal set of discrete frequenc-shift 
  values between the current minimum and maximum. This internal grid must
  contain a sufficient number of points, at least 5, to allow a reliable
  application of the layerwise causality procedure.

  The internal grid uses the same basic step size as before, but it is independent
  of the number of values that will later be used for the final selection. To avoid
  a degenerate or poorly resolved interval, we adjust the minimum value and the step
  size according to the width of the admissible interval.

  If the interval is extremely narrow, meaning that the minimum and maximum are
  almost identical, we slightly reduce the minimum value and significantly refine
  the step size. If the interval width is comparable to one step, we reduce the step
  size moderately. If the interval spans roughly two steps, we also refine the step
  size, but less aggressively. After this adjustment, the interval contains a
  suitable number of internal grid points with a reasonable spacing.

  These internal grid points are used exclusively for the layerwise causality
  procedure and are not reused in the final selection step.

- Step 4: Layerwise causality procedure and update of the minimum frequency shift
  --------------------------------------------------------------------------------
  We now apply the layerwise causality procedure using the internal frequency-shift
  grid.

  At each discrete value, the layerwise causality condition is evaluated at a predefined 
  subset of altitude levels. These levels are defined by
             zShift(k) = zmin + (k-1) * dZShift,   k = 1,...,NZShift, where
  where NZShift and dZShift are input parameters. If there exists an altitude level for 
  which the LC test fails, we stop the algorithm. Otherwise, the minimum frequency shift 
  AbsOmegaImMin is updated to the value determined by the procedure.

  After this update, every frequency shift between the new minimum AbsOmegaImMin and
  the maximum AbsOmegaImMax satisfies both the source-function reconstruction test
  and the layerwise causality condition at the selected altitude levels.

- Step 5: Final selection of the imaginary frequency shift
  ---------------------------------------------------------
  The final step is based on a new, independent set of NOmegaIm evenly spaced
  frequency-shift values covering the interval between the updated minimum
  AbsOmegaImMin and maximum AbsOmegaImMax. The number of values used at this stage
  is chosen independently of the internal grid employed earlier.

  For each value in this final set, the algorithm performs three tasks. First, it
  computes the corresponding wave parameters. Second, it evaluates the maximum
  values of the perturbed horizontal velocity, the perturbed vertical velocity, and
  the perturbed temperature. Third, it checks whether the layerwise causality
  condition is satisfied over the entire altitude range, rather than only at
  selected levels.

  Finally, among all frequency shifts for which the causality condition is satisfied
  over the full altitude range, the algorithm determines the center of mass in the
  space spanned by the maximum values of the perturbed horizontal velocity, vertical
  velocity, and temperature, and selects the solution whose maximum-amplitude vector
  is closest to this center of mass.

A graphical illustration of the algorithm steps is shown below.
Note that, as an input parameter, we specify the absolute value of the imaginary frequency 
shift, AbsOmegaIm, whereas in the code the imaginary frequency shift OmegaIm is negative, 
that is, OmegaIm = −AbsOmegaIm.

- Step 1. Check whether AbsOmegaImMin = 1.0e−6 satisfies the source function
  reconstruction test. In the figure below, this value is considered to be valid.
        ________________________________________________________________
       |                      ___________________                       |
       |                     |                   |                      |
       o---------------------X---------o---------X----------------------o
   OmegaImMin            OmegaImMax     0    AbsOmegaImMin         AbsOmegaImMax
      -1.e-4              -1.e-6                1.e-6                 1.e-4

- Step 2. Apply the source function reconstruction test to determine AbsOmegaImMax.
  In the figure below, the initial value of AbsOmegaImMax is reduced from 1.e-4
  to 5.e-5.
             _____________________________________________________
            |                 ___________________                 |
            |                |                   |                |
       -----X----------------o---------o---------o----------------X-----
        OmegaImMin      OmegaImMax     0    AbsOmegaImMin    AbsOmegaImMax
          -5.e-5          -1.e-6                1.e-6           5.e-5

- Step 4. Apply the layerwise causality condition at specified altitude levels
  to determine AbsOmegaImMin. In the figure below, the value accepted in Step 1
  is increased from 1.e-6 to 5.e-6.
             ______________________________________________________
            |           _______________________________            |
            |          |                               |           |
       -----o----------X---------------o---------------X-----------o-----
        OmegaImMin  OmegaImMax         0         AbsOmegaImMin AbsOmegaImMax
          -5.e-5     -5.e-6                          5.e-6       5.e-5

- Step 5. Loop over the frequency shift and apply the Layerwise Causality Condition
  at all Altitudes (LCCA). In the figure below, only the solutions marked with X
  satisfy LCCA.
            _______________________________________________________
           |            _______________________________            |
           |LCCA |     |                               |    | LCCA |
       ----oXXXXXo-----o---------------o---------------o----oXXXXXXo-----
        OmegaImMin  OmegaImMax         0         AbsOmegaImMin AbsOmegaImMax
          -5.e-5     -5.e-6                          5.e-6       5.e-5




III. THIRD-PARTY SOFTWARE
This project includes or derives from the following third-party software:
1. LAPACK — 3-Clause BSD License
   (see LicenseLAPACK.txt)

2. sami2py — 3-Clause BSD License
   Copyright © 2021, sami2py development team
   (see LicenseSAMI2PY.txt)

3. Ooura FFT package
   Copyright © 1996–2001 T. Ooura
   (see LicenseOOURA.txt)
