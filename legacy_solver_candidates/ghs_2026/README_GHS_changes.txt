GHS-enabled nozzle DSMC code package

Files:
- Viscous_Nozzle_GHS.for: modified Fortran source
- common.txt: modified common blocks with /GHSPAR/
- Property.txt: unchanged
- InputData.txt: modified input with IMOLMODEL=2 and N2 GHS parameters

GHS N2 parameters added from Hash, Moss & Hassan (1994):
 sigma = 3.798e-10 m
 epsilon/k = 71.4 K
 a1 = 2.32842
 a2 = 6.54333
 omega1 = 0.0280087
 omega2 = 0.777166

Implementation notes:
- IMOLMODEL=1: original VHS/VSS path
- IMOLMODEL=2: GHS total cross-section in SELECT via CVR = g*sigma_T(g)
- ELASTIC uses hard-sphere angular scattering for GHS, consistent with GHS model description
- PROPERTIES uses GHS effective cross-section for the mean-free-path/Kn diagnostic
- Initial CCG maximum collision-rate estimate uses a conservative 10x GHS estimate

Recommended first tests:
1) Run adiabatic with IMOLMODEL=1 and compare to old results to ensure no accidental changes when GHS is off.
2) Run adiabatic with IMOLMODEL=2 and compare history, mass flow, wall temperature, and Kn.
3) Then run Qw=2.5e4 and Qw=7.5e4 sensitivity cases.
