# Spectral Database Modifications

_Last updated: 2026-06-05_

## 1. Added Lines

**M5-O3 peaks (column BT):** Added energies for elements Z = 83 (Bi) through Z = 96. These had been set to zero but are clearly visible for Bi in the spectrum, so the values were included.

**La M-series lines:** Added energies for the M1-N3, M3-O1, M3-O4, and M3-O5 lines of La, as these are clearly present.

**Ll and Ln lines:** These lines were missing for Mg, Al, Si, and P; their values were added from the xraydb database. The Ll weight was set to 1 and the Ln weight kept at 0, consistent with NIST. This distinction has negligible practical effect, as the line energies are nearly identical.

## 2. Added or Adjusted Line Weights

**M2-N1 lines:** Assigned a weight of 0.0005 to the M2-N1 line of Rb, Sr, Y, and Zr, based on its presence in Kr and Nb. The peak is clearly visible in Y.

**Ca Lα lines:** Added weights for the Ca Lα1 and Lα2 lines, setting them equal to those of Sc, which works perfectly.

**Mζ1, Mζ2 lines:** Doubled the weights for Z = 75–83, as these were uncalibrated. More precise calibration would be desirable.

**Mg lines:** Quickly adjusted weights for Nb, Mo, Ag, Cd, and In, none of which were calibrated by NIST.

**Y, Nb (preliminary adjustments):** Adjusted the weight of the Y Ll line and the Nb Lβ1 line, both of which were off despite being calibrated. These changes were not subjected to rigorous testing.

## 3. Energy Shifts

**Mn:** Red-shifted the Kβ1 peak energy by 8.1 eV to match the measured spectra.

## 4. Code-Compatibility Adjustments

**Li, Be:** Modified the Kα1 line weight from 0 to 1 and the Kα2 weight vice versa, for compatibility with the code, which expects a Kα1 line to be present for all elements.
**Pr:** Modified the Mα1 line weight from 0.008 to 1 and the Mα2 weight vice versa, for compatibility with the code, which expects a Mα1 line to be the reference for M lines when present.

## 5. Recalibrated Weights

Recalibrated weights for Si, Bi, Pb, V, and In, along with the following:

| #  | Line      | Note                                              |
|----|-----------|---------------------------------------------------|
| 21 | Sc Lα     | Calibrated by NIST, but not corresponding to signal |
| 22 | Ti Lα     | Calibrated by NIST, but not corresponding to signal |
| 25 | Mn Kβ     | Calibrated by NIST, but not corresponding to signal |
| 26 | Fe Lα     | Calibrated by NIST, but not corresponding to signal |
| 28 | Ni Lα     | Calibrated by NIST, but slightly overestimated     |
| 34 | Se Lα     | Calibrated by NIST, but slightly underestimated    |
| 56 | Ba Mζ     | Previously calibrated, but at half the actual signal |
| 57 | La Mζ     | Previously calibrated, but at 60% of the actual signal |
| 59 | Pr Mα     | Calibrated by NIST, but not corresponding to signal |
| 72 | Hf Mα + N | Mostly uncalibrated from NIST                       |
| 72 | Hf Lα     | Partially uncalibrated by NIST                      |
| 73 | Ta Mα + N | Mostly uncalibrated from NIST                       |
| 73 | Ta Lα     | Partially uncalibrated by NIST                      |
| 74 | W Lα      | Partially uncalibrated by NIST                      |