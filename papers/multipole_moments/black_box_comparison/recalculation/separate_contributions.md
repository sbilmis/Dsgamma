# Separate contributions with lattice normalization

Each entry is the median of that contribution over the accepted windows. Perturbative includes heavy and strange lines. Twist 3 includes two- and three-particle terms; twist 4 includes two-particle, gluonic three-particle and electromagnetic three-particle terms. Contact is the Wilson-line term.

The +chi sign uses f_perp = -55.1 MeV; -chi flips only twist 2. Chi does not enter the two-point function. The condensate is unchanged elsewhere.

Magnetic contributions are in units of mu_N. Quadrupole contributions are to F_D-Q_H, where Q_H=1 for Ds and 0 for Bs. Medians of parts need not sum to the median total; per-window additivity is checked in the frozen output.

## Nominal mixing angles

### Magnetic moment: mu/mu_N

| State | Sign | Perturbative | Twist 2 | Twist 3 | Twist 4 | Contact |
|---|---|---:|---:|---:|---:|---:|
| Ds1(2460) |  +chi | 0.780010 | -0.568491 | -0.094003 | -0.069647 | -0.015260 |
| Ds1(2460) | -chi | 0.780010 | 0.568491 | -0.094003 | -0.069647 | -0.015260 |
| Ds1(2536) |  +chi | 0.638959 | 0.051564 | -0.337387 | -0.156162 | 0.097233 |
| Ds1(2536) | -chi | 0.638959 | -0.051564 | -0.337387 | -0.156162 | 0.097233 |
| Bs1(5750) |  +chi | 0.133431 | -0.314355 | -0.017663 | 0.017337 | 0.000000 |
| Bs1(5750) | -chi | 0.133431 | 0.314355 | -0.017663 | 0.017337 | 0.000000 |
| Bs1(5830) |  +chi | 0.480480 | 0.023798 | -0.329725 | 0.206843 | 0.000000 |
| Bs1(5830) | -chi | 0.480480 | -0.023798 | -0.329725 | 0.206843 | 0.000000 |

### Quadrupole form factor: F_D-Q_H

| State | Sign | Perturbative | Twist 2 | Twist 3 | Twist 4 | Contact |
|---|---|---:|---:|---:|---:|---:|
| Ds1(2460) |  +chi | -2.053195 | -2.643512 | 0.590327 | -0.154038 | 0.701903 |
| Ds1(2460) | -chi | -2.053195 | 2.643512 | 0.590327 | -0.154038 | 0.701903 |
| Ds1(2536) |  +chi | -2.216542 | 5.355061 | 1.791186 | -0.072485 | -1.020116 |
| Ds1(2536) | -chi | -2.216542 | -5.355061 | 1.791186 | -0.072485 | -1.020116 |
| Bs1(5750) |  +chi | -0.147135 | -1.432447 | -0.387724 | 0.266620 | 0.000000 |
| Bs1(5750) | -chi | -0.147135 | 1.432447 | -0.387724 | 0.266620 | 0.000000 |
| Bs1(5830) |  +chi | -5.728128 | 18.334977 | 8.001665 | -8.241479 | 0.000000 |
| Bs1(5830) | -chi | -5.728128 | -18.334977 | 8.001665 | -8.241479 | 0.000000 |

## 90 degrees, refitted pure-current windows

### Magnetic moment: mu/mu_N

| State | Sign | Perturbative | Twist 2 | Twist 3 | Twist 4 | Contact |
|---|---|---:|---:|---:|---:|---:|
| Ds1(2460) |  +chi | 0.839351 | -0.476423 | 0.026448 | -0.057200 | -0.002975 |
| Ds1(2460) | -chi | 0.839351 | 0.476423 | 0.026448 | -0.057200 | -0.002975 |
| Ds1(2536) |  +chi | 0.707351 | -0.470276 | -0.188986 | -0.089902 | 0.000371 |
| Ds1(2536) | -chi | 0.707351 | 0.470276 | -0.188986 | -0.089902 | 0.000371 |
| Bs1(5750) |  +chi | 0.160582 | -0.294854 | -0.019698 | 0.022629 | 0.000000 |
| Bs1(5750) | -chi | 0.160582 | 0.294854 | -0.019698 | 0.022629 | 0.000000 |
| Bs1(5830) |  +chi | 0.114196 | -0.260607 | -0.020594 | 0.012360 | 0.000000 |
| Bs1(5830) | -chi | 0.114196 | 0.260607 | -0.020594 | 0.012360 | 0.000000 |

### Quadrupole form factor: F_D-Q_H

| State | Sign | Perturbative | Twist 2 | Twist 3 | Twist 4 | Contact |
|---|---|---:|---:|---:|---:|---:|
| Ds1(2460) |  +chi | -0.283352 | 1.248852 | -1.329886 | 0.005127 | 0.026241 |
| Ds1(2460) | -chi | -0.283352 | -1.248852 | -1.329886 | 0.005127 | 0.026241 |
| Ds1(2536) |  +chi | -2.984153 | -3.452141 | 1.946491 | -0.256859 | 0.786012 |
| Ds1(2536) | -chi | -2.984153 | 3.452141 | 1.946491 | -0.256859 | 0.786012 |
| Bs1(5750) |  +chi | 1.032532 | 1.807016 | -0.828522 | -0.207140 | 0.000000 |
| Bs1(5750) | -chi | 1.032532 | -1.807016 | -0.828522 | -0.207140 | 0.000000 |
| Bs1(5830) |  +chi | -1.204771 | -2.942746 | 0.277191 | 0.343996 | 0.000000 |
| Bs1(5830) | -chi | -1.204771 | 2.942746 | 0.277191 | 0.343996 | 0.000000 |

Source: [committed summary CSV](blind_lattice_twist_split/twist_split_summary.csv), frozen in commit 54e6958.
