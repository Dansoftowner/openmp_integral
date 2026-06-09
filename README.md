# Numerical Approximation of the 2D Gaussian Integral

This project studies numerical integration of the two-dimensional Gaussian function
in C, using both sequential and OpenMP-parallel implementations. The goal is to
approximate the integral of

```text
f(x, y) = e^(-(x^2 + y^2))
```

with several numerical methods, and to compare sequential and parallel execution
times.

The known result over the full plane is:

```text
integral integral e^(-(x^2 + y^2)) dx dy = pi
```

The programs use this identity as a reference value. For the disk-limited version,
the corresponding closed-form analytical result is computed and used for error
measurement.

## Project Structure

```text
.
|-- gauss2d_integral_seq.c
|-- gauss2d_integral_par.c
|-- gauss2d_circle_integral_seq.c
|-- gauss2d_circle_integral_par.c
|-- gauss2d_montecarlo_integral_seq.c
|-- gauss2d_montecarlo_integral_par.c
|-- README.md
`-- README_hu.md
```

## Program Overview

| File | Type | Method | Parallelization |
| --- | --- | --- | --- |
| `gauss2d_integral_seq.c` | sequential | deterministic grid integration on a transformed domain | none |
| `gauss2d_integral_par.c` | parallel | deterministic grid integration on a transformed domain | OpenMP |
| `gauss2d_circle_integral_seq.c` | sequential | polar-coordinate integration over a disk | none |
| `gauss2d_circle_integral_par.c` | parallel | polar-coordinate integration over a disk | OpenMP |
| `gauss2d_montecarlo_integral_seq.c` | sequential | Monte Carlo sampling on a transformed domain | none |
| `gauss2d_montecarlo_integral_par.c` | parallel | Monte Carlo sampling on a transformed domain | OpenMP |

## Mathematical Background

### 1. Full-Plane Gaussian Integral

The two-dimensional Gaussian integral is:

```text
I = integral_{-inf}^{inf} integral_{-inf}^{inf} e^(-(x^2 + y^2)) dx dy
```

It can be separated into the product of two one-dimensional Gaussian integrals:

```text
I = (integral_{-inf}^{inf} e^(-x^2) dx)^2 = sqrt(pi)^2 = pi
```

Therefore, the full-plane programs compare the numerical result against `pi`.

### 2. Transforming the Infinite Domain

Directly applying a finite grid to the full plane is difficult because the domain is
infinite. The project uses the substitution:

```text
x = tan(u)
y = tan(v)
u, v in [-pi/2, pi/2]
```

The Jacobian terms are:

```text
dx/du = 1 / cos^2(u)
dy/dv = 1 / cos^2(v)
```

This gives the transformed integrand:

```text
g(u, v) = e^(-(tan^2(u) + tan^2(v))) / (cos^2(u) cos^2(v))
```

The `gauss2d_integral_*` and `gauss2d_montecarlo_integral_*` programs integrate this
transformed function over the finite domain
`[-pi/2, pi/2] x [-pi/2, pi/2]`.

### 3. Disk Integral in Polar Coordinates

The disk-based version integrates over a disk with radius `R`:

```text
I_R = integral_0^R integral_0^(2pi) e^(-r^2) r dphi dr
```

The `r` factor is the Jacobian term of the polar-coordinate transformation. The
analytical result is:

```text
I_R = pi * (1 - e^(-R^2))
```

The `gauss2d_circle_integral_*` programs print this analytical value and measure the
absolute numerical error against it.

## Numerical Methods

### Midpoint Rectangle Rule

The deterministic integration programs use the midpoint rectangle rule. The domain
is divided into equal parts, and the integrand is evaluated at the midpoint of each
cell.

For the transformed full-plane integral:

```text
dx = (b - a) / N
dy = (d - c) / M
```

where:

```text
a = c = -pi/2
b = d =  pi/2
```

The summation is:

```text
sum += g(a + (i + 0.5) * dx, c + (j + 0.5) * dy)
result = sum * dx * dy
```

For the disk version:

```text
dr = R / N_r
dphi = 2pi / N_phi
```

The summation is:

```text
sum += e^(-r_mid^2) * r_mid
result = sum * dr * dphi
```

### Monte Carlo Approximation

The Monte Carlo programs generate random points on the transformed
`[-pi/2, pi/2] x [-pi/2, pi/2]` domain. The integral is estimated as:

```text
result = area * average(g(u, v))
```

where the transformed domain area is:

```text
area = pi * pi
```

The advantage of Monte Carlo integration is that it is simple and naturally
parallelizable. Its drawback is statistical error: the result can vary between runs,
and the error usually decreases more slowly than with deterministic grid integration.

## Parallelization

The parallel programs use OpenMP.

### Grid Integration

The deterministic parallel programs parallelize the two nested loops with:

```c
#pragma omp parallel for collapse(2) reduction(+:sum)
```

Meaning:

- `parallel for`: distributes loop iterations among multiple threads;
- `collapse(2)`: flattens the two-dimensional loop nest into one larger iteration
  space;
- `reduction(+:sum)`: gives each thread a private partial sum and combines the
  partial sums safely at the end.

### Monte Carlo

The parallel Monte Carlo program starts a parallel region with:

```c
#pragma omp parallel reduction(+:sum)
```

Inside this region, the samples are distributed with:

```c
#pragma omp for
```

The sequential Monte Carlo program uses `rand()`, while the parallel version uses
`rand_r()`, because `rand()` has global state and is not a good choice for parallel
use.

## Requirements

Compilation requires a C compiler with support for:

- the C standard math library;
- `clock_gettime` and `CLOCK_MONOTONIC`;
- OpenMP for the parallel files.

On Linux or WSL, GCC is typically enough:

```bash
gcc --version
```

On Windows, WSL, MSYS2/MinGW, or another environment with OpenMP and the POSIX time
API is recommended.

## Compilation

### Sequential Programs

```bash
gcc -O2 gauss2d_integral_seq.c -lm -o gauss2d_integral_seq
gcc -O2 gauss2d_circle_integral_seq.c -lm -o gauss2d_circle_integral_seq
gcc -O2 gauss2d_montecarlo_integral_seq.c -lm -o gauss2d_montecarlo_integral_seq
```

### Parallel Programs

```bash
gcc -O2 -fopenmp gauss2d_integral_par.c -lm -o gauss2d_integral_par
gcc -O2 -fopenmp gauss2d_circle_integral_par.c -lm -o gauss2d_circle_integral_par
gcc -O2 -fopenmp gauss2d_montecarlo_integral_par.c -lm -o gauss2d_montecarlo_integral_par
```

Note: on some systems, the `-lm` flag should remain at the end of the command because
math-library linking can be order-dependent.

## Running

### Full-Plane Integral, Grid Method

Sequential:

```bash
./gauss2d_integral_seq <N> <M>
```

Parallel:

```bash
./gauss2d_integral_par <N> <M>
```

Example:

```bash
./gauss2d_integral_seq 1000 1000
./gauss2d_integral_par 1000 1000
```

Parameters:

- `N`: number of subdivisions in the first dimension;
- `M`: number of subdivisions in the second dimension.

### Disk Integral

Sequential:

```bash
./gauss2d_circle_integral_seq <R_radius> <N_r_subdivisions> <N_phi_subdivisions>
```

Parallel:

```bash
./gauss2d_circle_integral_par <R_radius> <N_r_subdivisions> <N_phi_subdivisions>
```

Example:

```bash
./gauss2d_circle_integral_seq 5 1000 1000
./gauss2d_circle_integral_par 5 1000 1000
```

Parameters:

- `R_radius`: radius of the disk;
- `N_r_subdivisions`: number of radial subdivisions;
- `N_phi_subdivisions`: number of angular subdivisions.

### Monte Carlo Integral

Sequential:

```bash
./gauss2d_montecarlo_integral_seq <N_samples>
```

Parallel:

```bash
./gauss2d_montecarlo_integral_par <N_samples>
```

Example:

```bash
./gauss2d_montecarlo_integral_seq 10000000
./gauss2d_montecarlo_integral_par 10000000
```

Parameter:

- `N_samples`: number of random sample points.

## Setting the Number of OpenMP Threads

The number of OpenMP threads can be controlled with the `OMP_NUM_THREADS`
environment variable.

Linux, WSL, or Git Bash:

```bash
export OMP_NUM_THREADS=4
./gauss2d_integral_par 2000 2000
```

PowerShell:

```powershell
$env:OMP_NUM_THREADS = "4"
.\gauss2d_integral_par.exe 2000 2000
```

## Output Interpretation

The programs usually print:

```text
Result: ...
Error: ...
Time: ...
```

The actual source files currently print Hungarian labels. The disk-based programs
also print:

```text
Analytical: ...
```

Meaning:

- `Result`: the numerically computed integral;
- `Analytical`: the closed-form reference value for the disk integral;
- `Error`: the absolute difference from the reference value;
- `Time`: measured execution time in seconds.

For the full-plane and Monte Carlo versions, the reference value is `pi`, so the
error is:

```text
|pi - result|
```

For the disk version:

```text
|pi * (1 - e^(-R^2)) - result|
```

## Performance Measurement Suggestion

To compare sequential and parallel versions, use the same input parameters.

Example measurement series for the full-plane grid integral:

```bash
./gauss2d_integral_seq 1000 1000
OMP_NUM_THREADS=2 ./gauss2d_integral_par 1000 1000
OMP_NUM_THREADS=4 ./gauss2d_integral_par 1000 1000
OMP_NUM_THREADS=8 ./gauss2d_integral_par 1000 1000
```

Useful metrics:

```text
speedup = T_sequential / T_parallel
efficiency = speedup / thread_count
```

For more reliable comparisons, run each measurement several times and use the
average execution time, because operating-system load and CPU state can affect the
results.

## Accuracy Considerations

- Larger `N`, `M`, `N_r`, and `N_phi` values usually improve deterministic
  integration accuracy, but increase runtime.
- In Monte Carlo integration, a larger sample count reduces statistical error, but
  the error can still vary from run to run.
- The transformed-domain method maps an infinite domain to a finite interval, but
  the `tan` and `cos` functions make the endpoints numerically sensitive. The
  programs use midpoint sampling, so they do not evaluate the singular endpoints
  directly.
- The programs use `long double` for summation and key floating-point calculations,
  which can improve numerical accuracy compared with plain `double`.

## Important Implementation Details

### Timing

All programs use `clock_gettime(CLOCK_MONOTONIC, ...)`. This uses a monotonic clock,
so it is not affected by manual or network-based system time changes.

### Mathematical Constants

The code defines its own `long double` precision pi constants:

```c
#define M_PIl   3.141592653589793238462643383279502884L
#define M_PI_2l 1.570796326794896619231321691639751442L
```

This is more portable than relying on whether the system math headers define the
`M_PI` macro.

### Reduction

In parallel summation, shared writes to `sum` would cause a race condition. The
OpenMP `reduction(+:sum)` clause solves this by giving each thread its own partial
sum and combining those partial sums at the end of the parallel region.

### Random Number Generation

The sequential Monte Carlo program uses:

```c
srand(time(NULL));
rand();
```

The parallel version initializes a per-thread seed:

```c
unsigned int seed = time(NULL) + omp_get_thread_num();
rand_r(&seed);
```

This reduces problems caused by global random-generator state. For scientific Monte
Carlo simulations, a higher-quality and reproducible random-number generator would
be a useful improvement.

## Known Limitations and Possible Improvements

- Input validation is minimal. Negative or zero subdivision counts can produce
  incorrect results or runtime errors.
- There is no shared build system, such as a `Makefile` or CMake configuration.
- Some Hungarian accented output strings in the source files currently appear to
  have incorrect character encoding.
- The parallel Monte Carlo run is not fully reproducible because it initializes
  random seeds from the current time.
- The `gauss2d` helper function remains in some files, but in the optimized disk
  version it is not actively used because `e^(-r^2)` is computed directly.

Possible improvements:

- add a `Makefile` for compiling all programs with one command;
- add input validation and clearer error messages;
- add a measurement script for testing multiple thread counts and problem sizes;
- add an optional seed parameter for reproducible Monte Carlo runs;
- write results to CSV for performance comparison;
- use a higher-quality random-number generator for the Monte Carlo versions.

## Suggested Report Structure

If a separate written report or presentation is prepared for this project, the
following structure can be used:

1. Problem statement
2. Mathematical background: 2D Gaussian integral and analytical results
3. Numerical methods: midpoint rectangle rule and Monte Carlo
4. Parallelization with OpenMP
5. Compilation and execution
6. Measurement results
7. Speedup and efficiency analysis
8. Accuracy and error analysis
9. Conclusions and possible improvements

