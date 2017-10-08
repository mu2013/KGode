## Test environments
* local OS X install, R 3.3.0
* ubuntu 12.04 (on travis-ci), R 3.3.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 3 note

* This is a new release.
* checking R code for possible problems ... [6s] NOTE
third: no visible global function definition for 'tail'
Undefined global functions or variables:
  tail
Consider adding
  importFrom("utils", "tail")
to your NAMESPACE file.
** running examples for arch 'i386' ... [457s] NOTE
Examples with CPU or elapsed time > 10s
              user system elapsed
warpfun     253.44   0.00  253.55
warpInitLen  84.10   0.00   84.15
third        46.20   0.01   46.22
crossv       44.15   0.00   44.21
rkg          20.25   0.00   20.19
** running examples for arch 'x64' ... [566s] NOTE
Examples with CPU or elapsed time > 10s
              user system elapsed
warpfun     309.80   0.00  309.93
warpInitLen 108.76   0.00  108.80
third        56.19   0.01   56.24
crossv       55.87   0.03   55.94
rkg          23.79   0.00   23.78
rkhs         10.41   0.00   10.40

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

