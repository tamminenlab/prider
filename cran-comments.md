## Resubmission of an archived package

prider was archived on 2022-10-10 because its dependency 'blaster' was
archived. This submission resolves that by removing the 'blaster'
dependency entirely: the only function used from it (read_fasta) has
been replaced with a small internal base-R implementation. No other
user-facing behavior has changed.

## Test environments

* local: macOS 26.2 (aarch64), R 4.4.1
* GitHub Actions:
  - ubuntu-latest (R release, R oldrel-1, R devel)
  - macos-latest (R release)
  - windows-latest (R release)

## R CMD check results

0 errors | 0 warnings | 1 note

* The expected "New submission / Package was archived on CRAN" note
  from CRAN incoming feasibility.
