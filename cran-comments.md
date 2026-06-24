## Submission

This is a new submission of the `clusterIV` package.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

(The win-builder result is pending; confirm the line above matches the emailed
R-devel and R-release reports before submitting.)

## Test environments

* local: macOS (darwin), R 4.3.2 -- 0 errors, 0 warnings, only the New
  submission note (the PDF-manual and HTML-validation notes seen locally are
  environmental: no LaTeX / outdated HTML tidy on this machine)
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macOS-latest (R-release), windows-latest (R-release)
* win-builder: R-devel and R-release (submitted)

## Notes for the reviewer

* The DESCRIPTION and help pages contain terms that the spell checker may flag
  as possibly misspelled. These are intentional and used consistently:
  proper nouns (Frandsen, Leslie, McIntyre, Woodbury), the estimator acronyms
  (CJIVE, JIVE), and British spellings (residualise, optimise, partialled).
* The Description cites Frandsen, Leslie & McIntyre (2025), the paper whose
  estimator the package implements.
* The package depends only on base R (`Imports: stats`); there are no reverse
  dependencies (new package).
* All exported functions have documented return values (`\value`) and runnable
  examples that complete in well under five seconds and require no internet.
