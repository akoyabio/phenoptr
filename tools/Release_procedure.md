---
editor_options: 
  markdown: 
    wrap: 72
---

# Release procedure for phenoptr and phenoptrReports

Development work for the public release should be done on the `9000`
branch. Development work for the ABS whole-slide workflow is on the
`roi` branch. These notes are for a public release based on the `9000`
branch.

### New functions

New, exported functions must have doxygen comments including the
`@export` tag. The function should also be added to `_pkgdown.yml` so it
will appear in the Reference section of the documentation.

Unexported functions may also have doxygen comments. These should
include `@keyword internal` which keeps the help out of the package
index.

### Spellcheck

Spell-check with the command `spelling::spell_check_package()`. Fix any
errors and add any new words to the `inst/WORDLIST` file.

### Build help files

Build the help files using the Document menu item in the More menu of
the Build tab. This also rebuilds the `NAMESPACE` file.

### Run CMD Check

Use the Check button in the build tab to run R CMD Check. This performs
many integrity checks on the package. It should run with no test
failures and no errors or warnings. There are a few notes that can be
ignored:

-   Large installed package size
-   Many imports

Any other issues noted by CMD Check should be fixed. The most common
errors are missing package name qualifiers (e.g. using plain `mutate()`
instead of `dplyr::mutate()` and missing Imports.

One oddity is that CMD Check will give an error for symbols used in
tidyverse non-standard evaluation, such as column names used in `select`
and `filter` calls. The solution to this is to add the names to a call
to `utils::globalVariables`. See existing usage for examples.

### Update version number and release date

The package version number and release date should be updated in several
places:

-   The version number and date appear in the `DESCRIPTION` file and
    `NEWS.md`.
-   The phenoptr `README.md` also contains the version number in the
    citation section.

### Update NEWS.md

Make sure that `NEWS.md` contains information about all user-visible
changes in the release.

### Build package documentation

Run `pkgdown::build_site()` to rebuild the site documentation. Pay
attention to any warnings, for example about missing Reference entries.
Check the output in your web browser.

### Check in and merge

Check in any changes in the `9000` branch. Merge the `9000` branch into
`main`. Push both branches to GitHub.

### Create a tagged release

On GitHub, create a new release for the version. Create a tag for the
version. Copy the `NEWS.md` entry for the release into the release
notes.

### Create a new development version

Add the suffix `.9000` to the version number in `DESCRIPTION`. Create a
new entry in `NEWS.md` for the new development version.

### Merge with `roi` branch

Merge the new `9000` branch into the `roi` branch so it is up-to-date
with the changes. Push the `9000` and `roi` branches to GitHub.
