# Climate Data Integration for Modeling Building Energy Demand

R package **climbed**, version **1.1.3**

[![CRAN status](https://www.r-pkg.org/badges/version/climbed)](https://cran.r-project.org/package=climbed) [![R build status](https://github.com/hagento/climbed/workflows/check/badge.svg)](https://github.com/hagento/climbed/actions) [![codecov](https://codecov.io/gh/hagento/climbed/branch/master/graph/badge.svg)](https://app.codecov.io/gh/hagento/climbed) [![r-universe](https://pik-piam.r-universe.dev/badges/climbed)](https://pik-piam.r-universe.dev/builds)

## Purpose and Functionality

Prepare climate data to be used by building energy demand models.


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("climbed")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Questions / Problems

In case of questions / problems please contact Hagen Tockhorn <hagento@pik-potsdam.de>.

## Citation

To cite package **climbed** in publications use:

Tockhorn H (2025). "climbed: Climate Data Integration for Modeling Building Energy Demand - Version 1.1.3."

A BibTeX entry for LaTeX users is

 ```latex
@Misc{,
  title = {climbed: Climate Data Integration for Modeling Building Energy Demand - Version 1.1.3},
  author = {Hagen Tockhorn},
  date = {2025-06-11},
  year = {2025},
}
```
