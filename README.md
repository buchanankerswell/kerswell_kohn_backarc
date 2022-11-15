# Kerswell & Kohn (2023, G3)

![](figs/diff/CentralAmericaDiffComp.png)

This repository provides all materials for the manuscript *A Comparison of Heat Flow Interpolations near Subduction Zones* (Kerswell & Kohn 2023, G3)

This repository includes:

- All datasets required to compile the study
- R scripts to reproduce all results and figures
- A Makefile to compile the study

*This repository is self-contained but requires the following software (all open-source).*

## Prerequisite software

### R

This study is written in R. Follow the instructions at [R's homepage](https://www.r-project.org) to download and install the latest release of R on your machine.

### GDAL, GEOS, and PROJ

Geographic operations require the geographic libraries [gdal](https://gdal.org), [geos](https://trac.osgeo.org/geos), and [proj](https://proj.org). On Mac, the easiest way to get gdal, geos, and proj is to use Homebrew. If your machine does not have Homebrew, you can install it by running

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Or follow instructions from [Hombrew's webpage](https://brew.sh).

Once Homebrew is installed, the following will install the latest gdal, geos, and proj libraries together:

```
brew install pkg-config
brew install gdal
```

### sf

The trickiest R package to install is [sf](https://r-spatial.github.io/sf/), since it needs to be compiled from source and depends on gdal, geos, and proj. Installation instructions for mac are found [here](https://github.com/r-spatial/sf/issues/1536#issuecomment-727342736). For all other systems please see [sf's webpage](https://r-spatial.github.io/sf/).

## Running the study

```
# Clone this repository
git clone https://github.com/buchanankerswell/kerswell_kohn_backarc.git

# Change into the directory
cd kerswell_kohn_backarc

# Use Makefile to compile
make
```

This will check for required R packages and try to install missing packages automatically.

If all packages are found and available it will proceed to run the study with some initial prompts from the user.

# Open Science Framework

This repository can also be found at the official [OSF repository](https://osf.io/ca6zu/).

# Funding

This project was supported by the NSF grant OIA1545903 to M. Kohn, S. Penniston-Dorland, and M. Feineman