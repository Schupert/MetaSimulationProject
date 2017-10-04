#### Test installing packages for linux

getDependencies <- function(packs){
  dependencyNames <- unlist(
    tools::package_dependencies(packages = packs, db = available.packages(), 
                                which = c("Depends", "Imports"),
                                recursive = TRUE))
  packageNames <- union(packs, dependencyNames)
  packageNames
}
# Calculate dependencies
packages <- getDependencies(c("copula", "doParallel", "data.table", "foreach", "doRNG", "compiler", "metafor"))

# library(data.table)
# library(doParallel)
# library(foreach)
# library(doRNG)
# library(copula)
# library(compiler)
# library(metafor)

### wd =  "/panfs/panasas01/sscm/pb17951"

# Warning in install.packages(pkgFilenames, repos = NULL, type = "source") :
#   'lib = "/cm/shared/languages/R-3.0.2/lib64/R/library"' is not writable
# Would you like to use a personal library instead?  (y/n)

# Download the packages to the working directory.
# Package names and filenames are returned in a matrix.
setwd("D:/my_usb/packages/")
pkgInfo <- download.packages(pkgs = packages, destdir = getwd(), type = "source")
# Save just the package file names (basename() strips off the full paths leaving just the filename)
write.csv(file = "pkgFilenames.csv", basename(pkgInfo[, 2]), row.names = FALSE)


#### On site

# Set working directory to the location of the package files
setwd("D:/my_usb/packages/")

# Read the package filenames and install
pkgFilenames <- read.csv("pkgFilenames.csv", stringsAsFactors = FALSE)[, 1]
install.packages(pkgFilenames, repos = NULL, type = "source")

############################# Testing for Windows

#### Test installing packages 

getDependencies <- function(packs){
  dependencyNames <- unlist(
    tools::package_dependencies(packages = packs, db = available.packages(), 
                                which = c("Depends", "Imports"),
                                recursive = TRUE))
  packageNames <- union(packs, dependencyNames)
  packageNames
}
# Calculate dependencies
packages <- getDependencies(c("copula", "doParallel", "data.table", "foreach", "doRNG", "compiler", "metafor"))

# library(data.table)
# library(doParallel)
# library(foreach)
# library(doRNG)
# library(copula)
# library(compiler)
# library(metafor)

# Download the packages to the working directory.
# Package names and filenames are returned in a matrix.
setwd("D:/my_usb/packages/")
pkgInfo <- download.packages(pkgs = packages, destdir = getwd(), type = "win.binary")
# Save just the package file names (basename() strips off the full paths leaving just the filename)
write.csv(file = "pkgFilenames.csv", basename(pkgInfo[, 2]), row.names = FALSE)


#### On site

# Set working directory to the location of the package files
setwd("D:/my_usb/packages/")

# Read the package filenames and install
pkgFilenames <- read.csv("pkgFilenames.csv", stringsAsFactors = FALSE)[, 1]
install.packages(pkgFilenames, repos = NULL, type = "win.binary")
