# Njoy2ChiConverter
Converts NJoy output to Chi-Tech Multigroup Transport cross-section format (.csx format).

## How to use this converter

### Step 1: Download and setup your flavor of ENDF

- The entire ENDF8.0 library (Preferred) [ENDFVIII.0](https://www.nndc.bnl.gov/endf/b8.0/zips/ENDF-B-VIII.0.zip) (~500Mb on download, ~2Gb extracted)
- The entire ENDF7.1 library [ENDFVII.1](https://ndclx4.bnl.gov/gf/download/frsrelease/138/2242/ENDF-B-VII.1.tar.gz)
- Older versions. Go to Brookhaven National Laboratory's site and download what you need. [National Nuclear Data Center](https://www.nndc.bnl.gov/exfor/endf00.jsp)

Extract this library in a folder of your choice which we shall just call `ENDF_FOLDER`.

### Step 2: Download and install NJOY2016

This converter processes only NJOY2016 output files for now so please stick to this version. The project is hosted on GitHub at [https://github.com/njoy/NJOY2016](https://github.com/njoy/NJOY2016), which has a link for installation instructions. However, the installation instructions are for NJOY2021 and therefore we give a short summary of the equivalent instructions for NJOY2016.

Make sure you meet the prerequisites:

- C++17 or higher (If you have GCC or clang then you should be good)
- Fortran 2003 or higher
- Python 3.4+
- CMake 3.2+

After this is good to go, do the following. Go to a folder where you want to install NJOY2016

```shell
# Download the source code
git clone https://github.com/njoy/NJOY2016

# Get the desired version of NJOY21 (1.1.0 in this example)
cd NJOY2016
wget https://raw.githubusercontent.com/njoy/signatures/master/NJOY21/1.1.0-NJOY21.json
./metaconfigure/fetch_subprojects.py 1.1.0-NJOY21.json
```
Yes we know the file points to NJOY21. We checked and it works for NJOY2016.

```shell
# Configure the build process
mkdir bin
cd bin
cmake -D fetched_subprojects=true ../

# Build NJOY16
make

# Test NJOY16
make test
```

