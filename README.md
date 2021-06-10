# Njoy2ChiConverter
Prepares NJOY inputs and converts NJOY output to Chi-Tech Multigroup Transport cross-section format (.csx format).

## How to use this converter

### Step-0: Python requirements

Runs with python3.

### Step 1: Download and setup your flavor of ENDF

- The entire ENDF8.0 library (Preferred) [ENDFVIII.0](https://www.nndc.bnl.gov/endf/b8.0/zips/ENDF-B-VIII.0.zip) (~500Mb on download, ~2Gb extracted)
- The entire ENDF7.1 library [ENDFVII.1](https://ndclx4.bnl.gov/gf/download/frsrelease/138/2242/ENDF-B-VII.1.tar.gz)
- Older versions. Go to Brookhaven National Laboratory's site and download what you need. [National Nuclear Data Center](https://www.nndc.bnl.gov/exfor/endf00.jsp)

Extract this library in a folder of your choice which we shall just call `ENDF_FOLDER`.

For example, on Linux machines, you can do this using: 
- ```shell 
  curl -O https://www.nndc.bnl.gov/endf/b8.0/zips/ENDF-B-VIII.0.zip
  ```
- unzip using the ```unzip``` command
- set the envrionement variable `ENDF_ROOT` to point to `ENDF_FOLDER`. In bash, this is
  ```shell
  export ENDF_ROOT=<path-to-ENDF_FOLDER>
  ```
  **Note:** this is the absolute path.
  
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

Next add `njoy` to your environment variables. This may be different depending on your system. For `bash` users, in order to add the current `bin`-directory to your `PATH`-environment variable, type `pwd` and copy the path (let us use the place-holder `<path-to-njoy>` for this path). Next add the following to your `~/.bashrc` file:

```shell
export PATH=<path-to-njoy>:$PATH
```

### Step 3: Clone the converter to a folder of your choice

```shell
cd <folder-of-choice>
git clone https://github.com/Naktakala/Njoy2ChiConverter
cd Njoy2ChiConverter
```

### Step 4: Adapt one of the example scripts to your specific application

The example scripts are a collection of neutron only examples:
- `example1a.sh`, U235 with an xmas172 neutron group structure
- `example1b.sh`, same as 1a but using graphite, showcasing how to include S(a,b) thermal scattering
- `example1c.sh`, same as 1b but now with a custom weighting spectrum
- `example1d.sh`, same as 1b but with completely custom group structure

And a collection of neutron-gamma examples:
- `example2a`, O16 with the lanl30 neutron group structure and the lanl12 gamma group structure
- `example2b`, O16 with the xmas172 neutron group structure and the lanl48 gamma group structure


The basic input for each example is:
```shell
CWD=$PWD

#================================= Set properties here
export ENDF_ROOT=/Users/janv4/Desktop/Projects/ENDF/ENDF-B-VII.1
neutron_file="n-092_U_235.endf"

output_directory="../output/ENDF-B-VII-1/xmas172/"
output_file_prefix="U235"

#================================= Run NJOY
cd njoy_runner || exit

python generate_njoy_mgxs.py \
--path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
--temperature=293.6 \
--neutron_group_structure=22 \
--output_directory=$output_directory \
--output_filename=$output_file_prefix.njoy

cd "$CWD" || exit

#================================= Run converter
cd njoy_processor || exit

python njoy_processor.py \
--path_to_njoy_output=$output_directory/$output_file_prefix.njoy \
--output_file_path=$output_directory/$output_file_prefix.csx

cd "$CWD" || exit
```

## FAQs:

### FAQ-1: General format of input examples
The two python scripts are `generate_njoy_mgxs.py` and `njoy_processor.py`.

The `generate_njoy_mgxs.py` script basically runs NJOY and has the following inputs (most of which are optional):
```
# --path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
# --path_to_sab=$ENDF_ROOT/thermal_scatt/$sab_file \
# --inelastic_thermal_number=229 \
# --inelastic_thermal_num_atoms=1 \
# --elastic_thermal_number=230 \
# --path_to_gamma_endf= \
# --temperature=296.0 \
# --neutron_group_structure=22 \
# --neutron_weight_function=8 \
# --output_directory=$output_directory \
# --output_filename=$output_file_prefix.njoy \
# --gamma_group_structure=0 \
# --gamma_weight_function=2 \
# --custom_neutron_gs_file="" \
# --custom_gamma_gs_file="" \
# --custom_neutron_wt_file="" \
# --custom_gamma_wt_file="" \
```

The `njoy_processor.py` script converts NJOY output to Chi-cross-section format. It only has two required inputs:
```
cd njoy_processor

python njoy_processor.py \
--path_to_njoy_output=$output_directory/$output_file_prefix.njoy \
--output_file_path=$output_directory/$output_file_prefix.csx 
```

### FAQ-2: How to specify S(a,b) thermal scattering
Simply run `generate_njoy_mgxs.py` with `--inelastic_thermal_number=` and `--inelastic_thermal_num_atoms` options set. This can sometimes be tricky because sometimes one needs `--elastic_thermal_number` as well. Hydrogen in water is tricky in the fact that it requires `--inelastic_thermal_num_atoms=2`. 

There aren't a lot of materials that have S(a,b) inelastic treatment so it is worth doing some homework on them and verifying a semi-infinite medium spectrum like done in the `tests` folder.

### FAQ-3: How to produce neutron-gamma cross sections?
The moment you supply the option `--path_to_gamma_endf` then the script will know to run with gamma production and make `--gamma_group_structure` required.

### FAQ-4: Format of a custom weighting spectrum file
The file is in ENDF TAB1-record format which can be confusing. An example spectrum file is supplied in `njoy_automate2/spectrum_file.txt` and is the same for custom neutron AND gamma spectrums.

**Note:** remember the `/` terminator and the blank line at the end of the file.

### FAQ-5: Format of a custom group structure file
The first line of a group structure file is the total number of groups `G`. Then followed by the lowest energy cutoff then a total of `G` upper bin boundaries (all in eV):

```
# Number of groups G
6
# G+1 boundaries, lower bound first then all upper bounds (eV)
1.0e-5
5.0e-2
5.0e-1
1.0e2
1.0e5
1.0e6
20.0e6/

```

**Note:** remember the `/` terminator and the blank line at the end of the file.