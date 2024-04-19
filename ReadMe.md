# FABM-VIMS

This is a collection of [FABM](https://fabm.net) ports of the ICM model and the CoSINE model.

## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=vims -DFABM_VIMS_BASE=</path/to/fabm-vims>`

Here, `</path/to/fabm-vims>` is the directory with the FABM-VIMS code (the same directory that contains this ReadMe file). Note that `-DFABM_INSTITUTES=vims` will make FABM compile ICM and CoSINE as the *only* available biogeochemical models. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="vims;iow"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

For instance, to use VIMS models with the latest stable release of the Semi-implicit Cross-scale Hydroscience Integrated System Model (SCHISM), do the following:

```
git clone --recurse-submodules --depth 1 https://github.com/schism-dev/schism.git clone --depth 1 https://github.com/fabm-model/fabm.git
git clone --depth 1 https://github.com/schism-dev/fabm-vims.git
wget https://raw.githubusercontent.com/platipodium/fabm/patch-1/src/drivers/schism/fabm_driver.h -O fabm/src/drivers/schism/fabm_driver.h
mkdir build
cd build
cmake ../schism/src -DBLD_STANDALONE=ON -DUSE_FABM=ON -DFABM_BASE=../fabm -DFABM_INSTITUTES=vims -DFABM_VIMS_BASE=../fabm-vims
make pschism
```

This will create the `pschism` executable with support for FABM ICM and COSINE.

## How to run a FABM ICM or FABM CoSINE simulation

Two sample `fabm.yaml` files are provided: `fabm-vims-cosine.yaml` and `fabm-vims-icm.yaml` under `</path/to/fabm-vims>/testcases`. You can drop this file (and link to or rename it `fabm.yaml`) in the working directory of a FABM-compatible model such as SCHISM to use it during simulation.

## Report when it is working/not working for you!

Build test with build chain suggested above

| System | Type     | Compiler | Status |
|--------|----------|----------|--------|
| quoll  | arm64   | gfortran |        |
| femto  | x86-64   | ifort    |   yes   |
| kuro   | x86-64   | intel/openmpi    | yes     |
| strand | x86-64   | ifort    |        |
