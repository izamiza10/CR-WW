# PYTHIA configuration file.
# Generated on Mon Oct 27 16:13:24 EDT 2025 with the following command:
# ./configure --prefix=/Users/darrenduan/work/pythia8315 --with-hepmc3=/Users/darrenduan/work/rivet/local

# Install directory prefixes.
PREFIX_BIN=/Users/darrenduan/work/pythia8315/bin
PREFIX_INCLUDE=/Users/darrenduan/work/pythia8315/include
PREFIX_LIB=/Users/darrenduan/work/pythia8315/lib
PREFIX_SHARE=/Users/darrenduan/work/pythia8315/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_DTAGS=
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=/Users/darrenduan/work/rivet/local/bin/
HEPMC3_INCLUDE=-I/Users/darrenduan/work/rivet/local/include
HEPMC3_LIB=-L/Users/darrenduan/work/rivet/local/lib -Wl,-rpath,/Users/darrenduan/work/rivet/local/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

YODA_USE=false
YODA_CONFIG=
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

ROOT_USE=false
ROOT_CONFIG=
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
