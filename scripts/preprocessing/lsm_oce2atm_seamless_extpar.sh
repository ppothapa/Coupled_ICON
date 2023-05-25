#!/usr/bin/env bash
#_____________________________________________________________________________
# For coupled ocean-atmosphere ICON-seamless runs: 
#  Create fractional land-sea-mask (LSM) for atmospheric component strictly based
#  on logical (0/1) land-sea-mask of ocean component of coupled model.
#
# The LSM from the ICON-O ocean that have added and deleted ocean points of
# LSM of the atmosphere is added to the extpar file in this script.
# The atmospheric LSM for ICON-NWP is created by interpolating the global
# ocean grid onto the atmospheric grid using conservative remapping.
# The logical choices to adopt the atmosphere to the ocean LSM are
# defined in ICON (mo_ext_data_init.f90).
#
# Info:
# cell_sea_land_mask:
#   long_name = "sea (-2 inner, -1 boundary) land (2 inner, 1 boundary) mask for the cell"
#
# prototyp2              Stephan Lorenz                 - 2022-07-19
# finalize, add to code  Martin Koehler                 - 2022-11-25
# version for levante    Stephan Lorenz                 - 2022-01-05
#_____________________________________________________________________________


#  Revision Number

revision=r0003

#  Names

WORKDIR=/work/mh0287/m211032/Icon/Coupled/Create_fractional_lsm_0012_0035

atmos_gridID="0012"
atmos_refinement="R02B04"

ocean_gridID="0035"
ocean_refinement="R02B06"

extpname=20161124_tiles

#_____________________________________________________________________________

#  load necessary and defined modules:

set +x
source /sw/etc/profile.levante
set -ex

module load cdo/2.0.6-gcc-11.2.0
module load nco/5.0.6-gcc-11.2.0
#module load nco


mkdir -p $WORKDIR
cd $WORKDIR

AGRIDIR=/pool/data/ICON/grids/public/edzw                  #  atmospheric NWP grid and extpar directory
OGRIDIR=/pool/data/ICON/grids/public/mpim/${ocean_gridID}  #  oceanic MPI grid directory

ln -sf ${AGRIDIR}/icon_grid_${atmos_gridID}_${atmos_refinement}_G.nc                        atmos_grid.nc
ln -sf ${OGRIDIR}/icon_mask_${ocean_gridID}_${ocean_refinement}_G.nc                        ocean_mask.nc   #  global mask incl. landpoints
ln -sf ${AGRIDIR}/icon_extpar_${atmos_gridID}_${atmos_refinement}_G_${extpname}.nc          extpar.nc
#ln -sf ${OGRIDIR}/icon_grid_${ocean_gridID}_${ocean_refinement}_O.nc                        ocean_grid.nc


# Interpolation of LSM:

# create temporal mask on global ocean grid (real numbers):  0. for ocean, 1. for land points
#  - the default ocean_grid.nc is defined on ocean and coastal land points only, grid is not defined on inner land
#  - for this operation the LSM is needed on the global mask (ocean_mask.nc), which must match the geometry and the LSM
#    of the default ocean grid
#  - convert cell_sea_land_mask from integer to real (32bit or 64bit output, cdo calculates in 64bit)
cdo -b F64 gtc,0 ocean_mask.nc temp_ocean_mask_01.nc

#  - if a LSM is only available on the default grid, it must be remapped to the global grid by remapnn:
# cdo -b F64 gtc,0 -remapnn,ocean_mask.nc  -selvar,cell_sea_land_mask ocean_grid.nc temp_ocean_mask_01.nc

# YAC first order conservative remapping from ocean to atmospheric fractional grid (64 bit)
cdo -b F64 remapcon,atmos_grid.nc temp_ocean_mask_01.nc fractional_lsm_atm${atmos_gridID}_oce${ocean_gridID}_${revision}.nc

# merge new LSM into original extpar file

cdo -O merge fractional_lsm_atm${atmos_gridID}_oce${ocean_gridID}_${revision}.nc extpar.nc temp.nc
ncatted -O -a rawdata,global,a,c,'GLOBCOVER2009, FAO DSMW, GLOBE, Lake Database' temp.nc

mv temp.nc icon_extpar_withOceanLSM_atm${atmos_gridID}_${atmos_refinement}_oce${ocean_gridID}_${ocean_refinement}_G_${extpname}_${revision}.nc

# cleanup
rm temp_ocean_mask_01.nc

#_____________________________________________________________________________

exit

