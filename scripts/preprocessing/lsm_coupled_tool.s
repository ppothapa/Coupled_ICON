#!/bin/ksh
#_____________________________________________________________________________
# Combine land-sea masks (LSM) for coupled ocean-atmosphere ICON runs:
#
# The LSM from the ocean is strictly logical (0/1) and has added and deleted ocean points.
# It is interpolated with conservative remapping onto the atmospheric grid and added
# to the extpar file as a fraction called cell_sea_land_mask.  The logical choices to adopt 
# the atmosphere to the ocean LSM are definded in ICON (mo_ext_data_init.f90).
#
# Info:
# - cell_sea_land_mask (in the ocean grid):
#   long_name = "sea (-2 inner, -1 boundary) land (2 inner, 1 boundary) mask for the cell"
# - flow chart of LSM in ICON:
#   * lsm_coupled_tool.s: interpolate ocean-grid LSM (cell_sea_land_mask) on atmo-grid
#                         and append to extparfile with same name.
#   * mo_ext_data_init/read_ext_data_atm: read cell_sea_land_mask from extpar file 
#                                         and put in variable ext_data%atmo%lsm_ctr_c
#   * mo_ext_data_init/lsm_ocean_atmo: convert with some rules lsm_ctr_c 
#                                      to fr_land and fr_lake
#   * mo_atmo_coupling_frame: give ext_data%atm%lsm_ctr_c to YAC representing the atmo LSM
#
# original version       Rene Redler and Martin Koehler - 2021-02
# prototyp2              Stephan Lorenz                 - 2022-07-19
# finalize, add to code  Martin Koehler                 - 2022-11-25
# using global grid      Stephan Lorenz                 - 2023-01-05
#_____________________________________________________________________________


set -ex

module load cdo/prerelease/2.1.0

WORKDIR=/hpc/uwork/mkoehler/run-icon/coupled/proto2/lsm-extpar

GRIDDIR=/hpc/rhome/routfox/routfox/icon/grids/public/edzw  # grid directory
EXTPDIR=/hpc/rhome/routfox/routfox/icon/grids/public/edzw  # external parameter directory
OCENDIR=/hpc/uwork/mkoehler/run-icon/coupled/proto2

#DWD path: /hpc/rhome/routfox/routfox/icon/grids/public/edzw
#MPI path: /pool/data/ICON/grids/public/mpim

extpname=20161124_tiles       # 20161124_tiles  20180625_tiles  20200917_tiles

atmos_gridID="0012"           # 0012            0030            0024
atmos_refinement="R02B04"     # R02B04          R02B05          R02B06

ocean_gridID="0043"           # 0043                            0035          
ocean_refinement="R02B04"     # R02B04                          R02B06


#_____________________________________________________________________________


mkdir -p $WORKDIR
cd $WORKDIR

ln -sf ${GRIDDIR}/icon_grid_${atmos_gridID}_${atmos_refinement}_G.nc               atmos_grid.nc
ln -sf ${OCENDIR}/icon_grid_${ocean_gridID}_${ocean_refinement}_O.nc               ocean_grid.nc #not needed
ln -sf ${OCENDIR}/icon_mask_${ocean_gridID}_${ocean_refinement}_G.nc               ocean_mask.nc
ln -sf ${EXTPDIR}/icon_extpar_${atmos_gridID}_${atmos_refinement}_G_${extpname}.nc extpar.nc

# note: ocean_grid only has coastal points and is not suitable for LSM interpolation
#       ocean_mask is global and used here for LSM calculation (NOT available for 0043 ?)

# conservative remapping from ocean to atmospheric grid with fractional LSM [0,1]
# - cdo calculates in 64bit

cdo selname,cell_sea_land_mask    ocean_mask.nc lsm_-2_2.nc          #  ocean mask -2. to 2.
cdo -b F64 gtc,0                  lsm_-2_2.nc   lsm_0-1.nc           #  ocean mask 0.0/1.0 (gtc,0: >0; 64bit))
cdo -b F64 remapcon,atmos_grid.nc lsm_0-1.nc    lsm_atmos.nc         #  atmos mask 0.0-1.0
#cdo setmisstoc,1.0               lsm_atmos.nc  lsm_atmos_no_miss.nc #  missing data -> land 
# - not needed for global ocean_mask 

# merge new LSM into original extpar file

cdo -O merge  lsm_atmos.nc  extpar.nc                                            lsm_temp.nc
ncatted -O -a rawdata,global,a,c,'GLOBCOVER2009, FAO DSMW, GLOBE, Lake Database' lsm_temp.nc

mv lsm_temp.nc icon_extpar_oceLSM_a${atmos_gridID}_${atmos_refinement}_o${ocean_gridID}_${ocean_refinement}_${extpname}.nc

\rm -f lsm_-2_2.nc lsm_0-1.nc lsm_atmos.nc lsm_atmos.nc lsm_atmos_no_miss.nc lsm_temp.nc


#_____________________________________________________________________________


# Stephan's comments: Interpolation of LSM

# create temporal mask on global ocean grid (real numbers):  0. for ocean, 1. for land points
#  - the default ocean_grid.nc is defined on ocean and coastal land points only, grid is not defined on inner land
#  - for this operation the LSM is needed on the global mask (ocean_mask.nc), which must match the geometry and the LSM
#    of the default ocean grid
#  - convert cell_sea_land_mask from integer to real (32bit or 64bit output, cdo calculates in 64bit)
#cdo -b F64 gtc,0 ocean_mask.nc temp_ocean_mask_01.nc

#  - if a LSM is only available on the default grid, it must be remapped to the global grid by remapnn:
# cdo -b F64 gtc,0 -remapnn,ocean_mask.nc  -selvar,cell_sea_land_mask ocean_grid.nc temp_ocean_mask_01.nc

# YAC first order conservative remapping from ocean to atmospheric fractional grid (64 bit)
#cdo -b F64 remapcon,atmos_grid.nc temp_ocean_mask_01.nc fractional_lsm_atm${atmos_gridID}_oce${ocean_gridID}_${revision}.nc
