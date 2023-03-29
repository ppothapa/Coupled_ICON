# Release notes for icon-2.6.6

Later than expected: this is another ICON release

- OpenACC port for Nvidia GPUs and enabling HPE-Cray OpenACC use with AMD GPUs
  - almost all of the ICON atmospheric components are ready for use on GPU equiped machines. This has been a big effort by Nvidia, HPE/Cray, MeteoSwiss, DWD, CSCS, CSC, and MPIM (the order is arbitrary).
  - there is further progress coming up
- CSCS
  - preparation of daint/dom transition for Nvidia compiler bumb to 22.5.
  - added Docker files for EXCLAIM project in scripts/docker/exclaim
  - more test cases
- DWD
  - still persisting problem with OpenMP (non-reproducability) in buildbot tests.
  - bufixes and improvements in the data assimilation part and physics
  - two-moment cloud microphysics has been replaced by a new version
  - initial steps for wave model implementation
  - added jsbach/vdiff from former ICON echam-physics
  - enable/fix more problems with respect to run icon-seamless 
  - some code moves and consolidation
  - OpenACC beautification application running to allow passing a required CI test
  - more test cases
- MPIM
  - bugfixes in the ocean
  - added additional diagnostics
  - improvements in JSBACH (HD mask, LCC, former dynveg)
  - updating of build on levante and macos
  - fixes for claw
  - fixes for building ICON with some data assimilation components
  - Nvidia 22.9 is needed for JSBACH
  - ocean updated to zstar vertical coordinate
  - bump rte-rrtmgp to new version
- DKRZ
  - fix many errors and enabled new features with respect to coupling external components
  - new and updated config-wrappers
  - improve buildbot testing

For details take a look at:

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

# Release notes for icon-2.6.5

It is time for yet another ICON release.

After 10 month we managed to provide a new ICON release with new
machine configurations, a new buildbot, significantly extended,
improved GPU capabilities, first steps into ICON seamless and a
cleaned up version of the MPIM AES physics (now available as aes
physics and not echam anymore.

- CSCS: Introduction of more and optimized OpenACC code parts and
        introducing work to allow later more modularization for
        improving and consolidating of the ICON code.

- DKRZ: Bugfixes and more bugfixes, cdi-1.8.4 incorporates nom (joint
        effort of DKRZ and CIMD) cdi-pio, new versions of all tools
        maintained by DKRZ, a new buildbot (Ralf Müller) and strong
        support for getting ICON running on levante.

- DWD: Introduction of many improvements in the parameterizations, in
       tuning and adding bugfixes. Setting up on the new ECMWF system
       and supporting the data assimilation system. First steps for
       ICON-seamless are done, especially the coupling to ICON ocean
       is technically available and most of the time varying external
       data for climate simulations have been added. However, the code
       variant has to be tuned and extensively tested and some parts are
       still missing. This activity is a joint effort of DWD and MPIM
       colleagues.

- MCH: Is actively providing more, and more model components of the
       NWP physics for the use on GPUs. Furthermore they are
       supporting all activities related to GPUs and testing of large
       parts of the ICON code.

- MPIM AES: To simplify and streamline further high resolution
       modelling only, the echam physics has been renamed as aes
       physics and some former components have been removed as are the
       SSO gravity drag parameterization, as well the gravity drag
       parameterization for the middle atmosphere, and most prominent
       the convection parameterization. As cloud paramterizations the
       DWD graupel and two-moment scheme are available. A few more,
       minor, changes have been applied. Newly added is the former
       large eddy model component the Smagorinsky 'vertical' diffusion
       scheme.

- MPIM OES: Added a new C-grid seaice model (C. Mehlmann) and extensive
       improvements of OpenMP for hamocc. 

- MPIM LES: JSBACH got more model components and has been adapted for
       to all changes above its top-level.  - KIT: The aerosoly
       microphysics and plume components have been added as well as a
       large number of bugfixes.

- MPIM CIMD: just a lot of things in the area of coupling, configuring
       and building and model infrastrucuture, new machine setups
       for lumi (CSC), juwels/booster (JSC) and levante (DKRZ) 

- All developers: many of the work has been done in cross
  institutional working groups.


A special thanks is going to Sergey Kosukhin (MPIM-CIMD) for his work
on the configuration and build system of ICON making life of
developers much easier.

The work on the GPU code version, implemented based on OpenACC has
been supported by NVIDIA. Thanks a lot, Dmitry Alexeev (and his
colleagues).

For detailed information, please, have a look at

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

# Release notes for icon-2.6.4

Another new  release of icon is available.

It consists of

- fix for the soil water budget and snow-tiles
- integrate YAC2, corresponding changes in jsbach and the yaxt version are included
- bundled libraries are no checked at DWD and DKRZ

- fix current CDI-PIO setup for coupled setups
- use plain netcdf intead of CDI for restart IO

- fixes for the building system regarding new compilers, CLAW compatibility, changes in RTTOV and building against eccodes/grib_api

- bugfix for hfbasin diagnostic

- changes for the DWD NEC Aurora

- adjustments for ruby-0 setup (clean up of lsm masks, coupling, add irad_aero=12, buildbot test, scripts)

- New ICON-Land/JSBACH version: jsbach:master@39758f03
  This merge contains many changes, bug fixes, and improvements/new features in ICON-Land/JSBACH.
  The most important new features are:
    Time handling in JSBACH now supports nested domains (with ECHAM physics), plus fix for the time albedo is calculated for radiation
    New formulation for computation of roughness length (ported from JSBACH3) Implementation of anthropogenic land cover change
    Diagnostic 1d global mean JSBACH output for monitoring
    Implementation of the standalone JSBACH model as an ICON model component


- Implementation of the two-moment microphysics scheme by Seifert and Beheng (2006). The original NWP routines are used 'as is'. Use encapsulation.
- Unify variable long names.
- Add new testcases and/or bug fixes: RCEMIP_analytical, RCEtorus, Tconst, bubble.
- Add buildbot tests (S. Rast).
- Make iqneg output distinguishable. Fix typos.
- Add templates for 2 moment in NWP and SP physics.
- Add minmax diagnostics for microphysics inside SP.

- make dom_start_time and dom_end_time relative to experiment start/end time (Hauke Schulz)

- Simple fix of integer range in interface function for CDI-PIO

- update to icon-oes-1.3.05
- optional seaice initialization
- optional calculation of windstress in uncoupled ocean runs
- coupled ruby/dyamond: combine ocean and ice velocities in the coupling interface to the atmosphere

- extended_N-cycle (c8d10814 to 421c584e)

- update of mkexp scripts and mtime

- add perp. month/day

- For coupled configurations, the land and runoff are created using a
  fractional mask which is generated from a selected pair of ocean and
  atmosphere grids and mask. For clarity, these data are stored in a
  directories with names made up of both - atmosphere and ocean - grid
  IDs.

- With the current physics coupling only one set of pressure variables
  is needed. Pressure variables in mo_echam_phy_memory are therefore
  replaced. Use add_ref instead of add_var for several more fields.
  Reintroduce an optimised reading of Kinne aerosols, only those month
  are read in which are need for the current job rather than reading
  always the full year.

- introduce standard emissivity to all test cases
- Fix for Cariolle scheme
- Fix a problem with 10 m wind diagnostics in the coupled model
- In icon-dev.checksuite: Change LC-CTYPE to LC-ALL. With this correction icon-dev.checksuite works again on MacOS Big Sur.
- Bugfix for graupel initialization in Sapphire physics.
- Bug fix for ocean surface coupling in VDIFF:
     Missing A_klev+1/2 term is added, which induce ocean surface momentum stress to the opposite direction.
     Missing alpha term in wind stress diagnostic is added.
     Timing of wind stress diagnostic calculation is fixed.

- new test script exp.atm_ape_mlo_test for mixed layer ocean together with aquaplanet
- adjust experiment test scripts, remove several outdated and unused templates
- changes paths in templates from /pool/data/ICON/grids/private/rene/mpim to /pool/data/ICON/grids/public/mpim
- remove unused and completely outdated icon-authors.txt file
- update simple plume input file for amip reference experiment.

 - New reference data, tolerance ranges and AMIP reference experiment.

- New option for adaptive pressure bias correction at lateral boundaries (limited-area mode only)
- Revised upper boundary condition for vertically nested grids: reduces spurious reflection of vertically propagating sound and gravity waves at the uppermost nest interface level.
- Improved process description for gravity wave emission in SSO scheme (NWP)
- Update of Emvorado radar forward operator, including new options for dual polarization radars
- Update of effective radius coupling with radiation (NWP)
- OpenACC port for LHN code (NWP)
- optional diagnostic of lightning potential index and lightning flash density (NWP)
     new output fields: lfd_con, lfd_con_max, lpi_con, lpi_con_max, mlpi_con, mlpi_con_max, koi
- close output stream if open: fixes rare and random model crashes at simulation end on NEC

- removed hydrostatic DyCore
- removed interface to PSRAD radiation scheme (for NWP)
- removed COSMO ifdefs
- removed unused options for snow-cover fraction diagnosis (for NWP)

# Release notes for icon-2.6.3

The new release of icon is available.

It consists of

- a large number of bugfixes, refactorings, and optimizations for the models infrastructure
- bugfixes and improvements in the build environment
- consolidation of all QUBICC based enhancements and bugfixes (for the sapphire physics)
- many add-ons for the sapphire physics port to GPUs based on OpenACC and many further steps on porting the NWP physics to GPUs
- ART has made its first step to an git submodule external (draft implementation not to be used yet- no warranty)
- fixes and improvements in the ocean code including hamocc
- improvements in the data assimilation NWP physics coupling
- tuning of data assimilation 
- added rrtm-gp as radiation scheme for the sapphire physics on GPU
- much progress on cdi-pio use in many icon components
- refactoring of mpi communication library (with a focus on GPU to GPU communication)

For many more details visit:

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

June 7th, 2021

