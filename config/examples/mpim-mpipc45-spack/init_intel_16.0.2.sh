BUILD_ENV='export ICCCFG="/sw/jessie-x64/intel/intel-16.0.2/bin/rpath.cfg"; export IFORTCFG="$ICCCFG";'

PLAIN_CC='/sw/jessie-x64/intel/intel-16.0.2/bin/icc'

SW_ROOT='/scratch/local1/icon-bootstrap/opt/intel-16.0.2'

MPICH_ROOT="${SW_ROOT}/mpich-3.2.1-gubxzw3"
MPICH_FC="${MPICH_ROOT}/bin/mpif90"
MPICH_CC="${MPICH_ROOT}/bin/mpicc"

NETCDF_ROOT="${SW_ROOT}/netcdf-4.6.1-gwvk55n"
NETCDF_CPPFLAGS="-I${NETCDF_ROOT}/include"
NETCDF_LDFLAGS="-L${NETCDF_ROOT}/lib"
NETCDF_LIBS='-lnetcdf'

NETCDFF_ROOT="${SW_ROOT}/netcdf-fortran-4.4.4-js46vtx"
NETCDFF_FCFLAGS="-I${NETCDFF_ROOT}/include"
NETCDFF_LDFLAGS="-L${NETCDFF_ROOT}/lib"
NETCDFF_LIBS='-lnetcdff'

SCT_ROOT="${SW_ROOT}/libsct-develop-qexve7s"
SCT_FCFLAGS="-I${SCT_ROOT}/include"
SCT_LDFLAGS="-L${SCT_ROOT}/lib"
SCT_LIBS='-lsct'

YAXT_ROOT="${SW_ROOT}/yaxt-0.7.0-lazd2ez"
YAXT_FCFLAGS="-I${YAXT_ROOT}/include"
YAXT_LDFLAGS="-L${YAXT_ROOT}/lib"
YAXT_LIBS='-lyaxt'

GRIBAPI_ROOT="${SW_ROOT}/grib-api-1.24.0-x7mc6ta"
GRIBAPI_CPPFLAGS="-I${GRIBAPI_ROOT}/include"
GRIBAPI_LDFLAGS="-L${GRIBAPI_ROOT}/lib"
GRIBAPI_LIBS='-lgrib_api'

CDI_ROOT="${SW_ROOT}/libcdi-1.8.2-rrdsdao"
CDI_FCFLAGS="-I${CDI_ROOT}/include"
CDI_LDFLAGS="-L${CDI_ROOT}/lib"
CDI_LIBS='-lcdi_f2003 -lcdi'

MTIME_ROOT="${SW_ROOT}/libmtime-1.0.8-pul5lv4"
MTIME_FCFLAGS="-I${MTIME_ROOT}/include"
MTIME_CPPFLAGS="-I${MTIME_ROOT}/include"
MTIME_LDFLAGS="-L${MTIME_ROOT}/lib"
MTIME_LIBS='-lmtime'

BLAS_ROOT="${SW_ROOT}/netlib-lapack-3.8.0-yuqhwfu"
BLAS_LDFLAGS="-L${BLAS_ROOT}/lib"
BLAS_LIBS='-lblas'

LAPACK_ROOT="${SW_ROOT}/netlib-lapack-3.8.0-yuqhwfu"
LAPACK_LDFLAGS="-L${LAPACK_ROOT}/lib"
LAPACK_LIBS='-llapack'

XML2_ROOT="${SW_ROOT}/libxml2-2.9.8-hg44wox"
XML2_CPPFLAGS="-I${XML2_ROOT}/include/libxml2"
XML2_LDFLAGS="-L${XML2_ROOT}/lib"
XML2_LIBS='-lxml2'

YAC_ROOT="${SW_ROOT}/yac-1.5.2-rc-h4gpa5y"
YAC_FCFLAGS="-I${YAC_ROOT}/include"
YAC_LDFLAGS="-L${YAC_ROOT}/lib"
YAC_LIBS='-lyac'

SELF_ROOT="${SW_ROOT}/libself-0.2-tixcibg"
SELF_FCFLAGS="-I${SELF_ROOT}/include"
SELF_LDFLAGS="-L${SELF_ROOT}/lib"
SELF_LIBS='-lself'

