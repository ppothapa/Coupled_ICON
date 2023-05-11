# ACX_SUBDIR_INIT_CONFIG(SUBDIR,
#                        [OPTIONS = recursive-help adjust-args run],
#                        [BUILD-SUBDIR = SUBDIR],
#                        [CONFIG-EXEC = configure])
# -----------------------------------------------------------------------------
# Initializes an Autoconf-based project residing in the directory SUBDIR (path
# relative to the top source directory of the top-level project) for the
# following configuration depending on the following OPTIONS (a space-separated
# list, defaults to 'recursive-help adjust-args run'):
#     [no-]recursive-help    whether the help message of the configure script
#                            in SUBDIR must be shown together with the help
#                            message of the top-level configure script when the
#                            latter is called with the argument
#                            '--help=recursive'. The positive variant of the
#                            option requires SUBDIR to be provided literally,
#                            without the shell variable indirection.
#     [no-]adjust-args       whether arguments provided for the top-level
#                            configure script must undergo a number of basic
#                            adjustments (including but not limited to
#                            '--prefix', '--srcdir' and '--cache-file' options)
#                            before they are passed to the configure script in
#                            SUBDIR.
#     [no-]run               whether the configure script in SUBDIR must be
#                            run by the top-level configure script.
#
# The configuration of the SUBDIR project is done by calling the CONFIG-EXEC
# (defaults to 'configure') script from the BUILD-SUBDIR (path relative to the
# top build directory of the top-level project, defaults to SUBDIR) directory.
#
# Sets variable extra_src_subdirs to the space-separated lists of all
# initialized SUBDIRs and BUILD-SUBDIRs (accounting for possible shell
# branching).
#
AC_DEFUN([ACX_SUBDIR_INIT_CONFIG],
  [m4_ifblank([$1], [m4_fatal([SUBDIR ('$1') cannot be blank])])dnl
   m4_pushdef([acx_subdir_opt_recursive_help], [recursive-help])dnl
   m4_pushdef([acx_subdir_opt_adjust_args], [adjust-args])dnl
   m4_pushdef([acx_subdir_opt_run], [run])dnl
   m4_foreach_w([opt], [$2],
     [m4_case(opt,
        [recursive-help], [m4_define([acx_subdir_opt_recursive_help], opt)],
        [no-recursive-help], [m4_define([acx_subdir_opt_recursive_help], opt)],
        [adjust-args], [m4_define([acx_subdir_opt_adjust_args], opt)],
        [no-adjust-args], [m4_define([acx_subdir_opt_adjust_args], opt)],
        [run], [m4_define([acx_subdir_opt_run], opt)],
        [no-run], [m4_define([acx_subdir_opt_run], opt)],
        [m4_fatal([unknown option ']opt['])])])dnl
   m4_cond([acx_subdir_opt_recursive_help], [recursive-help],
     [AS_LITERAL_IF([$1],
        [m4_append([_AC_LIST_SUBDIRS], [$1], [
])],
        [m4_fatal([option ']acx_subdir_opt_recursive_help[' requires ]dnl
[SUBDIR ('$1') to have a literal value])])])dnl
   ASX_SRCDIRS("m4_default([$3], [$1])")
   AS_VAR_SET([_ACX_SUBDIR_RUN_CMD_VAR([$1])],
     ["'m4_ifval([$3], [$ac_top_srcdir/$1], [$ac_srcdir])/]dnl
[m4_ifval([$4], ['$4], [configure'])"])
   AS_VAR_SET([_ACX_SUBDIR_RUN_DIR_VAR([$1])], ["m4_default([$3], [$1])"])
   AS_VAR_SET([_ACX_SUBDIR_BUILD_TYPE_VAR([$1])], ['config'])
   m4_cond([acx_subdir_opt_adjust_args], [adjust-args],
     [AC_REQUIRE_SHELL_FN([acx_subdir_pre_adjust_config_args], [],
        [AS_VAR_SET_IF([acx_subdir_pre_adjusted_config_args], [],
           [acx_subdir_pre_adjusted_config_args=$ac_configure_args
            _ACX_SUBDIR_REMOVE_ARGS([acx_subdir_pre_adjusted_config_args],
              [[ACX_SUBDIR_CONFIG_PATTERN_STDPOS([cache-file])| \
                ACX_SUBDIR_CONFIG_PATTERN_STDPOS([srcdir])| \
                ACX_SUBDIR_CONFIG_PATTERN_STDPOS([prefix])], [1]],
              [[ACX_SUBDIR_CONFIG_PATTERN_STDOPT([cache-file])| \
                ACX_SUBDIR_CONFIG_PATTERN_STDOPT([srcdir])| \
                ACX_SUBDIR_CONFIG_PATTERN_STDOPT([prefix])| \
                --config-cache|-C| \
                ACX_SUBDIR_CONFIG_PATTERN_ENABLE([option-checking])], [0]])
            AS_VAR_APPEND([acx_subdir_pre_adjusted_config_args],
              [" '--disable-option-checking'"])
            AS_VAR_IF([prefix], [NONE],
              [acx_tmp="--prefix=$ac_default_prefix"],
              [acx_tmp="--prefix=$prefix"])
            ASX_ESCAPE_SINGLE_QUOTE([acx_tmp])
            AS_VAR_APPEND([acx_subdir_pre_adjusted_config_args],
              [" '$acx_tmp'"])])])dnl
      acx_subdir_pre_adjust_config_args
      AS_VAR_SET([_ACX_SUBDIR_RUN_ARG_VAR([$1])],
        [$acx_subdir_pre_adjusted_config_args])
      AS_VAR_IF([cache_file], ['/dev/null'],
        [acx_tmp=$cache_file],
        [AS_CASE([$cache_file],
           [[[\\/]]* | ?:[[\\/]]*],
           [acx_tmp=$cache_file],
           [acx_tmp="$ac_top_build_prefix$cache_file"])
         acx_tmp="$acx_tmp.AS_LITERAL_IF([$1],
                             [m4_translit([$1], [/], [.])],
                             [`echo "$1" | tr / .`])"])
      _ACX_SUBDIR_APPEND_ARGS([_ACX_SUBDIR_RUN_ARG_VAR([$1])],
        ["--cache-file=$acx_tmp"],
        ["--srcdir=m4_ifval([$3], [$ac_top_srcdir/$1], [$ac_srcdir])"])],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_ARG_VAR([$1])], [$ac_configure_args])])
   m4_divert_once([DEFAULTS], [extra_src_subdirs=])dnl
   AS_VAR_APPEND([extra_src_subdirs], [" $1"])
   m4_cond([acx_subdir_opt_run], [run],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [yes])[]dnl
      _ACX_SUBDIR_COMMANDS_PRE],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [no])])[]dnl
   m4_popdef([acx_subdir_opt_recursive_help])dnl
   m4_popdef([acx_subdir_opt_adjust_args])dnl
   m4_popdef([acx_subdir_opt_run])])

# ACX_SUBDIR_INIT_CMAKE(SUBDIR,
#                       [OPTIONS = adjust-args run],
#                       [BUILD-SUBDIR = SUBDIR/build],
#                       [CMAKE-EXEC = cmake])
# -----------------------------------------------------------------------------
# Initializes a CMake-based project residing in the directory SUBDIR (path
# relative to the top source directory of the top-level project) for the
# following configuration depending on the following OPTIONS (a space-separated
# list, defaults to 'adjust-args run'):
#     [no-]adjust-args    whether arguments provided for the top-level
#                         configure script must be translated into CMake
#                         arguments (including but not limited to
#                         'CC -> -DCMAKE_C_COMPILER',
#                         '--prefix -> -DCMAKE_INSTALL_PREFIX', etc.).
#     [no-]run            whether CMAKE-EXEC must be run for SUBDIR by the
#                         top-level configure script.
#
# The configuration of the SUBDIR project is done by calling CMAKE-EXEC
# (defaults to '${CMAKE-cmake}') from the BUILD-SUBDIR (path relative to the
# top build directory of the top-level project, defaults to SUBDIR/build)
# directory.
#
# Sets variable extra_src_subdirs to the space-separated lists of all
# initialized SUBDIRs and BUILD-SUBDIRs (accounting for possible shell
# branching).
#
AC_DEFUN([ACX_SUBDIR_INIT_CMAKE],
  [m4_ifblank([$1], [m4_fatal([SUBDIR ('$1') cannot be blank])])dnl
   m4_pushdef([acx_subdir_opt_adjust_args], [adjust-args])dnl
   m4_pushdef([acx_subdir_opt_run], [run])dnl
   m4_foreach_w([opt], [$2],
     [m4_case(opt,
        [adjust-args], [m4_define([acx_subdir_opt_adjust_args], opt)],
        [no-adjust-args], [m4_define([acx_subdir_opt_adjust_args], opt)],
        [run], [m4_define([acx_subdir_opt_run], opt)],
        [no-run], [m4_define([acx_subdir_opt_run], opt)],
        [m4_fatal([unknown option ']opt['])])])dnl
   AS_VAR_SET([_ACX_SUBDIR_RUN_CMD_VAR([$1])], ["m4_default([$4], ['cmake'])"])
   AS_VAR_SET([_ACX_SUBDIR_RUN_DIR_VAR([$1])], ["m4_default([$3], [$1/build])"])
   AS_VAR_SET([_ACX_SUBDIR_BUILD_TYPE_VAR([$1])], ['cmake'])
   m4_cond([acx_subdir_opt_adjust_args], [adjust-args],
     [AC_REQUIRE_SHELL_FN([acx_subdir_pre_adjust_cmake_args], [],
        [AS_VAR_SET_IF([acx_subdir_pre_adjusted_cmake_args], [],
           [AS_VAR_SET([acx_subdir_pre_adjusted_cmake_args],
              ["'-Wno-dev' '--no-warn-unused-cli' '-GUnix Makefiles'"])
            eval "set dummy $ac_configure_args"; shift
dnl Transform standard precious (influential environment) variables:
            m4_pushdef([acx_subdir_known_args],
              [[AR],
               [RANLIB],
               [CC, [CMAKE_C_COMPILER]],
               [CFLAGS, [CMAKE_C_FLAGS]],
               [CXX, [CMAKE_CXX_COMPILER]],
               [CXXFLAGS, [CMAKE_CXX_FLAGS]],
               [CPPFLAGS],
               [FC, [CMAKE_Fortran_COMPILER]],
               [FCFLAGS, [CMAKE_Fortran_FLAGS]],
               [LDFLAGS, [CMAKE_EXE_LINKER_FLAGS
                          CMAKE_MODULE_LINKER_FLAGS
                          CMAKE_SHARED_LINKER_FLAGS]],
               [LIBS, [CMAKE_C_STANDARD_LIBRARIES
                       CMAKE_CXX_STANDARD_LIBRARIES
                       CMAKE_Fortran_STANDARD_LIBRARIES]]])dnl
            acx_subdir_cmake_vars_to_transform=
            for acx_tmp; do
              AS_CASE([$acx_tmp],
                [m4_join([|],
                   m4_foreach([pair],
                     [acx_subdir_known_args], [m4_car(pair)=*,]))],
                [acx_arg_name=`expr "x$acx_tmp" : 'x\(@<:@^=@:>@*\)='`
                 AS_VAR_APPEND([acx_subdir_cmake_vars_to_transform],
                   [" $acx_arg_name"])
                 acx_arg_cmd_value=`expr "x$acx_tmp" : '@<:@^=@:>@*=\(.*\)'`
                 AS_VAR_COPY([acx_arg_${acx_arg_name}], [acx_arg_cmd_value])])
            done
dnl CMake requires the archiver and the archive indexer commands to be set as
dnl absolute paths. Otherwise, it will try to find the executable in the build
dnl directory. Also, AR and RANLIB are not supposed to be paths to executables
dnl with arguments because it will cause CMake to choke:
            for acx_tmp in AR RANLIB; do
              AS_CASE([" $acx_subdir_cmake_vars_to_transform "],
                [*" $acx_tmp "*],
                [acx_prog_search_abspath=unknown
                 AS_VAR_COPY([acx_prog_exec], [acx_arg_${acx_tmp}])
                 AS_CASE([$acx_prog_exec],
                   [*[[\\/]]*],
                   [AS_IF([AS_EXECUTABLE_P([$acx_prog_exec])],
                      [acx_prog_search_abspath=$acx_prog_exec])],
                   [_AS_PATH_WALK([],
                      [AS_IF([AS_EXECUTABLE_P(["$as_dir/$acx_prog_exec"])],
                         [acx_prog_search_abspath="$as_dir/$acx_prog_exec"
                          break])])])
                 AS_VAR_IF([acx_prog_search_abspath], [unknown],
                   [AC_MSG_WARN([unable to convert argument $acx_tmp to its dnl
CMake equivalent(s): absolute path to "$acx_prog_exec" is not found])],
                   [AS_VAR_COPY([acx_arg_${acx_tmp}_ABSPATH],
                      [acx_prog_search_abspath])
                    AS_VAR_APPEND([acx_subdir_cmake_vars_to_transform],
                      [" ${acx_tmp}_ABSPATH"])])])
            done
            m4_append([acx_subdir_known_args],
              [[AR_ABSPATH, [CMAKE_AR
                             CMAKE_C_COMPILER_AR
                             CMAKE_CXX_COMPILER_AR
                             CMAKE_Fortran_COMPILER_AR]],
               [RANLIB_ABSPATH, [CMAKE_RANLIB
                                 CMAKE_C_COMPILER_RANLIB
                                 CMAKE_CXX_COMPILER_RANLIB
                                 CMAKE_Fortran_COMPILER_RANLIB]]], [,])dnl
dnl CMake has no explicit support for CPPFLAGS, therefore, we append them to
dnl CFLAGS and CXXFLAGS:
            AS_CASE([" $acx_subdir_cmake_vars_to_transform "],
              [*' CPPFLAGS '*],
              [for acx_tmp in CFLAGS CXXFLAGS; do
                 AS_CASE([" $acx_subdir_cmake_vars_to_transform "],
                   [*" $acx_tmp "*],
                   [AS_VAR_APPEND([acx_arg_${acx_tmp}],
                      [" $acx_arg_CPPFLAGS"])],
                   [AS_VAR_COPY([acx_arg_${acx_tmp}], [acx_arg_CPPFLAGS])
                    AS_VAR_APPEND([acx_subdir_cmake_vars_to_transform],
                      [" $acx_tmp"])])
               done])
dnl Also process the installation prefix:
            m4_append([acx_subdir_known_args],
              [[PREFIX, [CMAKE_INSTALL_PREFIX]]], [,])dnl
            AS_VAR_IF([prefix], [NONE],
              [acx_arg_PREFIX=$ac_default_prefix],
              [acx_arg_PREFIX=$prefix])
            AS_VAR_APPEND([acx_subdir_cmake_vars_to_transform], [' PREFIX'])
dnl Append the transformed arguments:
            for acx_tmp in $acx_subdir_cmake_vars_to_transform; do
              acx_subdir_cmake_vars_to_set=
              AS_CASE([$acx_tmp],
                m4_foreach([pair], [acx_subdir_known_args],
                  [m4_quote(m4_car(pair)),
                   m4_quote(acx_subdir_cmake_vars_to_set=dnl
'm4_normalize(m4_cdr(pair))'),]))
              AS_IF([test -n "$acx_subdir_cmake_vars_to_set"],
                [AS_VAR_COPY([acx_subdir_quoted_value], [acx_arg_${acx_tmp}])
                 ASX_ESCAPE_SINGLE_QUOTE([acx_subdir_quoted_value])
                 for acx_subdir_cmake_var in $acx_subdir_cmake_vars_to_set; do
                   AS_VAR_APPEND([acx_subdir_pre_adjusted_cmake_args],
                     [" '-D$acx_subdir_cmake_var=$acx_subdir_quoted_value'"])
                 done])
              AS_UNSET([acx_arg_${acx_tmp}])
            done
            m4_popdef([acx_subdir_known_args])])])dnl
      acx_subdir_pre_adjust_cmake_args
      AS_VAR_SET([_ACX_SUBDIR_RUN_ARG_VAR([$1])],
        [$acx_subdir_pre_adjusted_cmake_args])
      m4_ifval([$3],
        [ASX_SRCDIRS(["$3"])
         acx_tmp="-S$ac_top_srcdir/$1"],
        [ASX_SRCDIRS(["$1/build"])
         AS_CASE([$ac_srcdir],
           [.], [acx_tmp='-S..'],
           [acx_tmp="-S$ac_top_srcdir/$1"])])
      ASX_ESCAPE_SINGLE_QUOTE([acx_tmp])
      AS_VAR_APPEND([_ACX_SUBDIR_RUN_ARG_VAR([$1])], [" '$acx_tmp'"])],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_ARG_VAR([$1])])])
   m4_divert_once([DEFAULTS], [extra_src_subdirs=])dnl
   AS_VAR_APPEND([extra_src_subdirs], [" $1"])
   m4_cond([acx_subdir_opt_run], [run],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [yes])[]dnl
      _ACX_SUBDIR_COMMANDS_PRE],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [no])])[]dnl
   m4_popdef([acx_subdir_opt_adjust_args])dnl
   m4_popdef([acx_subdir_opt_run])])

# ACX_SUBDIR_INIT_IFELSE(SUBDIR,
#                        [ACTION-IF-INITIALIZED],
#                        [ACTION-IF-NOT-INITIALIZED])
# -----------------------------------------------------------------------------
# Checks whether SUBDIR was actually (accounting for possible shell branching)
# initialized and runs ACTION-IF-INITIALIZED if that is the case. Otherwise,
# runs ACTION-IF-NOT-INITIALIZED.
#
# For example, the top-level configure script needs to run 'some-command'
# inside AC_CONFIG_COMMANDS_PRE but only if 'subdir/configure' was actually
# run. A way to implement that is to expand the following:
#
# AS_IF([some_complex_condition],
#   [ACX_SUBDIR_INIT_CONFIG([subdir])], [run])
# AC_CONFIG_COMMANDS_PRE(
#   [ACX_SUBDIR_INIT_IFELSE([subdir], [some-command])])
#
AC_DEFUN([ACX_SUBDIR_INIT_IFELSE],
  [AS_CASE([" $extra_src_subdirs "], [*' $1 '*], [$2], [$3])])

# ACX_SUBDIR_REMOVE_ARGS(SUBDIR,
#                        [PATTERN...])
# -----------------------------------------------------------------------------
# Expands to the shell script that removes all argument of the configuration
# command for SUBDIR that match PATTERNs. A PATTERN can be either a shell case
# pattern or a comma-separated pair of a shell case pattern and a non-negative
# integer representing the number of arguments that must be dropped after the
# pattern match.
#
# Consider using ACX_SUBDIR_CONFIG_PATTERN_ENABLE,
# ACX_SUBDIR_CONFIG_PATTERN_WITH, ACX_SUBDIR_CONFIG_PATTERN_STDOPT and
# ACX_SUBDIR_CONFIG_PATTERN_STDPOS below to generate patterns that match
# standard Autoconf options.
#
AC_DEFUN([ACX_SUBDIR_REMOVE_ARGS],
  [_ACX_SUBDIR_REMOVE_ARGS(_ACX_SUBDIR_RUN_ARG_VAR([$1]),
     m4_unquote(m4_cdr($@)))])

# ACX_SUBDIR_CONFIG_PATTERN_ENABLE(FEATURE)
# -----------------------------------------------------------------------------
# Expands to a shell case pattern that matches all valid arguments introduced
# with the standard Autoconf macro AC_ARG_ENABLE([PACKAGE]).
#
AC_DEFUN([ACX_SUBDIR_CONFIG_PATTERN_ENABLE],
  [[-enable-$1|-enable-$1=*|--enable-$1|--enable-$1=*|-disable-$1|]dnl
[--disable-$1]])

# ACX_SUBDIR_CONFIG_PATTERN_WITH(PACKAGE)
# -----------------------------------------------------------------------------
# Expands to a shell case pattern that matches all valid arguments introduced
# with the standard Autoconf macro AC_ARG_WITH([PACKAGE]).
#
AC_DEFUN([ACX_SUBDIR_CONFIG_PATTERN_WITH],
  [[-with-$1|-with-$1=*|--with-$1|--with-$1=*|-without-$1|--without-$1]])

# ACX_SUBDIR_CONFIG_PATTERN_STDOPT(ARG)
# -----------------------------------------------------------------------------
# Expands to a shell case pattern that matches all possible ways to set a
# standard Autoconf argument ARG (set by _AC_INIT_PARSE_ARGS as part of
# AC_INIT) as a single command-line option (i.e. without an extra argument for
# the value).
#
# For example, ACX_SUBDIR_CONFIG_PATTERN_STDOPT([prefix]) expands to
#     -prefix=*|--prefix=*|--prefi=*|--pref=*|--pre=*|--pr=*|--p=*
#
AC_DEFUN([ACX_SUBDIR_CONFIG_PATTERN_STDOPT],
  [-$1=*|--$1=*[]dnl
m4_for([index], m4_decr(m4_len([$1])), m4_if([$1], [srcdir], [2], [1]), [-1],
     [|--m4_substr([$1], [0], index)=*])])

# ACX_SUBDIR_CONFIG_PATTERN_STDPOS(ARG)
# -----------------------------------------------------------------------------
# Expands to a shell case pattern that matches all possible ways to set a
# standard Autoconf argument ARG (set by _AC_INIT_PARSE_ARGS as part of
# AC_INIT) as a pair of command-line options (i.e. the option name and the
# option value).
#
# For example, ACX_SUBDIR_CONFIG_PATTERN_STDPOS([prefix]) expands to
#     -prefix|--prefix|--prefi|--pref|--pre|--pr|--p
#
AC_DEFUN([ACX_SUBDIR_CONFIG_PATTERN_STDPOS],
  [-$1|--$1[]dnl
m4_for([index], m4_decr(m4_len([$1])), m4_if([$1], [srcdir], [2], [1]), [-1],
     [|--m4_substr([$1], [0], index)])])

# ACX_SUBDIR_APPEND_ARGS(SUBDIR,
#                        [ARG...])
# -----------------------------------------------------------------------------
# Expands to a shell script that appends arguments ARGs for the command that
# configures the SUBDIR directory.
#
AC_DEFUN([ACX_SUBDIR_APPEND_ARGS],
  [_ACX_SUBDIR_APPEND_ARGS(_ACX_SUBDIR_RUN_ARG_VAR([$1]),
     m4_unquote(m4_cdr($@)))])

# ACX_SUBDIR_RUN_RESET(SUBDIR,
#                      [VALUE = SWAP])
# -----------------------------------------------------------------------------
# Expands to a shell script that adds/removes SUBDIR on/from the list of
# directories that must be configured by the top-level configure script.
#
# The macro provides a way to override the [no-]run option of the
# ACX_SUBDIR_INIT_CONFIG macro (see above) at the run-time of the top-level
# configure script. If the VALUE argument of the macro is provided, it must
# be a literal or a shell expressions that evaluates either to 'yes' (add
# SUBDIR on the list) or to 'no' (remove SUBDIR from the list). An omitted
# VALUE implies that SUBDIR is added on the list if it is not there and removed
# from it otherwise.
#
AC_DEFUN([ACX_SUBDIR_RUN_RESET],
  [m4_ifval([$2],
     [AS_VAR_SET([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [$2])],
     [acx_tmp=yes
      AS_VAR_IF([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [yes], [acx_tmp=no])
      AS_VAR_COPY([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [acx_tmp])])[]dnl
   _ACX_SUBDIR_COMMANDS_PRE])

# ACX_SUBDIR_RUN_EXTRA(SUBDIR,
#                      [BEFORE],
#                      [AFTER])
# -----------------------------------------------------------------------------
# Expands to a shell script that registers additional shell commands that must
# be run BEFORE or AFTER the configuration of SUBDIR done by the top-level
# configure (that is, the resulting script has no effect if SUBDIR is removed
# from the list of directories that must be configured by the top-level
# configure script).
#
# The additional commands are run in the top build directory of the top-level
# project.
#
AC_DEFUN([ACX_SUBDIR_RUN_EXTRA],
  [m4_ifval([$2],
     [_ACX_SUBDIR_RUN_EXTRA([_ACX_SUBDIR_RUN_BEFORE_VAR([$1])], [$2])])dnl
   m4_ifval([$3],
     [m4_ifval([$2],[
])_ACX_SUBDIR_RUN_EXTRA([_ACX_SUBDIR_RUN_AFTER_VAR([$1])], [$3])])])

# ACX_SUBDIR_RUN_IFELSE(SUBDIR,
#                       [ACTION-IF-RUN],
#                       [ACTION-IF-NO-RUN])
# -----------------------------------------------------------------------------
# Checks whether SUBDIR is on the list of directories that must be configured
# by the top-level configure script and runs ACTION-IF-RUN if that is the case.
# Otherwise, runs ACTION-IF-NO-RUN.
#
AC_DEFUN([ACX_SUBDIR_RUN_IFELSE],
  [AS_VAR_IF([_ACX_SUBDIR_RUN_YESNO_VAR([$1])], [yes], [$2], [$3])])

# ACX_SUBDIR_GET_BUILD_DIR(VARIABLE,
#                          SUBDIR)
# -----------------------------------------------------------------------------
# Expands to a shell script that sets the shell variable VARIABLE to the name
# of build directory that corresponds to the SUBDIR source directory (i.e. to
# the value of the BUILD-SUBDIR option provided to the respective expansion of
# the ACX_SUBDIR_INIT_CONFIG macro).
#
AC_DEFUN([ACX_SUBDIR_GET_BUILD_DIR],
  [AS_VAR_COPY([$1], [_ACX_SUBDIR_RUN_DIR_VAR([$2])])])

# ACX_SUBDIR_GET_BUILD_TYPE(VARIABLE,
#                           SUBDIR)
# -----------------------------------------------------------------------------
# Expands to a shell script that sets the shell variable VARIABLE to the type
# of the build system of the SUBDIR source directory.
#
# Possible out values are:
#   "config" - Autoconf-based build system;
#   "cmake"  - CMake-based build system.
#
AC_DEFUN([ACX_SUBDIR_GET_BUILD_TYPE],
  [AS_VAR_COPY([$1], [_ACX_SUBDIR_BUILD_TYPE_VAR([$2])])])

# ACX_SUBDIR_GET_RUN_CMD(VARIABLE,
#                        SUBDIR)
# -----------------------------------------------------------------------------
# Expands to a shell script that sets the shell variable VARIABLE to the full
# command (together with the arguments) that configures SUBDIR when run from
# the corresponding build directory (i.e. from BUILD-SUBDIR provided to the
# respective expansion of the ACX_SUBDIR_INIT_CONFIG macro).
#
# The command neither creates the build directory nor switches to it.
# Therefore, if SUBDIR is not on the list of directories that must be
# configured by the top-level configure script, it is the user's responsibility
# to extend the command accordingly.
#
AC_DEFUN([ACX_SUBDIR_GET_RUN_CMD],
  [AS_VAR_SET([$1],
     ["AS_VAR_GET(_ACX_SUBDIR_RUN_CMD_VAR([$2])) dnl
AS_VAR_GET(_ACX_SUBDIR_RUN_ARG_VAR([$2]))"])])

# ACX_SUBDIR_QUERY_CONFIG_STATUS(VARIABLE,
#                                SUBDIR,
#                                TEMPLATE)
# -----------------------------------------------------------------------------
# Expands to a shell script that sets the shell variable VARIABLE to the result
# of the variable substitution TEMPLATE done by the config.status script
# residing inside the build directory that corresponds to the SUBDIR source
# directory. The macro provides means of getting the results of the configure
# script from the subdirectory to the top-level configure script.
#
# For example, 'subdir/configure' sets an output variable 'LIBM' and its value
# needs to be known in the top-level configure script. A way to implement that
# is to expand the following:
#
# ACX_SUBDIR_INIT_CONFIG([subdir])
# AC_CONFIG_COMMANDS_PRE(
#   [ACX_SUBDIR_QUERY_CONFIG_STATUS([SUBDIR_VAR], [subdir], [@LIBM@])
#    AC_SUBST([SUBDIR_VAR])])
#
AC_DEFUN([ACX_SUBDIR_QUERY_CONFIG_STATUS],
  [acx_tmp=`AS_ECHO([$3]) | dnl
"AS_VAR_GET(_ACX_SUBDIR_RUN_DIR_VAR([$2]))/config.status" dnl
-q --file=- 2>/dev/null`
   AS_IF([test $? -eq 0],
     [AS_VAR_COPY([$1], [acx_tmp])],
     [AC_MSG_ERROR(
        [unable to run ]dnl
['AS_VAR_GET(_ACX_SUBDIR_RUN_DIR_VAR([$2]))/config.status'])])])

# _ACX_SUBDIR_COMMANDS_PRE()
# -----------------------------------------------------------------------------
# Registers a sequence of commands that are run by the top-level configure
# script (right before creating its own config.status) that configure all
# subdirectories that are added on the respective list (see the [no-]run option
# of the ACX_SUBDIR_INIT_CONFIG macro above).
#
m4_define([_ACX_SUBDIR_COMMANDS_PRE],
  [m4_ifndef([_ACX_SUBDIR_COMMANDS_PRE_DEFINED],
     [AC_CONFIG_COMMANDS_PRE(
        [AS_IF([test "x$no_recursion" != xyes],
           [AS_VAR_IF([silent], [yes],
              [acx_subdir_silent_arg="'--silent'"],
              [acx_subdir_silent_arg=])
            acx_subdir_run_any=no
            for acx_subdir_srcdir in $extra_src_subdirs; do
              AS_VAR_IF(
                [_ACX_SUBDIR_RUN_YESNO_VAR([$acx_subdir_srcdir])], [yes],
                [acx_subdir_run_any=yes
                 AS_VAR_COPY([acx_subdir_builddir],
                   [_ACX_SUBDIR_RUN_DIR_VAR([$acx_subdir_srcdir])])
                 AS_VAR_COPY([acx_subdir_fns],
                   [_ACX_SUBDIR_RUN_BEFORE_VAR([$acx_subdir_srcdir])])
                 AS_IF([test -n "$acx_subdir_fns"],
                   [acx_tmp=dnl
"=== running extra commands before configuring $acx_subdir_srcdir"
                    _AS_ECHO_LOG([$acx_tmp])
                    _AS_ECHO([$acx_tmp])
                    for acx_tmp in $acx_subdir_fns; do
                      eval "$acx_tmp"
                    done])
                 acx_tmp=dnl
"=== configuring $acx_subdir_srcdir (in $acx_subdir_builddir)"
                 _AS_ECHO_LOG([$acx_tmp])
                 _AS_ECHO([$acx_tmp])
                 AS_MKDIR_P(["$acx_subdir_builddir"])
                 ACX_SUBDIR_GET_RUN_CMD([acx_subdir_run_cmd],
                   [$acx_subdir_srcdir])
                 acx_subdir_run_cmd=dnl
"( cd '$acx_subdir_builddir' && $acx_subdir_run_cmd $acx_subdir_silent_arg)"
                 AC_MSG_NOTICE([running $acx_subdir_run_cmd])
                 eval "$acx_subdir_run_cmd" || dnl
AC_MSG_ERROR([configuration of $acx_subdir_srcdir failed])
                 AS_VAR_COPY([acx_subdir_fns],
                   [_ACX_SUBDIR_RUN_AFTER_VAR([$acx_subdir_srcdir])])
                 AS_IF([test -n "$acx_subdir_fns"],
                   [acx_tmp=dnl
"=== running extra commands after configuring $acx_subdir_srcdir"
                    _AS_ECHO_LOG([$acx_tmp])
                    _AS_ECHO([$acx_tmp])
                    for acx_tmp in $acx_subdir_fns; do
                      eval "$acx_tmp"
                    done])])
               done
               AS_VAR_IF([acx_subdir_run_any], [yes],
                 [_AS_ECHO([===])
                  _AS_ECHO_LOG([===])])])])dnl
      m4_define([_ACX_SUBDIR_COMMANDS_PRE_DEFINED])])])

# _ACX_SUBDIR_REMOVE_ARGS(VARIABLE,
#                         [PATTERN...])
# -----------------------------------------------------------------------------
# Expands to the shell script that removes all argument of the configuration
# command stored in the shell variable VARIABLE as a space-separated list of
# single-quoted elements that match PATTERNs. A PATTERN can be either a shell
# case pattern or a comma-separated pair of a shell case pattern and a
# non-negative integer representing the number of arguments that must be
# dropped after the pattern match.
#
m4_define([_ACX_SUBDIR_REMOVE_ARGS],
  [eval "set dummy AS_VAR_GET([$1])"; shift
   AS_VAR_SET([$1])
   while test $[]# != 0; do
     acx_tmp=$[]1
     AS_CASE([$acx_tmp],
       m4_foreach([pattern], m4_cdr($@),
         [m4_ifnblank(m4_car(pattern),
            [m4_quote(m4_car(pattern)),
             m4_quote(
               m4_cond([m4_default(m4_argn(2, pattern), 0)],
                 [0], [],
                 [m4_for([], [1], m4_argn(2, pattern), [1],
                    [test 2 -gt $[]@%:@ || shift; ])])),])])dnl
       [ASX_ESCAPE_SINGLE_QUOTE([acx_tmp])
        AS_VAR_APPEND([$1], [" '$acx_tmp'"])])
     shift
   done])

# _ACX_SUBDIR_APPEND_ARGS(VARIABLE,
#                         [ARG...])
# -----------------------------------------------------------------------------
# Expands to a shell script that appends arguments ARGs to the value of the
# shell variable VARIABLE as single-quoted elements.
#
m4_define([_ACX_SUBDIR_APPEND_ARGS],
  [set dummy[]m4_foreach([arg], m4_cdr($@), [ arg]); shift
   for acx_tmp; do
     ASX_ESCAPE_SINGLE_QUOTE([acx_tmp])
     AS_VAR_APPEND([$1], [" '$acx_tmp'"])
   done])

# _ACX_SUBDIR_RUN_EXTRA(VARIABLE,
#                       [CMD])
# -----------------------------------------------------------------------------
# Registers a shell function that runs the CMD shell code and appends its name
# to the shell variable VARIABLE.
#
m4_define([_ACX_SUBDIR_RUN_EXTRA],
  [m4_pushdef([acx_subdir_extra_run_idx])dnl
   m4_ifdef([_ACX_SUBDIR_RUN_EXTRA_COUNT],
     [m4_define([acx_subdir_extra_run_idx], _ACX_SUBDIR_RUN_EXTRA_COUNT)],
     [m4_define([acx_subdir_extra_run_idx], 0)])dnl
   m4_define([_ACX_SUBDIR_RUN_EXTRA_COUNT],
     m4_incr(acx_subdir_extra_run_idx))dnl
   m4_pushdef([acx_subdir_extra_run_name],
     [acx_subdir_extra_run_[]acx_subdir_extra_run_idx])dnl
   AC_REQUIRE_SHELL_FN(acx_subdir_extra_run_name, [], [$2])dnl
   AS_VAR_APPEND([$1], [' acx_subdir_extra_run_name'])dnl
   m4_popdef([acx_subdir_extra_run_name])dnl
   m4_popdef([acx_subdir_extra_run_idx])])

# _ACX_SUBDIR_RUN_YESNO_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds the condition ('yes' or
# 'no') of whether the configure script runs the configuration command for
# directory SUBDIR.
#
m4_define([_ACX_SUBDIR_RUN_YESNO_VAR],
  [acx_subdir_run_condition_[]AS_TR_SH([$1])])

# _ACX_SUBDIR_RUN_BEFORE_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds names of the functions that
# must be run before configuring in directory SUBDIR.
#
m4_define([_ACX_SUBDIR_RUN_BEFORE_VAR],
  [acx_subdir_run_before_[]AS_TR_SH([$1])])

# _ACX_SUBDIR_RUN_AFTER_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds names of the functions that
# must be run after configuring in directory SUBDIR.
#
m4_define([_ACX_SUBDIR_RUN_AFTER_VAR],
  [acx_subdir_run_after_[]AS_TR_SH([$1])])

# _ACX_SUBDIR_RUN_DIR_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds the name of the run (build)
# directory of the configuration command of the source directory SUBDIR.
#
m4_define([_ACX_SUBDIR_RUN_DIR_VAR],
  [acx_subdir_run_dir_[]AS_TR_SH([$1])])

# _ACX_SUBDIR_BUILD_TYPE_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds the type of the build system
# of the source directory SUBDIR.
#
m4_define([_ACX_SUBDIR_BUILD_TYPE_VAR],
  [acx_subdir_build_type_[]AS_TR_SH([$1])])

# _ACX_SUBDIR_RUN_CMD_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds the command (without the
# arguments) that configures directory SUBDIR.
#
m4_define([_ACX_SUBDIR_RUN_CMD_VAR],
  [acx_subdir_run_cmd_[]AS_TR_SH([$1])])

# _ACX_SUBDIR_RUN_ARG_VAR(SUBDIR)
# -----------------------------------------------------------------------------
# Expands to the name of shell variable that holds arguments of the command
# that configures directory SUBDIR.
#
m4_define([_ACX_SUBDIR_RUN_ARG_VAR],
  [acx_subdir_run_args_[]AS_TR_SH([$1])])
