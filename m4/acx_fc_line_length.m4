# ACX_FC_LINE_LENGTH([LENGTH = 132],
#                    [ACTION-IF-SUCCESS = APPEND-FCFLAGS],
#                    [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to accept lines of length LENGTH,
# where LENGTH may be 80, 132 (default), or 'unlimited' for longer lines. The
# result is either "unknown", or the actual compiler flag required to to accept
# lines of length LENGTH.
#
# If successful, runs ACTION-IF-SUCCESS (defaults to appending the result to
# FCFLAGS), otherwise runs ACTION-IF-FAILURE (defaults to failing with an error
# message).
#
# The flag is cached in the acx_cv_fc_line_length_[]LENGTH.
#
# The implementation is based on the standard Autoconf macro AC_FC_LINE_LENGTH,
# which is monkey-patched to reduce the number of LANG switches and to avoid
# false negative results with the GFortran '-fimplicit-none' flag.
#
AC_DEFUN([ACX_FC_LINE_LENGTH],
  [AC_LANG_ASSERT([Fortran])dnl
dnl Fail instead of warning:
   m4_bmatch(m4_default([$1], [132]),
     [unlimited\|132\|80], [],
     [m4_fatal([Invalid LENGTH argument for $0: '$1'])])dnl
dnl Monkey-patch AC_FC_LINE_LENGTH:
   m4_pushdef([acx_cache_var], [acx_cv_fc_line_length_$1])dnl
   m4_pushdef([ac_cv_fc_line_length], [acx_cache_var])dnl
   m4_pushdef([AC_FC_LINE_LENGTH],
     m4_bpatsubsts(m4_dquote(m4_defn([AC_FC_LINE_LENGTH])),
       [\$ac_fc_line_length_test$], [\&
        implicit integer (a)],
       [AC_LANG_P\(OP\|USH\)(\[?Fortran\]?)]))dnl
dnl This macro does not have a special meaning for the value 'none'
dnl but AC_FC_LINE_LENGTH does. To account for the difference, we need to know
dnl whether 'none' came from the cache variable and set the cache variable to
dnl 'none' if it is set to an empty string:
   acx_fc_line_length_none_in_cache=no
   AS_VAR_SET_IF([acx_cache_var],
     [AS_CASE([$acx_cache_var],
        [none], [acx_fc_line_length_none_in_cache=yes],
        [''], [acx_cache_var=none])])
dnl AC_FC_LINE_LENGTH changes the FCFLAGS, which we do not want:
   acx_save_FCFLAGS=$FCFLAGS
   m4_quote(AC_FC_LINE_LENGTH([$1], [], [:]))
   FCFLAGS=$acx_save_FCFLAGS
   m4_popdef([AC_FC_LINE_LENGTH])dnl
   m4_popdef([ac_cv_fc_line_length])dnl
dnl Set the cache variable to an empty string if the value 'none' was set by
dnl AC_FC_LINE_LENGTH but not if the variable had it before the expansion:
   AS_IF([test "x$acx_cache_var$acx_fc_line_length_none_in_cache" = xnoneno],
     [acx_cache_var=])
   AS_VAR_IF([acx_cache_var], [unknown],
     [m4_default([$3],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
accept m4_default([$1], [132]) column source lines])])],
     [m4_default([$2],
        [AS_IF([test -n "$acx_cache_var"],
           [AS_VAR_APPEND([FCFLAGS], [" $acx_cache_var"])])])])
   m4_popdef([acx_cache_var])])
