# ACX_COMPARE_VERSION(VERSION_A,
#                     OP,
#                     VERSION_B,
#                     [ACTION-IF-TRUE],
#                     [ACTION-IF-FALSE])
# -----------------------------------------------------------------------------
# Wraps AX_COMPARE_VERSION
# (https://www.gnu.org/software/autoconf-archive/ax_compare_version.html) in
# order to suppress the redundant expansion of AC_PROG_AWK, which is needed
# only if OP is either 'eq0' or 'ne0'.
#
# See the documentation on AX_COMPARE_VERSION for more details.
#
AC_DEFUN([ACX_COMPARE_VERSION],
  [AC_PROVIDE_IFELSE([AC_PROG_AWK], [],
     [m4_bmatch([m4_substr([$2],[2])],
        [0], [],
        [dnl AC_PROG_AWK is redundant:
dnl Leave the clean-up (undo) witness:
         m4_pushdef([acx_compare_version_undo])dnl
dnl Pretend that AC_PROG_AWK is provided:
         m4_pushdef([m4_provide(AC_PROG_AWK)])])])dnl
   AX_COMPARE_VERSION($@)
   m4_ifdef([acx_compare_version_undo],
     [dnl Clean up:
      m4_popdef([m4_provide(AC_PROG_AWK)])dnl
      m4_popdef([acx_compare_version_undo])])])
