<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/buildbot/dwd_nec_patches)
to see this file rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

This directory contains a patch that we have to apply to RTE-RRTMGP when
building with the Intel Fortran Compiler Classic 2021.5 on Levante. The patch
introduces a performance-hitting workaround for the vectorization problem
reported [here](https://github.com/earth-system-radiation/rte-rrtmgp/issues/159).
The workaround has also been submitted as a pull request to RTE-RRTMGP (see
[here](https://github.com/earth-system-radiation/rte-rrtmgp/pull/170)). The
known affected compilers are Intel Fortran Compiler Classic 2021.4, 2021.5 and
2022.1. The problem has not yet been confirmed for the Intel Fortran Compiler
oneAPI compilers (a.k.a `ifx`), which might be the case though. It is not clear
when the compiler bug will be fixed (see
[here](https://community.intel.com/t5/Intel-Fortran-Compiler/Compiler-vectorization-bug/m-p/1362591)).
