<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/buildbot/dwd_nec_patches)
to see this file rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

This directory contains a patch that we have to apply to YAXT on the NEC machine
at DWD. The patch modifies `test_xmap_all2all_fail_run` test of YAXT, which
expects a failure from the `mpirun` process with exit code `3`. The existing
workaround implemented in YAXT for such cases is unfortunately not enough.
Therefore, the test is modified to pass if the `mpirun` program terminates with
any non-zero exit code. The problem has been reported upstream (see
[here](https://gitlab.dkrz.de/dkrz-sw/yaxt/-/issues/2)) and seems to be
addressed in the `master` branch (see
[here](https://gitlab.dkrz.de/dkrz-sw/yaxt/-/commit/5a7e9e42ae85fec4d900efdc19690ffe0c772956)).
