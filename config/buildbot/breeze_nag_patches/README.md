<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/buildbot/dwd_nec_patches)
to see this file rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

This directory contains a patch that we have to apply to RTE-RRTMGP when
building with the NAG compiler and the DACE modules for data assimilation
enabled. The issue has been
[reported upstream](https://github.com/earth-system-radiation/rte-rrtmgp/issues/158)
and is already
[fixed in the `develop` branch](https://github.com/earth-system-radiation/rte-rrtmgp/pull/165).
This patch can be removed once the changes find their way into the
[`autoconf`](https://github.com/earth-system-radiation/rte-rrtmgp/tree/autoconf)
branch.
