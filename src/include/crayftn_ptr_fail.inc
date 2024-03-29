! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Cray ftn compilers 8.4 and 8.6 are known to misidentify argument INTENT
! of pointer components, i.e. will disallow changes to an array pointed to
! by a pointer component of a TYPE
! therefore this case needs to be handled specially
!
#if defined _CRAYFTN
#define PTR_INTENT(arg) INTENT(inout)
#else
#define PTR_INTENT(arg) INTENT(arg)
#endif
