!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Offers everything which is publicly available when dealing with dispersions.
!!
module dispersions
  use dispiface
  use dispuff_module
  use dispuffdata
  use dispslaterkirkw
  use disp_ddmc
  use disp_ddmc_data
#ifdef WITH_DFTD3  
  use dispdftd3_module
#endif
  implicit none
  public

  type :: DispersionInp
    type(DispUffInp), allocatable :: uff
    type(DispSlaKirkInp), allocatable :: slakirk
    type(DispDDMCInp), allocatable :: ddmc
#ifdef WITH_DFTD3    
    type(DispDftD3Inp), allocatable :: dftd3
#endif    
  end type DispersionInp
  
end module dispersions
