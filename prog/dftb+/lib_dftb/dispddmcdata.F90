!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains parameters for the dDMC dispersion correction
module disp_dDMC_data
  use Accuracy
  use Message, only : error
  use CharManip, only : tolower
  use Constants
  implicit none
  private
  
  public :: getDDMCValues
  
  ! Contains dDMC data (chemical symbol, polarizability, c6 for free atoms, vdW radii)
  type TDDMC
    character(2) :: symbol
    integer  :: Z                ! Atomic number
    integer  :: incharge         ! Charge in the innershells
    real(dp) :: p                ! Polarizability
    real(dp) :: c6_d3            ! c6 for free atom as provided by Grimme in the D3 correction
    real(dp) :: vdWr             ! vdW radii. If 1e+4 means that the data is missing
  end type TDDMC
  
  type(TDDMC), parameter :: database(54) = (/ &
      & TDDMC("h ",  1,  0,  0.66679_dp,    7.590_dp,   1.1_dp), &
      & TDDMC("he",  2,  0,  0.20505_dp,    1.560_dp,   1.4_dp), &
      & TDDMC("li",  3,  0, 24.33000_dp, 1163.450_dp, 1e+04_dp), &
      & TDDMC("be",  4,  2,  5.60000_dp,  257.490_dp, 1e+04_dp), &
      & TDDMC("b ",  5,  2,  3.03000_dp,  107.180_dp, 1e+04_dp), &
      & TDDMC("c ",  6,  2,  1.67000_dp,   49.110_dp,  1.70_dp), &
      & TDDMC("n ",  7,  2,  1.10000_dp,   25.270_dp,  1.55_dp), &
      & TDDMC("o ",  8,  2,  0.80200_dp,   15.510_dp,  1.52_dp), &
      & TDDMC("f ",  9,  2,  0.55700_dp,    9.690_dp, 1e+04_dp), &
      & TDDMC("ne", 10,  2,  0.39560_dp,    6.290_dp, 1e+04_dp), &
      & TDDMC("na", 11,  2, 24.11000_dp, 1608.030_dp, 1e+04_dp), &
      & TDDMC("mg", 12, 10, 10.60000_dp,  683.380_dp, 1e+04_dp), &
      & TDDMC("al", 13, 10,  6.80000_dp,  540.540_dp, 1e+04_dp), &
      & TDDMC("si", 14, 10,  5.38000_dp,  317.860_dp, 1e+04_dp), &
      & TDDMC("p ", 15, 10,  3.63000_dp,  191.690_dp,   1.8_dp), &
      & TDDMC("s ", 16, 10,  2.90000_dp,  134.010_dp,   1.8_dp), &
      & TDDMC("cl", 17, 10,  2.18000_dp,   92.350_dp, 1e+04_dp), &
      & TDDMC("ar", 18, 10,  1.64110_dp,   64.650_dp, 1e+04_dp), &
      & TDDMC("k ", 19, 10, 43.06000_dp, 4983.500_dp, 1e+04_dp), &
      & TDDMC("ca", 20, 18, 22.80000_dp, 2352.690_dp, 1e+04_dp), &
      & TDDMC("sc", 21, 18, 17.80000_dp, 1522.470_dp, 1e+04_dp), &
      & TDDMC("ti", 22, 18, 14.60000_dp, 1361.920_dp, 1e+04_dp), &
      & TDDMC("v ", 23, 18, 12.40000_dp, 1116.100_dp, 1e+04_dp), &
      & TDDMC("cr", 24, 18, 11.60000_dp,  690.740_dp, 1e+04_dp), &
      & TDDMC("mn", 25, 18,  9.40000_dp,  802.750_dp, 1e+04_dp), &
      & TDDMC("fe", 26, 18,  8.40000_dp,  491.330_dp, 1e+04_dp), &
      & TDDMC("co", 27, 18,  7.50000_dp,  532.780_dp, 1e+04_dp), &
      & TDDMC("ni", 28, 36,  6.80000_dp,  574.740_dp, 1e+04_dp), &
      & TDDMC("cu", 29, 36,  6.20000_dp,  337.180_dp, 1e+04_dp), &
      & TDDMC("zn", 30, 36,  5.75000_dp,  340.520_dp, 1e+04_dp), &
      & TDDMC("ga", 31, 36,  8.12000_dp,  483.750_dp, 1e+04_dp), &
      & TDDMC("ge", 32, 36,  6.07000_dp,  363.550_dp, 1e+04_dp), &
      & TDDMC("as", 33, 36,  4.31000_dp,  262.950_dp, 1e+04_dp), &
      & TDDMC("se", 34, 36,  3.77000_dp,  213.670_dp, 1e+04_dp), &
      & TDDMC("br", 35, 36,  3.05000_dp,  167.130_dp, 1e+04_dp), &
      & TDDMC("kr", 36, 36,  2.48400_dp,  130.400_dp, 1e+04_dp), &
      & TDDMC("rb", 37, 36, 47.30000_dp, 6138.780_dp, 1e+04_dp), &
      & TDDMC("sr", 38, 36, 27.60000_dp, 3381.370_dp, 1e+04_dp), &
      & TDDMC("y ", 39, 36, 22.70000_dp, 2365.890_dp, 1e+04_dp), &
      & TDDMC("zr", 40, 36, 17.90000_dp, 1822.720_dp, 1e+04_dp), &
      & TDDMC("nb", 41, 36, 15.70000_dp, 1475.250_dp, 1e+04_dp), &
      & TDDMC("mo", 42, 36, 12.80000_dp,  845.900_dp, 1e+04_dp), &
      & TDDMC("tc", 43, 36, 11.40000_dp, 1067.020_dp, 1e+04_dp), &
      & TDDMC("ru", 44, 36,  9.60000_dp,  598.200_dp, 1e+04_dp), &
      & TDDMC("rh", 45, 36,  8.60000_dp,  713.940_dp, 1e+04_dp), &
      & TDDMC("pd", 46, 36,  4.80000_dp,  608.500_dp, 1e+04_dp), &
      & TDDMC("ag", 47, 36,  7.20000_dp,  426.750_dp, 1e+04_dp), &
      & TDDMC("cd", 48, 36,  7.36000_dp,  468.190_dp, 1e+04_dp), &
      & TDDMC("in", 49, 36, 10.20000_dp,  757.740_dp, 1e+04_dp), &
      & TDDMC("sn", 50, 36,  7.70000_dp,  627.570_dp, 1e+04_dp), &
      & TDDMC("sb", 51, 36,  6.60000_dp,  492.940_dp, 1e+04_dp), &
      & TDDMC("te", 52, 36,  5.50000_dp,  425.540_dp, 1e+04_dp), &
      & TDDMC("i ", 53, 36,  5.35000_dp,  351.970_dp, 1e+04_dp), &
      & TDDMC("xe", 54, 36,  4.04400_dp,  290.220_dp, 1e+04_dp) /)

contains

  subroutine getDDMCValues(name, c6free, polar, vdWr, Z, incharge, found)
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: polar, vdWr, c6free
    integer, intent(out) :: Z
    integer, intent(out), optional :: incharge
    logical, intent(out), optional :: found

    character(2) :: symbol
    integer :: ii

    symbol = trim(tolower(name))
    do ii = 1, size(database)
      if (database(ii)%symbol == symbol) then
        polar = database(ii)%p * AA__Bohr**3
        vdWr = database(ii)%vdWr * AA__Bohr
        Z = database(ii)%Z
        if (present(incharge)) then
          incharge = database(ii)%incharge
        end if
        c6free = database(ii)%c6_d3
        if (vdWr > .9e+04_dp) then
          call error("DDMC database search for element '" // trim(name) &
               &// "' van der Waals radius not found!")
        end if
        if (present(found)) then
          found = .true.
        end if
        return
      end if
    end do
    if (present(found)) then
      found = .false.
    else
      call error("DDMC database search for element '" // trim(name) &
          &// "' failed")
    end if

  end subroutine getDDMCValues
    
end module disp_dDMC_data
