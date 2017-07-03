!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dispersion correction based on Mulliken charges
!!
!! \ref DOI: 10.1002/qua.24887
!! \note Original code at https://github.com/grhawk/dDMC/
!!
module disp_ddmc
  use assert
  use accuracy
  use simplealgebra, only : determinant33
  use lapackroutines, only : matinv
  use periodic, only: TNeighborList, getNrOfNeighborsForAll, getLatticePoints
  use constants, only : pi
  use dispiface
  use dispcommon
  use message, only: error
  implicit none
  private

  public :: DispDDMCInp, DispDDMC, DispDDMC_init


  !!* Contains the initialisation data for the dDMC module
  type :: DispDDMCInp
    real(dp), allocatable :: params(:) !* parameter a,b,s of the damping function
    real(dp), allocatable :: c6free(:)  !* c6 of free atoms from ddmc_data module
    real(dp), allocatable :: vdWr(:)    !* van der Waals radii from ddmc_data module
    real(dp), allocatable :: polar(:)   !* polarizabilities from ddmc_data module
    integer, allocatable :: incharge(:) !* inner electrons from ddmc_data module
    integer, allocatable :: Z(:)    !* atomic number from ddmc_data module
  end type DispDDMCInp


  !!* Data for the dDMC type dispersion
  type, extends(DispersionIface) :: DispDDMC
    private
    real(dp), allocatable :: c6(:,:)  !* array to be allocated to contain all the c6_ij. (It is a
                                      !* waste of memory since only half will be useful but I am
                                      !* following the dftb+ approach)
    real(dp), allocatable :: c6aim(:) !* c6 of aim from ddmc_data module
    real(dp), allocatable :: c6free(:)
    real(dp), allocatable :: Z(:)
    real(dp), allocatable :: incharge(:)
    real(dp), allocatable :: vdWr(:)
    real(dp), dimension(3):: params
    real(dp), allocatable :: bi(:)     ! needed to bring around the bi values
    real(dp), allocatable :: bi_free(:)                !* Defined as croot(1/polar)
    real(dp), allocatable :: mulcharge(:) ! Much more like mulliken population
    !    real(dp), allocatable :: c6_ij(:,:)
    integer :: nAtom                  !* Nr. of atoms (without images)
    real(dp), allocatable :: energies(:)  !* Energies
    real(dp), allocatable :: gradients(:,:) !* Gradients (3, nAtom)
    real(dp) :: stress(3,3)           !* stress tensor components
    logical :: tPeriodic              !* If system is periodic
    real(dp) :: rCutoff               !* Real space cutoff
    real(dp) :: gCutoff               !* Reciprocal space cutoff
    real(dp) :: dampCutoff            !* Cutoff, where damping function = 1
    real(dp) :: eta                   !* Periodic summation parameter
    real(dp) :: vol                   !* Volume of the unit cell
    real(dp) :: maxR                  !????
    real(dp) :: c6sum_init            !* The c6sum computed in init could be needed afterwords when
                                      !* using sockets
    real(dp), allocatable:: gLatPoint(:,:)  !* Temporary dirty solution
    logical :: coordsUpdated          !* If first coordinate update done
    logical :: chargsUpdated          !* If first Mulliken charges update done
    logical :: tInit = .false.
  contains
    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress
    procedure :: getRCutoff
  end type DispDDMC



contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!  Public routines
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Initializes a dDMC instance.
  !!* @param sf Initialised instance at return.
  !!* @param inp Input parameters for dDMC
  !!* @param latVecs Lattice vectors, if the system is periodic
  !!* @param recVecs Reciprocal vectors, if the system is periodic
  !!* @param vol Volume of the unit cell, if the system is periodic
  subroutine DispDDMC_init(this, inp, nAtom, latVecs)
    type(DispDDMC), intent(inout) :: this
    type(DispDDMCInp), intent(in) :: inp
    integer, intent(in) :: nAtom
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iAt1, iAt2
    real(dp) :: recVecs(3, 3), invRecVecs(3, 3)
    real(dp) :: tol, c6max
    real(dp), parameter :: c6extimScale = 1.1

    @:ASSERT(.not. this%tInit)
    @:ASSERT(size(inp%params) == 3)
    @:ASSERT(size(inp%c6free) == nAtom)
    @:ASSERT(size(inp%vdWr) == nAtom)
    @:ASSERT(size(inp%polar) == nAtom)
    @:ASSERT(size(inp%incharge) == nAtom)
    @:ASSERT(size(inp%Z) == nAtom)
  #:call ASSERT_CODE
    if (present(latVecs)) then)
      @:ASSERT(all(shape(latVecs) == (/ 3, 3 /)))
    end if
  #:endcall ASSERT_CODE

    allocate(this%mulcharge(nAtom))
    allocate(this%c6aim(nAtom))
    allocate(this%c6(nAtom, nAtom))
    allocate(this%bi_free(nAtom))
    allocate(this%bi(nAtom))
    allocate(this%c6free(nAtom))
    allocate(this%Z(nAtom))
    allocate(this%incharge(nAtom))
    allocate(this%vdWr(nAtom))

    this%c6 = 0.0_dp
    this%bi_free = (inp%polar)**(-1.0_dp/3.0_dp)
    this%c6free = inp%c6free
    this%Z = inp%Z
    this%incharge = inp%incharge
    this%vdWr = inp%vdWr
    this%params = inp%params

    this%nAtom = nAtom
    this%rCutoff = 0.0_dp
    this%c6aim = 0.0_dp
    this%maxR = 0.0_dp
    this%c6sum_init = 0.0_dp
    this%bi = this%bi_free
    tol = epsilon(1.0_dp)

    ! Since I do not have the charges yet, all the c6 dependent distances are
    ! estimated using the biggest c6free multiplied by 1.1
    c6max = c6mix_(maxval(inp%c6free)*c6extimScale, maxval(inp%c6free)*c6extimScale)
    this%rCutoff = (c6max/tolDispersion)**(1.0_dp/6.0_dp) ! This comes from UFF

    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      this%vol = abs(determinant33(latVecs))
      invRecVecs(:,:) = latVecs / (2.0_dp * pi)
      recVecs(:,:) = transpose(invRecVecs)
      call matinv(recVecs)
      !! Scaling down optimal eta (as suggested in the literature) is purely
      !! empirical, it reduces the real space summation, and seems to yield
      !! shorter execution times. (It doesn't influence the result.)
      this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
      !! Since it would not be possible to have the c6aim, I am estimating the
      !! c6sum using the c6free.
      do iAt1 = 1,nAtom-1
        do iAt2 = iAt1+1, nAtom
          this%c6sum_init = this%c6sum_init + c6mix_(inp%C6free(iAt1), inp%C6free(iAt2))
        end do
      end do

      this%rCutoff = getMaxRDispersion(this%eta, this%c6sum_init, this%vol, &
        &tolDispersion)
      !! Cutoff, beyond which dispersion is purely 1/r^6 without damping
      this%dampCutoff = getDampCutoff_(inp%params, maxval(inp%vdWr), maxval(this%bi), tolDispersion)

      this%rCutoff = max(this%rCutoff, this%dampCutoff)
      this%gCutoff = getMaxGDispersion(this%eta, this%c6sum_init, tolDispersion)
      call getLatticePoints(this%gLatPoint, recVecs, invRecVecs, &
        &this%gCutoff, onlyInside=.true., reduceByInversion=.true., &
        &withoutOrigin=.true.)
      this%gLatPoint = matmul(recVecs, this%gLatPoint)
    end if

    allocate(this%energies(this%nAtom))
    allocate(this%gradients(3, this%nAtom))
    this%coordsUpdated = .false.
    this%tInit = .true.

  end subroutine DispDDMC_init


  !!* Notifies the objects about changed coordinates
  !!* @param sf Object instance.
  !!* @param inp Input parameters for dDMC
  !!* @param neigh Current neighbor list mapping
  !!* @param img2CentCell Current mapping of periodic images to the central cell
  !!* @param coords Current coordinates of the atoms
  subroutine updateCoords(this, neigh, img2CentCell, coords, species0, charges)
    class(DispDDMC), intent(inout) :: this
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:), charges(:,:,:)
    integer, intent(in) :: species0(:)

    integer :: iAt1, iAt2
    integer, allocatable :: nNeighReal(:) ! Neighbors for real space summation
    integer, allocatable :: nNeighDamp(:) ! Nr. of neighbors with damping

    @:ASSERT(this%tInit)

    !! Charge dimensions are: orbital, nAtom, spin
    do iAt1 = 1, this%nAtom
      this%mulcharge(iAt1) = sum(charges(:,iAt1,:))
    end do
    !! Compute the c6 for the atoms-in-molecule
    this%c6aim = ((this%mulcharge + this%incharge) / this%Z)**2 * this%c6free
    !! Scale the basym values with the inverse of the N/Z**1/3
    this%bi(:) = this%bi_free(:) * ( this%Z(:) / ( this%mulcharge(:) + this%incharge(:) )) &
        & **(1.0_dp / 3.0_dp)

    allocate(nNeighReal(this%nAtom))
    call getNrOfNeighborsForAll(nNeighReal, neigh, this%rCutoff)
    this%energies(:) = 0.0_dp
    this%gradients(:,:) = 0.0_dp

    if (this%tPeriodic) then
      do iAt1 = 1, this%nAtom
        do iAt2 = 1, iAt1
          this%c6(iAt2, iAt1) = c6mix_(this%c6aim(iAt1), this%c6aim(iAt2))
          if( iAt1 /= iAt2 )then
            this%c6(iAt1, iAt2) = this%c6(iAt2, iAt1)
          end if
        end do
      end do
      !! Make Ewald summation for a pure 1/r^6 interaction
      call addDispEGr_per_atom(this%nAtom, coords, nNeighReal, &
           &neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6, this%eta, &
           &this%vol, this%gLatPoint, this%energies, this%gradients,this%stress)
      !! Correct those terms, where damping is important
      allocate(nNeighDamp(this%nAtom))
      call getNrOfNeighborsForAll(nNeighDamp, neigh, this%dampCutoff)
      call addDispEnergyAndGrad_cluster_(this%nAtom, coords, nNeighDamp, &
        &neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6aim, this%VdWr, &
        &this%energies, this%gradients, this%params, this%bi, dampCorrection=-1.0_dp)
    else
      call addDispEnergyAndGrad_cluster_(this%nAtom, coords, nNeighReal, &
        &neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6aim, this%VdWr, &
        &this%energies, this%gradients, this%params, this%bi, dampCorrection=0.0_dp)
    end if
    this%coordsUpdated = .true.

  end subroutine updateCoords


  !!* Notifies object about updated lattice vectors
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal vectors
  !!* @param vol  New unit cell volume
  subroutine updateLatVecs(this, latVecs)
    class(DispDDMC), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:)

    real(dp) :: invRecVecs(3,3), recVecs(3,3)
    real(dp) :: c6sum

    this%vol = abs(determinant33(latVecs))
    invRecVecs(:,:) = latVecs / (2.0_dp * pi)
    recVecs(:,:) = transpose(invRecVecs)
    call matinv(recVecs)
    this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
    ! When using socket, at the first cycle, c6 are setted to 0. In that case extimate the
    ! c6sum as done in the init (see above). This could introduce a numerical error in the first
    ! communication to the socket since the cutoff are extimated instead of computed exactly.
    c6sum = sum(abs(this%c6))
    if (c6sum < 1.E-8_dp) then
      c6sum = this%c6sum_init
    end if

    this%rCutoff = getMaxRDispersion(this%eta, c6sum, this%vol, tolDispersion)
    !! Cutoff, beyond which dispersion is purely 1/r^6 without damping
    this%dampCutoff = getDampCutoff_(this%params, maxval(this%vdWr), maxval(this%bi), tolDispersion)
    this%rCutoff = max(this%rCutoff, this%dampCutoff)
    this%gCutoff = getMaxGDispersion(this%eta, c6sum, tolDispersion)
      call getLatticePoints(this%gLatPoint, recVecs, invRecVecs, &
        &this%gCutoff, onlyInside=.true., reduceByInversion=.true., &
        &withoutOrigin=.true.)
    this%gLatPoint = matmul(recVecs, this%gLatPoint)
    this%coordsUpdated = .false.

  end subroutine updateLatVecs


  !!* Returns the atomic resolved energies due to the dispersion.
  !!* @param sf Object instance
  !!* @param energies Contains the atomic energy contributions on exit.
  subroutine getEnergies(this, energies)
    class(DispDDMC), intent(inout) :: this
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%tInit)
    @:ASSERT(this%coordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies(:)

  end subroutine getEnergies


  !!* Adds the atomic gradients to the provided vector.
  !!* @param sf Object instance.
  !!* @param gradients The vector to increase by the gradients.
  subroutine addGradients(this, gradients)
    class(DispDDMC), intent(inout) :: this
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%tInit)
    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(gradients) == [ 3, this%nAtom ]))

    gradients(:,:) = gradients(:,:) + this%gradients(:,:)

  end subroutine addGradients

  !> Returns the stress tensor.
  !!
  !! \note The stress tensor is not calculated for this dispersion model
  !!     so the program is stopped, if this method is called.
  !!
  !! \param stress tensor from the dispersion
  !!
  subroutine getStress(this, stress)
    class(DispDDMC), intent(inout) :: this
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(stress) == [3,3]))

    stress = this%stress
    ! This should be tested better...

  end subroutine getStress


  !!* Estimates the real space cutoff of the dispersion interaction.
  !!* @param sf Object instance
  !!* @return Cutoff
  function getRCutoff(this) result(cutoff)
    class(DispDDMC), intent(inout) :: this
    real(dp) :: cutoff

    @:ASSERT(this%tInit)
    cutoff = this%rCutoff

  end function getRCutoff


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Private routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Adds the energy per atom and the gradients for the cluster case
  !!* @param nAtom Nr. of atoms (without periodic images)
  !!* @param coords Coordinates of the atoms (including images)
  !!* @param nNeighbors Nr. of neighbors for each atom
  !!* @param iNeighbor Neighborlist.
  !!* @param neighDist2 Square distances of the neighbours.
  !!* @param img2CentCell Mapping into the central cell.
  !!* @param c6 Van der Waals coefficients (nAtom, nAtom)
  !!* @param rVdW2 Scaled inverse van der Waals radii (nAtom, nAtom)
  !!* @param energies Updated energy vector at return
  !!* @param gradients Updated gradient vector at return
  !!* @param dampCorrection Adds the provided value to the damping function
  !!*   (use -1.0 to sum up damped 1/r^6 terms and subtract pure 1/r^6 ones, in
  !!*   order to correct periodic Ewald sum for the short range damped terms.)
  subroutine addDispEnergyAndGrad_cluster_(nAtom, coords, nNeighbors, &
    &iNeighbor, neighDist2, img2CentCell, c6, rVdW, energies, gradients, &
    &prms, b, dampCorrection)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:)
    real(dp), intent(in) :: rVdW(:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(inout) :: energies(:)
    real(dp), intent(inout) :: gradients(:,:)
    real(dp), intent(in) :: prms(3)
    real(dp), intent(in), optional :: dampCorrection

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: dist2, dist, rTmp
    real(dp) :: diff(3), gr(3)
    real(dp) :: corr

    if( present(dampCorrection) )then
      corr = dampCorrection
    else
      corr = 0.0_dp
    end if

    !! Cluster case => explicit sum of the contributions
    !! NOTE: the cluster summation also (ab)used in the periodic case, neighbors
    !! may go over the cell boundary -> img2CentCell needed for folding back.
    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if(( c6(iAt2f) == 0.0_dp ) .and. ( c6(iAt1) == 0.0_dp ))then
          call error("C6 terms zero for both partners")
        end if
        dist2 = neighDist2(iNeigh, iAt1)
        if (dist2 > 0.0) then
          !! Energy
          dist = sqrt(dist2)
          rTmp = -0.5_dp * (getDampFactor_(dist, prms, rVdW(iAt1), rVdW(iAt2f), b(iAt1), b(iAt2f)) &
              &+ corr) * c6mix_(c6(iAt2f), c6(iAt1)) / dist**6
          energies(iAt1) = energies(iAt1) + rTmp
          if (iAt1 /= iAt2f) then
            energies(iAt2f) = energies(iAt2f) + rTmp
          end if
          ! Gradients
          ! The correction used by the periodic systems needs always to be the derivative
          ! of the energy correction. That is why we have the 6R**(-7) term.
          diff(:) = (coords(:,iAt1) - coords(:,iAt2))
          ! Comment from Alberto Fabrizio => Be careful: a minus sign is already
          ! included into getDampGrad_. Also the Pure Ewald sign must be positive
          ! to subtract (see comment before)
          gr(:) = c6mix_(c6(iAt2f), c6(iAt1)) * diff(:) * &
              & (getDampGrad_(dist, prms, rVdW(iAt1), rVdW(iAt2f), b(iAt1), b(iAt2f)) &
              & + 6 * dist**(-7) * corr) * 1/dist
          gradients(:,iAt1) = gradients(:,iAt1) + gr(:)
          gradients(:,iAt2f) = gradients(:,iAt2f) - gr(:)
        end if
      end do
    end do

  end subroutine addDispEnergyAndGrad_cluster_

  !!* Returns the distance, beyond that the damping function equals approx. 1.
  !!* Bruteforce method, otherwise too complicated
  !!* @param R0   longest vanDerWaals radius
  !!* @param prms damping function's parameters
  !!* @param b    highest value of b
  !!* @param tol  Tolerance value.
  !!* @return     cutoff
  function getDampCutoff_(prms, R0, b, tol) result(xx)
    real(dp), intent(in) :: R0, b, tol
    real(dp), intent(in) :: prms(3) ! b0, a, s
    real(dp), parameter :: increase = 0.05_dp
    real(dp) :: xx

    xx = 1.1_dp*prms(2)*R0
    do while (abs(getDampFactor_(xx, prms, R0, R0, b, b) - 1._dp) > tol)
      xx = xx + increase
    end do

  end function getDampCutoff_

  !!* Returns the value of the damping factor at a given distance x
  !!* @param x    distance between the atoms
  !!* @param R0i  vDW radius of the atom i
  !!* @param R0j  vDW radius of the atom j
  !!* @param prms three damping function's parameters
  !!* @param b_ii TT damping factor atom i
  !!* @param b_jj TT damping factor atom j
  function getDampFactor_(x, prms, R0i, R0j, b_ii, b_jj) result(df)
    real(dp), intent(IN) :: x, R0i, R0j, b_ii, b_jj
    real(dp), intent(IN) :: prms(3)
    real(dp) :: df, Fd, TT, b_ij, bx

    b_ij = bmix_(prms(1)*b_ii, prms(1)*b_jj)

    Fd = 0.5*( 1.d0 + tanh( prms(3) *  ( x / ( prms(2) * ( R0i + r0j )) - 1.d0 ) ) )

    bx = b_ij * x

    TT = 1.d0 - &
      & ( exp( -bx ) * (1.d0 + bx + (bx)**2.d0/2.d0 + (bx)**3.d0/6.d0 + &
      & (bx)**4.d0/24.d0 + (bx)**5.d0/120.d0 + (bx)**6.d0/720.d0) )

    df = TT*Fd

  end function getDampFactor_


  function getDampGrad_(R, prms, R0i, R0j, b_ii, b_jj) result(dfp)
    real(dp), intent(IN) :: R, R0i, R0j, b_ii, b_jj
    real(dp), intent(IN) :: prms(3)
    real(dp) :: dfp, R0, bij!, a, b, s, x

    bij = bmix_(prms(1)*b_ii, prms(1)*b_jj)
    R0 = R0i + R0j

    ! s = prms(3)
    ! a = prms(2)
    ! b = prms(1)
    ! x = R

    !First Derivative computed from Alberto Fabrizio's expression

       dfp = -0.5*prms(3)/(R0*prms(2))* &
           & (1-(tanh(prms(3)*(R/(R0*prms(2))-1)))**2)* &
           & (1-exp(-bij*R)*(1+bij*R+(bij**2*R**2)/2+(bij**3*R**3)/6+(bij**4*R**4)/24 &
           & +(bij**5*R**5)/120+ &
           & (bij**6*R**6)/720))/(R**6)-0.5*(1+tanh(prms(3)*(R/(R0* &
           & prms(2))-1)))*0.001388889*bij**7*R**6*exp(-bij*R)/(R**6)+ &
           & 0.5*(1+tanh(prms(3)*(R/(R0*prms(2))-1)))* &
           & (1-exp(-bij*R)*(1+bij*R+(bij**2*R**2)/2+(bij**3*R**3)/6+(bij**4*R**4)/24 &
           & +(bij**5*R**5)/120+ &
           & (bij**6*R**6)/720))*6/(R**7)

  end function getDampGrad_

  function bmix_(bi,bj) result(bmix)
    real(dp),intent(IN) :: bi, bj
    real(dp) :: bmix

    bmix = 2 * bi * bj / (bi + bj)

  end function bmix_

  function c6mix_(c6_i, c6_j) result(cmix)
    real(dp),intent(IN) :: c6_i, c6_j
    real(dp) :: cmix

    cmix = 2   * c6_i * c6_j / ( c6_i + c6_j )

  end function c6mix_

end module disp_ddmc
