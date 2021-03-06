!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains subroutines to calculate repulsive pair contributions to energy
!!* and forces
module repulsive
  use assert
  use accuracy, only : dp
  use repcont
  implicit none
  private

  public :: getERep, getERepDeriv

  interface getERep
     module procedure getERep_total
     module procedure getERep_atoms
  end interface

contains

  !!* Subroutine for calculating total energy contribution of the repulsives.
  !!* @param reslt Total energy contribution.
  !!* @param coords coordinates (x,y,z, all atoms including possible images)
  !!* @param nNeighbors Number of neighbors for atoms in the central cell
  !!* @param iNeighbors Index of neighbors for a given atom.
  !!* @param species Species of atoms in the central cell.
  !!* @param img2CentCell Index of each atom in the central cell which the atom
  !!*   is mapped on.
  !!* @param repCont Container for repulsive potentials.
  subroutine getERep_total(reslt, coords, nNeighbors, iNeighbors, species, &
      &img2CentCell, repCont)
    real(dp), intent(out) :: reslt
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbors(0:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: img2CentCell(:)
    type(ORepCont), intent(in) :: repCont

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), dist, intermed

    reslt = 0.0_dp
    do iAt1 = 1, size(nNeighbors)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbors(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
        if (iAt2f /= iAt1) then
          reslt = reslt + intermed
        else
          reslt = reslt + 0.5_dp * intermed
        end if
      end do
    end do

  end subroutine getERep_total



  !!* Subroutine for repulsive energy contributions for each atom
  !!* @param reslt Energy for each atom.
  !!* @param coords coordinates (x,y,z, all atoms including possible images)
  !!* @param nNeighbors Number of neighbors for atoms in the central cell
  !!* @param iNeighbors Index of neighbors for a given atom.
  !!* @param species Species of atoms in the central cell.
  !!* @param repCont Container for repulsive potentials.
  !!* @param img2CentCell Index of each atom in the central cell which the atom
  !!*   is mapped on.
  subroutine getERep_atoms(reslt, coords, nNeighbors, iNeighbors, species,&
      &repCont, img2CentCell)
    real(dp), intent(out) :: reslt(:)
    real(dp), intent(in)  :: coords(:,:)
    integer, intent(in)   :: nNeighbors(:)
    integer, intent(in)   :: iNeighbors(0:,:)
    integer, intent(in)   :: species(:)
    type(ORepCont), intent(in) :: repCont
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), dist, intermed

    @:ASSERT(size(reslt) == size(nNeighbors))

    reslt(:) = 0.0_dp
    do iAt1 = 1, size(nNeighbors)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbors(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
        reslt(iAt1) = reslt(iAt1) + 0.5_dp * intermed
        if (iAt2f /= iAt1) then
          reslt(iAt2f) = reslt(iAt2f) + 0.5_dp * intermed
        end if
      end do
    end do

  end subroutine getERep_atoms



  !!* Subroutine for force contributions of the repulsives.
  !!* @param reslt Energy for each atom.
  !!* @param coords coordinates (x,y,z, all atoms including possible images)
  !!* @param nNeighbors Number of neighbors for atoms in the central cell
  !!* @param iNeighbors Index of neighbors for a given atom.
  !!* @param species Species of atoms in the central cell.
  !!* @param repCont Container for repulsive potentials.
  !!* @param img2CentCell Index of each atom in the central cell which the atom
  !!*   is mapped on.
  subroutine getERepDeriv(reslt, coords, nNeighbors, iNeighbors, species, &
      &repCont, img2CentCell)
    real(dp), intent(out) :: reslt(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbors(0:,:)
    integer, intent(in) :: species(:)
    type(ORepCont), intent(in) :: repCont
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), intermed(3)

    @:ASSERT(size(reslt,dim=1) == 3)

    reslt(:,:) = 0.0_dp
    do iAt1 = 1, size(nNeighbors)
      lpNeigh: do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbors(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        if (iAt2f == iAt1) then
          cycle lpNeigh
        end if
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        call getEnergyDeriv(repCont, intermed, vect, species(iAt1), &
            &species(iAt2))
        reslt(:,iAt1) = reslt(:,iAt1) + intermed(:)
        reslt(:,iAt2f) = reslt(:,iAt2f) - intermed(:)
      end do lpNeigh
    end do

  end subroutine getERepDeriv


end module repulsive
