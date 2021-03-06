!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains routines to calculate contributions to typical DFTB Hamiltonian
!!* parts using various generalisations of H_mu,nu = 0.5*S_mu,nu*(V_mu + V_nu)
module shift
  use assert
  use accuracy
  use commontypes

  implicit none

  private
  public :: add_shift, total_shift

  ! add shifts to a given Hamiltonian
  interface add_shift
    module procedure add_shift_atom
    module procedure add_shift_lshell
    module procedure add_shift_block
  end interface

  !!* Totals together shifts to get composites
  interface total_shift
    module procedure addatom_shell
    module procedure addshell_block
  end interface

contains

  !!* Regular atomic shift (potential is only dependent on number of atom)
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at atom sites
  subroutine add_shift_atom(ham,over,nNeighbour,iNeighbour,species,orb,iPair, &
      & nAtom,img2CentCell,shift)
    real(dp), intent(inout)     :: ham(:,:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:)

    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer :: iNeigh, iSpin, nSpin

    @:ASSERT(size(ham,dim=1)==size(over))
    @:ASSERT(size(ham,dim=2)==size(shift,dim=2))
    @:ASSERT(size(nNeighbour)==nAtom)
    @:ASSERT(size(iNeighbour,dim=2)==nAtom)
    @:ASSERT(size(species)>=maxval(iNeighbour))
    @:ASSERT(size(species)<=size(img2CentCell))
    @:ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    @:ASSERT(size(iPair,dim=2)==nAtom)
    @:ASSERT(size(shift,dim=1)==nAtom)

    nSpin = size(shift,dim=2)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    do iSpin = 1, nSpin
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) = &
              & ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) + &
              & over(iOrig+1:iOrig+nOrb2*nOrb1) * 0.5_dp * &
              & ( shift(iAt1,iSpin) + shift(iAt2f,iSpin) )
        end do
      end do
    end do

  end subroutine add_shift_atom

  !!* l-dependent shift (potential is dependent on number of atom and l-shell)
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add for each l-shell on all atom sites,
  !!* (0:lmax,1:nAtom)
  subroutine add_shift_lshell( ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift )
    real(dp), intent(inout)     :: ham(:,:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:,:)

    integer  :: iAt1, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer  :: iSh1, iSh2, iNeigh, iSpin, nSpin
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), rTmp

    @:ASSERT(size(ham,dim=1)==size(over))
    @:ASSERT(size(nNeighbour)==nAtom)
    @:ASSERT(size(iNeighbour,dim=2)==nAtom)
    @:ASSERT(size(species)>=maxval(iNeighbour))
    @:ASSERT(size(species)<=size(img2CentCell))
    @:ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    @:ASSERT(size(iPair,dim=2)==nAtom)
    @:ASSERT(size(shift,dim=1)==orb%mShell)
    @:ASSERT(size(shift,dim=2)==nAtom)
    @:ASSERT(size(ham,dim=2)==size(shift,dim=3))

    nSpin = size(shift,dim=3)

    do iSpin = 1, nSpin
      do iAt1= 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2f = img2CentCell(iNeighbour(iNeigh, iAt1))
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          do iSh1 = 1, orb%nShell(iSp1)
            rTmp = shift(iSh1, iAt1, iSpin)
            do iSh2 = 1, orb%nShell(iSp2)
              tmpH(orb%posShell(iSh2,iSp2):orb%posShell(iSh2+1,iSp2)-1, &
                  & orb%posShell(iSh1,iSp1):orb%posShell(iSh1+1,iSp1)-1) = &
                  & rTmp + shift(iSh2, iAt2f,iSpin)
            end do
          end do
          ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) = &
              & ham(iOrig+1:iOrig+nOrb2*nOrb1,ispin) + &
              & 0.5_dp * over(iOrig+1:iOrig+nOrb2*nOrb1) &
              & * reshape(tmpH(1:nOrb2, 1:nOrb1), (/nOrb2*nOrb1/))
        end do
      end do
    end do

  end subroutine add_shift_lshell

  !!* shift depending on occupation-matrix like potentials. To use this for
  !!* lm-dependent potentials, use a diagonal shift matrix
  !!* @param ham          The resulting Hamiltonian contribution.
  !!* @param over         The overlap matrix.
  !!* @param nNeighbour   Number of neighbours surrounding each atom.
  !!* @param iNeighbour   List of neighbours for each atom.
  !!* @param species      List of the species of each atom.
  !!* @param orb Contains Information about the atomic orbitals in the system
  !!* @param iPair        Indexing array for the Hamiltonian.
  !!* @parma img2CentCell Index mapping atoms onto the central cell atoms.
  !!* @param shift        Shift to add at atom sites, listed as
  !!* (0:nOrb,0:nOrb,1:nAtom)
  subroutine add_shift_block( ham,over,nNeighbour,iNeighbour,species,orb, &
      & iPair,nAtom,img2CentCell,shift )
    real(dp), intent(inout)     :: ham(:,:)
    real(dp), intent(in)        :: over(:)
    integer, intent(in)         :: nNeighbour(:)
    integer, intent(in)         :: iNeighbour(0:,:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: iPair(0:,:)
    integer, intent(in)         :: nAtom
    integer, intent(in)         :: img2CentCell(:)
    real(dp), intent(in)        :: shift(:,:,:,:)

    integer  :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer :: iNeigh, iSpin, nSpin
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), tmpS(orb%mOrb,orb%mOrb)

    @:ASSERT(size(ham,dim=1)==size(over))
    @:ASSERT(size(nNeighbour)==nAtom)
    @:ASSERT(size(iNeighbour,dim=2)==nAtom)
    @:ASSERT(size(species)>=maxval(iNeighbour))
    @:ASSERT(size(species)<=size(img2CentCell))
    @:ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    @:ASSERT(size(iPair,dim=2)==nAtom)
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(size(ham,dim=2)==size(shift,dim=4))

    nSpin = size(shift,dim=4)

    do iSpin = 1, nSpin
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          tmpS(1:nOrb2,1:nOrb1) = reshape( &
              & over(iOrig+1:iOrig+nOrb2*nOrb1),(/nOrb2,nOrb1/) )
          tmpH(1:nOrb2,1:nOrb1) = 0.5_dp * ( &
              & matmul(tmpS(1:nOrb2,1:nOrb1), &
              & shift(1:nOrb1,1:nOrb1,iAt1,iSpin)) + &
              & matmul(shift(1:nOrb2,1:nOrb2,iAt2f,iSpin), &
              & tmpS(1:nOrb2,1:nOrb1)) )
          ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) = &
              & ham(iOrig+1:iOrig+nOrb2*nOrb1,iSpin) + &
              & reshape(tmpH(1:nOrb2, 1:nOrb1), (/nOrb2*nOrb1/))
        end do
      end do
    end do

  end subroutine add_shift_block


  subroutine addatom_shell(shiftshell, atom, orb, species)
    real(dp), intent(inout) :: shiftshell(:,:,:)
    real(dp), intent(in) :: atom(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)

    integer iAtom, iSpin, nAtom, nSpin

    nAtom = size(atom, dim=1)
    nSpin = size(atom, dim=2)

    @:ASSERT(size(shiftshell, dim=1) == orb%mShell)
    @:ASSERT(size(shiftshell, dim=2) == nAtom)
    @:ASSERT(size(shiftshell, dim=3) == nSpin)
    @:ASSERT(size(species) >= nAtom)

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        shiftshell(1:orb%nShell(species(iAtom)),iAtom,iSpin) = &
            & shiftshell(1:orb%nShell(species(iAtom)),iAtom,iSpin) &
            & + atom(iAtom,iSpin)
      end do
    end do

  end subroutine addatom_shell


  subroutine addshell_block(shiftblock, shell, orb, species)
    real(dp), intent(inout) :: shiftblock(:,:,:,:)
    real(dp), intent(in) :: shell(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)

    integer iAt, iSpin, nAtom, nSpin, iSh, iSp, iOrb

    nAtom = size(shiftblock, dim=3)
    nSpin = size(shiftblock, dim=4)

    @:ASSERT(size(shiftblock, dim=1) == orb%mOrb)
    @:ASSERT(size(shiftblock, dim=2) == orb%mOrb)
    @:ASSERT(size(shell, dim=1) == orb%mShell)
    @:ASSERT(size(shell, dim=2) == nAtom)
    @:ASSERT(size(shell, dim=3) == nSpin)
    @:ASSERT(size(species) >= nAtom)

    do iSpin = 1, nSpin
      do iAt = 1, nAtom
        iSp = species(iAt)
        do iSh = 1, orb%nShell(iSp)
          do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
            shiftblock(iOrb,iOrb,iAt,iSpin) = shiftblock(iOrb,iOrb,iAt,iSpin)&
                & + shell(iSh,iAt,iSpin)
          end do
        end do
      end do
    end do

  end subroutine addshell_block


end module shift
