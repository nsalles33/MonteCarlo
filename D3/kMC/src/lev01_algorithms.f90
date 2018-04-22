! =================================================================================================
!     Kinetic Monte Carlo Kernel v.1.0
! -------------------------------------------------------------------------------------------------
!  Write by N. Salles for Master HPC SISSA/ICTP 2017/18
!  Lecture in collabaration with S. Sorella and R. Innocente
! =================================================================================================
!
!   All Kinetic Monte Carlo Algorithms Available
!   - Residence Time Algorithms ( BKL )
!   - Gillepsie
! =================================================================================================

  subroutine Algorithm( struc )
    use iso_c_binding
    use derived_types
    implicit none
    !
    type( KMC_type ) :: struc
    character (len=3) :: temp1
    character (len=9) :: temp2
    !
    !
    temp1 = trim(struc% algorithm)
    !if ( trim(struc% algorithm) == 'BKL' ) call algo_bkl( struc )
    if ( temp1 == 'BKL' ) call algo_bkl( struc )
    !
    temp2 = trim(struc% algorithm)
    !if ( trim(struc% algorithm) == 'Gillepsie' ) call algo_Gillepsie( struc )
    if ( temp2 == 'Gillepsie' ) call algo_Gillepsie( struc )
    !
  end subroutine Algorithm
! =================================================================================================

  subroutine algo_bkl( struc )
    use derived_types
    use lib_hub
    implicit none

    type( KMC_type ) :: struc
    integer :: isite, ievent

    !print*, " -= BKL Algorithm =- "

    call rate_sum_calc( struc )

    call choose_event_hub( struc, isite, ievent )

    if ( struc% conv == 1 ) return

    call event_applied_hub( struc, isite, ievent )

    call Time_increment( struc )
    !stop " in algo_BKL..."

  end subroutine algo_bkl
! .................................................................................................

    subroutine time_increment( struc )
      use iso_c_binding
      use derived_types
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      real( c_double )                  :: rdn
      !
      call random_number( rdn )
      struc% rand_time = - log( rdn )/struc% sum_rate
      struc% time = struc% time + struc% rand_time
      !
      !stop " Time_increment..."
      !
    end subroutine
! =================================================================================================

  subroutine algo_gillepsie( struc )
    use derived_types
    use lib_hub
    implicit none
    !
    type( KMC_type ) :: struc
    integer :: ievent, isite
    !
    call rate_sum_calc( struc )
    !
    call choose_event_hub( struc, isite, ievent )
    !
    if ( struc% conv == 1 ) return
    !
    call event_applied_hub( struc, isite, ievent )
    !
    !call Time_increment( struc )
    struc% time = struc% time + struc% rand_time
    !
  end subroutine algo_gillepsie
! =================================================================================================
