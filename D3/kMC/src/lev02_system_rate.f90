! =================================================================================================
!     Kinetic Monte Carlo Kernel v.1.0
! -------------------------------------------------------------------------------------------------
!  Write by N. Salles for Master HPC SISSA/ICTP 2017/18
!  Lecture in collabaration with S. Sorella and R. Innocente
! =================================================================================================
!
!
!    Routine to calculate the rate of each event and each site  
!
!
    subroutine rate_sum_calc( struc )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use lib_hub, only : event_rate_calc_hub

      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int )                  :: i, j0, jv
      real( c_double )                  :: sum_rate

      integer( c_int ), dimension(:), pointer :: nneig, neig
      real( c_double ), dimension(:), pointer :: rate, event_rate
      !
      !write (*,*) " EVENT_RATE_CALC_HUB... "
      call event_rate_calc_hub( struc )
      !
      !  -- Sum over all the site
      !
      sum_rate = 0.0
      call link_int1_ptr( struc% ptr_nneig,       nneig,      struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig,        neig,       nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      !
      do i = 1,struc% tot_sites
         sum_rate = sum_rate + rate( i )
      enddo
      !
      !print*, "SUM RATE", sum_rate, sum2
      if ( sum_rate == 0.0 ) call error( " ERROR: sum_rate = 0.0" )
      struc% sum_rate = sum_rate
      !
      !write (*,*) " DONE "
      !
    end subroutine rate_sum_calc
! ..................................................................................................





