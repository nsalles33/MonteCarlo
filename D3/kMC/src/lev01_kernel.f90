! =================================================================================================
!     Kinetic Monte Carlo Kernel v.1.0
! -------------------------------------------------------------------------------------------------
!  Write by N. Salles for Master HPC SISSA/ICTP 2017/18
!  Lecture in collabaration with S. Sorella and R. Innocente
! =================================================================================================
!
Program Kernel
  !
  use iso_c_binding
  use derived_types
  use sub_new_types
  use errors
  use lib_hub  

#ifdef SHARED_LIB
  use dlopen_lib
#endif

  implicit none
  !
  integer          :: nstep, u0, ios
  real             :: t_start, t_stop
  !
  type( KMC_type ) :: sys 
  !
  ! ................................................
  !
  ! ---------------- Input_file Recuperation 
  if ( iargc() >= 1 ) then ; call getarg( 1, sys% input_file )
  else ; call error( " Execution : ./EXE.x input_file " ); 
  endif    
  sys% input_file= trim(sys% input_file)//c_null_char
  write (*,*) "Input_file : ", sys% input_file
  !
  !
  ! ---------------- SYSTEM INITIALIZATION 
  !
  call Init_system( sys )
  !
  call init_table( sys )
  !
  !
  ! ---------------- OPEN CONF FILE - - - - - - - - - - - - - - - - - - - - 
  !
  open( newunit=u0, file="state_file.xyz",status="replace", iostat=ios )
  if (ios /= 0) call warning( " CAN'T OPEN => state_file.xyz " )
  !
  ! ---------------- - - - - - - - - - - - - - - - - - - - - - - - - - - - 


  
  ! ---------------- SYSTEM EVOLUTION
  !
  call cpu_time( t_start )
  !call system_clock( t_start )
  !
  !
  nstep = 0
  sys% step = nstep
  if ( MODULO( nstep, sys% freq_write ) == 0 ) &
     call analysis( sys )
  !
  !
  call print_state( sys, nstep, u0 )
  !
  !
  do while ( nstep < sys% max_step )
     !
     !
     nstep = nstep + 1
     sys% step = nstep
     call algorithm( sys )
     !
     !
     if ( sys% conv == 1 ) exit 
     !stop " Kernel..."
     !
     !
     if ( MODULO( nstep, sys% freq_write ) == 0 ) &
        call analysis( sys )
     !
     !
     call print_state( sys, nstep, u0 )
     !
     !
  enddo
  !
  !
  call cpu_time( t_stop )
  !call system_clock( t_stop )
  !
  !
  if ( sys% conv == 1 ) &
     write (*,*) " Finit with Convergence criturium "
  !
  if ( nstep == sys% max_step ) &
     write (*,*) " Finit with max_step criturium "
  !
  write (*,'(a,1x,f15.6,1x,a)') " TIME elapse :", t_stop - t_start, " second "
  !
  close( u0 )
  !
  !
  ! ---------------- FINILIZATION 
#ifdef SHARED_LIB
    call close_shared_lib
#endif
  !
  call destructor_kmc_type
  !
end program Kernel

! =================================================================================================





