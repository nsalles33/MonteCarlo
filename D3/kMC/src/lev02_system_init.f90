! =================================================================================================
!     Kinetic Monte Carlo Kernel v.1.0
! -------------------------------------------------------------------------------------------------
!  Write by N. Salles for Master HPC SISSA/ICTP 2017/18
!  Lecture in collabaration with S. Sorella and R. Innocente
! =================================================================================================
!
! =================================================================================================
!    Routine to initialize the system
!  To do (or think) :
!  - Format file for input_file and event library
!    
! =================================================================================================

    subroutine Init_system( struc )
      use derived_types
      use sub_new_types
      use errors
      !use lib_hub
#ifdef SHARED_LIB 
      use dlopen_lib
#else
      use user_system
#endif
      implicit none
      !
      type( KMC_type ), intent( inout ) :: struc
      !
#ifdef SHARED_LIB
      procedure( read_event ), bind( C ), pointer :: read_event_proc
#endif
      !
      !
      !  ::: Read input_file
      call read_input( struc )
      !print*, "... sortie de read_input -> print_kmc_type"
         call print_kmc_type( struc )
      !
      !  ::: READ EVENT FILE
      !
#ifdef SHARED_LIB 
      !
      call open_shared_lib( struc )
      call c_f_procpointer( proc_read_event, read_event_proc )
      call read_event_proc( struc )
      !
#else
      !
      call read_event( struc )
      !
#endif
      !call read_event_hub( struc )
      !
      !  call print_event( struc% event )
      !stop "system_init..."
      !
      !  ::: Initialization of node State
      !      This step depend on system...
      call distribution_state ( struc ) ! This routine is specific for vacancies diffusion
      !
      call neig_list( struc )
      !
    end subroutine Init_system
! =================================================================================================

    subroutine read_input( struc )
      use derived_types
      use sub_new_types
      use errors
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer                           :: i
      integer                           ::  u0, ios, node_state !, nstep
!      real                              :: temp
      character (len=500)                :: string !, input_event, algo_txt
      character (len=1)                 :: px !,py,pz
      logical                           :: EOF

      character (len=1)  :: delims
      CHARACTER (len=100), dimension (50) :: args
      integer :: nargs
      !
      call init_kmc_type( struc )
      !
      !sys_dim = 2
      !nsites = 10
      !struc% nprop = 0
      !
      !input_file = "input_KMC.dat"
      open( newunit=u0, file=trim(struc% input_file), iostat=ios )
      if ( ios /= 0 ) call error( "input_file does't open!!" )
      !write (*,*) trim(input_file)," open...",u0,ios
      !
      EOF = .false.
      delims = " "
      do while ( .not.EOF )
         !
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         !print*, "Lu: ", trim(args(1))," ", trim(args(2))
         !
         if ( args(1) == "verbose" )          read ( args(2), '(i1)' ) struc% bavard
         if ( args(1) == "system_dimension" ) read ( args(2), '(i5)' ) struc% sys_dim
         if ( args(1) == "Nber_node_1dim".or.  &
              args(1) == "nsite_x" )          read ( args(2), '(i5)' ) struc% nsites(1)
         if ( args(1) == "nsite_y" )          read ( args(2), '(i5)' ) struc% nsites(2)
         if ( args(1) == "nsite_z" )          read ( args(2), '(i5)' ) struc% nsites(3)
         if ( args(1) == "species" )          read ( args(2), '(i5)' ) struc% nspec
         !
         if ( args(1) == "node_prop" ) then
            read ( args(2), '(i5)' ) struc% node_state
            struc% event% nbond = struc% node_state
         endif
         !
         if ( args(1) == "input_event" ) then
            read ( args(2), '(a)'  ) struc% input_event
            struc% input_event = trim(struc% input_event)//c_null_char
         endif
         !
         if ( args(1) == "shared_library" ) then
            read ( args(2), '(a)'  ) struc% libname
            struc% libname = trim( struc% libname )//c_null_char
            print*, args(1), args(2),trim(struc% libname)
         endif
         !
         if ( args(1) == "algorithm" ) then
            read ( args(2), '(a)'  ) struc% algorithm
            struc% algorithm = trim(struc% algorithm)//c_null_char
         endif
         !
         if ( args(1) == "scale_parameter" )  read ( args(2), '(f6.6)' ) struc% scale
         !
         if ( args(1) == "partial_pressure" ) read ( args(2), '(i5)' ) struc% npressure
         if ( args(1) == "temperature" )      read ( args(2), '(f6.6)' ) struc% temp
         if ( args(1) == "kt" )               read ( args(2), '(f6.6)' ) struc% kt
         if ( args(1) == "nstep" )            read ( args(2), '(I10)' ) struc% max_step
         if ( args(1) == "freq_write" )       read ( args(2), '(I6)' ) struc% freq_write
         if ( args(1) == "calc_properties" )  read ( args(2), '(I6)' ) struc% nprop
         !
         if ( args(1) == "init_config" ) then
            read ( args(2), '(a)' ) struc%init_mod    
            !
            if ( trim(struc% init_mod) == "random" )  &
               read ( args(3), '(f6.6)' ) struc% per100
         endif
         !
         if ( args(1) == "bound_condition" )  then
            do i = 1,3
               read ( args(i+1), '(a)' ) px
               if ( trim(px) == "p") then ; struc% period( i ) = 1
               else                       ; struc% period( i ) = 0
               endif
            enddo
         endif
         !
      enddo
      write (*,*) " LECTURE ",struc% sys_dim, struc% nsites(1), node_state, struc% input_event, struc% npressure, struc% temp
      close( u0 )
      !
      !if ( sys_dim == 2 )   &
      !  call builder_kmc_type( sys_dim, nsites, algorithm=algo_txt, temperature=temp, step=nstep )
      !if ( sys_dim == 3 )   &
      !  call builder_kmc_type( sys_dim, nsites, size_y=ny, size_z=nz, algorithm=algo_txt, temperature=temp, step=nstep )
      !
      call builder_kmc_type( struc ) 
      !
      print*, "On sort de read_input..."
      !
    end subroutine read_input
! =================================================================================================

   subroutine distribution_state ( struc )
     use iso_c_binding
     use derived_types
     use sub_new_types
     use random
     use errors
     implicit none 
     !
     type( KMC_type ), intent( inout ) :: struc
     integer( c_int )                  :: i, n
     real( c_double )                  :: rnd
     !
     integer( c_int ), dimension(:), pointer :: site
     call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
     !
     print*, "Distribution_state mode : ", struc% init_mod, struc% per100
     !
     !
     if ( trim(struc% init_mod) == "random" ) then
        !
        if ( struc% per100 == 0.0 ) call warning(" DIST_STATE: per100 is 0.0!!!")
        !
        call set_random_seed()
        !
        n = 0
        do i = 1,struc% tot_sites
           call random_number( rnd ) 
           if ( rnd <= 1 - struc% per100 ) then
              site( i ) = 1
           else
              site( i ) = 0
              n = n + 1
           endif
        enddo
        !
        print*, " Random '1' : ", n
        !
     else if ( trim(struc% init_mod) == "species" ) then
        !
        if ( struc% nspec == 0 ) call error( "DIST_STATE: No species declared!!!")
        !
        do i = 1, struc% tot_sites
           site( i ) = 1
           print*, " site: ", i, " activity :", site( i )
        enddo
        !
     else if ( trim(struc% init_mod) == "none" ) then
        !
        do i = 1, struc% tot_sites
           site( i ) = 0
           !print*, " site: ", i, " activity :", site( i )
        enddo
        !
     else if ( trim(struc% init_mod) == "one" ) then
        !
        n = struc% nsites(1)/2
        if ( struc% sys_dim == 1 ) n = 0
        do i = 1, struc% tot_sites
           site( i ) = 1
           if ( i == (struc% tot_sites/2 + n) )  site( i ) = 0
        enddo
        !
     endif
     !
     !stop " Distrib...ifin"
     !
   end subroutine distribution_state
! =================================================================================================


















