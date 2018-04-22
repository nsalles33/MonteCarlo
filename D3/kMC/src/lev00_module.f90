! =================================================================================================
!     Kinetic Monte Carlo Kernel v.1.0
! -------------------------------------------------------------------------------------------------
!  Write by N. Salles for Master HPC SISSA/ICTP 2017/18
!  Lecture in collabaration with S. Sorella and R. Innocente
! =================================================================================================
!
! =================================================================================================
  module errors
    use iso_c_binding
    implicit none
  contains
! ............................................................................
!      subroutine error( text ) bind( C )
!        character ( len=*, kind=c_char ), intent( in ) :: text
!        print*, " ERROR: ",trim(text)
!        stop " ===== FATAL ERROR!!"
!      end subroutine error
! ............................................................................
      subroutine error( text ) bind( C )
        character ( kind=c_char ), dimension(*), intent( in ) :: text
        integer( c_int ) :: len, i
        len = 0
        do
           if ( text(len+1) == C_NULL_CHAR ) exit
           len = len + 1
        enddo
        print*, " ERROR: ",(text(i),i=1,len)
        stop " ===== FATAL ERROR!!"
      end subroutine error
! ............................................................................
!      subroutine warning( text ) bind( C )
!        character ( len=*, kind=c_char ), intent( in ) :: text
!        print*, " WARNING: ",trim(text)
!      end subroutine warning
! ............................................................................
      subroutine warning( text ) bind( C )
        character ( kind=c_char ), dimension(*), intent( in ) :: text
        integer( c_int ) :: len, i
        len = 0
        do
           if ( text(len+1) == C_NULL_CHAR ) exit
           len = len + 1
        enddo
        print*, ' WARNING: ',(text(i),i=1,len)
      end subroutine warning
! ............................................................................
  end module errors
! =================================================================================================
  module random
    use iso_c_binding
    implicit none
  contains
! ............................................................................
    subroutine set_random_seed() bind( C )
     ! ----- setting a random seed, based on current time -----
     INTEGER( c_int ) :: i_seed
     INTEGER( c_int ), DIMENSION(:), ALLOCATABLE :: a_seed
     INTEGER( c_int ), DIMENSION(1:8) :: dt_seed

     ! ----- Set up random seed portably -----
     CALL RANDOM_SEED( size = i_seed )
     ALLOCATE( a_seed( 1: i_seed ) )
     CALL RANDOM_SEED( get = a_seed )
     CALL DATE_AND_TIME( values = dt_seed )
     a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
     CALL RANDOM_SEED( put = a_seed )
     DEALLOCATE( a_seed )
    end subroutine set_random_seed
! ............................................................................

  end module random
! =================================================================================================

  module sub_new_types
    use iso_c_binding
    use derived_types
    use errors
    implicit none

  contains

  !end module sub_new_types

    ! ............................................................................
! ............................ KMC_T BUILDER .................................
! ............................................................................

    subroutine init_kmc_type( this ) bind( C )
      implicit none
      type( KMC_type ) :: this

      this% bavard = 0
      this% conv = 0
      this% period(:) = 0

      this% algorithm = "BKL"
      this% input_event = "event_lib"
      this% libname = "$.so"

      this% time = 0
      this% max_time = 1.0
      this% rand_time = 0.0

      this% freq_write = 1
      this% max_step = 1
      this% sum_rate = 0.0
      this% rand_rate = 0.0
      this% nspec = 0
      this% npressure = 0.0
      this% scale = 1.0

      this% temp = 0.0
      this% kt = 0.0
      this% f0 = 0.0

      this% tot_sites =  1
      this% sys_dim = 0
      this% per100 = 0.0

      this% nsites(:) = 1

    end subroutine init_kmc_type
! ............................................................................
!
    subroutine builder_kmc_type( this ) bind( C )
      use errors
      implicit none
      type( KMC_type ) :: this

      integer( c_int ) :: n
      real( c_double ) :: mem
      print*, " Enter in KMC_type allocation..."

      if ( this% temp == 0.0.and.this% kt == 0.0 ) then
         this% kt = 1.0 
      else if ( this% temp == 0.0.and.this% kt /= 0.0 ) then
         this% temp = this% kt/kb
      else if ( this% temp /= 0.0.and.this% kt == 0.0 ) then
         this% kt = kb*this% temp
         print*, " Kb.T =", this% kt
      endif

      this% f0 = 1.0 ;! 1e-12
      !
      !  ::: System Size
      !print*, 'System Size...'
      if ( this% sys_dim > 3 .or. this% sys_dim == 0 ) &
         call error( " BAD SYSTEM DIMENSION " )
      !
      if ( this% nsites(1) == 0 ) &
         call error( " Lx is not declared " )
      !
      if ( this% sys_dim == 2.and.this% nsites(2) == 0.and.this% nsites(1) /= 0 ) then
         this% nsites(2) = this% nsites(1)
         call warning( " System 2D square Lx = Ly " )
      endif
      !
      if ( this% sys_dim == 3.and.this% nsites(2) == 0.and.this% nsites(3) == 0 )  &
        call error( " System 3D => Ly and Lz must be declared " )
      !
      this% tot_sites = this% nsites(1)*this% nsites(2)*this% nsites(3)
      n = this% tot_sites
      !
      !  ::: Table allocation 
      !
      allocate( h_site(n),              &
                h_rate(n),              &
                h_nneig(n),             &
                h_neig(nvois*n),        &
                h_event_rate(nvois*n),  &
                h_nevt(n),              &
                h_event_site(nvois*n) )
      !
      if ( .not.allocated(h_site)  .or. .not.allocated(h_rate) .or.       &
           .not.allocated(h_nneig) .or. .not.allocated(h_neig) .or.       &
           .not.allocated(h_event_rate).or. .not.allocated(h_nevt) .or.   &
           .not.allocated(h_event_site) )   &
         call error( " KMC_type => CONSTRUCTOR problem..." )
      !
      !
      if ( this% nspec /= 0 ) then
          !
          allocate( h_spec( this% nspec ) )
          !
          if ( .not.allocated( h_spec ) ) &
             call error( " KMC_type => CONSTRUCTOR h_spec problem..." )
          !
      endif
      !
      !
      if ( this% npressure /= 0 ) then
         !
         allocate( h_ppressure( this% npressure ) )
         allocate( h_masse( this% npressure ) )
         !
         if ( .not.allocated( h_ppressure ).or. .not.allocated( h_masse) ) &
            call error( " KMC_type => CONSTRUCTOR h_ppressure problem..." )
         !
      endif
      !
      mem = ( n + nvois*n )*sizeof( mem ) + ( 3*n + nvois*n*2 )*sizeof(n)
      print*, " MEMORY ALLOCATED :: ",mem,"Bytes"
      !
      !  ::: Connection between ptr_table -> h_table
      this% ptr_site       = c_loc( h_site(1) )
      this% ptr_rate       = c_loc( h_rate(1) )
      this% ptr_neig       = c_loc( h_neig(1) )
      this% ptr_nneig      = c_loc( h_nneig(1) )
      this% ptr_event_rate = c_loc( h_event_rate(1) )
      this% ptr_nevt       = c_loc( h_nevt(1) )
      this% ptr_event_site = c_loc( h_event_site(1) )
      if ( allocated( h_spec ) ) &
         this% ptr_spec    = c_loc( h_spec(1) )
      if ( allocated( h_ppressure ) ) then
         this% ptr_pressure = c_loc( h_ppressure(1) ) 
         this% ptr_masse = c_loc( h_masse(1) )
      endif
      !
      !  ::: Properties table
      !print*, 'Prop Table allocation...'
      if ( this% nprop /= 0 ) allocate( h_prop( this% nprop ) ) !, this% txtprop(this% nprop) )
      !
      this% ptr_prop = c_loc( h_prop(1) )
      !
      print*, " KMC_type Constructor DONE"
      !
    end subroutine builder_kmc_type
!!!!! ............................................................................
    !
    subroutine destructor_kmc_type()  bind( C )
      implicit none
      !type( KMC_type ) :: this
      !
      !  :: destroy event_type
      call destructor_event_type !( this% event ) 
      !
      deallocate( h_site, h_rate, h_nneig, h_neig, h_event_rate,   &
                  h_nevt, h_event_site )
      !
      if ( allocated(h_site).or.allocated(h_rate).or.       &
           allocated(h_nneig).or.allocated(h_neig).or.      &
           allocated(h_event_rate).or.allocated(h_nevt).or. &
           allocated(h_event_site).or.allocated(h_ebond) )  &
         print*, " KMC_type => DESTRUCTOR problem ..."
      !
      if ( allocated(h_ppressure) ) deallocate( h_ppressure )
      if ( allocated(h_masse) )    deallocate( h_masse )
      !
    end subroutine destructor_kmc_type
!!!!! ............................................................................
    !
    subroutine Init_table( this ) bind( C )
      implicit none
      type( KMC_type ) :: this
      !integer( c_int ) :: i, jn, nv

      integer( c_int ), dimension(:), pointer :: nneig
      real( c_double ), dimension(:), pointer :: rate, prop, event_rate

      print*, " INIT_TABLE..."

      call link_int1_ptr( this% ptr_nneig, nneig, this% tot_sites )
      call link_real1_ptr( this% ptr_rate, rate, this% tot_sites )
      call link_real1_ptr( this% ptr_prop, prop, this% nprop )
      call link_real1_ptr( this% ptr_event_rate, event_rate, 5*this% tot_sites )
!  ::: Rate Initialization at 0.
      rate = 0.0
      event_rate = 0.0
      prop = 0.0
!  ::: Properties Table initialization 
      print*, " INIT_TABLE... DONE"
    end subroutine Init_table
! ............................................................................
!
    subroutine builder_event_type( this, nevt ) bind( C )
      implicit none
      type( event_type )   :: this
      integer( c_int ), intent( in ) :: nevt
      integer( c_int ) :: n

      print*, " Enter in EVENT_type constructor..."
      this% nevent = nevt
      n = this% nevent 

      allocate( h_i_state(n), h_f_state(n), h_f0(n), h_Ebarrier(n), h_dE(n), h_ebond(this% nbond, this% nbond) )
      if ( .not.allocated(h_i_state) .or. .not.allocated(h_f_state) .or.   &
           .not.allocated(h_ebarrier) .or. .not.allocated(h_de).or.        &
           .not.allocated(h_ebond).or. .not.allocated( h_f0 ) )          &
         call error( " EVENT_type => CONSTRUCTOR problem " )

      this% ptr_i_state  = c_loc( h_i_state(1)  )
      this% ptr_f_state  = c_loc( h_f_state(1)  )
      this% ptr_f0       = c_loc( h_f0(1)       )
      this% ptr_ebarrier = c_loc( h_ebarrier(1) )
      this% ptr_de       = c_loc( h_de(1)       )
      this% ptr_ebond    = c_loc( h_ebond(1,1)  )

      print*, " EVENT_type constructor DONE "

    end subroutine builder_event_type
! ............................................................................
!
    subroutine destructor_event_type() bind( C )
      implicit none
      !type( event_type ) :: this

      !nullify( this% ptr_i_state, this% ptr_f_state, this% ptr_ebarrier, this% ptr_de )
      if ( allocated(h_i_state) ) deallocate( h_i_state )
      if ( allocated(h_f_state) ) deallocate( h_f_state )
      if ( allocated(h_ebarrier) ) deallocate( h_ebarrier )
      if ( allocated(h_de) ) deallocate( h_de )
      if ( allocated(h_ebond) ) deallocate( h_ebond )
      if ( allocated(h_f0) ) deallocate( h_f0 )
      !
      if ( allocated(h_i_state).or.allocated(h_f_state).or.  &
           allocated(h_ebarrier).or.allocated(h_de) )        &
         call error( "EVENT_type => DESTRUCTOR problem " )

    end subroutine destructor_event_type
! ............................................................................
!
    subroutine print_kmc_type ( this ) bind( C )
      implicit none
      type( KMC_type ) :: this
      integer( c_int ) :: i
      !
      !print*, " - To do :: print_kmc_type "
      write (*,*) " ======== SYSTEM PARAMETERS ======= "
      write (*,*) " Algorithm         : ", this% algorithm
      write (*,*) " MAX KMC steps     : ", this% max_step
      write (*,*) " MAX KMC time (s?) : ", this% max_time
      write (*,*) " System Dimension  : ", this% sys_dim
      do i = 1,this%sys_dim
         write (*,*) " Nber of Node : ", this% nsites(i)
      enddo
      write (*,*) " System Size (nodes) : ", this% tot_sites
      write (*,*) " Number of Species   : ", this% nspec
      write (*,*) " Freq Write          : ", this% freq_write
      write (*,*) " Bound Condition     : ", ( this% period( i ), i=1,3 )
      write (*,*) " Temperature  (Kb.T) : ", this% temp, this% kt
      write (*,*) " Partial Pressure    : ", this% npressure
      write (*,*) " ================================== "

    end subroutine print_kmc_type
! ............................................................................

   subroutine print_event( this )  bind( C )
      implicit none
      type( event_type ) :: this
      integer( c_int ) :: i, id!, nevt2

      integer( c_int ), dimension(:), pointer :: init_state, final_state
      real( c_double ), dimension(:), pointer :: ebarrier, de, f0
      call link_int1_ptr( this% ptr_i_state,   init_state,  this% nevent)
      call link_int1_ptr( this% ptr_f_state,   final_state, this% nevent)
      call link_real1_ptr( this% ptr_f0,       f0,          this% nevent)
      call link_real1_ptr( this% ptr_ebarrier, ebarrier,    this% nevent)
      call link_real1_ptr( this% ptr_de,       de,          this% nevent)
      !
      !nevt2 = this% nevent/2
      !if ( nevt2 == 0 ) nevt2 = this% nevent
      !
      write (*,*) " =========== EVENT LIB ========= ", this% nevent
      do i = 1,this% nevent
         id = i
         write (*,*) " Evt:",id,init_state( id ), final_state( id ), &
                                f0( id ),ebarrier( id )!, de( id )
      if ( id == 2 ) stop "print_event"
         !if ( this% nevent > 1 ) then
         !   id = i + nevt2
         !   write (*,*) " Inv:",id,init_state( id ), final_state( id ), &
         !                       ebarrier( id ), de( id )
         !endif
      enddo
      !
    end subroutine print_event
! ............................................................................
!
    subroutine print_state( this, step, u0 ) bind( C )
      implicit none
      type( KMC_type ) :: this
      integer( c_int ), intent ( in ) :: step, u0
      integer( c_int ) :: ios, i, x, y, z, nx, ny, nxy, f
      real( c_double ) :: t_start, t_stop
      !
      integer( c_int ), dimension(:), pointer :: site
      real( c_double ), dimension(:), pointer :: prop
      call link_int1_ptr( this% ptr_site, site, this% tot_sites )
      call link_real1_ptr( this% ptr_prop, prop, this% nprop )
      !
      if ( MODULO( step, this% freq_write ) /= 0 ) return
      !write (*,*) " ...PRINT_STATE... "
      !
      nx = this% nsites(1)
      ny = this% nsites(2)
      nxy = nx*ny
      !
      !  ::: Write Configuration
      !    -- 3D cubic
      write (u0,fmt='(1x,I6)',iostat=ios) this% tot_sites
        if (ios /=0) write (*,*) " Problem write => state_file.xyz  "
      !
      write (u0,*,iostat=ios) step, this% sys_dim,"DIM", (this% nsites(i),i=1,this% sys_dim)
        if (ios /=0) write (*,*) " Problem write => state_file.xyz  "
      !
      call cpu_time( t_start )
      !
      do i = 0,this%tot_sites - 1
         x = MODULO( i, nx )
         y = MODULO( i/nx, ny )
         z = i / nxy
         write (u0,'(1x,I2,3(2x,I4))') site(i+1), x, y, z
      enddo
      !
      call cpu_time( t_stop )
      !write (*,'(1x,a,1x,E12.3,1x,a)') " Write Conf. time", t_stop - t_start, "second "
      !
      !  ::: Write Energetic Statistic
      if ( step == 0 )  then
        write (*,*) "    Step   |    Time   | rand_rate |    rand_time    " !,( " | ",trim(this% txtprop(i)),i=1,this%nprop )
        write (*,*) " ----------------------------------------------------------------- "
      endif
      f = this% nprop
      write (*,'(1x,i10,*(2x,E10.4))') step, this% time, this% rand_rate, this% rand_time, (prop(i), i=1, this% nprop)
      write (100,'(1x,i10,*(2x,E10.4))') step, this% time, this% rand_rate, this% rand_time, (prop(i), i=1, this% nprop)
      !
    end subroutine print_state
! ............................................................................
  
    subroutine link_int1_ptr( ptr_c, ptr_f, n ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n
      integer( c_int ), dimension(:), pointer :: ptr_f
      !print*, "link_int1_ptr..."
      call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
    end subroutine link_int1_ptr
! ............................................................................

    subroutine link_int2_ptr( ptr_c, ptr_f, n, m ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n, m
      integer( c_int ), dimension(:,:), pointer :: ptr_f
      !print*, "link_int2_ptr..."
      call c_f_pointer( ptr_c, ptr_f, shape=[n,m] )
    end subroutine link_int2_ptr
! ............................................................................

    subroutine link_real1_ptr( ptr_c, ptr_f, n ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n
      real( c_double ), dimension(:), pointer :: ptr_f
      !print*, "link_real1_ptr..."
      call c_f_pointer( ptr_c, ptr_f, shape=[ n ] )
    end subroutine link_real1_ptr
! ............................................................................

    subroutine link_real2_ptr( ptr_c, ptr_f, n, m ) bind( C )
      implicit none
      type( c_ptr ) :: ptr_c
      integer( c_int ), intent( in ) :: n,m
      real( c_double ), dimension(:,:), pointer :: ptr_f
      !print*, "link_real2_ptr..."
      call c_f_pointer( ptr_c, ptr_f, shape=[ n,m ] )
    end subroutine link_real2_ptr
! ............................................................................


  end module sub_new_types
! =================================================================================================
