  module variables
    use iso_c_binding
    implicit none
    !
    integer( c_int ), dimension(:), pointer :: site, nneig, neig, nevt, spec, event_site
    !
    real( c_double ), dimension(:), pointer :: rate, prop, event_rate, pressure, masse
    !
    ! ::: EVENT
    integer( c_int ), dimension(:), pointer :: init_state, final_state
    real( c_double ), dimension(:), pointer :: ebarrier, de, f0
    real( c_double ), dimension(:,:), pointer :: ebond
    !
    integer( c_int ) :: save_site
    !
  end module variables

! .....	.........................................................................................

    subroutine read_event( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: obj
      character (len=500,kind=c_char)   :: string
      integer( c_int )                  :: u0, i, id, ios, nevent, npress
      logical                           :: EOF
      !
      character (len=1, kind=c_char)  :: delims
      CHARACTER (len=100, kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs
      !
      !  ::: Init save_site
      !
      save_site = 0
      !
      !  ::: Lecture of state and rate of each node
      !
      open( newunit=u0, file=trim(obj% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )
      !
      EOF = .false.
      do while ( .not.EOF )
         !
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         !write (*,*) "Lu : ", nargs, (i,trim(args(i)), i=1,nargs)
         !
         ! ::: Lecture Of Event
         !
         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( obj% event, nevent )
            !exit
            !
            call link_int1_ptr( obj% event% ptr_i_state,   init_state,  obj% event% nevent)
            call link_int1_ptr( obj% event% ptr_f_state,   final_state, obj% event% nevent)
            call link_real1_ptr( obj% event% ptr_f0,       f0,          obj% event% nevent )
            call link_real1_ptr( obj% event% ptr_ebarrier, ebarrier,    obj% event% nevent)
            call link_real1_ptr( obj% event% ptr_de,       de,          obj% event% nevent)
            !
            do i = 1,obj% event% nevent
               !
               read (u0,*) id, init_state(id), final_state(id), f0(id), ebarrier(id), de(id)
               write (*,*) id, init_state(id), final_state(id), f0(id), ebarrier(id), de(id)
               !
            enddo
            !
         endif
         !
         ! ::: Lecture Of Partial Pressure Parameters
         !
         if ( args(1) == "partial_pressure" ) then
            !
            read( args(2), '(i5)' ) npress
            if ( npress /= obj% npressure ) call warning( "WARNING npressure /= partial_pressure" )
            !
            call link_real1_ptr( obj% ptr_pressure, pressure, obj% npressure)
            call link_real1_ptr( obj% ptr_masse, masse, obj% npressure)
            !
            do i = 1,npress
               read( u0, * ) id, pressure( id ), masse( id ) 
               print*, "Pressure", id, pressure( id ), masse( id )
            enddo
            !
         endif
         !
      enddo
      !
      !
      close( u0 )
      !stop " read_event..." 
      !
    end subroutine read_event
! .................................................................................................

    subroutine event_rate_calc( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: obj
      integer :: i
      !
      call link_int1_ptr( obj% ptr_site,             site,        obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig,            nneig,       obj% tot_sites )
      call link_int1_ptr( obj% ptr_nevt,             nevt,        obj% tot_sites )
      call link_int1_ptr( obj% event% ptr_i_state,   init_state,  obj% event% nevent )
      call link_int1_ptr( obj% event% ptr_f_state,   final_state, obj% event% nevent )
      call link_real1_ptr( obj% event% ptr_f0,       f0,          obj% event% nevent )
      call link_real1_ptr( obj% event% ptr_ebarrier, ebarrier,    obj% event% nevent )
      call link_real1_ptr( obj% event% ptr_de,       de,          obj% event% nevent )
      call link_real1_ptr( obj% ptr_rate,            rate,        obj% tot_sites )
      !
      call link_int1_ptr( obj% ptr_neig,           neig,        nvois*obj% tot_sites )
      call link_int1_ptr( obj% ptr_event_site,     event_site,  nvois*obj% tot_sites )
      call link_real1_ptr( obj% ptr_event_rate,    event_rate,  nvois*obj% tot_sites )

      call link_real1_ptr( obj% ptr_masse, masse, obj%npressure )

      !
      ! ::: We have 1 event by site : 
      !    Adsorption or Desorption
      !
      do i = 1,obj% tot_sites

      enddo
      !
      !
    end subroutine event_rate_calc
! .................................................................................................

    subroutine choose_event( struc, isite, ievent ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use random
      use variables
      implicit none
      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent
      !
      integer( c_int ) :: i, jn, j0
      real( c_double ) :: rsum, rdn, rrdn
      !
      call link_int1_ptr( struc% ptr_nevt,       nevt,       struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig,     nneig,       struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      call link_int1_ptr( struc% ptr_event_site,   event_site, nvois*struc% tot_sites )
      !
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
      !      
      rsum = 0.0
      do i = 1,struc% tot_sites
         !
         j0 = nneig( i )
         do jn = 1,nevt( i )
         !  !
            ievent = jn
            rsum = rsum + event_rate( j0 + jn )
            if ( rsum > rrdn ) exit
         !  !
         enddo
         !
         isite = i
         if ( rsum > rrdn ) exit
      enddo
      if ( rsum <= rrdn ) write (*,*) " PB choose event...", rsum,struc% sum_rate
      !
      save_site = isite
      !
    end subroutine choose_event
! .................................................................................................

    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none
      type( KMC_type )               :: struc
      integer( c_int ), intent( in ) :: is,   &  ! Site selected
                                        jn       ! "ievent" of site "is" selected 
      integer( c_int )               :: jevt, j0

      call link_int1_ptr( struc% ptr_site,             site,             struc% tot_sites )
      call link_int1_ptr( struc% event% ptr_i_state,   init_state,       struc% tot_sites )
      call link_int1_ptr( struc% event% ptr_f_state,   final_state,      struc% tot_sites )
      call link_int1_ptr( struc% ptr_event_site,       event_site,       nvois*struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig,            nneig,            struc% tot_sites )
      !
      j0 = nneig( is ) 
      jevt = event_site( j0 + jn )
      !
      if ( site( is ) /= init_state( jevt ) ) then
         call warning( " Problem site state not correspond to initial state event ")
      endif
      !
      site( is ) = final_state( jevt )
      !
    end subroutine event_applied
! .................................................................................................

    subroutine analyse( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use variables
      implicit none
      !
      type( KMC_type ), intent( inout ) :: obj
      !
      !
    end subroutine analyse
! .................................................................................................
















