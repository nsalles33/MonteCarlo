
  module variables
    use iso_c_binding
    implicit none
    !
    integer( c_int ), dimension(:), pointer :: site, nneig, neig, nevt, spec, event_site
    !
    real( c_double ), dimension(:), pointer :: rate, prop, event_rate
    !
    ! ::: EVENT
    integer( c_int ), dimension(:), pointer :: init_state, final_state
    real( c_double ), dimension(:), pointer :: ebarrier, de, f0
    real( c_double ), dimension(:,:), pointer :: ebond
    !
    ! ::: TIME
    real( c_double ) :: t_choose
    !
    ! ::: Local update
    integer( c_int ) :: save_site, save_gp(25)
    !
    ! ::: Brownian Mvt
    integer( c_int ) :: navg
    integer( c_int ), dimension(3) :: rold, pos
    !
  end module variables

!
!  USER routine:
!  - subroutine read_event( obj ) 
!  - subroutine event_rate_calc( obj )
!  - subroutine 
!
! .................................................................................................

    subroutine read_event( struc ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none
      !
      type( KMC_type ), intent( inout ) :: struc
      character (len=500,kind=c_char)    :: string
      integer( c_int )                  :: u0, i, id, ibd, jbd, ios, nevent, nbond
      logical                           :: EOF
      !
      character (len=1,kind=c_char)  :: delims = " "
      CHARACTER (len=100,kind=c_char), dimension (50) :: args
      integer( c_int ) :: nargs
      !
      !  ::: INIT save site
      !
      save_site = 0 
      rold = 0
      navg = 0
      t_choose = 0.0
      !
      !  ::: Lecture of state and rate of each node
      !
      print*, " INPUT_EVENT: ",trim(struc% input_event)
      open( newunit=u0, file=trim(struc% input_event), iostat=ios )
      if ( ios /= 0 ) call error( "input_event does't open!!" )
      !
      EOF = .false.
      do while ( .not.EOF )
         !
         call read_line( u0, string, EOF )
         call parse( trim(string), delims, args, nargs )
         !
         write (*,*) "Lu :", nargs, (i,trim(args(i)), i=1,nargs)
         !
         if ( args(1) == "Number_of_event" ) then
            read( args(2), '(i5)' ) nevent
            call builder_event_type( struc% event, nevent )
            !
            call link_int1_ptr( struc% event% ptr_i_state, init_state, struc% event% nevent)
            call link_int1_ptr( struc% event% ptr_f_state, final_state, struc% event% nevent)
            call link_real1_ptr( struc% event% ptr_f0, f0, struc% event% nevent)
            call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier, struc% event% nevent)
            call link_real1_ptr( struc% event% ptr_de, de, struc% event% nevent)
            !
            do i = 1,struc% event% nevent
               read (u0,*) id, init_state(id), final_state(id), f0(id), ebarrier(id), de(id)
               write (*,*) id, init_state(id), final_state(id), f0(id), ebarrier(id), de(id)
            enddo
            !
         endif
         !
         !
         if ( args(1) == "Energy_Bond" ) then
            !
            read( args(2), '(i5)' ) nbond
            !
            if ( nbond /= 0 ) then
               !
               call link_real2_ptr( struc% event% ptr_ebond, ebond, struc% event% nbond, struc% event% nbond)
               !
               do i = 1,nbond
                  call read_line( u0, string, EOF )
                  call parse( trim(string), delims, args, nargs )
                  read( args(2), '(i5)') ibd
                  read( args(3), '(i5)') jbd
                  read( args(4), '(f10.5)') ebond(ibd,jbd)
                  write (*,*) i, ibd, jbd, ebond( ibd, jbd )
               enddo
            endif
            !
         endif
         !
      enddo
      !
      !
      close( u0 )
      !
    end subroutine read_event
! .................................................................................................

    subroutine event_rate_calc( struc ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use errors
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ) ::  is 

      call link_int1_ptr( struc% ptr_site,             site,       struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig,            nneig,      struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig,             neig,       nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate,            rate,       struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate,      event_rate, nvois*struc% tot_sites )
      !
      call link_real1_ptr( struc% event% ptr_f0,       f0,         struc% event% nevent )
      call link_real1_ptr( struc% event% ptr_ebarrier, ebarrier,   struc% event% nevent )
      call link_real1_ptr( struc% event% ptr_de,       de,         struc% event% nevent )
      call link_real2_ptr( struc% event% ptr_ebond,    ebond,      struc% event% nbond, struc% event% nbond )
      !
      do is = 1,obj% tot_sites   

      enddo
      !
      !
      !
    end subroutine event_rate_calc
! ..................................................................................................

    subroutine choose_event( struc, isite, ievent ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use random
      use variables
      implicit none

      type( KMC_type ), intent( inout ) :: struc
      integer( c_int ), intent( inout ) :: isite, ievent

      integer( c_int ) :: i, j0, jn
      real( c_double ) :: rsum, rdn, rrdn, t_start, t_stop

      call link_int1_ptr( struc% ptr_site,        site,       struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig,       nneig,      struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig,        neig,       nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate,       rate,       struc% tot_sites )
      !
      !
      isite = 0 ; ievent = 0
      call random_number( rdn )
      rrdn = rdn*struc% sum_rate
      struc% rand_rate = rrdn
      !
      call cpu_time( t_start )
      !      
      rsum = 0.0
      do i = 1,struc% tot_sites
         !
         j0 = nneig( i )
         do jn = 1,neig( j0 )
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
      !
      if ( rsum <= rrdn ) then
         write (*,*) " PB choose event...", rsum, rrdn, struc% sum_rate
         stop " EXIT..."
      endif
      !
      call cpu_time( t_stop )
      !
      t_choose = t_choose + t_stop - t_start
      !
    end subroutine choose_event
! ..................................................................................................
    
    subroutine event_applied( struc, is, jn ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use variables
      implicit none
      !
      type( KMC_type )               :: struc
      integer( c_int ), intent( in ) :: is, jn
      integer( c_int )               :: j, jv, j0, k, k0, kv,l,  &
                                        igp, see, nx, ny, nz, nxy
      !
      integer( c_int ), dimension(3) :: x
      !
      !
      call link_int1_ptr( struc% ptr_site, site, struc% tot_sites )
      call link_int1_ptr( struc% ptr_nneig, nneig, struc% tot_sites )
      call link_int1_ptr( struc% ptr_neig, neig, nvois*struc% tot_sites )
      call link_real1_ptr( struc% ptr_rate, rate, struc% tot_sites )
      call link_real1_ptr( struc% ptr_event_rate, event_rate, nvois*struc% tot_sites )
      !
      !  ::: Flip 0 -> 1
      !
      j0 = nneig( is )
      if ( site( is ) == 1 ) then
         call warning( " site Flip Problem 0 = 1 ")
         write (*,*) is, site( is ), ( neig(j0+j), site( neig(j0+j) ),j=1,neig(j0) )
         write (*,*) is, rate( is ), ( neig(j0+j), event_rate( j0+j ),j=1,neig(j0) )
      endif
      !
      site( is ) = 1
      !
      ! ::: Empty the table
      rate( is ) = 0.0
      do jv = 1,neig( j0 )
         event_rate( j0 + jv ) = 0.0   
      enddo
      !
      j = neig( j0 + jn )
      if ( site( j ) == 0 ) call warning( " site Flip Problem 1 = 0 ")
      !
      site( j ) = 0
      !
      !
    end subroutine event_applied    
! ..................................................................................................
!
    subroutine analyse( obj ) bind( C )
      use iso_c_binding
      use derived_types
      use sub_new_types
      use variables
      implicit none
      !
      type( KMC_type ) :: obj
      integer( c_int ) :: i, j, j0, jv, k, k0, kv,   &
                          ngp, gpv, nc, mixgp, nvac, &
                          maxsize, max_gp, nx, ny, nz, nxy
      real( c_double ) :: x, y, z, t_start, t_stop
      !
      integer( c_int ), dimension( obj% tot_sites ) :: gp, histo
      integer( c_int ), dimension( -1:int(obj% tot_sites/2) ) :: clster
      !
      call link_int1_ptr( obj% ptr_site, site, obj% tot_sites )
      call link_int1_ptr( obj% ptr_nneig, nneig, obj% tot_sites )
      call link_int1_ptr( obj% ptr_neig, neig, nvois*obj% tot_sites )
      call link_real1_ptr( obj% ptr_prop, prop, obj% nprop )
      !
      !
      call cpu_time( t_start )
      !
      ! ::: Cluster size :::
      !
      ngp = 0 ; gp = 0 ; max_gp = 0
      clster = 0 ; nvac = 0
      do i = 1,obj% tot_sites
         !
         if ( site( i ) == 1 ) cycle
         !
         nvac = nvac + 1
         nc = 0
         gp( i ) = -1
         gpv = 0 ; mixgp = 0
         j0 = nneig( i )
         do jv = 1,neig( j0 )
            j = neig( j0 + jv )
            !
            if ( site( j ) == 1 ) cycle
            !
            if ( gp( j ) /= 0.and.gpv == 0 ) then
               gpv = gp( j )
            elseif ( gp( j ) /= 0.and.gpv /= 0.and.gp(j) /= gpv ) then
               mixgp = gp( j )
            endif
            !
            nc = nc + 1
            !
         enddo
         !
         !write (*,*) i,nc,ngp,gpv,mixgp
         !
         if ( nc /= 0.and.gpv /= 0.and.mixgp == 0 ) then
            gp( i ) = gpv
         elseif ( nc /= 0.and.gpv == 0 ) then
            ngp = ngp + 1
            gp( i ) = ngp
            do jv = 1,neig( j0 )
               j = neig( j0 + jv )
               if ( site( j ) == 0 ) gp( j ) = gp( i )
            enddo
         elseif ( mixgp /= 0 ) then
            !
            gp( i ) = min( mixgp, gpv )
            do j = 1,i-1
               if ( site( j ) == 1 ) cycle
               if ( gp( j ) == mixgp.or.gp( j ) == gpv ) then
                  gp( j ) = gp( i )
                  k0 = nneig( j )
                  do kv = 1,neig( k0 )
                     k = neig( k0 + kv )
                     if ( site( k ) == 0 ) gp( k ) = gp( i )
                  enddo
               endif
            enddo
            !
         endif
         !
      enddo
      !
      clster = 0
      max_gp = 0
      do i = 1,obj% tot_sites
         clster( gp(i) ) = clster( gp(i) ) + 1
         max_gp = max( max_gp, gp(i) )
      enddo
      !
      nvac = clster(-1); nc = 0
      histo = 0 ; maxsize = 1
      histo( 1 ) = clster( -1 )
      !
      do i = 1,max_gp
         nvac = nvac + clster( i )
         if ( clster(i) /= 0 ) nc = nc + 1
         histo( clster(i) ) = histo( clster(i) ) + 1
         maxsize = max( maxsize, clster(i) )
      enddo
      !
      prop( 1 ) = 0
      do i = 1,maxsize
         prop( 1 ) = prop( 1 ) + real(i*histo( i ))/real( clster(-1) + nc )
      enddo
      !
      !
      call cpu_time( t_stop )
      !
      !
    end subroutine analyse
! ..................................................................................................
! ..................................................................................................








