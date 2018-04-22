#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* ------- VARIABLE DECLARATION... 
   ------------------------------- */

// ::: PHYSICAL PARAMETER
double kb = 1.38065e-23; 
double kbev = 8.61733e-5;
double na = 6.02214e23;
double pi = 4.0*atan(1.0);
double u = 1.66e-27;
double ang = 1e-10;

// ::: STRUCTURE
struct event_type {

  int nevent, nbond, nchem_react;

  int *ptr_i_state,   // array contains Initial state 
      *ptr_f_state;   // array contains Final state 

  double *ptr_f0,          // array contains prefactor rate 
         *ptr_ebarrier,    // array contains activation energies
         *ptr_de,          // array contains difference energy between two miniums
         **ptr_ebond;      // 2D array contains bond energies
};

struct kmc_type {

  int bavard, conv;          // Verbose & convergence
  int period[3], nsites[3];  // Bound periodic condition & number of site 

  char algorithm[50],        // algorithm name: BKL or Gillepsie
       input_file[50],       // input_file name
       input_event[50],      // event file name
       libname[50],          // shared lib name with the path
       init_mod[50];         // Initialization mode: species, random, none 

  int tot_sites,             // Total number of sites
      max_step,              // max step criterium
      sys_dim,               // system dimension, 1, 2 or 3 D
      freq_write,            // write stat frequency
      nprop,                 // number of properties compute in analyse routine
      node_state,            // number of state of one site
      nspec,                 // number of different sepcies
      npressure,             // number of partial pressure
      step;                  // actual step

  double sum_rate,           // Total rate of system 
         rand_rate,          // random rate 
         time,               // actual time
         rand_time,          // random time
         max_time,           // max time criterium
         temp,               // temperature value (Kelvin)
         kt,                 // Energy of temperature
         per100,             // Percentage of 0 state in initial system
         f0,                 // prefactor rate value
         scale;              // scale value 

  int     *ptr_site, *ptr_nneig, *ptr_nevt, *ptr_neig, *ptr_event_site, *ptr_spec;
  double  *ptr_rate, *ptr_prop, *ptr_event_rate, *ptr_pressure, *ptr_masse;

  struct event_type event;
};


extern int *__hidden_table_MOD_h_i_state, *__hidden_table_MOD_h_f_state;
extern double *__hidden_table_MOD_h_f0, *__hidden_table_MOD_h_ebarrier, *__hidden_table_MOD_h_de, **__hidden_table_MOD_h_ebond;


/* ------- FUNCTION DECLARATION...
   ------------------------------- */

void builder_event_type( struct event_type *, int * );


/* ------- FUNCTION DEFINITION... 
   ------------------------------ */

char ** parsing( char *string, int len ) {
  int i = 0;
  char **pch;
  pch = malloc( len*sizeof(char) );
  pch[i] = strtok( string, "\t " );
  while ( pch[i] != NULL ) {
   // printf( "sparse %d %s \n",i, pch[i] );
    i++;
    pch[i] = strtok( NULL, "\t " );
  }
  return pch;
}
// ...........................................................

// ......................................................................

double random_number( double min, double max ) {
   return min + (double)rand() / ( (double)(RAND_MAX + 1.0)*( max - min) );
}
// ......................................................................
// ......................................................................

void builder_event_type( struct event_type *event, int *nevent ) {
  int int_size = sizeof( int );
  int double_size = sizeof( double );

  printf( " Enter in EVENT_type constructor...\n" );
  event->nevent = *nevent;

  __hidden_table_MOD_h_i_state = malloc( *nevent*int_size );
  __hidden_table_MOD_h_f_state = malloc( *nevent*int_size );
  __hidden_table_MOD_h_f0 = malloc( *nevent*double_size );
  __hidden_table_MOD_h_ebarrier = malloc( *nevent*double_size );
  __hidden_table_MOD_h_de = malloc( *nevent*double_size );

  __hidden_table_MOD_h_ebond = malloc( event->nbond*sizeof( double* ) );
  for (int i = 0; i < event->nbond; i++ ) 
    __hidden_table_MOD_h_ebond[i] = malloc( event->nbond*sizeof( double** ) );

  event->ptr_i_state  = __hidden_table_MOD_h_i_state;
  event->ptr_f_state  = __hidden_table_MOD_h_f_state;
  event->ptr_ebarrier = __hidden_table_MOD_h_ebarrier;
  event->ptr_f0       = __hidden_table_MOD_h_f0; 
  event->ptr_de       = __hidden_table_MOD_h_de;
  event->ptr_ebond    = __hidden_table_MOD_h_ebond;

  printf( " EVENT_type constructor DONE \n" );

}
// ...........................................................

void read_event( struct kmc_type *struc ) {

  char string[100], **word;
  char *b = string;
  int i, id, nevent, nbond, ibd, jbd;
  size_t bufsize = 100;
  size_t nword;

  FILE *fp;

  fp = fopen( struc->input_event, "r" );
  if ( fp == NULL ) perror( " PB open input_event...\n" );

  nword = getline(&b, &bufsize, fp );

  while ( strcmp(&string[0],"#") >= 0 ) {

     word = parsing( string, bufsize );

     if ( !strcmp(word[0],"Number_of_event") ) {
        nevent = atoi( word[1] );
        printf( " nevent %d\n ",nevent );
        builder_event_type( &struc->event, &nevent );

        /* ==== read the event ==== */
        nword = getline(&b, &bufsize, fp );
        int n = struc->event.nevent ;
        for ( id = 0; id < struc->event.nevent; id ++ ) {
            word = parsing( string, bufsize );

            struc->event.ptr_i_state[ id ]  = atoi( word[1] );
            struc->event.ptr_f_state[ id ]  = atoi( word[2] );
            struc->event.ptr_f0[ id ]       = atof( word[3] );
            struc->event.ptr_ebarrier[ id ] = atof( word[4] );
            struc->event.ptr_de[ id ]       = atof( word[5] );
            printf( " read %d event %d %d %f %f \n", id, struc->event.ptr_i_state[ id ], struc->event.ptr_f_state[ id ],
                  struc->event.ptr_de[ id ], struc->event.ptr_ebarrier[ id ], struc->event.ptr_de[ id ] );

            nword = getline(&b, &bufsize, fp );
         }

     }

     if( !strcmp( word[0], "Energy_Bond") ) {
         nbond = atoi( word[1] );
         printf( " dans energy_bond %s %d \n",word[0],nbond );

         /* ==== read Energy Bond ==== */
         for ( id = 0; id < nbond; id++ ) {
             nword = getline(&b, &bufsize, fp );
             word = parsing( string, bufsize );
             ibd = atoi( word[1] );
             jbd = atoi( word[2] );
             struc->event.ptr_ebond[ibd][jbd] = atof( word[3] );
             //printf( " %d %d %f %f --\n", ibd, jbd, struc->event.ptr_ebond[ibd][jbd] );
         }

      }

     nword = getline(&b, &bufsize, fp );

  } // --- WHILE ----


  fclose( fp );

  srand( time(NULL) );
 // exit(1);

}
// .................................................................................................

void event_rate_calc( struct kmc_type *struc ) {


}
// ....................................................................................
// ..................................................................................................

void choose_event( struct kmc_type *struc, int *isite, int *ievent ) {
   
   int i, j0, jn, j, nvj;

   *isite = 0;
   *ievent = 0;
   
   double rdn = random_number( 0.0, 1.0 );
   double rrdn = rdn*struc->sum_rate;
   double *evt_rate = struc->ptr_event_rate;
   struc-> rand_rate = rrdn;
   double rsum = 0.0;

   for (i = 0; i < struc->tot_sites; i++) {
      
      j0 = struc->ptr_nneig[ i ] - 1;
      nvj = struc->ptr_neig[ j0 ];

      for ( jn = 1; jn <= nvj; jn++) {

          *ievent = jn;
          rsum += evt_rate[ j0 + jn ];
          if ( rsum > rrdn ) break;

      }

      *isite = i;
      if ( rsum > rrdn ) break;

   }

   if ( rsum <= rrdn ) 
     printf( " PB choose_event... %f %f %f\n", rrdn, rsum, struc->sum_rate);
}
// ..................................................................................................

void event_applied( struct kmc_type *struc, int *is, int *jn ) {

  int j, j0;
  
  if ( struc->ptr_site[ *is ] == 1 ) {
    printf( " site Flip Problem 0 = 1 : %d %d \n", *is, struc->ptr_site[ *is ]);
  }

  struc->ptr_site[ *is ] = 1;
  j0 = struc->ptr_nneig[ *is ] - 1;
  j = struc->ptr_neig[ j0 + *jn ] - 1;

  if ( struc->ptr_site[ j ] == 0 ) 
    printf( " site Flip Problem 1 = 0 : %d %d \n", j, struc->ptr_site[ j ]);

  struc->ptr_site[ j ] = 0  ;
}
// ..................................................................................................

void analyse( struct kmc_type *obj ) {

  int i,j, j0, jv, jn,k0, k, kv, nvj, nvk, ngp, gpv, nc, mixgp, nvac, maxsize, max_gp;

  int gp[obj->tot_sites], histo[obj->tot_sites];
  int clster[obj->tot_sites];

  for ( i = 0; i < obj->tot_sites; i++) {
    clster[ i ] = 0;
    gp[ i ] = 0;
  }
/* ----------------------------------------
   1) We identify the cluster:
      gp = 0        : no already treated
      gp = 1        : size = 1
      gp = {2,..,N} : size > 1
   2) We count the size of cluster 
   ---------------------------------------- */

  ngp = 1; max_gp = 0; nvac = 0;
  for ( i = 0; i < obj->tot_sites; i++ ) {

      if ( obj->ptr_site[ i ] == 1 ) continue;

      nvac += 1;
      nc = 0;
      gp[ i ] = 1;
      gpv = 0; mixgp = 0;

      j0 = obj->ptr_nneig[ i ] - 1;
      nvj = obj->ptr_neig[ j0 ];

      for ( jn = 1; jn <= nvj; jn++ ) {

          j = obj->ptr_neig[ j0 + jn ] - 1;
          
          if ( obj->ptr_site[ j ] == 1 ) continue; 

          if ( gp[ j ] != 0 && gpv == 0 ) {
             gpv = gp[ j ];
          } else if ( gp[ j ] != 0 && gpv != 0 && gp[ j ] != gpv ) {
             mixgp = gp[ j ];
          }

          nc += 1;

      } // ----

      if ( nc != 0 && gpv != 0 && mixgp == 0 ) {

         gp[ i ] = gpv;

      } else if ( nc != 0 && gpv == 0 ) {

         ngp += 1;
         gp[ i ] = ngp;

         j0 = obj->ptr_nneig[ i ] - 1;
         nvj = obj->ptr_neig[ j0 ];

         for ( jv = 1; jv <= nvj; jv++ ) {
             j = obj->ptr_neig[ j0 + jv ] - 1;
             if ( obj->ptr_site[ j ] == 0 ) gp[ j ] = gp[ i ];
         }

      } else if ( mixgp != 0 ) {

         gp[ i ] = (int)fmin( mixgp, gpv );

         for ( j = 0; j < i; j++ ) {

             if ( obj->ptr_site[ j ] == 1 ) continue;

             if ( gp[ j ] == mixgp || gp[ j ] == gpv ) {
                gp[ j ] = gp[ i ];

                k0 = obj->ptr_nneig[ j ] - 1;
		nvk = obj->ptr_neig[ k0 ];

                for ( kv = 1; kv <= nvk; kv++ ) {
		    k = obj->ptr_neig[ k0 + kv ] - 1;
                    if ( obj->ptr_site[ k ] == 0 ) gp[ k ] = gp[ i ];
                } 

             }
         } 

      } // --- else if

  } // ---- for "i"

  for ( i = 0; i < obj->tot_sites; i++ ) {
      clster[ gp[ i ] ] += 1;
      max_gp = fmax( max_gp, gp[ i ] );
      histo[ i ] = 0;
  }

  nvac = 0;  nc = 0;
  maxsize = 1;
  for ( i = 0; i < max_gp; i++ ) {
      nvac += clster[ i ];
      if ( clster[i] != 0 ) nc += 1;
      histo[ clster[i] ] += 1;
      maxsize = (int)fmax( maxsize, clster[i] ); 
  }
  obj->ptr_prop[ 0 ] = 0.;
  for ( i = 0; i < maxsize; i++ ) {
      obj->ptr_prop[ 0 ] += ((double)( (i+1)*histo[i] ))/((double)( clster[1] + nc ));
  }
}

//! ..................................................................................................






