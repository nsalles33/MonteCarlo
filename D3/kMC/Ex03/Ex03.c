#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* ------- VARIABLE DECLARATION... 
   ------------------------------- */

// ::: PHYSICAL PARAMETER
double kbev = 8.61733e-5;
double kb = 1.38e-23;
double na = 6.02214e23;
double pi = 4.0*atan(1.0);
double u = 1.66e-27;
double ang = 1e-10;

// ::: USER VARIABLE

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
extern double *__hidden_table_MOD_h_ebarrier, *__hidden_table_MOD_h_de, **__hidden_table_MOD_h_ebond, 
              *__hidden_table_MOD_h_f0;


/* ------- FUNCTION DECLARATION...
   ------------------------------- */

void builder_event_type( struct event_type *, int );


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
// .................................................................................
// ......................................................................

void alloc_chem_event_type( struct event_type *event, int nspec ) {
  int int_size = sizeof( int );
  int double_size = sizeof( double );
  int n = event-> nevent + event-> nchem_react*nspec;

  printf( " Enter in EVENT_type constructor... %d %d\n", event->nchem_react, event->nevent );

  __hidden_table_MOD_h_i_state  = malloc( n*int_size );
  __hidden_table_MOD_h_f_state  = malloc( n*int_size );
  __hidden_table_MOD_h_f0       = malloc( n*double_size );
  __hidden_table_MOD_h_ebarrier = malloc( n*double_size );
  __hidden_table_MOD_h_de       = malloc( n*double_size );

  if ( event->nbond != 0 ) {
     __hidden_table_MOD_h_ebond = malloc( event->nbond*sizeof( double* ) );
     for (int i = 0; i < event->nbond; i++ )
        __hidden_table_MOD_h_ebond[i] = malloc( event->nbond*sizeof( double** ) );
  }

  event->ptr_i_state  = __hidden_table_MOD_h_i_state;
  event->ptr_f_state  = __hidden_table_MOD_h_f_state;
  event->ptr_f0       = __hidden_table_MOD_h_f0;
  event->ptr_ebarrier = __hidden_table_MOD_h_ebarrier;
  event->ptr_de       = __hidden_table_MOD_h_de;
  if ( event->nbond != 0 ) 
     event->ptr_ebond    = __hidden_table_MOD_h_ebond;

  printf( " EVENT_type constructor DONE \n" );

}
// ...........................................................
//.................................................................................

void read_event( struct kmc_type *struc ) {

  char string[100], **word;
  char *b = string;
  int i, id, j, nevent, nbond, ibd, jbd, nspec;
  size_t bufsize = 100;
  size_t nword;

  FILE *fp;

  // ::: READ EVENT_FILE :::

  fp = fopen( struc->input_event, "r" );
  if ( fp == NULL ) perror( " PB open input_event...\n" );

  nword = getline(&b, &bufsize, fp );

  while ( strcmp(&string[0],"#") >= 0 ) {

     // SPARSE the line....
     word = parsing( string, bufsize );


     /* ==== read the event ==== */
     if ( !strcmp(word[0],"Number_of_reaction") ) {
        struc->event.nchem_react = atoi( word[1] );
        alloc_chem_event_type( &struc->event, nspec );
        //break;

        for ( i = 0; i < struc->event.nchem_react; i ++ ) {
            nword = getline(&b, &bufsize, fp );
            word = parsing( string, bufsize );

            id = i*struc->event.nchem_react;
            for ( j = 1; j <= struc->nspec; j++ ) {
                struc->event.ptr_i_state[ id + j ] = atoi( word[ j ] );
                struc->event.ptr_f_state[ id + j ] = atoi( word[ struc->nspec + j ] );
                printf( " %d  i %d f %d ", id + j, struc->event.ptr_i_state[ id + j ], struc->event.ptr_f_state[ id + j ] );
            }
            struc->event.ptr_ebarrier[ i ] = atof( word[ 2*struc->nspec + 1 ] );
         }
     } // -------


     /* ==== READ PARTIAL PRESSURE ==== */

     if ( !strcmp( word[0], "Init_species" ) ) {
        
	nspec = atoi( word[1] );
        if ( nspec != struc-> nspec ) {
            printf( "WARNING: Number of species %d is different to input species %d\n", nspec, struc-> nspec );
            exit(1);
         }
         
         for ( i = 0; i < nspec; i++ ) {
             nword = getline(&b, &bufsize, fp );
             word = parsing( string, bufsize );
             struc->ptr_spec[ i ] = atof( word[0] );
             printf( " species: %d %d\n", i, struc->ptr_spec[i] );
         }

     } // -------


     // Read new line...
     nword = getline(&b, &bufsize, fp );


  } // --- END WHILE ----


  fclose( fp );

  srand( time(NULL) );

  //exit(1);

}
// .................................................................................................

void event_rate_calc( struct kmc_type *obj ) {
  
  int is ; 
  double h ;


  for ( is = 0; is < obj->tot_sites; is++ ) {


  } // -- For "is"

}
// .................................................................................................

void choose_event( struct kmc_type *struc, int *isite, int *ievent ) {

  int i, j0, jn, j, nevt, nber;
  double rdn, tmin = 1000.0, t;

   *isite = 0;
   *ievent = 0;


}
// .................................................................................................

void event_applied( struct kmc_type *obj, int *is, int *jn ) {


}
// .................................................................................................

void analyse( struct kmc_type *obj ) {

  int i;

  if ( obj->nprop != obj->nspec ) {
     printf( " Pb nprop %d != nspec %d \n", obj->nprop, obj->nspec );
     return;
  }

  obj->ptr_prop[ 0 ] = 0.0;
  for ( i = 0; i < obj->nspec; i++ ) { 
      obj->ptr_prop[ i ] = obj->ptr_spec[ i ];
  }

}
// .................................................................................................
















