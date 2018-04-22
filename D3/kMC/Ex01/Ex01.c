#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ------- VARIABLE DECLARATION...
   ------------------------------- */

// ::: PHYSICAL PARAMETER
double kbev = 8.61733e-5;
double kb = 1.38e-23;
double na = 6.02214e23;
double pi = 4.0 * atan(1.0);
double u = 1.66e-27;
double ang = 1e-10;

// ::: USER VARIABLE
int save_site;

// ::: STRUCTURE
struct event_type {

  int nevent, nbond, nchem_react;

  int *ptr_i_state, // array contains Initial state
      *ptr_f_state; // array contains Final state

  double *ptr_f0,    // array contains prefactor rate
      *ptr_ebarrier, // array contains activation energies
      *ptr_de,       // array contains difference energy between two miniums
      **ptr_ebond;   // 2D array contains bond energies
};

struct kmc_type {

  int bavard, conv;         // Verbose & convergence
  int period[3], nsites[3]; // Bound periodic condition & number of site

  char algorithm[50],  // algorithm name: BKL or Gillepsie
      input_file[50],  // input_file name
      input_event[50], // event file name
      libname[50],     // shared lib name with the path
      init_mod[50];    // Initialization mode: species, random, none

  int tot_sites,  // Total number of sites
      max_step,   // max step criterium
      sys_dim,    // system dimension, 1, 2 or 3 D
      freq_write, // write stat frequency
      nprop,      // number of properties compute in analyse routine
      node_state, // number of state of one site
      nspec,      // number of different sepcies
      npressure,  // number of partial pressure
      step;       // actual step

  double sum_rate, // Total rate of system
      rand_rate,   // random rate
      time,        // actual time
      rand_time,   // random time
      max_time,    // max time criterium
      temp,        // temperature value (Kelvin)
      kt,          // Energy of temperature
      per100,      // Percentage of 0 state in initial system
      f0,          // prefactor rate value
      scale;       // scale value

  int *ptr_site, *ptr_nneig, *ptr_nevt, *ptr_neig, *ptr_event_site, *ptr_spec;
  double *ptr_rate, *ptr_prop, *ptr_event_rate, *ptr_pressure, *ptr_masse;

  struct event_type event;
};

extern int *__hidden_table_MOD_h_i_state, *__hidden_table_MOD_h_f_state;
extern double *__hidden_table_MOD_h_ebarrier, *__hidden_table_MOD_h_de,
    **__hidden_table_MOD_h_ebond, *__hidden_table_MOD_h_f0;

/* ------- FUNCTION DECLARATION...
   ------------------------------- */

void builder_event_type(struct event_type *, int);

/* ------- FUNCTION DEFINITION...
   ------------------------------ */

char **parsing(char *string, int len) {
  int i = 0;
  char **pch;
  pch = malloc(len * sizeof(char));
  pch[i] = strtok(string, "\t ");
  while (pch[i] != NULL) {
    // printf( "sparse %d %s \n",i, pch[i] );
    i++;
    pch[i] = strtok(NULL, "\t ");
  }
  return pch;
}
// ...........................................................

double random_number(double min, double max) {
  return min + (double)rand() / ((double)(RAND_MAX + 1.0) * (max - min));
}
// ......................................................................
// .................................................................................
// ......................................................................

void builder_event_type(struct event_type *event, int nevent) {
  int int_size = sizeof(int);
  int double_size = sizeof(double);

  printf(" Enter in EVENT_type constructor... %d\n", nevent);
  event->nevent = nevent;

  __hidden_table_MOD_h_i_state = malloc(nevent * int_size);
  __hidden_table_MOD_h_f_state = malloc(nevent * int_size);
  __hidden_table_MOD_h_f0 = malloc(nevent * double_size);
  __hidden_table_MOD_h_ebarrier = malloc(nevent * double_size);
  __hidden_table_MOD_h_de = malloc(nevent * double_size);

  if (event->nbond != 0) {
    __hidden_table_MOD_h_ebond = malloc(event->nbond * sizeof(double *));
    for (int i = 0; i < event->nbond; i++)
      __hidden_table_MOD_h_ebond[i] = malloc(event->nbond * sizeof(double **));
  }

  event->ptr_i_state = __hidden_table_MOD_h_i_state;
  event->ptr_f_state = __hidden_table_MOD_h_f_state;
  event->ptr_f0 = __hidden_table_MOD_h_f0;
  event->ptr_ebarrier = __hidden_table_MOD_h_ebarrier;
  event->ptr_de = __hidden_table_MOD_h_de;
  if (event->nbond != 0)
    event->ptr_ebond = __hidden_table_MOD_h_ebond;

  printf(" EVENT_type constructor DONE \n");
}
// ...........................................................
//.................................................................................

void read_event(struct kmc_type *struc) {

  char string[100], **word;
  char *b = string;
  int i, id, nevent, nbond, ibd, jbd, npress;
  size_t bufsize = 100;
  size_t nword;

  FILE *fp;

  printf(" === PRINT STRUCTURE KMC ===\n");
  printf(" Periodicity %d %d %d \n", struc->period[0], struc->period[1],
         struc->period[2]);
  printf(" nsites %d %d %d\n", struc->nsites[0], struc->nsites[1],
         struc->nsites[2]);
  printf(" INPUT EVENT %s \n", struc->input_event);
  printf(" ALGORITHM   %s \n", struc->algorithm);
  printf(" INPUT FILE  %s \n", struc->input_file);
  printf(" LIBNAME     %s \n", struc->libname);

  // ::: INIT SAVE_SITE
  save_site = 0;

  // ::: READ EVENT_FILE :::

  fp = fopen(struc->input_event, "r");
  if (fp == NULL)
    perror(" PB open input_event...\n");

  nword = getline(&b, &bufsize, fp);

  while (strcmp(&string[0], "#") >= 0) {

    // SPARSE the line....
    word = parsing(string, bufsize);

    /* ==== read the event ==== */

    if (!strcmp(word[0], "Number_of_event")) {
      nevent = atoi(word[1]);
      printf(" nevent %d\n ", nevent);
      builder_event_type(&struc->event, nevent);

      for (id = 0; id < struc->event.nevent; id++) {
        nword = getline(&b, &bufsize, fp);
        word = parsing(string, bufsize);

        struc->event.ptr_i_state[id] = atoi(word[1]);
        struc->event.ptr_f_state[id] = atoi(word[2]);
        struc->event.ptr_f0[id] = atof(word[3]);
        struc->event.ptr_ebarrier[id] = atof(word[4]);
        struc->event.ptr_de[id] = atof(word[5]);
        printf(" read %d event %d %d %e %f %f \n", id,
               struc->event.ptr_i_state[id], struc->event.ptr_f_state[id],
               struc->event.ptr_f0[id], struc->event.ptr_ebarrier[id],
               struc->event.ptr_de[id]);
      }
    } // -------

    /* ==== READ PARTIAL PRESSURE ==== */

    if (!strcmp(word[0], "partial_pressure")) {

      npress = atoi(word[1]);
      if (npress != struc->npressure) {
        printf("WARNING: Number of Partial_pressure %d is different to input "
               "pressure %d\n",
               npress, struc->npressure);
        exit(1);
      }

      for (i = 0; i < npress; i++) {
        nword = getline(&b, &bufsize, fp);
        word = parsing(string, bufsize);
        id = atoi(word[0]);
        struc->ptr_pressure[id] = atof(word[1]);
        struc->ptr_masse[id] = atof(word[2]);
        printf(" PRESSURE: %d %e %f \n", id, struc->ptr_pressure[id],
               struc->ptr_masse[id]);
      }

    } // -------

    /* ==== read the energy bond ==== */

    if (!strcmp(word[0], "Energy_Bond")) {

      nbond = atoi(word[1]);

      for (id = 0; id < nbond; id++) {

        nword = getline(&b, &bufsize, fp);
        word = parsing(string, bufsize);

        ibd = atoi(word[1]);
        jbd = atoi(word[2]);
        struc->event.ptr_ebond[ibd][jbd] = atof(word[3]);
        printf(" BOND: %d %d %f --\n", ibd, jbd,
               struc->event.ptr_ebond[ibd][jbd]);
      }
    }

    // Read new line...
    nword = getline(&b, &bufsize, fp);

  } // --- END WHILE ----

  fclose(fp);

  srand(time(NULL));
}
// .................................................................................................

void event_rate_calc(struct kmc_type *obj) {}
// .................................................................................................

void choose_event(struct kmc_type *struc, int *isite, int *ievent) {

  int i, j0, jn, nevt;

  *isite = 0;
  *ievent = 0;

  double rdn = random_number(0.0, 1.0);
  double rrdn = rdn * struc->sum_rate;
  double *evt_rate = struc->ptr_event_rate;
  struc->rand_rate = rrdn;

  double rsum = 0.0;
  for (i = 0; i < struc->tot_sites; i++) {

    j0 = struc->ptr_nneig[i] - 1;
    nevt = struc->ptr_nevt[i];

    for (jn = 1; jn <= nevt; jn++) {

      *ievent = jn;
      rsum += evt_rate[j0 + jn];

      if (rsum > rrdn)
        break;
    }

    *isite = i;
    if (rsum > rrdn)
      break;
  }
  if (rsum <= rrdn)
    printf(" PB choose_event... %f %f %f %f\n", rdn, rrdn, rsum,
           struc->sum_rate);

  save_site = *isite;
}
// .................................................................................................

void event_applied(struct kmc_type *obj, int *is, int *jn) {

  int jevt, j0;

  j0 = obj->ptr_nneig[*is] - 1;
  jevt = obj->ptr_event_site[j0 + *jn];

  if (obj->ptr_site[*is] != obj->event.ptr_i_state[jevt])
    printf(" %d %d Problem site state (%d) not correspond to initial state "
           "(%d) event %d | %d \n",
           *is, *jn, obj->ptr_site[*is], obj->event.ptr_i_state[jevt], jevt,
           obj->ptr_nevt[*is]);

  obj->ptr_site[*is] = obj->event.ptr_f_state[jevt];
}
// .................................................................................................

void analyse(struct kmc_type *obj) {}
// .................................................................................................
