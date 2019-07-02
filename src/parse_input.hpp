#include <string>
#include <random>
#include <iostream>

#ifndef PARSE_INPUT_D
#define PARSE_INPUT_D

using namespace std;

#define EMPTY 0
#define FISH 1
#define SHARK 2

typedef struct {

  int species; 
  double **genome;
}specimen;




typedef struct {
  int L;
  int N; // number of genes
  int n; // number of phenotypes

  double rho_f; // Initial density of fishes
  double rho_s; // Initial density of sharks
  int N_steps; // Thermalization steps
  int N_sim; // Number of configuration to extract after the thermalization
  int N_step_between; // Number of steps between two extraction
  double N_mutations_fishes; // The average number of mutation per turn.
  double N_mutations_sharks; // The average number of mutation per turn.

  string Dir_Name; // The directory on which the data are saved

  string EVOLUTION_TYPE; // Can be "uniform", "exp" or "power"
  double EVOLUTION_EXP_CUTOFF;
  double EVOLUTION_POWER_EXP;
  double genome_norm;
  
  // The entropic forces here
  double m_fm, m_ff, m_sm, m_sf, m_sd;

  // The starting phenotypes
  double p_ff, p_sf, p_sd, p_fm, p_sm;

} options;


/*
 * This function reads the input from the given filename and
 * initializes the variables accordingly, generating
 * an options structure
 */
options read_input(string filename);

/*
 * Print the options on the standard output
 */
void PrintOptions(options opt);

#endif
