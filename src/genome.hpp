#include "parse_input.hpp"

/*
 * This file contains the genome of the simulation
 */


void genomeInitialAssignment(int Species, options opt, double ** &genome);


// Creates the new genome for the two sons of specimen.
// The genome is created with an average N_mutations mutations .
void Mythosis(specimen parent, options opt, specimen * son1, specimen * son2);

double GetSinglePheno(specimen element, int phenotype, options opt);

double P_mut(int kind, int pheno, options opt) ;

void copyGenome(double ** original, double ** dest, options opt);
void copyPlanet(specimen ** original, specimen ** dest, options opt);
void cleanGenome(double ** genome, options opt); // Death, the result is null

void randomMove(int tmp_x, int tmp_y, int L, int * dest_x, int * dest_y) ;
