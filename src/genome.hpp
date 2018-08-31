#include "parse_input.hpp"

#ifndef GENOME_D
#define GENOME_D
/*
 * This file contains the genome of the simulation
 */


void genomeInitialAssignment(int Species, options opt, double ** &genome);


// Creates the new genome for the two sons of specimen.
// The genome is created with an average N_mutations mutations .
void Mythosis(specimen parent, options opt, specimen * son1, specimen * son2);

double GetSinglePheno(specimen element, int phenotype, options opt);

double P_mut(int kind, int pheno, options opt) ;

// Copy the genom into the dest that is also allocated
void copyGenome(double ** original, double ** &dest, options opt);


// Copy planet does not allocate the dest
void copyPlanet(specimen ** original, specimen ** dest, options opt);
void cleanGenome(double ** genome, options opt); // Death, the result is null

void savePlanet(char* filename, options opt, specimen ** planet);

void randomMove(int tmp_x, int tmp_y, int L, int * dest_x, int * dest_y) ;

void FreeGenome(double ** , options);

void DestroyAll(options opt, specimen** planet);
#endif
