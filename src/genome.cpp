#include "genome.hpp"
#include "math.h"
#define DEBUG_GEN 0


int get_poisson(double mean) {
  double L = exp(-mean);
  int k = 0;
  double p = 1;
  double u;

  do {
    k++;
    u = rand() / (double) RAND_MAX;
    p *= u;
  } while (p > L);

  return k -1;
}

void genomeInitialAssignment(int Species, options opt, double ** &genome) {
  // Initialize the array
  genome = (double**) malloc(sizeof(double*) * opt.n);
  for (int i = 0; i < opt.n; ++i)
    genome[i] = (double*) malloc(sizeof(double) * opt.N);

  if (Species == FISH) {
    for (int i = 0; i < opt.N; ++i) {
      genome[0][ i] = opt.p_ff;
      genome[1][i] = opt.p_fm;
      genome[2][i] = 0;
    }
  } else if (Species == SHARK) {
    for (int i = 0; i < opt.N; ++i) {
      genome[0][ i] = opt.p_sm;
      genome[1][ i] = opt.p_sf;
      genome[2][ i] = opt.p_sd;
    }
  }
}

double P_mut(int kind, int pheno, options opt) {

  double x_av;
  if (kind == FISH && pheno == 0) x_av = opt.m_fm;
  else if (kind == FISH && pheno == 1) x_av = opt.m_ff;
  else if (kind == SHARK && pheno == 0) x_av = opt.m_sm;
  else if (kind == SHARK && pheno == 1) x_av = opt.m_sf;
  else if (kind == SHARK && pheno == 2) x_av = opt.m_sd;
  else {
    cerr << "Error in P_mut: kind and pheno must be good." << endl;
    cerr << "KIND: " <<  kind << "  PHENO: " << pheno << endl;
    exit(EXIT_FAILURE);
  }
  
  double u = rand() / (double) RAND_MAX;
  double alpha, m;

  if (x_av > .5) {
    alpha = (2 * x_av - 1) / (1 - x_av);
    m = pow(u, 1. / (alpha+ 1));
  } else {
    alpha = (1 - 2 * x_av) / (1. - x_av);
    m = 1 - pow(u, 1./ (alpha + 1.));
  }
  return m;
}

void Mythosis(specimen parent, options opt, specimen * g1, specimen * g2) {
  int kind = parent.species;

  // Get the number of mutations
  //poisson_distribution<int> pdist(opt.N_mutations);

  g1->species = kind;
  g2->species = kind;
  g1->genome = parent.genome;
  
  if (g2->genome != NULL) {
    if (DEBUG_GEN) cout <<  "      killing and reallocating " << g2->genome << endl;
    FreeGenome(g2->genome, opt);
  }
  
  // Allocate the second genome and copy
  copyGenome(parent.genome, g2->genome, opt);

  int gene_x;
  int pheno;

  // Mutate the first genome
  int N_mut = get_poisson(opt.N_mutations);
  if (DEBUG_GEN)  cout << "Extracted " << N_mut << " mutations." <<endl;
  for (int i = 0; i < N_mut; ++i) {
    // Get the mutating gene
    gene_x = rand() % opt.N;
    if (kind == FISH)
      pheno = rand() % 2;
    else if (kind == SHARK)
      pheno = rand() % 3;

    g1->genome[pheno][ gene_x] = P_mut(kind, pheno, opt);
  }
  
  // Mutate the second genome
  N_mut = get_poisson(opt.N_mutations);
  for (int i = 0; i < N_mut; ++i) {
    // Get the mutating gene
    gene_x = rand() % opt.N;
    if (kind == FISH)
      pheno = rand() % 2;
    else if (kind == SHARK)
      pheno = rand() % 3;

    g2->genome[pheno][ gene_x] = P_mut(kind, pheno, opt);
  }
}


double GetSinglePheno(specimen element, int phenotype, options opt) {
  double ret;
  if (opt.EVOLUTION_TYPE.compare("uniform") == 0) {
    ret = 0;
    for (int i = 0; i < opt.N; ++i) {
      if (DEBUG_GEN) {
	//cout << "         [GET SINGLE PHENO ] i = " << i << endl;
	//cout << "                             gene = " << element.genome[phenotype][ i] << endl;
      }
      ret += element.genome[phenotype][ i];
    }
    ret /= opt.N;
  } else if (opt.EVOLUTION_TYPE.compare("power") == 0) {
    ret = 0;
    for (int i = 0; i < opt.N; ++i) 
      ret += element.genome[phenotype][i] * pow(i+1, opt.EVOLUTION_POWER_EXP); 
    ret /= opt.genome_norm;
  } else {
    cerr << "ERROR, EVOLUTION_TYPE = " << opt.EVOLUTION_TYPE << " not yet implemented." << endl;
    exit(EXIT_FAILURE);
  }
  return ret;
}



void randomMove(int tmp_x, int tmp_y, int L, int * dest_x, int * dest_y) {
  double r = rand() % 3;
  double x, y, newPos_x, newPos_y;
  if(r==0) {
    x=1;
    y=0;
  }
  else if ( r==1) {
    x=-1;
    y=0;
  }
  else if (r==2) {
    x=0;
    y=1;
  }
  else{
    x=0;
    y=-1;
  }
                        
  newPos_x=tmp_x+x;
  newPos_y=tmp_y+y;
    
  // Boundary conditions
  if(newPos_x == L)
    newPos_x=0;
  else if (newPos_x == -1)
    newPos_x=L-1;

  if(newPos_y == L)
    newPos_y=0;
  else if(newPos_y == -1)
    newPos_y=L-1;

  *dest_x = newPos_x;
  *dest_y = newPos_y;
}


void copyGenome(double ** original, double ** &dest, options opt) {
  // Allocate the genome
  dest = (double ** )malloc(sizeof(double*) * opt.n);
  if (DEBUG_GEN) cout << "    GENOME ALLOCATION : " << dest << endl;
  for (int i = 0; i < opt.n; ++i) {
    dest[i] = (double*) malloc(sizeof(double) * opt.N);
    if (DEBUG_GEN) cout << "                    : undergene i = " << i << " = " << dest[i] << endl;
    for (int j = 0; j < opt.N; ++j) {
      dest[i][j] = original[i][j];
    }
  }
}


void copyPlanet(specimen ** original, specimen ** dest, options opt) {
  for (int x = 0; x < opt.L; ++x) {
    for (int y = 0; y < opt.L; ++y) {
      dest[x][y].species = original[x][y].species;
      copyGenome(original[x][y].genome, dest[x][y].genome, opt);
    }
  } 
}

void FreeGenome(double ** genome, options opt) {
  for (int i = 0; i < opt.n; ++i) {
    if (DEBUG_GEN) cout << "killing " << i << " pointer " << genome[i] <<  endl;
    free(genome[i]);
  }
  free(genome);
}


void DestroyAll(options opt, specimen** planet) {
  for (int x = 0; x <  opt.L; ++x) {
    for (int y = 0; y < opt.L; ++y) {
      if (planet[x][y].species == EMPTY) {
	if (planet[x][y].genome != NULL) {
	  cerr << "ERROR, site (" << x << " " << y << ") is EMPTY but not NULL" << endl;
	  exit(EXIT_FAILURE);
	}
      } else {
	FreeGenome(planet[x][y].genome, opt);
      }
    }
  }

  cout << "DISTRUCTOR FOUND NO ERRORS" << endl;
  exit(0);
}
