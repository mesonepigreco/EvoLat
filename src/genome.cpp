#include "genome.hpp"



void genomeInitialAssignment(int Species, options opt, double ** &genome) {
  // Initialize the array
  genome = (double**) malloc(sizeof(double*) * opt.n);
  for (int i = 0; i < opt.n; ++i)
    genome[i] = (double*) malloc(sizeof(double) * opt.N);

  if (Species == FISH) {
    for (int i = 0; i < opt.N; ++i) {
      genome[0][ i] = opt.p_ff;
      genome[1][i] = opt.p_fm;
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
  poisson_distribution<int> pdist(opt.N_mutations);

  g1->species = kind;
  g1->genome = parent.genome;

  // Allocate the second genome
  copyGenome(parent.genome, g2->genome, opt);

  int gene_x;
  int pheno;

  // Mutate the first genome
  int N_mut = pdist(opt.gen);
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
  N_mut = pdist(opt.gen);
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
  if (opt.EVOLUTION_TYPE.compare("uniform")) {
    ret = 0;
    for (int i = 0; i < opt.N; ++i) {
      ret += element.genome[phenotype][ i];
    }
    ret /= opt.N;
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

