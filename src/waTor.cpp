#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<random>

#include "parse_input.hpp"
#include "genome.hpp"

srand(time(NULL));   // should only be called once

#define fish_filiation_inibitor = 0
#define fish_rw_inibitor = 0
#define shark_filiation_inibitor = 0
#define shark_rw_inibitor = 0
#define predation_inibitor = 0


typedef struct specimen{

  int species;
  double **genome;
}specimen;

struct planet** waTor(struct options, int N_fish0, int N_shark0,   int init_config**, bool initialize = true, bool progress_bar = false){        

  int N_step = options.N_step;
  int L = options.L;
  int N_mutations = options.N_mutations;
  
  //equating specimen numbers
  int  N_fish=N_fish0;
  int  N_shark=N_shark0;
  int  N_fish_t=N_fish0;
  int  N_shark_t=N_shark0;
  
  
  bool IsSpecimentSelected;
  int Site1_Species;
  double **Site1_Genes;
  double p_move_f, p_filiation_f;
  double p_move_s,  p_filiation_s, p_death_s; 
  int newPos_x, newPos_y;
  
  //creating lattice
  specimen **planet = (specimen**) calloc(L, sizeof(specimen*));
  for(int i = 0; i < L; i++){
    planet[i] = (specimen*)calloc(L, sizeof(specimen)); 
  }; 
  
  //  starting simulation  
  int count=0;
  int tmp_x, tmp_y;
  
  if(initialize){
    
    if(progress_bar){
      printf("Initializing fishes");
    }
    
    //setting initial fishes
    while(count <= (N_fish0-1)){
    
      tmp_x= (int) lrand48*(L-1.)/RAND_MAX;
      tmp_y= (int) lrand48*(L-1.)/RAND_MAX;
      
      if(planet[tmp_x][tmp_y].species == EMPTY){
	
	planet[tmp_x][tmp_y].species = FISH;
	planet[tmp_x][tmp_y].genome = genomeInitialAssigment(FISH, n,N);
	count++;
      };
    };
    //setting sharks
    count=0;
    
    if(progress_bar){
      printf("initializing sharks\n");
    };
    
    while(count < (N_shark)){
      
      tmp_x= (int) lrand48*(L-1.)/RAND_MAX;
      tmp_y= (int) lrand48*(L-1.)/RAND_MAX;
      
      if (planet[tmp_x][tmp_y].species == EMPTY){
	
	planet[tmp_x][tmp_y].species = SHARK;
	genomeInitialAssigment(SHARK, opt, planet[tmp_x][tmp_y].genome);
	count += 1;
      };	
    };
  }else{
    copyPlanet(init_config, planet, opt);
  };
  
  //print "running"
  // running time
  int print_all_counter = 0;
  
  for(int t = 0; t < N_step; t++){
    
    //Save the configuration
    if(progress_bar){
      // sys.stdout.write("\rProgress %% %2d" % (t * 100 /N_step));
      //sys.stdout.flush();
    }
 
    
    //running a while cicle, in order to be sure that at each temporal step all specimen  move (in mean).
    for(int  j = 0; j< (L*L); j++){

      
      // extrancting a position in the LxL lattice...
      tmp_x= (int) lrand48*(L-1.)/RAND_MAX;
      tmp_y= (int) lrand48*(L-1.)/RAND_MAX;
      
      IsSpecimentSelected = !(planet[tmp_x][tmp_y].species ==EMPTY);
      if(!IsSpecimentSelected){
	continue;
      }
      
      Site1_Species = planet[tmp_x][tmp_y].species;
      
      if( Site1_Species == FISH){
	
	Site1_Genes = planet[tmp_x][tmp_y].genome;
	
	p_move_f = GetSinglePheno(planet[tmp_x][tmp_y], 0, opt); //  np.mean(Site1_Genes[0,:])
	p_filiation_f = GetSinglePheno(planet[tmp_x][tmp_y], 1, opt); //np.mean(Site1_Genes[1,:])
	prob = lrand48/((double) RAND_MAX);
	
	if(prob <= (p_move_f+p_filiation_f)){
	  
	  randomMove(tmp_x,tmp_y,L, &newPos_x, &newPos_y);  
	  Site2_Species = planet[newPos_x][newPos_y].species;
	  
	  
	  if(Site2_Species == EMPTY && prob > p_move_f && prob <= (p_filiation_f+p_move_f) && fish_filiation_inibitor == 0){
	    //filiation
	    
	    Mithosys(planet[tmp_x, tmp_y], opt, &planet[tmp_x][tmp_y], &planet[newPos_x][newPos_y]);
	    N_fish_t += 1;
	    
	  }else if(Site2_Species == EMPTY && fish_rw_inibitor == 0){
	    //fish rw
	    planet[newPos_x][newPos_y].species = FISH;
	    planet[newPos_x][newPos_y].genome = Site1_Genes;
	    planet[tmp_x][tmp_y].species = EMPTY;
	    planet[tmp_x][tmp_y].genome = NULL;
	  };
	};
      }else if(Site1_Species == SHARK){
	
	Site1_Genes = planet[tmp_x][tmp_y].genome;
	
	p_move_s = GetSinglePheno(planet[tmp_x][tmp_y], 0, opt); //#np.mean(Site1_Genes[0,:])
	p_filiation_s = GetSinglePheno(planet[tmp_x][tmp_y], 1, opt); //#np.mean(Site1_Genes[1,:])
	p_death_s = GetSinglePheno(planet[tmp_x][tmp_y], 2, opt); //#np.mean(Site1_Genes[2,:])
	
	prob_m = lrand48/((double) RAND_MAX);
	prob_f = lrand48/((double) RAND_MAX);
	prob_d = lrand48/((double) RAND_MAX);
	
	if(prob_m >= p_move_s && prob_d < p_death_s){
	  //shark spontaneous death 
	  planet[tmp_x][tmp_y].species = EMPTY;
	  planet[tmp_x][tmp_y].genome = NULL;
	  
	  N_shark_t -= 1;
	}else if(prob_m < p_move_s){
	  
	  randomMove(tmp_x,tmp_y, L, &newPos_x, &newPos_y);
	  Site2_Species = planet[newPos_x][newPos_y].species;
	  
	  // PAY ATTENTION! sharks can NOT die when the meet other sharks!
	  
	  if(Site2_Species == FISH && prob_f <= p_filiation_s && shark_filiation_inibitor == 0  && predation_inibitor == 0){
	    //shark filiation with predation
	    Mithosys(planet[tmp_x][tmp_y], opt, &planet[tmp_x][tmp_y], &planet[newPos_x][newPos_y]);
	    N_shark_t += 1;
	    N_fish_t -= 1;
	    
	  }else if(Site2_Species == FISH && predation_inibitor == 0){
	    
	    //shark predation
	    planet[newPos_x][newPos_y].species = SHARK;
	    planet[tmp_x][tmp_y].species =  EMPTY;
	    
	    planet[newPos_x][newPos_y].genome = Site1_Genes;
	    planet[tmp_x][tmp_y].genome = NULL;
	    
	    N_fish_t -= 1;
	  }else if(Site2_Species == EMPTY && shark_rw_inibitor == 0){
	    
	    //shark rw
	    planet[newPos_x][newPos_y].species = SHARK;
	    planet[tmp_x][tmp_y].species = EMPTY;
	    planet[newPos_x][newPos_y].genome = Site1_Genes;
	    planet[tmp_x][tmp_y]genome = NULL;
	    
	    
	    if(prob_d < (p_death_s)){
	      //shark death 
	      planet[newPos_x][newPos_y].species = EMPTY;
	      planet[newPos_x][newPos_y].genome  = NULL;
	      N_shark_t -= 1;
	    };
	  };
	};
      };
    };
    
    N_fish = N_fish_t;
    N_shark= N_shark_t;
  }
  return planet;
}
