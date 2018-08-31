#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<random>

#include "parse_input.hpp"
#include "genome.hpp"

//srand(time(NULL));   // should only be called once

#define fish_filiation_inibitor  0
#define fish_rw_inibitor  0
#define shark_filiation_inibitor  0
#define shark_rw_inibitor  0
#define predation_inibitor  0

#define DW 1


specimen** waTor(options opt, int N_steps, int N_fish0, int N_shark0,   specimen **init_config, bool initialize = true, bool progress_bar = false){        

  int L = opt.L;
  int N_mutations = opt.N_mutations;
  
  //equating specimen numbers
  int  N_fish=N_fish0;
  int  N_shark=N_shark0;
  int  N_fish_t=N_fish0;
  int  N_shark_t=N_shark0;
  
  
  bool IsSpecimentSelected;
  int Site1_Species, Site2_Species;
  double **Site1_Genes;
  double p_move_f, p_filiation_f;
  double p_move_s,  p_filiation_s, p_death_s; 
  int newPos_x, newPos_y;

  double prob_m, prob_d, prob_f;
  
  //creating lattice
  specimen **planet = (specimen**) calloc(L, sizeof(specimen*));
  for(int i = 0; i < L; i++){
    planet[i] = (specimen*)calloc(L, sizeof(specimen)); 
  }; 
  
  //  starting simulation  
  int count=0;
  int tmp_x, tmp_y;
  
  if(initialize){
    if (DW) cout << "Initizzed" << endl;
    
    if(progress_bar){
      printf("Initializing fishes");
    }
    
    //setting initial fishes
    while(count <= (N_fish0-1)){
    
      tmp_x= rand() % L;
      tmp_y=  rand() % L;
      
      if(planet[tmp_x][tmp_y].species == EMPTY){
	
	planet[tmp_x][tmp_y].species = FISH;
	genomeInitialAssignment(FISH, opt, planet[tmp_x][tmp_y].genome);
	count++;
      };
    };

    if (DW) cout << "Fish initialized" << endl;
    //setting sharks
    count=0;
    
    if(progress_bar){
      printf("initializing sharks\n");
    };
    
    while(count < (N_shark)){
      
      tmp_x= rand() % L;
      tmp_y=  rand() % L;
      if (planet[tmp_x][tmp_y].species == EMPTY){
	
	planet[tmp_x][tmp_y].species = SHARK;
	genomeInitialAssignment(SHARK, opt, planet[tmp_x][tmp_y].genome);
	count += 1;
      };	
    };
    if (DW) cout << "Shark initialized" << endl;
  }else{
    planet = init_config;
    //copyPlanet(init_config, planet, opt);
  };

  if (DW) cout << "Planet is ready for the sim." << endl;
  
  //print "running"
  // running time
  int print_all_counter = 0;
  double prob;
  
  for(int t = 0; t < N_steps; t++){
    
    //Save the configuration
    if(progress_bar){
      // sys.stdout.write("\rProgress %% %2d" % (t * 100 /N_step));
      //sys.stdout.flush();
    }
 
    if (DW) cout << "Sim step " << t <<  endl;
    //running a while cicle, in order to be sure that at each temporal step all specimen  move (in mean).
    for(int  j = 0; j< (L*L); j++){
      if (DW) cout << "   -> parastep i = " << t << " j = " << j << " / " << L*L << endl; 
      
      // extrancting a position in the LxL lattice...

      tmp_x= rand() % L;
      tmp_y= rand() % L;


      
      IsSpecimentSelected = !(planet[tmp_x][tmp_y].species ==EMPTY);
      if(!IsSpecimentSelected){
	continue;
      }
      if (DW) cout << "Selected (" << tmp_x << " "<<tmp_y << ") : " << planet[tmp_x][tmp_y].species << endl << "   === >";
      Site1_Species = planet[tmp_x][tmp_y].species;
      
      if( Site1_Species == FISH){

	
	Site1_Genes = planet[tmp_x][tmp_y].genome;
	if (DW) cout << "ABACAT1 " << Site1_Genes << endl;
	p_move_f = GetSinglePheno(planet[tmp_x][tmp_y], 0, opt); //  np.mean(Site1_Genes[0,:])
	if (DW) cout << "ABACAT2" << endl;
	p_filiation_f = GetSinglePheno(planet[tmp_x][tmp_y], 1, opt); //np.mean(Site1_Genes[1,:])
	prob = rand()/((double) RAND_MAX);
	if(prob <= (p_move_f+p_filiation_f)){
	  
	  randomMove(tmp_x,tmp_y,L, &newPos_x, &newPos_y);  
	  Site2_Species = planet[newPos_x][newPos_y].species;
	  if (DW) cout << " Move in ("<< newPos_x << " " << newPos_y <<")" << endl;	  
	  
	  if(Site2_Species == EMPTY && prob > p_move_f && prob <= (p_filiation_f+p_move_f) && fish_filiation_inibitor == 0){
	    //filiation
	    if (DW) cout << "     ====> Filiation" << endl;;
	    Mythosis(planet[tmp_x][ tmp_y], opt, &planet[tmp_x][tmp_y], &planet[newPos_x][newPos_y]);
	    N_fish_t += 1;
	    
	  }else if(Site2_Species == EMPTY && fish_rw_inibitor == 0){
	    //fish rw
	    if (DW) cout << "     ====> RW" << endl;
	    planet[newPos_x][newPos_y].species = FISH;
	    planet[newPos_x][newPos_y].genome = Site1_Genes;
	    planet[tmp_x][tmp_y].species = EMPTY;
	    planet[tmp_x][tmp_y].genome = NULL;
	  }
	  else { if (DW) cout << "     =====> Occupied" <<  endl;}
	} else { if (DW) cout << endl;}
      }else if(Site1_Species == SHARK){

	if (DW) cout << "SITIAO   " << planet[tmp_x][tmp_y].genome << endl;
	Site1_Genes = planet[tmp_x][tmp_y].genome;
	if (DW) cout << "ABACAT" << endl;
	
	p_move_s = GetSinglePheno(planet[tmp_x][tmp_y], 0, opt); //#np.mean(Site1_Genes[0,:])
	if (DW) cout << "ABACAT1" << endl;
	p_filiation_s = GetSinglePheno(planet[tmp_x][tmp_y], 1, opt); //#np.mean(Site1_Genes[1,:])
	if (DW) cout << "ABACAT2" << endl;
	p_death_s = GetSinglePheno(planet[tmp_x][tmp_y], 2, opt); //#np.mean(Site1_Genes[2,:])
	if (DW) cout << "ABACAT3" << endl;
	prob_m = rand()/((double) RAND_MAX);
	prob_f = rand()/((double) RAND_MAX);
	prob_d = rand()/((double) RAND_MAX);

	if (DW) cout << " m:" << prob_m << " f:" << prob_f << " d:" << prob_d << endl;
	
	if(prob_m >= p_move_s && prob_d < p_death_s){
	  //shark spontaneous death
	  if (DW) cout << " Death killed " << planet[tmp_x][tmp_y].genome  <<  endl;
	  planet[tmp_x][tmp_y].species = EMPTY;
	  FreeGenome(planet[tmp_x][tmp_y].genome, opt);
	  planet[tmp_x][tmp_y].genome = NULL;
	  
	  N_shark_t -= 1;
	}else if(prob_m < p_move_s){
	  
	  randomMove(tmp_x,tmp_y, L, &newPos_x, &newPos_y);
	  Site2_Species = planet[newPos_x][newPos_y].species;
	  
	  // PAY ATTENTION! sharks can NOT die when the meet other sharks!
	  if (DW) cout << " move in ("<< newPos_x << " " << newPos_y <<") ";
	  if(Site2_Species == FISH && prob_f <= p_filiation_s && shark_filiation_inibitor == 0  && predation_inibitor == 0){
	    //shark filiation with predatio
	    if (DW) cout << " Filiation" << endl;
	    Mythosis(planet[tmp_x][tmp_y], opt, &planet[tmp_x][tmp_y], &planet[newPos_x][newPos_y]);
	    N_shark_t += 1;
	    N_fish_t -= 1;
	    
	  }else if(Site2_Species == FISH && predation_inibitor == 0){
	    if (DW) cout << " Predation killed " << planet[newPos_x][newPos_y].genome <<  endl;
	    //shark predation
	    planet[newPos_x][newPos_y].species = SHARK;
	    planet[tmp_x][tmp_y].species =  EMPTY;
	    FreeGenome(planet[newPos_x][newPos_y].genome, opt);
	    planet[newPos_x][newPos_y].genome = Site1_Genes;
	    planet[tmp_x][tmp_y].genome = NULL;
	    
	    N_fish_t -= 1;
	  }else if(Site2_Species == EMPTY && shark_rw_inibitor == 0){
	    if (DW) cout << " RW" << endl;
	    //shark rw
	    planet[newPos_x][newPos_y].species = SHARK;
	    
	    planet[tmp_x][tmp_y].species = EMPTY;
	    planet[newPos_x][newPos_y].genome = Site1_Genes;
	    
	    planet[tmp_x][tmp_y].genome = NULL;
	    
	    
	    if(prob_d < (p_death_s)){
	      if (DW) cout << " Death killed " << planet[newPos_x][newPos_y].genome <<  endl;
	      //shark death 
	      planet[newPos_x][newPos_y].species = EMPTY;
	      FreeGenome(planet[newPos_x][newPos_y].genome, opt);
	      planet[newPos_x][newPos_y].genome  = NULL;
	      N_shark_t -= 1;
	    };
	  } else if (DW) cout << endl;
	} else if (DW) cout << endl;
      };
    };
    
    N_fish = N_fish_t;
    N_shark= N_shark_t;
  }
  return planet;
}



void savePlanet(char * fname, options opt, specimen** planet) {
  FILE * fp;
  fp = fopen(fname, "w");
  if (!fp) {
    cerr << "ERROR, cannot access to " << fname << endl;
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "# Species; X; Y; Pm; Pf; Pd\n");
  for (int x = 0; x < opt.L; x++) {
    for (int y = 0; y < opt.L; ++y) {
      if (planet[x][y].species != EMPTY) {
	
	fprintf(fp, "%d %d %d %.4f %.4f %.4f\n",
		planet[x][y].species, x, y,
		GetSinglePheno(planet[x][y], 0, opt),
		GetSinglePheno(planet[x][y], 1, opt),
		GetSinglePheno(planet[x][y], 2, opt));
      }
    }
  }
}
