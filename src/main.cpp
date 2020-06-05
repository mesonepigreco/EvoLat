#include <iostream>
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "waTor.hpp"
#include "parse_input.hpp"
#include "genome.hpp"
#include <getopt.h>
#include <time.h>
#include <fstream>

#define DEBUG_M 0

void PrintUsage(void);
using namespace std;

int main(int argc, char * argv[]) {

  int cmd = 0;
  string input_file;

  // Initialize the seed for the random generator
  unsigned int random_seed;
  ifstream file ("/dev/urandom", ios::binary);
  if (file.is_open()) {
    char * memblock;
    int size = sizeof(int);
    memblock = new char [size];
    file.read (memblock, size);
    file.close();
    random_seed = *reinterpret_cast<int*>(memblock);
    delete[] memblock;
  } else{ random_seed = time(0);}

  srand(random_seed);

  // Parse the command line
  cmd = getopt(argc, argv, "i:h");
  while (cmd != -1) {
    switch (cmd) {
    case 'i':
      input_file = optarg;
      break;
    case 'h':
      PrintUsage();
      exit(0);
      break;
    default:
      PrintUsage();
      exit(EXIT_FAILURE);
      break;
    }
    cmd = getopt(argc, argv, "i:h");
  }


  // Read the configuration file
  options opt;
  opt = read_input(input_file);

  specimen ** configuration;

  // Get the initial fishes and sharks
  int N_fish0, N_shark0;

  N_fish0 = opt.rho_f * opt.L * opt.L;
  N_shark0 = opt.rho_s * opt.L * opt.L;

  cout << "Starting the termalization..." << endl;
  
  // Thermalize
  configuration = waTor(opt, opt.N_steps, N_fish0, N_shark0, configuration);

  cout << "End of the thermalization ..." << endl;
  
  char fname[256];
  for (int i = 0; i < opt.N_sim; ++i) {

    if (DEBUG_M) cout << "Cycle " << i << endl;
    // Proceed with the simulation
    configuration = waTor(opt, opt.N_step_between, N_fish0, N_shark0, configuration, false, false);
    
    cout << "Extraction step number " << i + 1 << endl;

    // Get the filename
    fname[0] = 0;
    sprintf(fname, "%s/conf%05d.dat", opt.Dir_Name.c_str(), i);

    // Save the configuration
    cout << "Saving planet in file " << fname << endl;
    savePlanet(fname, opt, configuration);
    
    // Save the average genome
    fname[0] = 0;
    sprintf(fname, "%s/genom%05d.dat", opt.Dir_Name.c_str(), i);
    cout << "Saving genome in file " << fname << endl;
    SaveGenome(fname, opt, configuration);
  }

  cout << "Destruction" << endl;
  DestroyAll(opt, configuration);
  cout << "Simulation done." << endl;
}



void PrintUsage(void) {
  cout << "Alleluia" << endl;
}
