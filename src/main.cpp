#include <iostream>
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "wator.hpp"
#include "parse_input.h"

void PrintUsage(void);
using namespace std;

int main(int argc, char * argv[]) {

  int opt = 0;
  string input_file;


  // Parse the command line
  opt = getopt(argc, argv, "i:h");
  while (opt != -1) {
    switch (opt) {
    case 'i':
      input_file.assing(optarg);
      break;
    case 'h':
      PrintUsage();
      exit(0)
      break;
    default:
      PrintUsage();
      exit(EXIT_FAILURE);
    }
    opt = getopt(argc, argv, "i:h");
  }


  // Read the configuration file
  options opt;
  opt = read_input(optarg);

  specimen ** configuration;

  // Get the initial fishes and sharks
  int N_fish0, N_shark0;

  N_fish0 = opt.rho_f * opt.L * opt.L;
  N_shark0 = opt.rho_s * opt.L * opt.L;

  // Thermalize
  configuration = waTor(opt, opt.N_steps, N_fish0, N_shark0, configuration);

  char fname[256];
  for (int i = 0; i < opt.N_sim; ++i) {

    // Proceed with the simulation
    configuration = waTor(opt, opt.N_step_between, N_fish0, N_shark0, configuration, false, false);
    
    cout << "Extraction step number " << i + 1 << endl;

    // Get the filename
    fname[0] = 0;
    sprintf(fname, "%s/conf%05d.dat", opt.Dir_Name.c_str(), i);

    // Save the configuration
    cout << "Saving in file " << fname << endl;
    savePlanet(fname, opt, configuration);
    
  }
  cout << "Simulation done." << endl;
}
