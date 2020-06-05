/*
 * This source reads the input file
 */

#include "parse_input.hpp"
#include <stdlib.h>
#include <libconfig.h++>

#define GEN_LEN 3

using namespace std;
using namespace libconfig;

options read_input(string filename) {
  // Load the configuration file
  Config cfg;

  // Load the file
  try {
    cfg.readFile(filename.c_str());
  } catch (const FileIOException &fioex) {
    cerr << "I/O error while reading " << filename << endl;
    exit(EXIT_FAILURE);
  } catch (const ParseException &pex) {
    cerr << "Error while parsing the input file:" << endl;
    cerr << "File: " << pex.getFile() << " => Line: " << pex.getLine() << endl;
    cerr << pex.getError() << endl;
    exit(EXIT_FAILURE);
  }
  options opt;

  // Get all the value from input
  // The following variables are mandatory
  opt.n = GEN_LEN;
  try {
    opt.L = cfg.lookup("L");
    opt.N = cfg.lookup("N");
    opt.rho_f = cfg.lookup("rho_f");
    opt.rho_s = cfg.lookup("rho_s");
    opt.N_steps = cfg.lookup("N_steps");
    opt.N_sim = cfg.lookup("N_sim");
    opt.N_step_between = cfg.lookup("N_step_between");

    string s = cfg.lookup("Dir_Name");
    opt.Dir_Name = s;
  } catch (const SettingNotFoundException &errorfound) {
    cerr << "Error, the setting " << errorfound.getPath() << " is not present in the input." << endl;
    exit(EXIT_FAILURE);
  } catch (const SettingException &errortype){
    cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl;
    exit(EXIT_FAILURE);
  }

  // Now read non mandatory options, and set them to their default value
  opt.EVOLUTION_TYPE = "uniform";
  try {
    string s = cfg.lookup("EVOLUTION_TYPE");
    opt.EVOLUTION_TYPE = s;
    if (! (opt.EVOLUTION_TYPE.compare("uniform") ||
	   opt.EVOLUTION_TYPE.compare("exp") ||
	   opt.EVOLUTION_TYPE.compare("power"))) {
      cerr << "ERROR, the EVOLUTION_TYPE can be one between 'uniform', 'exp' or 'power'." << endl;
      exit(EXIT_FAILURE);
    }
  } catch (const SettingNotFoundException &errorfound) {
  } catch (const SettingException &errortype){
    cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl;
    exit(EXIT_FAILURE);
  }

  opt.EVOLUTION_EXP_CUTOFF = 20;
  try {
    double cut = cfg.lookup("EVOLUTION_EXP_CUTOFF");
    opt.EVOLUTION_EXP_CUTOFF = cut;
    if (opt.EVOLUTION_TYPE.compare("exp") != 0) {
      cerr << "Error, EVOLUTION_TYPE does not match with EXP_CUTOFF" << endl;
      throw "";
    }
  } catch (const SettingNotFoundException &errorfound) {
  } catch (const SettingException &errortype){
    cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl;
    exit(EXIT_FAILURE);
  }
  
  opt.N_mutations_fishes = 0;
  opt.N_mutations_sharks = 0;
  try {
    if (cfg.exists("N_mutations")) {
      opt.N_mutations_fishes = cfg.lookup("N_mutations");
      opt.N_mutations_sharks = cfg.lookup("N_mutations");
    }

    if (cfg.exists("N_mut_sharks")) {
      opt.N_mutations_sharks = cfg.lookup("N_mut_sharks");
    }

    if (cfg.exists("N_mut_fishes")) {
      opt.N_mutations_fishes = cfg.lookup("N_mut_fishes");
    }
  } catch (const SettingNotFoundException &errorfound) {
  } catch (const SettingException &errortype){
    cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl;
    exit(EXIT_FAILURE);
  }

  
  opt.EVOLUTION_POWER_EXP = -2;
  try {
    double pe =  cfg.lookup("EVOLUTION_POWER_EXP");
    opt.EVOLUTION_POWER_EXP =pe;
    opt.genome_norm = 0;
    for (int i = 0; i < opt.N; ++i) opt.genome_norm += pow(i+1, pe);

    if (opt.EVOLUTION_TYPE.compare("power") != 0) {
      cerr << "Error, EVOLUTION_TYPE does not match with POWER_EXP" << endl;
      throw "";
    }
  } catch (const SettingNotFoundException &errorfound) {
  } catch (const SettingException &errortype){
    cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl;
    exit(EXIT_FAILURE);
  }

  // Entropic forces initialization
  opt.m_fm = 0.5;try {opt.m_fm = cfg.lookup("EntropicForces.m_fm"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.m_ff = 0.5;try {opt.m_ff = cfg.lookup("EntropicForces.m_ff"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.m_sf = 0.5;try {opt.m_sf = cfg.lookup("EntropicForces.m_sf"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.m_sm = 0.5;try {opt.m_sm = cfg.lookup("EntropicForces.m_sm"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.m_sd = 0.5;try {opt.m_sd = cfg.lookup("EntropicForces.m_sd"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}

  // Initial phenotypes initialization
  opt.p_ff = 0.3;try {opt.p_ff = cfg.lookup("InitialPhenotypes.p_ff"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.p_fm = 0.7;try {opt.p_fm = cfg.lookup("InitialPhenotypes.p_fm"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.p_sm = 0.4;try {opt.p_sm = cfg.lookup("InitialPhenotypes.p_sm"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.p_sf = 0.5;try {opt.p_sf = cfg.lookup("InitialPhenotypes.p_sf"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}
  opt.p_sd = 0.1;try {opt.p_sd = cfg.lookup("InitialPhenotypes.p_sd"); } catch (const SettingNotFoundException &errorfound) {}catch (const SettingException &errortype){ cerr << "Error, the setting " << errortype.getPath() << " is of the wrong type." << endl; exit(EXIT_FAILURE);}

  // Check consistency
  if (opt.p_ff + opt.p_fm > 1) {
    cerr << "Error, p_ff + p_fm is greater than" << endl;
    exit(EXIT_FAILURE);
  }
  return opt;
}


void PrintOptions(options opt) {
  // This function prints all the options in standard output.
  // In this way one can check if everithing has been readed correctly

  cout << std::fixed;
  cout << " ===== INPUT OPTIONS ===== " << endl;
  cout << " SYSTEM SIZE: " << opt.L << endl;
  cout << " GENOME SIZE: " << opt.N << endl;
  cout << " PHENOTYPE NUMBERS: " << opt.n << endl;
  cout << " INITIAL FISH DENSITY: " << opt.rho_f << endl;
  cout << " INITIAL SHARK DENSITY: " << opt.rho_f << endl;
  cout << " THERMALIZATION STEPS: " << opt.N_steps << endl;
  cout << " TOTAL NUMBER OF CONFIGURATIONS: " << opt.N_sim << endl;
  cout << " NUMBER OF STEPS BETWEEN TWO CONFIG: " << opt.N_step_between << endl;
  cout << endl << std::scientific;
  cout << " FISH MUTATION RATE: " << opt.N_mutations_fishes << endl;
  cout << " SHARK MUTATION RATE: " << opt.N_mutations_sharks << endl;
  cout << " EVOLUTION TYPE: " << opt.EVOLUTION_TYPE << endl;
  if (opt.EVOLUTION_TYPE == "exp") 
    cout << "   EXPONENT CUTOFF: " << opt.EVOLUTION_EXP_CUTOFF << endl;
  else if (opt.EVOLUTION_TYPE == "power")
    cout << "   POWER EXPONENT: " << opt.EVOLUTION_POWER_EXP << endl;
  
  cout << endl;
  cout << " INTIIAL CONDITIONS (p_f, p_m, p_d)" << endl;
  cout << " FISHES: " << opt.p_ff << "\t" << opt.p_fm << "\t" << 0 << endl;
  cout << " SHARKS: " << opt.p_sf << "\t" << opt.p_sm << "\t" << opt.p_sd << endl; 

  cout << " ENTROPIC FORCES (m_f, m_m, m_d)" << endl;
  cout << " FISHES: " << opt.m_ff << "\t" << opt.m_fm << "\t" << 0 << endl;
  cout << " SHARKS: " << opt.m_sf << "\t" << opt.m_sm << "\t" << opt.m_sd << endl; 
  cout << endl;
}
