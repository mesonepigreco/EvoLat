#include "parse_input.hpp"
#ifndef WATOR_D
#define WATOR_D
specimen** waTor(options opt, int steps, int N_fish0, int N_shark0,   specimen** init_config, bool initialize = true, bool progress_bar = false);

void savePlanet(char * fname, options opt, specimen** planet);
void SaveGenome(char * fname, options opt, specimen ** planet) ;
#endif
