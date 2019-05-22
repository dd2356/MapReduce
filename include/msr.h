#include <stdio.h>
#include <math.h>
#include <map>
#include <stdlib.h>


struct Pair; 
struct Word; 

void map(char *data, std::map<Word,long> *out_map);
void shuffle(std::map<Word,long> *map, int size, int *out_offsets, Pair *out_data); 
void reduce(Pair *data, std::map<Word,long> *out_map); 
