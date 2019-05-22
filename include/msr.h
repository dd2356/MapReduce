#include <stdio.h>
#include <math.h>
#include <map>
#include <stdlib.h>


struct Pair; 
struct Word; 

void map(char* data,std::map<Word,long>* map);
void shuffle(std::map<Word,long>* map, int size, Pair*); 
void reduce(Pair*, std::map<Word,long>* map); 

