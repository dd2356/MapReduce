#include "msr.h"

#define WORD_SIZE 20

struct Word { 
    char word[WORD_SIZE]; 
}; 

struct Pair { 
    long count;
    struct Word word; 
}; 


void map(char *data, std::map<Word,long> *map) { 

}

void shuffle(std::map<Word,long> *map, int size, int *out_offsets, Pair *out_data) {

}

void reduce(Pair *data, std::map<Word,long> *out_map) {

}
