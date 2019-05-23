#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <string.h>

#define WORD_SIZE 20

struct Word { 
	char word[WORD_SIZE];
	bool operator==(const Word &other) const {
		return strcmp(word, other.word) == 0;
	}
}; 

struct Pair { 
    long count;
    char word[WORD_SIZE]; 
}; 

namespace std {

	template <>
	struct hash<Word>
	{
		std::size_t operator()(const Word& word) const {

			unsigned long h = 5381;
			int c;
			const char *w = word.word;

			while ((c = *w++)) {
				// h = ((h << 5) + h) + c;  hash * 33 + c;
				h = ((h << 5) + h) + (c + 32 * (('A' <= c) & (c <= 'Z')));
			}
			return h;
		}
	};
}

void map(char *data, int size, std::unordered_map<Word, long> &out_map);
void shuffle(std::unordered_map<Word,long> &map, int size, 
	int *out_counts, int *out_offsets, Pair *out_data); 
void reduce(Pair *data, int n, std::unordered_map<Word,long> &out_map); 
