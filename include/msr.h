#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <string.h>

#define WORD_SIZE 20

struct Word { 
	char word[WORD_SIZE];
	bool operator==(const Word &other) const {
		bool equal = true;
		int i = 0;
		while (equal && word[i] != '\0') {
			equal &= (word[i] + 32 * (('A' <= word[i]) & (word[i] <= 'Z'))) 
				== (other.word[i] + 32 * (
					('A' <= other.word[i]) & (other.word[i] <= 'Z')
				)
			);
			i++;
		}
		return equal & (other.word[i] == '\0');
		// return strcmp(word, other.word) == 0;
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

void map(char *data, int size, int overlap, std::unordered_map<Word, long> &out_map);
void shuffle(std::unordered_map<Word,long> &map, int size, 
	int *out_counts, int *out_offsets, Pair *out_data); 
void reduce(Pair *data, int n, std::unordered_map<Word,long> &out_map); 
