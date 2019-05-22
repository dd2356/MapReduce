#include "msr.h"
#include <string>
#include <iostream>
#include <regex>
#include <unordered_map>

void map(char *data, int size, std::unordered_map<Word,long> &map) { 
	std::regex r("\\b\\w+\\b");
	// std::regex r("<p>.+</p>");
	std::cmatch sm;
	auto words_begin = std::cregex_iterator(data, data + size, r);
	auto words_end = std::cregex_iterator();
	// std::cout << "Found " 
		// << std::distance(words_begin, words_end) 
		// << " words:\n";
	for (std::cregex_iterator i = words_begin; i != words_end; ++i) {
		std::cmatch match = *i;
		// printf("%lu -> %lu\n", match.position(), match.length());
		if (match.length() < WORD_SIZE) {
			Word w;
			memcpy(w.word, data + match.position(), match.length());
			w.word[match.length()] = '\0';
			if (map.count(w) == 0) {
				map[w] = 0;
			}
			map[w] += 1;
			// printf("%s: %ld\n", w.word, map[w]);
		}
	}
}


void shuffle(std::unordered_map<Word,long> *map, int size, int *out_offsets, Pair *out_data) {

}

void reduce(Pair *data, std::unordered_map<Word,long> *out_map) {

}
