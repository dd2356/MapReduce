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


void shuffle(std::unordered_map<Word,long> &map, int size, 
	int *out_counts, int *out_offsets, Pair *out_data) {

	memset(out_counts, 0, size * sizeof(*out_counts));
    std::hash<Word> word_hasher;
	for (auto& it: map) {
		Word w = it.first;
		// long count = it.second;
		size_t target_process = word_hasher(w) % size;
	    // printf("%s: %ld, %lu\n", w.word, count, hash);
	    out_counts[target_process]++;
	}
	// return;
	out_offsets[0] = 0;
	for (int i = 1; i < size; i++) {
		out_offsets[i] = out_offsets[i-1] + out_counts[i-1];
	}

	int *temp_counts = (int*) calloc(size, sizeof(int));
	for (auto& it: map) {
		Word w = it.first;
		long count = it.second;
		size_t target_process = word_hasher(w) % size;
		Pair p;
		memcpy(p.word, w.word, WORD_SIZE);
		p.count = count;
	    int index = out_offsets[target_process] + temp_counts[target_process];
	    // printf("index: %d / %lu (%d, %d)\n", index, map.size(), out_counts[target_process], temp_counts[target_process]);
	    temp_counts[target_process]++;
	    out_data[index] = p;
	}
}

void reduce(Pair *data, std::unordered_map<Word,long> *out_map) {

}

