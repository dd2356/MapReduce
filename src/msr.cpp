#include "msr.h"
#include <string>
#include <iostream>
#include <regex>
#include <unordered_map>
#include <omp.h>


void map(char *data, int size, std::unordered_map<Word,long> &map) { 
	std::vector<int> newlines;
	std::vector<int> paragraph_starts;
	std::vector<int> paragraph_ends;

	for (int i = 0; i < size; i++) {
		if (data[i] == '\n') {
			newlines.push_back(i);
		}
	}
	for (unsigned int i = 0; i < newlines.size(); i++) {
		int ind = newlines[i]+1;
		bool is_paragraph = (data[ind] == '<') 
			&& (data[ind+1] == 'p') && (data[ind+2] == '>');
			// printf("paragraph_starts: %c%c%c, %d\n", data[ind+1], data, is_paragraph);
		if (i < newlines.size() - 1 && is_paragraph) {
			int next_ind = newlines[i+1];
			paragraph_starts.push_back(ind);
			paragraph_ends.push_back(next_ind);
		}
	}

	// std::regex r("\\b\\w+\\b");
	// std::cmatch sm;
	Word w;
	bool in_tag = false;
	bool valid_word;
	int temp_ind = 0;

	for (unsigned int i = 0; i < paragraph_starts.size(); i++) {
		int s_ind = paragraph_starts[i]+2;
		int e_ind = paragraph_ends[i];
		for (int i = s_ind; i < e_ind; i++) {
			if (data[i] == '<') {
				in_tag = true;
				continue;
			} else if (data[i] == '>') {
				in_tag = false;
				continue;
			}
			if (in_tag) {
				continue;
			} else {
				// printf("looking for word\n");
				temp_ind = 0;
				valid_word = true;
				while (data[i+temp_ind] != ' ' && data[i+temp_ind] != '<' 
					&& i + temp_ind < e_ind && temp_ind < WORD_SIZE) {
					if ((data[i+temp_ind] < 'a' || data[i+temp_ind] > 'z') 
						&& (data[i+temp_ind] < 'A' || data[i+temp_ind] > 'Z')) {
						valid_word = false;
					}
					temp_ind++;
				}
				if (data[i+temp_ind] == '<') {
					// printf("found tag\n");
					i += temp_ind - 1;
					continue;
				} else if (valid_word && temp_ind < WORD_SIZE && temp_ind > 1) {
					memcpy(w.word, &data[i], temp_ind);
					w.word[temp_ind] = '\0';
					// printf("(%s)\n", w.word);
					if (map.count(w) == 0) {
						map[w] = 0;
					}
					map[w] += 1;
					i += temp_ind;
				} else {
					// printf("word too long\n");
				}
			}
		}

	}

}

void map2(char *data, int size, std::unordered_map<Word,long> &map) { 

	std::regex r("\\b\\w+\\b");
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

void reduce(Pair *data, int n, std::unordered_map<Word,long> &out_map) {
	for (int i = 0; i < n; i++) {
		// printf("test\n");
		Word w;
		// printf("test: %s\n", data[i].word);
		memcpy(w.word, data[i].word, WORD_SIZE);
		// printf("test\n");
		// w.word = data[i].word;
		int count = data[i].count;
		if (out_map.count(w) == 0) {
			out_map[w] = 0;
		}
		out_map[w] += count;
	}
}

