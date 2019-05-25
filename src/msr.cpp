#include "msr.h"
#include <string>
#include <iostream>
#include <regex>
#include <unordered_map>
#include <omp.h>

// #define DEBUG 

void find_newlines_and_paragraphs(char *data, int size, int overlap,
	std::vector<int> &newlines, std::vector<int> &paragraph_starts, 
	std::vector<int> &paragraph_ends) {

	bool found_new_line = false;
	// unsigned int all_paragraphs = 0;

	// find all lines in original text
	for (int i = 0; i < size; i++) {
		if (data[i] == '\n') {
			newlines.push_back(i);
		}
	}

	// find first newline in overlap region
	for (int i = 0; i < overlap; i++) {
		if (data[size + i] == '\n') {
			newlines.push_back(size + i);
			found_new_line = true;
			// printf("found next newline at %d\n", i);
			break;
		}
	}

	// printf("Start of map: %ld\n", newlines.size());
	for (int i = 0; i < (int)newlines.size() - 1; i++) {
		int ind = newlines[i]+1;
		while (data[ind] == '\t' || data[ind] == ' ') {
			ind++;
		}

		bool is_paragraph = (data[ind] == '<') 
			&& (data[ind+1] == 'p') && (data[ind+2] == '>');

		if (is_paragraph) {
			int next_ind = newlines[i+1];
			paragraph_starts.push_back(ind);
			paragraph_ends.push_back(next_ind);
			// if (i == (int)newlines.size() - 2 && found_new_line) {
				// printf("paragraph found at border between %d and %d\n", ind, next_ind);
				// for (int j = ind; j < next_ind; j++) {
					// printf("%c", data[j]);
				// }
			// }
		}
		// for (int j = newlines[i]+1; j < newlines[i+1]; j++) {
			// if ((data[j] == '<') 
			// && (data[j+1] == 'p') && (data[j+2] == '>')) {
				// all_paragraphs++;
				// break;
			// }
		// }
	}
	// if (paragraph_starts.size() != all_paragraphs) {
		// printf("paragraphs: %ld/%d\n", paragraph_starts.size(), all_paragraphs);
	// }
}

inline bool is_separator(char a) {
	return a == ' ' || a == ',' || a == '.' || a == '!' 
		|| a == '?' || a == ':' || a == ';';
}

inline bool is_breakpoint(char a) {
	return is_separator(a) || a == '<';
}

inline bool is_letter(char a) {
	return (a >= 'a' && a <= 'z') || (a >= 'A' && a <= 'Z') 
		|| (a == '\'' || a == '-');
}

// Debug
// #ifdef DEBUG
#define NON_VALID_WORD 2
#define TOO_LONG   3
#define TOO_SHORT 4
#define SUCCESS 5
// #endif
// Debug end

int try_get_word(char *data, std::unordered_map<Word,long> &map, 
	int s_ind, int e_ind, int &i, Word &w, long &length_counter) {
	
	int temp_length = 0;
	bool valid_word = true;

	while (!is_breakpoint(data[i+temp_length])
		&& i + temp_length < e_ind && temp_length < WORD_SIZE) {
		if (!is_letter(data[i+temp_length])) {
			valid_word = false;
		}
		temp_length++;
	}
	if (valid_word && temp_length < WORD_SIZE && temp_length > 1) {
		memcpy(w.word, &data[i], temp_length);
		w.word[temp_length] = '\0';
		map[w]++;
		length_counter += temp_length - 1;
		if (data[i+temp_length] == '<') {
			temp_length--;
		}
		i += temp_length;
		return SUCCESS;
	}
	return valid_word ? (temp_length <= 1 ? TOO_SHORT : TOO_LONG) 
		: NON_VALID_WORD;
}

void map(char *data, int size, int overlap, std::unordered_map<Word,long> &map) { 
	std::vector<int> newlines, paragraph_starts, paragraph_ends;
	find_newlines_and_paragraphs(data, size, overlap, 
		newlines, paragraph_starts, paragraph_ends);
	long total_word_length = 0; // also debug variable
#ifdef DEBUG
	long total_text_length = 0;
	int non_valid_words = 0;
	int success_words = 0;
	int too_long = 0;
	int too_short = 0;
	for (unsigned int i = 0; i < paragraph_starts.size(); i++) {
		total_text_length += paragraph_ends[i] - paragraph_starts[i];
	}
#endif

	Word w;
	bool in_tag = false;

	for (unsigned int i = 0; i < paragraph_starts.size(); i++) {
		int s_ind = paragraph_starts[i]+2;
		int e_ind = paragraph_ends[i];
		for (int i = s_ind; i < e_ind; i++) {
			if (data[i] == '<') {
				in_tag = true;
			} else if (data[i] == '>') {
				in_tag = false;
				continue;
			} else if (!is_letter(data[i])) {
				continue;
			}
			if (!in_tag) {
#ifdef DEBUG
				int success = try_get_word(data, map, s_ind, e_ind, 
					i, w, total_word_length);
				non_valid_words += success == NON_VALID_WORD;
				too_long += success == TOO_LONG;
				too_short += success == TOO_SHORT;
				success_words += success == SUCCESS;
#else
				try_get_word(data, map, s_ind, e_ind, i, w, total_word_length);				
#endif
			}
		}
	}
#ifdef DEBUG
	int total_words = non_valid_words + too_long + too_short + success_words;

	printf("newlines: %ld, paragraphs: %ld, total length: %ld, total word "
		"length: %ld, nv: %.2f, tl: %.2f, ts: %.2f, success: %.2f\n", 
		newlines.size(), paragraph_starts.size(), 
		total_text_length, total_word_length, 
		non_valid_words / (double)(total_words),
		too_long / (double)(total_words),
		too_short / (double)(total_words),
		success_words / (double)(total_words)
	);
#endif
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
			map[w]++;
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
		out_map[w] += count;
	}
}

