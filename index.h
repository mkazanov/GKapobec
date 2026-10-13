#pragma once

#include <vector>
#include <string>
#include <filesystem>
#include <cstdint>

#include "common.h"

namespace fs = std::filesystem;


// structure for chromosome
struct Chromosome_index {
	std::string name; //chr1 - chr22, chrX, chrY
	unsigned int size = 0;
	unsigned int vector_size;

	std::vector<uint64_t> A;
	std::vector<uint64_t> T;
	std::vector<uint64_t> C;
	std::vector<uint64_t> G;

	std::vector<uint64_t> end3_hairpin;
	std::vector<uint64_t> end5_hairpin;

	std::vector<uint64_t> TC;
	std::vector<uint64_t> GA;
	std::vector<uint64_t> TC_end3_hairpin;
	std::vector<uint64_t> GA_end5_hairpin;

	Chromosome_index() = default;

	Chromosome_index(const std::string& chromosome_name, unsigned int chromosome_size) :
		name(chromosome_name),
		size(chromosome_size),
		vector_size((chromosome_size + 63) / 64),
		A(vector_size),
		T(vector_size),
		C(vector_size),
		G(vector_size),
		end3_hairpin(vector_size),
		end5_hairpin(vector_size),
		TC(vector_size),
		GA(vector_size),
		TC_end3_hairpin(vector_size),
		GA_end5_hairpin(vector_size)
	{
	}
};



//structure for hairpin
struct Hloop {
	unsigned int start0;
	unsigned int end0;
	unsigned short looppos;
	unsigned short looplen;
};



//functuions to work with bits
inline void set_bit(std::vector<uint64_t>& v, unsigned int pos) {

	unsigned int vector_pos = pos / 64;
	unsigned int num_pos = pos % 64;
	uint64_t bit = uint64_t(1) << num_pos;

	v[vector_pos] |= bit;
}


inline void clear_bit(std::vector<uint64_t>& v, unsigned int pos) {

	unsigned int vector_pos = pos / 64;
	unsigned int num_pos = pos % 64;
	uint64_t mask = ~(uint64_t(1) << num_pos);

	v[vector_pos] &= mask;

}


inline bool get_bit(const std::vector<uint64_t>& v, unsigned int pos) {

	unsigned int vector_pos = pos / 64;
	unsigned int num_pos = pos % 64;
	uint64_t shift = v[vector_pos] >> num_pos;

	bool bit = shift & uint64_t(1);
	return bit;
}


inline void bit_and(const std::vector<uint64_t>& v, const std::vector<uint64_t>& u, std::vector<uint64_t>& result) {

	for (int i = 0; i < v.size(); ++i)
		result[i] = v[i] & u[i];
}

inline void bit_and_inplace(const std::vector<uint64_t>& v, std::vector<uint64_t>& result) {

	for (int i = 0; i < v.size(); ++i)
		result[i] &= v[i];
}



// function to shift vector (right)
void shift_bit_right(const std::vector<uint64_t>& original, std::vector<uint64_t>& shifted);

//function to shift vector (left)
void shift_bit_left(const std::vector<uint64_t>& original, std::vector<uint64_t>& shifted);

//function for reading hairpin file
void read_hairpins(std::vector<Chromosome_index>& chromosomes, const fs::path& executable_dir, const std::string& reference);

//function for generating index
void generate_index(const ProgrammOptions& options);
	