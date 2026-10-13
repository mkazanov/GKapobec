#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <cctype>


#include "index.h"

#define CHR_CNT 24
#define HAIRPIN_FILE_H19 "aaw2872-data-s2.txt"
#define FASTA ".fa"
#define BUFFER_SIZE (1 << 20)

using namespace std;
namespace fs = filesystem;



//function to shift vector<uint64_t> (right)
void shift_bit_right(const vector<uint64_t>& original, vector<uint64_t>& shifted) {

	for (size_t i = original.size() - 1; i > 0; --i) 
	{
		uint64_t last_bit = original[i - 1] >> 63;
		shifted[i] = (original[i] << 1) | last_bit;
	}
	shifted[0] = original[0] << 1;
}



//function to shift vector<uint64_t> (left)
void shift_bit_left(const vector<uint64_t>& original, vector<uint64_t>& shifted) {

	size_t n = original.size();
	for (int i = 0; i < n - 1; ++i) 
	{
		uint64_t first_bit = original[i + 1] << 63;
		shifted[i] = (original[i] >> 1) | first_bit;
	}
	shifted[n - 1] = original[n - 1] >> 1;
}



//function for reading hairpin file
void read_hairpins(vector<Chromosome_index>& chromosomes, const fs::path& executable_dir, const string& reference) {
	
	fs::path hairpin_path;
	if (reference == "hg19")
		hairpin_path = executable_dir / "data" / HAIRPIN_FILE_H19;
	else
		throw runtime_error("unsupported reference");

	ifstream h(hairpin_path);
	if (!h)
		throw runtime_error("cannot open hairpin file: " + hairpin_path.string());

	string line;
	getline(h, line);

	while (getline(h, line)) {

		stringstream ss(line);

		string tmp;
		unsigned short chr_index;
		unsigned int pos;
		unsigned short looplen;
		unsigned short looppos;

		ss >> tmp;
		ss >> chr_index;
		--chr_index;
		ss >> pos;

		for (int i = 0; i < 9; i++)
			ss >> tmp;

		ss >> looplen;
		ss >> looppos;

		unsigned int pos0 = pos - 1;
		unsigned int chr_len = chromosomes[chr_index].size;
		if (pos0 >= chr_len)
			throw runtime_error("incorrect hairpin coordinates");

		bool reverse;

		if (get_bit(chromosomes[chr_index].A, pos0) || get_bit(chromosomes[chr_index].C, pos0))
			reverse = false;
		else if (get_bit(chromosomes[chr_index].G, pos0) || get_bit(chromosomes[chr_index].T, pos0))
			reverse = true;
		else
			continue;

		unsigned int start0;
		unsigned int end0;

		if (!reverse) {

			if (pos0 < looppos - 1)
				throw runtime_error("incorrect hairpin coordinates");

			start0 = pos0 - (looppos - 1);
			end0 = start0 + looplen - 1;

		}
		else {

			if (pos0 < looplen - looppos)
				throw runtime_error("incorrect hairpin coordinates");

			start0 = pos0 - (looplen - looppos);
			end0 = start0 + looplen - 1;
		}

		if (end0 >= chr_len)
			throw runtime_error("incorrect hairpin coordinates");

		set_bit(chromosomes[chr_index].end3_hairpin, end0);
		set_bit(chromosomes[chr_index].end5_hairpin, start0);
	}
}



//function for generating index 
void generate_index(const ProgrammOptions& options) {

	fs::path executable_dir = getExecutableDirectory();
	fs::path reference_path = executable_dir / "data" / (options.reference + FASTA);

	fs::path index_path = executable_dir / "data" / "index.bin";

	if (fs::exists(index_path))
	{
		std::cout << "Index already exists.\n";
		return;
	}

	unordered_map<string, unsigned int> chromosome_lengths;
	get_lengths(chromosome_lengths, executable_dir, options.reference);

	vector<Chromosome_index> chromosomes(CHR_CNT);
	for (unordered_map<string, unsigned int>::iterator it = chromosome_lengths.begin(); it != chromosome_lengths.end(); ++it) {
		string chr_name = it->first;
		unsigned int length = it->second;
		chromosomes[get_index(chr_name)] = Chromosome_index(chr_name, length);
	}

	ifstream genome(reference_path);
	if (!genome)
		throw runtime_error("cannot open reference file: " + reference_path.string());

	string line;
	int current_chr = -1;
	unsigned int current_pos = 0;

	while (getline(genome, line)) 
	{
		if (line.empty()) continue;

		if (line[0] == '>') {

			if (current_chr != -1)
			{
				if (current_pos != chromosomes[current_chr].size)
					throw runtime_error("chromosome length mismatch");
			}

			string chr_name = line.substr(1);

			if (chromosome_lengths.find(chr_name) != chromosome_lengths.end()) {
				current_chr = get_index(chr_name);
				current_pos = 0;
			}
			else {
				current_chr = -1;
			}
			continue;
		}

		if (current_chr == -1)
			continue;

		for (char nucleotide : line) {

			if (current_pos >= chromosomes[current_chr].size)
				throw runtime_error("chromosome length mismatch");

			nucleotide = toupper(static_cast<unsigned char>(nucleotide));

			switch (nucleotide)
			{
			case 'A':
				set_bit(chromosomes[current_chr].A, current_pos);
				break;

			case 'C':
				set_bit(chromosomes[current_chr].C, current_pos);
				break;

			case 'G':
				set_bit(chromosomes[current_chr].G, current_pos);
				break;

			case 'T':
				set_bit(chromosomes[current_chr].T, current_pos);
				break;

			case 'N':
				break;

			default:
				throw runtime_error("Unexpected nucleotide");
			}

			++current_pos;
		}
	}
	
	if (current_chr != -1)
	{
		if (current_pos != chromosomes[current_chr].size)
			throw runtime_error("chromosome length mismatch");
	}


	read_hairpins(chromosomes, executable_dir, options.reference);


	ofstream index(index_path, ios::binary);
	if (!index)
		throw std::runtime_error("cannot create index file: " + index_path.string());

	vector<uint8_t> buffer;
	buffer.reserve(BUFFER_SIZE + 16);

	uint8_t byte;
	uint8_t code;

	for (Chromosome_index& chr : chromosomes) 
	{
		shift_bit_right(chr.T, chr.TC);
		bit_and_inplace(chr.C, chr.TC);

		shift_bit_left(chr.A, chr.GA);
		bit_and_inplace(chr.G, chr.GA);

		bit_and(chr.TC, chr.end3_hairpin, chr.TC_end3_hairpin);
		bit_and(chr.GA, chr.end5_hairpin, chr.GA_end5_hairpin);

		byte = 0;
		for (int vector_pos = 0; vector_pos < chr.vector_size; ++vector_pos) {

			unsigned short boundary = 64;
			if (vector_pos == chr.vector_size - 1) {
				boundary = (chr.size - 1) % 64 + 1;
			}
			
			for (int bit = 0; bit < boundary; ++bit) {

				bool a = chr.A[vector_pos] & 1;
				bool c = chr.C[vector_pos] & 1;
				bool t = chr.T[vector_pos] & 1;
				bool g = chr.G[vector_pos] & 1;

				bool tc = chr.TC[vector_pos] & 1;
				bool ga = chr.GA[vector_pos] & 1;
				bool tc_end = chr.TC_end3_hairpin[vector_pos] & 1;
				bool ga_end = chr.GA_end5_hairpin[vector_pos] & 1;


				if (a)
					code = 0b00'00;
				else if (c)
					code = 0b00'01;
				else if (t)
					code = 0b00'11;
				else if (g)
					code = 0b00'10;
				else
					code = 0b11'11;


				if (c) {
					
					if (tc_end)
						code = 0b10'01;
					else if (tc)
						code = 0b01'01;
				}

				if (g) {

					if (ga_end)
						code = 0b10'10;
					else if (ga)
						code = 0b01'10;
				}



				if (bit % 2 == 0) {
					byte = code;
				}
				else {
					byte |= (code << 4);
					buffer.push_back(byte);
				}

				if (buffer.size() >= BUFFER_SIZE) {

					index.write(reinterpret_cast<const char*> (buffer.data()), static_cast<streamsize> (buffer.size()));
					buffer.clear();

					if (index.fail())
						throw std::runtime_error("error writing index file");
				}


				chr.A[vector_pos] >>= 1;
				chr.C[vector_pos] >>= 1;
				chr.T[vector_pos] >>= 1;
				chr.G[vector_pos] >>= 1;

				chr.GA_end5_hairpin[vector_pos] >>= 1;
				chr.TC_end3_hairpin[vector_pos] >>= 1;

				chr.TC[vector_pos] >>= 1;
				chr.GA[vector_pos] >>= 1;
			}

			if (boundary % 2 == 1)
				buffer.push_back(byte);
		}

		if (buffer.size() != 0) {
			index.write(reinterpret_cast<const char*> (buffer.data()), static_cast<streamsize> (buffer.size()));
			buffer.clear();

			if (index.fail())
				throw std::runtime_error("error writing index file");
		}
	}
}