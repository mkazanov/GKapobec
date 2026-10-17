#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

#define MAX_CHR_LENGTH 250'000'000
#define CHR_CNT 24
#define DEFAULT_WINDOW_LEN 1'000'000
#define DEFAULT_JUMP 1'000'000

// APOBEC activity measure
// hairpin enrichment of pure APOBEC mutagenesis (estimated on PCAWG samples with strong activity)
#define HAIRPIN_ENRICHMENT_APOBEC 8.0
// hairpin enrichment must be significant at this level, otherwise activity is 0
#define HAIRPIN_P_VALUE_THRESHOLD 0.001
// boundaries of activity categories
#define ACTIVITY_LOW 0.01
#define ACTIVITY_MODERATE 0.05
#define ACTIVITY_HIGH 0.20

namespace fs = std::filesystem;

// structure for program options
struct ProgrammOptions {
	std::string vcf_file;
	std::string vcf_list;
	std::string vcf_dir;
	std::string output_dir_path;
	std::string reference;
	std::string genome_file_path;
	std::string chromosome;
	std::string alg_type = "undefined";
	int window_len = -1;
	int jump = -1;
	int start_position = -1;
	int end_position = -1; //include
	bool use_boundaries = false;
	bool generate = false;
	bool pass_only = false; // use only variants with FILTER = PASS or .

	bool one_chr() const { return !chromosome.empty(); }
	bool use_vcf_file() const { return !vcf_file.empty(); }
	bool use_vcf_list() const { return !vcf_list.empty(); }
	bool use_vcf_dir() const { return !vcf_dir.empty(); }
	bool use_window() const { return window_len != -1; }
};


// function to get executable directory 
fs::path getExecutableDirectory();

//function to get chromosome lengths
void get_lengths(std::unordered_map<std::string, unsigned int>& chromosome_lengths, const fs::path& executable_dir, const std::string& reference);

//function to get chr index by name
int get_index(const std::string& chr_name);

//function to get chr name by index
std::string get_name(int chr_index);
