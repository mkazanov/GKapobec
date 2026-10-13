#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

#define MAX_CHR_LENGTH 250'000'000
#define CHR_CNT 24

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
