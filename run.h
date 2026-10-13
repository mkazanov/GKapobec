#pragma once

#include <string>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "common.h"

namespace fs = std::filesystem;


//structure for statistics
struct Stat {
	// coverage
	unsigned int total_len; 
	unsigned int tc_ds;
	unsigned int c_ds;
	unsigned int tc_ds_3_end_hairpin;
	// mutations
	unsigned int mutations_all;
	unsigned int mutations_tc_ds;
	unsigned int mutations_tc_ds_apobec;
	unsigned int mutations_c_ds;
	unsigned int mutations_tc_ds_apobec_hairpin_end;
	//pO - probability, t - targets coverage. p0 = 2/3 * t
	//first method: TC ds / total len - (only APOBEC mutagenesis)
	const double t_theor_first = 1.0 / 8.0;
	double t_real_first;
	const double p0_theor_first = (t_theor_first * 2) / 3;
	double p0_real_first;
	unsigned int k_first; // mutations in TC ds, c -> t or c -> g
	unsigned int n_first; // all mutations
	//second method: TC ds / C ds - (there is other mutagenesis not in C ds)
	const double t_theor_second = 1.0 / 4.0;
	double t_real_second;
	const double p0_theor_second = (t_theor_second * 2) / 3;
	double p0_real_second;
	unsigned int k_second; // mutations in TC ds, c -> t or c -> g
	unsigned int n_second; //mutations in C ds
	//third method: TC ds, c -> t or c -> g / TC ds
	const double t_theor_third = 1.0;
	const double p0_theor_third = 2.0 / 3.0;
	unsigned int k_third; // mutations in TC ds, c -> t or c -> g
	unsigned int n_third; // mutations in TC ds
	//fourth method: TC ds on 3' end of hairpin loop / TC ds
	double t_real_fourth;
	double p0_real_fourth;
	unsigned int k_fourth; // mutations in TC ds on 3' end of hairpin loop, c -> t or c -> g
	unsigned int n_fourth; // mutations in TC ds

	Stat(unsigned int total_len_, unsigned int tc_ds_, unsigned int c_ds_, unsigned int tc_ds_3_end_hairpin_, unsigned int mutations_all_, unsigned int mutations_tc_ds_, unsigned int mutations_tc_ds_apobec_, unsigned int mutations_c_ds_, unsigned int mutations_tc_ds_apobec_hairpin_end_):
		total_len(total_len_),
		tc_ds(tc_ds_),
		c_ds(c_ds_),
		tc_ds_3_end_hairpin(tc_ds_3_end_hairpin_),

		mutations_all(mutations_all_),
		mutations_tc_ds(mutations_tc_ds_),
		mutations_tc_ds_apobec(mutations_tc_ds_apobec_),
		mutations_c_ds(mutations_c_ds_),
		mutations_tc_ds_apobec_hairpin_end(mutations_tc_ds_apobec_hairpin_end_)

		{
			t_real_first = (double)tc_ds / total_len;
			p0_real_first = (t_real_first * 2.0) / 3.0;
			k_first = mutations_tc_ds_apobec;
			n_first = mutations_all;

			if (c_ds != 0)
				t_real_second = (double)tc_ds / c_ds;
			else
				throw std::runtime_error("there is no c ds in window, it must be too short");
			p0_real_second = (t_real_second * 2.0) / 3.0;
			k_second = mutations_tc_ds_apobec;
			n_second = mutations_c_ds;

			k_third = mutations_tc_ds_apobec;
			n_third = mutations_tc_ds;

			if (tc_ds != 0)
				t_real_fourth = (double)tc_ds_3_end_hairpin / tc_ds;
			else
				throw std::runtime_error("there is no tc ds in window, it must be too short");
			p0_real_fourth = (t_real_fourth * 2.0) / 3.0;
			k_fourth = mutations_tc_ds_apobec_hairpin_end;
			n_fourth = mutations_tc_ds;
		}
};


//struct for mutation
struct Mut {
	unsigned int position;
	char ref;
	char alt;

	bool operator<(const Mut& other) const { return position < other.position; }
};


//function to define nucleotide
inline char define_nucleotide(uint8_t byte) {

	uint8_t nucleotide = byte & 0b0000'0011;
	char nucl = 'N';

	if (byte == 0b0000'1111) {
		return nucl;
	}
	else {
		switch (nucleotide) {
		case 0b0000'0000:
			nucl = 'A';
			break;

		case 0b0000'0001:
			nucl = 'C';
			break;

		case 0b0000'0010:
			nucl = 'G';
			break;
		case 0b0000'0011:
			nucl = 'T';
			break;
		}
	}

	return nucl;
}


//functions for p-value
inline double binom_pmf(double p, unsigned k, unsigned n) {
    if (k > n)    return 0.0;
    if (p <= 0.0) return k == 0 ? 1.0 : 0.0;
    if (p >= 1.0) return k == n ? 1.0 : 0.0;
    return std::exp(std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0)
                    + k * std::log(p) + (n - k) * std::log1p(-p));
}

inline double binom_p_value(double p0, unsigned k, unsigned n) {
    if (k == 0)          return 1.0;  
    if (n == 0 || k > n) return 0.0;
    double s = 0.0;
    for (unsigned i = n; i >= k; --i) {   
        s += binom_pmf(p0, i, n);
        if (i == 0) break;                
    }
    return std::min(s, 1.0);
}


//function to unpack index
void unpack_index(std::ifstream& index, const std::string& chromosome, std::vector<std::uint8_t>& charset_chr, std::unordered_map<std::string, unsigned int>& chr_lens, const ProgrammOptions& options);

//function to define algorythm type;
void define_alg(ProgrammOptions& options);

//function to get array with vcf path
std::vector<fs::path> get_vcf_paths(const ProgrammOptions& options);

//function to get mutations from vcf file
void load_vcf(const fs::path& vcf_path, std::vector<std::vector<Mut>>& chr_mutations, const ProgrammOptions& options);

//function to wtire line in output file
void write_line(std::ofstream& output_file, const Stat& w_stat, unsigned int window_start, unsigned int window_end, const ProgrammOptions& options, const fs::path& vcf_path);

//function to process one window
Stat window_stat(unsigned int window_start, unsigned int window_end, const ProgrammOptions& options, const std::vector<std::uint8_t>& charset_chr, const std::vector<std::vector<Mut>>& chr_mutations);

//funstion for windows
void process_windows(const ProgrammOptions& options, const std::vector<std::uint8_t>& charset_chr, const std::vector<std::vector<Mut>>& chr_mutations, std::ofstream& output_file, const fs::path& vcf_path);

//slow algorythm
void slow_alg(ProgrammOptions& options);

//function to process VCF
void run(ProgrammOptions& options);

