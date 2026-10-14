#include <unordered_map>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>


#include "run.h"

using namespace std;

//function to define algorythm type
void define_alg(ProgrammOptions& options) {
	options.alg_type = "slow";
}


//function to unpack index
void unpack_index(ifstream& index, const string& chromosome, vector<uint8_t>& charset_chr, unordered_map<string, unsigned int>& chr_lens, const ProgrammOptions& options) {

	short chr_index = get_index(chromosome);
	unsigned int left = 0;

	for (int i = 0; i < chr_index; ++i) {
		string chr;
		if (i == 22)
			chr = "chrX";
		else if (i == 23)
			chr = "chrY";
		else
			chr = "chr" + to_string(i + 1);

		left += (chr_lens[chr] + 1) / 2;
	}

	unsigned int start;
	unsigned int end;

	if (options.use_boundaries) {
		start = options.start_position;
		end = options.end_position;
	}

	else {
		start = 0;
		end = chr_lens[chromosome] - 1;
	}

	if (end < start || end >= chr_lens[chromosome])
		throw runtime_error("invalid chromosome boundaries");

	unsigned int length = end - start + 1;
	charset_chr.resize(length);
	index.seekg(left + start / 2);

	unsigned int pos = 0;

	while (pos < length)
	{
		uint8_t byte;

		if (!index.read(reinterpret_cast<char*>(&byte), 1))
			throw runtime_error("cannot read index file");

		if (start % 2 == 1 && pos == 0) {
			charset_chr[pos++] = byte >> 4;
		}

		else {
			charset_chr[pos++] = byte & 0b00001111;
			if (pos < length) {
				charset_chr[pos++] = byte >> 4;
			}
		}
	}
}


//function to get array with vcf path
vector<fs::path> get_vcf_paths(const ProgrammOptions& options) 
{
	vector<fs::path> vcf_paths;

	if (options.use_vcf_file()) {
		vcf_paths.push_back(options.vcf_file);
	}
	else if (options.use_vcf_dir()) {
		for (const fs::directory_entry& p : fs::directory_iterator(options.vcf_dir)) {
			if (p.is_regular_file() && p.path().extension() == ".vcf")
				vcf_paths.push_back(p.path());
		}
	}
	else if (options.use_vcf_list())
	{
		ifstream list(options.vcf_list);

		if (!list)
			throw runtime_error("cannot open vcf list file");

		string line;
		while (getline(list, line)) {
			if (!line.empty())
				vcf_paths.push_back(fs::path(line));
		}
	}

	return vcf_paths;
}


//function to get mutations from vcf file
void load_vcf(const fs::path& vcf_path, vector<vector<Mut>>& chr_mutations, const ProgrammOptions& options) {
	
	ifstream vcf(vcf_path);

	if (!vcf)
		throw runtime_error("cannot open vcf file: " + vcf_path.string());

	string line;
	while (getline(vcf, line)) {

		if (line.empty() || line[0] == '#')
			continue;

		stringstream ss(line);
		
		string chr;
		unsigned int pos;
		string id;
		char ref;
		char alt;

		ss >> chr >> pos >> id >> ref >> alt;
		pos--;

		if (options.one_chr() && "chr" + chr != options.chromosome)
			continue;

		short chr_index = get_index("chr" + chr);
	
		Mut mutation;

		mutation.position = pos;
		mutation.alt = toupper(static_cast<unsigned char>(alt));
		mutation.ref = toupper(static_cast<unsigned char>(ref));

		chr_mutations[chr_index].push_back(mutation);
	}
}


//function to process one window
Stat window_stat(unsigned int window_start, unsigned int window_end, const ProgrammOptions& options, const vector<uint8_t>& charset_chr, const vector<vector<Mut>>& chr_mutations) {
	
	unsigned int local_start = window_start - options.start_position;
	unsigned int local_end = window_end - options.start_position;
	short chr_index = get_index(options.chromosome);

	//coverage
	unsigned int total_len = window_end - window_start + 1;
	unsigned int tc_ds = 0;
	unsigned int c_ds = 0;
	unsigned int tc_ds_3_end_hairpin = 0;

	for (unsigned int i = local_start; i <= local_end; ++i) {

		uint8_t byte = charset_chr[i];

		char nucl = define_nucleotide(byte);
		uint8_t status = (byte >> 2) & 0b0000'0011;

		// N positions are excluded from window length
		if (nucl == 'N') {
			total_len--;
			continue;
		}

		if (nucl == 'C' || nucl == 'G') {
			c_ds++;
			if (status != 0b0000'0000) {
				tc_ds++;
				if (status == 0b0000'0010) {
					tc_ds_3_end_hairpin++;
				}
			}
		}
	}

	//mutations
	unsigned int mutations_all = 0;
	unsigned int mutations_tc_ds = 0;
	unsigned int mutations_tc_ds_apobec = 0;
	unsigned int mutations_c_ds = 0;
	unsigned int mutations_tc_ds_apobec_hairpin_end = 0;
	
	for (int i = 0; i < chr_mutations[chr_index].size(); ++i) {
		const Mut& mutation = chr_mutations[chr_index][i];

		if (mutation.position < window_start)
			continue;
		if (mutation.position > window_end)
			break;

		unsigned int local_pos = mutation.position - options.start_position;
		uint8_t byte = charset_chr[local_pos];
		uint8_t status = (byte >> 2) & 0b0000'0011;
		char nucl = define_nucleotide(byte);

		if (nucl != mutation.ref)
			throw runtime_error("different nucleotides in vcf and reference genome");


		bool is_apobec_alt = false;
		if (nucl == 'C')
			is_apobec_alt = (mutation.alt == 'T') || (mutation.alt == 'G');
		if (nucl == 'G')
			is_apobec_alt = (mutation.alt == 'A') || (mutation.alt == 'C');
		
		mutations_all++;
		if (nucl == 'C' || nucl == 'G') {

			mutations_c_ds++;

			if (status != 0b0000'0000) {

				mutations_tc_ds++;

				if (is_apobec_alt)
					mutations_tc_ds_apobec++;

				if (status == 0b0000'0010) {

					if (is_apobec_alt)
						mutations_tc_ds_apobec_hairpin_end++;
				}
			}
		}
	}

	Stat w_stat(total_len, tc_ds, c_ds, tc_ds_3_end_hairpin, mutations_all, mutations_tc_ds, mutations_tc_ds_apobec, mutations_c_ds, mutations_tc_ds_apobec_hairpin_end);
	return w_stat;
}


//function to wtire line in output file
void write_line(ofstream& output_file, const Stat& w_stat, unsigned int window_start, unsigned int window_end, const ProgrammOptions& options, const fs::path& vcf_path) {

	double pv1 = binom_p_value(w_stat.p0_real_first, w_stat.k_first, w_stat.n_first);
	double pv2 = binom_p_value(w_stat.p0_real_second, w_stat.k_second, w_stat.n_second);
	double pv3 = binom_p_value(w_stat.p0_theor_third, w_stat.k_third, w_stat.n_third);
	double pv4 = binom_p_value(w_stat.p0_real_fourth, w_stat.k_fourth, w_stat.n_fourth);
	
	output_file
		<< vcf_path.string() << '\t'
		<< options.chromosome << '\t'
		<< (window_start + 1) << '\t'
		<< (window_end + 1) << '\t'

		<< w_stat.p0_real_first << '\t'
		<< w_stat.k_first << '\t'
		<< w_stat.n_first << '\t'
		<< pv1 << '\t'

		<< w_stat.p0_real_second << '\t'
		<< w_stat.k_second << '\t'
		<< w_stat.n_second << '\t'
		<< pv2 << '\t'

		<< w_stat.p0_theor_third << '\t'
		<< w_stat.k_third << '\t'
		<< w_stat.n_third << '\t'
		<< pv3 << '\t'

		<< w_stat.p0_real_fourth << '\t'
		<< w_stat.k_fourth << '\t'
		<< w_stat.n_fourth << '\t'
		<< pv4 << '\n';
}


//funstion for windows
void process_windows(const ProgrammOptions& options, const vector<uint8_t>& charset_chr, const vector<vector<Mut>>& chr_mutations, ofstream& output_file, const fs::path& vcf_path) {
	for (unsigned int window_start = options.start_position; window_start + options.window_len - 1 <= options.end_position; window_start += options.jump) {
		unsigned int window_end = window_start + options.window_len - 1;
		Stat w_stat = window_stat(window_start, window_end, options, charset_chr, chr_mutations);
		write_line(output_file, w_stat, window_start, window_end, options, vcf_path);
	}
}


//slow algorythm
void slow_alg(ProgrammOptions& options) {

	fs::path executable_dir = getExecutableDirectory();
	unordered_map<string, unsigned int> chr_lens;
	get_lengths(chr_lens, executable_dir, options.reference);

	fs::path index_path = executable_dir / "data" / "index.bin";
	ifstream index(index_path, ios::binary);

	if (!index)
		throw runtime_error("first generate index.bin file");

	if (options.one_chr()) 
	{
		short chr_index = get_index(options.chromosome);
		vector<uint8_t> charset_chr;

		if (options.one_chr() && !options.use_boundaries) {
			options.start_position = 0;
			options.end_position = chr_lens[options.chromosome] - 1;
		}

		charset_chr.resize(options.end_position - options.start_position + 1);

		unpack_index(index, options.chromosome, charset_chr, chr_lens, options);

		vector<fs::path> vcf_paths = get_vcf_paths(options);
		ofstream output_file(fs::path(options.output_dir_path) / "results.tsv");
		if (!output_file)
			throw runtime_error("cannot open output file");
		output_file << fixed << setprecision(10);

		for (fs::path vcf_path : vcf_paths)
		{
			vector<vector<Mut>> chr_mutations(24);
			load_vcf(vcf_path, chr_mutations, options);
			sort(chr_mutations[chr_index].begin(), chr_mutations[chr_index].end());

			if (options.use_window()) {
				process_windows(options, charset_chr, chr_mutations, output_file, vcf_path);
			}

			else {
				Stat w_stat = window_stat(options.start_position, options.end_position, options, charset_chr, chr_mutations);
				write_line(output_file, w_stat, options.start_position, options.end_position, options, vcf_path);
			}
		}
	}

	else 
	{
		throw runtime_error("in development");
		//vector<vector<uint8_t>> charset_full(CHR_CNT);

		//for (int i = 0; i < CHR_CNT; ++i) {
		//	string chr;
		//	if (i <= 21)
		//		chr = "chr" + to_string(i + 1);
		//	if (i == 22)
		//		chr = "chrX";
		//	if (i == 23)
		//		chr = "chrY";

		//	charset_full[i].resize(chr_lens[chr]);
		//	unpack_index(index, chr, charset_full[i], chr_lens);
	}
}


//function to process VCF
void run(ProgrammOptions& options) {

	if (options.alg_type == "undefined" || options.alg_type == "fast") 
		define_alg(options);

	if (options.alg_type == "slow") {
		slow_alg(options);
	}
}