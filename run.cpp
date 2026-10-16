#include <unordered_map>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <zlib.h>

#include "run.h"

using namespace std;


//class to read lines from plain or gzip (bgzip) file
//zlib reads not compressed files as is, so one class is used for .vcf and .vcf.gz
class LineReader {
public:
	explicit LineReader(const fs::path& path) {
		file = gzopen(path.string().c_str(), "rb");
		if (file)
			gzbuffer(file, 1 << 20);
	}

	~LineReader() {
		if (file)
			gzclose(file);
	}

	LineReader(const LineReader&) = delete;
	LineReader& operator=(const LineReader&) = delete;

	bool is_open() const { return file != nullptr; }

	//function to read next line (without \n and \r), returns false at the end of file
	bool getline(string& line) {

		line.clear();
		char buffer[4096];

		while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {

			line += buffer;

			if (line.back() == '\n') {
				line.pop_back();
				if (!line.empty() && line.back() == '\r')
					line.pop_back();
				return true;
			}
		}

		// end of file or error (e.g. truncated gzip file)
		int error;
		const char* message = gzerror(file, &error);
		if (error != Z_OK)
			throw runtime_error(string("cannot read file: ") + message);

		return !line.empty();
	}

private:
	gzFile file = nullptr;
};

//function to define algorythm type
void define_alg(ProgrammOptions& options) {
	options.alg_type = "slow";
}


//function to read packed index of region (region name, start and end must be set)
void read_region(ifstream& index, ChrRegion& region, unordered_map<string, unsigned int>& chr_lens) {

	int chr_index = get_index(region.name);
	uint64_t left = 0;

	for (int i = 0; i < chr_index; ++i)
		left += (chr_lens[get_name(i)] + 1) / 2;

	if (region.end < region.start || region.end >= chr_lens[region.name])
		throw runtime_error("invalid chromosome boundaries");

	unsigned int bytes_cnt = region.end / 2 - region.start / 2 + 1;

	region.packed.resize(bytes_cnt);
	index.clear();
	index.seekg(static_cast<streamoff>(left + region.start / 2));

	if (!index.read(reinterpret_cast<char*>(region.packed.data()), bytes_cnt))
		throw runtime_error("cannot read index file");
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
			string file_name = p.path().filename().string();
			bool is_vcf = file_name.size() > 4 && file_name.compare(file_name.size() - 4, 4, ".vcf") == 0;
			bool is_vcf_gz = file_name.size() > 7 && file_name.compare(file_name.size() - 7, 7, ".vcf.gz") == 0;

			if (p.is_regular_file() && (is_vcf || is_vcf_gz))
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


//function to get single nucleotide substitutions from ref and one alt allele
//alleles of equal length are compared position by position: CA -> TA is C -> T, CC -> TT (MNV) is two C -> T
//alleles of different length (indels) and symbolic alleles (<DEL>, breakends) give no substitutions
//positions with non ACGT bases (*, N, ...) are skipped
//position of substitution is offset from vcf POS
vector<Mut> get_substitutions(const string& ref, const string& alt) {

	vector<Mut> substitutions;

	if (ref.size() != alt.size())
		return substitutions;

	if (alt.find_first_of("<>[]") != string::npos)
		return substitutions;

	const string nucleotides = "ACGT";

	for (unsigned int i = 0; i < ref.size(); ++i) {

		if (ref[i] == alt[i])
			continue;

		if (nucleotides.find(ref[i]) == string::npos || nucleotides.find(alt[i]) == string::npos)
			continue;

		Mut substitution;
		substitution.position = i;
		substitution.ref = ref[i];
		substitution.alt = alt[i];
		substitutions.push_back(substitution);
	}

	return substitutions;
}


//function to get mutations from vcf file
void load_vcf(const fs::path& vcf_path, vector<vector<Mut>>& chr_mutations, const ProgrammOptions& options, const unordered_map<string, unsigned int>& chr_lens) {

	// .vcf or .vcf.gz
	LineReader vcf(vcf_path);

	if (!vcf.is_open())
		throw runtime_error("cannot open vcf file: " + vcf_path.string());

	string line;
	unsigned int line_number = 0;
	while (vcf.getline(line)) {

		line_number++;

		if (line.empty() || line[0] == '#')
			continue;

		stringstream ss(line);

		string chr;
		unsigned int pos;
		string id;
		string ref;
		string alts;

		ss >> chr >> pos >> id >> ref >> alts;

		if (!ss || pos == 0)
			throw runtime_error("cannot parse vcf line " + to_string(line_number));

		pos--;

		// chromosome name can be with or without "chr" prefix
		if (chr.rfind("chr", 0) != 0)
			chr = "chr" + chr;

		// skip contigs which are not in reference (MT, GL..., etc.)
		unordered_map<string, unsigned int>::const_iterator chr_len = chr_lens.find(chr);
		if (chr_len == chr_lens.end())
			continue;

		// position out of chromosome (maybe other reference)
		if (pos >= chr_len->second)
			throw runtime_error("position out of chromosome in vcf line " + to_string(line_number));

		if (options.one_chr() && chr != options.chromosome)
			continue;

		// --pass-only: skip variants which did not pass filters (FILTER is not PASS or .)
		if (options.pass_only) {
			string qual;
			string filter;
			ss >> qual >> filter;
			if (!ss)
				throw runtime_error("no FILTER column in vcf line " + to_string(line_number));
			if (filter != "PASS" && filter != ".")
				continue;
		}

		for (char& c : ref)
			c = toupper(static_cast<unsigned char>(c));

		// multiallelic position (ALT = T,G) gives several substitutions, MNV are split, indels are ignored
		stringstream alt_ss(alts);
		string alt;
		while (getline(alt_ss, alt, ',')) {

			for (char& c : alt)
				c = toupper(static_cast<unsigned char>(c));

			for (Mut mutation : get_substitutions(ref, alt)) {

				mutation.position += pos;

				if (mutation.position >= chr_len->second)
					throw runtime_error("position out of chromosome in vcf line " + to_string(line_number));

				chr_mutations[get_index(chr)].push_back(mutation);
			}
		}
	}

	for (vector<Mut>& mutations : chr_mutations)
		sort(mutations.begin(), mutations.end());
}


//function to count coverage of part of region (0-based chromosome coordinates, end include)
Coverage region_coverage(const ChrRegion& region, unsigned int start, unsigned int end) {

	Coverage coverage;

	for (unsigned int i = start; i <= end; ++i) {

		uint8_t byte = region.code_at(i);

		char nucl = define_nucleotide(byte);
		uint8_t status = (byte >> 2) & 0b0000'0011;

		// N positions are excluded from length
		if (nucl == 'N')
			continue;

		coverage.total_len++;

		if (nucl == 'C' || nucl == 'G') {
			coverage.c_ds++;
			if (status != 0b0000'0000) {
				coverage.tc_ds++;
				if (status == 0b0000'0010) {
					coverage.tc_ds_3_end_hairpin++;
				}
			}
		}
	}

	return coverage;
}


//function to split region into windows (last window can be shorter)
vector<Window> make_windows(const ChrRegion& region, const ProgrammOptions& options) {

	vector<Window> windows;

	for (unsigned int window_start = region.start; window_start <= region.end; window_start += options.jump) {

		unsigned int window_end = min(window_start + options.window_len - 1, region.end);

		Window window;
		window.start = window_start;
		window.end = window_end;
		window.coverage = region_coverage(region, window_start, window_end);
		windows.push_back(window);

		if (window_end == region.end)
			break;
	}

	return windows;
}


//function to check mutations with reference and save their index codes
void annotate_mutations(vector<Mut>& mutations, const ChrRegion& region, const fs::path& vcf_path) {

	for (Mut& mutation : mutations) {

		if (mutation.position < region.start || mutation.position > region.end)
			continue;

		uint8_t byte = region.code_at(mutation.position);

		if (define_nucleotide(byte) != mutation.ref)
			throw runtime_error("different nucleotides in vcf and reference genome: " + vcf_path.string() + ", position " + to_string(mutation.position + 1));

		mutation.code = byte;
	}
}


//function to count mutations in window (mutations must be sorted and annotated)
MutCounts count_mutations(const vector<Mut>& mutations, unsigned int window_start, unsigned int window_end) {

	MutCounts counts;

	Mut first;
	first.position = window_start;

	for (vector<Mut>::const_iterator it = lower_bound(mutations.begin(), mutations.end(), first); it != mutations.end() && it->position <= window_end; ++it) {

		const Mut& mutation = *it;
		uint8_t status = (mutation.code >> 2) & 0b0000'0011;
		char nucl = define_nucleotide(mutation.code);

		bool is_apobec_alt = false;
		if (nucl == 'C')
			is_apobec_alt = (mutation.alt == 'T') || (mutation.alt == 'G');
		if (nucl == 'G')
			is_apobec_alt = (mutation.alt == 'A') || (mutation.alt == 'C');

		counts.all++;
		if (nucl == 'C' || nucl == 'G') {

			counts.c_ds++;

			if (status != 0b0000'0000) {

				counts.tc_ds++;

				if (is_apobec_alt)
					counts.tc_ds_apobec++;

				if (status == 0b0000'0010) {

					if (is_apobec_alt)
						counts.tc_ds_apobec_hairpin_end++;
				}
			}
		}
	}

	return counts;
}


//function to write header of output file
void write_header(ostream& output_file) {

	output_file
		<< "sample" << '\t'
		<< "level" << '\t'
		<< "chromosome" << '\t'
		<< "start" << '\t'
		<< "end" << '\t';

	for (int i = 1; i <= 4; ++i) {
		string n = to_string(i);
		output_file
			<< "p0_" << n << '\t' << "k_" << n << '\t' << "n_" << n << '\t' << "p_value_" << n << '\t'
			<< "neg_log10_p_" << n << '\t' << "enrichment_" << n << '\t'
			<< "f_" << n << '\t' << "f_low_" << n << '\t' << "f_high_" << n << '\t';
	}

	output_file
		<< "len_without_N" << '\t'
		<< "c_ds" << '\t'
		<< "tc_ds" << '\t'
		<< "tc_ds_3_end_hairpin" << '\t'
		<< "apobec_activity" << '\t'
		<< "apobec_category" << '\n';
}


//function to wtire line in output file (start and end are 1-based)
void write_line(ostream& output_file, const Stat& w_stat, const string& sample, const string& level, const string& chromosome, uint64_t start, uint64_t end) {

	struct { double p0; unsigned int k; unsigned int n; } tests[4] = {
		{ w_stat.p0_real_first, w_stat.k_first, w_stat.n_first },
		{ w_stat.p0_real_second, w_stat.k_second, w_stat.n_second },
		{ w_stat.p0_real_third, w_stat.k_third, w_stat.n_third },
		{ w_stat.p0_real_fourth, w_stat.k_fourth, w_stat.n_fourth },
	};

	output_file
		<< sample << '\t'
		<< level << '\t'
		<< chromosome << '\t'
		<< start << '\t'
		<< end << '\t';

	TestResult results[4];
	for (int i = 0; i < 4; ++i)
		results[i] = binom_test(tests[i].p0, tests[i].k, tests[i].n);

	for (int i = 0; i < 4; ++i) {
		const TestResult& res = results[i];
		output_file
			<< NA_double{ tests[i].p0 } << '\t'
			<< tests[i].k << '\t'
			<< tests[i].n << '\t'
			<< NA_double{ res.p_value } << '\t'
			<< NA_double{ res.neg_log10_p } << '\t'
			<< NA_double{ res.enrichment } << '\t'
			<< NA_double{ res.f } << '\t'
			<< NA_double{ res.f_low } << '\t'
			<< NA_double{ res.f_high } << '\t';
	}

	Activity activity = apobec_activity(w_stat, results);

	output_file
		// coverage
		<< w_stat.total_len << '\t'
		<< w_stat.c_ds << '\t'
		<< w_stat.tc_ds << '\t'
		<< w_stat.tc_ds_3_end_hairpin << '\t'
		// APOBEC activity
		<< NA_double{ activity.value } << '\t'
		<< activity.category << '\n';
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

	// chromosomes to process: one (-c) or all (genome mode)
	vector<string> chromosomes;
	if (options.one_chr()) {
		if (chr_lens.find(options.chromosome) == chr_lens.end())
			throw runtime_error("unknown chromosome: " + options.chromosome);
		chromosomes.push_back(options.chromosome);
	}
	else {
		for (int i = 0; i < CHR_CNT; ++i) {
			if (chr_lens.find(get_name(i)) != chr_lens.end())
				chromosomes.push_back(get_name(i));
		}
	}

	// packed index and coverage of regions and windows (once for all samples)
	vector<ChrRegion> regions;
	for (const string& chromosome : chromosomes) {

		ChrRegion region;
		region.name = chromosome;

		if (options.use_boundaries) {
			region.level = "region";
			region.start = options.start_position;
			region.end = options.end_position;
		}
		else {
			region.level = "chromosome";
			region.start = 0;
			region.end = chr_lens[chromosome] - 1;
		}

		read_region(index, region, chr_lens);

		region.coverage = region_coverage(region, region.start, region.end);

		if (options.use_window())
			region.windows = make_windows(region, options);

		regions.push_back(std::move(region));
	}

	ofstream output_file(fs::path(options.output_dir_path) / "results.tsv");
	if (!output_file)
		throw runtime_error("cannot open output file");
	output_file << setprecision(8); // significant digits, small p-values are written as 1.2345678e-25

	write_header(output_file);

	// file with skipped samples and reasons
	ofstream skipped_file(fs::path(options.output_dir_path) / "skipped.tsv");
	if (!skipped_file)
		throw runtime_error("cannot open file for skipped samples");
	skipped_file << "sample" << '\t' << "reason" << '\n';

	unsigned int processed_cnt = 0;
	unsigned int skipped_cnt = 0;

	// samples are processed one by one (only one vcf in memory)
	// lines are grouped by sample: windows and total for each chromosome, then genome
	for (const fs::path& vcf_path : get_vcf_paths(options)) {

		string sample_name = vcf_path.string();

		// broken sample (cannot read, cannot parse, other reference) is skipped with warning
		// all checks are before writing, so there are no lines of skipped sample in results
		vector<vector<Mut>> chr_mutations(CHR_CNT);
		try {
			load_vcf(vcf_path, chr_mutations, options, chr_lens);

			for (const ChrRegion& region : regions)
				annotate_mutations(chr_mutations[get_index(region.name)], region, vcf_path);
		}
		catch (const exception& e) {
			cerr << "warning: sample skipped: " << sample_name << ": " << e.what() << endl;
			skipped_file << sample_name << '\t' << e.what() << '\n';
			skipped_cnt++;
			continue;
		}

		Coverage genome_coverage;
		MutCounts genome_mutations;
		uint64_t genome_len = 0;

		for (const ChrRegion& region : regions) {

			const vector<Mut>& mutations = chr_mutations[get_index(region.name)];

			for (const Window& window : region.windows) {
				Stat w_stat(window.coverage, count_mutations(mutations, window.start, window.end));
				write_line(output_file, w_stat, sample_name, "window", region.name, window.start + 1, window.end + 1);
			}

			MutCounts region_mutations = count_mutations(mutations, region.start, region.end);
			Stat r_stat(region.coverage, region_mutations);
			write_line(output_file, r_stat, sample_name, region.level, region.name, region.start + 1, region.end + 1);

			genome_coverage += region.coverage;
			genome_mutations += region_mutations;
			genome_len += region.end - region.start + 1;
		}

		if (!options.one_chr()) {
			Stat g_stat(genome_coverage, genome_mutations);
			write_line(output_file, g_stat, sample_name, "genome", "genome", 1, genome_len);
		}

		processed_cnt++;
	}

	cerr << "samples processed: " << processed_cnt << ", skipped: " << skipped_cnt << endl;

	if (processed_cnt == 0)
		throw runtime_error("no samples were processed");
}


//function to process VCF
void run(ProgrammOptions& options) {

	if (options.alg_type == "undefined" || options.alg_type == "fast")
		define_alg(options);

	if (options.alg_type == "slow") {
		slow_alg(options);
	}
}
