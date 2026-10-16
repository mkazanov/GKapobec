#pragma once

#include <string>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <ostream>
#include <vector>
#include <cstdint>

#include "common.h"

namespace fs = std::filesystem;


// value for undefined statistics (written as NA)
const double NA = std::numeric_limits<double>::quiet_NaN();

//structure to write double as NA if it is undefined (and infinity as Inf)
struct NA_double {
	double value;
};

inline std::ostream& operator<<(std::ostream& out, const NA_double& v) {
	if (std::isnan(v.value))
		return out << "NA";
	if (std::isinf(v.value))
		return out << (v.value > 0 ? "Inf" : "-Inf");
	return out << v.value;
}


//structure for coverage of region (targets)
struct Coverage {
	unsigned int total_len = 0; // without N
	unsigned int c_ds = 0;
	unsigned int tc_ds = 0;
	unsigned int tc_ds_3_end_hairpin = 0;

	Coverage& operator+=(const Coverage& other) {
		total_len += other.total_len;
		c_ds += other.c_ds;
		tc_ds += other.tc_ds;
		tc_ds_3_end_hairpin += other.tc_ds_3_end_hairpin;
		return *this;
	}
};


//structure for mutation counts in region
struct MutCounts {
	unsigned int all = 0;
	unsigned int c_ds = 0;
	unsigned int tc_ds = 0;
	unsigned int tc_ds_apobec = 0;
	unsigned int tc_ds_apobec_hairpin_end = 0;

	MutCounts& operator+=(const MutCounts& other) {
		all += other.all;
		c_ds += other.c_ds;
		tc_ds += other.tc_ds;
		tc_ds_apobec += other.tc_ds_apobec;
		tc_ds_apobec_hairpin_end += other.tc_ds_apobec_hairpin_end;
		return *this;
	}
};


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
	double p0_real_third; // p0_theor_third, or NA if there is no TC ds in window
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
			// if there are no targets in window, t and p0 are NA (NaN)
			if (total_len != 0)
				t_real_first = (double)tc_ds / total_len;
			else
				t_real_first = NA;
			p0_real_first = (t_real_first * 2.0) / 3.0;
			k_first = mutations_tc_ds_apobec;
			n_first = mutations_all;

			if (c_ds != 0)
				t_real_second = (double)tc_ds / c_ds;
			else
				t_real_second = NA;
			p0_real_second = (t_real_second * 2.0) / 3.0;
			k_second = mutations_tc_ds_apobec;
			n_second = mutations_c_ds;

			if (tc_ds != 0)
				p0_real_third = p0_theor_third;
			else
				p0_real_third = NA;
			k_third = mutations_tc_ds_apobec;
			n_third = mutations_tc_ds;

			if (tc_ds != 0)
				t_real_fourth = (double)tc_ds_3_end_hairpin / tc_ds;
			else
				t_real_fourth = NA;
			p0_real_fourth = (t_real_fourth * 2.0) / 3.0;
			k_fourth = mutations_tc_ds_apobec_hairpin_end;
			n_fourth = mutations_tc_ds;
		}

	Stat(const Coverage& coverage, const MutCounts& mutations) :
		Stat(coverage.total_len, coverage.tc_ds, coverage.c_ds, coverage.tc_ds_3_end_hairpin,
			mutations.all, mutations.tc_ds, mutations.tc_ds_apobec, mutations.c_ds, mutations.tc_ds_apobec_hairpin_end)
	{
	}
};


//struct for mutation
struct Mut {
	unsigned int position;
	char ref;
	char alt;
	uint8_t code = 0b1111'1111; // index code of reference position, set in annotate_mutations

	bool operator<(const Mut& other) const { return position < other.position; }
};


//struct for window
struct Window {
	unsigned int start;
	unsigned int end; //include
	Coverage coverage;
};


//struct for analysed region of chromosome (whole chromosome or -s/-e)
struct ChrRegion {
	std::string name;
	std::string level; // "chromosome" or "region"
	unsigned int start;
	unsigned int end; //include
	Coverage coverage;
	std::vector<Window> windows;
	std::vector<std::uint8_t> packed; // index bytes of region (2 positions per byte)

	//function to get index code of position (0-based chromosome coordinate)
	std::uint8_t code_at(unsigned int pos) const {
		std::uint8_t byte = packed[pos / 2 - start / 2];
		// even position - low 4 bits, odd position - high 4 bits
		return (pos % 2 == 0) ? (byte & 0b0000'1111) : (byte >> 4);
	}
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


//functions for statistics

//continued fraction for incomplete beta function (modified Lentz method)
inline double beta_cf(double a, double b, double x) {
	const double tiny = 1e-300;
	const double eps = 1e-15;
	double c = 1.0;
	double d = 1.0 - (a + b) * x / (a + 1.0);
	if (std::fabs(d) < tiny) d = tiny;
	d = 1.0 / d;
	double h = d;
	for (int m = 1; m < 1'000'000; ++m) {
		double m2 = 2.0 * m;
		// even step
		double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
		d = 1.0 + aa * d;
		if (std::fabs(d) < tiny) d = tiny;
		c = 1.0 + aa / c;
		if (std::fabs(c) < tiny) c = tiny;
		d = 1.0 / d;
		h *= d * c;
		// odd step
		aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
		d = 1.0 + aa * d;
		if (std::fabs(d) < tiny) d = tiny;
		c = 1.0 + aa / c;
		if (std::fabs(c) < tiny) c = tiny;
		d = 1.0 / d;
		double delta = d * c;
		h *= delta;
		if (std::fabs(delta - 1.0) < eps)
			break;
	}
	return h;
}

//function to get natural log of regularized incomplete beta function I_x(a, b), a > 0, b > 0
inline double log_beta_inc(double a, double b, double x) {
	if (x <= 0.0) return -std::numeric_limits<double>::infinity();
	if (x >= 1.0) return 0.0;
	double log_front = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log1p(-x);
	if (x < (a + 1.0) / (a + b + 2.0))
		return log_front + std::log(beta_cf(a, b, x)) - std::log(a);
	// I_x(a, b) = 1 - I_(1-x)(b, a), here I_x(a, b) is not small
	double other = std::exp(log_front + std::log(beta_cf(b, a, 1.0 - x)) - std::log(b));
	return std::log1p(-std::min(other, 1.0));
}

//function to get natural log of one-sided binomial p-value P(X >= k), X ~ Bin(n, p0)
inline double log_binom_p_value(double p0, unsigned k, unsigned n) {
	if (std::isnan(p0)) return NA;
	if (k == 0)         return 0.0;
	if (k > n)          return -std::numeric_limits<double>::infinity();
	if (p0 <= 0.0)      return -std::numeric_limits<double>::infinity();
	if (p0 >= 1.0)      return 0.0;
	// P(X >= k) = I_p0(k, n - k + 1)
	return log_beta_inc(k, n - k + 1.0, p0);
}

//function to get quantile of beta distribution (bisection)
inline double beta_quantile(double q, double a, double b) {
	double log_q = std::log(q);
	double lo = 0.0;
	double hi = 1.0;
	for (int i = 0; i < 64; ++i) {
		double mid = (lo + hi) / 2.0;
		if (log_beta_inc(a, b, mid) < log_q)
			lo = mid;
		else
			hi = mid;
	}
	return (lo + hi) / 2.0;
}

//confidence level for confidence intervals
const double CONF_LEVEL = 0.95;

//structure for result of one binomial test
struct TestResult {
	double p_value = NA;
	double neg_log10_p = NA; // -log10(p_value), computed without underflow
	double enrichment = NA;  // (k / n) / p0
	double f = NA;           // fraction of APOBEC mutations: (k / n - p0) / (1 - p0)
	double f_low = NA;       // Clopper-Pearson confidence interval for f
	double f_high = NA;
};

//function to compute binomial test, enrichment and fraction of APOBEC mutations
inline TestResult binom_test(double p0, unsigned k, unsigned n) {
	TestResult res;
	if (std::isnan(p0))
		return res;

	double log_p = log_binom_p_value(p0, k, n);
	res.p_value = std::exp(log_p);
	res.neg_log10_p = -log_p / std::log(10.0) + 0.0; // + 0.0 to avoid -0

	if (n == 0)
		return res;
	double ratio = (double)k / n;
	if (p0 > 0.0)
		res.enrichment = ratio / p0;
	if (p0 < 1.0) {
		double alpha = 1.0 - CONF_LEVEL;
		double low = (k == 0) ? 0.0 : beta_quantile(alpha / 2.0, k, n - k + 1.0);
		double high = (k == n) ? 1.0 : beta_quantile(1.0 - alpha / 2.0, k + 1.0, (double)(n - k));
		res.f = (ratio - p0) / (1.0 - p0);
		res.f_low = (low - p0) / (1.0 - p0);
		res.f_high = (high - p0) / (1.0 - p0);
	}
	return res;
}


//structure for final APOBEC activity measure
struct Activity {
	double value = NA;    // lower estimate of the fraction of all substitutions made by APOBEC, 0..1
	double high = NA;     // upper bound of the final interval
	const char* category = "NA";
};


//function to get activity category by value
inline const char* activity_category(double value) {
	if (std::isnan(value))         return "NA";
	if (value >= ACTIVITY_HIGH)     return "HIGH ACTIVITY";
	if (value >= ACTIVITY_MODERATE) return "MODERATE ACTIVITY";
	if (value >= ACTIVITY_LOW)      return "LOW ACTIVITY";
	return "NO ACTIVITY";
}


//function to compute APOBEC activity: start with test 4 (hairpins, immune to other mutagenesis in TC),
//then use tests 3, 2, 1 while they agree with it; the result never exceeds the hairpin estimate
inline Activity apobec_activity(const Stat& stat, const TestResult tests[4]) {

	Activity act;

	unsigned int n1 = stat.mutations_all;
	unsigned int k = stat.mutations_tc_ds_apobec;
	unsigned int k4 = stat.mutations_tc_ds_apobec_hairpin_end;

	// without targets or without mutations the measure is undefined
	if (stat.tc_ds == 0 || stat.tc_ds_3_end_hairpin == 0 || n1 == 0)
		return act;

	// t - fraction of TC/GA targets at the 3' end of a hairpin loop
	double t = (double)stat.tc_ds_3_end_hairpin / stat.tc_ds;

	// hairpin enrichment among APOBEC mutations must be significant, otherwise APOBEC is not confirmed
	// p-value is compared in log space, so a very small one does not underflow to 0
	if (!(log_binom_p_value(t, k4, k) < std::log(HAIRPIN_P_VALUE_THRESHOLD))) {
		act.value = 0.0;
		act.category = activity_category(act.value);
		return act;
	}

	// pi - fraction of true APOBEC among k mutations: observed enrichment r/t is a mix of
	// HAIRPIN_ENRICHMENT_APOBEC (APOBEC) and 1 (any other process)
	auto purity = [t](double r) {
		double pi = (r / t - 1.0) / (HAIRPIN_ENRICHMENT_APOBEC - 1.0);
		return std::min(std::max(pi, 0.0), 1.0);
	};

	double share = (double)k / n1;
	double r_low = (k4 == 0) ? 0.0 : beta_quantile((1.0 - CONF_LEVEL) / 2.0, k4, k - k4 + 1.0);
	double r_high = (k4 >= k) ? 1.0 : beta_quantile(1.0 - (1.0 - CONF_LEVEL) / 2.0, k4 + 1.0, (double)(k - k4));

	double a4 = purity((double)k4 / k) * share;
	double low = purity(r_low) * share;
	double high = purity(r_high) * share;

	// tests 3, 2, 1: each one estimates the same fraction, but assumes more about the background
	unsigned int n[3] = { stat.n_first, stat.n_second, stat.n_third };
	for (int i = 2; i >= 0; --i) {

		double m = (double)n[i] / n1; // from "fraction among n_i" to "fraction among all substitutions"
		double test_low = tests[i].f_low * m;
		double test_high = tests[i].f_high * m;

		if (std::isnan(test_low) || std::isnan(test_high))
			break;
		if (test_high < low) // test claims less than the current interval allows
			break;
		if (test_low > high) // test claims more: it is contaminated by other mutagenesis
			break;

		low = std::max(low, test_low);
		high = std::min(high, test_high);
	}

	// the measure never exceeds the hairpin estimate, which is the only one immune to other mutagenesis in TC
	act.value = std::min(low, a4);
	act.high = high;
	act.category = activity_category(act.value);
	return act;
}


//function to read packed index of region (region name, start and end must be set)
void read_region(std::ifstream& index, ChrRegion& region, std::unordered_map<std::string, unsigned int>& chr_lens);

//function to define algorythm type;
void define_alg(ProgrammOptions& options);

//function to get array with vcf path
std::vector<fs::path> get_vcf_paths(const ProgrammOptions& options);

//function to get single nucleotide substitutions from ref and one alt allele (MNV are split, indels give nothing)
std::vector<Mut> get_substitutions(const std::string& ref, const std::string& alt);

//function to get mutations from vcf file
void load_vcf(const fs::path& vcf_path, std::vector<std::vector<Mut>>& chr_mutations, const ProgrammOptions& options, const std::unordered_map<std::string, unsigned int>& chr_lens);

//function to count coverage of part of region (0-based chromosome coordinates, end include)
Coverage region_coverage(const ChrRegion& region, unsigned int start, unsigned int end);

//function to split region into windows (last window can be shorter)
std::vector<Window> make_windows(const ChrRegion& region, const ProgrammOptions& options);

//function to check mutations with reference and save their index codes
void annotate_mutations(std::vector<Mut>& mutations, const ChrRegion& region, const fs::path& vcf_path);

//function to count mutations in window (mutations must be sorted and annotated)
MutCounts count_mutations(const std::vector<Mut>& mutations, unsigned int window_start, unsigned int window_end);

//function to write header of output file
void write_header(std::ostream& output_file);

//function to wtire line in output file (start and end are 1-based)
void write_line(std::ostream& output_file, const Stat& w_stat, const std::string& sample, const std::string& level, const std::string& chromosome, std::uint64_t start, std::uint64_t end);

//slow algorythm
void slow_alg(ProgrammOptions& options);

//function to process VCF
void run(ProgrammOptions& options);

