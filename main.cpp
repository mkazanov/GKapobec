#include <iostream>
#include <algorithm>
#include <cctype>

#include "index.h"
#include "run.h"
#include "common.h"


using namespace std;


//exception for incorrect option value
struct OptionError : runtime_error {
    using runtime_error::runtime_error;
};


//function to check that argument is natural number (not option)
static bool is_number(const string& s) {
    return !s.empty() && all_of(s.begin(), s.end(), [](unsigned char c) { return isdigit(c); });
}


//function to convert option value to number with understandable error
static int to_number(const string& option, const string& value) {
    if (!is_number(value))
        throw OptionError("option " + option + " requires a non-negative integer, got '" + value + "'");
    try {
        return stoi(value);
    }
    catch (const out_of_range&) {
        throw OptionError("value of option " + option + " is too large: " + value);
    }
}


//function to print help
static void print_usage(ostream& out, const string& prog) {
    out << "GKapobec - APOBEC activity detector\n"
        << "\n"
        << "Usage:\n"
        << "  " << prog << " -g -r hg19\n"
        << "      build genome index data/index.bin (once)\n"
        << "  " << prog << " -r hg19 (-v FILE | -l LIST | -d DIR) -o OUT_DIR [-c CHR [-s START -e END]] [-w [LEN] [-j STEP]] [--pass-only]\n"
        << "      analyse VCF files\n"
        << "\n"
        << "Options:\n"
        << "  -r REF        reference (required), only hg19 is supported\n"
        << "  -g            build index from data/hg19.fa, data/hg19_chromosome_length.csv\n"
        << "                and data/aaw2872-data-s2.txt (data/ must be next to the executable)\n"
        << "  -v FILE       single VCF file (.vcf or .vcf.gz)\n"
        << "  -l LIST       text file with VCF paths, one per line\n"
        << "  -d DIR        directory with VCF files (*.vcf, *.vcf.gz)\n"
        << "                exactly one of -v, -l, -d must be given\n"
        << "  -o OUT_DIR    output directory (required, must exist)\n"
        << "  -c CHR        analyse one chromosome, e.g. chr1, chrX (default: whole genome)\n"
        << "  -s START      region start, 1-based (only with -c, together with -e)\n"
        << "  -e END        region end, 1-based, inclusive (only with -c, together with -s)\n"
        << "  -w [LEN]      sliding window of length LEN\n"
        << "                without LEN: window " << DEFAULT_WINDOW_LEN << ", step " << DEFAULT_JUMP << "\n"
        << "  -j STEP       window step (only with -w), default: LEN\n"
        << "  --pass-only   use only variants with FILTER = PASS or . (default: all variants)\n"
        << "  -h, --help    show this help\n"
        << "\n"
        << "Output: OUT_DIR/results.tsv (statistics), OUT_DIR/skipped.tsv (skipped samples)\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " -r hg19 -d vcfs/ -o results/ --pass-only\n"
        << "  " << prog << " -r hg19 -l samples.txt -o results/ -w\n"
        << "  " << prog << " -r hg19 -v sample.vcf.gz -o results/ -c chr3 -w 100000\n"
        << "  " << prog << " -r hg19 -v sample.vcf -o results/ -c chr1 -s 1000001 -e 2000000 -w 50000 -j 10000\n";
}


int main(int argc, char* argv[]) {

    string prog = fs::path(argv[0]).filename().string();

    //function to print error with hint about help
    auto usage_error = [&prog](const string& message) {
        cerr << prog << ": " << message << "\n"
             << "Try '" << prog << " -h' for more information." << endl;
        return 1;
    };

    if (argc == 1) {
        print_usage(cerr, prog);
        return 1;
    }

    try
    {
        ProgrammOptions options;
        bool default_window = false;
        for (int i = 1; i < argc; ++i) {
            string arg = argv[i];
            bool has_value = i + 1 < argc;
            if (arg == "-h" || arg == "--help") {
                print_usage(cout, prog);
                return 0;
            }
            else if (arg == "-g")
                options.generate = true;
            else if (arg == "--pass-only")
                options.pass_only = true;
            else if (arg == "-o" && has_value)
                options.output_dir_path = argv[++i];
            else if (arg == "-d" && has_value)
                options.vcf_dir = argv[++i];
            else if (arg == "-l" && has_value)
                options.vcf_list = argv[++i];
            else if (arg == "-v" && has_value)
                options.vcf_file = argv[++i];
            else if (arg == "-r" && has_value)
                options.reference = argv[++i];
            else if (arg == "-c" && has_value)
                options.chromosome = argv[++i];
            else if (arg == "-s" && has_value)
                options.start_position = to_number(arg, argv[++i]); //1-based
            else if (arg == "-e" && has_value)
                options.end_position = to_number(arg, argv[++i]);
            else if (arg == "-w") {
                // -w without value: default window and jump
                if (has_value && is_number(argv[i + 1])) {
                    options.window_len = to_number(arg, argv[++i]);
                }
                else {
                    options.window_len = DEFAULT_WINDOW_LEN;
                    default_window = true;
                }
            }
            else if (arg == "-j" && has_value)
                options.jump = to_number(arg, argv[++i]);
            else if (arg == "-a" && has_value)
                options.alg_type = argv[++i];
            else if (arg == "-o" || arg == "-d" || arg == "-l" || arg == "-v" || arg == "-r" ||
                     arg == "-c" || arg == "-s" || arg == "-e" || arg == "-j" || arg == "-a")
                return usage_error("option " + arg + " requires a value");
            else
                return usage_error("unknown option: " + arg);
        }


        if (options.reference.empty())
            return usage_error("reference must be provided (-r hg19)");

        if (options.generate) {
            generate_index(options);
            return 0;
        }


        if (options.alg_type != "slow" && options.alg_type != "fast" && options.alg_type != "undefined")
            return usage_error("incorrect algorythm type");

        if (!options.use_window() && options.jump != -1)
            return usage_error("jump len (-j) can be used only with window (-w)");

        // -w without -j: default jump for default window, otherwise jump = window len
        if (options.use_window() && options.jump == -1)
            options.jump = default_window ? DEFAULT_JUMP : options.window_len;

        if (options.use_window() && (options.window_len <= 0 || options.jump <= 0))
            return usage_error("window len and jump len must be positive");

        if (options.use_vcf_file() + options.use_vcf_list() + options.use_vcf_dir() > 1)
            return usage_error("only one of -v, -l, -d can be used");
        if (!(options.use_vcf_file() || options.use_vcf_list() || options.use_vcf_dir()))
            return usage_error("VCF input must be provided (-v FILE, -l LIST or -d DIR)");

        if (options.output_dir_path.empty())
            return usage_error("output directory must be provided (-o OUT_DIR)");

        if (options.start_position != -1 || options.end_position != -1) {

            options.use_boundaries = true;

            if (options.start_position == -1 || options.end_position == -1)
                return usage_error("both start (-s) and end (-e) positions must be provided");

            if (options.start_position < 1 || options.end_position < 1)
                return usage_error("positions must be natural numbers");

            if (options.end_position < options.start_position)
                return usage_error("end position must be more or equal than start position");

            if (options.chromosome.empty())
                return usage_error("if start and end positions are used, a chromosome must be specified (-c)");

            options.start_position--;
            options.end_position--;
        }
        run(options);
        return 0;
    }

    catch (const OptionError& e)
    {
        return usage_error(e.what());
    }
    catch (const exception& e)
    {
        cerr << "error: " << e.what() << endl;
        return 1;
    }
}
