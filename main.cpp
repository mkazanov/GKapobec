#include <iostream>

#include "index.h"
#include "run.h"
#include "common.h"


using namespace std;

int main(int argc, char* argv[]) {

    try 
    {
        ProgrammOptions options;
        for (int i = 1; i < argc; ++i) {
            string arg = argv[i];
            if (arg == "-g")
                options.generate = true;
            else if (arg == "-o" && i + 1 < argc)
                options.output_dir_path = argv[++i];
            else if (arg == "-d" && i + 1 < argc)
                options.vcf_dir = argv[++i];
            else if (arg == "-l" && i + 1 < argc)
                options.vcf_list = argv[++i];
            else if (arg == "-v" && i + 1 < argc)
                options.vcf_file = argv[++i];
            else if (arg == "-r" && i + 1 < argc)
                options.reference = argv[++i];
            else if (arg == "-c" && i + 1 < argc)
                options.chromosome = argv[++i];
            else if (arg == "-s" && i + 1 < argc)
                options.start_position = stoi(argv[++i]); //1-based
            else if (arg == "-e" && i + 1 < argc)
                options.end_position = stoi(argv[++i]);
            else if (arg == "-w" && i + 1 < argc)
                options.window_len = stoi(argv[++i]);
            else if (arg == "-j" && i + 1 < argc)
                options.jump = stoi(argv[++i]);
            else if (arg == "-a" && i + 1 < argc)
                options.alg_type = argv[++i];
            else {
                cerr << "Missing argument or unknown option: " << arg << endl;
                return 1;
            }
        }


        if (options.reference.empty()) {
            cerr << "reference must be provided" << endl;
            return 1;
        }

        if (options.generate) {
            generate_index(options);
            return 0;
        }


        if (options.alg_type != "slow" && options.alg_type != "fast" && options.alg_type != "undefined") {
            cerr << "incorrect algorythm type";
            return 1;
        }

        if ((options.window_len < 0 && options.jump > 0) || (options.window_len > 0 && options.jump < 0)) {
            cerr << "both window len and jump len must be provided" << endl;
            return 1;
        }

        if (options.use_vcf_file() + options.use_vcf_list() + options.use_vcf_dir() > 1) {
            cerr << "vcf directory, vcf path list or vcf file are used simultaneously" << endl;
            return 1;
        }
        if (!(options.use_vcf_file() || options.use_vcf_list() || options.use_vcf_dir())) {
            cerr << "vcf folder path, vcf list path or vcf directory path must be provided" << endl;
            return 1;
        }

        if (options.output_dir_path.empty()) {
            cerr << "output directory path must be provided" << endl;
            return 1;
        }

        if (options.start_position != -1 || options.end_position != -1) {

            options.use_boundaries = true;

            if (options.start_position == -1 || options.end_position == -1) {
                cerr << "both start and end positions must be provided" << endl;
                return 1;
            }

            if (options.start_position < 1 || options.end_position < 1) {
                cerr << "positions must be natural numbers" << endl;
                return 1;
            }

            if (options.end_position < options.start_position) {
                cerr << "end position must be more or equal than start position" << endl;
                return 1;
            }

            if (options.chromosome.empty()) {
                cerr << "if a window is used, a chromosome must be specified" << endl;
                return 1;
            }

            options.start_position--;
            options.end_position--;
        }
        run(options);
        return 0;
    }

    catch (const exception& e)
    {
        cerr << "error: " << e.what() << endl;
        return 1;
    }
}
 
