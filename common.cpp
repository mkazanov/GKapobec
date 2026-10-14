#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <stdexcept>

#include "common.h"
	
#ifdef _WIN32
#include <windows.h>
#endif 

#ifdef __linux__
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#define CHR_LEN_FILE "_chromosome_length.csv"


using namespace std;



// function to get executable directory 
fs::path getExecutableDirectory() {

#ifdef _WIN32

	char p[MAX_PATH];
	int len = GetModuleFileNameA(nullptr, p, MAX_PATH);

	if (len == 0)
		throw runtime_error("Cannot get executable path");

	return fs::path(p).parent_path();

#elif __linux__

	char p[4000];
	int len = readlink("/proc/self/exe", p, sizeof(p) - 1);

	if (len == -1)
		throw runtime_error("Cannot get executable path");

	p[len] = '\0';

	return fs::path(p).parent_path();

#elif __APPLE__

	char p[4000];
	uint32_t size = sizeof(p);

	if (_NSGetExecutablePath(p, &size) != 0)
		throw runtime_error("Cannot get executable path");

	return fs::path(p).parent_path();

#else

#error "Unsupported platform"

#endif 
}



//function to get chromosome lengths
void get_lengths(unordered_map<string, unsigned int>& chromosome_lengths, const fs::path& executable_dir, const string& reference) {

	fs::path lengths_path = executable_dir / "data" / (reference + CHR_LEN_FILE);
	ifstream file(lengths_path);

	if (!file)
		throw runtime_error("cannot open file with chromosome length: " + lengths_path.string());

	string line;
	getline(file, line);

	string chromosome;
	unsigned int length;

	while (getline(file, line)) {
		stringstream ss(line);
		getline(ss, chromosome, ',');
		ss >> length;
		chromosome_lengths["chr" + chromosome] = length;
	}
}



//function to get chr index by name
int get_index(const string& chr) { 

	if (chr == "chrX")
		return 22;

	if (chr == "chrY")
		return 23;

	return stoi(chr.substr(3)) - 1;
}



//function to get chr name by index
string get_name(int chr_index) {

	if (chr_index == 22)
		return "chrX";

	if (chr_index == 23)
		return "chrY";

	return "chr" + to_string(chr_index + 1);
}
