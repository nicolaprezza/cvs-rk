#include <rk_function.hpp>
#include <mod_int.hpp>
#include <matrix.hpp>
#include <fstream>
#include "internal/local_search.hpp"

using namespace cvs_rk;

int main(int argc,char** argv) {

	std::ifstream infile(argv[1]);

	vector<string> rows;

	std::string line;
	while (std::getline(infile, line)) rows.push_back(line);

	auto A = local_search(rows);

}
