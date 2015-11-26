#include <rk_function.hpp>
#include <mod_int.hpp>
#include <matrix.hpp>
#include <opt_algorithm.hpp>
#include <fstream>

using namespace cvs_rk;

int main(int argc,char** argv) {

	std::ifstream infile(argv[1]);

	vector<string> rows;

	std::string line;
	while (std::getline(infile, line)) rows.push_back(line);

	auto A = opt_algorithm(rows);
	A.run();

}
