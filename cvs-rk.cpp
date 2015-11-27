#include <rk_function.hpp>
#include <mod_int.hpp>
#include <matrix.hpp>
#include <fstream>
#include "internal/local_search.hpp"

using namespace cvs_rk;

int main(int argc,char** argv) {

	if(argc != 2){

		cout << "Usage:	" << endl;
		cout <<	"  cvs-rk file.txt" << endl;
		cout << "where 'file.txt' contains a M x L ASCII matrix" <<endl;

		exit(0);

	}

	std::ifstream infile(argv[1]);

	vector<string> rows;

	std::string line;
	int nr_rows=-1;

	while (std::getline(infile, line)){

		assert(nr_rows<0 or nr_rows == line.length());

		rows.push_back(line);
		nr_rows = line.length();

	}

	auto A = local_search(rows);

	srand(time(NULL));
	auto B = vector<bool>(nr_rows);

	for(ulint i=0;i<B.size();++i) B[i] = rand()%2;

	A.run(B);

}
