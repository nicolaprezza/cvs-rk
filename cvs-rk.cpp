#include <rk_function.hpp>
#include <mod_int.hpp>
#include <matrix.hpp>
#include <fstream>
#include "internal/local_search.hpp"

using namespace cvs_rk;

int main(int argc,char** argv) {

	if(argc != 3){

		cout << "Usage:	" << endl;
		cout <<	"  cvs-rk <file.txt> <P>" << endl;
		cout << "where 'file.txt' contains a M x L ASCII matrix and 0<P<1 is " << endl;
		cout <<	"the probability of selecting a column in the initial solution" <<endl;

		exit(0);

	}


	double P = atof(argv[2]);

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
	auto B = vector<bool>(nr_rows,true);

	for(ulint i=0;i<B.size();++i){

		B[i] = ((double)rand()/RAND_MAX) < P;

	}

	A.run_fixed_n(B);

}
