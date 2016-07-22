#include <rk_function.hpp>
#include <mod_int.hpp>
#include <matrix.hpp>
#include <fstream>
#include "internal/local_search.hpp"

using namespace cvs_rk;

string infile;
int n=0;
int rep=0;
bool exhaustive=false;
string fixed_cols;

/*
 * generate uniformly a random bit-vector of length L with
 * n bits set
 */
vector<bool> rand_vec(int L, int n){

	auto counts = vector<int>(n+1);

	for(int i=0;i<L-n;++i) counts[rand()%(n+1)]++;

	vector<bool> B;

	for(int i=0;i<counts[0];++i) B.push_back(false);

	for(int i=1;i<counts.size();++i){

		B.push_back(true);
		for(int j=0;j<counts[i];++j) B.push_back(false);

	}

	return B;

}

void help(){

	cout << "Usage:	" << endl;
	cout <<	"  cvs-rk [OPTIONS]" << endl;

	cout << "  MANDATORY OPTIONS:" << endl;
	cout << "    -i <file.txt> input file containing a M x L ASCII matrix "<<endl;
	cout << "                  stored row-wise" << endl;
	cout <<	"    -n <n> is the number of bits set in the initial solution" <<endl;
	cout <<	"    -r <rep> is the number of repetitions of the search algorithm." <<endl;
	cout <<	"       At each repetition, the search starts from a random" <<endl;
	cout <<	"       bitvector with n bits set." <<endl;

	cout << "  OTHER OPTIONS:" << endl;
	cout << "    -e exhaustive neighbor search: pick the best among all L(N-L)" << endl;
	cout <<	"       neighbors. Default: false (random search: pick the first" <<endl;
	cout <<	"       neighbor that improves the counts entropy). " <<endl;

	//cout <<	"    -f <c1,c2,c3,...,cm> string specifying which columns to keep" <<endl;
	//cout <<	"       fixed during search" <<endl;
	exit(0);

}

void parse_args(int & i,int tot_args, char** argv){

	if(i>=tot_args) help();

	//read option
	string op(argv[i++]);

	if(op.compare("-i")==0){

		//read argument
		infile = string(argv[i++]);

	}
	if(op.compare("-n")==0){

		//read argument
		n = atoi(argv[i++]);

	}
	if(op.compare("-r")==0){

		//read argument
		rep = atoi(argv[i++]);

	}
	if(op.compare("-e")==0){

		//read argument
		exhaustive = true;

	}
	/*if(op.compare("-f")==0){

		//read argument
		infile = string(argv[i++]);

	}*/


}

int main(int argc,char** argv) {

	srand(time(NULL));

	int i=1;

	if(argc<=1) help();

	while(i<argc) parse_args(i,argc,argv);

	if(infile.compare("")==0 or n == 0 or rep == 0) help();

	cout << "input file: " << infile << endl;

	vector<string> rows;

	std::string line;
	int nr_columns=-1;

	std::ifstream infile_f(infile);

	while (std::getline(infile_f, line)){

		assert(nr_columns<0 or nr_columns == line.length());

		rows.push_back(line);
		nr_columns = line.length();

	}

	assert(n<nr_columns);

	cout << "Input matrix size is " << rows.size() << " x " << nr_columns << endl;

	auto A = local_search(rows);

	//number of times a position is selected as critical
	auto critical_var_counts = vector<ulint>(nr_columns);

	double mean_HK=0;
	double mean_H=0;

	auto hks = vector<double>(rep);
	auto hs = vector<double>(rep);

	//repeat rep times
	for(ulint r = 0;r<rep;++r){

		auto B = rand_vec(nr_columns,n);

		assert(n==weight(B));

		local_search::search_result result;

		if(exhaustive)
			result = A.run_fixed_n_all_neighbors(B);
		else
			result = A.run_fixed_n_first_neighbor_eq(B);

		ulint j=0;
		for(auto b:result.B) critical_var_counts[j++] += b;

		mean_HK += result.HK;
		mean_H += result.H;

		hks[r] = result.HK;
		hs[r] = result.H;

		cout<<endl;
		for(auto b:result.B) cout << b;cout<<endl;
		cout << "Counts entropy: H(K) = " << result.HK << endl;
		cout << "Set entropy: H(S) = " << result.H << endl;
		cout << "distinct rows = " << result.distinct_rows << endl;

	}

	mean_HK /= rep;
	mean_H /= rep;

	double var_hk=0;
	for(auto h:hks) var_hk += pow(h-mean_HK,2);

	double var_h=0;
	for(auto h:hs) var_h += pow(h-mean_H,2);


	cout << endl << "Search terminated. Results:" << endl;
	cout << "Mean H(K) over " << rep << " runs: " << mean_HK << ". Variance = " << var_hk << endl;
	cout << "Mean H(S) over " << rep << " runs: " << mean_H << ". Variance = " << var_h << endl;

	cout << endl << "For each matrix column, number of times the column has been selected:" << endl;
	for(auto c:critical_var_counts) cout << c << "\t";

}
