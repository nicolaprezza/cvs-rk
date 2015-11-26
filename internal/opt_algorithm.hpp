/*
 * opt_algorithm.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: nico
 */

#ifndef INTERNAL_OPT_ALGORITHM_HPP_
#define INTERNAL_OPT_ALGORITHM_HPP_

#include <includes.hpp>
#include <matrix.hpp>

namespace cvs_rk{

class opt_algorithm{

public:

	typedef pair<vector<bool>, vector<matrix_t::mod_int_t> > candidate_t;

	opt_algorithm(){}

	opt_algorithm(vector<string>& rows){

		M = {rows};

	}

	void run(){

		//start point: all columns selected
		candidate_t C = {

				vector<bool>(M.n_columns(),true),
				M.hashed_rows()

		};

		explored.insert(C.first);

		//best H(K): entropy of current solution
		double best_HK = HK<matrix_t::mod_int_t>(C.second);
		double C_HK = best_HK;

		candidate_t best_C = C;

		bool local_max = false;

		//proceed until a local maximum is reached
		while(not local_max){

			//init as true, then check all neighbors
			//if at least 1 neighbor improves the score
			//of C, then local_max = false
			local_max = true;

			//try a flip in all column positions
			for(ulint j = 0;j<M.n_columns();++j){

				auto B = C.first;
				B[j] = B[j] ^ true; //flip bit

				//if new candidate has not been tried yet
				if(explored.find(B) == explored.end()){

					explored.insert(B);

					bool sum = (not C.first[j]) and B[j]; //0->1

					candidate_t C1;

					if(sum){

						C1 = {
								B,
								C.second + M.column(j)
						};

					}else{

						C1 = {
								B,
								C.second - M.column(j)
						};

					}

					double new_HK = HK<matrix_t::mod_int_t>(C1.second);

					if(new_HK>=best_HK){

						best_HK = new_HK;
						best_C = C1;

					}

					//if at least 1 neighbor strictly improves
					//C's entropy, then C is not a local maximum
					local_max = best_HK <= C_HK;

				}

			}

			C  = best_C;
			C_HK = best_HK;

			cout << "Current solution: " << endl;
			for(auto b:C.first) cout << b;cout<<endl;
			cout << "H(K) = " << C_HK << endl;

			if(local_max){

				cout << "Local maximum reached. Search terminated." << endl;

			}

		}

		cout << endl;
		cout << "Solution found:" << endl;
		for(auto b:C.first) cout << b;cout<<endl;
		cout << "H(K) = " << C_HK << endl;

	}

private:

	//matrix of characters
	matrix_t M;

	//elements of the search space that have already been explored
	set<vector<bool> > explored;


};

}

#endif /* INTERNAL_OPT_ALGORITHM_HPP_ */
