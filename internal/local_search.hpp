/*
 * local_search.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: nico
 */

#ifndef INTERNAL_LOCAL_SEARCH_HPP_
#define INTERNAL_LOCAL_SEARCH_HPP_

#include <includes.hpp>
#include <matrix.hpp>

namespace cvs_rk{

class local_search{

public:

	typedef pair<vector<bool>, vector<matrix_t::mod_int_t> > candidate_t;

	local_search(){}

	local_search(vector<string>& rows){

		M = {rows};

	}

	/*
	 * default: all columns (vector <111...1>)
	 */
	void run(vector<bool> initial_solution = {}){

		if(initial_solution.size()==0) initial_solution = vector<bool>(M.n_columns(),true);
		auto hashed_vec = M.hashed_rows();

		for(ulint i=0;i<M.n_columns();++i){

			//remove column from hashed vector
			if(not initial_solution[i]){

				hashed_vec = hashed_vec - M.column(i);

			}

		}

		candidate_t C = {

				initial_solution,
				hashed_vec

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

			//try all neighbors
			//the set of neighbors is built trying all possible O(L^2) combinations
			//of flips 0->1 and 1->0
			for(ulint j0 = 0;j0<=M.n_columns();++j0){//0->1. When j0=M.n_columns(), do not flip any 0-bit

				for(ulint j1 = 0;j1<=M.n_columns();++j1){//1->0. When j1=M.n_columns(), do not flip any 1-bit

					auto B = C.first;
					candidate_t C1;		//neighbor

					bool new_solution = false;

					//only one between a 0->1 flip and a 1->0 flip
					if(j0 == M.n_columns() xor j1 == M.n_columns()){

						if(j0 == M.n_columns() and B[j1]){

							B[j1] = false;

							if(explored.find(B) == explored.end()){

								new_solution = true;

								C1 = {
										B,
										C.second - M.column(j1)
								};

							}

						}

						if(j1 == M.n_columns() and not B[j0]){

							B[j0] = true;

							if(explored.find(B) == explored.end()){

								new_solution = true;

								C1 = {
										B,
										C.second + M.column(j0)
								};

							}

						}

					}

					//both flips
					if(j0 < M.n_columns() and j1 < M.n_columns()){

						if(not B[j0] and B[j1]){

							B[j0] = true;
							B[j1] = false;

							new_solution = true;

							if(explored.find(B) == explored.end()){

								C1 = {
										B,
										(C.second + M.column(j0))-M.column(j1)
								};

							}

						}

					}

					//if new candidate has not been tried yet
					if(new_solution){

						explored.insert(B);

						double new_HK = HK<matrix_t::mod_int_t>(C1.second);

						//cout << new_HK << endl;

						if(new_HK>=best_HK){

							best_HK = new_HK;
							best_C = C1;

						}

						//if at least 1 neighbor strictly improves
						//C's entropy, then C is not a local maximum
						local_max = best_HK <= C_HK;

					}

				}

			}//end neighbor search cycle

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

#endif /* INTERNAL_LOCAL_SEARCH_HPP_ */