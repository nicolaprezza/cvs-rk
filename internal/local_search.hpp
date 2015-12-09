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

	struct search_result{

		vector<bool> B;
		double HK;
		double H;
		int distinct_rows;

	};

	typedef pair<vector<bool>, vector<matrix_t::mod_int_t> > candidate_t;

	local_search(){}

	local_search(vector<string>& rows){

		M = {rows};

	}

	/*
	 * run local search keeping n fixed
	 * default: all columns (vector <111...1>)
	 */
	search_result run_fixed_n(vector<bool> initial_solution){

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
			//of both a flip 0->1 and a flip 1->0
			for(ulint j0 = 0;j0<M.n_columns();++j0){//0->1

				for(ulint j1 = 0;j1<M.n_columns();++j1){//1->0

					if(not C.first[j0] and C.first[j1]){

						auto B = C.first;
						candidate_t new_C;		//neighbor

						B[j0] = true;
						B[j1] = false;

						if(explored.find(B) == explored.end()){

							explored.insert(B);

							auto sum = C.second + M.column(j0);

							new_C = {
									B,
									sum-M.column(j1)
							};

							double new_HK = HK<matrix_t::mod_int_t>(new_C.second);

							//cout << new_HK << endl;

							if(new_HK>=best_HK){

								best_HK = new_HK;
								best_C = new_C;

							}

						}

					}

				}

			}//end neighbor search cycle

			//if at least 1 neighbor strictly improves
			//C's entropy, then C is not a local maximum
			local_max = best_HK <= C_HK;

			if(not local_max){
				C  = best_C;
				C_HK = best_HK;
			}

			cout << "Current best solution: H(K) = " << C_HK << endl;

			if(local_max){

				cout << "Local maximum reached. Search terminated." << endl;

			}

		}

		auto comp = [](matrix_t::mod_int_t x, matrix_t::mod_int_t y){ return x < y; };
		std::set<matrix_t::mod_int_t,decltype(comp)> distinct_rows(comp);

		for(auto c:C.second) distinct_rows.insert(c);

		return {C.first, C_HK, H<matrix_t::mod_int_t>(C.second), int(distinct_rows.size())};

	}

private:

	ulint psum(vector<bool>& B){
		ulint ps=0;
		for(auto b:B) ps += b;
		return ps;
	}

	//matrix of characters
	matrix_t M;

	//elements of the search space that have already been explored
	set<vector<bool> > explored;


};

}

#endif /* INTERNAL_LOCAL_SEARCH_HPP_ */
