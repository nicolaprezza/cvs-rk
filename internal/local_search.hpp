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
#include <map>

using namespace std;

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
	 *
	 * all O(L^2) neighbors are explored
	 *
	 */
	search_result run_fixed_n_all_neighbors(vector<bool> initial_solution){

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

							auto sum = C.second + M.column(j0);	//sum j0-th column
							auto new_vec = sum-M.column(j1);	//subtract j1-th column

							//new candidate
							new_C = {
									B,
									new_vec
							};

							//compute counts entropy
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

	/*
	 * run local search keeping n fixed
	 * default: all columns (vector <111...1>)
	 *
	 * we pick random neighbors until we find one that
	 * improves the counts entropy
	 *
	 */
	search_result run_fixed_n_first_neighbor(vector<bool> initial_solution){

		std::srand ( time(NULL) );

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

		//to avoid loops, in explored we keep all already explored bitvectors
		explored.insert(C.first);

		//best H(K): best entropy found up to now
		double C_HK = HK<matrix_t::mod_int_t>(C.second);

		//we stop search when eaither we found a neighbor that improves
		//the entropy OR we visited all neighbors (local maximum found)
		bool stop_search = false;

		while(not stop_search){

			stop_search = true;

			auto bitvector = C.first;

			//fill a vector with positions of 1's
			vector<int> ones;
			vector<int> zeros;
			{
				int i=0;
				for(auto b:bitvector)
					if(b)
						ones.push_back(i++);
					else
						zeros.push_back(i++);
			}

			//this map maps the position of a 1 to a
			//pair <N,vec>, where:
			//vec is a random permutation of the positions of the 0's
			//N is how many elements we have left to extract from vec
			//at the beginning, the map is empty. Every time we query
			//an empty bucket ones_to_zeros[x], we set
			//ones_to_zeros[x] = <zeros.size(),random_shuffle(zeros)>
			//if we extracted all possible elements from ones_to_zeros[x],
			//then we remove x from ones.
			std::map<int, pair<int,vector<int> > > ones_to_zeros;

			//do we have to continue visiting neighbors?
			bool continue_visit = true;

			int n_nb=0;//number of visited neighbors

			while(continue_visit){

				//pick a random 1 position
				int rand_pos = rand()%ones.size();
				int one_pos = ones[rand_pos];

				//we haven't yet created the permutation of 0's. create it.
				if(ones_to_zeros[one_pos].second.size()==0){

					std::random_shuffle(zeros.begin(),zeros.end());

					ones_to_zeros[one_pos] = { zeros.size(), zeros };

				}

				//now ones_to_zeros[one_pos] contains a proper pair <N, vec>

				assert(ones_to_zeros[one_pos].first>0);
				int zero_pos = ones_to_zeros[one_pos].second[ zeros.size()-ones_to_zeros[one_pos].first-1 ];

				ones_to_zeros[one_pos].first--;

				//if no more elements to extract, delete one_pos from ones
				if(ones_to_zeros[one_pos].first==0) ones.erase(ones.begin() + rand_pos);

				//now we have a pair <one_pos, zero_pos> of positions: compute new candidate

				auto B = C.first;		//bitvector of neighbor
				candidate_t neighbor_C;		//neighbor

				B[zero_pos] = true;		//flip bits
				B[one_pos] = false;

				if(explored.find(B) == explored.end()){

					explored.insert(B);

					auto sum = C.second + M.column(zero_pos);	//sum zero_pos-th column
					auto new_vec = sum-M.column(one_pos);	//subtract one_pos-th column

					//new candidate
					neighbor_C = {
							B,
							new_vec
					};

					//compute counts entropy
					double neighbor_HK = HK<matrix_t::mod_int_t>(neighbor_C.second);

					//cout << new_HK << endl;

					if(neighbor_HK>C_HK){

						C_HK = neighbor_HK;
						C = neighbor_C;

						continue_visit = false;

					}

				}

				//no more neighbors to visit
				if(ones.size()==0) continue_visit=false;

				n_nb++;

			}

			cout << n_nb <<  " visited neighbors"<<endl;

			//if we visited all neighbors, stop search.
			//C is the best solution found.
			stop_search = ones.size()==0;

			cout << "Current best solution: H(K) = " << C_HK << endl;

			if(stop_search){

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
