/*
 * matrix.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: nico
 */

#ifndef INTERNAL_MATRIX_HPP_
#define INTERNAL_MATRIX_HPP_


#include <includes.hpp>
#include <rk_function.hpp>

namespace cvs_rk{

template<class hash_function_t>
class matrix{

public:

	//modular integer type (defined inside the hash function
	using mod_int_t = typename hash_function_t::mod_int_t;

	//empty constructor
	matrix(){}

	//matrix constructor: takes as input a NONEMPTY string matrix (in format vector of rows)
	matrix(vector<string>& rows){

		assert(rows.size()>0);
		ulint c = rows[0].size();

		//all rows must have the same length
		for(auto r:rows){
			assert(r.size() == c);
		}

		this->rows = rows;

		columns = { c, vector<mod_int_t>(rows.size()) };

		//for each row, compute associated
		//weighted modular vector and transfer it
		//into columns
		for(ulint i=0;i<rows.size();++i){

			vector<mod_int_t> r_mod = H.to_mod_vector(rows[i]);

			for(ulint j = 0;j<c;++j){

				columns[j][i] = r_mod[j];

			}

		}

	}

	//i-th row as string
	string& row(ulint i){

		assert(i<rows.size());
		return rows[i];

	}

	/*
	* i-th column as vector of weighted integers modulo q
	*
	* equals i-th column in string format multiplied by s^i (mod q),
	* where s = 256 is the alphabet size
	*
	*/
	vector<mod_int_t>& column(ulint i){

		assert(i<columns.size());
		return columns[i];

	}

	ulint n_rows(){
		return rows.size();
	}

	ulint n_columns(){
		return columns.size();
	}

	/*
	 * returns the vector of rows hash values
	 */
	vector<mod_int_t> hashed_rows(){

		auto result = vector<mod_int_t>(n_rows());

		ulint i = 0;
		for(auto r : rows) result[i++] = H(r);

		return result;

	}

private:

	hash_function_t H;

	vector<string> rows;
	vector<vector<mod_int_t> > columns;

};

typedef matrix<rk_function_t> matrix_t;

}

#endif /* INTERNAL_MATRIX_HPP_ */
