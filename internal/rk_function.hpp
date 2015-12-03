/*
 * rk_function.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: nico
 *
 *  encodes a Rabin-Karp hash function
 *
 */

#ifndef INTERNAL_RK_FUNCTION_HPP_
#define INTERNAL_RK_FUNCTION_HPP_

#include <includes.hpp>
#include <mod_int.hpp>

namespace cvs_rk{

//the class is a template on the modulo q
template<ulint q>
class rk_function{

public:

	//integers modulo q
	typedef mod_int<q> mod_int_t;

	rk_function(){

		//we must be able to multiply
		//by 2 without overflows
		assert(q>1);
		assert(q<(ulint(1)<<63));

		sigma = 256;

	}


	/*
	 * input: string
	 * output: a vector <c0, c1*s^1, c2*s^2, ..., cm*s^m> mod q, where c0,..., cm are
	 * the input string's characters and s is the alphabet size (256)
	 *
	 */
	vector<mod_int_t> to_mod_vector(string& S){

		mod_int_t coefficient = 1;

		auto result = vector<mod_int_t>(S.length());
		ulint i = 0;

		for(auto c : S){

			result[i++] = mod_int_t(c) * coefficient;
			coefficient = coefficient * sigma;

		}

		return result;

	}

	/*
	 * hash value of input string S
	 */
	mod_int_t operator()(string& S){

		mod_int_t H = 0;

		for(ulint i=0;i<S.length();++i){

			H = H * sigma;
			H = H + mod_int_t(S[S.length()-i-1]);

		}

		assert(H == sum(to_mod_vector(S)));

		return H;

	}

private:

	mod_int_t sum(vector<mod_int_t>& v){

		mod_int_t s=0;

		for(auto x:v) s = s + x;

		return s;

	}

	//alphabet size
	mod_int_t sigma;


};

//default: q is a prime near 2^63
typedef rk_function<(ulint(1)<<63)-25> rk_function_t;

}


#endif /* INTERNAL_RK_FUNCTION_HPP_ */
