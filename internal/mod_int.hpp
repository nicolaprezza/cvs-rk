/*
 * mod_int.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: nico
 *
 *  Integers modulo q<2^63, with addition, subtraction and multiplication
 *
 */

#ifndef INTERNAL_MOD_INT_HPP_
#define INTERNAL_MOD_INT_HPP_

#include <includes.hpp>

namespace cvs_rk{

template<ulint q>
class mod_int{

public:

	//new modular integer initialized with value x
	mod_int(ulint x = 0){

		assert(q>1);
		assert(q<(ulint(1)<<63));

		this->val=x%q;

	}

	//assign an integer to this modular int
	mod_int operator=(ulint x){

		this->val=x%q;
		return *this;

	}

	operator ulint(){
		return val;
	}

	//sum two modular integers
	mod_int operator+(mod_int b){

		ulint l = this->val;
		ulint r = b.val;

		return {(l+r)%q};

	}

	//subtract two modular integers
	mod_int operator-(mod_int b){

		ulint l = this->val;
		ulint r = b.val;

		ulint res = (l>=r ? l-r : q-(r-l) );

		return {res};

	}

	//multiply two modular integers
	mod_int operator*(mod_int b){

		auto m = __uint128_t(this->val)*__uint128_t(b.val);

		return {ulint(m)%q};

	}

	bool operator<(mod_int& b){
		return this->val < b.val;
	}

private:

	ulint val;

};

/*
 * component-wise sum between vectors of modular integers
 */
template<ulint q>
std::vector<mod_int<q> > operator+(std::vector<mod_int<q> >& a, std::vector<mod_int<q> >& b){

	assert(a.size() == b.size());

    auto result = std::vector<mod_int<q> >(a.size());

    for(ulint i=0;i<a.size();++i) result[i] = a[i] + b[i];

    /*std::transform(	a.begin(),
    				a.end(),
					b.begin(),
					result.begin(),
					[](mod_int<q>& a, mod_int<q>& b) { return a+b; });*/

    return result;
}

/*
 * component-wise difference between vectors of modular integers
 */
template<ulint q>
std::vector<mod_int<q> > operator-(std::vector<mod_int<q> >& a, std::vector<mod_int<q> >& b){

	assert(a.size() == b.size());

    auto result = std::vector<mod_int<q> >(a.size());

    for(ulint i=0;i<a.size();++i) result[i] = a[i] - b[i];

    /*std::transform(	a.begin(),
    				a.end(),
					b.begin(),
					result.begin(),
					[](mod_int<q>& a, mod_int<q>& b) { return a-b; });*/

    return result;
}

}

#endif /* INTERNAL_MOD_INT_HPP_ */
