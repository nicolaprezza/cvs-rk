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

		ulint l1 = this->val >> 32;
		ulint l2 = (this->val << 32) >> 32;

		ulint r1 = b.val >> 32;
		ulint r2 = (b.val << 32) >> 32;

		ulint l1_r1 = mult_pow_2(l1*r1,64);
		ulint l1_r2 = mult_pow_2(l1*r2,32);
		ulint l2_r1 = mult_pow_2(l2*r1,32);
		ulint l2_r2 = (l2*r2)%q;

		ulint res = (l1_r1 + l1_r2)%q;
		res = (res + l2_r1)%q;
		res = (res + l2_r2)%q;

		return {res};

	}

	bool operator<(mod_int& b){
		return this->val < b.val;
	}

private:

	//return x*2^i mod q
	ulint mult_pow_2(ulint x,uint i){

		x = x%q;

		for(uint j=0;j<i;++j) x = (x<<1)%q;

		return x;

	}

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
