#ifndef INCLUDES_HPP_
#define INCLUDES_HPP_

#include <algorithm>
#include <functional>
#include "stdint.h"
#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <math.h>

using namespace std;

namespace cvs_rk{

#define WORD_SIZE 64;

typedef unsigned char uchar;
typedef uint64_t ulint;
typedef uint32_t uint;

using namespace std;


/*
 * input: vector S of objects having an operator < defined
 * output: max count k (i.e. number of times a sequence appears repeated)
 */
template<typename T>
ulint get_max_k(vector<T>& S){

	std::sort(S.begin(),S.end());

	ulint k = 0;
	ulint max_k=0;

	for(ulint i = 0;i<S.size();++i){

		if( i>0 && ( S[i-1] < S[i] ) ){

			if(k>max_k) max_k=k;

			k = 1;

		}else{

			k++;

		}

	}

	if(k>max_k) max_k=k;

	return max_k;

}

/*
 * input: vector S of objects having an operator < defined
 * output: vector [m_k], k=0,1, ... (m_0 = 0 always). m_k = number of
 * distinct S elements that appear exactly k times in S
 */
template<typename T>
vector<ulint> get_counts_mk(vector<T>& S){

	ulint max_k = get_max_k<T>(S);

	auto m_k = vector<ulint>(max_k+1,0);

	//std::sort(S.begin(),S.end());
	ulint k = 0;

	for(ulint i = 0;i<S.size();++i){

		if( i>0 && ( S[i-1] < S[i] ) ){

			m_k[k]++;
			k = 1;

		}else{

			k++;

		}

	}

	m_k[k]++;

	return m_k;

}

/*
 * input: vector S of objects having an operator < defined
 * output: counts entropy of vector S
 */
template<typename T>
double HK(vector<T> S){

	auto m_k = get_counts_mk<T>(S);

	double Hk = 0;
	double M = 0;

	for(ulint k=0;k<m_k.size();++k){

		M += k*m_k[k];

	}

	for(ulint k=0;k<m_k.size();++k){

		if(m_k[k]>0) Hk -= (double(k*m_k[k])/M) * log2(double(k*m_k[k])/M);

	}

	return Hk;

}


}

#endif /* INCLUDES_HPP_ */
