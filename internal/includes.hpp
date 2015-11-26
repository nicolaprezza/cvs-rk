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
 * output: counts entropy of vector S
 */
template<typename T>
double HK(vector<T> S){

	std::sort(S.begin(),S.end());

	vector<ulint> counts;
	ulint c_temp = 0;

	for(ulint i = 0;i<S.size();++i){

		if( i>0 && ( S[i-1] < S[i] ) ){

			for(ulint j=0;j<c_temp;++j) counts.push_back(c_temp);

			c_temp = 1;

		}else{

			c_temp++;

		}

	}

	for(ulint j=0;j<c_temp;++j) counts.push_back(c_temp);

	assert(counts.size()==S.size());

	std::sort(counts.begin(),counts.end());

	double Hk = 0;
	c_temp = 0;

	for(ulint i = 0;i<counts.size();++i){

		if( (i>0 && ( counts[i]>counts[i-1] )) or ( i == counts.size()-1 ) ){

			if(i == counts.size()-1) c_temp++;

			double f = double(c_temp)/counts.size();
			Hk -= f*log2(f);

			c_temp = 1;

		}else{

			c_temp++;

		}

	}

	return Hk;

}

}

#endif /* INCLUDES_HPP_ */
