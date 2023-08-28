//Gabriele Bernardino - Functions to implement asset allocation pseudocodes from Dan Stefanica's NLA Primer

#ifndef ASSETALLOC_HPP
#define ASSETALLOC_HPP

#include "Matrix.hpp"
#include "qbVector.hpp"
#include "Cholesky.hpp"
#include <tuple>
using namespace std;

//*********************************************************************************
// Minimum variance portfolio
//*********************************************************************************
template <typename Ty> using resTuple = tuple<qbVector<Ty>, Ty, Ty>; //return type
//function returning (1) portfolio weights, (2) weight of cash, (3) std dev of portfolio
template <typename T, int N>
resTuple<T> minimumVariance(const Matrix<T, N, N>& cov, const qbVector<T>& mu, const T& rf, const T& mu_req) {
	
	if (mu.size() != N) {
		throw invalid_argument("Size of expected return vector does not match.");
	}

	qbVector<T> excessRet = mu - rf;
	qbVector<T> x = Cholesky_solve(cov, excessRet);

	qbVector<T> w_min = ((mu_req - rf) / qbVector<T>::dot(excessRet, x)) * x;

	T w_min_cash = 1 - w_min.Sum();
	//for (int i = 0; i < N; ++i) w_min_cash -= w_min[i];
	T sig = sqrt(qbVector<T>::dot(w_min * cov, w_min));

	resTuple<T> res{ w_min,w_min_cash,sig };
	return res;
}

#endif // !ASSETALLOC_HPP
