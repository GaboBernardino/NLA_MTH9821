//Gabriele Bernardino - Functions to implement LU decomposition and solve linear systems (Chapter 2 of Dan Stefanica's NLA Primer)

#ifndef LU_DECOMPOSITION_HPP
#define LU_DECOMPOSITION_HPP

#include "UpperTriangular.hpp"
#include "LowerTriangular.hpp"
#include "BandedMatrix.hpp"
#include "qbVector.hpp"
#include "helpers.hpp"
#include <tuple>
#include <numeric>
#include <type_traits>
using namespace std;

//concepts for banded matrices:
//template <typename matr, typename Ty, int NR, int NC>
//concept IsBidiagLo = is_same_v<BandedMatrix<Ty, NR, NC, 0, 1>, matr>;

//this function computes the solution x to Lx=b, where L is lower triangular
template <typename T, int N>
qbVector<T> forward_subst(const LowerTriangular<T, N>& L, const qbVector<T>& b) {
	//size check
	if (N != b.size()) {
		throw invalid_argument("Sizes of inputs don't match.");
	}
	//defining variables
	qbVector<T> x(N);
	T runningSum;
	//algo:
	x.setElem(0, b[0] / L[{0, 0}]);
	for (int i = 1; i < N; ++i) {
		runningSum = T(0);
		for (int j = 0; j < i; ++j) {
			runningSum += L[{i, j}] * x[j];
		}
		if (_close_enough(L[{i, i}], T(0))) {
			throw invalid_argument("Matrix is singular");
		}
		x.setElem(i, (b[i] - runningSum) / L[{i, i}]);
	}

	return x;
}

//specialization for bidiagonal matrices
template <typename T, int N>
qbVector<T> forward_subst(const BandedMatrix<T, N, N, 0, 1>& L, const qbVector<T>& b) {
	//size check
	if (N != b.size()) {
		throw invalid_argument("Sizes of inputs don't match.");
	}
	//defining variables
	qbVector<T> x(N);
	//algo:
	x.setElem(0, b[0] / L[{0, 0}]);
	for (int i = 1; i < N; ++i) {
		if (_close_enough(L[{i, i}], T(0))) {
			throw invalid_argument("Matrix is singular");
		}
		x.setElem(i, (b[i] - L[{i, i - 1}] * x[i - 1]) / L[{i, i}]);
	}

	return x;
}


//this function computes the solution x to Ux=b, where U is upper triangular
template <typename T, int N>
qbVector<T> backward_subst(const UpperTriangular<T, N>& U, const qbVector<T>& b) {
	//size check
	if (N != b.size()) {
		throw invalid_argument("Sizes of inputs don't match.");
	}
	//defining variables
	qbVector<T> x(N);
	T runningSum;
	//algo:
	x.setElem(N - 1, b[N - 1] / U.at(N - 1, N - 1));
	for (int i = N - 1; i > 0; i--) { //size_t is an unsigned type!
		runningSum = T(0);
		for (int j = i; j < N; ++j) {
			runningSum += U.at(i - 1, j) * x[j];
		}
		x.setElem(i - 1, (b[i - 1] - runningSum) / U.at(i - 1, i - 1));
	}
	return x;
}

//specialization for bidiagonal matrices
template <typename T, int N>
qbVector<T> backward_subst(const BandedMatrix<T, N, N, 1, 0>& U, const qbVector<T>& b) {
	//size check
	if (N != b.size()) {
		throw invalid_argument("Sizes of inputs don't match.");
	}
	//defining variables
	qbVector<T> x(N);
	//algo:
	x.setElem(N - 1, b[N - 1] / U[{N - 1, N - 1}]);
	for (int i = N - 2; i >= 0; --i) {
		if (U[{i, i}] == 0) {
			throw invalid_argument("Matrix is singular");
		}
		x.setElem(i, (b[i] - U[{i, i + 1}] * x[i + 1]) / U[{i, i}]);
	}

	return x;
}

//**********************************************************************************
// FUNCTIONS - LU decomposition
//**********************************************************************************

template <typename Ty, int n> using LUpair = pair<LowerTriangular<Ty, n>, UpperTriangular<Ty, n>>;
//this function returns the LU decomposition without row pivoting (NB: L(i,i)=1 for all i)
template <typename T, int N>
LUpair<T, N> LU_no_pivot(const Matrix<T, N, N>& A_) {
	Matrix<T, N, N> A(A_); //create copy bc it gets modified
	LowerTriangular<T, N> L; UpperTriangular<T, N> U;

	for (int i = 0; i < N - 1; ++i) {
		for (int k = i; k < N; ++k) {
			U.at(i, k) = A.at(i, k); //i-th row of U
			L.at(k, i) = A.at(k, i) / U.at(i, i); //i-th column of l
		}
		for (int j = i; j < N; ++j) {
			for (int k = i; k < N; k++) {
				A.at(j, k) = A.at(j, k) - L.at(j, i) * U.at(i, k);
			}
		}
	}

	L.at(N - 1, N - 1) = 1; U.at(N - 1, N - 1) = A.at(N - 1, N - 1);

	LUpair<T, N> res(L, U);
	return res;
}

//linear system solver using LU decomposition without pivoting
template <typename T, int N>
qbVector<T> LU_solve_no_pivot(const Matrix<T, N, N>& A, const qbVector<T>& b) {
	LUpair<T, N> lu = LU_no_pivot(A);
	qbVector<T> y = forward_subst(lu.first, b);
	qbVector<T> x = backward_subst(lu.second, y);
	return x;
}

//specializations for tridiagonal matrices
template <typename Ty, int n> using LUpairDiag = pair<BandedMatrix<Ty, n, n, 0, 1>, BandedMatrix<Ty, n, n, 1, 0>>;
//this function returns the LU decomposition without row pivoting of a tridiagonal matrix
template <typename T, int N>
LUpairDiag<T, N> LU_no_pivot(const BandedMatrix<T, N, N, 1, 1>& A_) {
	Matrix<T, N, N> A(A_); //create copy bc it gets modified
	BandedMatrix<T, N, N, 0, 1> L; BandedMatrix<T, N, N, 1, 0> U;

	for (int i = 0; i < N - 1; ++i) {
		L.at(i, i) = T(1);
		L.at(i + 1, i) = A.at(i + 1, i) / A.at(i, i);

		U.at(i, i) = A.at(i, i);
		U.at(i, i + 1) = A.at(i, i + 1);

		A.at(i + 1, i + 1) = A.at(i + 1, i + 1) - L.at(i + 1, i) * U.at(i, i + 1);
	}
	L.at(N - 1, N - 1) = T(1); U.at(N - 1, N - 1) = A.at(N - 1, N - 1);

	LUpairDiag<T, N> res(L, U);
	return res;
}

template <typename Ty, int n, int band> using LUpairBand = pair<BandedMatrix<Ty, n, n, 0, band>, BandedMatrix<Ty, n, n, band, 0>>;
//this function returns the LU decomposition without row pivoting of a banded matrix of any band
template <typename T, int N, int BAND>
LUpairBand<T, N, BAND> LU_no_pivot_band(const BandedMatrix<T, N, N, BAND, BAND>& A_) {
	Matrix<T, N, N> A(A_); //create copy bc it gets modified
	BandedMatrix<T, N, N, 0, BAND> L; BandedMatrix<T, N, N, BAND, 0> U;

	for (int i = 0; i < N - 1; ++i) {
		//cout << i << endl;
		L.at(i, i) = T(1);
		U.at(i, i) = A.at(i, i);

		for (int k = i + 1; k < min(N, i + BAND + 1); ++k) { //only need to loop to the band
			U.at(i, k) = A.at(i, k); //i-th row of U
			L.at(k, i) = A.at(k, i) / U.at(i, i); //i-th column of L
		}
		//now update A
		for (int j = i + 1; j < min(N, i + BAND + 1); ++j) {
			for (int k = i + 1; k < min(N, i + BAND + 1); ++k) {
				A.at(j, k) = A.at(j, k) - L.at(j, i) * U.at(i, k);
			}
		}
	}
	L.at(N - 1, N - 1) = T(1); U.at(N - 1, N - 1) = A.at(N - 1, N - 1);

	LUpairBand<T, N, BAND> res(L, U);
	return res;
}

//linear system solver using LU decomposition without pivoting - specialization for banded matrices
template <typename T, int N>
qbVector<T> LU_solve_no_pivot(const BandedMatrix<T, N, N, 1, 1>& A, const qbVector<T>& b) {
	LUpairDiag<T, N> lu = LU_no_pivot(A);
	qbVector<T> y = forward_subst(lu.first, b);
	qbVector<T> x = backward_subst(lu.second, y);
	return x;
}


template <typename Ty, int n> using LUtuple = tuple<Matrix<Ty, n, n>, LowerTriangular<Ty, n>, UpperTriangular<Ty, n>>;
//LU decomposition WITH row pivoting
template <typename T, int N>
LUtuple<T, N> LU_row_pivot(const Matrix<T, N, N>& A_) { //NB: P is stored as vector of position of 1 on i-th line
	Matrix<T, N, N> A(A_); //create copy bc it gets modified
	LowerTriangular<T, N> L; UpperTriangular<T, N> U;
	L.SetToId(); //initialize to identity matrix
	vector<int> p(N); iota(p.begin(), p.end(), 0); //init to range from 0 to N-1


	for (int i = 0; i < N - 1; ++i) {
		int imax = i;
		for (int m = i + 1; m < N; ++m) {
			if (fabs(A.at(m, i)) > fabs(A.at(imax, i))) imax = m; //find index of largest entry of what's left of i-th column
		}
		//swap row i and imax of A
		A.SwapRows(i, imax, i, N);
		//swap i and imax of P
		std::swap(p[i], p[imax]);
		if (i > 0) {
			L.SwapRows(i, imax); //swap row i and imax of L
		}
		for (int j = i; j < N; ++j) {
			L.at(j, i) = A.at(j, i) / A.at(i, i); //i-th col of L
			U.at(i, j) = A.at(i, j); //i-th row of U
		}
		for (int j = i + 1; j < N; ++j) {
			for (int k = i + 1; k < N; k++) {
				A.at(j, k) = A.at(j, k) - L.at(j, i) * U.at(i, k); //update remaining matrix
			}
		}
	}
	L.at(N - 1, N - 1) = 1; U.at(N - 1, N - 1) = A.at(N - 1, N - 1);

	//construct permutation matrix from vector of entries:
	Matrix<T, N, N> P;
	for (int i = 0; i < N; ++i) P.at(p[i], i) = T(1);

	LUtuple<T, N> res{ P,L,U };
	return res;
}

//linear system solver using LU decomposition with row pivoting
template <typename T, int N>
qbVector<T> LU_solve(const Matrix<T, N, N>& A, const qbVector<T>& b) {
	LUtuple<T, N> lu = LU_row_pivot(A);
	qbVector<T> y = forward_subst(get<1>(lu), get<0>(lu) * b);
	qbVector<T> x = backward_subst(get<2>(lu), y);
	return x;
}

#endif // !LU_DECOMPOSITION_HPP