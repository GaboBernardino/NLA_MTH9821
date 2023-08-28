//Gabriele Bernardino - Cholesky decomposition from Dan Stefanica's NLA Primer

#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include "Matrix.hpp"
#include "UpperTriangular.hpp"
#include "qbVector.hpp"
#include "LU_decomposition.hpp"
#include <limits>
using namespace std;

//************************************************************
// cholesky decomposition - return upper tri s.t. U'U=A AND the line where decomposition failed
//************************************************************
template <typename T, int N>
pair<UpperTriangular<T, N>, int> Cholesky(const Matrix<T, N, N> A_) {

	Matrix<T, N, N> A(A_); //copy bc it gets modified

	UpperTriangular<T, N> U;

	for (int i = 0; i < N - 1; ++i) {
		if (A.at(i, i) < T(0) || _close_enough(A.at(i, i), T(0))) {
			return make_pair(UpperTriangular<T, N>(), i);
		}
		U.at(i, i) = sqrt(A.at(i, i));
		T num = U.at(i, i);

		//i-th row of U
		for (int k = i + 1; k < N; ++k) {
			U.at(i, k) = A.at(i, k) / num;
		}

		for (int j = i + 1; j < N; ++j) {
			for (int k = j; k < N; ++k) {
				A.at(j, k) = A.at(j, k) - U.at(i, j) * U.at(i, k);
			}
		}
	}

	if (A.at(N - 1, N - 1) <= T(0)) return make_pair(UpperTriangular<T, N>(), N - 1);
	U.at(N - 1, N - 1) = sqrt(A.at(N - 1, N - 1));

	return make_pair(U, N);
}

//specialization for banded matrices - MUST BE SYMMETRIC!!!
template <typename T, int N>
pair<BandedMatrix<T, N, N, 1, 0>, bool> Cholesky(const BandedMatrix<T, N, N, 1, 1> A_) {
	
	BandedMatrix<T, N, N, 1, 1> A(A_); //copy bc it gets modified

	BandedMatrix<T, N, N, 1, 0> U;
	bool isValid = true;

	for (int i = 0; i < N - 1; ++i) {
		if (A.at(i, i) < T(0) || _close_enough(A.at(i, i), T(0))) {
			isValid = false; //OR MAYBE =i
			break;
		}
		U.at(i, i) = sqrt(A.at(i, i));
		U.at(i, i + 1) = A.at(i, i + 1) / U.at(i, i);

		A.at(i + 1, i + 1) = A.at(i + 1, i + 1) - pow(U.at(i, i + 1), 2);
	}

	if (A.at(N - 1, N - 1) <= T(0)) isValid = false;
	else U.at(N - 1, N - 1) = sqrt(A.at(N - 1, N - 1));

	return make_pair(U, isValid);
}

//************************************************************
// cholesky decomposition to solve linear system
//************************************************************
template <typename T, int N>
qbVector<T> Cholesky_solve(const Matrix<T, N, N>& A, const qbVector<T>& b) {
	auto [U, check] = Cholesky(A); //cholesky decompose A

	if (check == N) {
		LowerTriangular<T, N> L = U.Transpose();

		qbVector<T> y = forward_subst(L, b);
		qbVector<T> x = backward_subst(U, y); //solution

		return x;
	}

	else {
		cout << "Decomposition failed." << endl;
		return qbVector<T>(N, numeric_limits<T>::infinity()); //maybe just (N)
	}
}

//specialization for banded matrix - MUST BE SYMMETRIC!!!
template <typename T, int N>
qbVector<T> Cholesky_solve(const BandedMatrix<T, N, N, 1, 1>& A, const qbVector<T>& b) {
	
	auto [U, check] = Cholesky(A); //cholesky decompose A

	if (check) {
		BandedMatrix<T, N, N, 0, 1>L = U.Transpose();

		qbVector<T> y = forward_subst(L, b);
		qbVector<T> x = backward_subst(U, y); //solution

		return x;
	}

	else {
		cout << "Decomposition failed." << endl;
		return qbVector<T>(N, numeric_limits<T>::infinity()); //maybe just (N)
	}
}

//************************************************************************
// OLS solver
//************************************************************************
template <typename T, int NR, int NC>
qbVector<T> OLS(const Matrix<T, NR, NC>& X, const qbVector<T>& y) {

	qbVector<T> b = Cholesky_solve((X.Transpose() * X), (X.Transpose() * y));
	return b;
}

#endif // !CHOLESKY_HPP
