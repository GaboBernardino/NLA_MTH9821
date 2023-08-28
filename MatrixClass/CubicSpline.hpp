//Gabriele Bernardino - Cubic spline interpolation from Dan Stefanica's NLA Primer

#ifndef CUBICSPLINE_HPP
#define CUBICSPLINE_HPP

#include "LU_decomposition.hpp"
using namespace std;

template <typename Ty> using CoefficientsTuple = tuple<qbVector<Ty>, qbVector<Ty>, qbVector<Ty>, qbVector<Ty>>;
//*****************************************************************************
// efficient cubic spline interpolation using tridiagonal solver
//*****************************************************************************
template <typename T, int N>
CoefficientsTuple<T> CubicSpline(qbVector<T> x, qbVector<T> v) {
	if (x.size() != N + 1 || v.size() != N + 1 || x.size() != v.size()) {
		throw invalid_argument("Dimensions do not match");
	}

	//z vector, rhs of linear system
	qbVector<T> z(N - 1);
	for (int i = 0; i < N - 1; ++i) {
		z.at(i) = 6 * ((v[i + 2] - v[i + 1]) / (x[i + 2] - x[i + 1]) - (v[i + 1] - v[i]) / (x[i + 1] - x[i]));
	}
	
	//matrix M of coefficients
	BandedMatrix<T, N - 1, N - 1, 1, 1> M;	
	for (int i = 0; i < N - 1; ++i) M.at(i, i) = 2 * (x[i + 2] - x[i]);
	for (int i = 0; i < N - 2; ++i) M.at(i, i + 1) = x[i + 2] - x[i + 1];
	for (int i = 1; i < N - 1; ++i) M.at(i, i - 1) = x[i + 1] - x[i];

	qbVector<double> w = LU_solve_no_pivot(M, z);
	qbVector<double> c(N), d(N);
	for (int i = 1; i < N - 1; ++i) {
		c[i] = (w[i - 1] * x[i + 1] - w[i] * x[i]) / (2 * (x[i + 1] - x[i]));
		d[i] = (w[i] - w[i - 1]) / (6 * (x[i + 1] - x[i]));
	}
	c[0] = -w[0] * x[0] / (2 * (x[1] - x[0])); c[N - 1] = w[N - 2] * x[N] / (2 * (x[N] - x[N - 1]));
	d[0] = w[0] / (6 * (x[1] - x[0])); d[N - 1] = -w[N - 2] / (6 * (x[N] - x[N - 1]));
	
	qbVector<double> q(N), r(N);
	for (int i = 0; i < N; ++i) {
		q[i] = v[i] - c[i] * x[i] * x[i] - d[i] * pow(x[i], 3);
		r[i] = v[i + 1] - c[i] * x[i + 1] * x[i + 1] - d[i] * pow(x[i + 1], 3);
	}

	qbVector<T> a(N), b(N);
	for (int i = 0; i < N; ++i) {
		a[i] = (q[i] * x[i + 1] - r[i] * x[i]) / (x[i + 1] - x[i]);
		b[i] = (r[i] - q[i]) / (x[i + 1] - x[i]);
	}

	CoefficientsTuple<T> res(a, b, c, d);
	return res;
}

//*****************************************************************************
// helper function to print coefficients
//*****************************************************************************
template <typename T>
void print_coefficients(const CoefficientsTuple<T>& t) {
	cout << "a: " << get<0>(t) << endl;
	cout << "b: " << get<1>(t) << endl;
	cout << "c: " << get<2>(t) << endl;
	cout << "d: " << get<3>(t) << endl;
	cout << endl;
}


#endif // !CUBICSPLINE_HPP
