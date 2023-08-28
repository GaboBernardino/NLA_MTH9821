//Gabriele Bernardino - derived class for Lower Triangular matrices

#ifndef UPPERTRIANGULAR_HPP
#define UPPERTRIANGULAR_HPP


#include "Matrix.hpp"
#include "LowerTriangular.hpp"
using namespace std;

template <class T, int N = 1>
class UpperTriangular : public Matrix<T, N, N> {
private:
	//using directives for more readable access to member data
	using Base = Matrix<T, N, N>;
public:
	//constructors and destructor
	UpperTriangular();
	UpperTriangular(const T& val);
	UpperTriangular(const T* inputData);
	UpperTriangular(const Matrix<T, N, N>& inputMatrix);
	UpperTriangular(const vector<T>* inputData);

	//accessing elements
	T& operator [](const pair<int, int>& idx);
	const T& operator [](const pair<int, int>& index) const;
	T& at(const int& row, const int& col);
	const T& at(const int& row, const int& col) const;

	// Override the setElem function to ensure elements outside the triangular region are set to zero
	bool setElem(const int& row, const int& col, const T& element);

	//Override SwapRows function to ensure matrix remains UpperTriangular
	void SwapRows(const int& r1, const int& r2, int start = 0, int end = N);

	//Override transpose to return a LowerTri
	LowerTriangular<T, N> Transpose() const;
};

/***********************************************************************************************************
Constructors and destructor
***********************************************************************************************************/
template <class T, int N> UpperTriangular<T, N>::UpperTriangular() : Base() {
	Base::mData = make_unique<T[]>(N * N);
	for (int i = 0; i < N * N; ++i) Base::mData[i] = T(0);
}

template <class T, int N> UpperTriangular<T, N>::UpperTriangular(const T& val) : Base() {
	Base::mData = make_unique<T[]>(N * N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (j >= i) Base::mData[i * N + j] = val;
			else Base::mData[i * N + j] = T(0);
		}
	}
}

template <class T, int N> UpperTriangular<T, N>::UpperTriangular(const T* inputData) : Base() {
	if ((sizeof(inputData) / sizeof(inputData[0])) != N * N) {
		throw invalid_argument("Input size does not match.");
	}
	Base::mData = make_unique<T[]>(N * N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (j >= i) Base::mData[i * N + j] = inputData[i * N + j];
			else Base::mData[i * N + j] = T(0);
		}
	}
}

template <class T, int N> UpperTriangular<T, N>::UpperTriangular(const vector<T>* inputData) : Base() {
	if (inputData->size() != (N * N)) {
		throw invalid_argument("Input size does not match.");
	}
	Base::mData = make_unique<T[]>(N * N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (j >= i) Base::mData[i * N + j] = inputData->at(i * N + j);
			else Base::mData[i * N + j] = T(0);
		}
	}
}

//copy ctor
template <class T, int N> UpperTriangular<T, N>::UpperTriangular(const Matrix<T, N, N>& inputMatrix) : Base() {
	Base::mData = make_unique<T[]>(N * N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (j >= i) Base::mData[i * N + j] = inputMatrix.mData[i * N + j];
			else Base::mData[i * N + j] = T(0);
		}
	}
}

//assignment operator
//..........

//*********************************************************
// element access
//*********************************************************
template <class T, int N>
bool UpperTriangular<T, N>::setElem(const int& row, const int& col, const T& element) {
	if (row > col) {
		// Elements below the main diagonal should be zero
		return false;
	}
	return Base::setElem(row, col, element);
}

template <class T, int N>
T& UpperTriangular<T, N>::operator [](const pair<int, int>& idx) {
	if (idx.second < idx.first) {
		throw out_of_range("Cannot modify elements below the diagonal.");
	}
	int index = Base::Sub2Ind(idx.first, idx.second);
	return Base::mData[index];
}

template <class T, int N>
const T& UpperTriangular<T, N>::operator [](const pair<int, int>& idx) const {
	return Base::operator[](idx);
}

//not sure if we need both at() and getElem()
template <class T, int N>
T& UpperTriangular<T, N>::at(const int& row, const int& col) {
	if (col < row) {
		throw out_of_range("Cannot modify elements below the diagonal.");
	}
	int index = Base::Sub2Ind(row, col);
	return Base::mData[index];
}

template <class T, int N>
const T& UpperTriangular<T, N>::at(const int& row, const int& col) const {
	return Base::at(row, col);
}

//*********************************************************
// row swapping - DOES NOT WORKKKKKKK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//*********************************************************
template <class T, int N>
void UpperTriangular<T, N>::SwapRows(const int& r1, const int& r2, int start, int end) {
	if (start < max(r1, r2)) start = max(r1, r2);
	Base::SwapRows(r1, r2, start, end);
}

//********************************************************************************
// transpose
//********************************************************************************
template <class T, int N>
LowerTriangular<T, N> UpperTriangular<T, N>::Transpose() const {
	LowerTriangular<T, N> res;
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			res.setElem(j, i, this->at(i, j));
		}
	}
	return res;
}

#endif // !UpperTriangular_HPP