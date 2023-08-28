//Gabriele Bernardino - derived class for banded matrices

#ifndef BANDEDMATRIX_HPP
#define BANDEDMATRIX_HPP


#include "Matrix.hpp"
using namespace std;

template <class T, int NR = 1, int NC = 1, int Upper = 0, int Lower = 0>
class BandedMatrix : public Matrix<T, NR, NC> {
private:
	//using directives for more readable access to member data
	using Base = Matrix<T, NR, NC>;
	
	//int up = Upper;
	//int lo = Lower;
	//maybe there's a better way
public:
	//constructors and destructor
	BandedMatrix();
	BandedMatrix(const T& val);
	BandedMatrix(const T* inputData);
	BandedMatrix(const BandedMatrix<T, NR, NC>& inputMatrix);
	BandedMatrix(const vector<T>* inputData);

	//accessing elements
	T& operator [](const pair<int, int>& idx);
	const T& operator [](const pair<int, int>& index) const;
	T& at(const int& row, const int& col);
	const T& at(const int& row, const int& col) const;

	// Override the setElem function to ensure elements outside the triangular region are set to zero
	bool setElem(const int& row, const int& col, const T& element);

	//Can't swap rows here
	void SwapRows(const int& r1, const int& r2, int start = 0, int end = NC) = delete;

	//Override transpose to invert Upper and Lower
	BandedMatrix<T, NC, NR, Lower, Upper> Transpose() const;
};

/***********************************************************************************************************
Constructors and destructor
***********************************************************************************************************/
template <class T, int NR, int NC, int Upper, int Lower>
BandedMatrix<T, NR, NC, Upper, Lower>::BandedMatrix() : Base() {
	Base::mData = make_unique<T[]>(NR * NC);
	for (int i = 0; i < NR * NC; ++i) Base::mData[i] = T(0);
}

template <class T, int NR, int NC, int Upper, int Lower>
BandedMatrix<T, NR, NC, Upper, Lower>::BandedMatrix(const T& val) : Base() {
	Base::mData = make_unique<T[]>(NR * NC);
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			if ((j >= i && j - i <= Upper) || (j < i && i - j <= Lower))
				Base::mData[i * NC + j] = val;
			else Base::mData[i * NC + j] = T(0);
		}
	}
}

template <class T, int NR, int NC, int Upper, int Lower>
BandedMatrix<T, NR, NC, Upper, Lower>::BandedMatrix(const T* inputData) : Base() {
	if ((sizeof(inputData) / sizeof(inputData[0])) != NR * NC) {
		throw invalid_argument("Input size does not match.");
	}
	Base::mData = make_unique<T[]>(NR * NC);
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			if ((j >= i && j - i <= Upper) || (j < i && i - j <= Lower))
				Base::mData[i * NC + j] = inputData[i * NC + j];
			else Base::mData[i * NC + j] = T(0);
		}
	}
}

template <class T, int NR, int NC, int Upper, int Lower>
BandedMatrix<T, NR, NC, Upper, Lower>::BandedMatrix(const vector<T>* inputData) : Base() {
	if (inputData->size() != (NR * NC)) {
		throw invalid_argument("Input size does not match.");
	}
	Base::mData = make_unique<T[]>(NR * NC);
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			if ((j >= i && j - i <= Upper) || (j < i && i - j <= Lower))
				Base::mData[i * NC + j] = inputData->at(i * NC + j);
			else Base::mData[i * NC + j] = T(0);
		}
	}
}

//copy ctor
template <class T, int NR, int NC, int Upper, int Lower>
BandedMatrix<T, NR, NC, Upper, Lower>::BandedMatrix(const BandedMatrix<T, NR, NC>& inputMatrix) : Base() {
	Base::mData = make_unique<T[]>(NR * NC);
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			if ((j >= i && j - i <= Upper) || (j < i && i - j <= Lower))
				Base::mData[i * NC + j] = inputMatrix.mData[i * NC + j];
			else Base::mData[i * NC + j] = T(0);
		}
	}
}

//assignment operator
//..........

//*********************************************************
// element access
//*********************************************************
template <class T, int NR, int NC, int Upper, int Lower>
bool BandedMatrix<T, NR, NC, Upper, Lower>::setElem(const int& row, const int& col, const T& element) {
	if ((row < col && col - row > Upper) || (row > col && row - col > Lower)) {
		// Elements outside band should be zero
		return false;
	}
	return Base::setElem(row, col, element);
}

template <class T, int NR, int NC, int Upper, int Lower>
T& BandedMatrix<T, NR, NC, Upper, Lower>::operator [](const pair<int, int>& idx) {
	if ((idx.first < idx.second && idx.second - idx.first > Upper) || (idx.first > idx.second && idx.first - idx.second > Lower)) {
		throw out_of_range("Cannot modify elements outside of the band.");
	}
	int index = Base::Sub2Ind(idx.first, idx.second);
	return Base::mData[index];
}

template <class T, int NR, int NC, int Upper, int Lower>
const T& BandedMatrix<T, NR, NC, Upper, Lower>::operator [](const pair<int, int>& idx) const {
	return Base::operator[](idx);
}

//not sure if we need both at() and getElem()
template <class T, int NR, int NC, int Upper, int Lower>
T& BandedMatrix<T, NR, NC, Upper, Lower>::at(const int& row, const int& col) {
	if ((row < col && col - row > Upper) || (row > col && row - col > Lower)) {
		throw out_of_range("Cannot modify elements below the diagonal.");
	}
	int index = Base::Sub2Ind(row, col);
	return Base::mData[index];
}

template <class T, int NR, int NC, int Upper, int Lower>
const T& BandedMatrix<T, NR, NC, Upper, Lower>::at(const int& row, const int& col) const {
	return Base::at(row, col);
}

//*********************************************************
// row swapping - DOES NOT WORKKKKKKK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//*********************************************************
/*template <class T, int NR, int NC>
void BandedMatrix<T, NR, NC>::SwapRows(const int& r1, const int& r2, int start, int end) {
	if (start < max(r1, r2)) start = max(r1, r2);
	Base::SwapRows(r1, r2, start, end);
}*/

//********************************************************************************
// transpose
//********************************************************************************
template <class T, int NR, int NC, int Upper, int Lower>
BandedMatrix<T, NC, NR, Lower, Upper> BandedMatrix<T, NR, NC, Upper, Lower>::Transpose() const {
	BandedMatrix<T, NC, NR, Lower, Upper> res;
	for (int i = 0; i < NC; ++i) {
		for (int j = 0; j < NR; ++j) {
			res.setElem(i, j, this->at(j, i));
		}
	}
	return res;
}

#endif // !BANDEDMATRIX_HPP