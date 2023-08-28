//Gabriele Bernardino - Matrix class for linear algebra stuff

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <iomanip>
#include <functional>
#include <memory>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <cassert>
#include "qbVector.hpp"
#include "helpers.hpp"
using namespace std;

//da aggiungere: assignment operator, matrix*vector

template <class T, int NR = 1, int NC = 1> class Matrix {
protected:
	unique_ptr<T[]> mData; //unique ptr to 1d built-in array
	int nElements = NR * NC; //useful in for loops

	//helpers
	int Sub2Ind(const int& row, const int& col) const;
	bool IsSquare() const;
public:
	//constructors and destructor
	Matrix();
	Matrix(const T& val);
	Matrix(const T* inputData);
	Matrix(const vector<T>* inputData);
	//Matrix(csvfile inputMatrix);
	Matrix(const Matrix<T, NR, NC>& inputMatrix);
	//~Matrix();

	//write assignment operator
	//Matrix<T, NR, NC>& operator = (const Matrix<T, NR, NC>& other);

	//Encapsulation
	T getElem(const int& row, const int& col) const;
	bool setElem(const int& row, const int& col, const T& element);
	int nRows() const;
	int nCols() const;

	//accessing elements
	T& operator [](const pair<int, int>& idx);
	const T& operator [](const pair<int, int>& index) const;
	T& at(const int& row, const int& col);
	const T& at(const int& row, const int& col) const;

	qbVector<T> row(const int& idx, int start = 0, int end = NC);
	//const vector<T>& row(const int& idx) const;
	//vector<T>& row(Matrix<T, NR, NC>& matrix, const int& idx);
	//const vector<T>& row(const Matrix<T, NR, NC>& matrix, const int& idx);
	qbVector<T> col(const int& idx, int start = 0, int end = NR);
	//const vector<T>& col(const int& idx) const;

	template <int nr, int nc>
	Matrix<T, nr, nc> subMatrix(const int& startRow, const int& startCol);
	//const Matrix<T>& subMatrix(const int& row, const int& idx) const;

	//computations
	T Det() const;

	//printing
	void Print(const size_t& prec = 4) const;

	//configuration methods
	//bool Resize(const int& numRows, const int& numCols); //should make it template and create and return a new matrix
	void SetToId();
	Matrix<T, NR - 1, NC - 1> remove(const int& row, const int& col) const; //this is to create a submatrix without a row and a col (useful for determinant)
	void SwapRows(const int& r1, const int& r2, int start = 0, int end = NC);
	Matrix<T, NC, NR> Transpose() const;

	//operator overloading (bool)
	bool operator== (const Matrix<T, NR, NC>& rhs);
	bool Compare(const Matrix<T, NR, NC>& matrix1, double tol) const;

	//operator overloading (arithmetics)
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator+ (const Matrix<U, NR, NC>& lhs, const Matrix<U, NR, NC>& rhs);
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator+ (const U& lhs, const Matrix<U, NR, NC>& rhs);
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator+ (const Matrix<U, NR, NC>& lhs, const U& rhs);

	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator- (const Matrix<U, NR, NC>& lhs, const Matrix<U, NR, NC>& rhs);
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator- (const U& lhs, const Matrix<U, NR, NC>& rhs);
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator- (const Matrix<U, NR, NC>& lhs, const U& rhs);

	template <class U, int NR, int NC, int NC2> friend Matrix<U, NR, NC2> operator* (const Matrix<U, NR, NC>& lhs, const Matrix<U, NC, NC2>& rhs);
	template <class U, int NR, int NC> friend qbVector<U> operator* (const qbVector<U>& lhs, const Matrix<U, NR, NC>& rhs);
	template <class U, int NR, int NC> friend qbVector<U> operator* (const Matrix<U, NR, NC>& lhs, const qbVector<U>& rhs);
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator* (const U& lhs, const Matrix<U, NR, NC>& rhs);
	template <class U, int NR, int NC> friend Matrix<U, NR, NC> operator* (const Matrix<U, NR, NC>& lhs, const U& rhs);
};

/***********************************************************************************************************
Constructors and destructor
***********************************************************************************************************/
template <class T, int NR, int NC> Matrix<T, NR, NC>::Matrix() {
	mData = make_unique<T[]>(nElements);
	for (int i = 0; i < nElements; ++i) mData[i] = T(0);
}

template <class T, int NR, int NC> Matrix<T, NR, NC>::Matrix(const T& val) {
	mData = make_unique<T[]>(nElements);
	for (int i = 0; i < nElements; ++i) mData[i] = val;
}

template <class T, int NR, int NC> Matrix<T, NR, NC>::Matrix(const T* inputData) {
	if ((sizeof(inputData) / sizeof(inputData[0])) != nElements) {
		throw invalid_argument("Input size does not match.");
	}
	mData = make_unique<T[]>(nElements);
	for (int i = 0; i < nElements; ++i) mData[i] = inputData[i];
}

template <class T, int NR, int NC> Matrix<T, NR, NC>::Matrix(const vector<T>* inputData) {
	if (inputData->size() != (NR * NC)) {
		throw invalid_argument("Input size does not match.");
	}
	mData = make_unique<T[]>(nElements);
	for (int i = 0; i < nElements; ++i) mData[i] = inputData->at(i);
}

//copy ctor
template <class T, int NR, int NC> Matrix<T, NR, NC>::Matrix(const Matrix<T, NR, NC>& inputMatrix) {
	mData = make_unique<T[]>(nElements);
	for (int i = 0; i < nElements; ++i) mData[i] = inputMatrix.mData[i];
}

//assignment operator - GOTTA FIX
/*template <class T, int NR, int NC>
Matrix<T, NR, NC>& Matrix<T, NR, NC>::operator = (const Matrix<T, NR, NC>& other) {
	if (this == &) return *this;
	mData = other.mData;
	nElements = other.nElements;
	return *this;
}*/

/***********************************************************************************************************
Getters and setters
***********************************************************************************************************/
template <class T, int NR, int NC>
T Matrix<T, NR, NC>::getElem(const int& row, const int& col) const {
	//dont need size check here because Sub2Ind does it
	int idx = Sub2Ind(row, col);
	return mData[idx];
}

template <class T, int NR, int NC>
bool Matrix<T, NR, NC>::setElem(const int& row, const int& col, const T& element) {
	int idx = Sub2Ind(row, col);
	mData[idx] = element;
	return true;
}

template <class T, int NR, int NC>
int Matrix<T, NR, NC>::nRows() const {
	return NR;
}

template <class T, int NR, int NC>
int Matrix<T, NR, NC>::nCols() const {
	return NC;
}

/***********************************************************************************************************
accessing elements
***********************************************************************************************************/
template <class T, int NR, int NC>
T& Matrix<T, NR, NC>::operator [](const pair<int, int>& idx) {
	int index = Sub2Ind(idx.first, idx.second);
	return mData[index];
}

template <class T, int NR, int NC>
const T& Matrix<T, NR, NC>::operator [](const pair<int, int>& idx) const {
	int index = Sub2Ind(idx.first, idx.second);
	return mData[index];
}

//not sure if we need both at() and getElem()
template <class T, int NR, int NC>
T& Matrix<T, NR, NC>::at(const int& row, const int& col) {
	int index = Sub2Ind(row, col);
	return mData[index];
}

template <class T, int NR, int NC>
const T& Matrix<T, NR, NC>::at(const int& row, const int& col) const {
	int index = Sub2Ind(row, col);
	return mData[index];
}

//rows, columns, submatrices:
/*template <class T, int NR, int NC>
vector<T>& Matrix<T, NR, NC>::row(const int& idx) {
	if (idx >= NR || idx < 0) {
		throw out_of_range("Index must be non-negative");
	}
	return *reinterpret_cast<const vector<T>*>(mData.get() + idx * NC);
}

template <class T, int NR, int NC>
const vector<T>& Matrix<T, NR, NC>::row(const int& idx) const {
	if (idx >= NR || idx < 0) {
		throw out_of_range("Index must be non-negative");
	}
	return *reinterpret_cast<const vector<T>*>(mData.get() + idx * NC);
}

template <class T, int NR, int NC>
vector<T>& Matrix<T, NR, NC>::row(Matrix<T, NR, NC>& matrix, const int& idx) {
	if (idx < 0 || idx >= NR) {
		throw std::out_of_range("Invalid row index");
	}
	return matrix.row(idx);
}

template <class T, int NR, int NC>
const vector<T>& Matrix<T, NR, NC>::row(const Matrix<T, NR, NC>& matrix, const int& idx) {
	if (idx < 0 || idx >= NR) {
		throw std::out_of_range("Invalid row index");
	}
	return matrix.row(idx);
}*/
template <class T, int NR, int NC>
qbVector<T> Matrix<T, NR, NC>::row(const int& idx, int start, int end) {
	if (idx >= NR || idx < 0) {
		throw out_of_range("Index must be non-negative");
	}
	if (end > NC) end = NC; if (start < 0) start = 0;

	qbVector<T> res(end - start);
	for (int i = start; i < end; ++i) res[i - start] = mData[idx * NC + i];
	return res;
}

//funzione per assegnare un vettore a una riga?

template <class T, int NR, int NC>
qbVector<T> Matrix<T, NR, NC>::col(const int& idx, int start, int end) {
	if (idx >= NC || idx < 0) {
		throw out_of_range("Index must be non-negative");
	}
	if (end > NR) end = NR; if (start < 0) start = 0;

	qbVector<T> res(end - start);
	for (int i = start; i < end; ++i) res[i] = mData[idx + NR * i];
	return res;
}

//funzione per assegnare un vettore a una colonna?

template <class T, int NR, int NC> template <int nr, int nc>
Matrix<T, nr, nc> Matrix<T, NR, NC>::subMatrix(const int& startRow, const int& startCol) {
	if ((startRow + nr >= NR) || (startRow < 0) || (startCol + nc >= NC) || (startCol < 0)) {
		throw out_of_range("Indexes must be non-negative");
	}
	Matrix<T, nr, nc> res;
	for (int i = 0; i < nr; ++i)
		for (int j = 0; j < nc; ++j)
			res[{i, j}] = this->at(startRow + i, startCol + j);
	return res;
}

/***********************************************************************************************************
Operator overloading
***********************************************************************************************************/

template <class T, int NR, int NC>
bool Matrix<T, NR, NC>::operator==(const Matrix<T, NR, NC>& rhs) {
	bool flag = true;
	for (int i = 0; i < this->nElements; ++i) {
		if (!_close_enough(this->mData[i], rhs.mData[i]))
			flag = false;
		break;
	}
	return flag;
}

//matrix+matrix; assuming user uses same-sized matrices
template <class T, int NR, int NC>
Matrix<T, NR, NC> operator+(const Matrix<T, NR, NC>& lhs, const Matrix<T, NR, NC>& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs.mData[i] + rhs.mData[i];
	}

	return result;
}

//scalar + matrix
template <class T, int NR, int NC>
Matrix<T, NR, NC> operator+(const T& lhs, const Matrix<T, NR, NC>& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; ++i) {
		result.mData[i] = lhs + rhs.mData[i];
	}

	return result;
}

//matrix + scalar
template <class T, int NR, int NC>
Matrix<T, NR, NC> operator+(const Matrix<T, NR, NC>& lhs, const T& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs.mData[i] + rhs;
	}

	return result;
}

//matrix-matrix; assuming user uses same-sized matrices
template <class T, int NR, int NC> Matrix<T, NR, NC> operator-(const Matrix<T, NR, NC>& lhs, const Matrix<T, NR, NC>& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs.mData[i] - rhs.mData[i];
	}

	return result;
}

//scalar - matrix
template <class T, int NR, int NC> Matrix<T, NR, NC> operator-(const T& lhs, const Matrix<T, NR, NC>& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs - rhs.mData[i];
	}

	return result;
}

//matrix - scalar
template <class T, int NR, int NC> Matrix<T, NR, NC> operator-(const Matrix<T, NR, NC>& lhs, const T& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs.mData[i] - rhs;
	}

	return result;
}

//matrix*matrix; assuming user uses same-sized matrices
//************* SHOULD PROBABLY MAKE MORE EFFICIENT *******************
template <class T, int NR, int NC, int NC2>
Matrix<T, NR, NC2> operator*(const Matrix<T, NR, NC>& lhs, const Matrix<T, NC, NC2>& rhs) {
	Matrix<T, NR, NC2> res;
	T cumSum;
	//loop through each lhs row
	for (int lhsR = 0; lhsR < NR; ++lhsR) {
		//loop through each rhs column
		for (int rhsC = 0; rhsC < NC2; ++rhsC) {
			cumSum = T(0);
			//loop through each element of this lhs row
			for (int lhsC = 0; lhsC < NC; ++lhsC) {
				cumSum += lhs.at(lhsR, lhsC) * rhs.at(lhsC, rhsC);
			}
			res.at(lhsR, rhsC) = cumSum;
		}
	}

	return res;
}

//vector*matrix
template <class T, int NR, int NC> qbVector<T> operator*(const qbVector<T>& lhs, const Matrix<T, NR, NC>& rhs) {
	//static_assert(lhs.size() == NR, "Dimensions do not match.");

	qbVector<T> res(NC);
	T cumSum;

	//loop thru matrix columns
	for (int i = 0; i < NC; ++i) {
		cumSum = T(0);
		//loop thru each column's rows
		for (int j = 0; j < NR; ++j) {
			cumSum += lhs.at(j) * rhs.at(j, i);
		}
		res.at(i) = cumSum;
	}
	return res;
}

//matrix*vector
template <class T, int NR, int NC> qbVector<T> operator*(const Matrix<T, NR, NC>& lhs, const qbVector<T>& rhs) {
	//static_assert(rhs.size() == NC, "Dimensions do not match.");

	qbVector<T> res(NR);
	T cumSum;

	//loop thru matrix rows
	for (int i = 0; i < NR; ++i) {
		cumSum = T(0);
		//loop thru each row's columns
		for (int j = 0; j < NC; ++j) {
			cumSum += lhs.at(i, j) * rhs.at(j);
		}
		res.at(i) = cumSum;
	}
	return res;
}

//scalar * matrix
template <class T, int NR, int NC> Matrix<T, NR, NC> operator*(const T& lhs, const Matrix<T, NR, NC>& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs * rhs.mData[i];
	}

	return result;
}

//matrix * scalar
template <class T, int NR, int NC> Matrix<T, NR, NC> operator*(const Matrix<T, NR, NC>& lhs, const T& rhs) {
	Matrix<T, NR, NC> result;
	int nElements = NR * NC;
	for (int i = 0; i < nElements; i++) {
		result.mData[i] = lhs.mData[i] * rhs;
	}

	return result;
}

/***********************************************************************************************************
Computations
***********************************************************************************************************/
/*template <class T, int NR, int NC>
T Matrix<T, NR, NC>::Det() const {
	if (!(this->IsSquare())) {
		throw invalid_argument("Matrix must be square");
	}
	if (NR == 1) return mData[0];
	if (NR == 2) return this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0);

	else {
		T det = T(0);
		T sign = T(1);
		for (int j = 0; j < NC; ++j) {
			Matrix<T, NR, NC> sub = this->remove(0, j);
			det += (sign * this->at(0, j) * sub.Det());
			sign *= -1;
		}
		return det;
	}
}*/
//DOES NOT WORK, MAYBE TRY USING GAUSSIAN ELIMINATION
//Is it really needed for numerical stuff tho?

/***********************************************************************************************************
Configuration
***********************************************************************************************************/
template <class T, int NR, int NC>
void Matrix<T, NR, NC>::SetToId() {
	if (!(this->IsSquare())) {
		throw invalid_argument("Matrix must be square");
	}
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			if (i == j) this->setElem(i, j, T(1));
			else this->setElem(i, j, T(0));
		}
	}
}

template <class T, int NR, int NC>
Matrix<T, NR - 1, NC - 1> Matrix<T, NR, NC>::remove(const int& row, const int& col) const {
	//static_assert(nr >= 0 && nr <= NR - 1, "Invalid number of rows for submatrix");
	//static_assert(nc >= 0 && nc <= NC - 1, "Invalid number of columns for submatrix");

	Matrix<T, NR - 1, NC - 1> sub;
	int rowNum = 0, colNum = 0;
	for (int i = 0; i < NR; ++i) {
		if (i == row) continue; //ignore the row to eliminate
		colNum = 0;
		for (int j = 0; j < NC; ++j) {
			if (j == col) continue; //ignore the column to eliminate
			sub.setElem(rowNum, colNum, this->at(i, j));
			++colNum;
		}
		++rowNum;
	}
	return sub;
}

template <class T, int NR, int NC>
void Matrix<T, NR, NC>::SwapRows(const int& r1, const int& r2, int start, int end) {
	//start is included, end is not
	if (r1 == r2) return;
	//if (end > NC) end = NC; if (start < 0) start = 0;
	qbVector<T> temp = this->row(r1, start, end);
	for (int i = start; i < end; ++i) {
		this->setElem(r1, i, this->at(r2, i));
		this->setElem(r2, i, temp[i - start]);
	}
}

template <class T, int NR, int NC>
Matrix<T, NC, NR> Matrix<T, NR, NC>::Transpose() const {
	Matrix<T, NC, NR> res;
	for (int i = 0; i < NC; ++i) {
		for (int j = 0; j < NR; ++j) {
			res.setElem(i, j, this->at(j, i));
		}
	}
	return res;
}

/***********************************************************************************************************
Printing
***********************************************************************************************************/
template <class T, int NR, int NC>
void Matrix<T, NR, NC>::Print(const size_t& prec) const {
	for (int i = 0; i < NR; ++i) {
		for (int j = 0; j < NC; ++j) {
			cout.precision(prec);
			cout << setw(prec * 2) << this->getElem(i, j);
		}
		cout << endl;
	}
}

//also overload ostream operator <<

/***********************************************************************************************************
Helpers
***********************************************************************************************************/
template <class T, int NR, int NC>
int Matrix<T, NR, NC>::Sub2Ind(const int& row, const int& col) const {
	if ((row >= NR) || (row < 0) || (col >= NC) || (col < 0)) {
		throw out_of_range("Indexes must be non-negative");
	}
	return (row * NC) + col; //turns the subscript into the index of the linear array
}

template <class T, int NR, int NC>
bool Matrix<T, NR, NC>::IsSquare() const {
	return (NR == NC);
}

#endif // !MATRIX_HPP