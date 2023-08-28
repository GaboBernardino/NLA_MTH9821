/*
Template class for vector specific operations

From Linear Algebra in C++ - Part 4 - A simple Vector class by QuantitativeBytes
*/

#ifndef QBVECTOR_HPP
#define QBVECTOR_HPP

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <initializer_list>
#include <vector>
using namespace std;

template <class T> class qbVector {
public:
    //constructors and destructor
    qbVector();
    qbVector(const size_t& s);
    qbVector(const size_t& s, const T& val);
    qbVector(const vector<T>& inputData);
    qbVector(initializer_list<T> values);
    ~qbVector();

    //Encapsulation
    T getElem(int index) const;
    bool setElem(int index, T elem);
    size_t size() const;

    //accessing elements
    T& operator [](const int& idx);
    const T& operator [](const int& index) const;
    T& at(const int& idx);
    const T& at(const int& idx) const;

    //overload operators
    qbVector<T> operator+ (const qbVector<T>& rhs) const;
    qbVector<T> operator- (const qbVector<T>& rhs) const;
    qbVector<T> operator* (const T& rhs) const; //vector*scalar product

    template <class U> friend qbVector<U> operator* (const U& lhs, const qbVector<U>& rhs); //scalar*vector product

    template <class U> friend qbVector<U> operator+ (const U& lhs, const qbVector<U>& rhs);
    template <class U> friend qbVector<U> operator+ (const qbVector<U>& lhs, const U& rhs);
    template <class U> friend qbVector<U> operator- (const U& lhs, const qbVector<U>& rhs);
    template <class U> friend qbVector<U> operator- (const qbVector<U>& lhs, const U& rhs);

    //static functions (wtf do they do?) - can be instantiated w/o reference to specific instance of class
    static T dot(const qbVector<T>& a, const qbVector<T>& b);
    static qbVector<T> cross(const qbVector<T>& a, const qbVector<T>& b);

    void Print();
    template <class U> friend ostream& operator<<(ostream& stream, const qbVector<T>& v);
    T Norm();
    T Sum();

private:
    vector<T> m_vectorData;
    size_t m_nDims;

};

/***********************************************************************************************************
Constructor and destructor
***********************************************************************************************************/
template <class T> qbVector<T>::qbVector() : m_vectorData(vector<T>()), m_nDims(0) {}

template <class T> qbVector<T>::qbVector(const size_t& s) : m_vectorData(vector<T>(s)), m_nDims(s) {}

template <class T> qbVector<T>::qbVector(const size_t& s, const T& val) : m_vectorData(vector<T>(s, val)), m_nDims(s) {}

template <class T> qbVector<T>::qbVector(const vector<T>& inputData) : m_vectorData(inputData), m_nDims(inputData.size()) {}

template <class T> qbVector<T>::qbVector(initializer_list<T> values) : m_vectorData(values), m_nDims(values.size()) {}

template <class T> qbVector<T>::~qbVector() {
    //stl vector class already cleans up
}

/***********************************************************************************************************
Element functions
***********************************************************************************************************/
template <class T> T qbVector<T>::getElem(int index) const {
    return m_vectorData[index];
}

template <class T> bool qbVector<T>::setElem(int index, T elem) {
    if (index >= 0 && index < m_nDims) {
        m_vectorData.at(index) = elem;
        return true;
    }
    else return false;
}

template <class T> size_t qbVector<T>::size() const {
    return m_nDims;
}

/***********************************************************************************************************
accessing elements
***********************************************************************************************************/
template <class T>
T& qbVector<T>::operator [](const int& idx) {
    return m_vectorData[idx];
}

template <class T>
const T& qbVector<T>::operator [](const int& idx) const {
    return m_vectorData[idx];
}

template <class T>
T& qbVector<T>::at(const int& idx) {
    return m_vectorData[idx];
}

template <class T>
const T& qbVector<T>::at(const int& idx) const {
    return m_vectorData[idx];
}


/***********************************************************************************************************
printing
***********************************************************************************************************/
template <class T>
ostream& operator<<(ostream& stream, const qbVector<T>& v) {
    stream << "[";
    for (int i = 0; i < v.size(); i++) {
        stream << v.getElem(i);
        if (i < v.size() - 1) stream << ", ";
    }
    stream << "]";
    return stream;
}

template <class T> void qbVector<T>::Print() {
    int nDim = this->size();
    cout << "[";
    for (int i = 0; i < nDim; i++) {
        //cout.precision(3);
        cout << this->getElem(i);
        if (i < nDim - 1) cout << ", ";
    }
    cout << "]\n";
}

/***********************************************************************************************************
computations
***********************************************************************************************************/
template <class T> T qbVector<T>::Norm() {
    int nDim = this->size();
    T cumSum = T(0);
    for (int i = 0; i < nDim; i++) {
        cumSum += pow(this->getElem(i), 2);
    }
    return sqrt(cumSum);
}

template <class T> T qbVector<T>::Sum() {
    int nDim = this->size();
    T cumSum = T(0);
    for (int i = 0; i < nDim; i++) {
        cumSum += this->getElem(i);
    }
    return cumSum;
}

/***********************************************************************************************************
Overload operators
***********************************************************************************************************/
template <class T>
qbVector<T> qbVector<T>::operator+(const qbVector<T>& rhs) const {
    //check same size
    if (m_nDims != rhs.m_nDims)
        throw invalid_argument("Vector dimensions must match.");

    vector<T> result;
    for (int i = 0; i < m_nDims; ++i)
        result.push_back(m_vectorData.at(i) + rhs.m_vectorData.at(i));
    return qbVector<T>(result);
}

template <class T>
qbVector<T> qbVector<T>::operator-(const qbVector<T>& rhs) const {
    //check same size
    if (m_nDims != rhs.m_nDims)
        throw invalid_argument("Vector dimensions must match.");

    vector<T> result;
    for (int i = 0; i < m_nDims; ++i)
        result.push_back(m_vectorData.at(i) - rhs.m_vectorData.at(i));
    return qbVector<T>(result);
}

template <class T>
qbVector<T> qbVector<T>::operator*(const T& rhs) const {
    vector<T> result;
    for (int i = 0; i < m_nDims; ++i)
        result.push_back(m_vectorData.at(i) * rhs);
    return qbVector<T>(result);
}

template <class T>
qbVector<T> operator*(const T& lhs, const qbVector<T>& rhs) {
    vector<T> result;
    for (int i = 0; i < rhs.m_nDims; ++i)
        result.push_back(lhs * rhs.m_vectorData.at(i));
    return qbVector<T>(result);
}

template <class T>
qbVector<T> operator+(const T& lhs, const qbVector<T>& rhs) {
    vector<T> result;
    for (int i = 0; i < rhs.m_nDims; ++i)
        result.push_back(lhs + rhs.m_vectorData.at(i));
    return qbVector<T>(result);
}

template <class T>
qbVector<T> operator+(const qbVector<T>& lhs, const T& rhs) {
    vector<T> result;
    for (int i = 0; i < lhs.m_nDims; ++i)
        result.push_back(lhs.m_vectorData.at(i) + rhs);
    return qbVector<T>(result);
}

template <class T>
qbVector<T> operator-(const T& lhs, const qbVector<T>& rhs) {
    vector<T> result;
    for (int i = 0; i < rhs.m_nDims; ++i)
        result.push_back(lhs - rhs.m_vectorData.at(i));
    return qbVector<T>(result);
}

template <class T>
qbVector<T> operator-(const qbVector<T>& lhs, const T& rhs) {
    vector<T> result;
    for (int i = 0; i < lhs.m_nDims; ++i)
        result.push_back(lhs.m_vectorData.at(i) - rhs);
    return qbVector<T>(result);
}


/***********************************************************************************************************
Static functions
***********************************************************************************************************/
template <class T> T qbVector<T>::dot(const qbVector<T>& a, const qbVector<T>& b) {
    //check same size
    if (a.m_nDims != b.m_nDims)
        throw invalid_argument("Vector dimensions must match.");

    T cumSum = T(0.0); //no. specified as a double but must be whatever type T is
    for (int i = 0; i < a.m_nDims; ++i)
        cumSum += a.m_vectorData.at(i) * b.m_vectorData.at(i);

    return cumSum;
}

template <class T> qbVector<T> qbVector<T>::cross(const qbVector<T>& a, const qbVector<T>& b) {
    //check same size
    if (a.m_nDims != b.m_nDims)
        throw invalid_argument("Vector dimensions must match.");

    //also, cross product can only be computed for 3 (or 7) dim vectors
    if (a.m_nDims != 3)
        throw invalid_argument("Cross product can only be computed for 3D vectors.");

    vector<T> result;
    result.push_back(a.m_vectorData.at(1) * b.m_vectorData.at(2) - a.m_vectorData.at(2) * b.m_vectorData.at(1));
    result.push_back(-(a.m_vectorData.at(0) * b.m_vectorData.at(2) - a.m_vectorData.at(2) * b.m_vectorData.at(0))); //watch out for the -
    result.push_back(a.m_vectorData.at(0) * b.m_vectorData.at(1) - a.m_vectorData.at(1) * b.m_vectorData.at(0));

    return result;
}

#endif // QBVECTOR_HPP