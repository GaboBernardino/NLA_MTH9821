//Gabriele Bernardino - some helper functions for the NLA class

#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <cmath>
using namespace std;

//*****************************************************************************
// function to check whether floating point numbers are the same
//*****************************************************************************
template <typename T>
bool _close_enough(const T& a, const T& b, const T& tol = 1e-9) {
	return fabs(a - b) < tol;
}

#endif // !HELPERS_HPP