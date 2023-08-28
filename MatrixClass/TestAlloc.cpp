//Gabriele Bernardino - Code to test Asset Allocation functions

#include "AssetAlloc.hpp"
#include "CubicSpline.hpp"
using namespace std;

int main() {

	cout << "****************** Minimum variance portfolio ******************" << endl;
	vector<double> sigData{ .0105,.0028,.0039,.002,.0025,.0028,.008,.0034,.002,.0018,.0039,.0034,.0175,.003,.0028,.002,.002,.003,.0084,.002,.0025,.0018,.0028,.002,.0074 };
	Matrix<double, 5, 5> sig(&sigData);
	qbVector<double> ret{ .2219,-.0925,.5672,.1327,.2811 };
	auto [w, cash, std] = minimumVariance(sig, ret, 0., .1784);
	cout << w << endl;

	//********************************************************************
	// Testing cubic spline
	//********************************************************************
	cout << "\n****************** Cubic spline interpolation ******************" << endl;
	qbVector<double> x{ 0., 2. / 12., 0.5, 1., 20. / 12. };
	qbVector<double> v{ 0.005, 0.0065, 0.0085, 0.0105, 0.012 };
	auto coeff = CubicSpline<double, 4>(x, v);
	print_coefficients(coeff);
	

	system("PAUSE");
	return 0;
}