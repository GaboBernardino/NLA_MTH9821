//Gabriele Bernardino - Code to test my Matrix class

#include "Matrix.hpp"
#include "UpperTriangular.hpp"
#include "LowerTriangular.hpp"
#include "qbVector.hpp"
using namespace std;

int main() {

	cout << "********************* Testing constructors and printing *********************" << endl;
	Matrix<double, 3, 2> m1(1.5);
	cout << "m1:\n"; m1.Print();

	vector<int> v{ 3,2,4,6,11,-7,-1,10 };
	Matrix<int, 2, 4> m2(&v);
	cout << "m2:\n"; m2.Print();

	cout << "m2[1,3] = " << m2.getElem(1, 3) << endl;

	m1.setElem(0, 0, -100.);
	cout << "m1 after modifying:\n"; m1.Print();

	cout << "Using [] operator:\n";
	cout << "m2[0,2] = " << m2[{0, 2}] << endl;
	m1[{1, 1}] = -100;
	cout << "m1 after modifying:\n"; m1.Print();

	//cout << m2[{100, 10}];
	Matrix<int, 1, 2> sub2 = m2.subMatrix<1, 2>(0, 1);
	cout << "Submatrix m2[0,1:2]:\n"; sub2.Print();

	qbVector<double> r1 = m1.row(0);
	cout << "First row of m1: " << r1 << endl;

	qbVector<int> r2 = m2.row(1, 1, 3);
	cout << "Second row of m2, elems at [1, 3): " << r2 << endl;

	cout << "\n********************* Testing operators *********************" << endl;

	Matrix<double, 2, 2> m3(-1.);
	Matrix<double, 3, 2> prod = m1 * m3;
	cout << "m3:\n"; m3.Print();
	cout << "m1*m3:\n"; prod.Print();

	cout << "5+m2:\n"; (5 + m2).Print();

	qbVector<double> v1(3, 1.5);
	cout << "v1 = "; v1.Print();
	cout << "-2*v1 = " << (-2. * v1) << endl;

	cout << "\n********************* Testing some functions *********************" << endl;
	//determinant:
	vector<double> sqData{ 1,0,0,3,5,1,-7,2,3 };
	Matrix<double, 3, 3> sq(&sqData);
	cout << "Matrix:\n"; sq.Print();
	//not workingcout << "Determinant of the matrix = " << sq.Det() << endl; //not working

	sq.SwapRows(0, 2, 0, 2);
	cout << "After swapping rows 1 and 3 from [0, 2):\n"; sq.Print();

	vector<double> triData{7.,8.,9.,10.,0, 4., 5., 6.,0,0, 2., 3.,0, 0, 0,1.};
	//vector<double> triData{ 1,0,0,0,2,3,0,0,4,5,6,0,7,8,9,10 };
	UpperTriangular<double, 4> tri(&triData);
	cout << "Swapping rows of Triangular matrix:\n";
	tri.SwapRows(0, 2);
	tri.Print();


	cout << endl;
	system("PAUSE");
	return 0;
}