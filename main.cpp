//Gabriele Bernardino - Code to test Cholesky decomposition based on Dan Stefanica's NLA Primer

#include "Cholesky.hpp"
#include "StopWatch.hpp"
#include "BandedMatrix.hpp"
#include "CubicSpline.hpp"
using namespace std;

int main() {

	StopWatch patek; //to measure time

	cout << "******************* Testin' the transpose function *******************" << endl;
	vector<int> data1{ 1,2,3,4,5,6 };
	Matrix<int, 2, 3> m1(&data1);
	m1.Print();
	cout << "Transposed:\n"; m1.Transpose().Print();

	cout << "Triangular matrices:\n";
	vector<int> data2{ 1,0,0,2,3,0,4,5,6 };
	LowerTriangular<int, 3> m2(&data2);
	m2.Print();
	cout << "Transposed:\n"; m2.Transpose().Print();

	cout << "\n******************* Testin' Cholesky decomposition *******************" << endl;
	vector<double> data3{ 4,12,-16,12,37,-43,-16,-43,98 };
	Matrix<double, 3, 3> A(&data3);
	cout << "A:\n"; A.Print();

	patek.Start();
	auto [U, check] = Cholesky(A); //UpperTriangular<double, 3>, int
	if (check == 3) {
		cout << "U (s.t. U'U=A):\n"; U.Print();
	}
	else cout << "Decomposition failed at row " << check << endl;
	patek.Stop();

	cout << "\n(Elapsed time: " << patek.GetTime() << "s)\n\n"; patek.Reset();


	cout << "Check: U'U = \n"; (U.Transpose() * U).Print();

	cout << "\n******************* Testin' Cholesky linear solver *******************" << endl;
	vector<double> data4{ 2.,3.,-1.,4.,-2.,1.,-3.,2.,-2. };
	Matrix<double, 3, 3> A4(&data4);
	cout << "A:\n"; A4.Print();
	qbVector<double> b4{ 5.,1.,-3. };
	cout << "b = " << b4 << endl;

	patek.Start();
	qbVector<double> solve4 = Cholesky_solve(A, b4);
	cout << "x (s.t. Ax=b) = " << solve4 << endl;
	patek.Stop();

	cout << "\n(Elapsed time: " << patek.GetTime() << "s)"; patek.Reset();
	//compare with LU:
	patek.Start();
	qbVector<double> solveLU = LU_solve(A4, b4);
	patek.Stop();
	cout << "\n(Elapsed time, LU solver: " << patek.GetTime() << "s)\n\n"; patek.Reset();

	cout << "Check: Ax = " << (A * solve4) << endl;

	cout << "\n******************* Testin' Banded matrices *******************" << endl;
	vector<int> data5{ 1,2,0,3,4,5,0,6,7,0,0,8 };
	BandedMatrix<int, 4, 3, 1, 1> m5(&data5);
	m5.Print();
	cout << "Transposed:\n"; m5.Transpose().Print();

	vector<double> bandData{ 5.,-1.,0.,-1.,4.,-1.,0.,-1.,4. };
	BandedMatrix<double, 3, 3, 1, 1> Aband(&bandData);
	cout << "A:\n"; Aband.Print();

	patek.Start();
	auto [Uband, checkBand] = Cholesky(Aband);
	if (checkBand) {
		cout << "U (s.t. U'U=A):\n"; Uband.Print();
	}
	else cout << "Decomposition failed." << endl;
	patek.Stop();

	cout << "\n(Elapsed time: " << patek.GetTime() << "s)\n\n"; patek.Reset();
	cout << "Check: U'U = \n"; (Uband.Transpose() * Uband).Print();

	cout << "\n**************** NLA summer HW - Exercise 21 ****************" << endl;
	vector<double> HWdata{ 1, 0.2, 0.3, 0.2, 1, -0.2, 0.3, -0.2, 1 };
	Matrix<double, 3, 3> HWm(&HWdata);
	auto [cholHW, HWcheck] = Cholesky(HWm);
	if (HWcheck == 3) {
		cout << "U =\n";
		cholHW.Print();
		cout << "Check U'U =\n";
		(cholHW.Transpose() * cholHW).Print();
	}
	
	cout << endl;
	system("PAUSE");
	return 0;
}