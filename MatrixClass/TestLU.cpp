//Gabriele Bernardino - Code to test LU decomposition (chapter 2 of Dan Stefanica's NLA primer)

#include "LU_decomposition.hpp"
#include "StopWatch.hpp"
using namespace std;

int main() {

	StopWatch ap;

	cout << "********************* Testing Bwd and Fwd substitution *********************" << endl;
	vector<double> vFwd{2,0,0,1,2,0,1,1,2};
	LowerTriangular<double, 3> L(&vFwd);
	qbVector<double> bFwd{ 2,-1,0 };
	
	ap.Start();
	qbVector<double> xFwd = forward_subst(L, bFwd);
	ap.Stop();
	cout << "L:\n"; L.Print();
	cout << "b: " << bFwd << endl;
	cout << "x s.t. Lx=b: " << xFwd << endl;
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();

	vector<double> vBwd{ 2.,1.,1.,0.,2.,1.,0.,0.,2. };
	UpperTriangular<double, 3> U(&vBwd);
	qbVector<double> bBwd{ 9.0,4.0,4.0 };

	ap.Start();
	qbVector<double> xBwd = backward_subst(U, bBwd);
	ap.Stop();
	cout << "U:\n"; U.Print();
	cout << "b: " << bBwd << endl;
	cout << "x s.t. Ux=b: " << xBwd << endl;
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();

	cout << "BIDIAGONAL MATRICES:\n";
	BandedMatrix<double, 4, 4, 1, 0> Uband(3.);
	cout << "U:\n"; Uband.Print();
	qbVector<double> bBand{ 10.5,15,21.,12. };
	cout << "b: " << bBand << endl;

	ap.Start();
	qbVector<double> xBand = backward_subst(Uband, bBand);
	ap.Stop();
	cout << "x s.t. Ux=b: " << xBand << endl;
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();

	cout << "\n********************** LU decomposition - no pivoting **********************" << endl;
	vector<double> v1{ 4.,3.,6.,3 };
	Matrix<double, 2, 2> A1(&v1);
	cout << "A:\n"; A1.Print();

	ap.Start();
	pair<LowerTriangular<double, 2>, UpperTriangular<double, 2>> LU1 = LU_no_pivot(A1);
	cout << "L:\n"; LU1.first.Print();
	cout << "U:\n"; LU1.second.Print();
	ap.Stop();
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();

	cout << "LINEAR SOLVER:" << endl;
	vector<double> v2{ 2.,1.,1.,1.,2.,2.,1.,1.,2. };
	Matrix<double, 3, 3> A2(&v2);
	cout << "A:\n"; A2.Print();
	qbVector<double> b2{ 2.,4.,6. };
	cout << "b = " << b2 << endl;
	
	ap.Start();
	qbVector<double> solve1 = LU_solve_no_pivot(A2, b2);
	cout << "x (s.t. Ax=b) = " << solve1 << endl;
	ap.Stop();
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();

	cout << "BANDED MATRICES:\n";
	vector<double> bandData{ 1.,2.,0.,3.,4.,5.,0.,6.,7. };
	BandedMatrix<double, 3, 3, 1, 1> Aband(&bandData);
	cout << "A:\n"; Aband.Print();

	ap.Start();
	//LUpairDiag<double, 3> LUband = LU_no_pivot(Aband);
	LUpairBand<double, 3, 1> LUband = LU_no_pivot_band(Aband);
	cout << "L:\n"; LUband.first.Print();
	cout << "U:\n"; LUband.second.Print();
	ap.Stop();
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();
	cout << "Check: L*U =\n";
	(LUband.first * LUband.second).Print();
	

	cout << "\n********************** LU decomposition **********************" << endl;

	vector<double> v3{ 2.,1.,1.,0.,4.,3.,3.,1.,8.,7.,9.,5.,6.,7.,9.,8. };
	Matrix<double, 4, 4> m(&v3);
	cout << "A:\n"; m.Print();

	ap.Start();
	auto PLU = LU_row_pivot(m);
	Matrix<double, 4, 4> P = get<0>(PLU);
	LowerTriangular<double, 4> Ldec = get<1>(PLU);
	UpperTriangular<double, 4> Udec = get<2>(PLU);
	cout << "P:\n"; P.Print();
	cout << "L:\n"; Ldec.Print();
	cout << "U:\n"; Udec.Print();
	ap.Stop();
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();

	cout << "Check: P*L*U:\n"; (P * Ldec * Udec).Print();

	cout << "\nLINEAR SOLVER:\n";

	vector<double> lin3Data{ 2.,3.,-1.,4.,-2.,1.,-3.,2.,-2. };
	Matrix<double, 3, 3> A3(&lin3Data);
	cout << "A:\n"; A3.Print();
	qbVector<double> b3{ 5.,1.,-3. };
	cout << "b = " << b3 << endl;

	ap.Start();
	qbVector<double> solve2 = LU_solve(A3, b3);
	cout << "x (s.t. Ax=b) = " << solve2 << endl;
	ap.Stop();
	cout << "(Elapsed time: " << ap.GetTime() << "s)\n\n"; ap.Reset();
	
	cout << "Check: A*x = " << (A3 * solve2) << endl;

	cout << endl;
	system("PAUSE");
	return 0;
}