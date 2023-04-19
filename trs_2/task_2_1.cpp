#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

constexpr auto pi = 3.14159265358979323846;

using namespace std;

void tridiagonal_matrix(const int n, const double* A, const double* B, const double* C, const double* D, double* X) {
	double* P = new double[n - 1];
	double* Q = new double[n];

	// forward
	P[0] = C[0] / B[0];
	Q[0] = D[0] / B[0];
	for (int i = 1; i < n; ++i) {
		if (i < n - 1) P[i] = C[i] / (B[i] - A[i] * P[i - 1]);
		Q[i] = (D[i] - A[i] * Q[i - 1]) / (B[i] - A[i] * P[i - 1]);
	}

	//backward
	X[n - 1] = Q[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		X[i] = Q[i] - P[i] * X[i + 1];
	}
	delete[] P;
	delete[] Q;
}

double sol(double t, double x)
{
	return (1 - exp(-t)) * sin(2 * pi * x) + x * x;
}

int main()
{
	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_2\\out\\u1_out.csv");

	vector <int> m = { 10};
	vector <int> n = { 540 };


	for (int k = 0; k < m.size(); k++)
	{
		double M = 0, diff;
		vector < vector <double> > u(n[k], vector <double>(m[k]));

		double tau = 1. / (n[k]);
		double h = 1. / (m[k]);

		for (int i = 0; i < m[k]; i++)
		{
			u[0][i] = (i * i * h * h);
		}

		double* A = new double[m[k]];
		double* B = new double[m[k]];
		double* C = new double[m[k]];
		double* F = new double[m[k]];
		double* U_next = new double[m[k]];

		for (int i = 1; i < m[k] - 1; i++)
		{
			A[i] = -1. / (4 * pi * pi * h * h);
			B[i] = 1. / tau + 1. / (2 * pi * pi * h * h);
			C[i] = -1. / (4 * pi * pi * h * h);
		}

		A[0] = 0;
		B[0] = -1. / h;
		C[0] = 1. / h;
		A[m[k] - 1] = -1. / h;
		B[m[k] - 1] = 1. + 1. / h;
		C[m[k] - 1] = 0;
		
		for (int j = 1; j < n[k]; j++)
		{	
			double t = j * tau;
			F[0] = 2 * pi * (1 - exp(-t));
			for (int i = 1; i < m[k] - 1; i++)
			{
				double x = i * h;
				F[i] = u[j - 1][i] / tau + sin(2*pi*x) - 1. / 2. / pi / pi;
			}
			F[m[k] - 1] = 3 + 2 * pi * (1 - exp(-t));

			//vector<double> U_next;
			tridiagonal_matrix(m[k], A, B, C, F, U_next);

			for (int i = 0; i < m[k]; i++)
			{
				u[j][i] = U_next[i];
			}
		}

		for (int j = 0; j < n[k]; j++)
		{
			for (int i = 0; i < m[k]; i++)
			{
				//u_1 << j * tau << " " << i * h << " " << u[j][i] << endl;
				diff = abs(u[j][i] - sol(j * tau, i * h));
				if (diff > M) {
					M = diff;
					//cout << j * tau << " " << i * h << " " << M << endl;
				}
			}
		}
		cout << "Step h(x): " << h << " Step T(t): " << tau << " Error: " << M << endl;
	}
	u_1.close();
	return 0;
}
