#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

constexpr auto pi = 3.14159265358979323846;

using namespace std;

double k(double u)
{
	return u;
	//return 1 / pi / pi / 4;
}

double F(double u)
{
	return u * u * u;
	//return 1;
}

double f(double x)
{
	return sin(2 * pi * x) - 1. / 2. / pi / pi;
}

void tridiagonal_matrix(const int n, const vector<double>& A, const vector<double>& B, const vector<double>& C, const vector<double>& D, vector<double>& X) {
	vector<double> P(n - 1);
	vector<double> Q(n);

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
}

double sol(double t, double x)
{
	return (1 - exp(-t)) * sin(2 * pi * x) + x * x;
}

int main()
{
	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_2\\out\\u3_out.txt");

	vector <int> m = { 10 };
	vector <int> n = { 540 };


	for (int z = 0; z < m.size(); z++)
	{
		double M = 0, diff = 0;
		vector < vector <double> > u(n[z], vector <double>(m[z]));

		double tau = 1. / (n[z]);
		double h = 1. / (m[z]);

		for (int i = 0; i < m[z]; i++)
		{
			u[0][i] = (i * i * h * h);
		}

		vector<double> A(m[z]);
		vector<double> B(m[z]);
		vector<double> C(m[z]);
		vector<double> D(m[z]);


		A[0] = 0;
		B[0] = -1. / h;
		C[0] = 1. / h;
		A[m[z] - 1] = -1. / h;
		B[m[z] - 1] = 1. + 1. / h;
		C[m[z] - 1] = 0;

		for (int j = 1; j < n[z]; j++)
		{
			double t = j * tau;
			D[0] = 2 * pi * (1 - exp(-t));
			for (int i = 1; i < m[z] - 1; i++)
			{
				double x = i * h;
				A[i] = -k(u[j - 1][i] / 2. + u[j - 1][i - 1] / 2.) / h / h;
				B[i] = 1. / tau + (k(u[j - 1][i] / 2. + u[j - 1][i - 1] / 2.) + k(u[j - 1][i + 1] / 2. + u[j - 1][i] / 2.)) / h / h;
				C[i] = -k(u[j - 1][i + 1] / 2. + u[j - 1][i] / 2.)  / h / h;
				D[i] = u[j - 1][i] / tau + f(x) * F(u[j - 1][i]);
			}
			D[m[z] - 1] = 3 + 2 * pi * (1 - exp(-t));

			tridiagonal_matrix(m[z], A, B, C, D, u[j]);
		}

		for (int j = 0; j < n[z]; j++)
		{
			for (int i = 0; i < m[z]; i++)
			{
				u_1 << j * tau << " " << i * h << " " << u[j][i] << endl;
				/*diff = abs(u[j][i] - sol(j * tau, i * h));
				if (diff > M) {
					M = diff;
					cout << j * tau << " " << i * h << " " << M << endl;
				}*/
			}
		}
		//cout << "Step h(x): " << h << " Step T(t): " << tau << " Error: " << M << endl;
	}
	u_1.close();
	return 0;
}
