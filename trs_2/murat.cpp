#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include<iterator>

constexpr auto pi = 3.14159265358979323846;

using namespace std;

double k(double u)
{
	return sin(u);
	//return 1 / pi / pi / 4;
}

double F(double u)
{
	return u;
	//return 1;
}

double f(double x, double t)
{
	return t * x * x - t * t;
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
	ofstream u_4("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_2\\out\\murat.csv");

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
			u[0][i] = (i * h);
			u_4 << 0 << " " << i * h << " " << u[0][i] << endl;
			cout << 0 << " " << i * h << " " << u[0][i] << endl;
		}

		vector<double> A(m[z]);
		vector<double> B(m[z]);
		vector<double> C(m[z]);
		vector<double> D(m[z]);
		vector<double> prev;

		A[0] = 0;
		B[0] = -1. / h;
		C[0] = 1. / h;
		A[m[z] - 1] = -1. / h;
		B[m[z] - 1] = 1. + 1. / h;
		C[m[z] - 1] = 0;

		for (int j = 1; j < n[z]; j++)
		{
			double t = j * tau;
			D[0] = 1;
			for (int i = 1; i < m[z] - 1; i++)
			{
				double x = i * h;
				A[i] = -k(u[j - 1][i] / 2. + u[j - 1][i - 1] / 2.) / h / h;
				B[i] = 1. / tau + (k(u[j - 1][i] / 2. + u[j - 1][i - 1] / 2.) + k(u[j - 1][i + 1] / 2. + u[j - 1][i] / 2.)) / h / h;
				C[i] = -k(u[j - 1][i + 1] / 2. + u[j - 1][i] / 2.) / h / h;
				D[i] = u[j - 1][i] / tau + f(x,t) * F(u[j - 1][i]);
			}
			D[m[z] - 1] = 3. / 2. * (t * t) + 2;

			tridiagonal_matrix(m[z], A, B, C, D, u[j]);

			//auto max_elemenT = std::max_element(u[j].begin(), u[j].end());
			//max_value = k(*max_elemenT);

			int Q = 0;
			do
			{
				M = 0;
				Q++;
				prev.assign(u[j].begin(), u[j].end());
				for (int i = 1; i < m[z] - 1; i++)
				{
					double x = i * h;
					A[i] = -k(u[j][i] / 2. + u[j][i - 1] / 2.) / h / h;
					B[i] = 1. / tau + (k(u[j][i] / 2. + u[j][i - 1] / 2.) + k(u[j][i + 1] / 2. + u[j][i] / 2.)) / h / h;
					C[i] = -k(u[j][i + 1] / 2. + u[j][i] / 2.) / h / h;
					D[i] = u[j][i] / tau + f(x, t) * F(u[j][i]);
				}
				tridiagonal_matrix(m[z], A, B, C, D, u[j]);

				for (int s = 0; s < m[z]; s++)
				{
					diff = abs(k(u[j][s]) - k(prev[s]));
					if (diff > M) {
						M = diff;
					}
				}
				prev.clear();
				//cout << Q << " ";
			} while (M > 1e-6);

			cout << endl;
			for (int i = 0; i < m[z]; i++)
			{
				u_4 << j * tau << " " << i * h << " " << u[j][i] << endl;
				cout << j * tau << " " << i * h << " " << u[j][i] << endl;
			}
		}

		//for (int j = 0; j < n[z]; j++)
		//{
		//	for (int i = 0; i < m[z]; i++)
		//	{
		//		u_4 << j * tau << " " << i * h << " " << u[j][i] << endl << endl;
		//		/*diff = abs(u[j][i] - sol(j * tau, i * h));
		//		if (diff > M) {
		//			M = diff;
		//			cout << j * tau << " " << i * h << " " << M << endl;
		//		}*/
		//	}
		//	u_4 << endl << endl;
		//}
		////cout << "Step h(x): " << h << " Step T(t): " << tau << " Error: " << M << endl;
	}
	u_4.close();
	return 0;
}
