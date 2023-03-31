#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

constexpr auto pi = 3.14159265358979323846;

using namespace std;

double sol(double t, double x)
{
	return (1 - exp(-t)) * sin(2 * pi * x) + x * x;
}

int main()
{
	ofstream u_1("C:\\Users\\PETA4\\Desktop\\3\\6_sem\\trs\\labs\\trs_2\\out\\u1_out.csv");

	const int m = 1000;
	const int n = 54000;
	double M = 0, diff;

	vector < vector <double> > u(n, vector <double>(m));

	double tau = 1. / (n - 1);
	double h = 1. / (m - 1);

	for (int i = 0; i < m; i++)
	{
		u[0][i] = (i * i * h * h);
	}

	for (int j = 1; j < n; j++)
	{
		for (int i = 1; i < m - 1; i++)
		{
			u[j][i] = tau / (4 * pi * pi * h * h) * (u[j - 1][i - 1] - 2 * u[j - 1][i] + u[j - 1][i + 1])
				+ tau * sin(2 * pi * i * h) - tau / (2 * pi * pi) + u[j - 1][i];
		}
		u[j][0] = - 2 * pi * h * (1 - exp(-tau * j)) + u[j][1];
		u[j][m - 1] = (h * (3 + 2 * pi * (1 - exp(-j * tau))) + u[j][m - 2]) / (1 + h);
	}

	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			//u_1 << j * tau << " " << i * h << " " << u[j][i] << endl;
			diff = abs(u[j][i] - sol(j * tau, i * h));
			if (diff > M) {
				M = diff;
				cout << j * tau << " " << i * h << " " << M << endl;
			}
		}
	}

	cout << M;
	u_1.close();
	return 0;
}
