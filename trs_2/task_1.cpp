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

	vector <int> m = { 10, 100, 1000, 10000 };
	vector <int> n = { 540, 5400, 54000, 54000 };

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

		for (int j = 1; j < n[k]; j++)
		{
			for (int i = 1; i < m[k] - 1; i++)
			{
				u[j][i] = tau / (4 * pi * pi * h * h) * (u[j - 1][i - 1] - 2 * u[j - 1][i] + u[j - 1][i + 1])
					+ tau * sin(2 * pi * i * h) - tau / (2 * pi * pi) + u[j - 1][i];
			}
			u[j][0] = - 2 * pi * h * (1 - exp(-tau * j)) + u[j][1];
			u[j][m[k] - 1] = (h * (3 + 2 * pi * (1 - exp(-j * tau))) + u[j][m[k] - 2]) / (1 + h);
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
