#include<iostream>
#include<vector>

using namespace std;

double func_1(double x1, double x2) {
	return ((2 * pow(x1, 3)) - pow(x2, 2) - 1);
}

double func_2(double x1, double x2) {
	return ((x1 * pow(x2, 3)) - x2 - 4);
}

void solution(vector<vector<double>> matrix, vector<double>& xx, vector<double> F)
{
	int N = matrix.size();
	for (int i = 0; i < N; i++)
	{
		double tmp = matrix[i][i];
		for (int j = N; j >= i; j--) {
			if (j == N) {
				F[i] /= tmp;
			}
			else {
				matrix[i][j] /= tmp;
			}
		}
		for (int j = i + 1; j < N; j++)
		{
			tmp = matrix[j][i];
			for (int k = N; k >= i; k--) {
				if (k == N) {
					F[j] /= tmp;
				}
				else {
					matrix[j][k] -= tmp * matrix[i][k];
				}
			}
		}
	}

	xx[N - 1] = -F[N - 1];
	for (int i = N - 2; i >= 0; i--)
	{
		xx[i] = -F[i];
		for (int j = i + 1; j < N; j++)
			xx[i] -= matrix[i][j] * xx[j];
	}

}

void makeVectorOfDiscrepancy(vector<double(*)(double, double)> M, vector<double> xx, vector<double>& vect)
{
	for (int i = 0; i < M.size(); i++)
	{
		// F(i) = 0
		double temp = 0;
		temp += M[i](xx[0],xx[1]);
		vect[i] = temp;
	}

}
void makeVectorOfJacoby(vector<double(*)(double, double)> M, vector<double> xx, vector<vector<double>>& vect) {
	for (int i = 0; i < vect.size(); i++) {
		for (int j = 0; j < vect[i].size(); j++) {
			vect[i][j] = 0;
			double temp = 0.01 * xx[j];
			xx[j] += temp;
			vect[i][j] += M[i](xx[0], xx[1]);
			xx[j] -= temp;
			vect[i][j] -= M[i](xx[0], xx[1]);
			vect[i][j] /= temp;
		}
	}
}
double sigma_1(vector<double> F) {
	double max = abs(F[0]);
	for (int i = 1; i < F.size(); i++) {
		if (abs(F[i]) > max) {
			max = abs(F[i]);
		}
	}
	return max;
}
double sigma_2(vector<double> temp, vector<double> x) {
	double max;
	if (abs(temp[0]) < 1) {
		max = abs(temp[0] - x[0]);
	}
	else {
		max = abs((temp[0] - x[0]) / temp[0]);
	}
	for (int i = 1; i < temp.size(); i++) {
		double t;
		if (abs(temp[i]) < 1) {
			t = temp[i] - x[i];
		}
		else {
			t = (temp[i] - x[i]) / temp[i];
		}
		if (abs(t) > max) {
			max = t;
		}
	}
	return max;
}

int main() {
	int n, NIT = 500;
	double E1 = pow(10, -9), E2 = E1;
	cout << "Enter amount of var:\n";
	cin >> n;
	vector<double> x(n);
	vector<double> dx(n);
	cout << "Enter first approximation\n";
	for (int i = 0; i < n; i++) {
		cin >> x[i];
	}
	vector<double(*)(double, double)> Functions(n);
	Functions[0] = func_1;
	Functions[1] = func_2;
	int iter = 0;
	vector<double> vect_F(n);
	vector<vector<double>> vect_J(n);
	for (int i = 0; i < n; i++) {
		vect_J[i].resize(n);
	}
	while (true) {
		double s1;
		double s2;
		vector<double> temp(n);
		makeVectorOfDiscrepancy(Functions, x, vect_F);
		makeVectorOfJacoby(Functions, x, vect_J);
		solution(vect_J, dx, vect_F);
		for (int i = 0; i < n; i++) {
			temp[i] = x[i] + dx[i];
		}
		iter++;
		s1 = sigma_1(vect_F);
		s2 = sigma_2(temp, x);
		cout << "iteration " << iter << ", sigma 1 = " << s1 << ", sigma 2 = " << s2 << endl;
		if (s1 <= E1 && s2 <= E2) {
			for (int i = 0; i < x.size(); i++) {
				cout << "x" << i + 1 << " = " << x[i] << endl;
			}
			break;
		}
		else if (iter == NIT) {
			for (int i = 0; i < x.size(); i++) {
				cout << "x" << i + 1 << " = " << x[i] << endl;
			}
			cout << "IER = 2";
			return 2;
		}
		for (int i = 0; i < n; i++) {
			x[i] = temp[i];
		}
	}
	return 0;
}
