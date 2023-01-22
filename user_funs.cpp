#include"user_funs.h"
#include<vector>
#include<math.h>
#include"matrix.h"
#include"opt_alg.h"

# define M_PI 3.14159265358979323846

double f(double x) {
	return -cos(0.1 * x) * exp(-pow(0.1 * x - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x, 2);
}

matrix fun(matrix x, matrix ud1, matrix ud2) {
	return -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
}

matrix fun2(matrix x0, matrix ud1, matrix ud2) {
	return pow(x0[0](0), 2) + pow(x0[0](1),2) - cos(2.5 * 3.14 * x0[0](0)) - cos(2.5 * 3.14 * x0[0](1)) + 2;
}


double MetodaFibonacciego(double ap, double bp, double epsilonp) {
	vector<double> a;
	vector<double> b;
	vector<double> c;
	vector<double> d;
	int k = 0;

	while (true) {
		if (FiN(k) > (bp - ap) / epsilonp) {		
			break;
		}
		else {
			k++;
		}
	}

	a.push_back(ap); // 50
	b.push_back(bp); // 70
	c.push_back(b[0] - FiN(k-1) / FiN(k) * (b[0] - a[0])); // 70 - 1597 / 2584 * (70 - 50) = 57,639
	d.push_back(a[0] + b[0] - c[0]); // 50 + 70 - 57,639 = 62,361
	int fcalls = 4;

	for (int i = 0; i < k - 3; i++) {
		if ( f(c[i]) < f(d[i]) ) {
			a.push_back(a[i]);
			b.push_back(d[i]);
		}
		else {
			b.push_back(b[i]); // b[1] = b[0]
			a.push_back(c[i]); // a[1] = c[0]
		}
		c.push_back(b[i + 1] - 1.0 * FiN(k - i - 2) / FiN(k - i - 1) * (b[i + 1] - a[i + 1])); // 70 - 987 / 1597 * (70 - 57,639) // 62,36048
		d.push_back(a[i + 1] + b[i + 1] - c[i + 1]); // 57,639 + 70 - 62,36048 = 65,27852
		fcalls += 2;
	}

	int pom = c.size();
	cout << c[pom - 1] << endl;
	cout << fcalls << endl;

	return c[pom-1];
}

double FiN(int n) {
	return (1 / sqrt(5) * pow( ( 1+sqrt(5) ) / 2, n ) - 1 / sqrt(5) * pow((1 - sqrt(5)) / 2, n));
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	double a = 0.98, b = 0.63, g = 9.81;
	double Pa = 0.75, Va0 = 5, Ta0 = 90;
	double Pb = 1, Vb0 = 1, Tb0 = 10;
	double Tbin = 10, Fbin = 0.01;
	double Db = 0.00365665;
	
	matrix dY(3, 1);
	double Faout = a * b * m2d(ud2) * sqrt(2 * g * Y(0) / Pa);
	double Fbout = a * b * Db * sqrt(2 * g * Y(1) / Pb);

	dY(0) = -Faout;
	dY(1) = Faout + Fbin - Fbout;
	dY(2) = Fbin / Y(1) * (Tbin - Y(2)) + Faout / Y(1) * (Ta0 - Y(2));
	return dY;
}

matrix fR(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix YO = matrix(3, new double[3]{ 5, 1, 10 });

	matrix* Y = solve_ode(df1, 0, 1, 1000, YO, ud1, x);

	// Y[0] -> czas, Y[1] -> Va, Vb, Tb

	double maxT = Y[1](0, 2);
	int n = get_len(Y[0]);

	for (int i = 1; i < n; i++)
		if (maxT < Y[1](i, 2)) {
			maxT = Y[1](i, 2);
		}

	/*
	for (int i = 0; i < n; i++) {
		cout << Y[0](i) << ": " << Y[1](i, 0) << " " << Y[1](i, 1) << " "<<  Y[1](i, 2) << endl;
	}
	*/
	y = matrix(abs(m2d(maxT) - 50)); // Tu zmieniac temperature w przykladach, teraz 50

	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	double mr = 1;
	double l = 0.5;
	double mc = 9; // wyniki podane przez prowadzaca sa dla 10 xd
	double b = 0.5;

	double I = mr * l * l / 3 + mc * l * l;
	matrix dY(2, 1);

	dY(0) = Y(1);
	dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;

	return dY;
}

matrix fR2(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = 0;
	matrix Y0(2, 1);
	matrix Yref(2, new double[2]{ 3.14, 0 });
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, x);
	int n = get_len(Y[0]);
	for (int i = 0; i < n; i++) {
		y = y + 10 * pow(Yref(0) - Y[1](i, 0), 2) + pow(Yref(1) - Y[1](i, 1), 2) 
			+ pow(x(0) * (Yref(0) - Y[1](i, 0)) + x(1) * (Yref(1) - Y[1](i, 1)), 2);
		cout << (double)i/10  << ";" << Y[1](i, 0) << ";" << Y[1](i, 1) << endl;
	}
	y = y * 0.1;
	//cout << x(0) << "  " << x(1) << endl;
	return y;
}

matrix Fun3(matrix x1, matrix x2, matrix ud1) {
	return (sin(M_PI * sqrt(pow(x1() / M_PI, 2) + pow(x2() / M_PI, 2)))) / (M_PI * sqrt(pow(x1() / M_PI, 2) + pow(x2() / M_PI, 2)));
}

bool g1(matrix x1) {
	if (-x1() + 1 <= 0)
		return true;
	else
		return false;
}
bool g2(matrix x2) {
	if (-x2() + 1 <= 0)
		return true;
	else
		return false;
}
bool g3(matrix x1, matrix x2, double alpha) {
	if (sqrt(pow(x1(), 2)+pow(x2(),2)) - alpha <= 0)
		return true;
	else
		return false;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
	double C = 0.47;
	double r = 0.12; // 12 cm
	double m = 0.6; // 600 g
	double ro = 1.2; // 1.2 kg/m3
	double g = 9.81; // 9.81 m/s2

	double S = M_PI * pow(r, 2);

	double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));

	double Fmx = M_PI * ro * Y(3) * m2d(ud2) * pow(r, 3);
	double Fmy = M_PI * ro * Y(1) * m2d(ud2) * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = Y(1);
	dY(1) = (-Dx - Fmx) / m;
	dY(2) = Y(3);
	dY(3) = (-m * g - Dy - Fmy) / m;

	return dY;
}

matrix fR3(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0(4, new double[4]{ 0, x(0), 100, 0 });
	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, ud1, x(1)); // to ud1?
	int n = get_len(Y[0]);
	int i0 = 0, i50 = 0;

	//cout << x(0) << " " << x(1) << endl;

	for (int i = 0; i < n; i++) {
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
			i0 = i;

		cout << Y[1](i, 0) << ";";
		//cout << Y[1](i, 1) << ";";
		cout << Y[1](i, 2) << ";";
		//cout << Y[1](i, 3) << ";";
		cout << endl;
	}

	y = -Y[1](i0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 25 > 0)
		y = y + ud2 * pow(abs(x(1)) - 25, 2);
	if (abs(Y[1](i50, 0) - 5) - 1 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5) - 1, 2);

	//cout << Y0(0) << " " << Y0(1) << " " << Y0(2) << " " << Y0(3) << endl;
	//cout << x(0) << " " << x(1) << endl; // z optymalnego dla tego uruchomic symulacje

	//cout << "y = " << y << endl;

	return y;
}

matrix fun3(matrix x, matrix ud1, matrix ud2) {
	double arg = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
	matrix y = sin(arg) / arg;
	// y = pow(x(0),2) + pow(x(1),2);

	if (ud2(1) > 1) { //kara zew
		if (-x(0) + 1 > 0)		// g1
			y = y + (ud2)(0) * pow(-x(0) + 1, 2);
		if (-x(1) + 1 > 0)		// g2
			y = y + (ud2)(0) * pow(-x(1) + 1, 2);
		if (norm(x) - (ud1)(0) > 0)  // g3
			y = y + (ud2)(0) * pow(norm(x) - (ud1)(0), 2);
	}
	else { //kara wew
		if (-x(0) + 1 > 0)
			y = 1e10;
		else
			y = y - (ud2)(0) / (-x(0) + 1);

		if (-x(1) + 1 > 0)
			y = 1e10;
		else
			y = y - (ud2)(0) / (-x(1) + 1);

		if (norm(x) - (ud1)(0) > 0)
			y = 1e10;
		else
			y = y - (ud2)(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (ud1)(0));
	}
	return y;
}

matrix fun4(matrix x, matrix ud1, matrix ud2) {
	return pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
}

matrix grad4(matrix x, matrix ud1, matrix ud2) {
	matrix g(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	return g;
}

matrix hesj4(matrix x, matrix ud1, matrix ud2) {
	matrix H(2, 2);
	H(0, 0) = 10;
	H(0, 1) = 8;
	H(1, 0) = 8;
	H(1, 1) = 10;
	return H;
}

matrix fT4(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	if (isnan(ud2(0, 0)))
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	else
		y = fT4(ud2[0] + x * ud2[1], 0, ud1);
	return y;
}

matrix fR4(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	int m = 100;
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if(solution::f_calls == 1) {
	ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}
	int P = 0;
	double h;
	y = 0;
	for (int i = 0; i < m; i++) {
		h = m2d(trans(x) * X[i]);
		h = 1.0 / (1.0 + exp(-h));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	}
	y = y / m;
	return y;
}

matrix gf(matrix x, matrix ud1, matrix ud2) {
	int m = 100;
	int n = get_len(x);
	matrix g(n, 1);
	static matrix X(n, m), Y(1, m);
	if (solution::g_calls == 1) {
		ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}

	double h;
	for (int j = 0; j < n; ++j) {
		for (int i = 0; i < m; ++i) {
			h = m2d(trans(x) * X[i]);
			h = 1 / (1 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
	}
	return g;
}

matrix fT5(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	/**/
	if (isnan(ud2(0, 0))) {
		y = matrix(2, 1);
		y(0) = ud1(1) * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2)); // ud1(1) to a
		y(1) = (1.0 / ud1(1)) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	} else {
		matrix yt;
		yt = fT5(ud2[0] + x * ud2[1], ud1, NAN); // ud2[0] to xi a ud2[1] to di x tutaj to h
		y = ud1(0) * yt(0) + (1 - ud1(0)) * yt(1); // ud1(0) to w
	}
	return y;
	
	 /*
	if (ud2 == NULL) {
		y = matrix(2, 1);
		y(0) = ud1[0]() * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2)); // ud1(1) to a
		y(1) = 1.0 / ud1[0]() * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	}
	else {
		solution T;
		T.x = ud2[0] + x * ud2[1];
		T.fit_fun(fT5, ud1);
		y = ud1[1]() * T.y(0) + (1 - ud1[1]()) * T.y(1);
	}
	return y;
	*/
}

matrix fR5(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	if (isnan(ud2(0, 0))) {
		y = matrix(3, 1);
		double ro = 7800, P = 1e3, E = 207e9;
		y(0) = ro * x(0) * 3.14 * pow(x(1), 2) / 4; // ro * l * pi * r^2 - masa
		y(1) = 64 * P * pow(x(0), 3) / (3 * E * 3.14 * pow(x(1), 4)); // ugiecie 64 * P * l^3 / 3 * E * pi * d^4
		y(2) = 32 * P * x(0) / (3.14 * pow(x(1), 3)); // naprezenie 31 * P * l / pi * d^3
	} else {
		matrix yt, xt = ud2[0] + x * ud2[1]; // ud2[0] = xi, ud2[1] = di
		yt = fR5(xt, ud1, NAN);
		y = ud1 * (yt(0) - 0.06) / (1.53 - 0.06) + (1 - ud1) * (yt(1) - 5.25e-6) / (0.0032 - 5.25e-6);
		double c = 1e10;
		if (xt(0) < 0.1) // jesli l mniejsze niz 100mm
			y = y + c * (pow(0.1 - xt(0), 2));
		if (xt(0) > 1) // jesli l > lm
			y = y + c * (pow(xt(0) - 1, 2));		
		if (xt(1) < 0.01) // jesli srednica <10mm
			y = y + c * (pow(0.01 - xt(1), 2));		
		if (xt(1) > 0.05) // jesli d>50mm
			y = y + c * (pow(xt(1) - 0.05, 2));		
		if (yt(1) > 0.005) // jesli ugiecie wieksze niz 5mm
			y = y + c * (pow(yt(1) - 0.005, 2));
		if (yt(2) > 300e6) // jesli naprezenie wieksze niz 300MPa
			y = y + c * (pow(yt(2) - 300e6, 2));
	}
	return y;
}

matrix f5_1(double a, matrix x, matrix ud1, matrix ud2) {
	return (a * (pow(x(0) - 2), 2) + (pow(x(1) - 2), 2));
}
matrix f5_2(double a, matrix x, matrix ud1, matrix ud2) {
	return ((1 / a) * (pow(x(0) + 2), 2) + (pow(x(1) + 2), 2));
}

matrix f5(matrix x, matrix ud1, matrix ud2) {
	if (isnan(ud2(0, 0))) {
		matrix y;
		y = matrix(2, new double[2]{ 0., 0. });
		y(0) = ud1(1) * (pow(x(0) - 2, 2) + (pow(x(1) - 2, 2)));
		y(1) = (1 / ud1(1)) * (pow(x(0) + 2, 2) + (pow(x(1) + 2, 2)));
		return y;
	}
	else {
		return ud1(0) * f5_1(ud1(1), ud2[0] + x * ud2[1], ud1, NAN) + (1 - ud1(0)) * f5_2(ud1(1), ud2[0] + x * ud2[1], ud1, NAN);
	}
}

matrix fT6(matrix x, matrix ud1, matrix ud2) {
	return pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
}

matrix df6(double t, matrix Y, matrix ud1, matrix ud2) {
	double b1 = ud2(0);
	double b2 = ud2(1);

	double k1 = 1;
	double k2 = 1;

	double m1 = 5;
	double m2 = 5;

	double F = 1;

	matrix dY(4, 1);
	dY(0) = Y(1);
	dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
	dY(2) = Y(3);
	dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;

	return dY;
}

matrix fR6(matrix x, matrix ud1, matrix ud2) {
	
	matrix y;
	int N = 1001;
	matrix X(N, 2);

	ifstream in("polozenia.txt");
	in >> X;
	in.close();

	matrix YO(4, new double[4]{ 0, 0, 0, 0 });
	matrix* Y = solve_ode(df6, 0, 0.1, 100, YO, ud1, x);

	y = 0;

	for (int i = 0; i < N; i++) {
		y = y + abs(X(i, 0) - Y[1](i, 0)) + abs(X(i, 1) - Y[1](i, 2));
		cout << Y[1](i, 0) << ";" << Y[1](i, 2) << endl;
	}

	y = y / (2 * N);

	return y;
}