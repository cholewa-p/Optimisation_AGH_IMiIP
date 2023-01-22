/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
*********************************************/

#include"opt_alg.h"
#include <fstream>

void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		//lab1();
		//lab2();
		//lab3();
		//lab4();

		//lab5();

		lab6();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}



	system("pause");
	return 0;
}

void lab1()
{
	//cout << fib(fun, -10, 10, 1e-5) << endl;
	//cout << fib(fun, 10, 100, 1e-5) << endl; 
	//cout << lag(fun, -100, 100, 0.001, 1e-7, 1000) << endl;
	//cout << lag(fun, -10, 1, 0.0001, 1e-7, 1000) << endl;
	//double* p = expansion(fun, 10, 1.0, 1.5, 1000);	 
	//cout << p[0] << " " << p[1] << endl;
	double epsilon = 1e-10, gamma = 1e-200, Nmax = 1000;
	cout << fib(fR, 1e-4, 1e-2, epsilon) << endl;
	//cout << lag(fR, 1e-4, 1e-2, epsilon, gamma, Nmax) << endl;

	ofstream out("proba1.csv");

	srand(time(NULL));
	double alpha = rand() % 100 + 101;
	double d = 1;
	double x0 = rand() % 200 - 100;
	alpha /= 100;
	double* p = new double[2];

	/*
	for (int i = 0; i < 100; i++) {
		alpha = 2.5;
		x0 = rand() % 200 - 100;
		//alpha /= 100;

		p = expansion(fun, x0, d, alpha, 10000);
		out << x0 << ";" << p[0] << ";" << p[1] << ";" << solution::f_calls <<";";
		solution::f_calls = 0;

		solution wynik_fib = fib(fun, p[0], p[1], 1e-5);
		out << wynik_fib.x(0) << ";" << wynik_fib.y(0) << ";" << solution::f_calls << ";";
		string lok = "";
		solution::f_calls = 0;

		if (wynik_fib.y(0) < -0.75)
			lok = "globalne";
		else
			lok = "lokalne";
		out << lok << ";";

		solution wynik_lag = lag(fun, p[0], p[1], 0.0001, 1e-7, 1000);
		out << wynik_lag.x(0) << ";" << wynik_lag.y(0) << ";"<< solution::f_calls << ";";
		solution::f_calls = 0;

		if (wynik_lag.y(0) < -0.75)
			lok = "globalne";
		else
			lok = "lokalne";

		out << lok << ";\n";
		cout << i << endl;
	}
	*/
	out.close();

	//cout << fib(fun, -100, 100, 1e-5) << endl;
	//cout << lag(fun, -100, 100, 1e-5, 1e-7, 1000) << endl;
}

void lab2()
{
	double s = 1.5;
	double alphaHJ = 0.5;
	double alphaR = 2;
	double beta = 0.5;
	double epsilon = 1e-3;
	int Nmax = 1000;
	matrix x0 = matrix(2, 1, 1.0);
	matrix x1 = matrix(2, 1, 1);

	//x0(0) = -0.86;
	//x0(1) = -0.75;

	//cout << HJ(fun2, x0, s, alphaHJ, epsilon, Nmax) << endl;
	//cout << Rosen(fun2, x0, s, alphaR, beta, epsilon, Nmax) << endl;
	solution::clear_calls();
	//cout << HJ(fR2, x1, s, alphaHJ, epsilon, Nmax) << endl;
	cout << Rosen(fR2, x1, matrix(2, 1, s), alphaR, beta, epsilon, Nmax) << endl;

	//x1(0) = 2.76953; // dla mc = 10
	//x1(1) = 3.18359;

	x1(0) = 2.75195; // dla mc = 9
	x1(1) = 2.99219;
	/*
	//cout << fR2(x1, NULL, NULL) << endl;
	ifstream inx("x15.txt");
	ifstream iny("y15.txt");

	for (int i = 0; i < 100; i++) {

		inx >> x0(0);
		iny >> x0(1);
		//cout << x0(0) << ";" << x0(1) << ";";
		solution wynik_R = Rosen(fun2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax);
		cout << wynik_R.x(0) << ";" << wynik_R.x(1) << ";" << wynik_R.y << solution::f_calls << endl;
		
		solution::clear_calls();
	}*/
	
	x0(0) = 0.88;
	x0(1) = 0.84;
	//cout << HJ(fun2, x0, s, alphaHJ, epsilon, Nmax) << endl;
	//matrix s(2, 1, 0.5);
	//cout << Rosen(fun2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax) << endl;
}

void lab3()
{
	matrix x0 = matrix(2, 1, 1.0);

	double c0 = 1;
	matrix a[3] = { 4, 4.4934, 5 };

	const double epsilon = 1e-3;
	const int Nmax = 10000;

	const int ktory = 2;

	// Generacja danych do pierwszej czesci zadania
	double pkt1[100] = { 0 };
	double pkt2[100] = { 0 };

	//cout.precision(20);

	for (int i = 0; i < 100; i++) {
		do
			x0 = 5 * rand_mat(2, 1) + 1; // generowanie punktow znajdujacych sie w obszarze funkcji
		while (norm(x0) > a[ktory]);

		pkt1[i] = x0(0);
		pkt2[i] = x0(1);
	}
	// Wygenerowane wartosci sa dluzsze niz wyswietla
	// Prowadzi to do niedokladnosci - wyswietlany wynik moze sie roznic od tego dla dokladnie takiej wartosci

	/**/
	for (int i = 0; i < 100; i++) {
		x0(0) = pkt1[i];
		x0(1) = pkt2[i];

		cout << x0(0) << ";" << x0(1) << ";";

		solution zew_opt = pen(fun3, x0, c0, 2, epsilon, Nmax, a[ktory]);
		cout << zew_opt.x(0) << ";" << zew_opt.x(1) << ";" << norm(zew_opt.x) << ";" << zew_opt.y[0] << solution::f_calls << ";";
		solution::clear_calls();

		solution wew_opt = pen(fun3, x0, c0, 0.5, epsilon, Nmax, a[ktory]);
		cout << wew_opt.x(0) << ";" << wew_opt.x(1) << ";" << norm(wew_opt.x) << ";" << wew_opt.y[0] << solution::f_calls << endl;
		solution::clear_calls();
	}
	
	
	// Problem rzeczywisty
	/*
	matrix x1 = matrix(2, 1);
	x1(0) = 0.;
	x1(1) = 0.;
	
	//cout << pen(fR3, x1, c0, 2, epsilon, Nmax) << endl;

	//-3.3473 25.0006

	x1(0) = -3.3473;
	x1(1) = 25;

	fR3(x1, c0, 2);
	*/
}

void lab4()
{
	matrix x0 = matrix(2, 1, 0.0);
	const double epsilon = 1e-3;
	const int Nmax = 10000;

	double h = 0.05;

	double pkt1[100] = { 0 };
	double pkt2[100] = { 0 };

	/*
	cout << SD(fT4, grad4, x0, h, epsilon, Nmax) << endl;
	cout << CG(fT4, grad4, x0, h, epsilon, Nmax) << endl;
	cout << Newton(fT4, grad4, hesj4, x0, h, epsilon, Nmax) << endl;
	*/

	for (int i = 0; i < 100; i++) {

		x0 = 20 * rand_mat(2, 1) -10; // generowanie punktow znajdujacych sie w obszarze funkcji
		pkt1[i] = x0(0);
		pkt2[i] = x0(1);
	}
	/*
	for(int i = 0; i < 100; i++) {
		x0(0) = pkt1[i];
		x0(1) = pkt2[i];

		cout << x0(0) << ";" << x0(1) << ";";

		solution SD_sol = SD(fT4, grad4, x0, h, epsilon, Nmax);
		cout << SD_sol.x(0) << ";" << SD_sol.x(1) << ";" << SD_sol.y[0] << SD_sol.f_calls << ";" << SD_sol.g_calls << ";";
		solution::clear_calls();
		solution CG_sol = CG(fT4, grad4, x0, h, epsilon, Nmax);
		cout << CG_sol.x(0) << ";" << CG_sol.x(1) << ";" << CG_sol.y[0] << CG_sol.f_calls << ";" << CG_sol.g_calls << ";";
		solution::clear_calls();
		solution Nw_sol = Newton(fT4, grad4, hesj4, x0, h, epsilon, Nmax);
		cout << Nw_sol.x(0) << ";" << Nw_sol.x(1) << ";" << Nw_sol.y[0] << Nw_sol.f_calls << ";" << Nw_sol.g_calls << ";" << Nw_sol.H_calls << ";";
		solution::clear_calls();
		
		cout << endl;
	}
	*/
	/*
	x0(0) = 6.2231;
	x0(1) = -0.238002;

	Newton(fT4, grad4, hesj4, x0, -1, epsilon, Nmax);	
	*/
	/**/
	matrix x1(3, 1, 0.0);

	solution Real_sol = CG(fR4, gf, x1, 0.001, 0.000001, Nmax);

	int m = 100;
	static matrix X(3, m), Y(1, m);
	ifstream in("XData.txt");
	in >> X;
	in.close();
	in.open("YData.txt");
	in >> Y;
	in.close();

	double P = 0.0;

	for (int i = 0; i < 100; i++) {
		double h = 1.0 / (1 + exp(-(trans(Real_sol.x) * X[i])()));
		if (lroundf(h) == Y(0, i))
			h = 1;
		else 
			h = 0;

		P += h;
	}
	P /= m;

	cout << Real_sol.x(0,0) << " " << Real_sol.x(1, 0) << " " << Real_sol.x(2, 0) << " " << Real_sol.y(0, 0) << " " << P << " " << Real_sol.g_calls << endl;
	
}

void lab5()
{
	/*
	matrix x0 = 20 * rand_mat(2, 1) - 10, * ud = new matrix[2];
	double epsilon = 1e-3, w = 0;
	int Nmax = 5000;
	double a = 1;
	ud[0] = a;
	ud[1] = w;

	solution opt = Powell(fT5, x0, epsilon, Nmax, *ud, NAN);
	*/
	/**/

	// Problem Testowy
	/**/
	solution solPow;
	double a = 10;
	double w = 0.0;
	matrix ud1(2, new double[2] {w, a});

	const double epsilon = 1e-3;
	const int Nmax = 10000;

	for (int i = 0; i < 101; i++) {
		matrix x0 = 20 * rand_mat(2, 1) - 10;

		ud1(1) = 1.;
		solPow = Powell(f5, x0, epsilon, Nmax, ud1, NAN);
		cout <<  x0(0) << ";" << x0(1) << ";" << solPow.x(0) << ";" << solPow.x(1) << ";" << solPow.y(0) << ";" << solPow.y(1) << ";" << solPow.f_calls << ";";
		solution::clear_calls();

		ud1(1) = 10.;
		solPow = Powell(f5, x0, epsilon, Nmax, ud1, NAN);
		cout << solPow.x(0) << ";" << solPow.x(1) << ";" << solPow.y(0) << ";" << solPow.y(1) << ";" << solPow.f_calls << ";";
		solution::clear_calls();

		ud1(1) = 100.;
		solPow = Powell(f5, x0, epsilon, Nmax, ud1, NAN);
		cout << solPow.x(0) << ";" << solPow.x(1) << ";" << solPow.y(0) << ";" << solPow.y(1) << ";" << solPow.f_calls << endl;
		solution::clear_calls();

		ud1(0) += 0.01;
	}
	

	// Problem rzeczywisty
	/*
	solution naszSolution;
	double waga = 0.0;
	matrix ud1(waga);
	for (int i = 0; i < 101; i++) {
		double l, d;
		l = (900. * ((double)rand() / (double)RAND_MAX)) + 100.;
		d = (40. * ((double)rand() / (double)RAND_MAX)) - 10.;
		double pom[2] = { l,d };
		matrix x0(2, pom);
		naszSolution = Powell(fR5, x0, 0.001, 1000, ud1, 0);
		cout << l << " " << d << " " << naszSolution.x(0) << " " << naszSolution.x(1) << " " << naszSolution.y(0) << " " << naszSolution.y(1) << " " << " " << solution::f_calls << endl;
		solution::clear_calls();

		ud1(0) += 0.01;
	}
	*/
	
}

void lab6()
{
	int N = 2, Nmax = 10000, mi = 20, lambda = 40;
	double epsilon = 1e-3;

	// Zadanie Testowe
	/*
	matrix limits(2, 2), sigma0(2, 1);
	limits(0, 0) = limits(1, 0) = -5;
	limits(0, 1) = limits(1, 1) = 5;
	sigma0(0) = sigma0(1) = 1;

	solution solEvo = EA(fT6, N, limits, mi, lambda, sigma0, epsilon, Nmax);

	string s = "NIE";
	sigma0(0) = sigma0(1) = 100; // sigma - 0.01 / 0.1 / 1 / 10 / 100
	for (int i = 0; i < 100; ++i) {
		solution::clear_calls();
		solEvo = EA(fT6, N, limits, mi, lambda, sigma0, epsilon, Nmax);

		if (solEvo.y < 0.001) 
			s = "TAK";

		cout << solEvo.x(0) << ";" << solEvo.x(1) << ";" << solEvo.y << solEvo.f_calls << ";"  << s << endl;
		s = "NIE";
	}
	*/

	// Problem rzeczywisty
	matrix limits(2, 2), sigma1(2, 1);
	limits(0, 0) = limits(1, 0) = 0.1;
	limits(0, 1) = limits(1, 1) = 3;

	sigma1(0) = sigma1(1) = 1;

	matrix simulation = matrix(0);
	solution test;

	// rzeczywisty przeprowadzono dla epsilon = 1e-4
	//test = EA(fR6, N, limits, mi, lambda, sigma1, epsilon, Nmax, simulation);

	//cout << test << endl;
	matrix xopt(2, 1);
	xopt(0) = 2.10424;
	xopt(1) = 0.0142826;

	fR6(xopt, NAN, NAN);


}