#include"opt_alg.h"
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		solution XO(x0), XI(x0 + d);
		XO.fit_fun(ff, ud1, ud2);
		XI.fit_fun(ff, ud1, ud2);

		if (XO.y == XI.y)
		{
			p[0] = XO.x(0);
			p[1] = XI.x(0);
			return p;
		}
		if (XO.y < XI.y)
		{
			d *= -1;
			XI.x = XO.x + d;
			XI.fit_fun(ff, ud1, ud2);
			if (XI.y >= XO.y)
			{
				p[0] = XI.x(0);
				p[1] = -XI.x(0);
				return p;
			}
		}
		solution X2;
		int i = 1;
		while (true)
		{
			X2.x = x0 + pow(alpha, i) * d;
			X2.fit_fun(ff, ud1, ud2);
			if (X2.y >= XI.y || solution::f_calls > Nmax)
				break;
			XO = XI;
			XI = X2;
			++i;
		}
		if (d > 0)
		{
			p[0] = XO.x(0);
			p[1] = X2.x(0);
		}
		else
		{
			p[0] = X2.x(0);
			p[1] = XO.x(0);
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{ //(*ff)(x, ud1, ud2) ud1, ud2 - dodatkowe parametry, nie potrzebne teraz
	try
	{
		solution Xopt;
		Xopt.ud = b - a;
		int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
		int* F = new int[n] {1, 1};
		for (int i = 2; i < n; ++i)
			F[i] = F[i - 2] + F[i - 1];
		solution A(a), B(b), C, D;
		C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);
		cout << B.x - A.x << endl;
		
		for (int i = 0; i <= n - 3; ++i) 
		{
			if (C.y < D.y) 
			{
				A = A;
				B = D;
			}
			else 
			{
				B = B;
				A = C;
			}
			C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
			D.x = A.x + B.x - C.x;
			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);
			cout << B.x - A.x << endl;
		}

		Xopt = C;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.ud = b - a;
		solution A(a), B(b), C, D, D_old(a);
		C.x = (a + b) / 2;
		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);
		double l, m;
		cout << B.x - A.x << endl;
		while (true)
		{
			l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y * (pow(C.x) - pow(A.x)) + C.y * (pow(A.x) - pow(B.x)));
			m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));
			if (m <= 0)
			{
				Xopt = D_old;
				Xopt.flag = 2;
				return Xopt;
			}
			D.x = 0.5 * l / m;
			D.fit_fun(ff, ud1, ud2);
			if (A.x <= D.x && D.x <= C.x)
			{
				if (D.y <= C.y)
				{
					A = A;
					B = C;
					C = D;					
				}
				else
				{
					A = D;
					C = C;
					B = B;
				}
			}
			else if (C.x <= D.x && D.x <= B.x)
			{
				if (D.y <= C.y)
				{
					A = C;
					C = D;
					B = B;
				}
				else
				{
					A = A;
					C = C;
					B = D;
				}
			}
			else
			{
				Xopt = D_old;
				Xopt.flag = 2;
				return Xopt;
			}
			cout << B.x - A.x << endl;
			Xopt.ud.add_row((B.x - A.x)());
			if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma)
			{
				Xopt = D;
				Xopt.flag = 0;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = D;
				Xopt.flag = 1;
				break;
			}
			D_old = D;

		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution XB, XB_old, X;
		XB.x = x0;
		XB.fit_fun(ff);
		while (true)
		{
			//cout << XB.x(0) << endl;
			//cout << " ";
			cout << XB.x(1) << endl;
			X = HJ_trial(ff, XB, s);
			if (X.y < XB.y)
			{
				while (true)
				{
					XB_old = XB;
					XB = X;
					X.x = 2.0 * XB.x - XB_old.x;
					X.fit_fun(ff);
					X = HJ_trial(ff, X, s);
					if (X.y >= XB.y)
						break;
					if (solution::f_calls > Nmax)
						return XB;
				}
			}
			else
				s *= alpha;
			if (s < epsilon || Nmax < solution::f_calls) {
				XB.flag = 0;
				//cout << XB.x(0) << endl;
				//cout << " ";
				cout << XB.x(1) << endl;
				return XB;
			}
		}
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int* n = get_size(XB.x);
		//matrix D = ident_mat(n[0]);
		matrix D(n[0], n[0]);
		for (int i = 0; i < n[0]; ++i)
			D(i, i) = 1;
		solution X;
		for (int i = 0; i < n[0]; ++i)
		{
			X.x = XB.x + s * D[i];
			X.fit_fun(ff);
			if (X.y < XB.y)
				XB = X;
			else
			{
				X.x = XB.x - s * D[i];
				X.fit_fun(ff);
				if (X.y < XB.y)
					XB = X;
			}
		}
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_dim(x0);
		matrix l(n, 1), p(n, 1), s(s0), D(n, n);
		for (int i = 0; i < n; ++i)
			D(i, i) = 1;
		solution X(x0), Xt;
		X.fit_fun(ff, ud1, ud2);
		while (true)
		{
			//cout << X.x(0) << endl;
			//cout << " ";
			//cout << X.x(1) << endl;
			for (int i = 0; i < n; ++i)
			{
				Xt.x = X.x + s(i) * D[i];
				Xt.fit_fun(ff, ud1, ud2);
				if (Xt.y < X.y)
				{
					X = Xt;
					l(i) += s(i);
					s(i) *= alpha;
				}
				else
				{
					p(i)++;
					s(i) *= -beta;
				}
			}
			bool change = true;
			for (int i = 0; i < n; ++i)
				if (p(i) == 0 || l(i) == 0)
				{
					change = false;
					break;
				}
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = l(i);
				Q = D * Q;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + trans(Q[i]) * D[j] * D[j];
					v = (Q[i] - temp) / norm(Q[i] - temp);
					D.set_col(v, i);
				}
				s = s0;
				l = matrix(n, 1);
				p = matrix(n, 1);
			}
			double max_s = abs(s(0));
			for (int i = 1; i < n; ++i)
				if (max_s < abs(s(i)))
					max_s = abs(s(i));
			if (max_s < epsilon || solution::f_calls > Nmax) {
				X.flag = 0;
				return X;
			}
		}
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		//solution Xopt;
		double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
		solution X(x0), X1;
		matrix c0(2, new double[2]{ c,dc });
		while (true) {
			X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c0);
			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax) {
				X1.flag = 0;
				return X1;
			}
			X = X1;
			c0(0) = c0(0) * dc;
		}
		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		int n = get_len(x0);
		matrix D = ident_mat(n);
		int N = n + 1;
		solution* S = new solution[N];
		S[0].x = x0;		
		S[0].fit_fun(ff, ud1, ud2);

		for (int i = 1; i < N; ++i) {
			S[i].x = S[0].x + s * D[i - 1]; 
			S[i].fit_fun(ff, ud1, ud2);
		}

		solution PR, PE, PN;
		matrix pc;
		int i_min, i_max;

		while (true) {
			i_min = i_max = 0;
			for (int i = 1; i < N; ++i) {
				if (S[i].y(0) < S[i_min].y(0))
					i_min = i;
				if (S[i].y(0) > S[i_max].y(0))
					i_max = i;
			}

			pc = matrix(n, 1);

			for (int i = 0; i < N; ++i)
				if (i != i_max)
					pc = pc + S[i].x;

			pc = pc / (N - 1);
			PR.x = pc + alpha * (pc - S[i_max].x);
			PR.fit_fun(ff, ud1, ud2);

			if (PR.y(0) < S[i_max].y(0) && S[i_min].y(0) <= PR.y(0))
				S[i_max] = PR;				
			else if (PR.y(0) < S[i_min].y(0))	{
				PE.x = pc + gamma * (PR.x - pc);
				PE.fit_fun(ff, ud1, ud2);
				
				if (PE.y(0) < PR.y(0))
					S[i_max] = PE;
				else
					S[i_max] = PR;
			}
			else {
				PN.x = pc + beta * (S[i_max].x - pc);
				PN.fit_fun(ff, ud1, ud2);
				if (PN.y(0) < S[i_max].y(0))
					S[i_max] = PN;
				else {
					for (int i = 0; i < N; ++i)
						if (i != i_min) {
							S[i].x = delta * (S[i].x + S[i_min].x);
							S[i].fit_fun(ff, ud1, ud2);
						}
				}
			}
			double max_s = norm(S[i_min].x - S[0].x);

			for (int i = 1; i < N; ++i)
				if (max_s < norm(S[i_min].x - S[i].x))
					max_s = norm(S[i_min].x - S[i].x);

			if (max_s < epsilon)
				return S[i_min];
		}
		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		solution X, X1;
		X.x = x0;
		int n = get_len(x0);
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab;
		while (true) {
			//X.grad(gf);
			//d = -X.g;
			d = -X.grad(gf, ud1, ud2);
			if (h0 < 0) {
				//P[0] = X.x;
				//P[1] = d;
				P.set_col(X.x, 0);
				P.set_col(d, 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
				//cout << X1.x(0) << ";" << X1.x(1) << endl;
			}
			else {
				X1.x = X.x + h0 * d;
				//cout << X1.x(0) << ";" << X1.x(1) << endl;
			}
			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1, ud2);
				X1.flag = 0;
				return X1;
			}
			X = X1;
		}

		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab, beta;
		//X.grad(gf);
		//d = -X.g;
		d = -X.grad(gf, ud1, ud2);
		while (true) {
			if (h0 < 0)	{
				P.set_col(X.x, 0);
				P.set_col(d, 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
				//cout << X1.x(0) << ";" << X1.x(1) << endl;
			}
			else {
				X1.x = X.x + h0 * d;
				//cout << X1.x(0) << ";" << X1.x(1) << endl;
			}
			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1);
				X1.flag = 0;
				return X1;
			}
			X1.grad(gf);
			beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
			d = -X1.g + beta * d;
			X = X1;
		}
		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab;
		while (true) {
			X.grad(gf);
			X.hess(Hf);
			d = -inv(X.H) * X.g;
			if (h0 < 0) {
				P[0] = X.x;
				P[1] = d;
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
				//cout << X1.x(0) << ";" << X1.x(1) << endl;
			}
			else {
				X1.x = X.x + h0 * d;
				//cout << X1.x(0) << ";" << X1.x(1) << endl;
			}
			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1);
				X1.flag = 0;
				return X1;
			}
			X = X1;
		}

		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		double alfa = (sqrt(5) - 1) / 2;
		solution A, B, C, D;
		A.x = a;
		B.x = b;
		C.x = B.x - alfa * (B.x - A.x);
		C.fit_fun(ff, ud1, ud2);
		D.x = A.x + alfa * (B.x - A.x);
		D.fit_fun(ff, ud1, ud2);

		while (true) {
			if (C.y < D.y) {
				B = D;
				D = C;
				C.x = B.x - alfa * (B.x - A.x);
				C.fit_fun(ff, ud1, ud2);
			}
			else {
				A = C;
				C = D;
				D.x = A.x + alfa * (B.x - A.x);
				D.fit_fun(ff, ud1, ud2);
			}
			if (B.x - A.x < epsilon || solution::f_calls > Nmax) {
				A.x = (A.x + B.x) / 2;
				A.fit_fun(ff, ud1, ud2);
				A.flag = 0;
				return A;
			}
		}

		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		/*
		int n = get_len(x0);
		matrix D = ident_mat(n), * A = new matrix[2];
		solution X, P, h;
		X.x = x0;
		double* ab;

		while (true) {
			P = X;
			for (int i = 0; i < n; ++i) {
				A[0] = P.x;
				A[1] = D[i];
				ab = find_ab(ff, 0, 1, 1.2, Nmax, ud1, *A);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, *A);
				P.x = P.x + h.x * D[i];
			}
			if (norm(P.x - X.x) < epsilon || solution::f_calls > Nmax) {
				P.fit_fun(ff, ud1);
				return P;
			}
			for (int i = 0; i < n - 1; ++i)
				D.set_col(D[i + 1], i);

			D.set_col(P.x - X.x, n - 1);
			A[0] = P.x;
			A[1] = D[n - 1];
			ab = find_ab(ff, 0, 1, 1.2, Nmax, ud1, *A);
			h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, *A);
			X.x = P.x + h.x * D[n - 1];
		}
		
		//return Xopt;
		*/
		/**/
		solution Xopt;
		//Tu wpisz kod funkcji
		double V1[] = { 1,0 };
		double V2[] = { 0,1 };
		matrix e1(2, V1);
		matrix e2(2, V2);
		matrix e[2] = { e1,e2 };

		int n = get_len(x0);

		matrix d[2];
		matrix d1(n, 1);
		matrix d2(n, 1);

		solution x;
		solution p0, p, h[2];

		x.x = x0;

		p0.x = x0;
		p = p0;
		int i = 0;

		matrix P(n, 2);
		double* ab;

		d1 = e[0];
		d2 = e[1];
		d[0] = d1;
		d[1] = d2;

		do {
			p0 = x;
			for (int j = 1; j <= n; j++) {

				P.set_col(p.x, 0);
				P.set_col(d[j - 1], 1);

				ab = find_ab(ff, 0, 1, 1.2, Nmax, ud1, P);
				h[j - 1] = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);

				p.x = p.x + h[j - 1].x * d[j - 1];
			}

			if (norm(p.x - x.x) < epsilon) {
				x.fit_fun(ff, ud1);
				return Xopt = x;
			}

			for (int j = 1; j <= n - 1; j++)
				d[0] = d[1];

			d[1] = p.x - p0.x;

			P.set_col(p.x, 0);
			P.set_col(d[1], 1);

			ab = find_ab(ff, 0, 1, 1.2, Nmax, ud1, P);
			h[0] = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);

			p.x = p.x + h[0].x * d[1];
			x = p;

			++i;
		} while (solution::f_calls < Nmax);
		
		
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

double* find_ab(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0, 0 };

		int i = 0;
		solution X0, X1, X[3]; // X{ i-1 ; i ; i+1 }

		for (int j = 0; j < 3; j++) {
			matrix temp(0.);
			X[j].x = temp;
			X[j].fit_fun(ff, ud1, ud2);
		}

		matrix mx0(x0);
		X0.x = mx0;
		X1.x = X0.x + d;
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		if (X1.y == X0.y) {
			p[0] = X1.x(0);
			p[1] = X0.x(0) - d;
			return p;
		}

		if (X1.y > X0.y) {
			d = -d;
			X1.x = X0.x + d;
			X0.fit_fun(ff, ud1, ud2);
			X1.fit_fun(ff, ud1, ud2);

			if (X1.y >= X0.y) {
				p[0] = X1.x(0);
				p[1] = X0.x(0) - d;
				return p;
			}
		}

		do {
			if (solution::f_calls > Nmax) {
				throw ("Przekroczono maksymalna liczbe iteracji");
			}
			i++;
			X[0] = X[1];
			X[1] = X[2];
			X[2].x = X0.x + pow(alpha, i) * d;
			X[2].fit_fun(ff, ud1, ud2);
		} while (X[1].y <= X[2].y || i < 2);

		if (d > 0) {
			p[0] = X[0].x(0);
			p[1] = X[2].x(0);
		}
		else {
			p[0] = X[2].x(0);
			p[1] = X[0].x(0);
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion2(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		//Tu wpisz kod funkcji
		//mi liczebnosc populacji
		solution* P = new solution[mi + lambda];
		solution* Pm = new solution[mi];
		default_random_engine gen;
		gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
		normal_distribution<double> distr(0.0, 1.0);
		matrix IFF(mi, 1), temp(N, 2); //IFF macierz z przystosowaniami //temp - kopia osobnikia
		double r, s, s_IFF;
		double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5); //tau, tau1 - mutacja
		int j_min; // najlepsze rozwiazanie
		for (int i = 0; i < mi; ++i)
		{
			P[i].x = matrix(N, 2);
			for (int j = 0; j < N; ++j)
			{
				P[i].x(j, 0) = (limits(j, 1) - limits(j, 0)) * rand_mat(1, 1)() + limits(j, 0);
				P[i].x(j, 1) = sigma0(j);
			}
			P[i].fit_fun(ff, ud1, ud2);
			if (P[i].y < epsilon)
				return P[i];
			
		}
		while (true)
		{
			s_IFF = 0;
			for (int i = 0; i < mi; ++i)
			{
				IFF(i) = 1 / P[i].y();
				s_IFF += IFF(i);
			}
			for (int i = 0; i < lambda; ++i)
			{
				r = s_IFF * rand_mat(1, 1)();
				s = 0;
				for (int j = 0; j < mi; ++j)
				{
					s += IFF(j);
					if (r <= s)
					{
						P[mi + i] = P[j]; // j - wylosowany osobnik
						break;
					}
				}
			}
			//mutajca
			for (int i = 0; i < lambda; ++i)
			{
				r = distr(gen);
				for (int j = 0; j < N; ++j)
				{
					P[mi + i].x(j, 1) *= exp(tau1 * r + tau * distr(gen));
					P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * distr(gen);
				}
			}
			//krzyzowanie
			for (int i = 0; i < lambda; i += 2)
			{
				r = rand_mat(1, 1)();
				temp = P[mi + i].x;  //jeden z rodzicow
				P[mi + i].x = r * P[mi + i].x + (1 - r) * P[mi + i + 1].x;  //pierwszy potomek
				P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r) * temp;  //drugi potomek
			}
			//ocena osobnikow
			for (int i = 0; i < lambda; ++i)
			{
				P[mi + i].fit_fun(ff, ud1, ud2);
				if (P[mi + i].y < epsilon) 	//ocena rozwiazania
					return P[mi + i];

			}
			//wskazanie najelpszych osobnikow
			for (int i = 0; i < mi; ++i)
			{
				j_min = 0;
				for (int j = 1; j < mi + lambda; ++j)
					if (P[j_min].y > P[j].y)
						j_min = j;
				Pm[i] = P[j_min];
				P[j_min].y = 1e10;
			}
			for (int i = 0; i < mi; ++i)
				P[i] = Pm[i];  //P[i] najlepsza populacja
			if (solution::f_calls > Nmax)
				return P[0];
			
		}
		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
