#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std;
#include "Burgers.h"

#define AREA 4.0
#define Plus 1
#define Minus 0
#define Forward 1
#define Backward 0
#define e 10.0e-8

// public functions
Burgers::Burgers()
{
	ReadInput();
	NFX = AREA / dx;
	
	s_f.resize(NFX + 4);
	s_b.resize(NFX + 4);
	Q_p.resize(NFX + 4);
	X.resize(NFX + 4); // X[0] & X[1] & X[NFX + 2] & X[NFX + 3] for ghost cell
	if ((InitCond == 1) || (InitCond == 2)) X[0] = -0.5*AREA - 1.5*dx;
	else if (InitCond == 3) X[0] = -1.5*dx;
	else cout << "!!!Wrong initial Condition!!!\n";
	for (int i = 1; i < NFX + 4; ++i)
	{
		X[i] = X[0] + dx*(double)i;
	}
	Q = Initial();
	cout << "----------Input Parameters----------\n";
	cout << "	Initial Condition : " << InitCond << "\n";
	cout << "	Limiter Condition : " << LimitCond << "\n";
	cout << "	Target Time       : " << TargetT << "\n";
	cout << "	Grid size         : " << dx << "\n";
	cout << "	Number of grids   : " << NFX << "\n";
	cout << "------------------------------------\n";
}

Burgers::~Burgers()
{

}

std::vector<double> Burgers::RHS()
{
	BC();
	SaveQ();
	s_f = ShockSpeed(Forward);
	s_b = ShockSpeed(Backward);
	dt = TimeStep();
	std::vector<double> rhs;
	rhs.resize(NFX + 4);
	for (int i = 2; i < NFX + 2; ++i)
	{
		rhs[i] = - (1.0 / dx)*(AdvSpeed(Plus, i)*(Q_p[i] - Q_p[i - 1]) + AdvSpeed(Minus, i)*(Q_p[i + 1] - Q_p[i])) - (1.0 / dx)*(Flux(Forward, i, dt) - Flux(Backward, i, dt));
	}

	return rhs;
}

std::vector<double> Burgers::Progress()
{
	BC();
	SaveQ();
	s_f = ShockSpeed(Forward);
	s_b = ShockSpeed(Backward);
	dt = TimeStep();
	for (int i = 2; i < NFX + 2; ++i)
	{
		Q[i] = Q_p[i] - (dt / dx)*(AdvSpeed(Plus, i)*(Q_p[i] - Q_p[i - 1]) + AdvSpeed(Minus, i)*(Q_p[i + 1] - Q_p[i])) - (dt / dx)*(Flux(Forward, i, dt) - Flux(Backward, i, dt));
	}

	return Q;
}

std::vector<double> Burgers::Progress(const double dT)
{
	BC();
	SaveQ();
	s_f = ShockSpeed(Forward);
	s_b = ShockSpeed(Backward);
	dt = dT;
	for (int i = 2; i < NFX + 2; ++i)
	{
		Q[i] = Q_p[i] - (dt / dx)*(AdvSpeed(Plus, i)*(Q_p[i] - Q_p[i - 1]) + AdvSpeed(Minus, i)*(Q_p[i + 1] - Q_p[i])) - (dt / dx)*(Flux(Forward, i, dt) - Flux(Backward, i, dt));
	}

	return Q;
}

void Burgers::WriteOutput(string index) const
{
	ofstream file;
	string FileName = "./output/";
	FileName += "Burgers_";
	string Scheme = "";
	switch (InitCond)
	{
	case 1:
		FileName += "square_";
		break;
	case 2:
		FileName += "expansion_";
		break;
	case 3:
		FileName += "sine_";
		break;
	default:
		FileName += "TypeError_";
	}
	switch (LimitCond)
	{
	case 1:
		FileName += "upwind";
		Scheme += "upwind";
		break;
	case 2:
		FileName += "minmod";
		Scheme += "minmod";
		break;
	case 3:
		FileName += "superbee";
		Scheme += "superbee";
		break;
	case 4:
		FileName += "MC";
		Scheme += "MC";
		break;
	case 5:
		FileName += "vanLeer";
		Scheme += "van Leer";
		break;
	default:
		FileName += "TypeError";
	}
	FileName += "_";
	FileName += index;
	FileName += ".plt";
	file.open(FileName, ios::trunc);
	if (file.is_open())
	{
		cout << "output file open\n";
		file << "Variables = X, Q\n";
		file << "Zone t = \"" << Scheme << "\", i=" << NFX << ", f=point\n";
		for (int i = 2; i < NFX + 2; ++i)
		{
			file << X[i] << "\t" << Q[i] << "\n";
		}
	}
	else cout << "output file error\n";

	file.close();
}

// protected functions

void Burgers::ReadInput()
{
	ifstream file;
	string buff;
	string FileName = "./input.inp";
	file.open(FileName);
	if (file.is_open())
	{
		cout << "----------Input File Open-----------\n";
		file >> buff >> buff >> buff >> buff >> buff >> InitCond 
			>> LimitCond >> TargetT >> dx >> CFL;
	}
	else cout << "-------Cannot Find Input File-------\n";
	file.close();
}

std::vector<double> Burgers::Initial()
{
	std::vector<double> U;
	U.resize(NFX + 4);
	switch (InitCond)
	{
	case 1:
		for (int i = 0; i < NFX + 4; ++i)
		{
			if (X[i] <= 0.0) U[i] = 2.0;
			else U[i] = 0.0;
		}
		break;
	case 2:
		for (int i = 0; i < NFX + 4; ++i)
		{
			if (X[i] <= 0.0) U[i] = -1.0;
			else U[i] = 1.0;
		}
		break;
	case 3:
		for (int i = 0; i < NFX + 4; ++i)
		{
			//if ((X[i] <= 3.0) && (X[i] >= 1.0)) U[i] = sin(M_PI*(X[i] - 1.0));
			if (X[i] <= 2.0) U[i] = sin(M_PI*X[i]);
			else U[i] = 0.0;
		}
		break;
	}

	return U;
}

void Burgers::SaveQ()
{
	for (int i = 0; i < NFX + 4; ++i)
	{
		Q_p[i] = Q[i];
	}
}

std::vector<double> Burgers::ShockSpeed(bool Direc)
{
	std::vector<double> s;
	s.resize(NFX + 4);
	if (Direc == Forward)
	{
		for (int i = 0; i < NFX + 3; ++i)
		{
			if (abs(Q_p[i + 1] - Q_p[i]) < e) s[i] = Q_p[i];
			else s[i] = 0.5*(Q_p[i + 1] + Q_p[i]);
		}
	}
	else if (Direc == Backward)
	{
		for (int i = 1; i < NFX + 4; ++i)
		{
			if (abs(Q_p[i] - Q_p[i - 1]) < e) s[i] = Q_p[i];
			else s[i] = 0.5*(Q_p[i] + Q_p[i - 1]);
		}
	}
	else cout << "error in ShockSpeed\n";

	return s;
}

double Burgers::PiFunction(bool Direc, int n) const
{
	double Pi = 0.0;
	if (LimitCond == 1) return Pi;

	double theta = 0.0;
	if (Direc == Forward)
	{
		if (abs(Q_p[n + 1] - Q_p[n]) < e)
		{
			if (s_f[n] >= 0.0) return PiExcept(Q_p[n], Q_p[n - 1]);
			else return PiExcept(Q_p[n + 2], Q_p[n + 1]);
		}
		else
		{
			if (s_f[n] >= 0.0) theta = (Q_p[n] - Q_p[n - 1]) / (Q_p[n + 1] - Q_p[n]);
			else theta = (Q_p[n + 2] - Q_p[n + 1]) / (Q_p[n + 1] - Q_p[n]);
		}		
	}
	else if (Direc == Backward)
	{
		if (abs(Q_p[n] - Q_p[n - 1]) < e)
		{
			if (s_b[n] >= 0.0) return PiExcept(Q_p[n - 1], Q_p[n - 2]);
			else return PiExcept(Q_p[n + 1], Q_p[n]);
		}
		if (s_b[n] >= 0.0) theta = (Q_p[n - 1] - Q[n - 2]) / (Q_p[n] - Q_p[n - 1]);
		else theta = (Q_p[n + 1] - Q_p[n]) / (Q_p[n] - Q_p[n - 1]);
	}
	else cout << "error in PiFunction\n";

	std::vector<double> List;
	switch (LimitCond)
	{
	case 2: // minmod
		Pi = Minmod(1.0, theta);
		break;
	case 3: // superbee
		List.clear();
		List.push_back(0.0);
		List.push_back(min(1.0, 2.0*theta));
		List.push_back(min(2.0, theta));
		Pi = *std::max_element(List.begin(),List.end());
		break;
	case 4: // MC
		List.clear();
		List.push_back(0.5*(1.0 + theta));
		List.push_back(2.0);
		List.push_back(2.0*theta);
		Pi = max(0.0, *min_element(List.begin(), List.end()));
		break;
	case 5: // van Leer
		Pi = (theta + abs(theta)) / (1.0 + abs(theta));
		break;
	default:
		cout << "wrong limiter type\n";
	}

	return Pi;
}

double Burgers::AdvSpeed(bool sign, int n) const
{
	if (sign == Plus)
	{
		if (s_b[n] > 0.0) return s_b[n];
		else return 0.0;
	}
	else if (sign == Minus)
	{
		if (s_f[n] < 0.0) return s_f[n];
		else return 0.0;
	}
	else cout << "error in AdvSpeed\n";
}

double Burgers::TimeStep() const
{
	double maxQ = 0.0;
	for (int i = 2; i < NFX + 2; ++i)
	{
		maxQ = max(maxQ, abs(Q[i]));
	}

	double timestep = CFL * dx / maxQ;
	
	return timestep;
}

double Burgers::Flux(bool Direct, int n, double dt) const
{
	if (Direct == Forward)
	{
		return 0.5*abs(s_f[n])*(1.0 - (dt / dx)*abs(s_f[n]))*(Q[n + 1] - Q[n])*PiFunction(Forward, n);
	}
	else if (Direct == Backward)
	{
		return 0.5*abs(s_b[n])*(1.0 - (dt / dx)*abs(s_b[n]))*(Q[n] - Q[n - 1])*PiFunction(Backward, n);
	}
	else cout << "error in Flux\n";
}

double Burgers::PiExcept(double Q1, double Q2) const
{
	if (abs(Q1 - Q2) < e) return 1.0;
	else if ((Q1 - Q2) >= e)
	{
		switch (LimitCond)
		{
		case 2:
			return 1.0;
			break;
		case 3:
			return 2.0;
			break;
		case 4:
			return 2.0;
			break;
		case 5:
			return 2.0;
			break;
		default:
			return 0.0;
		}
	}
	else return 0.0;
}

double Burgers::Minmod(double a, double b) const
{
	if ((a*b > 0) && (abs(a) < abs(b))) return a;
	if ((a*b > 0) && (abs(a) > abs(b))) return b;
	else return 0;
}

void Burgers::BC()
{
	/*Q[0] = Q[NFX];
	Q[1] = Q[NFX + 1];
	Q[NFX + 2] = Q[2];
	Q[NFX + 3] = Q[3];*/
	switch (InitCond)
	{
	case 1:
		Q[0] = 2.0;
		Q[1] = 2.0;
		Q[NFX + 2] = 0.0;
		Q[NFX + 3] = 0.0;
		break;
	case 2:
		Q[0] = -1.0;
		Q[1] = -1.0;
		Q[NFX + 2] = 1.0;
		Q[NFX + 3] = 1.0;
		break;
	case 3:
		Q[0] = 0.0;
		Q[1] = 0.0;
		Q[NFX + 2] = 0.0;
		Q[NFX + 3] = 0.0;
		break;
	}
}


