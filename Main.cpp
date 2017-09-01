#include <iostream>
using namespace std;
#include "Burgers.h"

#define e 10.0e-8

int main()
{
	double time;
	double dt = 0.0;
	Burgers solver;
	std::vector<double> L;
	std::vector<double> Q1;
	std::vector<double> Q2;
	std::vector<double> Q;

	int NFX = solver.getNFX();
	double TargetTime = solver.get_TargetTime();
	solver.WriteOutput("Initial");

	Q = solver.getQ();
	
	for (time = 0.0; (time - TargetTime + dt) < e; time += dt)
	{
		dt = 0.0;
		Q = solver.getQ();
		L = solver.RHS();
		dt = solver.getDT();

		Q1.clear();
		Q1 = Q;
		for (int i = 2; i < NFX + 2; ++i)
		{
			Q1[i] = Q[i] + dt*L[i];
		}

		solver.setQ(Q1);
		L = solver.RHS();
		//dt += solver.getDT();

		Q2.clear();
		Q2 = Q;
		for (int i = 2; i < NFX + 2; ++i)
		{
			Q2[i] = 0.75*Q[i] + 0.25*(Q1[i] + dt*L[i]);
		}

		solver.setQ(Q2);
		L = solver.RHS();
		//dt += solver.getDT();

		for (int i = 2; i < NFX + 2; ++i)
		{
			Q[i] = (Q[i] + 2.0*Q2[i] + 2.0*dt*L[i]) / 3.0;
		}

		solver.setQ(Q);
	}
	// last step
	dt = TargetTime - time;
	if (dt > e)
	{
		Q = solver.getQ();
		L = solver.RHS();

		Q1.clear();
		Q1 = Q;
		for (int i = 2; i < NFX + 2; ++i)
		{
			Q1[i] = Q[i] + dt*L[i];
		}

		solver.setQ(Q1);
		L = solver.RHS();
		//dt += solver.getDT();

		Q2.clear();
		Q2 = Q;
		for (int i = 2; i < NFX + 2; ++i)
		{
			Q2[i] = 0.75*Q[i] + 0.25*(Q1[i] + dt*L[i]);
		}

		solver.setQ(Q2);
		L = solver.RHS();
		//dt += solver.getDT();

		for (int i = 2; i < NFX + 2; ++i)
		{
			Q[i] = (Q[i] + 2.0*Q2[i] + 2.0*dt*L[i]) / 3.0;
		}

		solver.setQ(Q);
		time += dt;
	}
	
	cout << "current time = " << time << "\n";
	solver.WriteOutput("result");
}