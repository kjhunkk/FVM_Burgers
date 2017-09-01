#ifndef BURGERS_H
#define BURGERS_H

#include <vector>
#include <string>

class Burgers
{
public:
	Burgers();
	~Burgers();

	// get target time from input variables
	inline double get_TargetTime() const { return TargetT; };

	inline double getDT() const { return dt; };

	inline double getNFX() const { return NFX; };

	inline std::vector<double> getQ() const { return Q; };

	inline void setQ(std::vector<double> input) { Q = input; };

	// calculate RHS
	std::vector<double> RHS();

	// update solution variables with timestep()
	std::vector<double> Progress();

	// update solution variables with const dt / timestep
	std::vector<double> Progress(const double);

	// export solution variables / file index
	void WriteOutput(string) const;

protected:
	// variables
	int InitCond;
	int LimitCond;
	int NFX;
	double TargetT;
	double CFL;
	double dx;
	double dt;

	// cell quantities
	std::vector<double> Q;
	std::vector<double> Q_p;
	std::vector<double> s_f;
	std::vector<double> s_b;
	std::vector<double> X;

	// functions
	// read input file
	void ReadInput();

	// initialize solution variables
	std::vector<double> Initial();

	// save current solution variable to previous solution
	void SaveQ();

	// calculate shock speed / direction(forward, backward)
	std::vector<double> ShockSpeed(bool);

	// calculate Pi function / direction, index
	double PiFunction(bool, int) const;

	// calculate advection speed / sign(puls, minus), index
	double AdvSpeed(bool, int) const;

	// calculate time step from CFL & grid size
	double TimeStep() const;

	// calculate modified flux function with limiter / direction(forward, backward), index, time step
	double Flux(bool, int, double) const;

	// return when exception occured / e.g INF, NAN / Q1, Q2
	double PiExcept(double, double) const;

	// calculate minmod function / Q1, Q2
	double Minmod(double, double) const;

	// apply boundary condition
	void BC();
};

#endif BURGERS_H