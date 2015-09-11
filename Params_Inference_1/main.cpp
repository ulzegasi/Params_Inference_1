#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <float.h>
#include <array>
#include <random>
#include <algorithm>
#include <functional>
#include "mdfun.h"
using namespace std;

// Useful directories
const string dir = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/";
const string dir3 = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/input_data/";

// Read time interval from file
bool tdat( vector<double> & t_limits )
{
	ifstream ifs(dir+"t.dat");
	if (! ifs)
	{
		cerr << "\nInvalid file name or location "
			<< dir+"t.dat" << " ... Aborting...\n\n";
		return false;
	}
	else
	{	
		int range = 5002;
		string line;
		int pos = 0, ix = 0;
		string col1; double col2;
		
		while (getline(ifs, line) && ix < range) // read one line from ifs
		{
			if (ix == 1 || ix == (range-1))
			{
				istringstream iss(line); // access line as a stream
				iss >> col1 >> col2;
				t_limits[pos++]=col2;
			}
			++ix;
		}
		return true;
	}
}

// Read input data from file
bool indat( vector<double> & y )
{
	ifstream ifs( dir3 + "y_n10_K50_g02_s10_sinr.dat" );
	if (! ifs)
	{
		cerr << "\nInvalid file name or location "
			<< dir3 + "y_n10_K50_g02_s10_sinr.dat" << " ... Aborting...\n\n";
		return false;
	}
	else
	{	
		string line;
		double yvalue;
		while (getline(ifs, line))   // read one line from ifs
		{
			istringstream iss(line); // access line as a stream
			while (iss >> yvalue)
				y.push_back(yvalue);
		}
		return true;
	}
}

// Read "true" system realization from file
bool Sdat( vector<double> & S )
{
	ifstream ifs( dir3 + "St_n10_K50_g02_s10_sinr.dat" );
	if (! ifs)
	{
		cerr << "\nInvalid file name or location "
			<< dir3 + "St_n10_K50_g02_s10_sinr.dat" << " ... Aborting...\n\n";
		return false;
	}
	else
	{	
		string line;
		double Svalue;
		while (getline(ifs, line)) // read one line from ifs
		{
			istringstream iss(line); // access line as a stream
			iss >> Svalue;
			S.push_back(Svalue);
		}
		return true;
	}
}

// Generate a vector of normally distributed random numbers with mean mu (default 0) and standard deviation std (default 1)
void nrand(vector<double> & vec, double mu = 0.0, double std = 1.0)
{
	static mt19937_64 engine(1807201122102013);
	// static mt19937_64 engine(time(NULL));
	static normal_distribution<double> norm_dist(mu,std);
	for (int ix = 0; ix < vec.size(); ++ix)
		vec[ix] = norm_dist(engine);
}

// Generate a random number from uniform ditribution between a and b (default 0 and 1) 
double urand(double a = 0.0, double b = 1.0)
{
	static mt19937_64 engine(1807201122102013);
	// static mt19937_64 engine(time(NULL));
	static uniform_real_distribution<double> unif_dist(a,b);
	double res = unif_dist(engine);
	return res;
}

// Sum of elements in a vector
template<typename elemType> elemType vsum(const vector<elemType> & vec)
{
	elemType sum_of_elems = 0;
	for_each(vec.begin(), vec.end(), [&](elemType elem){sum_of_elems += elem;});
	return sum_of_elems;
}

// Sums two vectors element-wise
template<typename elemType> vector<elemType> vsum(const vector<elemType> & vec1, const vector<elemType> & vec2)
{
	vector<elemType> vec_plus_vec(vec1.size());
	transform(vec1.begin(), vec1.end(), vec2.begin(), vec_plus_vec.begin(), plus<elemType>());
	return vec_plus_vec;
}

// Square root of elements of a vector
template<typename elemType> vector<double> vsqrt(const vector<elemType> & vec)
{
	vector<double> sqrt_vec(vec.size());
	vector<double>::iterator iter = sqrt_vec.begin();
	for_each(vec.begin(), vec.end(), [&](elemType elem){
		*iter = sqrt(elem); ++iter;});
	return sqrt_vec;
}

// Square of elements of a vector
template<typename elemType> vector<elemType> vsquare(const vector<elemType> & vec)
{
	vector<elemType> vec2(vec.size());
	vector<elemType>::iterator iter = vec2.begin();
	for_each(vec.begin(), vec.end(), [&](elemType elem){
		*iter = pow(elem,2.0); ++iter;});
	return vec2;
}

// Multiplies a vector by a constant
template<typename elemType> vector<elemType> vtimes(elemType val, const vector<elemType> & vec)
{
	vector<elemType> vec_times_val(vec.size());
	vector<elemType>::iterator iter = vec_times_val.begin();
	for_each(vec.begin(), vec.end(), [&](elemType elem){
		*iter = val*elem; ++iter;});
	return vec_times_val;
}

// Multiplies two vectors element-wise
template<typename elemType> vector<elemType> vtimes(const vector<elemType> & vec1, const vector<elemType> & vec2)
{
	vector<elemType> vec_times_vec(vec1.size());
	transform(vec1.begin(), vec1.end(), vec2.begin(), vec_times_vec.begin(), multiplies<elemType>());
	return vec_times_vec;
}

// Divides a vector by a constant
template<typename elemType> vector<elemType> vdiv(elemType val, const vector<elemType> & vec)
{
	vector<elemType> vec_div_val(vec.size());
	vector<elemType>::iterator iter = vec_div_val.begin();
	for_each(vec.begin(), vec.end(), [&](elemType elem){
		*iter = elem/val; ++iter;});
	return vec_div_val;
}

// Divides two vectors element-wise
template<typename elemType> vector<elemType> vdiv(const vector<elemType> & vec1, const vector<elemType> & vec2)
{
	vector<elemType> vec_div_vec(vec1.size());
	transform(vec1.begin(), vec1.end(), vec2.begin(), vec_div_vec.begin(), divides<elemType>());
	return vec_div_vec;
}



int main(){

	// ******************************************************* //
	// *******************  Preliminaries ******************** //
	// ******************************************************* //

	// **************  System parameters *************** //
	const int nparams = 2;
	const int n = 10;
	const int j = 30;
	const int N = int (n*j+1);

	// **************  Time parameters *************** //
	double T, dt = 0;
	vector<double> t_limits(2);
	vector<double> t(N);
	vector<int> ty(n+1);

	if (!tdat(t_limits))   // Read time points from file
		return -1;

	T  = t_limits.back()-t_limits.front();  // Total time
	dt = T/(double)(N-1);					// Time step

	t[0] = t_limits.front();				// Time points (N)
	for (int ix = 1; ix < N; ++ix)
		t[ix] = t[ix-1] + dt;

	for (int ix = 0; ix < n+1; ++ix)	   // Indexes of measurement points (= n+1 boundary beads)
		ty[ix] = ix*j+1;

	// **************  HMC parameters *************** //
	const int nsample_burnin = 0;         // Number of points in the MCMC
	const int nsample_eff = 1000;
	const int nsample = nsample_eff + nsample_burnin;

	const double dtau = 0.25;  // MD time step
	const int n_napa = 3;      // Number of NAPA steps

	const double true_K = 50.0;   // Retention time
	const double true_gam = 0.2;  // Dimensionless noise parameter
	const double true_bet = sqrt(T*true_gam/true_K);
	const double sigma = 0.10;     // Measurement noise
	vector<double> true_theta(nparams);  // Parameters to be inferred (beta, tau)
	true_theta[0] = true_bet;
	true_theta[1] = true_gam;


	double K = 200.0;    // Initial state
	double gam = 0.5;
	double bet = sqrt(T*gam/K);                                   
	vector<double> theta(nparams);
	vector<double> * theta_pt = &theta;
	theta[0] = bet;
	theta[1] = gam;

	double K_min   = 0.0;   // Parameter limits
	double gam_min = 0.0;
	double K_max   = 1000.0;
	double gam_max = 5.0;

	// **************  Generate data  *************** // 
	vector<double> r(N);         // Input signal (rain)
	vector<double> lnr_der(N,0); // Log-derivative of the rain
	vector<double> y;
	vector<double> S;
	vector<double> long_t;
	vector<double> q_init(n+1);
	vector<double> bq(n+1);
	vector<double> q(N,0);
	
	for (int ix = 0; ix < N; ++ix)
		r[ix] = pow(sin(t[ix]/100),2.0) + 0.1;

	for (int ix = 0; ix < N-1; ++ix)
		lnr_der[ix] = (log(r[ix+1])-log(r[ix]))/dt;

	if (!indat(y))   // Read data points from file
		return -1;

	if (!Sdat(S))   // Read "true" system realization from file
		return -1;

	double tiny_dt = T/((S.size())-1);	// Time step for the "true" system realization 
	long_t.push_back(t[0]);				    // Time points (N)
	for (int ix = 1; ix < S.size(); ++ix)
		long_t.push_back(long_t[ix-1] + tiny_dt);

	for (int ix = 0; ix < n+1; ++ix)
	{
		q_init[ix] = (1/true_bet)*log(y[ix]/r[ty[ix]-1]);
		bq[ix] = true_bet*q_init[ix];
	}

	for (int s = 0; s < n; ++s)
	{
		double step = (q_init[s+1]-q_init[s])/j;
		q[ty[s]-1] = q_init[s];
		q[ty[s+1]-1] = q_init[s+1];
		for (int ix = 1; ix < j; ++ix)
			q[ty[s]-1+ix] = q[ty[s]-1] + ix*step;
	}

	// **************************************************************** //
	// ******************  Hamiltonian Monte Carlo  ******************* //
	// **************************************************************** //
	
	// **************  Transformations q -> u  *************** // 
    // (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)

	vector<double> u(N);

	for (int s = 0; s < 1; ++s)
	{
		u[s*j] = q[s*j];
		for (int k = 2; k < j+1; ++k)
			u[s*j+k-1] = q[s*j+k-1] - ( (k-1)*q[s*j+k] + q[s*j] )/k;
	}
	u[N-1] = q[N-1];

	// **************  Init chains storage  *************** //
	vector< vector<double> > theta_sample;
	vector< vector<double> > u_sample;
	vector<double> energies;

	theta_sample.push_back(theta);
	u_sample.push_back(u);

	// **************  Init containers for temp values  *************** //
	vector<double> theta_save(nparams);
	vector<double> u_save(N);
	vector<double> * theta_save_pt = &theta_save;
	vector<double> * u_save_pt = &u_save;

	// **************  Init container for masses and momenta  *************** //
	vector<double> mp(nparams+N);
	vector<double> sqrt_mp(nparams+N);
	vector<double> p(nparams+N);
	vector<double> * mp_pt = &mp;
	vector<double> * sqrt_mp_pt = &sqrt_mp;
	vector<double> * p_pt = &p;

	// **************  Init vector container for normal random numbers  *************** //
	vector<double> randvec(nparams+N); 

	// **************  Energies  *************** //
	double H_old, H_new;
	double accept_prob;

	// **************  Containers for execution times  *************** //
	vector<double> time_respa(nsample_eff);
	vector<double> time_respa_s(nsample_eff);
	vector<double> time_respa_f(nsample_eff);

	// **************  Masses (burn-in)  *************** //
	double m_bdy = 10.0;     // m = m_q / dt
	double m_stg = 1.0;      // we assume m_q prop. to dt ==> m = costant     
	array<double,nparams> m_theta;
	m_theta[0]=1.0;
	m_theta[1]=1.0;

	for (int ix = 0; ix < nparams; ++ix)
		(*mp_pt)[ix] = m_theta[ix];
	for (int s = 1; s < n+1; ++s)
	{
		(*mp_pt)[nparams+(s-1)*j] = m_bdy;                       
		for (int k = 2; k < j+1; ++k) 
			(*mp_pt)[nparams+(s-1)*j+k-1] = m_stg;                   
	}
	(*mp_pt)[nparams+N-1] = m_bdy;  

	*sqrt_mp_pt = vsqrt(*mp_pt);

	// ************** Loading potentials, derivatives and Napa *************** //
	cout << "\nLoading potentials, derivatives and napa...\n-------------------------------\n";
	cout << "\nTo be implemented...\n";
	// reload("$dir/ParInf_Fun_AD_napa.jl")

	// ************** HMC loops *************** //
	int reject_counter_burnin = 0;

	cout << "\nStarting HMC loops (burn-in)...\n---------------------------------\n";
	clock_t tinit = clock();

	for (int counter = 1; counter <= nsample_burnin; ++counter)
		cout << "\nHMC loops here ... to be implemented ...\n";

	// ************** Redefinition of masses (effective HMC loops) *************** //
	double m_bdy_burnin = m_bdy;
	double m_stg_burnin = m_stg;
	array<double,nparams> m_theta_burnin = m_theta;

	m_bdy = 720.0;       // m = m_q / dt
	m_stg = 130.0;       // we assume m_q prop. to dt ==> m = costant     
	m_theta[0] = 150.0;  // refers to beta
	m_theta[1] = 150.0;  // refers to gamma

	for (int ix = 0; ix < nparams; ++ix)
		(*mp_pt)[ix] = m_theta[ix];
	for (int s = 1; s < n+1; ++s)
	{
		(*mp_pt)[nparams+(s-1)*j] = m_bdy;                       
		for (int k = 2; k < j+1; ++k) 
			(*mp_pt)[nparams+(s-1)*j+k-1] = m_stg;                   
	}
	(*mp_pt)[nparams+N-1] = m_bdy;

	*sqrt_mp_pt = vsqrt(*mp_pt);

	// ************** Effective HMC loops *************** //
	int reject_counter = 0;
	cout << "\nStarting effective HMC loops...\n---------------------------------\n";

	for (int counter = (nsample_burnin + 1); counter <= nsample; ++counter)
	{
		clock_t t0 = clock();
		// Sample momenta
		nrand(randvec);
		*p_pt = vtimes(*sqrt_mp_pt,randvec);

		// Calculate energy
		// H_old = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
		H_old = vsum(vdiv(vsquare(*p_pt),vtimes(2.0,*mp_pt))); // + V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u)
		energies.push_back(H_old);

		if ( _isnan(H_old) != 0) {
			cerr << "\n\nIteration " << counter << " --> energy values diverged...\n\n";
			return -1;}

		// Save current state
		(*theta_save_pt) = theta;
		(*u_save_pt) = u;

		// MD Integration:
		clock_t t1 = clock();
		for (int counter_napa = 1; counter_napa <= n_napa; ++counter_napa)
		{
			napa(theta, u, p, mp, counter, n, j, N, nparams, T, dt, dtau, m_stg, m_bdy, 
				time_respa_f, time_respa_s, nsample_burnin);
		}
		time_respa[counter-nsample_burnin-1] = ((float)(clock()-t1)/CLOCKS_PER_SEC);

		// Caluclate energy of proposal state => Metropolis accept/reject
		if ( theta[1] <= gam_min || theta[1] >= gam_max || (T*theta[1]/pow(theta[0],2) <= K_min) || (T*theta[1]/pow(theta[0],2) >= K_max) )
		{
			theta = (*theta_save_pt); u = (*u_save_pt);
			reject_counter += 1;
		}
		else if ( find_if(u.begin(), u.end(), [&](double el){return (el > 10.0);}) != u.end() || 
			find_if(u.begin(), u.end(), [&](double el){return (el < -10.0);}) != u.end() )
		{
			theta = (*theta_save_pt); u = (*u_save_pt);
			reject_counter += 1;
		}
		else
		{
			H_new = vsum(vdiv(vsquare(*p_pt),vtimes(2.0,*mp_pt))); // + V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u);
			accept_prob = min(1.0,exp(H_old-H_new));
			if (urand() > accept_prob)
			{
				theta = (*theta_save_pt); u = (*u_save_pt);
				reject_counter += 1;
			}
		}

		theta_sample.push_back(theta);
		u_sample.push_back(u);

		if (counter%100 == 0)
			cout << "\n" << counter << " loops completed in " << ((float)(clock()-tinit)/CLOCKS_PER_SEC) << " seconds\n";
	}

	// ***************************************************************** //
	// **********************  End of HMC loops  *********************** //
	// ***************************************************************** //
	cout << "\n\n";
	cout << "Run completed in " << ((float)(clock()-tinit)/CLOCKS_PER_SEC) << " seconds\n";
	cout << "NAPA cycles in " << vsum(time_respa) << " seconds\n";
	cout << "Slow NAPA in " << vsum(time_respa_s) << " seconds\n";
	cout << "Fast NAPA in " << vsum(time_respa_f) << " seconds\n";
	
	// **************  Transformations q -> u  *************** // 
    // (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)

	cout << "\n\n";
	return 0;
}