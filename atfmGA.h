#pragma once
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <numeric>
# include <vector>


using namespace std;
// 
//  Change any of these parameters to match your needs 
//
# define MAXGENS 1			   //Number of generations/iterations
# define TOURNAMENTSIZE 3      //Number of contestants in tournament selection______Recomendaded values: 2,3
# define PXOVER 0.99           //Probability of crossover___________________________Choose value inside this interval [0,1]
# define PMUTATION 0.05        //Probability of mutation____________________________Choose value inside this interval [0,1]
# define INITIALSOLSSIZE 1	   //Size of the initial population_____________________INITIALSOLSSIZE must be equal or bigger than POPSIZE
# define POPSIZE 1			   //Size of the generational population
# define ELITRATIO 0.2		   //Ratio of elistim___________________________________Choose value inside this interval [0,1]

//
//  Each GENOTYPE is a member of the population, with
//  gene: a vector of variables,
//  fitness: the fitness
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//

//Population data structure
struct genotype
{
	vector<int> gene; // integer representation
	//double infeasibility;
	double fitness;
	double rfitness;
	double cfitness;
};

//Random value generator functions
int i4_uniform_ab(
	int					a,
	int					b,
	int					&seed
	);

double r8_uniform_ab(
	double				a, 
	double				b, 
	int					&seed
	);

//Selection strategies
void selector(
	int					&seed,
	vector<genotype>	&population,
	vector<genotype>	&parentsPopulation,
	int					&nSols
	);

void tournament(
	int					&seed,
	vector<genotype>	&population,
	vector<genotype>	parentsPopulation,
	int					&nSols
	);

void elitist(
	vector<genotype>	&population,
	int					&nSols, 
	int					&nvars
	);

//Mutation Operators
void mutate(
	int					&seed,
	vector<genotype>	&population
	);

void randomResetting(
	int					&seed,
	vector<genotype>	&population
	);

void biasResetting(
	int					&seed,
	vector<genotype>	&population, 
	vector<double>		cdfStartSlot,
	int					gen
	);

void creepMutation(
	int					&seed,
	vector<genotype>	&population,
	vector<double>		cdfStartSlot, 
	int					gen
	);

//Crossover operators
void crossover(
	int					&seed,
	vector<genotype>	&population,
	int					&nVar, 
	int					&nSols,
	int					gen
	);

void ShuffleX(
	int					one,
	int					two,
	double				pGeneShuffle,
	vector<genotype>	&population,
	int					&nVar
	);

void Xover1point(
	int					father,
	int					mother,
	int					&seed,
	vector<genotype>	&population, 
	int					&nVar
	);

void Xover2point(
	int					father,
	int					mother,
	int					&seed,
	vector<genotype>	&population, 
	int					&nVar
	);

//Other functions
int myRandom(
	int					range
	);

vector<int> ShuffleIndx(
	int					N
	);

void initializePopulation(
	vector<genotype>	&population,
	int					&nVar
	);
