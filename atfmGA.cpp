# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <chrono>
# include <cstring>
# include <numeric>
# include <vector>
# include <algorithm>

#include "atfmGA.h"


using namespace std;


//****************************************************************************80

// random generator function:
int myRandom(int range)
{
	//Returns a random value
	return rand() % range;
}

vector<int> ShuffleIndx(int N)
//****************************************************************************80
// 
//  Purpose:
//
//    ShuffleIndx performs a shuffling of elements in a vector. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    André Serrano
//
//  Parameters:
//
//    Input:	N			number of elements to shuffle
//
//    Output:	indexes		vector of shuffled indexes
//
{
	//Generates a random number between 0 and N correspondent to their own indexes
	srand(unsigned(time(0)));

	//vector of indexes to shuffle 
	vector<int> indexes;

	// set the indexes:
	for (int i = 0; i < N; ++i)
	{
		indexes.push_back(i); 
	}

	//Shuffle the indexes
	random_shuffle(indexes.begin(), indexes.end(), myRandom);

	return indexes;
}

void ShuffleX(int one, int two, double pGeneShuffle, vector<genotype> &population, int &nVar)
//****************************************************************************80
// 
//  Purpose:
//
//    ShuffleX performs uniform crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    André Serrano
//
//  Parameters:
//
//    Input: int one, two		 Indices of the two parents.
//			 pGeneShuffle		 Probability of gene shuffling inside the interval [0,1]
//			 population			 Genotype/GA space
//			 nVar				 Number of variables
//
//    Output:
//
{
	
	int geneIndxShuffle;
	int t;

	//Generate a vector with indexes shuffled randomnly 
	vector<int> geneindexes = ShuffleIndx(nVar);

	//number of genes to shuffle
	int nGeneShuffle = pGeneShuffle * nVar;

	//  Swap genes in random positions in the geneindexes vector.
	for (int i = 0; i < nVar; i++)
	{
		//Get the gene index to swap
		geneIndxShuffle = geneindexes[i];

		//Get temporary value of the first individual's gene
		t = population[one].gene[geneIndxShuffle];
		//Perform swapping
		population[one].gene[geneIndxShuffle] = population[two].gene[geneIndxShuffle];
		population[two].gene[geneIndxShuffle] = t;
	}

	return;
}

void crossover(int &seed, vector<genotype> &population, int &nVar, int &nSols, int gen)
//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//	  Modified by André Serrano
//
//  Parameters:
//
//    Input:	 seed			seed for the random number generator,
//				 population		genotype/GA space,
//				 nVar			number of variables,
//				 nSols			number of feasible solutions
//				 gen			generation index
//
{
	const double a = 0.0;
	const double b = 1.0;
	int mem;
	int one;
	int first = 0;
	double x;
	double randPshuffle;
	double Pxover;
	int maxGens = MAXGENS;

	//Dynamic Parameter Update (DPU)
	if (gen <= 0.25*(double)maxGens)
	{
		Pxover = 0.99;
	}
	else if (gen > 0.25*(double)maxGens && gen <= 0.75*(double)maxGens)
	{
		Pxover = 1.435 - 1.78*(double)gen / (double)maxGens;
	}
	else
	{
		Pxover = 0.1;
	}

	for (mem = 0; mem < POPSIZE; ++mem)
	{
		//Compute random number between 0 and 1
		x = r8_uniform_ab(a, b, seed);

		//Perform crossover if the probability of crossover [0,1] is superior than x
		if (x < Pxover)
		{
			++first;
			//Randomize the probability of gene's number to shuffle
			randPshuffle = r8_uniform_ab(0, 1.0, seed);

			if (first % 2 == 0) 
			{

				//one and mem are indexes of the parents who are going to mate 
				//ShuffleX(one, mem, randPshuffle, population, nVar); //UNCOMMENT TO USE

				//one and mem are indexes of the parents who are going to mate 
				//Xover1point(one, mem, seed, population, nVar); //UNCOMMENT TO USE

				Xover2point(one, mem, seed, population, nVar);

			}
			else
			{
				one = mem;
				
				//In case of populations with odd individuals
				if (mem == POPSIZE - 1)
				{
					//get the index of a random individual to mate with the last individual
					int last = i4_uniform_ab(0, POPSIZE - 2, seed);

					//ShuffleX(one, last, randPshuffle, population, nVar); //UNCOMMENT TO USE

					//Xover1point(one, last, seed, population, nVar); //UNCOMMENT TO USE

					Xover2point(one, last, seed, population, nVar);

				}
			}

		}
	}
	return;
}
//****************************************************************************80

void elitist(vector<genotype> &population, int &nSols, int &nvars)

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//	  Modified by André Serrano
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
	int i;
	double best; //lower fitness
	int best_mem;
	double worst; //higher fitness
	int worst_mem;

	printf("elistist\n");

	best = population[0].fitness;
	worst = population[0].fitness;

	for (i = 0; i < nSols - 1; ++i)
	{
		if (population[i + 1].fitness < population[i].fitness)
		{

			if (best <= population[i].fitness)
			{
				best = population[i].fitness;
				best_mem = i;
			}

			if (population[i + 1].fitness <= worst)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}

		}
		else
		{

			if (population[i].fitness <= best)
			{
				best = population[i].fitness;
				best_mem = i;
			}

			if (worst <= population[i + 1].fitness)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}

		}

	}
	// 
	//  If the best individual from the new population is better than 
	//  the best individual from the previous population, then 
	//  copy the best from the new population; else replace the 
	//  worst individual from the current population with the 
	//  best one from the previous generation                     
	//
	if (population[nSols-1].fitness >= best)
	{
		for (i = 0; i < nvars; i++)
		{
			population[nSols-1].gene[i] = population[best_mem].gene[i];
		}
		population[nSols-1].fitness = population[best_mem].fitness;
	}
	else
	{
		for (i = 0; i < nvars; i++)
		{
			population[worst_mem].gene[i] = population[nSols-1].gene[i];
		}
		population[worst_mem].fitness = population[nSols-1].fitness;
	}

	return;
}
//****************************************************************************80


int i4_uniform_ab(int a, int b, int &seed)
//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
	int c;
	const int i4_huge = 2147483647;
	int k;
	float r;
	int value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "I4_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}
	//
	//  Guarantee A <= B.
	//
	if (b < a)
	{
		c = a;
		a = b;
		b = c;
	}

	k = seed / 127773;

	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	r = (float)(seed)* 4.656612875E-10;
	//
	//  Scale R to lie between A-0.5 and B+0.5.
	//
	r = (1.0 - r) * ((float)a - 0.5)
		+ r   * ((float)b + 0.5);
	//
	//  Use rounding to convert R to an integer between A and B.
	//
	value = round(r);
	//
	//  Guarantee A <= VALUE <= B.
	//
	if (value < a)
	{
		value = a;
	}
	if (b < value)
	{
		value = b;
	}

	return value;
}
//****************************************************************************80

void initializePopulation(vector<genotype> &population, int &nVar)
//****************************************************************************80
// 
//  Purpose:
//
//    Initializes the genes to zero and its attributes 
//
//  Discussion:
//
//    It also initializes (to zero) all fitness values for each
//    member of the population. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//   Original version by Dennis Cormier and Sita Raghavan.
//   This C++ version by John Burkardt.
//	 Modified by André Serrano
//
//  Parameters:
//
//    Input: population		all individuals,
//		     nVar			number of variables.
//
//
{
	for (auto& p : population) //(p is an individual)
	{
		//p.infeasibility = 0.0;
		p.fitness = 0.0;
		p.rfitness = 0.0;
		p.cfitness = 0.0;

		//each individual/chromossome has size equal to the number of variables
		p.gene.resize(nVar);
			
		//Set all the genes to zero
		for (auto& gene : p.gene)
		{
			gene = 0;
		}
	}

	return;
}
//****************************************************************************80

void randomResetting(int &seed, vector<genotype> &population)
//****************************************************************************80
// 
//  Purpose:
//
//    ´Random Resetting mutation 
//
//  Discussion:
//
//    Randomly allocates departure slots to variables
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//	 André Serrano
//
//  Parameters:
//
//    Input: population		genotype/GA space,
//		     seed			random number.
//
//
{
	double pr;
	double x;
	double Pmutation;
	int maxGens = MAXGENS;

	for (auto& p : population)
	{
		for (auto& gene : p.gene)
		{
			x = r8_uniform_ab(0.0, 1.0, seed);

			if (x < PMUTATION)
			{
				//random integer number generated from 0 to 11
				gene = i4_uniform_ab(0, 11, seed);
			}
		}
	}

	return;
}


void mutate(int &seed, vector<genotype> &population)
//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a shifted biased mutation. 
//
//  Discussion:
//
//    NOT FINISHED
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    André Serrano
//
//  Parameters:
//
//    Input:	seed			 a seed for the random number generator,
//				population		 genotype/ GA space
//
{
	const double a = 0.0;
	const double b = 1.0;

	const double delayLB = 0.0;
	const double delayUB = 11.0;
	double x;
	double randPmut;
	double pBias;

	for (auto& p : population)
	{
		//pbias = [0.0, 0.8] will increase the probability of shifting to lower values
		pBias = r8_uniform_ab(0.0, 1.0, seed);
		
		randPmut = r8_uniform_ab(0.0001, PMUTATION, seed);
		for (auto& gene : p.gene)
		{
			x = r8_uniform_ab(a, b, seed);
			if (x < randPmut)
			{
				if (pBias < r8_uniform_ab(0.0, 1.0, seed))
				{
					gene = i4_uniform_ab(0, gene, seed);
				}
				else
				{
					gene = i4_uniform_ab(gene, 11, seed);
				}
			}
			////Compute a random number between 0 and 1
			//x = r8_uniform_ab(a, b, seed);
			////Randomize probability of mutation from 0.001 and an user defined value
			//randPmut = r8_uniform_ab(0.001, PMUTATION, seed);
			//if (x < randPmut)
			//{
			//	//Computes a random number in the delay start slots interval [0, 11]
			//	gene = round(r8_uniform_ab(delayLB, delayUB, seed)); //i4_uniform_ab(0, 11, seed)
			//}
		}
	}

	return;
}
//****************************************************************************80

void biasResetting(int &seed, vector<genotype> &population, vector<double> cdfStartSlot, int gen)
//****************************************************************************80
// 
//  Purpose:
//
//    Performs a biased resetting mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random biased value favouring startslots with lower delays
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    André Serrano
//
//  Parameters:
//
//    Input: seed			 a seed for the random number generator,
//			 population		 genotype/GA space,
//			 cdfStartSlot	 probability distribution for the startslots
//			 gen			 generation index
//
{
	double pr;
	double x;
	double Pmutation;
	int maxGens = MAXGENS;

	//Dynamic Parameter Update (DPU)
	//Pmutation = 0.1*(1 - (double)gen / (0.75*(double)maxGens)) + 0.001;
	if (gen <= 0.25*(double)maxGens)
	{
		Pmutation = 0.1;
	}
	else if (gen > 0.25*(double)maxGens && gen <= 0.75*(double)maxGens)
	{
		Pmutation = 0.1495 - 0.198*(double)gen / (double)maxGens;
	}
	else
	{
		Pmutation = 0.001;
	}

	for (auto& p : population)
	{
		for (auto& gene : p.gene)
		{
			x = r8_uniform_ab(0.0, 1.0, seed);

			if (x < Pmutation)
			{
				//random number generated from 1 to 0
				pr = r8_uniform_ab(0.0, 1.0, seed);

				//pr is in the first interval
				if (pr < cdfStartSlot[0])
				{
					gene = 0;
				}
				else
				{
					//pr is not in the first interval
					for (int k = 0; k <= 11; k++)
					{
						if (cdfStartSlot[k] <= pr && pr < cdfStartSlot[k + 1])
							gene = k + 1;  //gene = k + 1; Bias Resetting || gene = gene -/+ (k + 1); Creep
					}
				}
			}
		}
	}

	return;
}
//****************************************************************************80

void creepMutation(int &seed, vector<genotype> &population, vector<double> cdfStartSlot, int gen)
//****************************************************************************80
// 
//  Purpose:
//
//    Performs a biased resetting mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is shifted (added or subtrated) by a biased value favouring smaller changes/shifts
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    André Serrano
//
//  Parameters:
//
//    Input: seed			 a seed for the random number generator,
//			 population		 genotype/GA space,
//			 cdfStartSlot	 probability distribution for the startslots
//			 gen			 generation index
//
{
	double pr;
	double x;
	double addOrsub;
	double Pmutation;
	int maxGens = MAXGENS;

	//Dynamic Parameter Update (DPU)
	//Pmutation = 0.1*(1 - (double)gen / (0.75*(double)maxGens)) + 0.001;
	if (gen <= 0.25*(double)maxGens)
	{
		Pmutation = 0.1;
	}
	else if (gen > 0.25*(double)maxGens && gen <= 0.75*(double)maxGens)
	{
		Pmutation = 0.1495 - 0.198*(double)gen / (double)maxGens;
	}
	else
	{
		Pmutation = 0.001;
	}

	for (auto& p : population)
	{
		for (auto& gene : p.gene)
		{
			if (Pmutation < PMUTATION)
			{
				//random number generated from 1 to 0
				pr = r8_uniform_ab(0.0, 1.0, seed);

				//probability of adding or subtracting 50/50 chance
				addOrsub = r8_uniform_ab(0.0, 1.0, seed);

				//pr is in the first interval
				if (pr < cdfStartSlot[0])
				{
					//50 / 50 chance
					if (addOrsub < 0.01)
					{
						gene = gene + 1;  //Add
					}
					else
					{
						gene = gene - 1;  //Subtract
					}
				}
				else
				{
					//pr is not in the first interval
					for (int k = 0; k <= 10; k++)
					{
						if (cdfStartSlot[k] <= pr && pr < cdfStartSlot[k + 1])
						{
							//50 / 50 chance
							if (addOrsub < 0.5)
							{
								gene = gene + (k + 1);  //Add
							}
							else
							{
								gene = gene - (k + 1);  //Subtract
							}
						}
					}
				}
				//Make sure that the gene doesn't go out of boundaries [0;11]
				if (gene < 0)
				{
					gene = 0;
				}
				if (gene > 11)
				{
					gene = 11;
				}
			}
		}
	}

	return;
}
//****************************************************************************80

double r8_uniform_ab(double a, double b, int &seed)
//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
	int i4_huge = 2147483647;
	int k;
	double value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "R8_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	k = seed / 127773;

	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	value = (double)(seed)* 4.656612875E-10;

	value = a + (b - a) * value;

	return value;
}
//****************************************************************************80


void selector(int &seed, vector<genotype> &population, vector<genotype> &parentsPopulation, int &nSols)
//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is a Rank-Based Selection operator.
//
//  Discussion:
//
//    This selection strategy named Roulette Whell selection favours individuals with better fitnesses. The latter 
//	  have an higher chance to be picked as a parent and thus higher chance to pass its genetic material to its children
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//	  Modified by André Serrano
//
//  Parameters:
//
//    Input:	seed				 a seed for the random number generator,
//				population			 genotype/GA space,
//				parentsPopulation	 parents population which will procriate,
//				nSols				 number of feasible solutions.
//
{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	int mem;
	double p;
	double sum;

	//
	//  Find the total fitness of the population.
	//
	sum = 0.0;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		if (population[mem].fitness == 0)
		{
			sum = sum;
		}
		else
		{
			sum = sum + 1/population[mem].fitness; //MAXIMIZE: population[mem].fitness, MINIMIZE: 1/population[mem].fitness
		}
	}

	//
	//  Calculate the relative fitness of each member.
	//
	for (mem = 0; mem < POPSIZE; mem++)
	{
		if (population[mem].fitness == 0)
		{
			population[mem].rfitness = 0;
		}
		else
		{
			population[mem].rfitness = (1/population[mem].fitness) / sum; //MAXIMIZE: population[mem].fitness / sum, MINIMIZE: (1/population[mem].fitness) / sum
		}
	}

	// 
	//  Calculate the cumulative fitness.
	//
	population[0].cfitness = population[0].rfitness;
	for (mem = 1; mem < POPSIZE; mem++)
	{
		population[mem].cfitness = population[mem - 1].cfitness +
			population[mem].rfitness;
	}


	// 
	//  Select survivors using cumulative fitness. 
	//
	for (i = 0; i < POPSIZE; i++)
	{
		//spin the wheel by computing a random number from 0 to 1
		p = r8_uniform_ab(a, b, seed);

		//The wheel stopped before the first individual
		if (p < population[0].cfitness) //if(population[0].cfitness==0)
		{
			parentsPopulation[i] = population[0];
		}
		else
		{
			for (j = 0; j < POPSIZE; j++)
			{
				//The wheel stopped between two individuals
				if (population[j].cfitness <= p && p < population[j + 1].cfitness)
				{
					parentsPopulation[i] = population[j + 1];
				}
			}
		}
	}

	// 
	//  Overwrite the old population with the new one. Should I overwrite all the population?
	//
	for (i = 0; i < POPSIZE; i++)
	{
		population[i] = parentsPopulation[i];
	}

	return;
}
//****************************************************************************80

void tournament(int &seed, vector<genotype> &population, vector<genotype> parentsPopulation, int &nSols)
//****************************************************************************80
// 
//  Purpose:
//
//    TOURNAMENT is a Standard Proportional Selection operator.
//
//  Discussion:
//
//    This selection strategy named Tournament selection favours individuals with better fitnesses. The latter 
//	  have an higher chance to be picked as a parent and thus higher chance to pass its genetic material to its children
//	  Select the parents with the highest fitness value from the previous population. 
//	  In the tournament selection, there is no arithmetical computation based on the fitness value, but only comparison between individuals by fitness value.
//	  The number of the individuals taking part in the tournament is called tournament size.The selection performs following the steps as below :
//	  1) Randomly select several individuals from the population to take part in the tournament.Choose the individual that has the highest fitness value
//	  from the individuals selected above by comparing the fitness value of each individual.Then the chosen one is copied into the next generation of the
//	  population.
//	  2) Repeat step1 n times where n is the number of individuals of the population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    André Serrano
//
//  Parameters:
//
//    Input:	seed				 a seed for the random number generator,
//				population			 genotype/GA space,
//				parentsPopulation	 parents population which will procriate,
//				nSols				 number of feasible solutions.
//
{	
	//Initializing tournament population
	vector<genotype> tournamentpop(TOURNAMENTSIZE);

	for (int i = 0; i < POPSIZE; i++)
	{		
		//Initializing tournament population
		for (auto& t : tournamentpop)
		{
			t.fitness = 0.0;

			for (auto& gene : t.gene)
			{
				gene = 0;
			}
		}

		//Create the tournament population with random individuals from the current population
		for (int j = 0; j < TOURNAMENTSIZE; j++)
		{
			//Generates a random index of the current population
			int p = i4_uniform_ab(0, POPSIZE-1, seed); //100 - > nSols

			//Fill tournament population with random individuals
			tournamentpop[j] = population[p];
		}

		double winner = tournamentpop[0].fitness;
		int winner_ind = 0;

		//Get the index of the tournament's winner
		for (int k = 0; k < TOURNAMENTSIZE - 1; k++)
		{
			if (tournamentpop[k + 1].fitness < tournamentpop[k].fitness)
			{
				if (tournamentpop[k + 1].fitness <= winner)
				{
					winner_ind = k + 1;
				}
			}
			else
			{
				if (tournamentpop[k].fitness <= winner)
				{
					winner_ind = k;
				}
			}
		}

		//Populate the parents with the winners of the tournament 
		parentsPopulation[i] = tournamentpop[winner_ind];

		//Probability of win? like:
		/*
		z = random
		if (z < PWIN){
		   newpopulation[i] = tournamentpop[winner_ind];
		else
		   ix = round(random(1, TOURNAMENTSIZE)
		   newpopulation[i] = tournamentpop[ix]
		}
		*/	
	}

	//Overwrites the previous population with the selected parents. These individuals will proceed into the recombination.
	population = parentsPopulation;


	return;
}
//

void Xover1point(int father, int mother, int &seed, vector<genotype> &population, int &nVar)
//****************************************************************************80
// 
//  Purpose:
//
//    Xover1point performs one point crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//	  Modified by André Serrano
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input:  father, mother	 the indices of the two parents,
//			  population		 genotype/GA space,
//			  nVar				 number of feasible solutions.
//
{
	int i;
	int point;
	double t;

	// 
	//  Select the crossover point.
	//
	point = i4_uniform_ab(0, nVar - 1, seed);
	//
	//  Swap genes in positions 0 through POINT-1.
	//
	for (i = 0; i < point; i++)
	{
		t = population[father].gene[i];
		population[father].gene[i] = population[mother].gene[i];
		population[mother].gene[i] = t;
	}

	return;
}

void Xover2point(int father, int mother, int &seed, vector<genotype> &population, int &nVar)
//****************************************************************************80
// 
//  Purpose:
//
//    Xover2point performs two point crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2018
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//	  Modified by André Serrano
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input:  father, mother	 the indices of the two parents,
//			  population		 genotype/GA space,
//			  nVar				 number of feasible solutions.
//
{

	int point1;
	int point2;
	
	// 
	//  Select the crossover points.
	//
	point1 = i4_uniform_ab(0, nVar - 1, seed);
	point2 = i4_uniform_ab(0, nVar - 1, seed);

	//In case point1 is to the left of point2
	if (point1 > point2)
	{
		swap(point1, point2);
	}

	//Perform swaping
	swap_ranges(begin(population[father].gene) + point1,
		begin(population[father].gene) + point2,
		begin(population[mother].gene) + point1);

	return;
}
