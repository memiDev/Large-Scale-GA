# Large-Scale_GA
A genetic algorithm developed for a large scale binary problems


Genetic Algorithms (GAs) have the ability to deliver a ’good-enough’ solution ’fast-enough’. This makes GAs attractive for use in solving optimization problems. The reasons why GAs are needed are as follows. In NP-Hard problem the most powerful computing systems take an immense amount of time to solve the problem to optimality or to near-optimality.
In such a scenario, GAs prove to be an efficient tool to provide usable near-optimal solutions in a short amount of time.

Traditional calculus based methods work by starting at a random point and by moving in the direction of the gradient, until the top of the hill is reached. This technique is efficient and works very well for single-peaked objective functions like the cost function in linear regression. But, in most real-world situations, a very complex problem arises called as landscapes, which are made of many peaks and many valleys, which causes such methods to fail, as they suffer from an inherent tendency of getting stuck at the local optima, see figure bellow.

GAs have several advantages over other metaheuristics such as:
• Search efficiently in problem with large spaces and large number of parameters involved;
• Robust with respect to the complexity of the search problem;
• Use a population of solution instead of searching only one solution at a time.
However, they also present some limitations such as:
• Fitness values are calculated repeatedly which might be computationally expensive for some problems;
• Due to their stochastic nature, there are no guarantees on the optimality or the quality of the solution;
• A good understanding of the problem is required in order to properly implement, otherwise GA may not converge to the optimal or near-optimal solution.
