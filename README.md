## SpinMCFork.jl

This is a fork of the package [SpinMC.jl](https://github.com/fbuessen/SpinMC.jl). It adds a recursive tuning of the parallel tempering temperature distribution that converges to a set of temperatures where exchange acceptance rates between replicas are approximately the same (as described in ["Markov Chain Monte Carlo Simulations and Their Statistical Analysis"](https://www.worldscientific.com/worldscibooks/10.1142/5602), Section 6.3).

Usage of the package remains completly the same (refer to the original package), to tune parallel tempering temperature set add the following kwargs to the MonteCarlo constructor:
* `recursionIterations`: number of iterations the recursive tuning algorithm will perform, ideally, the point where the temperature set has converged (default: 0). 
* `recursionRate`: Rate at which the temperatures will be tuned, in units of sweeps (default: 1000).

Experimentation with these parameters is required for every specific problem, as these kind of algorithms can easily become unstable.
