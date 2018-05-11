# HPC-Monte Carlo problems

## Problem 1: Infinite variance distributions

![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/1.jpg)
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/2.jpg)
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/3.jpg)

## Problem 2: Sum of uniformly distributed random variables

![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/4.jpg)
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/5.jpg)



### The histogram obtained from counting the number occurances.
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/histogram.png)

### The probability for n=4 compared with the actual gaussian probability.
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/prob_4.png)


### The probability for n=6 compared with the actual gaussian probability.
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/prob_6.png)


### The probability for n=4 and n=6 compared with the actual gaussian probability.
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/4_6_probability.png)

### The probabilities for different n values
![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/probability.png)

### The error in the probability 
```math_def
The relative error % = 100*[prob_n(i) - gaussian(i)]/gaussian(i)
```
* Where the prob_n represents the probability of n sums of the random variable.
The value of maximun realative error percentage starts with a really high value but decrese very sharply and then get flatened at the value of n is incresed. The minimum value error percentage is obtained for the of n of 210.

![](https://github.com/rjtkp/MonteCarlo/blob/master/D1/error_percent.png)


## Problem 3: Ising model with arbitrary range

* I have written Monte code of Carlo algorithm which generates samples of spin configurations distributed according to the Boltzmann weight `exp(−H(~σ )/T )` and by using a finite number M of independent walkers with independent seed initialization whcich can be found [here](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/montecarlo.c). 
### The varying enery of the system with temperature
* Here we can clearly observe a transition in the average energy in the system at temp `~ 9.5` where the system goes from one stable state to another state. After temp. 15 the energy satuarates and stops increasing.
![energy](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/energy.png)
### The varying magnetization of the system with temperature
* Here a simila rtrend is seen but in opposite direction. The average square of magnetization statys constant for the temperature closer to zero then as the temparature inceases there is more fluctuation among the spins of the system. There complete chaos when the temperature is too high and as a result the magnatization goes to zero. 
![magnetization](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/magnetization.png)

### Variation of the intended values through the sweeps
* We see the values are far away from the equalibrium in the biginning but as the number of sweeps incease, it travels towards the statble equalibrium for the given parameters.
![magnetization](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/for_sweep.png)
![magnetization](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/forsweeps_mag.png)


![speed up](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/speed_up.png)
* I have done different studies on the systems with changing the temparatue, changing the size  

* For larger N values it makes sense to parallelize the code using the hybrid OpenMP-mpi paradigm so that proccesses between different nodes can talk with each other using the message passing interface which is restricted for the OpenMp thread as indeprendent node dont have shared memory among them.

