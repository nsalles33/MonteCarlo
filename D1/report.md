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
* To study the system with varying temparature I am taking 100 spins and 500 sweeps for the temperature starting from 0 incresed to 40 at a step of 0.8 each time. I am taking the mean from the last 10 sweeps for energy and magnetization of teh system then calculating the mean and standard devaition for the values of last 10 sweeps and plotting them with respect to temperature.
### The varying enery of the system with temperature
* Here we can clearly observe a transition in the average energy in the system at temp `~ 9.5` where the system goes from one stable state to another state. After temp. 15 the energy satuarates and stops increasing.
![energy](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/energy.png)
### The varying magnetization of the system with temperature
* Here a simila rtrend is seen but in opposite direction. The average square of magnetization statys constant for the temperature closer to zero then as the temparature inceases there is more fluctuation among the spins of the system. There complete chaos when the temperature is too high and as a result the magnatization goes to zero. 
![magnetization](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/magnetization.png)

### Variation of the intended values through the sweeps
* Here my system has 1000 spins and the temp.is set 10 and the system goes through 200 sweeps and at the end of each sweep I am calculating the average energy and the average square of magnatization of the system and plotting them against the sweep number.
* From all data through the sweeps the average energy is calculated to be ` -5600.205175` with a standard deviation of ` 198.793850`.
* The average square of magnetization is calculated to be ` 0.720993` with a standard deviation of `0.030410`.
* We see the values are far away from the equalibrium in the biginning but as the number of sweeps incease, it travels towards the statble equalibrium for the given parameters.
![magnetization](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/for_sweep.png)
![magnetization](https://github.com/rjtkp/MonteCarlo/blob/master/D1/ising/forsweeps_mag.png)

### Speed Up

* For larger N values it makes sense to parallelize the code using the hybrid OpenMP-mpi paradigm so that proccesses between different nodes can talk with each other using the message passing interface which is restricted for the OpenMp thread as indeprendent node dont have shared memory among them.

* With changing the numbe rof threds there was no improvement in the performence as openblas already utilizes the multitreading atrchitecture underneath. The code speands most of the time in the loop of sweeps and inside that loop the multiplication for the energy calculation and update takes the most amont of time  where blas uses the most number threads available for the system. The contribution from all the other parts is negligible. So, no visible speed up is observed.
* When I use `export OMP_NUM_THREADS=1` and run I get average time per sweeps as `0.597459 sec` but with `export OMP_NUM_THREADS=20`, I get `0.578473 sec` per sweed on average for 1000 spins and 100 walker.
#### PS: I have done all tests in Ulysses with specifications of
```
processor	: 19
vendor_id	: GenuineIntel
cpu family	: 6
model		: 62
model name	: Intel(R) Xeon(R) CPU E5-2680 v2 @ 2.80GHz
stepping	: 4
microcode	: 1064
cpu MHz		: 2799.939
cache size	: 25600 KB
physical id	: 1
siblings	: 10
core id		: 12
cpu cores	: 10
apicid		: 56
initial apicid	: 56
fpu		: yes
fpu_exception	: yes
```


