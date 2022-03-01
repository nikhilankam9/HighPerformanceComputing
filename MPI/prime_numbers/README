Performance Evaluation and Optimising the Seive of Eratosthenes to calculate prime numbers in the order of 10^10 using the OpenMPI library in C.

Parallel Algorithm:
Decomposition: The input size is decomposed into equal parts based on the number of processes available and the current prime calculated is broadcasted to all the other processes.
First divisible value in each process: If the first value is not divisible by the prime then the remainder(r = firstValue%prime) is used to calculate the first divisible prime number using the formula startIndex = prime - remainder. The meaning of this formula is to say how many places should the startValue move to be divisible by prime.

Eliminating even integers:
Fact: Leverages the fact that 2 is the only even prime number.
Decomposition: The input size is reduced by half by the above fact and the prime number can be deduced from the index of the array using the formula p = 2*i + 1
Prime calculation: The same approach is followed as in the first algorithm. 
First divisible value in each process: The distance between the lowValue and the next divisible prime is calculated using the formula startVal = lowVal + prime - remainder since only odd integers are considered if the startVal from the previous formula turned out to be even then prime is again added to generate the first divisible value(since even + odd = odd).

Elimination broadcast:
The main communication between all the processes in the previous methods is to broadcast the current sieving prime, to eliminate that step all the primes <sqrt(n) are precalculated and stored in each of the processes.

Cache hit:
To leverage the cache, a block of the input array in each process is processed completely before stepping onto the next block.

![Benchmark](https://github.com/nikhilankam9/high_perf_computing/blob/prime_numbers/MPI/prime_numbers/benchmark.png)

By using various techniques we could see an improvement of around 7 fold.
Credit: Execution tested using the [Foundry](https://itrss.mst.edu/cluster/foundry/) cluster at MST Rolla, and the experimentation is done solely for educational purposes.