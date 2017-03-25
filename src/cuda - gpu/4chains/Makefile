make: 
	mv MCMC_4D.c MCMC_4D.cu
	mv MCMC_8D.c MCMC_8D.cu
	nvcc MCMC_4D.cu -o MCMC_4D.out
	nvcc MCMC_8D.cu -o MCMC_8D.out
clean:
	rm -rf MCMC*.out
