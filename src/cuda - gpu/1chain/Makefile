make: 
	mv MCMC_4D.c MCMC_4D.cu
	mv MCMC_8D.c MCMC_8D.cu
	mv MCMC_10D.c MCMC_10D.cu
	nvcc MCMC_4D.cu -o MCMC_4D.out
	nvcc MCMC_8D.cu -o MCMC_8D.out
	nvcc MCMC_10D.cu -o MCMC_10D.out
clean:
	rm -rf MCMC*.out
