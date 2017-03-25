
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include<curand.h>
#include<curand_kernel.h>
#include<time.h>


const double pi = 3.1415926535897;
const double TUNED_VAR = 0.755;

double sampleNormal();
double expectancy(double x[], int size);

__global__ void addKernel1(double *y_c1, double *y_now1, double *y_c2, double *y_now2 ,double *y_c3, double *y_now3, double *y_c4, double *y_now4 ,double *temp1, double *temp2, double *temp3, double *temp4)
{
  int idx = blockDim.x * threadIdx.y + threadIdx.x;
	temp1[idx] = ((exp(-(pow(((y_c1[idx]-1) /1),2))/2))/1/sqrt(2*pi))/((exp(-(	pow(((y_now1[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	temp2[idx] = ((exp(-(pow(((y_c2[idx]-1) /1),2))/2))/1/sqrt(2*pi))/((exp(-(pow(((y_now2[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	temp3[idx] = ((exp(-(pow(((y_c3[idx]-1) /1),2))/2))/1/sqrt(2*pi))/((exp(-(	pow(((y_now3[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	temp4[idx] = ((exp(-(pow(((y_c4[idx]-1) /1),2))/2))/1/sqrt(2*pi))/((exp(-(pow(((y_now4[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	__syncthreads();
}

__global__ void addKernel2(double *y_now1, double *y_now2, double *y_now3, double *y_now4, double *temp1, double *temp2, double *temp3, double *temp4)
{
  int idx = blockDim.x * threadIdx.y + threadIdx.x;
	temp1[idx] = ((exp(-(pow(((y_now1[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	temp2[idx] = ((exp(-(pow(((y_now2[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	temp3[idx] = ((exp(-(pow(((y_now3[idx]-1) /1),2))/2))/1/sqrt(2*pi));
	temp4[idx] = ((exp(-(pow(((y_now4[idx]-1) /1),2))/2))/1/sqrt(2*pi));
  __syncthreads();
}

double calcMult(double a[], int size){
	double result = 1;
	int i =0;
	for(i = 0;i<size;i++){
		result = a[i]*result;
	}
	return result;
}

void print(double a[], int size){
	int i;
	for(i=0;i<size;i++){
		printf("%lf  ", a[i]);
	}
	printf("\n");
}

int main()
{
  int size = 8;
  double y_c1[8] = {1,1,1,1.1,1,1,1};
  double y_c2[8] = {1,1,1,1,1,1,1,1};
  double y_c3[8] = {1,1,1,1,1,1,1,1};
  double y_c4[8] = {1,1,1,1,1,1,1,1};
  double y_now1[8] = {1,1,1,1,1,1,1,1};
  double y_now2[8] = {1,1,1,1,1,1,1,1};
  double y_now3[8] = {1,1,1,1,1,1,1,1};
  double y_now4[8] = {1,1,1,1,1,1,1,1};
  double temp1[8] = { 0 };
  double temp2[8] = { 0 };
  double temp3[8] = { 0 };
  double temp4[8] = { 0 };
  double exchange[8] = { 0 };
  double exchanged = 0;
  double EX_MCMC = 0;
  double TOSS_A_DICE = 0.0;
  double * y_c1_cuda = NULL;
  double * y_now1_cuda = NULL;
  double * y_c2_cuda = NULL;
  double * y_now2_cuda = NULL;
  double * y_c3_cuda = NULL;
  double * y_now3_cuda = NULL;
  double * y_c4_cuda = NULL;
  double * y_now4_cuda = NULL;
  double * temp1_cuda = NULL;
  double * temp2_cuda = NULL;
  double * temp3_cuda = NULL;
  double * temp4_cuda = NULL;
  double criteria = 0;
  int accepted1 = 0;
  int rejected1 = 0;
  int accepted2 = 0;
  int rejected2 = 0;
  int accepted3 = 0;
  int rejected3 = 0;
  int accepted4 = 0;
  int rejected4 = 0;
  int i = 0;
  int j = 0;
  float t1, t2;
  t1=clock();
  
  //memory alloc for gpu
  cudaMalloc((void**)&y_c1_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_c2_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_now1_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_now2_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_c3_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_c4_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_now3_cuda, size * sizeof(double));
  cudaMalloc((void**)&y_now4_cuda, size * sizeof(double));
  cudaMalloc((void**)&temp1_cuda, size * sizeof(double));
  cudaMalloc((void**)&temp2_cuda, size * sizeof(double));
  cudaMalloc((void**)&temp3_cuda, size * sizeof(double));
  cudaMalloc((void**)&temp4_cuda, size * sizeof(double));
  
  srand (time(NULL));
  
  while(criteria<10000){
	for(j=0;j<size;j++){
		y_c1[j] = y_now1[j]+(TUNED_VAR*sampleNormal());
		y_c2[j] = y_now2[j]+(TUNED_VAR*sampleNormal());
		y_c3[j] = y_now1[j]+(TUNED_VAR*sampleNormal());
		y_c4[j] = y_now2[j]+(TUNED_VAR*sampleNormal());
	}
	
	cudaMemcpy(y_c1_cuda, y_c1, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_c2_cuda, y_c2, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_now2_cuda, y_now2, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_now1_cuda, y_now1, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_c3_cuda, y_c3, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_c4_cuda, y_c4, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_now4_cuda, y_now4, size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(y_now3_cuda, y_now3, size * sizeof(double), cudaMemcpyHostToDevice);
	
	
	addKernel1 <<<1, size >>> (y_c1_cuda, y_now1_cuda, y_c2_cuda, y_now2_cuda,y_c3_cuda, y_now3_cuda, y_c4_cuda, y_now4_cuda,temp1_cuda, temp2_cuda, temp3_cuda, temp4_cuda);
	
	cudaMemcpy(temp1, temp1_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp2, temp2_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp3, temp3_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp4, temp4_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if(TOSS_A_DICE < calcMult(temp1, size))
	{
		cudaMemcpy(y_now1_cuda, y_c1_cuda, size * sizeof(double), cudaMemcpyDeviceToDevice);
		accepted1=accepted1+1;
	} 
	else 
	{  
		rejected1=rejected1+1;
	}
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if(TOSS_A_DICE < calcMult(temp2, size))
	{
		cudaMemcpy(y_now2_cuda, y_c2_cuda, size * sizeof(double), cudaMemcpyDeviceToDevice);
		accepted2=accepted2+1;
	} 
	else 
	{  
		rejected2=rejected2+1;
	}
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if(TOSS_A_DICE < calcMult(temp3, size))
	{
		cudaMemcpy(y_now3_cuda, y_c3_cuda, size * sizeof(double), cudaMemcpyDeviceToDevice);
		accepted3=accepted3+1;
	} 
	else 
	{  
		rejected3=rejected3+1;
	}
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if(TOSS_A_DICE < calcMult(temp4, size))
	{
		cudaMemcpy(y_now4_cuda, y_c4_cuda, size * sizeof(double), cudaMemcpyDeviceToDevice);
		accepted4=accepted4+1;
	} 
	else 
	{  
		rejected4=rejected4+1;
	}
	
	
	addKernel2 <<<1, size >>> (y_now1_cuda, y_now2_cuda, y_now3_cuda, y_now4_cuda, temp1_cuda, temp2_cuda, temp3_cuda, temp4_cuda);
	
	cudaMemcpy(temp1, temp1_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp2, temp2_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp3, temp3_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp4, temp4_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	
	double p1 = calcMult(temp1, size);
	double p2 = calcMult(temp2, size);
	double p3 = calcMult(temp3, size);
	double p4 = calcMult(temp4, size);
	
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if (TOSS_A_DICE < (p2*pow(p2,0.5625)/p1*pow(p2,0.5625))){
		for(j=0;j<size;j++){
			exchange[j] = y_now2[j];
			y_now2[j] = y_now1[j];
			y_now1[j] = exchange[j];
			exchanged++;
		}
	}
	
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if (TOSS_A_DICE < (pow(p3,0.0625)*pow(p4,0.25)/pow(p3,0.25)*pow(p4,0.0625))){
		for(j=0;j<size;j++){
			exchange[j] = y_now2[j];
			y_now2[j] = y_now1[j];
			y_now1[j] = exchange[j];
			exchanged++;
		}
	}
	
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	if (TOSS_A_DICE < (pow(p1,0.25)*p3/p1*pow(p3,0.0625))){
		for(j=0;j<size;j++){
			exchange[j] = y_now2[j];
			y_now2[j] = y_now1[j];
			y_now1[j] = exchange[j];
			exchanged++;
		}
	}
	
	EX_MCMC=EX_MCMC + expectancy(y_now1, size);
	criteria=criteria+1;
	i++;
  }
  
  t2 = clock();
  printf("accepted1: %d\nrejected1: %d\n",accepted1, rejected1);
  printf("accepted2: %d\nrejected2: %d\n",accepted2, rejected2);
  printf("accepted3: %d\nrejected3: %d\n",accepted3, rejected3);
  printf("accepted4: %d\nrejected4: %d\n",accepted4, rejected4);
  printf("exchanged: %lf\n",exchanged);
  printf("Converged Value:%lf\n", EX_MCMC/i);
  printf("Exec time: %lf\n",(t2-t1)/double(CLOCKS_PER_SEC));

  cudaFree(y_c1_cuda);
  cudaFree(y_now1_cuda);
  cudaFree(y_c2_cuda);
  cudaFree(y_now2_cuda);
  cudaFree(y_c3_cuda);
  cudaFree(y_now3_cuda);
  cudaFree(y_now4_cuda);
  cudaFree(y_c4_cuda);
  cudaFree(temp1_cuda);
  cudaFree(temp2_cuda);
  cudaFree(temp3_cuda);
  cudaFree(temp4_cuda);
  
  return 0;
}

double sampleNormal() {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r <= 0 || r > 1) return sampleNormal();
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

double expectancy (double x[], int size) 
{
    double result = 1;
    int i=0;
    for(i=0; i<size;i++){
        result = result*x[i];
    }
    return result;
}