
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include<curand.h>
#include<curand_kernel.h>
#include<time.h>

const double pi = 3.1415926535897;
const double TUNED_VAR = 0.447;

double sampleNormal();

__global__ void addKernel1(double *y_c, double *y_now, double *temp1)
{
  int idx = blockDim.x * threadIdx.y + threadIdx.x;
  int c;
  temp1[idx] = ((exp(-(pow(((y_c[idx]-1) /1),2))/2))/1/sqrt(2*pi))/((exp(-(pow(((y_now[idx]-1) /1),2))/2))/1/sqrt(2*pi));
  for ( c = 1 ; c <= 32767 ; c++ )
    {}
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
  double y_c[10] = {1,1,1,1,1,1,1,1,1,1};
  double y_now[10] = {1,1,1,1,1,1,1,1,1,1};
  double temp1[10] = { 0 };
  double TOSS_A_DICE = 0.0;
  double * y_c_cuda = NULL;
  double * y_now_cuda = NULL;
  double * temp1_cuda = NULL;
  int accepted = 0;
  int rejected = 0;
  int size =10;
  int i = 0;
  int j = 0;
  float t1, t2;
  t1=clock();
  
  //memory alloc for gpu
  cudaMalloc((void**)&y_c_cuda, size * sizeof(double));
  cudaMemcpy(y_c_cuda, y_c, size * sizeof(double), cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&y_now_cuda, size * sizeof(double));
  cudaMemcpy(y_now_cuda, y_now, size * sizeof(double), cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&temp1_cuda, size * sizeof(double));
  
  srand (time(NULL));
  
  while(i<10000){
	//addKernel <<<1, size >>> (y_c_cuda, y_now_cuda,time(NULL));
	for(j=0;j<size;j++){
		y_c[j] = y_now[j]+(TUNED_VAR*sampleNormal());
	}
	cudaMemcpy(y_c_cuda, y_c, size * sizeof(double), cudaMemcpyHostToDevice);
	
	TOSS_A_DICE = (rand() % 10000) / 10000.0;
	addKernel1 <<<1, size >>> (y_c_cuda, y_now_cuda, temp1_cuda);
	
	cudaMemcpy(temp1, temp1_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	
	if(TOSS_A_DICE < calcMult(temp1, size))
	{
		cudaMemcpy(y_now_cuda, y_c_cuda, size * sizeof(double), cudaMemcpyDeviceToDevice);
		accepted=accepted+1;
	} 
	else 
	{  
		rejected=rejected+1;
	}
	cudaMemcpy(y_c, y_c_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(y_now, y_now_cuda, size * sizeof(double), cudaMemcpyDeviceToHost);
	//print(y_now, size);
	i++;
  }
  
  printf("accepted: %d\n rejected: %d\n",accepted, rejected);
  t2 = clock();
  printf("Exec time: %lf\n",(t2-t1)/double(CLOCKS_PER_SEC));
   
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