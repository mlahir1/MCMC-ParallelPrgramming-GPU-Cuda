#include <iostream>
#include <stdlib.h>     //srand, rand 
#include <time.h>       // time 
#include "math.h"       // for RAND, and rand, log and sqrt
#include <ctime>
using namespace std;

double sampleNormal();  // function declaration
double mypdf(double x[], int size);
double expectancy(double x[], int size);

const double pi = 3.1415926535897;

int main()
{

float exe_time=0;
float t1=clock();
double TOSS_A_DICE;
double VAR_TUNE[4] = {0.458,0.458,0.458,0.458};

double p1, p2, p1_2, p2_2, p3 ,p3_2, p4, p4_2;
int accepted1=0;
int rejected1=0;
int accepted2=0;
int rejected2=0;
int accepted3=0;
int rejected3=0;
int accepted4=0;
int rejected4=0;
int exchanged = 0;
double EX_MCMC = 0;
int c_size = 4;
int d_size =8;
double y_now[c_size][d_size];
double y_c[c_size][d_size];
double exchange[c_size];
int i=0;
int j=0;
int criteria=0;
//initialize random seed
srand (time(NULL));

for(i=0;i<c_size;i++){
    for(j=0;j<d_size;j++){
        y_c[i][j] = 1;
        y_now[i][j] = 1;
    }
}
while (criteria<100000) {
    // chain 1 and the main 1    
    for (i=0;i<c_size;i++){
		for(j=0;j<d_size;j++){
			y_c[i][j] = y_now[i][j] + VAR_TUNE[j]*sampleNormal();
		}
	}
    TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
	//cout << TOSS_A_DICE <<"  "<< mypdf(y_c[0], d_size)/mypdf(y_now[0], d_size)<< endl;
    if (TOSS_A_DICE < (mypdf(y_c[0], d_size)/mypdf(y_now[0], d_size))) 
    {
        for (i=0;i<d_size;i++){
            y_now[0][i] = y_c[0][i];
        }
        accepted1=accepted1+1;
    } 
    else 
    {  
        rejected1=rejected1+1;
    }
    
    TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
    if (TOSS_A_DICE < (pow(mypdf(y_c[1], d_size),0.5625)/pow(mypdf(y_now[1], d_size),0.5625))) 
    {
        for (i=0;i<d_size;i++){
            y_now[1][i] = y_c[1][i];
        }
        accepted2=accepted2+1;

    } 
    else 
    {  
        rejected2=rejected2+1;
    }
    TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
    if (TOSS_A_DICE < (pow(mypdf(y_c[2], d_size),0.25)/pow(mypdf(y_now[2], d_size),0.25))) 
    {
        for (i=0;i<d_size;i++){
            y_now[2][i] = y_c[2][i];
        }
        accepted3=accepted3+1;

    } 
    else 
    {  
        rejected3=rejected3+1;
    }
    TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
    if (TOSS_A_DICE < (pow(mypdf(y_c[3], d_size),0.0625)/pow(mypdf(y_now[3], d_size),0.0625))) 
    {
        for (i=0;i<d_size;i++){
            y_now[3][i] = y_c[3][i];
        }
        accepted4=accepted4+1;
    } 
    else 
    {  
        rejected4=rejected4+1;
    }

    if(accepted1>rejected1){VAR_TUNE[0] = VAR_TUNE[0]*1.1;}
    else{VAR_TUNE[0]=VAR_TUNE[0]*0.9;}
    if(accepted2>rejected2){VAR_TUNE[1] = VAR_TUNE[1]*1.1;}
    else{VAR_TUNE[1]=VAR_TUNE[1]*0.9;}
    if(accepted3>rejected3){VAR_TUNE[2] = VAR_TUNE[2]*1.1;}
    else{VAR_TUNE[2]=VAR_TUNE[2]*0.9;}
    if(accepted4>rejected4){VAR_TUNE[3] = VAR_TUNE[3]*1.1;}
    else{VAR_TUNE[3]=VAR_TUNE[3]*0.9;}
    // exchange --- done till here

    p1 = mypdf(y_now[0],d_size);
    p2 = pow(mypdf(y_now[1],d_size),0.5625);
    p1_2 = mypdf(y_now[1],d_size);
    p2_2 = pow(mypdf(y_now[0],d_size),0.5625);

    TOSS_A_DICE = (rand() % 10000) / 10000.0;
    if(TOSS_A_DICE < (p2_2*p1_2)/(p1*p2) ){
        for(int i =0;i<d_size;i++){
            exchange[i] = y_now[0][i];
            y_now[1][i] = y_now[0][i];
            y_now[0][i] = exchange[i];
        }
    }
    p3 = pow(mypdf(y_now[2],d_size),0.25);
    p4 = pow(mypdf(y_now[3],d_size),0.0625);
    p3_2 = pow(mypdf(y_now[3],d_size),0.25);
    p4_2 = pow(mypdf(y_now[2],d_size),0.625);
    
    TOSS_A_DICE = (rand() % 10000) / 10000.0;
    if(TOSS_A_DICE < (p4_2*p3_2)/(p3*p4) ){
        for(int i =0;i<d_size;i++){
            exchange[i] = y_now[0][i];
            y_now[1][i] = y_now[0][i];
            y_now[0][i] = exchange[i];
        }
    }
    
    p3 = pow(mypdf(y_now[2],d_size),0.25);
    p1 = mypdf(y_now[0],d_size);
    p3_2 = pow(mypdf(y_now[0],d_size),0.25);
    p1_2 = mypdf(y_now[2],d_size);
    
    TOSS_A_DICE = (rand() % 10000) / 10000.0;
    if(TOSS_A_DICE < (p1_2*p3_2)/(p1*p3)){
        for(int i =0;i<d_size;i++){
            exchange[i] = y_now[0][i];
            y_now[1][i] = y_now[0][i];
            y_now[0][i] = exchange[i];
        }
		exchanged++;
    }
    
    EX_MCMC=EX_MCMC+expectancy(y_now[0], d_size);
    criteria++;
}
float t2 = clock();
exe_time = (t2-t1)/double(CLOCKS_PER_SEC);

cout << "Accepted Chain1: " <<   accepted1 << endl;
cout << "Rejected Chain1: " <<   rejected1 << endl;
cout << "Accepted Chain2: " <<   accepted2 << endl;
cout << "Rejected Chain2: " <<   rejected2 << endl;
cout << "Accepted Chain1: " <<   accepted3 << endl;
cout << "Rejected Chain1: " <<   rejected3 << endl;
cout << "Accepted Chain2: " <<   accepted4 << endl;
cout << "Rejected Chain2: " <<   rejected4 << endl;
cout << "Exchanged b/w 2 Chains: " <<   exchanged << endl;
cout << "Converged Value: " <<  EX_MCMC/i << endl;
cout <<"Execution time: " <<  exe_time << endl;

}

double sampleNormal() {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return sampleNormal();
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

double mypdf (double x[], int size) 
{ 
    double result = 1;
    int i=0;
    for(int i=0; i<size;i++){
        result = result*((exp(-(pow(((x[i]-1) /1),2))/2))/1/sqrt(2*pi));
    }
    return result;
}

double expectancy (double x[], int size) 
{
    double result = 1;
    int i=0;
    for(int i=0; i<size;i++){
        result = result*x[i];
    }
    return result;
}