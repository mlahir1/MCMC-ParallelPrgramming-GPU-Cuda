#include <iostream>
#include <stdlib.h>     //srand, rand 
#include <time.h>       // time 
#include "math.h"       // for RAND, and rand, log and sqrt
#include <ctime>
using namespace std;    // for cout

double sampleNormal();  // function declaration
double mypdf     (double x1,double x2);
double expectancy(double x1,double x2);

const double pi = 3.1415926535897;
const double TUNED_VAR = 1.129;

int main()
{

double ITERATION=10;
double AVERAGE_I=0;
float TOTAL_TIME=0;

for (int z=0; z<ITERATION; ++z) 
{
	

	
float exe_time=0;
float t1=clock();
double TOSS_A_DICE;

double accepted=0;
double rejected=0;
double accepted2=0;
double rejected2=0;
double exchanged=0;

double EX_MCMC = 0;
double y1_c=1;
double y2_c=1;

double y1_c2=1;
double y2_c2=1;

// initialization
double y1_now=1;
double y2_now=1;

double y1_now2=1;
double y2_now2=1;

double exchange1=0;
double exchange2=0;
double exchange3=0;

double i=1;
double criteria=0;


//initialize random seed
srand (time(NULL));


while (criteria<100000) 
{
// chain 1 and the main 1	
y1_c = y1_now + TUNED_VAR*sampleNormal();
y2_c = y2_now + TUNED_VAR*sampleNormal();

TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
if (TOSS_A_DICE < mypdf(y1_c,y2_c) /mypdf(y1_now,y2_now) ) 
{
	y1_now = y1_c;
    y2_now = y2_c;

    accepted=accepted+1;

} 
else 
{  
    rejected=rejected+1;
}



//chain 2 and use of pow(base,exponent)
y1_c2 = y1_now2 + TUNED_VAR*sampleNormal();
y2_c2 = y2_now2 + TUNED_VAR*sampleNormal();

TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
if (TOSS_A_DICE < pow(mypdf(y1_c2,y2_c2),0.25) /pow(mypdf(y1_now2,y2_now2),0.25) )
{
	y1_now2 = y1_c2;
    y2_now2 = y2_c2;
 
    accepted2=accepted2+1;

} 
else 
{  
    rejected2=rejected2+1;
}

// exchange
TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
   if (TOSS_A_DICE < ((mypdf(y1_now2,y2_now2)*(pow(mypdf(y1_now,y2_now),0.25))) / (mypdf(y1_now,y2_now)*( pow(mypdf(y1_now2,y2_now2),0.25)))))
       {
	   exchange1 = y1_now2;
	   exchange2 = y2_now2;

       y1_now2  = y1_now;
       y2_now2  = y2_now;

       y1_now = exchange1;
       y2_now = exchange2;
 	   
	   exchanged = exchanged + 1;
	   }



EX_MCMC=EX_MCMC+expectancy(y1_now,y2_now);
    
if (((EX_MCMC/i)<1.05) && ((EX_MCMC/i)>0.95)) {criteria=criteria+1;} else {criteria=0;}

i=i+1;


}
float t2 = clock();
exe_time = (t2-t1)/double(CLOCKS_PER_SEC);

cout << "Accepted Chain1: " <<   accepted << endl;
cout << "Rejected Chain1: " <<   rejected << endl;
cout << "Accepted Chain2: " <<   accepted2 << endl;
cout << "Rejected Chain2: " <<   rejected2 << endl;
cout << "Exchanged b/w 2 Chains: " <<   exchanged << endl;
cout << "Converged Value: " <<  EX_MCMC/i << endl;
cout<<"Execution time: "    <<  exe_time << endl;

AVERAGE_I = AVERAGE_I + i;
TOTAL_TIME = TOTAL_TIME + exe_time;


}

cout<<"Average Number of Iterations: "    <<  AVERAGE_I/ITERATION << endl;
cout<<"Average Time of Execution: "    <<  TOTAL_TIME/ITERATION << endl;
return 0;

}


double sampleNormal() {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return sampleNormal();
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

double mypdf (double x1, double x2) 
{ 
double pdf;

double f1 = (exp(-(pow(((x1-1) /1),2))/2))/1/sqrt(2*pi);
double f2 = (exp(-(pow(((x2-1) /1),2))/2))/1/sqrt(2*pi);


pdf = f1*f2;

return pdf;   
}

double expectancy (double x1, double x2) 
{ 
double expected;

expected = x1*x2;

return expected;   
}
