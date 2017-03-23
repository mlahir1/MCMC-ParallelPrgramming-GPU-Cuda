#include <iostream>
#include <stdlib.h>     //srand, rand 
#include <time.h>       // time 
#include "math.h"       // for RAND, and rand, log and sqrt
#include <ctime>
using namespace std;    // for cout

double sampleNormal();  // function declaration
double mypdf     (double x1,double x2,double x3,double x4,double x5,double x6,double x7);
double expectancy(double x1,double x2,double x3,double x4,double x5,double x6,double x7);

const double pi = 3.1415926535897;
const double TUNED_VAR = 0.510;

int main()
{
float t1=clock();
double TOSS_A_DICE;

double accepted=0;
double rejected=0;
double EX_MCMC = 0;
double y1_c=1;
double y2_c=1;
double y3_c=1;
double y4_c=1;
double y5_c=1;
double y6_c=1;
double y7_c=1;

// initialization
double y1_now=1;
double y2_now=1;
double y3_now=1;
double y4_now=1;
double y5_now=1;
double y6_now=1;
double y7_now=1;

double i=1;
double criteria=0;


//initialize random seed
srand (time(NULL));


while (criteria<100000) 
{
y1_c = y1_now + TUNED_VAR*sampleNormal();
y2_c = y2_now + TUNED_VAR*sampleNormal();
y3_c = y3_now + TUNED_VAR*sampleNormal();
y4_c = y4_now + TUNED_VAR*sampleNormal();
y5_c = y5_now + TUNED_VAR*sampleNormal();
y6_c = y6_now + TUNED_VAR*sampleNormal();
y7_c = y7_now + TUNED_VAR*sampleNormal();

TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
if (TOSS_A_DICE < mypdf(y1_c,y2_c,y3_c,y4_c,y5_c,y6_c,y7_c) /mypdf(y1_now,y2_now,y3_now,y4_now,y5_now,y6_now,y7_now) ) 
{
	y1_now = y1_c;
    y2_now = y2_c;
    y3_now = y3_c;
    y4_now = y4_c;
    y5_now = y5_c;
    y6_now = y6_c;
    y7_now = y7_c;
   
    accepted=accepted+1;

} 
else 
{  
    rejected=rejected+1;
}

EX_MCMC=EX_MCMC+expectancy(y1_now,y2_now,y3_now,y4_now,y5_now,y6_now,y7_now);
    
if (((EX_MCMC/i)<1.05) && ((EX_MCMC/i)>0.95)) {criteria=criteria+1;} else {criteria=0;}

i=i+1;


}

cout <<   accepted << endl;
cout <<   rejected << endl;
cout <<   EX_MCMC/i << endl;

float t2 = clock();
cout<<"Exe time: " << (t2-t1)/double(CLOCKS_PER_SEC) << endl;

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

double mypdf (double x1, double x2, double x3, double x4, double x5, double x6, double x7 ) 
{ 
double pdf;

double f1 = (exp(-(pow(((x1-1) /1),2))/2))/1/sqrt(2*pi);
double f2 = (exp(-(pow(((x2-1) /1),2))/2))/1/sqrt(2*pi);
double f3 = (exp(-(pow(((x3-1) /1),2))/2))/1/sqrt(2*pi);
double f4 = (exp(-(pow(((x4-1) /1),2))/2))/1/sqrt(2*pi);
double f5 = (exp(-(pow(((x5-1) /1),2))/2))/1/sqrt(2*pi);
double f6 = (exp(-(pow(((x6-1) /1),2))/2))/1/sqrt(2*pi);
double f7 = (exp(-(pow(((x7-1) /1),2))/2))/1/sqrt(2*pi);

pdf = f1*f2*f3*f4*f5*f6*f7;

return pdf;   
}

double expectancy (double x1, double x2, double x3, double x4, double x5, double x6, double x7) 
{ 
double expected;

expected = x1*x2*x3*x4*x5*x6*x7;

return expected;   
}
