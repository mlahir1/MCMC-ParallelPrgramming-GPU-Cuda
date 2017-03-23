#include <iostream>
#include <stdlib.h>     //srand, rand 
#include <time.h>       // time 
#include "math.h"       // for RAND, and rand, log and sqrt
#include <ctime>
using namespace std;    // for cout

double sampleNormal();  // function declaration
double mypdf     (double x1);
double expectancy(double x1);

const double pi = 3.1415926535897;
const double TUNED_VAR = 2.000;

int main()
{
float t1=clock();
double TOSS_A_DICE;

double accepted=0;
double rejected=0;
double EX_MCMC = 0;
double y1_c=1;

// initialization
double y1_now=1;

double i=1;
double criteria=0;


//initialize random seed
srand (time(NULL));


while (i<100000) 
{
	y1_c = y1_now + TUNED_VAR*sampleNormal();

	TOSS_A_DICE = (rand() % 10000) / 10000.0; // random between 0 and 1
	if (TOSS_A_DICE < mypdf(y1_c) /mypdf(y1_now) ) 
	{
		y1_now = y1_c;

		
		accepted=accepted+1;

	} 
	else 
	{  
		rejected=rejected+1;
	}

	EX_MCMC=EX_MCMC+expectancy(y1_now);
		
	//if (((EX_MCMC/i)<1.05) && ((EX_MCMC/i)>0.95)) {criteria=criteria+1;} else {criteria=0;}
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

double mypdf (double x1) 
{ 
double pdf;

double f1 = (exp(-(pow(((x1-1) /1),2))/2))/1/sqrt(2*pi);


pdf = f1;

return pdf;   
}

double expectancy (double x1 ) 
{ 
double expected;

expected = x1;

return expected;   
}
