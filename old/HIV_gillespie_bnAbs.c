/*
 Gillespie Algorithm (Direct Method) for HIV simulation in C
 ref: http://pubs.acs.org/doi/abs/10.1021/j100540a008 

	Solves the model:

	0 -> S		(aS)
	S -> 0		(dS)
	S -> I		(Bt(t)*V)
	I -> p*V	(dI)
	V -> 0		(g)

 How to compile:
 > gcc HIV_gillespie_bnAbs.c -lm -o stochastic_hiv
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

//change simulation parameters
#define ENDTIME   6*7.0		//end of time [days]
#define TIMESTEP  0.05		//interval of output
#define N         5			//number of reaction equations
#define M         3			//number of state variables
int x[M];				//population of chemical species
double c[N];			//transition rates
double p[N];			//propencity vectors
int T[N][M];			//transition matrix for updating x[]
	
void init(int x[], double c[], int T[5][3]){

	//variables for model
	float aS = 50;		//birth rate susceptible cells [cells/uL/day]
	float dS = 0.5;		//susceptible cell death rate [1/day]
	float Bt = 0.0001;	//viral infectivity [cells/virus/day]
	float dI = 1.0;		//death rate of infected cells [1/day]
	int p  = 1000; 		//viral production rate [virus/cell/day] must be integer
	float g  = 23.0;	//viral clearance rate [1/day]

	float R0 = aS/dS*Bt/dI*10000.0/23.0;
	printf("%s","R_0 = ");
	printf("%9.6f",R0);
	printf("%s","\n");
		
	float V0 = 0.03; //initial concentration of virus
		
	//initial conditions
	x[0] = aS/dS; 				//S_0 [cells per uL]
	x[1] = 23.0/10000.0*V0;  	//I_0 [cells per uL]
	x[2] = V0;  				//V_0 [copies per uL]
				 
	//transition rates
	c[0] = aS;
	c[1] = dS;
	c[2] = Bt;
	c[3] = dI;
	c[4] = g;
			
	// transition matrix T[i][j]: i=transition number, j=index of state that transitions
	T[0][0] =  1;
	T[0][1] =  0;
	T[0][2] =  0;
	T[1][0] = -1;
	T[1][1] =  0;
	T[1][2] =  0;
	T[2][0] = -1;
	T[2][1] =  1;
	T[2][2] = -1;
	T[3][0] =  0;
	T[3][1] = -1;
	T[3][2] =  p;
	T[4][0] =  0;
	T[4][1] =  0;
	T[4][2] = -1;

}

//function that chooses which transition happens based on propensity vector
int select_reaction(double p[], int pn, double sum_propencity, double r){
	int reaction = -1; 	//initialize the reaction selection integer
	double sp = 0.0;	//keep track of each window in terms of reaction probability
	r = r * sum_propencity; //calculate the total length of windows
	//loop through each process, if sp in that reaction window choose it
	int i;
	for(i=0; i<pn; i++){
		sp += p[i];
		if(r < sp){
			reaction = i;
			break;
		}
	}
	return reaction;
}

//function that updates the propensity vector given state vector
void update_p(double p[], double c[], int x[]){
	p[0] = c[0];	    	//constant production 
	p[1] = c[1]*x[0];		//density dependent susceptible death
	p[2] = c[2]*x[0]*x[2];	//infection
	p[3] = c[3]*x[1];		//viral production from density dependent infected cell death
	p[4] = c[4]*x[2];		//density dependent viral clearance
}

//function that updates the state variables based on which reaction occured 
//(selected by select_reaction)
void update_x(int x[], int T[5][3], int reaction){
	int i;
	for(i=0; i<3; i++){
		x[i] += T[reaction][i];
	}
}

//print the output to a file
void output(FILE *out, double t, int x[], int xn){
	static double output_t = 0.0;
	int i;
	if(output_t <= t){
		fprintf(out, "%f", t);
		for(i=0; i<xn; i++){
			fprintf(out, "\t%d", x[i]); 
		}
		fprintf(out, "\n");
		output_t += TIMESTEP;
	}
}

//helper function that sums all the propensities to normalize
double sum(double a[], int n){
	int i;
	double s=0.0;
	for(i=0; i<n; i++) 
		s += a[i];
	return(s);
}

//int main(void){
int main(int argc, char **argv){

	// initialization
	double sum_propencity = 0.0;	//intialize sum of propencities
	double tau=0.0;					//step of time
	double t=0.0;					//time
	double r;						//random number
	int reaction;					//reaction number selected

	init(x, c, T);

	//printf("argv[2]"); 
	
	//int SEED;
    FILE *out;
	
    if (argc >= 2){
    	out = fopen(argv[1], "w");
		//SEED = argv[2];
		}
	else { 
		out = fopen("out.txt", "w");
		//SEED=123;
		}
	
		
	time_t tt;
   
	/* Intializes random number generator */
	srand((unsigned) time(&tt));
	//srand(SEED);
			
	//main loop over time
	while(t < ENDTIME){
	
		//print output
		output(out, t, x, M);
	
		//update propencity vector
		update_p(p, c, x);
		sum_propencity = sum(p, N);
	
		//get a sample tau, just to be sure check >0 probability
		if(sum_propencity > 0){
			tau = -log((double)rand()/INT_MAX) / sum_propencity;
		}else{
			break;
		}
	
		// select transition
		r = (double)rand()/INT_MAX;
		reaction = select_reaction(p, N, sum_propencity, r);
	
		// update chemical species
		update_x(x, T, reaction);	
	
		// time
		t += tau;
	
	}
	fclose(out);
	

	return(0);
}


