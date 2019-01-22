/*
 Gillespie Algorithm (Direct Method) for HIV simulation in C
 ref: http://pubs.acs.org/doi/abs/10.1021/j100540a008 

	Solves the model:

	0 ->(aS) 	S
	S ->(dS) 	0
	S ->(Bt*V) 	I
	I ->(dI) 	p*V
	V ->(g) 	0

 How to compile:
 > gcc HIV_gillespie.c -lm -o stochastic_hiv
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

//change simulation parameters
#define ENDTIME   10.0		//end of time [days]
#define TIMESTEP  0.05		//interval of output
#define N         6			//number of reaction equations
#define M         3			//number of state variables
int x[M];				//population of chemical species
double c[N];			//transition rates
double prop[N];			//propencity vectors
int T[N][M];			//transition matrix for updating x[]
	
void init(int x[], double c[], int T[5][3]){

	//initial conditions
	x[0] = (int) c[0]/c[1]; //S_0
//	x[1] = 1;  				//I_0
//	x[2] = 0;//(int) c[4];  	//V_0
	x[1] = 0;  				//I_0
	x[2] = 1;  	//V_0
				 
	// transition matrix T[i][j]: i=transition number, j=index of state that transitions
	T[0][0] =  1; T[0][1] =  0; T[0][2] =  0;
	T[1][0] = -1; T[1][1] =  0; T[1][2] =  0;
	T[2][0] = -1; T[2][1] =  1; T[2][2] = -1;
	T[3][0] =  0; T[3][1] = -1; T[3][2] =  (int) c[4];
	T[4][0] =  0; T[4][1] =  0; T[4][2] = -1;

}

//function that chooses which transition happens based on propensity vector
int select_reaction(double prop[], int pn, double sum_propencity, double r){
	int reaction = -1; 	//initialize the reaction selection integer
	double sp = 0.0;	//keep track of each window in terms of reaction probability
	r = r * sum_propencity; //calculate the total length of windows
	//loop through each process, if sp in that reaction window choose it
	int i;
	for(i=0; i<pn; i++){
		sp += prop[i];
		if(r < sp){
			reaction = i;
			break;
		}
	}
	return reaction;
}

//function that updates the propensity vector given state vector
void update_p(double prop[], double c[], int x[]){
	prop[0] = c[0];	    	//constant production 
	prop[1] = c[1]*x[0];		//density dependent susceptible death
	prop[2] = c[2]*x[0]*x[2];	//infection
	prop[3] = c[3]*x[1];		//infected cell burst
	prop[4] = c[5]*x[2];		//density dependent viral clearance
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
	double sum_propencity=0.0;	//intialize sum of propencities
	double dt=0.0;				//step of time
	double t=0.0;				//time
	double r;					//random number
	int reaction;				//reaction number selected

	//variables for model
	float aS=1e8, dS=0.5, Bt=1e-10, dI=1.0, p=10000.0, g=23.0;

	float R0 = aS/dS*Bt/dI*p/g;
	
	printf("%s","R_0 = ");
	printf("%9.6f",R0);
	//printf("%s","\n");


    FILE *out;
	
	//int a;
	//a = strtol(argv[1], NULL, 0);
	
    if (argc >= 2){
    	out = fopen(argv[1], "w");
	}
	else { 
		out = fopen("out.txt", "w");
	}
	
   
	//intializes random number generator
	time_t tt;
	srand((unsigned) time(&tt));		
		
	//transition rates
	c[0]=aS; c[1]=dS; c[2]=Bt; c[3]=dI; c[4]=p; c[5]=g;
	//bnAb decay rate
	double a1  = 0.01;
	float bnAb = 1.0;

	init(x, c, T);

	//main loop over time
	while(t < ENDTIME){
	
		bnAb = 0;//0.5*exp(-a1*t); //proportion of infections prevented 
			
		c[2] = c[2] * (1-bnAb);
		
		
		//print output
		output(out, t, x, M);
		
		//printf("%9.6f",c[2]);
	
		//update propencity vector
		update_p(prop, c, x);
		sum_propencity = sum(prop, N);
	
		//get a sample dt, just to be sure check >0 probability
		if(sum_propencity > 0){
			dt = -log((double)rand()/INT_MAX) / sum_propencity;
		}
		else{
			break;
		}
	
		// select transition
		r = (double) rand()/INT_MAX;
		reaction = select_reaction(prop, N, sum_propencity, r);
	
		// update chemical species
		update_x(x, T, reaction);	
	
		// time
		t += dt;
	
	}
	fclose(out);
	

	return(0);
}


