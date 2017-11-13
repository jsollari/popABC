/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		1st May 2009
	
*/
#ifndef ABC_H_
#define ABC_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mylib.h"

#define MAXSSTATS 15		//number of possible used summary statistics
#define MAXSSTATS_M 6		//number of possible used summary statistics for Microssatelites
#define MAXSSTATS_S 9		//number of possible used summary statistics for Sequence Data
#define COALEVNT 3			//different types of events in the Genetic Trees
#define MAXDATA 250000000	//maximum number of data (nparam+nstats)*nsim to be analysed
 				 
// struct used to save the parameters of a prior distribution 
struct prior{
	int type;				//Distribution of the prior
	double *p;				//parameters of the distribution
};

// struct used to save the migration weights matrix 
struct migweights{
	int type;				//Distribution of the prior
	double ***m;			//migweights per npop per ntev per npop
};

// struct used to save the information about a gene tree
struct data{
	int npop;				//number of populations
	int nloc;				//number of loci
	int **nsamp;			//pointer to a 2d-array(population, number of samples)
	int *tsamp;				//total number of samples per loci
	int ***freq;			//frequency of the diferent haplotypes by loci by pop
	int *ldna;				//array that gives the number of all the diferent haplotypes by loci
	int *Nmax;				//maximum number of different alleles/haplotypes
	//only used in microssatelites analysis
	int **valM;				//all the diferent haplotypes by loci
	//only used in dna sequence analysis
	int *lsites;			//array(number of sites of the dna sequency (number of mutations), loci)
	char ***valS;			//2d-array(all the diferent haplotypes, loci)
};

// struct used to save the parameters to build a gene tree 
struct params{
	int npop;				//number of populations
	int niter;				//number of iterations
	int nevt;				//number of events
	int nrSTR;				//number of linked STR's loci
	double gent;			//number of generation time
	int ntop;				//number of possible topologies
	double *ploidy;			//heriditage scalar per locus
	char *type;				//DNA type
	int nloc;				//number of loci
	int **nsamp;			//pointer to 2d array: population vs number of samples
	double *mu;				//pointer to array with mutation rate in diferente loci
	double *rec;			//pointer to array with recombination rate in different loci
	int topol;				//topology
	int *seq;				//pointer to array with information about split population events
	double *tev;			//pointer to array with events' times
	double **psize;			//2d array: populations' Ne size vs number of events + 1 
	double **mig;			//2d array: migration rates vs number of events
};

/*
	This function samples values from the priors and stores them
	Defined in samplePriors.c
		
	@param pm 		- pointer to the parameters to be used in the simulation
  	@param outline1 - place to store all the first part of the data
 	@param out_mut  - pointer to output file .mut
  	@param printMut - indicate if the mutation rates are going to be printed or not
	@param foundSTR - check if STR's are present in the study
	@param foundSNP - check if SNP's are present in the study
	@param ltype    - DNA type per loci
*/
void sampPriors(struct params *pm,char *outline1,FILE *out_mut,FILE *out_rec,int printMut,int printRec,int foundSTR,int foundSNP,char *ltype);

/*
	This function calulates the summary statistics from simulated data and stores them
	Defined in summStats.c
	
	@param data - struture with the informations about a genealogical tree
	@param lsstats - list of the used sstats (0-absent;1-present)
	@param outp - string that will store the summary statistics
	@param foundSTR - check if STR's are present in the study
	@param foundSNP - check if SNP's are present in the study
	@param ltype - DNA type per loci
*/
void summStats(struct data *data,int *lsstats,char *outp,int foundSTR,int foundSNP, char*ltype);

/*
	This function simulats the genetic data
	Defined in geneticTree.c
	
	@param pm 	   - parameters used to build a geneology tree
	@param data    - data from a geneology tree
	@param printit - use to print or not the informations to a separate file
	@param citer   - current iteration
	@param path	   - path where the output files are going to be store in
*/
void sampLikelihood(struct params *pm, struct data *data,int printit,long citer, char *path);

/*
	Frees the memory allocated in the first iteration
	Defined in abc.c
	
	@param nloc number of loci
	@param npop number of populations
	@param nevt number of events 	
*/
void freetree(int nloc,int npop,int nevt);

static const int RecIter = 50000,	//stepsize of iteration when the program records information
 				 DnaSizeM = 100,	//size of microssatellite dna
 				 DnaSizeS = 10;		//size of dna sequence
struct prior P_top,					//prior of topology
 			 P_mutM,				//hiperprior of mutation rate Microsatelites				
 			 P_mutS,				//hiperprior of mutation rate Sequence data
 			 P_recM,				//hiperprior of recombination rate Microsatelites				
 			 P_recS,				//hiperprior of recombination rate Sequence data
 			 *P_t,					//priors of time events
			 **P_mig,				//prior of migration
			 **P_psize;				//priors of Ne by events				 
struct migweights M_migw;			//migration weights matrix

#endif /*ABC_H_*/
