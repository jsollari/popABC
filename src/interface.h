/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		12th May 2009
*/
#ifndef INTERFACE_H_
#define INTERFACE_H_

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
#define NPOP 2
#define LOCMAXCHAR 10
#define MAXMODEL 5
#define NPLOIDY 2
#define MAXCHAR 256     //max number of char of the pops names
#define MAXNPRIOR 10000		//maximum number of parameter values to be stored
 				 
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
	Start of the ABC program. Simulates coordenated points (summary_statistics, parameters) to be
	used in an Aproximate Bayesian Computation model to estimate the true values of a genealogical tree
	(sequence data).
	
	@param input_prs - .prs input filename 
	@param input_ssz - .ssz input filename
	@param input_sst - .sst input filename
	@param output    - .dat output filename
	@param printIt   - print or not the .len and .frq files of every genetic tree (0-don't print; 1-print)
	@param printMut  - print or not the .mut (0-don't print; 1-print)
	@param printRec  - print or not the .rec (0-don't print; 1-print)
*/
int abc(char *input_prs,char *input_ssz,char *input_sst, char* output,int printIt,int printMut,int printRec);

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
int sampPriors(struct params *pm,char *outline1,FILE *out_mut,FILE *out_rec,int printMut,int printRec,int foundSTR,int foundSNP,char *ltype);

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
int sampLikelihood(struct params *pm, struct data *data,int printit,long citer, char *path);

/*
	Frees the memory allocated in the first iteration
	Defined in abc.c
	
	@param nloc number of loci
	@param npop number of populations
	@param nevt number of events 	
*/
void freetree(int nloc,int npop,int nevt);

/*
	This function gets an sample file (.pop) and creates a table file (.len)
	
	@param input  - input filename
	@param output - output filename
*/	 
int createFreqTab(char *intput,char *output);

/*
	This function gets an IMa input file and creates a table file (.len)
	
	@param input  - input filename
	@param output - output filename
*/	 
int createFreqTab2(char *intput,char *output);

/*
	This function gets a GenePop input file and creates a table file (.len)
	
	@param input  - input filename
	@param output - output filename
*/	 
int createFreqTab3(char *intput,char *output);

/*
	This function gets a Nexus file (.nex) and creates a table file (.len)
	
	@param input  - input filename
	@param output - output filename
*/	 
int createFreqTab4(char *intput,char *output);

/*
	This funtion uses a freq_tab_length file and creates a output file which contains
	the summary statistics of the given file and its sample sizes
	
	@param input filename (*.len)
	@param input filename (*.sst)
	@param output filename (*.trg AND *.szz)
	@return an integer to check for errors
*/
int maketarget(char *input_len,char *input_sst,char *output);

/*
	This program prints out a proportion, specified by tolerance, of the simulated data with an ABC approach
	which are closer in Euclidian distance to the target summary statistics. It will also create a file
	which will contain the simulated populational tree parameters (only) in the first 10000 lines of the data.
	These can then be used to build a posterior distributions with which we can compare prior distributions

	@arg filename with simulated data
	@arg filename with target summary statistics (summary statistics from our "real" data)
	@arg filename of the output of the program
	@arg number of parameters
	@arg number of summstats
	@arg tolerance of analysis (between 0 and 1)
*/
int firstpass(char *input,char *target,char *output,int nparam,int nsstas,double tol);

/*
	This file creates a .prs file.
	
	@param outname	- name of the output file
	@param niter	- number of simulations to run
	@param genet	- generation time of the populations species
	@param npop		- number of population
	@param nloc		- number of loci
	@param lplo		- herederiter scalar per locus
	@param ltype	- type of DNA data per locus
	@param pr_top	- prior struct - topology
	@param pr_Ne	- list of prior struct (size=npop*2-1) - pop size
	@param pr_tev	- list of prior struct (size=npop-1) - tev
	@param pr_mig	- list of prior struct (size=2*npop-2) - mig
	@param pr_mut	- prior struct - mutSTR
	@param pr_mut	- prior struct - mutSNP
	@param pr_rec	- prior struct - recSTR
	@param pr_rec	- prior struct - recSNP
*/
int makeprior(char *output,
			  int niter,int genet,int npop,int nloc,		//1sr line
			  double *lplo,									//2nd line
			  char *ltype,									//3th line
			  struct prior pr_top,							//4th line
			  struct prior *pr_Ne,							//5th line
			  struct prior *pr_tev,							//6th line
			  struct prior *pr_mig,							//7th line
			  struct prior pr_mutSTR,						//8th line
			  struct prior pr_mutSNP,						//9th line
			  struct prior pr_recSTR,						//10th line
			  struct prior pr_recSNP,						//11th line
			  struct migweights migw);						//12th line

/*
	This file creates a .sst file.
	
	@param outname	- name of the output file
	@param lsstats	- list of the used sstats (0 - don't use; 1 - use)
*/
int makestats(char *output,int *lsstats);

/*
	This funtion takes a .len file and creates a sample populations for it
	@param filename of .len file
	@param filename of the output of the function
	@return an integer to check for errors
*/
int makepop(char *input, char *output);

/*
	It creates an ABC .len file given an IM input file.
	
	@param input filename
	@param output filename
	@return a integer to check for errors
*/
int convertToABC(char *input,char *output);

/*
	This function joins files together

	@param ninp number of files to join
	@param linp list of the input filenames to join
	@param out output filename
	@return an integer to check for errors
*/
int joindata(int ninp, char *linp[],char *outp);

static const int RecIter = 50000,	//stepsize of iteration when the program records information
 				 RecLine = 500000,  //record step per lines analysed
				 DnaSizeM = 100,	//size of microssatellite dna
 				 DnaSizeS = 10;		//size of dna sequence
int Exepathsize;						//size of the path to the executable
char *Exepath;							//path to the executable
struct prior P_top,					//prior of topology
 			 P_mutM,				//hiperprior of mutation rate Microsatelites				
 			 P_mutS,				//hiperprior of mutation rate Sequence data
 			 P_recM,				//hiperprior of recombination rate Microsatelites				
 			 P_recS,				//hiperprior of recombination rate Sequence data
 			 *P_t,					//priors of time events
			 **P_mig,				//prior of migration
			 **P_psize;				//priors of Ne by events				 
struct migweights M_migw;			//migration weights matrix

#endif /*INTERFACE_H_*/
