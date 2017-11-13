/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		1th May 2009
	
*/
#ifndef MAKETARGET_H_
#define MAKETARGET_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mylib.h"

#define MAXSSTATS 15		//number of possible used summary statistics
#define MAXSSTATS_M 6		//number of possible used summary statistics for Microssatelites
#define MAXSSTATS_S 9		//number of possible used summary statistics for Sequence Data

// struct used to save the information about a gene tree
struct data{
	int npop;				//number of populations
	int nloc;				//number of loci
	int **nsamp;			//pointer to a 2d-array(population, number of samples)
	int *tsamp;				//total number of samples per loci
	int ***freq;			//frequency of the diferent haplotypes by loci by pop
	int *ldna;				//array that gives the number of all the diferent haplotypes by loci
	int *Nmax;				//maximum number of the diferent haplotypes
	//only used in microssatelites analysis
	int **valM;			//all the diferent haplotypes by loci
	//only used in dna sequence analysis
	int *lsites;			//array(number of sites of the dna sequency (number of mutations), loci)
	char ***valS;		//2d-array(all the diferent haplotypes, loci)
};

/*
	Prints the summary statistics to a char array (Microsatellites)
	Defined in summStats.c
	
	@param data - struture with the informations about a genealogical tree
	@param lsstats - list of the used sstats (0-absent;1-present)
	@param outp - string that will store the summary statistics
	@param foundSTR - check if STR's are present in the study
	@param foundSNP - check if SNP's are present in the study
	@param ltype - DNA type per loci
*/
void summStats(struct data *data,int *lsstats,char *outp,int foundSTR,int foundSNP,char *ltype);

#endif /*MAKETARGET_H_*/
