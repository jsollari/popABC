/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		13th November 2008
*/
#include "interface.h"

/*
	This function calculates the allele number in one population

	@param nhap  - total number of alleles
	@param freq  - frequency of an allele
	@return the allele number in one population
*/
int calc_nhap(int nhap, int freq[]);

/*
	This function calculates the Shannon's index of one population

	@param nhap  - number of all the haplotypes
	@param npop  - number of population
	@param nsamp - number of samples of a population
	@param freq  - frequency of haplotypes of a population in a locus
	@param cpop  - current populations
	@result number of haplotypes of one population
*/
double calc_shanon(int nhap, int nsamp, int *freq, int cpop);

/* only used in Microsatellite analysis ********************************************************************/

/*
	This function calculates the hetetozygosity in a population

	@param sum   - size of the sample
	@param freq  - frequency of a microsatellite size (in number of occurrencies
	@param nhap  - number of diferent microsatellite sizes
	@return the hetetozygosity in one population
*/
double calc_het(int sum,int freq[], int nhap);

/*
	This function calculates the variance of allele length in a population

	@param sum   - size of the sample
	@param freq  - frequency of a microsatellite size (in number of occurrencies)
	@param val   - diferent sizes of a microsatellite
	@param nhap  - number of diferent microsatellite sizes
	@return variance of allele length in one population
*/
double calc_var(int sum,int freq[],int val[],int nhap);

/*
	This function calculates the curtosis of allele length in a population

	@param nhap  - total number of microsatellite sizes
	@param freq  - frequency of a microsatellite size (in number of occurrencies)
	@return the allele number in one population
*/
double calc_curt(int nhap, int *freq);


/*
	This function calculates the Nm estimator using the heterozygosity

	@param Hw - heterozygosity within pop
	@param Ha - heterozygosity among pop
	@result Nm_H
*/
double calc_Nm_H(double Hw,double Ha);

/* only used in Sequence Data analysis *********************************************************************/
/*
	This function is used to caluclate the pairwise difference for one population

	@param sum     - number of samples
	@param freq    - frequency of haplotypes of a population in a locus
	@param lsites  - size of a haplotype in a locus
	@param val     - all the haplotypes in a locus
	@param nhap    - number of all the haplotypes
	@result pairwise difference of one population
*/
double calc_pi(int sum,int freq[], int lsites, char **val,int nhap);

/*
	Number of diferent sites in two diferent haplotypes

	@param n  - number of sites in the haplotype
	@param vi - one haplotype
	@param v2 - other haplotype
	@result number of diferent sites in the two haplotypes
*/
int numdiff(int n, char *v1, char *v2);

/*
	This function is used to caluclate the number of segregate sites of one population

	@param sum     - number of samples
	@param freq    - frequency of haplotypes of a population in a locus
	@param lsites  - size of a haplotype in a locus
	@param val     - all the haplotypes in a locus
	@param nhap    - number of all the haplotypes
	@param lsnp	   - SNP in a locus in a population
	@result number of segregate sites of one population
*/
double calc_segsites(int sum,int freq[],int lsites, char **val,int nhap, int *lsnp);

/*
	This function is used to caluclate the number of segregate sites of populations pooled together

	@param sum     - number of samples
	@param freq    - frequency of haplotypes of a population in a locus
	@param lsites  - size of a haplotype in a locus
	@param val     - all the haplotypes in a locus
	@param nhap    - number of all the haplotypes
	@result number of segregate sites of one population
*/
double calc_segsites2(int sum,int freq[],int lsites, char **val,int nhap);

/*
	This function calculates the Mutation Frequency Spectrum (MFS)

	@param val    - array with all the haplotypes
	@param freq   - array with the frequency of all the haplotypes
	@param nhap   - number of different haplotypes
	@param nsites - number of sites of a locus
	@param cloc   - considered locus
	@param cpop   - considered population
	@param mfs    - Mutation Frequency Spectrum
*/
void fillMFS(char **val,int *freq,int nhap,int nsites,int nsamp,int cloc,int cpop,int***mfs);

/*
	This function changes a char '1' or '0' to a int 1 or 0

	@param c1 - char '1' or '0'
	@result int 1 or 0
*/
int segtoint(char c1);

/*
	This function calculates the mean of the Mutation Frequency Spectrum (MFS)

	@param nsites - number of sites of a locus
	@param mfs	  - MFS
	@param snp	  - number of SNP
	@result mean of MFS
*/
double calc_meanMFS(int nsites,int *mfs, int snp);

/*
	This function calculates the stdev of the Mutation Frequency Spectrum (MFS)

	@param nsites - number of sites of a locus
	@param mfs	  - MFS
	@param snp	  - number of SNP
	@param lsnp	  - list the SNPs
	@param mMFS	  - mean of MFS
	@result stdev of MFS
*/
double calc_stdevMFS(int nsites,int *mfs, int snp, int *lsnp, double mMFS);

/*
	This function calculates the Nm estimator using the segregating sites

	@param Sw - no segregating sites within pop
	@param Sa - no segregating sites among pop
	@result Nm_S
*/
double calc_Nm_S(double Sw,double Sa);

/*
	This function calculates the number of private segregating sites

	@param nsites   - number of sites of the current locus
	@param npop     - total number of populations
	@param pop      - current population
	@param loc		- current locus
	@param lsnp		- SNP per site per loci per pop
	@result privateS
*/
double calc_privS(int nsites, int npop, int pop, int loc,int*** lsnp);

/*
	This function calculates the average frequency of private segregating sites

	@param nsites   - number of sites of the current locus
	@param npop     - total number of populations
	@param nhap		- 
	@param pop      - current population
	@param loc      - current locus
	@param lsnp     - SNP per locus per pop per site
	@param freq     - array with the frequency of all the haplotypes
	@param val      - array with all the haplotypes
	@result S_1
*/
double calc_S_1(int nsites, int npop, int nhap, int pop,int loc,int ***lsnp,int **freq,char **val);

void summStats(struct data *data,int *lsstats,char *outp,int foundSTR,int foundSNP,char *ltype){
	int cloc,cdna,cpop,cpop2,csite,ic,	//iterators
		maxDna,							//maximum number of different haplotypes of all loci
		npop,							//number of population
		nloc,							//number of loci
		tpop,							//iterator for two populations
		tsamp,							//size of sample of two populations
		skip,							//number of loci with 1 or less samples
		npair,							//number of population pairs
		*tfreq,							//sum of the frequencies of diferent microsatellites size of two populations
		*ldna,							//number of diferent alleles sizes by loci
		**nsamp,						//size of sample per locus per population
		***freq;						//all the frequencies of the diferent microsatellite sizes by pop and loci
	/*only used in microssatellites data*/
	int	**valM;							//all the diferent microsatellite sizes by loci
	double sstats1_t,					//heterozygosity for all populations pooled together
		   *sstats1,					//heterozygosity difference in each population
	   	   *sstats2,					//variance of allele length in each population
		   *sstats3,					//allele no. in each population
	       *sstats4,					//curtosis of allele lengths in each population
	       *sstats5,					//Shanon's index in each population
		   *sstats6,					//Nm_H in each population
		   *sstats1_2,					//heterozygosity for each pair of population pooled together
	       *sstats2_2,					//variance of allele length for each pair of population pooled together
	       *sstats3_2,					//allele no.for each pair of population pooled togther
	       *sstats4_2,					//curtosis of allele lengths for each pair of population pooled together
	       *sstats5_2;					//Shanon's index for each pair of population pooled togther
	/*only used in sequence data*/
	int *lsites,						//size of the haplotypes by loci
		***mfs,							//mutation frequency spectrum (MFS) by loci by pop
		**snp,							//number of SNPs by loci by pop
		***lsnp;						//SNP per sites per loci per pop		
	double sstats8_t,					//segregate sites no. for all populations pooled together
		   *sstats7,					//pairwise difference in each population
	   	   *sstats8,					//segregate sites no. in each population
		   *sstats9,					//haplotypes no. in each population
		   *sstats10,					//Shanon's index in each population
	       *sstats11,					//mean of mutation frequency spectrum (MFS) in each population
	   	   *sstats12,					//stdev of mutation frequency spectrum (MFS) in each population
	   	   *sstats13,					//Nm_S in each population
	   	   *sstats14,					//private S in each population
	       *sstats15,					//S(1) in each population
	       *sstats7_2,					//pairwise difference for each pair of population pooled together
	       *sstats8_2,					//segregate sites no. for each pair of population pooled together
	       *sstats9_2,					//haplotypes no. for each pair of population pooled togther
	       *sstats10_2,					//Shanon's index for each pair of population pooled togther
		   **mMFS;						//mean of MFS by loci by pop
	char ***valS;						//all the diferent haplotypes by loci

	npop = data->npop;
	nloc = data->nloc;
	ldna = data->ldna;
	nsamp = data->nsamp;
	freq = data->freq;
	npair = combinations(npop,2);
	if(foundSTR)
		valM = data->valM;
	if(foundSNP){
		valS = data->valS;
		lsites = data->lsites;
	}

	/*allocate memory*/
	maxDna = ldna[0];
	for(cloc=1 ; cloc<nloc ; cloc++){
		if(maxDna<ldna[cloc])
			maxDna = ldna[cloc];
	}
	tfreq = (int *)malloc(maxDna*sizeof(int));
	if(foundSTR){
		sstats1 = (double *)malloc(npop*sizeof(double));
		sstats2 = (double *)malloc(npop*sizeof(double));
		sstats3 = (double *)malloc(npop*sizeof(double));
		sstats4 = (double *)malloc(npop*sizeof(double));
		sstats5 = (double *)malloc(npop*sizeof(double));
		sstats6 = (double *)malloc(npop*sizeof(double));
		sstats1_2 = (double *)malloc(npair*sizeof(double));
		sstats2_2 = (double *)malloc(npair*sizeof(double));
		sstats3_2 = (double *)malloc(npair*sizeof(double));
		sstats4_2 = (double *)malloc(npair*sizeof(double));
		sstats5_2 = (double *)malloc(npair*sizeof(double));
	}
	if(foundSNP){
		snp = (int **)malloc(npop*sizeof(int*));
		lsnp = (int ***)malloc(npop*sizeof(int**));
		mfs = (int ***)malloc(npop*sizeof(int**));
		mMFS = (double **)malloc(npop*sizeof(double*));
		for(cpop=0; cpop<npop; cpop++){
			snp[cpop] = (int *)malloc(nloc*sizeof(int));
			lsnp[cpop] = (int **)malloc(nloc*sizeof(int*));
			mMFS[cpop] = (double *)malloc(nloc*sizeof(double));
			mfs[cpop] = (int **)malloc(nloc*sizeof(int*));
			for(cloc=0; cloc<nloc; cloc++){
				mfs[cpop][cloc] = (int *)malloc(lsites[cloc]*sizeof(int));
				lsnp[cpop][cloc] = (int *)malloc(lsites[cloc]*sizeof(int));
			}
		}
		sstats7 = (double *)malloc(npop*sizeof(double));
		sstats8 = (double *)malloc(npop*sizeof(double));
		sstats9 = (double *)malloc(npop*sizeof(double));
		sstats10 = (double *)malloc(npop*sizeof(double));
		sstats11 = (double *)malloc(npop*sizeof(double));
		sstats12 = (double *)malloc(npop*sizeof(double));
		sstats13 = (double *)malloc(npop*sizeof(double));
		sstats14 = (double *)malloc(npop*sizeof(double));
		sstats15 = (double *)malloc(npop*sizeof(double));
		sstats7_2 = (double *)malloc(npair*sizeof(double));
		sstats8_2 = (double *)malloc(npair*sizeof(double));
		sstats9_2 = (double *)malloc(npair*sizeof(double));
		sstats10_2 = (double *)malloc(npair*sizeof(double));
	}
	/*start calculating the summary statistics*/
	memset(outp,'\0',sizeof(outp));
	if(foundSTR){
		/*hetetozygosity in each population */
		if(lsstats[0]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats1[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
						++skip;
						continue;
					}
					sstats1[cpop] += calc_het(nsamp[cloc][cpop],freq[cloc][cpop],ldna[cloc]);
				}
				if(nloc-skip==0)
					sstats1[cpop]=0;
				else
					sstats1[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats1[cpop]);
			}
		}
		/*variance of allele length in each population */
		if(lsstats[1]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats2[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
						++skip;
						continue;
					}
					sstats2[cpop] += calc_var(nsamp[cloc][cpop],freq[cloc][cpop],valM[cloc],ldna[cloc]);
				}
				if(nloc-skip==0)
					sstats2[cpop]=0;
				else
					sstats2[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats2[cpop]);
			}
		}
		/*allele number in each population */
		if(lsstats[2]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats3[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
						++skip;
						continue;
					}
					sstats3[cpop] += calc_nhap(ldna[cloc],freq[cloc][cpop]);
				}
				if(nloc-skip==0)
					sstats3[cpop]=0;
				else
					sstats3[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats3[cpop]);
			}
		}
		/*curtosis of allele length in each population */
		if(lsstats[3]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats4[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
						++skip;
						continue;
					}
					sstats4[cpop] += calc_curt(ldna[cloc],freq[cloc][cpop]);
				}
				if(nloc-skip==0)
					sstats4[cpop]=0;
				else
					sstats4[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats4[cpop]);
			}
		}
		/*Shanon's index for each population */
		if(lsstats[4]==1){
			for(cpop=0 ; cpop<npop ; ++cpop){
				sstats5[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
						++skip;
						continue;
					}
					sstats5[cpop] += calc_shanon(ldna[cloc],nsamp[cloc][cpop],freq[cloc][cpop],cpop);
				}
				if(nloc-skip==0)
					sstats5[cpop]=0;
				else
					sstats5[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats5[cpop]);
			}
		}
		/*heterozygosity for all populations pooled together */
		if(lsstats[5]==1 && npop>1){
			tsamp=0;
			for(cdna=0;cdna<maxDna;++cdna)
				tfreq[cdna]=0;
			for(cpop = 0;cpop <npop;++cpop)	{
				for(cloc=0;cloc<nloc;++cloc){
					tsamp += nsamp[cloc][cpop];
					for(cdna=0;cdna<ldna[cloc];++cdna)
						tfreq[cdna] += freq[cloc][cpop][cdna];
				}
			}
			sstats1_t = 0;
			skip = 0;
			for(cloc=0;cloc<nloc;++cloc){
				if(tsamp <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
					++skip;
					continue;
				}
				sstats1_t += calc_het(tsamp,tfreq,ldna[cloc]);
			}
			if(nloc-skip==0)
				sstats1_t=0;
			else
				sstats1_t /= nloc-skip;
		}
		/*Nm_H in each population */
		if(lsstats[5]==1 && npop>1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats6[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
						++skip;
						continue;
					}
					sstats6[cpop] += calc_Nm_H(sstats1[cpop],sstats1_t);
				}
				if(nloc-skip==0)
					sstats6[cpop]=0;
				else
					sstats6[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats6[cpop]);
			}
		}
		/*heterozygosity for each pair of population pooled together */
		if(lsstats[0]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats1_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's')	{
							++skip;
							continue;
						}
						tsamp = nsamp[cloc][cpop]+nsamp[cloc][cpop2];
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats1_2[ic] += calc_het(tsamp,tfreq,ldna[cloc]);
					}
					if(nloc-skip==0)
						sstats1_2[ic]=0;
					else
						sstats1_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats1_2[ic]);
					++ic;
				}
			}
		}
		/*variance of allele lengths for each pair of population pooled together */
		if(lsstats[1]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats2_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
							++skip;
							continue;
						}
						tsamp = nsamp[cloc][cpop]+nsamp[cloc][cpop2];
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats2_2[ic] += calc_var(tsamp,tfreq,valM[cloc],ldna[cloc]);
					}
					if(nloc-skip==0)
						sstats2_2[ic]=0;
					else
						sstats2_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats2_2[ic]);
					++ic;
				}
			}
		}
		/*number of alleles for each pair of population pooled together */
		if(lsstats[2]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats3_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] < 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
							++skip;
							continue;
						}
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats3_2[ic] += calc_nhap(ldna[cloc],tfreq);
					}
					if(nloc-skip==0)
						sstats3_2[ic]=0;
					else
						sstats3_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats3_2[ic]);
					++ic;
				}
			}
		}
		/*curtosis of allele lengths for each pair of population pooled together */
		if(lsstats[3]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats4_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] <= 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
							++skip;
							continue;
						}
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats4_2[ic] += calc_curt(ldna[cloc],tfreq);
					}
					if(nloc-skip==0)
						sstats4_2[ic]=0;
					else
						sstats4_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats4_2[ic]);
					++ic;
				}
			}
		}
		/*Shanon's index for each pair of populations pooled together */
		if(lsstats[4]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats5_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] < 1 || ltype[cloc] == 'S' || ltype[cloc] == 's'){
							++skip;
							continue;
						}
						tsamp = nsamp[cloc][cpop]+nsamp[cloc][cpop2];
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats5_2[ic] += calc_shanon(ldna[cloc],tsamp,tfreq,cpop);
					}
					if(nloc-skip==0)
						sstats5_2[ic]=0;
					else
						sstats5_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats5_2[ic]);
					++ic;
				}
			}
		}
	}
	if(foundSNP){
		/*pairwise difference for each population */
		if(lsstats[6]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats7[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats7[cpop] += calc_pi(nsamp[cloc][cpop],freq[cloc][cpop],lsites[cloc],valS[cloc],ldna[cloc]);
				}
				if(nloc-skip==0)
					sstats7[cpop]=0;
				else
					sstats7[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats7[cpop]);
			}
		}
		//fill lsnp
		if(lsstats[7]==1){
			for(cpop=0; cpop<npop; cpop++){
				for(cloc=0; cloc<nloc; cloc++){
					for(csite=0;csite<lsites[cloc];csite++)
						lsnp[cpop][cloc][csite]=0;
				}
			}
		}
		/*number of segregating sites for each population */
		if(lsstats[7]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats8[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					snp[cpop][cloc]=calc_segsites(nsamp[cloc][cpop],freq[cloc][cpop],lsites[cloc],
												   valS[cloc],ldna[cloc],lsnp[cpop][cloc]);
					sstats8[cpop]+=snp[cpop][cloc];
				}
				if(nloc-skip==0)
					sstats8[cpop]=0;
				else
					sstats8[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats8[cpop]);
			}
		}
		/*number of haplotypes for each population */
		if(lsstats[8]==1){
			for(cpop = 0;cpop <npop;++cpop)	{
				sstats9[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc)
				{
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats9[cpop] += calc_nhap(ldna[cloc],freq[cloc][cpop]);
				}
				if(nloc-skip==0)
					sstats9[cpop]=0;
				else
					sstats9[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats9[cpop]);
			}
		}
		/*Shanon's index for each population */
		if(lsstats[9]==1){
			for(cpop=0 ; cpop<npop ; ++cpop){
				sstats10[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats10[cpop] += calc_shanon(ldna[cloc],nsamp[cloc][cpop],freq[cloc][cpop],cpop);
				}
				if(nloc-skip==0)
					sstats10[cpop]=0;
				else
					sstats10[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats10[cpop]);
			}
		}
		/*build a mutation frequency spectrum (MFS)*/
		if(lsstats[10]==1 || lsstats[11]==1){
			for(cpop = 0;cpop <npop;++cpop)	{
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						continue;
					}
					fillMFS(valS[cloc],freq[cloc][cpop],ldna[cloc],lsites[cloc],nsamp[cloc][cpop],cloc,cpop,mfs);
				}
			}
		}
		/*mean of mutation frequency spectrum (MFS) in each population*/
		if(lsstats[10]==1){
			for(cpop = 0;cpop <npop;++cpop)	{
				sstats11[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					mMFS[cpop][cloc]=calc_meanMFS(lsites[cloc],mfs[cpop][cloc],snp[cpop][cloc]);
					sstats11[cpop]+=mMFS[cpop][cloc]; 
				}
				if(nloc-skip==0)
					sstats11[cpop]=0;
				else
					sstats11[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats11[cpop]);
			}
		}
		/*stdev of mutation frequency spectrum (MFS) in each population*/
		if(lsstats[11]==1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats12[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats12[cpop]+= calc_stdevMFS(lsites[cloc],mfs[cpop][cloc],snp[cpop][cloc],lsnp[cpop][cloc],mMFS[cpop][cloc]);
				}
				if(nloc-skip==0)
					sstats12[cpop]=0;
				else
					sstats12[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats12[cpop]);
			}
		}
		/*number of segregate sites for all the populations pooled together */
		if(lsstats[12]==1 && npop>1){
			tsamp=0;
			for(cdna=0;cdna<maxDna;++cdna)
				tfreq[cdna]=0;
			for(cpop = 0;cpop <npop;++cpop){
				for(cloc=0;cloc<nloc;++cloc){
					tsamp += nsamp[cloc][cpop];
					for(cdna=0;cdna<ldna[cloc];++cdna)
						tfreq[cdna] += freq[cloc][cpop][cdna];
				}
			}
			sstats8_t = 0;
			skip = 0;
			for(cloc=0;cloc<nloc;++cloc)
			{
				if(tsamp <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm')
				{
					++skip;
					continue;
				}
				sstats8_t += calc_segsites2(tsamp,tfreq,lsites[cloc],valS[cloc],ldna[cloc]);
			}
			if(nloc-skip==0)
				sstats8_t=0;
			else
				sstats8_t /= nloc-skip;
		}
		/*Nm_S in each population */
		if(lsstats[12]==1 && npop>1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats13[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats13[cpop] += calc_Nm_S(sstats8[cpop],sstats8_t);
				}
				if(nloc-skip==0)
					sstats13[cpop]=0;
				else
					sstats13[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats13[cpop]);
			}
		}
		/*private S in each population */
		if(lsstats[13]==1 && npop>1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats14[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats14[cpop]+=calc_privS(lsites[cloc],npop,cpop,cloc,lsnp);
				}
				if(nloc-skip==0)
					sstats14[cpop]=0;
				else
					sstats14[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats14[cpop]);
			}
		}
		/*S(1) in each population */
		if(lsstats[14]==1 && npop>1){
			for(cpop = 0;cpop <npop;++cpop){
				sstats15[cpop] = 0;
				skip = 0;
				for(cloc=0;cloc<nloc;++cloc){
					if(nsamp[cloc][cpop] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
						++skip;
						continue;
					}
					sstats15[cpop] += calc_S_1(lsites[cloc],npop,ldna[cloc],cpop,cloc,lsnp,freq[cloc],valS[cloc]);
				}
				if(nloc-skip==0)
					sstats15[cpop]=0;
				else
					sstats15[cpop] /= nloc-skip;
				sprintf(outp+strlen(outp),"%g ",sstats15[cpop]);
			}
		}
		/*pairwise difference for each pair of populations pooled together */
		if(lsstats[6]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats7_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
							++skip;
							continue;
						}
						tsamp = nsamp[cloc][cpop]+nsamp[cloc][cpop2];
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats7_2[ic] += calc_pi(tsamp,tfreq,lsites[cloc],valS[cloc],ldna[cloc]);
					}
					if(nloc-skip==0)
						sstats7_2[ic]=0;
					else
						sstats7_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats7_2[ic]);
					++ic;
				}
			}
		}
		/*number of segregate sites for each pair of populations pooled together */
		if(lsstats[7]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats8_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] <= 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
							++skip;
							continue;
						}
						tsamp = nsamp[cloc][cpop]+nsamp[cloc][cpop2];
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats8_2[ic] += calc_segsites2(tsamp,tfreq,lsites[cloc],valS[cloc],ldna[cloc]);
					}
					if(nloc-skip==0)
						sstats8_2[ic]=0;
					else
						sstats8_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats8_2[ic]);
					++ic;
				}
			}
		}
		/*number of haplotypes for each pair of populations pooled together */
		if(lsstats[8]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats9_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
							++skip;
							continue;
						}
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats9_2[ic] += calc_nhap(ldna[cloc],tfreq);
					}
					if(nloc-skip==0)
						sstats9_2[ic]=0;
					else
						sstats9_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats9_2[ic]);
					++ic;
				}
			}
		}
		/*Shanon's index for each pair of populations pooled together */
		if(lsstats[9]==1){
			for(cpop = 0,ic=0;cpop <npop;++cpop){
				for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
					sstats10_2[ic] = 0;
					skip = 0;
					for(cloc=0;cloc<nloc;++cloc){
						if(nsamp[cloc][cpop] + nsamp[cloc][cpop2] < 1 || ltype[cloc] == 'M' || ltype[cloc] == 'm'){
							++skip;
							continue;
						}
						tsamp = nsamp[cloc][cpop]+nsamp[cloc][cpop2];
						for(cdna=0;cdna<ldna[cloc];++cdna)
							tfreq[cdna] = freq[cloc][cpop][cdna] + freq[cloc][cpop2][cdna];
						sstats10_2[ic] += calc_shanon(ldna[cloc],tsamp,tfreq,cpop);
					}
					if(nloc-skip==0)
						sstats10_2[ic]=0;
					else
						sstats10_2[ic] /= nloc-skip;
					sprintf(outp+strlen(outp),"%g ",sstats10_2[ic]);
					++ic;
				}
			}
		}
	}
	//write informations to string outp
	sprintf(outp+strlen(outp),"\n");

	/*free stuff*/
	free(tfreq);
	if(foundSTR){
		free(sstats1);
		free(sstats2);
		free(sstats3);
		free(sstats4);
		free(sstats5);
		free(sstats6);
		free(sstats1_2);
		free(sstats2_2);
		free(sstats3_2);
		free(sstats4_2);
		free(sstats5_2);
	}
	if(foundSNP){
		for(cpop=0; cpop<npop; cpop++){
			for(cloc=0; cloc<nloc; cloc++){
				free(mfs[cpop][cloc]);
				free(lsnp[cpop][cloc]);
			}
			free(snp[cpop]);
			free(lsnp[cpop]);
			free(mMFS[cpop]);
			free(mfs[cpop]);
		}
		free(snp);
		free(lsnp);
		free(mMFS);
		free(mfs);
		free(sstats7);
		free(sstats8);
		free(sstats9);
		free(sstats10);
		free(sstats11);
		free(sstats12);
		free(sstats13);
		free(sstats14);
		free(sstats15);
		free(sstats7_2);
		free(sstats8_2);
		free(sstats9_2);
		free(sstats10_2);
	}

} //end of get_summstats

int calc_nhap(int nhap, int freq[]){
	int chap,	//iterator
		ndna;

	ndna = 0;
	for(chap=0;chap<nhap;chap++)
		if(freq[chap] > 0)
			++ndna;
	return ndna;

} //end of calc_nhap

double calc_shanon(int nhap, int nsamp, int *freq, int cpop){
	int chap;		//iterators
	double shan,
		   p;

	shan = 0.0;
	for(chap=0;chap<nhap;++chap){
		if(freq[chap]!=0&&nsamp!=0){
			p = freq[chap]/(double)nsamp;
			shan -= p*(log(p)/log(2.0));
		}
	}
	return(shan);

} //end of calc_shanon

double calc_var(int sum,int freq[],int val[],int nhap){
	int chap;	//iterator
	double x,
		   x2;

	x = 0;
	for(chap=0;chap<nhap;++chap){
		x += val[chap]*freq[chap];
	}
	x /= sum;
	
	x2 = 0;
	for(chap=0;chap<nhap;++chap){
		x2 += (val[chap]-x)*(val[chap]-x)*freq[chap];
	}
	x2 /= sum-1;
	return(x2);

} //end of calc_var

double calc_curt(int nhap, int *freq){
	int chap;	//iterator
	double curt,
		   *x1,
		   *x2,
		   *x3,
		   *x4,
		   *min,
		   *max,
		   *dfreq;

	dfreq=(double*)malloc(nhap*sizeof(double));
	x1=(double*)malloc(nhap*sizeof(double));
	x2=(double*)malloc(nhap*sizeof(double));
	x3=(double*)malloc(nhap*sizeof(double));
	x4=(double*)malloc(nhap*sizeof(double));
	min=(double*)malloc(nhap*sizeof(double));
	max=(double*)malloc(nhap*sizeof(double));
	for(chap=0 ; chap<nhap ; chap++)
		dfreq[chap]=freq[chap];
	curt=0.0;
	mom(dfreq,nhap,x1,x2,x3,x4,min,max);
	curt=*x4;

	free(dfreq);
	free(x1);
	free(x2);
	free(x3);
	free(x4);
	free(min);
	free(max);
	return curt;

} //end of calc_curt

double calc_het(int sum,int freq[], int nhap){
	int chap;	//iterator
	double het;
	
	het = 0;
	for(chap=0;chap<nhap;++chap){
		het += freq[chap]*freq[chap];
	}
	het = sum/(sum-1.0)*(1-het/(sum*sum));
	return(het);

} //end of calc_het

double calc_Nm_H(double Hw,double Ha){
	if(Ha-Hw<=0)
		return Hw;
	else
		return Hw/(1 + Ha - Hw);

} //end of calc_Nm_H

double calc_pi(int sum,int freq[], int nsites, char **val,int nhap){
	int chap,chap2;		//iterators
	double het;
	
	het = 0;
	for(chap=0;chap<nhap;++chap){
		for(chap2=chap+1;chap2<nhap;++chap2){
			het += freq[chap]*freq[chap2]*numdiff(nsites,val[chap],val[chap2]);
		}
	}
	het /= sum*(sum-1)/2;
	return(het);

} //end of calc_pi

int numdiff(int nsites, char *hap1, char *hap2){
	int csite,
		sum;

	sum = 0;
	for(csite=0;csite<nsites;++csite)
		if(hap1[csite] != hap2[csite])
			++sum;
	return sum;

} //end of numdiff

double calc_segsites(int sum,int freq[],int nsites, char **val,int nhap,int *lsnp){
	int csite,chap,		//iterators
		ip;
	double nS;

	nS = 0;
	/* find first haplotype that has freq > 0 in this pop */
	for(chap=0;chap<nhap;++chap)
		if(freq[chap] > 0)
			break;
	ip = chap;
	for(csite=0;csite<nsites;++csite){
		for(chap=1;chap<nhap;++chap){
			if(freq[chap] == 0)continue; /* ignore haplotype with 0 freq */
			if(val[chap][csite] != val[ip][csite]){
				nS++;
				lsnp[csite]=1;
				break;
			}
		}
	}
	return nS;

} //end of calc_segsites

double calc_segsites2(int sum,int freq[],int nsites, char **val,int nhap){
	int csite,chap,		//iterators
		ip;
	double x;

	x = 0;
	/* find first haplotype that has freq > 0 in this pop */
	for(chap=0;chap<nhap;++chap)
		if(freq[chap] > 0)
			break;
	ip = chap;
	for(csite=0;csite<nsites;++csite){
		for(chap=1;chap<nhap;++chap){
			if(freq[chap] == 0)continue; /* ignore haplotype with 0 freq */
			if(val[chap][csite] != val[ip][csite]){
				++x;
				break;
			}
		}
	}
	return x;

} //end of calc_segsites

void fillMFS(char **val,int *freq,int nhap,int nsites,int nsamp,int cloc,int cpop,int ***mfs){
	int chap,csite,		//iterators
		aux;
		
	//fill mfs[][][]
	for(csite=0; csite<nsites; csite++){
		aux=0;
		for(chap=0; chap<nhap; chap++){
			aux+=freq[chap]*segtoint(val[chap][csite]);
		}
		if(aux>nsamp/2)
			aux=nsamp-aux;
		mfs[cpop][cloc][csite]=aux;
	}

} //end of fillMFS

int segtoint(char c1){
	if(c1=='0')
		return 0;
	else
		return 1;

} //end of segtoint

double calc_meanMFS(int nsites,int *mfs,int snp){
	int csite;			//iterators
	double aver;		//average of MFS

	aver=0;
	for(csite=0; csite<nsites; csite++){
		aver+=mfs[csite];
	}
	if(snp>0)
		return aver/(double)snp;
	else
		return 0;

} //end of calc_meanMFS

double calc_stdevMFS(int nsites,int *mfs,int snp,int *lsnp,double mMFS){
	int csite;			//iterators
	double stdev;		//stdev of MFS

	if(snp>1){
		stdev = 0;
		for(csite=0; csite<nsites; csite++){
			if(lsnp[csite])
				stdev += pow(mfs[csite]-mMFS,2);
		}
		stdev = stdev/((double)snp-1);
		return pow(stdev,0.5);
	}
	else
		return 0;

} //end of calc_stdevMFS

double calc_Nm_S(double Sw,double Sa){
	if(Sa-Sw<=0)
		return Sw;
	else
		return Sw/(Sa - Sw);

} //end of calc_Nm_S

double calc_privS(int nsites, int npop, int pop, int loc,int ***lsnp){
	int cpop,csite,		//iterators
		priv,
		nS;
	
	nS = 0;
	for(csite=0;csite<nsites;++csite){
		if(lsnp[pop][loc][csite]==1){
			priv = 1;
			for(cpop=0; cpop<npop; cpop++){
				if(cpop==pop)
					continue;
				if(lsnp[cpop][loc][csite]==1){
					priv=0;
					break;
				}
			}
			if(priv){
				++nS;
			}
		}
	}
	return nS;

} //end of calc_privS

double calc_S_1(int nsites, int npop, int nhap, int pop, int loc, int ***lsnp, int **freq,char **val){
	int cpop,chap,csite,	//iterators
		gotS,
		priv,
		nS;
	double S_1;
	char otherS;

	nS = 0;
	S_1 = 0.0;
	for(csite=0;csite<nsites;++csite){
		if(lsnp[pop][loc][csite]==1){
			gotS = 0;
			priv = 1;
			for(cpop=0; cpop<npop; cpop++){
				if(cpop==pop)
					continue;
				if(lsnp[cpop][loc][csite]==1){
					priv=0;
					break;
				}
				else{
					for(chap=0; chap<nhap && !gotS; chap++){
						if(freq[cpop][chap] == 0)continue; /* ignore haplotype with 0 freq */
						otherS = val[chap][csite];
						gotS = 1;
					}
				}
			}
			if(priv){
				for(chap=0; chap<nhap; chap++){
					if(freq[pop][chap] == 0)continue; /* ignore haplotype with 0 freq */
					if(val[chap][csite]!=otherS)
						S_1 +=freq[pop][chap];
				}
				++nS;
			}
		}
	}
	if(nS == 0)
		S_1 = 0;
	else
		S_1 /= nS;
	return S_1;

} //end of calc_S_1
