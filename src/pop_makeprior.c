/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		12th May 2009
*/
#include "interface.h"
					 
int makeprior(char *output,int niter,int genet,int npop,int nloc,double *lplo,
			  char *ltype,struct prior pr_top,struct prior *pr_Ne,struct prior *pr_tev,
			  struct prior *pr_mig,struct prior pr_mutSTR,struct prior pr_mutSNP,
			  struct prior pr_recSTR,struct prior pr_recSNP,struct migweights migw){
	int i,cparam,cloc,cpop,ctev,cpop2,	//iterators
		allTev,							//check if Tev are defined with only one prior
		nNe,							//no. of population (modern and ancient)
		nmig,							//no. of migration rates (modern and ancient
		ntev;							//no. of time events
	char *prs_name;						//name of the output file
	FILE *out_prs;						//pntr to .prs

	prs_name = malloc(5 + strlen(output)*sizeof(char));
	strcpy(prs_name,output);
	out_prs = fopen(strcat(prs_name,".prs"),"w");
	if(out_prs == NULL)
		return 1; //couldn't create .prs file

	/*write the first line to .prs file*/
	fprintf(out_prs,"%d %d %d %d \n\n",niter,genet,npop,nloc);
	/*write heredity scalars per locus to .prs file*/
	for(cloc=0 ; cloc<nloc ; cloc++)
		fprintf(out_prs,"%.2lf ",lplo[cloc]);
	fprintf(out_prs,"\n\n");
	/*write DNA type per locus to .prs file*/
	for(cloc=0 ; cloc<nloc ; cloc++)
		fprintf(out_prs,"%c ",ltype[cloc]);
	/*write the top priors to .prs file*/
	ntev = npop-1;
	if(pr_top.type==0){
		fprintf(out_prs,"\n\n%d\n",pr_top.type);
	}
	else if(pr_top.type==1){
		fprintf(out_prs,"\n\n%d",pr_top.type);
		fprintf(out_prs," %.0lf\n",pr_top.p[0]);
	}
	else if(pr_top.type==2){
		fprintf(out_prs,"\n\n%d",pr_top.type);
		for(i=0; i<2*ntev; i++)
			fprintf(out_prs," %.0lf",pr_top.p[i]);
		fprintf(out_prs,"\n");
	}
	else if(pr_top.type==3){
		fprintf(out_prs,"\n\n%d %.0lf\n",pr_top.type,pr_top.p[0]);
	}
	else if(pr_top.type==4){
		fprintf(out_prs,"\n\n%d",pr_top.type);
		fprintf(out_prs," %.0lf %.0lf\n",pr_top.p[0],pr_top.p[1]);
	}
	else{
		fprintf(out_prs,"\n\n%d",pr_top.type);
		for(i=0; i<2*ntev+1; i++)
			fprintf(out_prs," %.0lf",pr_top.p[i]);
		fprintf(out_prs,"\n");
	}
	/*write the Ne priors to .prs file*/
	nNe = npop*2-1;
	for(i=0 ; i<nNe ; i++){
		fprintf(out_prs,"\n%d",pr_Ne[i].type);
		if(pr_Ne[i].type==1){
			for(cparam=0 ; cparam<2 ; cparam++)
				fprintf(out_prs," %g",pr_Ne[i].p[cparam]);
		}
		else if(pr_Ne[i].type==2){
			for(cparam=0 ; cparam<4 ; cparam++)
				fprintf(out_prs," %g",pr_Ne[i].p[cparam]);
		}
		if((i+1)==npop)
			fprintf(out_prs,"\n");
	}
	fprintf(out_prs,"\n");

	if(npop>1){
		/*write the tev priors to .prs file*/
		allTev=0;
		for(i=0 ; i<ntev && !allTev; i++){
			fprintf(out_prs,"\n%d",pr_tev[i].type);
			if(pr_tev[i].type==3||pr_tev[i].type==4){
				allTev = 1;
			}
			if(pr_tev[i].type==1||pr_tev[i].type==3){
				for(cparam=0 ; cparam<2 ; cparam++)
					fprintf(out_prs," %g",pr_tev[i].p[cparam]);
			}
			else if(pr_tev[i].type==2||pr_tev[i].type==4){
				for(cparam=0 ; cparam<4 ; cparam++)
					fprintf(out_prs," %g",pr_tev[i].p[cparam]);
			}
		}
		fprintf(out_prs,"\n");
		/*write the mig prior to .prs file*/
		nmig = npop*2-2;
		for(i=0 ; i<nmig ; i++){
			fprintf(out_prs,"\n%d",pr_mig[i].type);
			if(pr_mig[i].type==1||pr_mig[i].type==3){
				for(cparam=0 ; cparam<2 ; cparam++)
					fprintf(out_prs," %g",pr_mig[i].p[cparam]);
			}
			else if(pr_mig[i].type==2||pr_mig[i].type==4){
				for(cparam=0 ; cparam<4 ; cparam++)
					fprintf(out_prs," %g",pr_mig[i].p[cparam]);
			}
			if((i+1)==npop)
				fprintf(out_prs,"\n");
		}
		fprintf(out_prs,"\n");
	}
	/*write the mut prior to .prs file*/
	fprintf(out_prs,"\n%d ",pr_mutSTR.type);
	if(pr_mutSTR.type==1||pr_mutSTR.type==2){
		for(cparam=0 ; cparam<4 ; cparam++)
			fprintf(out_prs,"%g ",pr_mutSTR.p[cparam]);
	}
	fprintf(out_prs,"\n%d ",pr_mutSNP.type);
	if(pr_mutSNP.type==1||pr_mutSNP.type==2){
		for(cparam=0 ; cparam<4 ; cparam++)
			fprintf(out_prs,"%g ",pr_mutSNP.p[cparam]);
	}
	fprintf(out_prs,"\n");
	/*write the rec prior to .prs file*/
	fprintf(out_prs,"\n%d ",pr_recSTR.type);
	if(pr_recSTR.type==1||pr_recSTR.type==2){
		for(cparam=0 ; cparam<4 ; cparam++)
			fprintf(out_prs,"%g ",pr_recSTR.p[cparam]);
	}
	fprintf(out_prs,"\n%d ",pr_recSNP.type);
	if(pr_recSNP.type==1||pr_recSNP.type==2){
		for(cparam=0 ; cparam<4 ; cparam++)
			fprintf(out_prs,"%g ",pr_recSNP.p[cparam]);
	}
	/*write the mig weights to .prs file*/
	fprintf(out_prs,"\n");
	if(npop>2){
		if(migw.type==0)
			fprintf(out_prs,"\n0");
		else{
			fprintf(out_prs,"\n1");
			for(cpop=0;cpop<npop;cpop++){
				fprintf(out_prs,"\n");			
				for(ctev=0;ctev<ntev;ctev++){
					for(cpop2=0;cpop2<npop;cpop2++){
						fprintf(out_prs,"%g ",migw.m[cpop][ctev][cpop2]);
					}
					fprintf(out_prs,"\n");			
				}
			}
		}
	}
	/*write the legend*/
	if(npop==1){
		fprintf(out_prs,"----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"PopABC - Mark Beaumont & Joao Lopes                                             01/05/09\n\n");
		fprintf(out_prs,">no_iterations, generation_time, no_populations, no_loci\n\n");
		fprintf(out_prs,">escalar per locus (autosome - 1; X-linked - 0.75; Y-linked or mitDNA - 0.25)\n\n");
		fprintf(out_prs,">type of DNA data (s - sequence; m - microssatelites)\n\n");
		fprintf(out_prs,">topology:       0 - uniform distribution;\n");
		fprintf(out_prs,"                 3 - uniform distribution (and choose a Model marker).\n\n");
		fprintf(out_prs,">ne1 params:     1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n\n");
		fprintf(out_prs,">mutM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">mutS params\n\n");
		fprintf(out_prs,">recM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">recS params\n\n");
	}
	else if(npop==2){
		fprintf(out_prs,"\n----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"PopABC - Mark Beaumont & Joao Lopes                                             01/05/09\n\n");
		fprintf(out_prs,">no_iterations, generation_time, no_populations, no_loci\n\n");
		fprintf(out_prs,">escalar per locus (autosome - 1; X-linked - 0.75; Y-linked or mitDNA - 0.25)\n\n");
		fprintf(out_prs,">type of DNA data (s - sequence; m - microssatelites)\n\n");
		fprintf(out_prs,">topology:       0 - uniform distribution;\n");
		fprintf(out_prs,"                 3 - uniform distribution (and choose a Model marker).\n\n");
		fprintf(out_prs,">ne1 params:     1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,">ne2 params\n\n");
		fprintf(out_prs,">neanc1 params\n\n");
		fprintf(out_prs,">t1 params:      1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n\n");
		fprintf(out_prs,">mig1 params:    0 - zero migration;\n");
		fprintf(out_prs,"                 1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution;\n");
		fprintf(out_prs,"                 3 - uniform distribution (on number of migrations);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (on number of migrations).\n");
		fprintf(out_prs,"                 [for 3 and 4 real mig rate is calculated as nmig/Ne]\n");
		fprintf(out_prs,">mig2 params\n\n");
		fprintf(out_prs,">mutM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">mutS params\n\n");
		fprintf(out_prs,">recM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">recS params\n");
		fprintf(out_prs,"----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"Tree topology:\n\n");
		fprintf(out_prs,"      ||       PopA1\n");
		fprintf(out_prs,"      ||         |\n");
		fprintf(out_prs,"      ||         |\n");
		fprintf(out_prs,"    t1||     ---------\n");
		fprintf(out_prs,"      ||     |       |\n");
		fprintf(out_prs,"      ||     |       |\n");
		fprintf(out_prs,"      \\/   Pop1    Pop2\n\n");
	}
	else if(npop==3){
		fprintf(out_prs,"\n----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"PopABC - Mark Beaumont & Joao Lopes                                             01/05/09\n\n");
		fprintf(out_prs,">no_iterations, generation_time, no_populations, no_loci\n\n");
		fprintf(out_prs,">escalar per locus (autosome - 1; X-linked - 0.75; Y-linked or mitDNA - 0.25)\n\n");
		fprintf(out_prs,">type of DNA data (s - sequence; m - microssatelites)\n\n");
		fprintf(out_prs,">topology:       0 - uniform distribution;\n");
		fprintf(out_prs,"                 1 - choose topology from a list;\n");
		fprintf(out_prs,"                 2 - specify topology manually [e.g. ((Pop1,Pop2)Pop3) -> 1 2 2 3];\n");
		fprintf(out_prs,"                 3 - uniform distribution (and choose a Model marker);\n");
		fprintf(out_prs,"                 4 - choose topology from a list (and choose a Model marker);\n");
		fprintf(out_prs,"                 5 - specify topology manually (and choose a Model marker).\n\n");
		fprintf(out_prs,">ne1 params:     1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,">ne2 params\n");
		fprintf(out_prs,">ne3 params\n\n");
		fprintf(out_prs,">neanc1 params\n\n");
		fprintf(out_prs,">neanc2 params\n");
		fprintf(out_prs,">t1 params:      1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,"                 3 - uniform distribtuion (for all time events);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (for all time events).\n");
		fprintf(out_prs,"                 [for 1 and 2 t(n) is added to t(n+1)]\n");
		fprintf(out_prs,"                 [for 3 and 4 set only one priors for all t(n)]\n");
		fprintf(out_prs,">t2 params\n\n");
		fprintf(out_prs,">mig1 params:    0 - zero migration;\n");
		fprintf(out_prs,"                 1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution\n");
		fprintf(out_prs,"                 3 - uniform distribution (on number of migrations);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (on number of migrations).\n");
		fprintf(out_prs,"                 [for 3 and 4 real mig rate is calculated as nmig/Ne]\n");
		fprintf(out_prs,">mig2 params\n");
		fprintf(out_prs,">mig3 params\n\n");
		fprintf(out_prs,">miganc1 params\n\n");
		fprintf(out_prs,">mutM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">mutS params\n\n");
		fprintf(out_prs,">recM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">recS params\n\n");
		fprintf(out_prs,">migweight:      0 - do not use migweights matrix;\n");
		fprintf(out_prs,"                 1 - use migweights matrix as following:\n\n");
		fprintf(out_prs,"                 0     mw112 mw113\n");
		fprintf(out_prs,"                 0     mw122 mw123\n\n");
		fprintf(out_prs,"                 mw211 0     mw213\n");
		fprintf(out_prs,"                 mw221 0     mw223\n\n");
		fprintf(out_prs,"                 mw311 mw312 0\n");
		fprintf(out_prs,"                 mw321 mw322 0\n\n");
		fprintf(out_prs,"                 , where mwitj is the prob that the fraction of migrantes in pop i comes\n");
		fprintf(out_prs,"                 from pop j at a period of time before time event t. Sum of prob should\n");
		fprintf(out_prs,"                 be equal to 1.\n");
		fprintf(out_prs,"                 [only use migweight if the topology is specified (option 1,2,4 or 5)]\n");
		fprintf(out_prs,"----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"Tree topology:\n\n");
		fprintf(out_prs,"      ||             PopA2\n");
		fprintf(out_prs,"      ||               |\n");
		fprintf(out_prs,"    t2||         ------------\n");
		fprintf(out_prs,"      ||         |          |\n");
		fprintf(out_prs,"      ||       PopA1        |\n");
		fprintf(out_prs,"      ||         |          |\n");
		fprintf(out_prs,"      ||         |          |\n");
		fprintf(out_prs,"    t1||     ---------      |\n");
		fprintf(out_prs,"      ||     |       |      |\n");
		fprintf(out_prs,"      ||     |       |      |\n");
		fprintf(out_prs,"      \\/   Pop     Pop    Pop\n\n");
	}
	else if(npop==4){
		fprintf(out_prs,"\n----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"PopABC - Mark Beaumont & Joao Lopes                                             01/05/09\n\n");
		fprintf(out_prs,">no_iterations, generation_time, no_populations, no_loci\n\n");
		fprintf(out_prs,">escalar per locus (autosome - 1; X-linked - 0.75; Y-linked or mitDNA - 0.25)\n\n");
		fprintf(out_prs,">type of DNA data (s - sequence; m - microssatelites)\n\n");
		fprintf(out_prs,">topology:       0 - uniform distribution;\n");
		fprintf(out_prs,"                 1 - choose topology from a list;\n");
		fprintf(out_prs,"                 2 - specify topology manually [e.g. ((Pop1,Pop2)Pop3) -> 1 2 2 3];\n");
		fprintf(out_prs,"                 3 - uniform distribution (and choose a Model marker);\n");
		fprintf(out_prs,"                 4 - choose topology from a list (and choose a Model marker);\n");
		fprintf(out_prs,"                 5 - specify topology manually (and choose a Model marker).\n\n");
		fprintf(out_prs,">ne1 params:     1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,">ne2 params\n");
		fprintf(out_prs,">ne3 params\n");
		fprintf(out_prs,">ne4 params\n\n");
		fprintf(out_prs,">neanc1 params\n");
		fprintf(out_prs,">neanc2 params\n");
		fprintf(out_prs,">neanc3 params\n\n");
		fprintf(out_prs,">t1 params:      1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,"                 3 - uniform distribtuion (for all time events);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (for all time events).\n");
		fprintf(out_prs,"                 [for 1 and 2 t(n) is added to t(n+1)]\n");
		fprintf(out_prs,"                 [for 3 and 4 set only one priors for all t(n)]\n");
		fprintf(out_prs,">t2 params\n");
		fprintf(out_prs,">t3 params\n\n");
		fprintf(out_prs,">mig1 params:    0 - zero migration;\n");
		fprintf(out_prs,"                 1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution\n");
		fprintf(out_prs,"                 3 - uniform distribution (on number of migrations);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (on number of migrations).\n");
		fprintf(out_prs,"                 [for 3 and 4 real mig rate is calculated as nmig/Ne]\n");
		fprintf(out_prs,">mig2 params\n");
		fprintf(out_prs,">mig3 params\n");
		fprintf(out_prs,">mig4 params\n\n");
		fprintf(out_prs,">miganc1 params\n");
		fprintf(out_prs,">miganc2 params\n\n");
		fprintf(out_prs,">mutM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">mutS params\n\n");
		fprintf(out_prs,">recM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">recS params\n\n");
		fprintf(out_prs,">migweight:      0 - do not use migweights matrix;\n");
		fprintf(out_prs,"                 1 - use migweights matrix as following:\n\n");
		fprintf(out_prs,"                 0     mw112 mw113 mw114\n");
		fprintf(out_prs,"                 0     mw122 mw123 mw124\n");
		fprintf(out_prs,"                 0     mw132 mw133 mw134\n\n");
		fprintf(out_prs,"                 mw211 0     mw213 mw214\n");
		fprintf(out_prs,"                 mw221 0     mw223 mw224\n");
		fprintf(out_prs,"                 mw231 0     mw233 mw234\n\n");
		fprintf(out_prs,"                 mw311 mw312 0     mw314\n");
		fprintf(out_prs,"                 mw321 mw322 0     mw324\n");
		fprintf(out_prs,"                 mw331 mw332 0     mw334\n\n");
		fprintf(out_prs,"                 , where mwitj is the prob that the fraction of migrantes in pop i comes\n");
		fprintf(out_prs,"                 from pop j at a period of time before time event t. Sum of prob should\n");
		fprintf(out_prs,"                 be equal to 1.\n");
		fprintf(out_prs,"                 [only use migweight if the topology is specified (option 1,2,4 or 5)]\n");
		fprintf(out_prs,"----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"Tree topology:\n\n");
		fprintf(out_prs,"    ||                   PopA3                          PopA3\n");
		fprintf(out_prs,"    ||                     |                              |\n");
		fprintf(out_prs,"  t3||              -------------                ----------------\n");
		fprintf(out_prs,"    ||              |           |                |              |\n");
		fprintf(out_prs,"    ||            PopA2         |              PopA2            |\n");
		fprintf(out_prs,"    ||              |           |                |              |\n");
		fprintf(out_prs,"  t2||         -----------      |     OR      --------          |\n");
		fprintf(out_prs,"    ||         |         |      |             |      |          |\n");
		fprintf(out_prs,"    ||       PopA1       |      |             |      |        PopA1\n");
		fprintf(out_prs,"    ||         |         |      |             |      |          |\n");
		fprintf(out_prs,"  t1||     --------      |      |             |      |      --------\n");
		fprintf(out_prs,"    ||     |      |      |      |             |      |      |      |\n");
		fprintf(out_prs,"    \\/   Pop    Pop    Pop    Pop           Pop    Pop    Pop    Pop\n\n");
	}
	else if(npop==5){
		fprintf(out_prs,"\n----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"PopABC - Mark Beaumont & Joao Lopes                                             01/05/09\n\n");
		fprintf(out_prs,">no_iterations, generation_time, no_populations, no_loci\n\n");
		fprintf(out_prs,">escalar per locus (autosome - 1; X-linked - 0.75; Y-linked or mitDNA - 0.25)\n\n");
		fprintf(out_prs,">type of DNA data (s - sequence; m - microssatelites)\n\n");
		fprintf(out_prs,">topology:       0 - uniform distribution;\n");
		fprintf(out_prs,"                 1 - choose topology from a list;\n");
		fprintf(out_prs,"                 2 - specify topology manually [e.g. ((Pop1,Pop2)Pop3) -> 1 2 2 3];\n");
		fprintf(out_prs,"                 3 - uniform distribution (and choose a Model marker);\n");
		fprintf(out_prs,"                 4 - choose topology from a list (and choose a Model marker);\n");
		fprintf(out_prs,"                 5 - specify topology manually (and choose a Model marker).\n\n");
		fprintf(out_prs,">ne1 params:     1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,">ne2 params\n");
		fprintf(out_prs,">ne3 params\n");
		fprintf(out_prs,">ne4 params\n");
		fprintf(out_prs,">ne5 params\n\n");
		fprintf(out_prs,">neanc1 params\n");
		fprintf(out_prs,">neanc2 params\n");
		fprintf(out_prs,">neanc3 params\n");
		fprintf(out_prs,">neanc4 params\n\n");
		fprintf(out_prs,">t1 params:      1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,"                 3 - uniform distribtuion (for all time events);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (for all time events).\n");
		fprintf(out_prs,"                 [for 1 and 2 t(n) is added to t(n+1)]\n");
		fprintf(out_prs,"                 [for 3 and 4 set only one priors for all t(n)]\n");
		fprintf(out_prs,">t2 params\n");
		fprintf(out_prs,">t3 params\n");
		fprintf(out_prs,">t4 params\n\n");
		fprintf(out_prs,">mig1 params:    0 - zero migration;\n");
		fprintf(out_prs,"                 1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution\n");
		fprintf(out_prs,"                 3 - uniform distribution (on number of migrations);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (on number of migrations).\n");
		fprintf(out_prs,"                 [for 3 and 4 real mig rate is calculated as nmig/Ne]\n");
		fprintf(out_prs,">mig2 params\n");
		fprintf(out_prs,">mig3 params\n");
		fprintf(out_prs,">mig4 params\n");
		fprintf(out_prs,">mig5 params\n\n");
		fprintf(out_prs,">miganc1 params\n");
		fprintf(out_prs,">miganc2 params\n");
		fprintf(out_prs,">miganc3 params\n\n");
		fprintf(out_prs,">mutM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">mutS params\n\n");
		fprintf(out_prs,">recM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">recS params\n\n");
		fprintf(out_prs,">migweight:      0 - do not use migweights matrix;\n");
		fprintf(out_prs,"                 1 - use migweights matrix as following:\n\n");
		fprintf(out_prs,"                 0     mw112 mw113 mw114 mw115\n");
		fprintf(out_prs,"                 0     mw122 mw123 mw124 mw125\n");
		fprintf(out_prs,"                 0     mw132 mw133 mw134 mw135\n");
		fprintf(out_prs,"                 0     mw142 mw143 mw144 mw145\n\n");
		fprintf(out_prs,"                 mw211 0     mw213 mw214 mw215\n");
		fprintf(out_prs,"                 mw221 0     mw223 mw224 mw225\n");
		fprintf(out_prs,"                 mw231 0     mw233 mw234 mw235\n");
		fprintf(out_prs,"                 mw241 0     mw243 mw244 mw245\n\n");
		fprintf(out_prs,"                 mw311 mw312 0     mw314 mw315\n");
		fprintf(out_prs,"                 mw321 mw322 0     mw324 mw325\n");
		fprintf(out_prs,"                 mw331 mw332 0     mw334 mw335\n");
		fprintf(out_prs,"                 mw341 mw342 0     mw344 mw345\n\n");
		fprintf(out_prs,"                 mw411 mw412 mw413 0     mw415\n");
		fprintf(out_prs,"                 mw421 mw422 mw423 0     mw425\n");
		fprintf(out_prs,"                 mw431 mw432 mw433 0     mw435\n");
		fprintf(out_prs,"                 mw441 mw442 mw443 0     mw445\n\n");
		fprintf(out_prs,"                 mw511 mw512 mw513 mw514 0\n");
		fprintf(out_prs,"                 mw521 mw522 mw523 mw524 0\n");
		fprintf(out_prs,"                 mw531 mw532 mw533 mw534 0\n");
		fprintf(out_prs,"                 mw541 mw542 mw543 mw544 0\n\n");
		fprintf(out_prs,"                 , where mwitj is the prob that the fraction of migrantes in pop i comes\n");
		fprintf(out_prs,"                 from pop j at a period of time before time event t. Sum of prob should\n");
		fprintf(out_prs,"                 be equal to 1.\n");
		fprintf(out_prs,"                 [only use migweight if the topology is specified (option 1,2,4 or 5)]\n");
		fprintf(out_prs,"----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"Tree topology:\n\n");
		fprintf(out_prs,"    ||               Neanc                 Neanc                           Neanc\n");
		fprintf(out_prs,"    ||                 |                     |                               |\n");
		fprintf(out_prs,"  t4||             --------             ----------                      -----------\n");
		fprintf(out_prs,"    ||             |      |             |        |                      |         |\n");
		fprintf(out_prs,"    ||           Neanc    |             |      Neanc                  Neanc       |\n");
		fprintf(out_prs,"    ||             |      |             |        |                      |         |\n");
		fprintf(out_prs,"  t3||         --------   |             |     -------               ---------     |\n");
		fprintf(out_prs,"    ||         |      |   |             |     |     |               |       |     |\n");
		fprintf(out_prs,"    ||       Neanc    |   |     OR      |     |   Neanc     OR      |     Neanc   |\n");
		fprintf(out_prs,"    ||         |      |   |             |     |     |               |       |     |\n");
		fprintf(out_prs,"  t2||      -------   |   |             |     |   -----             |     -----   |\n");
		fprintf(out_prs,"    ||      |     |   |   |             |     |   |   |             |     |   |   |\n");
		fprintf(out_prs,"    ||    Neanc   |   |   |             |     |   |   |           Neanc   |   |   |\n");
		fprintf(out_prs,"    ||      |     |   |   |             |     |   |   |             |     |   |   |\n");
		fprintf(out_prs,"  t1||    -----   |   |   |           -----   |   |   |           -----   |   |   |\n");
		fprintf(out_prs,"    ||    |   |   |   |   |           |   |   |   |   |           |   |   |   |   |\n");
		fprintf(out_prs,"    \\/    Ne  Ne  Ne  Ne  Ne          Ne  Ne  Ne  Ne  Ne          Ne  Ne  Ne  Ne  Ne\n\n");
	}
	else{
		fprintf(out_prs,"\n----------------------------------------------------------------------------------------\n");
		fprintf(out_prs,"PopABC - Mark Beaumont & Joao Lopes                                             01/05/09\n\n");
		fprintf(out_prs,">no_iterations, generation_time, no_populations, no_loci\n\n");
		fprintf(out_prs,">escalar per locus (autosome - 1; X-linked - 0.75; Y-linked or mitDNA - 0.25)\n\n");
		fprintf(out_prs,">type of DNA data (s - sequence; m - microssatelites)\n\n");
		fprintf(out_prs,">topology:       2 - specify topology manually [e.g. ((Pop1,Pop2)Pop3) -> 1 2 2 3];\n");
		fprintf(out_prs,"                 5 - specify topology manually (and choose a Model marker).\n\n");
		fprintf(out_prs,">ne1 params:     1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,">ne2 params\n");
		fprintf(out_prs,">ne3 params\n");
		fprintf(out_prs,">ne4 params\n");
		fprintf(out_prs,">ne5 params\n");
		fprintf(out_prs,">(...)\n\n");
		fprintf(out_prs,">neanc1 params\n");
		fprintf(out_prs,">neanc2 params\n");
		fprintf(out_prs,">neanc3 params\n");
		fprintf(out_prs,">neanc4 params\n");
		fprintf(out_prs,">(...)\n\n");
		fprintf(out_prs,">t1 params:      1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution.\n");
		fprintf(out_prs,"                 3 - uniform distribtuion (for all time events);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (for all time events).\n");
		fprintf(out_prs,"                 [for 1 and 2 t(n) is added to t(n+1)]\n");
		fprintf(out_prs,"                 [for 3 and 4 set only one priors for all t(n)]\n");
		fprintf(out_prs,">t2 params\n");
		fprintf(out_prs,">t3 params\n");
		fprintf(out_prs,">t4 params\n");
		fprintf(out_prs,">(...)\n\n");
		fprintf(out_prs,">mig1 params:    0 - zero migration;\n");
		fprintf(out_prs,"                 1 - uniform distribtuion;\n");
		fprintf(out_prs,"                 2 - generalized gamma distribution\n");
		fprintf(out_prs,"                 3 - uniform distribution (on number of migrations);\n");
		fprintf(out_prs,"                 4 - generalized gamma distribution (on number of migrations).\n");
		fprintf(out_prs,"                 [for 3 and 4 real mig rate is calculated as nmig/Ne]\n");
		fprintf(out_prs,">mig2 params\n");
		fprintf(out_prs,">mig3 params\n");
		fprintf(out_prs,">mig4 params\n");
		fprintf(out_prs,">mig5 params\n");
		fprintf(out_prs,">(...)\n\n");
		fprintf(out_prs,">miganc1 params\n");
		fprintf(out_prs,">miganc2 params\n");
		fprintf(out_prs,">miganc3 params\n");
		fprintf(out_prs,">(...)\n\n");
		fprintf(out_prs,">mutM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">mutS params\n\n");
		fprintf(out_prs,">recM params:    0 - zero mutation;\n");
		fprintf(out_prs,"                 1 - lognormal distribution: (mean of mean(log10); stdev of mean(log10);\n");
		fprintf(out_prs,"                 mean of Sdev(log10); stdev of stdev(log10). Stdev truncated at 0.\n");
		fprintf(out_prs,"                 2 - normal distribution: (mean of mean; stdev of mean; mean of Sdev;\n");
		fprintf(out_prs,"                 stdev of stdev. Stdev truncated at 0.\n");
		fprintf(out_prs,">recS params\n\n");
		fprintf(out_prs,">migweight:      0 - do not use migweights matrix;\n");
		fprintf(out_prs,"                 1 - use migweights matrix as following:\n\n");
		fprintf(out_prs,"                 0     mw112 mw113 mw114 mw115 ...\n");
		fprintf(out_prs,"                 0     mw122 mw123 mw124 mw125 ...\n");
		fprintf(out_prs,"                 0     mw132 mw133 mw134 mw135 ...\n");
		fprintf(out_prs,"                 0     mw142 mw143 mw144 mw145 ...\n");
		fprintf(out_prs,"                 ...   ...   ...   ...   ...   ...\n\n");
		fprintf(out_prs,"                 mw211 0     mw213 mw214 mw215 ...\n");
		fprintf(out_prs,"                 mw221 0     mw223 mw224 mw225 ...\n");
		fprintf(out_prs,"                 mw231 0     mw233 mw234 mw235 ...\n");
		fprintf(out_prs,"                 mw241 0     mw243 mw244 mw245 ...\n");
		fprintf(out_prs,"                 ...   ...   ...   ...   ...   ...\n\n");
		fprintf(out_prs,"                 mw311 mw312 0     mw314 mw315 ...\n");
		fprintf(out_prs,"                 mw321 mw322 0     mw324 mw325 ...\n");
		fprintf(out_prs,"                 mw331 mw332 0     mw334 mw335 ...\n");
		fprintf(out_prs,"                 mw341 mw342 0     mw344 mw345 ...\n");
		fprintf(out_prs,"                 ...   ...   ...   ...   ...   ...\n\n");
		fprintf(out_prs,"                 mw411 mw412 mw413 0     mw415 ...\n");
		fprintf(out_prs,"                 mw421 mw422 mw423 0     mw425 ...\n");
		fprintf(out_prs,"                 mw431 mw432 mw433 0     mw435 ...\n");
		fprintf(out_prs,"                 mw441 mw442 mw443 0     mw445 ...\n");
		fprintf(out_prs,"                 ...   ...   ...   ...   ...   ...\n\n");
		fprintf(out_prs,"                 mw511 mw512 mw513 mw514 0     ...\n");
		fprintf(out_prs,"                 mw521 mw522 mw523 mw524 0     ...\n");
		fprintf(out_prs,"                 mw531 mw532 mw533 mw534 0     ...\n");
		fprintf(out_prs,"                 mw541 mw542 mw543 mw544 0     ...\n");
		fprintf(out_prs,"                 ...   ...   ...   ...   ...   ...\n\n");
		fprintf(out_prs,"                 ...\n\n");
		fprintf(out_prs,"                 , where mwitj is the prob that the fraction of migrantes in pop i comes\n");
		fprintf(out_prs,"                 from pop j at a period of time before time event t. Sum of prob should\n");
		fprintf(out_prs,"                 be equal to 1.\n");
		fprintf(out_prs,"                 [only use migweight if the topology is specified (option 1,2,4 or 5)]\n");
	}
	
	/*free stuff*/
	free(prs_name);
	free(lplo);
	free(ltype);
	if(pr_top.type!=0)
		free(pr_top.p);
	for(cpop=0; cpop<npop*2-1; cpop++)
		free(pr_Ne[cpop].p);
	free(pr_Ne);	
	if(npop>1){
		if(pr_tev[0].type!=3 && pr_tev[0].type!=4){
			for(ctev=0; ctev<ntev; ctev++)
				free(pr_tev[ctev].p);
			free(pr_tev);
		}
		else{
			free(pr_tev[0].p);
			free(pr_tev);
		}
		for(cpop=0; cpop<npop*2-2; cpop++)
			if(pr_mig[cpop].type!=0)
				free(pr_mig[cpop].p);
		free(pr_mig);
	}
	if(pr_mutSTR.type!=0)
		free(pr_mutSTR.p);
	if(pr_mutSNP.type!=0)
		free(pr_mutSNP.p);
	if(pr_recSTR.type!=0)
		free(pr_recSTR.p);
	if(pr_recSNP.type!=0)
		free(pr_recSNP.p);
	if(migw.type==1){
		for(cpop=0;cpop<npop;cpop++){
			for(ctev=0;ctev<ntev;ctev++)
				free(migw.m[cpop][ctev]);
			free(migw.m[cpop]);
		}
		free(migw.m);
	}
	
	
	//close files
	fclose(out_prs);

	//no errors
	return 0;

} //end of main
