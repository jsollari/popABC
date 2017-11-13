/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		1st May 2009
	
	NBB - based on Mark's genparmf.c and genparmfS.c
*/

#include "abc.h"

/*
	This function calculates the number of topologies according to the npop
	
	@param npop - numer of populations
	@result number of topologies
*/
int notop(int npop);

/*
	This function frees the memory used for the DNA sequence data simulations 

	@param nloc - number of loci
  	@param na	- number of diferente sequencies
  	@param val	- all the diferente sequencies by loci
  	@param ltype - DNA type per loci
*/
void freeMem(int nloc, int *na,char ***val,char * ltype);

/*
	Start of the ABC program. Simulates coordenated points (summary_statistics, parameters) to be
	used in an Aproximate Bayesian Computation model to estimate the true values of a genealogical tree.
	
	@arg .prs input filename 
	@arg .ssz input filename
	@arg .sst input filename
	@arg .dat output filename
	@arg print or not the .len and .frq files of every genetic tree (0-don't print; 1-print)
	@arg print or not the .mut (0-don't print; 1-print)
	@arg print or not the .rec (0-don't print; 1-print)
*/
int main(int argc, char* argv[]){
	int csamp,cpop,cpop2,cloc,cevt,i,cSTR,	//iterators
		foundSTR,				//check if STR's are present
		foundSNP,				//check if SNP's are present
		nrSTR,					//number of linked STR's
		npop,					//nr. populations
		nloc,					//nr. loci
		nevt,					//nr. time events
		nstats_aux,				//auxiliar to get the nr. of summary statistics used
		nstats,					//nr. of summary statistics used
		nparams,				//nr. of parameters used
		printIt,				//print or not the information about frequency of haplotyes, sizes of haplotypes, number of diferent sizes, number of loci
		printMut,				//print or not to a file the mutation rate of each locus
		printRec,				//print or not to a file the recombination rate of each locus
		cut,					//auxiliar to work around strings
		homesize,				//size of the path from where the program is being run
		outsize,				//size of output filename
		started,				//auxiliar to help print lsstats2 into the file
		lsstats[MAXSSTATS],		//list of the sstats (0 - absent; 1 - present)
		*sampleSTR;				//sample sizes of the linked STR's
	double	ploidySTR;			//ploidy of the linked STR's
	long citer,					//iterator
		 niter;					//nr. iterations
	char type,					//type of analysis
		 *lsstats2[MAXSSTATS] = {"H","varL","k_M","curL","sH_M","NmH","pi","S","k_S","sH_S","avMFS","sdMFS","NmS","privS","S(1)"},
		 *home,					//path where the program is being run
		 *path,					//path to which the output files are going to be store in
		 *outline1,				//stores parameters informations
		 *outline2, 			//stores summstats informations
		 *name_dat,				//output filename
		 *name_inf,				//output filename
		 *name_mut,				//output filename
		 *name_rec;				//output filename
	FILE *inp_prs,				//pntr to .prs
		 *inp_ssz,				//pntr to .ssz
		 *inp_sst,				//pntr to .ssz
		 *out_mut,				//pntr to .mut
		 *out_rec,				//pntr to .rec
		 *out_dat,				//pntr to .dat
		 *out_inf;				//pntr to .txt
	time_t startClock,			//time_t when the program starts
		   endClock;	 		//time_t when the program ends
	struct params pm;			//struct that stores the parameter set values
	struct data data;			//struct that stores the data of a genetic tree 
	const struct tm *startTime,	//struct time when the program starts
					*endTime;	//struct time when the program ends

	if(!(argc == 8))
		printerr("needs .prs file, .ssz file, .sst file, output filename (no extension), printIt, printMut, printRec");

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 

	inp_prs = fopen(argv[1],"r");			//input .prs
	if(inp_prs == NULL)
		printerr("cannot open .prs file");
	
	inp_ssz = fopen(argv[2],"r");			//input .ssz
	if(inp_ssz == NULL)
		printerr("cannot open .ssz file");

	inp_sst = fopen(argv[3],"r");			//input .sst
	if(inp_sst == NULL)
		printerr("cannot open .sst file");
	
	/*get the path from where the progam is being run*/
	homesize = strlen(argv[0])+ 5;
	home = (char *)malloc(homesize*sizeof(char));
	strcpy(home,argv[0]);
	for(cut=-1,i=homesize-1 ; home[i]!='/' && home[i]!='\\' && i >= 0 ; cut++,i--){
		home[i]='\0';
	}
	home = realloc(home,(homesize-cut)*sizeof(char));
	
	/*get the path to output files*/
	outsize = strlen(argv[4]) + 5;
	path = (char *)malloc(outsize*sizeof(char));
	strcpy(path,argv[4]);
	for(cut=-1,i=outsize-1 ; path[i]!='/' && path[i]!='\\' && i >= 0 ; cut++,i--){
		path[i]='\0';
	}
	path = realloc(path,(outsize-cut)*sizeof(char));
	
	name_dat = (char *)malloc(outsize*sizeof(char));
	strcpy(name_dat,argv[4]);
	out_dat = fopen(strcat(name_dat,".dat"),"w");		//output .dat
	free(name_dat);
	if(out_dat == NULL)
		printerr("cannot create .dat file");
	
	name_inf = (char *)malloc(outsize*sizeof(char));
	strcpy(name_inf,argv[4]);
	out_inf = fopen(strcat(name_inf,".txt"),"w");		//output info.txt
	if(out_inf == NULL)
		printerr("cannot create .txt file");
	free(name_inf);	
	
	printIt = atoi(argv[5]);							//value of printIt
	
	printMut = atoi(argv[6]);							//value of printMut

	printRec = atoi(argv[7]);							//value of printRec
		
	fscanf(inp_prs,"%d %lf %d %d",&pm.niter,&pm.gent,&pm.npop,&pm.nloc);
	pm.nevt = pm.npop - 1;
	pm.ntop = notop(pm.npop);

	npop = pm.npop;
	nloc = pm.nloc;
	niter = pm.niter;
	nevt = pm.nevt;
 	
	if(npop>1){
		pm.mig = (double **)malloc(nevt*sizeof(double*));	
		pm.tev = (double *)malloc(nevt*sizeof(double));
	}
	pm.mu = (double *)malloc(nloc*sizeof(double));
	pm.rec = (double *)malloc(nloc*sizeof(double));
	pm.ploidy = (double *)malloc(nloc*sizeof(double));
	pm.type = (char *)malloc((nloc+1)*sizeof(char));
	pm.nsamp = (int **)malloc(nloc*sizeof(int *));
	for(cloc=0;cloc<nloc;++cloc)
		pm.nsamp[cloc] = (int *)malloc(npop*sizeof(int));
	for(cevt=0;cevt<nevt;cevt++)
		pm.mig[cevt] = (double *)malloc(npop*sizeof(double));
	pm.psize = (double **)malloc((nevt+1)*sizeof(double*));
	for(cevt=0;cevt<=nevt;cevt++)
		pm.psize[cevt] = (double *)malloc(npop*sizeof(double));

	for(cloc=0 ; cloc<nloc ; cloc++)
		fscanf(inp_prs,"%lf",&pm.ploidy[cloc]);
	foundSTR = 0;
	foundSNP = 0;
	for(cloc=0 ; cloc<nloc ; cloc++){
		fscanf(inp_prs,"%s",&pm.type[cloc]);
		if(pm.type[cloc]=='m'||pm.type[cloc]=='M')
			foundSTR = 1;
		if(pm.type[cloc]=='s'||pm.type[cloc]=='S')
			foundSNP = 1;
	}
	for(cloc=0;cloc<nloc;++cloc){
		for(cpop=0;cpop<npop;++cpop){
			if(EOF == fscanf(inp_ssz,"%d",&(pm.nsamp[cloc][cpop])))
				printerr("reading .ssz file");		
		}
	}	
	fclose(inp_ssz);
	for(i=0; i<MAXSSTATS; i++){
		while(getc(inp_sst)!='#');	
		if(EOF == fscanf(inp_sst,"%d",&lsstats[i]))
			printerr("reading .sst file");
	}
	fclose(inp_sst);

	/*check for errors when choosing summary statistics*/
	if(npop==1 && (lsstats[5] || lsstats[12] || lsstats[13] || lsstats[14]))
		printerr("following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)");
	//if we want to use the Nm_H we have to calculate H
	if(lsstats[5]== 1)
		lsstats[0] = 1;
	//if we want to use the sdMFS we have to calculate mMFS
	if(lsstats[11]== 1)
		lsstats[10] = 1;
	//if we want to use the mMFS, Nm_S, privateS or S(1) we have to calculate S
	if(lsstats[10]== 1 || lsstats[12]== 1 || lsstats[13]== 1 || lsstats[14]== 1 )
		lsstats[7] = 1;
	//counting the number of summary statistics
	if(npop>1){
		//applied to each population and to the populations pooled together
		nstats_aux = 0;
		if(foundSTR)
			for(i=0 ; i<MAXSSTATS_M ; i++){
				if((i==0||i==1||i==2||i==3||i==4) && lsstats[i])
					nstats_aux++;
			}
		if(foundSNP)
			for(i=MAXSSTATS_M ; i<MAXSSTATS ; i++){
				if((i==6||i==7||i==8||i==9) && lsstats[i])
					nstats_aux++;
			}
		nstats = npop*nstats_aux+(combinations(npop,2))*nstats_aux;	
		//applied only to each population
		nstats_aux = 0;
		if(foundSTR)
			for(i=0 ; i<MAXSSTATS ; i++){
				if((i==5) && lsstats[i])
					nstats_aux++;
			}
		if(foundSNP)
			for(i=MAXSSTATS_M ; i<MAXSSTATS ; i++){
				if((i==10||i==11||i==12||i==13||i==14) && lsstats[i])
					nstats_aux++;
			}
		nstats += npop*nstats_aux;
	}
	else{
		nstats = 0;
		if(foundSTR)
			for(i=0 ; i<MAXSSTATS_M ; i++){
				if(lsstats[i])
					nstats++;
			}
		if(foundSNP)
			for(i=MAXSSTATS_M ; i<MAXSSTATS ; i++){
				if(lsstats[i])
					nstats++;
			}
	}

	/*writting first part of output to .txt file*/
	fprintf(out_inf,"Simulate - Mark Beaumont & Joao Lopes\n\n");
	fprintf(out_inf," last updated 01/05/09\n\n");
	fprintf(out_inf,"INPUT AND STARTING INFORMATION\n");
	fprintf(out_inf,"-------------------------------\n");
	fprintf(out_inf,"Priors file :           %s\n",argv[1]);
	fprintf(out_inf,"Sample Sizes file:      %s\n",argv[2]);
	fprintf(out_inf,"SummaryStatistics file: %s\n\n",argv[3]);
	if(foundSTR && !foundSNP){
		fprintf(out_inf,"Type of genetic data analysed: microssatelites\n\n");
		fprintf(out_inf,"Number of summary statistics:  %d\n",nstats);
		if(lsstats[0])
			fprintf(out_inf," -heterozygosity\n"); 
		if(lsstats[1])
			fprintf(out_inf," -variance of alleles length\n"); 
		if(lsstats[2])
			fprintf(out_inf," -number of alleles\n"); 
		if(lsstats[3])
			fprintf(out_inf," -curtosis of alleles length\n"); 
		if(lsstats[4])
			fprintf(out_inf," -Shanon's index\n"); 
		if(lsstats[5])
			fprintf(out_inf," -Nm estimator based on H\n"); 
		fprintf(out_inf,"\nMutation model used: single-step model\n\n");
	}
	else if(foundSNP && !foundSTR){
		fprintf(out_inf,"Type of genetic data analysed: dna sequency\n\n");
		fprintf(out_inf,"Number of summary statistics:  %d\n",nstats);
		if(lsstats[6])
			fprintf(out_inf," -mean of pairwise differences\n"); 
		if(lsstats[7])
			fprintf(out_inf," -number of segregating sites\n"); 
		if(lsstats[8])
			fprintf(out_inf," -number of haplotypes\n"); 
		if(lsstats[9])
			fprintf(out_inf," -Shanon's index\n"); 
		if(lsstats[10])
			fprintf(out_inf," -mean of MFS\n"); 
		if(lsstats[11])
			fprintf(out_inf," -stdev of MFS\n"); 
		if(lsstats[12])
			fprintf(out_inf," -Nm estimator based on S\n"); 
		if(lsstats[13])
			fprintf(out_inf," -private S\n"); 
		if(lsstats[14])
			fprintf(out_inf," -S(1)\n"); 
		fprintf(out_inf,"\nMutation model used: infinite-segregated-sites model\n\n");
	}
	else{
		fprintf(out_inf,"Type of genetic data analysed: microssatelites\n");
		fprintf(out_inf,"                               dna sequency\n\n");
		fprintf(out_inf,"Number of summary statistics:  %d\n",nstats);
		fprintf(out_inf," Microssatelites:\n");
		if(lsstats[0])
			fprintf(out_inf," -heterozygosity\n"); 
		if(lsstats[1])
			fprintf(out_inf," -variance of alleles length\n"); 
		if(lsstats[2])
			fprintf(out_inf," -number of alleles\n"); 
		if(lsstats[3])
			fprintf(out_inf," -curtosis of alleles length\n"); 
		if(lsstats[4])
			fprintf(out_inf," -Shanon's index\n"); 
		if(lsstats[5])
			fprintf(out_inf," -Nm estimator based on H\n"); 
		fprintf(out_inf," Sequence data:\n");
		if(lsstats[6])
			fprintf(out_inf," -mean of pairwise differences\n"); 
		if(lsstats[7])
			fprintf(out_inf," -number of segregating sites\n"); 
		if(lsstats[8])
			fprintf(out_inf," -number of haplotypes\n"); 
		if(lsstats[9])
			fprintf(out_inf," -Shanon's index\n"); 
		if(lsstats[10])
			fprintf(out_inf," -mean of MFS\n"); 
		if(lsstats[11])
			fprintf(out_inf," -stdev of MFS\n"); 
		if(lsstats[12])
			fprintf(out_inf," -Nm estimator based on S\n"); 
		if(lsstats[13])
			fprintf(out_inf," -private S\n"); 
		if(lsstats[14])
			fprintf(out_inf," -S(1)\n"); 
		fprintf(out_inf,"\nMutation model used: \n");
		fprintf(out_inf,"  Microssatelites: single-step model\n");
		fprintf(out_inf,"  Sequence data: infinite-segregated-sites model\n\n");	
	}
	fprintf(out_inf,"Number of Modern Populations:  %d\n",npop);
	fprintf(out_inf,"Number of Loci analysed:       %d\n",nloc);
	fprintf(out_inf,"Type of the loci analysed:     ");
	for(cloc=0 ; cloc<nloc ; cloc++)
		fprintf(out_inf,"%c ",pm.type[cloc]);
	fprintf(out_inf,"\nScalar of the loci analysed:   ");
	for(cloc=0 ; cloc<nloc ; cloc++)
		fprintf(out_inf,"%.2lf ",pm.ploidy[cloc]);
	fprintf(out_inf,"\nGeneration time (in years):    %g\n\n",pm.gent);
	fprintf(out_inf,"RUN INFORMATION\n");
	fprintf(out_inf,"---------------------------\n");
	fprintf(out_inf,"Number of geneological trees created: %d\n\n",niter);
			
	/*Reads the topology type*/
	fscanf(inp_prs,"%d",&P_top.type);
	if(npop > 2){
		fprintf(out_inf,"Populational tree topology:\n");
		//uniform distribution along the number of different topologies
		if(P_top.type == 0)
			fprintf(out_inf,"- Prior distribution (Topol):  uniform(0,%d)\n",pm.ntop);
		//read values from file
		else if(P_top.type == 1){
			P_top.p = (double*)malloc(sizeof(double));
			fscanf(inp_prs,"%lf ",&P_top.p[0]);
			fprintf(out_inf,"- Fixed Topology: %.0lf\n",P_top.p[0]);
		}
		//discriminate the join events
		else if(P_top.type == 2){
			P_top.p = (double*)malloc(2*nevt*sizeof(double));
			for(i=0 ; i<2*nevt ; i++){
				fscanf(inp_prs,"%lf ",&P_top.p[i]);
			}
			fprintf(out_inf,"- Fixed Topology: (");
			for(i=0 ; i<2*nevt ; i++){
				fprintf(out_inf,"%.0lf",P_top.p[i]);
				if(i%2==0)
					fprintf(out_inf,",");
				else if(i!=2*nevt-1)
					fprintf(out_inf,")(");
			}
			fprintf(out_inf,")\n");
		}
		//uniform distribution along the number of different topologies (with specific marker)
		else if(P_top.type == 3){
			P_top.p = (double*)malloc(sizeof(double));
			fscanf(inp_prs,"%lf ",&P_top.p[0]);
			fprintf(out_inf,"- Prior distribution (Topol):  uniform(0,%d) [with Model marker:%.0lf]\n",pm.ntop,P_top.p[0]);
		}
		//read values from file (with specific marker)
		else if(P_top.type == 4){
			P_top.p = (double*)malloc(2*sizeof(double));
			fscanf(inp_prs,"%lf ",&P_top.p[0]);
			fscanf(inp_prs,"%lf ",&P_top.p[1]);
			fprintf(out_inf,"- Fixed Topology: %.0lf [with Model marker:%.0lf]\n",P_top.p[0],P_top.p[1]);
		}
		//discriminate the join events (with specific marker)
		else if(P_top.type == 5){
			P_top.p = (double*)malloc((1+(2*nevt))*sizeof(double));
			for(i=0 ; i<2*nevt ; i++){
				fscanf(inp_prs,"%lf ",&P_top.p[i]);
			}
			fscanf(inp_prs,"%lf ",&P_top.p[2*nevt]);
			fprintf(out_inf,"- Fixed Topology: (");
			for(i=0 ; i<2*nevt ; i++){
				fprintf(out_inf,"%.0lf",P_top.p[i]);
				if(i%2==0)
					fprintf(out_inf,",");
				else if(i!=2*nevt-1)
					fprintf(out_inf,")(");
			}
			fprintf(out_inf,") [with Model marker:%.0lf]\n",P_top.p[2*nevt]);
		}
		else
			printerr("wrong topology prior type");
	}
	if(npop<=2){
		fprintf(out_inf,"Populational tree topology:\n");
		if(P_top.type == 0)
			fprintf(out_inf,"- Single topology\n");
		else if(P_top.type == 3){
			P_top.p = (double*)malloc(sizeof(double));
			fscanf(inp_prs,"%lf ",&P_top.p[0]);
			fprintf(out_inf,"- Single topology [with Model marker:%.0lf]\n",P_top.p[0]);
		}
		else if(P_top.type == 1 || P_top.type == 2 || P_top.type == 4 || P_top.type == 5)
			printerr("with 1 or 2 pops only topology prior type 0 or 3 can be used");
		else
			printerr("wrong topology prior type");
	}
	/*Reads the Populations'Ne priors*/
	P_psize = (struct prior **)malloc((nevt+1)*sizeof(struct prior *));
	fprintf(out_inf,"Efective population size:\n");
	for(cevt=0;cevt<=nevt;++cevt){
		P_psize[cevt] = (struct prior *)malloc(npop*sizeof(struct prior));
		for(cpop=0;cpop<npop;++cpop){
			fscanf(inp_prs,"%d",&P_psize[cevt][cpop].type);
			//uniformly between 2 limits
			if(P_psize[cevt][cpop].type == 1){
				P_psize[cevt][cpop].p = (double*)malloc(2*sizeof(double));
				if(cevt==0)
					fprintf(out_inf,"- Prior distribution (Ne%d):    uniform(",cpop+1);
				else
					fprintf(out_inf,"- Prior distribution (NeAnc%d): uniform(",cevt);
				for(i=0;i<2;++i){
					fscanf(inp_prs,"%lf",&P_psize[cevt][cpop].p[i]);
					if(i==1)
						fprintf(out_inf,"%.0lf)\n",P_psize[cevt][cpop].p[i]);
					else	
						fprintf(out_inf,"%.0lf,",P_psize[cevt][cpop].p[i]);
				}
			}
			//use generalized gamma
			else if(P_psize[cevt][cpop].type == 2){
				P_psize[cevt][cpop].p = (double*)malloc(4*sizeof(double));
				if(cevt==0)
					fprintf(out_inf,"- Prior distribution (Ne%d):    generalized gamma(",cpop+1);
				else
					fprintf(out_inf,"- Prior distribution (NeAnc%d): generalized gamma(",cevt);
				for(i=0;i<4;++i){
					fscanf(inp_prs,"%lf",&P_psize[cevt][cpop].p[i]);
					if(i==3)
						fprintf(out_inf,"%.0lf)\n",P_psize[cevt][cpop].p[i]);
					else	
						fprintf(out_inf,"%.0lf,",P_psize[cevt][cpop].p[i]);
				}
			}
			else
				printerr("wrong Ne priors");
			//each anc pop is just done once
			if(cevt > 0)	
				break; 
		}
	}
	if(npop>1){
		/*Reads the events'time priors*/
		P_t = (struct prior *)malloc(nevt*sizeof(struct prior));
		fprintf(out_inf,"Time events:\n");
		for(cevt=0;cevt<nevt;++cevt){
			fscanf(inp_prs,"%d",&P_t[cevt].type);
			//use uniform distribution
			if(P_t[cevt].type == 1){
				P_t[cevt].p = (double*)malloc(2*sizeof(double));
				if(cevt>0){
					fprintf(out_inf,"- Prior distribution (tev%d):   tev%d + uniform(",cevt+1,cevt);
					for(i=0;i<2;++i){
						fscanf(inp_prs,"%lf",&P_t[cevt].p[i]);
						if(i==1)
							fprintf(out_inf,"%.0lf)\n",P_t[cevt].p[i]);
						else	
							fprintf(out_inf,"%.0lf,",P_t[cevt].p[i]);
					}
				}
                else{
                    fprintf(out_inf,"- Prior distribution (tev%d):   uniform(",cevt+1);
                    for(i=0;i<2;++i){
                        fscanf(inp_prs,"%lf",&P_t[cevt].p[i]);
                        if(i==1)
                            fprintf(out_inf,"%.0lf)\n",P_t[cevt].p[i]);
                        else    
                            fprintf(out_inf,"%.0lf,",P_t[cevt].p[i]);
                    }
                }
			}
			//use generalized gamma
			else if(P_t[cevt].type == 2){
				P_t[cevt].p = (double*)malloc(4*sizeof(double));
				if(cevt>0){
					fprintf(out_inf,"- Prior distribution (tev%d):   tev%d + generalized gamma(",cevt+1,cevt);
					for(i=0;i<4;++i){
						fscanf(inp_prs,"%lf",&P_t[cevt].p[i]);
						if(i==3)
							fprintf(out_inf,"%g)\n",P_t[cevt].p[i]);
						else	
							fprintf(out_inf,"%g,",P_t[cevt].p[i]);
					}
				}
                else{
                    fprintf(out_inf,"- Prior distribution (tev%d):   generalized gamma(",cevt+1);
                    for(i=0;i<4;++i){
                        fscanf(inp_prs,"%lf",&P_t[cevt].p[i]);
                        if(i==3)
                            fprintf(out_inf,"%g)\n",P_t[cevt].p[i]);
                        else    
                            fprintf(out_inf,"%g,",P_t[cevt].p[i]);
                    }
                }
			}
			//use uniform distribution (all tev)
			else if(P_t[cevt].type == 3){
				P_t[cevt].p = (double*)malloc(2*sizeof(double));
				fprintf(out_inf,"- Prior distribution (all tev): uniform(");
				for(i=0;i<2;++i){
					fscanf(inp_prs,"%lf",&P_t[cevt].p[i]);
					if(i==1)
						fprintf(out_inf,"%.0lf)\n",P_t[cevt].p[i]);
					else	
						fprintf(out_inf,"%.0lf,",P_t[cevt].p[i]);
				}
				if(cevt!=0)
					printerr("if using tev type 3 prior, define prior only once");
				break;
			}
			//use generalized gamma (all tev)
			else if(P_t[cevt].type == 4){
				P_t[cevt].p = (double*)malloc(4*sizeof(double));
				fprintf(out_inf,"- Prior distribution (all tev): generalized gamma(");
				for(i=0;i<4;++i){
					fscanf(inp_prs,"%lf",&P_t[cevt].p[i]);
					if(i==3)
						fprintf(out_inf,"%g)\n",P_t[cevt].p[i]);
					else	
						fprintf(out_inf,"%g,",P_t[cevt].p[i]);
				}
				if(cevt!=0)
					printerr("if using tev type 4 prior, define prior only once");
				break;
			}
			else 
				printerr("wrong tev priors");
		}
		/*Reads the migration priors*/
		P_mig = (struct prior **)malloc(nevt*sizeof(struct prior *));
		fprintf(out_inf,"Migration rates:\n");
		for(cevt=0;cevt<nevt;++cevt){
			P_mig[cevt] = (struct prior *)malloc(npop*sizeof(struct prior));
			for(cpop=0;cpop<npop;++cpop){
				fscanf(inp_prs,"%d",&P_mig[cevt][cpop].type);
				//no migration
				if(P_mig[cevt][cpop].type == 0){
					P_mig[cevt][cpop].p = (double*)malloc(2*sizeof(double));
					if(cevt==0)
						fprintf(out_inf,"- mig%d: No migration\n",cpop+1);
					else	
						fprintf(out_inf,"- migAnc%d: No migration\n",cevt);
				}
				//uniform	
				else if(P_mig[cevt][cpop].type == 1){
					P_mig[cevt][cpop].p = (double*)malloc(2*sizeof(double));
					if(cevt==0)
						fprintf(out_inf,"- Prior distribution (mig%d):    uniform(",cpop+1);
					else
						fprintf(out_inf,"- Prior distribution (migAnc%d): uniform(",cevt);
					for(i=0;i<2;++i){
						fscanf(inp_prs,"%lf",&P_mig[cevt][cpop].p[i]);
						if(i==1)
							fprintf(out_inf,"%g)\n",P_mig[cevt][cpop].p[i]);
						else	
							fprintf(out_inf,"%g,",P_mig[cevt][cpop].p[i]);
					}
				}
				//generalized gamma
				else if(P_mig[cevt][cpop].type == 2){
					P_mig[cevt][cpop].p = (double*)malloc(4*sizeof(double));
					if(cevt==0)
						fprintf(out_inf,"- Prior distribution (mig%d):    generalized gamma(",cpop+1);
					else
						fprintf(out_inf,"- Prior distribution (migAnc%d): generalized gamma(",cevt);
					for(i=0;i<4;++i){
						fscanf(inp_prs,"%lf",&P_mig[cevt][cpop].p[i]);
						if(i==3)
							fprintf(out_inf,"%g)\n",P_mig[cevt][cpop].p[i]);
						else	
							fprintf(out_inf,"%g,",P_mig[cevt][cpop].p[i]);
					}
				}
				//uniform on number of migrants
				else if(P_mig[cevt][cpop].type == 3){
					P_mig[cevt][cpop].p = (double*)malloc(2*sizeof(double));
					if(cevt==0)
						fprintf(out_inf,"- Prior distribution (mig%d):    uniform(",cpop+1);
					else
						fprintf(out_inf,"- Prior distribution (migAnc%d): uniform(",cevt);
					for(i=0;i<2;++i){
						fscanf(inp_prs,"%lf",&P_mig[cevt][cpop].p[i]);
						if(i==1)
							fprintf(out_inf,"%g) on emigrants number\n",P_mig[cevt][cpop].p[i]);
						else	
							fprintf(out_inf,"%g,",P_mig[cevt][cpop].p[i]);
					}
				}
				//generalized gamma on number of migrants
				else if(P_mig[cevt][cpop].type == 4){
					P_mig[cevt][cpop].p = (double*)malloc(4*sizeof(double));
					if(cevt==0)
						fprintf(out_inf,"- Prior distribution (mig%d):    generalized gamma(",cpop+1);
					else
						fprintf(out_inf,"- Prior distribution (migAnc%d): generalized gamma(",cevt);					
					for(i=0;i<4;++i){
						fscanf(inp_prs,"%lf",&P_mig[cevt][cpop].p[i]);
						if(i==3)
							fprintf(out_inf,"%g) on emigrants number\n",P_mig[cevt][cpop].p[i]);
						else	
							fprintf(out_inf,"%g,",P_mig[cevt][cpop].p[i]);
					}
				}
				else				 
					printerr("wrong migration prior type");
				//each anc pop is just done once
				if(cevt > 0)	
					break; 
				
			}
		}
	}
	else{
		fprintf(out_inf,"Time events:\n");
		fprintf(out_inf,"- No time events -\n");
		fprintf(out_inf,"Migration rates:\n");
		fprintf(out_inf,"- No migration -\n");
	}
	/*Reads the mutation rates (Microsatelites)*/
	fprintf(out_inf,"Mutation rates (Microsatelites):\n");
	fscanf(inp_prs,"%d",&P_mutM.type);
	//no mutation present
	if(P_mutM.type == 0){
		P_mutM.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- No mutation -\n");
	}
	//lognormal distribution
	else if(P_mutM.type == 1){
		P_mutM.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- HiperPrior distribution:     lognormal(");
		for(i=0;i<4;i++){
			fscanf(inp_prs,"%lf ",&P_mutM.p[i]);
			if(i==3)
				fprintf(out_inf,"%g)\n",P_mutM.p[i]);			
			else
				fprintf(out_inf,"%g,",P_mutM.p[i]);
		}
	}
    //normal distribution
    else if(P_mutM.type == 2){
        P_mutM.p = (double*)malloc(4*sizeof(double));
        fprintf(out_inf,"- HiperPrior distribution:     normal(");
        for(i=0;i<4;i++){
            fscanf(inp_prs,"%lf ",&P_mutM.p[i]);
            if(i==3)
                fprintf(out_inf,"%g)\n",P_mutM.p[i]);            
            else
                fprintf(out_inf,"%g,",P_mutM.p[i]);
        }
    }
	else{
		printerr("wrong mutation prior type (Microsatelite)");
	}
	/*Reads the mutation rates (Sequence data)*/
	fprintf(out_inf,"Mutation rates (Sequence data):\n");
	fscanf(inp_prs,"%d",&P_mutS.type);
	//no mutation present
	if(P_mutS.type == 0){
		P_mutS.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- No mutation -\n");
	}
	//lognormal distribution
	else if(P_mutS.type == 1){
		P_mutS.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- HiperPrior distribution:     lognormal(");
		for(i=0;i<4;i++){
			fscanf(inp_prs,"%lf ",&P_mutS.p[i]);
			if(i==3)
				fprintf(out_inf,"%g)\n",P_mutS.p[i]);			
			else
				fprintf(out_inf,"%g,",P_mutS.p[i]);
		}
	}
    //normal distribution
    else if(P_mutS.type == 2){
        P_mutS.p = (double*)malloc(4*sizeof(double));
        fprintf(out_inf,"- HiperPrior distribution:     normal(");
        for(i=0;i<4;i++){
            fscanf(inp_prs,"%lf ",&P_mutS.p[i]);
            if(i==3)
                fprintf(out_inf,"%g)\n",P_mutS.p[i]);            
            else
                fprintf(out_inf,"%g,",P_mutS.p[i]);
        }
    }
	else{
		printerr("wrong mutation prior type (Sequence data)");
	}
	/*Reads the recombination rates (Microsatelites)*/
	fprintf(out_inf,"Recombination rates (Microsatelites):\n");
	fscanf(inp_prs,"%d",&P_recM.type);
	//no recombination present
	if(P_recM.type == 0){
		P_recM.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- No recombination -\n");
	}
	//lognormal distribution
	else if(P_recM.type == 1){
		P_recM.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- HiperPrior distribution:     lognormal(");
		for(i=0;i<4;i++){
			fscanf(inp_prs,"%lf ",&P_recM.p[i]);
			if(i==3)
				fprintf(out_inf,"%g)\n",P_recM.p[i]);			
			else
				fprintf(out_inf,"%g,",P_recM.p[i]);
		}
	}
    //normal distribution
    else if(P_recM.type == 2){
        P_recM.p = (double*)malloc(4*sizeof(double));
        fprintf(out_inf,"- HiperPrior distribution:     normal(");
        for(i=0;i<4;i++){
            fscanf(inp_prs,"%lf ",&P_recM.p[i]);
            if(i==3)
                fprintf(out_inf,"%g)\n",P_recM.p[i]);            
            else
                fprintf(out_inf,"%g,",P_recM.p[i]);
        }
    }
	else{
		printerr("wrong recombination prior type (Microsatelite");
	}
	/*Reads the recombination rates (Sequence data)*/
	fprintf(out_inf,"Recombination rates (Sequence data):\n");
	fscanf(inp_prs,"%d",&P_recS.type);
	//no recombination present
	if(P_recS.type == 0){
		P_recS.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- No recombination -\n");
	}
	//lognormal distribution
	else if(P_recS.type == 1){
		P_recS.p = (double*)malloc(4*sizeof(double));
		fprintf(out_inf,"- HiperPrior distribution:     lognormal(");
		for(i=0;i<4;i++){
			fscanf(inp_prs,"%lf ",&P_recS.p[i]);
			if(i==3)
				fprintf(out_inf,"%g)\n",P_recS.p[i]);			
			else
				fprintf(out_inf,"%g,",P_recS.p[i]);
		}
	}
    //normal distribution
    else if(P_recS.type == 2){
        P_recS.p = (double*)malloc(4*sizeof(double));
        fprintf(out_inf,"- HiperPrior distribution:     normal(");
        for(i=0;i<4;i++){
            fscanf(inp_prs,"%lf ",&P_recS.p[i]);
            if(i==3)
                fprintf(out_inf,"%g)\n",P_recS.p[i]);            
            else
                fprintf(out_inf,"%g,",P_recS.p[i]);
        }
    }
	else{
		printerr("wrong recombination prior type (Sequence data)");
	}
	/*Reads the migration weights*/
	if(npop>2){
		fprintf(out_inf,"Migration weights matrix:\n");
		fscanf(inp_prs,"%d",&M_migw.type);
		//no use of migration weights
		if(M_migw.type == 0){
			fprintf(out_inf,"- No migration weights in use -\n");
		}
		//use of migration weights
		else if(M_migw.type == 1){
			if(P_top.type==0 || P_top.type==3)
				printerr("to use migration weights matrix, topology must be defined (option 1,2,4 or 5)");
			
			M_migw.m = (double***)malloc(npop*sizeof(double**));
			for(cpop=0; cpop<npop; cpop++){
				fprintf(out_inf,"- %d:\n", cpop);			
				M_migw.m[cpop] = (double**)malloc(nevt*sizeof(double*));
				for(cevt=0; cevt<nevt; cevt++){
					M_migw.m[cpop][cevt] = (double*)malloc(npop*sizeof(double));
					for(cpop2=0; cpop2<npop; cpop2++){
						fscanf(inp_prs,"%lf",&M_migw.m[cpop][cevt][cpop2]);
						if(cpop2==0)
							fprintf(out_inf,"  %6.4g ",M_migw.m[cpop][cevt][cpop2]);			
						else
							fprintf(out_inf,"%6.4g ",M_migw.m[cpop][cevt][cpop2]);			
					}
					fprintf(out_inf,"\n");			
				}
				fprintf(out_inf,"\n");			
			}
		}
		else
			printerr("wrong migration weights matrix related option");
	}
	
	fclose(inp_prs);

	/*counting the number of parameters*/
	nparams=0;
	if(npop>1){
		if(P_top.type==3||P_top.type==4||P_top.type==5)
			nparams++;
		if(foundSTR)
			nparams+=2;
		if(foundSNP)
			nparams+=2;
		if(foundSTR)
			nparams+=2;
		if(foundSNP)
			nparams+=2;
		nparams+= 1+(npop-1)+(2*npop-1)+(2*npop-2);
	}
	else{
		if(P_top.type==3)
			nparams++;
		if(foundSTR)
			nparams+=2;
		if(foundSNP)
			nparams+=2;
		if(foundSTR)
			nparams+=2;
		if(foundSNP)
			nparams+=2;
		nparams+= 1+(2*npop-1);
	}

	//check if there are too many nparams+nstats
	if((nparams+nstats)*niter>MAXDATA)
		printerr(".dat file would be too big to be analyzed");

	/*check if the STR's loci are linked*/
	if(foundSTR && P_recM.type != 0){
		//join all the STR loci
		nrSTR = 0;
		sampleSTR = (int*)malloc(nloc*sizeof(int));
		for(cloc=0,i=0;cloc<nloc;cloc++){
			if(pm.type[cloc]=='m'||pm.type[cloc]=='M'){
				if(nrSTR==0){
					ploidySTR = pm.ploidy[cloc];
					for(cpop=0; cpop<npop; cpop++)
						sampleSTR[cpop] = pm.nsamp[cloc][cpop];
				}
				else if(ploidySTR != pm.ploidy[cloc])
					printerr("when using rec with STR all loci must have the same scalar");
				else
					for(cpop=0;cpop<npop; cpop++)
						if(sampleSTR[cpop] != pm.nsamp[cloc][cpop])
							printerr("when using rec with STR all loci must have the same samples per population");	
				nrSTR++;
			}
			else{
				pm.type[i] = pm.type[cloc];
				pm.ploidy[i] = pm.ploidy[cloc];
				for(cpop=0;cpop<npop;cpop++){
					pm.nsamp[i][cpop] = pm.nsamp[cloc][cpop];
				}
				i++;
			}
		}
		for(cSTR=0; cSTR<nrSTR; cSTR++){
			pm.type[i+cSTR] = 'm';
			pm.ploidy[i+cSTR] = ploidySTR;
			for(cpop=0;cpop<npop;cpop++){
				pm.nsamp[i+cSTR][cpop] = sampleSTR[cpop];
			}
		}
	}
	else
		nrSTR = 1;
	pm.nrSTR=nrSTR;

	/*write second part of the .txt output file*/
	fprintf(out_inf,"\nOUTPUT FILES:\n");
	fprintf(out_inf,"-------------------------------\n");
	if(printIt)
		for(citer=0;citer<niter;citer++)
			fprintf(out_inf,"Simulation data file %d:     %sdata%d.len\n",citer+1,path,citer);
	if(printMut)
		fprintf(out_inf,"Mutation rates file:         %s.mut\n",argv[4]);
	if(printRec)
		fprintf(out_inf,"Recombination rates file:    %s.rec\n",argv[4]);
	fprintf(out_inf,"Data file:                   %s.dat\n",argv[4]);
	fprintf(out_inf,"%d parameters [top",nparams);
	if(P_top.type==3||P_top.type==4||P_top.type==5)
		fprintf(out_inf,"|marker");
	if(foundSTR)
		fprintf(out_inf,"|avMutM|sdMutM");
	if(foundSNP)
		fprintf(out_inf,"|avMutS|sdMutS");
	if(foundSTR)
		fprintf(out_inf,"|avRecM|sdRecM");
	if(foundSNP)
		fprintf(out_inf,"|avRecS|sdRecS");
	if(npop>1){
		for(cevt=0;cevt<nevt;++cevt)
			fprintf(out_inf,"|t%d",cevt+1);
		for(cevt=0;cevt<=nevt;++cevt){
			for(cpop=0;cpop<npop;cpop++){
				if(cevt==0)
					fprintf(out_inf,"|Ne%d",cpop+1);
				else
					fprintf(out_inf,"|NeA%d",cevt);
				if(cevt > 0)	
				break; 
			}
		}
		for(cevt=0;cevt<nevt;++cevt){
			for(cpop=0;cpop<npop;cpop++){
				if(cevt==0)
					fprintf(out_inf,"|mig%d",cpop+1);
				else
					fprintf(out_inf,"|migA%d",cevt);
				if(cevt > 0)	
				break; 
			}
		}
	}
	else
		fprintf(out_inf,"|Ne1");	
	fprintf(out_inf,"]\n%d summstats [",nstats);
	started=0;
	if(foundSTR)
		for(i=0;i<MAXSSTATS_M;i++){
			if(lsstats[i])
				for(cpop=0; cpop<npop; cpop++)
					if(started)
						fprintf(out_inf,"|%s%d",lsstats2[i],cpop+1);
					else{
						fprintf(out_inf,"%s%d",lsstats2[i],cpop+1);
						started=1;
					}
					
		}
	if(foundSTR)
		for(i=0;i<MAXSSTATS_M;i++){
			if((i==0||i==1||i==2||i==3||i==4||i==6||i==7||i==8||i==9) && lsstats[i])
				for(cpop = 0;cpop <npop;++cpop){
					for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
						if(started)
							fprintf(out_inf,"|%s%d-%d",lsstats2[i],cpop+1,cpop2+1);
						else{
							fprintf(out_inf,"%s%d-%d",lsstats2[i],cpop+1,cpop2+1);
							started=1;
						}
					
					}
				}
		}
	if(foundSNP)
		for(i=MAXSSTATS_M;i<MAXSSTATS;i++){
			if(lsstats[i])
				for(cpop=0; cpop<npop; cpop++)
					if(started)
						fprintf(out_inf,"|%s%d",lsstats2[i],cpop+1);
					else{
						fprintf(out_inf,"%s%d",lsstats2[i],cpop+1);
						started=1;
					}
					
		}
	if(foundSNP)
		for(i=MAXSSTATS_M;i<MAXSSTATS;i++){
			if((i==0||i==1||i==2||i==3||i==4||i==6||i==7||i==8||i==9) && lsstats[i])
				for(cpop = 0;cpop <npop;++cpop){
					for(cpop2 = cpop+1;cpop2 < npop;++cpop2){
						if(started)
							fprintf(out_inf,"|%s%d-%d",lsstats2[i],cpop+1,cpop2+1);
						else{
							fprintf(out_inf,"%s%d-%d",lsstats2[i],cpop+1,cpop2+1);
							started=1;
						}
					
					}
				}
		}
		
	fprintf(out_inf,"]\n\nRUNTIME:\n");
	fprintf(out_inf,"-------------------------------\n");
	fprintf(out_inf,"Starting date: %s",asctime(startTime));
	
	data.nsamp = pm.nsamp;
	data.npop = npop;
	data.nloc = nloc;
	data.tsamp = (int *)malloc(nloc*sizeof(int));		
	for(cloc=0; cloc<nloc; cloc++){													
		data.tsamp[cloc]=0;								
		for(cpop=0; cpop<npop; cpop++)					
			data.tsamp[cloc]+=data.nsamp[cloc][cpop];	
	}													
	data.Nmax = (int *)malloc(nloc*sizeof(int *));
	data.ldna = (int *)malloc(nloc*sizeof(int *));
	data.freq = (int ***)malloc(nloc*sizeof(int **));
	if(foundSTR){
		data.valM = (int **)malloc(nloc*sizeof(int *));
	}
	if(foundSNP){
		data.valS = (char ***)malloc(nloc*sizeof(char **));
		data.lsites = (int *)malloc(nloc*sizeof(int *));
	}
	for(cloc = 0;cloc<nloc;++cloc){
		data.Nmax[cloc] = 0;
		data.freq[cloc] = (int **)malloc(npop*sizeof(int *));
		for(cpop=0;cpop<npop;++cpop){
			data.freq[cloc][cpop] = (int *)malloc(data.tsamp[cloc]*sizeof(int));	
		}
		if(foundSTR)
			data.valM[cloc] = (int *)malloc(sizeof(int));
		if(foundSNP){
			data.lsites[cloc] = 0;
			data.valS[cloc] = (char **)malloc(data.tsamp[cloc]*sizeof(char *));
		}
	}	

    /*check if the mutation rates are going to be recorded*/
    if(printMut){
    	name_mut = (char *)malloc(outsize*sizeof(char));
		strcpy(name_mut,argv[4]);
		out_mut = fopen(strcat(name_mut,".mut"),"w");	//output .mut
		if(out_mut == NULL)
			printerr("cannot create .mut file");
		free(name_mut);
    }

    /*check if the recombination rates are going to be recorded*/
    if(printRec){
    	name_rec = (char *)malloc(outsize*sizeof(char));
		strcpy(name_rec,argv[4]);
		out_rec = fopen(strcat(name_rec,".rec"),"w");	//output .rec
		if(out_rec == NULL)
			printerr("cannot create .rec file");
		free(name_rec);
    }
	
	outline1 = (char*)malloc(nparams*15*sizeof(char));
	if(nstats!=0)
		outline2 = (char*)malloc(nstats*15*sizeof(char));
	else
		outline2 = (char*)malloc(15*sizeof(char));

	/*Generate all the genetic trees*/
	opengfsr(home);
	printf("\n--\nNumber of Genetic trees to be simulated: %d\n",niter);
	printf("\nStarting of ABC program...\n ");
	for(citer=0;citer<niter;++citer){
		sampPriors(&pm,outline1,out_mut,out_rec,printMut,printRec,foundSTR,foundSNP,pm.type);
		sampLikelihood(&pm,&data,printIt,citer,path);
		summStats(&data,lsstats,outline2,foundSTR,foundSNP,pm.type);
		freeMem(nloc,data.tsamp,data.valS,pm.type);
		fprintf(out_dat,"%s%s",outline1,outline2);
		if(citer!=0 && citer%RecIter==0){
			printf("\nGenetic trees simulated so far: %d\n",citer);
		}
	}
	printf("\n...simulations ending\n");
	printf("\nGenetic trees simulated: %d\n--\n",citer);
	closegfsr(home);
	
	/*get ending time and save it to .txt file*/
	time( &endClock );   				// Get time in seconds
	endTime = localtime( &endClock );  	// Convert time to struct tm form 
	fprintf(out_inf,"Ending date:   %s\n",asctime(endTime));
	fprintf(out_inf,"\nEND OF OUTPUT\n");

	/*free stuff*/
	if(foundSTR && P_recM.type != 0)
		free(sampleSTR);
	freetree(nloc - pm.nrSTR + 1,npop,nevt);
	free(home);
	free(path);
	if(npop>1){
		if(P_t[0].type==1||P_t[0].type==2){
			for(cevt=0;cevt<nevt;cevt++){
				free(pm.mig[cevt]);
				free(P_t[cevt].p);
			}
		}
		else{
			for(cevt=0;cevt<nevt;cevt++){
				free(pm.mig[cevt]);
			}
			free(P_t[0].p);
		}
		for(cevt=0;cevt<nevt;++cevt){
			for(cpop=0;cpop<npop;++cpop){
					free(P_mig[cevt][cpop].p);
				if(cevt>0)
					break;
			}
			free(P_mig[cevt]);
		}
		free(pm.mig);
		free(pm.tev);
		free(P_t);
		free(P_mig);
	}
	if(npop>2){
		if(M_migw.type==1){
			for(cpop=0; cpop<npop; cpop++){
				for(cevt=0; cevt<nevt; cevt++)
					free(M_migw.m[cpop][cevt]);
				free(M_migw.m[cpop]);
			}
			free(M_migw.m);
		}
	}
	free(P_mutM.p);
	free(P_mutS.p);
	free(P_recM.p);
	free(P_recS.p);
	free(pm.mu);
	free(pm.rec);
	free(pm.ploidy);
	free(pm.type);
	for(cloc=0;cloc<nloc;++cloc)
		free(pm.nsamp[cloc]);
	free(pm.nsamp);
	
	for(cevt=0;cevt<=nevt;++cevt){
		for(cpop=0;cpop<npop;++cpop){
			free(P_psize[cevt][cpop].p);
			if(cevt>0)
				break;
		}
		free(P_psize[cevt]);
		free(pm.psize[cevt]);
	}
	if(P_top.type!=0)
		free(P_top.p);
	free(P_psize);
	free(pm.psize);
	free(outline2);
	free(outline1);
	free(data.Nmax);
	free(data.ldna);
	free(data.tsamp);					
	for(cloc = 0;cloc<nloc;++cloc){
		for(cpop=0;cpop<npop;++cpop)
			free(data.freq[cloc][cpop]);
		free(data.freq[cloc]);
	}
	free(data.freq);
	if(foundSNP){
		for(cloc = 0;cloc<nloc;++cloc)
			free(data.valS[cloc]);
		free(data.valS);
		free(data.lsites);
	}
	if(foundSTR){
		for(cloc = 0;cloc<nloc;++cloc)
			free(data.valM[cloc]);
		free(data.valM);
	}
	/*close files*/
	fclose(out_dat);
	fclose(out_inf);	
	if(printMut)
		fclose(out_mut);
	if(printRec)
		fclose(out_rec);

}	//end of main

int notop(int n){
	int i,
		prod;

	prod=1;
	for(i=2;i<=n;i++)
		prod*=combinations(i,2);
	return prod;

}// end of notop

void freeMem(int nloc, int *tsamp,char ***val,char *ltype){
	int cloc,cdna;	//iterators
	
		for(cloc=0 ; cloc<nloc ; cloc++){
			if(ltype[cloc]=='s'||ltype[cloc]=='S'){
				for(cdna=0;cdna<tsamp[cloc];++cdna){
					if(val[cloc][cdna]!=NULL){
						free(val[cloc][cdna]);
						val[cloc][cdna] = NULL;
					}
				}
			}
		}

} //end of freeMem
