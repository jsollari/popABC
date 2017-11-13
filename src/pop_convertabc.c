/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		8th May 2009
	
	NBB - based on Mark's freqtab2.c and freqtabS.c
*/

#include "interface.h"

/*
	This function gets an input file with individuals data and calculates the frequence
	of their populations. Can be used with 1 or more overall populations.
	
	@arg input filename
	@arg output filename
*/	 
int createFreqTab(char *input,char *output){
	int cloc,cdna,cindiv,csites,cpop,	//iterators
		nindiv,							//number of individuals
		nloc,							//number of loci
		npop,							//total number of populations
		isum,							//check if the population as any haplotype information on a specific locus
		foundSTR,						//check if STR's are present
		foundSNP,						//check if SNP's are present
		outsize,						//size of output file name 
  		equalDna,						//check if the haplotype is diferent
		*pop = NULL,					//list of the belonging population of all individuals
		*ldna = NULL,					//number of diferent haplotypes per loci
		***alln;						//frequency of dna data by no_pop by nr_diff_dna by loci
	char c1,							//auxiliar to read from file
  		 *ltype,						//DNA type per loci
  		 **popc = NULL,					//list with each individual's population identifiers (A,B,...)
		 **lpop = NULL,					//list of population identifiers (A,B,...)
		 *outname;						//name of the output.len file
	FILE *inp,							//pntr to input file
	 	 *outp;							//pntr to output_length file
	time_t startClock;					//time when the program starts
	const struct tm *startTime;			//struct time when the program starts
/* only used in sequence dna analysis */
	int	nsites,							//number of sites in an allele of a loci
		*lsites;						//number of sites in the dna sequence per loci	 
	char ***valS,						//list of all different haplotypes ordered according to alln[][] by pop by loc
		 ***genotS = NULL;				//list of each individual's haplotype per loci, per allele
/* only used in microsatellites analysis */
	int	i1,								//auxiliar to read from file
		**IDsort,						//sorted valM[]
    	**valM,							//list of all different haplotypes ordered according to alln[][] by pop by loc
		**genotM = NULL;				//list of each individual's haplotype per loci, per allele

	inp = fopen(input,"r");
	if(inp == NULL)
		return 1; //cannot open .pop file

	outsize = strlen(output) + 5;
	outname = (char *)malloc(outsize*sizeof(char));
	strcpy(outname,output);
	outp = fopen(strcat(outname,".len"),"w");	
	if(outp == NULL)
		return 2; //cannot create .len file

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 


	/*if # (commentary) run to end of line*/
	c1 = getc(inp);
	while(c1 == '#'){
		while(!isendline(c1 = getc(inp)));
		while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
	}
	ungetc(c1,inp);

	fscanf(inp,"%d",&nloc);
	
 	ltype = (char*)malloc(1+nloc*sizeof(char));
	foundSTR = 0;
	foundSNP = 0;	
	for(cloc=0 ; cloc<nloc ; cloc++){
		fscanf(inp,"%s",&ltype[cloc]);
		if(ltype[cloc]=='m'||ltype[cloc]=='M')
			foundSTR = 1;
		if(ltype[cloc]=='s'||ltype[cloc]=='S')
			foundSNP = 1;
	}

	if(foundSNP){
        lsites = (int *)malloc(nloc*sizeof(int));
    }
    
    /*run through all the individuals*/
	npop = 0;	
	for(cindiv=0 ; ; cindiv++){
		while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');

		//1st for(;;) STOP condition: reached EOF
		if(c1== EOF)
			break;
		//if # (commentary) run to end of line
		while(c1 == '#'){
			while(!isendline(c1 = getc(inp)));
			while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
		}
		ungetc(c1,inp);
        
		pop = (int *)myAlloc(pop,(cindiv+1)*sizeof(int));	
		popc = (char**)myAlloc(popc,(cindiv+1)*sizeof(char*));
        popc[cindiv] = NULL;
        popc[cindiv] = (char *)myAlloc(popc[cindiv],MAXCHAR*sizeof(char));
        
        fscanf(inp,"%s",popc[cindiv]);
        
		while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t'); 
		if(!isdigit(c1))
			return 3; //individual has no information
		ungetc(c1,inp);
		
		//choose type of analysis
		if(foundSTR){
			genotM = (int **)myAlloc(genotM,(cindiv+1)*sizeof(int*));
			genotM[cindiv] = NULL;
		    genotM[cindiv] = (int *)malloc(nloc*sizeof(int));
        }
		if(foundSNP){
			genotS = (char ***)myAlloc(genotS,(cindiv+1)*sizeof(char**));
			genotS[cindiv] = NULL;
            genotS[cindiv] = (char **)malloc(nloc*sizeof(char*));
            for(cloc=0 ; cloc<nloc; cloc++){
                genotS[cindiv][cloc] = NULL;
            }
 		}
        			
		//run through all the loci of one individual
		for(cloc=0 ; cloc<nloc; cloc++){
			//choose type of analysis
			if(ltype[cloc]=='M' || ltype[cloc]=='m'){
				fscanf (inp, "%d", &i1);
				genotM[cindiv][cloc] = i1;								
			}
			else{
				//run through all the sites of a allele of a locus of an individual
				for(csites=0 ; isdigit(c1 = getc(inp)); csites++){
					genotS[cindiv][cloc] = (char *)myAlloc(genotS[cindiv][cloc],(csites+1)*sizeof(char));
					genotS[cindiv][cloc][csites] = c1;
				}
				//dealing with the first individual, check number of sites of current locus
				if(cindiv==0){
					lsites[cloc]=csites;
				}
				//not dealing with first individual, check if number of sites are the same
				else{
					if(lsites[cloc]!= csites)
						return 4; //different number of loci or number of sites in a loci
				}		
			}
			while(isspace(c1=getc(inp))||c1=='\t');
			ungetc(c1,inp);
		}	//end of for(cloc=0; ; )
		
		//get the number of populations, assign numeric identifiers to populations
		for(cpop=0 ; cpop<npop ; ++cpop){
			// check if the current alphabetic identifier is already listed
			if(strcmp(lpop[cpop],popc[cindiv])==0){
				pop[cindiv] = cpop+1;
				break;
			}
		}
		// the population identifier is not listed yet
		if(cpop == npop){
			lpop = (char**)myAlloc(lpop,(cpop+1)*sizeof(char*));
            lpop[cpop] = NULL;
            lpop[cpop] = (char*)myAlloc(lpop[cpop],MAXCHAR*sizeof(char));
			strcpy(lpop[cpop],popc[cindiv]); //add pop identifier to lpop[]
			pop[cindiv] = cpop+1;		     //add numeric population identifier to pop[]
			++npop;						     //increase population number by one
		}			
	}	//end of for(cindiv=0; ; )
	
	nindiv = cindiv;

	/*allocate memory to valM, IDsort, valS and alln*/
	if(foundSTR){
		valM = (int **)malloc(nloc*sizeof(int *));
		IDsort = (int **)malloc(nloc*sizeof(int *));	
	}
	if(foundSNP)
		valS = (char ***)malloc(nloc*sizeof(char **));
	for(cloc=0; cloc<nloc ;++cloc){
		if(foundSTR){
			valM[cloc] = NULL;
			IDsort[cloc] = NULL;	
		}
		if(foundSNP){
			valS[cloc] = NULL;
		}
	}
	alln = (int ***)malloc(npop*sizeof(int **));
	for(cpop=0;cpop<npop;++cpop){
		alln[cpop] = (int **)malloc(nloc*sizeof(int *));
		
		for(cloc=0; cloc<nloc ;++cloc){
			alln[cpop][cloc] = (int *)malloc(nindiv*sizeof(int));			
		}
	}	

	/*run through every loci*/
	for(cloc=0;cloc<nloc;++cloc){
		ldna = (int *)myAlloc(ldna, (cloc+1)*sizeof(int));
		ldna[cloc]=0;
		
		//run through every individual		
		for(cindiv=0;cindiv<nindiv;++cindiv){
			//choose type of analysis
			if(ltype[cloc]=='M' || ltype[cloc]=='m'){
				valM[cloc] = (int *)myAlloc(valM[cloc], (ldna[cloc]+1)*sizeof(int));	
				IDsort[cloc] = (int *)myAlloc(IDsort[cloc], (ldna[cloc]+1)*sizeof(int));
						
				//check if there is information in the current individual			
				if(genotM[cindiv][cloc] == 0){
					continue;
				}
			}
			else{
				valS[cloc] = (char **)myAlloc(valS[cloc], (ldna[cloc]+1)*sizeof(char *));	
			
				//check if there is information in the current individual			
				if(genotS[cindiv][cloc][0] == '9'){
					continue;
				}
			}
				
			//check if the current haplotypeis already listed
			for(cdna=0;cdna<ldna[cloc];++cdna){
				equalDna = 1;						
					
				//choose type of analysis
				if(ltype[cloc]=='M' || ltype[cloc]=='m'){
					if(genotM[cindiv][cloc] != valM[cloc][cdna])
						equalDna = 0;
				}
				else{
					for(csites=0 ; csites<lsites[cloc] ; csites++){
						if(genotS[cindiv][cloc][csites] != valS[cloc][cdna][csites]){
							equalDna = 0;
							break;
						}
					}
				}
				if(equalDna){
					++alln[pop[cindiv]-1][cloc][cdna];	//increase freq of current haplotype
					break;
				}
			}
					
			//current haplotype is not listed yet
			if(cdna>=ldna[cloc]){
				for(cpop=0;cpop<npop;++cpop){
					//get current haplotype's population identifier 
					if(cpop == pop[cindiv]-1)
						continue;
					alln[cpop][cloc][ldna[cloc]] = 0;		//initiate as 0 the freq of current haplotype in other populations
				}
				alln[pop[cindiv]-1][cloc][ldna[cloc]] = 1;	//initiate as 1 the freq of current haplotype in its population
					
				//choose type of analysis
				if(ltype[cloc]=='M' || ltype[cloc]=='m'){
					valM[cloc][cdna] = genotM[cindiv][cloc];		//adds current haplotype to valM[][]
				}
				else{
					valS[cloc][cdna] = (char *)malloc((lsites[cloc]+1)*sizeof(char));
					for(csites=0; csites<lsites[cloc] ; csites++){
						valS[cloc][cdna][csites] = genotS[cindiv][cloc][csites];		//adds current haplotype to valS[][][]	
					}
					valS[cloc][cdna][csites] = '\0';
				}
				++ldna[cloc];	//increase by one the number of different haplotypes of a particular locus
			}
		}
	}
	
	/*write to file (analysis of sequence data)*/
	fprintf(outp,"# PopABC file converted to .len file with pop2table1.0\n");
	fprintf(outp,"# input file:  %s\n",input);
	fprintf(outp,"# output file: %s.len\n",output);
	fprintf(outp,"# date: %s\n",asctime(startTime));
	for(cpop=0;cpop<npop;++cpop){
		fprintf(outp,"#Population %d - pop%d\n",cpop+1,cpop+1);
	}
	for(cloc=0;cloc<nloc;++cloc){
		fprintf(outp,"#Locus %d - loc%d\n",cloc+1,cloc+1);
	}	
    fprintf(outp,"\n%d\n%d\n",npop,nloc);										//outp: npop, nloc
	
	for(cloc=0 ; cloc<nloc ; cloc++){
		fprintf(outp,"%c ",ltype[cloc]);
	}
	fprintf(outp,"\n\n",ltype[cloc]);											//outp: ltype[cloc]
	
	for(cloc=0 ; cloc<nloc ; cloc++){
		if(ltype[cloc]=='M' || ltype[cloc]=='m'){
 			isorti('a',ldna[cloc],valM[cloc],IDsort[cloc]);						//defined in myutil.c
 		}
		fprintf(outp,"%d\n",ldna[cloc]);										//outp: ndna
		
		for(cpop=0;cpop<npop;++cpop){
			for(cdna=0;cdna<ldna[cloc];++cdna){
				//choose type of analysis
				if(ltype[cloc]=='M' || ltype[cloc]=='m')
					fprintf(outp,"%4d ",alln[cpop][cloc][IDsort[cloc][cdna]]);	//(M)outp: alln[][]
				else
					fprintf(outp,"%4d ",alln[cpop][cloc][cdna]);				//(S)outp: alln[][]
			}
			fprintf(outp,"\n");
		}
		fprintf(outp,"\n");
		
		if(ltype[cloc]=='S' || ltype[cloc]=='s'){
			for(cdna=0;cdna<ldna[cloc];++cdna)
				fprintf(outp,"%4d ",cdna);										//(S)outp: cdna
			fprintf(outp,"\n\n");
			fprintf(outp,"%d \n",lsites[cloc]);									//(S)outp: lsites[]
		}	
		for(cdna=0;cdna<ldna[cloc];++cdna){
			//choose type of analysis
			if(ltype[cloc]=='M' || ltype[cloc]=='m')
				fprintf(outp,"%4d ",valM[cloc][IDsort[cloc][cdna]]);			//(M)outp: val[]
			else
				fprintf(outp,"%d %s\n",cdna,valS[cloc][cdna]);					//(S)outp: cdna, val[]
		}
		if(ltype[cloc]=='M' || ltype[cloc]=='m')
			fprintf(outp,"\n\n");
		else	
			fprintf(outp,"\n");
	}
	
	/*free stuff*/
	for(cindiv=0 ; cindiv<nindiv ; cindiv++){
		for(cloc=0 ; cloc<nloc ; cloc++){
			if(ltype[cloc]=='s'||ltype[cloc]=='S')
				free(genotS[cindiv][cloc]);
		}	
		if(foundSTR)
			free(genotM[cindiv]);
		if(foundSNP)
			free(genotS[cindiv]);
	}	
	for(cpop=0 ; cpop<npop ; cpop++){
		for(cloc=0; cloc<nloc ; cloc++)
			free(alln[cpop][cloc]);
		free(alln[cpop]);
	}
	for(cloc=0 ; cloc<nloc ; cloc++){
		if(ltype[cloc]=='s'||ltype[cloc]=='S')
			for(cdna=0 ; cdna<ldna[cloc] ; cdna++)
				free(valS[cloc][cdna]);
		if(ltype[cloc]=='m'||ltype[cloc]=='M'){
			free(IDsort[cloc]);
    		free(valM[cloc]);
    	}
		else
			free(valS[cloc]);
	}	
	if(foundSTR){
		free(valM);
		free(genotM);
		free(IDsort);		
	}
	if(foundSNP){
    	free(lsites);
		free(valS);
		free(genotS);
	}
    for(cindiv=0 ; cindiv<nindiv ; cindiv++){
        free(popc[cindiv]);
    }
    for(cpop=0 ; cpop<npop ; cpop++){
        free(lpop[cpop]);
    }    
	free(ldna);
	free(alln);
	free(ltype);	
	free(pop);
	free(popc);
	free(lpop);
	free(outname);
    	
	/*close files	*/
	fclose(inp);
	fclose(outp);

	return 0;
	
}	//end of main
