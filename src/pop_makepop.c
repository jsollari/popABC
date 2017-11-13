/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		12th May 2009
*/
#include "interface.h"

int makepop(char *input, char *output){
	int cpop,cloc,cdna,cindiv,csites,i,	//iterators
		foundSTR,						//check if STR's are present
		foundSNP,						//check if SNP's are present
		npop,							//number of population
		outsize,						//size of output file name 
		nloc,							//number of loci
		sum,							//sum of all individuals of a sample population
		*maxsamp,						//maximum number of sample in all loci of all the populations
		*ldna,							//number of diferent dna data
		***freq;						//freq of dna data by pop by diff_dna_data
	char c1,							//auxiliar to read file
		 *outname,						//name of the output.len file
		 *ltype;						//DNA type per loci
	FILE *inp_len,						//pntr to the input file
		 *outp;							//pntr to the output file
	time_t startClock;					//time when the program starts
	const struct tm *startTime;			//struct time when the program starts
/* only used in microsatellites analysis */		 
	int **valM,							//list of all the different dna sequence by loci
		***indivM;						//matrix with all the individuals by pop, by max_samp, by loci
/* only used in sequence dna analysis */
	char ***valS,						//list of all the different microsatellite size by loci
		 ****indivS;					//matrix with all the individuals by pop, by max_samp, be loci
	int dump,							//dump to put unwanted elements of input file
		*lsites;						//lsites of each sequency

	inp_len = fopen(input,"r");
	if(inp_len==NULL)
		return 1; //error: cannot open the input file

	outsize = strlen(output) + 5;
	outname = (char *)malloc(outsize*sizeof(char));
	strcpy(outname,output);
	outp = fopen(strcat(outname,".pop"),"w");	
	if(outp == NULL)
		return 2; //error: cannot create .pop file

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 
	
   	c1=getc(inp_len);
	while( c1 =='#'){
		while(!isendline(c1=getc(inp_len)));
		ungetc(c1,inp_len);
		while(isendline(c1=getc(inp_len)));
	}	
	ungetc(c1,inp_len);

	fscanf(inp_len,"%d",&npop);
	fscanf(inp_len,"%d",&nloc);
	
	ltype = (char*)malloc(1+nloc*sizeof(char));
	foundSTR = 0;
	foundSNP = 0;	
	for(cloc=0 ; cloc<nloc ; cloc++){
		fscanf(inp_len,"%s",&ltype[cloc]);
		if(ltype[cloc]=='m'||ltype[cloc]=='M')
			foundSTR = 1;
		if(ltype[cloc]=='s'||ltype[cloc]=='S')
			foundSNP = 1;
	}

	ldna = (int *)malloc(nloc*sizeof(int));
	maxsamp = (int *)malloc(npop*sizeof(int));
	freq = (int ***)malloc(nloc*sizeof(int **));
	
	if(foundSTR)
		valM = (int **)malloc(nloc*sizeof(int *));
	if(foundSNP){
		valS = (char ***)malloc(nloc*sizeof(char **));
		lsites = (int *)malloc(nloc*sizeof(int));
	}
	
	for(cpop=0;cpop<npop;++cpop)
		maxsamp [cpop]= 0;
	//fills freq and val	
	for(cloc=0;cloc<nloc;++cloc){
		fscanf(inp_len,"%d",&ldna[cloc]);
		freq[cloc] = (int **)malloc(npop*sizeof(int *));
		
		if(ltype[cloc] == 'm' || ltype[cloc] == 'M')
			valM[cloc] = (int *)malloc(ldna[cloc]*sizeof(int));
		else
			valS[cloc] = (char **)malloc(ldna[cloc]*sizeof(char*));
			
		for(cpop=0;cpop<npop;++cpop){
			freq[cloc][cpop] = (int *)malloc(ldna[cloc]*sizeof(int));
			sum = 0;
			for(cdna=0;cdna<ldna[cloc];++cdna){
				fscanf(inp_len,"%d",&(freq[cloc][cpop][cdna]));
				sum += freq[cloc][cpop][cdna];
			}
			if (maxsamp[cpop] < sum)
				maxsamp[cpop] = sum;
		}
				
		if(ltype[cloc] == 'm' || ltype[cloc] == 'M'){
			for(cdna=0;cdna<ldna[cloc];++cdna)
				fscanf(inp_len,"%d",&(valM[cloc][cdna]));
		}
		else{
			for(cdna=0;cdna<ldna[cloc];++cdna){
				fscanf(inp_len,"%d",&dump);
			}
			fscanf(inp_len,"%d",&lsites[cloc]);
			for(cdna=0;cdna<ldna[cloc];++cdna){
				valS[cloc][cdna] = (char *)malloc(lsites[cloc]*sizeof(char));
				fscanf(inp_len,"%d ",&dump);
				for(csites=0;csites<lsites[cloc];csites++){
					fscanf(inp_len,"%c",&(valS[cloc][cdna][csites]));
				}
			}
		}
	}
	
	//allocate memory to indiv
	if(foundSTR){ 	
		indivM = (int ***)malloc(npop*sizeof(int **));
		for(cpop=0 ; cpop<npop ; ++cpop){
			indivM[cpop]= (int **)malloc(maxsamp[cpop]*sizeof(int *));
			for(cindiv=0;cindiv<maxsamp[cpop];cindiv++){
				indivM[cpop][cindiv] = (int *)malloc(nloc*sizeof(int));
				for(cloc=0;cloc<nloc;cloc++){
					indivM[cpop][cindiv][cloc] = 0;				
				}	
			}
		}
	}
	if(foundSNP){
		indivS = (char ****)malloc(npop*sizeof(char ***));
		for(cpop=0;cpop<npop;++cpop){
			indivS[cpop]= (char ***)malloc(maxsamp[cpop]*sizeof(char **));
			for(cindiv=0;cindiv<maxsamp[cpop];cindiv++){
				indivS[cpop][cindiv] = (char **)malloc(nloc*sizeof(char *));
				for(cloc=0;cloc<nloc;cloc++){
					if(ltype[cloc]=='s'||ltype[cloc]=='S'){
						indivS[cpop][cindiv][cloc] = (char *)malloc(lsites[cloc]*sizeof(char));
						for(csites=0;csites<lsites[cloc];csites++){
							indivS[cpop][cindiv][cloc][csites] = '9';
						}
					}
				}	
			}
		}
	}
	
	//fills indiv
	for(cloc=0;cloc<nloc;cloc++){
		for(cpop=0;cpop<npop;++cpop){
			cindiv=0;
			for(cdna=0;cdna<ldna[cloc];++cdna){
				for(i=0;i<freq[cloc][cpop][cdna];i++){
					if(ltype[cloc] == 'm' || ltype[cloc] == 'M')
						indivM[cpop][cindiv][cloc] = valM[cloc][cdna];
					else
						for(csites=0;csites<lsites[cloc];csites++)
							indivS[cpop][cindiv][cloc][csites] = valS[cloc][cdna][csites];
					cindiv++;
				}			
			}
		}
	}
	
	//writes the output file
	fprintf(outp,"# Table file converted to .pop file with table2pop1.0\n");
	fprintf(outp,"# input file:  %s\n",input);
	fprintf(outp,"# output file: %s.pop\n",output);
	fprintf(outp,"# date: %s\n",asctime(startTime));
	fprintf(outp,"%d\n",nloc);													//outp:nloc
	for(cloc=0; cloc<nloc; cloc++)
		fprintf(outp,"%c ",ltype[cloc]);										//outp:ltype[cloc]
	fprintf(outp,"\n\n");
	for(cpop=0;cpop<npop;++cpop){
		for(cindiv=0;cindiv < maxsamp[cpop];++cindiv){
			fprintf(outp,"pop%d ",cpop+1);										//outp:pop
			for(cloc=0;cloc<nloc;++cloc){
				if(ltype[cloc] == 'm' || ltype[cloc] == 'M')
					fprintf(outp,"%d",indivM[cpop][cindiv][cloc]);				//outp: indivM[][]
				else
					for(csites=0;csites<lsites[cloc];csites++)
						fprintf(outp,"%c",indivS[cpop][cindiv][cloc][csites]);	//outp: indivS[][]
				fprintf(outp,"  ");
			}
		fprintf(outp,"\n");
		}	
	}

	//frees allocated memory
	if(foundSTR){
		for(cloc=0 ; cloc<nloc ; cloc++){
			if(ltype[cloc]=='m'||ltype[cloc]=='M')
				free(valM[cloc]);	
		}
		free(valM);
		for(cpop=0 ; cpop<npop ; cpop++){
			for(cindiv=0 ; cindiv<maxsamp[cpop] ; cindiv++)
				free(indivM[cpop][cindiv]);
			free(indivM[cpop]);
		}
		free(indivM);
	}
	if(foundSNP){
		for(cloc=0 ; cloc<nloc ; cloc++){
			if(ltype[cloc]=='s'||ltype[cloc]=='S'){
				for(cdna=0 ; cdna<ldna[cloc] ; cdna++)
					free(valS[cloc][cdna]);
				free(valS[cloc]);
			}
		}
		free(valS);
		for(cpop=0 ; cpop<npop ; cpop++){
			for(cindiv=0 ; cindiv<maxsamp[cpop] ; cindiv++){
				for(cloc=0 ; cloc<nloc ; cloc++){
					if(ltype[cloc]=='s'||ltype[cloc]=='S')
						free(indivS[cpop][cindiv][cloc]);
				}
				free(indivS[cpop][cindiv]);
			}
			free(indivS[cpop]);
		}
		free(indivS);	
		free(lsites);
	}
	for(cloc=0 ; cloc<nloc ; cloc++){
		for(cpop=0 ; cpop<npop ; cpop++)
			free(freq[cloc][cpop]);
		free(freq[cloc]);
	}
	free(freq);
	free(maxsamp);
	free(ldna);
	free(ltype);
	free(outname);
	
	//closes the files opened
	fclose(inp_len);
	fclose(outp);
	
	return 0;
			
}	//end of main
