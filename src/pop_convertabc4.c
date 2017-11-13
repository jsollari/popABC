/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		8th May 2009

*/

#include "interface.h"

/*
	It creates an ABC .len file given a Nexus file (.nex)
	
	@arg length format filename
	@arg output filename
*/
int createFreqTab4(char *input,char *output){
	int cloc, cpop, cdna, csamp, csite, i, j,	//iterators
		nr1,					//retrieve info from input file
		nr2,					//retrieve info from input file
		reachEOF,				//check if EOF has been reached
		endSets,				//check if end of Sets has been reached
		endPop,					//check if end of pop has been reached
		firstInd,				//check if using '-' to define populations
		equalDna,				//auxiliar to build the .len file
		outsize,				//size of output file name 
  		npop,					//number of populations
  		nsites,					//number of bases
  		ndna,					//number of different dna data
		nsamp,					//total number of samples
		*lsamp,					//size of sample per pop	
		*id,					//identify the belonging population of each used individual
		 **freq;				//freq of dna data by no_pop by no. of diff dna data
	char c1,					//gets the value of a char temporarly
		 info1[MAXCHAR],		//retrieve info from the input file
		 info2[MAXCHAR],		//retrieve info from the input file
		 aux[MAXCHAR],			//auxiliar
		 *outname,				//name of the output file 
		 *mainseq,				//main sequence of a loci
		 **lpop=NULL,			//list of the populations names
		 **indiv,				//list of samples
		 **valS,				//list of all the diff dna sequencies
		 ***genotS;				//list of all the individual's haplotype per pop per loci 
	FILE *inp,					//pntr to the input file
		 *outp;					//pntr to the output file
	time_t startClock;			//time when the program starts
	const struct tm *startTime;	//struct time when the program starts

	inp = fopen(input,"r");
	if(inp == NULL)
		return 1; //cannot open .nex file
	
	outsize = strlen(output) + 5;
	outname = (char *)malloc(outsize*sizeof(char));
	strcpy(outname,output);
	outp = fopen(strcat(outname,".len"),"w");
	if(outp == NULL)
		return 2; //cannot create output file

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 

	fgets(aux,MAXCHAR,inp);					//#NEXUS
	
	reachEOF=0;
	while(!reachEOF){
		//reads input file
		while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
		if(c1 == '[')
			while((c1 = getc(inp))!=']');
		else if(c1==EOF)
			reachEOF=1;
		else
			ungetc(c1,inp);
		if(!reachEOF){
			while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
			ungetc(c1,inp);
			fscanf(inp,"%s",aux);
			if(strcmp(aux,"begin")==0||strcmp(aux,"BEGIN")==0||strcmp(aux,"Begin")==0){
			   	fscanf(inp,"%s",aux);
				if(strcmp(aux,"data;")==0||strcmp(aux,"DATA;")==0||strcmp(aux,"Data;")==0){
				   	fscanf(inp,"%s",aux);		//DIMENSIONS
				   	fscanf(inp,"%s",aux);
					i=0; j=0;
				   	while(!isdigit(aux[i]))
				   		i++;
				   	while(isdigit(aux[i+j])){
				   		info1[j]=aux[i+j];
				   		j++;
				   	}
				   	nsamp=atoi(info1);
				   	fscanf(inp,"%s",aux);
					i=0; j=0;
				   	while(!isdigit(aux[i]))
				   		i++;
				   	while(isdigit(aux[i+j])){
				   		info1[j]=aux[i+j];
				   		j++;
				   	}
				   	nsites=atoi(info1);
				   	fscanf(inp,"%s",aux);		//FORMAT
				   	fscanf(inp,"%s",aux);
					i=0; j=0;
				   	while(aux[i]!='=')
				   		i++;
				   	while(isalpha(aux[1+i+j])){
				   		info1[j]=aux[1+i+j];
				   		j++;
				   	}
				   	info1[j]='\0';
				   	if(!(strcmp(info1,"DNA")==0||strcmp(info1,"dna")==0||strcmp(info1,"Dna")==0))
				   		return 3; //program can't deal with specified data
					
					//MATRIX: go through all the individuals
				   	fscanf(inp,"%s",aux);		//MATRIX
					indiv=(char**)malloc(nsamp*sizeof(char*));
				   	for(csamp=0; csamp<nsamp; csamp++)
						indiv[csamp]=(char*)malloc((nsites+1)*sizeof(char));
					csamp=0;
					while(csamp<nsamp){
						while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
						if(c1 == '[')
							while((c1 = getc(inp))!=']');
						else	
							ungetc(c1,inp);
						while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
						if(c1 == '[')
							ungetc(c1,inp);
						else{
							ungetc(c1,inp);
							fscanf(inp,"%s",aux);
							while(isspace(c1=getc(inp))||c1=='\t');
							ungetc(c1,inp);
							csite=0;
							while(csite<nsites){
								if((c1=getc(inp))=='T'||c1=='A'||c1=='C'||c1=='G'||
								               c1=='t'||c1=='a'||c1=='c'||c1=='g'||c1=='N'||c1=='-'){
									indiv[csamp][csite]=c1;
									csite++;
								}
								else if(isspace(c1)||isendline(c1)||c1=='\t'){
									while(isspace(c1=getc(inp))||c1=='\t');
									ungetc(c1,inp);
								}
								else
									return 4; //problem reading MATRIX
							}
							csamp++;
						}
					}
					while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
					ungetc(c1,inp);
					fscanf(inp,"%s",aux);		//;
					fscanf(inp,"%s",aux);		//END;	
				}
				else if(strcmp(aux,"sets;")==0||strcmp(aux,"SETS;")==0||strcmp(aux,"Sets;")==0){
				   	id=(int*)malloc(nsamp*sizeof(int));
				   	for(csamp=0; csamp<nsamp; csamp++)
						id[csamp]=-1;
					//SETS: go through all the populations
				    cpop=0;
				   	endSets=0;
				   	while(!endSets){
				   		fscanf(inp,"%s",aux);
				   		if(strcmp(aux,"end;")==0||strcmp(aux,"END;")==0||strcmp(aux,"End;")==0){
				   			endSets=1;
				   		}
				   		else if(strcmp(aux,"TAXSET")==0||strcmp(aux,"taxset")==0||strcmp(aux,"Taxset")==0){
							lpop=(char**)myAlloc(lpop,(cpop+1)*sizeof(char*));
					   		lpop[cpop]=(char*)malloc(MAXCHAR*sizeof(char));
					   		fscanf(inp,"%s",aux);
							strcpy(lpop[cpop],aux);
							fscanf(inp,"%s",aux);	//=
							//get individuals belonging to pop
							endPop=0;
							while(!endPop){
								fscanf(inp,"%s",aux);
								i=0;
								firstInd=1;
								while(i<strlen(aux)){
									if(aux[i]==';'){
										endPop=1;
										i++;
									}							
									else if(aux[i]=='-'){
										firstInd=0;
										i++;
									}							
									else if(isdigit(aux[i])&&firstInd){
										j=0;
										while(isdigit(aux[i])){
											info1[j]=aux[i];
											i++;
											j++;
										}
										info1[j]='\0';
										nr1=atoi(info1);
										
									}
									else if(isdigit(aux[i])&&(!firstInd)){
										j=0;
										while(isdigit(aux[i])){
											info2[j]=aux[i];
											i++;
											j++;
										}
										info2[j]='\0';
										nr2=atoi(info2);
									}
									else
										return 5; //problem reading TAXSET
								}
								if(firstInd)
									id[nr1-1]=cpop;
								else{
									if(nr1>=nr2)
										return 6; //problem defining TAXSET: n1>=n2
									for(i=nr1;i<=nr2;i++)
										id[i-1]=cpop;
								}
							}
					   		cpop++;
				   		}
				   		else
				   			return 7; //unrecognized specifier after BEGIN SETS (use TAXSET)
				   	}
				   	npop=cpop;
				}
				else
					return 8; //unrecognized specifier after BEGIN (use DATA or SETS)
			}
			else if(aux[0]=='['){
				for(i=strlen(aux)-1; i>=0; i--)
					ungetc(aux[i],inp);
			}
			else
				return 9; //unrecognized command (only BEGIN accepted)
		}
	}
	
	
	//fill lsamp and genotS
	genotS = (char***)malloc(npop*sizeof(char**));
	lsamp = (int *)malloc(npop*sizeof(int));
	for(cpop=0; cpop<npop; cpop++){
		lsamp[cpop]=0;
		genotS[cpop] = (char**)malloc(nsamp*sizeof(char*));
		for(csamp=0; csamp<nsamp; csamp++){
			genotS[cpop][csamp] = (char *)malloc(nsites*sizeof(char));
		}
	}
	for(csamp=0; csamp<nsamp; csamp++){
		if(id[csamp]!=-1){
			for(csite=0; csite<nsites; csite++){
				genotS[id[csamp]][lsamp[id[csamp]]][csite] = indiv[csamp][csite];
			}
			lsamp[id[csamp]]++;
		}
	}
		
	//allocate memory to freq, IDsort, valS and ldna
	freq = (int **)malloc(npop*sizeof(int *));
	valS = (char **)malloc(nsamp*sizeof(char *));
	for(cpop=0;cpop<npop;++cpop){
		freq[cpop] = (int *)malloc(nsamp*sizeof(int));
	}	

	//fill freq[cpop][cdna] and valS[cdna][csite]
	ndna=0;
	for(cpop=0;cpop<npop;cpop++){
		//run through every individual		
		for(csamp=0;csamp<lsamp[cpop];++csamp){
			//check if the current haplotype is already listed
			for(cdna=0;cdna<ndna;++cdna){
				equalDna = 1;						
				for(csite=0 ; csite<nsites ; csite++){
					if(genotS[cpop][csamp][csite] != valS[cdna][csite]){
						if(valS[cdna][csite]=='N'||valS[cdna][csite]=='-'){
							valS[cdna][csite] = genotS[cpop][csamp][csite];
						}
						else if(genotS[cpop][csamp][csite]=='A'||genotS[cpop][csamp][csite]=='a'||
								genotS[cpop][csamp][csite]=='T'||genotS[cpop][csamp][csite]=='t'||
						  		genotS[cpop][csamp][csite]=='G'||genotS[cpop][csamp][csite]=='g'||
						  		genotS[cpop][csamp][csite]=='C'||genotS[cpop][csamp][csite]=='c'){
						  	equalDna = 0;
							break;
						}
					}
				}
				if(equalDna){
					++freq[cpop][cdna];	//increase freq of current haplotype
					break;
				}
			}
			//current haplotype is not listed yet
			if(cdna>=ndna){
				for(j=0;j<npop;++j){
					freq[j][ndna] = 0;//initiate as 0 the freq of current haplotype in other populations
				}
				freq[cpop][ndna] = 1;	//initiate as 1 the freq of current haplotype in its population
				
				valS[cdna] = (char *)malloc((nsites+1)*sizeof(char));
				for(csite=0; csite<nsites ; csite++){
					valS[cdna][csite] = genotS[cpop][csamp][csite];		//adds current haplotype to valS[][][]	
				}
				valS[cdna][csite] = '\0';
				++ndna;	//increase by one the number of different haplotypes of a particular locus
			}
		}
	}
	//transforms dna sequences into an array of segregating sites	
	mainseq = (char *)malloc(nsites*sizeof(char));
    for(csite=0; csite<nsites; csite++)
		mainseq[csite] = valS[0][csite];
	for(cdna=0;cdna<ndna;++cdna){
		for(csite=0;csite<nsites;++csite){
			if(valS[cdna][csite]!=mainseq[csite]&&(valS[cdna][csite]=='A'||
												   valS[cdna][csite]=='a'||
												   valS[cdna][csite]=='G'||
												   valS[cdna][csite]=='g'||
												   valS[cdna][csite]=='T'||
												   valS[cdna][csite]=='t'||
												   valS[cdna][csite]=='C'||
												   valS[cdna][csite]=='c')){
				valS[cdna][csite] = '1';
			}
			else{
				valS[cdna][csite] = '0';
			}
		}
	}

	//writes the output file
	fprintf(outp,"# Nexus file converted to .len file with nexus2table1.0\n");
	fprintf(outp,"# input file:  %s\n",input);
	fprintf(outp,"# output file: %s.len\n",output);
	fprintf(outp,"# date: %s#\n",asctime(startTime));
	for(cpop=0;cpop<npop;++cpop){
		fprintf(outp,"#Population %d - %s\n",cpop+1,lpop[cpop]);
	}
	fprintf(outp,"#Locus 1 - locus1\n");
	fprintf(outp,"\n%d\n",npop);												//outp: npop
	fprintf(outp,"1\ns\n\n");	 												//outp: nloc, ltype		
	fprintf(outp,"%d\n",ndna);													//outp: ndna
	for(cpop=0;cpop<npop;++cpop){
		for(cdna=0;cdna<ndna;++cdna){
			fprintf(outp,"%4d ",freq[cpop][cdna]);								//(S)outp: freq[][]
		}
		fprintf(outp,"\n");
	}
	fprintf(outp,"\n");
	for(cdna=0;cdna<ndna;++cdna)
		fprintf(outp,"%4d ",cdna);												//(S)outp: cdna
	fprintf(outp,"\n\n");
	fprintf(outp,"%d \n",nsites);												//(S)outp: nsites
	
	for(cdna=0;cdna<ndna;++cdna)
		fprintf(outp,"%d %s\n",cdna,valS[cdna]);								//(S)outp: cdna, valS[]
	fprintf(outp,"\n");
	
	//free stuff
	for(csamp=0 ; csamp<nsamp ; csamp++){
		free(indiv[csamp]);	
	}
	free(indiv);
	for(cpop=0 ; cpop<npop ; cpop++){
		for(csamp=0 ; csamp<nsamp ; csamp++){
			free(genotS[cpop][csamp]);	
		}
		free(genotS[cpop]);
	}	
	free(genotS);
	for(cpop=0; cpop<npop ; cpop++){
		free(freq[cpop]);
		free(lpop[cpop]);
	}
	free(freq);
	free(lpop);
	for(cdna=0 ; cdna<ndna; cdna++)
		free(valS[cdna]);
	free(valS);
	free(mainseq);
	free(lsamp);
	free(outname);
	free(id);
		
	//close files
	fclose(inp);
	fclose(outp);

	return 0;

} //end of main
