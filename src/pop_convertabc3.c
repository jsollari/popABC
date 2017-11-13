/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		8th May 2009

*/

#include "interface.h"

/*
	It creates an ABC .len file given an GenePop input file.
	
	@arg length format filename
	@arg output filename
*/
int createFreqTab3(char *input,char *output){
	int cloc, cpop, cdna, csamp, i, j, cplo,	//iterators
		nsamp,					//total nr. of samples
		loctype,				//0 = locus divided by ','; 1 = locus divided by lines
		searchInd,				//check if ',' is found in the sample file
		reachEND,				//check if EOF have been reached
		searchPop,				//check if the word 'pop' is found in the sample file
		npop,					//number of populations
		nloc,					//number of loci
		equalDna,				//auxiliar to build the .len file
		outsize,				//size of output file name 
  		*ldna,					//number of different dna data by loci
		*lplo,					//nr of markers per loci
		*lsamp = NULL,			//total number of samples per loci
		***freq,				//freq of dna data by no_pop by no. of diff dna data
		**valM,					//list of all the diff microsatellite sizes
		**IDsort,				//list of the order of the diff microsatellite sizes
		***genotM2,				//list of all the individual's haplotype per pop per loci
		****genotM;				//list of all the individual's haplotype per pop per loci
	char c1,					//gets the value of a char temporarly
		 info[LOCMAXCHAR+1],	//information of a locus
		 info1[LOCMAXCHAR+1],	//information of a locus
		 info2[LOCMAXCHAR+1],	//information of a locus
		 aux[MAXCHAR],			//auxiliar
		 *outname,				//name of the output file 
		 **locname;				//name of the loci
	FILE *inp,					//pntr to the input file
		 *outp;					//pntr to the output file
	time_t startClock;			//time when the program starts
	const struct tm *startTime;	//struct time when the program starts

	inp = fopen(input,"r");
	if(inp == NULL)
		return 1; //cannot open sample file
	
	outsize = strlen(output) + 5;
	outname = (char *)malloc(outsize*sizeof(char));
	strcpy(outname,output);
	outp = fopen(strcat(outname,".len"),"w");
	if(outp == NULL)
		return 2; //cannot create .len file

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 

	//get nloc
	fgets(aux,MAXCHAR,inp);
	nloc=1;
	loctype=1;
	fgets(aux,MAXCHAR,inp);
	for(i=0;i<strlen(aux);i++){
		if(aux[i]==','){
			loctype=0;
			nloc+=1;
		}
	}
	if(loctype){
		searchPop=1;
		i=0;		
		while(searchPop){
			fgets(aux,MAXCHAR,inp);
			if(feof(inp))
				return 3; //file ended before the display of data
			j=0;
			while(isspace(aux[j])||aux[j]=='\t')
				j++;
			if(aux[j]=='P'||aux[j]=='p'){
				if(aux[j+1]=='O'||aux[j+1]=='o')
					if(aux[j+2]=='P'||aux[j+2]=='p')
						searchPop=0;
			}
			else
				i++;
		}
		nloc+=i;
	}
	else
		fgets(aux,MAXCHAR,inp);

	//get npop, lsamp
	reachEND=0;		
	cpop=0;
	while(!reachEND){
		lsamp = (int *)myAlloc(lsamp,(cpop+1)*sizeof(int));
		csamp=0;
		searchPop=1;
		while(searchPop && !reachEND){
			fgets(aux,MAXCHAR,inp);
			if(feof(inp))
				reachEND=1;
			j=0;
			while(isspace(aux[j])||aux[j]=='\t')
				j++;
			if(aux[j]=='P'||aux[j]=='p')
				if(aux[j+1]=='O'||aux[j+1]=='o')
					if(aux[j+2]=='P'||aux[j+2]=='p')
						searchPop=0;
			searchInd=1;
			for( ; j<strlen(aux) && searchInd; j++){
				if(aux[j]==','){
					csamp++;
					searchInd=0;
				}
			}
		}
		lsamp[cpop]=csamp;
		cpop++;
	}
	npop=cpop;

	rewind(inp);
	//print first line to output file
	fprintf(outp,"# GenePop file converted to .len file with genepop2table1.0\n");
	fprintf(outp,"# input file:  %s\n",input);
	fprintf(outp,"# output file: %s.len\n",output);
	fprintf(outp,"# date: %s#\n#",asctime(startTime));
	fgets(aux,MAXCHAR,inp);
	for(i=0;!isendline(aux[i]);i++)
		fprintf(outp,"%c",aux[i]);

	fprintf(outp,"\n");

	//get locus names
	locname=(char**)malloc(nloc*sizeof(char*));
	if(loctype){
		for(cloc=0;cloc<nloc;cloc++){
			locname[cloc] = malloc(MAXCHAR*sizeof(char));
			fgets(aux,MAXCHAR,inp);
			for(i=0;!isendline(aux[i]);i++)
				locname[cloc][i]=aux[i];
			locname[cloc][i]='\0';
		}	
	}
	else{
		fgets(aux,MAXCHAR,inp);
		for(j=0,cloc=0; cloc<nloc; cloc++){
			i=0;
			locname[cloc] = malloc(MAXCHAR*sizeof(char));
			while(aux[j]!=','&& !isendline(aux[j])){
				locname[cloc][i++]=aux[j++];
			}
			j++;
			locname[cloc][i]='\0';
		}
	}

	//allocate memory to lsamp, lplo and genotM
	lplo = (int *)malloc(nloc*sizeof(int));
	genotM = (int****)malloc(nloc*sizeof(int***));
	for(cloc=0;cloc<nloc;cloc++){
		genotM[cloc] = (int***)malloc(npop*sizeof(int**));
		for(cpop=0; cpop<npop; cpop++){
			genotM[cloc][cpop] = (int**)malloc(lsamp[cpop]*sizeof(int*));
				for(csamp=0; csamp<lsamp[cpop]; csamp++)
					genotM[cloc][cpop][csamp]=(int*)malloc(NPLOIDY*sizeof(int));
		}
	}
	
	//fill genotM[cloc][cpop][csamp][cplo]
	for(cpop=0 ; cpop<npop ; cpop++){
		fgets(aux,MAXCHAR,inp);
		for(csamp=0; csamp<lsamp[cpop] ; csamp++){
			fgets(aux,MAXCHAR,inp);
			j=0;
			while(aux[j]!=',')
				j++;
			for(cloc=0; cloc<nloc; cloc++){
				for( ; !isdigit(aux[j]); j++){
					if(isendline(aux[j])){
						fgets(aux,MAXCHAR,inp);
						j=0;
					}
				}
				for(i=0;isdigit(aux[j]);i++,j++){
					info[i]=aux[j];
				}
				info[i]='\0';
				if(strlen(info)==2){
					genotM[cloc][cpop][csamp][0]=atoi(info);
					lplo[cloc]=1;
				}
				else if(strlen(info)==3){
					genotM[cloc][cpop][csamp][0]=atoi(info);
					lplo[cloc]=1;
				}
				else if(strlen(info)==4){
					for(i=0;i<2;i++){
						info1[i]=info[i];
						info2[i]=info[i+2];
					}
					info1[i]='\0';
					info2[i]='\0';
					genotM[cloc][cpop][csamp][0]=atoi(info1);
					genotM[cloc][cpop][csamp][1]=atoi(info2);					
					lplo[cloc]=2;
				}
				else if(strlen(info)==6){
					for(i=0;i<3;i++){
						info1[i]=info[i];
						info2[i]=info[i+3];
					}
					info1[i]='\0';
					info2[i]='\0';
					genotM[cloc][cpop][csamp][0]=atoi(info1);
					genotM[cloc][cpop][csamp][1]=atoi(info2);					
					lplo[cloc]=2;
				}
				else
					return 4; //incorrect number of digits in locus

			}
		}
	}
	genotM2=(int***)malloc(nloc*sizeof(int**));
	for(cloc=0; cloc<nloc; cloc++){
		genotM2[cloc] = (int**)malloc(npop*sizeof(int*));
		for(cpop=0; cpop<npop; cpop++){
			genotM2[cloc][cpop] = (int*)malloc((lsamp[cpop]*NPLOIDY)*sizeof(int));
			for(i=0, csamp=0; csamp<lsamp[cpop]; csamp++){
				for(cplo=0; cplo<lplo[cloc]; cplo++){
					genotM2[cloc][cpop][i]=genotM[cloc][cpop][csamp][cplo];
					i++;
				}
			}
		}
	}
	
	nsamp=0;
	for(cpop=0; cpop<npop; cpop++){
		nsamp+=lsamp[cpop];
	}
	//allocate memory to freq, valM, IDsort, valS and ldna
	valM = (int **)malloc(nloc*sizeof(int *));
	IDsort = (int **)malloc(nloc*sizeof(int *));	
	ldna = (int *)malloc((nloc+1)*sizeof(int));
	for(cloc=0; cloc<nloc; cloc++){
		ldna[cloc]=0;
		valM[cloc] = malloc(nsamp*sizeof(int));
		IDsort[cloc] = malloc(nsamp*sizeof(int));
	}
	freq = (int ***)malloc(npop*sizeof(int **));
	for(cpop=0; cpop<npop; cpop++){
		freq[cpop] = (int **)malloc(nloc*sizeof(int *));	
		for(cloc=0; cloc<nloc ; cloc++){
			freq[cpop][cloc] = (int *)malloc(nsamp*sizeof(int));
		}
	}	

	//fill freq[cpop][cloc][cdna], valM[cloc][cdna] and IDsort[cloc][cdna]
	for(cloc=0;cloc<nloc;++cloc){
		for(cpop=0;cpop<npop;cpop++){
			for(i=0,csamp=0;csamp<lsamp[cpop];++csamp){
				for(cplo=0; cplo<lplo[cloc]; cplo++){
					//check if there is information in the current individual			
					if(genotM2[cloc][cpop][i] == 0){
						i++;
						continue;
					}
				
					//check if the current haplotype is already listed
					for(cdna=0;cdna<ldna[cloc];++cdna){
						equalDna = 1;						
					
						if(genotM2[cloc][cpop][i] != valM[cloc][cdna])
							equalDna = 0;
						if(equalDna){
							++freq[cpop][cloc][cdna];	//increase freq of current haplotype
							break;
						}
					}
					//current haplotype is not listed yet
					if(cdna>=ldna[cloc]){
						for(j=0;j<npop;++j){
							freq[j][cloc][ldna[cloc]] = 0;	//initiate as 0 the freq of current haplotype in other populations
						}
						freq[cpop][cloc][ldna[cloc]] = 1;	//initiate as 1 the freq of current haplotype in its population
						
						valM[cloc][cdna] = genotM2[cloc][cpop][i];		//adds current haplotype to valM[][]
						++ldna[cloc];	//increase by one the number of different haplotypes of a particular locus
					}
					i++;
				}
			}
		}
	}

	//writes the output file
	for(cpop=0; cpop<npop; cpop++)
		fprintf(outp,"#Population %d - pop%d\n",cpop+1,cpop+1);
	for(cloc=0;cloc<nloc;++cloc){
		fprintf(outp,"#Locus %d - %s\n",cloc+1,locname[cloc]);
	}
	fprintf(outp,"\n%d\n",npop);											//outp: npop
	fprintf(outp,"%d\n",nloc);	 											//outp: nloc		
	
	for(cloc=0 ; cloc<nloc ; cloc++){
		fprintf(outp,"m ");												//outp: ltype
	}
	fprintf(outp,"\n\n"); 														
	
	for(cloc=0; cloc<nloc ; cloc++){
		fprintf(outp,"%d\n",ldna[cloc]);									//outp: ndna

		isorti('a',ldna[cloc],valM[cloc],IDsort[cloc]);

		for(cpop=0;cpop<npop;++cpop){
			for(cdna=0;cdna<ldna[cloc];++cdna){
				fprintf(outp,"%4d ",freq[cpop][cloc][IDsort[cloc][cdna]]);	//outp: freq[][]
			}
			fprintf(outp,"\n");
		}
		fprintf(outp,"\n");
			
		for(cdna=0;cdna<ldna[cloc];++cdna){
			fprintf(outp,"%4d ",valM[cloc][IDsort[cloc][cdna]]);				//outp: valM[]
		}
		fprintf(outp,"\n\n");
	}
	
	//free stuff
	for(cloc=0 ; cloc<nloc ; cloc++){
		for(cpop=0 ; cpop<npop ; cpop++){
			for(csamp=0 ; csamp<lsamp[cpop] ; csamp++){
				free(genotM[cloc][cpop][csamp]);
			}
				free(genotM[cloc][cpop]);
		}	
		free(genotM[cloc]);
	}
	free(genotM);
	for(cloc=0; cloc<nloc ; cloc++){
		for(cpop=0 ; cpop<npop ; cpop++){
			free(genotM2[cloc][cpop]);
		}	
		free(genotM2[cloc]);
	}
	free(genotM2);
	for(cpop=0; cpop<npop ; cpop++){
		for(cloc=0; cloc<nloc ; cloc++){
			free(freq[cpop][cloc]);
		}
	    free(freq[cpop]);
	}
	free(freq);		
	for(cloc=0; cloc<nloc ; cloc++){
		free(valM[cloc]);
		free(IDsort[cloc]);
	}
	free(valM);
	free(IDsort);
	for(cloc=0 ; cloc<nloc ; cloc++){
		free(locname[cloc]);
	}
	free(locname);			
	free(lplo);
	free(ldna);		
	free(lsamp);
	free(outname);
		
	//close files
	fclose(inp);
	fclose(outp);

	return 0;

} //end of main
