/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		6th May 2009

*/

#include "interface.h"

/*
	It creates an ABC .len file given an IM input file.
	
	@arg length format filename
	@arg output filename
*/
int createFreqTab2(char *input,char *output){
	int cloc, cpop, cdna, cindiv, csite, i, j, cmark,	//iterators
		nmark,					//nr. of total markers used
		nloc,					//number of loci
		equalDna,				//auxiliar to build the .len file
		outsize,				//size of output file name 
  		*ldna,					//number of different dna data by loci
		*lmark,					//nr of markers per loci
		*nsamp,					//total number of samples per loci
		**lsamp,				//size of sample per pop per loci	
		***freq;				//freq of dna data by no_pop by no. of diff dna data
	char type,					//type of analysis choosen
		 c1,					//gets the value of a char temporarly
		 aux[MAXCHAR],			//auxiliar
		 pop1[MAXCHAR],			//name of population 1
		 pop2[MAXCHAR],			//name of population 2
		 DNAtype[MAXCHAR],		//current locus mutation model used in IM
		 *ltype,				//DNA type per loci
		 *outname,				//name of the output file 
		 **locname;				//name of the loci
	FILE *inp,					//pntr to the input file
		 *outp;					//pntr to the output file
	time_t startClock;			//time when the program starts
	const struct tm *startTime;	//struct time when the program starts
// only used in microssatelite dna analysis	
	int **valM,					//list of all the diff microsatellite sizes
		**IDsort,				//list of the order of the diff microsatellite sizes
		***genotM2,				//list of all the individual's haplotype per pop per loci
		****genotM;				//list of all the individual's haplotype per pop per loci
// only used in sequence dna analysis
	char **mainseq,				//main sequence of a loci
		 ***valS,				//list of all the diff dna sequencies
		 ****genotS2,			//list of all the individual's haplotype per pop per loci 
		 ****genotS;			//list of all the individual's haplotype per pop per loci 
	int *lsites;				//number of bases per loci

	inp = fopen(input,"r");
	if(inp == NULL)
		return 1; //cannot open sample file
	
	outsize = strlen(output) + 5;
	outname = (char *)malloc(outsize*sizeof(char));
	strcpy(outname,output);
	outp = fopen(strcat(outname,".len"),"w");
	if(outp == NULL)
		return 2; //cannot create output file

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 

	//print first line to output file
	fprintf(outp,"# IMa file converted to .len file with im2table1.0\n");
	fprintf(outp,"# input file:  %s\n",input);
	fprintf(outp,"# output file: %s.len\n",output);
	fprintf(outp,"# date: %s#\n#",asctime(startTime));
	while(!isendline(c1=getc(inp))){
		fprintf(outp,"%c",c1);
	}
	fprintf(outp,"\n");

	//if # (commentary) run to end of line
	c1 = getc(inp);
	while(c1 == '#'){
		while(!isendline(c1 = getc(inp)));
		while(isspace(c1=getc(inp))||isendline(c1)||c1=='\t');
	}
	ungetc(c1,inp);

	fscanf(inp,"%s",&pop1);
	fscanf(inp,"%s",&pop2);
	fscanf(inp,"%d",&nloc);

	//allocate memory to DNAtype, lsamp, locname, genotM or genotS and lsites
	ltype = (char *)malloc(nloc*sizeof(char));
	lmark = (int *)malloc(nloc*sizeof(int));
	lsamp = (int **)malloc(nloc*sizeof(int*));
	nsamp = (int *)malloc(nloc*sizeof(int));
	locname = (char **)malloc(nloc*sizeof(char*));
	lsites = (int *)malloc(nloc*sizeof(int));
	genotM = (int****)malloc(nloc*sizeof(int***));
	genotS = (char****)malloc(nloc*sizeof(char***));
	for(cloc=0;cloc<nloc;cloc++){
		lsamp[cloc] = (int*)malloc(NPOP*sizeof(int));
		locname[cloc] = (char *)malloc(MAXCHAR*sizeof(char));
		genotM[cloc] = (int***)malloc(NPOP*sizeof(int**));
		genotS[cloc] = (char***)malloc(NPOP*sizeof(char**));
	}

	//fill genotM[cloc][cpop][cindiv] or genotS[cloc][cpop][cindiv][csite]
	for(cloc=0;cloc<nloc;cloc++){
		while(isendline(c1=getc(inp)));
		ungetc(c1,inp);
		for(i=0 ; !isspace(c1=getc(inp)) ; i++)
			locname[cloc][i]=c1;
		locname[cloc][i]='\0';
		fscanf(inp,"%d",&lsamp[cloc][0]);
		fscanf(inp,"%d",&lsamp[cloc][1]);
		fscanf(inp,"%d",&lsites[cloc]);
		fscanf(inp,"%s",&DNAtype);
		if(DNAtype[0]=='I'){
			ltype[cloc]='s';
			lmark[cloc]=1;
		}
		else if(DNAtype[0]=='S'){
			ltype[cloc]='m';
			for(i=1; isdigit(DNAtype[i]); i++)
				aux[i-1]=DNAtype[i];
			lmark[cloc]=atoi(aux);
		}
		else
			return 3; //program can't deal with specified data
		while(!isendline(c1=getc(inp)));
		ungetc(c1,inp);
		while(isendline(c1=getc(inp)));
		ungetc(c1,inp);
		
		for(cpop=0 ; cpop<NPOP ; cpop++){
			//choose between diferent analysis		
			if(ltype[cloc] == 'm'){
				genotM[cloc][cpop] = (int **)malloc((lsamp[cloc][cpop])*sizeof(int*));
				for(cindiv=0 ; cindiv<lsamp[cloc][cpop] ; cindiv++){
					genotM[cloc][cpop][cindiv] = (int *)malloc(lmark[cloc]*sizeof(int));
					for(i=0 ; i<=LOCMAXCHAR ; i++)
						c1=getc(inp);
					if(isspace(c1)){
						while(isspace(c1=getc(inp)));
					}
					ungetc(c1,inp);
					for(cmark=0; cmark<lmark[cloc]; cmark++)
						fscanf(inp,"%d",&genotM[cloc][cpop][cindiv][cmark]);
					while(!isendline(c1=getc(inp)));
					ungetc(c1,inp);
					while(isendline(c1=getc(inp)));
					ungetc(c1,inp);
				}
			}
			else{
				genotS[cloc][cpop] = (char **)malloc((lsamp[cloc][cpop])*sizeof(char *));
				for(cindiv=0 ; cindiv<lsamp[cloc][cpop] ; cindiv++){
					genotS[cloc][cpop][cindiv] = (char *)malloc(lsites[cloc]*sizeof(char));
					for(i=0 ; i<=LOCMAXCHAR ; i++)
						c1=getc(inp);
					if(isspace(c1)){
						while(isspace(c1=getc(inp)));
					}
					ungetc(c1,inp);
					for(csite=0 ; csite < lsites[cloc] ; csite++)
						fscanf(inp,"%c",&genotS[cloc][cpop][cindiv][csite]);
					while(!isendline(c1=getc(inp)));
					ungetc(c1,inp);
					while(isendline(c1=getc(inp)));
					ungetc(c1,inp);
				}	
			}
		}
	}
	for(cloc=0; cloc<nloc; cloc++){
		nsamp[cloc]=0;
		for(cpop=0; cpop<NPOP; cpop++){
			nsamp[cloc]+=lsamp[cloc][cpop];
		}
	}
	
	nmark=0;
	for(cloc=0; cloc<nloc; cloc++)
		nmark+=lmark[cloc];

	genotM2=(int***)malloc(nmark*sizeof(int**));
	for(i=0,cloc=0; cloc<nloc; cloc++){
		for(cmark=0; cmark<lmark[cloc]; cmark++){
			if(ltype[cloc]=='m'){
				genotM2[i] = (int**)malloc(NPOP*sizeof(int*));
				for(cpop=0; cpop<NPOP; cpop++){
					genotM2[i][cpop] = (int*)malloc((lsamp[cloc][cpop])*sizeof(int));
					for(cindiv=0; cindiv<lsamp[cloc][cpop]; cindiv++)
						genotM2[i][cpop][cindiv]=genotM[cloc][cpop][cindiv][cmark];
				}
			}
			i++;
		}
	}

	genotS2=(char****)malloc(nmark*sizeof(char***));
	for(i=0,cloc=0; cloc<nloc; cloc++){
		for(cmark=0; cmark<lmark[cloc]; cmark++){
			if(ltype[cloc]=='s'){
				genotS2[i] = (char***)malloc(NPOP*sizeof(char**));
				for(cpop=0; cpop<NPOP; cpop++){
					genotS2[i][cpop] = (char**)malloc((lsamp[cloc][cpop])*sizeof(char*));
					for(cindiv=0; cindiv<lsamp[cloc][cpop]; cindiv++){
						genotS2[i][cpop][cindiv] = (char *)malloc(lsites[cloc]*sizeof(char));
						for(csite=0 ; csite < lsites[cloc] ; csite++)
							genotS2[i][cpop][cindiv][csite]=genotS[cloc][cpop][cindiv][csite];
					}
				}
			}
			i++;
		}
	}
		
	//allocate memory to freq, valM, IDsort, valS and ldna
	valM = (int **)malloc(nmark*sizeof(int *));
	IDsort = (int **)malloc(nmark*sizeof(int *));	
	valS = (char ***)malloc(nmark*sizeof(char **));
	ldna = (int *)malloc((nmark+1)*sizeof(int));
	for(i=0,cloc=0; cloc<nloc ;++cloc){
		for(cmark=0; cmark<lmark[cloc];cmark++){
			ldna[i]=0;
			if(ltype[cloc]=='m'){
				valM[i] = malloc(nsamp[cloc]*sizeof(int));
				IDsort[i] = malloc(nsamp[cloc]*sizeof(int));
			}
			else{
				valS[i] = (char **)malloc((nsamp[cloc]+1)*sizeof(char *));
			}
			i++;	
		}
	}
	freq = (int ***)malloc(NPOP*sizeof(int **));
	for(cpop=0;cpop<NPOP;++cpop){
		freq[cpop] = (int **)malloc(nmark*sizeof(int *));	
		for(i=0,cloc=0; cloc<nloc ;++cloc){
			for(cmark=0; cmark<lmark[cloc];cmark++){
				freq[cpop][i] = (int *)malloc(nsamp[cloc]*sizeof(int));
				i++;
			}
		}
	}	

	//fill freq[cpop][cloc][cdna], valM[cloc][cdna] and IDsort[cloc][cdna]
	for(i=0,cloc=0;cloc<nloc;++cloc){
		for(cmark=0; cmark<lmark[cloc];cmark++){
			if(ltype[cloc]=='m'){
				for(cpop=0;cpop<NPOP;cpop++){
					//run through every individual		
					for(cindiv=0;cindiv<lsamp[cloc][cpop];++cindiv){
						//check if there is information in the current individual			
						if(genotM2[i][cpop][cindiv] == 0){
							continue;
						}
				
						//check if the current haplotype is already listed
						for(cdna=0;cdna<ldna[i];++cdna){
							equalDna = 1;						
						
							if(genotM2[i][cpop][cindiv] != valM[i][cdna])
								equalDna = 0;
							if(equalDna){
								++freq[cpop][i][cdna];	//increase freq of current haplotype
								break;
							}
						}
						//current haplotype is not listed yet
						if(cdna>=ldna[i]){
							for(j=0;j<NPOP;++j){
								freq[j][i][ldna[i]] = 0;	//initiate as 0 the freq of current haplotype in other populations
							}
							freq[cpop][i][ldna[i]] = 1;	//initiate as 1 the freq of current haplotype in its population
							
							valM[i][cdna] = genotM2[i][cpop][cindiv];		//adds current haplotype to valM[][]
							++ldna[i];	//increase by one the number of different haplotypes of a particular locus
						}
					}
				}
			}
			i++;
		}
	}

	//fill freq[cpop][cloc][cdna] and valS[cloc][cdna][csite]
	for(i=0,cloc=0;cloc<nloc;++cloc){
		for(cmark=0; cmark<lmark[cloc];cmark++){
			if(ltype[cloc]=='s'){
				for(cpop=0;cpop<NPOP;cpop++){
					//run through every individual		
					for(cindiv=0;cindiv<lsamp[cloc][cpop];++cindiv){
						//check if the current haplotype is already listed
						for(cdna=0;cdna<ldna[i];++cdna){
							equalDna = 1;						
							for(csite=0 ; csite<lsites[cloc] ; csite++){
								if(genotS2[i][cpop][cindiv][csite] != valS[i][cdna][csite]){
									if(valS[i][cdna][csite]=='N'){
										valS[i][cdna][csite] = genotS2[i][cpop][cindiv][csite];
									}
									else if(genotS2[i][cpop][cindiv][csite]=='A'||genotS2[i][cpop][cindiv][csite]=='a'||
											genotS2[i][cpop][cindiv][csite]=='T'||genotS2[i][cpop][cindiv][csite]=='t'||
									  		genotS2[i][cpop][cindiv][csite]=='G'||genotS2[i][cpop][cindiv][csite]=='g'||
									  		genotS2[i][cpop][cindiv][csite]=='C'||genotS2[i][cpop][cindiv][csite]=='c'){
										equalDna = 0;
										break;
									}
								}
							}
							if(equalDna){
								++freq[cpop][i][cdna];	//increase freq of current haplotype
								break;
							}
						}
						//current haplotype is not listed yet
						if(cdna>=ldna[i]){
							for(j=0;j<NPOP;++j){
								freq[j][i][ldna[i]] = 0;//initiate as 0 the freq of current haplotype in other populations
							}
							freq[cpop][i][ldna[i]] = 1;	//initiate as 1 the freq of current haplotype in its population
							
							valS[i][cdna] = (char *)malloc((lsites[cloc]+1)*sizeof(char));
							for(csite=0; csite<lsites[cloc] ; csite++){
								valS[i][cdna][csite] = genotS2[i][cpop][cindiv][csite];		//adds current haplotype to valS[][][]	
							}
							valS[i][cdna][csite] = '\0';
							++ldna[i];	//increase by one the number of different haplotypes of a particular locus
						}
					}
				}
			}
			i++;
		}
	}
	//transforms dna sequences into an array of segregating sites	
	mainseq = (char **)malloc(nmark*sizeof(char *));
    for(i=0,cloc=0 ; cloc<nloc ; cloc++){
    	for(cmark=0; cmark<lmark[cloc]; cmark++){
			if(ltype[cloc]=='s'){					
		        mainseq[i] = (char *)malloc(lsites[cloc]*sizeof(char));
	            for(csite=0; csite<lsites[cloc]; csite++){
					mainseq[i][csite] = valS[i][0][csite];
		        }
				for(cdna=0;cdna<ldna[i];++cdna){
					for(csite=0;csite<lsites[cloc];++csite){
						if(valS[i][cdna][csite]!=mainseq[i][csite]&&(valS[i][cdna][csite]=='A'||
																	 valS[i][cdna][csite]=='a'||
																	 valS[i][cdna][csite]=='G'||
																	 valS[i][cdna][csite]=='g'||
																	 valS[i][cdna][csite]=='T'||
																	 valS[i][cdna][csite]=='t'||
																	 valS[i][cdna][csite]=='C'||
																	 valS[i][cdna][csite]=='c')){
							valS[i][cdna][csite] = '1';
						}
						else{
							valS[i][cdna][csite] = '0';
						}
					}
				} 
			}
			i++;
    	}
    }

	//writes the output file
	fprintf(outp,"#Population 1 - %s\n",pop1);
	fprintf(outp,"#Population 2 - %s\n",pop2);	

	for(cloc=0,i=0;cloc<nloc;++cloc){
	   	for(cmark=0; cmark<lmark[cloc]; cmark++){
			if(lmark[cloc]>1)
				fprintf(outp,"#Locus %d - %s_%d\n",i+1,locname[cloc],cmark+1);
			else
				fprintf(outp,"#Locus %d - %s\n",i+1,locname[cloc]);
			i++;
	   	}
	}
	fprintf(outp,"\n2\n"); 														//outp: npop
	fprintf(outp,"%d\n",nmark);	 												//outp: nloc		
	
	for(cloc=0 ; cloc<nloc ; cloc++){
	   	for(cmark=0; cmark<lmark[cloc]; cmark++)
			fprintf(outp,"%c ",ltype[cloc]);									//outp: ltype
	}
	fprintf(outp,"\n\n"); 														
	
	for(cloc=0,i=0 ; cloc<nloc ; cloc++){
	   	for(cmark=0; cmark<lmark[cloc]; cmark++){
			fprintf(outp,"%d\n",ldna[i]);										//outp: ndna
		
			if(ltype[cloc]=='m')
				isorti('a',ldna[i],valM[i],IDsort[i]);
			
			for(cpop=0;cpop<NPOP;++cpop){
				for(cdna=0;cdna<ldna[i];++cdna){
					if(ltype[cloc]=='m')
						fprintf(outp,"%4d ",freq[cpop][i][IDsort[i][cdna]]);	//(M)outp: freq[][]
					else
						fprintf(outp,"%4d ",freq[cpop][i][cdna]);				//(S)outp: freq[][]
				}
				fprintf(outp,"\n");
			}
			fprintf(outp,"\n");
			
			if(ltype[cloc]=='s'){
				for(cdna=0;cdna<ldna[i];++cdna)
					fprintf(outp,"%4d ",cdna);										//(S)outp: cdna
				fprintf(outp,"\n\n");
				fprintf(outp,"%d \n",lsites[cloc]);									//(S)outp: lsites[]
			}	
			
			for(cdna=0;cdna<ldna[i];++cdna){
				if(ltype[cloc]=='m')
					fprintf(outp,"%4d ",valM[i][IDsort[i][cdna]]);			//(M)outp: valM[]
				else
					fprintf(outp,"%d %s\n",cdna,valS[i][cdna]);					//(S)outp: cdna, valS[]
			}
			if(ltype[cloc]=='m')
				fprintf(outp,"\n\n");
			else
				fprintf(outp,"\n");
			i++;
		}
	}
	
	//free stuff
	for(cloc=0 ; cloc<nloc ; cloc++){
		if(ltype[cloc]=='s'){
			for(cpop=0 ; cpop<NPOP ; cpop++){
				for(cindiv=0 ; cindiv<lsamp[cloc][cpop] ; cindiv++){
						free(genotS[cloc][cpop][cindiv]);	
				}
					free(genotS[cloc][cpop]);
			}	
		}
		free(genotS[cloc]);
	}
	free(genotS);
	for(cloc=0 ; cloc<nloc ; cloc++){
		if(ltype[cloc]=='m'){
			for(cpop=0 ; cpop<NPOP ; cpop++){
				for(cindiv=0 ; cindiv<lsamp[cloc][cpop] ; cindiv++){
					free(genotM[cloc][cpop][cindiv]);
				}
					free(genotM[cloc][cpop]);
			}	
		}
		free(genotM[cloc]);
	}
	free(genotM);
	for(cloc=0, i=0; cloc<nloc ; cloc++){
		for(cmark=0;cmark<lmark[cloc];cmark++){
			if(ltype[cloc]=='s'){
				for(cpop=0 ; cpop<NPOP ; cpop++){
					for(cindiv=0 ; cindiv<lsamp[cloc][cpop] ; cindiv++){
							free(genotS2[i][cpop][cindiv]);	
					}
						free(genotS2[i][cpop]);
				}	
				free(genotS2[i]);
			}
			i++;
		}
	}
	free(genotS2);
	for(cloc=0, i=0; cloc<nloc ; cloc++){
		for(cmark=0;cmark<lmark[cloc];cmark++){
			if(ltype[cloc]=='m'){
				for(cpop=0 ; cpop<NPOP ; cpop++){
					free(genotM2[i][cpop]);
				}	
				free(genotM2[i]);
			}
			i++;
		}
	}
	free(genotM2);
	for(cpop=0; cpop<NPOP ; cpop++){
		for(cloc=0,i=0; cloc<nloc ; cloc++){
			for(cmark=0;cmark<lmark[cloc];cmark++){
				free(freq[cpop][i]);
				i++;
			}
		}
	    free(freq[cpop]);
	}
	free(freq);		
	for(cloc=0, i=0; cloc<nloc ; cloc++){
		for(cmark=0; cmark<lmark[cloc];cmark++){
			if(ltype[cloc]=='s'){
				for(cdna=0 ; cdna<ldna[i] ; cdna++)
					free(valS[i][cdna]);
				free(valS[i]);
				free(mainseq[i]);
			}
			else{
				free(valM[i]);
				free(IDsort[i]);
			}
			i++;
		}
	}
	free(mainseq);
	free(valS);
	free(valM);
	free(IDsort);
	for(cloc=0 ; cloc<nloc ; cloc++){
		free(lsamp[cloc]);
		free(locname[cloc]);
	}
	free(lsamp);
	free(locname);			
	free(lmark);
	free(ldna);		
	free(nsamp);
	free(lsites);
	free(ltype);
	free(outname);
		
	//close files
	fclose(inp);
	fclose(outp);

	return 0;

} //end of main
