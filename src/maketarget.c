/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		1th May 2009
	
	NBB - based on Mark's make_target.c and make_targetS.c
*/

#include "maketarget.h"

/*
	This funtion uses a freq_tab_length file and creates a output file which contains
	the summary statistics of the given file and its sample sizes
	
	@arg input filename (*.len)
	@arg input filename (*.sst)
	@arg output filename (*.trg AND *.szz)
*/
int main(int argc, char* argv[]){
	int cloc, cpop, cpop2, cdna, i,	//iterators
		dump,					//variable just to put the information from data that isn't necessary, dump
		started,				//auxiliar to help print sstats to the .txt file
		foundSTR,				//check if STR's are present
		foundSNP,				//check if SNP's are present
		nstats_aux,				//auxiliar to calculate the number of summ stats
		nstats,					//number of summ stats
		outsize,				//size of output file name
		sum,					//sample size in one population
		nloc,					//number of loci
		npop,					//number of populations
		ndna,					//number of all the diferent haplotypes
		lsstats[MAXSSTATS];	//list of used summary statistics
	char c1,					//auxiliar
		 *ltype,				//DNA type per loci
		 *outline,				//store the target data summary statistic
		 *name_inf,             //filename of the report file
 		 *name_trg,				//name of the output.tfq file
		 *lsstats2[MAXSSTATS] = {"H","varL","k_M","curL","sH_M","NmH","pi","S","k_S","sH_S","avMFS","sdMFS","NmS","privS","S(1)"},
		 *name_ssz;				//name of the output.out file
	FILE *inp_len,				//pntr to data file
		 *inp_sst,				//pntr to data file
         *out_inf,              //pntr to the report file
		 *out_trg,				//pntr to target output file
		 *out_ssz;				//pntr to samp_size output file
	time_t startClock,			//time_t when the program starts
		   endClock;	 		//time_t when the program ends
	struct data block,			//target data
				*data;			//pntr to target data	
	const struct tm *startTime,	//struct time when the program starts
					*endTime;	//struct time when the program ends

	if(argc != 4)
		printerr("needs .len file, .sst file, output filename (no extension)");
		
	inp_len = fopen(argv[1],"r");						//input filename
	if(inp_len == NULL)
		printerr("cannot open .len file");
	inp_sst = fopen(argv[2],"r");			//input .sst
	if(inp_sst == NULL)
		printerr("cannot open .sst file");
	outsize = strlen(argv[3]) + 5;
	name_trg = (char *)malloc(outsize*sizeof(char));
	name_ssz = (char *)malloc(outsize*sizeof(char));
    name_inf = (char *)malloc((outsize+4)*sizeof(char));
	strcpy(name_trg,argv[3]);
	strcpy(name_ssz,argv[3]);
    strcpy(name_inf,argv[3]);
	out_trg = fopen(strcat(name_trg,".trg"),"w");	//out_trg filename
	if(out_trg == NULL)
		printerr("cannot create .trg file");
	out_ssz = fopen(strcat(name_ssz,".ssz"),"w");	//out_sps filename
	if(out_ssz == NULL)
		printerr("cannot create .ssz file");
    out_inf = fopen(strcat(name_inf,"_trg.txt"),"w");       //open out_inf
    if(out_inf == NULL)
        printerr("cannot create .txt file");

	for(i=0; i<MAXSSTATS; i++){
		while(getc(inp_sst)!='#');	
		if(EOF == fscanf(inp_sst,"%d",&lsstats[i]))
			printerr("error reading .sst file");
	}
	fclose(inp_sst);
	/*check for errors when choosing summary statistics*/
	if(npop==1 && (lsstats[5] || lsstats[12] || lsstats[13] || lsstats[14]))
		printerr("Following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)");
	//if we want to use the Nm_H we have to calculate H
	if(lsstats[5]== 1)
		lsstats[0] = 1;
	//if we want to use the sdMFS we have to calculate mMFS
	if(lsstats[11]== 1)
		lsstats[10] = 1;
	//if we want to use the mMFS, Nm_S, privateS or S(1) we have to calculate S
	if(lsstats[10]== 1 || lsstats[12]== 1 || lsstats[13]== 1 || lsstats[14]== 1 )
		lsstats[7] = 1;

	time( &startClock );   					// Get time in seconds
	startTime = localtime( &startClock );  	// Convert time to struct tm form 
	
	/*if # (commentary) run to end of line*/
	c1=getc(inp_len);
	while(c1 == '#'){
		while(!isendline(c1 = getc(inp_len)));
		while(isspace(c1=getc(inp_len))||isendline(c1)||c1=='\t');
	}
	ungetc(c1,inp_len);
	
	fscanf(inp_len,"%d %d",&npop,&nloc);

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

	data = &block;			
	data->npop = npop;
	data->nloc = nloc;
	data->ldna = (int *)malloc(nloc*sizeof(int));
	data->tsamp = (int *)malloc(nloc*sizeof(int));		
	data->nsamp = (int **)malloc(nloc*sizeof(int *));
	data->freq = (int ***)malloc(nloc*sizeof(int **));

	if(foundSTR)
		data->valM = (int **)malloc(nloc*sizeof(int *));
	if(foundSNP){
		data->lsites = (int *)malloc(nloc*sizeof(int));
		data->valS = (char ***)malloc(nloc*sizeof(char **));
	}
		
	/*run through every locus*/
	for(cloc=0;cloc<nloc;++cloc){
		data->freq[cloc] = (int **)malloc(npop*sizeof(int *));
		data->nsamp[cloc] = (int *)malloc(npop*sizeof(int));
		data->tsamp[cloc]=0;
		if(foundSNP)
			data->lsites[cloc]=0;
		if(ltype[cloc]=='m'||ltype[cloc]=='M'){
			fscanf(inp_len,"%d",&ndna);
			data->ldna[cloc] = ndna;
			data->valM[cloc] = (int *)malloc(ndna*sizeof(int));
			
			//run through every population
			for(cpop=0;cpop<npop;++cpop){
				sum = 0;
				data->freq[cloc][cpop] = (int *)malloc(ndna*sizeof(int));
				for(cdna=0;cdna<ndna;++cdna){
					fscanf(inp_len,"%d",&(data->freq[cloc][cpop][cdna]));
					sum += data->freq[cloc][cpop][cdna];
				}
				data->nsamp[cloc][cpop] = sum;
				data->tsamp[cloc]+=sum;	
				fprintf(out_ssz,"%d ",sum);		//out_sps: sum						
			}
			fprintf(out_ssz,"\n");

			for(cdna=0;cdna<ndna;++cdna)
				fscanf(inp_len,"%d",&(data->valM[cloc][cdna]));
		}
		else{
			fscanf(inp_len,"%d",&ndna);
			data->ldna[cloc] = ndna;
			data->valS[cloc] = (char **)malloc(ndna*sizeof(char *));
		
			//run through every population
			for(cpop=0;cpop<npop;++cpop){
				sum = 0;
				data->freq[cloc][cpop] = (int *)malloc(ndna*sizeof(int));
				for(cdna=0;cdna<ndna;++cdna){
					fscanf(inp_len,"%d",&(data->freq[cloc][cpop][cdna]));
					sum += data->freq[cloc][cpop][cdna];
				}
				data->nsamp[cloc][cpop] = sum;
				data->tsamp[cloc]+=sum;	
				fprintf(out_ssz,"%d ",sum);		//NT- out_sps: sum						
			}
			fprintf(out_ssz,"\n");
						
			for(cdna=0;cdna<ndna;++cdna)
				fscanf(inp_len,"%d",&dump); 
	
			fscanf(inp_len,"%d",&(data->lsites[cloc]));
			if(data->lsites[cloc]!=0)
				for(cdna=0;cdna<ndna;++cdna){
					data->valS[cloc][cdna] = (char *)malloc((data->lsites[cloc]+1)*sizeof(char));
					fscanf(inp_len,"%d %s ",&dump,data->valS[cloc][cdna]);
				}
			else
				printerr("use of Sequence data locus without information, this should be discarded");
		}
	}

	/*counting the number of summary statistics*/
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
			for(i=0 ; i<MAXSSTATS_M ; i++){
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
		if(foundSTR){
			for(i=0 ; i<MAXSSTATS_M ; i++){
				if(lsstats[i])
					nstats++;
			}
		}
		if(foundSNP){
			for(i=MAXSSTATS_M ; i<MAXSSTATS ; i++){
				if(lsstats[i])
					nstats++;
			}
		}
	}

	/*writting first part of output to .txt file*/
	fprintf(out_inf,"SummData - Mark Beaumont & Joao Lopes\n\n");
	fprintf(out_inf," last updated 01/05/09\n\n");
	fprintf(out_inf,"INPUT AND STARTING INFORMATION\n");
	fprintf(out_inf,"-------------------------------\n");
	fprintf(out_inf,"Table file :            %s\n",argv[1]);
	fprintf(out_inf,"SummaryStatistics file: %s\n\n",argv[2]);
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
		fprintf(out_inf,"%c ",ltype[cloc]);
	fprintf(out_inf,"\n\nOUTPUT FILES:\n");
	fprintf(out_inf,"-------------------------------\n");
	fprintf(out_inf,"Sample size file:       %s.ssz\n",argv[3]);
	fprintf(out_inf,"Target file:            %s.trg\n",argv[3]);
	fprintf(out_inf,"%d summstats [",nstats);
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

	
	outline = malloc(nstats*15*sizeof(char));
	
	summStats(data,lsstats,outline,foundSTR,foundSNP,ltype);
	
	fprintf(out_trg,"%s",outline);			//out_trg: outline

	/*get ending time and save it to .txt file*/
	time( &endClock );   				// Get time in seconds
	endTime = localtime( &endClock );  	// Convert time to struct tm form 
	fprintf(out_inf,"Ending date:   %s\n",asctime(endTime));
	fprintf(out_inf,"\nEND OF OUTPUT\n");

	/*free stuff*/
	free(outline);
	free(name_trg);
	free(name_ssz);
    free(name_inf); 
	
	for(cloc=0 ; cloc<nloc ; cloc++){
		if(ltype[cloc]=='M'||ltype[cloc]=='m')
			free(data->valM[cloc]);
		else{
			for(cdna=0; cdna<data->ldna[cloc]; cdna++)
					free(data->valS[cloc][cdna]);
				free(data->valS[cloc]);
		}
	}
	
	if(foundSTR)
		free(data->valM);
	if(foundSNP){
		free(data->valS);
		free(data->lsites);
	}
	
	for(cloc=0 ; cloc<nloc ; cloc++){
		for(cpop=0 ; cpop<npop ; cpop++)
			free(data->freq[cloc][cpop]);
		free(data->freq[cloc]);
		free(data->nsamp[cloc]);
	}
	free(data->freq);
	free(data->nsamp);
	free(data->tsamp);
	free(data->ldna);
	free(ltype);
	
	/*close files*/
	fclose(inp_len);
    fclose(out_inf);
	fclose(out_trg);
	fclose(out_ssz);
	
}//end of main
