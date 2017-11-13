/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		1th May 2009
	 
	NBB - based on Mark's firstpass.c
*/

#include "firstpass.h"

/*
	This program prints out a proportion, specified by tolerance, of the simulated data with an ABC approach
	which are closer in Euclidian distance to the target summary statistics. It will also create a file
	which will contain the simulated populational tree parameters (only) in the first 10000 lines of the data.
	These can then be used to build a posterior distributions with which we can compare prior distributions

	@arg filename with simulated data
	@arg filename with target summary statistics (summary statistics from our "real" data)
	@arg filename of the output of the program
	@arg number of parameters
	@arg number of summstats
	@arg tolerance of analysis (between 0 and 1)
*/
int main(int argc, char *argv[]){
	int cstat,cparam,i,			//iterators
		outsize,				//size of the output filename
		nclosest,				//tolerance*number of simulations
		npriors,				//no. sets to choose from the random prior set values (0-10000)
		nstats,					//no. of total summary statistics calculated
		nparam,					//no. of parameters of the populacional model
		*index,					//stores the original position of sorted elements of the summary statistics
		gap,					//auxiliar to build .pri
		isim;					//index of a simulation
	long cline,csim,			//iterators
		 nsim;					//no. of simulations of genetic trees from the same populational tree
	double  d1,					//auxiliar to read text
			tol,				//tolerance value
			*trgStat,			//list of the target summary statistics
			**iparams,			//stores the value of each parameter for the first points (0-10000)
			**stats,			//stores the value of each summary statistics for each tree
			**params,			//stores the value of each parameter for each tree
			**fstats,			//stores the value of each summary statistics for the closest points
			**fparams,			//stores the value of each parameter for the closest points
			*ave,				//stores the average of each list of different summary statistic
			*sdev,				//stores the standard deviation of each list of different summary statistic
			d3,d4,d5,d6,		//auxiliar values to get the average and the variance of each list of each summary statistic
			*dist;				//stores the distance between each simulated data and the target data
	char *name_dat,				//filename with the simulated data
		 *name_trg,				//filename with the target summary statistics
		 *name_inf,             //filename of the report file
         *name_rej,				//filename of the first output (nearest simulations in a given tolerance)
		 *name_pri;				//filename of the second output (first 10000 simulations)
	FILE *out_dat,				//pntr to the simulated data
		 *out_trg,				//pntr to the target summary statitstics
         *out_inf,              //pntr to the report file
		 *out_rej,				//pntr to the file with first output (nearest simulations in a given tolerance)
		 *out_pri;				//pntr to the file with second output (first 10000 simulations)
    time_t startClock,          //time_t when the program starts
           endClock;            //time_t when the program ends
    const struct tm *startTime, //struct time when the program starts
                    *endTime;   //struct time when the program ends

 	if(!(argc==7))
		printerr("needs .dat file, .trg file, output filename (no extension), nparams, nsstats, tolerance.");

	name_dat = argv[1];
	name_trg = argv[2];
	sscanf(argv[4],"%d",&nparam);
	sscanf(argv[5],"%d",&nstats);
	sscanf(argv[6],"%lf",&tol);

	out_dat = fopen(name_dat,"r");						//open out_dat
	if(out_dat == NULL)
		printerr("cannot open .dat file");

	out_trg = fopen(name_trg,"r");						//open out_trg
	if(out_trg == NULL)
		printerr("cannot open .trg file");

	name_rej = malloc((strlen(argv[3])+5)*sizeof(char));
	strcpy(name_rej,argv[3]);
	out_rej = fopen(strcat(name_rej,".rej"),"w");		//open out_rej
	if(out_rej == NULL)
		printerr("cannot create .rej file");
	free(name_rej);

	name_pri = malloc((strlen(argv[3])+5)*sizeof(char));
	strcpy(name_pri,argv[3]);
	out_pri = fopen(strcat(name_pri,".pri"),"w");		//open out_pri
	if(out_trg == NULL)
		printerr("cannot create .pri file");
	free(name_pri);

    name_inf = malloc((strlen(argv[3])+9)*sizeof(char));
    strcpy(name_inf,argv[3]);
    out_inf = fopen(strcat(name_inf,"_rej.txt"),"w");       //open out_inf
    if(out_inf == NULL)
        printerr("cannot create .txt file");
    free(name_inf); 
    
    /*fill trgStat[]*/
	trgStat = (double *)malloc(nstats*sizeof(double));
	for(cstat=0;cstat<nstats;++cstat)
		fscanf(out_trg,"%lf",&trgStat[cstat]);

    /*get the starting time*/
    time( &startClock );                    // Get time in seconds
    startTime = localtime( &startClock );   // Convert time to struct tm form 

	printf("\n\n--\nCheck integrity of .dat:\n");
	printf("\nStart checking...\n");

	/*Verify if the data file is well build and get the number of geneological trees simulated*/
	cline=0;
	while(1){
		for(cparam=0 ; cparam<nparam ; ++cparam){
			i = fscanf(out_dat,"%lf",&d3);
			if(i == EOF){
				if(cparam==0)
					break;
				else
					printerr(".dat file ended prematurely");
			}
		}
		if(cparam==0)
			break;
		for(cstat=0;cstat<nstats;++cstat){
			i = fscanf(out_dat,"%lf",&d3);
			if(i == EOF)
				printerr(".dat file ended prematurely");
		}
		++cline;
		if(cline!=0 && cline%RecLine==0){
			printf("\nLines analysed so far: %d\n",cline);
		}
	}
	printf("\n...total lines analysed: %d\n\n--\n",cline);
	nsim = cline;

	if(nsim*(nparam+nstats)>MAXDATA)
		printerr(".dat file is too big to be analysed");

    nclosest = (int) nsim*tol;

    /*writting first part of output to .txt file*/
    fprintf(out_inf,"Rejection - Mark Beaumont & Joao Lopes\n\n");
    fprintf(out_inf," last updated 01/05/09\n\n");
    fprintf(out_inf,"INPUT AND STARTING INFORMATION\n");
    fprintf(out_inf,"-------------------------------\n");
    fprintf(out_inf,"Data file:                %s\n",argv[1]);
    fprintf(out_inf,"Target file:              %s\n",argv[2]);
    fprintf(out_inf,"\nNumber of used parameters:  %d\n",nparam);
    fprintf(out_inf,"Number of used summstats: %d\n",nstats);
    fprintf(out_inf,"\nOUTPUT FILES:\n");
    fprintf(out_inf,"-------------------------------\n");
    fprintf(out_inf,"Output file:            %s.rej\n",argv[3]);
    fprintf(out_inf,"Priors file:            %s.pri\n",argv[3]);
    fprintf(out_inf,"\nRUN INFORMATION\n");
    fprintf(out_inf,"---------------------------\n");
    fprintf(out_inf,"Number of points analysed: %d\n",nsim);
    fprintf(out_inf,"relative tolerance: %g\n",tol);
    fprintf(out_inf,"absolute tolerance: %d\n\n",nclosest);
    fprintf(out_inf,"             average      stdev        target \n\n");

	printf("Getting the average and the stdev of the sumStats\n");
    printf("           average    stdev      target \n\n");

	/*fill stats[][]*/
	stats = (double **)malloc(nstats*sizeof(double *));
	for(cstat=0;cstat<nstats;++cstat)
		stats[cstat] = (double *)malloc(nsim*sizeof(double));
	fseek(out_dat,0,SEEK_SET);
	for(csim=0;csim<nsim;++csim){
		for(cparam=0 ; cparam<nparam ; ++cparam){
			fscanf(out_dat,"%lf",&d1);
		}
		for(cstat=0 ; cstat<nstats ; ++cstat){
			fscanf(out_dat,"%lf",&stats[cstat][csim]);
		}
	}

	/*get the average and the variance of each summStat and normalize it ((X-average)/sdev)*/
	/*normalize target summary statistics*/
	ave = (double *)malloc(nstats*sizeof(double));
	sdev = (double *)malloc(nstats*sizeof(double));
	for(cstat=0;cstat<nstats;++cstat){
		mom(stats[cstat],nsim,&ave[cstat],&sdev[cstat],&d3,&d4,&d5,&d6);
		printf("\nsumstat%d  %lf %lf    %lf\n",cstat+1,ave[cstat],sdev[cstat],trgStat[cstat]);
		fprintf(out_inf,"\nsumstat%d  %10.4g %10.4g    %10.4g\n",cstat+1,ave[cstat],sdev[cstat],trgStat[cstat]);
        for(csim=0;csim<nsim;++csim){
			if(sdev[cstat]==0)
				stats[cstat][csim]=0;
			else
				stats[cstat][csim] = (stats[cstat][csim] - ave[cstat])/sdev[cstat];
		}
		if(sdev[cstat]==0)
			trgStat[cstat]=0;
		else
			trgStat[cstat] = (trgStat[cstat] - ave[cstat])/sdev[cstat];
	}

	printf("\n--\n");
	printf("Ordering and choosing the closest simulated sumStats (%d) to the \"real\" sumStats:\n",nclosest);
	printf("\nObtaining the points distances...");
	fflush(NULL);

	/*get the Euclidian distance between each simulated summary statistics and the target(empirical) summary statistics */
	dist = (double *)malloc(nsim*sizeof(double));
	for(csim=0;csim<nsim;++csim){
		dist[csim] = 0;
		for(cstat=0;cstat<nstats;++cstat){
			dist[csim] += (stats[cstat][csim] - trgStat[cstat])*(stats[cstat][csim] - trgStat[cstat]);
		}
		dist[csim] = sqrt(dist[csim]);
	}

	/*sort the elements position of dist to index*/
	index = (int *)malloc(nsim*sizeof(int));
	dsorti('a',nsim,dist,index);

    /*free memory 1*/
    free(trgStat);
    free(dist);

	printf(" done");
	printf("\nAllocate the points information in the memory...");
	fflush(NULL);

	/*get the sumstats of the closest points*/
	fstats = (double **)malloc(nstats*sizeof(double *));
	for(cstat=0;cstat<nstats;++cstat)
		fstats[cstat] = (double *)malloc(nclosest*sizeof(double));
	for(csim=0 ; csim<nclosest ; ++csim){
		isim = index[csim];
		for(cstat=0;cstat<nstats;++cstat)
			fstats[cstat][csim] = stats[cstat][isim]*sdev[cstat] + ave[cstat];
	}

	/*free memory 2*/
	for(cstat=0;cstat<nstats;++cstat)
		free(stats[cstat]);
	free(stats);

	/*fill params[][]*/
	params = (double **)malloc(nparam*sizeof(double *));
	for(cparam=0;cparam<nparam;++cparam)
		params[cparam] = (double *)malloc(nsim*sizeof(double));
	fseek(out_dat,0,SEEK_SET);
	for(csim=0;csim<nsim;++csim){
		for(cparam=0 ; cparam<nparam ; ++cparam){
			fscanf(out_dat,"%lf",&params[cparam][csim]);
		}
		for(cstat=0 ; cstat<nstats ; ++cstat){
			fscanf(out_dat,"%lf",&d1);
		}
	}

	/*fill iparam[][]*/
	if(MAXNPRIOR > nsim)
		npriors = nsim;
	else
		npriors = MAXNPRIOR;
	iparams = (double **)malloc(nparam*sizeof(double *));
	for(cparam=0;cparam<nparam;++cparam)
		iparams[cparam] = (double *)malloc(npriors*sizeof(double));
	fseek(out_dat,0,SEEK_SET);

	gap=nsim/npriors;
	for(i=0,csim=0; csim<npriors; ++csim){
		i = csim*gap;
		for(cparam=0 ; cparam<nparam ; ++cparam){
			iparams[cparam][csim] = params[cparam][i];
		}
	}

	/*prints a random sample of the parameters to a file*/
	for(csim=0;csim<npriors;++csim){
		for(cparam=0;cparam<nparam;++cparam)
			fprintf(out_pri,"%g ",iparams[cparam][csim]);
		fprintf(out_pri,"\n");
	}

	/*free memory 3*/
	for(cparam=0;cparam<nparam;++cparam)
		free(iparams[cparam]);
	free(iparams);

    printf(" done");
    printf("\nGetting the closest points to the target file...");
    fflush(NULL);

	/*get the params of the closest points*/
	fparams = (double **)malloc(nparam*sizeof(double *));
	for(cparam=0;cparam<nparam;++cparam)
		fparams[cparam] = (double *)malloc(nclosest*sizeof(double));
	for(csim=0 ; csim<nclosest ; ++csim){
		isim = index[csim];
		for(cparam=0;cparam<nparam;++cparam)
			fparams[cparam][csim] = params[cparam][isim];
	}

	/*free memory 4*/
	free(ave);
	free(sdev);
	free(index);
	for(cparam=0;cparam<nparam;++cparam)
		free(params[cparam]);
	free(params);

	printf(" done");
	printf("\nPrinting the points to a file...");
	fflush(NULL);

	/*prints the nearest series of simulated parameters to the empirical data to a file*/
	for(csim=0 ; csim<nclosest ; ++csim){
		for(cparam=0;cparam<nparam;++cparam){
			fprintf(out_rej,"%g ", fparams[cparam][csim]);
		}
		for(cstat=0;cstat<nstats;++cstat){
			fprintf(out_rej,"%g ",fstats[cstat][csim]);
		}
		fprintf(out_rej,"\n");
	}

	printf(" done\n\n--\n\n");

	/*free memory 5*/
	for(cstat=0;cstat<nstats;++cstat)
		free(fstats[cstat]);
	free(fstats);
	for(cparam=0;cparam<nparam;++cparam)
		free(fparams[cparam]);
	free(fparams);

    /*get ending time and save it to .txt file*/
    fprintf(out_inf,"\nStarting date: %s",asctime(startTime));
    time( &endClock );                  // Get time in seconds
    endTime = localtime( &endClock );   // Convert time to struct tm form 
    fprintf(out_inf,"Ending date:   %s\n",asctime(endTime));
    fprintf(out_inf,"\nEND OF OUTPUT\n");

	/*close files*/
	fclose(out_dat);
	fclose(out_trg);
    fclose(out_inf);
	fclose(out_rej);
	fclose(out_pri);

}	//end of main
