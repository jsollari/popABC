/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		13th November 2008
*/
#include "interface.h"
					 
int makestats(char *output,int *lsstats){
	FILE *out_sst;	//pntr to .sst
	int cstat;		//iterators
	char *sst_name;				//name of the output file

	sst_name = malloc((5+strlen(output))*sizeof(char));
	strcpy(sst_name,output);
	out_sst = fopen(strcat(sst_name,".sst"),"w");
	if(out_sst == NULL)
		return 1; //couldn't create .sst file

	/*write the .sst file*/
	for(cstat=0; cstat<MAXSSTATS; cstat++){
		switch(cstat){
		case 0:
		fprintf(out_sst," Microsatellites data:\n");
		fprintf(out_sst," -heterozygosity                #%d\n",lsstats[cstat]);
		break;
		case 1:
		fprintf(out_sst," -variance of alleles length    #%d\n",lsstats[cstat]);
		break;
		case 2:
		fprintf(out_sst," -number of alleles             #%d\n",lsstats[cstat]);
		break;
		case 3:
		fprintf(out_sst," -curtosis of alleles length    #%d\n",lsstats[cstat]);
		break;
		case 4:
		fprintf(out_sst," -Shanon's index                #%d\n",lsstats[cstat]);
		break;
		case 5:
		fprintf(out_sst," -Nm estimator based on H       #%d\n",lsstats[cstat]);
		break;
		case 6:
		fprintf(out_sst,"\n Sequence data:\n");
		fprintf(out_sst," -mean of pairwise differences  #%d\n",lsstats[cstat]);
		break;
		case 7:
		fprintf(out_sst," -number of segregating sites   #%d\n",lsstats[cstat]);
		break;
		case 8:
		fprintf(out_sst," -number of haplotypes          #%d\n",lsstats[cstat]);
		break;
		case 9:
		fprintf(out_sst," -Shanon's index                #%d\n",lsstats[cstat]);
		break;
		case 10:
		fprintf(out_sst," -mean of MFS                   #%d\n",lsstats[cstat]);
		break;
		case 11:
		fprintf(out_sst," -stdev of MFS                  #%d\n",lsstats[cstat]);
		break;
		case 12:
		fprintf(out_sst," -Nm estimator based on S       #%d\n",lsstats[cstat]);
		break;
		case 13:
		fprintf(out_sst," -private segregating sites     #%d\n",lsstats[cstat]);
		break;
		case 14:
		fprintf(out_sst," -S(1)                          #%d\n",lsstats[cstat]);
		break;
		}
	}

	/*free stuff*/
	free(sst_name);
	
	/*close files*/
	fclose(out_sst);

	return 0;

} //end of main
