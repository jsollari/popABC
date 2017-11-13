/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		13th November 2008
*/
#include "interface.h"

/*
	This function fills the topology related variables
	
	@param pm   - structer with the parameters values
	@param nevt - number of time events
	@param ntop - type of topology chosen
	@param npop - number of populations
*/
int choosetop(struct params *pm,int nevt,int ntop,int npop);

int sampPriors(struct params *pm,char *outline1,FILE *out_mut,FILE *out_rec,int printMut,int printRec,int foundSTR,int foundSNP,char *ltype){
	int cpop,cevt,cloc,i,			//iterators
		error,						//gets the error type
		npop,						//number of populations
		nevt,						//number of events
		ntop,						//number of possible topologies
		nloc;						//number of loci
	double nmig,					//number of migrants for each tevent
		   mut_meanM,				//mutation prior parameter
		   mut_devM,				//mutation prior parameter
		   mut_meanS,				//mutation prior parameter
		   mut_devS,				//mutation prior parameter
		   rec_meanM,				//mutation prior parameter
		   rec_devM,				//mutation prior parameter
		   rec_meanS,				//mutation prior parameter
		   rec_devS,				//mutation prior parameter
		   *unsorted_tev;			//unsorted times of event (prior type 3 and 4)

	npop = pm->npop;
	nevt = pm->nevt;
	ntop = pm->ntop;
	nloc = pm->nloc;

	//Set topology parameters from priors
	if(npop>1){
		error = choosetop(pm,nevt,ntop,npop);
		if(error)
			return error;
	}
	else
		pm->topol = 0;
					
	/*Set Ne parameters from priors*/
	for(cevt=0;cevt<=nevt;++cevt){
		for(cpop=0;cpop<npop;++cpop){
			//uniformly between 2 limits
			if(P_psize[cevt][cpop].type==1){
				pm->psize[cevt][cpop] =  gfsr8()*(P_psize[cevt][cpop].p[1] - P_psize[cevt][cpop].p[0]) + P_psize[cevt][cpop].p[0];
			}
			//use generalized gamma
			else if(P_psize[cevt][cpop].type == 2){
				pm->psize[cevt][cpop] =  rgengamma(P_psize[cevt][cpop].p[0], P_psize[cevt][cpop].p[1], P_psize[cevt][cpop].p[2], P_psize[cevt][cpop].p[3]);
			}
			// each anc pop is just done once
			if(cevt != 0)
				break;
		}
	}
	if(npop>1){
		/*Set tev parameters from priors*/
		if(P_t[0].type == 1 || P_t[0].type == 2){
			for(cevt=0;cevt<nevt;++cevt){
				//uniformly distributed
				if(P_t[cevt].type == 1)
					pm->tev[cevt] = gfsr8()*(P_t[cevt].p[1]-P_t[cevt].p[0])+P_t[cevt].p[0];	
				//generalized gamma distribution
				else if(P_t[cevt].type == 2)
					pm->tev[cevt] = rgengamma(P_t[cevt].p[0],P_t[cevt].p[1],P_t[cevt].p[2],P_t[cevt].p[3]);
				//to define next tev we allways sum the last one)
				if(cevt > 0) 
					pm->tev[cevt] = pm->tev[cevt-1] + pm->tev[cevt];
			}
		}
		else{
			unsorted_tev = (double*)malloc(nevt*sizeof(double));
			for(cevt=0;cevt<nevt;++cevt){
				//uniformly distributed (all tev)
				if(P_t[0].type == 3)
					unsorted_tev[cevt] = gfsr8()*(P_t[0].p[1]-P_t[0].p[0])+P_t[0].p[0];	
				//generalized gamma distribution (all tev)
				else if(P_t[0].type == 4)
					unsorted_tev[cevt] = rgengamma(P_t[0].p[0],P_t[0].p[1],P_t[0].p[2],P_t[0].p[3]);
			}
			//sort the sampled tev (type 3 and 4)
			shell_sort_double(unsorted_tev,nevt);
			for(cevt=0;cevt<nevt;++cevt){
				pm->tev[cevt] = unsorted_tev[cevt];
			}
		}
		//transform time(years) into time(generation)
		for(cevt=0;cevt<nevt;++cevt)
			pm->tev[cevt] /= pm->gent;
		//check for errors in the time values
		for(cevt=1 ; cevt<nevt ; ++cevt)
			if(pm->tev[cevt] <= pm->tev[cevt-1])
				return 28; //times not monotonic increasing
		/*Set migration parameters from priors*/
		for(cevt=0;cevt<nevt;++cevt){
			for(cpop=0;cpop<npop;++cpop){
				//no migration
				if (P_mig[cevt][cpop].type == 0)
					pm->mig[cevt][cpop] = 0;
				//uniform between 2 limits
				else if(P_mig[cevt][cpop].type == 1) 
					pm->mig[cevt][cpop] = gfsr8()*(P_mig[cevt][cpop].p[1]-P_mig[cevt][cpop].p[0]) + P_mig[cevt][cpop].p[0];
				//generalized gamma
				else if(P_mig[cevt][cpop].type == 2)
					pm->mig[cevt][cpop] = rgengamma(P_mig[cevt][cpop].p[0],P_mig[cevt][cpop].p[1],P_mig[cevt][cpop].p[2],P_mig[cevt][cpop].p[3]);
				//uniform between 2 limits on nmig
				else if(P_mig[cevt][cpop].type == 3){
					nmig = P_mig[cevt][cpop].p[0] + gfsr8()*(P_mig[cevt][cpop].p[1] - P_mig[cevt][cpop].p[0]);
					pm->mig[cevt][cpop] = nmig/pm->psize[cevt][cpop];
				}
				//generalized gamma on nmig
				else if(P_mig[cevt][cpop].type == 4){
					nmig = rgengamma(P_mig[cevt][cpop].p[0],P_mig[cevt][cpop].p[1],P_mig[cevt][cpop].p[2],P_mig[cevt][cpop].p[3]);
					pm->mig[cevt][cpop] = nmig/pm->psize[cevt][cpop];
				}
				// each anc pop is just done once
				if(cevt>0)
					break;
			}
		}
	}
	/*Set mutation parameters*/
	mut_meanM = 0.0;
	mut_devM = 0.0;
	mut_meanS = 0.0;
	mut_devS = 0.0;
	if(P_mutM.type!=0){
		mut_meanM = norm8()*P_mutM.p[1]+P_mutM.p[0];				
		mut_devM = fabs(norm8()*P_mutM.p[3]+P_mutM.p[2]);
	}
	if(P_mutS.type!=0){
		mut_meanS = norm8()*P_mutS.p[1]+P_mutS.p[0];				
		mut_devS = fabs(norm8()*P_mutS.p[3]+P_mutS.p[2]);		
	}
	for(cloc=0; cloc<nloc; ++cloc){
		if(ltype[cloc]=='m'||ltype[cloc]=='M'){
			if(P_mutM.type==0){
				pm->mu[cloc]=0.0;
			}
			else if(P_mutM.type==1){
				pm->mu[cloc] = pow(10.0,norm8()*mut_devM + mut_meanM);	
			}
            else if(P_mutM.type==2){
                pm->mu[cloc] = fabs(norm8()*mut_devM + mut_meanM);  
            }
		}
		else{
			if(P_mutS.type==0){
				pm->mu[cloc]=0.0;
			}
            else if(P_mutS.type==1){
                pm->mu[cloc] = pow(10.0,norm8()*mut_devS + mut_meanS);  
            }
            else if(P_mutS.type==2){
                pm->mu[cloc] = fabs(norm8()*mut_devS + mut_meanS);  
            }
		}
	}
	//print mutation rates into .mut file
	if(printMut){
		for(cloc=0;cloc<nloc;++cloc)
			fprintf(out_mut,"%e ",pm->mu[cloc]);
		fprintf(out_mut,"\n");
	}
	/*Set recombination parameters*/
	rec_meanM = 0.0;
	rec_devM = 0.0;
	rec_meanS = 0.0;
	rec_devS = 0.0;
	if(P_recM.type!=0){
		rec_meanM = norm8()*P_recM.p[1]+P_recM.p[0];				
		rec_devM = fabs(norm8()*P_recM.p[3]+P_recM.p[2]);
	}
	if(P_recS.type!=0){
		rec_meanS = norm8()*P_recS.p[1]+P_recS.p[0];				
		rec_devS = fabs(norm8()*P_recS.p[3]+P_recS.p[2]);		
	}
	for(cloc=0; cloc<nloc; ++cloc){
		if(ltype[cloc]=='m'||ltype[cloc]=='M'){
			if(P_recM.type==0){
				pm->rec[cloc]=0.0;
			}
			else if(P_recM.type==1){
				pm->rec[cloc] = pow(10.0,norm8()*rec_devM + rec_meanM);	
			}
            else if(P_recM.type==2){
                pm->rec[cloc] = fabs(norm8()*rec_devM + rec_meanM);  
            }
		}
		else{
			if(P_recS.type==0){
				pm->rec[cloc]=0.0;
			}
            else if(P_recS.type==1){
                pm->rec[cloc] = pow(10.0,norm8()*rec_devS + rec_meanS);  
            }
            else if(P_recS.type==2){
                pm->rec[cloc] = fabs(norm8()*rec_devS + rec_meanS);  
            }
		}
	}
	//print recombination rates into .mut file
	if(printRec){
		for(cloc=0;cloc<nloc;++cloc)
			fprintf(out_rec,"%e ",pm->rec[cloc]);
		fprintf(out_rec,"\n");
	}

	/*write parameters to output string*/	
	sprintf(outline1,"%d ",pm->topol);											//topology 
	if(P_top.type==3)
		sprintf(outline1+strlen(outline1),"%g ",P_top.p[0]);					//marker 
	else if(P_top.type==4)
		sprintf(outline1+strlen(outline1),"%g ",P_top.p[1]);					//marker 
	else if(P_top.type==5)
		sprintf(outline1+strlen(outline1),"%g ",P_top.p[2*nevt]);			    //marker 
	if(foundSTR)
		sprintf(outline1+strlen(outline1),"%g %g ",mut_meanM,mut_devM);	        //logmu_meanM logmu_devM
	if(foundSNP)
		sprintf(outline1+strlen(outline1),"%g %g ",mut_meanS,mut_devS);	        //logmu_meanS logmu_devS
	if(foundSTR)
		sprintf(outline1+strlen(outline1),"%g %g ",rec_meanM,rec_devM);	        //logmu_meanM logmu_devM
	if(foundSNP)
		sprintf(outline1+strlen(outline1),"%g %g ",rec_meanS,rec_devS);	        //logmu_meanS logmu_devS
	for(cevt=0;cevt<nevt;++cevt)
		sprintf(outline1+strlen(outline1),"%g ",pm->tev[cevt]*pm->gent);		//tev1*gent ...
	for(cevt=0;cevt<=nevt;++cevt){
		for(cpop=0;cpop<npop;++cpop){
			sprintf(outline1+strlen(outline1),"%g ",pm->psize[cevt][cpop]);		//ne1 ne2 ... nanc
			if(cevt!= 0)
				break;
		}
	}
	for(cevt=0;cevt<nevt;++cevt){
		for(cpop=0;cpop<npop;++cpop){
			sprintf(outline1+strlen(outline1),"%g ",pm->mig[cevt][cpop]);		// mig1 mig2 ... miganc
			if(cevt!= 0)
				break;
		}
	}

	/*free stuff*/
	if(npop>1){
		if(P_t[0].type == 3 || P_t[0].type == 4)
			free(unsorted_tev);
	}

	return 0;
	
}	// end of sample_from_priors

int choosetop(struct params *pm,int nevt,int ntop,int npop){
	int i;
	static int seq1[]  ={1,0,2,0,3,0,4,0},	//sequence of splitting no.1		
			   seq2[]  ={2,0,1,0,3,0,4,0},	//sequence of splitting no.2		
			   seq3[]  ={2,1,0,1,3,1,4,1},	//sequence of splitting no.3		
			   seq4[]  ={1,0,3,0,2,0,4,0},	//sequence of splitting no.4		
			   seq5[]  ={2,0,3,0,1,0,4,0},	//sequence of splitting no.5		
			   seq6[]  ={3,0,2,0,1,0,4,0},	//sequence of splitting no.6		
			   seq7[]  ={3,0,1,0,2,0,4,0},	//sequence of splitting no.7		
			   seq8[]  ={2,1,3,1,0,1,4,1},	//sequence of splitting no.8		
			   seq9[]  ={3,1,2,1,0,1,4,1},	//sequence of splitting no.9		
			   seq10[] ={3,1,0,1,2,1,4,1},	//sequence of splitting no.10		
			   seq11[] ={3,2,0,2,1,2,4,2},	//sequence of splitting no.11		
			   seq12[] ={3,2,1,2,0,2,4,2},	//sequence of splitting no.12		
			   seq13[] ={1,0,3,2,2,0,4,0},	//sequence of splitting no.13		
			   seq14[] ={2,0,3,1,1,0,4,0},	//sequence of splitting no.14		
			   seq15[] ={3,0,2,1,1,0,4,0},	//sequence of splitting no.15		
			   seq16[] ={2,3,0,1,1,3,4,3},	//sequence of splitting no.16		
			   seq17[] ={1,3,0,2,2,3,4,3},	//sequence of splitting no.17		
			   seq18[] ={1,2,0,3,2,3,4,3},	//sequence of splitting no.18		
			   seq19[] ={4,0,1,0,2,0,3,0},	//sequence of splitting no.19		
			   seq20[] ={4,1,1,0,2,0,3,0},	//sequence of splitting no.20		
			   seq21[] ={4,2,1,0,2,0,3,0},	//sequence of splitting no.21		
			   seq22[] ={4,3,1,0,2,0,3,0},	//sequence of splitting no.22		
			   seq23[] ={1,0,4,0,2,0,3,0},	//sequence of splitting no.23		
			   seq24[] ={1,0,4,2,2,0,3,0},	//sequence of splitting no.24		
			   seq25[] ={1,0,4,3,2,0,3,0},	//sequence of splitting no.25		
			   seq26[] ={1,0,2,0,4,0,3,0},	//sequence of splitting no.26		
			   seq27[] ={1,0,2,0,4,3,3,0},	//sequence of splitting no.27		
			   seq28[] ={4,0,2,0,1,0,3,0},	//sequence of splitting no.28		
			   seq29[] ={4,1,2,0,1,0,3,0},	//sequence of splitting no.29		
			   seq30[] ={4,2,2,0,1,0,3,0},	//sequence of splitting no.30		
			   seq31[] ={4,3,2,0,1,0,3,0},	//sequence of splitting no.31		
			   seq32[] ={2,0,4,0,1,0,3,0},	//sequence of splitting no.32		
			   seq33[] ={2,0,4,1,1,0,3,0},	//sequence of splitting no.33		
			   seq34[] ={2,0,4,3,1,0,3,0},	//sequence of splitting no.34		
			   seq35[] ={2,0,1,0,4,0,3,0},	//sequence of splitting no.35		
			   seq36[] ={2,0,1,0,4,3,3,0},	//sequence of splitting no.36		
			   seq37[] ={4,0,2,1,0,1,3,1},	//sequence of splitting no.37		
			   seq38[] ={4,1,2,1,0,1,3,1},	//sequence of splitting no.38		
			   seq39[] ={4,2,2,1,0,1,3,1},	//sequence of splitting no.39		
			   seq40[] ={4,3,2,1,0,1,3,1},	//sequence of splitting no.40		
			   seq41[] ={2,1,4,0,0,1,3,1},	//sequence of splitting no.41		
			   seq42[] ={2,1,4,1,0,1,3,1},	//sequence of splitting no.42		
			   seq43[] ={2,1,4,3,0,1,3,1},	//sequence of splitting no.43		
			   seq44[] ={2,1,0,1,4,1,3,1},	//sequence of splitting no.44		
			   seq45[] ={2,1,0,1,4,3,3,1},	//sequence of splitting no.45		
			   seq46[] ={4,0,1,0,3,0,2,0},	//sequence of splitting no.46		
			   seq47[] ={4,1,1,0,3,0,2,0},	//sequence of splitting no.47		
			   seq48[] ={4,2,1,0,3,0,2,0},	//sequence of splitting no.48		
			   seq49[] ={4,3,1,0,3,0,2,0},	//sequence of splitting no.49		
			   seq50[] ={1,0,4,0,3,0,2,0},	//sequence of splitting no.50		
			   seq51[] ={1,0,4,2,3,0,2,0},	//sequence of splitting no.51		
			   seq52[] ={1,0,4,3,3,0,2,0},	//sequence of splitting no.52		
			   seq53[] ={1,0,3,0,4,0,2,0},	//sequence of splitting no.53		
			   seq54[] ={1,0,3,0,4,2,2,0},	//sequence of splitting no.54		
			   seq55[] ={4,0,2,0,3,0,1,0},	//sequence of splitting no.55		
			   seq56[] ={4,1,2,0,3,0,1,0},	//sequence of splitting no.56		
			   seq57[] ={4,2,2,0,3,0,1,0},	//sequence of splitting no.57		
			   seq58[] ={4,3,2,0,3,0,1,0},	//sequence of splitting no.58		
			   seq59[] ={2,0,4,0,3,0,1,0},	//sequence of splitting no.59		
			   seq60[] ={2,0,4,1,3,0,1,0},	//sequence of splitting no.60		
			   seq61[] ={2,0,4,3,3,0,1,0},	//sequence of splitting no.61		
			   seq62[] ={2,0,3,0,4,0,1,0},	//sequence of splitting no.62		
			   seq63[] ={2,0,3,0,4,1,1,0},	//sequence of splitting no.63		
			   seq64[] ={4,0,3,0,2,0,1,0},	//sequence of splitting no.64		
			   seq65[] ={4,1,3,0,2,0,1,0},	//sequence of splitting no.65		
			   seq66[] ={4,2,3,0,2,0,1,0},	//sequence of splitting no.66		
			   seq67[] ={4,3,3,0,2,0,1,0},	//sequence of splitting no.67		
			   seq68[] ={3,0,4,0,2,0,1,0},	//sequence of splitting no.68		
			   seq69[] ={3,0,4,1,2,0,1,0},	//sequence of splitting no.69		
			   seq70[] ={3,0,4,2,2,0,1,0},	//sequence of splitting no.70		
			   seq71[] ={3,0,2,0,4,0,1,0},	//sequence of splitting no.71		
			   seq72[] ={3,0,2,0,4,1,1,0},	//sequence of splitting no.72		
			   seq73[] ={4,0,3,0,1,0,2,0},	//sequence of splitting no.73		
			   seq74[] ={4,1,3,0,1,0,2,0},	//sequence of splitting no.74		
			   seq75[] ={4,2,3,0,1,0,2,0},	//sequence of splitting no.75		
			   seq76[] ={4,3,3,0,1,0,2,0},	//sequence of splitting no.76		
			   seq77[] ={3,0,4,0,1,0,2,0},	//sequence of splitting no.77		
			   seq78[] ={3,0,4,1,1,0,2,0},	//sequence of splitting no.78		
			   seq79[] ={3,0,4,2,1,0,2,0},	//sequence of splitting no.79		
			   seq80[] ={3,0,1,0,4,0,2,0},	//sequence of splitting no.80		
			   seq81[] ={3,0,1,0,4,2,2,0},	//sequence of splitting no.81		
			   seq82[] ={4,0,2,1,3,1,0,1},	//sequence of splitting no.82		
			   seq83[] ={4,1,2,1,3,1,0,1},	//sequence of splitting no.83		
			   seq84[] ={4,2,2,1,3,1,0,1},	//sequence of splitting no.84		
			   seq85[] ={4,3,2,1,3,1,0,1},	//sequence of splitting no.85		
			   seq86[] ={2,1,4,0,3,1,0,1},	//sequence of splitting no.86		
			   seq87[] ={2,1,4,1,3,1,0,1},	//sequence of splitting no.87		
			   seq88[] ={2,1,4,3,3,1,0,1},	//sequence of splitting no.88		
			   seq89[] ={2,1,3,1,4,0,0,1},	//sequence of splitting no.89		
			   seq90[] ={2,1,3,1,4,1,0,1},	//sequence of splitting no.90		
			   seq91[] ={4,0,3,1,2,1,0,1},	//sequence of splitting no.91		
			   seq92[] ={4,1,3,1,2,1,0,1},	//sequence of splitting no.92		
			   seq93[] ={4,2,3,1,2,1,0,1},	//sequence of splitting no.93		
			   seq94[] ={4,3,3,1,2,1,0,1},	//sequence of splitting no.94		
			   seq95[] ={3,1,4,0,2,1,0,1},	//sequence of splitting no.95		
			   seq96[] ={3,1,4,1,2,1,0,1},	//sequence of splitting no.96		
			   seq97[] ={3,1,4,2,2,1,0,1},	//sequence of splitting no.97		
			   seq98[] ={3,1,2,1,4,0,0,1},	//sequence of splitting no.98		
			   seq99[] ={3,1,2,1,4,1,0,1},	//sequence of splitting no.99		
			   seq100[]={4,0,3,1,0,1,2,1},	//sequence of splitting no.100		
			   seq101[]={4,1,3,1,0,1,2,1},	//sequence of splitting no.101		
			   seq102[]={4,2,3,1,0,1,2,1},	//sequence of splitting no.102		
			   seq103[]={4,3,3,1,0,1,2,1},	//sequence of splitting no.103		
			   seq104[]={3,1,4,0,0,1,2,1},	//sequence of splitting no.104		
			   seq105[]={3,1,4,1,0,1,2,1},	//sequence of splitting no.105		
			   seq106[]={3,1,4,2,0,1,2,1},	//sequence of splitting no.106		
			   seq107[]={3,1,0,1,4,1,2,1},	//sequence of splitting no.107		
			   seq108[]={3,1,0,1,4,2,2,1},	//sequence of splitting no.108		
			   seq109[]={4,0,3,2,0,2,1,2},	//sequence of splitting no.109		
			   seq110[]={4,1,3,2,0,2,1,2},	//sequence of splitting no.110		
			   seq111[]={4,2,3,2,0,2,1,2},	//sequence of splitting no.111		
			   seq112[]={4,3,3,2,0,2,1,2},	//sequence of splitting no.112		
			   seq113[]={3,2,4,0,0,2,1,2},	//sequence of splitting no.113		
			   seq114[]={3,2,4,1,0,2,1,2},	//sequence of splitting no.114		
			   seq115[]={3,2,4,2,0,2,1,2},	//sequence of splitting no.115		
			   seq116[]={3,2,0,2,4,1,1,2},	//sequence of splitting no.116		
			   seq117[]={3,2,0,2,4,2,1,2},	//sequence of splitting no.117		
			   seq118[]={4,0,3,2,1,2,0,2},	//sequence of splitting no.118		
			   seq119[]={4,1,3,2,1,2,0,2},	//sequence of splitting no.119		
			   seq120[]={4,2,3,2,1,2,0,2},	//sequence of splitting no.120		
			   seq121[]={4,3,3,2,1,2,0,2},	//sequence of splitting no.121		
			   seq122[]={3,2,4,0,1,2,0,2},	//sequence of splitting no.122		
			   seq123[]={3,2,4,1,1,2,0,2},	//sequence of splitting no.123		
			   seq124[]={3,2,4,2,1,2,0,2},	//sequence of splitting no.124		
			   seq125[]={3,2,1,2,4,0,0,2},	//sequence of splitting no.125		
			   seq126[]={3,2,1,2,4,2,0,2},	//sequence of splitting no.126		
			   seq127[]={4,0,1,0,3,2,2,0},	//sequence of splitting no.127		
			   seq128[]={4,1,1,0,3,2,2,0},	//sequence of splitting no.128		
			   seq129[]={4,2,1,0,3,2,2,0},	//sequence of splitting no.129		
			   seq130[]={4,3,1,0,3,2,2,0},	//sequence of splitting no.130		
			   seq131[]={1,0,4,0,3,2,2,0},	//sequence of splitting no.131		
			   seq132[]={1,0,4,2,3,2,2,0},	//sequence of splitting no.132		
			   seq133[]={1,0,4,3,3,2,2,0},	//sequence of splitting no.133		
			   seq134[]={1,0,3,2,4,0,2,0},	//sequence of splitting no.134		
			   seq135[]={1,0,3,2,4,2,2,0},	//sequence of splitting no.135		
			   seq136[]={4,0,2,0,3,1,1,0},	//sequence of splitting no.136		
			   seq137[]={4,1,2,0,3,1,1,0},	//sequence of splitting no.137		
			   seq138[]={4,2,2,0,3,1,1,0},	//sequence of splitting no.138		
			   seq139[]={4,3,2,0,3,1,1,0},	//sequence of splitting no.139		
			   seq140[]={2,0,4,0,3,1,1,0},	//sequence of splitting no.140		
			   seq141[]={2,0,4,1,3,1,1,0},	//sequence of splitting no.141		
			   seq142[]={2,0,4,3,3,1,1,0},	//sequence of splitting no.142		
			   seq143[]={2,0,3,1,4,0,1,0},	//sequence of splitting no.143		
			   seq144[]={2,0,3,1,4,1,1,0},	//sequence of splitting no.144		
			   seq145[]={4,0,3,0,2,1,1,0},	//sequence of splitting no.145		
			   seq146[]={4,1,3,0,2,1,1,0},	//sequence of splitting no.146		
			   seq147[]={4,2,3,0,2,1,1,0},	//sequence of splitting no.147		
			   seq148[]={4,3,3,0,2,1,1,0},	//sequence of splitting no.148		
			   seq149[]={3,0,4,0,2,1,1,0},	//sequence of splitting no.149		
			   seq150[]={3,0,4,1,2,1,1,0},	//sequence of splitting no.150		
			   seq151[]={3,0,4,2,2,1,1,0},	//sequence of splitting no.151		
			   seq152[]={3,0,2,1,4,0,1,0},	//sequence of splitting no.152		
			   seq153[]={3,0,2,1,4,1,1,0},	//sequence of splitting no.153		
			   seq154[]={4,0,2,3,0,1,1,3},	//sequence of splitting no.154		
			   seq155[]={4,1,2,3,0,1,1,3},	//sequence of splitting no.155		
			   seq156[]={4,2,2,3,0,1,1,3},	//sequence of splitting no.156		
			   seq157[]={4,3,2,3,0,1,1,3},	//sequence of splitting no.157		
			   seq158[]={2,3,4,0,0,1,1,3},	//sequence of splitting no.158		
			   seq159[]={2,3,4,1,0,1,1,3},	//sequence of splitting no.159		
			   seq160[]={2,3,4,3,0,1,1,3},	//sequence of splitting no.160		
			   seq161[]={2,3,0,1,4,1,1,3},	//sequence of splitting no.161		
			   seq162[]={2,3,0,1,4,3,1,3},	//sequence of splitting no.162		
			   seq163[]={4,0,1,3,0,2,2,3},	//sequence of splitting no.163		
			   seq164[]={4,1,1,3,0,2,2,3},	//sequence of splitting no.164		
			   seq165[]={4,2,1,3,0,2,2,3},	//sequence of splitting no.165		
			   seq166[]={4,3,1,3,0,2,2,3},	//sequence of splitting no.166		
			   seq167[]={1,3,4,0,0,2,2,3},	//sequence of splitting no.167		
			   seq168[]={1,3,4,2,0,2,2,3},	//sequence of splitting no.168		
			   seq169[]={1,3,4,3,0,2,2,3},	//sequence of splitting no.169		
			   seq170[]={1,3,0,2,4,2,2,3},	//sequence of splitting no.170		
			   seq171[]={1,3,0,2,4,3,2,3},	//sequence of splitting no.171		
			   seq172[]={4,0,1,2,0,3,2,3},	//sequence of splitting no.172		
			   seq173[]={4,1,1,2,0,3,2,3},	//sequence of splitting no.173		
			   seq174[]={4,2,1,2,0,3,2,3},	//sequence of splitting no.174		
			   seq175[]={4,3,1,2,0,3,2,3},	//sequence of splitting no.175		
			   seq176[]={1,2,4,0,0,3,2,3},	//sequence of splitting no.176		
			   seq177[]={1,2,4,2,0,3,2,3},	//sequence of splitting no.177		
			   seq178[]={1,2,4,3,0,3,2,3},	//sequence of splitting no.178		
			   seq179[]={1,2,0,3,4,2,2,3},	//sequence of splitting no.179		
			   seq180[]={1,2,0,3,4,3,2,3};	//sequence of splitting no.180		

	/*Set topology parameters from priors*/
	if((P_top.type==0 || P_top.type==1 ||P_top.type==3 || P_top.type==4) && npop < 6){
		if(P_top.type == 0 || P_top.type == 3)
			i = disrand(1,ntop);
		else if(P_top.type == 1 || P_top.type == 4)
			i = P_top.p[0];
		switch(i){
			case 1:   pm->topol = 1;   pm->seq = seq1;   break;
			case 2:   pm->topol = 2;   pm->seq = seq2;   break;
			case 3:   pm->topol = 3;   pm->seq = seq3;   break;
			case 4:	  pm->topol = 4;   pm->seq = seq4;   break;
			case 5:   pm->topol = 5;   pm->seq = seq5;   break;
			case 6:   pm->topol = 6;   pm->seq = seq6;   break;
			case 7:   pm->topol = 7;   pm->seq = seq7;   break;
			case 8:   pm->topol = 8;   pm->seq = seq8;   break;
			case 9:   pm->topol = 9;   pm->seq = seq9;   break;
			case 10:  pm->topol = 10;  pm->seq = seq10;  break;
			case 11:  pm->topol = 11;  pm->seq = seq11;  break;
			case 12:  pm->topol = 12;  pm->seq = seq12;  break;
			case 13:  pm->topol = 13;  pm->seq = seq13;  break;
			case 14:  pm->topol = 14;  pm->seq = seq14;  break;
			case 15:  pm->topol = 15;  pm->seq = seq15;  break;
			case 16:  pm->topol = 16;  pm->seq = seq16;  break;
			case 17:  pm->topol = 17;  pm->seq = seq17;  break;
			case 18:  pm->topol = 18;  pm->seq = seq18;  break;
			case 19:  pm->topol = 19;  pm->seq = seq19;  break;
			case 20:  pm->topol = 20;  pm->seq = seq20;  break;
			case 21:  pm->topol = 21;  pm->seq = seq21;  break;
			case 22:  pm->topol = 22;  pm->seq = seq22;  break;
			case 23:  pm->topol = 23;  pm->seq = seq23;  break;
			case 24:  pm->topol = 24;  pm->seq = seq24;  break;
			case 25:  pm->topol = 25;  pm->seq = seq25;  break;
			case 26:  pm->topol = 26;  pm->seq = seq26;  break;
			case 27:  pm->topol = 27;  pm->seq = seq27;  break;
			case 28:  pm->topol = 28;  pm->seq = seq28;  break;
			case 29:  pm->topol = 29;  pm->seq = seq29;  break;
			case 30:  pm->topol = 30;  pm->seq = seq30;  break;
			case 31:  pm->topol = 31;  pm->seq = seq31;  break;
			case 32:  pm->topol = 32;  pm->seq = seq32;  break;
			case 33:  pm->topol = 33;  pm->seq = seq33;  break;
			case 34:  pm->topol = 34;  pm->seq = seq34;  break;
			case 35:  pm->topol = 35;  pm->seq = seq35;  break;
			case 36:  pm->topol = 36;  pm->seq = seq36;  break;
			case 37:  pm->topol = 37;  pm->seq = seq37;  break;
			case 38:  pm->topol = 38;  pm->seq = seq38;  break;
			case 39:  pm->topol = 39;  pm->seq = seq39;  break;
			case 40:  pm->topol = 40;  pm->seq = seq40;  break;
			case 41:  pm->topol = 41;  pm->seq = seq41;  break;
			case 42:  pm->topol = 42;  pm->seq = seq42;  break;
			case 43:  pm->topol = 43;  pm->seq = seq43;  break;
			case 44:  pm->topol = 44;  pm->seq = seq44;  break;
			case 45:  pm->topol = 45;  pm->seq = seq45;  break;
			case 46:  pm->topol = 46;  pm->seq = seq46;  break;
			case 47:  pm->topol = 47;  pm->seq = seq47;  break;
			case 48:  pm->topol = 48;  pm->seq = seq48;  break;
			case 49:  pm->topol = 49;  pm->seq = seq49;  break;
			case 50:  pm->topol = 50;  pm->seq = seq50;  break;
			case 51:  pm->topol = 51;  pm->seq = seq51;  break;
			case 52:  pm->topol = 52;  pm->seq = seq52;  break;
			case 53:  pm->topol = 53;  pm->seq = seq53;  break;
			case 54:  pm->topol = 54;  pm->seq = seq54;  break;
			case 55:  pm->topol = 55;  pm->seq = seq55;  break;
			case 56:  pm->topol = 56;  pm->seq = seq56;  break;
			case 57:  pm->topol = 57;  pm->seq = seq57;  break;
			case 58:  pm->topol = 58;  pm->seq = seq58;  break;
			case 59:  pm->topol = 59;  pm->seq = seq59;  break;
			case 60:  pm->topol = 60;  pm->seq = seq60;  break;
			case 61:  pm->topol = 61;  pm->seq = seq61;  break;
			case 62:  pm->topol = 62;  pm->seq = seq62;  break;
			case 63:  pm->topol = 63;  pm->seq = seq63;  break;
			case 64:  pm->topol = 64;  pm->seq = seq64;  break;
			case 65:  pm->topol = 65;  pm->seq = seq65;  break;
			case 66:  pm->topol = 66;  pm->seq = seq66;  break;
			case 67:  pm->topol = 67;  pm->seq = seq67;  break;
			case 68:  pm->topol = 68;  pm->seq = seq68;  break;
			case 69:  pm->topol = 69;  pm->seq = seq69;  break;
			case 70:  pm->topol = 70;  pm->seq = seq70;  break;
			case 71:  pm->topol = 71;  pm->seq = seq71;  break;
			case 72:  pm->topol = 72;  pm->seq = seq72;  break;
			case 73:  pm->topol = 73;  pm->seq = seq73;  break;
			case 74:  pm->topol = 74;  pm->seq = seq74;  break;
			case 75:  pm->topol = 75;  pm->seq = seq75;  break;
			case 76:  pm->topol = 76;  pm->seq = seq76;  break;
			case 77:  pm->topol = 77;  pm->seq = seq77;  break;
			case 78:  pm->topol = 78;  pm->seq = seq78;  break;
			case 79:  pm->topol = 79;  pm->seq = seq79;  break;
			case 80:  pm->topol = 80;  pm->seq = seq80;  break;
			case 81:  pm->topol = 81;  pm->seq = seq81;  break;
			case 82:  pm->topol = 82;  pm->seq = seq82;  break;
			case 83:  pm->topol = 83;  pm->seq = seq83;  break;
			case 84:  pm->topol = 84;  pm->seq = seq84;  break;
			case 85:  pm->topol = 85;  pm->seq = seq85;  break;
			case 86:  pm->topol = 86;  pm->seq = seq86;  break;
			case 87:  pm->topol = 87;  pm->seq = seq87;  break;
			case 88:  pm->topol = 88;  pm->seq = seq88;  break;
			case 89:  pm->topol = 89;  pm->seq = seq89;  break;
			case 90:  pm->topol = 90;  pm->seq = seq90;  break;
			case 91:  pm->topol = 91;  pm->seq = seq91;  break;
			case 92:  pm->topol = 92;  pm->seq = seq92;  break;
			case 93:  pm->topol = 93;  pm->seq = seq93;  break;
			case 94:  pm->topol = 94;  pm->seq = seq94;  break;
			case 95:  pm->topol = 95;  pm->seq = seq95;  break;
			case 96:  pm->topol = 96;  pm->seq = seq96;  break;
			case 97:  pm->topol = 97;  pm->seq = seq97;  break;
			case 98:  pm->topol = 98;  pm->seq = seq98;  break;
			case 99:  pm->topol = 99;  pm->seq = seq99;  break;
			case 100: pm->topol = 100; pm->seq = seq100; break;
			case 101: pm->topol = 101; pm->seq = seq101; break;
			case 102: pm->topol = 102; pm->seq = seq102; break;
			case 103: pm->topol = 103; pm->seq = seq103; break;
			case 104: pm->topol = 104; pm->seq = seq104; break;
			case 105: pm->topol = 105; pm->seq = seq105; break;
			case 106: pm->topol = 106; pm->seq = seq106; break;
			case 107: pm->topol = 107; pm->seq = seq107; break;
			case 108: pm->topol = 108; pm->seq = seq108; break;
			case 109: pm->topol = 109; pm->seq = seq109; break;
			case 110: pm->topol = 110; pm->seq = seq110; break;
			case 111: pm->topol = 111; pm->seq = seq111; break;
			case 112: pm->topol = 112; pm->seq = seq112; break;
			case 113: pm->topol = 113; pm->seq = seq113; break;
			case 114: pm->topol = 114; pm->seq = seq114; break;
			case 115: pm->topol = 115; pm->seq = seq115; break;
			case 116: pm->topol = 116; pm->seq = seq116; break;
			case 117: pm->topol = 117; pm->seq = seq117; break;
			case 118: pm->topol = 118; pm->seq = seq118; break;
			case 119: pm->topol = 119; pm->seq = seq119; break;
			case 120: pm->topol = 120; pm->seq = seq120; break;
			case 121: pm->topol = 121; pm->seq = seq121; break;
			case 122: pm->topol = 122; pm->seq = seq122; break;
			case 123: pm->topol = 123; pm->seq = seq123; break;
			case 124: pm->topol = 124; pm->seq = seq124; break;
			case 125: pm->topol = 125; pm->seq = seq125; break;
			case 126: pm->topol = 126; pm->seq = seq126; break;
			case 127: pm->topol = 127; pm->seq = seq127; break;
			case 128: pm->topol = 128; pm->seq = seq128; break;
			case 129: pm->topol = 129; pm->seq = seq129; break;
			case 130: pm->topol = 130; pm->seq = seq130; break;
			case 131: pm->topol = 131; pm->seq = seq131; break;
			case 132: pm->topol = 132; pm->seq = seq132; break;
			case 133: pm->topol = 133; pm->seq = seq133; break;
			case 134: pm->topol = 134; pm->seq = seq134; break;
			case 135: pm->topol = 135; pm->seq = seq135; break;
			case 136: pm->topol = 136; pm->seq = seq136; break;
			case 137: pm->topol = 137; pm->seq = seq137; break;
			case 138: pm->topol = 138; pm->seq = seq138; break;
			case 139: pm->topol = 139; pm->seq = seq139; break;
			case 140: pm->topol = 140; pm->seq = seq140; break;
			case 141: pm->topol = 141; pm->seq = seq141; break;
			case 142: pm->topol = 142; pm->seq = seq142; break;
			case 143: pm->topol = 143; pm->seq = seq143; break;
			case 144: pm->topol = 144; pm->seq = seq144; break;
			case 145: pm->topol = 145; pm->seq = seq145; break;
			case 146: pm->topol = 146; pm->seq = seq146; break;
			case 147: pm->topol = 147; pm->seq = seq147; break;
			case 148: pm->topol = 148; pm->seq = seq148; break;
			case 149: pm->topol = 149; pm->seq = seq149; break;
			case 150: pm->topol = 150; pm->seq = seq150; break;
			case 151: pm->topol = 151; pm->seq = seq151; break;
			case 152: pm->topol = 152; pm->seq = seq152; break;
			case 153: pm->topol = 153; pm->seq = seq153; break;
			case 154: pm->topol = 154; pm->seq = seq154; break;
			case 155: pm->topol = 155; pm->seq = seq155; break;
			case 156: pm->topol = 156; pm->seq = seq156; break;
			case 157: pm->topol = 157; pm->seq = seq157; break;
			case 158: pm->topol = 158; pm->seq = seq158; break;
			case 159: pm->topol = 159; pm->seq = seq159; break;
			case 160: pm->topol = 160; pm->seq = seq160; break;
			case 161: pm->topol = 161; pm->seq = seq161; break;
			case 162: pm->topol = 162; pm->seq = seq162; break;
			case 163: pm->topol = 163; pm->seq = seq163; break;
			case 164: pm->topol = 164; pm->seq = seq164; break;
			case 165: pm->topol = 165; pm->seq = seq165; break;
			case 166: pm->topol = 166; pm->seq = seq166; break;
			case 167: pm->topol = 167; pm->seq = seq167; break;
			case 168: pm->topol = 168; pm->seq = seq168; break;
			case 169: pm->topol = 169; pm->seq = seq169; break;
			case 170: pm->topol = 170; pm->seq = seq170; break;
			case 171: pm->topol = 171; pm->seq = seq171; break;
			case 172: pm->topol = 172; pm->seq = seq172; break;
			case 173: pm->topol = 173; pm->seq = seq173; break;
			case 174: pm->topol = 174; pm->seq = seq174; break;
			case 175: pm->topol = 175; pm->seq = seq175; break;
			case 176: pm->topol = 176; pm->seq = seq176; break;
			case 177: pm->topol = 177; pm->seq = seq177; break;
			case 178: pm->topol = 178; pm->seq = seq178; break;
			case 179: pm->topol = 179; pm->seq = seq179; break;
			case 180: pm->topol = 180; pm->seq = seq180; break;
			default: return 29; //choosing branching sequence
		}
		return 0;
	}
	else if (P_top.type == 2 || P_top.type == 5){
		pm->topol = 0;
		pm->seq = (int*)malloc(2*nevt*sizeof(int));
		for(i=0 ; i<2*nevt ; i++)
			pm->seq[i]=P_top.p[i];
		return 0;
	}
	else
		return 30; //more than 5 pop, top as to be fixed using option 2 or 5
		
}
