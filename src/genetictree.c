/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		17th October 2008
	
	NBB - based on Mark's splitmig.c and splitmigS.c
*/

#include "abc.h"

//struct that contains information about a node of a genetic tree
struct node{
	int spop;					//subpopulation no.
	int ospop;					//original pop. in sample node (meaningless in other nodes)
	double time;				//time of the coalescente event
	struct node *a[2];			//2 direct ancient node (for coalescent use a[0])
	struct node *d[2];			//2 direct descendents nodes
	/*only used when recombination is present*/
	int cut;					//place where the recombination occurs
	int cbeg; 					//leftmost ancestral site
	int cend; 					//rightmost ancestral site
	int *desc; 					//number of descendents in the sample in a locus
	/*only used with microsatellite data*/
	int *dnaM;					//size of microssatellites
	/*only used with sequence data*/
	double *dnaS;					//dna sequence with info about mut
	int len;					//current no. of mut on the branch of this node
	int lenMax;					//current max mut no. on the branch of this node	
};

/*
	Performs the simulation of a genetic tree.
	
	@param lind   - no. indiv. of a samp. pop. for a given locus
	@param nind   - no. indiv. of a samp. pop. for all loci
	@param mut 	  - mut. rate for a given loci
	@param rec 	  - rec. rate for a given loci
	@param data   - structure to store the data from the genetic tree
	@param cloc   - current locus
	@param type   - type of analysis
	@param maxRec - max no. sites when considering rec events
*/
void simulation(int lind[],int nind,double mut,double rec,struct data *data,int cloc,char type,int maxRec);

/*
	Join two pop. originated in a spliting event
	
	@param from_pop - pop. that is going to join the other
	@param to_pop   - pop. that is going to be joined
	@param currPop - pointer to the identifier of the population being evaluated
*/
void colonize(int from_pop,int to_pop,int *currPop);

/*
	Function that creates a recombination event

	@param cpop - pop. where the coalescente event happens
	@param type 	 - type of analysis
	@param nind - number of samples in the population
	@param time - time of the recombination event
	@param maxRec - max no. sites when considering rec events
*/
void recombination(int cpop,char type, int nind, double time,int *currNode,int*currIndiv,int maxRec);

/*
	Function that creates a coalescent event
	
	@param cpop 	 - pop. where the coalescente event happens
	@param type 	 - type of analysis
	@param nind 	 - number of samples in the population
	@param time		 - expected time for next event
	@param currNode  - pointer to the identifier of the node being evaluated
	@param currIndiv - pointer to the identifier of the individual being evaluated
	@param maxRec - max no. sites when considering rec events
*/
void coalescence(int cpop,char type,int nind,double time,int *currNode,int*currIndiv,int maxRec);

/*
	Function that creates a migration event
	
	@param ito       - pop. from where the migration occurs (as seen back in time)
	@param npop      - number of populations
	@param currPop   - pointer to the identifier of the population being evaluated
	@param currEvent - identifier of the splitting time event being evaluated
*/
void migration(int ito,int npop,int *currPop,int currEvent);

/*
	Adds information (mutations) to the created genetic trees
	
	@param p    - MRCA (most recent common ancester)
	@param data - place where the information about the genetic tree is stored
	@param cloc - current locus
	@param mu	- mut. rate of the given locus
	@param nind - number of sample on a population for the locus
	@param type - type of analysis 
	@param nmut - pointer to the number of mutations
	@param maxRec - max no. sites when considering rec events
	@param lmut - list of the mutations
*/
void treefill(struct node *p,struct data *data, int cloc, double mut,int nind,char type,int *nmut,int maxRec,double **lmut);

/*
	Fills the val[] and the freq[] of the struct data of a genetic tree (Microsatelites)
	
	@param p     - node from the sample
	@param data  - data struct related to the current genetic tree
	@param cloc  - current locus
	@param nrSTR - nr of linked STR's loci
*/
void makeSampM(struct node *p,struct data *data, int cloc, int nrSTR);

/*
	Fills the val[] and the freq[] of the struct data of a genetic tree (Segregating sites)
	
	@param p      - node from the sample
	@param data   - data struct related to the current genetic tree
	@param cloc   - current locus
	@param nmut   - number of mutations
	@param lmut   - list of the mutations (only used if recombination is present
	@param maxRec - check if recombination is present
*/
void makeSampS(struct node *p,struct data *data, int cloc, int nmut,double *lmut,int maxRec);

/*
	Sort the freq[] and the val[] by microsatellite length sizes
	
	@param data  - struture where the freq_arr and val_arr are stored
	@param cloc  - loci no. 
	@param nrSTR - nr of linked STR's loci
*/
void reorder(struct data *data,int cloc,int nrSTR);

/*
	Prints the information to a file
	
	@param data   - data struct related to the current tree
	@param ltype  - list of the DNA type per locus (m or s)
	@param path   - path to save the output into
*/
void printdata(struct data* data,char* ltype, char* path);

/*
	Gets the list of mut in lineage and returns an ordered list of mut

	@param dna    - sequence of lineology mut and its ordinal no. in the overall mut
	@param lmut	  - list of the mutations
	@param len	  - no. of mutations in sequence
	@param nsites - no. of the overall mut
	@return - an ordered mutations array
*/
int *posToOrd(double *dna,double *lmut,int len,int nsites);

int *nInd1,					//current size of sample per pop
	*listInd,				//current size of lnode[][][]
	*listIndex,				//index of indiv per sample pop
	*recLength,				//list of the rec. rates in the modern pop
	*lind1,					//original size of samples per loci
	**lind2,				//original size of samples per loci per pop.
	**listSplit,			//stores the way split events will happen (from_pop,to_pop)
	**nInd2;				//current size of sample per locus per pop
double  *coalRate,			//list of the coal. rates in the modern pop
		*migRate,			//list of the mig. rates in the modern pop
		*tEvent,			//list of the time of events in generations (last value is close to infinite)
		**den,				//list of the prob of happening an event (migration OR coalescent) for each pop
		**Nevent,	  		//pop. sizes during the splitting events in time
		**iEvent;	  		//migration rates during the spitting events in time
struct node **allnode,		//list of nodes of a genetic tree (sample indiv. + coal. indiv.)
			***lnode;		//list of nodes of a genetic tree per pop, per sample indiv.
 		
void sampLikelihood(struct params *pm, struct data *data,int printit,long int citer,char *path){
	int	cevt, cpop, cloc, cind, i, cdna,//iterators
		maxRec,							//max no. sites when considering rec events
		maxSamp,						//maximum sample size
		maxNe,							//maximum population size
		nR,								//expected number of recombination events
		npop,							//number of populations
		nevt,							//number of events
		nloc;							//number of loci
	double sSamp,						//auxiliar to calculate MaxRec
		   *mut,						//mut. rates in different loci
		   *rec;						//rec. rates in different loci
 
	mut = pm->mu;
	rec = pm->rec;
	npop = pm->npop;
	nevt = pm->nevt;
	nloc = pm->nloc - pm->nrSTR + 1;

 	/*first iteration*/
	if(citer == 0){
		listIndex = (int *)malloc(npop*sizeof(int));	
		listInd = (int *)malloc(npop*sizeof(int));
		lnode = (struct node ***)malloc(npop*sizeof(struct node **));
		den = (double **)malloc(npop*sizeof(double*));
		migRate = (double *)malloc(npop*sizeof(double));
		coalRate = (double *)malloc(npop*sizeof(double));
		recLength = (int *)malloc(npop*sizeof(int));
		nInd1 = (int *)malloc(npop*sizeof(int));
		nInd2 = (int **)malloc(nloc*sizeof(int *));
		lind1 = (int *)malloc(nloc*sizeof(int));
		lind2 = (int **)malloc(nloc*sizeof(int *));
		tEvent = (double*)malloc((nevt+1)*sizeof(double));
 		Nevent = (double**)malloc((nevt+1)*sizeof(double*));
 		iEvent = (double**)malloc((nevt+1)*sizeof(double*));
 		listSplit = (int **)malloc((nevt+1)*sizeof(int*));
		for(cpop=0;cpop<npop;++cpop){
			den[cpop] = (double *)malloc(COALEVNT*sizeof(double ));
		}
		for(cloc=0;cloc<nloc;++cloc){
			nInd2[cloc] = (int *) malloc(npop*sizeof(int));
			lind2[cloc] = (int *) malloc(npop*sizeof(int));
		}
		for(cevt=0; cevt<=nevt ; ++cevt){
			Nevent[cevt] = (double*)malloc(cpop*sizeof(double));
 			iEvent[cevt] = (double*)malloc(cpop*sizeof(double));
 			if(cevt < nevt){
				listSplit[cevt] = (int *)malloc(2*sizeof(int));
			}	
		}
	}
	
	/*fill listSplit[][] and tEvent[]*/
	if(npop>1){
		for(cevt=0, i=0 ; cevt<nevt ; ++cevt){
			listSplit[cevt][0] = pm->seq[i++];						//from pop
			listSplit[cevt][1] = pm->seq[i++];						//to pop
			if(cevt<nevt)
				tEvent[cevt] = pm->tev[cevt];						//time events
		}
		tEvent[cevt] = 1.e100;										//big no.
	}
	else{
		tEvent[0] = 1.e100;											//big no.	
	}		
	
	for(cloc=0;cloc<nloc;++cloc){
	 	/*fill iEvent[][]*/
		if(npop>1){
		 	for(cevt=0 ; cevt<=nevt ; ++cevt){
				for(cpop=0;cpop<npop;++cpop)
					if(P_mig[0][cpop].type == 4 || P_mig[0][cpop].type == 5)
						iEvent[cevt][cpop] = pm->mig[0][cpop]/2*pm->ploidy[cloc];	//mig rates of modern pop. to all the columns of iEvent	
					else
						iEvent[cevt][cpop] = pm->mig[0][cpop];
				//Choose the pop. who joins an ancient pop. during the time events
				if(cevt>0 && cevt<nevt){
					for(i=0; i < cevt ; i++){
						if(P_mig[i+1][0].type == 4 || P_mig[i+1][0].type == 5)
							iEvent[cevt][pm->seq[2*i+1]] = pm->mig[i+1][0]/2*pm->ploidy[cloc];
						else
							iEvent[cevt][pm->seq[2*i+1]] = pm->mig[i+1][0];
						iEvent[cevt][pm->seq[2*i]] = -1;
					}
				}
				if(cevt==nevt){
					for(i=0; i < cevt ; i++){
						iEvent[cevt][pm->seq[2*i+1]] = 0;
						iEvent[cevt][pm->seq[2*i]] = -1;
					}
				}
			}		
		}
		else{
			iEvent[0][0] = 0;
		}		
		/*fill listIndex*/
		for(cpop=0;cpop<npop;++cpop){
			listIndex[cpop]=cpop;
		}	
		/*fill Nevent[][]*/
	 	for(cevt=0 ; cevt<=nevt ; ++cevt){
			for(cpop=0;cpop<npop;++cpop)
				Nevent[cevt][cpop] = 2*pm->ploidy[cloc]*pm->psize[0][cpop];		//size of modern pop. to all the columns of Nevent
			
			//Choose the pop. who joins an ancient pop. during the time events
			if(cevt>0){
				for(i=0; i < cevt ; i++){
					Nevent[cevt][pm->seq[2*i+1]] = 2*pm->ploidy[cloc]*pm->psize[i+1][0];	//gives the size of the first ancient pop. to one cell of all the columns of Nevent (except the first one)	 
					Nevent[cevt][pm->seq[2*i]] = -1;								//sets the size of the pop. who coalesce with the other to -1 
				}
			}
		}		
		/*fill nInd2[][], lind2[][] and lind1[]*/
		lind1[cloc] = 0;
		for(cpop=0;cpop<npop;++cpop){
			nInd2[cloc][cpop] = pm->nsamp[cloc][cpop];
			lind2[cloc][cpop] = nInd2[cloc][cpop];
			lind1[cloc] += lind2[cloc][cpop];
		}
		/*calculate MaxRec*/
		if(pm->type[cloc]=='m'||pm->type[cloc]=='M'){
			maxRec = pm->nrSTR;
		}
		else{
			if(rec[cloc]!=0){
				maxSamp=0;
				for(cpop=0;cpop<npop;cpop++){
					if(maxSamp < nInd2[cloc][cpop])
						maxSamp = nInd2[cloc][cpop];
				}
				sSamp=0;
				for(i=1; i<maxSamp; i++){
					sSamp+=1/(double)i;
				}
				maxNe = 0;
				for(cevt=0 ; cevt<=nevt ; ++cevt){
					for(cpop=0;cpop<npop;++cpop){
						if(maxNe < Nevent[cevt][cpop])
							maxNe = Nevent[cevt][cpop];
					}
				}
				nR = (int)(rec[cloc]*maxNe*4*sSamp+0.5);
				if(nR!=0)
					maxRec = nR*100;
				else
					maxRec = 100;		
			}
			else
				maxRec = 1;
		}
		//genetic tree simulation
		simulation(lind2[cloc],lind1[cloc],mut[cloc],rec[cloc],data,cloc,pm->type[cloc],maxRec);
	}
	if(printit)
		printdata(data,pm->type,path);
	
	/*free stuff*/
	if(P_top.type == 2 || P_top.type == 5)
		free(pm->seq);		
		
} // end of sampLikelihood

void simulation(int lind[],int nind,double mut,double rec,struct data *data,int cloc,char type,int maxRec){
	int csamp,cpop,cind,cnode,i,j,cSTR,	//iterators
		currEvent,						//splitting time event being evaluated
		currPop,						//identifier of the population being evaluated
		currIndiv,						//identifier of the individual being evaluated					
		currNode,						//identifier of the node being evaluated
		currNpop,						//number of populations being evaluated
		nloc,							//number of loci
		npop,							//number of populations	
		maxnode,						//current length of allnode
		from_pop, to_pop;				//joining of from_pop to to_pop during ancester pop. listSplit
	double timeevent,					//expected time to next event
		   dtop,						//rate of all events in all the pop (1/dtop is the time for the next event (mig. OR coal.) for every pop)
		   rand;						//random no. uniformely distributed between 0 and 1
	/*only used in DNA sequence data*/
	int nmut;							//number of mutations
	double *lmut = NULL;				//list of the mutations
	
	nloc = data->nloc;
	npop = data->npop;
	currPop = npop;
	
	if(nind<=1)
		printerr("the locus have less than 2 samples, it should be cleared");
	
	/*fill listInd[] and nInd1[]*/
	for(cpop=0;cpop<npop;++cpop){
		nInd1[cpop] = lind[cpop];
		listInd[cpop] = nind;
	}
	//allocate memory to allnode
	maxnode = 2*(nind+5);
	allnode = (struct node **)malloc(maxnode*sizeof(struct node *));	
	//allocate memory to lnode[]
	for(cpop=0,i = 0 ; cpop<npop ; ++cpop){
		lnode[cpop] = (struct node **)malloc(listInd[cpop]*sizeof(struct node *));
		for(cind=0;cind<nInd1[cpop];++cind){ 
			lnode[cpop][cind] = (struct node *)malloc(sizeof(struct node));
			lnode[cpop][cind]->d[0] = lnode[cpop][cind]->d[1]= NULL;
			lnode[cpop][cind]->a[0] = lnode[cpop][cind]->a[1] = NULL;
			lnode[cpop][cind]->time = 0.0;
			lnode[cpop][cind]->spop = listIndex[cpop];
			lnode[cpop][cind]->ospop = listIndex[cpop];
			lnode[cpop][cind]->cut = -1;
			lnode[cpop][cind]->cbeg = 0;
			lnode[cpop][cind]->cend = maxRec-1;
			lnode[cpop][cind]->desc = (int *)malloc(maxRec*sizeof(int));
			for(j=0;j<maxRec;j++){
				lnode[cpop][cind]->desc[j] = 1;
			}
			if(type=='M'||type=='m'){
				lnode[cpop][cind]->dnaM=NULL;
			}
			else{
				lnode[cpop][cind]->dnaS = NULL;
				lnode[cpop][cind]->len = 0;
				lnode[cpop][cind]->lenMax = 0;
			}
			allnode[i] = lnode[cpop][cind];
			++i;
		}
		recLength[cpop] = (maxRec-1)*nInd1[cpop];
	}
	/*rearrange populations in case of 0 population samples */
	for(i=npop-1 ; i>=0 ; i--){	
		if(nInd1[i] == 0){
			for(cpop=i;cpop<currPop-1;++cpop){
				listIndex[cpop] = listIndex[cpop+1];
				nInd1[cpop] = nInd1[cpop+1];  
				
				//increase memory to lnode[]
				if(nInd1[cpop] >= listInd[cpop]-5){
					listInd[cpop] *= 2;
					lnode[cpop] = (struct node **)realloc(lnode[cpop],listInd[cpop]*sizeof(struct node *));
				}
				
				for(cind=0;cind<nInd1[cpop];++cind)
					lnode[cpop][cind] = lnode[cpop+1][cind];
			}
			--currPop;					//the origin pop. is erased
		}
	}
	/*set the coalescente and migration rate*/
	currEvent = 0;
	for(cpop=0;cpop<npop;++cpop){
		coalRate[cpop] = 1.0/Nevent[currEvent][cpop];
		migRate[cpop] = iEvent[currEvent][cpop];
	}
	
	/*Build the genetic tree (stops when there is only one node left)*/
	currIndiv = nind;
	currNpop = npop;
	currNode = nind;
	timeevent = 0.0;
	while(1){
		//increase memory to lnode
		for(cpop=0;cpop<currPop;++cpop){ 
			if(nInd1[cpop] >= listInd[cpop]-5){
				listInd[cpop] *= 2;
				lnode[cpop] = (struct node **)realloc(lnode[cpop],listInd[cpop]*sizeof(struct node *));
			}
		}
		//increase memory to allnode
		if(currNode >= maxnode-5){ 
			maxnode=2*(maxnode+5);
			allnode = (struct node **)realloc(allnode,maxnode*sizeof(struct node *));
		}
		if(maxRec == 1)
			den[0][0] = 0.0;
		else	
			den[0][0] = ((double)recLength[0]/((double)maxRec-1.0))*rec;			//prob of a recombination event in pop 0
		den[0][1] = den[0][0] + coalRate[listIndex[0]]*nInd1[0]*(nInd1[0]-1.0)*0.5;	//prob of a coalescent event in pop 0
		den[0][2] = den[0][1] + migRate[listIndex[0]]*nInd1[0];						//prob of a migration event in pop 0
		for(cpop=1;cpop<currPop;++cpop){
			if(maxRec == 1)
				den[cpop][0] = den[cpop-1][2] + 0.0;
			else
				den[cpop][0] = den[cpop-1][2] + ((double)recLength[cpop]/((double)maxRec-1.0))*rec;										//prob of a recombination event in pop cpop
			den[cpop][1] = den[cpop][0] + coalRate[listIndex[cpop]]*nInd1[cpop]*(nInd1[cpop]-1.0)*0.5;//prob of a coalescent event in pop cpop
			den[cpop][2] = den[cpop][1] + migRate[listIndex[cpop]]*nInd1[cpop];							//prob of a migration event in pop cpop
		}
		dtop = den[currPop-1][2];	//rate of all events combined
		timeevent +=  expdev()/dtop; //expected time to next event	
		
		/*splitting event reached, join descendent pops*/
		if(timeevent >= tEvent[currEvent]){
			from_pop = listSplit[currEvent][0];	//assigned the pop. that will join the other
			to_pop = listSplit[currEvent][1];	//assigned the pop. that will be joined

			colonize(from_pop,to_pop,&currPop);
			timeevent = tEvent[currEvent];
			++currEvent;						//choose the next event
			--currNpop;							//the no. of pop. are decreased

			//update Ne values
			for(cpop=0;cpop<npop;++cpop){
				if(Nevent[currEvent][cpop] > 0){
					coalRate[cpop] = 1.0/Nevent[currEvent][cpop];
					if(currNpop == 1)
						migRate[cpop] = 0;
					else
						migRate[cpop] = iEvent[currEvent][cpop];
				}
				//for non-existent pops
				else{
					coalRate[cpop] = 0;
					migRate[cpop] = 0;
				}
			}
		}
		/*choose next event(coalescence, migration or recombination)*/
		else{
			rand = gfsr8();

			//choose between all the pop. and type of events which one is going to happen
			//i = 0(recombination event), 1(coalescent event), 2(migration event)
			for(cpop=0;cpop < currPop;++cpop){
				for(i=0 ; i<COALEVNT ; ++i){
					if(rand < den[cpop][i]/dtop)
						break;
				}
				if(i<COALEVNT)
					break;
			}
			//trap possibility that den[cpop][i]/dtop not quite 1
			if(cpop == currPop){
				cpop = currPop-1;
				i = 2;
			}
			/*recombination event*/
			if(i==0){
				recombination(cpop,type,lind[cpop],timeevent,&currNode,&currIndiv,maxRec);
			}
			/*coalescent event*/
			else if(i==1){
				coalescence(cpop,type,lind[cpop],timeevent,&currNode,&currIndiv,maxRec);
				
				//STOP condition: while(1)
				if(currIndiv == 1)
					break;
			}
			/*migration event*/
			else{
				migration(cpop,npop,&currPop,currEvent);
			}
		}
	} // end of while(1)

	/*initialize MRCA dna data*/
	if(type=='M'||type=='m'){
		allnode[currNode-1]->dnaM = (int *)malloc(maxRec*sizeof(int));
		if(maxRec!=1){
			for(cSTR=0; cSTR<maxRec; ++cSTR){
				if(!(allnode[currNode-1]->d[0]->desc[cSTR] == nind || allnode[currNode-1]->d[1]->desc[cSTR] == nind))
					allnode[currNode-1]->dnaM[cSTR] = DnaSizeM;	//ancestral type - set to what we want
			}
		}
		else{
			allnode[currNode-1]->dnaM[0] = DnaSizeM;
		}
	}
	else{
		allnode[currNode-1]->len = 0;			//ancestral type - length of dna sequence starts in 0
		allnode[currNode-1]->lenMax = DnaSizeS;	//ancestral type - set to what we want
		allnode[currNode-1]->dnaS = (double *)malloc(allnode[currNode-1]->lenMax*sizeof(double));				   			  

		nmut = 0;
	}
	
	/*fill all the nodes of the tree*/
	for(cnode=currNode-1;cnode>=0;--cnode){
		//not dealing with the MRCA
		if(allnode[cnode]->a[0]!=NULL || allnode[cnode]->a[1]!=NULL)		
			treefill(allnode[cnode],data,cloc,mut,lind[cpop],type,&nmut,maxRec,&lmut);
	}

	/*calculate freq() and val() from the samples nodes*/
	if(type=='M'||type=='m'){
		for(cSTR=0; cSTR<maxRec; cSTR++)
			data->ldna[cloc+cSTR] = 0;
	}
	else{
		data->ldna[cloc] = 0;
		data->lsites[cloc] = nmut;
		for(csamp=0; csamp<data->tsamp[cloc]; csamp++)
			data->valS[cloc][csamp] = (char*)malloc((nmut+1)*sizeof(char));
	}
	for(cind=0;cind<nind;++cind){
		if(type=='M'||type=='m')
			makeSampM(allnode[cind],data,cloc,maxRec);
		else
			makeSampS(allnode[cind],data,cloc,nmut,lmut,maxRec);
	}
	if(type=='M'||type=='m')
		reorder(data,cloc,maxRec);
		
	/*free memory*/
	for(cnode=0;cnode<currNode;++cnode){
		free(allnode[cnode]->desc);
		if(type=='M'||type=='m'){	
			free(allnode[cnode]->dnaM);
		}
		else{
			free(allnode[cnode]->dnaS);
		}
		free(allnode[cnode]);
	}
	if((type=='S'||type=='s')&&lmut!=NULL){	
		free(lmut);
	}	
	free(allnode);
	for(cpop=0;cpop<npop;++cpop)
		free(lnode[cpop]);
	
}	//end of simulation

void colonize(int from_pop,int to_pop,int *currPop){
	int cind,cpop,		//iterator
		ifrom,			//index of the pop. that is going to join the other
		ito,			//index of the pop. that is going to be joined
		nindiv;			//(size-1) of the pop. that is going to join the other
	struct node *temp;	//temporary pntr to the individual that is going to join the other pop.

	//get value of ifrom
	for(ifrom=0;ifrom<(*currPop);++ifrom){
		if(listIndex[ifrom] == from_pop)
			break;
	}
	//The from_pop doesn't have any individuals (it is not present in listIndex[])
	if(ifrom==(*currPop))
		return;

	//get value of ito	
	for(ito=0;ito<(*currPop);++ito){
		if(listIndex[ito] == to_pop)
			break;
	}

	//start of stuff
	nindiv = nInd1[ifrom]-1;
	while(1){
		//STOP condition - while(1) if the size of origin reaches 0 end the cycle
		if(nindiv < 0)
			break;
		
		lnode[ifrom][nindiv]->spop = to_pop;						//change the type of pop. the individual belongs to
		temp = lnode[ifrom][nindiv];								//point the temp pntr to the individal
		lnode[ifrom][nindiv] = lnode[ifrom][nInd1[ifrom]-1];		//points the old pointer to the last individual of the pop.

		//no individuals in the destination pop
		if(ito == (*currPop)){
			//only one individual in the origin
			if(nInd1[ifrom] == 1){
				ito = ifrom;
				listIndex[ito] = to_pop;	//the origin becames the joined pop.
				nInd1[ito] = 0;				//no individuals at the destination pop momentanly
			}
			//more than one individual in origin
			else{
				++(*currPop);				//the destination pop. is added
				listIndex[ito] = to_pop;	//the destiation pop. is created	
				nInd1[ito] = 0;				//no individuals in the destination momentanly
				--nInd1[ifrom];				//the origin decreases its size by one
				recLength[ifrom]-= (temp->cend - temp->cbeg);	//remove number of recombination places in the lineage removed from donor deme
				recLength[ito] = (temp->cend - temp->cbeg); 	//add these to those in the recipient deme
			}
		}
		//at least one individual in the destination pop
		else{
			//only one individual in the origin
			if(nInd1[ifrom] == 1){
				nInd1[ifrom]--;
				recLength[ifrom] -= (temp->cend - temp->cbeg);
				recLength[ito] += (temp->cend - temp->cbeg);
				for(cpop=ifrom;cpop<(*currPop)-1;++cpop){
					listIndex[cpop] = listIndex[cpop+1];
					nInd1[cpop] = nInd1[cpop+1];
					recLength[cpop] = recLength[cpop+1];
					
					//increase memory to lnode[]
					if(nInd1[cpop] >= listInd[cpop]-5){
						listInd[cpop] *= 2;
						lnode[cpop] = (struct node **)realloc(lnode[cpop],listInd[cpop]*sizeof(struct node *));
					}
					
					for(cind=0;cind<nInd1[cpop];++cind)
						lnode[cpop][cind] = lnode[cpop+1][cind];
				}
				nInd1[cpop]=0;
				recLength[cpop]=0;
				--(*currPop);					//the origin pop. is erased
				if(ito>ifrom)
					--ito;					//update the index of the created destination pop.
			}
			//more than an individual in origin
			else{
				--nInd1[ifrom];			//decrease the size of the origin by one individual
				recLength[ifrom] -= (temp->cend - temp->cbeg);
				recLength[ito] += (temp->cend - temp->cbeg);
			}
		}
		lnode[ito][nInd1[ito]] = temp;	//point the last individual of the destination
		++nInd1[ito];					//increase the size of the destination
		
		//increase memory to lnode[]
		if(nInd1[ito] >= listInd[ito]){
			listInd[ito] *= 2;
			lnode[ito] = (struct node **)realloc(lnode[ito],listInd[ito]*sizeof(struct node *));
		}
		--nindiv;								//decrease the total size of origin
	} //end of while(1)

} // end of colonize

void recombination(int cpop,char type,int nind,double time,int *currNode,int *currIndiv,int maxRec){
	int cind,crec,i,	//iterators
		ind1,			//node that suffers recombination
		cut,			//place where the recombination happens
		spc,			//space from recLength travelled so far
		ind1_spc;		//space from recLength that suffers rec
	struct node *p1,	//ancester to the left of the breakpoint
				*p2;	//ancester to the right of the breakpoint

	//choose an individual to apply recombination on it
	ind1_spc = disrand(1,recLength[cpop]);
	spc=0;
	for(cind=0;cind<nInd1[cpop];++cind){
		spc+= (lnode[cpop][cind]->cend - lnode[cpop][cind]->cbeg);
		if(spc>=ind1_spc){
			ind1=cind;
			break;
		}
	}
	//choose the place where the breakpoint ocurres
	cut = disrand(lnode[cpop][ind1]->cbeg,lnode[cpop][ind1]->cend-1);
	
	//create ancester node to the left of break point
	p1 = (struct node *)malloc(sizeof(struct node));
	p1->d[0] = lnode[cpop][ind1];						//descendent
	p1->d[1] = p1->a[0] = p1->a[1] = NULL;				//ancesters
	p1->time = time;									//time of event
	p1->spop = listIndex[cpop];     			        //population identifier
	p1->ospop = listIndex[cpop];						//population identifier
	if(type=='m'||type=='M'){
		p1->dnaM = NULL;								//iniciate dnaM
	}
	else{
		p1->dnaS = NULL;                    			//DNA sequence
		p1->len = 0;									//length of DNA sequence
		p1->lenMax = 0;									//maximum length of DNA sequence
	}
	p1->desc = (int *)malloc(maxRec*sizeof(int));		//list to keep track of number of descendents
	p1->cbeg = p1->d[0]->cbeg;              			//first site with info
	p1->cend = cut;										//last site with info
	for(i=0;i<=cut;++i)
		p1->desc[i] = p1->d[0]->desc[i];
	for(i=cut+1;i<maxRec;++i)
		p1->desc[i] = 0;
	//readjust cend
	for(i=p1->cend;i >= 0;--i){
		if(p1->desc[i] > 0 && p1->desc[i] < nind){
			p1->cend = i;
			break;
		}
	}

	//create ancester node to the right of break point
	p2 = (struct node *)malloc(sizeof(struct node));
	p2->d[0] = lnode[cpop][ind1];						//descendent
	p2->d[1] = p2->a[0]= p2->a[1] = NULL;				//ancesters
	p2->time = time;									//time of event
	p2->spop = listIndex[cpop];							//population identifier
	p2->ospop = listIndex[cpop];						//population identifier
	if(type=='m'||type=='M'){
		p2->dnaM = NULL;								//iniciate dnaM
	}
	else{
		p2->dnaS = NULL;                    		    //DNA sequence
		p2->len = 0;									//length of DNA sequence
		p2->lenMax = 0;									//maximum length of DNA sequence
	}
	p2->desc = (int *)malloc(maxRec*sizeof(int));		//list to keep track of number of descendents
  	p2->cbeg = cut+1;										//first site with info
	p2->cend = p2->d[0]->cend;							//last site with info
	for(i=0;i<=cut;++i)
		p2->desc[i] = 0; 
	for(i=cut+1;i<maxRec;++i)
		p2->desc[i] = p2->d[0]->desc[i];
	//readjust cbeg
	for(i=p2->cbeg;i <maxRec;++i){
		if(p2->desc[i] > 0 && p2->desc[i] < nind){
			p2->cbeg = i;
			break;
		}
	}
	recLength[cpop] -= (lnode[cpop][ind1]->cend - lnode[cpop][ind1]->cbeg);
	recLength[cpop] += (p1->cend - p1->cbeg) + (p2->cend - p2->cbeg);

	//update the list of existing nodes (lnode[cpop][cindiv])
	lnode[cpop][ind1]->a[0] = p1;
	lnode[cpop][ind1]->a[1] = p2;
	lnode[cpop][ind1]->cut = cut;
	lnode[cpop][ind1] = p1;
	lnode[cpop][nInd1[cpop]] = p2;
	++nInd1[cpop];
	++(*currIndiv);

	//update the list of all nodes (allnode[cnode])
	allnode[(*currNode)++] = p1;
	allnode[(*currNode)++] = p2;

} //end of recombination

void coalescence(int cpop,char type,int nind,double timeevent,int *currNode,int *currIndiv,int maxRec){
	int crec,i,			//iterator
		ind1,			//individual that coalesces
		ind2,			//individual that coalesces
		jstart,			//auxiliar
		jend,			//auxiliar
		temp;			//auxiliar
	struct node *p1;	//node product of a coalescent event

	while(1){
		ind1 = disrand(0,nInd1[cpop]-1);
		ind2 = disrand(0,nInd1[cpop]-1);
		if(ind2 != ind1)
			break;
	}
	if(ind1 > ind2){
		temp = ind1;
		ind1 = ind2;
		ind2 = temp;
	}
	p1 = (struct node *) malloc(sizeof(struct node));
	p1->time = timeevent;						//assigning the time event to the ancester of the coalescence
	p1->d[0] = lnode[cpop][ind1];				//assigning the first descendent of the ancester of the coalescence
	p1->d[1] = lnode[cpop][ind2];				//assigning the second descendent of the ancester of the coalescence
	p1 -> a[0] = p1 -> a[1] =  NULL;				//no ancester yet
	p1->spop = listIndex[cpop];					//assingning the belonging pop. to the ancester of the coalesence
	p1->desc = (int *)malloc(maxRec*sizeof(int));
	for(i=0;i<maxRec;++i)
		p1->desc[i] = p1->d[0]->desc[i]+p1->d[1]->desc[i];
	p1->cbeg = p1->cend = 0;
	if(type=='M'||type=='m')
		p1->dnaM = NULL;						//iniciate dnaM
	else{
		p1->lenMax = 0;							//maximum length of DNA sequence
		p1->len = 0;							//iniciate len
		p1->dnaS = NULL;						//iniciate dnaS
	}
 	//choose the leftmost cbeg (guaranteed that desc is 0 or nind to the left of this
	if(p1->d[0]->cbeg < p1->d[1]->cbeg)
		jstart = p1->d[0]->cbeg;	
	else
		jstart = p1->d[1]->cbeg;
	for(i=jstart;i<maxRec;++i){
		if(p1->desc[i] > 0 && p1->desc[i] < nind){
			p1->cbeg = i;
			break;
		}
	} 
	//choose the rightmost cend (guaranteed that desc is 0 or nind to the right of this)
	if(p1->d[0]->cend > p1->d[1]->cend)
		jend = p1->d[0]->cend;
	else
		jend = p1->d[1]->cend;
	for(i=jend; i>= 0; --i){
		if(p1->desc[i] > 0 && p1->desc[i] < nind){
			p1->cend = i;
			break;
		}
	}
	//update the overall length of the population that can have a recombination event
	recLength[cpop] -= (p1->d[0]->cend - p1->d[0]->cbeg);
	recLength[cpop] -= (p1->d[1]->cend - p1->d[1]->cbeg);
	recLength[cpop] += (p1->cend - p1->cbeg);
	
	lnode[cpop][ind1]->a[0] = p1;	//assigning the ancester of the first descendent of the coalescence
	lnode[cpop][ind2]->a[0] = p1;	//assigning the ancester of the first descendent of the coalescence
	lnode[cpop][ind1]->a[1] = lnode[cpop][ind2]->a[1] = NULL;
	lnode[cpop][ind1] = p1; 		//the first descendent is replaced with his ancester of the coalescence
	lnode[cpop][ind2] = lnode[cpop][nInd1[cpop]-1];	//the pointer is now pointing to the last individual of the pop.
	--nInd1[cpop];					//no. of individuals in a sample that are still electable to be chosen to migrate or coalesce
	--(*currIndiv);					//total no. of individuals that are still electable to be chosen to migrate or coalesce
	allnode[(*currNode)] = p1; 		//list of the individuals that are the product of a coalescent event
	++(*currNode);					//iterator for the allnode

} // end of coalescence

void migration(int ifrom,int npop,int *currPop,int currEvent){
	int cpop, cind,		//iterators
		from_ind,		//migrating individual from the origin pop. (as seen back in time)
		to_pop,			//pop. that is the destination of migration (as seen back in time)
		getito,			//boolean to keep track if ito was chosen
		ito;			//index of the listIndex for the pop. that is the destination (as seen back in time)
	double sumweights,	//check for error on setting the migration weights matrix
		   rand;		//random number for the migration weights
	struct node *temp;	//individual that is going to join the other pop.

	from_ind = disrand(0,nInd1[ifrom]-1);

	//choose the destiny pop. (as seen back in time)
	if(npop<=2 || M_migw.type==0){
		while(1){
			to_pop = disrand(0,npop-1);	
			if(to_pop == listIndex[ifrom])
				continue;
			if(Nevent[currEvent][to_pop] <= 0)
				continue;
			break;
		}
	}
	else{
		//check for errors in the migration matrix
		if(M_migw.m[listIndex[ifrom]][currEvent][listIndex[ifrom]]!=0)
			printerr("migration weights for the considered population have to be 0");
		sumweights=0.0;
		for(cpop=0; cpop<npop; cpop++)
			sumweights+=M_migw.m[listIndex[ifrom]][currEvent][cpop];
		if(sumweights!=1.0)
			printerr("the sum of the prob for a given migration as to be 1.0");
		
		//choose the destiny pop. (as seen back in time)
		getito = 1;
		while(getito){
			rand = gfsr8();
			sumweights=0.0;
			for(cpop=0; cpop<npop; cpop++){
				sumweights+=M_migw.m[listIndex[ifrom]][currEvent][cpop];
				if(sumweights>=rand){
					to_pop = cpop;
					if(Nevent[currEvent][to_pop]>0)
						getito=0;
					break;
				}	
			}
		}
	}
	//find the index of the destiny pop.
	for(ito=0;ito<(*currPop);++ito){
		if(listIndex[ito] == to_pop)
			break;
	}

	lnode[ifrom][from_ind]->spop = to_pop;	//new belonging pop. assigned to the migrating individual
	temp = lnode[ifrom][from_ind];			//pointer temp grabs the migrate individual
	lnode[ifrom][from_ind] = lnode[ifrom][nInd1[ifrom]-1];	//the pntr will point to the last individual of the pop.

	//destin pop. doesn't have any individuals
	if(ito == (*currPop)){
		//only one individual in the origin pop.
		if(nInd1[ifrom] == 1){
			ito = ifrom;					
			listIndex[ito] = to_pop;	//the origin pop. disapears and the destiny pop. appears 
			nInd1[ito] = 0;				//the size of the origin pop. is zero
		}
		//more than one individual in the origin pop.
		else{
			++(*currPop);									//a new pop. appears
			listIndex[ito] = to_pop;					//the new pop. is registered in the listIndex	
			nInd1[ito] = 0;
			--nInd1[ifrom];								//the origin pop. decrease the no. of individuals
			recLength[ifrom]-= (temp->cend - temp->cbeg);	//remove number of recombination places in the lineage removed from donor deme
			recLength[ito] = (temp->cend - temp->cbeg); 	//add these to those in the recipient deme
		}
	}
	//destination pop. still have some individuals
	else{
		//only one individual in the origin specie
		if(nInd1[ifrom] == 1){
			nInd1[ifrom]--;
			recLength[ifrom] -= (temp->cend - temp->cbeg);
			recLength[ito] += (temp->cend - temp->cbeg);
			for(cpop=ifrom;cpop<(*currPop)-1;++cpop){
				listIndex[cpop] = listIndex[cpop+1];
				nInd1[cpop] = nInd1[cpop+1];
				recLength[cpop] = recLength[cpop+1];
				
				//increase memory to lnode[]
				if(nInd1[cpop] >= listInd[cpop]-5){
					listInd[cpop] *= 2;			
					lnode[cpop] = (struct node **)realloc(lnode[cpop],listInd[cpop]*sizeof(struct node *));
				}
	
				for(cind=0;cind<nInd1[cpop];++cind)
					lnode[cpop][cind] = lnode[cpop+1][cind];
			}
			nInd1[cpop]=0;
			recLength[cpop]=0;
			--(*currPop);			//the origin pop. disapears		
			if(ito>ifrom)
				--ito;			//the index of the destiny pop. is reallocated
		}
		//more than one individual in the origin pop.
		else{
			--nInd1[ifrom];			//one less individual in the origin pop.
			recLength[ifrom] -= (temp->cend - temp->cbeg);
			recLength[ito] += (temp->cend - temp->cbeg);
		}
	}
	lnode[ito][nInd1[ito]] = temp;	//the destiny pop. gets the migrate individual
	++nInd1[ito];					//the size of the destiny pop. increases by one
	
	//increase memory to lnode[]
	if(nInd1[ito] >= listInd[ito]){
		listInd[ito] *= 2;
		lnode[ito] = (struct node **)realloc(lnode[ito],listInd[ito]*sizeof(struct node *));
	}
	
} // end of migration

void treefill(struct node *p,struct data *data, int cloc,double mut,int nind,char type,int *nmut,int maxRec,double **lmut){
	int cmut,clen,cSTR,	//iterator
		pos,			//current mutation in a particular sample
		mutInt, 		//position of the mutation within a 0-1 segment
		len,			//number of mutations of a sample
		mutno;			//random mut. no. of mut in the node
	double time,		//time between nodes
		   mutDouble;	//position of the mutation within a 0-1 segment
	
	/*recombination - join fragments*/
	if(p->a[0] != NULL && p->a[1] != NULL){
		if(type=='m'||type=='M'){
			p->dnaM = (int *)malloc(maxRec*sizeof(int));
			for(cSTR=0; cSTR<=p->cut; ++cSTR){
				if(p->d[1] == NULL)
					p->dnaM[cSTR] = p->a[0]->dnaM[cSTR];
				else{
					if(p->desc[cSTR] == nind && !(p->d[0]->desc[cSTR] == nind || p->d[1]->desc[cSTR] == nind))
						p->dnaM[cSTR] = DnaSizeM;
					else
						p->dnaM[cSTR] = p->a[0]->dnaM[cSTR];
				}
			}
			for(cSTR=p->cut+1; cSTR<maxRec; ++cSTR){
				if(p->d[1] == NULL)
					p->dnaM[cSTR] = p->a[1]->dnaM[cSTR];
				else{
					if(p->desc[cSTR] == nind && !(p->d[0]->desc[cSTR] == nind || p->d[1]->desc[cSTR] == nind))
						p->dnaM[cSTR] = DnaSizeM;
					else
						p->dnaM[cSTR] = p->a[1]->dnaM[cSTR];
				}
			}	
		}
		else{
			p->dnaS = (double*)malloc((p->a[0]->lenMax + p->a[1]->lenMax)*sizeof(double));
		
			len = 0;
			//get all the mutations from ind1
			for(clen=0; clen<p->a[0]->len ;clen++){
				p->dnaS[len] = p->a[0]->dnaS[clen];
				len++;
			}
			//get all the mutations from ind2 that are diferent from ind1
			for(clen=0; clen<p->a[1]->len; clen++){
				p->dnaS[len]=p->a[1]->dnaS[clen];
				len++;
			}
			p->len = len;
			p->lenMax = len;
		}
	}
	/*coalescence - point node to its ancester*/
	else{
		//choose analysis type
		if(type=='M'||type=='m'){
			p->dnaM = (int *)malloc(maxRec*sizeof(int));
			if(maxRec!=1){
				for(cSTR=0;cSTR<maxRec;++cSTR){
					if(p->d[1] == NULL)
						p->dnaM[cSTR] = p->a[0]->dnaM[cSTR];
					else{
						if(p->desc[cSTR] == nind && !(p->d[0]->desc[cSTR] == nind || p->d[1]->desc[cSTR] == nind))
							p->dnaM[cSTR] = DnaSizeM;
						else
							p->dnaM[cSTR] = p->a[0]->dnaM[cSTR];
					}
				}
			}
			else{
				p->dnaM[0] = p->a[0]->dnaM[0];
			}
		}
		else if(type=='S'||type=='s'){
			if(maxRec!=1){
				p->lenMax = p->a[0]->lenMax;
				p->dnaS = (double *)malloc(p->lenMax*sizeof(double));
				pos=0;
				for(clen=0 ; clen<p->a[0]->len ; ++clen){
					mutInt = (int)((p->a[0]->dnaS[clen]*((double)maxRec-1.0))+0.5);
					if(p->desc[mutInt] != 0 && p->desc[mutInt] != nind){
						p->dnaS[pos] = p->a[0]->dnaS[clen];	//adds the new mutation position
						pos++;
					}
				}
				p->len = pos;
			}
			else{
				p->lenMax = p->a[0]->lenMax;
				p->len = p->a[0]->len;
				p->dnaS = (double *)malloc(p->lenMax*sizeof(double));
				for(clen=0; clen<p->a[0]->len; clen++)
					p->dnaS[clen]=p->a[0]->dnaS[clen];
			}
		}
	}
	/*add mutation events*/
	time = p->a[0]->time - p->time;
	if(type=='M'||type=='m'){
		if(maxRec!=1){
			for(cSTR=0;cSTR<maxRec;++cSTR){
				if(p->desc[cSTR] == 0)
					continue;
				if(p->desc[cSTR] == nind)
					continue;
				mutno = poidev(time*mut);
				for(cmut=0;cmut<mutno;++cmut)
					p->dnaM[cSTR] = p->dnaM[cSTR] + disrand(0,1)*2-1;	//adds or not a mut (change length)	
			}
		}
		else{
			mutno = poidev(time*mut);
			for(cmut=0;cmut<mutno;++cmut)
				p->dnaM[0] = p->dnaM[0] + disrand(0,1)*2-1;	//adds or not a mut (change length)	
		}
	}
	else if(type=='S'||type=='s'){
		mutno = poidev(time*mut);
		if(maxRec!=1){
			for(cmut=0;cmut<mutno;++cmut){
				if(p->len == p->lenMax){
					p->dnaS = (double *)realloc(p->dnaS,(p->lenMax+DnaSizeS)*sizeof(double));
					p->lenMax += DnaSizeS;
				}
				(*lmut) = (double*)myAlloc((*lmut),((*nmut)+1)*sizeof(double));
				
				mutDouble = gfsr8();
				mutInt = (int)((mutDouble*((double)maxRec-1.0))+0.5);
				if(p->desc[mutInt] != 0){
					(*lmut)[(*nmut)] = mutDouble;	//adds the new mutation to the list of mutations
					p->dnaS[p->len] = mutDouble;	//adds a mut. to the current site of dna sequency
					p->len++;					//branch mut. no. is increased by one 
					(*nmut)++;					//overall mut. no. is increased by one
				}
			}
		}
		else{
			for(cmut=0;cmut<mutno;++cmut){
				if(p->len == p->lenMax){
					p->dnaS = (double *)realloc(p->dnaS,(p->lenMax+DnaSizeS)*sizeof(double));
					p->lenMax += DnaSizeS;
				}
				p->dnaS[p->len] = (double)(*nmut);	//adds a mut. to the current site of dna sequency
				p->len++;					//branch mut. no. is increased by one 
				(*nmut)++;					//overall mut. no. is increased by one		
			}
		}
	}
	
} //end of treefill

void makeSampM(struct node *p,struct data *data, int cloc, int nrSTR){
	int call,cSTR,		//iterator
		npop,			//number of pop
		spop,			//origin pop
		ispop;			//index for a particular pop
	
	npop = data->npop;
	
	//this individual belongs to the sample (last generation)
	if(p->d[0] == NULL && p->d[1] == NULL){ 
		for(cSTR=0; cSTR<nrSTR; cSTR++){
			//gets the place where the current allele is stored
			for(call=0 ; call < data->ldna[cloc+cSTR] ; ++call){
				if(data->valM[cloc+cSTR][call] == p->dnaM[cSTR])
					break;
			}
			spop = p->ospop;								
			//current allele is present, increase its freq 
			if(call < data->ldna[cloc+cSTR]){
				++data->freq[cloc+cSTR][spop][call]; 	//add one more individual
			}	
			//current dna isn't listed yet 
			else{
				//adds the allele to the list of alleles		
				data->valM[cloc+cSTR] = realloc(data->valM[cloc+cSTR],(data->Nmax[cloc+cSTR]+1)*sizeof(int));
				data->valM[cloc+cSTR][call] = p->dnaM[cSTR]; 		//get dna size
				//declares the freq of the current dna in all the pop.
				for(ispop=0;ispop<npop;++ispop){
					data->freq[cloc+cSTR][ispop] = realloc(data->freq[cloc+cSTR][ispop],(data->Nmax[cloc+cSTR]+1)*sizeof(int));
					data->freq[cloc+cSTR][ispop][call] = 0;	//iniciate freq of the dna size to all pop.
				}
				data->freq[cloc+cSTR][spop][call] = 1;		//increase the freq of the current allele in the current pop.
				++data->ldna[cloc+cSTR];						//increase the no. of all alleles
				
				if(data->ldna[cloc+cSTR] > data->Nmax[cloc+cSTR])
					data->Nmax[cloc+cSTR] ++;
			}
		}
	}

} //end of makeSampM

void makeSampS(struct node *p,struct data *data, int cloc,int nsites,double *lmut,int maxRec){
	int i, call, cpop, cmut, csites,	//iterator
		npop,							//number of populations
		spop,							//belonging pop. of the current node
		*orderdna;		//sequence of mutations ordered
	char *dna;							//dna sequence ([0,0,1,0,1]) 
	
	npop = data->npop;
	spop = p->ospop;

	dna = (char *)malloc((nsites+1)*sizeof(char));
	if(maxRec!=1){
		//put the mutation list in order
		shell_sort_double(lmut,nsites);
		//get the ordered positions of the dna samples 
		orderdna = posToOrd(p->dnaS,lmut,p->len,nsites);
		shell_sort_int(orderdna,p->len);
	
		cmut = 0;
		for(csites=0;csites<nsites;++csites){
			//no more mut in the linealogy OR linealogy mut. is bigger than the overall mut.
			if(cmut >= p->len || orderdna[cmut] > csites)
				dna[csites] = '0';	//add a no mut. mark
			//linealogy mut. is the same order as the overall mut.
			else{
				dna[csites] = '1';	//add a mut. mark
				++cmut;				//next linealogy mut.	
			}
		}
		dna[csites] = '\0';
	}
	else{
		cmut = 0;
		for(csites=0; csites<nsites; ++csites){
			//no more mut in the linealogy OR linealogy mut. is bigger than the overall mut.
			if(cmut >= p->len || p->dnaS[cmut] > csites)
				dna[csites] = '0';	//add a no mut. mark
			//linealogy mut. is the same order as the overall mut.
			else{
				dna[csites] = '1';	//add a mut. mark
				++cmut;				//next linealogy mut.	
			}
		}
		dna[csites] = '\0';	
	}
			
	//gets the place where the current dna is stored
	for(call=0; call<data->ldna[cloc]; ++call){
		if(strcmp(data->valS[cloc][call],dna)==0){
			break;
		}
	}
	//current haplotype is present, increase its freq 
	if(call<data->ldna[cloc]){
		++data->freq[cloc][spop][call];
	}	
	//current dna isn't listed yet 
	else{
		//adds the haplotype to the list of haplotypes		
		strcpy(data->valS[cloc][call],dna);	
		//declares the freq of the current dna in all the pop.
		for(cpop=0;cpop<npop;++cpop){
			data->freq[cloc][cpop][call] = 0;
		}
		data->freq[cloc][spop][call] = 1;	//increase the freq of the current dna in the current pop.
		++data->ldna[cloc];					//increase the no. of all haplotypes
		if(data->ldna[cloc] > data->Nmax[cloc])		
			data->Nmax[cloc]++;
	}
	
	/*free stuff*/
	free(dna);
	if(maxRec!=1)
		free(orderdna);

}	//end of makeSampS

void reorder(struct data *data,int cloc,int nrSTR){
	int call,cpop,cSTR,					//iterators
		npop,						//number of populations
		*sorted,					//list of the frequency sorted by microsatellite sizes
		*indx;						//indexes of the sorted list of microsatellite sizes
	
	npop = data->npop;
	
	//get the indexes of the sorted list of microssatelite sizes
	for(cSTR=0; cSTR<nrSTR; cSTR++){
		sorted = (int *)malloc(data->ldna[cloc+cSTR]*sizeof(int));
		indx = (int *)malloc(data->ldna[cloc+cSTR]*sizeof(int));
		isorti('a',data->ldna[cloc+cSTR],data->valM[cloc+cSTR],indx);	
		
		for(cpop=0;cpop<npop;++cpop){
			for(call=0 ; call < data->ldna[cloc+cSTR] ; ++call){
				sorted[call] = data->freq[cloc+cSTR][cpop][indx[call]];
			}
			for(call=0 ; call < data->ldna[cloc+cSTR] ; ++call){
				data->freq[cloc+cSTR][cpop][call] = sorted[call];
			}
		}
		
		for(call=0;call<data->ldna[cloc+cSTR];++call){
			sorted[call] = data->valM[cloc+cSTR][indx[call]];
		}
		for(call=0;call<data->ldna[cloc+cSTR];++call){
			data->valM[cloc+cSTR][call] = sorted[call];
		}
		//free stuff
		free(sorted);
		free(indx);
	}
		
	
} // end of reorder

void printdata(struct data* data, char* ltype, char *path){
	static int iter=0;		//keeps track of the iterations to give its name to the output files
	int cloc,call,cpop,		//iterators
		npop,				//number of populations
		nloc,				//number of loci
		outsize;			//size of the output file
	char *outname4;			//name of the output filename

	FILE *out_len;			//pntr to .len
	
	npop = data->npop;
	nloc = data->nloc;
	
	outsize=strlen(path) + 10;
	outname4 = (char *)malloc(outsize*sizeof(char));
	strcpy(outname4,path);
	sprintf(outname4+strlen(outname4),"data");
	sprintf(outname4+strlen(outname4),"%d",iter);
	out_len = fopen(strcat(outname4,".len"),"w");	//output .len
	if(out_len == NULL)
		printerr("cannot create .len file");
	
	++iter;

	fprintf(out_len,"%d\n%d\n\n",npop, nloc);					//out_len: npop, nloc
	for(cloc=0;cloc<nloc;++cloc){
		fprintf(out_len,"%c ",ltype[cloc]);						//out_len: ltype[cloc]
	}
	fprintf(out_len,"\n\n");
			
	for(cloc=0;cloc<nloc;++cloc){
		fprintf(out_len,"%d\n",data->ldna[cloc]);				//out_len: noall[cloc]
		for(cpop=0;cpop<npop;++cpop){
			for(call=0;call<data->ldna[cloc];++call){
				fprintf(out_len,"%d\t",data->freq[cloc][cpop][call]);	//out_len: freq[cloc][cpop][call]
			}
			fprintf(out_len,"\n");
		}
		fprintf(out_len,"\n");
		for(call=0;call<data->ldna[cloc];++call){
			if(ltype[cloc]=='M' || ltype[cloc]=='m')
				fprintf(out_len,"%d\t",data->valM[cloc][call]);	//out_len: vals[cloc][call]
			else
				fprintf(out_len,"%d\t",call);					//out_len: call
		}
		fprintf(out_len,"\n\n");

		if(ltype[cloc]=='S' || ltype[cloc]=='s'){
			fprintf(out_len,"%d\n",data->lsites[cloc]);
			for(call=0;call<data->ldna[cloc];++call)
				fprintf(out_len,"%d %s\n",call,data->valS[cloc][call]);	//out_len: call, vals[cloc][call]
			fprintf(out_len,"\n\n");
		}
	}
	
	/*free stuff*/
	free(outname4);

	/*close files*/
	fclose(out_len);
	
}	//end of printdata

int *posToOrd(double *dna,double *lmut,int len,int nsites){
	int cmut,pos,	//iterator
		*out;		//array with the postion of the mutations

	out = (int *)malloc(len*sizeof(int));

	for(pos=0; pos<len; pos++){
		for(cmut=0; cmut<nsites; cmut++){
			if(dna[pos]==lmut[cmut]){
				out[pos] = cmut;
				break;
			}
		}
	}

	return out;

} //end of posToOrd

void freetree(int nloc,int npop,int nevt){
	int cpop, cloc, cevt;	//iterators
	for(cpop=0;cpop<npop;++cpop){
		free(den[cpop]);
	}
	for(cloc=0;cloc<nloc;++cloc){
		free(nInd2[cloc]);
		free(lind2[cloc]);
	}
	for(cevt=0 ; cevt<=nevt ; ++cevt){
		free(Nevent[cevt]);
 		free(iEvent[cevt]);
		if(cevt < nevt){
			free(listSplit[cevt]);
		}	
	}
	free(listIndex);
	free(listInd);
	free(lnode);
	free(den);
	free(recLength);
	free(migRate);
	free(coalRate);
	free(nInd1);
	free(nInd2);
	free(lind1);
	free(lind2);
	free(tEvent);
	free(Nevent);
	free(iEvent);
	free(listSplit);
		
} //end of freetree
