/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		12th May 2009
*/
#include "interface.h"

/*
	Presents the main menu to user

*/
void mainMenu();

/*
	Simulates data sets

	@return errors
*/
int simulate();

/*
	Run the ABC algorithm

	@return errors
*/
int runAbc();

/*
	Presents the customize ABC menu to user

*/
void customMenu();

/*
	Run the simulation step of ABC algorithm

	@return errors
*/
int runabc1();

/*
	Run the rejection step of ABC algorithm

	@return errors
*/
int runabc2();

/*
	Check for suitability of performing the regression step of ABC Algorithm

	@return errors
*/
int runabc3();

/*
	Presents the tools menu to user

*/
void toolsMenu();

/*
	Build a prior file (.prs)

	@return errors
*/
int buildprior();

/*
	Check if a population still exists
	
	@arg avlpop list of all available populations
	@pop population which existence should be checked
	@size number of existing populations
	
	@return int (0 - if pop doesn't exist; 1 - if pop exists) 
*/
int exists(int *avlpop,int pop, int size);

/*
	Take out a pop from the list of available pop
	
	@arg avlpop list of all available populations
	@pop population to be taken out
	@size number of existing populations
	
	@return updated list of the available pop 
*/
int *takeout(int *avlpop,int pop,int size);

/*
	Build a summary statistics set file (.sst)

	@return errors
*/
int buildsstats();

/*
	Create a target (.trg) and a sample size (.ssz) files from a table file (.len)

	@return errors
*/
int buildtarget();

/*
	Convert a sample file (.pop) into a table file (.len)

	@return errors
*/
int createtab();

/*
	Convert a IMa sample file into a table file (.len)

	@return errors
*/
int createtab2();

/*
	Convert a GenePop file into a table file (.len)

	@return errors
*/
int createtab3();

/*
	Convert a Nexus file (.nex) into a table file (.len)

	@return errors
*/
int createtab4();

/*
	Convert a table file (.len) into a sample file (.pop)

	@return errors
*/
int createpop();

/*
	Join files together

	@return errors
*/
int joinfiles();

/*
	Find segregating sites (.sgs)

	@return errors
*/
int segsites();

/*
	It creates an ASCII interface to run ABC approachs to Populacional Biology problems
*/
int main(int argc, char *argv[]){
	int i,		//iterator
		cut;	//auxiliar to work around strings
		
	/*get the path from where the progam is being run*/
	Exepathsize = strlen(argv[0])+ 5;
	Exepath = (char *)malloc(Exepathsize*sizeof(char));
	strcpy(Exepath,argv[0]);
	for(cut=-1,i=Exepathsize-1 ; Exepath[i]!='/' && Exepath[i]!='\\' && i >= 0 ; cut++,i--){
		Exepath[i]='\0';
	}
	Exepath = realloc(Exepath,(Exepathsize-cut)*sizeof(char));

	printf("\n\n  ----------------------------------\n");
	printf("   popABC - Mark Beaumont & Joao Lopes\n\n");
	printf("   version: 1.0               12/05/09\n");
	printf("\n  ----------------------------------\n");

	mainMenu();
	
	//free stuff
	free(Exepath);

} //end of main

void mainMenu(){
	char str[MAXCHAR];
	int option,
		r;

	do {
	    printf("\n\n  *Main Menu*\n\n");
	    printf("  1. Create a simulated \"real\" data\n\n");
	    printf("  2. Run the ABC algorithm on \"real\" data\n\n");
	    printf("  3. Customize the ABC algorithm\n\n");
	    printf("  4. Tools\n\n");
	    printf("  5. Exit\n\n");
	    printf("  choose an option: ");
		scanf("%s",&str);
		option = atoi(str);
	    switch (option){
	    case 1:
		printf("\n\n\n  *Simulate Data with ABC*\n\n");
		r = simulate();
		if(r==0)
			printf("\n\n>simulation completed (report created).");
		else if(r==1)
			printf("\n\n>error::cannot open the .prs file!");
		else if(r==2)
			printf("\n\n>error::cannot open the .ssz file!");
		else if(r==3)
			printf("\n\n>error::cannot open the .sst file!");
		else if(r==4)
			printf("\n\n>error::cannot create .dat file!");
		else if(r==5)
			printf("\n\n>error::cannot create info file!");
		else if(r==6)
			printf("\n\n>error::while reading .ssz file!");
		else if(r==7)
			printf("\n\n>error::while reading .sst file!");
		else if(r==8)
			printf("\n\n>error::following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)!");
		else if(r==9)
			printf("\n\n>error::wrong topology prior type!");
		else if(r==10)
			printf("\n\n>error::with 1 or 2 pops only topology prior type 0 or 3 can be used!");
		else if(r==11)
			printf("\n\n>error::wrong Ne prior type!");
		else if(r==12)
			printf("\n\n>error::type-2 tev is followed by another type-2 tev!");
		else if(r==13)
			printf("\n\n>error::if using tev type 3 prior, define prior only once!");
		else if(r==14)
			printf("\n\n>error::if using tev type 4 prior, define prior only once!");
		else if(r==15)
			printf("\n\n>error::wrong Tev prior type!");
		else if(r==16)
			printf("\n\n>error::wrong Mig prior type!");
		else if(r==17)
			printf("\n\n>error::wrong Mut prior type (Microsatellites)!");
		else if(r==18)
			printf("\n\n>error::wrong Mut prior type (Sequence Data)!");
		else if(r==19)
			printf("\n\n>error::wrong Rec prior type (Microsatellites)!");
		else if(r==20)
			printf("\n\n>error::wrong Rec prior type (Sequence Data)!");
		else if(r==21)
			printf("\n\n>error::to use migration weights matrix, topology must be defined (option 1,2,4 or 5)!");
		else if(r==22)
			printf("\n\n>error::wrong migration weights matrix option!");
		else if(r==23)
			printf("\n\n>error::.dat file would be too big to be analyzed!");
		else if(r==24)
			printf("\n\n>error::when using rec with STR all STR loci must have the same scalar!");
		else if(r==25)
			printf("\n\n>error::when using rec with STR all STR loci must have the same samples per population!");
		else if(r==26)
			printf("\n\n>error::cannot create .mut file!");
		else if(r==27)
			printf("\n\n>error::cannot create .rec file!");
		else if(r==28)
			printf("\n\n>error::times are not monotonic increasing!");
		else if(r==29)
			printf("\n\n>error::choosing branching sequence!");
		else if(r==30)
			printf("\n\n>error::with more than 5 populations the topology as to be fixed using option 2 or 5!");
		else if(r==31)
			printf("\n\n>error::the locus have less than samples, it should be cleared!");
		else if(r==32)
			printf("\n\n>error::migration weights for the considered population have to be 0!");
		else if(r==33)
			printf("\n\n>error::the sum of the prob for a given migration as to be 1.0!");
		else if(r==34)
			printf("\n\n>error::cannot create .len file!");
		break;

	    case 2:
		printf("\n\n\n  *Run ABC algorithm*\n\n");
		r = runAbc();
		if(r==0)
			printf("\n\n>ABC run completed\n");
		else if(r==1)
			printf("\n\n>error::cannot open .len file!");	
		else if(r==2)
			printf("\n\n>error::cannot open .sst file!");	
		else if(r==3)
			printf("\n\n>error::cannot create .trg file!");		
		else if(r==4)
			printf("\n\n>error::cannot create .ssz file!");
		else if(r==5)
			printf("\n\n>error::cannot create info file!");
		else if(r==6)
			printf("\n\n>error::reading .sst file!");
		else if(r==7)
			printf("\n\n>error::following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)!");
		else if(r==8)
			printf("\n\n>error::sequence data locus without information, discard it!");
		else if(r==101)
			printf("\n\n>error::cannot open the .prs file!");
		else if(r==102)
			printf("\n\n>error::cannot open the .ssz file!");
		else if(r==103)
			printf("\n\n>error::cannot open the .sst file!");
		else if(r==104)
			printf("\n\n>error::cannot create .dat file!");
		else if(r==105)
			printf("\n\n>error::cannot create info file!");
		else if(r==106)
			printf("\n\n>error::while reading .ssz file!");
		else if(r==107)
			printf("\n\n>error::while reading .sst file!");
		else if(r==108)
			printf("\n\n>error::following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)!");
		else if(r==109)
			printf("\n\n>error::wrong topology prior type!");
		else if(r==110)
			printf("\n\n>error::with 1 or 2 pops only topology prior type 0 or 3 can be used!");
		else if(r==111)
			printf("\n\n>error::wrong Ne prior type!");
		else if(r==112)
			printf("\n\n>error::type-2 tev is followed by another type-2 tev!");
		else if(r==113)
			printf("\n\n>error::if using tev type 3 prior, define prior only once!");
		else if(r==114)
			printf("\n\n>error::if using tev type 4 prior, define prior only once!");
		else if(r==115)
			printf("\n\n>error::wrong Tev prior type!");
		else if(r==116)
			printf("\n\n>error::wrong Mig prior type!");
		else if(r==117)
			printf("\n\n>error::wrong Mut prior type (Microsatellites)!");
		else if(r==118)
			printf("\n\n>error::wrong Mut prior type (Sequence Data)!");
		else if(r==119)
			printf("\n\n>error::wrong Rec prior type (Microsatellites)!");
		else if(r==120)
			printf("\n\n>error::wrong Rec prior type (Sequence Data)!");
		else if(r==121)
			printf("\n\n>error::to use migration weights matrix, topology must be defined (option 1,2,4 or 5)!");
		else if(r==122)
			printf("\n\n>error::wrong migration weights matrix option!");
		else if(r==123)
			printf("\n\n>error::.dat file would be too big to be analyzed!");
		else if(r==124)
			printf("\n\n>error::when using rec with STR all STR loci must have the same scalar!");
		else if(r==125)
			printf("\n\n>error::when using rec with STR all STR loci must have the same samples per population!");
		else if(r==126)
			printf("\n\n>error::cannot create .mut file!");
		else if(r==127)
			printf("\n\n>error::cannot create .rec file!");
		else if(r==128)
			printf("\n\n>error::times are not monotonic increasing!");
		else if(r==129)
			printf("\n\n>error::choosing branching sequence!");
		else if(r==130)
			printf("\n\n>error::with more than 5 populations the topology as to be fixed using option 2 or 5!");
		else if(r==131)
			printf("\n\n>error::the locus have less than samples, it should be cleared!");
		else if(r==132)
			printf("\n\n>error::migration weights for the considered population have to be 0!");
		else if(r==133)
			printf("\n\n>error::the sum of the prob for a given migration as to be 1.0!");
		else if(r==134)
			printf("\n\n>error::cannot create .len file!");
		else if(r==201)
			printf("\n\n>error::could't open .dat file!");
		else if(r==202)
			printf("\n\n>error::cannot open .trg file!");
		else if(r==203)
			printf("\n\n>error::cannot create .rej file!");
		else if(r==204)
			printf("\n\n>error::cannot create .pri file!");
		else if(r==205)
			printf("\n\n>error::cannot create .txt file!");
		else if(r==206)
			printf("\n\n>error::file .dat ended prematurely!");
		else if(r==207)
			printf("\n\n>error::.dat file is too big to be analysed!");
		break;
		
	    case 3:
		customMenu();
		break;

	    case 4:
		toolsMenu();
		break;

	    case 5:
		printf("\n\n...exiting popABC\n\n");
		break;

	    default:
		printf("\n\n>error::option not available\n");
		break;
	    }
	}
	while (option != 5);

} // end of mainMenu

int simulate(){
	char str[MAXCHAR],			//auxiliar to get the input from the user
		 input_prs[MAXCHAR],	//input file name
		 input_ssz[MAXCHAR],	//input file name
		 input_sst[MAXCHAR],	//input file name
		 output[MAXCHAR];		//input file name

	//input filename
	printf("\n  choose prior file (.prs) (max %d character)",MAXCHAR-1);
	printf("\n  SHOULD HAVE ONLY ONE ITERATION: ");
	scanf("%s",input_prs);

	//input filename
	printf("\n  choose sample size file (.ssz) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_ssz);

	//input filename
	printf("\n  choose summary statistics set file (.sst) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_sst);

	//output filename
	printf("\n  insert name for output file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	return abc(input_prs,input_ssz,input_sst,output,1,1,1);

}  //ends of simulate

int runAbc(){
	int	lsstats[MAXSSTATS],		//list of the summary statistics used
	 	printMut,				//choose if the mutation rates should be print
		printRec,				//choose if the recombination rates should be print
		nparam,					//number of parameters
		nstats,					//number of summary statistics
		error;					//gets any error
	double tol;					//tolerance for the rejection method
	char str[MAXCHAR],			//auxiliar to get the input from the user
		 input_len[MAXCHAR],	//.len file name
		 input_trg[MAXCHAR],	//.trg file name
		 input_prs[MAXCHAR],	//.prs file name
		 input_ssz[MAXCHAR],	//.ssz file name
		 input_sst[MAXCHAR],	//input file name
		 input_dat[MAXCHAR],	//.dat file name
		 input_rej[MAXCHAR];	//.rej file name

	//input filename
	printf("\n  choose table file (.len) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_len);

	//input filename
	printf("\n  choose summary statistics set file (.sst) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_sst);

	//input filename
	printf("\n  insert name for target file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_trg);

	error = maketarget(input_len,input_sst,input_trg);
	if(error)
		return error;
	else
		printf("\n\n>'real' data summarization completed (report created).\n\n");
			
	//.prs filename
	printf("\n  choose prior file (.prs) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_prs);

	//.dat filename
	printf("\n  insert name for data file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_dat);
	
	//choose if file with mutation rate will be created
	printf("\n");
	do{
		printf("  choose if file with mutation rate will be created (0 - file not created;");
		printf("\n  1 - file created): ");
		scanf("%s",&str);
		printMut = atoi(str);
	}
	while(printMut !=0 && printMut !=1);

	//choose if file with recombination rate will be created
	printf("\n");
	do{
		printf("  choose if file with recombination rate will be created (0 - file not created;");
		printf("\n  1 - file created): ");
		scanf("%s",&str);
		printRec = atoi(str);
	}
	while(printRec !=0 && printRec !=1);

	strcpy(input_ssz,input_trg);
	strcat(input_ssz,".ssz"); 

	error = abc(input_prs,input_ssz,input_sst,input_dat,0,printMut,printRec);
	if(error)
		return error+100;
	else printf("\n\n>simulations completed (report created).\n\n");

	//.rej filename
	printf("\n  insert name for rejection file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_rej);

	//choose number of parameters
	printf("\n");
	do{
		printf("  specify the number of parameters used (>0): ");
		scanf("%s",&str);
		nparam = atof(str);
	}
	while(nparam <=0);

	//choose number of summary statistics
	printf("\n");
	do{
		printf("  specify the number of summary statistics used (>0): ");
		scanf("%s",&str);
		nstats = atof(str);
	}
	while(nstats <=0);

	//choose tolerance
	printf("\n");
	do{
		printf("  choose tolerance for rejection step (0 < tol < 1): ");
		scanf("%s",&str);
		tol = atof(str);
	}
	while(tol <=0 || tol >=1);

	strcat(input_dat,".dat");
	strcat(input_trg,".trg");

	error = firstpass(input_dat,input_trg,input_rej,nparam,nstats,tol);
	if(error)
		return error+200;
	else
		printf("\n\n>rejection step completed (report created).\n\n");

	return 0;
	
}  //end of runABC

void customMenu(){
	char str[MAXCHAR];
	int option,
		r;

	do {
	    printf("\n\n\n  *Customize ABC Menu*\n\n");
	    printf("  1. Summarize 'real' data with summstats set\n\n");
	    printf("  2. Run simulations\n\n");
	    printf("  3. Perform rejection step\n\n");
	    printf("  4. Exit to Main Menu\n\n");
	    printf("  choose an option: ");
		scanf("%s",&str);
		option = atoi(str);
	    switch (option){
	    case 1:
		printf("\n\n\n  *Run ABC - 'real' data summarization step*\n\n");
		r = buildtarget();
		if(r==0)
			printf("\n\n>'real' data summarization completed  (report created).");
		else if(r==1)
			printf("\n\n>error::cannot open the .len file!");	
		else if(r==2)
			printf("\n\n>error::cannot open the .sst file!");	
		else if(r==3)
			printf("\n\n>error::cannot create .trg file!");		
		else if(r==4)
			printf("\n\n>error::cannot create .ssz file!");
		else if(r==5)
			printf("\n\n>error::cannot create info file!");
		else if(r==6)
			printf("\n\n>error::reading .sst file!");
		else if(r==7)
			printf("\n\n>error::following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)!");
		else if(r==8)
			printf("\n\n>error::sequence data locus without information, discard it!");
		break;

	    case 2:
		printf("\n\n\n  *Run ABC - simulation step*\n\n");
		r = runabc1();
		if(r==0)
			printf("\n\n>simulation step completed (report created).");
		else if(r==1)
			printf("\n\n>error::cannot open the .prs file!");
		else if(r==2)
			printf("\n\n>error::cannot open the .ssz file!");
		else if(r==3)
			printf("\n\n>error::cannot open the .sst file!");
		else if(r==4)
			printf("\n\n>error::cannot create .dat file!");
		else if(r==5)
			printf("\n\n>error::cannot create info file!");
		else if(r==6)
			printf("\n\n>error::while reading .ssz file!");
		else if(r==7)
			printf("\n\n>error::while reading .sst file!");
		else if(r==8)
			printf("\n\n>error::following sstats can only be chosen with 2 or more pops: Nm_H; Nm_S; privS; S(1)!");
		else if(r==9)
			printf("\n\n>error::wrong topology prior type!");
		else if(r==10)
			printf("\n\n>error::with 1 or 2 pops only topology prior type 0 or 3 can be used!");
		else if(r==11)
			printf("\n\n>error::wrong Ne prior type!");
		else if(r==12)
			printf("\n\n>error::type-2 tev is followed by another type-2 tev!");
		else if(r==13)
			printf("\n\n>error::if using tev type 3 prior, define prior only once!");
		else if(r==14)
			printf("\n\n>error::if using tev type 4 prior, define prior only once!");
		else if(r==15)
			printf("\n\n>error::wrong Tev prior type!");
		else if(r==16)
			printf("\n\n>error::wrong Mig prior type!");
		else if(r==17)
			printf("\n\n>error::wrong Mut prior type (Microsatellites)!");
		else if(r==18)
			printf("\n\n>error::wrong Mut prior type (Sequence Data)!");
		else if(r==19)
			printf("\n\n>error::wrong Rec prior type (Microsatellites)!");
		else if(r==20)
			printf("\n\n>error::wrong Rec prior type (Sequence Data)!");
		else if(r==21)
			printf("\n\n>error::to use migration weights matrix, topology must be defined (option 1,2,4 or 5)!");
		else if(r==22)
			printf("\n\n>error::wrong migration weights matrix option!");
		else if(r==23)
			printf("\n\n>error::.dat file would be too big to be analyzed!");
		else if(r==24)
			printf("\n\n>error::when using rec with STR all STR loci must have the same scalar!");
		else if(r==25)
			printf("\n\n>error::when using rec with STR all STR loci must have the same samples per population!");
		else if(r==26)
			printf("\n\n>error::cannot create .mut file!");
		else if(r==27)
			printf("\n\n>error::cannot create .rec file!");
		else if(r==28)
			printf("\n\n>error::times are not monotonic increasing!");
		else if(r==29)
			printf("\n\n>error::choosing branching sequence!");
		else if(r==30)
			printf("\n\n>error::with more than 5 populations the topology as to be fixed using option 2 or 5!");
		else if(r==31)
			printf("\n\n>error::the locus have less than samples, it should be cleared!");
		else if(r==32)
			printf("\n\n>error::migration weights for the considered population have to be 0!");
		else if(r==33)
			printf("\n\n>error::the sum of the prob for a given migration as to be 1.0!");
		else if(r==34)
			printf("\n\n>error::cannot create .len file!");
		break;

	    case 3:
		printf("\n\n\n  *Run ABC - rejection step*\n\n");
		r = runabc2();
		if(r==0)
			printf("\n\n>rejection step completed (report created).");
		else if(r==1)
			printf("\n\n>error::cannot open .dat file!");
		else if(r==2)
			printf("\n\n>error::cannot open .trg file!");
		else if(r==3)
			printf("\n\n>error::cannot create .rej file!");
		else if(r==4)
			printf("\n\n>error::cannot create .pri file!");
		else if(r==5)
			printf("\n\n>error::cannot create .txt file!");
		else if(r==6)
			printf("\n\n>error::file .dat ended prematurely!");
		else if(r==7)
			printf("\n\n>error::.dat file is too big to be analysed!");
		break;

	    case 4:
		break;

	    default:
		printf("\n\n>error::option not available\n");
		break;
	    }
	}
	while (option != 4);
	
}  //end of customMenu

int buildtarget(){
	int cstat,i,			//iterators
		lsstats[MAXSSTATS],	//list of the summary statistics used
	 	nmarks,				//param number of markers used (0-no marks; 2-two marks; 3-three marks; 4- four marks)
	 	cmark,				//auxiliar to choose the sample populations
	 	mark1,				//marker of first sample population (only if npop >= 3 AND nsample == 2)
	 	mark2,				//marker of second sample population (only if npop >= 3 AND nsample == 2 OR 3 OR 4)
	 	mark3;				//marker of third sample population (only if npop >= 3 AND nsample == 3 OR 4)
	char str[MAXCHAR],		//auxiliar to get the input from the user
		 input_len[MAXCHAR],//input file name
		 input_sst[MAXCHAR],//input file name
		 output[MAXCHAR],	//input file name
		 type;				//analysis type (M - microsatellites; S - sequence data; B - both)

	
	//input filename
	printf("\n  choose a table file (.len) (max %d characters): ",MAXCHAR-1);
	scanf("%s",input_len);

	//input filename
	printf("\n  choose a summary statistics set file (.sst) (max %d characters): ",MAXCHAR-1);
	scanf("%s",input_sst);

	//output filename
	printf("\n  choose a name for the target file (no extension) (max %d",MAXCHAR-1);
	printf("\n  characters): ");
	scanf("%s",output);

	return maketarget(input_len,input_sst,output);
	
} //end of buildtarget

int runabc1(){
	int printIt,				//choose if the table files of the simulations should be printed
		printMut,				//choose if the mutation rates should be printed
		printRec;				//choose if the recombination rates should be printed
	char str[MAXCHAR],			//auxiliar to get the input from the user
		 input_prs[MAXCHAR],	//input file name
		 input_ssz[MAXCHAR],	//input file name
		 input_sst[MAXCHAR],	//input file name
		 output[MAXCHAR];		//input file name

	//input filename
	printf("\n  choose prior file (.prs) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_prs);

	//input filename
	printf("\n  choose sample size file (.ssz) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_ssz);

	//input filename
	printf("\n  choose summary statistics set file (.sst) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_sst);

	//output filename
	printf("\n  insert name for output file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	//choose if table files of the simulations should be printed
	printf("\n");
	do{
		printf("  choose if table files (.len) for each simulation will be created (0 - file");
		printf("\n  not created;1 - file created)");
		printf("\n  THIS OPTION MEANS ONE .LEN FILE PER SIMULATION: ");
		scanf("%s",&str);
		printIt = atoi(str);
	}
	while(printIt !=0 && printIt !=1);

	//choose if file with mutation rate will be created
	printf("\n");
	do{
		printf("  choose if file with mutation rate will be created (0 - file not created;");
		printf("\n  1 - file created): ");
		scanf("%s",&str);
		printMut = atoi(str);
	}
	while(printMut !=0 && printMut !=1);

	//choose if file with recombination rate will be created
	printf("\n");
	do{
		printf("  choose if file with recombination rate will be created (0 - file not created;");
		printf("\n  1 - file created): ");
		scanf("%s",&str);
		printRec = atoi(str);
	}
	while(printRec !=0 && printRec !=1);
		
	return abc(input_prs,input_ssz,input_sst,output,printIt,printMut,printRec);
	
} //end of runabc1

int runabc2(){
	int	nparam,					//number of parameters
	 	nstats;					//number of summary statistics
	double tol;					//tolerance for the rejection method
	char str[MAXCHAR],			//auxiliar to get the input from the user
		 input_dat[MAXCHAR],	//.dat file name
		 input_trg[MAXCHAR],	//.trg file name
		 input_rej[MAXCHAR];	//.rej file name

	//.dat filename
	printf("\n  choose data file (.dat) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_dat);

	//.trg filename
	printf("\n  choose target file (.trg) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_trg);

	//.rej filename
	printf("\n  insert name for rejection file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",input_rej);

	//choose number of parameters
	printf("\n");
	do{
		printf("  specify the number of parameters used (>0): ");
		scanf("%s",&str);
		nparam = atof(str);
	}
	while(nparam <=0);

	//choose number of summary statistics
	printf("\n");
	do{
		printf("  specify the number of summary statistics used (>0): ");
		scanf("%s",&str);
		nstats = atof(str);
	}
	while(nstats <=0);

	//choose tolerance
	printf("\n");
	do{
		printf("  choose tolerance for rejection step (0 < tol < 1): ");
		scanf("%s",&str);
		tol = atof(str);
	}
	while(tol <=0 || tol >=1);

	return firstpass(input_dat,input_trg,input_rej,nparam,nstats,tol);
	
} //end of runabc2

void toolsMenu(){
	char str[MAXCHAR];
	int option,
		r;

	do {
	    printf("\n\n\n  *Tools Menu*\n\n");
	    printf("  1. Build prior file (.prs)\n\n");
	    printf("  2. Build summary statistics set file (.sst)\n\n");
	    printf("  3. Convert sample file (.pop) to table file (.len)\n\n");
	    printf("  4. Convert IMa input file to table file (.len)\n\n");
	    printf("  5. Convert GenePop input file to table file (.len)\n\n");
	    printf("  6. Convert Nexus file (.nex) to table file (.len)\n\n");
	    printf("  7. Create sample file (.pop) from table file (.len)\n\n");
	    printf("  8. Join data files together (.dat)\n\n");
	    printf("  9. Exit to Main Menu\n\n");
	    printf("  choose an option: ");
		scanf("%s",&str);
		option = atoi(str);
	    switch (option)
	    {
	    case 1:
		printf("\n\n\n  *Build prior file (.prs)*\n\n");
		r = buildprior();
		if(r==0)
			printf("\n\n>prior file (.prs) created.");
		else if(r==1)
			printf("\n\n>error::cannot create .prs file!");	
		break;

	    case 2:
		printf("\n\n\n  *Build summary statistics set file (.prs)*\n\n");
		r = buildsstats();
		if(r==0)
			printf("\n\n>summary statistics set file (.sst) created.");
		else if(r==1)
			printf("\n\n>error::cannot create .sst file!");	
		break;

	    case 3:
		printf("\n\n\n  *Convert sample file (.pop) to table file (.len)*\n\n");
		r = createtab();
		if(r==0)
			printf("\n\n>table file (.len) build.");
		else if(r==1)
			printf("\n\n>error::cannot open .pop file!");	
		else if(r==2)
			printf("\n\n>error::cannot create .len file!");	
		else if(r==3)
			printf("\n\n>error::individual has no information!");
		else if(r==4)
			printf("\n\n>error::different number of ploidy or number of sites in alleles of a loci!");
		break;

	    case 4:
		printf("\n\n\n  *Convert IMa file to table file (.len)*\n\n");
		r = createtab2();
		if(r==0)
			printf("\n\n>table file (.len) build.");
		else if(r==1)
			printf("\n\n>error::cannot open sample file!");	
		else if(r==2)
			printf("\n\n>error::cannot create .len file!");	
		else if(r==3)
			printf("\n\n>error::cannot deal with specified data type!");	
		break;

	    case 5:
		printf("\n\n\n  *Convert GenePop file to table file (.len)*\n\n");
		r = createtab3();
		if(r==0)
			printf("\n\n>table file (.len) build.");
		else if(r==1)
			printf("\n\n>error::cannot open sample file!");	
		else if(r==2)
			printf("\n\n>error::cannot create .len file!");	
		else if(r==3)
			printf("\n\n>error::file ended before the display of data!");	
		else if(r==4)
			printf("\n\n>error::incorrect number of digits in locus!");	
		break;

	    case 6:
		printf("\n\n\n  *Convert Nexus file (.nex) to table file (.len)*\n\n");
		r = createtab4();
		if(r==0)
			printf("\n\n>table file (.len) build.");
		else if(r==1)
			printf("\n\n>error::cannot open .nex file!");	
		else if(r==2)
			printf("\n\n>error::cannot create .len file!");	
		else if(r==3)
			printf("\n\n>error::cannot deal with specified data type!");	
		else if(r==4)
			printf("\n\n>error::problems reading MATRIX!");	
		else if(r==5)
			printf("\n\n>error::problems reading TAXSET!");	
		else if(r==6)
			printf("\n\n>error::problem defining TAXSET: n1>=n2!");	
		else if(r==7)
			printf("\n\n>error::unrecognized specifier after BEGIN SETS (use TAXSET)!");	
		else if(r==8)
			printf("\n\n>error::unrecognized specifier after BEGIN (use DATA or SETS)!");	
		else if(r==9)
			printf("\n\n>error::unrecognized command (only BEGIN accepted)!");	
		break;
		
	    case 7:
		printf("\n\n\n  *Create sample file (.pop) from table file (.len)*\n\n");
		r = createpop();
		if(r==0)
			printf("\n\n>sample file (.pop) converted.");
		else if(r==1)
			printf("\n\n>error::cannot open .pop files!");	
		else if(r==2)
			printf("\n\n>error::cannot create .len file!");		
		break;

	    case 8:
		printf("\n\n\n  *Join data files (.dat) together*\n\n");
		r = joinfiles();
		if(r==0)
			printf("\n\n>files joined.");
		else if(r==1)
			printf("\n\n>error::cannot create output file!");	
		else if(r==2)
			printf("\n\n>error::cannot open an input files!");		
		break;

	    case 9:
		break;

	    default:
		printf("\n\n>error::option not available\n");
		break;
	    }
	}
	while (option != 9);

}  //end of toolsMenu

int buildprior(){
	int cloc,cpop,ctev,i,j,cpop2,	//iterator
		allTev,						//check if only on prior for all tevs is used
		error,						//saves the error type
		d1,							//auxiliar to get info from the user
		from,to,					//auxiliars to construct the fixed topology when the option 2 is chosen
		niter,						//number of iterations
		genet,						//generation time
		npop,						//number of populations
		ntev,						//number of time events
		psize,						//number of populations (used in topology prior type 2)
		nloc,						//number of loci
		*avlpop;					//populations available to chose from when fixing the topology with option 2
	double *lplo;					//list of the ploidy per loci
	char str[MAXCHAR],				//auxiliar to get the input from the user
		 output[MAXCHAR],			//output filename
		 *ltype;					//list of the DNA types analysed per loci
	struct prior pr_top,			//topology prior parameters
				 pr_mutSTR,			//STR mutation prior parameters
				 pr_mutSNP,			//SNP mutation prior parameters
				 pr_recSTR,			//STR recombination prior parameters
				 pr_recSNP,			//SNP recombination prior parameters
			 	 *pr_mig,			//list of the mig prior parameters
				 *pr_Ne,			//list of the Ne prior parameters
			 	 *pr_tev;			//list of the tev prior paramters
	struct migweights migw; 		//matrix with the migration weigths
	
	/*number of iterations*/
	do{
		printf("  choose the number of iterations (>0): ");
		scanf("%s",&str);
		niter = atoi(str);
	}
	while(niter <=0);
	/*generation time*/
	printf("\n");
	do{
		printf("  choose the generation time (>0): ");
		scanf("%s",&str);
		genet = atoi(str);
	}
	while(genet <=0);
	/*number of pop*/
	printf("\n");
	do{
		printf("  choose the number of populations (>0): ");
		scanf("%s",&str);
		npop = atoi(str);
	}
	while(npop <=0);
	/*number of loci*/
	printf("\n");
	do{
		printf("  choose the number of loci (>0): ");
		scanf("%s",&str);
		nloc = atoi(str);
	}
	while(nloc <=0);
	/*inherety scalar*/
	lplo = (double*)malloc(nloc*sizeof(double));
	if(nloc==1){
		printf("\n");
		do{
			printf("  choose the hereditary scalar of the locus [autossome = 1.00; X-linked = 0.75;\n");
			printf("  haploid = 0.5; mitochodrial DNA or Y-linked = 0.25]: ");
			scanf("%s",&str);
			lplo[0] = atof(str);
		}
		while(lplo[0]!=1.0 && lplo[0]!=0.75 && lplo[0]!=0.5 && lplo[0]!=0.25);
	}
	else{
		printf("\n  choose %d scalars, one for each loci [autossome = 1.00; X-linked = 0.75;\n",nloc);
		printf("  haploid = 0.5; mitochodrial DNA or Y-linked = 0.25]:\n");
		for(cloc=0; cloc<nloc; cloc++){
			do{
				printf("    locus %d: ",cloc+1);
				scanf("%s",&str);
				lplo[cloc] = atof(str);
			}
		while(lplo[0]!=1.0 && lplo[0]!=0.75 && lplo[0]!=0.5 && lplo[0]!=0.25);
		}
	}
	/*DNA types*/
	ltype = (char*)malloc(nloc*sizeof(char));
	if(nloc==1){
		printf("\n");
		do{	
			printf("  choose the DNA type of the locus [s - DNA sequence or m - microsatellites]: ");
			scanf("%s",&str);
			ltype[0] = str[0];
		}
		while(ltype[0]!='m' && ltype[0]!='s' && ltype[0]!='M' && ltype[0]!='S');
	}
	else{
		printf("\n  choose the DNA type of each locus [s - DNA sequence or m - microsatellites]:\n");
		for(cloc=0; cloc<nloc; cloc++){
			do{
				printf("    locus %d: ",cloc+1);
				scanf("%s",&str);
				ltype[cloc] = str[0];
			}
			while(ltype[cloc]!='m' && ltype[cloc]!='s' && ltype[cloc]!='M' && ltype[cloc]!='S');
		}
	}
	/*topology priors*/
	if(npop<=2){
		printf("\n");
		do{
			printf("  choose to add a marker to the simulations [0 - no marker; 3 - marker]: ");
			scanf("%s",&str);
			pr_top.type = atoi(str);
		}
		while(pr_top.type!=0 && pr_top.type!=3);
	}
	else if(npop<=5){
		printf("\n");
		do{
			printf("  choose the topology prior and to add or not a marker to the simulations   \n");
			printf("  [No marker: 0 - random; 1 - predefined topology; 2 - write topology]      \n");
			printf("  [With marker: 3 - random; 4 - predefined topology; 5 - write topology]: ");
			scanf("%s",&str);
			pr_top.type = atoi(str);
		}
		while(pr_top.type!=0 && pr_top.type!=1 && pr_top.type!=2 && pr_top.type!=3 && pr_top.type!=4 && pr_top.type!=5);
	}
	else{
		printf("\n");
		do{
			printf("  choose to add a marker to the simulations [2 - no marker; 5 - marker]: ");
			scanf("%s",&str);
			pr_top.type = atoi(str);
		}
		while(pr_top.type!=2 && pr_top.type!=5);
	}

	if(pr_top.type==1){
		pr_top.p = (double *)malloc(sizeof(double));
		if(npop == 3){
			printf("\n");
			do{
				printf("  choose the predifined topology [1 to 3 (See user-guide)]: ");
				scanf("%s",&str);
				pr_top.p[0] = atof(str);
			}
			while(pr_top.p[0]<1 || pr_top.p[0]>3);
		}
		else if(npop == 4){
			printf("\n");
			do{
				printf("  choose the predifined topology [1 to 18 (See user-guide)]: ");
				scanf("%s",&str);
				pr_top.p[0] = atoi(str);
			}
			while(pr_top.p[0]<1 || pr_top.p[0]>18);
		}
		else if(npop == 5){
			printf("\n");
			do{
				printf("  choose the predifined topology [1 to 180 (See user-guide)]: ");
				scanf("%s",&str);
				pr_top.p[0] = atoi(str);
			}
			while(pr_top.p[0]<1 || pr_top.p[0]>180);
		}
	}
	else if(pr_top.type==2){
		pr_top.p = (double *)malloc(2*(npop+1)*sizeof(double));
		avlpop = (int*)malloc(npop*sizeof(int));
		for(i=0; i<npop; i++)
			avlpop[i] = i;
		psize = npop;
		printf("\n  choose the topology to be used (%d numbers)",2*(npop-1));
		for(i=0; i<2*(npop-1) ;i=i+2){
			printf("\n  available populations:(");
			for(j=0; j<psize-1; j++)
				printf("%d ",avlpop[j]);
			printf("%d)\n",avlpop[j]);
			do{
				printf("    from - ");
				scanf("%s",&str);
				from = atoi(str);
				printf("    to - ");
				scanf("%s",&str);
				to = atoi(str);
			}
			while(!exists(avlpop,from,npop) || !exists(avlpop,to,npop) || from == to);
		 	pr_top.p[i] = from;
		 	pr_top.p[i+1] = to;
		 	avlpop = takeout(avlpop,from,psize);
		 	psize --;
		 }
	}
	else if(pr_top.type==3){
		pr_top.p = (double *)malloc(sizeof(double));
		if(npop <= 4){
			printf("\n");
			do{
				printf("  choose the marker to be used (>=0): ");
				scanf("%s",&str);
				pr_top.p[0] = atof(str);
			}
			while(pr_top.p[0]<0);
		}
	}
	else if(pr_top.type==4){
		pr_top.p = (double *)malloc(2*sizeof(double));
		if(npop == 3){
			printf("\n");
			do{
				printf("  choose the predifined topology [1 to 3 (See user-guide)]: ");
				scanf("%s",&str);
				pr_top.p[0] = atof(str);
			}
			while(pr_top.p[0]<1 || pr_top.p[0]>3);
			printf("\n");
			do{
				printf("  choose the marker to be used (>=0): ");
				scanf("%s",&str);
				pr_top.p[1] = atof(str);
			}
			while(pr_top.p[1]<0);
		}
		else if(npop==4){
			printf("\n");
			do{
				printf("  choose the predifined topology [1 to 18 (See user-guide)]: ");
				scanf("%s",&str);
				pr_top.p[0] = atoi(str);
			}
			while(pr_top.p[0]<1 || pr_top.p[0]>18);
			printf("\n");
			do{
				printf("  choose the marker to be used (>=0): ");
				scanf("%s",&str);
				pr_top.p[1] = atof(str);
			}
			while(pr_top.p[1]<0);
		}
		else if(npop==5){
			printf("\n");
			do{
				printf("  choose the predifined topology [1 to 180 (See user-guide)]: ");
				scanf("%s",&str);
				pr_top.p[0] = atoi(str);
			}
			while(pr_top.p[0]<1 || pr_top.p[0]>180);
			printf("\n");
			do{
				printf("  choose the marker to be used (>=0): ");
				scanf("%s",&str);
				pr_top.p[1] = atof(str);
			}
			while(pr_top.p[1]<0);
		}
	}
	else if(pr_top.type==5){
		pr_top.p = (double *)malloc((1+2*(npop+1))*sizeof(double));
		avlpop = (int*)malloc(npop*sizeof(int));
		for(i=0; i<npop; i++)
			avlpop[i] = i;
		psize = npop;
		printf("\n  choose the topology to be used (%d numbers)",2*(npop-1));
		for(i=0; i<2*(npop-1) ;i=i+2){
			printf("\n  available populations:(");
			for(j=0; j<psize-1; j++)
				printf("%d ",avlpop[j]);
			printf("%d)\n",avlpop[j]);
			do{
				printf("    from - ");
				scanf("%s",&str);
				from = atoi(str);
				printf("    to - ");
				scanf("%s",&str);
				to = atoi(str);
			}
			while(!exists(avlpop,from,npop) || !exists(avlpop,to,npop) || from == to);
		 	pr_top.p[i] = from;
		 	pr_top.p[i+1] = to;
		 	avlpop = takeout(avlpop,from,psize);
		 	psize --;
		 }
		printf("\n");
		do{
			printf("  choose the marker to be used (>=0): ");
			scanf("%s",&str);
			pr_top.p[i] = atof(str);
		}
		while(pr_top.p[i]<0);
	}
	/*Ne priors*/
	pr_Ne = (struct prior *)malloc((npop*2-1)*sizeof(struct prior));
	for(cpop=0; cpop<npop*2-1; cpop++){
		if(cpop<npop){
			printf("\n");
			do{
				printf("  choose the Ne%d prior [1 - uniform distribution; 2 - generalized gamma\n",cpop+1);
		        printf("  distribution]: ");
				scanf("%s",&str);
				pr_Ne[cpop].type = atoi(str);
			}
			while(pr_Ne[cpop].type!=1 && pr_Ne[cpop].type!=2 && pr_Ne[cpop].type!=3);
		}
		else{
			printf("\n");
			do{
				printf("  choose the NeAnc%d prior [1 - uniform distribution; 2 - generalized gamma\n",cpop+1-npop);
				printf("  distribution]: ");
				scanf("%s",&str);
				pr_Ne[cpop].type = atoi(str);
			}
			while(pr_Ne[cpop].type!=1 && pr_Ne[cpop].type!=2 && pr_Ne[cpop].type!=3);
		}
		if(pr_Ne[cpop].type==1){
			pr_Ne[cpop].p = (double *)malloc(2*sizeof(double));
			do{
				printf("  choose lower boundary of the prior(>=0): ");
				scanf("%s",&str);
				pr_Ne[cpop].p[0] = atof(str);
				printf("  choose upper boundary of the prior(>0): ");
				scanf("%s",&str);
				pr_Ne[cpop].p[1] = atof(str);
			}
			while(pr_Ne[cpop].p[0]>pr_Ne[cpop].p[1] || pr_Ne[cpop].p[0]<0 || pr_Ne[cpop].p[1]<=0);

		}
		else{
			pr_Ne[cpop].p = (double *)malloc(4*sizeof(double));
			do{
				printf("  choose location parameter (=>0): ");
				scanf("%s",&str);
				pr_Ne[cpop].p[0] = atof(str);
			}
			while(pr_Ne[cpop].p[0]<0);
			do{
				printf("  choose scale parameter (>0): ");
				scanf("%s",&str);
				pr_Ne[cpop].p[1] = atof(str);
			}
			while(pr_Ne[cpop].p[1]<=0);
			do{
				printf("  choose shape parameter (>0): ");
				scanf("%s",&str);
				pr_Ne[cpop].p[2] = atof(str);
			}
			while(pr_Ne[cpop].p[2]<=0);
			do{
				printf("  choose 2nd shape parameter (>0): ");
				scanf("%s",&str);
				pr_Ne[cpop].p[3] = atof(str);
			}
			while(pr_Ne[cpop].p[3]<=0);
		}
	}
	if(npop>1){
		/*tev priors*/
		ntev = npop-1;
		pr_tev = (struct prior *)malloc((ntev)*sizeof(struct prior));
		allTev = 0;
		for(ctev=0; ctev<ntev; ctev++){
			if(!allTev){
				printf("\n");
				if(ctev==0 && npop>2){
					do{
						printf("  choose the Tev%d prior [1 - uniform distribution; 2 - generalized gamma\n",ctev+1);
			    	    printf("  distribution; 3 - uniform distribution (for all tev) 4 - generalized\n");
			    	    printf("  gamma distribution (for all tev)]: ");
						scanf("%s",&str);
						pr_tev[ctev].type = atoi(str);
					}
					while(pr_tev[ctev].type!=1 && pr_tev[ctev].type!=2 && pr_tev[ctev].type!=3 && pr_tev[ctev].type!=4);
				}
				else{
					do{
						printf("  choose the Tev%d prior [1 - uniform distribution; 2 - generalized gamma\n",ctev+1);
			    	    printf("  distribution]: ");
						scanf("%s",&str);
						d1 = atoi(str);
		           		pr_tev[ctev].type= d1;
					}
					while(pr_tev[ctev].type!=1 && pr_tev[ctev].type!=2);
				}
				if(pr_tev[ctev].type==3||pr_tev[ctev].type==4){
					allTev=1;
				}
				if(pr_tev[ctev].type==1||pr_tev[ctev].type==3){
					pr_tev[ctev].p = (double *)malloc(2*sizeof(double));
					do{
						printf("  choose lower boundary of the prior (>=0): ");
						scanf("%s",&str);
						pr_tev[ctev].p[0] = atof(str);
						printf("  choose upper boundary of the prior (>0): ");
						scanf("%s",&str);
						pr_tev[ctev].p[1] = atof(str);
					}
					while(pr_tev[ctev].p[0]>pr_tev[ctev].p[1] || pr_tev[ctev].p[0]<0 || pr_tev[ctev].p[1]<=0);
				}
				else if(pr_tev[ctev].type==2||pr_tev[ctev].type==4){
					pr_tev[ctev].p = (double *)malloc(4*sizeof(double));
					do{
						printf("  choose location parameter (>=0): ");
						scanf("%s",&str);
						pr_tev[ctev].p[0] = atof(str);
					}
					while(pr_tev[ctev].p[0]<0);
					do{
						printf("  choose scale parameter (>0): ");
						scanf("%s",&str);
						pr_tev[ctev].p[1] = atof(str);
					}
					while(pr_tev[ctev].p[1]<=0);
					do{
						printf("  choose shape parameter (>0): ");
						scanf("%s",&str);
						pr_tev[ctev].p[2] = atof(str);
					}
					while(pr_tev[ctev].p[2]<=0);
					do{
						printf("  choose 2nd shape parameter (>0): ");
						scanf("%s",&str);
						pr_tev[ctev].p[3] = atof(str);
					}
					while(pr_tev[ctev].p[3]<=0);
				}
			}
		}
		/*migration priors*/
		pr_mig = (struct prior *)malloc((npop*2-2)*sizeof(struct prior));
		for(cpop=0; cpop<npop*2-2; cpop++){
			if(cpop<npop){
				printf("\n");
				do{
					printf("  choose the mig%d rate prior [0 - zero migration; 1 - uniform\n",cpop+1);
			        printf("  distribution; 2 - generalized gamma distribution; 3 - uniform\n");
			        printf("  distribution (on number of migrations); 4 - generalized gamma\n");
			        printf("  distribution(on number of migrations)]: ");
					scanf("%s",&str);
					pr_mig[cpop].type = atoi(str);
				}
				while(pr_mig[cpop].type!=0 && pr_mig[cpop].type!=1 && pr_mig[cpop].type!=2 && pr_mig[cpop].type!=3 && pr_mig[cpop].type!=4 && pr_mig[cpop].type!=5);
			}
			else{
				printf("\n");
				do{
					printf("  choose the migA%d rate prior [0 - zero migration; 1 - uniform\n",cpop+1-npop);
			        printf("  distribution; 2 - generalized gamma distribution; 3 - uniform\n");
			        printf("  distribution (on number of migrations); 4 - generalized gamma\n");
			        printf("  distribution(on number of migrations)]: ");
					scanf("%s",&str);
					pr_mig[cpop].type = atoi(str);
				}
				while(pr_mig[cpop].type!=0 && pr_mig[cpop].type!=1 && pr_mig[cpop].type!=2 && pr_mig[cpop].type!=3 && pr_mig[cpop].type!=4 && pr_mig[cpop].type!=5);
			}
	
			if(pr_mig[cpop].type==1||pr_mig[cpop].type==3){
				pr_mig[cpop].p = (double *)malloc(2*sizeof(double));
				do{
					printf("  choose lower boundary of the prior (>=0): ");
					scanf("%s",&str);
					pr_mig[cpop].p[0] = atof(str);
					printf("  choose upper boundary of the prior (>0): ");
					scanf("%s",&str);
					pr_mig[cpop].p[1] = atof(str);
				}
				while(pr_mig[cpop].p[0]>pr_mig[cpop].p[1] || pr_mig[cpop].p[0]<0 || pr_mig[cpop].p[1]<=0);
			}
			else if(pr_mig[cpop].type==2||pr_mig[cpop].type==4){
				pr_mig[cpop].p = (double *)malloc(4*sizeof(double));
				do{
					printf("  choose location parameter (=>0): ");
					scanf("%s",&str);
					pr_mig[cpop].p[0] = atof(str);
				}
				while(pr_mig[cpop].p[0]<0);
				do{
					printf("  choose scale parameter (>0): ");
					scanf("%s",&str);
					pr_mig[cpop].p[1] = atof(str);
				}
				while(pr_mig[cpop].p[1]<=0);
				do{
					printf("  choose shape parameter (>0): ");
					scanf("%s",&str);
					pr_mig[cpop].p[2] = atof(str);
				}
				while(pr_mig[cpop].p[2]<=0);
				do{
					printf("  choose 2nd shape parameter (>0): ");
					scanf("%s",&str);
					pr_mig[cpop].p[3] = atof(str);
				}
				while(pr_mig[cpop].p[3]<=0);
			}
		}
	}
	/*mutation priors*/
	printf("\n");
	if(nloc==1){
		do{
			printf("  choose STR mutation rate prior [0 - no mutation; 1 - lognormal distribution];\n");
			printf("  2 - normal distribution]: ");
			scanf("%s",&str);
			pr_mutSTR.type = atoi(str);
		}
		while(pr_mutSTR.type!=0 && pr_mutSTR.type!=1 && pr_mutSTR.type!=2);
	}
	else{
		do{
			printf("  choose STR mutation rate hiperprior [0 - no mutation; 1 - lognormal\n");
			printf("  distribution; 2 - normal distribution]: ");
			scanf("%s",&str);
			pr_mutSTR.type = atoi(str);
		}
		while(pr_mutSTR.type!=0 && pr_mutSTR.type!=1 && pr_mutSTR.type!=2);
	}
	if(pr_mutSTR.type==1){
		pr_mutSTR.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[0] = atof(str);
			}
			while(pr_mutSTR.p[0]>=0);

			pr_mutSTR.p[1] = 0;

			do{
				printf("  choose the stdev of the prior [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[2] = atof(str);
			}
			while(pr_mutSTR.p[2]<0);

			pr_mutSTR.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[0] = atof(str);
			}
			while(pr_mutSTR.p[0]>=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_mutSTR.p[1] = atof(str);
			}
			while(pr_mutSTR.p[1]<0);
			do{
				printf("  choose the mean of all stdevs [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[2] = atof(str);
			}
			while(pr_mutSTR.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_mutSTR.p[3] = atof(str);
			}
			while(pr_mutSTR.p[3]<0);
		}
	}
	if(pr_mutSTR.type==2){
		pr_mutSTR.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior (>0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[0] = atof(str);
			}
			while(pr_mutSTR.p[0]<=0);

			pr_mutSTR.p[1] = 0;

			do{
				printf("  choose the stdev of the prior (>=0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[2] = atof(str);
			}
			while(pr_mutSTR.p[2]<0);

			pr_mutSTR.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means (>0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[0] = atof(str);
			}
			while(pr_mutSTR.p[0]<=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_mutSTR.p[1] = atof(str);
			}
			while(pr_mutSTR.p[1]<0);
			do{
				printf("  choose the mean of all stdevs (>=0)]: ");
				scanf("%s",&str);
				pr_mutSTR.p[2] = atof(str);
			}
			while(pr_mutSTR.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_mutSTR.p[3] = atof(str);
			}
			while(pr_mutSTR.p[3]<0);
		}
	}
	printf("\n");
	if(nloc==1){
		do{
			printf("  choose SNP mutation rate prior [0 - no mutation; 1 - lognormal distribution];\n");
			printf("  2 - normal distribution]: ");
			scanf("%s",&str);
			pr_mutSNP.type = atoi(str);
		}
		while(pr_mutSNP.type!=0 && pr_mutSNP.type!=1 && pr_mutSNP.type!=2);
	}
	else{
		do{
			printf("  choose SNP mutation rate hiperprior [0 - no mutation; 1 - lognormal\n");
			printf("  distribution; 2 - normal distribution]: ");
			scanf("%s",&str);
			pr_mutSNP.type = atoi(str);
		}
		while(pr_mutSNP.type!=0 && pr_mutSNP.type!=1 && pr_mutSNP.type!=2);
	}
	if(pr_mutSNP.type==1){
		pr_mutSNP.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[0] = atof(str);
			}
			while(pr_mutSNP.p[0]>=0);

			pr_mutSNP.p[1] = 0;

			do{
				printf("  choose the stdev of the prior [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[2] = atof(str);
			}
			while(pr_mutSNP.p[2]<0);

			pr_mutSNP.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[0] = atof(str);
			}
			while(pr_mutSNP.p[0]>=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_mutSNP.p[1] = atof(str);
			}
			while(pr_mutSNP.p[1]<0);
			do{
				printf("  choose the mean of all stdevs [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[2] = atof(str);
			}
			while(pr_mutSNP.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_mutSNP.p[3] = atof(str);
			}
			while(pr_mutSNP.p[3]<0);
		}
	}
	if(pr_mutSNP.type==2){
		pr_mutSNP.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior (>0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[0] = atof(str);
			}
			while(pr_mutSNP.p[0]<=0);

			pr_mutSNP.p[1] = 0;

			do{
				printf("  choose the stdev of the prior (>=0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[2] = atof(str);
			}
			while(pr_mutSNP.p[2]<0);

			pr_mutSNP.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means (>0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[0] = atof(str);
			}
			while(pr_mutSNP.p[0]<=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_mutSNP.p[1] = atof(str);
			}
			while(pr_mutSNP.p[1]<0);
			do{
				printf("  choose the mean of all stdevs (>=0)]: ");
				scanf("%s",&str);
				pr_mutSNP.p[2] = atof(str);
			}
			while(pr_mutSNP.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_mutSNP.p[3] = atof(str);
			}
			while(pr_mutSNP.p[3]<0);
		}
	}
	/*recombination priors*/
	printf("\n");
	if(nloc==1){
		do{
			printf("  choose STR recombination prior [0 - no recombination; 1 - lognormal\n");
			printf("  distribution]; 2 - normal distribution]: ");
			scanf("%s",&str);
			pr_recSTR.type = atoi(str);
		}
		while(pr_recSTR.type!=0 && pr_recSTR.type!=1 && pr_recSTR.type!=2);
	}
	else{
		do{
			printf("  choose STR recombination rate hiperprior [0 - no recombination; 1 - lognormal\n");
			printf("  distribution; 2 - normal distribution]: ");
			scanf("%s",&str);
			pr_recSTR.type = atoi(str);
		}
		while(pr_recSTR.type!=0 && pr_recSTR.type!=1 && pr_recSTR.type!=2);
	}
	if(pr_recSTR.type==1){
		pr_recSTR.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[0] = atof(str);
			}
			while(pr_recSTR.p[0]>=0);

			pr_recSTR.p[1] = 0;

			do{
				printf("  choose the stdev of the prior [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[2] = atof(str);
			}
			while(pr_recSTR.p[2]<0);

			pr_recSTR.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[0] = atof(str);
			}
			while(pr_recSTR.p[0]>=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_recSTR.p[1] = atof(str);
			}
			while(pr_recSTR.p[1]<0);
			do{
				printf("  choose the mean of all stdevs [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[2] = atof(str);
			}
			while(pr_recSTR.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_recSTR.p[3] = atof(str);
			}
			while(pr_recSTR.p[3]<0);
		}
	}
	if(pr_recSTR.type==2){
		pr_recSTR.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior (>0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[0] = atof(str);
			}
			while(pr_recSTR.p[0]<=0);

			pr_recSTR.p[1] = 0;

			do{
				printf("  choose the stdev of the prior (>=0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[2] = atof(str);
			}
			while(pr_recSTR.p[2]<0);

			pr_recSTR.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means (>0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[0] = atof(str);
			}
			while(pr_recSTR.p[0]<=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_recSTR.p[1] = atof(str);
			}
			while(pr_recSTR.p[1]<0);
			do{
				printf("  choose the mean of all stdevs (>=0)]: ");
				scanf("%s",&str);
				pr_recSTR.p[2] = atof(str);
			}
			while(pr_recSTR.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_recSTR.p[3] = atof(str);
			}
			while(pr_recSTR.p[3]<0);
		}
	}
	printf("\n");
	if(nloc==1){
		do{
			printf("  choose SNP recombination rate prior [0 - no recombination; 1 - lognormal\n");
			printf("  distribution]; 2 - normal distribution]: ");
			scanf("%s",&str);
			pr_recSNP.type = atoi(str);
		}
		while(pr_recSNP.type!=0 && pr_recSNP.type!=1 && pr_recSNP.type!=2);
	}
	else{
		do{
			printf("  choose SNP recombination rate hiperprior [0 - no recombination; 1 - lognormal\n");
			printf("  distribution; 2 - normal distribution]: ");
			scanf("%s",&str);
			pr_recSNP.type = atoi(str);
		}
		while(pr_recSNP.type!=0 && pr_recSNP.type!=1 && pr_recSNP.type!=2);
	}
	if(pr_recSNP.type==1){
		pr_recSNP.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[0] = atof(str);
			}
			while(pr_recSNP.p[0]>=0);

			pr_recSNP.p[1] = 0;

			do{
				printf("  choose the stdev of the prior [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[2] = atof(str);
			}
			while(pr_recSNP.p[2]<0);

			pr_recSNP.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means [lognormal scale (<0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[0] = atof(str);
			}
			while(pr_recSNP.p[0]>=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_recSNP.p[1] = atof(str);
			}
			while(pr_recSNP.p[1]<0);
			do{
				printf("  choose the mean of all stdevs [lognormal scale(>=0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[2] = atof(str);
			}
			while(pr_recSNP.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_recSNP.p[3] = atof(str);
			}
			while(pr_recSNP.p[3]<0);
		}
	}
	if(pr_recSNP.type==2){
		pr_recSNP.p = (double *)malloc(4*sizeof(double));

		if(nloc==1){
			do{
				printf("  choose the mean of the prior (>0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[0] = atof(str);
			}
			while(pr_recSNP.p[0]<=0);

			pr_recSNP.p[1] = 0;

			do{
				printf("  choose the stdev of the prior (>=0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[2] = atof(str);
			}
			while(pr_recSNP.p[2]<0);

			pr_recSNP.p[3] = 0;
		}
		else{
			do{
				printf("  choose the mean of all means (>0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[0] = atof(str);
			}
			while(pr_recSNP.p[0]<=0);
			do{
				printf("  choose the stdev of all means (>=0): ");
				scanf("%s",&str);
				pr_recSNP.p[1] = atof(str);
			}
			while(pr_recSNP.p[1]<0);
			do{
				printf("  choose the mean of all stdevs (>=0)]: ");
				scanf("%s",&str);
				pr_recSNP.p[2] = atof(str);
			}
			while(pr_recSNP.p[2]<0);
			do{
				printf("  choose the stdev of all stdevs (>=0): ");
				scanf("%s",&str);
				pr_recSNP.p[3] = atof(str);
			}
			while(pr_recSNP.p[3]<0);
		}
	}
	/*migration matrix*/
	if(npop>=3 && (pr_top.type==1 || pr_top.type==2 || pr_top.type==4 || pr_top.type==5)){
		printf("\n");
		do{
			printf("  choose to set or not migration weights [0 - do not use; 1 - use]: ");
			scanf("%s",&str);
			migw.type = atoi(str);
		}
		while(migw.type!=0 && migw.type!=1);
	
		if(migw.type==1){
			migw.m = (double ***)malloc(npop*sizeof(double**));
			for(cpop=0;cpop<npop;cpop++){
				migw.m[cpop] = (double **)malloc(ntev*sizeof(double*));
				for(ctev=0;ctev<ntev;ctev++)
					migw.m[cpop][ctev] = (double *)malloc(npop*sizeof(double));
			}
			printf("\n  choose the percentage of migrants:\n");
			for(cpop=0;cpop<npop;cpop++){
				for(ctev=0;ctev<ntev;ctev++){
					for(cpop2=0;cpop2<npop;cpop2++){
						if(cpop==cpop2)
							migw.m[cpop][ctev][cpop2] = 0;
						else{
							do{
								printf("    in pop %d from pop %d at time %d (0 <= per <= 1): ",cpop,cpop2,ctev);
								scanf("%s",&str);
								migw.m[cpop][ctev][cpop2] = atof(str);
							}
							while(migw.m[cpop][ctev][cpop2]<0 || migw.m[cpop][ctev][cpop2]>1);
						}
					}
				}
			}
		}
	}
	else
		migw.type=0;
	/*output filename*/
	printf("\n  insert name for prior file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);
	
	/*run makeprior()*/
	error = makeprior(output,niter,genet,npop,nloc,lplo,ltype,pr_top,pr_Ne,pr_tev,pr_mig,
					 pr_mutSTR,pr_mutSNP,pr_recSTR,pr_recSNP,migw);
	
	/*free stuff*/
	if(pr_top.type==2||pr_top.type==5){
		free(avlpop);
	}
	return error;

}  //ends of buildprior

int exists(int *avlpop,int pop,int size){
	int i,
		foundit;
	
	foundit = 0;	
	for(i=0; i<size && !foundit; i++){
		if(avlpop[i] == pop)
			foundit = 1;
	}
	
	return foundit;
		
} //end of exists

int *takeout(int *avlpop,int pop, int size){
	int i,			//iterator
		foundit,	//check if the pop was found
		*aftlist;	//updated list
	
	aftlist = (int *)malloc((size-1)*sizeof(int));
	
	foundit = 0;
	for(i=0; i<size; i++){
		if(foundit)
			aftlist[i-1] = avlpop[i];
		else{
			if(avlpop[i]==pop)
				foundit = 1;	
			else
				aftlist[i] = avlpop[i];
		}
	}
	
	//free stuff
	free(avlpop);
	
	return aftlist;
	
} //end of takeout

int buildsstats(){
	int cstat,					//iterator
		lsstats[MAXSSTATS];		//list of the summary statistics used
	char str[MAXCHAR],			//auxiliar to get the input from the user
		 output[MAXCHAR];		//output filename
		 	
	//summary statistics
	for(cstat=0; cstat<MAXSSTATS; cstat++)
		lsstats[cstat]=0;
	printf("  choose the summary statistics to be used (0 - absent; 1 - present):\n");
	for(cstat=0; cstat<MAXSSTATS; cstat++){
		do{
			switch(cstat){
			case 0:
			printf("\n  microsatellites data\n");
			printf("    heterozygosity: ");
			break;
			case 1:
			printf("    variance of alleles length: ");
			break;
			case 2:
			printf("    number of alleles: ");
			break;
			case 3:
			printf("    curtosis of alleles length: ");
			break;
			case 4:
			printf("    Shanon's index: ");
			break;
			case 5:
			printf("    Nm estimator based on H: ");
			break;
			case 6:
			printf("\n  sequence data\n");
			printf("    mean of pairwise differences: ");
			break;
			case 7:
			printf("    number of segregating sites: ");
			break;
			case 8:
			printf("    number of haplotypes: ");
			break;
			case 9:
			printf("    Shanon's index: ");
			break;
			case 10:
			printf("    mean of MFS: ");
			break;
			case 11:
			printf("    stdev of MFS: ");
			break;
			case 12:
			printf("    Nm estimator based on S: ");
			break;
			case 13:
			printf("    private segregating sites: ");
			break;
			case 14:
			printf("    S(1): ");
			break;
			}
			scanf("%s",&str);
			lsstats[cstat] = atoi(str);
		}
		while(lsstats[cstat]!=1 && lsstats[cstat]!=0);
	}
	/*output filename*/
	printf("\n  insert name for summary statistics set file (no extension)");
	printf("\n  (max %d character): ",MAXCHAR-1);
	scanf("%s",output);
	
	/*run makestats()*/
	return makestats(output,lsstats);

}  //ends of buildsstats

int createtab(){
	char str[MAXCHAR],		//auxiliar to get the input from the user
		 input[MAXCHAR],	//input file name
		 output[MAXCHAR];	//input file name

	//input filename
	printf("\n  choose sample file (.pop) (max %d character): ",MAXCHAR-1);
	scanf("%s",input);

	//output filename
	printf("\n  insert name for table file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	return createFreqTab(input,output);
	
} //end of createtab1

int createtab2(){
	char str[MAXCHAR],		//auxiliar to get the input from the user
		 input[MAXCHAR],	//input file name
		 output[MAXCHAR];	//input file name

	//input filename
	printf("\n  choose IMa sample file (max %d character): ",MAXCHAR-1);
	scanf("%s",input);

	//output filename
	printf("\n  insert name for table file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	return createFreqTab2(input,output);
	
} //end of createtab2

int createtab3(){
	char str[MAXCHAR],		//auxiliar to get the input from the user
		 input[MAXCHAR],	//input file name
		 output[MAXCHAR];	//input file name

	//input filename
	printf("\n  choose GenePop sample file (max %d character): ",MAXCHAR-1);
	scanf("%s",input);

	//output filename
	printf("\n  insert name for table file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	return createFreqTab3(input,output);
	
} //end of createtab3

int createtab4(){
	char str[MAXCHAR],		//auxiliar to get the input from the user
		 input[MAXCHAR],	//input file name
		 output[MAXCHAR];	//input file name

	//input filename
	printf("\n  choose Nexus file (.nex) (max %d character): ",MAXCHAR-1);
	scanf("%s",input);

	//output filename
	printf("\n  insert name for table file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	return createFreqTab4(input,output);
	
} //end of createtab4

int createpop(){
	char str[MAXCHAR],		//auxiliar to get the input from the user
		 input[MAXCHAR],	//input file name
		 output[MAXCHAR];	//input file name

	//input filename
	printf("\n  choose table file (.len) (max %d character): ",MAXCHAR-1);
	scanf("%s",input);

	//output filename
	printf("\n  insert name for sample file (no extension) (max %d character): ",MAXCHAR-1);
	scanf("%s",output);

	return makepop(input,output);
	
} //end of convertTab

int joinfiles(){
	int i,				//iterator
		error,			//saves the error type
		ninp;			//number of files to join
	char outp[MAXCHAR],	//output file name
		 str[MAXCHAR],	//auxiliar to get the input from the user
		 **linp;		//list of the files to join together
		 
	//number of files to join
	do{
		printf("  choose the number of files to join (>1): ");
		scanf("%s",&str);
		ninp = atoi(str);
	}
	while(ninp <=1);

	//path of the files to join
	linp = (char**)malloc(ninp*sizeof(char*));
	printf("\n  choose %d files to be joint (.dat)\n",ninp);
	for(i=0; i<ninp; i++){
		linp[i] = (char*)malloc(MAXCHAR*sizeof(char));
		printf("    file %d (max %d character): ",i+1,MAXCHAR-1);
		scanf("%s",linp[i]);
	}
	
	//output filename
	printf("\n  insert name for created file (.dat) (max %d character): ",MAXCHAR-1);
	scanf("%s",outp);
	
	//run joindata()
	error = joindata(ninp,linp,outp);
	
	//free stuff
	for(i=0; i<ninp; i++){
		free(linp[i]);
	}
	free(linp);
	
	return error;
	
} //end of joinfiles
