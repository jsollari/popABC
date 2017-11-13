/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 1th May 2009

	NBB - based on Mark’s turnit
	
	@update: - the INTFILE is used are on the same folder as turnit.exe
*/

#include "turnit.h"

/*
	Rewrite the random numbers table

	@arg - number of iterations to run
*/
int main(int argc, char *argv[])
{
	int iter,i,		//iterator
		cut,		//auxiliar to help on get the path from where turnit.exe is running
		niter,		//number of iterator
		homesize;	//size of the path from where the program is being run
	char *home;		//path where the program is being run

	//check error in arguments
	if(argc != 2)
		printerr("needs niter");
	for(iter=0 ; iter<strlen(argv[1]) ; iter++)
		if(!(isdigit(argv[1][iter])))
			printerr("argument is not a number");
	
	//get the path from where the progame is being run
	homesize = strlen(argv[0])+ 5;
	home = (char *)malloc(homesize*sizeof(char));
	strcpy(home,argv[0]);
	for(cut=-1,i=homesize-1 ; home[i]!='/' && home[i]!='\\' && i >= 0 ; cut++,i--)
	{
		home[i]='\0';
	}
	home = realloc(home,(homesize-cut)*sizeof(char));
	
	niter = atoi(argv[1]);

	//randomization
	printf("\n\n  start of randomization of INTFILE...\n");
	opengfsr(home);
	for(iter=0 ; iter < niter ; iter++)
	{
		intrand();
	}
	closegfsr(home);
	printf("\n  ...INTFILE randomized\n\n");
	
	//free stuff
	free(home);

} //end of main
