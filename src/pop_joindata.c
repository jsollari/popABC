/*
	@author:	joao lopes
	@workplace: Reading University
	@date: 		13th November 2008
*/
#include "interface.h"

int joindata(int ninp, char *linp[],char *outp){
	int i;			//iterator
	char c1;
	FILE *out,		//pntr to the output file
		 *in;		//pntr to the output files

	out = fopen(outp,"w");					
	if(out == NULL)
		return 1; //can't create output file

	printf("\n");
	for(i=0 ; i<ninp ; i++){
		in = fopen(linp[i],"r");					
		printf(">copying %s\n",linp[i]);
		if(in == NULL)
			return 2; //error::no input file
		while(1){
			c1=getc(in);

			//STOP condition
			if(c1 == EOF)
				break;

			putc(c1,out);
		}
		close(in);

	}

	//close file
	fclose(out);
	
	return 0;

} //end of main
