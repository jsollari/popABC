/*
	@author:	  joao lopes
	@workplace:   Reading University
	@date: 		  29th May 2009

	NBB - based on Mark's myutil.h	
*/
#ifndef _MYLIB_H_
#define _MYLIB_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI  3.141592653589793

/*
	It creates the random table used to create random numbers
*/
void opengfsr(char* path);

/*
	It saves random information back to INTFILE and closes it
*/
void closegfsr(char* path);

/*
	It prints the error to a ERRFILE and ends the program
	
	@param s - information about the type of error
*/
void printerr(char *s);
	
/*
	It returns 0 or 1
	
	@return the value 0 or 1 with equal probability
*/
int intrand(void);

/*
	It returns a random int in a uniform distribution between l(exclusive) and t(inclusive)
	
	@param l - smallest value of a uniform distribution
	@param t - gratest value of a uniform distribution
	@return a value from a uniform distribution
*/
int disrand(int l, int t);

/*
	It returns a random float variable uniformely distributed between 0.0,1.0
	
	@return a value from a uniform distribution between 0.0 and 1.0
*/
float gfsr4(void);

/*
	It returns a random double uniformely distributed between 0.0,1.0
	
	@return a value from a uniform distribution between 0.0 and 1.0
*/
double gfsr8(void);

/*
	It returns a random value from a rgamma distribution with parameters a and scale
	
	@param a 	 -  parameter of a gama distribution
	@param scale - parameter of a gama distribution
	@return a value of a gama distribution
*/
double rgamma(double a, double scale);

/*
	It returns a random int value from a poisson distribution with parameter xm
	
	@param - parameter of a poisson distribution
	@return a value of a poisson distribution
*/
int poidev(float xm);

/*
	It returns a random value from a lgamma distribution with parameter arg

	@param - parameter of a gamma distribution
	@return a value of a gamma distribution
 	
 */
double lgamma(double arg);

/*
	Exponential distribution
	
	@return a value of a exponential destribution normalized
*/
double expdev(void);

/*
	It returns a random value of a standard normal distribution
	
	@return value from a standard normal distribution
*/
double norm8(void);

/*
	This function gets a serie of values and gets some statistics about them
	
	@param x   - serie of values to be analised
	@param n   - number of values of the series
	@return x1  - mean
	@return x2  - variance 
	@return x3  - skew 
	@return x4  - curtosis 
	@return min - minimum value of the serie of values analised
	@return max - maximum value of the serie of values analised 
*/
void mom(double x[],int n,double *x1,double *x2,double *x3,double *x4, double *min,double *max);

/*
	Get a sorted list from list x and store the original position of the sorted elements in indx
	
	@param char dir   - ('a' or 'A') or anything else
	@param int n 	  - number of elements
	@param double *x  - list from wich the elements will be sorted
	@param int *indx  - list with the original positions of the sorted elements
*/
void dsorti(char dir,int n,double * x,int *indx);

/*
	Sort the list x and put it sorted in indx
	
	@param char dir  - ('a' or 'A') or 
	@param int n 	 - number of elements
	@param int *x 	 - list to be sorted
	@param int *indx - list sorted
*/
void isorti(char dir, int n, int * x, int *indx);

/*
	This fuction returns a value from a beta distribution with parameters a[0] and a[1]
	
	@param a - list with the parameters used in the simulation
	@return a value from a beta distribution
*/
double betasim(double a[]);

/*
	This function returns a value from a gama distribution with parameters a, b, c and k 
	
	@param a - parameter of a gama distribution
	@param b - parameter of a gama distribution
	@param c - parameter of a gama distribution
	@param k - parameter of a gama distribution
	@return a value from a gama distribution 
*/
double rgengamma(double a, double b, double c, double k);

/*
	This function calculates combinations 
	Defined in summStats.c
	
	@param n - total number of cases
	@param k - size of combinations
	@return combinations of n choose k
*/
int combinations(int n,int k);

/*
	This function calculates permutations
	Defined in summStats.c
	
	@param n - number of permutations
	@return the permutations of n
*/
int permutation(int n);

/*
	Replace malloc and realloc with just one function
	
	@param p - pntr to set the memory to
	@param n - size of the memory we want to allocate
	@return malloc(n) or realloc(p,n)
*/
void *myAlloc(void *p, size_t n);

/*
	This function gets a ordinal number and returns the correspondent prefix

	@param cpop - ordinal number
	@return the number correspondent prefix
*/
char *intToPrefix(int n);

/*
	This function gets the ordinal number in the alphabet of a letter an returns the correspondent char

	@param cpop - ordinal number in the alphabet of a letter
	@return the number correspondent char
*/
char intToChar(int cpop);

/*
	Check if a character is a change-line character (only used in sequence data analysis)
	
	@param c - given character
	@return a boolean (true if c is a change line symbol, false if it's not) 
*/
int isendline(char c);

/*
	Shell-sort - sort a list of type double variables

	@param A   	- list to be sorted
	@param size - size of the list to be sorted
*/
void shell_sort_double(double A[],int size);

/*
	Shell-sort - sort a list of type int variables

	@param A   	- list to be sorted
	@param size - size of the list to be sorted
*/
void shell_sort_int(int A[],int size);

#endif //_MYLIB_H_
