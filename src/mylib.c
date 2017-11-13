/*
	@author:	  joao lopes
	@workplace:   Reading University
	@date: 		  29th May 2009

	NBB - based on Mark's myutil.h	
*/

#include "mylib.h"

#define repeat for(;;)

int jindic,
	rand_table[98];

void opengfsr(char*path){
	static int rand[] =	{
		-899909286,-1806404725,1582468931,-1678321118,97185498,1145178031,-814937526,256730278,-1691631069,
		1718777475,-453536501,1857826690,893882730,57758276,834697662,-803668054,-591757859,-1657313141,
		88619219,1998479650,-701262517,-1471231107,-1431743148,2132644097,1359033410,-657319434,284575260,
		720437758,979592561,1796725011,-1118945874,899985249,-1324873621,259372545,2035744715,-8273512,
		275256135,289288160,843894067,2061106821,371848653,-1679272953,1522134473,-989926898,-1648586077,
		79513443,-1408033903,811767484,-2087304794,1636025030,-1529047970,-1616471647,358603871,-553830325,
		-764029463,277775098,-529497046,84613338,-1122707738,1574486947,2118784361,-895839015,1689888366,
		1173044288,562273445,-1594728793,2140441606,1892789422,-2051052488,95435543,-1707679106,-546465673,
		-150225236,970702808,1691166307,-1814393727,586361783,-914226192,-553100016,-996795202,-522545274,
		-1975984216,1038272682,427900780,559898624,-1952252899,-1770534366,-1930252922,674961545,814541338,
		1834206839,-1117614972,-1650957533,-1723516185,1757781308,870369684,1911065966,-2130395123,91
	};
	FILE *rt;
	char* aux; 
	int j,
		randInt;

	aux = (char*)malloc((strlen(path)+8)*sizeof(char));
	strcpy(aux,path);
	aux = strcat(aux,"INTFILE");
	rt = fopen(aux,"r");
	if(rt==NULL){
		for(j=0;j<98;++j){
			rand_table[j] = rand[j];
		}
		jindic = rand[j];
	}
	else{
		for(j=0;j<98;++j){
			fscanf(rt,"%d",&randInt);
			if(randInt == 0)
				printerr("0 in INTFILE");
			rand_table[j] = randInt;
		}
		fscanf(rt,"%d",&jindic);
		fclose(rt);
	}
	free(aux);

}	//end of opengfsr

void closegfsr(char*path){
	FILE 	*rt;
	char* aux; 
	int 	j;

	aux = (char*)malloc((strlen(path)+8)*sizeof(char));
	strcpy(aux,path);
	aux = strcat(aux,"INTFILE");
	rt = fopen(aux,"w");
	if(rt==NULL){
		printerr("INTFILE couldn't be created\n");
	}
	for(j=0;j<98;++j)
		fprintf(rt,"%d\n",rand_table[j]);
	fprintf(rt,"%d\n",jindic);
	free(aux);
	fclose(rt);
	
}	//end of closegfsr

void printerr(char *s)
{
	FILE *erf;
	erf = fopen("ERRFILE","w");
	printf("error:: %s\n",s);
	fprintf(erf,"error:: %s\n",s);
	fclose(erf);
	exit(1);
}

int intrand(void)
{
      int k;
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
      return(rand_table[jindic]);
}

int disrand(int l,int t)
{
      int k;
      if(t<l){
      	printerr("error in disrand\n");
      }
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return((unsigned)rand_table[jindic]%(t-l+1)+l);
}

float gfsr4(void)
{
      int k;
      ++jindic;
      if(jindic>97)
      	jindic=0;
      k = jindic+27;
      if(k>97)
      	k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return(((unsigned)rand_table[jindic] + 1.0)/4294967298.0);
}

double gfsr8(void)
{
    int k;
    ++jindic;
    if(jindic>97)
		jindic=0;
    k = jindic+27;
    if(k>97)
		k=k-98;
    rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return(((unsigned)rand_table[jindic] + 1.0)/4294967298.0);
}

double rgamma(double a, double scale)
{
	static double a1 = 0.3333333;
	static double a2 = -0.250003;
	static double a3 = 0.2000062;
	static double a4 = -0.1662921;
	static double a5 = 0.1423657;
	static double a6 = -0.1367177;
	static double a7 = 0.1233795;
	static double e1 = 1.0;
	static double e2 = 0.4999897;
	static double e3 = 0.166829;
	static double e4 = 0.0407753;
	static double e5 = 0.010293;
	static double q1 = 0.04166669;
	static double q2 = 0.02083148;
	static double q3 = 0.00801191;
	static double q4 = 0.00144121;
	static double q5 = -7.388e-5;
	static double q6 = 2.4511e-4;
	static double q7 = 2.424e-4;
	static double sqrt32 = 5.656854;

	static double aa = 0.;
	static double aaa = 0.;

	static double b, c, d, e, p, q, r, s, t, u, v, w, x;
	static double q0, s2, si;
	double ret_val;

	if (a < 1.0)
	{
		/* alternate method for parameters a below 1 */
		/* 0.36787944117144232159 = exp(-1) */
		aa = 0.0;
		b = 1.0 + 0.36787944117144232159 * a;
		repeat 
		{
			p = b * gfsr8();
			if (p >= 1.0) 
			{
				ret_val = -log((b - p) / a);
				if (expdev() >= (1.0 - a) * log(ret_val))
					break;
			}
			else
			{
				ret_val = exp(log(p) / a);
				if (expdev() >= ret_val)
					break;
			}
		}
	return scale * ret_val;
	}
	/* Step 1: Recalculations of s2, s, d if a has changed */
	if (a != aa) 
	{
		aa = a;
		s2 = a - 0.5;
		s = sqrt(s2);
		d = sqrt32 - s * 12.0;
	}
	
	/* Step 2: t = standard normal deviate, */
	/* x = (s,1/2)-normal deviate. */
	/* immediate acceptance (i) */
	t = norm8();
	x = s + 0.5 * t;
	ret_val = x * x;
	if (t >= 0.0)
		return scale * ret_val;

	/* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
	u = gfsr8();
	if (d * u <= t * t * t) 
	{
		return scale * ret_val;
	}

	/* Step 4: recalculations of q0, b, si, c if necessary */
	if (a != aaa) 
	{
		aaa = a;
		r = 1.0 / a;
		q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r + q2) * r + q1) * r;

		/* Approximation depending on size of parameter a */
		/* The constants in the expressions for b, si and */
		/* c were established by numerical experiments */
		if (a <= 3.686) 
		{
			b = 0.463 + s + 0.178 * s2;
			si = 1.235;
			c = 0.195 / s - 0.079 + 0.16 * s;
		}
		else if (a <= 13.022) 
		{
			b = 1.654 + 0.0076 * s2;
			si = 1.68 / s + 0.275;
			c = 0.062 / s + 0.024;
		}
		else 
		{
			b = 1.77;
			si = 0.75;
			c = 0.1515 / s;
		}
	}
	
	/* Step 5: no quotient test if x not positive */
	if (x > 0.0)
	{

		/* Step 6: calculation of v and quotient q */
		v = t / (s + s);
		if (fabs(v) <= 0.25)
			q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
		else
			q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);

		/* Step 7: quotient acceptance (q) */
		if (log(1.0 - u) <= q)
			return scale * ret_val;
	}

	/* Step 8: e = standard exponential deviate */
	/* u= 0,1 -uniform deviate */
	/* t=(b,si)-double exponential (laplace) sample */
	repeat
	{
		e = expdev();
		u = gfsr8();
		u = u + u - 1.0;
		if (u < 0.0)
			t = b - si * e;
		else
			t = b + si * e;
		/* Step  9:  rejection if t < tau(1) = -0.71874483771719 */
		if (t >= -0.71874483771719) 
		{
			/* Step 10:  calculation of v and quotient q */
			v = t / (s + s);
			if (fabs(v) <= 0.25)
				q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
			else
				q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);

			/* Step 11:  hat acceptance (h) */
			/* (if q not positive go to step 8) */
			if (q > 0.0)
			{
				if (q <= 0.5)
					w = ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q;
				else
					w = exp(q) - 1.0;

				/* if t is rejected */
				/* sample again at step 8 */
				if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
					break;
			}
		}
	}
	x = s + 0.5 * t;
	return scale * x * x;
}

int poidev(float xm)
{
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			em += 1.0;
			t *= gfsr4();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-lgamma(xm+1.0);
		}
		do {
			do {
				y=tan(PI*gfsr4());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-lgamma(em+1.0)-g);
		} while (gfsr4() > t);
	}
	return (int) em+0.5;
	
}	//end of poidev

/* Inserted from myutil.c::lgamma() ***************************************************************************/
/* R0.16::lgamma.c modificated by MAB (20.3.97) */

#undef M_PI
#define M_PI		3.141592653589793238462643383276

/*
	Auxiliar to lgamma()

	@param double
	@return double
*/
static double posarg(double);

/*
	Auxiliar to lgamma()

	@param double
	@return double
*/
static double negarg(double);

/* 
	Equation 6.1.41 Abramowitz and Stegun, auxiliar to lgamma()
	(See also ACM algorithm 291)

	@param double
	@return double
*/
static double asform(double);

int signgamR1 = 0;
static double hl2pi = 0.9189385332046727417803297; //log(2*pi)/2 and pi
static double xpi = M_PI; 						   //pi
static int M = 6, N = 8; 						   //coefficients from Cheney and Hart
static double p1[] =
{
	0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2,
},
			  p2[] =
{
	-.42353689509744089647e5,
	-.20886861789269887364e5,
	-.87627102978521489560e4,
	-.20085274013072791214e4,
	-.43933044406002567613e3,
	-.50108693752970953015e2,
	-.67449507245925289918e1,
	0.0,
},
			  q2[] =
{
	-.42353689509744090010e5,
	-.29803853309256649932e4,
	0.99403074150827709015e4,
	-.15286072737795220248e4,
	-.49902852662143904834e3,
	0.18949823415702801641e3,
	-.23081551524580124562e2,
	0.10000000000000000000e1,
};

double lgamma(double arg)
{
	signgamR1 = 1.0;
	if (arg <= 0.0)
		return (negarg(arg));
	if (arg > 8.0)
		return (asform(arg));
	return (log(posarg(arg)));
}

static double asform(double arg)
{
	double log();
	double n, argsq;
	int i;

	argsq = 1. / (arg * arg);
	for (n = 0, i = M - 1; i >= 0; i--) {
		n = n * argsq + p1[i];
	}
	return ((arg - .5) * log(arg) - arg + hl2pi + n / arg);
}

static double negarg(double arg)
{
	double temp;
	double log(), sin(), posarg();

	arg = -arg;
	temp = sin(xpi * arg);
	if (temp == 0.0)
		printerr("negarg: temp == 0.0");
	if (temp < 0.0)
		temp = -temp;
	else
		signgamR1 = -1;
	return (-log(arg * posarg(arg) * temp / xpi));
}

static double posarg(double arg)
{
	double n, d, s;
	register i;

	if (arg < 2.0)
		return (posarg(arg + 1.0) / arg);
	if (arg > 3.0)
		return ((arg - 1.0) * posarg(arg - 1.0));

	s = arg - 2.;
	for (n = 0, d = 0, i = N - 1; i >= 0; i--) {
		n = n * s + p2[i];
		d = d * s + q2[i];
	}
	return (n / d);
}
/* end of lgamma insertion ************************************************************************************/

double expdev(void)
{

      int k;
      ++jindic;
      if(jindic>97)jindic=0;
      k = jindic+27;
      if(k>97)k=k-98;
      rand_table[jindic] = rand_table[jindic]^rand_table[k];
	return(-log(((unsigned)rand_table[jindic] + 1.0)/4294967298.0));

}	//end of expdev

double norm8(void)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if(iset == 0)
	{
		do {
			v1=2.0*gfsr8()-1.0;
			v2=2.0*gfsr8()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} 
	else
	{
		iset=0;
		return gset;
	}
	
}	//end of norm8

void mom(double x[],int n,double *x1,double *x2,double *x3,double *x4,double *min,double *max)
{
	int i;
	double s1,		//1st momment
		   s2,		//2nd momment [(x.-x)^2]
		   s3,		//3rd momment [(x.-x)^3]
		   s4,		//4th momment [(x.-x)^4]
		   an,an1,dx,dx2,xi,var;

	s1 = x[0];
	s2 = 0.0;
	s3 = 0.0;
	s4 = 0.0;
	*min = s1;
	*max = s1;
	for(i=1;i<n;++i)
	{
		xi = x[i];
		an = i+1;
		an1 = i;
		dx = (xi-s1)/an;
		dx2 = dx*dx;
		s4 -= dx*(4.0*s3-dx*(6.0*s2+an1*(1.0+pow(an1,3.0))*dx2));
		s3 -= dx*(3.0*s2-an*an1*(an-2.0)*dx2);
		s2 += an*an1*dx2;
		s1 += dx;
		if(xi<*min)*min=xi;
		if(xi>*max)*max=xi;
	}
	*x1 = s1;											//getting the mean
	var = n>1 ? s2/(n-1) : 0.0;
	*x2 = sqrt(var);									//getting the variance
	*x3 = var>0.0 ? 1.0/(n-1)*s3/pow(var,1.5) : 0.0;	//getting the skew out of the 3rd mommment
	*x4 = var>0.0 ? 1.0/(n-1)*s4/pow(var,2.0)-3.0 : 0.0;//getting the curtosis out of the 4th momment

}	//end of mom

void dsorti(char dir,int n,double * x,int *indx)  /* This is adapted from R 0.16 */
{
	int i, j, h, asc,indtmp;
	double xtmp,*priv;

	priv = (double *)malloc(n*sizeof(double));
	for(j=0;j<n;++j)priv[j] = x[j];

	if(dir == 'a' || dir == 'A')asc = 1;
	else asc = 0;

	h = 1;
	for(j=0;j<n;++j)indx[j] = j;
	do {
		h = 3 * h + 1;
	}
	while (h <= n);

	do {
		h = h / 3;
		for (i = h; i < n; i++) {
			xtmp = priv[i];
			indtmp = indx[i];
			j = i;
			if(asc){
				while (priv[j - h] > xtmp) {
					priv[j] = priv[j - h];
					indx[j] = indx[j-h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}
			else{
				while (priv[j - h] < xtmp) {
					priv[j] = priv[j - h];
					indx[j] = indx[j-h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}

		end:	priv[j] = xtmp;indx[j] = indtmp;;
		}
	} while (h != 1);
	free(priv);
	
}	//end of dsorti

void isorti(char dir, int n, int * x, int *indx){
	int i, j, h,
		asc,
		indtmp,
		xtmp,
		*priv;		 

	priv = (int *)malloc(n*sizeof(int));
	for(j=0 ; j<n ; ++j)
		priv[j] = x[j];
	if(dir == 'a' || dir == 'A')
		asc = 1;
	else
		asc = 0;
	for(j=0;j<n;++j)
		indx[j] = j;
	h = 1;
	do {
		h = 3 * h + 1;
	}
	while (h <= n);

	do {
		h = h / 3;
		for (i = h; i < n; i++){
			xtmp = priv[i];
			indtmp = indx[i];
			j = i;
			if(asc){
				while (priv[j - h] > xtmp){
					priv[j] = priv[j - h];
					indx[j] = indx[j-h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}
			else{
				while (priv[j - h] < xtmp){
					priv[j] = priv[j - h];
					indx[j] = indx[j-h];
					j = j - h;
					if (j < h)
						goto end;
				}
			}
		end:	priv[j] = xtmp;indx[j] = indtmp;
		}
	} 
	while (h != 1);
	free(priv);
	
} //end isorti

double betasim(double a[])
{
    double x1,x2;
    x1 = rgamma(a[0],1.0);
    x2 = rgamma(a[1],1.0);
    return x1/(x1+x2);

} //end of betasim

double rgengamma(double a, double b, double c, double k)
{
	return(pow(rgamma(c,1),1.0/k)*b + a);

} //end of rgengamma

int combinations(int n,int k)
{
	return permutation(n)/(permutation(k)*permutation(n-k));

} //end of combinations

int permutation(int n)
{
	int i,
		result = 1;
	
	for(i=n ; i>0 ; i--)
	{
		result *= i;
	}
	return result;
	
} //end of permutation

void *myAlloc(void *p, size_t n)
{
	if (n == 0)
	{
		n = 1;
		free (p);
		p = NULL;
	}
	if (p == NULL)
		return malloc (n);
	return realloc (p, n);
	
} //end of myAlloc

char *intToPrefix(int n)
{
	switch(n)
	{
		break;
		case 0:
			return("no-");		
		case 1:
			return("mono");
        break;
        case 2:
			return("di");
        break;
        case 3:
			return("tri");
        break;
        case 4:
			return("tetra");
        break;
        case 5:
			return("penta");
        break;
        case 6:
			return("hexa");
        break;
        case 7:
			return("hepta");
        break;
        case 8:
			return("octa");
        break;
        case 9:
			return("ennea");
        break;
        case 10:
			return("hexa");
        break;
        case 11:
			return("hendeca");
        break;
        case 12:
			return("dodeca");
        break;
        case 13:
			return("trideca");
        break;
        case 14:
			return("tetradeca");
        break;
        case 15:
			return("pentadeca");
        break;
        case 16:
			return("hexadeca");
        break;
        case 17:
			return("heptadeca");
        break;
        case 18:
			return("octadeca");
        break;
        case 19:
			return("enneadeca");
        break;
        case 20:
			return("icosa");
        break;
        default:
        	return ("poly");
	}
	
}	//end of intToPrefix

char intToChar(int cpop)
{
	switch(cpop)
	{
		case 0:
			return('A');
		break;		
		case 1:
			return('B');
        break;
        case 2:
			return('C');
        break;
        case 3:
			return('D');
        break;
        case 4:
			return('E');
        break;
        case 5:
			return('F');
        break;
        case 6:
			return('G');
        break;
        case 7:
			return('H');
        break;
        case 8:
			return('I');
        break;
        case 9:
			return('J');
        break;
        case 10:
			return('K');
        break;
        case 11:
			return('L');
        break;
        case 12:
			return('M');
        break;
        case 13:
			return('N');
        break;
        case 14:
			return('O');
        break;
        case 15:
			return('P');
        break;
        case 16:
			return('Q');
        break;
        case 17:
			return('R');
        break;
        case 18:
			return('S');
        break;
        case 19:
			return('T');
        break;
        case 20:
			return('U');
        break;
        case 21:
			return('V');
        break;
        case 22:
			return('W');
        break;
        case 23:
			return('X');
        break;
        case 24:
			return('Y');
        break;
        case 25:
			return('Z');
        break;
        default:
        	return ('\0');
	}
	
}	//end of intToChar

int isendline(char c)
{
	if(c=='\f')			//form feed
		return 1;
	if(c=='\n')			//new line
		return 1;
	if(c=='\r')			//carriage return
		return 1;
	if(c=='\v')			//vertical tab
		return 1;
	return 0;
	
}	//end of isendline

void shell_sort_double(double A[], int size)
{
    int i, j, incrmnt;
    double temp;

    incrmnt = size/2;
    while (incrmnt > 0)
    {
        for (i = incrmnt; i < size; i++)
        {
            j = i;
            temp = A[i];
            while ((j >= incrmnt) && (A[j - incrmnt] > temp))
            {
                A[j] = A[(j -= incrmnt)];
            }
            A[j] = temp;
        }
        incrmnt /= 2;
    }
    
} //end of shell_sort

void shell_sort_int(int A[], int size)
{
    int i, j, incrmnt,temp;

    incrmnt = size/2;
    while (incrmnt > 0)
    {
        for (i = incrmnt; i < size; i++)
        {
            j = i;
            temp = A[i];
            while ((j >= incrmnt) && (A[j - incrmnt] > temp))
            {
                A[j] = A[(j -= incrmnt)];
            }
            A[j] = temp;
        }
        incrmnt /= 2;
    }
    
} //end of shell_sort
