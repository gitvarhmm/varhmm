//Copyright 2016 Thomas M. Poulsen

#include <cstdlib>
#include <string.h> 

#include <math.h>
#include <time.h>

#include <algorithm>
#include <vector>

#include <limits>
#include <fstream>
#include <iostream>

#include <sstream>

using namespace std;


struct vfwdbackpair
{
	double ptotalfwd;
	double ptotalback;
	double vmax;
	double vmaxratio;
	double fratio;
	string vouts1;
	string vouts2;
	int matchlength;
	int nrepeats;
	int ninsert;
};


struct vfwdback
{
	double ptotalfwd;
	double ptotalback;
};



//const double DNEGINF = -(numeric_limits<double>::infinity());
const int AMIN = -(int)(numeric_limits<double>::infinity()); //-100000;

const int NOSTATETYPES = 6;

const int BEGINSTATE = 10;
const int ENDSTATE = 11;

const int MATCHSTATE = 1000;
const int OUTPUTSTATEX = 2000;
const int OUTPUTSTATEY = 3000;	
const int MATCHSTATEPH = 1001;
const int MATCHSTATEHO = 1002;
const int OUTPUTSTATEXPH = 2001;
const int OUTPUTSTATEXHO = 2002;
const int OUTPUTSTATEYPH = 3001;	
const int OUTPUTSTATEYHO = 3002;	
const int SILENTSTATE = 4000;	

const int RX1STATE = 5000;	
const int RY1STATE = 6000;	

const int RX2STATE = 7000;	
const int RY2STATE = 8000;	

const double PMIN = 0.0000001;

const double  MAXMISMAP              = 0.001;

const int     MINTRAINLENGTH         = 150;
const string sMINTRAINLENGTHCMD      = "l";
const string sMINTRAINLENGTHDESC     = "  Only use alignments with this minimum length              ";

const int     MAXTRAINITERATIONS     = 5;
const string sMAXTRAINITERATIONSCMD  = "i";
const string sMAXTRAINITERATIONSDESC = "  Maximum number of training iterations (if -1, no maximum) ";

const double  TRAINSTOPTHRESH        = 0.8;
const string sTRAINSTOPTHRESHCMD     = "d";
const string sTRAINSTOPTHRESHDESC    = "  Stop training if total prob. of all align. improves less  ";
const string sTRAINSTOPTHRESHDESC2   = "    than |log(d)| (if -1, then d->inifinity)";

const double VERSION = 1.01;

#define hmmv(i,j,k) hmmv[(k+j*nstates+i*ndatatotal2*nstates)]
#define hmmptrs(i,j,k) hmmptrs[(k+j*nstates+i*ndatatotal2*nstates)]
#define hmmx(i,j,k) hmmx[(k+j*nstates+i*ndatatotal2*nstates)]
#define hmmy(i,j,k) hmmy[(k+j*nstates+i*ndatatotal2*nstates)]

#define hmmprobs(i,j,k) hmmprobs[(i+j*noutputs+k*noutputs*noutputs)]
#define hmmprobsho(i,j,k) hmmprobsho[(i+j*noutputsho+k*noutputsho*noutputsho)]
#define hmmprobsho1(i,j,k) hmmprobsho1[(i+j*noutputsho1+k*noutputsho1*noutputsho1)]
#define hmmprobsho2(i,j,k) hmmprobsho2[(i+j*noutputsho2+k*noutputsho2*noutputsho2)]
#define hmmprobsho3(i,j,k) hmmprobsho3[(i+j*noutputsho3+k*noutputsho3*noutputsho3)]
#define hmmesum(i,j,k) hmmesum[(i+j*noutputs+k*noutputs*noutputs)]
#define hmmesums(i,j,k) hmmesums[(i+j*noutputs+k*noutputs*noutputs)]


int GetRandomInt(int minimum, int maximum)
{
     double r = ( (double)rand() / (double)(RAND_MAX) ); 
	 return (int)(r*((maximum+1)-minimum)) + minimum;
}

double GetRandomDouble(double minimum, double maximum)
{	 
     double r = ( (double)rand() / (double)(RAND_MAX) ); 
	 return r*(maximum-minimum) + minimum;
}


int** CreateIntArray(int rows, int cols)
{
	int i,j;
		
	int** a = new int*[rows];
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new int[cols];
		
		for(j = 0; j < cols; j++)
		{
			a[i][j] = 0;
		}
	}
	
	return a;
}

char** CreateCharArray(int rows, int cols)
{
	int i,j;
		
	char** a = new char*[rows];
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new char[cols];
	}
	
	return a;
}


short int** CreateShortIntArray(int rows, int cols)
{
	int i,j;
		
	short int** a = new short int*[rows];
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new short int[cols];
		
		for(j = 0; j < cols; j++)
		{
			a[i][j] = 0;
		}
	}
	
	return a;
}



string ReverseString(string s)
{
	return string ( s.rbegin(), s.rend() );
}

double** CreateDoubleArray(int rows, int cols)
{
	int i,j;
	
	double** a = new double*[rows];	
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new double[cols];
		
		for(j = 0; j < cols; j++)
		{
			a[i][j] = 0;
		}
	}
	
	return a;
}


float** CreateFloatArray(int rows, int cols)
{
	int i,j;
	
	float** a = new float*[rows];	
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new float[cols];
		
		for(j = 0; j < cols; j++)
		{
			a[i][j] = 0;
		}
	}
	
	return a;
}


string** CreateStringArray(int rows, int cols)
{
	int i,j;
		
	string** a = new string*[rows];
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new string[cols];
		
		for(j = 0; j < cols; j++)
		{
			a[i][j] = "";
		}
	}
	
	return a;
}


void Clear2DimArray(int** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}

void Clear2DimArray(short int** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}

void Clear2DimArray(double** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}


void Clear2DimArray(float** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}


void Clear2DimArray(string** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}
	
int ConvertNucToInt(char c)
{
	switch(c)
	{
		case 'A':
		return 0;		
		case 'C':
		return 1;
		case 'G':
		return 2;
		case 'T':
		return 3;	
		default:
		return 0;				
	}
}

char ConvertIntToNuc(int i)
{
	switch(i)
	{
		case 0:
		return 'A';		
		case 1:
		return 'C';
		case 2:
		return 'G';
		case 3:
		return 'T';	
		default:
		return ' ';
	}
}


int ConvertCharToIndex2(char c)
{
	
	int n;
	
	n = -1;
	
	switch(toupper(c))
	{
		case 'A': n = 0; break;
		case 'C': n = 1; break;
		case 'G': n = 2; break;
		case 'T': n = 3; break;		
	}	
	
	return n;								
}


int ConvertCharToIndex(char c, int nbase)
{
	
	int n;
	
	if(nbase == 4)
	{		
		switch(toupper(c))
		{
			case 'A': n = 0; break;
			case 'C': n = 1; break;
			case 'G': n = 2; break;
			case 'T': n = 3; break;
		}									
	}
	else if(nbase == 20)
	{
		switch(c)
		{
			
			
			case 'a': n = 0; break;
			case 'r': n = 1; break;
			case 'n': n = 2; break;
			case 'd': n = 3; break;
			case 'c': n = 4; break;
			case 'q': n = 5; break;
			case 'e': n = 6; break;
			case 'g': n = 7; break;
			case 'h': n = 8; break;
			case 'i': n = 9; break;
			case 'l': n = 10; break;
			case 'k': n = 11; break;
			case 'm': n = 12; break;
			case 'f': n = 13; break;
			case 'p': n = 14; break;
			case 's': n = 15; break;
			case 't': n = 16; break;
			case 'w': n = 17; break;
			case 'y': n = 18; break;
			case 'v': n = 19; break;
		}						
	}
	
	return n;
	
}


int ConvertOrder1StringToIndex2(int n1, int n2)
{
	return 4*n2+n1;
}

int ConvertOrder2StringToIndex2(int n1, int n2, int n3)
{
	return 16*n3+4*n2+n1;
}

int ConvertOrder3StringToIndex2(int n1, int n2, int n3, int n4)
{
	return 64*n4+16*n3+4*n2+n1;
}

int ConvertOrder4StringToIndex2(int n1, int n2, int n3, int n4, int n5)
{
	return n5*256+64*n4+16*n3+4*n2+n1;
}

int ConvertOrder5StringToIndex2(int n1, int n2, int n3, int n4, int n5, int n6)
{
	return n6*1024+n5*256+64*n4+16*n3+4*n2+n1;
}

int ConvertOrder6StringToIndex2(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{
	return n7*4096+n6*1024+n5*256+64*n4+16*n3+4*n2+n1;
}

int ConvertOrder7StringToIndex2(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	return n8*16384+n7*4096+n6*1024+n5*256+64*n4+16*n3+4*n2+n1;
}

int ConvertOrder8StringToIndex2(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9)
{
	return n9*65536+n8*16384+n7*4096+n6*1024+n5*256+64*n4+16*n3+4*n2+n1;
}

int ConvertOrder9StringToIndex2(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10)
{
	return n10*262144+n9*65536+n8*16384+n7*4096+n6*1024+n5*256+64*n4+16*n3+4*n2+n1;
}



int ConvertOrder2StringToIndex(int* seq, int n)
{
	return 16*seq[n]+4*seq[n-1]+seq[n-2];
}

int ConvertOrder3StringToIndex(int* seq, int n)
{
	return 64*seq[n]+16*seq[n-1]+4*seq[n-2]+seq[n-3];
}

int ConvertOrderStringToIndex(int* seq, int n, int norder)
{
	switch(norder)
	{
		case 0: return seq[n];			
		case 1: return 4*seq[n]+seq[n-1];
		case 2: return 16*seq[n]+4*seq[n-1]+seq[n-2];
		case 3: return 64*seq[n]+16*seq[n-1]+4*seq[n-2]+seq[n-3];				
		case 4: return 256*seq[n]+64*seq[n-1]+16*seq[n-2]+4*seq[n-3]+seq[n-4];
		case 5: return 1024*seq[n]+256*seq[n-1]+64*seq[n-2]+16*seq[n-3]+4*seq[n-4]+seq[n-5];
		case 6: return 4096*seq[n]+1024*seq[n-1]+256*seq[n-2]+64*seq[n-3]+16*seq[n-4]+4*seq[n-5]+seq[n-6];
		case 7: return 16384*seq[n]+4096*seq[n-1]+1024*seq[n-2]+256*seq[n-3]+64*seq[n-4]+16*seq[n-5]+4*seq[n-6]+seq[n-7];
		case 8: return 65536*seq[n]+16384*seq[n-1]+4096*seq[n-2]+1024*seq[n-3]+256*seq[n-4]+64*seq[n-5]+16*seq[n-6]+4*seq[n-7]+seq[n-8];
		default: return 0;
	}
}

int CheckIfOrderOk(int* seq, int n, int norder)
{
	int i;
	int bok;
	
	bok = 1;
	for(i = 0; i <= norder; i++)
	{
		if(seq[n-i] == -1) 
		{
			bok = 0;
			break;
		}
	}
	
	return bok;
}

int ConvertOrder0StringToIndex2(int* seq, int n)
{
	return seq[n];
}

int ConvertOrder1StringToIndex2(int* seq, int n)
{
	if(seq[n] == -1 || seq[n-1] == -1) return -1;
	else return 4*seq[n]+seq[n-1];
}


int ConvertOrder2StringToIndex2(int* seq, int n)
{
	if(seq[n] == -1 || seq[n-1] == -1 || seq[n-2] == -1) return -1;
	else return 16*seq[n]+4*seq[n-1]+seq[n-2];
}

int ConvertOrder3StringToIndex2(int* seq, int n)
{
	if(seq[n] == -1 || seq[n-1] == -1 || seq[n-2] == -1 || seq[n-3] == -1) return -1;
	else return 64*seq[n]+16*seq[n-1]+4*seq[n-2]+seq[n-3];
}



int ConvertOrder0StringToIndex2ForTrain(int* seq, int n)
{
	return seq[n];
}

int ConvertOrder1StringToIndex2ForTrain(int* seq, int n)
{
	int n1 = n;
	int n2 = n-1;
	
	while(n1 >= 0 && seq[n1] == -1)
	{
		n1--;
		n2--;		
	}
	
	if(n2 < 0) return -1;
	
	while(n2 >= 0 && seq[n2] == -1)
	{
		n2--;		
	}
	
	if(n2 < 0) return -1;
	
	
	return 4*seq[n1]+seq[n2];
}


int ConvertOrder2StringToIndex2ForTrain(int* seq, int n)
{
	int n1 = n;
	int n2 = n-1;
	int n3 = n-2;
	
	while(n1 >= 0 && seq[n1] == -1)
	{
		n1--;
		n2--;
		n3--;
	}
	
	if(n3 < 0) return -1;
	
	while(n2 >= 0 && seq[n2] == -1)
	{
		n2--;
		n3--;
	}
	
	if(n3 < 0) return -1;
	
	while(n3 >= 0 && seq[n3] == -1)
	{		
		n3--;		
	}
	
	if(n3 < 0) return -1;
		
	return 16*seq[n1]+4*seq[n2]+seq[n3];
}


int ConvertOrder3StringToIndex2ForTrain(int* seq, int n)
{
	int n1 = n;
	int n2 = n-1;
	int n3 = n-2;
	int n4 = n-3;
	
	if(seq[n1] == -1) return -1;
	
	while(n2 >= 0 && seq[n2] == -1)
	{
		n2--;
		n3--;
		n4--;
	}
	
	if(n4 < 0) return -1;
	
	while(n3 >= 0 && seq[n3] == -1)
	{		
		n3--;
		n4--;
	}
	
	if(n4 < 0) return -1;
	
	while(n4 >= 0 && seq[n4] == -1)
	{
		n4--;
	}
	
	if(n4 < 0) return -1;
	
	return 64*seq[n1]+16*seq[n2]+4*seq[n3]+seq[n4];
}

int ConvertOrderStringToIndex2(int* seq, int n, int norder)
{
	switch(norder)
	{
		case 1: return ConvertOrder0StringToIndex2(seq, n);
		case 2: return ConvertOrder1StringToIndex2(seq, n);
		case 3: return ConvertOrder2StringToIndex2(seq, n);
		case 4: return ConvertOrder3StringToIndex2(seq, n);
		default: return 0;
	}
}



int ConvertOrderStringToIndex2ForTrain(int* seq, int n, int norder)
{
	switch(norder)
	{
		case 1: return ConvertOrder0StringToIndex2ForTrain(seq, n);
		case 2: return ConvertOrder1StringToIndex2ForTrain(seq, n);
		case 3: return ConvertOrder2StringToIndex2ForTrain(seq, n);
		case 4: return ConvertOrder3StringToIndex2ForTrain(seq, n);
		default: return 0;
	}
}

int StringToIndex(string str)
{
	int i, k, index;
	index = 0;
	
	int nlen = str.length();
	
	for(i = 0; i < nlen; i++)
	{
		if(str[i] == 'A') k = 0;
		else if(str[i] == 'C') k = 1;
		else if(str[i] == 'G') k = 2;
		else if(str[i] == 'T') k = 3;															
							
		index += (int)pow(4, i)*k;
	}
	
	return index;
}



int StringToIndex(string str, int nbases)
{
	char c;
	int i, k, index;
	index = 0;
	
	int nlen = str.length();		
		
	if(nbases == 4)
	{
		for(i = 0; i < nlen; i++)
		{
			/*
			if(str[i] == 'A') k = 0;
			else if(str[i] == 'G') k = 1;
			else if(str[i] == 'C') k = 2;
			else if(str[i] == 'T') k = 3;															
			*/
			
			switch(toupper(str[i]))
			{
				case 'A': k = 0; break;
				case 'C': k = 1; break;
				case 'G': k = 2; break;
				case 'T': k = 3; break;
			}
							
			index += (int)pow(4, i)*k;
		}
	}
	else if(nbases == 20)
	{
		
		if(nlen == 1)
		{	
			c = str[0];	
			switch(c)
			{				
				case 'a': k = 0; break;
				case 'r': k = 1; break;
				case 'n': k = 2; break;
				case 'd': k = 3; break;
				case 'c': k = 4; break;
				case 'q': k = 5; break;
				case 'e': k = 6; break;
				case 'g': k = 7; break;
				case 'h': k = 8; break;
				case 'i': k = 9; break;
				case 'l': k = 10; break;
				case 'k': k = 11; break;
				case 'm': k = 12; break;
				case 'f': k = 13; break;
				case 'p': k = 14; break;
				case 's': k = 15; break;
				case 't': k = 16; break;
				case 'w': k = 17; break;
				case 'y': k = 18; break;
				case 'v': k = 19; break;
			}
			
			return k;	
		}
		else if(nlen == 2)
		{		
			c = str[0];				
			switch(c)
			{								
				case 'a': k = 0; break;
				case 'r': k = 1; break;
				case 'n': k = 2; break;
				case 'd': k = 3; break;
				case 'c': k = 4; break;
				case 'q': k = 5; break;
				case 'e': k = 6; break;
				case 'g': k = 7; break;
				case 'h': k = 8; break;
				case 'i': k = 9; break;
				case 'l': k = 10; break;
				case 'k': k = 11; break;
				case 'm': k = 12; break;
				case 'f': k = 13; break;
				case 'p': k = 14; break;
				case 's': k = 15; break;
				case 't': k = 16; break;
				case 'w': k = 17; break;
				case 'y': k = 18; break;
				case 'v': k = 19; break;
			}
			
			index = k;
			
			c = str[1];				
			switch(c)
			{								
				case 'a': k = 0; break;
				case 'r': k = 1; break;
				case 'n': k = 2; break;
				case 'd': k = 3; break;
				case 'c': k = 4; break;
				case 'q': k = 5; break;
				case 'e': k = 6; break;
				case 'g': k = 7; break;
				case 'h': k = 8; break;
				case 'i': k = 9; break;
				case 'l': k = 10; break;
				case 'k': k = 11; break;
				case 'm': k = 12; break;
				case 'f': k = 13; break;
				case 'p': k = 14; break;
				case 's': k = 15; break;
				case 't': k = 16; break;
				case 'w': k = 17; break;
				case 'y': k = 18; break;
				case 'v': k = 19; break;
			}
			
			index += 20*k;
			
			return index;
			
		}
		else
		{
			for(i = 0; i < nlen; i++)
			{
				switch(str[i])
				{								
					case 'a': k = 0; break;
					case 'r': k = 1; break;
					case 'n': k = 2; break;
					case 'd': k = 3; break;
					case 'c': k = 4; break;
					case 'q': k = 5; break;
					case 'e': k = 6; break;
					case 'g': k = 7; break;
					case 'h': k = 8; break;
					case 'i': k = 9; break;
					case 'l': k = 10; break;
					case 'k': k = 11; break;
					case 'm': k = 12; break;
					case 'f': k = 13; break;
					case 'p': k = 14; break;
					case 's': k = 15; break;
					case 't': k = 16; break;
					case 'w': k = 17; break;
					case 'y': k = 18; break;
					case 'v': k = 19; break;
				}
			}			
			index += (int)pow(20, i)*k;
		}	
	}
	
	return index;
}



void ConvertStringToInts(string s, int* data, int nbase)
{
	int i;
	
	int nlen = s.length();
	
	for(i = 0; i < nlen; i++) data[i] = ConvertCharToIndex(s[i], nbase);
	
	return;
}


void ConvertStringToInts2(string s, int* data)
{
	int i;
	
	int nlen = s.length();
	
	for(i = 0; i < nlen; i++) data[i] = ConvertCharToIndex2(s[i]);
	
	return;
}





void OutputToFileAlignProbs(ofstream &savehmmfile, double* hmmprobsho, int norder, int narrayindex, int noutputsho)
{	
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int n;
	int n1;
	int n2;	
	int nseq1;
	int nseq2;
	int nbases;		
	double d;
	string s;
	
	nbases = 4;
	
	savehmmfile << "p" << norder << ";";
	switch(norder)
	{
		case 0:			
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
				{									
					nseq1 = i1;																		
					nseq2 = i2;		
								
					savehmmfile << hmmprobsho(nseq1, nseq2, narrayindex) << ",";																																						
				}
		break;
		case 1:
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
						{									
							nseq1 = ConvertOrder1StringToIndex2(i1, i2);																		
							nseq2 = ConvertOrder1StringToIndex2(i3, i4);									
									
							savehmmfile << hmmprobsho(nseq1, nseq2, narrayindex) << ",";																																		
						}
		break;
		case 2:
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
							for(i5 = 0; i5 < nbases; i5++)					
								for(i6 = 0; i6 < nbases; i6++)
								{									
									nseq1 = ConvertOrder2StringToIndex2(i1, i2, i3);																		
									nseq2 = ConvertOrder2StringToIndex2(i4, i5, i6);									
									
									savehmmfile << hmmprobsho(nseq1, nseq2, narrayindex) << ",";																																
								}
		break;
		case 3:			
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
							for(i5 = 0; i5 < nbases; i5++)					
								for(i6 = 0; i6 < nbases; i6++)
									for(i7 = 0; i7 < nbases; i7++)
										for(i8 = 0; i8 < nbases; i8++)		
										{									
											nseq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);																		
											nseq2 = ConvertOrder3StringToIndex2(i5, i6, i7, i8);									
									
											savehmmfile << hmmprobsho(nseq1, nseq2, narrayindex) << ",";																																
										}
		break;
	}
	
	savehmmfile << "\n";
	
}




vfwdbackpair ViterbiTrainVariedOrderHMM(int* data1, int* data2, string sdata1, string sdata2, int ndatalen1, 
									   int ndatalen2,int noutputs, int nstates, double* hmmprobs, double* hmmprobsho1, 
									   double* hmmprobsho2, double* hmmprobsho3, double* as[],  int* statetypes, 
									   double* hmmesum, double* asum[], int nbeginstate, int nendstate, int norder, 
									   int noutputsho1, int noutputsho2, int noutputsho3)
{
	
	int i;
	int j;
	int k;
	int l;
	int n1;
	int n2;
	int indexi;
	int indexj; 
	int ndata1; 
	int ndata2; 
	int index1; 
	int index2; 
	int nstatetype; 
	int imax; 
	int vstate; 	
	int nstart;
	int ndatatotal1;
	int ndatatotal2;
	int ndataproduct; 
	int noutputproduct; 
	int nindex;
	int norderuse;
	int xcount;
	int ycount;
	int xsum;
	int ysum; 
	int xmax; 
	int ymax;
	
	string s1; 
	string s2;
	string vouts1;
	string vouts2;
	
	double el; 
	double el1; 
	double el2; 
	double dmax;  
	double vmax; 
	double vsum;  
	
	vfwdbackpair vf;
	
	nstart = 0;
	vouts1 = ""; 
	vouts2 = "";
	
	ndatatotal1 = ndatalen1+1; 
	ndatatotal2 = ndatalen2+1;
	ndataproduct = ndatatotal1*ndatatotal2;
	noutputproduct = noutputs*noutputs;
	
	double* hmmv = new double[(ndatatotal1)*(ndatatotal2)*nstates];	
	int* hmmptrs = new int[(ndatatotal1)*(ndatatotal2)*nstates];	
	
	int* hmmx  = new int[(ndatatotal1)*(ndatatotal2)*nstates];
	int* hmmy  = new int[(ndatatotal1)*(ndatatotal2)*nstates];


	for(i = 0; i < nstates; i++)
	{		
		for(j = 0; j < ndatatotal1; j++)
		{
			for(k = 0; k < ndatatotal2; k++)
			{	
				hmmptrs(j, k, i) = -1;	
				
				hmmx(j,k,i) = 0;		
				hmmy(j,k,i) = 0;
			}			
		}		
		
		
		hmmv(0, 0, i) = AMIN;	
		hmmptrs(0, 0, i) = -1;
		hmmptrs(0, 1, i) = -1;
		hmmptrs(1, 0, i) = -1;
		hmmptrs(1, 1, i) = -1;
				
		if(i == nbeginstate)
		{
			hmmv(0, 0, i) = 0;	
		}
	}


	for(i = nstart; i < ndatatotal1; i++)
	{			
		for(j = nstart; j < ndatatotal2; j++)
		{					
			for(l = 0; l < nstates; l++)
			{	
			
				imax = -1;
				dmax = AMIN;
				xmax = 0;
				ymax = 0;

				el = 0;				
			
				
				for(k = 0; k < nstates; k++)
				{																								
					nstatetype = statetypes[l];																			
					vsum = AMIN;	
														
					if(as[k][l] > AMIN)
					{			
							
						switch(nstatetype) 
						{
							
							case MATCHSTATEHO:
							
								indexi = i-1;
								indexj = j-1;						
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;

								if(indexi >= 0 && indexj >= 0)
								{					
									xcount = hmmx[nindex];																								
									ycount = hmmy[nindex];			
															
									norderuse = min(xcount, ycount);
									if(norderuse > norder) norderuse = norder;
									
									if(indexi >= norderuse && indexj >= norderuse && norderuse > 0) 
									{
										ndata1 = ConvertOrderStringToIndex(data1,indexi,norderuse);
										ndata2 = ConvertOrderStringToIndex(data2,indexj,norderuse);
										
										if(norderuse == 3) el = hmmprobsho3(ndata1, ndata2, l);
										else if(norderuse == 2) el = hmmprobsho2(ndata1, ndata2, l);																		
										else el = hmmprobsho1(ndata1, ndata2, l);	

										if(el == 1)
										{
											ndata1 = data1[indexi];
											ndata2 = data2[indexj];
											
											el = hmmprobs(ndata1, ndata2, l);		
										}
										else if(el == 2)
										{
											ndata1 = ConvertOrderStringToIndex(data1,indexi,1);
											ndata2 = ConvertOrderStringToIndex(data2,indexj,1);
											
											el = hmmprobsho1(ndata1, ndata2, l);	
										}										
										else if(el == 3)
										{
											ndata1 = ConvertOrderStringToIndex(data1,indexi,2);
											ndata2 = ConvertOrderStringToIndex(data2,indexj,2);
											
											el = hmmprobsho2(ndata1, ndata2, l);	
										}	
									}
									else 
									{
										ndata1 = data1[indexi];
										ndata2 = data2[indexj];	
									
										el = hmmprobs(ndata1, ndata2, l);
									}									
									
									vsum = hmmv[nindex] + as[k][l];																		
										
									xsum = xcount+1;
									ysum = ycount+1;										
								}	
								
							break;

							case OUTPUTSTATEXHO:
							
								indexi = i-1;
								indexj = j;
							
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
							
								if(indexi >= 0 && indexj >= 0)
								{												
									xcount = hmmx[nindex];																										
															
									norderuse = xcount;
									if(norderuse > norder) norderuse = norder;
																																																			
									if(indexi >= norderuse && norderuse > 0) 
									{
										ndata1 = ConvertOrderStringToIndex(data1,indexi,norderuse);										

										if(norderuse == 3) el = hmmprobsho3(ndata1, 0, l); 																
										else if(norderuse == 2) el = hmmprobsho2(ndata1, 0, l); 																
										else el = hmmprobsho1(ndata1, 0, l); 

										if(el == 1)
										{
											ndata1 = data1[indexi];																		
											el = hmmprobs(ndata1, 0, l);		
										}
										else if(el == 2)
										{
											ndata1 = ConvertOrderStringToIndex(data1,indexi,1);																					
											el = hmmprobsho1(ndata1, 0, l);		
										}	
										else if(el == 3)
										{
											ndata1 = ConvertOrderStringToIndex(data1,indexi,2);																					
											el = hmmprobsho2(ndata1, 0, l);		
										}
										
									}
									else
									{
										ndata1 = data1[indexi];	
										el = hmmprobs(ndata1, 0, l); 
									} 
									
									vsum = hmmv[nindex] + as[k][l];
									
									xsum = xcount+1;
									ysum = 0;
																																
								}
								else vsum = AMIN;														
								
							break;
							
							
							case OUTPUTSTATEYHO:
							
								indexi = i;
								indexj = j-1;
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
						
								if(indexi >= 0 && indexj >= 0)
								{												
									ycount = hmmy[nindex];			
															
									norderuse = ycount;
									if(norderuse > norder) norderuse = norder;
																																																							
									if(indexj >= norderuse && norderuse > 0) 
									{																				
										ndata2 = ConvertOrderStringToIndex(data2,indexj,norderuse);	
										
										if(norderuse == 3) el = hmmprobsho3(ndata2, 0, l);
										else if(norderuse == 2) el = hmmprobsho2(ndata2, 0, l);
										else el = hmmprobsho1(ndata2, 0, l); 	
										
										if(el == 1)
										{
											ndata2 = data2[indexj];																	
											el = hmmprobs(ndata2, 0, l);		
										}
										else if(el == 2)
										{
											ndata2 = ConvertOrderStringToIndex(data2,indexj,1);																			
											el = hmmprobsho1(ndata2, 0, l);		
										}
										else if(el == 3)
										{
											ndata2 = ConvertOrderStringToIndex(data2,indexj,2);																					
											el = hmmprobsho2(ndata2, 0, l);		
										}	

									}
									else
									{
										ndata2 = data2[indexj];	
										el = hmmprobs(ndata2, 0, l); 
									} 
							
									vsum = hmmv[nindex] + as[k][l];									
																			
									xsum = 0;
									ysum = ycount+1;																								
																																
								}
								else vsum = AMIN;														
								
							break;
							
							
						} // switch																																	
																				
					} //as[k][l] > AMIN
										
					if(vsum > dmax)
					{
						dmax = vsum;
						imax = k;	
						
						xmax = xsum;
						ymax = ysum;
					}
					
				} //for(k = 0; k < nstates; k++)
				
		
				
				if(i > nstart || j > nstart)
				{	
					nindex = l+j*nstates+i*ndatatotal2*nstates;	
					
					hmmv[nindex] = el+dmax;																											
					hmmptrs[nindex] = imax;												
					hmmx[nindex] = xmax;													
					hmmy[nindex] = ymax;		
				}
				
				
						
			} //for(l = 0; l < nstates; l++)
			

		}			
	}
			

	vmax = AMIN;								
	for(k = 0; k < nstates; k++)
	{	
		nstatetype = statetypes[k];
		
		hmmv((ndatalen1),(ndatalen2),k) += as[k][nendstate];
				
		if(hmmv((ndatalen1),(ndatalen2),k) > vmax) 
		{
			vmax = hmmv((ndatalen1),(ndatalen2),k); 		
			vstate = k;						
		}		
	}		

	
	vf.vmax = vmax;

		
	n1 = ndatalen1;
	n2 = ndatalen2;
			
	vector<int> v;
	vector<int>::iterator it;

	it = v.begin();
	it = v.insert(it, vstate);
		
	string states = "";
		
		
		

	l = vstate;
	int ntotal = 0;
	
	while(n1 > 1 || n2 > 1)
	{		
		nstatetype = statetypes[l];	
		
			
		k = hmmptrs(n1, n2, l);			
		asum[k][l]++;		
		ntotal++;					
						
		it = v.insert(it, k);
		
		if(nstatetype == OUTPUTSTATEX || nstatetype == OUTPUTSTATEXHO) 		n1--;	
		else if(nstatetype == OUTPUTSTATEY || nstatetype == OUTPUTSTATEYHO) n2--;
		else if(nstatetype == MATCHSTATE || nstatetype == MATCHSTATEHO)
		{
			n1--;
			n2--;
		}
		
		l = k;		
	}
	
	asum[0][l]++;

		
	
	j = 0;
	k = 0;
		
	for (i=0; i<v.size(); i++)
	{		
		vstate = v.at(i);
						
		if(statetypes[vstate] == OUTPUTSTATEX || statetypes[vstate] == OUTPUTSTATEXHO) 
		{							
			vouts1 += sdata1[j];			
			vouts2 += "-";							
			states += "x";			
			j++;				
				
		}
		else if(statetypes[vstate] == OUTPUTSTATEY || statetypes[vstate] == OUTPUTSTATEYHO) 
		{					
			vouts1 += "-";			
			vouts2 += sdata2[k];						
			states += "y";								
			k++;
		}
		else if(statetypes[vstate] == MATCHSTATE || statetypes[vstate] == MATCHSTATEHO) 
		{
			vouts1 += sdata1[j];
			vouts2 += sdata2[k];				
			j++;
			k++;						
			states += "m";				
		}
	}
	
	vf.vouts1 = vouts1;
	vf.vouts2 = vouts2;
	


	delete[] hmmv;
	delete[] hmmptrs;
	delete[] hmmx;
	delete[] hmmy;
	
	
	
	return vf;
}





vfwdbackpair ViterbiTrain(int* data1, int* data2, string sdata1, string sdata2, int ndatalen1, 
									   int ndatalen2,int noutputs, int nstates, double* hmmprobs, double* hmmprobsho, 
									   double* as[],  int* statetypes, double* hmmesum, double* asum[],
									   int nbeginstate, int nendstate, int norder, int noutputsho)
{
	
	int i;
	int j;
	int k;
	int l;
	int n1;
	int n2;
	int indexi;
	int indexj; 
	int ndata1; 
	int ndata2; 
	int index1; 
	int index2; 
	int nstatetype; 
	int imax; 
	int vstate; 	
	int nstart;
	int ndatatotal1;
	int ndatatotal2;
	int ndataproduct; 
	int noutputproduct; 
	int nindex;
	int norderuse;
	int xcount;
	int ycount;
	int xsum;
	int ysum; 
	int xmax; 
	int ymax;
	
	string s1; 
	string s2;
	string vouts1;
	string vouts2;
	
	double el; 
	double el1; 
	double el2; 
	double dmax;  
	double vmax; 
	double vsum;  
	
	vfwdbackpair vf;
	
	nstart = 0;
	vouts1 = ""; 
	vouts2 = "";
	
	ndatatotal1 = ndatalen1+1; 
	ndatatotal2 = ndatalen2+1;
	ndataproduct = ndatatotal1*ndatatotal2;
	noutputproduct = noutputs*noutputs;
	
	double* hmmv = new double[(ndatatotal1)*(ndatatotal2)*nstates];	
	int* hmmptrs = new int[(ndatatotal1)*(ndatatotal2)*nstates];	
	
	int* hmmx  = new int[(ndatatotal1)*(ndatatotal2)*nstates];
	int* hmmy  = new int[(ndatatotal1)*(ndatatotal2)*nstates];

	for(i = 0; i < nstates; i++)
	{		
		for(j = 0; j < ndatatotal1; j++)
		{
			for(k = 0; k < ndatatotal2; k++)
			{	
				hmmptrs(j, k, i) = -1;					
				hmmx(j,k,i) = 0;		
				hmmy(j,k,i) = 0;
			}			
		}		
		
		
		hmmv(0, 0, i) = AMIN;	
		hmmptrs(0, 0, i) = -1;
		hmmptrs(0, 1, i) = -1;
		hmmptrs(1, 0, i) = -1;
		hmmptrs(1, 1, i) = -1;
		
		
		if(i == nbeginstate)
		{
			hmmv(0, 0, i) = 0;	
		}
	}


	for(i = nstart; i < ndatatotal1; i++)
	{	
		for(j = nstart; j < ndatatotal2; j++)
		{					
			for(l = 0; l < nstates; l++)
			{	
			
				imax = -1;
				dmax = AMIN;
				xmax = 0;
				ymax = 0;

				el = 0;				
			
				
				for(k = 0; k < nstates; k++)
				{																								
					nstatetype = statetypes[l];																			
					vsum = AMIN;	
														
					if(as[k][l] > AMIN)
					{			
							
						switch(nstatetype) 
						{
							
							case OUTPUTSTATEX:
							
								indexi = i-1;
								indexj = j;
							
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
							
								if(indexi >= 0 && indexj >= 0)
								{																																	
									ndata1 = data1[indexi];																												
									el = hmmprobs(ndata1, 0, l);													
									vsum = hmmv[nindex] + as[k][l];
								}
								else vsum = AMIN;														
								
							break;
						
							case OUTPUTSTATEY:
						
								indexi = i;
								indexj = j-1;
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
						
								if(indexi >= 0 && indexj >= 0)
								{																																		
									ndata2 = data2[indexj];																			
									el = hmmprobs(ndata2, 0, l);										
									vsum = hmmv[nindex] + as[k][l];																																																																														
								}
								else vsum = AMIN;
	
							break;
							
							case MATCHSTATE:
							
								indexi = i-1;
								indexj = j-1;						
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
						
								if(indexi >= 0 && indexj >= 0)
								{																																			
									ndata1 = data1[indexi];
									ndata2 = data2[indexj];	
																		
									el = hmmprobs(ndata1, ndata2, l); 																										
									vsum = hmmv[nindex] + as[k][l];																																																			
								}	
								
							break;
							
							
							
							case MATCHSTATEHO:
							
								indexi = i-1;
								indexj = j-1;						
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;

								if(indexi >= 0 && indexj >= 0)
								{					
									xcount = hmmx[nindex];																								
									ycount = hmmy[nindex];			
															
									norderuse = min(xcount, ycount);
									if(norderuse > norder) norderuse = norder;
									
									if(indexi >= norderuse && indexj >= norderuse && norderuse > 0) 
									{
										ndata1 = ConvertOrderStringToIndex(data1,indexi,norderuse);
										ndata2 = ConvertOrderStringToIndex(data2,indexj,norderuse);

										el = hmmprobsho(ndata1, ndata2, l);																		
									}
									else 
									{
										ndata1 = data1[indexi];
										ndata2 = data2[indexj];	
									
										el = hmmprobs(ndata1, ndata2, l);
									}									
									
									vsum = hmmv[nindex] + as[k][l];																		
										
									xsum = xcount+1;
									ysum = ycount+1;										
								}	
								
							break;

							case OUTPUTSTATEXHO:
							
								indexi = i-1;
								indexj = j;
							
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
							
								if(indexi >= 0 && indexj >= 0)
								{												
									xcount = hmmx[nindex];																										
															
									norderuse = xcount;
									if(norderuse > norder) norderuse = norder;
																																																			
									if(indexi >= norderuse && norderuse > 0) 
									{
										ndata1 = ConvertOrderStringToIndex(data1,indexi,norderuse);										
										el = hmmprobsho(ndata1, 0, l); 
									}
									else
									{
										ndata1 = data1[indexi];	
										el = hmmprobs(ndata1, 0, l); 
									} 
									
									vsum = hmmv[nindex] + as[k][l];
									
									xsum = xcount+1;
									ysum = 0;
																																
								}
								else vsum = AMIN;														
								
							break;
							
							
							case OUTPUTSTATEYHO:
							
								indexi = i;
								indexj = j-1;
								nindex = k+indexj*nstates+indexi*ndatatotal2*nstates;
								
								if(indexi >= 0 && indexj >= 0)
								{												
									ycount = hmmy[nindex];			
															
									norderuse = ycount;
									if(norderuse > norder) norderuse = norder;
																																																							
									if(indexj >= norderuse && norderuse > 0) 
									{																				
										ndata2 = ConvertOrderStringToIndex(data2,indexj,norderuse);	
										el = hmmprobsho(ndata2, 0, l); 
									}
									else
									{
										ndata2 = data2[indexj];	
										el = hmmprobs(ndata2, 0, l); 
									} 
							
									vsum = hmmv[nindex] + as[k][l];									
																			
									xsum = 0;
									ysum = ycount+1;																								
																																
								}
								else vsum = AMIN;														
								
							break;
							
							
						} // switch																																	
																				
					} //as[k][l] > AMIN
					
											
											
					if(vsum > dmax)
					{
						dmax = vsum;					
						imax = k;	
						
						xmax = xsum;
						ymax = ysum;
						
					}
					
				} //for(k = 0; k < nstates; k++)
				
		
				
				if(i > nstart || j > nstart)
				{	
					nindex = l+j*nstates+i*ndatatotal2*nstates;	
					
					hmmv[nindex] = el+dmax;																		
									
					hmmptrs[nindex] = imax;		
									
					hmmx[nindex] = xmax;													
					hmmy[nindex] = ymax;
						
				}
				
				
						
			} //for(l = 0; l < nstates; l++)
			

		}			
	}

	vmax = AMIN;								
	for(k = 0; k < nstates; k++)
	{	
		nstatetype = statetypes[k];
		
		hmmv((ndatalen1),(ndatalen2),k) += as[k][nendstate];
				
		if(hmmv((ndatalen1),(ndatalen2),k) > vmax) 
		{
			vmax = hmmv((ndatalen1),(ndatalen2),k); 		
			vstate = k;						
		}		
	}		

	
	vf.vmax = vmax;
		
	n1 = ndatalen1;
	n2 = ndatalen2;
			
	vector<int> v;
	vector<int>::iterator it;

	it = v.begin();
	it = v.insert(it, vstate);
		
	string states = "";	

	l = vstate;

	while(n1 > 1 || n2 > 1)
	{		
		nstatetype = statetypes[l];	
		
			
		k = hmmptrs(n1, n2, l);			
		asum[k][l]++;					
						
		it = v.insert(it, k);
		
		if(nstatetype == OUTPUTSTATEX || nstatetype == OUTPUTSTATEXHO) 		n1--;	
		else if(nstatetype == OUTPUTSTATEY || nstatetype == OUTPUTSTATEYHO) n2--;
		else if(nstatetype == MATCHSTATE || nstatetype == MATCHSTATEHO)
		{
			n1--;
			n2--;
		}
		
		l = k;		
	}
	
	asum[0][l]++;

	j = 0;
	k = 0;
		
	for (i=0; i<v.size(); i++)
	{		
		vstate = v.at(i);
						
		if(statetypes[vstate] == OUTPUTSTATEX || statetypes[vstate] == OUTPUTSTATEXHO) 
		{							
			vouts1 += sdata1[j];			
			vouts2 += "-";							
			states += "x";			
			j++;				
				
		}
		else if(statetypes[vstate] == OUTPUTSTATEY || statetypes[vstate] == OUTPUTSTATEYHO) 
		{					
			vouts1 += "-";			
			vouts2 += sdata2[k];						
			states += "y";								
			k++;
		}
		else if(statetypes[vstate] == MATCHSTATE || statetypes[vstate] == MATCHSTATEHO) 
		{
			vouts1 += sdata1[j];
			vouts2 += sdata2[k];				
			j++;
			k++;						
			states += "m";				
		}
	}
	
	vf.vouts1 = vouts1;
	vf.vouts2 = vouts2;
	

	delete[] hmmv;
	delete[] hmmptrs;
	delete[] hmmx;
	delete[] hmmy;
	
	
	
	return vf;
}








void OutputToFileOutputProbs(ofstream &savehmmfile, double* hmmprobsho, int norder, int narrayindex, int noutputsho)
{	
	int i1;
	int i2;
	int i3;
	int i4;
	int n;
	int n1;
	int n2;	
	int nbases;		
	double d;
	string s;
	
	nbases = 4;
	
	if(narrayindex == 1) savehmmfile << "q" << norder << "x;";
	else if(narrayindex == 3) savehmmfile << "q" << norder << "y;";
	
	switch(norder)
	{
		case 0:			
			for(i1 = 0; i1 < nbases; i1++)
			{	
				n = i1;				
				savehmmfile << hmmprobsho(n, 0, narrayindex) << ",";									
			}		
		break;
		case 1:
			for(i1 = 0; i1 < nbases; i1++)
			{	
				for(i2 = 0; i2 < nbases; i2++)
				{
					n = ConvertOrder1StringToIndex2(i1, i2);
					savehmmfile << hmmprobsho(n, 0, narrayindex) << ",";		
				}
			}
		break;
		case 2:
			for(i1 = 0; i1 < nbases; i1++)
			{	
				for(i2 = 0; i2 < nbases; i2++)
				{
					for(i3 = 0; i3 < nbases; i3++)
					{			
						n = ConvertOrder2StringToIndex2(i1, i2, i3);
						savehmmfile << hmmprobsho(n, 0, narrayindex) << ",";									
					}			
				}
			}
		break;
		case 3:
			for(i1 = 0; i1 < nbases; i1++)
			{	
				for(i2 = 0; i2 < nbases; i2++)
				{
					for(i3 = 0; i3 < nbases; i3++)
					{
						for(i4 = 0; i4 < nbases; i4++)
						{
							n = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
							savehmmfile << hmmprobsho(n, 0, narrayindex) << ",";				
						}				
					}			
				}
			}
		break;
	}
	
	savehmmfile << "\n";
}


	
	



void ReadAlignProbs(string line, double* hmmprobsho, int norder, int narrayindex, int noutputsho)
{	
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int n;
	int n1;
	int n2;	
	int nseq1;
	int nseq2;
	int nbases;		
	double d;
	string s;
	
	nbases = 4;
	
	n1 = 3;
	n2 = 0;
	switch(norder)
	{
		case 0:			
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
				{									
					nseq1 = i1;																		
					nseq2 = i2;		
									
					n2 = line.find(",", n1);	
					
					s = line.substr(n1,n2-n1);
					n1 = n2+1;
					
					d = atof(s.c_str());
					hmmprobsho(nseq1, nseq2, narrayindex) = d;																																			
				}
		break;
		case 1:
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
						{									
							nseq1 = ConvertOrder1StringToIndex2(i1, i2);																		
							nseq2 = ConvertOrder1StringToIndex2(i3, i4);									
									
							n2 = line.find(",", n1);	
					
							s = line.substr(n1,n2-n1);
							n1 = n2+1;
					
							d = atof(s.c_str());
							hmmprobsho(nseq1, nseq2, narrayindex) = d;																																			
						}
		break;
		case 2:
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
							for(i5 = 0; i5 < nbases; i5++)					
								for(i6 = 0; i6 < nbases; i6++)
								{									
									nseq1 = ConvertOrder2StringToIndex2(i1, i2, i3);																		
									nseq2 = ConvertOrder2StringToIndex2(i4, i5, i6);									
									
									n2 = line.find(",", n1);	
					
									s = line.substr(n1,n2-n1);
									n1 = n2+1;
					
									d = atof(s.c_str());
									hmmprobsho(nseq1, nseq2, narrayindex) = d;																																			
								}
		break;
		case 3:			
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
							for(i5 = 0; i5 < nbases; i5++)					
								for(i6 = 0; i6 < nbases; i6++)
									for(i7 = 0; i7 < nbases; i7++)
										for(i8 = 0; i8 < nbases; i8++)		
										{									
											nseq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);																		
											nseq2 = ConvertOrder3StringToIndex2(i5, i6, i7, i8);									
									
											n2 = line.find(",", n1);	
					
											s = line.substr(n1,n2-n1);
											n1 = n2+1;
					
											d = atof(s.c_str());
											hmmprobsho(nseq1, nseq2, narrayindex) = d;																																			
										}
		break;
	}
	
}





void ReadOutputProbs(string line, double* hmmprobsho, int norder, int narrayindex, int noutputsho)
{	
	int i1;
	int i2;
	int i3;
	int i4;
	int n;
	int n1;
	int n2;	
	int nbases;		
	double d;
	string s;
	
	nbases = 4;
	
	n1 = 4;
	n2 = 0;
	switch(norder)
	{
		case 0:			
			for(i1 = 0; i1 < nbases; i1++)
			{	
				n = i1;
				n2 = line.find(",", n1);	
					
				s = line.substr(n1,n2-n1);
				n1 = n2+1;
					
				d = atof(s.c_str());				
				hmmprobsho(n, 0, narrayindex) = d;							
			}		
		break;
		case 1:
			for(i1 = 0; i1 < nbases; i1++)
			{	
				for(i2 = 0; i2 < nbases; i2++)
				{
					n = ConvertOrder1StringToIndex2(i1, i2);
					n2 = line.find(",", n1);	
					
					s = line.substr(n1,n2-n1);
					n1 = n2+1;
					
					d = atof(s.c_str());				
					hmmprobsho(n, 0, narrayindex) = d;			
				}
			}
		break;
		case 2:
			for(i1 = 0; i1 < nbases; i1++)
			{	
				for(i2 = 0; i2 < nbases; i2++)
				{
					for(i3 = 0; i3 < nbases; i3++)
					{			
						n = ConvertOrder2StringToIndex2(i1, i2, i3);
						n2 = line.find(",", n1);	
					
						s = line.substr(n1,n2-n1);
						n1 = n2+1;
					
						d = atof(s.c_str());				
						hmmprobsho(n, 0, narrayindex) = d;										
					}			
				}
			}
		break;
		case 3:
			for(i1 = 0; i1 < nbases; i1++)
			{	
				for(i2 = 0; i2 < nbases; i2++)
				{
					for(i3 = 0; i3 < nbases; i3++)
					{
						for(i4 = 0; i4 < nbases; i4++)
						{
							n = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
							n2 = line.find(",", n1);	
					
							s = line.substr(n1,n2-n1);
							n1 = n2+1;
					
							d = atof(s.c_str());				
							hmmprobsho(n, 0, narrayindex) = d;					
						}				
					}			
				}
			}
		break;
	}
	
}
	


























	
int CreateProbMatrixFromTrainDataVariedOrder(string sfilein1, string sfilein2, double* hmmprobs, double* hmmprobsho1, 
											  double* hmmprobsho2, double* hmmprobsho3, int noutputs, int noutputsho1, 
											  int noutputsho2, int noutputsho3)
{	
	int i;
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	
	int nbases;	
	int totalcounts;
	int totalxcounts;
	int totalycounts;
	int n1seq1;
	int n1seq2;
	int n2seq1;
	int n2seq2;
	int n3seq1;
	int n3seq2;
	int n4seq1;
	int n4seq2;
	int nlengtho1;
	int nlengtho2;
	int nlengtho3;
	int nlengtho4;
	int nordermax;
	int nmincount;
	
	int nseq1;
	int nseq2;
	int nSeqLength1;
	int nSeqLength2;
	
	double d;
	double dminval;
	
	string line;
	string s1;
	string s2;
	string sref;
	string sread;
	
	int* seq1;
	int* seq2;
	
	
	ifstream finref(sfilein1.c_str());
	ifstream finread(sfilein2.c_str());
	
	if(finref.good() == 0)
	{
		cout << "Problem accessing " << sfilein1 << ". Please check file exists.\n\n";
		return 0;
	}
	else if(finread.good() == 0)
	{
		cout << "Problem accessing " << sfilein2 << ". Please check file exists.\n\n";
		return 0;
	}
		
	nordermax = 4;
	
	nlengtho1 = nordermax;
	nlengtho2 = (int)pow(4, nordermax-2);
	nlengtho3 = (int)pow(4, nordermax-1);
	nlengtho4 = (int)pow(4, nordermax);
	

	//1st order
	int** countso1 = CreateIntArray(nlengtho1, nlengtho1);
	int* countsxq1 = new int[nlengtho1];
	int* countsyq1 = new int[nlengtho1];

	//2nd order
	int** countso2o1 = CreateIntArray(nlengtho1, nlengtho1);
	int* countsxq2q1 = new int[nlengtho1];
	int* countsyq2q1 = new int[nlengtho1];
	
	int** countso2 = CreateIntArray(nlengtho2, nlengtho2);	
	int* countsxq2 = new int[nlengtho2];
	int* countsyq2 = new int[nlengtho2];


	//3rd order
	int** countso3o2 = CreateIntArray(nlengtho2, nlengtho2);
	int* countsxq3q2 = new int[nlengtho2];
	int* countsyq3q2 = new int[nlengtho2];
	
	int** countso3 = CreateIntArray(nlengtho3, nlengtho3);	
	int* countsxq3 = new int[nlengtho3];
	int* countsyq3 = new int[nlengtho3];
		
	
	//4th order
	int** countso4o3 = CreateIntArray(nlengtho3, nlengtho3);
	int* countsxq4q3 = new int[nlengtho3];
	int* countsyq4q3 = new int[nlengtho3];
		
	int** countso4 = CreateIntArray(nlengtho4, nlengtho4);
	int* countsxq4 = new int[nlengtho4];
	int* countsyq4 = new int[nlengtho4];
	
	double** probso1 = CreateDoubleArray(nlengtho1, nlengtho1);
	double** probso2 = CreateDoubleArray(nlengtho2, nlengtho2);
	double** probso3 = CreateDoubleArray(nlengtho3, nlengtho3);
	double** probso4 = CreateDoubleArray(nlengtho4, nlengtho4);
	
	
	
	dminval = 0.0000000001; //0.00001;
	
	nmincount = 100;
	
	for(i = 0; i < nlengtho1; i++) 
	{
		countsxq1[i] = 0;
		countsyq1[i] = 0;
		
		countsxq2q1[i] = 0;
		countsyq2q1[i] = 0;
	}
	
	for(i = 0; i < nlengtho2; i++) 
	{
		countsxq3q2[i] = 0;
		countsyq3q2[i] = 0;
		
		countsxq2[i] = 0;
		countsyq2[i] = 0;
	}
	
	for(i = 0; i < nlengtho3; i++) 
	{
		countsxq4q3[i] = 0;
		countsyq4q3[i] = 0;
		
		countsxq3[i] = 0;
		countsyq3[i] = 0;
	}
	
	for(i = 0; i < nlengtho4; i++) 
	{
		countsxq4[i] = 0;
		countsyq4[i] = 0;
	}
	
	totalcounts = 0;
	totalxcounts = 0;
	totalycounts = 0;
	
	sref = "";
	sread = "";
	getline(finref, line);	
	getline(finread, line);	
	
	
	while(finref.eof() == 0)
	{
		getline(finref, line);		
		sref += line;
		
		getline(finread, line);
		sread += line;
	}
	
	
	nSeqLength1 = sref.length();
	nSeqLength2 = sread.length();
			
	seq1 = new int[nSeqLength1];
	seq2 = new int[nSeqLength2];
	
	ConvertStringToInts2(sref, seq1);
	ConvertStringToInts2(sread, seq2);
			
	for(i = 0; i < nSeqLength1; i++)
	{
			
		if(i >= 1)
		{
			n1seq1 = ConvertOrderStringToIndex2(seq1, i-1, 1);
			n2seq1 = ConvertOrderStringToIndex2(seq1, i, 2);
				
			n1seq2 = ConvertOrderStringToIndex2(seq2, i-1, 1);
			n2seq2 = ConvertOrderStringToIndex2(seq2, i, 2);
				
			if(n1seq1 != -1 && n2seq1 != -1)
			{
				countsxq2q1[n1seq1]++;
				countsxq2[n2seq1]++;
			}
				
			if(n1seq2 != -1 && n2seq2 != -1)
			{
				countsyq2q1[n1seq2]++;
				countsyq2[n2seq2]++;
			}
				
			if(n1seq1 != -1 && n2seq1 != -1 && n1seq2 != -1 && n2seq2 != -1)
			{
				countso2o1[n1seq1][n1seq2]++;
				countso2[n2seq1][n2seq2]++;					
			}
		}
			
		if(i >= 2)
		{
			n2seq1 = ConvertOrderStringToIndex2(seq1, i-1, 2);
			n3seq1 = ConvertOrderStringToIndex2(seq1, i, 3);
				
			n2seq2 = ConvertOrderStringToIndex2(seq2, i-1, 2);
			n3seq2 = ConvertOrderStringToIndex2(seq2, i, 3);
				
			if(n2seq1 != -1 && n3seq1 != -1)
			{
				countsxq3q2[n2seq1]++;
				countsxq3[n3seq1]++;
			}
				
			if(n2seq2 != -1 && n3seq2 != -1)
			{
				countsyq3q2[n2seq2]++;
				countsyq3[n3seq2]++;
			}
				
			if(n2seq1 != -1 && n3seq1 != -1 && n2seq2 != -1 && n3seq2 != -1)
			{
				countso3o2[n2seq1][n2seq2]++;
				countso3[n3seq1][n3seq2]++;					
			}
		}
			
		if(i >= 3)
		{
			n3seq1 = ConvertOrder2StringToIndex2(seq1, i-1);
			n4seq1 = ConvertOrder3StringToIndex2(seq1, i);
				
			n3seq2 = ConvertOrder2StringToIndex2(seq2, i-1);
			n4seq2 = ConvertOrder3StringToIndex2(seq2, i);
			
			if(n3seq1 != -1 && n4seq1 != -1)
			{
				countsxq4q3[n3seq1]++;
				countsxq4[n4seq1]++;
			}
					
			if(n3seq2 != -1 && n4seq2 != -1)
			{
				countsyq4q3[n3seq2]++;
				countsyq4[n4seq2]++;
			}
				
			if(n3seq1 != -1 && n4seq1 != -1 && n3seq2 != -1 && n4seq2 != -1)
			{
				//cout << n3seq1 << "," << n3seq2 << "," << nlengtho3 << "\n";
				countso4o3[n3seq1][n3seq2]++;
				countso4[n4seq1][n4seq2]++;					
			}
		}
			
		nseq1 = seq1[i];
		nseq2 = seq2[i];
			
		if(nseq1 != -1 && nseq2 != -1)
		{
			countso1[nseq1][nseq2]++;
			totalcounts++;
		}
			
		if(nseq1 != -1)
		{
			countsxq1[nseq1]++;
			totalxcounts++;
		}
			
		if(nseq2 != -1)
		{
			countsyq1[nseq2]++;
			totalycounts++;
		}
		
	}
	
	delete[] seq1;
	delete[] seq2;

	/*
	 	for(i1 = 0; i1 < nbases; i1++)
		for(i2 = 0; i2 < nbases; i2++)
		{									
			nseq1 = i1;																		
			nseq2 = i2;		
									
			n2 = line.find(",", n1);	
					
			s = line.substr(n1,n2-n1);
			n1 = n2+1;
					
			d = atof(s.c_str());
			hmmprobs(nseq1, nseq2, narrayindex) = d;
		}

	*/
	
	/*
	int npseudocount = 1;
	int npseudocountalign = 16;
	int npseudocountoutput = 4;*/
	
	int npseudocount = 0;
	int npseudocountalign = 0;
	int npseudocountoutput = 0;
	
	nbases = 4;
	
	for(i1 = 0; i1 < nbases; i1++)
	{
		for(i2 = 0; i2 < nbases; i2++)
		{
			d = (double)(countso1[i1][i2]+npseudocount)/(double)(totalcounts+npseudocountalign);	
			if(d < dminval) d = dminval;
			
			hmmprobs(i1, i2, 2) = log(d);
	
		}				
	}
	
		
	for(i1 = 0; i1 < nbases; i1++)
	{
		d = (double)(countsxq1[i1]+npseudocount)/(double)(totalxcounts+npseudocountoutput);	
		if(d < dminval) d = dminval;
		
		hmmprobs(i1, 0, 1) = log(d);
		
		//cout << "i1: " << i1 << "," << log(d) << "\n";
	}
	

	for(i1 = 0; i1 < nbases; i1++)
	{
		 d = (double)(countsyq1[i1]+npseudocount)/(double)(totalycounts+npseudocountoutput);
		 if(d < dminval) d = dminval;
		 
		 hmmprobs(i1, 0, 3) = log(d);
		 
		// cout << "i2: " << i1 << "," << log(d) << "\n";
	}	
	
	
	
	//Order 1
		line = "p1:";			
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{						
						n1seq1 = i1;
						n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
									
						n1seq2 = i3;
						n2seq2 = ConvertOrder1StringToIndex2(i3, i4);
					
						probso2[n2seq1][n2seq2] = (double)(countso2[n2seq1][n2seq2]+npseudocount)/(double)(countso2o1[n1seq1][n1seq2]+npseudocountalign);
									
						if(probso2[n2seq1][n2seq2] == 0) probso2[n2seq1][n2seq2] = dminval;					
					
						if(countso2o1[n1seq1][n1seq2] >= nmincount) hmmprobsho1(n2seq1, n2seq2, 2) = log(probso2[n2seq1][n2seq2]);	
						else hmmprobsho1(n2seq1, n2seq2, 2) = 1;
						
						//hmmprobsho1(n2seq1, n2seq2, 2) = 1;
						
																											
					}
				}						
			}
		}
		//ReadAlignProbs(line, hmmprobsho, norder, 2, noutputsho);
		
		line = "q1x:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				n1seq1 = i1;
				n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
						
				d = (double)(countsxq2[n2seq1]+npseudocount)/(double)(countsxq2q1[n1seq1]+npseudocountoutput);
				if(d < dminval) d = dminval;	
				
				if(countsxq2q1[n1seq1] >= nmincount) hmmprobsho1(n2seq1, 0, 1) = log(d);
				else hmmprobsho1(n2seq1, 0, 1) = 1;
				
				//hmmprobsho1(n2seq1, 0, 1) = 1;
			}
		}
		
	
		line = "q1y:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				n1seq1 = i1;
				n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
									
				d = (double)(countsyq2[n2seq1]+npseudocount)/(double)(countsyq2q1[n1seq1]+npseudocountoutput);
				if(d < dminval) d = dminval;
				
				if(countsyq2q1[n1seq1] >= nmincount) hmmprobsho1(n2seq1, 0, 3) = log(d);
				else hmmprobsho1(n2seq1, 0, 3) = 1;
				
				//hmmprobsho1(n2seq1, 0, 3) = 1;
			}
		}
		

	
	
	//order 2

		line ="p2:";			
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						for(i5 = 0; i5 < nbases; i5++)
						{
							for(i6 = 0; i6 < nbases; i6++)
							{								
								n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
								n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
									
								n2seq2 = ConvertOrder1StringToIndex2(i4, i5);
								n3seq2 = ConvertOrder2StringToIndex2(i4, i5, i6);
					
								probso3[n3seq1][n3seq2] = (double)(countso3[n3seq1][n3seq2]+npseudocount)/(double)(countso3o2[n2seq1][n2seq2]+npseudocountalign);
								if(probso3[n3seq1][n3seq2] < dminval) probso3[n3seq1][n3seq2] = dminval;
								
								if(countso3o2[n2seq1][n2seq2] >= nmincount) hmmprobsho2(n3seq1,n3seq2,2) = log(probso3[n3seq1][n3seq2]);
								else 
								{
									if(countso2o1[i2][i5] > nmincount) hmmprobsho2(n3seq1,n3seq2,2) = 2;
									else hmmprobsho2(n3seq1,n3seq2,2) = 1;																													
								}
								
								if(i1 == 1 && i2 == 1 && i3 == 1 && i4 == 1 && i5 == 1 && i6 == 1)
								{
									cout << countso3[n3seq1][n3seq2] << "\n";
									cout << countso3o2[n2seq1][n2seq2] << "," <<  nmincount << "\n";
									cout << probso3[n3seq1][n3seq2] << "," << hmmprobsho2(n3seq1,n3seq2,2) << "\n";
									
									cout << "\n";
									
									n3seq1 = ConvertOrder2StringToIndex2(i1, i2, 0);
									n3seq2 = ConvertOrder2StringToIndex2(i4, i5, 0);
									cout << "A:" << countso3[n3seq1][n3seq2] << "\n";
									
									n3seq1 = ConvertOrder2StringToIndex2(i1, i2, 1);
									n3seq2 = ConvertOrder2StringToIndex2(i4, i5, 1);
									cout << "C:" << countso3[n3seq1][n3seq2] << "\n";
									
									n3seq1 = ConvertOrder2StringToIndex2(i1, i2, 2);
									n3seq2 = ConvertOrder2StringToIndex2(i4, i5, 2);
									cout << "G:" << countso3[n3seq1][n3seq2] << "\n";
									
									n3seq1 = ConvertOrder2StringToIndex2(i1, i2, 3);
									n3seq2 = ConvertOrder2StringToIndex2(i4, i5, 4);
									cout << "T:" << countso3[n3seq1][n3seq2] << "\n";
									
									cout << "total counts:\n";
									cout << countso3o2[n2seq1][n2seq2] << "\n";
								}
								
																		
								//hmmprobsho2(n3seq1,n3seq2,2) = 1;
							}									
						}
					}
				}						
			}
		}
		//ReadAlignProbs(line, hmmprobsho, norder, 2, noutputsho);
	
		line = "q2x:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
					n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
						
					d = (double)(countsxq3[n3seq1]+npseudocount)/(double)(countsxq3q2[n2seq1]+npseudocountoutput);
					if(d < dminval) d = dminval;	
					
					if(countsxq3q2[n2seq1] >= nmincount) hmmprobsho2(n3seq1,0,1) = log(d);
					else
					{
						if(countsxq2q1[i2] > nmincount) hmmprobsho2(n3seq1,0,1) = 2;
						else hmmprobsho2(n3seq1,0,1) = 1;	
					}
					
					//hmmprobsho2(n3seq1,0,1) = 1;	

				}			
			}
		}		

		line = "q2y:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
					n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
									
					d = (double)(countsyq3[n3seq1]+npseudocount)/(double)(countsyq3q2[n2seq1]+npseudocountoutput);
					if(d < dminval) d = dminval;		
					
					if(countsyq3q2[n2seq1] >= nmincount) hmmprobsho2(n3seq1,0,3) = log(d);
					else 
					{
						if(countsyq2q1[i2] > nmincount) hmmprobsho2(n3seq1,0,3) = 2;
						else hmmprobsho2(n3seq1,0,3) = 1;
					}
					
					//hmmprobsho2(n3seq1,0,3) = 1;
				}			
			}
		}			
	
	
	
	
	
	
	
	//order 3
	
	line = "p3:";	
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						for(i5 = 0; i5 < nbases; i5++)
						{
							for(i6 = 0; i6 < nbases; i6++)
							{								
								for(i7 = 0; i7 < nbases; i7++)
								{
									for(i8 = 0; i8 < nbases; i8++)		
									{
										n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
										n4seq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
									
										n3seq2 = ConvertOrder2StringToIndex2(i5, i6, i7);
										n4seq2 = ConvertOrder3StringToIndex2(i5, i6, i7, i8);
					
										probso4[n4seq1][n4seq2] = (double)(countso4[n4seq1][n4seq2]+npseudocount)/(double)(countso4o3[n3seq1][n3seq2]+npseudocountalign);									
										if(probso4[n4seq1][n4seq2] < dminval) probso4[n4seq1][n4seq2] = dminval;									
									
										//hmmprobsho(n4seq1,n4seq2,2) = log(probso4[n4seq1][n4seq2]);							
										
										if(countso4o3[n3seq1][n3seq2] >= nmincount) hmmprobsho3(n4seq1,n4seq2,2) = log(probso4[n4seq1][n4seq2]);
										else 
										{
											n2seq1 = ConvertOrder1StringToIndex2(i2, i3);																
											n2seq2 = ConvertOrder1StringToIndex2(i6, i7);
								
											if(countso3o2[n2seq1][n2seq2] > nmincount) hmmprobsho3(n4seq1,n4seq2,2) = 3;
											else if(countso2o1[i3][i7] > nmincount) hmmprobsho3(n4seq1,n4seq2,2) = 2;
											else hmmprobsho3(n4seq1,n4seq2,2) = 1;																													
										}
								
									}
																	
								}								
							}													
						}
					}
				}												
			}
		}
		//ReadAlignProbs(line, hmmprobsho, norder, 2, noutputsho);
		
		line = "q3x:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
						n4seq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
						
						d = (double)(countsxq4[n4seq1]+npseudocount)/(double)(countsxq4q3[n3seq1]+npseudocountoutput);
						if(d < dminval) d = dminval;								
						
						//hmmprobsho(n4seq1,0,1) = log(d);									
						
						if(countsxq4q3[n3seq1] >= nmincount) hmmprobsho3(n4seq1,0,1) = log(d);
						else
						{	
							n2seq1 = ConvertOrder1StringToIndex2(i2, i3);
							
							if(countsxq3q2[n2seq1] > nmincount) hmmprobsho3(n4seq1,0,1) = 3;
							else if(countsxq2q1[i3] > nmincount) hmmprobsho3(n4seq1,0,1) = 2;
							else hmmprobsho3(n4seq1,0,1) = 1;	
						}
					
					}				
				}			
			}
		}
		
		//ReadOutputProbs(line, hmmprobsho, norder, 1, noutputsho);
		
		line = "q3y:";
		for(i1 = 0; i1 < nbases; i1++)
		{	
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
						n4seq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
									
						d = (double)(countsyq4[n4seq1]+npseudocount)/(double)(countsyq4q3[n3seq1]+npseudocountoutput);
						if(d < dminval) d = dminval;				
					
						//hmmprobsho(n4seq1,0,3) = log(d);						
						
						if(countsyq4q3[n3seq1] >= nmincount) hmmprobsho3(n4seq1,0,3) = log(d);
						else
						{	
							n2seq1 = ConvertOrder1StringToIndex2(i2, i3);
							
							if(countsyq3q2[nseq1] > nmincount) hmmprobsho3(n4seq1,0,3) = 3;
							else if(countsyq2q1[i3] > nmincount) hmmprobsho3(n4seq1,0,3) = 2;
							else hmmprobsho3(n4seq1,0,3) = 1;	
						}
					}
				}			
			}
		}
		
	
	
	
	finref.close();
	finread.close();
		
	
	delete[] countsxq2q1;
	delete[] countsyq2q1;
	delete[] countsxq3q2;
	delete[] countsyq3q2;
	delete[] countsxq4q3;
	delete[] countsyq4q3;
	
	delete[] countsxq1;
	delete[] countsyq1;
	delete[] countsxq2;
	delete[] countsyq2;
	delete[] countsxq3;
	delete[] countsyq3;
	delete[] countsxq4;
	delete[] countsyq4;
		
	Clear2DimArray(countso2o1, nlengtho1);	
	Clear2DimArray(countso3o2, nlengtho2);
	Clear2DimArray(countso4o3, nlengtho3);
	Clear2DimArray(countso1, nlengtho1);
	Clear2DimArray(countso2, nlengtho2);	
	Clear2DimArray(countso3, nlengtho3);
	Clear2DimArray(countso4, nlengtho4);
	Clear2DimArray(probso1, nlengtho1);
	Clear2DimArray(probso2, nlengtho2);
	Clear2DimArray(probso3, nlengtho3);
	Clear2DimArray(probso4, nlengtho4);
	
	return 1;
}





	
int CreateProbMatrixFromTrainData(string sfilein1, string sfilein2, double* hmmprobs, double* hmmprobsho, 
								   int norder, int noutputs, int noutputsho)
{	
	int i;
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	
	int nbases;	
	int totalcounts;
	int totalxcounts;
	int totalycounts;
	int n1seq1;
	int n1seq2;
	int n2seq1;
	int n2seq2;
	int n3seq1;
	int n3seq2;
	int n4seq1;
	int n4seq2;
	int nlengtho1;
	int nlengtho2;
	int nlengtho3;
	int nlengtho4;
	
	int nseq1;
	int nseq2;
	int nSeqLength1;
	int nSeqLength2;
	
	double d;
	
	string line;
	string s1;
	string s2;
	string sref;
	string sread;
	
	int* seq1;
	int* seq2;
	
	
	ifstream finref(sfilein1.c_str());
	ifstream finread(sfilein2.c_str());
	
	if(finref.good() == 0)
	{
		cout << "Problem accessing " << sfilein1 << ". Please check file exists.\n\n";
		return 0;
	}
	else if(finread.good() == 0)
	{
		cout << "Problem accessing " << sfilein2 << ". Please check file exists.\n\n";
		return 0;
	}
	
	
	int nordermax = 4;
	
	nlengtho1 = nordermax;
	nlengtho2 = (int)pow(4, nordermax-2);
	nlengtho3 = (int)pow(4, nordermax-1);
	nlengtho4 = (int)pow(4, nordermax);
	
	

	//1st order
	int** countso1 = CreateIntArray(nlengtho1, nlengtho1);
	int* countsxq1 = new int[nlengtho1];
	int* countsyq1 = new int[nlengtho1];

	//2nd order
	int** countso2o1 = CreateIntArray(nlengtho1, nlengtho1);
	int* countsxq2q1 = new int[nlengtho1];
	int* countsyq2q1 = new int[nlengtho1];
	
	int** countso2 = CreateIntArray(nlengtho2, nlengtho2);	
	int* countsxq2 = new int[nlengtho2];
	int* countsyq2 = new int[nlengtho2];


	//3rd order
	int** countso3o2 = CreateIntArray(nlengtho2, nlengtho2);
	int* countsxq3q2 = new int[nlengtho2];
	int* countsyq3q2 = new int[nlengtho2];
	
	int** countso3 = CreateIntArray(nlengtho3, nlengtho3);	
	int* countsxq3 = new int[nlengtho3];
	int* countsyq3 = new int[nlengtho3];
		
	
	//4th order
	int** countso4o3 = CreateIntArray(nlengtho3, nlengtho3);
	int* countsxq4q3 = new int[nlengtho3];
	int* countsyq4q3 = new int[nlengtho3];
		
	int** countso4 = CreateIntArray(nlengtho4, nlengtho4);
	int* countsxq4 = new int[nlengtho4];
	int* countsyq4 = new int[nlengtho4];
	
	double** probso1 = CreateDoubleArray(nlengtho1, nlengtho1);
	double** probso2 = CreateDoubleArray(nlengtho2, nlengtho2);
	double** probso3 = CreateDoubleArray(nlengtho3, nlengtho3);
	double** probso4 = CreateDoubleArray(nlengtho4, nlengtho4);
	
	double dminval;
	
	//dminval = exp(-12.3856166667); //0.0000000001; //0.00001;
	dminval = 0.0000000001; //0.00001;

	
	//cout << log(dminval) << "\n";
	
	for(i = 0; i < nlengtho1; i++) 
	{
		countsxq1[i] = 0;
		countsyq1[i] = 0;
		
		countsxq2q1[i] = 0;
		countsyq2q1[i] = 0;
	}
	
	for(i = 0; i < nlengtho2; i++) 
	{
		countsxq3q2[i] = 0;
		countsyq3q2[i] = 0;
		
		countsxq2[i] = 0;
		countsyq2[i] = 0;
	}
	
	for(i = 0; i < nlengtho3; i++) 
	{
		countsxq4q3[i] = 0;
		countsyq4q3[i] = 0;
		
		countsxq3[i] = 0;
		countsyq3[i] = 0;
	}
	
	for(i = 0; i < nlengtho4; i++) 
	{
		countsxq4[i] = 0;
		countsyq4[i] = 0;
	}
	
	totalcounts = 0;
	totalxcounts = 0;
	totalycounts = 0;
	
	sref = "";
	sread = "";
	getline(finref, line);	
	getline(finread, line);	
	
	
	while(finref.eof() == 0)
	{
		getline(finref, line);		
		sref += line;
		
		getline(finread, line);
		sread += line;
	}
	
	
	nSeqLength1 = sref.length();
	nSeqLength2 = sread.length();
			
	seq1 = new int[nSeqLength1];
	seq2 = new int[nSeqLength2];
	
	ConvertStringToInts2(sref, seq1);
	ConvertStringToInts2(sread, seq2);
			
	for(i = 0; i < nSeqLength1; i++)
	{
			
		if(i >= 1)
		{
			n1seq1 = ConvertOrderStringToIndex2(seq1, i-1, 1);
			n2seq1 = ConvertOrderStringToIndex2(seq1, i, 2);
				
			n1seq2 = ConvertOrderStringToIndex2(seq2, i-1, 1);
			n2seq2 = ConvertOrderStringToIndex2(seq2, i, 2);
				
			if(n1seq1 != -1 && n2seq1 != -1)
			{
				countsxq2q1[n1seq1]++;
				countsxq2[n2seq1]++;
			}
				
			if(n1seq2 != -1 && n2seq2 != -1)
			{
				countsyq2q1[n1seq2]++;
				countsyq2[n2seq2]++;
			}
				
			if(n1seq1 != -1 && n2seq1 != -1 && n1seq2 != -1 && n2seq2 != -1)
			{
				countso2o1[n1seq1][n1seq2]++;
				countso2[n2seq1][n2seq2]++;					
			}
		}
			
		if(i >= 2)
		{
			n2seq1 = ConvertOrderStringToIndex2(seq1, i-1, 2);
			n3seq1 = ConvertOrderStringToIndex2(seq1, i, 3);
				
			n2seq2 = ConvertOrderStringToIndex2(seq2, i-1, 2);
			n3seq2 = ConvertOrderStringToIndex2(seq2, i, 3);
				
			if(n2seq1 != -1 && n3seq1 != -1)
			{
				countsxq3q2[n2seq1]++;
				countsxq3[n3seq1]++;
			}
				
			if(n2seq2 != -1 && n3seq2 != -1)
			{
				countsyq3q2[n2seq2]++;
				countsyq3[n3seq2]++;
			}
				
			if(n2seq1 != -1 && n3seq1 != -1 && n2seq2 != -1 && n3seq2 != -1)
			{
				countso3o2[n2seq1][n2seq2]++;
				countso3[n3seq1][n3seq2]++;					
			}
		}
			
		if(i >= 3)
		{
			n3seq1 = ConvertOrder2StringToIndex2(seq1, i-1);
			n4seq1 = ConvertOrder3StringToIndex2(seq1, i);
				
			n3seq2 = ConvertOrder2StringToIndex2(seq2, i-1);
			n4seq2 = ConvertOrder3StringToIndex2(seq2, i);
			
			if(n3seq1 != -1 && n4seq1 != -1)
			{
				countsxq4q3[n3seq1]++;
				countsxq4[n4seq1]++;
			}
					
			if(n3seq2 != -1 && n4seq2 != -1)
			{
				countsyq4q3[n3seq2]++;
				countsyq4[n4seq2]++;
			}
				
			if(n3seq1 != -1 && n4seq1 != -1 && n3seq2 != -1 && n4seq2 != -1)
			{
				countso4o3[n3seq1][n3seq2]++;
				countso4[n4seq1][n4seq2]++;					
			}
		}
			
		nseq1 = seq1[i];
		nseq2 = seq2[i];
			
		if(nseq1 != -1 && nseq2 != -1)
		{
			countso1[nseq1][nseq2]++;
			totalcounts++;
		}
			
		if(nseq1 != -1)
		{
			countsxq1[nseq1]++;
			totalxcounts++;
		}
			
		if(nseq2 != -1)
		{
			countsyq1[nseq2]++;
			totalycounts++;
		}
		
	}
	
	delete[] seq1;
	delete[] seq2;
	
	int npseudocount = 1;
	int npseudocountalign = 16;
	int npseudocountoutput = 4;

	nbases = 4;
	
	for(i1 = 0; i1 < nbases; i1++)
	{
		for(i2 = 0; i2 < nbases; i2++)
		{
			d = (double)(countso1[i1][i2]+npseudocount)/(double)(totalcounts+npseudocountalign);	
			
			hmmprobs(i1, i2, 2) = log(d);
			
		}				
	}
		
	for(i1 = 0; i1 < nbases; i1++)
	{
		d = (double)(countsxq1[i1]+npseudocount)/(double)(totalxcounts+npseudocountoutput);	
		
		hmmprobs(i1, 0, 1) = log(d);
	}			

	for(i1 = 0; i1 < nbases; i1++)
	{
		 d = (double)(countsyq1[i1]+npseudocount)/(double)(totalycounts+npseudocountoutput);
		 
		 hmmprobs(i1, 0, 3) = log(d);
	}	
	
	if(norder == 1)
	{
		line = "p1:";			
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{						
						n1seq1 = i1;
						n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
									
						n1seq2 = i3;
						n2seq2 = ConvertOrder1StringToIndex2(i3, i4);
					
						probso2[n2seq1][n2seq2] = (double)(countso2[n2seq1][n2seq2]+npseudocount)/(double)(countso2o1[n1seq1][n1seq2]+npseudocountalign);				
										
						hmmprobsho(n2seq1, n2seq2, 2) = log(probso2[n2seq1][n2seq2]);																						
					}
				}						
			}
		}

		line = "q1x:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				n1seq1 = i1;
				n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
						
				d = (double)(countsxq2[n2seq1]+npseudocount)/(double)(countsxq2q1[n1seq1]+npseudocountoutput);
				
				hmmprobsho(n2seq1, 0, 1) = log(d);
			}
		}
		
	
		line = "q1y:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				n1seq1 = i1;
				n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
									
				d = (double)(countsyq2[n2seq1]+npseudocount)/(double)(countsyq2q1[n1seq1]+npseudocountoutput);
				
				hmmprobsho(n2seq1, 0, 3) = log(d);																					
			}
		}
	}
	
	
	if(norder == 2)
	{
		line ="p2:";			
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						for(i5 = 0; i5 < nbases; i5++)
						{
							for(i6 = 0; i6 < nbases; i6++)
							{								
								n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
								n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
									
								n2seq2 = ConvertOrder1StringToIndex2(i4, i5);
								n3seq2 = ConvertOrder2StringToIndex2(i4, i5, i6);
								
								probso3[n3seq1][n3seq2] = (double)(countso3[n3seq1][n3seq2]+npseudocount)/(double)(countso3o2[n2seq1][n2seq2]+npseudocountalign);
								
								hmmprobsho(n3seq1,n3seq2,2) = log(probso3[n3seq1][n3seq2]);
																		
							}									
						}
					}
				}						
			}
		}
	
		line = "q2x:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
					n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
					
					d = (double)(countsxq3[n3seq1]+npseudocount)/(double)(countsxq3q2[n2seq1]+npseudocountoutput);
					
					hmmprobsho(n3seq1,0,1) = log(d);							

				}			
			}
		}		

		line = "q2y:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					n2seq1 = ConvertOrder1StringToIndex2(i1, i2);
					n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
					
					d = (double)(countsyq3[n3seq1]+npseudocount)/(double)(countsyq3q2[n2seq1]+npseudocountoutput);
					
					hmmprobsho(n3seq1,0,3) = log(d);
				}			
			}
		}	
	}
	
	
	
	
	
	if(norder == 3)
	{
		
		line = "p3:";	
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						for(i5 = 0; i5 < nbases; i5++)
						{
							for(i6 = 0; i6 < nbases; i6++)
							{								
								for(i7 = 0; i7 < nbases; i7++)
								{
									for(i8 = 0; i8 < nbases; i8++)		
									{
										n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
										n4seq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
									
										n3seq2 = ConvertOrder2StringToIndex2(i5, i6, i7);
										n4seq2 = ConvertOrder3StringToIndex2(i5, i6, i7, i8);
					
										probso4[n4seq1][n4seq2] = (double)(countso4[n4seq1][n4seq2]+npseudocount)/(double)(countso4o3[n3seq1][n3seq2]+npseudocountalign);									
										if(probso4[n4seq1][n4seq2] < dminval) probso4[n4seq1][n4seq2] = dminval;									
									
										hmmprobsho(n4seq1,n4seq2,2) = log(probso4[n4seq1][n4seq2]);							
									}
																	
								}								
							}													
						}
					}
				}												
			}
		}
		
		line = "q3x:";
		for(i1 = 0; i1 < nbases; i1++)
		{
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
						n4seq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
						
						d = (double)(countsxq4[n4seq1]+npseudocount)/(double)(countsxq4q3[n3seq1]+npseudocountoutput);
						if(d < dminval) d = dminval;								
						
						hmmprobsho(n4seq1,0,1) = log(d);									
					}				
				}			
			}
		}
		
		
		line = "q3y:";
		for(i1 = 0; i1 < nbases; i1++)
		{	
			for(i2 = 0; i2 < nbases; i2++)
			{
				for(i3 = 0; i3 < nbases; i3++)
				{
					for(i4 = 0; i4 < nbases; i4++)
					{
						n3seq1 = ConvertOrder2StringToIndex2(i1, i2, i3);
						n4seq1 = ConvertOrder3StringToIndex2(i1, i2, i3, i4);
									
						d = (double)(countsyq4[n4seq1]+npseudocount)/(double)(countsyq4q3[n3seq1]+npseudocountoutput);
						if(d < dminval) d = dminval;				
					
						hmmprobsho(n4seq1,0,3) = log(d);						
					}
				}			
			}
		}
	
	}
	

	finref.close();
	finread.close();
		
	
	delete[] countsxq2q1;
	delete[] countsyq2q1;
	delete[] countsxq3q2;
	delete[] countsyq3q2;
	delete[] countsxq4q3;
	delete[] countsyq4q3;
	
	delete[] countsxq1;
	delete[] countsyq1;
	delete[] countsxq2;
	delete[] countsyq2;
	delete[] countsxq3;
	delete[] countsyq3;
	delete[] countsxq4;
	delete[] countsyq4;
		
	Clear2DimArray(countso2o1, nlengtho1);	
	Clear2DimArray(countso3o2, nlengtho2);
	Clear2DimArray(countso4o3, nlengtho3);
	Clear2DimArray(countso1, nlengtho1);
	Clear2DimArray(countso2, nlengtho2);	
	Clear2DimArray(countso3, nlengtho3);
	Clear2DimArray(countso4, nlengtho4);
	Clear2DimArray(probso1, nlengtho1);
	Clear2DimArray(probso2, nlengtho2);
	Clear2DimArray(probso3, nlengtho3);
	Clear2DimArray(probso4, nlengtho4);
	
	return 1;
}










int InitializeHMMVariables(int nstates, int norder, int noutputs, int noutputsho, int  noutputsho1, int noutputsho2, 
							int noutputsho3, int bVariableOrder, double dtau, int* statetypes, double* hmmesum, 
							double* hmmprobs, double* hmmprobsho, double* hmmprobsho1, double* hmmprobsho2, 
							double* hmmprobsho3, double* as[], string sinputfiletrain1b, string sinputfiletrain2b)
{
	int i;
	int j;
	int k;
	int l;
	int binit;
	int statetypefrom;
	int statetypeto;
	
	double d;
	
	int* nconnectioncounts = new int[nstates];
	
	for(i = 0; i < noutputs; i++)	
		for(j = 0; j < noutputs; j++)	
			for(k = 0; k < nstates; k++) hmmesum(i,j,k) = 0;

	
	for(i = 0; i < nstates; i++)	
	{			
		statetypefrom = statetypes[i];
				
		for(j = 0; j < nstates; j++)	
		{			
			statetypeto = statetypes[j];
			
			if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEX) nconnectioncounts[i]++;			
			else if(statetypefrom == BEGINSTATE && statetypeto == MATCHSTATE) nconnectioncounts[i]++;			
			else if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEY) nconnectioncounts[i]++;	
			else if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEXHO) nconnectioncounts[i]++;			
			else if(statetypefrom == BEGINSTATE && statetypeto == MATCHSTATEHO) nconnectioncounts[i]++;			
			else if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEYHO) nconnectioncounts[i]++;					
			else if(statetypefrom == BEGINSTATE && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++; 
			else if(statetypefrom == BEGINSTATE && statetypeto == ENDSTATE)
			{
			}			
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == OUTPUTSTATEX) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == MATCHSTATE) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == ENDSTATE)
			{
			}
			else if(statetypefrom == OUTPUTSTATEXHO && statetypeto == OUTPUTSTATEXHO) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEXHO && statetypeto == MATCHSTATEHO) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEXHO && statetypeto == ENDSTATE)
			{
			}
			else if(statetypefrom == MATCHSTATE && statetypeto == MATCHSTATE) nconnectioncounts[i]++;			
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEX) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEY) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATE && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEXPH) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEYPH) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATE && statetypeto == ENDSTATE)	
			{				
			}			
			else if(statetypefrom == MATCHSTATEHO && statetypeto == MATCHSTATEHO) nconnectioncounts[i]++;			
			else if(statetypefrom == MATCHSTATEHO && statetypeto == OUTPUTSTATEXHO) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATEHO && statetypeto == OUTPUTSTATEYHO) nconnectioncounts[i]++;			
			else if(statetypefrom == MATCHSTATEHO && statetypeto == ENDSTATE)	
			{				
			}			
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == OUTPUTSTATEY) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == MATCHSTATE) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == ENDSTATE)
			{
			}
			else if(statetypefrom == OUTPUTSTATEYHO && statetypeto == OUTPUTSTATEYHO) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEYHO && statetypeto == MATCHSTATEHO) nconnectioncounts[i]++;			
			else if(statetypefrom == OUTPUTSTATEYHO && statetypeto == ENDSTATE)
			{
			}
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == OUTPUTSTATEXPH) nconnectioncounts[i]++;			
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == MATCHSTATE) nconnectioncounts[i]++;
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == ENDSTATE) 
			{
			}
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == OUTPUTSTATEYPH) nconnectioncounts[i]++;	
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++;		
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == MATCHSTATE)  nconnectioncounts[i]++;	
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == ENDSTATE) 
			{
			}			
			else if(statetypefrom == MATCHSTATEPH && statetypeto == MATCHSTATE) nconnectioncounts[i]++; 
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEX) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEY) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEXPH) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEYPH) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == MATCHSTATEPH) nconnectioncounts[i]++;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == ENDSTATE) 
			{
			}   							
			 
		}
	}
	

		
	for(i = 0; i < nstates; i++)	
	{			
		statetypefrom = statetypes[i];
				
		for(j = 0; j < nstates; j++)	
		{
			as[i][j] = AMIN;									
			
			statetypeto = statetypes[j];
			
			d = log(1.0/(double)nconnectioncounts[i]);

			if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEX) as[i][j] = d;
			else if(statetypefrom == BEGINSTATE && statetypeto == MATCHSTATE) as[i][j] = d;
			else if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEY) as[i][j] = d;			
			else if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEXHO) as[i][j] = d;
			else if(statetypefrom == BEGINSTATE && statetypeto == MATCHSTATEHO) as[i][j] = d;
			else if(statetypefrom == BEGINSTATE && statetypeto == OUTPUTSTATEYHO) as[i][j] = d;			
			else if(statetypefrom == BEGINSTATE && statetypeto == ENDSTATE) as[i][j] = dtau;
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == OUTPUTSTATEX) as[i][j] = d;			
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == MATCHSTATE) as[i][j] = d;			
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == MATCHSTATEPH) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == ENDSTATE) as[i][j] = dtau;			
			else if(statetypefrom == OUTPUTSTATEXHO && statetypeto == OUTPUTSTATEXHO) as[i][j] = d;			
			else if(statetypefrom == OUTPUTSTATEXHO && statetypeto == MATCHSTATEHO) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEXHO && statetypeto == ENDSTATE) as[i][j] = dtau;
			else if(statetypefrom == MATCHSTATE && statetypeto == MATCHSTATE) as[i][j] = d; 			   
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEX) as[i][j] = d;			
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEY) as[i][j] = d;				
			else if(statetypefrom == MATCHSTATEHO && statetypeto == MATCHSTATEHO) as[i][j] = d; 			   
			else if(statetypefrom == MATCHSTATEHO && statetypeto == OUTPUTSTATEXHO) as[i][j] = d;			
			else if(statetypefrom == MATCHSTATEHO && statetypeto == OUTPUTSTATEYHO) as[i][j] = d;				
			else if(statetypefrom == MATCHSTATEHO && statetypeto == ENDSTATE)	as[i][j] = dtau;
			else if(statetypefrom == MATCHSTATE && statetypeto == MATCHSTATEPH) as[i][j] = d;
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEXPH) as[i][j] = d;			
			else if(statetypefrom == MATCHSTATE && statetypeto == OUTPUTSTATEYPH) as[i][j] = d;			
			else if(statetypefrom == MATCHSTATE && statetypeto == ENDSTATE)	as[i][j] = dtau;		
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == OUTPUTSTATEY) as[i][j] = d;					
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == MATCHSTATE) as[i][j] = d;	
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == MATCHSTATEPH) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == ENDSTATE) as[i][j] = dtau;
			else if(statetypefrom == OUTPUTSTATEYHO && statetypeto == OUTPUTSTATEYHO) as[i][j] = d;					
			else if(statetypefrom == OUTPUTSTATEYHO && statetypeto == MATCHSTATEHO) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEYHO && statetypeto == ENDSTATE) as[i][j] = dtau;
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == OUTPUTSTATEXPH) as[i][j] = d; 
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == ENDSTATE) as[i][j] = dtau;
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == MATCHSTATE) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEXPH && statetypeto == MATCHSTATEPH) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == OUTPUTSTATEYPH) as[i][j] = d; 
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == ENDSTATE) as[i][j] = dtau;
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == MATCHSTATE) as[i][j] = d;
			else if(statetypefrom == OUTPUTSTATEYPH && statetypeto == MATCHSTATEPH) as[i][j] = d;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == ENDSTATE) as[i][j] = dtau;		   				
			else if(statetypefrom == MATCHSTATEPH && statetypeto == MATCHSTATE) as[i][j] = d;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEX) as[i][j] = d;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEY) as[i][j] = d;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEXPH) as[i][j] = d;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == OUTPUTSTATEYPH) as[i][j] = d;
			else if(statetypefrom == MATCHSTATEPH && statetypeto == MATCHSTATEPH) as[i][j] = d;
			else if(statetypefrom == BEGINSTATE && statetypeto == MATCHSTATEPH) as[i][j] = d; 		
			else if(statetypefrom == MATCHSTATE && statetypeto == MATCHSTATEPH) as[i][j] = d; 
			else if(statetypefrom == OUTPUTSTATEX && statetypeto == MATCHSTATEPH) as[i][j] = d; 
			else if(statetypefrom == OUTPUTSTATEY && statetypeto == MATCHSTATEPH) as[i][j] = d; 
			
		}
	}	
	
	delete[] nconnectioncounts;
	
	if(bVariableOrder == 0) 
	{
		binit = CreateProbMatrixFromTrainData(sinputfiletrain1b, sinputfiletrain2b, hmmprobs, 
									   hmmprobsho, norder, noutputs, noutputsho);
	}
	else
	{
		binit = CreateProbMatrixFromTrainDataVariedOrder(sinputfiletrain1b, sinputfiletrain2b, hmmprobs, 
												 hmmprobsho1, hmmprobsho2, hmmprobsho3, noutputs, 
												 noutputsho1, noutputsho2, noutputsho3);
	
	}
	
	if(binit == 0)
	{
		return 0;
	}
	
	/////////////////////////////
	
	
	cout << "as:\n";
	for(k = 0; k < nstates; k++)
	{
		for(l = 0; l < nstates; l++)
		{
			cout << exp(as[k][l]) << ",";
		}
		cout << "\n";
	}
	cout << "\n";
	
	
	cout << "hmmprobs - match state:\n";
	for(k = 0; k < noutputs; k++)
	{
		for(l = 0; l < noutputs; l++)
		{
			cout << exp(hmmprobs(k,l,2)) << ",";
		}
			
		cout << "\n";
	}
		
	cout << "hmmprobs - q1:\n";
	for(k = 0; k < noutputs; k++)
	{
		for(l = 0; l < noutputs; l++)
		{
			cout << exp(hmmprobs(k,l,1)) << ",";
		}
			
		cout << "\n";
	}
		
	cout << "hmmprobs - q2:\n";
	for(k = 0; k < noutputs; k++)
	{
		for(l = 0; l < noutputs; l++)
		{
			cout << exp(hmmprobs(k,l,3)) << ",";
		}
			
		cout << "\n";
	}
	cout << "\n";
	
	return binit;
	
}

/*
void TrainVariedOrder(int norder, string soutresults, int nminlen)
{	
	
	int i; 
	int j; 
	int k; 
	int l;
	int nbeginstate; 
	int nendstate; 
	int noutputs; 
	int noutputsho; 
	int noutputsho1;
	int noutputsho2;
	int noutputsho3;
	int nstates;
	int statetypefrom, statetypeto;
	int nbases;
	int nSeqLength1;
	int nSeqLength2;	
	int maxcount;
	int mincount;
	int ntrainseq;
	int nmaxlength;	
	int btrain;
	int boutput;
	int bVariableOrder;
	int ntotaltrainseq;
	int ncount;	 
	int boverfit;

	double d;	
	double pallseqtotal;
	double pallseqtotalold;
	double pallseqinitial;
	double pdiff;
	double dtau;			
	double dthresh; 
	double asmin;
	double atemp;
	double dscoremin;
	double dstd;
	double dscoremean;
	
	string line;
	string sinputfiletrain1;
	string sinputfiletrain2;
	string sinputfiletrain1b;
	string sinputfiletrain2b;
	string sinputfiletest;
	string soutputfile1;
	string soutputfile2;	
	string sseq1;
	string sseq2;
	
	vfwdbackpair vf;

	nstates = 5;
	nbeginstate = 0;
	nendstate = nstates-1;	
			
	int* statetypes = new int[nstates];
	

	statetypes[0] = BEGINSTATE;	
	statetypes[1] = OUTPUTSTATEXHO;
	statetypes[2] = MATCHSTATEHO;
	statetypes[3] = OUTPUTSTATEYHO;		
	statetypes[4] = ENDSTATE;		
	
	norder = 3;
	nbases = 4;
	noutputs = 4;		
	
	noutputsho = -1;
	noutputsho1 = (int)pow(4, 2);
	noutputsho2 = (int)pow(4, 3);
	noutputsho3 = (int)pow(4, 4);
	
	
	
	vfwdbackpair vf;
	
	pallseqtotal = 0;
	pallseqtotalold = 0;
	pdiff = 0;
	
	dthresh = fabs(log(0.8)); 
	maxcount = 5; //8; 		//max number of training steps
	mincount = 0; // 2;   	//min number of training steps

	ntrainseq = -1; //304479; //160000;//10000; //200000; //20000; //1000
	nmaxlength = 500; //200
	
	btrain = 1;
	boutput = 0;
	asmin = log(PMIN);
	dtau = log(0.01);	
	bVariableOrder = 0;
	
	sinputfiletrain1 = "lasttrainseqnew1";
	sinputfiletrain2 = "lasttrainseqnew2";	
	sinputfiletrain1b = "lasttrainseqnew1b";
	sinputfiletrain2b = "lasttrainseqnew2b";
	
	
	soutputfile1 = soutresults;
	soutputfile2 = soutresults;
	soutputfile1 += "_out1";
	soutputfile2 += "_out2";
	
	int* seq1;
	int* seq2;	

	double* dscores;
	
	int* seq1ho;
	int* seq2ho;	


	double* hmmprobsho1 = new double[noutputsho1*noutputsho1*nstates];
	double* hmmprobsho2 = new double[noutputsho2*noutputsho2*nstates];
	double* hmmprobsho3 = new double[noutputsho3*noutputsho3*nstates];
	
	double* hmmprobs = new double[noutputs*noutputs*nstates];	
	double* hmmesum = new double[noutputs*noutputs*nstates];
	
	double** asum = CreateDoubleArray(nstates, nstates);
	double** asums = CreateDoubleArray(nstates, nstates);

	double* hmmesums = new double[noutputs*noutputs*nstates];
	double* hmmasums = new double[noutputs*noutputs*nstates];
	
	double* esumsums = new double[nstates];
	double* asumsums = new double[nstates];
	
	double** matchprobs = CreateDoubleArray(noutputs,noutputs);	
	double** as = CreateDoubleArray(nstates, nstates);
	string** alabels = CreateStringArray(nstates, nstates);
	
	
	
		
	InitializeHMMVariables(nstates, norder, noutputs, noutputsho, noutputsho1, noutputsho2,
							noutputsho3, bVariableOrder, dtau, statetypes, hmmesums, hmmprobs, 
							hmmprobsho, as, sinputfiletrain1b, sinputfiletrain2b);
	
	
	for(i = 0; i < 4; i++)
	{
		hmmprobs(i, 0, 1) = qxs[i];	
		hmmprobs(i, 0, 3) = qys[i];	
	}
	
	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++) hmmprobs(i, j, 2) = ps2[i][j];	
	}
		

	/////////////////////////////////////////////////////////////
	//Start of training routines
	/////////////////////////////////////////////////////////////

	ifstream infiletrain1(sinputfiletrain1.c_str());	
	ifstream infiletrain2(sinputfiletrain2.c_str());	
		
	ncount = 0;
	ntotaltrainseq = 0;

	while(infiletrain1.eof() == 0)
	{
		getline(infiletrain1,line);
		
		if(line.length() > 0) ntotaltrainseq++;
	}
	
	
	infiletrain1.close();
	infiletrain2.close();
		
	
	int ncountseq;
		
	if(ntrainseq == -1) ntrainseq = ntotaltrainseq;
	double* dscores = new double[ntrainseq];
	
	double dscoremin;
	
	//ofs.open("test.txt", std::ofstream::out | std::ofstream::trunc);
	
	ofstream out1;
	ofstream out2;

	while(btrain == 1)
	{					
		//infiletrain1.seekg(0);
		//infiletrain2.seekg(0);
		
		infiletrain1.open(sinputfiletrain1.c_str());	
		infiletrain2.open(sinputfiletrain2.c_str());
		
		out1.open(soutputfile1.c_str(), std::ofstream::out | std::ofstream::trunc);
		out2.open(soutputfile2.c_str(), std::ofstream::out | std::ofstream::trunc);

		//infiletrain.open("matches_gappedb_chr1_chr1", ifstream::in); 	
		
		//cout << "count = " << ncount << "\n";
		
		/////////////////Initialization
		
		pallseqtotal = 0;
		
		dscoremin = 0;
		
		//Initialize E and A sums (E represents the probability output matrix and A the transfer matrix)
		//Esum contains the total number of outputs for each output type in each state
		//Esums sums Esum across all training sequences
		//Asum contains the total number of transfers from one state to another		
		//Asums sums Asum across all training sequences		
		for(k = 0; k < nstates; k++)
		{	
			esumsums[k] = 0;	//Used to sum across outputs for esums
			asumsums[k] = 0;	//Used to sum across states for asums
					
			for(j = 0; j < noutputs; j++)
			{
				for(l = 0; l < noutputs; l++) hmmesums(j,l,k) = 0;						
			}			
			
			for(l = 0; l < nstates; l++) asums[k][l] = 0;			
		}
	
		////////////////Begin training
	

		i = 0;
		
		ncountseq = 0;

		//while(infiletrain.eof() == 0 && i < ntrainseq)
		
		while(infiletrain1.eof() == 0 && i < ntrainseq) //&& i < 10000) // && i < ntrainseq)
		{					
			//Asum and Esum are summed across each training sequence so initialize at start of each training sequence
			//as opposed to esums and asums which are across all the training sequences
			for(k = 0; k < nstates; k++)
			{	
				for(j = 0; j < noutputs; j++)
				{
					for(l = 0; l < noutputs; l++) hmmesum(j,l,k) = 0;				
				}			
			
				for(l = 0; l < nstates; l++)  asum[k][l] = 0;		
							
			}
		
			
			getline(infiletrain1,line);
			sseq1 = line;		
			getline(infiletrain2,line);
			sseq2 = line;
			
			nSeqLength1 = sseq1.length();
			nSeqLength2 = sseq2.length();
			
			seq1 = new int[nSeqLength1];
			seq2 = new int[nSeqLength2];
	
			ConvertStringToInts(sseq1, seq1, nbases);
			ConvertStringToInts(sseq2, seq2, nbases);
		

			if(nSeqLength1 > nminlen && nSeqLength2 > nminlen && nSeqLength1 < nmaxlength && nSeqLength2 < nmaxlength) // && abs(nSeqLength1-nSeqLength2) <= 0)
			//if(nSeqLength1 > 40 && nSeqLength2 > 40 && nSeqLength1 < nmaxlength && nSeqLength2 < nmaxlength) // && abs(nSeqLength1-nSeqLength2) <= 0)
			{											 
						
				
				
	
				vf = ViterbiTrainVariedOrderHMM(seq1, seq2, sseq1, sseq2, nSeqLength1, nSeqLength2,
	 					    noutputs, nstates, hmmprobs, hmmprobsho1, hmmprobsho2, hmmprobsho3, as, statetypes, hmmesum, asum,						   
						    nbeginstate, nendstate, qxs, qys,
						    ps2, norder, noutputsho1, noutputsho2, noutputsho3);
	
				out1 << vf.vouts1 << "\n";
				out2 << vf.vouts2 << "\n";
				out1.flush();
				out2.flush();
				

				
				//cout << vf.vmax/(double)(nSeqLength2) << "\n";
				//return;
				dscores[ncountseq] = vf.vmax/(double)(nSeqLength2);
				
				if(dscores[ncountseq] < dscoremin) dscoremin = dscores[ncountseq];
				
				//cout << vf.vouts1 << "\n";
				//cout << vf.vouts2 << "\n";
				//cout << vf.vmax << "\n\n";
				//return;	
				
				pallseqtotal += exp(vf.vmax);

				for(k = 0; k < nstates; k++)
				{					
					for(l = 0; l < nstates; l++)
					{										
						asums[k][l] += asum[k][l];																		
					}			
				}	
				
				ncountseq++;

			}
			
			if(ncount == 0) pallseqinitial = pallseqtotal;
				
			delete[] seq1;
			delete[] seq2;
			
			i++;
			
		}	//while(!fin.IsEnd() && i < ntrainseq)
		
		out1.close();
		out2.close();
		
		cout << "pallseqtotal = " << log(pallseqtotal) << "\n";
		//return;	
		
		if(ncount > 0)
		{
			pdiff = log(pallseqtotal) - log(pallseqtotalold);
		
			cout << "\n" << ncount << "," << pallseqtotal << "," << pallseqtotalold << "\n";
			cout << log(pallseqtotalold) << "," << log(pallseqtotal) << "," << pdiff << "\n\n";
			
			if(pdiff <= dthresh && ncount > mincount) btrain = 0;		
			else if(ncount >= maxcount) btrain = 0;						
			
			//cout << ncount << "," << maxcount << "," << btrain << "\n";
		}
				
				
		pallseqtotalold = pallseqtotal;			
		asmin = log(0.0000001); //log(0.01); //log(0.0000001); //log(0.00000000000001);
		
		//cout << "hmmp: " << hmmprobs(0,0,2) << "," << hmmprobs(0,1,2) << "\n";
		
		if(btrain == 1)
		{		
			
			//CalcProbMatrixFromTrainDataNew(soutputfile1, soutputfile2, sinputprobsfile);	
			//ReadinAllOutputProbs(hmmprobs, hmmprobsho, norder, noutputs, noutputsho, sinputprobsfile);
			
			CreateProbMatrixFromTrainDataVariedOrder(soutputfile1, soutputfile2, hmmprobs, hmmprobsho1, hmmprobsho2, hmmprobsho3, noutputs, noutputsho1, noutputsho2, noutputsho3);
			
			cout << "\nhmmprobs - match state:\n";
			for(k = 0; k < noutputs; k++)
			{
				for(l = 0; l < noutputs; l++)
				{
					cout << exp(hmmprobs(k,l,2)) << ",";
				}
			
				cout << "\n";
			}
		
			cout << "hmmprobs - q1:\n";
			for(k = 0; k < noutputs; k++)
			{
				for(l = 0; l < noutputs; l++)
				{
					cout << exp(hmmprobs(k,l,1)) << ",";
				}
				
				cout << "\n";
			}
		
			cout << "hmmprobs - q2:\n";
			for(k = 0; k < noutputs; k++)
			{
				for(l = 0; l < noutputs; l++)
				{
					cout << exp(hmmprobs(k,l,3)) << ",";
				}
			
				cout << "\n";
			}
			cout << "\n";
			
						
			//Calculate the total number of outputs and total number of transfers across the training sequences (esumsums and asumsums)		
			for(k = 0; k < nstates; k++)
			{												
				asumsums[k] = 0;					
				for(l = 0; l < nstates; l++) asumsums[k] += asums[k][l];	
			}
		
			boverfit = 0;
				
			for(k = 0; k < nstates; k++)
			{										
						
				for(l = 0; l < nstates; l++)
				{				
					double aold  = exp(as[k][l]);
									
					as[k][l] = AMIN;
					if(asumsums[k] > 0 && asums[k][l] > 0) 
					{						
						atemp = log(asums[k][l]/asumsums[k]);
						
						if(atemp < asmin) 
						{
							atemp = asmin;
							boverfit = 1;
						}
						as[k][l] = atemp;			
					}
					else
					{
						//cout << k << "," << l << ": " << asumsums[k] << "," << asums[k][l] <<"\n";
					}
					
					if(statetypes[k] != ENDSTATE && statetypes[l] == ENDSTATE) as[k][l] = log(0.001); 
				}					
																		
			}
			

			if(boverfit == 1) cout << "\n----------Warning: overfitting may have occured----------\n";
			
		}
		//cout << "hmmp: " << hmmprobs(0,0,2) << "," << hmmprobs(0,1,2) << "\n";		
				
		ncount++;	
		
		
		infiletrain1.close();	
		infiletrain2.close();

		//btrain = 0;	
		
	} // end train	

	infiletrain1.close();
	infiletrain2.close();	
	//infiletrain.close();
	
	
	double dscoremean = 0;
	for(i = 0; i < ncountseq; i++)
	{
		dscoremean += dscores[i];
	}
	
	dscoremean = dscoremean/(double)ncountseq;
	
	double dstd = 0;
	for(i = 0; i < ncountseq; i++)
	{
		dstd += (dscores[i]-dscoremean)*(dscores[i]-dscoremean);
	}
	
	dstd = sqrt(dstd/(double)ncountseq);
	
	cout << "dscoremean = " << dscoremean << ", dstd = " << dstd << ", dscoremin = " << dscoremin << "\n";
	
	
		cout << "as:\n";
		for(k = 0; k < nstates; k++)
		{
			for(l = 0; l < nstates; l++)
			{
				cout << exp(as[k][l]) << ",";
			}
			cout << "\n";
		}
		cout << "\n";
		
		cout << "hmmprobs - match state:\n";
		for(k = 0; k < noutputs; k++)
		{
			for(l = 0; l < noutputs; l++)
			{
				cout << exp(hmmprobs(k,l,2)) << ",";
			}
			
			cout << "\n";
		}
		
		cout << "hmmprobs - q1:\n";
		for(k = 0; k < noutputs; k++)
		{
			for(l = 0; l < noutputs; l++)
			{
				cout << exp(hmmprobs(k,l,1)) << ",";
			}
			
			cout << "\n";
		}
		
		cout << "hmmprobs - q2:\n";
		for(k = 0; k < noutputs; k++)
		{
			for(l = 0; l < noutputs; l++)
			{
				cout << exp(hmmprobs(k,l,3)) << ",";
			}
			
			cout << "\n";
		}
		
		cout << "\n";
	
	
	

	ofstream savehmmfile(soutresults.c_str());
	
	//output number of states, outputs, etc.
	savehmmfile << "n1;" << nstates << ";\n"; 
	savehmmfile << "n2;" << noutputs << ";\n"; 
	savehmmfile << "n3;" << nbeginstate << ";\n";
	savehmmfile << "n4;" << nendstate << ";\n";
	savehmmfile << "n5;" << nbases << ";\n";
	savehmmfile << "n6;" << norder << ";\n";
		
	
	
	//output statetypes
	savehmmfile << "st;";
	for(i = 0; i < nstates; i++)
	{
		savehmmfile << statetypes[i] << ",";
	}	
	savehmmfile << "\n";
	
	
	//output transitions
	savehmmfile << "as;";	
	for(i = 0; i < nstates; i++)
	{
		for(j = 0; j < nstates; j++)
		{
			savehmmfile << as[i][j] << ",";
		}									
	}	
	savehmmfile << "\n";
	
		
	savehmmfile << "q0x;";	
	for(k = 0; k < noutputs; k++)
	{		
		savehmmfile << hmmprobs(k,0,1) << ",";		
	}
	savehmmfile << "\n";
		
	savehmmfile << "q0y;";	
	for(k = 0; k < noutputs; k++)
	{		
		savehmmfile << hmmprobs(k,0,3) << ",";		
	}
	savehmmfile << "\n";
	
	
	savehmmfile << "p0;";	
	for(k = 0; k < noutputs; k++)
	{
		for(l = 0; l < noutputs; l++)
		{
			savehmmfile << hmmprobs(k,l,2) << ",";
		}
	}
	savehmmfile << "\n";
	
	OutputToFileOutputProbs(savehmmfile, hmmprobsho1, 1, 1, noutputsho1);
	OutputToFileOutputProbs(savehmmfile, hmmprobsho1, 1, 3, noutputsho1);
	OutputToFileAlignProbs(savehmmfile, hmmprobsho1, 1, 2, noutputsho1);
	
	OutputToFileOutputProbs(savehmmfile, hmmprobsho2, 2, 1, noutputsho2);
	OutputToFileOutputProbs(savehmmfile, hmmprobsho2, 2, 3, noutputsho2);
	OutputToFileAlignProbs(savehmmfile, hmmprobsho2, 2, 2, noutputsho2);
	
	OutputToFileOutputProbs(savehmmfile, hmmprobsho3, 3, 1, noutputsho3);
	OutputToFileOutputProbs(savehmmfile, hmmprobsho3, 3, 3, noutputsho3);
	OutputToFileAlignProbs(savehmmfile, hmmprobsho3, 3, 2, noutputsho3);


	
	savehmmfile.close();			
		
	delete[] qxs;
	delete[] qys;
	//delete[] ps1[];
	
	delete[] hmmprobs;		
	delete[] hmmprobsho1;	
	delete[] hmmprobsho2;	
	delete[] hmmprobsho3;
	delete[] hmmesum;		
	delete[] hmmesums;		
	delete[] hmmasums;	
	delete[] esumsums;
	delete[] asumsums;		
	delete[] statetypes;	
	
	delete[] dscores;

	Clear2DimArray(matchprobs, noutputs);
	Clear2DimArray(as, nstates);
	
	Clear2DimArray(asum, nstates);
	Clear2DimArray(asums, nstates);
	
	Clear2DimArray(alabels, nstates);				
	Clear2DimArray(ps2, nbases);	
	
}

*/


//Train(sfiletrain, sfileout, norder, bVariableOrder, nminlen, nmintrainlength, nmaxiterations, dTrainStopThresh);


void Train(string sinputfiletrain, string soutresults, int norder, int bVariableOrder, 
		  int nminlen, int maxcount, double dthresh)
{	
	
	int i; 
	int j; 
	int k; 
	int l;
	int nbeginstate; 
	int nendstate; 
	int noutputs; 
	int noutputsho; 
	int noutputsho1;
	int noutputsho2;
	int noutputsho3;
	int nstates;
	int statetypefrom, statetypeto;
	int nbases;
	int nSeqLength1;
	int nSeqLength2;	
	int mincount;
	int ntrainseq;
	int nmaxlength;	
	int btrain;
	int binit;
	int boutput;
	int ntotaltrainseq;
	int ncount;	 
	int boverfit;
	int ncountseq;

	double d;	
	double pallseqtotal;
	double pallseqtotalold;
	double pallseqinitial;
	double pdiff;
	double dtau;			
	double asmin;
	double atemp;
	double dscoremin;
	double dstd;
	double dscoremean;
	
	string line;
	string sinputfiletrain1;
	string sinputfiletrain2;
	string sinputfiletrain1b;
	string sinputfiletrain2b;
	string sinputfiletest;
	string soutputfile1;
	string soutputfile2;	
	string sseq1;
	string sseq2;
	
	double* hmmprobsho;
	double* hmmprobsho1;
	double* hmmprobsho2;
	double* hmmprobsho3;
	
	vfwdbackpair vf;

	nstates = 5;
	nbeginstate = 0;
	nendstate = nstates-1;
	
	cout << "\nInitializing variables and arrays...\n";
	
	if(bVariableOrder == 0)
	{
		noutputsho = (int)pow(4, norder+1);
		noutputsho1 = -1;
		noutputsho2 = -1;
		noutputsho3 = -1;
		
		hmmprobsho  = new double[noutputsho*noutputsho*nstates]; 
		hmmprobsho1 = NULL;
		hmmprobsho2 = NULL;
		hmmprobsho3 = NULL;
	}
	else
	{
		noutputsho = -1;
		noutputsho1 = (int)pow(4, 2);
		noutputsho2 = (int)pow(4, 3);
		noutputsho3 = (int)pow(4, 4);	
		
		hmmprobsho  = NULL;
		hmmprobsho1 = new double[noutputsho1*noutputsho1*nstates];
		hmmprobsho2 = new double[noutputsho2*noutputsho2*nstates];
		hmmprobsho3 = new double[noutputsho3*noutputsho3*nstates];
	}
		
	int* statetypes = new int[nstates];

	if(norder == 0)
	{
		statetypes[0] = BEGINSTATE;	
		statetypes[1] = OUTPUTSTATEX;
		statetypes[2] = MATCHSTATE;
		statetypes[3] = OUTPUTSTATEY;		
		statetypes[4] = ENDSTATE;	
	}
	else
	{
		statetypes[0] = BEGINSTATE;	
		statetypes[1] = OUTPUTSTATEXHO;
		statetypes[2] = MATCHSTATEHO;
		statetypes[3] = OUTPUTSTATEYHO;		
		statetypes[4] = ENDSTATE;		
	}
	
	nbases = 4;
	noutputs = 4;		
			
	pallseqtotal = 0;
	pallseqtotalold = 0;
	pdiff = 0;
	
	dthresh = fabs(log(dthresh)); 
	mincount = 0; //min number of training steps

	ntrainseq = -1; 
	nmaxlength = 500; 
	
	btrain = 1;
	boutput = 0;
	asmin = log(PMIN);
	dtau = log(0.01);	
	
	sinputfiletrain1  = sinputfiletrain;
	sinputfiletrain2  = sinputfiletrain;
	sinputfiletrain1b = sinputfiletrain;
	sinputfiletrain2b = sinputfiletrain;
	
	sinputfiletrain1  += "1";
	sinputfiletrain2  += "2";	
	sinputfiletrain1b += "1b";
	sinputfiletrain2b += "2b";
	
	soutputfile1 = soutresults;
	soutputfile2 = soutresults;
	soutputfile1 += "_out1";
	soutputfile2 += "_out2";
	
	int* seq1;
	int* seq2;	

	double* dscores;
	
	double* hmmprobs = new double[noutputs*noutputs*nstates];	
	double* hmmesum = new double[noutputs*noutputs*nstates];
	
	double** asum = CreateDoubleArray(nstates, nstates);
	double** asums = CreateDoubleArray(nstates, nstates);

	double* hmmesums = new double[noutputs*noutputs*nstates];
	double* hmmasums = new double[noutputs*noutputs*nstates];
	
	double* esumsums = new double[nstates];
	double* asumsums = new double[nstates];
	
	double** matchprobs = CreateDoubleArray(noutputs,noutputs);	
	double** as = CreateDoubleArray(nstates, nstates);
	string** alabels = CreateStringArray(nstates, nstates);
	
	binit = InitializeHMMVariables(nstates, norder, noutputs, noutputsho, noutputsho1, noutputsho2,
							noutputsho3, bVariableOrder, dtau, statetypes, hmmesums, hmmprobs, 
							hmmprobsho, hmmprobsho1, hmmprobsho2, hmmprobsho3, as, 
							sinputfiletrain1b, sinputfiletrain2b);
	
	if(binit == 0) 
	{
		cout << "Failed to initialize all of variables required for training\n";
		return;
	}

	/////////////////////////////////////////////////////////////
	//Start of training routines
	/////////////////////////////////////////////////////////////

	ifstream infiletrain1(sinputfiletrain1.c_str());	
	ifstream infiletrain2(sinputfiletrain2.c_str());
		
	
	if(infiletrain1.good() == 0 )
	{
		cout << "Problem accessing " << sinputfiletrain1 << ". Please check file exists.\n\n";
		return;
	}
	else if(infiletrain2.good() == 0)
	{
		cout << "Problem accessing " << sinputfiletrain2 << ". Please check file exists.\n\n";
		return;
	}
	
	ncount = 0;
	ntotaltrainseq = 0;

	while(infiletrain1.eof() == 0)
	{
		getline(infiletrain1,line);
		
		if(line.length() > 0) ntotaltrainseq++;
	}
	
	
	infiletrain1.close();	
	infiletrain2.close();	
	
	
	
		
	if(ntrainseq == -1) ntrainseq = ntotaltrainseq;
	dscores = new double[ntrainseq];
	
	ofstream out1;
	ofstream out2;
	
	cout << "Training...\n\n";

	while(btrain == 1)
	{					
		infiletrain1.open(sinputfiletrain1.c_str());	
		infiletrain2.open(sinputfiletrain2.c_str());	
		
		out1.open(soutputfile1.c_str(), std::ofstream::out | std::ofstream::trunc);
		out2.open(soutputfile2.c_str(), std::ofstream::out | std::ofstream::trunc);
		
		/////////////////Initialization
		
		pallseqtotal = 0;		
		dscoremin = 0;
		
		//Initialize E and A sums (E represents the probability output matrix and A the transfer matrix)
		//Esum contains the total number of outputs for each output type in each state
		//Esums sums Esum across all training sequences
		//Asum contains the total number of transfers from one state to another		
		//Asums sums Asum across all training sequences		
		for(k = 0; k < nstates; k++)
		{	
			esumsums[k] = 0;	//Used to sum across outputs for esums
			asumsums[k] = 0;	//Used to sum across states for asums
					
			for(j = 0; j < noutputs; j++)
			{
				for(l = 0; l < noutputs; l++) hmmesums(j,l,k) = 0;						
			}			
			
			for(l = 0; l < nstates; l++) asums[k][l] = 0;			
		}
	
		////////////////Begin training

		i = 0;		
		ncountseq = 0;

		while(infiletrain1.eof() == 0 && i < ntrainseq) 
		{					
			//Asum and Esum are summed across each training sequence so initialize at start of each training sequence
			//as opposed to esums and asums which are across all the training sequences
			for(k = 0; k < nstates; k++)
			{	
				for(j = 0; j < noutputs; j++)
				{
					for(l = 0; l < noutputs; l++) hmmesum(j,l,k) = 0;				
				}			
			
				for(l = 0; l < nstates; l++)  asum[k][l] = 0;		
							
			}
		
			getline(infiletrain1,line);
			sseq1 = line;		
			getline(infiletrain2,line);
			sseq2 = line;
			
			nSeqLength1 = sseq1.length();
			nSeqLength2 = sseq2.length();
			
			seq1 = new int[nSeqLength1];
			seq2 = new int[nSeqLength2];
	
			ConvertStringToInts(sseq1, seq1, nbases);
			ConvertStringToInts(sseq2, seq2, nbases);
		

			if(nSeqLength1 > nminlen && nSeqLength2 > nminlen && nSeqLength1 < nmaxlength && nSeqLength2 < nmaxlength) 			
			{		
				if(bVariableOrder == 0)
				{									 
					vf = ViterbiTrain(seq1, seq2, sseq1, sseq2, nSeqLength1, nSeqLength2,
	 					    noutputs, nstates, hmmprobs, hmmprobsho, as, statetypes, hmmesum, asum,						   
						    nbeginstate, nendstate, norder, noutputsho);
				}
				else
				{
					vf = ViterbiTrainVariedOrderHMM(seq1, seq2, sseq1, sseq2, nSeqLength1, nSeqLength2,
	 					    noutputs, nstates, hmmprobs, hmmprobsho1, hmmprobsho2, hmmprobsho3, as, statetypes, hmmesum, asum,						   
						    nbeginstate, nendstate, norder, noutputsho1, noutputsho2, noutputsho3);

				}
						    
				out1 << vf.vouts1 << "\n";
				out2 << vf.vouts2 << "\n";
				out1.flush();
				out2.flush();
				
				dscores[ncountseq] = vf.vmax/(double)(nSeqLength2);
				
				if(dscores[ncountseq] < dscoremin) dscoremin = dscores[ncountseq];
				
				pallseqtotal += exp(vf.vmax);

				for(k = 0; k < nstates; k++)
				{					
					for(l = 0; l < nstates; l++)
					{										
						asums[k][l] += asum[k][l];																		
					}			
				}	
				
				ncountseq++;

			}
				
			delete[] seq1;
			delete[] seq2;
			
			i++;
			
		}	//while(!fin.IsEnd() && i < ntrainseq)
		
		out1.close();
		out2.close();
		
		cout << "pallseqtotal = " << log(pallseqtotal) << "\n";	
		
		if(ncount > 0)
		{
			pdiff = log(pallseqtotal) - log(pallseqtotalold);
		
			cout << "\n" << ncount << "," << pallseqtotal << "," << pallseqtotalold << "\n";
			cout << log(pallseqtotalold) << "," << log(pallseqtotal) << "," << pdiff << "\n\n";
			
			if(dthresh != -1 && pdiff <= dthresh && ncount > mincount) btrain = 0;		
			else if(maxcount != -1 && ncount >= maxcount) btrain = 0;												
		}
				
				
		pallseqtotalold = pallseqtotal;			
		 
		
		if(btrain == 1)
		{		
			
			if(bVariableOrder == 0)
			{
				CreateProbMatrixFromTrainData(soutputfile1, soutputfile2, hmmprobs, hmmprobsho, 
											  norder, noutputs, noutputsho);
			}
			else
			{
				CreateProbMatrixFromTrainDataVariedOrder(soutputfile1, soutputfile2, hmmprobs,
														 hmmprobsho1, hmmprobsho2, hmmprobsho3, noutputs, 
														 noutputsho1, noutputsho2, noutputsho3);
	
			}
			
			cout << "\nhmmprobs - match state:\n";
			for(k = 0; k < noutputs; k++)
			{
				for(l = 0; l < noutputs; l++)
				{
					cout << exp(hmmprobs(k,l,2)) << ",";
				}
			
				cout << "\n";
			}
		
			cout << "hmmprobs - q1:\n";
			for(k = 0; k < noutputs; k++)
			{
				for(l = 0; l < noutputs; l++)
				{
					cout << exp(hmmprobs(k,l,1)) << ",";
				}
				
				cout << "\n";
			}
		
			cout << "hmmprobs - q2:\n";
			for(k = 0; k < noutputs; k++)
			{
				for(l = 0; l < noutputs; l++)
				{
					cout << exp(hmmprobs(k,l,3)) << ",";
				}
			
				cout << "\n";
			}
			cout << "\n";
			
						
			//Calculate the total number of outputs and total number of transfers across the training sequences (esumsums and asumsums)		
			for(k = 0; k < nstates; k++)
			{												
				asumsums[k] = 0;					
				for(l = 0; l < nstates; l++) asumsums[k] += asums[k][l];	
			}
		
			boverfit = 0;
				
			for(k = 0; k < nstates; k++)
			{										
						
				for(l = 0; l < nstates; l++)
				{				
					double aold  = exp(as[k][l]);
									
					as[k][l] = AMIN;
					if(asumsums[k] > 0 && asums[k][l] > 0) 
					{						
						atemp = log(asums[k][l]/asumsums[k]);
						
						if(atemp < asmin) 
						{
							atemp = asmin;
							boverfit = 1;
						}
						as[k][l] = atemp;			
					}
					else
					{
						//cout << k << "," << l << ": " << asumsums[k] << "," << asums[k][l] <<"\n";
					}
					
					if(statetypes[k] != ENDSTATE && statetypes[l] == ENDSTATE) as[k][l] = log(0.001); 
				}					
																		
			}
			
			
			
			cout << "as:\n";
			for(k = 0; k < nstates; k++)
			{
				for(l = 0; l < nstates; l++)
				{
					cout << exp(as[k][l]) << ",";
				}
				cout << "\n";
			}
			cout << "\n";
			

			if(boverfit == 1) cout << "\n----------Warning: overfitting may have occured----------\n";
			
		}			
				
		ncount++;	
				
		infiletrain1.close();	
		infiletrain2.close();	
		
	} // end train	

	infiletrain1.close();
	infiletrain2.close();	
	
	
	
	dscoremean = 0;
	for(i = 0; i < ncountseq; i++)
	{
		dscoremean += dscores[i];
	}
	
	dscoremean = dscoremean/(double)ncountseq;
	
	dstd = 0;
	for(i = 0; i < ncountseq; i++)
	{
		dstd += (dscores[i]-dscoremean)*(dscores[i]-dscoremean);
	}
	
	dstd = sqrt(dstd/(double)ncountseq);
	
	cout << "dscoremean = " << dscoremean << ", dstd = " << dstd << ", dscoremin = " << dscoremin << "\n";
	
	
		cout << "as:\n";
		for(k = 0; k < nstates; k++)
		{
			for(l = 0; l < nstates; l++)
			{
				cout << exp(as[k][l]) << ",";
			}
			cout << "\n";
		}
		cout << "\n";
		
		cout << "hmmprobs - match state:\n";
		for(k = 0; k < noutputs; k++)
		{
			for(l = 0; l < noutputs; l++)
			{
				cout << exp(hmmprobs(k,l,2)) << ",";
			}
			
			cout << "\n";
		}
		
		cout << "hmmprobs - q1:\n";
		for(k = 0; k < noutputs; k++)
		{
			for(l = 0; l < noutputs; l++)
			{
				cout << exp(hmmprobs(k,l,1)) << ",";
			}
			
			cout << "\n";
		}
		
		cout << "hmmprobs - q2:\n";
		for(k = 0; k < noutputs; k++)
		{
			for(l = 0; l < noutputs; l++)
			{
				cout << exp(hmmprobs(k,l,3)) << ",";
			}
			
			cout << "\n";
		}
		
		cout << "\n";
	
	
	

	ofstream savehmmfile(soutresults.c_str());
	
	//output number of states, outputs, etc.
	savehmmfile << "n1;" << nstates << ";\n"; 
	savehmmfile << "n2;" << noutputs << ";\n"; 
	savehmmfile << "n3;" << nbeginstate << ";\n";
	savehmmfile << "n4;" << nendstate << ";\n";
	savehmmfile << "n5;" << nbases << ";\n";
	savehmmfile << "n6;" << norder << ";\n";
	savehmmfile << "n7;" << bVariableOrder << ";\n";

		
	//output statetypes
	savehmmfile << "st;";
	for(i = 0; i < nstates; i++)
	{
		savehmmfile << statetypes[i] << ",";
	}	
	savehmmfile << "\n";
	
	
	//output transitions
	savehmmfile << "as;";	
	for(i = 0; i < nstates; i++)
	{
		for(j = 0; j < nstates; j++)
		{
			savehmmfile << as[i][j] << ",";
		}									
	}	
	savehmmfile << "\n";
	
		
	savehmmfile << "q0x;";	
	for(k = 0; k < noutputs; k++)
	{		
		savehmmfile << hmmprobs(k,0,1) << ",";		
	}
	savehmmfile << "\n";
		
	savehmmfile << "q0y;";	
	for(k = 0; k < noutputs; k++)
	{		
		savehmmfile << hmmprobs(k,0,3) << ",";		
	}
	savehmmfile << "\n";
	
	
	savehmmfile << "p0;";	
	for(k = 0; k < noutputs; k++)
	{
		for(l = 0; l < noutputs; l++)
		{
			savehmmfile << hmmprobs(k,l,2) << ",";
		}
	}
	savehmmfile << "\n";
	
	
	if(bVariableOrder == 0)
	{
		OutputToFileOutputProbs(savehmmfile, hmmprobsho, norder, 1, noutputsho);
		OutputToFileOutputProbs(savehmmfile, hmmprobsho, norder, 3, noutputsho);
		OutputToFileAlignProbs(savehmmfile, hmmprobsho, norder, 2, noutputsho);
		
		delete[] hmmprobsho;
	}
	else
	{
		OutputToFileOutputProbs(savehmmfile, hmmprobsho1, 1, 1, noutputsho1);
		OutputToFileOutputProbs(savehmmfile, hmmprobsho1, 1, 3, noutputsho1);
		OutputToFileAlignProbs(savehmmfile, hmmprobsho1, 1, 2, noutputsho1);
	
		OutputToFileOutputProbs(savehmmfile, hmmprobsho2, 2, 1, noutputsho2);
		OutputToFileOutputProbs(savehmmfile, hmmprobsho2, 2, 3, noutputsho2);
		OutputToFileAlignProbs(savehmmfile, hmmprobsho2, 2, 2, noutputsho2);
	
		OutputToFileOutputProbs(savehmmfile, hmmprobsho3, 3, 1, noutputsho3);
		OutputToFileOutputProbs(savehmmfile, hmmprobsho3, 3, 3, noutputsho3);
		OutputToFileAlignProbs(savehmmfile, hmmprobsho3, 3, 2, noutputsho3);
		
		delete[] hmmprobsho1;
		delete[] hmmprobsho2;
		delete[] hmmprobsho3;
	}
	
	savehmmfile.close();			

	delete[] hmmprobs;		
	delete[] hmmesum;		
	delete[] hmmesums;		
	delete[] hmmasums;	
	delete[] esumsums;
	delete[] asumsums;		
	delete[] statetypes;	
	
	delete[] dscores;
	
	Clear2DimArray(matchprobs, noutputs);
	Clear2DimArray(as, nstates);
	
	Clear2DimArray(asum, nstates);
	Clear2DimArray(asums, nstates);
	
	Clear2DimArray(alabels, nstates);	
	
	cout << "Finished!\n\n";			
	
}






void CreateTrainingDataFromMafFile(string sfilemaf, string sfileout, double dmaxmismap)
{
	int i;
	int n1;
	int n2;
	double d;
	
	string s1;
	string s2;
	string sout1;
	string sout2;
	string line1;
	string line2;
	string line3;
	string line4;
	string sfileout1;
	string sfileout2;
	string sfileout1b;
	string sfileout2b;
	
	sfileout1 = sfileout;
	sfileout1 += "1";
	
	sfileout2 = sfileout;
	sfileout2 += "2";	
	
	sfileout1b = sfileout;
	sfileout1b += "1b";
	
	sfileout2b = sfileout;
	sfileout2b += "2b";	
	
	ifstream finmaf(sfilemaf.c_str());
	
	if(finmaf.good() == 0) 
	{
		cout << "\nFailed to find .maf input file, please check file name and that file input and output arguments were passed correctly.\n\n";
		return;
	}
	
	
	ofstream fout(sfileout.c_str());
	ofstream fout1(sfileout1.c_str());
	ofstream fout2(sfileout2.c_str());
	ofstream fout1b(sfileout1b.c_str());
	ofstream fout2b(sfileout2b.c_str());
	

	cout << "\nCreating training data files...\n\n";
	
	while(finmaf.eof() == 0)
	{
		getline(finmaf, line1);								
		
		if(line1.length() > 1 && line1[0] == 'a')
		{			
			n1 = line1.find("mismap=",0)+7;
			
			s1 = line1.substr(n1, line1.length()-n1);
			d = atof(s1.c_str());
			
			getline(finmaf, line2);
			getline(finmaf, line3);
			getline(finmaf, line4);
			
			if(d < dmaxmismap)
			{
				n1 = line2.find("+", 0);
				if(n1 < 0) n1 = line2.find("-", 0);			
				n1++;
				while(line2[n1] == ' ') n1++;
				n1 = line2.find(" ", n1)+1;
			
				n2 = line3.find("+", 0);
				if(n2 < 0) n2 = line3.find("-", 0);
				n2++;
				while(line3[n2] == ' ') n2++;			
				n2 = line3.find(" ", n2)+1;
		
				s1 = line2.substr(n1, line2.length()-n1);
				s2 = line3.substr(n2, line3.length()-n2);
			
				sout1 = "";
				sout2 = "";
			
				for(i = 0; i < s1.length(); i++) 			
					if(s1[i] != '-') sout1 += s1[i];
			
				for(i = 0; i < s2.length(); i++) 			
					if(s2[i] != '-') sout2 += s2[i];

				fout << sout1 << "\n";
				fout << sout2 << "\n\n";					
			
				fout1 << sout1 << "\n";
				fout2 << sout2 << "\n";	
				
				fout1b << s1 << "\n";
				fout2b << s2 << "\n";			
			}
		}
	}
	
	cout << "Finished!\n\n";
}


void PrintHelp(int bPrintVersion)
{
	cout << "\n"; 
	cout << "How to run (training):\n";
	cout << "----------------------\n";
	cout << "varhmmtrain [training data file] [output file] [norder] (optional arguments)\n";
	cout << "\n";	
	cout << "Important: If an order of 3 is specified, a variable ordered HMM is generated and\n";
	cout << "trained. For 0th and higher order models, 0th to 2nd order models can be generated.\n";
	cout << "Greater values than 3 for [norder] are not supported.\n\n";

	cout << "Optional arguments:\n";
	cout << "-------------------\n";
	cout << "-" << sMINTRAINLENGTHCMD     << sMINTRAINLENGTHDESC      << "[" << MINTRAINLENGTH     << "]" << "\n";
	cout << "-" << sMAXTRAINITERATIONSCMD << sMAXTRAINITERATIONSDESC  << "[" << MAXTRAINITERATIONS << "]" << "\n";
	cout << "-" << sTRAINSTOPTHRESHCMD    << sTRAINSTOPTHRESHDESC     << "[" << TRAINSTOPTHRESH    << "]" << "\n";
	cout <<                                  sTRAINSTOPTHRESHDESC2    <<"\n";
	cout << "\n";
	
	cout << "How to create training data from maf file:\n";
	cout << "------------------------------------------\n";
	cout << "varhmmtrain -C [maf file] [output file] (alignments with mismap prob. < 0.001 used)\n";
	cout << "varhmmtrain -C [maf file] [output file] -m [max mismap probability]\n\n";
	cout << "Note: This generates 5 seperate files that use [output file] as the base filename. Use\n";
	cout << "the base filename as given by [output file] for varhmmtrain [training data file].\n\n";
	
	cout << "Examples:\n";
	cout << "---------\n";
	cout << "varhmmtrain mytrainset hmm2_mytrainset 2\n";
	cout << "varhmmtrain mytrainset hmm3_mytrainset 3 -i 10 -d 0.5 -l 150\n";
	cout << "varhmmtrain -C alignments.maf mytrainset\n";
	
	cout << "\n";

	if(bPrintVersion == 1)
	{
		cout << "Version: " << VERSION << "\n"; 
		cout << "\n";
	}
}

int main(int argc,char *argv[])
{
	int i;
	int narg;
	int norder;
	int bCreateTrainDataset;
	int bVariableOrder;
	int nmintrainlength;				   
	int nmaxiterations;
						
	double dmaxmismap;
	double dtrainstopthresh;
	
	string s1;
	string s2;
	string sarg;
	string sfiletrain;
	string sfilemaf;
	string sfileout;
		
	bCreateTrainDataset = 0;	
	dmaxmismap = MAXMISMAP;
	nmintrainlength = MINTRAINLENGTH;	
	nmaxiterations = MAXTRAINITERATIONS;
	dtrainstopthresh = TRAINSTOPTHRESH;
	
	narg = argc;
	
	
	if(narg < 4)
	{
		PrintHelp(1);
		return 0;
	}
	
	
	i = 1;	
	while(i < narg)
	{
		s1 = argv[i];
		
		//Check that an optional argument isn't passed before required ones
		
		if(i == 1)
		{
			if(s1.length() > 0 && s1[0] == '-' && s1[1] != 'C')
			{
				PrintHelp(0);
				cout << "Error: an optional argument appears to have been placed before the first 3 required file arguments. Please see above syntax.\n\n";
				return 0;
			}			
		}
		else if(i <= 3)
		{
			if(s1[0] == '-' && bCreateTrainDataset == 0)
			{
				PrintHelp(0);
				cout << "Error: an optional argument appears to have been placed before the first 3 required file arguments. Please see above syntax.\n\n";
				return 0;
			}	
		}		
		

		switch(i)
		{
			case 1: 
				if(s1.length() > 0 && s1[0] == '-' && s1[1] == 'C') bCreateTrainDataset = 1;
				else sfiletrain = s1;
			break;
			case 2: 
				if(bCreateTrainDataset == 1) sfilemaf = s1;
				else sfileout = s1;					
			break;
			case 3:
				if(bCreateTrainDataset == 1) sfileout = s1;
				else 
				{
					if(s1[0] == '0' || s1[0] == '1' || s1[0] == '2' || s1[0] == '3')
					{
						norder = atoi(s1.c_str());
					}
					else
					{
						PrintHelp(0);
						cout << "Error: only order values between 0 and 3 supported. Please see help above.\n";
					}
				}		
			break;
			default:	
				if(bCreateTrainDataset == 1)
				{					
					if(i == 4)
					{
						if(s1.length() > 1 && s1[0] == '-' && s1[1] == 'm')
						{							
						}
						else
						{
							PrintHelp(0);
							cout << "Error: error in argument(s) passed after -C [hmm training file]. See help for details\n";
						}
					}
					else if(i == 5)
					{
						dmaxmismap = atof(s1.c_str());
					}
				}
				else
				{
					if(s1[0] != '-')
					{
						PrintHelp(0);
						cout << "Error: there appears to be a '-' missing in the optional argument " << s1 << ". Please see above syntax.\n\n";
						return 0;
					}
					sarg = s1.substr(1, s1.length()-1);
				
					i++;
					if(i > narg-1)
					{
						PrintHelp(0);
						cout << "Error: there appears to be an argument missing after " << s1 << ". Please see above syntax.\n\n";
						return 0;
					}
					s2 = argv[i];		
										
					if(sarg.compare(sMINTRAINLENGTHCMD) == 0)
					{
						nmintrainlength = atoi(s2.c_str());				   
				    
						if(nmintrainlength < 0) 
					    {
							cout << "\nError: Invalid value for " << "-" << sMINTRAINLENGTHCMD << "n";
							return 0;
						}
					}
					else if(sarg.compare(sMAXTRAINITERATIONSCMD) == 0)
					{
						nmaxiterations = atoi(s2.c_str());
						
						if(nmaxiterations == 0 || nmaxiterations < -1) 
					    {
							cout << "\nError: Invalid value for " << "-" << MAXTRAINITERATIONS << "n";
							return 0;
						}
					}
					else if(sarg.compare(sTRAINSTOPTHRESHCMD) == 0)
					{
						dtrainstopthresh =  atof(s2.c_str());
						
						if(dtrainstopthresh < -1)
						{
							cout << "\nError: Invalid value for " << "-" << sTRAINSTOPTHRESHCMD << "n";
							return 0;
						}
					}
				}	
					
			break;	
		}
		
		i++;
	}
	
	
	if(bCreateTrainDataset == 0)
	{
		if(norder <= 2) bVariableOrder = 0;
		else if(norder == 3) bVariableOrder = 1;
		else	
		{
	
			PrintHelp(0);
			cout << "Error: Only order values between 0 and 3 supported. Please see help above.\n";
		}
	
		if(dtrainstopthresh == -1 && nmaxiterations == -1)
		{
			
			PrintHelp(0);
			cout << "Error: Both " <<  "-" << MAXTRAINITERATIONS << " and " <<  "-" << sTRAINSTOPTHRESHCMD << " cannot be set to 1. Please see help above.\n";
		}
	}
	
	if(bCreateTrainDataset==1)
	{
		CreateTrainingDataFromMafFile(sfilemaf, sfileout, dmaxmismap);
	}
	else
	{
		Train(sfiletrain, sfileout, norder, bVariableOrder, nmintrainlength, nmaxiterations, dtrainstopthresh);	
	}
	
	
}

