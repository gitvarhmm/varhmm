//Copyright 2016 Thomas M. Poulsen

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>   

#include <fstream>
#include <sstream>
#include <limits>
#include <vector>

#include <math.h> 
#include <string.h> 

#include <thread>
using namespace std;


struct charcountstruct
{
	int nA;
	int nC;
	int nG;
	int nT;
};

			
struct vitstruct
{
	double ptotalfwd;
	double ptotalback;
	double vmax;
	double vmaxdiv;
	string vouts1;
	string vouts2;
	string voutinfo;
	int nreflength;
	int ninsert;
};


struct vfwdback
{
	double ptotalfwd;
	double ptotalback;
};


struct numpos
{
	int number;
	int pos;
};


typedef int BOOL;
const int AMIN = -(int)(numeric_limits<double>::infinity()); 




const double VERSION = 1.01;

const int NOSTATETYPES = 6;

const int BEGINSTATE = 10;
const int ENDSTATE = 11;

const int MATCHSTATE = 1000;
const int OUTPUTSTATEX = 2000;
const int OUTPUTSTATEY = 3000;	

const int MATCHSTATEHO = 1002;
const int OUTPUTSTATEXHO = 2002;
const int OUTPUTSTATEYHO = 3002;
	
const int SILENTSTATE = 4000;	

const int RX1STATE = 5000;	
const int RY1STATE = 6000;	

const int RX2STATE = 7000;	
const int RY2STATE = 8000;


const int MAXCOUNTDIFF = 7;
const int COUNTWINDOWLENGTH = 40;
const int MAXREALIGN = 5;
const double MINSCORE = -2.5;
const double MINSCORE2 = -3.2;
const double MINALIGNSCORE = -2.5; 
const int THREADS = 0;
const int BFASTQ = 1;
const string PRECALCMATRIXFILE = "";  
const string OUTPUTFILE = ""; 

const string sMAXCOUNTDIFFARGCMD   = "d";
const string sCOUNTWINDOWLENGTHCMD = "l";
const string sMAXREALIGNCMD        = "n";
const string sMINSCORECMD          = "s8";
const string sMINSCORE2CMD         = "s24";
const string sMINALIGNSCORECMD     = "s";
const string sBFASTQCMD       	   = "Q";
const string sTHREADSCMD           = "T";
const string sPRECALCMATRIXFILECMD = "M";  
const string sOUTPUTFILECMD        = "O";

const string sMAXCOUNTDIFFARGDESC   = "   Max. count difference between 2 sequences   ";
const string sCOUNTWINDOWLENGTHDESC = "   Count window length used with count diff.   ";
const string sMAXREALIGNDESC        = "   Max. number of alignment searches           ";
const string sMINSCOREDESC          = "  Min. score required for 8-length alignment  ";
const string sMINSCORE2DESC         = " Min. score required for 24-length alignment ";
const string sMINALIGNSCOREDESC     = "   Min. score required for full alignment      ";
const string sBFASTQDESC            = "   Input format for reads (0=fasta, 1=fastq)   ";
const string sTHREADSDESC           = "   Number of  threads                          ";
const string sPRECALCMATRIXFILEDESC = "   Filename: 8-length score matrix file (if not given,\n     only exact alignments used in seeding step)";  
const string sOUTPUTFILEDESC        = "   Filename: Output file (if not given, then generated\n     based on reference,reads,and training file names)"; 

								   


#define hmmx(i,j,k) hmmx[(k+j*nstates+i*nindexfactor)]
#define hmmy(i,j,k) hmmy[(k+j*nstates+i*nindexfactor)]

#define hmmv(i,j,k) hmmv[(k+j*nstates+i*nindexfactor)]
#define hmmptrs(i,j,k) hmmptrs[(k+j*nstates+i*nindexfactor)]

#define hmma(i,j,k) hmma[(k+j*nstates+i*nindexfactor)]
#define hmmf(i,j,k) hmmf[(k+j*nstates+i*nindexfactor)]
#define hmmb(i,j,k) hmmb[(k+j*nstates+i*nindexfactor)]

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


unsigned char** CreateUnsignedCharArray(int rows, int cols)
{
	int i,j;
	
	unsigned char** a = new unsigned char*[rows];	
	
	for(i = 0; i < rows; i++)
	{
		a[i] = new unsigned char[cols];
		
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
		
		for(j = 0; j < cols; j++)
		{
			a[i][j] = 0;
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

void Clear2DimArray(string** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}

void Clear2DimArray(unsigned char** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}

void Clear2DimArray(char** v, int len)
{
	int i; 
	
	for(i = 0; i < len; i++) delete[] v[i];
	
	delete[] v;
}
	
int ConvertFourthOrderStringToIndex(int* seq, int n)
{
	return 64*seq[n]+16*seq[n-1]+4*seq[n-2]+seq[n-3];
}

int ConvertFourthOrderStringToIndex(int n1, int n2, int n3, int n4)
{
	return 64*n4+16*n3+4*n2+n1;
}

int ConvertOrderStringToIndex(int* seq, int n, int norder)
{
	switch(norder)
	{
		case 0: return seq[n];			
		case 1: return 4*seq[n]+seq[n-1];
		case 2: return 16*seq[n]+4*seq[n-1]+seq[n-2];
		case 3: return 64*seq[n]+16*seq[n-1]+4*seq[n-2]+seq[n-3];				
	}
}

int ConvertSecondOrderStringToIndex2(int n1, int n2)
{
	return 4*n2+n1;
}

int ConvertThirdOrderStringToIndex2(int n1, int n2, int n3)
{
	return 16*n3+4*n2+n1;
}

int ConvertFourthOrderStringToIndex2(int n1, int n2, int n3, int n4)
{
	return 64*n4+16*n3+4*n2+n1;
}
	
char ReverseStrandNuc(char c)
{
	switch(c)
	{
		case 'A':
		return 'T';		
		case 'C':
		return 'G';
		case 'G':
		return 'C';
		case 'T':
		return 'A';			
	}
}

short int ConvertNucToInt(char c)
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

void ConvertStringToInts(string s, int* data, int nbase)
{
	int i;
	
	int nlen = s.length();
	
	for(i = 0; i < nlen; i++) data[i] = ConvertCharToIndex(s[i], nbase);
	
	return;
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
							nseq1 = ConvertSecondOrderStringToIndex2(i1, i2);																		
							nseq2 = ConvertSecondOrderStringToIndex2(i3, i4);									
									
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
									nseq1 = ConvertThirdOrderStringToIndex2(i1, i2, i3);																		
									nseq2 = ConvertThirdOrderStringToIndex2(i4, i5, i6);									
									
									n2 = line.find(",", n1);	
					
									s = line.substr(n1,n2-n1);
									n1 = n2+1;
					
									d = atof(s.c_str());
									hmmprobsho(nseq1, nseq2, narrayindex) = d;																																			
								}
		break;
		case 3:	
		int nbla = 0;		
			for(i1 = 0; i1 < nbases; i1++)
				for(i2 = 0; i2 < nbases; i2++)
					for(i3 = 0; i3 < nbases; i3++)
						for(i4 = 0; i4 < nbases; i4++)				
							for(i5 = 0; i5 < nbases; i5++)					
								for(i6 = 0; i6 < nbases; i6++)
									for(i7 = 0; i7 < nbases; i7++)
										for(i8 = 0; i8 < nbases; i8++)		
										{									
											nseq1 = ConvertFourthOrderStringToIndex2(i1, i2, i3, i4);																		
											nseq2 = ConvertFourthOrderStringToIndex2(i5, i6, i7, i8);									
									
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
					n = ConvertSecondOrderStringToIndex2(i1, i2);
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
						n = ConvertThirdOrderStringToIndex2(i1, i2, i3);
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
							n = ConvertFourthOrderStringToIndex2(i1, i2, i3, i4);
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
	






vitstruct ViterbiVariedOrderHMM(int* data1, int* data2, string sdata1, string sdata2, string sseqinfo, int ndatalen1, 
									   int ndatalen2, int nstates, double* as[], double* aslogs[], int* statetypes, 
									   int nbeginstate, int nendstate, double* qxs, double* qys, double** p1s, double* hmmprobsho1,
									   double* hmmprobsho2, double* hmmprobsho3, int nbases, int boutput, 
									   double* hmmv, int* hmmptrs, int nstartfrom, 
									   int nindexfactor, int norder, int noutputsho1, int noutputsho2, 
									   int noutputsho3, int* hmmx, int* hmmy, int nnearendbufferlen)			
{
	
	int i;
	int j; 
	int k; 
	int l;
	int indexi; 
	int indexj; 
	int ndata1; 
	int ndata2;
	int nstatetype; 
	int imax; 
	int vstate;

	int nlen;
	int nstart;
	int indexm; 
	int indexn;
	int ndatatotal1; 
	int ndatatotal2;
	int nindex;
	int mindex;
	int ninsert;
	int nmatches;
	
	int nif;
	int nifm1; 
	int njf; 
	int njfm1;	
	int norderuse;
	int nnearend;
	int xcount;
	int ycount;
	int xsum; 
	int ysum; 
	int xmax; 
	int ymax;
	
	
	double dscore;
	double dmaxscore;
	double a;
	double el; 
	double dmax; 
	double vsum; 
	double vmax;		
		
	string s1;
	string s2; 
	string vouts1; 
	string vouts2; 
	string voutinfo;
	
	vitstruct vf;

	vouts1 = ""; 
	vouts2 = "";
	voutinfo = "";

	nstart = 0;
		
	//Add one instance before the start of the data (for the begin state)	
	ndatatotal1 = ndatalen1+1; 
	ndatatotal2 = ndatalen2+1;		

	nnearend = ndatalen1 - nnearendbufferlen;
	
	xcount = 0;
	ycount = 0;

	if(nstartfrom == 0)
	{
		for(i = 0; i < nstates; i++)
		{		
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
	}
		

	//Forward algorithm and viterbi, at i = 0 and j = 0 we are in the begin state
	for(i = 0; i < nstartfrom; i++)
	{							
		nif = i*nindexfactor;
		nifm1 = (i-1)*nindexfactor;
			
		for(j = nstartfrom; j < ndatatotal2; j++)
		{		
			njf = j*nstates;
			njfm1 = (j-1)*nstates;
			
			for(l = 0; l < nstates; l++)	 //To state
			{		
				hmmptrs(i, j, l) = -1;
				hmmx(i,j,l) = 0;		
				hmmy(i,j,l) = 0;
															
				imax = -1;
				dmax = AMIN;
				dmaxscore = AMIN;
				
				xmax = 0;
				ymax = 0;
			
				el = 0;
				
				nstatetype = statetypes[l];									
				
				for(k = 0; k < nstates; k++) // From state
				{							
					
					vsum = AMIN;
					a = as[k][l];
					
					if(a > AMIN)
					{																			
						switch(nstatetype) 
						{		

							case MATCHSTATEHO:
							
								indexi = i-1;
								indexj = j-1;						
								
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njfm1+nifm1;					
										
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
											
											el = p1s[ndata1][ndata2];		
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
									
										el = p1s[ndata1][ndata2];
									}

									vsum = hmmv[nindex] + a + el;
									
									xsum = xcount+1;
									ysum = ycount+1;																											
								}	
								
							break;

							case OUTPUTSTATEXHO:
							
								indexi = i-1;
								indexj = j;
							
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njf+nifm1;
																																																						
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
											el = qxs[ndata1];		
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
										el = qxs[ndata1]; 
									} 
	
									vsum = hmmv[nindex] + a + el;	
									
									xsum = xcount+1;
									ysum = 0;																																	
								}
								else vsum = AMIN;														
								
							break;
							
							
							case OUTPUTSTATEYHO:
							
								indexi = i;
								indexj = j-1;

								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nif;	
																																																							
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
											el = qys[ndata2];		
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
										el = qys[ndata2];
									} 
								
									vsum = hmmv[nindex] + a + el;												
									
									xsum = 0;
									ysum = ycount+1;	
								}
								else vsum = AMIN;
																				
								
							break;
				
							
						} // switch																																	
																				
					} //a > AMIN


					if(vsum > dmax) 
					{
						dmax = vsum;						
						imax = k;
						
						xmax = xsum;
						ymax = ysum;
																																	
					}					
				} 
	
				if(i > nstart || j > nstart)
				{																							
					mindex = l+j*nstates+i*nindexfactor;				
					
					hmmv[mindex] = dmax;						
					hmmptrs[mindex] = imax;	
					
					hmmx[mindex] = xmax;													
					hmmy[mindex] = ymax;									
				}	
															
			} //for(l = 0; l < nstates; l++)
	
		}			
	}

	for(i = nstartfrom; i < ndatatotal1; i++)
	{		
		nif = i*nindexfactor;
		nifm1 = (i-1)*nindexfactor;
					
		for(j = 0; j < ndatatotal2; j++)
		{	
			njf = j*nstates;
			njfm1 = (j-1)*nstates;
						
														
			for(l = 0; l < nstates; l++)
			{		
				
				hmmptrs(i, j, l) = -1;
							
				hmmx(i,j,l) = 0;		
				hmmy(i,j,l) = 0;
																		
				imax = -1;
				dmax = AMIN;
				dmaxscore = AMIN;
				xmax = 0;
				ymax = 0;

				el = 0;
				
				nstatetype = statetypes[l];									
				
				for(k = 0; k < nstates; k++)
				{																							
				
					vsum = AMIN;
					a = as[k][l];
					
					if(a > AMIN)
					{													
						
						switch(nstatetype) 
						{																		
							case MATCHSTATEHO:
							
								indexi = i-1;
								indexj = j-1;						

								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njfm1+nifm1;					
										
									xcount = hmmx[nindex];														
									ycount = hmmy[nindex];																		
									norderuse = min(xcount, ycount);
									
									if(norderuse > norder) norderuse = norder;
									
									if(indexj > nnearend) norderuse = 0;
									
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
											
											el = p1s[ndata1][ndata2];		
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
									
										el = p1s[ndata1][ndata2];
									}

									vsum = hmmv[nindex] + a + el;
									
									
									
									xsum = xcount+1;
									ysum = ycount+1;																											
								}	
								
							break;

							case OUTPUTSTATEXHO:
							
								indexi = i-1;
								indexj = j;
							
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njf+nifm1;
																																																						
									xcount = hmmx[nindex];																										
															
									norderuse = xcount;
									if(norderuse > norder) norderuse = norder;
										
									if(indexj > nnearend) norderuse = 0;
																																																			
									if(indexi >= norderuse && norderuse > 0) 
									{
										ndata1 = ConvertOrderStringToIndex(data1,indexi,norderuse);	
										
										if(norderuse == 3) el = hmmprobsho3(ndata1, 0, l); 																
										else if(norderuse == 2) el = hmmprobsho2(ndata1, 0, l); 																
										else el = hmmprobsho1(ndata1, 0, l); 	
	
										
										if(el == 1)
										{
											ndata1 = data1[indexi];																		
											el = qxs[ndata1];		
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
										el = qxs[ndata1]; 
									} 
	
									vsum = hmmv[nindex] + a + el;	
									
									xsum = xcount+1;
									ysum = 0;																																	
								}
								else vsum = AMIN;														
								
							break;
							
							
							case OUTPUTSTATEYHO:
							
								indexi = i;
								indexj = j-1;

								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nif;	
																																																							
									ycount = hmmy[nindex];			
															
									norderuse = ycount;
									if(norderuse > norder) norderuse = norder;
										
									if(indexj > nnearend) norderuse = 0;	
																																																						
									if(indexj >= norderuse && norderuse > 0) 
									{																		
										ndata2 = ConvertOrderStringToIndex(data2,indexj,norderuse);	 
										
										if(norderuse == 3) el = hmmprobsho3(ndata2, 0, l);
										else if(norderuse == 2) el = hmmprobsho2(ndata2, 0, l);
										else el = hmmprobsho1(ndata2, 0, l); 	
										
										if(el == 1)
										{
											ndata2 = data2[indexj];																		
											el = qys[ndata2];		
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
										el = qys[ndata2];
									} 
								
									vsum = hmmv[nindex] + a + el;			
									
									
									xsum = 0;
									ysum = ycount+1;	
								}
								else vsum = AMIN;
																				
								
							break;

					
						
						} // switch																																	
																				
					} //a > AMIN
					
					
					if(vsum > dmax)// || (vsum == dmax && k == 2))
					{
						dmax = vsum;
						imax = k;
						xmax = xsum;
						ymax = ysum;	
						
						
												
					}					
				} 
		
				if(i > nstart || j > nstart)
				{																							
					mindex = l+j*nstates+i*nindexfactor;											
															
					hmmv[mindex] = dmax;						
					hmmptrs[mindex] = imax;		
					hmmx[mindex] = xmax;													
					hmmy[mindex] = ymax;
																				
				}	
															
			} //for(l = 0; l < nstates; l++)
		
		}			
	}

	
	vmax = AMIN;								
	for(k = 0; k < nstates; k++)
	{			
		hmmv((ndatalen1),(ndatalen2),k) += as[k][nendstate];
		
		if(hmmv((ndatalen1),(ndatalen2),k) > vmax) 
		{
			vmax = hmmv((ndatalen1),(ndatalen2),k); 
			vstate = k;						
		}	
		
	}	
	
		
	int vstatefinal = vstate;
	
	
	vf.vmax = vmax;
	vf.vmaxdiv = vmax/(double)ndatalen2;
	
	if(boutput == 0) return vf;	
	
	int n1, n2; 
		
	n1 = ndatalen1;
	n2 = ndatalen2;			
		
	vector<int> v;
	vector<int>::iterator it;
	
	vector<double> vvals;
	vector<double>::iterator vit;

	it = v.begin();
	it = v.insert(it, vstate);	
	
	vit = vvals.begin();
	vit = vvals.insert(vit, hmmv((ndatalen1),(ndatalen2),vstate));	
		
	string states = "";				
	string statesoff = "";
		
	int nc = 0;
	
	while(n1 > 1 || n2 > 1 && (n1 > 0 || n2 > 0))
	{
		nstatetype = statetypes[vstate];
			
		vstate = hmmptrs(n1, n2, vstate);
						
		it = v.insert(it, vstate);
		
			
		if(nstatetype == OUTPUTSTATEXHO || nstatetype == OUTPUTSTATEX)
		{	
			n1--;
		}
		else if(nstatetype == OUTPUTSTATEYHO || nstatetype == OUTPUTSTATEY)
		{
			n2--;
		}
		else 
		{
			n1--;
			n2--;
		}
	
		vit = vvals.insert(vit, hmmv(n1,n2,vstate));
		
	}
	
	ninsert = 0;			
	nmatches = 0;
		
	j = 0;
	k = 0;
				

	for (i=0; i<v.size(); i++)
	{		
		vstate = v.at(i);
		
		if(statetypes[vstate] == OUTPUTSTATEXHO || statetypes[vstate] == OUTPUTSTATEX) 
		{			
			vouts1 += sdata1[j];			
			vouts2 += "-";		
			voutinfo += "-";
														
			j++;	
				
			states += "x";				
			statesoff += '=';				
				
			ninsert++;	
								
		}
		else if(statetypes[vstate] == OUTPUTSTATEYHO || statetypes[vstate] == OUTPUTSTATEY) 
		{	
			vouts1 += "-";			
			vouts2 += sdata2[k];	
			voutinfo += sseqinfo[k];						
				
			k++;			
				
			states += "y";						
			statesoff += '=';
				
			ninsert++;											
				
		}
		else if(statetypes[vstate] == MATCHSTATEHO || statetypes[vstate] == MATCHSTATE) 
		{
			vouts1 += sdata1[j];
			vouts2 += sdata2[k];	
			voutinfo += sseqinfo[k];	
				
			j++;
			k++;	
								
			states += "m";
			statesoff += '=';
					
			nmatches++;
				
		}
	}
		
	
	vstate = vstatefinal;
	i = vvals.size()-1;
	j = vouts1.length()-1;
	k = vouts2.length()-1;

	
	while(k > 0 && vouts2[k] == '-')
	{
		
		vouts1.erase(j,1);
		vouts2.erase(k,1);
		voutinfo.erase(j,1);
		
		j--;
		k--;
		i--;
		
		vf.vmax = vvals.at(i) + as[v.at(i)][nendstate];
		ndatalen1--;								
	}
	
	
	k = vouts2.length()-1;
		
	vf.vmaxdiv = vf.vmax/(double)ndatalen2;
	vf.vouts1 = vouts1;
	vf.vouts2 = vouts2;
	vf.voutinfo = voutinfo;
	
	vf.nreflength = ndatalen1;			
	vf.ninsert = ninsert;
	
	return vf;
}



vitstruct Viterbi(int* data1, int* data2, string sdata1, string sdata2, string sseqinfo, int ndatalen1, 
									   int ndatalen2, int nstates, double* as[], double* aslogs[], int* statetypes, 
									   int nbeginstate, int nendstate, double* qxs, double* qys, double** p1s, double* hmmprobsho,
									   int nbases, int boutput, double* hmmv, int* hmmptrs, int nstartfrom, 
									   int nindexfactor, int norder, int noutputsho, int* hmmx, int* hmmy,
									   int nnearendbufferlen)			
{
	
	int i;
	int j; 
	int k; 
	int l;
	int indexi; 
	int indexj; 
	int ndata1; 
	int ndata2;
	int nstatetype; 
	int imax; 
	int vstate;

	int nlen; 
	int nstart;
	int indexm; 
	int indexn;
	int ndatatotal1; 
	int ndatatotal2; 
	int nindex; 
	int mindex;
	int ninsert;
	int nmatches;
	
	int nif; 
	int nifm1; 
	int njf; 
	int njfm1;	
	int norderuse;
	int nnearend;
	int xcount; 
	int ycount; 
	int xsum; 
	int ysum; 
	int xmax; 
	int ymax;
	
	
	
	double dscore; 
	double dmaxscore;
	double a; 
	double el; 
	double dmax; 
	double vsum;  
	double vmax;		
		
	string s1; 
	string s2; 
	string vouts1;
	string vouts2; 
	string voutinfo;
	
	vitstruct vf;

	vouts1 = ""; 
	vouts2 = "";
	voutinfo = "";

	nstart = 0;
		
	//Add one instance before the start of the data (for the begin state)	
	ndatatotal1 = ndatalen1+1; 
	ndatatotal2 = ndatalen2+1;		

	nnearend = ndatalen1 - nnearendbufferlen;
	

	xcount = 0;
	ycount = 0;
	

	if(nstartfrom == 0)
	{
		for(i = 0; i < nstates; i++)
		{		
			hmmv(0, 0, i) = AMIN;	

			hmmptrs(0, 0, i) = -1;
			hmmptrs(0, 1, i) = -1;
			hmmptrs(1, 0, i) = -1;
			hmmptrs(1, 1, i) = -1;

			// Last index in b array is the last data index so add from 				
			//hmmb(ndatatotal1-1, ndatatotal2-1, i) = as[i][nendstate];
			if(i == nbeginstate)
			{	
				hmmv(0, 0, i) = 0;
			}
		}
	}
		

	for(i = 0; i < nstartfrom; i++)
	{							
		nif = i*nindexfactor;
		nifm1 = (i-1)*nindexfactor;
			
		for(j = nstartfrom; j < ndatatotal2; j++)
		{		
			njf = j*nstates;
			njfm1 = (j-1)*nstates;
			
			for(l = 0; l < nstates; l++)   //To state
			{		
				hmmptrs(i, j, l) = -1;
				hmmx(i,j,l) = 0;		
				hmmy(i,j,l) = 0;
															
				imax = -1;
				dmax = AMIN;
				dmaxscore = AMIN;

				
				xmax = 0;
				ymax = 0;
			
				el = 0;
				
				nstatetype = statetypes[l];									
				
				for(k = 0; k < nstates; k++) //From state
				{												
					vsum = AMIN;
					a = as[k][l];					
					
					if(a > AMIN)
					{																			
						switch(nstatetype) 
						{		
						
							case OUTPUTSTATEX:
							
								indexi = i-1;
								indexj = j;			
																											
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njf+nifm1;																										
									ndata1 = data1[indexi];																																					
									el = qxs[ndata1];																																																										
									vsum = hmmv[nindex] + a + el;												
								
								}
								else vsum = AMIN;														
								
							break;
						
							case OUTPUTSTATEY:
																						
								indexi = i;
								indexj = j-1;
																			
								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nif;													
									ndata2 = data2[indexj];																											
									el = qys[ndata2];																																																
									vsum = hmmv[nindex] + a + el;																																																																																											
								}
								else vsum = AMIN;

							break;
							
							case MATCHSTATE:
														
								indexi = i-1;
								indexj = j-1;	
																	
								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nifm1;
									ndata1 = data1[indexi];
									ndata2 = data2[indexj];																																						
									el = p1s[ndata1][ndata2];																																																				
									vsum = hmmv[nindex] + a + el;

								}	
							
							break;
							
							
							case MATCHSTATEHO:
							
								indexi = i-1;
								indexj = j-1;						
						
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njfm1+nifm1;					
										
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
									
										el = p1s[ndata1][ndata2];
									}

									vsum = hmmv[nindex] + a + el;
									
									xsum = xcount+1;
									ysum = ycount+1;																											
								}	
								
							break;

							case OUTPUTSTATEXHO:
							
								indexi = i-1;
								indexj = j;

								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njf+nifm1;
																																																						
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
										el = qxs[ndata1]; 
									} 
	
									vsum = hmmv[nindex] + a + el;	
									
									xsum = xcount+1;
									ysum = 0;																																	
								}
								else vsum = AMIN;														
								
							break;
							
							
							case OUTPUTSTATEYHO:
							
								indexi = i;
								indexj = j-1;								
						
								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nif;	
																																																							
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
										el = qys[ndata2];
									} 
								
									vsum = hmmv[nindex] + a + el;												
									
									xsum = 0;
									ysum = ycount+1;	
								}
								else vsum = AMIN;
																				
								
							break;
							
							
						} // switch																																	
																				
					} //a > AMIN



					if(vsum > dmax)
					{
						dmax = vsum;
						imax = k;
						
						xmax = xsum;
						ymax = ysum;
																																	
					}					
				} 
	
				if(i > nstart || j > nstart)
				{																							
					mindex = l+j*nstates+i*nindexfactor;				
				
					hmmv[mindex] = dmax;						
					hmmptrs[mindex] = imax;	
					
					hmmx[mindex] = xmax;													
					hmmy[mindex] = ymax;									
				}	
															
			} //for(l = 0; l < nstates; l++)
	
		}			
	}


	for(i = nstartfrom; i < ndatatotal1; i++)
	{		
		nif = i*nindexfactor;
		nifm1 = (i-1)*nindexfactor;
					
		for(j = 0; j < ndatatotal2; j++)
		{	
			njf = j*nstates;
			njfm1 = (j-1)*nstates;
						
														
			for(l = 0; l < nstates; l++)
			{		
				
				hmmptrs(i, j, l) = -1;
							
				hmmx(i,j,l) = 0;		
				hmmy(i,j,l) = 0;
																		
				imax = -1;
				dmax = AMIN;
				dmaxscore = AMIN;
				xmax = 0;
				ymax = 0;

				el = 0;
				
				nstatetype = statetypes[l];									
				
				for(k = 0; k < nstates; k++)
				{																							
				
					vsum = AMIN;
					a = as[k][l];
					
					if(a > AMIN)
					{													
						
						switch(nstatetype) // state we are going to 
						{			
							
							case OUTPUTSTATEX:
							
								indexi = i-1;
								indexj = j;			
																
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njf+nifm1;																										
									ndata1 = data1[indexi];																																					
									el = qxs[ndata1];																																																										
									vsum = hmmv[nindex] + a + el;												
								
								}
								else vsum = AMIN;														
								
							break;
						
							case OUTPUTSTATEY:
																						
								indexi = i;
								indexj = j-1;

								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nif;													
									ndata2 = data2[indexj];																											
									el = qys[ndata2];																																																
									vsum = hmmv[nindex] + a + el;																																																																																											
								}
								else vsum = AMIN;

							break;
							
							case MATCHSTATE:
														
								indexi = i-1;
								indexj = j-1;	
																	
								if(indexi >= 0 && indexj >= 0)
								{		
									nindex = k+njfm1+nifm1;
									ndata1 = data1[indexi];
									ndata2 = data2[indexj];																																						
									el = p1s[ndata1][ndata2];																																																				
									vsum = hmmv[nindex] + a + el;
								}	
							
							break;
									
											
							case MATCHSTATEHO:
							
								indexi = i-1;
								indexj = j-1;						
	
								if(indexi >= 0 && indexj >= 0)
								{				
									nindex = k+njfm1+nifm1;	
										
									xcount = hmmx[nindex];																								
									ycount = hmmy[nindex];																				
															
									norderuse = min(xcount, ycount);
									if(norderuse > norder) norderuse = norder;
				
									if(indexj > nnearend) norderuse = 0;
									
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
																				
										el = p1s[ndata1][ndata2];										
									}

									vsum = hmmv[nindex] + a + el;
									
									xsum = xcount+1;
									ysum = ycount+1;	
																																				
								}	
								
							break;
							
							case OUTPUTSTATEXHO:
							
								indexi = i-1;
								indexj = j;
							
								if(indexi >= 0 && indexj >= 0)
								{
									nindex = k+njf+nifm1;																																													
									xcount = hmmx[nindex];																										
															
									norderuse = xcount;
									if(norderuse > norder) norderuse = norder;
									
									if(indexj > nnearend) norderuse = 0;
																																																		
									if(indexi >= norderuse && norderuse > 0) 
									{
										ndata1 = ConvertOrderStringToIndex(data1,indexi,norderuse);										
										el = hmmprobsho(ndata1, 0, l); 
									}
									else
									{
										ndata1 = data1[indexi];	
										el = qxs[ndata1];
									} 
	
									vsum = hmmv[nindex] + a + el;	
									
									xsum = xcount+1;
									ysum = 0;																																
								}
								else vsum = AMIN;														
								
							break;
							
							
							case OUTPUTSTATEYHO:
							
								indexi = i;
								indexj = j-1;

								if(indexi >= 0 && indexj >= 0)
								{	
									nindex = k+njfm1+nif;			
																																																						
									ycount = hmmy[nindex];																		
									norderuse = ycount;
									if(norderuse > norder) norderuse = norder;
	
									if(indexj > nnearend) norderuse = 0;
																																																							
									if(indexj >= norderuse && norderuse > 0) 
									{																				
										ndata2 = ConvertOrderStringToIndex(data2,indexj,norderuse);	
										el = hmmprobsho(ndata2, 0, l); 
									}
									else
									{
										ndata2 = data2[indexj];	
										el = qys[ndata2];
									} 
								
									vsum = hmmv[nindex] + a + el;
									
									xsum = 0;
									ysum = ycount+1;												
								}
								else vsum = AMIN;														
								
							break;
					
						
						} // switch																																	
																				
					} //a > AMIN


					if(vsum > dmax)
					{
						dmax = vsum;						
						imax = k;
						xmax = xsum;
						ymax = ysum;															
												
					}					
				} 
		
				if(i > nstart || j > nstart)
				{																							
					mindex = l+j*nstates+i*nindexfactor;											
															
					hmmv[mindex] = dmax;						
					hmmptrs[mindex] = imax;		
					hmmx[mindex] = xmax;													
					hmmy[mindex] = ymax;															
				}	
															
			} //for(l = 0; l < nstates; l++)
		
		}			
	}
	
	
	// Get best state at end of viterbi 
	vmax = AMIN;								
	for(k = 0; k < nstates; k++)
	{			
		hmmv((ndatalen1),(ndatalen2),k) += as[k][nendstate];
		
		if(hmmv((ndatalen1),(ndatalen2),k) > vmax) 
		{
			vmax = hmmv((ndatalen1),(ndatalen2),k); 
			vstate = k;						
		}	

	}	
		
	int vstatefinal = vstate;
	
	vf.vmax = vmax;
	vf.vmaxdiv = vmax/(double)ndatalen2;
	
	if(boutput == 0) return vf;	
	
	int n1, n2; 
		
	n1 = ndatalen1;
	n2 = ndatalen2;			
		
	vector<int> v;
	vector<int>::iterator it;
	
	vector<double> vvals;
	vector<double>::iterator vit;

	it = v.begin();
	it = v.insert(it, vstate);	
	
	vit = vvals.begin();
	vit = vvals.insert(vit, hmmv((ndatalen1),(ndatalen2),vstate));	
		
	string states = "";				
	string statesoff = "";
		
	int nc = 0;
	
	while(n1 > 1 || n2 > 1 && (n1 > 0 || n2 > 0))
	{
		nstatetype = statetypes[vstate];
		
		vstate = hmmptrs(n1, n2, vstate);
						
		it = v.insert(it, vstate);
		
			
		if(nstatetype == OUTPUTSTATEXHO || nstatetype == OUTPUTSTATEX)
		{	
			n1--;
		}
		else if(nstatetype == OUTPUTSTATEYHO || nstatetype == OUTPUTSTATEY)
		{
			n2--;
		}
		else 
		{
			n1--;
			n2--;
		}
	
		vit = vvals.insert(vit, hmmv(n1,n2,vstate));
		
	}
	
	ninsert = 0;	
	nmatches = 0;
		
	j = 0;
	k = 0;
				
		
	for (i=0; i<v.size(); i++)
	{		
		vstate = v.at(i);
		
		if(statetypes[vstate] == OUTPUTSTATEXHO || statetypes[vstate] == OUTPUTSTATEX) 
		{			
			vouts1 += sdata1[j];			
			vouts2 += "-";		
			voutinfo += "-";
														
			j++;	
				
			states += "x";				
			statesoff += '=';				
				
			ninsert++;	
								
		}
		else if(statetypes[vstate] == OUTPUTSTATEYHO || statetypes[vstate] == OUTPUTSTATEY) 
		{	
			vouts1 += "-";			
			vouts2 += sdata2[k];	
			voutinfo += sseqinfo[k];						
				
			k++;			
				
			states += "y";						
			statesoff += '=';
				
			ninsert++;											
				
		}
		else if(statetypes[vstate] == MATCHSTATEHO || statetypes[vstate] == MATCHSTATE) 
		{
			vouts1 += sdata1[j];
			vouts2 += sdata2[k];	
			voutinfo += sseqinfo[k];	
				
			j++;
			k++;	
								
			states += "m";
			statesoff += '=';
					
			nmatches++;
				
		}
	}
		
	vstate = vstatefinal;
	i = vvals.size()-1;
	j = vouts1.length()-1;
	k = vouts2.length()-1;
	
	while(k > 0 && vouts2[k] == '-')
	{
		
		vouts1.erase(j,1);
		vouts2.erase(k,1);
		voutinfo.erase(j,1);
		
		j--;
		k--;
		i--;
		
		vf.vmax = vvals.at(i) + as[v.at(i)][nendstate];
		ndatalen1--;									
	}	
		
	vf.vmaxdiv = vmax/(double)ndatalen2;
			
	vf.vouts1 = vouts1;
	vf.vouts2 = vouts2;
	vf.voutinfo = voutinfo;
		
	vf.nreflength = ndatalen1;	
	vf.ninsert = ninsert;

	return vf;
}









charcountstruct CountChars(string s)
{
	int i;
	charcountstruct cs;
	
	cs.nA = 0;
	cs.nC = 0;
	cs.nG = 0;
	cs.nT = 0;
	
	for(i = 0; i < s.length(); i++)
	{
		switch(s[i])
		{
			case 'A':
				cs.nA++;
			break;
			case 'C':
				cs.nC++;
			break;
			case 'G':
				cs.nG++;
			break;
			case 'T':
				cs.nT++;
			break;			
		}
	}
	
	return cs;
}



void quickSort(unsigned short int arr[], unsigned int arr2[], int left, int right) 
{
      int i = left, j = right;
      int tmp;
      int n = (int)(((float)left + (float)right) / 2.0);
      int pivot = arr[n];
 
      /* partition */
      while (i <= j) 
      {
            while (arr[i] < pivot) i++;
            while (arr[j] > pivot) j--;
            if (i <= j) 
            {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  
                  tmp = arr2[i];
                  arr2[i] = arr2[j];
                  arr2[j] = tmp;
                  
                  i++;
                  j--;
            }
      }
 
      /* recursion */
      if (left < j)
      {
         quickSort(arr, arr2, left, j);
	  }
      if (i < right)
      {
         quickSort(arr, arr2, i, right);
	  }
}



void quickSortStrings(string* items, unsigned short int arr[], unsigned int arr2[], int left, int right)
{
  int i, j;
  int tmp;
  //char *x;
  //char temp[24];
  string x;
  string temp;

  i = left;
  j = right;
  x = items[(left+right)/2];

  do {
    while((items[i].compare(x) < 0) && (i < right)) {
       i++;
    }
    while((items[j].compare(x) > 0) && (j > left)) {
        j--;
    }
    if(i <= j) {	  	  
	  temp = items[i];
	  items[i] = items[j];
	  items[j] = temp;
	  /*
      strcpy(temp, items[i]);
      strcpy(items[i], items[j]);
      strcpy(items[j], temp);*/
      
      tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
      
      tmp = arr2[i];
      arr2[i] = arr2[j];
      arr2[j] = tmp;
      
      i++;
      j--;
   }
  } while(i <= j);

  if(left < j) {
     quickSortStrings(items, arr, arr2, left, j);
  }
  if(i < right) {
     quickSortStrings(items, arr, arr2, i, right);
  }
}



unsigned short int ConvertIntsToID(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
{
	return 16384*i1+4096*i2+1024*i3+256*i4+64*i5+16*i6+4*i7+i8;
}


unsigned int ConvertIntsToID(unsigned int* nvals, int nindex, int nidlength)
{
	int i;
	int nid;
	
	nid = 0;
	
	for(i = 0; i < nidlength; i++)
	{
		nid += (int)pow(4,nidlength-1-i)*nvals[nindex+i];
	}
	
	return nid;
}


unsigned int ConvertIntsToID(unsigned short int* nvals, int nindex, int nidlength)
{
	int i;
	int nid;
	
	nid = 0;
	
	for(i = 0; i < nidlength; i++)
	{
		nid += (int)pow(4,nidlength-1-i)*nvals[nindex+i];
	}
	
	return nid;
}


int ConvertDoubleScoreCharScore(double d)
{
	int i;
	
	//i = (int)(((d+1)/5.0)*256);
	//if(i > 256) i = 256;
	
	i = (int)(((d+1)/5.0)*126);
	if(i > 126) i = 126;
	return i;
}




void InitializeHMM(int* statetypes, double** as, double** aslogs, double** p1s, double* qxs, double* qys, 
				   ifstream &infilehmm, int nbases, int nstates, int noutputsho, double* hmmprobsho, int norder)
{
	int i;
	int j;
	int n;
	int n1;
	int n2;
	
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int np4;
	int n4seq1;
	int n4seq2;
	
	string s;
	string line;
	
	int nmaxoffset = 0;
	
	/////////////////////////////////////////////////////////////////////////HMM initialization
	
	//nstatetypes
	getline(infilehmm,line);	
	n1 = line.find("st;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: st not found";
		return;
	}
	n1 += 3;
						
	for(i = 0; i < nstates; i++)
	{
		n2 = line.find(",", n1);
		s = line.substr(n1, n2-n1);
		
		statetypes[i] = atoi(s.c_str());					
		
		n1 = n2+1;
	}

	//as
	getline(infilehmm,line);	
	n1 = line.find("as;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: as not found";
		return;
	}
	n1 += 3;
						
	for(i = 0; i < nstates; i++)
	{
		for(j = 0; j < nstates; j++)
		{
			n2 = line.find(",", n1);
			s = line.substr(n1, n2-n1);
		
			as[i][j] = atof(s.c_str());					
		
			n1 = n2+1;
		}
	}
	
	
	//load order: qxs, qys, qx2s, qy2s, qoffs, p1s, p2s, poffs, indexes	
	
	
	//qxs
	getline(infilehmm,line);	
	n1 = line.find("x;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: qxs not found";
		return;
	}
	n1 += 2; 
	
	for(i = 0; i < nbases; i++)
	{
		n2 = line.find(",", n1);
		s = line.substr(n1, n2-n1);
				
		qxs[i] = atof(s.c_str());

		n1 = n2+1;
	}
	
	
	//qys
	getline(infilehmm,line);	
	n1 = line.find("y;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: qys not found";
		return;
	}
	n1 += 2; 
	
	for(i = 0; i < nbases; i++)
	{
		n2 = line.find(",", n1);
		s = line.substr(n1, n2-n1);
				
		qys[i] = atof(s.c_str());

		n1 = n2+1;
	}
	

	//p1s
	getline(infilehmm,line);	
	n1 = line.find("p0;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: p0s not found";
		return;
	}
	n1 += 3; 
						
	for(i = 0; i < nbases; i++)
	{
		for(j = 0; j < nbases; j++)
		{
			n2 = line.find(",", n1);
			s = line.substr(n1, n2-n1);
				
			p1s[i][j] = atof(s.c_str());			

			n1 = n2+1;
		}
	}
	
	if(norder > 0)
	{
		getline(infilehmm,line);	
		ReadOutputProbs(line, hmmprobsho, norder, 1, noutputsho);
	
		getline(infilehmm,line);
		ReadOutputProbs(line, hmmprobsho, norder, 3, noutputsho);
	
		getline(infilehmm,line);	
		ReadAlignProbs(line, hmmprobsho, norder, 2, noutputsho);
	}
}






void InitializeVariedOrderHMM(int* statetypes, double** as, double** aslogs, double** p1s, double* qxs, double* qys, ifstream &infilehmm, int nbases, int nstates, int noutputsho1, int noutputsho2, int noutputsho3, double* hmmprobsho1, double* hmmprobsho2, double* hmmprobsho3)
{
	int i;
	int j;
	int n;
	int n1;
	int n2;
	
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;
	int np4;
	int n4seq1;
	int n4seq2;
	
	string s;
	string line;
	
	int nmaxoffset = 0;
	
	/////////////////////////////////////////////////////////////////////////HMM initialization


	//nstatetypes
	getline(infilehmm,line);	
	n1 = line.find("st;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: st not found";
		return;
	}
	n1 += 3;
						
	for(i = 0; i < nstates; i++)
	{
		n2 = line.find(",", n1);
		s = line.substr(n1, n2-n1);
		
		statetypes[i] = atoi(s.c_str());					
		
		n1 = n2+1;
	}

	//as
	getline(infilehmm,line);	
	n1 = line.find("as;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: as not found";
		return;
	}
	n1 += 3;
						
	for(i = 0; i < nstates; i++)
	{
		for(j = 0; j < nstates; j++)
		{
			n2 = line.find(",", n1);
			s = line.substr(n1, n2-n1);
		
			as[i][j] = atof(s.c_str());					
		
			n1 = n2+1;
		}
	}
	
	
	//load order: qxs, qys, qx2s, qy2s, qoffs, p1s, p2s, poffs, indexes	
	
	
	//qxs
	getline(infilehmm,line);	
	n1 = line.find("x;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: qxs not found";
		return;
	}
	n1 += 2; 
	
	for(i = 0; i < nbases; i++)
	{
		n2 = line.find(",", n1);
		s = line.substr(n1, n2-n1);
				
		qxs[i] = atof(s.c_str());

		n1 = n2+1;
	}
	
	
	//qys
	getline(infilehmm,line);	
	n1 = line.find("y;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: qys not found";
		return;
	}
	n1 += 2; 
	
	for(i = 0; i < nbases; i++)
	{
		n2 = line.find(",", n1);
		s = line.substr(n1, n2-n1);
				
		qys[i] = atof(s.c_str());

		n1 = n2+1;
	}
	

	//p1s
	getline(infilehmm,line);	
	n1 = line.find("p0;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: p0s not found";
		return;
	}
	n1 += 3; 
						
	for(i = 0; i < nbases; i++)
	{
		for(j = 0; j < nbases; j++)
		{
			n2 = line.find(",", n1);
			s = line.substr(n1, n2-n1);
				
			p1s[i][j] = atof(s.c_str());			

			n1 = n2+1;
		}
	}
	

	getline(infilehmm,line);		
	ReadOutputProbs(line, hmmprobsho1, 1, 1, noutputsho1);
	
	getline(infilehmm,line);	
	ReadOutputProbs(line, hmmprobsho1, 1, 3, noutputsho1);
	
	getline(infilehmm,line);		
	ReadAlignProbs(line, hmmprobsho1, 1, 2, noutputsho1);
	
	
	getline(infilehmm,line);	
	ReadOutputProbs(line, hmmprobsho2, 2, 1, noutputsho2);
	
	getline(infilehmm,line);
	ReadOutputProbs(line, hmmprobsho2, 2, 3, noutputsho2);
	
	getline(infilehmm,line);	
	ReadAlignProbs(line, hmmprobsho2, 2, 2, noutputsho2);
	
	
	getline(infilehmm,line);	
	ReadOutputProbs(line, hmmprobsho3, 3, 1, noutputsho3);
	
	getline(infilehmm,line);
	ReadOutputProbs(line, hmmprobsho3, 3, 3, noutputsho3);
	
	getline(infilehmm,line);	
	ReadAlignProbs(line, hmmprobsho3, 3, 2, noutputsho3);
	
	
}





void LoadHMMBaseValues(ifstream &infilehmm, int &nstates, int &noutputs, 
					   int &nbeginstate, int &nendstate, int &nbases,
					   int &norder, int &bVariableOrder)
{
	int n1;
	int n2;
	
	string s;
	string line;
	
	// nstates
	getline(infilehmm,line);
	n1 = line.find("n1;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n1 not found\n";
		return;
	}
	n1 += 3; //add length of "n1;" to get to next starting point
	n2 = line.find(";", n1+1);			
	s = line.substr(n1, n2-n1);
	nstates = atoi(s.c_str());		
	
	//noutputs
	getline(infilehmm,line);	
	n1 = line.find("n2;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n2 not found\n";
		return;
	}
	n1 += 3; 
	n2 = line.find(";", n1+1);		
	s = line.substr(n1, n2-n1);
	noutputs = atoi(s.c_str());
	
	
	//nbeginstate
	getline(infilehmm,line);	
	n1 = line.find("n3;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n3 not found\n";
		return;
	}
	n1 += 3; 
	n2 = line.find(";", n1+1);		
	s = line.substr(n1, n2-n1);
	nbeginstate = atoi(s.c_str());
	
	
	//nendstate
	getline(infilehmm,line);	
	n1 = line.find("n4;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n4 not found\n";
		return;
	}
	n1 += 3; 
	n2 = line.find(";", n1+1);		
	s = line.substr(n1, n2-n1);
	nendstate = atoi(s.c_str());
	
	
	//nbases
	getline(infilehmm,line);	
	n1 = line.find("n5;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n5 not found\n";
		return;
	}
	n1 += 3; 
	n2 = line.find(";", n1+1);		
	s = line.substr(n1, n2-n1);
	nbases = atoi(s.c_str());
	
	//norder
	getline(infilehmm,line);	
	n1 = line.find("n6;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n6 not found\n";
		return;
	}
	n1 += 3; 
	n2 = line.find(";", n1+1);		
	s = line.substr(n1, n2-n1);
	norder = atoi(s.c_str());
	
	//bvariableorder
	getline(infilehmm,line);	
	n1 = line.find("n7;", 0);
	if(n1 == -1) 
	{
		cout << "Error loading hmm variables: n7 not found\n";
		return;
	}
	n1 += 3; 
	n2 = line.find(";", n1+1);		
	s = line.substr(n1, n2-n1);
	bVariableOrder = atoi(s.c_str());
}



void CalcAlignScores(string shmmfile, int binfo, int i1from, int i2from, int i3from, int ncurrentthread)	
{
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int i7;
	int i8;		
	int j1;
	int j2;
	int j3;
	int j4;
	int j5;
	int j6;
	int j7;
	int j8;		
	int n;	
	int n1;
	int n2;	
	int i1to;
	int i2to;
	int i3to;
	
	int nid1;
	int nid2;
	int ncharscore;
	int nlinemod;	
	int nsum;	
	int nmaxchar;
	int nreadlength;
	int noutputsho;
	int noutputsho1;
	int noutputsho2;
	int noutputsho3;
	
	int nbases;
	int nstates; 
	int noutputs; 
	int nbeginstate; 
	int nendstate;
	int norder;
	int bVariableOrder;	
	int nindexfactorfull;	
	
	double d;
	double dscore;
	
	string s;
	string line;
	string soutalignscores;	
	
	ofstream foutalignscores;
	
	vitstruct vf;	

	/////////////////Initialize HMM		
	ifstream infilehmm;
	infilehmm.open(shmmfile.c_str(), ifstream::in);	

	
	LoadHMMBaseValues(infilehmm, nstates, noutputs, 
					   nbeginstate, nendstate, nbases,
					   norder, bVariableOrder);

	//cout << "HMM order = " << norder-1 << "\n";
	//if(bVariableOrder == 0) cout << "Variable order = No\n\n";
	//else cout << "Variable order = Yes\n\n";
	
	int* statetypes = new int[nstates];
	double** as = CreateDoubleArray(nstates, nstates);
	double** aslogs = CreateDoubleArray(nstates, nstates);		
	double** p1s = CreateDoubleArray(nbases, nbases);
	double*  qxs = new double[nbases];
	double*  qys = new double[nbases];	
	
	noutputsho = (int)pow(4, norder+1);
	double* hmmprobsho = new double[noutputsho*noutputsho*nstates];

	noutputsho1 = (int)pow(4, 2);
	noutputsho2 = (int)pow(4, 3);
	noutputsho3 = (int)pow(4, 4);
	
	double* hmmprobsho1 = new double[noutputsho1*noutputsho1*nstates];
	double* hmmprobsho2 = new double[noutputsho2*noutputsho2*nstates];
	double* hmmprobsho3 = new double[noutputsho3*noutputsho3*nstates];
	
	
	if(bVariableOrder == 1) 
	{
		InitializeVariedOrderHMM(statetypes, as, aslogs, p1s, qxs, qys, infilehmm, nbases, nstates, noutputsho1, noutputsho2, noutputsho3, hmmprobsho1, hmmprobsho2, hmmprobsho3);
	}
	else 
	{
		InitializeHMM(statetypes, as, aslogs, p1s, qxs, qys, infilehmm, nbases, nstates, noutputsho, hmmprobsho, norder);
	}
	
	infilehmm.close();
	
	
	//////////////Initialize variables for viterbi
	
	nreadlength = 8;
	
	int* readshort = new int[nreadlength];	
	int* refshort = new int[nreadlength];
	
	int hmmnlen = (nreadlength+1)*(nreadlength+1)*nstates;	
		
	int* hmmx = new int[hmmnlen];	
	int* hmmy = new int[hmmnlen];	
		
	double* hmmv = new double[hmmnlen];		
	int* hmmptrs = new int[hmmnlen];			
	
	
	//////////////Initialize variables alignemnt score calculations
			 	
		
	nindexfactorfull = (nreadlength+1)*nstates;
	
	nmaxchar = 4;		
	nlinemod = 0;
	
	if(i1from == -1)
	{
		i1from = 0;
		i1to = 4;
		i2from = 0;
		i2to = 4;
		i3from = 0;
		i3to = 4;
	}
	else
	{
		i1to    = i1from+1;
		i2to    = i2from+1;
		i3to    = i3from+1;
	}
	
	//////////////Calc HMM alignments
	
	
	if(binfo == 1)
	{
		cout << "Examples of alignments and their scores using " << shmmfile << " paramters:\n";
		
		for(i1 = i1from; i1 < i1to; i1++) 
			for(i2 = i2from; i2 < i2to; i2++)
				for(i3 = i3from; i3 < i3to; i3++)
					for(i4 = 0; i4 < nmaxchar; i4++)
						for(i5 = 0; i5 < nmaxchar; i5++)
							for(i6 = 0; i6 < nmaxchar; i6++)
								for(i7 = 0; i7 < nmaxchar; i7++)
									for(i8 = 0; i8 < nmaxchar; i8++)
										for(j1 = 0; j1 < nmaxchar; j1++) 
											for(j2 = 0; j2 < nmaxchar; j2++)
												for(j3 = 0; j3 < nmaxchar; j3++)
													for(j4 = 0; j4 < nmaxchar; j4++)
														for(j5 = 0; j5 < nmaxchar; j5++)
															for(j6 = 0; j6 < nmaxchar; j6++)
																for(j7 = 0; j7 < nmaxchar; j7++)
																	for(j8 = 0; j8 < nmaxchar; j8++)
																	{
																		nsum = i2+i3+i4+i5+i6+i7+i8+j2+j3+j4+j5+j6+j7+j8;
																	
																		if(nsum == 0 || (i2 == i1 && i3 == i1 && i4 == i1 && i5 == i1 && i6 == i1 && i7 == i1 && i8 == i1 && j1 == i1 && j2 == i1 && j3 == i1 && j4 == i1 && j5 == i1 && j6 == i1 && j7 == i1 && j8 == i1) || (i1 == 0 && i2 == 1 && i3 == 2 && i4 == 3 && i5 == 0 && i6 == 1 && i7 == 2 && i8 == 3
																		&& j1 == i1 && j2 == i2 && j3 == i3 && j4 == i4 && j5 == i5 && j6 == i6 && j7 == i7 && j8 == i8))
																		{
																			refshort[0] = i1;
																			refshort[1] = i2;
																			refshort[2] = i3;
																			refshort[3] = i4;
																			refshort[4] = i5;
																			refshort[5] = i6;
																			refshort[6] = i7;
																			refshort[7] = i8;
																			
																			readshort[0] = j1;
																			readshort[1] = j2;
																			readshort[2] = j3;
																			readshort[3] = j4;
																			readshort[4] = j5;
																			readshort[5] = j6;
																			readshort[6] = j7;
																			readshort[7] = j8;
																			
																			/*
																			if(bVariableOrder == 1)
																			{
																				vf = ViterbiVariedOrderHMM(refshort, readshort, "", "", "", nreadlength, nreadlength,
																											   nstates, as, aslogs, statetypes, 
																											   nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
																											   nbases, 0, 0, 0,
																											   hmmv, hmmvratio, hmmptrs, 0, nindexfactorfull,-1000, 
																											   norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, 0);	
																			}
																			else
																			{
																				vf = Viterbi(refshort, readshort, "", "", "", nreadlength, nreadlength,
																										 nstates, as, aslogs, statetypes, 
																										 nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
																										 nbases, 0, 0, 0,
																										 hmmv, hmmvratio, hmmptrs, 0, nindexfactorfull,-1000, norder, noutputsho, hmmx, hmmy);	
																			}
																			*/
																			
																			dscore = -vf.vmax/(float)nreadlength;
																																																									
																			ncharscore = ConvertDoubleScoreCharScore(dscore);
																			

																			if(i1 == 0 && j1 == 0 && nsum == 0)
																			{
																				cout << "AAAAAAAA" << "\n";
																				cout << "AAAAAAAA" << "\n";
																				cout << dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 1 && i2 == i1 && i3 == i1 && i4 == i1 && i5 == i1 && i6 == i1 && i7 == i1 && i8 == i1  && j1 == i1 && j2 == i1 && j3 == i1 && j4 == i1 && j5 == i1 && j6 == i1 && j7 == i1 && j8 == i1)
																			{
																				cout << "CCCCCCCC" << "\n";
																				cout << "CCCCCCCC" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 2 && i2 == i1 && i3 == i1 && i4 == i1 && i5 == i1 && i6 == i1 && i7 == i1 && i8 == i1  && j1 == i1 && j2 == i1 && j3 == i1 && j4 == i1 && j5 == i1 && j6 == i1 && j7 == i1 && j8 == i1)
																			{	
																				cout << "GGGGGGGG" << "\n";
																				cout << "GGGGGGGG" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 3 && j1 == 3 && i2 == i1 && i3 == i1 && i4 == i1 && i5 == i1 && i6 == i1 && i7 == i1 && i8 == i1  && j1 == i1 && j2 == i1 && j3 == i1 && j4 == i1 && j5 == i1 && j6 == i1 && j7 == i1 && j8 == i1)
																			{
																				cout << "TTTTTTTT" << "\n";
																				cout << "TTTTTTTT" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 0 && i2 == 1 && i3 == 2 && i4 == 3 && i5 == 0 && i6 == 1 && i7 == 2 && i8 == 3 && j1 == i1 && j2 == i2 && j3 == i3 && j4 == i4 && j5 == i5 && j6 == i6 && j7 == i7 && j8 == i8)
																			{
																				cout << "ACGTACGT" << "\n";
																				cout << "ACGTACGT" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			}
																			else if(i1 == 1 && j1 == 0 && nsum == 0)
																			{
																				cout << "CAAAAAAA" << "\n";
																				cout << "AAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 2 && j1 == 0 && nsum == 0)
																			{
																				cout << "GAAAAAAA" << "\n";
																				cout << "AAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 3 && j1 == 0 && nsum == 0)
																			{
																				cout << "TAAAAAAA" << "\n";
																				cout << "AAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 0 && j1 == 1 && nsum == 0)
																			{
																				cout << "AAAAAAAA" << "\n";
																				cout << "CAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 0 && j1 == 2 && nsum == 0)
																			{
																				cout << "AAAAAAAA" << "\n";
																				cout << "GAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 0 && j1 == 3 && nsum == 0)
																			{
																				cout << "AAAAAAAA" << "\n";
																				cout << "TAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 1 && j1 == 2 && nsum == 0)
																			{
																				cout << "CAAAAAAA" << "\n";
																				cout << "GAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i2 == 2 && j1 == 1 && nsum == 0)
																			{
																				cout << "GAAAAAAA" << "\n";
																				cout << "CAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 1 && j1 == 3 && nsum == 0)
																			{
																				cout << "CAAAAAAA" << "\n";
																				cout << "TAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i2 == 3 && j1 == 1 && nsum == 0)
																			{
																				cout << "TAAAAAAA" << "\n";
																				cout << "CAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i1 == 2 && j1 == 3 && nsum == 0)
																			{
																				cout << "GAAAAAAA" << "\n";
																				cout << "TAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																			else if(i2 == 3 && j1 == 2 && nsum == 0)
																			{
																				cout << "TAAAAAAA" << "\n";
																				cout << "GAAAAAAA" << "\n";
																				cout << -dscore << "," << ncharscore << "\n\n";
																			} 
																		}																						
																	}
	}		
	else
	{			
		for(i1 = i1from; i1 < i1to; i1++) 
			for(i2 = i2from; i2 < i2to; i2++)
				for(i3 = i3from; i3 < i3to; i3++)
				{
								
					if(i1 > 0 || i2 > 0 || i3 > 0) foutalignscores.close();
					
					soutalignscores = shmmfile;
					soutalignscores += "_alignscores_8length_";	
					std::ostringstream ot;
					ot << i1;
					soutalignscores += ot.str();
	
					std::ostringstream ot2;
					ot2 << i2;
					soutalignscores += "_";
					soutalignscores += ot2.str();
	
					std::ostringstream ot3;
					ot3 << i3;
					soutalignscores += "_";
					soutalignscores += ot3.str();
					
					foutalignscores.open(soutalignscores.c_str());
					
					if(ncurrentthread == 0)
					{
						d = 25*i1+6.25*i2+1.5625*i3;
						n = round(d);
						cout << "\r" << n << "% processed" << std::flush;						
					}
					
					for(i4 = 0; i4 < nmaxchar; i4++)
						for(i5 = 0; i5 < nmaxchar; i5++)
							for(i6 = 0; i6 < nmaxchar; i6++)
								for(i7 = 0; i7 < nmaxchar; i7++)
									for(i8 = 0; i8 < nmaxchar; i8++)
										for(j1 = 0; j1 < nmaxchar; j1++) 
											for(j2 = 0; j2 < nmaxchar; j2++)
												for(j3 = 0; j3 < nmaxchar; j3++)
													for(j4 = 0; j4 < nmaxchar; j4++)
														for(j5 = 0; j5 < nmaxchar; j5++)
															for(j6 = 0; j6 < nmaxchar; j6++)
																for(j7 = 0; j7 < nmaxchar; j7++)
																	for(j8 = 0; j8 < nmaxchar; j8++)
																	{																		
																		refshort[0] = i1;
																		refshort[1] = i2;
																		refshort[2] = i3;
																		refshort[3] = i4;
																		refshort[4] = i5;
																		refshort[5] = i6;
																		refshort[6] = i7;
																		refshort[7] = i8;
																			
																		readshort[0] = j1;
																		readshort[1] = j2;
																		readshort[2] = j3;
																		readshort[3] = j4;
																		readshort[4] = j5;
																		readshort[5] = j6;
																		readshort[6] = j7;
																		readshort[7] = j8;
																			
																			
																		if(bVariableOrder == 1)
																		{
																			vf = ViterbiVariedOrderHMM(refshort, readshort, "", "", "", nreadlength, nreadlength,
																											   nstates, as, aslogs, statetypes, 
																											   nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
																											   nbases, 0, hmmv, hmmptrs, 0, nindexfactorfull, 
																											   norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, 0);	
																		}
																		else
																		{
																			vf = Viterbi(refshort, readshort, "", "", "", nreadlength, nreadlength,
																										 nstates, as, aslogs, statetypes, 
																										 nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
																										 nbases, 0, hmmv, hmmptrs, 0, nindexfactorfull, 
																										 norder, noutputsho, hmmx, hmmy, 0);	
																		}																																				
					
																		
																		dscore = -vf.vmax/(float)nreadlength;
																													
																																							
																		nid1 = ConvertIntsToID(i1, i2, i3, i4, i5, i6, i7, i8);
																		nid2 = ConvertIntsToID(j1, j2, j3, j4, j5, j6, j7, j8);
																																										  
																		ncharscore = ConvertDoubleScoreCharScore(dscore);
																			
																		foutalignscores << (char)ncharscore;
																			
																		nlinemod++;
																		if(nlinemod % 100 == 0) foutalignscores << "\n";																																							
		
				}													}
	}
	
	foutalignscores.close();	
	
	delete[] statetypes;		
	Clear2DimArray(as, nstates);
	Clear2DimArray(aslogs, nstates);
	Clear2DimArray(p1s, nbases);	
	delete[] qxs;
	delete[] qys;
	
	
	delete[] readshort;
	delete[] refshort;
	
	delete[] hmmprobsho;
	delete[] hmmprobsho1;
	delete[] hmmprobsho2;
	delete[] hmmprobsho3;		
		
	delete[] hmmx;
	delete[] hmmy;
	delete[] hmmv;	

	delete[] hmmptrs;													
}




void CreateScoreMatrix(string sfilehmm, int nThreads)
{
	int i;
	int j;
	int k;
	int l;
	int binit;
	int nThreadCount;
	
	if(nThreads <= 0)
	{
		cout << "\nCreating alignment matrix...\n\n";
		CalcAlignScores(sfilehmm, 0, -1, -1, -1, 0);
		
		cout << "\n\nFinished creating alignment matrix\n\n";		 
		CalcAlignScores(sfilehmm, 1, -1, -1, -1, 0);
		cout << "\n";
	}
	else
	{
		std::thread* Threads = new std::thread[nThreads]; 
		
		nThreadCount = 0;
		binit = 0;
		
		cout << "\nCreating alignment matrix (threading used)...\n\n";
		
		for(i = 0; i < 4; i++)
			for(j = 0; j < 4; j++)
				for(k = 0; k < 4; k++) 
				{
					if(nThreadCount < nThreads)
					{
						Threads[nThreadCount]=std::thread(CalcAlignScores,sfilehmm, 0, i, j, k, nThreadCount);
						nThreadCount++;
					}
					else
					{						
						for(l = 0; l < nThreadCount; l++) Threads[l].join();						
						nThreadCount = 0;					
					}
					
				}
				
			
		for(l = 0; l < nThreadCount; l++) Threads[l].join();		
		delete[] Threads;
		
		cout << "\r" << "100% processed" << std::flush;
		
		cout << "\n\nFinished creating alignment matrix\n\n";		 
		CalcAlignScores(sfilehmm, 1, -1, -1, -1, 0);
		cout << "\n";
	
		
	}
	
}


string GenerateOutputString(int nreadlength, int nalignedposition, int bfastq, double dscore,
							string sreadid, string salignedtogenome, string sstrand, vitstruct vf)
{
	int i;
	int nreadidtextlength;

	string s;
	//string stemp;
	string stemp2;
	string sout;
	string srefalignout;
	string sreadalignout;
	string squalityscoresout;
	
	sout = "";
	
	srefalignout = vf.vouts1;
	sreadalignout = vf.vouts2;
	squalityscoresout = vf.voutinfo;
			
	std::ostringstream outpos;
	std::ostringstream outscore;
	std::ostringstream outreadlength;
	std::ostringstream outreflength;
	
	outscore.precision(10);
	
	outpos                << nalignedposition;									
	outscore              << dscore;	
	outreadlength         << nreadlength;	
	outreflength          << vf.nreflength;				
	
	s = salignedtogenome;
	s += " ";
	s += outpos.str();
																						
	sout += "a score = ";
	sout += outscore.str();
	sout += "\n";			
	
	sout += "s ";								
	sout += s;
	sout += " ";
	sout += outreflength.str();					
	sout += " + ";
	sout += outreflength.str();
	sout += " ";
	sout += srefalignout;
	sout += "\n";		
					
	nreadidtextlength = sreadid.length()+3;								
	if(s.length() <= nreadidtextlength) nreadidtextlength = 1;
	
	string stemp(s.length()+(outreadlength.str()).length()-sreadid.length()-4, ' ');
		
	stemp2 = "s ";
	stemp2 += sreadid;
	stemp2 += stemp;
	stemp2 += "0 ";
	stemp2 += outreadlength.str();
	stemp2 += " ";
	stemp2 += sstrand;
	stemp2 += " ";
	stemp2 += outreadlength.str();
	stemp2 += " ";
					
	sout += stemp2;
	sout += sreadalignout;
	sout += "\n";
	
	string stemp3(stemp2.length()-nreadidtextlength+1, ' ');
			
	if(bfastq == 1)
	{										
		sout += "q ";				
		sout += sreadid;
		sout += stemp3;
		sout += squalityscoresout;
		sout += "\n\n";
	}
	else sout += "\n";
	
	return sout;
}


void HMMAlign(int nfile, int nfrom, int nto, int nreads, int nmaxreadlength, string sref, string srefrev,
					int nCountWindowLen, int nstringidlen, int nstringidlen2, int nMinReadLength,  
					int nmaxrealign, double dminscore, double dminscore2, double dminalignscore, int nMaxCountDiff, int nminscore,					
					int nbases, int nstates, int nbeginstate, int nendstate, int* statetypes, double** as, double** aslogs, double** p1s, double* qxs, double* qys,
					unsigned int* nrefids, unsigned short int* nrefs, unsigned int* nrefrevids, unsigned short int* nrefrevs, char** allidscores,
					char* refcountsA, char* refcountsC, char* refcountsG, char* refcountsT,
					char* refrevcountsA, char* refrevcountsC, char* refrevcountsG, char* refrevcountsT,
					unsigned short int* nreadids, unsigned int* nreadindexes, string* sreads, string* sreadids, 
					string* sshortreads, string* sreadphredscores, int norder, int noutputsho, double* hmmprobsho, 
					string shmmfile, string salignedtogenome, string sreffile, int noutputsho1, int noutputsho2, 
					int noutputsho3, double* hmmprobsho1, double* hmmprobsho2, double* hmmprobsho3, int bVariableOrder, 
					int nnearendbufferlen, int bUsePreCalcMatrix, int bfastq, string sfileout)	
{
	int i;
	int j;
	int k;	
	int l;	
	int ncount;
	int nlen;
	int n;

	int nreadid;
	int nreadidreversestrand;
	int nrefid;
	int nrefrevid;
	int nprevreadid;
	int nreadlength;
	int nreflength;	
	int nSeqLength1;
	int nSeqLength2;	
	int nreflengthalign;
	int nalignedpositionrev;
	int nindexfactorshort;
	int nindexfactorfull;
	int nalignedindex;
	
	int ninitialalignsrev;
				
	int nreadindex;
	int nalignedposition;
	int nrealign;	
	int nstartfromindex;
	int nstartfromindex2;

	int ninitialaligns;
	int ninitialaligns2;
	
	int ndiff;
	int ntotaldiff;
	int ndiffcount;	
	int bprevalign;

	int nreadlengthnew;
	int nreflengthnew;
	
	double d;
	double dmax;
	double dmax2;
	double dAvgDiff;	
	double dscore;				
	
	string s;
	string s1;
	string s2;
	string s3;
	string sout;
	string line;
	string strand;		
	string srefaligned;
	string sfile;
	string sread;
	string sreadinfo;
	string sseq1; 
	string sseq2;
	string sseqinfo;
	string sreadreversestrand;	
	string sreadalign;
	string srefalign;	
	string sreadshort1;
	string sreadshort2;
	string srefshort1;
	string srefshort2;
	string sreadid;
	string sshortread;
	string sprevread;
	string sreadnew;
	string srefalignnew;
	
	charcountstruct cs;
	vitstruct vf;
			
	nreflength = sref.length();
	
	ofstream fout;					
	fout.open(sfileout.c_str(), ofstream::out);		
	
	int* seq1;
	int* seq2; 
	
	unsigned short int* readcountAs = new unsigned short int[nmaxreadlength];
	unsigned short int* readcountCs = new unsigned short int[nmaxreadlength];
	unsigned short int* readcountGs = new unsigned short int[nmaxreadlength];
	unsigned short int* readcountTs = new unsigned short int[nmaxreadlength];
	
	int nidscoreslength = (int)pow(4,nstringidlen);			
	long int hmmnlen = (nmaxreadlength+1)*(nmaxreadlength+1)*nstates;		
	
	int* hmmx = new int[hmmnlen];	
	int* hmmy = new int[hmmnlen];	
		
	double* hmmv = new double[hmmnlen];	
	int* hmmptrs = new int[hmmnlen];	
	
	double* idscores = new double[nidscoreslength];
	double* idscoresrev = new double[nidscoreslength];
	
	int* readshort1 = new int[nstringidlen];
	int* readshort2 = new int[nstringidlen2];
	
	int* refshort1 = new int[nstringidlen];
	int* refshort2 = new int[nstringidlen2];
	
	int* readalign = new int[nmaxreadlength];
	int* refalign = new int[nmaxreadlength];
		
	unsigned int* nlines = new unsigned int[nstringidlen];
	
	int* nalignedpos = new int[nreflength];	
	int* nalignedposrev = new int[nreflength];

	const double dmaxinit = -10000;
	
	nindexfactorshort = (nstringidlen+1)*nstates;
	nprevreadid = -1;

	for(i = nfrom; i < nto; i++)		
	{				
		nreadid = nreadids[i];
		nreadindex = nreadindexes[i];	
		
		sreadid = sreadids[nreadindex];	
		sread = sreads[nreadindex];
		sreadinfo = sreadphredscores[nreadindex];	

		nreadlength = sread.length();
		nindexfactorfull = (nreadlength+1)*nstates;
		
		bprevalign = 0;
		dmax = dmaxinit;	
		dmax2 = dmaxinit;
		nalignedposition = -1000;
		nrealign = 0;

	
		while((dmax < dminalignscore || (dmax2 == dmaxinit && nrealign <= 1)) && nreadlength >= nMinReadLength && nrealign < nmaxrealign)
		{	
			if(nreadid != nprevreadid)
			{				
				//Foward strand 
				ninitialaligns = 0;
				nalignedposition = -1000;
			
				for(j = 0; j < nidscoreslength; j++)
				{
					idscores[j] = -1000;
				}
				
				//Reverse strand 							
				ninitialalignsrev = 0;
			
				for(j = 0; j < nidscoreslength; j++)
				{
					idscoresrev[j] = -1000;
				}
				
				if(bUsePreCalcMatrix == 1)
				{
					//Foward strand 					
					n = 0;
					nlen = nreflength-nreadlength;
					while(n < nlen)
					{	
						nrefid = nrefids[n];				
				
						if(idscores[nrefid] == -1000) idscores[nrefid] = (int)allidscores[nrefid][nreadid];												
					
						if(idscores[nrefid] < nminscore)  
						{																							
							nalignedpos[ninitialaligns] = n;
							ninitialaligns++;						
						}
										
						n++;																
					}
				
					//Reverse strand 												
					n = 0;				
					while(n < nlen)
					{		
						nrefid = nrefrevids[n];				
				
						if(idscoresrev[nrefid] == -1000) idscoresrev[nrefid] = (int)allidscores[nrefid][nreadid];												
					
						if(idscoresrev[nrefid] < nminscore)  
						{																							
							nalignedposrev[ninitialalignsrev] = n;
							ninitialalignsrev++;					
						}
										
						n++;																
					}								
				}
				else //Only exact matches
				{
					
					//Foward strand 								
					n = 0;
					nlen = nreflength-nreadlength;
					while(n < nlen)
					{	
						nrefid = nrefids[n];				
				
						if(idscores[nrefid] == -1000) idscores[nrefid] = (int)allidscores[nrefid][nreadid];												
					
						if(idscores[nrefid] > 0)  
						{																							
							nalignedpos[ninitialaligns] = n;
							ninitialaligns++;						
						}
										
						n++;																
					}
				
					//Reverse strand 												
					n = 0;				
					while(n < nlen)
					{		
						nrefid = nrefrevids[n];				
				
						if(idscoresrev[nrefid] == -1000) idscoresrev[nrefid] = (int)allidscores[nrefid][nreadid];												
					
						if(idscoresrev[nrefid] > 0)  
						{																							
							nalignedposrev[ninitialalignsrev] = n;
							ninitialalignsrev++;					
						}
										
						n++;																
					}
				
					
				}	
				
				nprevreadid = nreadid;				
										
			}					
					
		
			
			//////////////Forward strand
			if(ninitialaligns > 0)
			{
				for(k = nstringidlen; k < nreadlength-nCountWindowLen; k += nCountWindowLen)
				{	
					cs = CountChars(sread.substr(k,nCountWindowLen));
					
					readcountAs[k] = cs.nA;
					readcountCs[k] = cs.nC;
					readcountGs[k] = cs.nG;
					readcountTs[k] = cs.nT;
				}
			}			
			
			//Initial alignments made (using first basepars, e.g. length 8). Now extend alignment 
			//further for those matches above cut of score (e.g. extend to length 24)
			for(j = 0; j < ninitialaligns; j++)
			{
								
				n = nalignedpos[j];
				nrefid = nrefids[n];
								
				ntotaldiff = 0;
				ndiffcount = 0;
								
				for(k = nstringidlen; k < nreadlength-nCountWindowLen; k += nCountWindowLen)
				{							
					ndiff = abs(readcountAs[k]-(int)refcountsA[n+k])+abs(readcountCs[k]-(int)refcountsC[n+k])+abs(readcountGs[k]-(int)refcountsG[n+k])+abs(readcountTs[k]-(int)refcountsT[n+k]);
					ntotaldiff += ndiff;
							
					ndiffcount++;							
				}
										
				if(ndiffcount > 0) dAvgDiff = (double)ntotaldiff/(double)ndiffcount;
				else dAvgDiff = 0;
				
												
				if(dAvgDiff < nMaxCountDiff)
				{
					for(k = 0; k < nstringidlen2; k++) 
					{
						readshort2[k] = ConvertNucToInt(sread[k]);
						refshort2[k] = nrefs[n+k];	
					}
			
					sreadshort2 = sread.substr(0, nstringidlen2);		
					srefshort2 = sref.substr(n, nstringidlen2);													
				
					
					nstartfromindex = 0;


					if(bVariableOrder == 1)
					{
						vf = ViterbiVariedOrderHMM(refshort2, readshort2, srefshort2, sreadshort2, "", nstringidlen2, 
										nstringidlen2, nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
										nbases, 0, hmmv, hmmptrs, nstartfromindex, nindexfactorfull,
										norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, nnearendbufferlen);	
					}
					else
					{
						vf = Viterbi(refshort2, readshort2, srefshort2, sreadshort2, "", nstringidlen2, nstringidlen2,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
										nbases, 0, hmmv, hmmptrs, nstartfromindex, nindexfactorfull, 
										norder, noutputsho, hmmx, hmmy, nnearendbufferlen);	
					}
					
					dscore = vf.vmaxdiv;

				}
				else dscore = dminscore2-1;


				//Full alignments
				if(dscore > dminscore2)  
				{		
					
					for(k = 0; k < nreadlength; k++) 
					{
						readalign[k] = ConvertNucToInt(sread[k]);
						refalign[k] = nrefs[n+k];	
					}
								
					srefalign = sref.substr(n, nreadlength);	
																								
					if(bVariableOrder == 1)
					{
						vf = ViterbiVariedOrderHMM(refalign, readalign, srefalign, sread, sreadinfo, nreadlength, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
										nbases, 1, hmmv, hmmptrs, nstringidlen2, nindexfactorfull, 
										norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, nnearendbufferlen);	
								
					}
					else
					{
						vf = Viterbi(refalign, readalign, srefalign, sread, sreadinfo, nreadlength, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
										nbases, 1, hmmv, hmmptrs, nstringidlen2, nindexfactorfull, 
										norder, noutputsho, hmmx, hmmy, nnearendbufferlen);					
					}

					dscore = vf.vmaxdiv;

					
					//Extend reference string length if alignment ends with insertions in the reference string
					s = vf.vouts1;															
					ncount = 0;
					l = s.length()-1;
					while(l > 0 && s[l] == '-')
					{
						l--;
						ncount++;
					}				
												
					
					if(ncount > 0 && n+ncount+nreadlength < nreflength && ncount+nreadlength < nmaxreadlength)
					{	
									
						for(l = 0; l < ncount; l++)
						{
							refalign[nreadlength+l] = nrefs[n+nreadlength+l];
							srefalign += sref[n+nreadlength+l];												
						}																		
						
						
						if(bVariableOrder == 1)
						{
							vf = ViterbiVariedOrderHMM(refalign, readalign, srefalign, sread, sreadinfo, nreadlength+ncount, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
										nbases, 1, hmmv, hmmptrs, nreadlength-1, nindexfactorfull, 
										norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, nnearendbufferlen);	
						}
						else
						{
							vf = Viterbi(refalign, readalign, srefalign, sread, sreadinfo, nreadlength+ncount, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
										nbases, 1, hmmv, hmmptrs, nreadlength-1, nindexfactorfull, 
										norder, noutputsho, hmmx, hmmy, nnearendbufferlen);					

						}	
						
						
						dscore = vf.vmaxdiv;
						
						
					}
										
					if(dscore > dmax) 
					{
						if(dmax > dmaxinit) dmax2 = dmax;
						dmax = dscore;						
					}
					else if(dscore > dmax2) 
					{
						dmax2 = dmax;
					}
					
					if(dscore >= dminalignscore) 
					{						
						nalignedposition = n;
						strand = "+";
												
						sout = GenerateOutputString(nreadlength, nalignedposition, bfastq, dscore,
													sreadid, salignedtogenome, strand, vf);
											
						fout << sout;
						
					}
					
					
																			
				}				   
			}			


			//////////////Rev strand
			if(ninitialalignsrev > 0)
			{
				for(k = nstringidlen; k < nreadlength-nCountWindowLen; k += nCountWindowLen)
				{	
					cs = CountChars(sread.substr(k,nCountWindowLen));
					
					readcountAs[k] = cs.nA;
					readcountCs[k] = cs.nC;
					readcountGs[k] = cs.nG;
					readcountTs[k] = cs.nT;
				}
			}		
			
			
			//Initial alignments made (using first basepars, e.g. length 8). Now extend alignment 
			//further for those matches above cut of score (e.g. extend to length 24)
			for(j = 0; j < ninitialalignsrev; j++)
			{								
				n = nalignedposrev[j];
								
				ntotaldiff = 0;
				ndiffcount = 0;
								
				for(k = nstringidlen; k < nreadlength-nCountWindowLen; k += nCountWindowLen)
				{							
					ndiff = abs(readcountAs[k]-(int)refrevcountsA[n+k])+abs(readcountCs[k]-(int)refrevcountsC[n+k])+abs(readcountGs[k]-(int)refrevcountsG[n+k])+abs(readcountTs[k]-(int)refrevcountsT[n+k]);
					ntotaldiff += ndiff;
							
					ndiffcount++;							
				}
										
				if(ndiffcount > 0) dAvgDiff = (double)ntotaldiff/(double)ndiffcount;
				else dAvgDiff = 0;
				
				
				if(dAvgDiff < nMaxCountDiff)
				{
					for(k = 0; k < nstringidlen2; k++) 
					{
						readshort2[k] = ConvertNucToInt(sread[k]);
						refshort2[k] = nrefrevs[n+k];	
					}
			
					sreadshort2 = sread.substr(0, nstringidlen2);		
					srefshort2 = srefrev.substr(n, nstringidlen2);													
				
					
					nstartfromindex = 0;


					
					if(bVariableOrder == 1)
					{
						vf = ViterbiVariedOrderHMM(refshort2, readshort2, srefshort2, sreadshort2, "", nstringidlen2, nstringidlen2,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
										nbases, 0, hmmv, hmmptrs, nstartfromindex, nindexfactorfull, norder, 
										noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, nnearendbufferlen);	
					}
					else
					{
						vf = Viterbi(refshort2, readshort2, srefshort2, sreadshort2, "", nstringidlen2, nstringidlen2,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
										nbases, 0, hmmv, hmmptrs, nstartfromindex, nindexfactorfull, 
										norder, noutputsho, hmmx, hmmy, nnearendbufferlen);	
					}
					
					dscore = vf.vmaxdiv;		

				}
				else dscore = dminscore2-1;


				//Full alignments
				if(dscore > dminscore2)  
				{			
							
					for(k = 0; k < nreadlength; k++) 
					{
						readalign[k] = ConvertNucToInt(sread[k]);
						refalign[k] = nrefrevs[n+k];	
					}
								
					srefalign = srefrev.substr(n, nreadlength);	
												
					if(bVariableOrder == 1)
					{
						vf = ViterbiVariedOrderHMM(refalign, readalign, srefalign, sread, sreadinfo, nreadlength, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
										nbases, 1, hmmv, hmmptrs, nstringidlen2, nindexfactorfull,
										norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, nnearendbufferlen);	
								
					}
					else
					{
						vf = Viterbi(refalign, readalign, srefalign, sread, sreadinfo, nreadlength, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
										nbases, 1, hmmv, hmmptrs, nstringidlen2, nindexfactorfull, 
										norder, noutputsho, hmmx, hmmy, nnearendbufferlen);					
					}

					dscore = vf.vmaxdiv;

					
					
					//Extend reference string length if alignment ends with insertions in the reference string
					s = vf.vouts1;															
					ncount = 0;
					l = s.length()-1;
					while(l > 0 && s[l] == '-')
					{
						l--;
						ncount++;
					}														
					
					if(ncount > 0 && n+ncount+nreadlength < nreflength && ncount+nreadlength < nmaxreadlength)
					{	
									
						for(l = 0; l < ncount; l++)
						{
							refalign[nreadlength+l] = nrefrevs[n+nreadlength+l];
							srefalign += srefrev[n+nreadlength+l];												
						}																		
						
						
						if(bVariableOrder == 1)
						{
							vf = ViterbiVariedOrderHMM(refalign, readalign, srefalign, sread, sreadinfo, nreadlength+ncount, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho1, hmmprobsho2, hmmprobsho3,
										nbases, 1, hmmv, hmmptrs, nreadlength-1, nindexfactorfull, 
										norder, noutputsho1, noutputsho2, noutputsho3, hmmx, hmmy, nnearendbufferlen);	
						}
						else
						{
							vf = Viterbi(refalign, readalign, srefalign, sread, sreadinfo, nreadlength+ncount, nreadlength,
										nstates, as, aslogs, statetypes, 
										nbeginstate, nendstate, qxs, qys, p1s, hmmprobsho,
										nbases, 1, hmmv, hmmptrs, nreadlength-1, nindexfactorfull,  
										norder, noutputsho, hmmx, hmmy, nnearendbufferlen);					

						}	
						
						
						dscore = vf.vmaxdiv;
												
						
					}
		
					
					if(dscore > dmax) 
					{
						if(dmax > dmaxinit) dmax2 = dmax;
						dmax = dscore;						
					}
					else if(dscore > dmax2) 
					{
						dmax2 = dmax;
					}
					
					if(dscore >= dminalignscore) 
					{
						nalignedpositionrev = srefrev.length()-n;
						strand = "-";
						
						sout = GenerateOutputString(nreadlength, nalignedpositionrev, bfastq, dscore,
													sreadid, salignedtogenome, strand, vf);																			
				
						fout << sout;
						
					}
					
					
					
																			
				}				   
			}			
													
			nrealign++;
			
			if(dmax < dminalignscore && nreadlength >= (nstringidlen+nMinReadLength) && nrealign < nmaxrealign)
			{
				sread.erase(0,nstringidlen);
				sreadinfo.erase(0,nstringidlen);
				nreadlength = sread.length();
				nindexfactorfull = (nreadlength+1)*nstates;
				
				for(j = 0; j < nstringidlen; j++) nlines[j] = ConvertNucToInt(sread[j]);				
				nreadid = ConvertIntsToID(nlines, 0, nstringidlen);	
			}
			else if((dmax < dminalignscore || (dmax2 == dmaxinit && nrealign <= 1)) && nreadlength >= (nstringidlen+nMinReadLength) && nrealign < nmaxrealign)			
			{			
		
				sread.erase(0,2);
				sreadinfo.erase(0,2);
				
				nreadlength = sread.length();
				nindexfactorfull = (nreadlength+1)*nstates;
				
				for(j = 0; j < nstringidlen; j++) nlines[j] = ConvertNucToInt(sread[j]);				
				nreadid = ConvertIntsToID(nlines, 0, nstringidlen);	
			}
			else 
			{							
				nrealign = nmaxrealign;
			}
			
		}				
		
	}		
	

	
	delete[] nlines;
	
	delete[] readcountAs;
	delete[] readcountCs;
	delete[] readcountGs;
	delete[] readcountTs;
	
	delete[] idscores;
	delete[] idscoresrev;
		
	delete[] readshort1;
	delete[] readshort2;	
	delete[] refshort1;
	delete[] refshort2;

	delete[] readalign;
	delete[] refalign;
	delete[] nalignedpos;
	delete[] nalignedposrev;
	
	delete[] hmmx;
	delete[] hmmy;
	delete[] hmmv;	
	delete[] hmmptrs;
}



void RunHMMAlignments(string sfileref, string sfilereads, string sfilehmm, string sfilescorematrix,
					  string sfileout, int bfastq, int nNumberOfThreads, double dminscore, double dminscore2,
					  double dminalignscore, int nmaxrealign, int nCountWindowLen, int nMaxCountDiff)
{
	int i;
	int j;
	int n;
	int nlen;
	int nfile;
	int nfrom;
	int nto;
	int nreflength;	
	int nstringidlen;
	int nstringidlen2;
	int nMinReadLength;
	int nreads;
	int nlinepos;
	int nmaxreadlength;
	int nidscoreslength;
	int nbufmaxreadlength;
	int bUsePreCalcMatrix;
	
	int nMaxChar;
	int n1;
	int n2;
	int nminscore;	
	int nid1;
	int nid2;
	
	int norder; 
	int noutputsho; 
	int noutputsho1;
	int noutputsho2;
	int noutputsho3;
	int bVariableOrder;
	int nnearendbufferlen;
	
	string s;
	string sref;
	string srefrev;
	string sinalignscores;
	string sfileoutthread; 
	
	string line;
	string line1;
	string line2;
	string line3;
	string line4;
	string salignedtogenome;
	
	charcountstruct cs;
	
	//////////////Initalization of variables
	nstringidlen = 8;
	nstringidlen2 = 24;
	nMinReadLength = 20;
	nbufmaxreadlength = 10;
	nnearendbufferlen = 10;
	nMaxChar = 4; 
	
	nminscore = ConvertDoubleScoreCharScore(-dminscore);	
	
	bUsePreCalcMatrix = 0;
	if(sfilescorematrix.length() > 0) bUsePreCalcMatrix = 1;	
	/////////////////////////////////////////////////////						

	cout << "\nLoading reference and read files, initalizing arrays...\n\n";
	
	//////////Load reference file
	if(sfileout.length() == 0)
	{
		string sfilereadsedit = sfilereads;
		string sfilerefedit = sfileref;
		string sfilehmmedit = sfilehmm;
		
		n1 = sfilereads.find(".",0);		
		if(n1 > 0) sfilereadsedit = sfilereads.substr(0, n1);
				
		n1 = sfileref.find(".",0);		
		if(n1 > 0) sfilerefedit = sfileref.substr(0, n1);						
		
		n1 = sfilehmm.find(".",0);		
		if(n1 > 0) sfilehmmedit = sfilehmm.substr(0, n1);
		
		sfileout = sfilereadsedit;
		sfileout += "_";
		sfileout += sfilerefedit;
		sfileout += "_";
		sfileout += sfilehmmedit;				
	}	
				
	nidscoreslength = (int)pow(4,nstringidlen);
		
	ifstream finref(sfileref.c_str());					
	getline(finref, line);
	
	n1 = line.find(" ", 2);
	salignedtogenome = line.substr(1, n1-1);			
	while(finref.eof() == 0)
	{
		getline(finref, line);
		sref += line;
	}	
	nreflength = sref.length();

	srefrev = ReverseString(sref);
	
	nlen = srefrev.length();
	for(i = 0; i < nlen; i++)
	{
		srefrev[i] = ReverseStrandNuc(srefrev[i]);					
	}				
		
	
	//////////Initialize arrays
		
	unsigned int* nrefids = new unsigned int[nreflength];		
	unsigned short int* nrefs = new unsigned short int[nreflength];
	
	unsigned int* nrefrevids = new unsigned int[nreflength];		
	unsigned short int* nrefrevs = new unsigned short int[nreflength];

	char* refcountsA = new char[nreflength];
	char* refcountsC = new char[nreflength];
	char* refcountsG = new char[nreflength];
	char* refcountsT = new char[nreflength];
	
	char* refrevcountsA = new char[nreflength];
	char* refrevcountsC = new char[nreflength];
	char* refrevcountsG = new char[nreflength];
	char* refrevcountsT = new char[nreflength];

	nlen = nreflength-nstringidlen;
	
	//Arrays for + strand
	for(i = 0; i < nstringidlen; i++)
	{	
		nrefs[i] = ConvertNucToInt(sref[i]);
		
		cs = CountChars(sref.substr(i,nCountWindowLen));		
		refcountsA[i] = (char)cs.nA;
		refcountsC[i] = (char)cs.nC;
		refcountsG[i] = (char)cs.nG;
		refcountsT[i] = (char)cs.nT;
		
	}
	for(i = 0; i < nlen; i++)		 		
	{	
		nrefs[i+nstringidlen] = ConvertNucToInt(sref[i+nstringidlen]);
		nrefids[i] = ConvertIntsToID(nrefs, i, nstringidlen);
		
		cs = CountChars(sref.substr(i,nCountWindowLen));		
		refcountsA[i] = (char)cs.nA;
		refcountsC[i] = (char)cs.nC;
		refcountsG[i] = (char)cs.nG;			
		refcountsT[i] = (char)cs.nT;
	}
	
	for(j = i; j < nreflength; j++)
	{	
		cs = CountChars(sref.substr(j,nCountWindowLen));		
		refcountsA[j] = (char)cs.nA;
		refcountsC[j] = (char)cs.nC;
		refcountsG[j] = (char)cs.nG;			
		refcountsT[j] = (char)cs.nT;
	}

	//Arrays for - strand
	for(i = 0; i < nstringidlen; i++)
	{	
		nrefrevs[i] = ConvertNucToInt(srefrev[i]);
		
		cs = CountChars(srefrev.substr(i,nCountWindowLen));		
		refrevcountsA[i] = (char)cs.nA;
		refrevcountsC[i] = (char)cs.nC;
		refrevcountsG[i] = (char)cs.nG;
		refrevcountsT[i] = (char)cs.nT;
		
	}
	for(i = 0; i < nlen; i++)		 		
	{	
		nrefrevs[i+nstringidlen] = ConvertNucToInt(srefrev[i+nstringidlen]);
		nrefrevids[i] = ConvertIntsToID(nrefrevs, i, nstringidlen);
		
		cs = CountChars(srefrev.substr(i,nCountWindowLen));		
		refrevcountsA[i] = (char)cs.nA;
		refrevcountsC[i] = (char)cs.nC;
		refrevcountsG[i] = (char)cs.nG;			
		refrevcountsT[i] = (char)cs.nT;
	}
	
	for(j = i; j < nreflength; j++)
	{	
		cs = CountChars(srefrev.substr(j,nCountWindowLen));		
		refrevcountsA[j] = (char)cs.nA;
		refrevcountsC[j] = (char)cs.nC;
		refrevcountsG[j] = (char)cs.nG;			
		refrevcountsT[j] = (char)cs.nT;
	}

	ifstream finreads(sfilereads.c_str());


	//////////////Get number of reads
	nreads = 0;
	
	if(bfastq == 1) //fastq
	{
		while(finreads.eof() == 0)
		{
			getline(finreads, line);
		
			nlen = line.length();
				
			if(nlen > 1)
			{
				getline(finreads, line);
				
				if(line.length() >= nMinReadLength)
				{
					nreads++;
				}
				
				getline(finreads, line);
				getline(finreads, line);
			}
		}
	}
	else //fasta
	{
		while(finreads.eof() == 0)
		{
			getline(finreads, line);		
			nlen = line.length();
				
			if(nlen > 1)
			{
				getline(finreads, line);
				
				if(line.length() >= nMinReadLength)
				{
					nreads++;
				}
			}
		}		
	}
	
	finreads.close();
	finreads.open(sfilereads.c_str(), ifstream::in);
	
	cout << "Number of reads = " << nreads << "\n";
	
	unsigned short int* nreadids = new unsigned short int[nreads];	
	unsigned int* nreadindexes = new unsigned int[nreads];
	string* sreads = new string[nreads];
	string* sreadids = new string[nreads];
	
	string* sshortreads = new string[nreads];
	string* sreadphredscores = new string[nreads];
	unsigned int* nlines = new unsigned int[nstringidlen];	

	//////////////Parse reads file		
	nmaxreadlength = 0;

	n = 0;
	
	if(bfastq == 1) //fastq
	{	
		while(finreads.eof() == 0)
		{
			getline(finreads, line1);			
			nlen = line1.length();
				
			if(nlen > 1)
			{
				getline(finreads, line2);
				getline(finreads, line3);
				getline(finreads, line4);
			
			
				if(line2.length() >= nMinReadLength)
				{
					nreadindexes[n] = n;
							
					for(i = 0; i < nstringidlen; i++) nlines[i] = ConvertNucToInt(line2[i]);
					nreadids[n] = ConvertIntsToID(nlines, 0, nstringidlen);
				
					line1.erase(0,1); 
					sreadids[n] = line1;
					sreads[n] = line2;
					sreadphredscores[n] = line4;					

					sshortreads[n] = line2.substr(0, nstringidlen2);
				
					if(line2.length() > nmaxreadlength) nmaxreadlength = line2.length();		
				
					n++;
				}
			}				
		}	
	}
	else //fasta
	{
		while(finreads.eof() == 0)
		{
			getline(finreads, line1);			
			nlen = line1.length();
				
			if(nlen > 1)
			{
				getline(finreads, line2);
			
				if(line2.length() >= nMinReadLength)
				{
					nreadindexes[n] = n;
							
					for(i = 0; i < nstringidlen; i++) nlines[i] = ConvertNucToInt(line2[i]);
					nreadids[n] = ConvertIntsToID(nlines, 0, nstringidlen);
				
					if(line1[0] == '>') line1.erase(0,1);
					sreadids[n] = line1;
					sreads[n] = line2;	
					sreadphredscores[n] = line2; //used as a dummy line for quality scores for fasta			

					sshortreads[n] = line2.substr(0, nstringidlen2);
				
					if(line2.length() > nmaxreadlength) nmaxreadlength = line2.length();		
				
					n++;
				}
			}				
		}	
	}
	
	delete[] nlines;
	
	nmaxreadlength += nbufmaxreadlength;	
	nreads = n;

	cout << "Number of reads used (too short reads removed) = " << nreads << "\n\n";

	//////////////Sort reads
	quickSortStrings(sshortreads, nreadids, nreadindexes, 0, nreads-1);


	//////////////Load pre-calculated score matrix if matrix file specified by user 
	char** allidscores;

	if(bUsePreCalcMatrix == 1)
	{
		int i1, i2, i3, i4, i5, i6, i7, i8;
		int j1, j2, j3, j4, j5, j6, j7, j8;
	
		allidscores = CreateCharArray(nidscoreslength, nidscoreslength);	
		ifstream finalignscores;
	
		cout << "Loading pre-calculated score matrix\n";
		
		n = 0;
	
		for(i1 = 0; i1 < nMaxChar; i1++) 
		{
			for(i2 = 0; i2 < nMaxChar; i2++)
			{								
				for(i3 = 0; i3 < nMaxChar; i3++)
				{
					std::ostringstream os1;
					std::ostringstream os2;
					std::ostringstream os3;
	
					os1 << i1;
					os2 << i2;
					os3 << i3;
		
					sinalignscores = sfilescorematrix;
					sinalignscores += "_alignscores_8length_";
			
					sinalignscores += os1.str();
					sinalignscores += "_";
					sinalignscores += os2.str();
					sinalignscores += "_";
					sinalignscores += os3.str();
			
					if(i1 > 0 || i2 > 0 || i3 > 0) finalignscores.close();
			
					finalignscores.open(sinalignscores.c_str());
			
					getline(finalignscores, line);
			
					nlinepos = 0;
			
					cout << i1 << "," << i2 << "," << i3 << "\n";
			
					for(i4 = 0; i4 < nMaxChar; i4++)
						for(i5 = 0; i5 < nMaxChar; i5++)
							for(i6 = 0; i6 < nMaxChar; i6++)
								for(i7 = 0; i7 < nMaxChar; i7++)
									for(i8 = 0; i8 < nMaxChar; i8++)
										for(j1 = 0; j1 < nMaxChar; j1++) 
											for(j2 = 0; j2 < nMaxChar; j2++)
												for(j3 = 0; j3 < nMaxChar; j3++)
													for(j4 = 0; j4 < nMaxChar; j4++)
														for(j5 = 0; j5 < nMaxChar; j5++)
															for(j6 = 0; j6 < nMaxChar; j6++)
																for(j7 = 0; j7 < nMaxChar; j7++)
																	for(j8 = 0; j8 < nMaxChar; j8++)
																	{		
																		nid1 = ConvertIntsToID(i1, i2, i3, i4, i5, i6, i7, i8);
																		nid2 = ConvertIntsToID(j1, j2, j3, j4, j5, j6, j7, j8);
																	
																		allidscores[nid1][nid2] = line[nlinepos];
																	
																		nlinepos++;
																	
																		if(nlinepos >= line.length()) 
																		{
																			getline(finalignscores, line);
																			nlinepos = 0;
																		}
																	}
				}
			}
		}
		
		finalignscores.close();
		
		cout << "\n";
	}
	else
	{
		int i1, i2, i3, i4, i5, i6, i7, i8;
		int j1, j2, j3, j4, j5, j6, j7, j8;
	
		allidscores = CreateCharArray(nidscoreslength, nidscoreslength);	
		ifstream finalignscores;
	
		cout << "No pre-calculated score matrix specified. Setting up matrix for exact matches only.\n";
	
		n = 0;
	
		for(i1 = 0; i1 < nMaxChar; i1++) 
		{
			for(i2 = 0; i2 < nMaxChar; i2++)
			{								
				for(i3 = 0; i3 < nMaxChar; i3++)
				{
								
					//cout << i1 << "," << i2 << "," << i3 << "\n";
			
					for(i4 = 0; i4 < nMaxChar; i4++)
						for(i5 = 0; i5 < nMaxChar; i5++)
							for(i6 = 0; i6 < nMaxChar; i6++)
								for(i7 = 0; i7 < nMaxChar; i7++)
									for(i8 = 0; i8 < nMaxChar; i8++)
										for(j1 = 0; j1 < nMaxChar; j1++) 
											for(j2 = 0; j2 < nMaxChar; j2++)
												for(j3 = 0; j3 < nMaxChar; j3++)
													for(j4 = 0; j4 < nMaxChar; j4++)
														for(j5 = 0; j5 < nMaxChar; j5++)
															for(j6 = 0; j6 < nMaxChar; j6++)
																for(j7 = 0; j7 < nMaxChar; j7++)
																	for(j8 = 0; j8 < nMaxChar; j8++)
																	{		
																		nid1 = ConvertIntsToID(i1, i2, i3, i4, i5, i6, i7, i8);
																		nid2 = ConvertIntsToID(j1, j2, j3, j4, j5, j6, j7, j8);
																	
																		if(nid1 == nid2) allidscores[nid1][nid2] = 1;
																		else allidscores[nid1][nid2] = 0;																																			
																	}
				}
			}
		}
		
		finalignscores.close();
	}
	
	
	
	
	
	/////////////////Initialize HMM
	ifstream infilehmm;
	infilehmm.open(sfilehmm.c_str(), ifstream::in);	

	int nbases;
	int nstates; 
	int noutputs; 
	int nbeginstate; 
	int nendstate;	
	
	LoadHMMBaseValues(infilehmm, nstates, noutputs, 
					   nbeginstate, nendstate, nbases,
					   norder, bVariableOrder);			

	/////////////////Finished HMM initilization

	double* hmmprobsho;
	double* hmmprobsho1;
	double* hmmprobsho2;
	double* hmmprobsho3;
	
	cout << "HMM order = " << norder << "\n";
	if(bVariableOrder == 0) cout << "Variable order = No\n\n";
	else cout << "Variable order = Yes\n\n";
	
	int* statetypes = new int[nstates];
	double** as = CreateDoubleArray(nstates, nstates);
	double** aslogs = CreateDoubleArray(nstates, nstates);		
	double** p1s = CreateDoubleArray(nbases, nbases);
	double*  qxs = new double[nbases];
	double*  qys = new double[nbases];	
	
	noutputsho = (int)pow(4, norder+1);	

	
	if(bVariableOrder == 1) 
	{
		noutputsho1 = (int)pow(4, 2);
		noutputsho2 = (int)pow(4, 3);
		noutputsho3 = (int)pow(4, 4);
	
		hmmprobsho  = NULL;
		hmmprobsho1 = new double[noutputsho1*noutputsho1*nstates];
		hmmprobsho2 = new double[noutputsho2*noutputsho2*nstates];
		hmmprobsho3 = new double[noutputsho3*noutputsho3*nstates];
	
		InitializeVariedOrderHMM(statetypes, as, aslogs, p1s, qxs, qys, infilehmm, nbases, nstates, noutputsho1, noutputsho2, noutputsho3, hmmprobsho1, hmmprobsho2, hmmprobsho3);
	}
	else 
	{
		hmmprobsho = new double[noutputsho*noutputsho*nstates];
		hmmprobsho1 = NULL;
		hmmprobsho2 = NULL;
		hmmprobsho3 = NULL;
		
		InitializeHMM(statetypes, as, aslogs, p1s, qxs, qys, infilehmm, nbases, nstates, noutputsho, hmmprobsho, norder);
	}
	infilehmm.close();
	
	
	if(nNumberOfThreads > 0)
	{
		cout << "Aligning reads (threaded)...\n";

		int nReadsPerThread = (int)ceil((double)nreads/(double)nNumberOfThreads);
		std::vector<std::thread> Threads;
		Threads.reserve(nNumberOfThreads);
					
		for (i = 0; i < nNumberOfThreads; ++i) 
		{
			nfile = i;
			nfrom = i*nReadsPerThread;
			if(i > 0) nfrom++;
			nto = i*nReadsPerThread+nReadsPerThread;
		
			if(nto > nreads) nto = nreads;
			
			
			std::ostringstream os;
			os << i;
					
			sfileoutthread = sfileout;			
			sfileoutthread += "_t";
			sfileoutthread += os.str();
			sfileoutthread += ".maf";
		
			cout << i << "," << nfrom << "," << nto << "\n";

			Threads.push_back(std::thread(HMMAlign, nfile, nfrom, nto, nreads, nmaxreadlength, sref, srefrev,
					nCountWindowLen, nstringidlen, nstringidlen2, nMinReadLength,  
					nmaxrealign, dminscore, dminscore2, dminalignscore, nMaxCountDiff, nminscore,					
					nbases, nstates, nbeginstate, nendstate, statetypes, as, aslogs, p1s, qxs, qys,
					nrefids, nrefs, nrefrevids, nrefrevs, allidscores,
					refcountsA, refcountsC, refcountsG, refcountsT,
					refrevcountsA, refrevcountsC, refrevcountsG, refrevcountsT,
					nreadids, nreadindexes, sreads, sreadids,
					sshortreads, sreadphredscores, norder, noutputsho, hmmprobsho, 
					sfilehmm, salignedtogenome, sfileref, noutputsho1, noutputsho2, noutputsho3, 
					hmmprobsho1, hmmprobsho2, hmmprobsho3, bVariableOrder, nnearendbufferlen, 
					bUsePreCalcMatrix, bfastq, sfileoutthread));				
					
		}
	
		for (auto &i: Threads) 
		{
			i.join();
		}

	}
	else
	{	
		sfileout += ".maf";
		
		cout << "Aligning reads...\n";
		
		nfile = 0;
		nfrom = 0;
		nto = nreads;
	
		HMMAlign(nfile, nfrom, nto, nreads, nmaxreadlength, sref, srefrev,
					nCountWindowLen, nstringidlen, nstringidlen2, nMinReadLength,  
					nmaxrealign, dminscore, dminscore2, dminalignscore, nMaxCountDiff, nminscore,					
					nbases, nstates, nbeginstate, nendstate, statetypes, as, aslogs, p1s, qxs, qys,
					nrefids, nrefs, nrefrevids, nrefrevs, allidscores,
					refcountsA, refcountsC, refcountsG, refcountsT,
					refrevcountsA, refrevcountsC, refrevcountsG, refrevcountsT,
					nreadids, nreadindexes, sreads, sreadids,
					sshortreads, sreadphredscores, norder, noutputsho, hmmprobsho, 
					sfilehmm, salignedtogenome, sfileref, noutputsho1, noutputsho2, noutputsho3, 
					hmmprobsho1, hmmprobsho2, hmmprobsho3, bVariableOrder, nnearendbufferlen, 
					bUsePreCalcMatrix, bfastq, sfileout);	
						
	}
  
					
	if(bUsePreCalcMatrix == 1) Clear2DimArray(allidscores,nidscoreslength);
	
		
	delete[] refcountsA;
	delete[] refcountsC;
	delete[] refcountsG;
	delete[] refcountsT;
	
	delete[] refrevcountsA;
	delete[] refrevcountsC;
	delete[] refrevcountsG;
	delete[] refrevcountsT;
	
	delete[] nrefs;
	delete[] nrefids;		
	delete[] nrefrevs;
	delete[] nrefrevids;		
	
	delete[] nreadindexes;
	delete[] nreadids;
	delete[] sreads;
	delete[] sreadids;
	delete[] sshortreads;
	delete[] sreadphredscores;
		

	delete[] statetypes;		
	Clear2DimArray(as, nstates);
	Clear2DimArray(aslogs, nstates);
	Clear2DimArray(p1s, nbases);	
	delete[] qxs;
	delete[] qys;	
	
	if(bVariableOrder == 1) 
	{
		delete[] hmmprobsho1;
		delete[] hmmprobsho2;
		delete[] hmmprobsho3;
	}
	else
	{
		delete[] hmmprobsho;
	}
	
}


void PrintHelp(int bPrintVersion)
{
	cout << "\n"; 
	cout << "How to run (alignment):\n";
	cout << "---------------------- \n";
	cout << "varhmm [reference file] [reads file] [hmm training file] (optional arguments)\n";
	cout << "\n";	

	cout << "Optional arguments:\n";
	cout << "-------------------\n";
	cout << "-" << sBFASTQCMD            << sBFASTQDESC            << "[" << BFASTQ            << "]" << "\n";
	cout << "-" << sTHREADSCMD           << sTHREADSDESC           << "[" << THREADS           << "]" << "\n";
	cout << "\n";
	cout << "-" << sMAXCOUNTDIFFARGCMD   << sMAXCOUNTDIFFARGDESC   << "[" << MAXCOUNTDIFF      << "]" << "\n";
	cout << "-" << sCOUNTWINDOWLENGTHCMD << sCOUNTWINDOWLENGTHDESC << "[" << COUNTWINDOWLENGTH << "]" << "\n";
	cout << "-" << sMAXREALIGNCMD        << sMAXREALIGNDESC        << "[" << MAXREALIGN       << "]" << "\n";
	cout << "-" << sMINSCORECMD          << sMINSCOREDESC          << "[" << MINSCORE          << "]" << "\n";
	cout << "-" << sMINSCORE2CMD         << sMINSCORE2DESC         << "[" << MINSCORE2         << "]" << "\n";
	cout << "-" << sMINALIGNSCORECMD     << sMINALIGNSCOREDESC     << "[" << MINALIGNSCORE     << "]" << "\n";
	cout << "\n";	
	cout << "-" << sPRECALCMATRIXFILECMD << sPRECALCMATRIXFILEDESC << "\n";
	cout << "-" << sOUTPUTFILECMD        << sOUTPUTFILEDESC        << "\n";
	cout << "\n";
	
	cout << "How to create a score matrix (optional):\n";
	cout << "----------------------------------------\n";
	cout << "varhmm -C [hmm training file]\n";
	cout << "varhmm -C [hmm training file] -T [# of threads] (threaded version)\n\n";
	
	cout << "Output format:\n";
	cout << "--------------\n";
	cout << "Alignment output is made using the multiple alignment format (MAF) as\n"; 
	cout << "described in the UCSC Genome FAQ, and can be converted to other formats\n";
	cout << "such as SAM using for example the maf-convert command in the alignment\n";
	cout << "software Last.\n";  
	cout << "\n";
	
	cout << "Examples:\n";
	cout << "---------\n";
	cout << "varhmm ecoli_dh10b.fasta reads.fastq train_param\n";
	cout << "varhmm ecoli_dh10b.fasta reads.fasta train_param -Q 0\n";
	cout << "varhmm ecoli_dh10b.fasta reads.fastq train_param -T 4 -n 3 -s -3.0\n";
	cout << "varhmm ecoli_dh10b.fasta reads.fastq train_param -M smatrix -s8 -3.5\n";	
	cout << "varhmm -C train_param -T 4\n";

	if(bPrintVersion == 1)
	{
		cout << "Version: " << VERSION << "\n"; 
		cout << "\n";
	}
}

int main(int argc,char *argv[])
{
	int i;
	int n;
	int narg;
	int nmaxcountdiff;
	int nCountWindowLen;
	int nmaxrealign;
	int nthreads;
	int bfastq;
	int bUsePreCalcMatrix;
	int bCreateScoreMatrix;
	int bThreadCreateScoreMatrix;
	int nThreadsCreateScoreMatrix;
	
	double dminscore;
	double dminscore2;
	double dminalignscore;
	
	string s1;
	string s2;
	string sarg;
	string sfileref;
	string sfilereads;
	string sfilehmm;	
	string sfilescorematrix;		
	string sfileout;

	narg = argc;
	
	if(narg < 3)
	{
		PrintHelp(1);
		return 0;
	}
	
	bCreateScoreMatrix = 0;
	bThreadCreateScoreMatrix = 0;
	nThreadsCreateScoreMatrix = 0;
	
	nmaxcountdiff = MAXCOUNTDIFF;
	nCountWindowLen = COUNTWINDOWLENGTH;	
	nmaxrealign = MAXREALIGN; 		
	dminscore = MINSCORE; 	
	dminscore2 = MINSCORE2;
	dminalignscore = MINALIGNSCORE;
	bfastq = BFASTQ;
	nthreads = THREADS;
	sfilescorematrix = "";
	
	//////////////Set variables as given by user
	i = 1;	
	while(i < narg)
	{
		//cout << i << "\n";
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
			if(s1[0] == '-' && bCreateScoreMatrix == 0)
			{
				PrintHelp(0);
				cout << "Error: an optional argument appears to have been placed before the first 3 required file arguments. Please see above syntax.\n\n";
				return 0;
			}	
			else if(i == 2 && s1[0] == '-')
			{
				PrintHelp(0);
				cout << "Error: an optional argument appears to have been placed before the required file arguments. Please see above syntax.\n\n";
				return 0;
			}			
		} 
		
		switch(i)
		{
			case 1: 
				if(s1.length() > 0 && s1[0] == '-' && s1[1] == 'C') bCreateScoreMatrix = 1;
				else sfileref = s1;
			break;
			case 2: 
				if(bCreateScoreMatrix == 1) sfilehmm = s1;
				else sfilereads = s1;					
			break;
			case 3:
				if(bCreateScoreMatrix == 1)
				{
					if(s1.length() > 1 && s1[0] == '-' && s1[1] == 'T') bThreadCreateScoreMatrix = 1;
					else 
					{
						PrintHelp(0);
						cout << "Error: only argument accepted after -C [hmm training file] is -T. See help for details\n";
					}
				}
				else sfilehmm = s1;
			break;
			default:	
				if(bCreateScoreMatrix == 1)
				{
					if(s1.length() > 0 && bThreadCreateScoreMatrix == 1 && i == 4 && s1[0] != '-')
					{
						nThreadsCreateScoreMatrix = atoi(s1.c_str());
					}
					else
					{
						PrintHelp(0);
						cout << "Error: error in argument(s) passed after -C [hmm training file]. See help for details\n";
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
										
					if(sarg.compare(sMAXCOUNTDIFFARGCMD) == 0)
					{
						nmaxcountdiff = atoi(s2.c_str());				   
				    
						if(nmaxcountdiff < 0) 
					    {
							cout << "\nError: Invalid value for " << "-" << sMAXCOUNTDIFFARGCMD << "n";
							return 0;
						}
					}	
					else if(sarg.compare(sCOUNTWINDOWLENGTHCMD) == 0)
					{
						nCountWindowLen = atoi(s2.c_str());
				   
						if(nCountWindowLen <= 0) 
						{
							cout << "\nError: Invalid value for " << "-" << sCOUNTWINDOWLENGTHCMD << "n";
							return 0;
						}
					}	
					else if(sarg.compare(sMAXREALIGNCMD) == 0)
					{
						nmaxrealign = atoi(s2.c_str());
				   
						if(nmaxrealign  < 0) 
						{
							cout << "\nError: Invalid value for " << "-" << sMAXREALIGNCMD << "n";
							return 0;
						}
					}	
					else if(sarg.compare(sMINSCORECMD) == 0)
					{
						dminscore = atoi(s2.c_str());
					}	
					else if(sarg.compare(sMINSCORE2CMD) == 0)
					{
						dminscore2 = atof(s2.c_str());
					}	
					else if(sarg.compare(sMINALIGNSCORECMD) == 0)
					{
						dminalignscore = atof(s2.c_str());
					}					
					else if(sarg.compare(sTHREADSCMD) == 0)
					{
						nthreads = atoi(s2.c_str());
					}					
					else if(sarg.compare(sBFASTQCMD) == 0)
					{
						bfastq = atoi(s2.c_str());
					}	
					else if(sarg.compare(sPRECALCMATRIXFILECMD) == 0)
					{
						sfilescorematrix = s2;
					}	
					else if(sarg.compare(sOUTPUTFILECMD) == 0)
					{
						sfileout = s2;
					}	
					else
					{
						PrintHelp(0);
						cout << "Error: optional argument not found, please see above syntax.\n\n";
						return 0;
					}															
				}	
			break;
		}
		
		i++;
		
	}
	
	n = sfilereads.find(".", 0);
	
	if(n > 0) 
	{
		n++;
		s1 = sfilereads.substr(n, sfilereads.length()-n);
		for(i = 0; i < s1.length(); i++) s1[i] = toupper(s1[i]);
		
		
		if(s1.compare("FASTQ") == 0 && bfastq == 0)
		{	
			cout << "\n";
			cout << "The read file ends in .fastq and so appears to be in fastq\nformat but the read file input was set to fasta. Please\nchange the read file or change the input format using:\n";
			cout << "\n";
			cout << "-" << sBFASTQCMD            << sBFASTQDESC            << "[" << BFASTQ            << "]" << "\n";
			cout << "\n";
			return 0;
		}
		else if(s1.compare("FASTA") == 0 && bfastq == 1)
		{	
			cout << "\n";
			cout << "The read file ends in .fasta and so appears to be in fasta\nformat but the read file input is set to fastq. Please\nchange the read file or change the input format using:\n";
			cout << "\n";
			cout << "-" << sBFASTQCMD            << sBFASTQDESC            << "[" << BFASTQ            << "]" << "\n";
			cout << "\n";
			return 0;
		}
		else if(s1.compare("FA") == 0 && bfastq == 1)
		{	
			cout << "The read file ends in .fa and so appears to be in fasta\nformat but the read file input is set to fastq. Please\nchange the read file or change the input format using:\n";
			cout << "\n";
			cout << "-" << sBFASTQCMD            << sBFASTQDESC            << "[" << BFASTQ            << "]" << "\n";
			cout << "\n";
			return 0;
		}
	}
	
	if(bCreateScoreMatrix == 0)
	{
		RunHMMAlignments(sfileref, sfilereads, sfilehmm, sfilescorematrix,
					  sfileout, bfastq, nthreads, dminscore, dminscore2,
					  dminalignscore, nmaxrealign, nCountWindowLen, nmaxcountdiff);		
	}
	else
	{
		CreateScoreMatrix(sfilehmm, nThreadsCreateScoreMatrix);
	}
	
	return 0;
	
}
