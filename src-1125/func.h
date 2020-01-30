#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mrc.h"
#include "mrcz.h"
#include "mrct.h"


#include <complex>
#include <vector>
#include <string>
using namespace std;

struct CPLX
{
	float x;
	float y;
};

CPLX cXc(CPLX a, CPLX b);

void fft2d(float* buf, int nsam);
void ifft2d(float* buf, int nsam);


void SetFastFFT(float *buf, int nsam);
void ReleaseFastFFT();
void fft2d_fast();
void ifft2d_fast();


void buf2fft(float *buf, float *fft, int nsam);
void fft2buf(float *buf, float *fft, int nsam);
void buf2mrc(const char *filename, float* buf, int nx, int ny, int nz);


int verifyCropSize(int nsamin, int offsetx, int offsety, int nsamout, float FutureBin);
int crop(float *bufin, int nsamin, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin=1.0);  //return new nsam, may different from nsamout
int crop2fft(float *bufin, int nsamin, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin=1.0);

int verifyCropSize(int nx, int ny, int offsetx, int offsety, int nsamout, float FutureBin);
int crop(float *bufin, int nx, int ny, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin=1.0);  //return new nsam, may different from nsamout
int crop2fft(float *bufin, int nx, int ny, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin=1.0);

void MinMaxMean(float *buf, size_t size, float &min, float &max, float &mean);
float STD(float *buf, size_t size, float mean);

void shift2d_phase(float *buf, int nsam, float sx, float sy);
void CPLXAdd(CPLX *dst, CPLX *src, int sizec);

void FFTDataToLogDisp(float *pIn, float *pOut, int dim);

void rmNoiseCC(float* hboxmap, int box, int nframe, int peakR);
void cosmask2d(complex<float> *pfft, int nsam);

void histMinMax(float *buf, int size, float &min, float &max, float threshold=0.005);
void buf2Disp(char *disp, float *buf, int size);
void buf2DispShort(short *disp, float *buf, int size);

void FFTModulusToDispBuf(float *pIn, float *pOut, int nsam);
void BinFFTDispBufToChar(char *pDisp, int dispdim, float *pfft, int nsam);
void BinFFTDispBufToChar(short *pDisp, int dispdim, float *pfft, int nsam);

bool isMRC(const char * m_fn, int verb=0);
bool isMRCZ(const char * m_fn, int verb=0);
bool isTIFF(const char * m_fn, int verb=0);
