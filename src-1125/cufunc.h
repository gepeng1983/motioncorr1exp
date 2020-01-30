#pragma once

#include <cuda.h>
#include <cufft.h>
#include <string>
#include <vector>
using namespace std;

struct MASK
{
	int x;
	int y;
	float z;
};

bool initGPU(int GPUNum);
bool ResetGPU();
void siginthandler(int param);
int getGPUList(vector<string> &namelist);
void GPUMemCheck(size_t &theFree, size_t &theTotal);

bool GPUMemAlloc(void **buf, int size);
bool GPUMemZero(void **buf, int size);
bool GPUMemFree(void **buf);
bool GPUMemH2D(void *dst, void *src, int size);
bool GPUMemD2H(void *dst, void *src, int size);
bool GPUMemD2D(void *dst, void *src, int size);
bool GPUMemBinD2H(float *dst, float *src, int dst_nsam, int src_nsam);
bool GPUMemBinD2D(float *dst, float *src, int dst_nsam, int src_nsam);
cufftHandle GPUFFTPlan(int nsam);
cufftHandle GPUIFFTPlan(int nsam);
bool GPUFFT2d(float* dfft, cufftHandle plan);
bool GPUIFFT2d(float* dfft, cufftHandle plan);
void GPUFFTDestroy(cufftHandle &plan);
bool GPUSync();
void GPUAdd(float *dst, float *src, int size);
void GPUMultiplyNum(float *dst, float num, int size);

void GPUFFTLogModulus(float *dMod, float *dfft, int nsam, float scale);
void GPUFFTModulus(float *dMod, float *dfft, int nsam);

void MkPosList(int3 *list, int nsam, float inner_r, float outer_r);
void MkPosList(MASK *list, int nsam, float bfactor);
void GPUShiftCC(float *dfft, float *dsum, MASK *dposlist, float sx, float sy, int nsam);
void GPUShift(float *dfft, MASK *dposlist, float sx, float sy, int nsam);
float FindShift(float *dsrc,int nsam, float* hboxmap, int box, float &sx, float &sy, int wNoise=1); //wNoise==-1 to disable it
float FindShift(float* hboxmap, int box, float &sx, float &sy);

void testCUFFT();
