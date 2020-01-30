#pragma once
#include <complex>
#include <vector>
#include "mrc.h"
#include "mrcz.h"
#include "mrct.h"
#include "cufunc.h"
#include <string>
#include <pthread.h> 
using namespace std;


#define DISPDIM 400


struct APARA
{
	int crop_offsetx;  
	int crop_offsety;  
	int crop_nsam;   

	float bin;

	float xscale, yscale; // X and Y scale factor

	int nStart;  //first frame(0-base)
	int nEnd;    //last frame(0-base)
	int nStartSum;  //first frame to sum(0-base)
	int nEndSum;    //last frame to sum(0-base)
	int nEndSum2;    //last frame to sum(0-base) for second sum image


	int GPUNum;  // GPU device ID 

	float bfactor;  // in pix^2
	int CCPeakSearchDim;//search peak in this box
	int FrameDistOffset;
	int NoisePeakSize;
	float kiThresh;

	bool bSaveSubAreaCorrSum;
	bool bSaveRawSum;
	bool bSaveStackCorr;
	bool bSaveCCmap;
	bool bSaveLog;

	bool bAlignToMid;

	

	//diplay para
	float fftscale;
	bool bDispSumCorr;
	bool bDispFFTRaw;
	bool bDispFFTCorr;
	bool bDispCCMap;
	bool bDispFSC;
	bool bLogFSC;


	//reserved parameters for Dialog window
	float fscMax;
};

class CDFAlign
{
public:
	CDFAlign(void);
	~CDFAlign(void);

	bool m_bRun;
	//HWND m_dlgwnd;
	
	//Iterative alignment
	int m_iter;
	float m_iterthres;

	//output filename
	char m_fnStack[512];
	char m_fnRawsum[512];
	char m_fnAlignsum[512];
	char m_fnAlignsum2[512];
	char m_fnStackCorr[512];
	char m_fnCCmap[512];
	char m_fnLog[512];
	// dark and gain filename
	char m_fnDark[512];
	char m_fnNorm[512];
	
	//quick display filename
	bool m_bSaveDisp;
	char m_dispRawFFT[512];
	char m_dispCorrSum[512];
	char m_dispCorrFFT[512];
	
	
	APARA m_para;

	string m_log;
	int m_nsam;
	int m_nsamRaw;
	
	vector<pthread_t> m_tids;

public:
	static void Message(const char *str);
	
	void UpdateDisplay();
	void TextOutput(const char *str);

	void ImageOutput(float *buf);
	float *m_bufIm;
	short *m_dispIm;
	static void* ImageOutputThread(void *p);

	void FFTOutputCorr(float *buf);
	float *m_bufFFTCorr;
	short *m_dispFFTCorr;
	static void* FFTOutputCorrThread(void *p);

	void FFTOutputRaw(float *buf);
	float *m_bufFFTRaw;
	short *m_dispFFTRaw;
	static void* FFTOutputRawThread(void *p);

	void CCMapOutput(float *buf, void *pki);
	float *m_bufCCMap;
	vector<double> m_kiCCMap;
	static void* CCMapOutputThread(void *p);

	void PlotOutput(vector<complex<double> > &xy);
	vector<complex<double> > m_curve;

	void PlotFSC(float2* hRaw0, float2 *hRaw1, float2 *hCorr0, float2 *hCorr1,
					MASK *pPosList, int nsam, complex<double> direction);
	vector<complex<double> > m_fscRaw0;
	vector<complex<double> > m_fscRaw1;
	vector<complex<double> > m_fscCorr0;
	vector<complex<double> > m_fscCorr1;


	void Done();

public:
	void RunAlign();
	static void* ThreadFunc_cuAlign(void* p);

	int getNFrame();
	int getNFrame(const char *filename);
	MRCHeader getMrcHeader(const char *filename);
};

