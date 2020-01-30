#pragma once
#include "stdio.h"
#include "zlib.h"
#include "mrcheader.h"
#include "mrc.h"


class MRCZ : public MRC
{
public:
	MRCZ();
	MRCZ(const char *filename, const char *mode);
	~MRCZ();
	
public:	
	
	void printInfo();

	int open(const char *filename, const char *mode);
	void close();
	
	int read2DIm(void *buf, int n);
	int readLine(void *buf, int nimage, int nline);
	int readPixel(void *buf, int nimage, int nline, int npixel);
        int readnPixels(void *buf, int nimage, int nline, int npixel,int n);
	int write2DIm(void *buf, int n);
	int writeLine(void *buf, int nimage, int nline);
	int writePixel(void *buf, int nimage, int nline, int npixel);
	void setHeader(const MRCHeader *pheader);
	void updateHeader();
	
	int createMRC(float *data, int nx, int ny, int nz);
	int createMRC(short *data, int nx, int ny, int nz);

	int getWordLength();   //get word length in byte


private:
	
	gzFile m_fpz;

};











