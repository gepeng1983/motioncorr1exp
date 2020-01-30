#pragma once
#include "stdio.h"
#include "mrcheader.h"
#include "mrc.h"
#include <tiffio.h>


class MRCT : public MRC
{
public:
	MRCT();
	MRCT(const char *filename, const char *mode);
	~MRCT();

private:
	void fakeReadHeader();
	
public:	
	int open(const char *filename, const char *mode);
	void close();
	
	void printInfo();
	
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

protected:
	TIFF * m_fpt;
	

};











