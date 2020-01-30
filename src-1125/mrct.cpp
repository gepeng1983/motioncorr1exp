#include "mrct.h"
#include "stdio.h"
#include <cstring>

MRCT::MRCT()
{
	memset((void *)&m_header,0,sizeof(MRCHeader));
	m_header.mode=2;
	m_header.cellb[0]=90;
	m_header.cellb[1]=90;
	m_header.cellb[2]=90;
	m_header.mapc=1;
	m_header.mapr=2;
	m_header.maps=3;
	m_fp=NULL;
}

MRCT::MRCT(const char * fn, const char * mode)
{
	MRCT();
	open(fn,mode);
}

MRCT::~MRCT()
{
}

int MRCT::open(const char *filename, const char *mode)
{
	close();
	m_fpt=TIFFOpen(filename, mode);
	if (m_fpt==NULL) return 0;

	fakeReadHeader();
	return 1;
	
}

void MRCT::close()
{
	if (m_fpt!=NULL) TIFFClose(m_fpt);
	m_fpt=NULL;
}

int MRCT::write2DIm(void *buf, int n)
{
	printf("Writing to TIFF not supported\n");
}
	
int MRCT::writeLine(void *buf, int nimage, int nline)
{
	printf("Writing to TIFF not supported\n");
}
	
int MRCT::writePixel(void *buf, int nimage, int nline, int npixel)
{
	printf("Writing to TIFF not supported\n");
}

void MRCT::setHeader(const MRCHeader *pheader)
{
	printf("Writing to TIFF not supported\n");
}

void MRCT::updateHeader()
{
	printf("Writing to TIFF not supported\n");
}
	
int MRCT::createMRC(float *data, int nx, int ny, int nz)
{
	printf("Writing to TIFF not supported\n");
}

int MRCT::createMRC(short *data, int nx, int ny, int nz)
{
	printf("Writing to TIFF not supported\n");
}

void MRCT::fakeReadHeader()
{
	m_header.mode=7199;

	int tiffdir=0;
	int zcount=0;
	int nx,ny;
	do {
	TIFFGetField(m_fpt, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(m_fpt, TIFFTAG_IMAGELENGTH, &ny);
	zcount++;
	} while (tiffdir=TIFFReadDirectory(m_fpt));

	m_header.nx=nx;
	m_header.ny=ny;
	m_header.nz=zcount;

	TIFFSetDirectory(m_fpt, 0);
}

int MRCT::read2DIm(void *buf, int n)
{
	//TIFFPrintDirectory(m_fpt, stdout);

	uint16 bitspersample=0;
	TIFFGetField(m_fpt, TIFFTAG_BITSPERSAMPLE,&bitspersample);
	
	tsize_t stripSize=TIFFStripSize(m_fpt);
	tsize_t stripMax=TIFFNumberOfStrips(m_fpt);

	tsize_t bufferSize=stripMax*stripSize;
	tsize_t imageOffset=0;

	unsigned char * cdata;

	if((cdata=(unsigned char*)_TIFFmalloc(bufferSize))==NULL){
                fprintf(stderr,"Error: Could not allocate enough memory\n");
                return(-1);
        }

	//reading image as multiple strips
        //(***) CAUTION: the last strips may have empty values
        for (int stripCount=0;stripCount<stripMax;stripCount++)
	{
		tsize_t result;
                if((result=TIFFReadEncodedStrip(m_fpt,stripCount,cdata+imageOffset,stripSize))==-1)
		{
	                fprintf(stderr,"Error reading stripped image\n");return(-1);
		}
                imageOffset+=result;
        }
	
	float * fbuf = (float *)buf;
	for (int j=0; j<getNy(); j++)
		for (int i=0; i<getNx(); i++)
		{
			if (bitspersample==8)
			{
				fbuf[i+j*getNx()]=(float) ((unsigned char *)cdata)[i+j*getNx()];
			}
			else if (bitspersample==16)
			{
				fbuf[i+j*getNx()]=(float) ((unsigned short *)cdata)[i+j*getNx()];
			}
			else
			{
				printf("Only 8-bit or 16-bit TIFF is supported\n");
			}
		}
	
	TIFFReadDirectory(m_fpt);

	return getNx()*getNy()*4;
}

int MRCT::readLine(void *buf, int nimage, int nline)
{
	printf("To be implemented\n");
}

int MRCT::readPixel(void *buf, int nimage, int nline, int npixel)
{
	printf("To be implemented\n");
}

int MRCT::readnPixels(void *buf, int nimage, int nline, int npixel,int n)
{
	printf("To be implemented\n");
}

int MRCT::getWordLength()
{
	return 4;
}
