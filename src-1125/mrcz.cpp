#include "mrcz.h"
#include "stdio.h"
#include <cstring>
#include <zlib.h>

MRCZ::MRCZ()
{
	memset((void *)&m_header,0,sizeof(MRCHeader));
	m_header.mode=2;
	m_header.cellb[0]=90;
	m_header.cellb[1]=90;
	m_header.cellb[2]=90;
	m_header.mapc=1;
	m_header.mapr=2;
	m_header.maps=3;
	m_fpz=NULL;
}

MRCZ::MRCZ(const char *filename, const char *mode)
{
	MRCZ();
	int ret=open(filename,mode);
	if (ret > 0) return;
	else printf("Error opening file %s: %d\n", filename, ret);
}

MRCZ::~MRCZ()
{

}

int MRCZ::open(const char *filename, const char *mode)
{
	close();
	m_fpz=gzopen(filename,mode);
	if(m_fpz==NULL) return 0;
	
	//read file header
	gzrewind(m_fpz);
	if(gzread(m_fpz,&m_header,1024)<1024) return -1;
	
	return 1;
}

void MRCZ::close()
{
	if(m_fpz!=NULL) gzclose(m_fpz);
	m_fpz=NULL;
}

void MRCZ::printInfo()
{
	if(m_fpz==NULL) 
	{
		printf("No MRCZ file was opened!");
		return;
	}
	
	printf("\tMRC Header size:                   %12d\n",sizeof(MRCHeader));
	printf("\tNum of columns, rows,sections:     %12d %12d %12d\n",m_header.nx,m_header.ny,m_header.nz);
	printf("\tMode:                              %12d\n",m_header.mode);
	printf("\tNum of First column, row, section: %12d %12d %12d\n",m_header.nxstart, m_header.nystart, m_header.nzstart);
	printf("\tNum of intervals along x, y, z:    %12d %12d %12d\n",m_header.mx,m_header.my,m_header.mz);
	printf("\tCell dimensions in angstroms:      %12.3f %12.3f %12.3f\n",m_header.cella[0], m_header.cella[1],m_header.cella[2]);
	printf("\tCell angles in degrees:            %12.3f %12.3f %12.3f\n",m_header.cellb[0], m_header.cellb[1],m_header.cellb[2]);
	printf("\tAxis for cols, rows, sections:     %12d %12d %12d\n",m_header.mapc, m_header.mapr, m_header.maps);
	printf("\tMin, max, mean density value:      %12.3f %12.3f %12.3f\n",m_header.dmin, m_header.dmax, m_header.dmean);
	printf("\tSpace group number:                %12d\n",m_header.ispg);
	printf("\tNum of bytes for symmetry data:    %12d\n",m_header.nsymbt);
	//printf("\tExtra:                             %s\n",m_header.extra);
	printf("\tOrigin in X,Y,Z:                   %12.3f %12.3f %12.3f\n",m_header.origin[0], m_header.origin[1], m_header.origin[2]);
	printf("\tFile type:                         %c%c%c%c\n",m_header.map[0],m_header.map[1],m_header.map[2],m_header.map[3]);
	printf("\tMachine stamp:                     %12d\n",m_header.machst);
	printf("\trms deviationfrom mean density:    %12.3f\n",m_header.rms);
	printf("\tNum of labels being used:          %12d\n",m_header.nlabels);
	for(int i=0;i<m_header.nlabels;i++)
	{
		printf("\t\t%s\n",m_header.label[i]);
	}
}

int MRCZ::read2DIm(void *buf, int n)
{
        size_t ImSize=getImSize();
        size_t offset=1024+getSymdatasize()+(size_t)n*ImSize;
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzread(m_fpz, buf, ImSize);
}

int MRCZ::readLine(void *buf, int nimage, int nline)
{
        size_t ImSize=getImSize();
        size_t LineLength=getNy()*getWordLength();
        size_t offset=1024+getSymdatasize()+(size_t)nimage*ImSize+nline*LineLength;
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzread(m_fpz, buf, LineLength);
}

int MRCZ::readPixel(void *buf, int nimage, int nline, int npixel)
{
        size_t ImSize=getImSize();
        size_t LineLength=getNy()*getWordLength();
	if(npixel>=LineLength) return 0;
	size_t offset=1024+getSymdatasize()+nimage*ImSize+nline*LineLength+npixel*getWordLength();
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzread(m_fpz, buf, getWordLength());
}

int MRCZ::readnPixels(void *buf, int nimage, int nline, int npixel, int n)
{
        size_t ImSize=getImSize();
        size_t LineLength=getNy()*getWordLength();
        if(npixel>=LineLength) return 0;
        size_t offset=1024+getSymdatasize()+nimage*ImSize+nline*LineLength+npixel*getWordLength();
        if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
        return gzread(m_fpz, buf, n*getWordLength());
}

int MRCZ::write2DIm(void *buf, int n)
{
        size_t ImSize=(size_t)getImSize();
        size_t offset=1024+getSymdatasize()+(size_t)n*ImSize;
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzwrite(m_fpz, buf, ImSize);
}

int MRCZ::writeLine(void *buf, int nimage, int nline)
{
        size_t ImSize=getImSize();
        size_t LineLength=getNy()*getWordLength();
	size_t offset=1024+getSymdatasize()+nimage*ImSize+nline*LineLength;
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzwrite(m_fpz, buf, LineLength);
}

int MRCZ::writePixel(void *buf, int nimage, int nline, int npixel)
{
        size_t ImSize=getImSize();
        size_t LineLength=getNy()*getWordLength();
	if(npixel>=LineLength) return 0;
	size_t offset=1024+getSymdatasize()+nimage*ImSize+nline*LineLength+npixel*getWordLength();
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzwrite(m_fpz, buf, getWordLength());
}

void MRCZ::setHeader(const MRCHeader *pheader)
{
	memcpy((void *)&m_header, (void *)pheader, 1024);
	gzrewind(m_fpz);
	gzwrite(m_fpz, &m_header,1024);
}

void MRCZ::updateHeader()
{
	gzrewind(m_fpz);
	gzwrite(m_fpz, &m_header,1024);
}

int MRCZ::createMRC(float *data, int nx, int ny, int nz)
{
	m_header.nx=nx;
	m_header.ny=ny;
	m_header.nz=nz;
	m_header.mode=2;
	m_header.cellb[0]=90;
	m_header.cellb[1]=90;
	m_header.cellb[2]=90;
	m_header.mapc=1;
	m_header.mapr=2;
	m_header.maps=3;
	size_t i,size=nx*ny*nz;
	float min=data[0],max=data[0];
	double mean=0;
	for(i=0;i<size;i++)
	{
		if(min>data[i]) min=data[i];
		if(max<data[i]) max=data[i];
		mean+=data[i];
	}
	mean/=size;
	m_header.dmin=min;
	m_header.dmax=max;
	m_header.dmean=mean;
	updateHeader();
	
	size_t ImSize=size*sizeof(float);
	size_t offset=1024+getSymdatasize();
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzwrite(m_fpz, data, ImSize);
}

int MRCZ::createMRC(short *data, int nx, int ny, int nz)
{
	m_header.nx=nx;
	m_header.ny=ny;
	m_header.nz=nz;
	m_header.mode=1;
	m_header.cellb[0]=90;
	m_header.cellb[1]=90;
	m_header.cellb[2]=90;
	m_header.mapc=1;
	m_header.mapr=2;
	m_header.maps=3;
	size_t i,size=nx*ny*nz;
	short min=data[0],max=data[0];
	double mean=0;
	for(i=0;i<size;i++)
	{
		if(min>data[i]) min=data[i];
		if(max<data[i]) max=data[i];
		mean+=data[i];
	}
	mean/=size;
	m_header.dmin=min;
	m_header.dmax=max;
	m_header.dmean=mean;
	updateHeader();
	
	size_t ImSize=size*sizeof(short);
	size_t offset=1024+getSymdatasize();
	if(gzseek(m_fpz, offset, SEEK_SET)!=offset) return 0;
	return gzwrite(m_fpz, data, ImSize);
}

int MRCZ::getWordLength()
{
	return MRC::getWordLength();
}
