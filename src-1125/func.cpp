#include "func.h"
#include "fftw3.h"
#include <pthread.h>

#pragma comment(lib,"libfftw3f-3.lib")
fftwf_plan plan_fft_fast;
fftwf_plan plan_ifft_fast;


void fft2d(float* buf, int nsam)
{
	fftwf_plan plan_fft=fftwf_plan_dft_r2c_2d(nsam,nsam,buf,reinterpret_cast<fftwf_complex *>(buf),FFTW_ESTIMATE);  
	fftwf_execute(plan_fft);
	fftwf_destroy_plan(plan_fft);
}


void ifft2d(float* buf, int nsam)
{
	fftwf_plan plan_fft=fftwf_plan_dft_c2r_2d(nsam,nsam,reinterpret_cast<fftwf_complex *>(buf),buf,FFTW_ESTIMATE); 
	fftwf_execute(plan_fft);	
	fftwf_destroy_plan(plan_fft);
}



void SetFastFFT(float *buf, int nsam)
{
	plan_fft_fast=fftwf_plan_dft_r2c_2d(nsam,nsam,buf,reinterpret_cast<fftwf_complex *>(buf),FFTW_ESTIMATE); 
	plan_ifft_fast=fftwf_plan_dft_c2r_2d(nsam,nsam,reinterpret_cast<fftwf_complex *>(buf),buf,FFTW_ESTIMATE); 
}
void ReleaseFastFFT()
{
	fftwf_destroy_plan(plan_fft_fast);
	fftwf_destroy_plan(plan_ifft_fast);
}
void fft2d_fast()
{
	fftwf_execute(plan_fft_fast);
}
void ifft2d_fast()
{
	fftwf_execute(plan_ifft_fast);
}

CPLX cXc(CPLX a, CPLX b)
{
	CPLX c;
	c.x=a.x*b.x-a.y*b.y;
	c.y=a.x*b.y+a.y*b.x;
	return c;
}


int verifyCropSize(int nsamin, int offsetx, int offsety, int nsamout, float FutureBin)
{
	if(offsetx<0 || offsetx>=nsamin || offsety<0 || offsety>=nsamin || FutureBin<=0) return 0;

	if(nsamout<=0) nsamout=min(nsamin-offsetx,nsamin-offsety);

	if((nsamin-offsetx)<nsamout) nsamout=nsamin-offsetx;
	if((nsamin-offsety)<nsamout) nsamout=nsamin-offsety;
	if (fabs(FutureBin-1.0)<0.0001)
	{
		nsamout=nsamout/2*2;
	}
	else if (fabs(FutureBin-2.0)<0.0001)
	{
		nsamout=nsamout/4*4;
	}
	else
	{
		nsamout=nsamout/2*2;
	}
//	nsamout=nsamout/(2*FutureBin)*(2*FutureBin);

	return nsamout;
}

int verifyCropSize(int nx, int ny, int offsetx, int offsety, int nsamout, float FutureBin)
{
	if(offsetx<0 || offsetx>=nx || offsety<0 || offsety>=ny || FutureBin<=0) return 0;

	if(nsamout<=0) nsamout=min(nx-offsetx,ny-offsety);

	if((nx-offsetx)<nsamout) nsamout=nx-offsetx;
	if((ny-offsety)<nsamout) nsamout=ny-offsety;
	if (fabs(FutureBin-1.0)<0.0001)
	{
		nsamout=nsamout/2*2;
	}
	else if (fabs(FutureBin-2.0)<0.0001)
	{
		nsamout=nsamout/4*4;
	}
	else
	{
		nsamout=nsamout/2*2;
	}
//	nsamout=nsamout/(2*FutureBin)*(2*FutureBin);

	return nsamout;
}

int crop(float *bufin, int nsamin, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin)
{
	if(bufin==0 || bufout==0) return 0;

	nsamout=verifyCropSize(nsamin, offsetx, offsety, nsamout,FutureBin);
	if(nsamout==0) return 0;

	int i;
	for(i=0;i<nsamout;i++)
	{
		memcpy(bufout+i*nsamout, bufin+(i+offsety)*nsamin+offsetx, sizeof(float)*nsamout);
	}

	return nsamout;

}

int crop(float *bufin, int nx, int ny, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin)
{
	if(bufin==0 || bufout==0) return 0;

	nsamout=verifyCropSize(nx, ny, offsetx, offsety, nsamout,FutureBin);
	if(nsamout==0) return 0;

	int i;
	for(i=0;i<nsamout;i++)
	{
		memcpy(bufout+i*nsamout, bufin+(i+offsety)*nx+offsetx, sizeof(float)*nsamout);
	}

	return nsamout;

}

int crop2fft(float *bufin, int nsamin, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin)
{
	if(bufin==0 || bufout==0) return 0;

	nsamout=verifyCropSize(nsamin, offsetx, offsety, nsamout, FutureBin);
	//memset(bufout,0,(nsamout+2)*nsamout*sizeof(float));
	if(nsamout==0) return 0;

	int i;
	for(i=0;i<nsamout;i++)
	{
		memcpy(bufout+i*(nsamout+2), bufin+(i+offsety)*nsamin+offsetx, sizeof(float)*nsamout);
	}

	return nsamout;

}

int crop2fft(float *bufin, int nx, int ny, float *bufout, int offsetx, int offsety, int nsamout, float FutureBin)
{
	if(bufin==0 || bufout==0) return 0;

	nsamout=verifyCropSize(nx, ny, offsetx, offsety, nsamout, FutureBin);
	//memset(bufout,0,(nsamout+2)*nsamout*sizeof(float));
	if(nsamout==0) return 0;

	int i;
	for(i=0;i<nsamout;i++)
	{
		memcpy(bufout+i*(nsamout+2), bufin+(i+offsety)*nx+offsetx, sizeof(float)*nsamout);
	}

	return nsamout;

}

void buf2mrc(const char *filename, float* buf, int nx, int ny, int nz)
{
	MRC mrc;
	mrc.open(filename,"wb");
	mrc.createMRC(buf,nx,ny,nz);
	mrc.close();
}

void buf2fft(float *buf, float *fft, int nsam)
{
	int nsamb=nsam+2;
	int i;
	for(i=0;i<nsam;i++)
	{
		memcpy(fft+i*nsamb,buf+i*nsam,sizeof(float)*nsam);
	}
}

void fft2buf(float *buf, float *fft, int nsam)
{
	int nsamb=nsam+2;
	int i;
	for(i=0;i<nsam;i++)
	{
		memcpy(buf+i*nsam,fft+i*nsamb,sizeof(float)*nsam);
	}
}


void MinMaxMean(float *buf, size_t size, float &min, float &max, float &mean)
{
	if(buf==NULL || size<=0) return; 
	min=1e20;
	max=-1e20;
	double dmean=0.0;
	
	size_t i;
	for(i=0;i<size;i++)
	{
		if(min>buf[i]) min=buf[i];
		if(max<buf[i]) max=buf[i];
		dmean+=buf[i];
	}
	mean=dmean/size;
}

float STD(float *buf, size_t size, float mean)
{
	if(buf==NULL || size<=0) return 0; 
	
	double std=0.0;
	double val;
	size_t i;
	for(i=0;i<size;i++)
	{
		val=buf[i]-mean;
		std+=val*val;
	}
	
	return sqrt(std/size);
}

void shift2d_phase(float *buf, int nsam, float sx, float sy)
{
	CPLX *bufc=(CPLX *)buf;
	
	CPLX pshift;
	float shift;
	int hnsam=nsam/2;
	int hnsamb=hnsam+1;

	int i,j,jj,id,is;
	float shx=sx*6.2831852/nsam;
	float shy=sy*6.2831852/nsam;
	for(j=0;j<hnsam;j++)
	{
		id=j*hnsamb;
		is=(j+hnsam)*hnsamb;
		jj=j-hnsam;
		for(i=0;i<hnsamb;i++)
		{
			shift=i*shx+j*shy;
			pshift.x=cos(shift);
			pshift.y=sin(shift);
			bufc[id+i]=cXc(bufc[id+i],pshift);
			
			shift=i*shx+jj*shy;
			pshift.x=cos(shift);
			pshift.y=sin(shift);
			bufc[is+i]=cXc(bufc[is+i],pshift);
		}
	}
}


void CPLXAdd(CPLX *dst, CPLX *src, int sizec)
{
	for(int i=0;i<sizec;i++)
	{
		dst[i].x+=src[i].x;
		dst[i].y+=src[i].y;
	}
}


void FFTDataToLogDisp(float *pIn, float *pOut, int dim)
{
	//convert fft_complex to Intensity
	int fftsize=(dim+2)*dim/2;
	int hdim=dim/2;
	int width=(dim+2)/2;
	int id0=dim*dim/2+dim/2;
	int id1=dim*dim/2+dim/2;
	int id2=dim/2-dim*dim/2;
	int id3=3*dim*dim/2+dim/2;
	//pIn[0]=0;
	float absval;
	int i,j,is;

	int index;
	float val;
	float fs=0.00005;

	for(i=0;i<dim;i++)
	{
		is=i*width;
		for(j=0;j<hdim;j++)
		{
			absval=pIn[is+j];
			val=log(1+fs*absval);

			index=i*dim+j;
			if(i<hdim) pOut[id0+index]=pOut[id1-index]=val;
			else if(i>hdim) pOut[id2+index]=pOut[id3-index]=val;
		}
		pOut[i*dim]=pOut[i*dim+1];
	}

}

void rmNoiseCC(float* hboxmap, int box, int nframe, int peakR)
{
	int peakD=peakR*2-1;
	float *peak=new float[peakD*peakD];
	float *buf=0;
	memset(peak,0,sizeof(float)*peakD*peakD);

	int i,j,k;
	int offset=box/2-peakR+1;
	for(k=0;k<nframe;k++)
	{
		buf=hboxmap+box*box*k;
		for(j=0;j<peakD;j++)
			for(i=0;i<peakD;i++)
			{
				peak[j*peakD+i]+=buf[(j+offset)*box+i+offset]/nframe;
			}
	}

	for(k=0;k<nframe;k++)
	{
		buf=hboxmap+box*box*k;
		for(j=0;j<peakD;j++)
			for(i=0;i<peakD;i++)
			{
				buf[(j+offset)*box+i+offset]-=peak[j*peakD+i];
			}
	}

	delete [] peak;

}

void cosmask2d(complex<float> *pfft, int nsam)
{
	int i,j,ii,jj,id;
	
	int hnsam=nsam/2;
	int hnsamb=hnsam+1;
	float step=3.141592653589793/hnsam;
	int r;
	for(j=0;j<hnsam;j++)
		for(i=0;i<hnsamb;i++)
		{
			id=j*hnsamb+i;
			r=sqrt(float(i*i+j*j));
			if(r<hnsam)	pfft[id]*=cos(r*step)/2+0.5;
			else pfft[id]=0;

			id=(j+hnsam)*hnsamb+i;
			jj=j-hnsam;
			r=sqrt(float(i*i+jj*jj));
			if(r<hnsam)	pfft[id]*=cos(r*step)/2+0.5;
			else pfft[id]=0;
		}
}

void histMinMax(float *buf, int size, float &min, float &max, float threshold)
{
	int nbin=400;
	int *hist=new int[nbin];
	memset(hist,0,sizeof(int)*nbin);

	//find real min and max
	min=1e20;
	max=-1e20;
	int i,id;
	for(i=0;i<size;i++)
	{
		if(buf[i]<min) min=buf[i];
		if(buf[i]>max) max=buf[i];
	}
	if(min>=max) return;

	//get histogram
	float step=(max-min)/nbin;
	for(i=0;i<size;i++)
	{
		id=int((buf[i]-min)/step);
		if(id>=nbin) id=nbin-1;
		hist[id]++;
	}

	//remove 1% 
	int percent=int(size*threshold);
	int count=0;
	for(i=0;i<nbin;i++)
	{
		count+=hist[i];
		if(count>percent)
		{
			min+=i*step;
			break;
		}
	}
	count=0;
	for(i=nbin-1;i>0;i--)
	{
		count+=hist[i];
		if(count>percent)
		{
			max-=(nbin-1-i)*step;
			break;
		}
	}


	delete [] hist;
}

void buf2Disp(char *disp, float *buf, int size)
{
	float min,max;
	int i;
	histMinMax(buf,size, min, max);
	if(min>=max) return;

	float scale=255.0/(max-min);
	for(i=0;i<size;i++)
	{
		if(buf[i]<=min) disp[i]=0;
		else if(buf[i]>=max) disp[i]=255;
		else disp[i]=(buf[i]-min)*scale;
	}
}

void buf2DispShort(short *disp, float *buf, int size)
{
	float min,max;
	int i;
	histMinMax(buf,size, min, max);
	if(min>=max) return;

	float scale=255.0/(max-min);
	for(i=0;i<size;i++)
	{
		if(buf[i]<=min) disp[i]=0;
		else if(buf[i]>=max) disp[i]=255;
		else disp[i]=(buf[i]-min)*scale;
	}
}

//sizeof(pIn)=(nsam/2+1)*nsam;
//sizeof(pOut)=(nsam+2)*nsam;
void FFTModulusToDispBuf(float *pIn, float *pOut, int nsam)
{
	int hnsam=nsam/2;
	int hnsamb=nsam/2+1;
	int nsamb=nsam+2;

	int i,j,ii,jj,id0,id1;
	//bottom-right
	for(i=0;i<hnsam;i++)
	{
		memcpy(pOut+(hnsam+i)*nsamb+hnsam,pIn+i*hnsamb,sizeof(float)*hnsamb);
	}
	//top-right
	for(i=0;i<hnsam;i++)
	{
		memcpy(pOut+i*nsamb+hnsam,pIn+(hnsam+i)*hnsamb,sizeof(float)*hnsamb);
	}

	//get another half
	for(j=1;j<nsam;j++)
	{
		id0=j*nsamb;
		id1=(nsam-j)*nsamb;
		for(i=0;i<hnsam;i++)
		{
			pOut[id0+i]=pOut[id1+nsam-i];
		}
	}
	//set first half row which doesn't have corresping line
	memcpy(pOut,pOut+nsamb,sizeof(float)*hnsam);
}

void BinFFTDispBufToChar(char *pDisp, int dispdim, float *pfft, int nsam)
{
	int i;
	int hdispdim=dispdim/2;
	int dispdimb=dispdim+2;
	int nsamb=nsam+2;
	int dispsize=dispdim*dispdim;

	//fft
	fft2d(pfft,nsam);

	//bin by FFT
	float *buf=new float[(dispdim+2)*dispdim];
	float *disp=new float[dispdim*dispdim];
	for(i=0;i<hdispdim;i++)
	{
		memcpy(buf+i*dispdimb,pfft+i*nsamb,sizeof(float)*dispdimb);
		memcpy(buf+(dispdim-1-i)*dispdimb,pfft+(nsam-1-i)*nsamb,sizeof(float)*dispdimb);
	}

	//ifft disp
	ifft2d(buf,dispdim);

	//convert to char display
	fft2buf(disp,buf,dispdim);
	for(i=0;i<dispsize;i++) disp[i]/=dispsize;
	buf2Disp(pDisp,disp,dispdim*dispdim);
	

	delete [] buf;
	delete [] disp;
}

void BinFFTDispBufToChar(short *pDisp, int dispdim, float *pfft, int nsam)
{
	int i;
	int hdispdim=dispdim/2;
	int dispdimb=dispdim+2;
	int nsamb=nsam+2;
	int dispsize=dispdim*dispdim;

	//fft
	fft2d(pfft,nsam);

	//bin by FFT
	float *buf=new float[(dispdim+2)*dispdim];
	float *disp=new float[dispdim*dispdim];
	for(i=0;i<hdispdim;i++)
	{
		memcpy(buf+i*dispdimb,pfft+i*nsamb,sizeof(float)*dispdimb);
		memcpy(buf+(dispdim-1-i)*dispdimb,pfft+(nsam-1-i)*nsamb,sizeof(float)*dispdimb);
	}

	//ifft disp
	ifft2d(buf,dispdim);

	//convert to char display
	fft2buf(disp,buf,dispdim);
	for(i=0;i<dispsize;i++) disp[i]/=dispsize;
	buf2DispShort(pDisp,disp,dispdim*dispdim);
	

	delete [] buf;
	delete [] disp;
}

bool isMRC(const char * m_fp, int verb)
{
	MRC stack;
	int ret=stack.open(m_fp, "rb");
	if (ret>0) 
	{
		if (stack.m_header.map[0]=='M' && stack.m_header.map[1]=='A' && 
			stack.m_header.map[2]=='P' && stack.m_header.map[3]==' ')
		{
			printf("Auto-detected file format: MRC2000\n");
			return true;
		}
		else
			if (strstr(m_fp, "mrc") != NULL)
			{
				if (verb) printf("Auto-detected file format: old MRC with mrc extension.\n");
				return true;
			}
			else return false;
	}
	return false;
}
	
bool isMRCZ(const char * m_fp, int verb)
{
	FILE * fp=fopen(m_fp, "rb");
	if (fp==NULL) 
	{
		return false;
	}
	rewind(fp);
	int a=fgetc(fp);
	int b=fgetc(fp);
	fclose(fp);
	if (a==0x1f && b==0x8b)
	{
		MRCZ stack;
		int ret=stack.open(m_fp, "rb");
		if (ret > 0)
		{
			if (stack.m_header.map[0]=='M' && stack.m_header.map[1]=='A' && 
				stack.m_header.map[2]=='P' && stack.m_header.map[3]==' ')
			{	
				if (verb) printf("Auto-detected file format: MRC-gz\n");
				return true;
			}
		}
	}
	return false;
}
	
bool isTIFF(const char * m_fp, int verb)
{
	FILE * fp=fopen(m_fp, "rb");
	if (fp==NULL) 
	{
		return false;
	}
	rewind(fp);
	int a=fgetc(fp);
	int b=fgetc(fp);
	fclose(fp);
	if (a=='I' && b=='I' || a=='M' && b=='M')
	{
		if (verb) printf("Auto-detected file format: TIFF\n");
		return true;
	}
	return false;
}
	
