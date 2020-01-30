#include "DFAlign.h"
#include "func.h"
#include "mfunc.h"
#include "safefft.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


CDFAlign::CDFAlign(void)
{
	m_nsam=0;
	m_nsamRaw=0;
	m_iter=0;
	m_iterthres=0.0;
	//For image output
	m_bufIm=new float[(DISPDIM+2)*DISPDIM];
	m_dispIm=new short[DISPDIM*DISPDIM];

	m_bufFFTCorr=0;
	m_dispFFTCorr=new short[DISPDIM*DISPDIM];

	m_bufFFTRaw=0;
	m_dispFFTRaw=new short[DISPDIM*DISPDIM];

	m_bufCCMap=0;;

	initFFTWLock();

}


CDFAlign::~CDFAlign(void)
{
	if(m_bufIm!=0) delete [] m_bufIm;
	if(m_dispIm!=0) delete [] m_dispIm;

	if(m_bufFFTCorr!=0) delete [] m_bufFFTCorr;
	if(m_dispFFTCorr!=0) delete [] m_dispFFTCorr;

	if(m_bufFFTRaw!=0) delete [] m_bufFFTRaw;
	if(m_dispFFTRaw!=0) delete [] m_dispFFTRaw;

	if(m_bufCCMap!=0) delete [] m_bufCCMap;

	freeFFTWLock();
}


void CDFAlign::Message(const char *str)
{
	printf("%s\n",str);
}

void CDFAlign::UpdateDisplay()
{}

void CDFAlign::TextOutput(const char *str)
{
	m_log=str;
	printf("%s",m_log.c_str());
	//SendMessage(m_dlgwnd, WM_TSHOWLOG, 0,0);
	
	if(m_para.bSaveLog)
	{
		FILE *fp=fopen(m_fnLog,"a");
		fprintf(fp,"%s",m_log.c_str());
		fclose(fp);
	}
}

void* CDFAlign::ImageOutputThread(void *p)
{
	CDFAlign *pThis=(CDFAlign *)p;

	ifft2d(pThis->m_bufIm,DISPDIM);

	int size=DISPDIM*DISPDIM;
	float *buf=new float[size];
	fft2buf(buf,pThis->m_bufIm,DISPDIM);
	buf2DispShort(pThis->m_dispIm, buf, size);

	//add output code here
	//SendMessage(pThis->m_dlgwnd, WM_TSHOWIMAGE, 0,0);
	MRC mrc;
	mrc.open(pThis->m_dispCorrSum,"wb");
	mrc.createMRC(pThis->m_dispIm,DISPDIM,DISPDIM,1);
	mrc.close();
	
	delete [] buf;
	
	return (void *)0;
}
void CDFAlign::ImageOutput(float *buf)
{
	memcpy(m_bufIm,buf,sizeof(float)*(DISPDIM+2)*DISPDIM);
	
	pthread_t tid;
	int terror;
	terror=pthread_create(&tid,NULL,ImageOutputThread,(void *)this);
  	if(terror!=0)
  	{
		TextOutput("Error: Failed to create pthread: Image Output\n");;
   	return;
   }
   m_tids.push_back(tid);

}

void* CDFAlign::FFTOutputCorrThread(void *p)
{
	CDFAlign *pThis=(CDFAlign *)p;
	int nsam=pThis->m_nsam;
	float *buf=new float[(nsam+2)*nsam];
	FFTModulusToDispBuf(pThis->m_bufFFTCorr, buf, nsam);
	BinFFTDispBufToChar(pThis->m_dispFFTCorr, DISPDIM, buf, nsam);
	delete [] buf;

	//add output code here
	//SendMessage(pThis->m_dlgwnd, WM_TSHOWFFTCORR, 0,0);
	MRC mrc;
	mrc.open(pThis->m_dispCorrFFT,"wb");
	mrc.createMRC(pThis->m_dispFFTCorr,DISPDIM,DISPDIM,1);
	mrc.close();
	
	return (void *)0;
}
void CDFAlign::FFTOutputCorr(float *buf)
{
	if(m_bufFFTCorr==0) m_bufFFTCorr=new float[(m_nsam/2+1)*m_nsam];
	memcpy(m_bufFFTCorr,buf,sizeof(float)*(m_nsam/2+1)*m_nsam);
	
	pthread_t tid;
	int terror;
	terror=pthread_create(&tid,NULL,FFTOutputCorrThread,(void *)this);
  	if(terror!=0)
  	{
		TextOutput("Error: Failed to create pthread: FFT Output Corr\n");;
   	return;
   }
   m_tids.push_back(tid);
}

void* CDFAlign::FFTOutputRawThread(void *p)
{
	CDFAlign *pThis=(CDFAlign *)p;
	int nsam=pThis->m_nsam;
	float *buf=new float[(nsam+2)*nsam];
	FFTModulusToDispBuf(pThis->m_bufFFTRaw, buf, nsam);
	BinFFTDispBufToChar(pThis->m_dispFFTRaw, DISPDIM, buf, nsam);
	
	//add output code here
	//SendMessage(pThis->m_dlgwnd, WM_TSHOWFFTRAW, 0,0);
	MRC mrc;
	mrc.open(pThis->m_dispRawFFT,"wb");
	mrc.createMRC(pThis->m_dispFFTRaw,DISPDIM,DISPDIM,1);
	mrc.close();

	delete [] buf;
	
	return (void *)0;
}
void CDFAlign::FFTOutputRaw(float *buf)
{
	if(m_bufFFTRaw==0) m_bufFFTRaw=new float[(m_nsam/2+1)*m_nsam];
	memcpy(m_bufFFTRaw,buf,sizeof(float)*(m_nsam/2+1)*m_nsam);
	
	//pthread_t tid;
	//int terror;
	//terror=pthread_create(&tid,NULL,FFTOutputRawThread,(void *)this);
  	//if(terror!=0)
  	//{
	//	TextOutput("Error: Failed to create pthread: FFT Output Raw\n");;
   	//return;
        //}
        //m_tids.push_back(tid);
        FFTOutputRawThread( (void *) this);
}

void* CDFAlign::CCMapOutputThread(void *p)
{
	CDFAlign *pThis=(CDFAlign *)p;
	
	//add output code here
	//SendMessage(pThis->m_dlgwnd, WM_TSHOWCCMAP, 0,0);
	
	return (void *)0;
}
void CDFAlign::CCMapOutput(float *buf, void *pki)
{
	Vector<double> &ki=*(Vector<double> *)pki;
	if(m_bufCCMap!=0) delete [] m_bufCCMap;
	int size=m_para.CCPeakSearchDim*m_para.CCPeakSearchDim*ki.size();
	m_bufCCMap=new float[size];
	memcpy(m_bufCCMap,buf,sizeof(float)*size);
	m_kiCCMap.clear();
	for(int i=0;i<ki.size();i++) m_kiCCMap.push_back(ki[i]);

	pthread_t tid;
	int terror;
	terror=pthread_create(&tid,NULL,CCMapOutputThread,(void *)this);
  	if(terror!=0)
  	{
		TextOutput("Error: Failed to create pthread: CC Map Output\n");;
   	return;
   }
   m_tids.push_back(tid);

}

void CDFAlign::PlotFSC(float2* hRaw0, float2 *hRaw1, float2 *hCorr0, float2 *hCorr1,
					MASK *pPosList, int nsam, complex<double> direction)
{
	const int step=nsam/400;
	

	int nsamc=nsam/2+1;
	int sizec=nsamc*nsam;
	int nbox=nsamc/step+1;

	int i,id;
	float r,angle;
	float edge=cos(PI/4);

	float3 *fRaw0=new float3[nbox];   // x:cos(phase) y:amp^2 z:amp^2, along drift
	float3 *fRaw1=new float3[nbox];   // x:cos(phase) y:amp^2 z:amp^2, perpendicular to drift
	float3 *fCorr0=new float3[nbox];   // x:cos(phase) y:amp^2 z:amp^2, along drift
	float3 *fCorr1=new float3[nbox];   // x:cos(phase) y:amp^2 z:amp^2, perpendicular to drift
	memset(fRaw0,0,sizeof(float3)*nbox);
	memset(fRaw1,0,sizeof(float3)*nbox);
	memset(fCorr0,0,sizeof(float3)*nbox);
	memset(fCorr1,0,sizeof(float3)*nbox);

	if(abs(direction)>0.1) direction/=abs(direction);
	else direction=1.0;

	cuComplex a,b;

	
	for(i=0;i<sizec;i++)
	{
		r=sqrt(float(pPosList[i].x*pPosList[i].x+pPosList[i].y*pPosList[i].y));
		if(int(r)<=0 || r>=nsamc) continue;
		id=int(r/step);

		angle=fabs(pPosList[i].x*direction.real()+pPosList[i].y*direction.imag())/r;
		
		if(angle>edge) //along drift
		{
			a=hRaw0[i];
			b=hRaw1[i];
			fRaw0[id].x+=a.x*b.x+a.y*b.y;
			fRaw0[id].y+=a.x*a.x+a.y*a.y;
			fRaw0[id].z+=b.x*b.x+b.y*b.y;

			a=hCorr0[i];
			b=hCorr1[i];
			fCorr0[id].x+=a.x*b.x+a.y*b.y;
			fCorr0[id].y+=a.x*a.x+a.y*a.y;
			fCorr0[id].z+=b.x*b.x+b.y*b.y;
		}
		else
		{
			a=hRaw0[i];
			b=hRaw1[i];
			fRaw1[id].x+=a.x*b.x+a.y*b.y;
			fRaw1[id].y+=a.x*a.x+a.y*a.y;
			fRaw1[id].z+=b.x*b.x+b.y*b.y;

			a=hCorr0[i];
			b=hCorr1[i];
			fCorr1[id].x+=a.x*b.x+a.y*b.y;
			fCorr1[id].y+=a.x*a.x+a.y*a.y;
			fCorr1[id].z+=b.x*b.x+b.y*b.y;
		}
	}

	m_fscRaw0.resize(nbox,0.0);
	m_fscRaw1.resize(nbox,0.0);
	m_fscCorr0.resize(nbox,0.0);
	m_fscCorr1.resize(nbox,0.0);


	float t;
	for(i=0;i<nbox;i++)
	{
		t=fRaw0[i].y*fRaw0[i].z;
		if(t>0.00001) m_fscRaw0[i]=complex<double>(double(i)/nbox,fRaw0[i].x/sqrt(t));
		else m_fscRaw0[i]=complex<double>(double(i)/nbox,0.0);

		t=fRaw1[i].y*fRaw1[i].z;
		if(t>0.00001) m_fscRaw1[i]=complex<double>(double(i)/nbox,fRaw1[i].x/sqrt(t));
		else m_fscRaw1[i]=complex<double>(double(i)/nbox,0.0);

		t=fCorr0[i].y*fCorr0[i].z;
		if(t>0.00001) m_fscCorr0[i]=complex<double>(double(i)/nbox,fCorr0[i].x/sqrt(t));
		else m_fscCorr0[i]=complex<double>(double(i)/nbox,0.0);

		t=fCorr1[i].y*fCorr1[i].z;
		if(t>0.00001) m_fscCorr1[i]=complex<double>(double(i)/nbox,fCorr1[i].x/sqrt(t));
		else m_fscCorr1[i]=complex<double>(double(i)/nbox,0.0);
	}

	//add output code here
	//SendMessage(m_dlgwnd, WM_TSHOWFSC, 0,0);

	//output
	
	char str[512];
	TextOutput("\nFSC parallel(D) and perpendicular(U) to drift direction:\n");
	TextOutput("   Nq%       D_Raw       U_Raw       D_Corr       U_Corr\n");
	for(i=0;i<nbox;i++)
	{
		sprintf(str,"%7.2f %12.4f %12.4f %12.4f %12.4f\n",i*100.0/nbox,
			m_fscRaw0[i].imag(),m_fscRaw1[i].imag(),m_fscCorr0[i].imag(),m_fscCorr1[i].imag());
		TextOutput(str);
	}
	TextOutput("\n");


	delete [] fRaw0;
	delete [] fRaw1;
	delete [] fCorr0;
	delete [] fCorr1;
}

void CDFAlign::Done()
{
	   
	//add output code here
	//SendMessage(m_dlgwnd, WM_TDONE, 0,0);
	
}

void CDFAlign::PlotOutput(vector<complex<double> > &xy)
{
	m_curve=xy;
	//add output code here
	//SendMessage(m_dlgwnd, WM_TSHOWCURVE, 0,0);

}

int CDFAlign::getNFrame()
{
	MRC mrc;
	MRCZ mrcz;
	MRCT mrct;
	if(isMRC(m_fnStack) && mrc.open(m_fnStack,"rb") > 0) 
	{
		int n=mrc.getNz();
		mrc.close();
		return n;
	}
	else if (isMRCZ(m_fnStack) && mrcz.open(m_fnStack,"rb") > 0) 
        {
                int n=mrcz.getNz();
                mrcz.close();
                return n;
        }
	else if (isTIFF(m_fnStack) && mrct.open(m_fnStack,"rb") > 0) 
        {
                int n=mrct.getNz();
                mrct.close();
                return n;
        }
	return 0;
}
int CDFAlign::getNFrame(const char *filename)
{
	MRC mrc;
	MRCZ mrcz;
	MRCT mrct;
	if(isMRC(filename) && mrc.open(filename,"rb") > 0) 
	{
		int n=mrc.getNz();
		mrc.close();
		return n;
	}
	else if (isMRCZ(filename) && mrcz.open(filename,"rb") > 0)
        {
                int n=mrcz.getNz();
                mrcz.close();
                return n;
        }
	else if (isTIFF(filename) && mrct.open(filename,"rb") > 0)
        {
                int n=mrct.getNz();
                mrct.close();
                return n;
        }
	return 0;
}

MRCHeader CDFAlign::getMrcHeader(const char *filename)
{
	MRC mrc;
	MRCZ mrcz;
	MRCT mrct;
	MRCHeader header;
	memset(&header,0,sizeof(MRCHeader));
	if(isMRC(filename) && mrc.open(filename,"rb") > 0)
	{
		mrc.getHeader(&header);;
		mrc.close();
		return header;
	}
	else if (isMRCZ(filename) && mrcz.open(filename,"rb") > 0) 
        {
                mrcz.getHeader(&header);;
                mrcz.close();
                return header;
        }
	else if (isTIFF(filename) && mrct.open(filename,"rb") > 0) 
        {
                mrct.getHeader(&header);;
                mrct.close();
                return header;
        }
	return header;
}

void CDFAlign::RunAlign()
{
	// TODO: Add your control notification handler code here
	//UpdateData(true);
	
	m_tids.clear();
	
	pthread_t tid;
	int terror;
	terror=pthread_create(&tid,NULL,ThreadFunc_cuAlign,(void *)this);
  	if(terror!=0)
  	{
		TextOutput("Error: Failed to create pthread: Align\n");
		return;
   }
   m_tids.push_back(tid);
   
   //wait for finish
   void *TReturn;
   int i;
	for(i=0;i<m_tids.size();i++)
	{
   	terror=pthread_join(m_tids[i],&TReturn);
   	if(terror!=0)
   	{
      	TextOutput("Warnning: Thread doesn't exit. Something may be wrong.\n");
   	}
   	
   }
   m_tids.clear();
}


void* CDFAlign::ThreadFunc_cuAlign(void* p)
{
	CDFAlign *pThis=(CDFAlign *)p;
	APARA &para=pThis->m_para;
	pThis->m_bRun=true;

	char str[512];
	int j;

	//open stack file
	//test file type
	//
	MRC *  stack;
	
	if (isTIFF(pThis->m_fnStack, 1))
		stack = new MRCT(pThis->m_fnStack,"rb");
	else if (isMRCZ(pThis->m_fnStack, 1))
		stack = new MRCZ(pThis->m_fnStack,"rb");
	else if (isMRC(pThis->m_fnStack, 1))
		stack = new MRC(pThis->m_fnStack,"rb");
	else 
	{
                Message("Only MRC, MRC-gz, TIFF file formats are supported.");
                pThis->m_bRun=false;
                return (void *)0;
        }


	if (stack == NULL)
	{
		sprintf(str,"Error: Failed to open stack %s .",pThis->m_fnStack);
		Message(str);
		pThis->m_bRun=false;
		return (void *)0;
	}

	//get image size
	int nx=stack->getNx();
	int ny=stack->getNy();
	int nz=stack->getNz();
	sprintf(str,"Input Stack: Nx(%d) Ny(%d) Nz(%d) Mode(%d)\n\n",nx,ny,nz,stack->getMode());
	pThis->TextOutput(str);

	
	float bin=para.bin;
	if(bin<=0) return (void *)0;

	float xscale=1.0/para.xscale;
	float yscale=1.0/para.yscale;

	int offsetx=para.crop_offsetx;
	int offsety=para.crop_offsety;
	int nsamUnbin=verifyCropSize(nx,ny, offsetx, offsety, para.crop_nsam,bin);
	if(nsamUnbin<=0)
	{
		Message("Error: Wrong image Size.");
		pThis->m_bRun=false;
		return (void *)0;
	}
	int nsam=(int)((float)nsamUnbin/bin)/2*2;
	int nsamb=nsam+2;
//	if(bin==1) sprintf(str,"Crop Image: Offset(%d %d) Dim(%d)\n",offsetx,offsety,nsamUnbin);
//	else sprintf(str,"Crop Image: Offset(%d %d) RawDim(%d) BinnedDim(%d)\n",offsetx,offsety,nsamUnbin,nsam);
	sprintf(str,"Crop Image: Offset(%d %d) RawDim(%d) BinnedDim(%d)\n",offsetx,offsety,nsamUnbin,nsam);
	pThis->TextOutput(str);
	pThis->m_nsam=nsam;
	pThis->m_nsamRaw=(int)((float)nx/bin);

	//allocate memeory	
	size_t size=nsam*nsam;
	size_t sizeb=nsamb*nsam;
	int sizebUnbin=(nsamUnbin+2)*nsamUnbin;
	if(para.nStart<0) para.nStart=0;
	if(para.nEnd>=nz) para.nEnd=nz-1;
	int nframe=para.nEnd-para.nStart+1;
	pThis->UpdateDisplay();

	//host memory
	float *bufmrc=new float[nx*ny];
	float *bufmrc2=new float[nx*ny]; // dark/gain corrected
//	float *bufscale=new float[nx*ny]; // temp buffer for image stretching
//	float *bufscale2=new float[nx*ny]; // temp buffer for image stretching
	float *bufmrcfft=new float[sizebUnbin];
	float *bufdark=new float[nx*ny];
	float *bufnorm=new float[nx*ny];
	float *bufdark2=new float[nx*ny];
	float *bufnorm2=new float[nx*ny];

	float *htmp=new float[sizeb];
	float *hbuf=new float[sizeb*nframe];  //host memory for entir stack
	float *hbuf2=new float[sizeb*nframe];  //host memory for entire stack, dark/gain corrected
	float *hdisp=new float[sizeb];
	float *hFSCRaw0=new float[sizeb];  //even number
	float *hFSCRaw1=new float[sizeb];  //odd
	float *hFSCCorr0=new float[sizeb];  //even number
	float *hFSCCorr1=new float[sizeb];   //odd
	sprintf(str,"Allocate host memory: %f Gb\n",(6*nx*ny+sizeb*(2*nframe+6)+sizebUnbin)/256.0/1024.0/1024.0);
	pThis->TextOutput(str);
	if(hbuf==0 || hbuf2==0)
	{
		if(bufmrc!=NULL) delete [] bufmrc;
		if(bufmrc2!=NULL) delete [] bufmrc2;
		Message("Failed to allocate host memeory.");
		pThis->m_bRun=false;
		return (void *)0;
	}


	//device memory
	bool success=initGPU(para.GPUNum);
	if(!success)
	{
		sprintf(str,"Failed to initialize GPU #%d.",para.GPUNum);
		Message(str);
		delete [] bufmrc;
		delete [] hbuf;
		delete [] bufmrc2;
//		delete [] bufscale;
//		delete [] bufscale2;
		delete [] hbuf2;
		delete [] bufdark;
		delete [] bufnorm;
		delete [] bufdark2;
		delete [] bufnorm2;

		pThis->m_bRun=false;
		return (void *)0;
	}
	

	float *dsum=0;
	float *dsumcorr=0;
	float *dfft=0;
	float *dtmp=0;
	GPUMemAlloc((void **)&dsum,sizeof(float)*sizeb);	
	GPUMemAlloc((void **)&dsumcorr,sizeof(float)*sizeb);	
	GPUMemAlloc((void **)&dtmp,sizeof(float)*sizeb);
	cufftHandle fft_plan,ifft_plan;
	
	//prepare fft for unbinned image
	fft_plan=GPUFFTPlan(nsamUnbin);
	GPUSync();
	GPUMemAlloc((void **)&dfft,sizeof(float)*sizebUnbin);


	//make a list 
	int sizec=(nsam/2+1)*nsam;
	MASK *hPosList=new MASK[sizec];
	MASK *dPosList=0;
	MkPosList(hPosList,nsam,para.bfactor);
	GPUMemAlloc((void **)&dPosList,sizeof(MASK)*sizec);
	GPUMemH2D((void **)dPosList,(void **)hPosList,sizeof(MASK)*sizec);

	size_t theFree, theTotal;
	GPUMemCheck(theFree,theTotal);
	sprintf(str,"GPU memory:  free:%.0fMb    total:%.0fMb\n", theFree/1024.0/1024.0, theTotal/1024.0/1024.0);
	pThis->TextOutput(str);

//Read dark/gain reference
	MRC darkref;
	if(strlen(pThis->m_fnDark))
	{
		if(darkref.open(pThis->m_fnDark,"rb")<=0)
		{
			sprintf(str,"Error: Failed to open dark reference %s .",pThis->m_fnDark);
			Message(str);
			pThis->m_bRun=false;
			return (void *)0;
		}
		if (nx!=darkref.getNx() || ny!=darkref.getNy())
		{
			sprintf(str,"Error: Image dimension of dark reference %s differs from image stack .",pThis->m_fnDark);
			Message(str);
			pThis->m_bRun=false;
			return (void *)0;
		}
		if(darkref.read2DIm_32bit(bufdark,0)!=darkref.getImSize())
		{
			sprintf(str,"Error: Failed to read dark reference %s .",pThis->m_fnDark);
			Message(str);
			pThis->m_bRun=false;
			return (void *)0;
		}
		if (darkref.getNz() > 1)
		{
			Message("Processing two dark references.");
			darkref.read2DIm_32bit(bufdark2,1);
		}
		else
		{
			Message("Processing single dark reference.");
			//darkref.read2DIm_32bit(bufdark2,0);
			for(int k=0;k<nx*ny;k++)
			{
				bufdark2[k] = 0.0;
			}
		}
		
	}
	MRC normref;
	if(strlen(pThis->m_fnNorm))
	{
		if(normref.open(pThis->m_fnNorm,"rb")<=0)
		{
			sprintf(str,"Error: Failed to open gain reference %s .",pThis->m_fnNorm);
			Message(str);
			pThis->m_bRun=false;
			return (void *)0;
		}
		if (nx!=normref.getNx() || ny!=normref.getNy())
		{
			sprintf(str,"Error: Image dimension of gain reference %s differs from image stack .",pThis->m_fnNorm);
			Message(str);
			pThis->m_bRun=false;
			return (void *)0;
		}
		if(normref.read2DIm_32bit(bufnorm,0)!=normref.getImSize())
		{
			sprintf(str,"Error: Failed to read gain reference %s .",pThis->m_fnNorm);
			Message(str);
			pThis->m_bRun=false;
			return (void *)0;
		}
		if (normref.getNz() > 1)
		{
			Message("Processing two norm references.");
			normref.read2DIm_32bit(bufnorm2,1);
		}
		else
		{
			Message("Processing single norm reference.");
			//normref.read2DIm_32bit(bufnorm2,0);
			for(int k=0;k<nx*ny;k++)
			{
				bufnorm2[k] = 1.0;
			}
		}

	}
//End reading dark/gain reference
	//Read stack
	pThis->TextOutput("\nRead stack:\n");
	
	float sx=0;
	float sy=0;
	float shiftx,shifty,cc;
	float avgcc=0.0;
	bool bFSCEven=true;

	//1. calculate sum
	GPUMemZero((void **)&dsum,sizeof(float)*sizeb);
	GPUSync();
	GPUMemZero((void **)&dsumcorr,sizeof(float)*sizeb);
	GPUSync();
	for(j=para.nStart;j<=para.nEnd;j++)
	{
		//read from file and crop
		if(stack->read2DIm_32bit(bufmrc,j)!=stack->getImSize())
		{
			sprintf(str,"Error when reading #%03d\n",j);
			pThis->TextOutput(str);
		}
		// Do dark/gain correction
		if(strlen(pThis->m_fnNorm) && strlen(pThis->m_fnDark)) 
		{
			pThis->TextOutput("Correcting dark/gain...\n");
			for(int k=0;k<nx*ny;k++)
			{
				bufmrc2[k]=(bufmrc[k]-bufdark2[k])*bufnorm2[k];
				bufmrc[k]=(bufmrc[k]-bufdark[k])*bufnorm[k];
			}
		}
		else if (strlen(pThis->m_fnNorm))
		{
			pThis->TextOutput("Correcting gain...\n");
			for(int k=0;k<nx*ny;k++)
			{
				bufmrc2[k]=bufmrc[k]*bufnorm2[k];
				bufmrc[k]=bufmrc[k]*bufnorm[k];
			}
		}
		else if (strlen(pThis->m_fnDark))
		{
			pThis->TextOutput("Correcting dark...\n");
			for(int k=0;k<nx*ny;k++)
			{
				bufmrc2[k]=bufmrc[k]-bufdark2[k];
				bufmrc[k]=bufmrc[k]-bufdark[k];
			}
		}
		else
		{
//			for(int k=0;k<nx*ny;k++) bufmrc2[k]=bufmrc[k];
			memcpy(bufmrc2, bufmrc, sizeof(float)*nx*ny);
		}
// Stretch image
		if (xscale<1.0 || yscale<1.0)
		{
			int addr_orig, addr_new;
			for (int ky=ny-1; ky>0; ky--)
			{
				for (int kx=nx-1; kx>0; kx--)
				{
					// Nearest neighbor interpolation.
					addr_orig = int(ky*yscale+0.5)*nx+int(kx*xscale+0.5);
					addr_new = ky*nx+kx;
					bufmrc[addr_new]=bufmrc[addr_orig];
					bufmrc2[addr_new]=bufmrc2[addr_orig];
				}
			}
		}
// FFT the dark/gain corrected image.
		crop2fft(bufmrc2,nx,ny,bufmrcfft,offsetx,offsety,nsamUnbin,bin);
		//copy to GPU
		GPUMemH2D((void *)dfft,(void *)bufmrcfft,sizeof(float)*sizebUnbin);
		//do fft
		GPUFFT2d(dfft,fft_plan);
		GPUSync();
		//do binning
		if(bin>1.0001)
		{
			GPUMemBinD2D(dtmp, dfft, nsam, nsamUnbin);
			GPUMemD2D(dfft, dtmp, sizeof(float)*sizeb);
		}
		//Sum
//		if(j>=para.nStartSum && j<=para.nEndSum)
//		{
//			if(bFSCEven) GPUAdd(dsum,dfft,sizeb);
//			else GPUAdd(dsumcorr,dfft,sizeb);
//			bFSCEven=!bFSCEven;
//		}
		//copy ffted image to host
		GPUMemD2H((void *)(hbuf2+(j-para.nStart)*sizeb),(void *)dfft,sizeof(float)*sizeb);
		GPUSync();
// 

		crop2fft(bufmrc,nx,ny,bufmrcfft,offsetx,offsety,nsamUnbin,bin);
		
		//copy to GPU
		GPUMemH2D((void *)dfft,(void *)bufmrcfft,sizeof(float)*sizebUnbin);
		//do fft
		GPUFFT2d(dfft,fft_plan);
		GPUSync();

		//do binning
		if(bin>1.0001)
		{
			GPUMemBinD2D(dtmp, dfft, nsam, nsamUnbin);
			GPUMemD2D(dfft, dtmp, sizeof(float)*sizeb);
		}

		//Sum
		if(j>=para.nStartSum && j<=para.nEndSum)
		{
			if(bFSCEven) GPUAdd(dsum,dfft,sizeb);
			else GPUAdd(dsumcorr,dfft,sizeb);
			bFSCEven=!bFSCEven;
		}
		//copy ffted image to host
		GPUMemD2H((void *)(hbuf+(j-para.nStart)*sizeb),(void *)dfft,sizeof(float)*sizeb);
		GPUSync();

		sprintf(str,"......Read and sum frame #%03d   mean:%f\n",j,(hbuf+(j-para.nStart)*sizeb)[0]/nsam/nsam);
		pThis->TextOutput(str);
	}
	GPUMemD2H((void *)hFSCRaw0,(void *)dsum,sizeof(float)*sizeb);
	GPUMemD2H((void *)hFSCRaw1,(void *)dsumcorr,sizeof(float)*sizeb);
	GPUAdd(dsum,dsumcorr,sizeb);
	GPUSync();
	

	//free memory for unbined image
	delete [] bufmrcfft;
	delete stack;
	bufmrcfft=0;
	GPUMemFree((void **)&dfft);
	GPUFFTDestroy(fft_plan);
	fft_plan=0;	
	//finish GPU memory allocate
	GPUMemAlloc((void **)&dfft,sizeof(float)*sizeb);
	GPUMemZero((void **)&dsumcorr,sizeof(float)*sizeb);
	GPUSync();
	ifft_plan=GPUIFFTPlan(nsam);
	GPUSync();

	//Make fft modulus for display
	if(para.bDispFFTRaw)
	{
		GPUSync();
		GPUFFTLogModulus(dfft, dsum, nsam, para.fftscale);
		GPUSync();
		GPUMemD2H((void *)hdisp,(void *)dfft,sizeof(float)*(nsam/2+1)*nsam);
		//pThis->FFTOutputRaw(hdisp);   //has been move to below
	}
	//copy sum image to host for save and display
	if(para.bDispFFTRaw || para.bSaveRawSum)
	{
		GPUIFFT2d(dsum,ifft_plan);
		GPUSync();
		GPUMultiplyNum(dsum,1.0/size,sizeb);
		GPUMemD2H((void *)htmp,(void *)dsum,sizeof(float)*sizeb);
		fft2buf(bufmrc,htmp,nsam);
	}
	//save
	MRC mrcraw;
	if(para.bSaveRawSum)
	{
		//write to file
		mrcraw.open(pThis->m_fnRawsum,"wb");
		mrcraw.createMRC(bufmrc,nsam,nsam,1);
		//stats
		sprintf(str,"Mean=%f   Min=%f   Max=%f\n",mrcraw.m_header.dmean,mrcraw.m_header.dmin,mrcraw.m_header.dmax);
		pThis->TextOutput(str);
		mrcraw.close();
		sprintf(str,"Save Unaligned Sum to: %s\n",pThis->m_fnRawsum);
		pThis->TextOutput(str);
	}
	

	//2. frame to frame shift
	pThis->TextOutput("\nCalculate relative drift between frames\n");
	Matrix<complex<double> > A;
	vector<complex<int> > compList;
	int ncomp=OD_SetEquation_All(A,compList, nframe, para.FrameDistOffset);
	Vector<complex<double> > b=Vector<complex<double> >(ncomp);
	int box=para.CCPeakSearchDim;
	float *hboxmap=new float[box*box*ncomp];
	int par0,par1;
	for(j=0;j<ncomp;j++)
	{
		par0=compList[j].real();
		par1=compList[j].imag();
		//copy to GPU
		GPUMemH2D((void *)dsum,(void *)(hbuf+par0*sizeb),sizeof(float)*sizeb);
		GPUMemH2D((void *)dfft,(void *)(hbuf2+par1*sizeb),sizeof(float)*sizeb);
		//shift and cc
		sx=0;
		sy=0;
		GPUShiftCC(dfft, dsum, dPosList,sx, sy, nsam);
		GPUSync();
		//do ifft
		GPUIFFT2d(dfft,ifft_plan);
		GPUSync();
		//find shift
		cc=FindShift(dfft,nsam, hboxmap+j*box*box, box, shiftx, shifty, para.NoisePeakSize-1);
		b[j]=complex<double>(shiftx,shifty);
		avgcc+=cc;
		sprintf(str,"......%03d Frame #%03d VS #%03d xy-shift: %8.4f %8.4f      CC:%f\n",j,par0+para.nStart,par1+para.nStart,shiftx,shifty,cc);
		pThis->TextOutput(str);
	}

	//display the RawImageFFT here due to FFTW thread safety issue
	if(para.bDispFFTRaw)
	{
		pThis->FFTOutputRaw(hdisp);
	}


	//3. sovle overdetermined equation
	Vector<complex<double> > shift=lsSolver(A,b);
	Vector<double> ki=abs(A*shift-b);
	sprintf(str,"\n......ki: First round \n");
	pThis->TextOutput(str);
	for(j=0;j<ki.size();j++)
	{
		par0=compList[j].real();
		par1=compList[j].imag();
		sprintf(str,"......ki #%03d of Frame #%03d VS #%03d: %8.4lf \n",j+para.nStart,par0+para.nStart,par1+para.nStart,ki[j]);
		pThis->TextOutput(str);
	}
	sprintf(str,"................................Average ki: %8.4lf \n\n",sum(ki)/ki.size());
	pThis->TextOutput(str);
	//display CCMap
	if(para.bDispCCMap)
	{
		pThis->CCMapOutput(hboxmap,(void *)&ki);
	}
	//3.1 re-sovle overdetermined equation after removing large ki elments
	double kiThresh=para.kiThresh;
	vector<int> goodlist=OD_Threshold(A, b, ki, kiThresh);
	shift=lsSolver(A,b);
	ki=abs(A*shift-b);
	sprintf(str,"......ki: Second round \n");
	pThis->TextOutput(str);
	for(j=0;j<ki.size();j++)
	{
		par0=compList[goodlist[j] ].real();
		par1=compList[goodlist[j] ].imag();
		sprintf(str,"......ki #%03d of Frame #%03d VS #%03d: %8.4f \n",j+para.nStart,par0+para.nStart,par1+para.nStart,ki[j]);
		pThis->TextOutput(str);
	}
	sprintf(str,"................................Average ki: %8.4lf \n\n",sum(ki)/ki.size());
	pThis->TextOutput(str);

	//4. Do the iterative alignment
	//
	pThis->TextOutput("Begin Iterative Alignment... \n");
	bool converged=false;
	complex<double> zeroshift=0;
	complex<double> allshift=0;
	vector<complex<double> > newshift;
	newshift.push_back(zeroshift);
	for (j=0; j<shift.size();j++)
	{
		newshift.push_back(allshift+shift[j]);
		allshift+=shift[j];
	}

	vector<complex<double> > oldshift;
	oldshift=newshift;
	
	double rmsd=0;
	for (j=0; (j<pThis->m_iter) && (! converged); j++)
	{
		sprintf(str,"Iteration [%d/%d]\n", j, pThis->m_iter);
		pThis->TextOutput(str);
		//4.1 Do the sum
		//
		GPUMemH2D((void *)dsum,(void *)(hbuf),sizeof(float)*sizeb);
		GPUShift(dsum,dPosList,-oldshift[0].real(),-oldshift[0].imag(), nsam);
		sprintf(str,"......Add Frame #%03d with xy shift: %8.4lf %8.4lf\n",para.nStart,-oldshift[0].real(),-oldshift[0].imag());
		pThis->TextOutput(str);
		for (int k=para.nStart+1; k<=para.nEnd; k++)
		{
			int position=k-para.nStart;
			GPUMemH2D((void *)dfft,(void *)(hbuf+position*sizeb),sizeof(float)*sizeb);
			GPUSync();
			GPUShift(dfft,dPosList,-oldshift[position].real(),-oldshift[position].imag(), nsam);
			GPUSync();
			GPUAdd(dsum, dfft, sizeb);
			sprintf(str,"......Add Frame #%03d with xy shift: %8.4lf %8.4lf\n",k,-oldshift[position].real(),-oldshift[position].imag());
			pThis->TextOutput(str);
			GPUSync();
		}
			
		pThis->TextOutput("Sum complete... \n");

		//4.2 align to the sum
		//
		float *hboxmap2=new float[box*box*nframe];
		avgcc=0;
		for (int k=para.nStart; k<=para.nEnd; k++)
		{
			int position=k-para.nStart;
			GPUMemH2D((void *)dfft,(void *)(hbuf+position*sizeb),sizeof(float)*sizeb);
			GPUSync();
	//		sx=oldshift[position].real();
	//		sy=oldshift[position].imag();
			sx=0;
			sy=0;
			GPUShiftCC(dfft, dsum, dPosList,sx, sy, nsam);
			GPUSync();
			GPUIFFT2d(dfft,ifft_plan);
			GPUSync();
			cc=FindShift(dfft,nsam, hboxmap2+k*box*box, box, shiftx, shifty, para.NoisePeakSize-1);
			avgcc+=cc;
			newshift[position]=complex<double>(shiftx,shifty);
			complex<double> diffshift=oldshift[position]-newshift[position];
			sprintf(str,"......%03d Frame #%03d VS <sum> xy-shift: %8.4f %8.4f (%8.4f %8.4f) CC:%f\n",position,k,shiftx,shifty,diffshift.real(),diffshift.imag(),cc);
			pThis->TextOutput(str);
		}

		if (para.bSaveCCmap) 
		{
			sprintf(str, "iterCC%03d.mrc", j);
			buf2mrc(str,hboxmap2,box,box,(para.nStart-para.nEnd+1));
			sprintf(str,"Save CC map to: %s\n", str);
			pThis->TextOutput(str);
		}

		delete [] hboxmap2;
		
		//4.3 update shift and track changes
		//
		complex<double> refshift=newshift[0];
		for (int k=para.nStart; k<=para.nEnd; k++)
		{
			int position=k-para.nStart;
			complex<double> diffshift=oldshift[position]-newshift[position];
			oldshift[position]=newshift[position]-refshift;
			rmsd+=diffshift.real()*diffshift.real()+diffshift.imag()*diffshift.imag();
		}
		
		sprintf(str,"Average CC : %8.4f\n", avgcc/nframe);
		pThis->TextOutput(str);
		rmsd/=nframe;
		if (rmsd < pThis->m_iterthres*pThis->m_iterthres) 
		{
			converged=true;
			sprintf(str,"Converged at iteration [%d/%d]\n", j, pThis->m_iter);
			pThis->TextOutput(str);
		}
			
		sprintf(str,"RMSD to <sum> : %8.4f  (%8.4f for convergence)\n", sqrt(rmsd), pThis->m_iterthres);
		pThis->TextOutput(str);

	}
	sprintf(str, "Final RMSD to <sum> : %8.4f at iteration [%d/%d]\n\nShifts (delta) after iterative aligment\n", sqrt(rmsd), j-1, pThis->m_iter);
	pThis->TextOutput(str);

	//Convert oldshift to shift
	//
	sprintf(str, "...Frame #%03d xy-shift %8.4f %8.4f (%8.4f %8.4f)\n", 0,0,0,0,0);
	pThis->TextOutput(str);
	allshift=0;
	for (int k=0; k<shift.size(); k++)
	{
		allshift+=shift[k];
		complex<double> tempshift=0;
		tempshift=oldshift[k+1]-oldshift[k];
		complex<double> tempdiffshift=0;
		tempdiffshift=allshift-oldshift[k+1];
		sprintf(str, "...Frame #%03d xy-shift %8.4f %8.4f (%8.4f %8.4f)\n", k+1, oldshift[k+1].real(), oldshift[k+1].imag(), tempdiffshift.real(), tempdiffshift.imag());
		pThis->TextOutput(str);
		shift[k]=oldshift[k+1]-oldshift[k];
	}


	//output final shift
	pThis->TextOutput("Final shift:\n");
	vector<complex<double> > shiftlist;
	complex<double> totalshift=0;
	sprintf(str,"......Shift of Frame #%03d : %8.4f %8.4f\n",para.nStart,totalshift.real(),totalshift.imag());
//	sprintf(str,"......Shift of Frame #%03d : %8.4f %8.4f\n",para.nStart,oldshift[0].real(),oldshift[0].imag());
	pThis->TextOutput(str);
	shiftlist.push_back(totalshift);
//	shiftlist.push_back(oldshift[0]);
	for(j=0;j<shift.size();j++)
	{
		totalshift=totalshift+shift[j];
		sprintf(str,"......Shift of Frame #%03d : %8.4f %8.4f\n",j+para.nStart+1,totalshift.real(),totalshift.imag());
//		sprintf(str,"......Shift of Frame #%03d : %8.4f %8.4f\n",j+para.nStart+1,oldshift[j+1].real(),oldshift[j+1].imag());
		pThis->TextOutput(str);
		shiftlist.push_back(totalshift);
//		shiftlist.push_back(oldshift[j+1]);
	}
	pThis->PlotOutput(shiftlist);

	//save CCMap image
	if(para.bSaveCCmap) 
	{
		buf2mrc(pThis->m_fnCCmap,hboxmap,box,box,ncomp);
		sprintf(str,"Save CC map to: %s\n",pThis->m_fnCCmap);
		pThis->TextOutput(str);
	}
		

	
	MRC stackCorr;
	if(para.bSaveStackCorr)
	{
		stackCorr.open(pThis->m_fnStackCorr,"wb");
		stackCorr.m_header.nx=nsam;
		stackCorr.m_header.ny=nsam;
		stackCorr.m_header.nz=para.nEndSum-para.nStartSum+1;
		stackCorr.updateHeader();
	}

	//3. correct xy-shift
	int nStartSum=para.nStartSum-para.nStart;
	int nEndSum=para.nEndSum-para.nStart;
	int nEndSum2=para.nEndSum2-para.nStart;

	if(nStartSum<=0) nStartSum=0;
	if(nEndSum<=nStartSum || nEndSum>=nframe) nEndSum=nframe-1;
	if(nEndSum2<=nStartSum || nEndSum2>=nframe) nEndSum2=nframe-1;
	if(nEndSum2>nEndSum) nEndSum2=nEndSum;

	sprintf(str,"\nSum Frame #%03d - #%03d\n",nStartSum+para.nStart,nEndSum+para.nStart);
	pThis->TextOutput(str);
	//reset memory
	GPUMemZero((void **)&dsum,sizeof(float)*sizeb);
	GPUSync();
	GPUMemZero((void **)&dsumcorr,sizeof(float)*sizeb);
	GPUSync();
	//calculate middle frame shift
	complex<double> midshift=0.0;
	int RefFrame=nz/2+1;
	if(para.bAlignToMid) 
	{
		if(RefFrame<para.nStart) RefFrame=para.nStart;
		if(para.nStartSum>para.nEnd) para.nStartSum=para.nEnd;
		for(j=0;j<RefFrame-para.nStart;j++) midshift+=shift[j];
	}

	//Add(copy) first frame to GPU
	totalshift=0;
	for(j=1;j<nStartSum+1;j++)
	{
		totalshift+=shift[j-1];
	}
	GPUMemH2D((void *)dsumcorr,(void *)(hbuf+nStartSum*sizeb),sizeof(float)*sizeb);
	if(para.bAlignToMid) GPUShift(dsumcorr,dPosList,-totalshift.real()+midshift.real(),-totalshift.imag()+midshift.imag(), nsam);
	GPUSync();
	bFSCEven=false;
	sprintf(str,"......Add Frame #%03d with xy shift: %8.4lf %8.4lf\n",nStartSum+para.nStart,-totalshift.real()+midshift.real(),-totalshift.imag()+midshift.imag());
	pThis->TextOutput(str);
	//Save stack
	if(para.bSaveStackCorr)
	{
		GPUMemD2D((void *)dfft,(void *)dsumcorr,sizeof(float)*sizeb);
		GPUIFFT2d(dfft,ifft_plan);
		GPUSync();
		GPUMultiplyNum(dfft,1.0/size,sizeb);
		GPUSync();
		GPUMemD2H((void *)htmp,(void *)dfft,sizeof(float)*sizeb);
		fft2buf(bufmrc,htmp,nsam);
		stackCorr.write2DIm(bufmrc,0);
	}
	//*******
	//sum other frame
	for(j=nStartSum+1;j<=nEndSum2;j++)
	{
		totalshift+=shift[j-1];
		
		//copy to GPU
		GPUMemH2D((void *)dfft,(void *)(hbuf+j*sizeb),sizeof(float)*sizeb);
		//shift
		GPUShift(dfft,dPosList,-totalshift.real()+midshift.real(),-totalshift.imag()+midshift.imag(), nsam);
		GPUSync();
		//Sum
		if(bFSCEven) GPUAdd(dsumcorr,dfft,sizeb);
		else GPUAdd(dsum,dfft,sizeb);
		bFSCEven=!bFSCEven;

		sprintf(str,"......Add Frame #%03d with xy shift: %8.4lf %8.4lf\n",j+para.nStart,-totalshift.real()+midshift.real(),-totalshift.imag()+midshift.imag());
		pThis->TextOutput(str);

		//save stack
		if(para.bSaveStackCorr)
		{
			GPUIFFT2d(dfft,ifft_plan);
			GPUSync();
			GPUMultiplyNum(dfft,1.0/size,sizeb);
			GPUSync();
			GPUMemD2H((void *)htmp,(void *)dfft,sizeof(float)*sizeb);
			fft2buf(bufmrc,htmp,nsam);
			stackCorr.write2DIm(bufmrc,j-nStartSum);
		}
	}

//output the second sum image (less frames).
	if (para.nEndSum2>0 && para.nEndSum2<para.nEndSum)
	{
		GPUMemZero((void **)&dtmp,sizeof(float)*sizeb);
		GPUMemD2D(dtmp, dsumcorr, sizeof(float)*sizeb);
		GPUAdd(dtmp,dsum,sizeb);
		GPUSync();
		//copy sum image to host
		float *tsum=dtmp;
		GPUIFFT2d(tsum,ifft_plan);
		GPUMultiplyNum(tsum,1.0/size,sizeb);
		GPUMemD2H((void *)htmp,(void *)tsum,sizeof(float)*sizeb);
		fft2buf(bufmrc,htmp,nsam);
		//save
		char tempname[512]="";
		char subfilename[512]="";
		MRC mrc;
		if(para.bSaveSubAreaCorrSum) 
		{
			strncat(tempname,pThis->m_fnAlignsum2,strlen(pThis->m_fnAlignsum2)-4);
			sprintf(subfilename,"%s_ox%d_oy%d_dim%d.mrc",tempname,offsetx,offsety,nsamUnbin);
			mrc.open(subfilename,"wb");
		}
		else mrc.open(pThis->m_fnAlignsum2,"wb");
		mrc.createMRC(bufmrc,nsam,nsam,1);
		//stats
		sprintf(str,"Mean=%f   Min=%f   Max=%f\n",mrc.m_header.dmean,mrc.m_header.dmin,mrc.m_header.dmax);
		pThis->TextOutput(str);
		mrc.close();
		if(para.bSaveSubAreaCorrSum) sprintf(str,"Save Sum to: %s\n",subfilename);
		else sprintf(str,"Save Sum to: %s\n",pThis->m_fnAlignsum2);
		pThis->TextOutput(str);
	}

//add the rest of the frames
	for(j=nEndSum2+1;j<=nEndSum;j++)
	{
		totalshift+=shift[j-1];
		
		//copy to GPU
		GPUMemH2D((void *)dfft,(void *)(hbuf+j*sizeb),sizeof(float)*sizeb);
		//shift
		GPUShift(dfft,dPosList,-totalshift.real()+midshift.real(),-totalshift.imag()+midshift.imag(), nsam);
		GPUSync();
		//Sum
		if(bFSCEven) GPUAdd(dsumcorr,dfft,sizeb);
		else GPUAdd(dsum,dfft,sizeb);
		bFSCEven=!bFSCEven;

		sprintf(str,"......Add Frame #%03d with xy shift: %8.4lf %8.4lf\n",j+para.nStart,-totalshift.real()+midshift.real(),-totalshift.imag()+midshift.imag());
		pThis->TextOutput(str);

		//save stack
		if(para.bSaveStackCorr)
		{
			GPUIFFT2d(dfft,ifft_plan);
			GPUSync();
			GPUMultiplyNum(dfft,1.0/size,sizeb);
			GPUSync();
			GPUMemD2H((void *)htmp,(void *)dfft,sizeof(float)*sizeb);
			fft2buf(bufmrc,htmp,nsam);
			stackCorr.write2DIm(bufmrc,j-nStartSum);
		}
	}

	//close save stack
	if(para.bSaveStackCorr) stackCorr.close();
	//final sum
	GPUMemD2H((void *)hFSCCorr0,(void *)dsumcorr,sizeof(float)*sizeb);
	GPUMemD2H((void *)hFSCCorr1,(void *)dsum,sizeof(float)*sizeb);
	GPUAdd(dsumcorr,dsum,sizeb);
	GPUSync();

	//Make fft modulus for display
	if(para.bDispFFTCorr)
	{
		GPUSync();
		GPUFFTLogModulus(dfft, dsumcorr, nsam, para.fftscale);
		GPUSync();
		GPUMemD2H((void *)hdisp,(void *)dfft,sizeof(float)*(nsam/2+1)*nsam);
		pThis->FFTOutputCorr(hdisp);   //has been move to below
	}

	//copy binned sum to display
	if(para.bDispSumCorr)
	{
		GPUMemBinD2H(hdisp, dsumcorr, DISPDIM, nsam);
		pThis->ImageOutput(hdisp);
	}

	//copy sum image to host
	float *tsum=dsumcorr;
	GPUIFFT2d(tsum,ifft_plan);
	GPUMultiplyNum(tsum,1.0/size,sizeb);
	GPUMemD2H((void *)htmp,(void *)tsum,sizeof(float)*sizeb);
	fft2buf(bufmrc,htmp,nsam);
	
	//save
	char tempname[512]="";
	char subfilename[512]="";
	MRC mrc;
	if(para.bSaveSubAreaCorrSum) 
	{
		strncat(tempname,pThis->m_fnAlignsum,strlen(pThis->m_fnAlignsum)-4);
		sprintf(subfilename,"%s_ox%d_oy%d_dim%d.mrc",tempname,offsetx,offsety,nsamUnbin);
		mrc.open(subfilename,"wb");
	}
	else mrc.open(pThis->m_fnAlignsum,"wb");
	mrc.createMRC(bufmrc,nsam,nsam,1);
	//stats
	sprintf(str,"Mean=%f   Min=%f   Max=%f\n",mrc.m_header.dmean,mrc.m_header.dmin,mrc.m_header.dmax);
	pThis->TextOutput(str);
	mrc.close();
	if(para.bSaveSubAreaCorrSum) sprintf(str,"Save Sum to: %s\n",subfilename);
	else sprintf(str,"Save Sum to: %s\n",pThis->m_fnAlignsum);
	pThis->TextOutput(str);

	if(para.bLogFSC)
	{
		pThis->PlotFSC((cuComplex *)hFSCRaw0, (cuComplex *)hFSCRaw1, (cuComplex *)hFSCCorr0, 
					(cuComplex *)hFSCCorr1,hPosList,nsam,totalshift);
	}

	sprintf(str,"Done.\n");
	pThis->TextOutput(str);

	delete [] bufmrc;
	delete [] hbuf;
	delete [] bufmrc2;
//	delete [] bufscale;
//	delete [] bufscale2;
	delete [] hbuf2;
	delete [] hPosList;
	delete [] bufdark;
	delete [] bufnorm;
	delete [] bufdark2;
	delete [] bufnorm2;

	GPUMemFree((void **)&dPosList);
	GPUMemFree((void **)&dsum);
	GPUMemFree((void **)&dsumcorr);
	GPUMemFree((void **)&dfft);
	GPUMemFree((void **)&dtmp);
	//GPUFFTDestroy(fft_plan);
	GPUFFTDestroy(ifft_plan);

	delete [] htmp;
	delete [] hboxmap;
	delete [] hdisp;
	delete [] hFSCRaw0;
	delete [] hFSCRaw1;
	delete [] hFSCCorr0;
	delete [] hFSCCorr1;

	ResetGPU();
	pThis->Done();

	return (void *)0;
}
