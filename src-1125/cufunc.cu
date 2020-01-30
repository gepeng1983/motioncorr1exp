#include "cufunc.h"
#include <string.h>
#include "mrc.h"
#include "func.h"
#include <signal.h> 
#define PI 3.141592653589793
#define BLOCKSIZE 1024


static __device__ cuComplex conj(cuComplex f)
{
	f.y*=-1.0;
	return f;
}


static __device__ cuComplex cXc(cuComplex a, cuComplex b) // a*b
{
	cuComplex c;
	c.x=a.x*b.x-a.y*b.y;
	c.y=a.x*b.y+a.y*b.x;
	return c;
}

static __device__ float cabs(cuComplex a) // a*b
{
	return sqrt(a.x*a.x+a.y*a.y);
}

bool initGPU(int GPUNum)
{
	//initional CUDA device
	int ngpu;
	cudaGetDeviceCount(&ngpu);
	if(ngpu <= 0)
	{
		return false;
	}
	if(GPUNum>=ngpu)
	{
		printf("GPU ID %d is out of range(%d). Abort.\n",GPUNum,ngpu);
		return false;
	}

	cudaDeviceProp prop;
	if(cudaGetDeviceProperties(&prop, GPUNum) == cudaSuccess) 
	{
		printf("Use GPU: #%d %s\n",GPUNum,prop.name);
		if(prop.kernelExecTimeoutEnabled)
		{
			printf("Warnning: This GPU is also used for display, may not stable.\n");
		}
	}


	if(cudaSetDevice(GPUNum)!=cudaSuccess)
	{
		printf("Error: Failed to set CUDA Device #%d. Abort.\n",GPUNum);
		return false;
	}
	
	signal(SIGINT, siginthandler);

	return true;
}

bool ResetGPU()
{
	if(cudaDeviceReset()!=cudaSuccess) return false;
	return true;
}

void siginthandler(int param) 
{   
	if(ResetGPU()) printf(" GPU was reset successfully after process was killed.\n"); 
	else printf(" Error: Failed to reset GPU.\n"); 
	exit(1); 
}

int getGPUList(vector<string> &namelist)
{
	int ngpu;
	cudaGetDeviceCount(&ngpu);
	if(ngpu <= 0)
	{
		return 0;
	}

	namelist.clear();
	int i;
	string str;
	for(i = 0; i < ngpu; i++) 
	{
		cudaDeviceProp prop;
		if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) 
		{
			str=prop.name;
			namelist.push_back(str);
		}
	}
	
	return ngpu;
}

void GPUMemCheck(size_t &theFree, size_t &theTotal)
{
	cuMemGetInfo( &theFree, &theTotal );  
	//printf( "CARD returns:  free:%d  total:%d\n", theFree, theTotal);
}

bool GPUMemAlloc(void **buf, int size)
{
	if(cudaMalloc((void **)buf,size)!=cudaSuccess) return false;
	return true;
}

bool GPUMemZero(void **buf, int size)
{
	if(cudaMemset(*buf,0,size)!=cudaSuccess) return false;
	return true;
}

bool GPUMemFree(void **buf)
{
	if(cudaFree(*buf)!=cudaSuccess) return false;
	*buf=0;
	return true;
}

bool GPUMemH2D(void *dst, void *src, int size)
{
	if(cudaMemcpy(dst,src,size,cudaMemcpyHostToDevice)!=cudaSuccess) return false;
	return true;
}

bool GPUMemD2H(void *dst, void *src, int size)
{
	if(cudaMemcpy(dst,src,size,cudaMemcpyDeviceToHost)!=cudaSuccess) return false;
	return true;
}
bool GPUMemD2D(void *dst, void *src, int size)
{
	if(cudaMemcpy(dst,src,size,cudaMemcpyDeviceToDevice)!=cudaSuccess) return false;
	return true;
}

bool GPUMemBinD2H(float *dst, float *src, int dst_nsam, int src_nsam)
{
	int i;
	int size=sizeof(float)*(dst_nsam+2);
	for(i=0;i<dst_nsam/2;i++)
	{
		//up half
		if(cudaMemcpy(dst+i*(dst_nsam+2),src+i*(src_nsam+2),size,cudaMemcpyDeviceToHost)!=cudaSuccess) return false;

		//down half
		if(cudaMemcpy(dst+(dst_nsam-1-i)*(dst_nsam+2),src+(src_nsam-1-i)*(src_nsam+2),size,cudaMemcpyDeviceToHost)!=cudaSuccess) return false;
	}

	return true;
}
bool GPUMemBinD2D(float *dst, float *src, int dst_nsam, int src_nsam)
{
	int i;
	int size=sizeof(float)*(dst_nsam+2);
	for(i=0;i<dst_nsam/2;i++)
	{
		//up half
		if(cudaMemcpy(dst+i*(dst_nsam+2),src+i*(src_nsam+2),size,cudaMemcpyDeviceToDevice)!=cudaSuccess) return false;

		//down half
		if(cudaMemcpy(dst+(dst_nsam-1-i)*(dst_nsam+2),src+(src_nsam-1-i)*(src_nsam+2),size,cudaMemcpyDeviceToDevice)!=cudaSuccess) return false;
	}

	return true;
}

__global__ void cuFFTLogModulus(float *dMod, cuComplex *dfft, int size, float scale)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=size) return;
	dMod[id]=log(1+cabs(dfft[id])*scale);
}
void GPUFFTLogModulus(float *dMod, float *dfft, int nsam, float scale)
{
	int size=(nsam/2+1)*nsam;
	cuFFTLogModulus<<<size/BLOCKSIZE+1,BLOCKSIZE>>>(dMod,(cuComplex *)dfft,size,scale);
}

__global__ void cuFFTModulus(float *dMod, cuComplex *dfft, int size)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=size) return;
	dMod[id]=cabs(dfft[id]);
}
void GPUFFTModulus(float *dMod, float *dfft, int nsam)
{
	int size=(nsam/2+1)*nsam;
	cuFFTModulus<<<size/BLOCKSIZE+1,BLOCKSIZE>>>(dMod,(cuComplex *)dfft,size);
}

cufftHandle GPUFFTPlan(int nsam)
{
	cufftHandle plan;
	cufftResult r=cufftPlan2d(&plan,nsam,nsam,CUFFT_R2C);
	//cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_PADDING);
	return plan;
}
cufftHandle GPUIFFTPlan(int nsam)
{
	cufftHandle plan;
	cufftResult r=cufftPlan2d(&plan,nsam,nsam,CUFFT_C2R);
	//cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_FFTW_PADDING);
	return plan;
}

void GPUFFTDestroy(cufftHandle &plan)
{
	cufftDestroy(plan);
	plan=0;
}

bool GPUFFT2d(float* dfft, cufftHandle plan)
{
	if(cufftExecR2C(plan,(cufftReal*)dfft,(cufftComplex *)dfft)!=CUFFT_SUCCESS) return false;
	return true;
}

bool GPUIFFT2d(float* dfft, cufftHandle plan)
{
	if(cufftExecC2R(plan,(cufftComplex *)dfft,(cufftReal*)dfft)!=CUFFT_SUCCESS) return false;
	return true;
}

bool GPUSync()
{
	if(cudaThreadSynchronize()!=cudaSuccess) return false;
	if(cudaGetLastError()!=cudaSuccess) return false;
	return true;
}


__global__ void cuAdd(float *dst, float *src, int size)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=size) return;
	dst[id]+=src[id];
}
void GPUAdd(float *dst, float *src, int size)
{
	cuAdd<<<size/BLOCKSIZE+1,BLOCKSIZE>>>(dst,src,size);
}

__global__ void cuMultiplyNum(float *dst, float num, int size)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=size) return;
	dst[id]*=num;
}
void GPUMultiplyNum(float *dst, float num, int size)
{
	cuMultiplyNum<<<size/BLOCKSIZE+1,BLOCKSIZE>>>(dst,num,size);
}

//int3[x,y,sign]
void MkPosList(int3 *list, int nsam, float inner_r, float outer_r)
{
	int hnsamb=nsam/2+1;
	int hnsam=nsam/2;
	int i,j;
	int count=0;
	int r2;
	int ri2=int(inner_r*inner_r);
	int ro2=int(outer_r*outer_r);
	for(j=0;j<hnsam;j++)
		for(i=0;i<hnsamb;i++)
		{
			list[count].x=i;
			list[count].y=j;

			r2=list[count].x*list[count].x+list[count].y*list[count].y;
			if(r2<ri2 || r2>ro2) list[count].z=0;
			else 
			{
				if((list[count].x+list[count].y)%2==0) list[count].z=1;
				else list[count].z=-1;
			}

			count++;
		}
	for(j=hnsam;j<nsam;j++)
		for(i=0;i<hnsamb;i++)
		{
			list[count].x=i;
			list[count].y=j-nsam;
			r2=list[count].x*list[count].x+list[count].y*list[count].y;
			if(r2<ri2 || r2>ro2) list[count].z=0;
			else 
			{
				if((list[count].x+list[count].y)%2==0) list[count].z=1;
				else list[count].z=-1;
			}
			count++;
		}

}

//MASK[x,y,sign*bfactor]
void MkPosList(MASK *list, int nsam, float bfactor)
{
	int hnsamb=nsam/2+1;
	int hnsam=nsam/2;
	int i,j;
	int count=0;
	int r2;
	float m=-0.5*bfactor/nsam/nsam;
	for(j=0;j<hnsam;j++)
		for(i=0;i<hnsamb;i++)
		{
			list[count].x=i;
			list[count].y=j;

			r2=list[count].x*list[count].x+list[count].y*list[count].y;
			if((list[count].x+list[count].y)%2==0) list[count].z=exp(m*r2);
			else list[count].z=-exp(m*r2);

			count++;
		}
	for(j=hnsam;j<nsam;j++)
		for(i=0;i<hnsamb;i++)
		{
			list[count].x=i;
			list[count].y=j-nsam;
			r2=list[count].x*list[count].x+list[count].y*list[count].y;
			if((list[count].x+list[count].y)%2==0) list[count].z=exp(m*r2);
			else list[count].z=-exp(m*r2);

			count++;
		}

}

__global__ void cuShiftCC(cuComplex *dfft, cuComplex *dsum, MASK *dposlist,float shx, float shy, int size)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=size) return;

	MASK pos=dposlist[id];
	float shift=shx*pos.x+shy*pos.y;
	cuComplex phase;
	phase.x=cos(shift);
	phase.y=sin(shift);

	cuComplex val=cXc(dsum[id],conj(cXc(dfft[id],phase)));
	val.x/=size;
	val.x*=pos.z;
	val.y/=size;
	val.y*=pos.z;
	dfft[id]=val;

}
void GPUShiftCC(float *dfft, float *dsum, MASK *dposlist, float sx, float sy, int nsam)
{
	float shx=sx*2.0*PI/nsam;
	float shy=sy*2.0*PI/nsam;

	int size=(nsam/2+1)*nsam;
	cuShiftCC<<<size/BLOCKSIZE+1,BLOCKSIZE>>>((cuComplex *)dfft, (cuComplex *)dsum, dposlist,shx, shy, size);

}

__global__ void cuShift(cuComplex *dfft,MASK *dposlist,float shx, float shy, int size)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=size) return;

	MASK pos=dposlist[id];
	float shift=shx*pos.x+shy*pos.y;
	cuComplex phase;
	phase.x=cos(shift);
	phase.y=sin(shift);

	dfft[id]=cXc(dfft[id],phase);

}
void GPUShift(float *dfft, MASK *dposlist, float sx, float sy, int nsam)
{
	float shx=sx*2.0*PI/nsam;
	float shy=sy*2.0*PI/nsam;

	int size=(nsam/2+1)*nsam;
	cuShift<<<size/BLOCKSIZE+1,BLOCKSIZE>>>((cuComplex *)dfft,dposlist,shx, shy, size);

}

float FindShift(float *dsrc,int nsam, float* hboxmap, int box, float &sx, float &sy, int wNoise)
{
	int ori=(nsam-box)/2;
	int nsamb=nsam+2;
	int i,j,id,is;
	float bestcc=-1e9;

	//float *dst=new float[box*box];
	for(i=0;i<box;i++)
	{
		cudaMemcpy(hboxmap+i*box,dsrc+(i+ori)*nsamb+ori,sizeof(float)*box,cudaMemcpyDeviceToHost);
	}

	sx=0;
	sy=0;
	for(j=0;j<box;j++)
		for(i=0;i<box;i++)
		{
			id=j*box+i;

			if(abs(i-box/2)<=wNoise && abs(j-box/2)<=wNoise) continue;

			if(hboxmap[id]>bestcc)
			{
				bestcc=hboxmap[id];
				sx=i;
				sy=j;
			}
		}


	//Fourier interpolation
	int subbox=16; //box/4;
	float *hsubboxmap=new float[(subbox+2)*subbox];
	int offsetx=sx-subbox/2;
	int offsety=sy-subbox/2;
	if((offsetx+subbox)>box) offsetx=box-subbox;
	else if(offsetx<0) offsetx=0;
	if((offsety+subbox)>box) offsety=box-subbox;
	else if(offsety<0) offsety=0;
	//crop and fft
	crop2fft(hboxmap,box,hsubboxmap,offsetx,offsety,subbox);
	fft2d(hsubboxmap,subbox);
	cosmask2d((complex<float> *)hsubboxmap,subbox);
	//pad
	int scale=32;
	int wNoiseScaled=scale*wNoise;
	int pad=subbox*scale;
	float *hpadmap=new float[(pad+2)*pad];
	memset(hpadmap,0,sizeof(float)*(pad+2)*pad);
	for(i=0;i<subbox/2;i++)
	{
		memcpy(hpadmap+i*(pad+2),hsubboxmap+i*(subbox+2),sizeof(float)*(subbox+2));
	}
	for(i=0;i<subbox/2;i++)
	{
		memcpy(hpadmap+(pad-1-i)*(pad+2),hsubboxmap+(subbox-1-i)*(subbox+2),sizeof(float)*(subbox+2));
	}
	//ifft
	ifft2d(hpadmap,pad);
	//find shift
	int ox=(box/2-offsetx)*scale; //in order to avoid noise peak at box/2
	int oy=(box/2-offsety)*scale; //in order to avoid noise peak at box/2
	int sxp=0,syp=0;
	bestcc=hpadmap[0];
	for(j=0;j<pad;j++)
	{
		is=j*(pad+2);
		for(i=0;i<pad;i++)
		{
			if(abs(i-ox)<=wNoiseScaled && abs(j-oy)<=wNoiseScaled) continue;

			id=is+i;
			if(hpadmap[id]>bestcc)
			{
				bestcc=hpadmap[id];
				sxp=i;
				syp=j;
			}
		}
	}

	sx=offsetx+sxp/double(scale);
	sy=offsety+syp/double(scale);
	sx-=box/2;
	sy-=box/2;

	/*char filename[256];
	sprintf(filename,"D:\\UCSFImage\\DoseFragProcess\\data\\temp.mrc");
	MRC mrc;
	mrc.open(filename,"wb");
	mrc.createMRC(hpadmap,pad+2,pad,1);
	mrc.close();*/



	delete [] hsubboxmap;
	delete [] hpadmap;
	

	return bestcc/nsam/nsam/subbox/subbox;
}

float FindShift(float* hboxmap, int box, float &sx, float &sy)
{
	int i,j,id;
	float bestcc=-1e9;

	sx=0;
	sy=0;
	for(j=0;j<box;j++)
		for(i=0;i<box;i++)
		{
			id=j*box+i;
			if(hboxmap[id]>bestcc)
			{
				bestcc=hboxmap[id];
				sx=i;
				sy=j;
			}
		}
	sx-=box/2;
	sy-=box/2;

	return bestcc;
}

void testCUFFT()
{
	int i,j;
	int nsam=26;
	int size=(nsam+2)*nsam;
	float *h=new float[size];
	float *r=new float[size];
	memset(h,0,size*sizeof(float));
	memset(r,0,size*sizeof(float));
	for(j=0;j<nsam;j++)
		for(i=0;i<nsam;i++)
		{
			h[j*(nsam+2)+i]=i+1;
		}

	float *d=0;
	
	GPUMemAlloc((void **)&d,size*sizeof(float));
	GPUMemH2D(d,h,size*sizeof(float));

	cufftHandle fft_plan=GPUFFTPlan(nsam);
	cufftHandle ifft_plan=GPUIFFTPlan(nsam);

	GPUFFT2d(d,fft_plan);
	GPUSync();
	GPUIFFT2d(d,ifft_plan);
	GPUSync();

	GPUMultiplyNum(d,1.0/nsam/nsam,size);
	GPUMemD2H(r,d,sizeof(float)*size);

	char hstr[65536]="";
	char rstr[65536]="";
	char str[16]="";
	for(j=0;j<nsam;j++)
	{
		strcat(hstr,"\n");
		strcat(rstr,"\n");
		for(i=0;i<(nsam+2);i++)
		{
			sprintf(str,"%6.3f ",h[j*(nsam+2)+i]);
			strcat(hstr,str);
			sprintf(str,"%6.3f ",r[j*(nsam+2)+i]);
			strcat(rstr,str);
		}
		
	}
	
	GPUMemFree((void **)&d);
	GPUFFTDestroy(fft_plan);
	GPUFFTDestroy(ifft_plan);
	delete [] h;
	delete [] r;
	
/*GPUFFT2d(d,fft_plan);
	GPUSync();
	GPUIFFT2d(d,ifft_plan);
	GPUSync();*/
	return;
}
