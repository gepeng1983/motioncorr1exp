#include "DFAlign.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

APARA getDefaultPara()
{
	APARA p;
	p.crop_offsetx=0;  
	p.crop_offsety=0;  
	p.crop_nsam=0;   

	p.bin=1.0;

	p.xscale=1.0;
	p.yscale=1.0;

	p.nStart=0;  //first frame(0-base)
	p.nEnd=0;    //last frame(0-base)
	p.nStartSum=0;  //first frame to sum(0-base)
	p.nEndSum=0;    //last frame to sum(0-base)
	p.nEndSum2=0;    //last frame to sum(0-base)

	p.GPUNum=0;  // GPU device ID 

	p.bfactor=150;  // in pix^2
	p.CCPeakSearchDim=96;//search peak in this box
	p.FrameDistOffset=2;
	p.NoisePeakSize=0;
	p.kiThresh=1.0; //alignment error threshold in pixel

	p.bSaveSubAreaCorrSum=false;
	p.bSaveRawSum=false;
	p.bSaveStackCorr=false;
	p.bSaveCCmap=false;
	p.bSaveLog=true;

	p.bAlignToMid=1;

	

	//diplay para
	p.fftscale=0.0001;
	p.bDispSumCorr=true;
	p.bDispFFTRaw=true;
	p.bDispFFTCorr=true;
	p.bDispCCMap=true;
	p.bDispFSC=false;
	p.bLogFSC=false;


	//reserved parameters for Dialog window
	p.fscMax=0.25;
	
	return p;
}

bool VerifyPara(CDFAlign &align)
{
	//check number of frame
	int nFrame=align.getNFrame();
	if(nFrame<=1)
	{
		printf("Not a stack or wrong file. Abort\n");
		return false;
	}

	if(align.m_para.nStart<0)
	{
		align.m_para.nStart=0;
	}
	if(align.m_para.nEnd>=nFrame) 
	{
		align.m_para.nEnd=nFrame-1;
	}
	if(align.m_para.nEnd<=align.m_para.nStart)
	{
		align.m_para.nEnd=nFrame-1;
	}

	if(align.m_para.nStartSum<align.m_para.nStart)
	{
		align.m_para.nStartSum=align.m_para.nStart;
	}
	if(align.m_para.nEndSum>align.m_para.nEnd)
	{
		align.m_para.nEndSum=align.m_para.nEnd;
	}
	if(align.m_para.nEndSum<=align.m_para.nStartSum)
	{
		align.m_para.nEndSum=align.m_para.nEnd;
	}
	if(align.m_para.nEndSum2>align.m_para.nEnd)
	{
		align.m_para.nEndSum2=align.m_para.nEnd;
	}
	if(align.m_para.nEndSum2<=align.m_para.nStartSum)
	{
//		align.m_para.nEndSum2=align.m_para.nEnd;
		align.m_para.nEndSum2=0;
	}

	if(align.m_para.FrameDistOffset>=nFrame/2 || align.m_para.FrameDistOffset<0)
	{
		printf("Frame distance offset is too big or less than 0. Abort\n");
		return false;
	}

	if (align.m_para.nEndSum2 >= align.m_para.nEndSum)
	{
		printf("Last frame # for the second sum image cannot be larger than or equal to that for the first sum image.\n");
		return false;
	}
	
	return true;
}


bool getPara(int narg, char* argc[], CDFAlign &align)
{
	APARA &p=align.m_para;
	
	p=getDefaultPara();

	align.m_fnNorm[0]=0;
	align.m_fnDark[0]=0;

	align.m_bSaveDisp=true;
	strcpy(align.m_fnStack,argc[1]);
	
	if(narg%2==1)
	{
		printf("Wrong number of arguments. Abort.\n");
		return false;
	}
	
	//whether use default name
	bool bfnRawsum=true;
	bool bfnAlignsum=true;
	bool bfnAlignsum2=true;
	bool bfnStackCorr=true;
	bool bfnCCmap=true;
	bool bfnLog=true;
	

	int i;
	string option;
	for(i=2;i<narg;i++)
	{
		option=argc[i];
		if(option.compare("-crx")==0)
		{
			i++;
			p.crop_offsetx=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-cry")==0)
		{
			i++;
			p.crop_offsety=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-crd")==0)
		{
			i++;
			p.crop_nsam=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-bin")==0)
		{
			i++;
			p.bin=atof(argc[i]);
			if(p.bin>2.0001)
			{
				printf("Larger than maximum binning 2. Set to 2.\n");
				p.bin=2.0;
			}
			else if(p.bin<0.9999)
			{
				printf("Smaller than minimum binning 1. Set to 1.\n");
				p.bin=1.0;
			}

			continue;
		}
	
		if(option.compare("-scx")==0)
		{
			i++;
			p.xscale=atof(argc[i]);
			if (p.xscale<1.0)
			{
				printf("Stretching factor along X-axis needs to be >=1.0.\n");
				return false;
			}
			continue;
		}

		if(option.compare("-scy")==0)
		{
			i++;
			p.yscale=atof(argc[i]);
			if (p.yscale<1.0)
			{
				printf("Stretching factor along Y-axis needs to be >=1.0.\n");
				return false;
			}
			continue;
		}

		if(option.compare("-nst")==0)
		{
			i++;
			p.nStart=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-ned")==0)
		{
			i++;
			p.nEnd=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-nss")==0)
		{
			i++;
			p.nStartSum=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-nes")==0)
		{
			i++;
			p.nEndSum=atoi(argc[i]);
			continue;
		}

		if(option.compare("-2es")==0)
		{
			i++;
			p.nEndSum2=atoi(argc[i]);
			continue;
		}
	
		if(option.compare("-gpu")==0)
		{
			i++;
			p.GPUNum=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-bft")==0)
		{
			i++;
			p.bfactor=atof(argc[i]);
			continue;
		}
		
		if(option.compare("-pbx")==0)
		{
			i++;
			p.CCPeakSearchDim=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-fod")==0)
		{
			i++;
			p.FrameDistOffset=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-nps")==0)
		{
			i++;
			p.NoisePeakSize=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-kit")==0)
		{
			i++;
			p.kiThresh=atof(argc[i]);
			continue;
		}
		
		if(option.compare("-sub")==0)
		{
			i++;
			p.bSaveSubAreaCorrSum=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-srs")==0)
		{
			i++;
			p.bSaveRawSum=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-ssc")==0)
		{
			i++;
			p.bSaveStackCorr=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-scc")==0)
		{
			i++;
			p.bSaveCCmap=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-slg")==0)
		{
			i++;
			p.bSaveLog=atoi(argc[i]);
			continue;
		}
		if(option.compare("-atm")==0)
		{
			i++;
			p.bAlignToMid=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-dsp")==0)
		{
			i++;
			if(atoi(argc[i])==0)
			{
				p.bDispSumCorr=false;
				p.bDispFFTRaw=false;
				p.bDispFFTCorr=false;
				p.bDispCCMap=false;
				align.m_bSaveDisp=false;
			}
			continue;
		}
		
		if(option.compare("-fsc")==0)
		{
			i++;
			p.bLogFSC=atoi(argc[i]);
			continue;
		}
		
		if(option.compare("-frs")==0)
		{
			i++;
			strcpy(align.m_fnRawsum, argc[i]);
			bfnRawsum=false;
			continue;
		}
		
		if(option.compare("-fcs")==0)
		{
			i++;
			strcpy(align.m_fnAlignsum, argc[i]);
			bfnAlignsum=false;
			continue;
		}
			
		if(option.compare("-2cs")==0)
		{
			i++;
			strcpy(align.m_fnAlignsum2, argc[i]);
			bfnAlignsum2=false;
			continue;
		}
	
		if(option.compare("-fct")==0)
		{
			i++;
			strcpy(align.m_fnStackCorr, argc[i]);
			bfnStackCorr=false;
			continue;
		}
		
		if(option.compare("-fcm")==0)
		{
			i++;
			strcpy(align.m_fnCCmap, argc[i]);
			bfnCCmap=false;
			continue;
		}
		
		if(option.compare("-flg")==0)
		{
			i++;
			strcpy(align.m_fnLog, argc[i]);
			bfnLog=false;
			continue;
		}
		if(option.compare("-fnm")==0)
		{
			i++;
			strcpy(align.m_fnNorm, argc[i]);
			continue;
		}
		if(option.compare("-fdk")==0)
		{
			i++;
			strcpy(align.m_fnDark, argc[i]);
			continue;
		}

		if(option.compare("-itr")==0)
		{
			i++;
			align.m_iter=atoi(argc[i]);
			continue;
		}
	
		if(option.compare("-ith")==0)
		{
			i++;
			align.m_iterthres=atof(argc[i]);
			continue;
		}
	
		printf("Undefined option: %s .Abort.\n",argc[i]);
		return false;
	}
	
	//Default file name
	string fnp=argc[1];
	if(fnp.rfind(".")!=string::npos) fnp=fnp.substr(0,fnp.rfind("."));
//	if(p.bin==2) fnp+="_2x";
	if(p.bin>1.0001) fnp+="_binned";
	string path,prefix;
	if(fnp.rfind("/")==string::npos)
	{
		path=".";
		prefix=fnp;
	}
	else
	{
//		path=fnp.substr(0,fnp.rfind("/"));
		path=".";
		prefix=fnp.substr(fnp.rfind("/")+1);
	}

	string fn;
	if(bfnRawsum)
	{
		fn=prefix+"_SumRaw.mrc";
		strcpy(align.m_fnRawsum,fn.c_str());
	}
	if(bfnAlignsum)
	{
		fn=prefix+"_SumCorr.mrc";
		strcpy(align.m_fnAlignsum,fn.c_str());
	}
	if(bfnAlignsum2)
	{
		fn=prefix+"_SumCorr2.mrc";
		strcpy(align.m_fnAlignsum2,fn.c_str());
	}
	if(bfnStackCorr)
	{
		fn=prefix+"_Corrected.mrc";
		strcpy(align.m_fnStackCorr,fn.c_str());
	}
	if(bfnCCmap)
	{
		fn=prefix+"_CCMap.mrc";
		strcpy(align.m_fnCCmap,fn.c_str());
	}
	if(bfnLog)
	{
		fn=prefix+"_Log.txt";
		strcpy(align.m_fnLog,fn.c_str());
	}
	
	struct stat st;
	if(align.m_bSaveDisp) 
	{
		//make a folder
		if(stat("dosef_quick",&st) != 0) system("mkdir dosef_quick");
		
		//add this folder into filename
//		string path,prefix;
//		if(fnp.rfind("/")==string::npos)
//		{
//			path=".";
//			prefix=fnp;
//		}
//		else
//		{
//			path=fnp.substr(0,fnp.rfind("/"));
//			prefix=fnp.substr(fnp.rfind("/")+1);
//		}
		string fndsp=path+"/dosef_quick/"+prefix;
		
		//finish filename
		fn=fndsp+"_RawFFT.mrc";
		strcpy(align.m_dispRawFFT,fn.c_str());
	
		fn=fndsp+"_CorrSum.mrc";
		strcpy(align.m_dispCorrSum,fn.c_str());
	
		fn=fndsp+"_CorrFFT.mrc";
		strcpy(align.m_dispCorrFFT,fn.c_str());
	}
	
	return VerifyPara(align);
}

void ShowPara(CDFAlign &align)
{
	APARA &p=align.m_para;
	
	if(!p.bSaveLog) return;
	FILE *fp=fopen(align.m_fnLog,"w");
	
	fprintf(fp,"******Parameter List******\n");
	fprintf(fp,"Input:  %s\n",align.m_fnStack);
	fprintf(fp,"  -crx  %d\n",p.crop_offsetx);
	fprintf(fp,"  -cry  %d\n",p.crop_offsety);
	fprintf(fp,"  -crd  %d\n",p.crop_nsam);
	fprintf(fp,"  -bin  %f\n",p.bin);
	fprintf(fp,"  -scx  %f\n",p.xscale);
	fprintf(fp,"  -scy  %f\n",p.yscale);
	fprintf(fp,"  -nst  %d\n",p.nStart);		
	fprintf(fp,"  -ned  %d\n",p.nEnd);
	fprintf(fp,"  -nss  %d\n",p.nStartSum);
	fprintf(fp,"  -nes  %d\n",p.nEndSum);
	fprintf(fp,"  -2es  %d\n",p.nEndSum2);
	fprintf(fp,"  -gpu  %d\n",p.GPUNum);
	fprintf(fp,"  -bft  %f\n",p.bfactor);
	fprintf(fp,"  -pbx  %d\n",p.CCPeakSearchDim);
	fprintf(fp,"  -fod  %d\n",p.FrameDistOffset);
	fprintf(fp,"  -nps  %d\n",p.NoisePeakSize);	
	fprintf(fp,"  -kit  %f\n",p.kiThresh);	
	fprintf(fp,"  -sub  %d\n",p.bSaveSubAreaCorrSum);
									
	fprintf(fp,"  -srs  %d\n",p.bSaveRawSum);
	fprintf(fp,"  -ssc  %d\n",p.bSaveStackCorr);
	fprintf(fp,"  -scc  %d\n",p.bSaveCCmap);
	fprintf(fp,"  -slg  %d\n",p.bSaveLog);

	fprintf(fp,"  -atm  %d\n",p.bAlignToMid);
	fprintf(fp,"  -dsp  %d\n",align.m_bSaveDisp);
	fprintf(fp,"  -fsc  %d\n",p.bLogFSC);
	
	fprintf(fp,"  -fcs  %s\n",align.m_fnAlignsum);
	fprintf(fp,"  -2cs  %s\n",align.m_fnAlignsum2);

	if(p.bSaveRawSum)    fprintf(fp,"  -frs  %s\n",align.m_fnRawsum);	
	if(p.bSaveStackCorr) fprintf(fp,"  -fct  %s\n",align.m_fnStackCorr);
	if(p.bSaveCCmap)     fprintf(fp,"  -fcm  %s\n",align.m_fnCCmap);
	if(p.bSaveLog)       fprintf(fp,"  -flg  %s\n",align.m_fnLog);
	
	if(align.m_bSaveDisp)
	{
		fprintf(fp,"\nQuick results in dosef_quick folder:\n");
		fprintf(fp,"UnCorrected FFT: %s \n",align.m_dispRawFFT);
		fprintf(fp,"Corrected Sum:   %s \n",align.m_dispCorrSum);
		fprintf(fp,"Corrected FFT:   %s \n",align.m_dispCorrFFT);
	}
	if (strlen(align.m_fnDark)) fprintf(fp,"  -fdk  %s\n",align.m_fnDark);	
	if (strlen(align.m_fnNorm)) fprintf(fp,"  -fnm  %s\n",align.m_fnNorm);	

	fprintf(fp,"  -itr  %d\n",align.m_iter);
	fprintf(fp,"  -ith  %f\n",align.m_iterthres);

	fprintf(fp,"******Parameter End******\n\n");
	
	fclose(fp);

}

//Input: 1InputStack.mrc [OPTION VALUE]
int main(int narg, char* argc[])
{
	if(narg==1) 
	{
		printf("Dose Fractionation Tool:\n");
		printf("Drift correction\n\n");
		printf("    Input: InputStack.mrc [OPTION VALUE] ...\n");
		printf("          *Note: If not specify OPTION, the default value will be used.\n\n");
		printf("           OPTION     VALUE(Default)    Introduction  \n");
		printf("           -crx       0                 Image crop offset X\n");
		printf("           -cry       0                 Image crop offset Y\n");
		printf("           -crd       0                 Image crop dimension. If 0, use maximum size.\n");
		printf("           -bin       1.0               Float number between 1x and 2x. Bin stack before processing.\n");
		printf("           -scx       1.0               Float number (>=1.0). Stretch image along X axis.\n");
		printf("           -scy       1.0               Float number (>=1.0). Stretch image along Y axis.\n");

		printf("           -nst       0                 First frame (0-base) used in alignment.\n");		
		printf("           -ned       0                 Last frame (0-base) used in alignment. If 0, use maximum value.\n");
		printf("           -nss       0                 First frame (0-base) used for final sum.\n");
		printf("           -nes       0                 Last frame (0-base) used for final sum. If 0,use maximum value.\n");
		printf("           -2es       0                 Last frame (0-base) used for the second sum.\n");

		printf("           -gpu       0                 GPU device ID.\n");
		printf("           -bft       150               BFactor in pix^2.\n");
		printf("           -pbx       96                Box dimension for searching CC peak.\n");
		printf("           -fod       2                 Number of frame offset for frame comparision.\n");
		printf("           -nps       0                 Radius of noise peak.\n");	
		printf("           -kit       1.0               Threshold of alignment error in pixel.\n");	
		printf("           -sub       0                 1: Save as sub-area corrected sum. 0: Not. \n");								
		printf("           -srs       0                 1: Save uncorrected sum. 0: Not.\n");
		printf("           -ssc       0                 1: Save aligned stack. 0: Not.\n");
		printf("           -scc       0                 1: Save CC Map. 0: Not.\n");
		printf("           -slg       1                 1: Save Log. 0: Not.\n");
		printf("           -atm       1                 1: Align to middle frame. 0: Not.\n");
		printf("           -dsp       1                 1: Save quick results. 0: Not.\n");
		printf("           -fsc       0                 1: Calculate and log FSC. 0: Not.\n");
		printf("           -frs       FileName.mrc      Uncorrected sum\n");
		printf("           -fcs       FileName.mrc      Corrected sum\n");
		printf("           -2cs       FileName.mrc      Corrected second sum\n");

		printf("           -fct       FileName.mrc      Corrected stack\n");
		printf("           -fcm       FileName.mrc      CC map\n");
		printf("           -flg       FileName.txt      Log file\n");

		printf("           -fdk       FileName.mrc      Dark reference\n");
		printf("           -fnm       FileName.mrc      Gain reference\n");

		printf("           -itr       0                 Number of iteration for iterative alignment\n");
		printf("           -ith       0.0               Convergence threshold for iterative alignment (RMSD shift)\n");
		
		printf("\n  Wrote by Xueming Li @ Yifan Cheng Lab, UCSF\n");
		printf("\n  Modified by Jiansen Jiang @ Z. Hong Zhou Lab, UCLA\n");

		return 0;
	}
	
	CDFAlign align;
	if(!getPara(narg,argc,align))
	{
		printf("Error: Have wrong parameters. Abort!\n");
		return 1;
	}
	
	ShowPara(align);

	align.RunAlign();
	
	return 0;
}
