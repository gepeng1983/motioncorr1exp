
CXXFLAGS = -static -O3 -I./SP++3/include
CUFLAG = -O3 -lcuda -lcufft -arch=sm_60  -Xptxas -dlcm=ca  -lpthread -I./SP++3/include -I/usr/local/cuda-5.5/include
CULIBFLAG = -lfftw3f -lm -ltiff -lz
NVCC =  /usr/local/cuda/bin/nvcc
CXX = g++

all: DoseFProject

dosefobjs = func.o mrcz.o mrct.o mrc.o cufunc.o DFAlign.o safefft.o 
dosefbins = gpuinfo dosefgpu_driftcorr


DoseFProject: ${dosefbins}

gpuinfo: gpuinfo.cu
	${NVCC} -O3 gpuinfo.cu -o gpuinfo	
	cp -f gpuinfo ../bin

dosefgpu_driftcorr: ${dosefobjs}
	${NVCC} ${CUFLAG} ${CULIBFLAG} ${dosefobjs} dosefgpu_driftcorr.cpp -o dosefgpu_driftcorr
	cp -f dosefgpu_driftcorr ../bin


cufunc.o: cufunc.cu cufunc.h

DFAlign.o: DFAlign.cu DFAlign.h mrc.h mrcz.h mrct.h cufunc.h mrcheader.h

func.o: func.cpp func.h mrc.h mrcz.h mrct.h mrcheader.h

mrc.o: mrc.cpp mrc.h mrcheader.h

mrcz.o: mrcz.cpp mrc.h mrcz.h mrcheader.h

mrct.o: mrct.cpp mrct.h mrc.h mrcheader.h

savefft.o: savefft.cpp savefft.h






%.o: %.cu
	${NVCC} ${CUFLAG} -c $< -o $@

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	rm -rf *.o
	rm -rf ${dosefbins}
