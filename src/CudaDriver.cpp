#include "CudaDriver.h"
#include <cuda_runtime_api.h>
#include <cuda.h>

/*
 * Constructor. Creates new Cuda class.
*/ 
CudaDriver::CudaDriver() {
	deviceNum = 0;	
	availMem = 0;
	totalMem = 0;
}

CudaDriver::~CudaDriver() {

}
 
bool CudaDriver::initDevice() {
	int numDevices = 0;
	cudaGetDeviceCount(&numDevices);
	printCudaErrorGPU("Getting device count.");

	if(!numDevices) {
		fprintf(stderr, "There are no CUDA devices found.\n");
		return false;
	}
	else {
		if(numDevices == 1) {
			fprintf(stderr, "There is %d Cuda device found.\n", numDevices);
		}
		else {
			fprintf(stderr, "There are %d CUDA devices found.\n", numDevices);
		}
		return true;
	}
}

void CudaDriver::printCudaErrorGPU(const char* message) {
	cudaError_t cudaErr = cudaGetLastError();

	if(cudaSuccess != cudaErr) {
		fprintf(stderr, "CUDA Error: %s\n", message);
		fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(cudaErr));
	}
}

void CudaDriver::getAvailableMem() {
	cudaMemGetInfo(&availMem, &totalMem);
	printCudaErrorGPU("Checking the available memory.");
	fprintf(stderr, "CUDA device has %lu available bytes.\n", availMem);
	
}

void CudaDriver::copyBwtToDevice(bwt_t* bwt, const int availMem) {

}
