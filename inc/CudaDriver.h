#ifndef CUDA_DRIVER_H
#define CUDA_DRIVER_H

#include "GenomeBwt.h"

/*
 * The CudaDriver class is used to run all of the CUDA code, and to be an interface between the
 * GPU implementation and the rest of GNUMAP.
 */

class CudaDriver {
	public:
		/*
		 * init_device() will identify how many CUDA devices are available, and choose one. 
		 * @return true- there is at least one CUDA device available
		 * @return false- there are no available CUDA devices available
		 */
		static bool initDevice();

		/*
		 * print_cuda_error_GPU(const char *message) will print the last error message from CUDA
		 * to stderr.
		 * @input const char *message- the message to be printed
		 */
		static void printCudaErrorGPU(const char* message);

		/*
		 * getAvailableMem() will get the available memory of the CUDA device.
		 * @return the memory available on the CUDA device
		 */
		int getAvailableMem();

		/*
		 * copyBwtToDevice(bwt_t* bwt, const int availMem) will copy the BWT object from memory 
		 * to the CUDA device.
		 * @input bwt_t* bwt- a pointer to the BWT object
		 * @input int availMem- the available memory on the CUDA device
		 */
		void copyBwtToDevice(bwt_t* bwt, const int availMem);
	private:
		/*
		 * deviceNum is the CUDA device that will be utilized.
		 */
		int deviceNum;
};

#endif