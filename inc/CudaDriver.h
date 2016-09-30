#ifndef CUDA_DRIVER_H
#define CUDA_DRIVER_H

#include "GenomeBwt.h"

/*
 * The CudaDriver class is used to run all of the CUDA code, and to be an interface between the
 * GPU implementation and the rest of GNUMAP.
 */

class CudaDriver {
	public:
		CudaDriver();
		~CudaDriver();
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
		 */
		void getAvailableMem();

		/*
		 * copyBwtToDevice(bwaidx_t* bwt) will copy the BWT object from memory 
		 * to the CUDA device.
		 * @input bwaidx_t* index-- a pointer to the BWT object
		 * @return true if the index is copied, false if the index is not copied
		 */
		bool copyBwtToDevice(bwaidx_t* index);
	private:
		/*
		 * deviceNum is the CUDA device that will be utilized.
		 */
		int deviceNum;
		/*
		 * availMem is the memory available on the CUDA device.
		 */
		size_t availMem;
		/*
		 * totalMem is the total memory of the CUDA device.
		 */
		size_t totalMem;
};

#endif
