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

size_t getSizeOfbwt_t(bwt_t* bwt) {
	size_t size = 0;
	// the type bwtint_t is equivalent to a uint64_t
	size += sizeof(bwt->primary);
	size += sizeof(bwt->L2);
	size += sizeof(bwt->seq_len);
	size += sizeof(bwt->bwt_size);
	size += (size_t) (sizeof(bwt->bwt) * bwt->bwt_size); 
	size += sizeof(bwt->cnt_table);
	size += sizeof(bwt->sa_intv);
	size += sizeof(bwt->n_sa);
	size += (size_t) (sizeof(bwt->sa) * bwt->n_sa); // TODO check that bwt->n_sa is actually the length of the array bwt->sa
//	the members of struct bwt_t
//	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
//	bwtint_t L2[5]; // C(), cumulative count
//	bwtint_t seq_len; // sequence length
//	bwtint_t bwt_size; // size of bwt, about seq_len/4
//	uint32_t *bwt; // BWT
//	// occurance array, separated to two parts
//	uint32_t cnt_table[256];
//	// suffix array
//	int sa_intv;
//	bwtint_t n_sa;
//	bwtint_t *sa;
	return size;
}

size_t getSizeOfbntann1_t(bntann1_t* anns) {
	size_t size = 0;
	size += sizeof(anns->offset);
	size += sizeof(anns->len);
	size += sizeof(anns->n_ambs);
	size += sizeof(anns->gi);
	size += sizeof(anns->is_alt);
	size += sizeof(anns->name);
	size += sizeof(anns->anno);
//	the members of struct bntann1_t
//	int64_t offset;
//	int32_t len;
//	int32_t n_ambs;
//	uint32_t gi;
//	int32_t is_alt;
//	char *name, *anno;
	return size;
}

size_t getSizeOfbntamb1_t(bntamb1_t* ambs) {
	size_t size = 0;
	size += sizeof(ambs->offset);
	size += sizeof(ambs->len);
	size += sizeof(ambs->amb);
//	the members of struct bntamb1_t
//	int64_t offset;
//	int32_t len;
//	char amb;
	return size;
}

size_t getSizeOfFILE(FILE* fp_pac) {
	size_t size = 0;
	// TODO implement the sizeof type FILE
	return size;
}

size_t getSizeOfbntseq_t(bntseq_t* bns) {
	size_t size = 0;
	size += sizeof(bns->l_pac);
	size += sizeof(bns->n_seqs);
	size += sizeof(bns->seed);
	size += getSizeOfbntann1_t(bns->anns); 
	size += sizeof(bns->n_holes);
	size += getSizeOfbntamb1_t(bns->ambs);
	size += getSizeOfFILE(bns->fp_pac); 
//	the members of struct bntseq_t
//	int64_t l_pac;
//	int32_t n_seqs;
//	uint32_t seed;
//	bntann1_t *anns; // n_seqs elements
//	int32_t n_holes;
//	bntamb1_t *ambs; // n_holes elements
//	FILE *fp_pac;
	return size;
}

size_t getSizeOfSeq(uint8_t* pac) {
	size_t size = 0;
	// TODO figure out how to calculate the size of this memory
	return size;
}

size_t getSizeOfMem(uint8_t* mem) {
	size_t size = 0;
	// TODO figure out how to calculate the size ofthis memory
	return size;
}

size_t getSizeOfbwaidx_t(bwaidx_t* bwtIndex) {
	size_t size = 0;
	size += getSizeOfbwt_t(bwtIndex->bwt);
	size += getSizeOfbntseq_t(bwtIndex->bns);
	size += getSizeOfSeq(bwtIndex->pac);
	size += sizeof(bwtIndex->is_shm);
	size += sizeof(bwtIndex->l_mem);
	size += getSizeOfMem(bwtIndex->mem); 
	// the members of struct bwaidx_t
	// bwt_t    *bwt; // FM-index
	// bntseq_t *bns; // information on the reference sequences
	// uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
	// int    is_shm;
	// int64_t l_mem;
	// uint8_t  *mem; 
	return size;
}

bool CudaDriver::copyBwtToDevice(bwaidx_t* index) {
	if(getSizeOfbwaidx_t(index) > availMem) {
		return false;
	}

	return true;
}
