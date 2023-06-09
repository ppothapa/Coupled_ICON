#pragma once

#ifdef __HIP__
#include <iostream>
#include <hip/hip_runtime.h>
using gpuStream_t = hipStream_t;
#else
// CUDA
using gpuStream_t = cudaStream_t;
#endif

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C"
{
#endif

	void c_generate_index_list_gpu_single(
			const void* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int* nvalid, 
			int data_size, bool copy_to_host,
			gpuStream_t stream);

	void c_generate_index_list_gpu_batched(
			const int batch_size,
			const void* dev_conditions, const int stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, int data_size,
			gpuStream_t stream);

#ifdef __cplusplus
}
#endif
