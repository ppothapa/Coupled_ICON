#pragma once

#ifdef __HIP__
#include <iostream>
#include <hip/hip_runtime.h>
#endif

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C"
{
#endif
#ifndef __HIP__

	void c_generate_index_list_i1(
			const char* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int& nvalid, cudaStream_t stream);

	void c_generate_index_list_i4(
			const int* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int& nvalid, cudaStream_t stream);

	void c_generate_index_list_batched_i1(
			const int batch_size,
			const char* dev_conditions, const int stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, cudaStream_t stream);

	void c_generate_index_list_batched_i4(
			const int batch_size,
			const int* dev_conditions, const int stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, cudaStream_t stream);
#else

// HIP implementation

	void c_generate_index_list_i1(
			const char* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int& nvalid, hipStream_t stream);

	void c_generate_index_list_i4(
			const int* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int& nvalid, hipStream_t stream);

	void c_generate_index_list_batched_i1(
			const int batch_size,
			const char* dev_conditions, const int stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, hipStream_t stream);

	void c_generate_index_list_batched_i4(
			const int batch_size,
			const int* dev_conditions, const int stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, hipStream_t stream);

#endif


#ifdef __cplusplus
}
#endif
