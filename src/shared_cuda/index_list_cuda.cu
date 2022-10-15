#include "index_list_cuda.h"

#ifndef __HIP__
#include <cuda.h>
#include <cub/device/device_select.cuh>
#include <cub/iterator/counting_input_iterator.cuh>

template<typename T>
struct ZeroCmp
{
        const T* conditions;
        const int startid;

        ZeroCmp(const int startid, const T* conditions) :
                startid(startid), conditions(conditions)
        { }

        __device__ __host__ __forceinline__
        bool operator() (const int &id)
        {
          return (conditions[ id - startid ] != 0);
        }
};

template <typename T>
static
void c_generate_index_list_generic_device(
                        const T* dev_conditions,
                        const int startid, const int endid,
                        int* dev_indices,
                        int* dev_nvalid, cudaStream_t stream)
{
        static size_t storageSize = 0;
        static char* storage = nullptr;

        const int n = endid - startid + 1;

        // Argument is the offset of the first element
        cub::CountingInputIterator<int> iterator(startid);

        // Determine temporary device storage requirements
        size_t storageRequirement;
        cub::DeviceSelect::Flagged(nullptr, storageRequirement,
                        iterator, dev_conditions, dev_indices,
                        dev_nvalid, n, stream);

        // Allocate temporary storage (only if not enough)
        if (storageRequirement > storageSize)
        {
                cudaFree(storage);
                cudaMalloc(&storage, storageRequirement);
                storageSize = storageRequirement;
        }

        ZeroCmp<T> select(startid, dev_conditions);
        cub::DeviceSelect::If(storage, storageRequirement,
                        iterator, dev_indices,
                        dev_nvalid, n,
                        select, stream);
}


template <typename T>
static
void c_generate_index_list_batched_generic(
                        const int batch_size,
                        const T* dev_conditions, const int cond_stride,
                        const int startid, const int endid,
                        int* dev_indices, const int idx_stride,
                        int* dev_nvalid, cudaStream_t stream)
{
        for (int i = 0; i < batch_size; i++)
                c_generate_index_list_generic_device(
                                dev_conditions + cond_stride*i,
                                startid, endid,
                                dev_indices + idx_stride*i,
                                dev_nvalid + i, stream);
}

template <typename T>
static
void c_generate_index_list_generic(
                        const T* dev_conditions,
                        const int startid, const int endid,
                        int* dev_indices,
                        int& nvalid, cudaStream_t stream)
{
        static int* dev_nvalid = nullptr;
        if (dev_nvalid == nullptr)
                        cudaMalloc(&dev_nvalid, sizeof(int));

        c_generate_index_list_generic_device(
                        dev_conditions, startid, endid, dev_indices, dev_nvalid, stream);

        cudaMemcpyAsync(&nvalid, dev_nvalid, sizeof(int), cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
}

///
/// Exposed functions
/// 
/// Non-batched first
/// 

void c_generate_index_list_i1(
                        const char* dev_conditions,
                        const int startid, const int endid,
                        int* dev_indices,
                        int& nvalid, cudaStream_t stream)
{
        c_generate_index_list_generic(dev_conditions, startid, endid, dev_indices, nvalid, stream);
}

void c_generate_index_list_i4(
                        const int* dev_conditions,
                        const int startid, const int endid,
                        int* dev_indices,
                        int& nvalid, cudaStream_t stream)
{
        c_generate_index_list_generic(dev_conditions, startid, endid, dev_indices, nvalid, stream);
}

/// 
/// And now batched
/// 

void c_generate_index_list_batched_i1(
        const int batch_size,
        const char* dev_conditions, const int cond_stride,
        const int startid, const int endid,
        int* dev_indices, const int idx_stride,
        int* dev_nvalid, cudaStream_t stream)
{
c_generate_index_list_batched_generic(
                batch_size,
                dev_conditions, cond_stride,
                startid, endid,
                dev_indices, idx_stride,
                dev_nvalid, stream);
}

void c_generate_index_list_batched_i4(
                const int batch_size,
                const int* dev_conditions, const int cond_stride,
                const int startid, const int endid,
                int* dev_indices, const int idx_stride,
                int* dev_nvalid, cudaStream_t stream)
{
        c_generate_index_list_batched_generic(
                        batch_size,
                        dev_conditions, cond_stride,
                        startid, endid,
                        dev_indices, idx_stride,
                        dev_nvalid, stream);
}

#else

// HIP implementation

#include <hip/hip_runtime.h>
#include <hipcub/device/device_select.hpp>
#include <hipcub/iterator/counting_input_iterator.hpp>

template<typename T>
struct ZeroCmp
{
	const T* conditions;
	const int startid;

	ZeroCmp(const int startid, const T* conditions) :
		startid(startid), conditions(conditions)
	{ }

	__device__ __host__ __forceinline__
	bool operator() (const int &id)
	{
	  return (conditions[ id - startid ] != 0);
	}
};

template <typename T>
static
void c_generate_index_list_generic_device(
			const T* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int* dev_nvalid, hipStream_t stream)
{
	static size_t storageSize = 0;
	static char* storage = nullptr;

	const int n = endid - startid + 1;

	// Argument is the offset of the first element
	hipcub::CountingInputIterator<int> iterator(startid);

	// Determine temporary device storage requirements
	size_t storageRequirement;
	hipcub::DeviceSelect::Flagged(nullptr, storageRequirement,
			iterator, dev_conditions, dev_indices,
			dev_nvalid, n, 0);

	// Allocate temporary storage (only if not enough)
	if (storageRequirement > storageSize)
	{
		hipFree(storage);
		hipMalloc(&storage, storageRequirement);
		storageSize = storageRequirement;
	}

	ZeroCmp<T> select(startid, dev_conditions);
	hipcub::DeviceSelect::If(storage, storageRequirement,
			iterator, dev_indices,
			dev_nvalid, n,
			select, 0);
}


template <typename T>
static
void c_generate_index_list_batched_generic(
			const int batch_size,
			const T* dev_conditions, const int cond_stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, hipStream_t stream)
{
	for (int i = 0; i < batch_size; i++)
		c_generate_index_list_generic_device(
				dev_conditions + cond_stride*i,
				startid, endid,
				dev_indices + idx_stride*i,
				dev_nvalid + i, 0);
}

template <typename T>
static
void c_generate_index_list_generic(
			const T* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, hipStream_t stream)
{
	static int* dev_nvalid = nullptr;
	if (dev_nvalid == nullptr)
			hipMalloc(&dev_nvalid, sizeof(int));

	c_generate_index_list_generic_device(
			dev_conditions, startid, endid, dev_indices, dev_nvalid, 0);

	hipMemcpyAsync(&nvalid, dev_nvalid, sizeof(int), hipMemcpyDeviceToHost, 0);
	hipStreamSynchronize(0);
}

///
/// Exposed functions
/// 
/// Non-batched first
/// 

void c_generate_index_list_i1(
			const char* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, hipStream_t stream)
{
	c_generate_index_list_generic(dev_conditions, startid, endid, dev_indices, nvalid, 0);
}

void c_generate_index_list_i4(
			const int* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, hipStream_t stream)
{
	c_generate_index_list_generic(dev_conditions, startid, endid, dev_indices, nvalid, 0);
}

/// 
/// And now batched
/// 

void c_generate_index_list_batched_i1(
	const int batch_size,
	const char* dev_conditions, const int cond_stride,
	const int startid, const int endid,
	int* dev_indices, const int idx_stride,
	int* dev_nvalid, hipStream_t stream)
{
c_generate_index_list_batched_generic(
		batch_size,
		dev_conditions, cond_stride,
		startid, endid,
		dev_indices, idx_stride,
		dev_nvalid, 0);
}

void c_generate_index_list_batched_i4(
		const int batch_size,
		const int* dev_conditions, const int cond_stride,
		const int startid, const int endid,
		int* dev_indices, const int idx_stride,
		int* dev_nvalid, hipStream_t stream)
{
	c_generate_index_list_batched_generic(
			batch_size,
			dev_conditions, cond_stride,
			startid, endid,
			dev_indices, idx_stride,
			dev_nvalid, 0);
}

void initHIP(int deviceNum){
	hipSetDevice(deviceNum);
}
#endif 
