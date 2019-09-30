#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>

using namespace std;

cudaError_t cu_add(int* a, int* b, int* c, int arr_size);
cudaError_t cu_multi(float* c_m, float* a_m, float* b_m, int ha, int n, int wb);

__global__ void kernel_add(int* a, int* b, int* c) {
	int i = threadIdx.x;
	c[i] = a[i] + b[i];
}


__global__ void kernel_multi(float* c_m, float* a_m, float* b_m, int ha, int n, int wb) {
	int tid_x = threadIdx.x;
	int tid_y = threadIdx.y;

	int t_id = (blockIdx.y*blockDim.y + threadIdx.y)*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	int row = blockIdx.y*blockDim.y + threadIdx.y;
	int column = blockIdx.x*blockDim.x + threadIdx.x;

	float sum = 0;
	for (int i = 0; i != n; ++i) {
		//printf("(%f,%f,%d,%d,%d)\n", a_m[tid_y * n + i], b_m[i * n + tid_x],tid_x, tid_y, i);
		//sum += a_m[tid_y * n + i] * b_m[i * wb + tid_x];
		sum += a_m[row*n + i] * b_m[i*wb + column];
	}
	//c_m[tid_y*wb + tid_x] = sum;
	c_m[t_id] = sum;
}


cudaError_t cu_add(int* a, int* b, int* c, int arr_size) {
	int* dev_a = 0;
	int* dev_b = 0;
	int* dev_c = 0;

	cudaError_t cuda_stat;
	cuda_stat = cudaSetDevice(0);

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	cout << "device name: " << deviceProp.name << endl;
	cout << "the warpsize: " << deviceProp.warpSize << endl;
	cout << "the maxThreadPerBlock: " << deviceProp.maxThreadsPerBlock << endl;
	cout << "the maxThreadsDim: " << deviceProp.maxThreadsDim[0] << endl;
	cout << "the maxGridSize: " << deviceProp.maxGridSize[3] << endl;
	cout << "the compute capability: " << deviceProp.major << '.' << deviceProp.minor << endl;

	if (cuda_stat != cudaSuccess) {
		cout << "error in cudaSetDevice" << endl;
		goto Error;
	}

	cudaMalloc((void**)&dev_a, arr_size * sizeof(int));
	cudaMalloc((void**)&dev_b, arr_size * sizeof(int));
	cudaMalloc((void**)&dev_c, arr_size * sizeof(int));

	cuda_stat = cudaMemcpy(dev_a, a, arr_size * sizeof(int), cudaMemcpyHostToDevice);
	cuda_stat = cudaMemcpy(dev_b, b, arr_size * sizeof(int), cudaMemcpyHostToDevice);
	cuda_stat = cudaMemcpy(dev_c, c, arr_size * sizeof(int), cudaMemcpyHostToDevice);

	kernel_add << <1, arr_size >> > (dev_a, dev_b, dev_c);

	cuda_stat = cudaGetLastError();
	if (cuda_stat != cudaSuccess) {
		cout << "error in calling the kernel_add" << endl;
		goto Error;
	}

	cuda_stat = cudaDeviceSynchronize();
	cuda_stat = cudaMemcpy(c, dev_c, arr_size * sizeof(int), cudaMemcpyDeviceToHost);
	if (cuda_stat != cudaSuccess) {
		cout << "error in copy mem from D to H" << endl;
		goto Error;
	}

Error:
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);

	return cuda_stat;
}

cudaError_t cu_multi(float* c_m, float* a_m, float* b_m, int ha, int n, int wb) {
	float* dev_a = 0;
	float* dev_b = 0;
	float* dev_c = 0;

	dim3 threads(ha, wb);

	cudaError_t cuda_stat;
	cuda_stat = cudaSetDevice(0);
	if (cuda_stat != cudaSuccess) {
		cout << "error in setting device" << endl;
		goto Error;
	}

	cudaMalloc((void**)&dev_a, ha * n * sizeof(float));
	cudaMalloc((void**)&dev_b, n * wb * sizeof(float));
	cudaMalloc((void**)&dev_c, ha * wb * sizeof(float));

	cudaMemcpy(dev_a, a_m, ha * n * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b_m, n * wb * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_c, c_m, ha * wb * sizeof(float), cudaMemcpyHostToDevice);

	kernel_multi << <1, threads >> > (dev_c, dev_a, dev_b, ha, n, wb);
	//matrixmultiplicationFunc << <1, threads >> > (dev_c, dev_a, dev_b, ha, n, wb);

	cuda_stat = cudaGetLastError();
	if (cuda_stat != cudaSuccess) {
		cout << "error in calling kernel_multi" << endl;
		goto Error;
	}

	cuda_stat = cudaDeviceSynchronize();
	cuda_stat = cudaMemcpy(c_m, dev_c, ha * wb * sizeof(float), cudaMemcpyDeviceToHost);
	if (cuda_stat != cudaSuccess) {
		cout << "error in copy from D to H";
		goto Error;
	}

Error:
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);

	return cuda_stat;
}