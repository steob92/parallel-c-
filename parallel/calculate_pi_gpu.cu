#include <iostream>
#include <chrono>
#include <omp.h>

// Function to calculate pi
__global__
void assemble_pi (int num_steps, float *pi){

    float step = 1.0/(float) num_steps;
    float x;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < num_steps; i += stride){
        x = (i+0.5)*step;
        // Sum over each's own sum counter
        // sum[i] = 4.0/(1.0+x*x) * step;
        atomicAdd(&pi[0], 4.0/(1.0+x*x) * step);
    }
}

float calculate_pi (int num_steps){

    float *_pi;
    float pi;
    float step = 1.0/(float) num_steps;

    int blockSize = 1<<10;
    int numBlocks = (num_steps + blockSize - 1) / blockSize;
    // cudaMallocManaged(&sum, num_steps*sizeof(float));
    cudaMallocManaged(&_pi, sizeof(float));
    _pi[0] = 0;
    assemble_pi<<<numBlocks, blockSize>>>(num_steps, _pi);

    cudaDeviceSynchronize();
    pi = _pi[0];
    cudaFree(_pi);

    return pi;
}



int main(){

    // 1 million steps
    int num_steps = 1000000;
    int n_repeat = 500;

    // mean pi value
    float serial_pi = 0;

    // Run timer
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) serial_pi += calculate_pi(num_steps);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<float> serial_duration = end_time - start_time;

    // Print average results
    std::cout << "Serial Calculation of Pi: " << serial_pi / n_repeat 
              << std::endl << "Duration: " << serial_duration.count() / n_repeat<< " seconds" << std::endl;

}