#include <iostream>
#include <chrono>
#include <omp.h>

#define NUM_THREADS 16
#define PAD 64

using namespace std;

// Function to calculate pi
double calculate_pi (int num_steps){
    int i;
    double x, pi, sum = 0;
    double step = 1./ (double) num_steps;

    for (i = 0; i < num_steps; i++){
        x = (i+0.5)*step;
        sum = sum + 4.0/(1.0+x*x);
    }
    pi = sum * step;
    return pi;
}


// Parallel function to calculate pi
double calculate_pi_wrong (int num_steps){
    double x, pi, sum = 0;
    double step = 1./ (double) num_steps;
    int chunks = num_steps / NUM_THREADS;

    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int i;
        int id = omp_get_thread_num();
        for (i = id*chunks; i < (id+1)*chunks; i++){
            x = (i+0.5)*step;
            sum = sum + 4.0/(1.0+x*x);
        }
    }
    
    pi = sum * step;
    return pi;
}

// Parallel function to calculate pi
double calculate_pi_correct (int num_steps){

    double pi = 0.0;
    double sum [NUM_THREADS];
    double step = 1.0/(double) num_steps;
    int chunks = num_steps / NUM_THREADS;

    // Set up the parallel loop
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        // get the thread id
        int i, id;
        double x;
        id = omp_get_thread_num();
        // Each thread sums over its own chunk
        for ( i=id*chunks, sum[id]=0; i < (id+1)*chunks; i++){
            x = (i+0.5)*step;
            // Sum over each's own sum counter
            sum[id] += 4.0/(1.0+x*x);
        }
    }

    // Combine in serial
    for (int i = 0; i < NUM_THREADS; i++) pi += step * sum[i];
    return pi;

}    


// Padded parallel function to calculate pi
double calculate_pi_padded (int num_steps){

    double pi = 0.0;
    double sum [NUM_THREADS][PAD];
    double step = 1.0/(double) num_steps;
    int chunks = num_steps / NUM_THREADS;

    // Set up the parallel loop
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        // get the thread id
        int i, id;
        double x;
        id = omp_get_thread_num();
        // Each thread sums over its own chunk
        for ( i=id*chunks, sum[id][0]=0; i < (id+1)*chunks; i++){
            x = (i+0.5)*step;
            // Sum over each's own sum counter
            sum[id][0] += 4.0/(1.0+x*x);
        }
    }

    // Combine in serial
    for (int i = 0; i < NUM_THREADS; i++) pi += step * sum[i][0];
    return pi;

}



// parallel function to calculate pi using critical barrier
double calculate_pi_critical (int num_steps){

    double pi = 0.0;
    double step = 1.0/(double) num_steps;
    int chunks = num_steps / NUM_THREADS;
    double sum, x;

    // Set up the parallel loop
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel private(sum, x) shared (pi)
    {
        // get the thread id
        int i, id;
        id = omp_get_thread_num();
        // Each thread sums over its own chunk
        for ( i=id*chunks, sum=0; i < (id+1)*chunks; i++){
            x = (i+0.5)*step;
            // Sum over each's own sum counter
            sum += 4.0/(1.0+x*x);
        }
        
        // Use a critical barrier
        // Only one thread at a time
        #pragma omp critical
        pi += step * sum;

    }
    return pi;

}    



// parallel function to calculate pi using atomic barrier
double calculate_pi_atomic (int num_steps){

    double pi = 0.0;
    double step = 1.0/(double) num_steps;
    int chunks = num_steps / NUM_THREADS;
    double sum, x;

    // Set up the parallel loop
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel private(sum, x) shared(pi)
    {
        // get the thread id
        int i, id;
        id = omp_get_thread_num();
        // Each thread sums over its own chunk
        for ( i=id*chunks, sum=0; i < (id+1)*chunks; i++){
            x = (i+0.5)*step;
            // Sum over each's own sum counter
            sum += 4.0/(1.0+x*x);
        }
        
        // Use an atomic barrier
        // Only one thread at a time can modify
        #pragma omp atomic
        pi += step * sum;

    }
    return pi;

}    




// Padded parallel function to calculate pi
double calculate_pi_reduction (int num_steps){

    double pi = 0.0;
    double step = 1.0/(double) num_steps;
    int i;
    double sum = 0;
    double x;

    // Set up the parallel loop
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction (+:sum) private (x)
    for ( i=0; i < num_steps; i++){
        x = (i+0.5)*step;
        // Sum over each's own sum counter
        sum += 4.0/(1.0+x*x);
    }

    pi = step * sum;

    return pi;

}    

int main(){

    // 1 million steps
    int num_steps = 1000000 ;
    int n_repeat = 500;

    // mean pi value
    double serial_pi = 0;
    double parallel_pi_wrong = 0;    
    double parallel_pi_correct = 0;
    double parallel_pi_padded = 0;
    double parallel_pi_critical = 0;
    double parallel_pi_atomic = 0;
    double parallel_pi_reduction = 0;
    

    // Run timer
    auto start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) serial_pi += calculate_pi(num_steps);
    auto end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> serial_duration = end_time - start_time;

    start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) parallel_pi_wrong += calculate_pi_wrong(num_steps);
    end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> parallel_duration_wrong = end_time - start_time;


    start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) parallel_pi_correct += calculate_pi_correct(num_steps);
    end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> parallel_duration_correct = end_time - start_time;

    start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) parallel_pi_padded += calculate_pi_padded(num_steps);
    end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> parallel_duration_padded = end_time - start_time;


    start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) parallel_pi_critical += calculate_pi_critical(num_steps);
    end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> parallel_duration_critical = end_time - start_time;



    start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) parallel_pi_atomic += calculate_pi_atomic(num_steps);
    end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> parallel_duration_atomic = end_time - start_time;

    start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) parallel_pi_reduction += calculate_pi_reduction(num_steps);
    end_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double> parallel_duration_reduction = end_time - start_time;


    // Print average results
    cout << "Serial Calculation of Pi: " << serial_pi / n_repeat 
              << endl << "Duration: " << serial_duration.count() / n_repeat << " seconds" << endl;
    cout << "Parallel (Wrong) Calculation of Pi: " << parallel_pi_wrong / n_repeat 
              << endl << "Duration: " << parallel_duration_wrong.count() / n_repeat << " seconds (" << serial_duration.count() / parallel_duration_wrong.count() << ") "<< endl;
    cout << "Parallel (Correct) Calculation of Pi: " << parallel_pi_correct / n_repeat 
              << endl << "Duration: " << parallel_duration_correct.count() / n_repeat << " seconds (" << serial_duration.count() / parallel_duration_correct.count() << ") "<< endl;
    cout << "Parallel (Padded) Calculation of Pi: " << parallel_pi_padded / n_repeat 
              << endl << "Duration: " << parallel_duration_padded.count() / n_repeat << " seconds (" << serial_duration.count() / parallel_duration_padded.count() << ") "<< endl;
    cout << "Parallel (Critical) Calculation of Pi: " << parallel_pi_critical / n_repeat 
              << endl << "Duration: " << parallel_duration_critical.count() / n_repeat << " seconds (" << serial_duration.count() / parallel_duration_critical.count() << ") "<< endl;
    cout << "Parallel (Atomic) Calculation of Pi: " << parallel_pi_atomic / n_repeat 
              << endl << "Duration: " << parallel_duration_atomic.count() / n_repeat << " seconds (" << serial_duration.count() / parallel_duration_atomic.count() << ") "<< endl;
        cout << "Parallel (Reduction) Calculation of Pi: " << parallel_pi_reduction / n_repeat 
              << endl << "Duration: " << parallel_duration_reduction.count()  / n_repeat<< " seconds (" << serial_duration.count() / parallel_duration_reduction.count() << ") "<< endl;


}