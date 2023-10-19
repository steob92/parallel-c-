#include <iostream>
#include <chrono>

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



int main(){

    // 1 million steps
    int num_steps = 1000000;
    int n_repeat = 50;

    // mean pi value
    double serial_pi = 0;

    // Run timer
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) serial_pi += calculate_pi(num_steps);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> serial_duration = end_time - start_time;

    // Print average results
    std::cout << "Serial Calculation of Pi: " << serial_pi / n_repeat 
              << std::endl << "Duration: " << serial_duration.count() /n_repeat << " seconds" << std::endl;

}