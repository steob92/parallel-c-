#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 16

struct Point{
    double x;
    double y;
    double z;
};

std::vector <double>  distance_to_point_vector(
    std::vector<double> x, 
    std::vector<double> y, 
    std::vector<double> z
    )
{
    std::vector<double> c(x.size());
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for
    for (int i=0; i < x.size(); i++) {
        c[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
    }
    return c;
}


std::vector <double>  distance_to_point_struct(
    std::vector<Point> p
    )
{
    std::vector<double> c(p.size());
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for
    for (int i=0; i < p.size(); i++){
        c[i] = sqrt(p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z);
    } 
    return c;
}

int main(){

    // 1 million steps
    int num_entries = 1000000;
    int n_repeat = 500;
    int i = 0;

    std::vector <double> a(num_entries);
    std::vector <double> b(num_entries);
    std::vector <double> c(num_entries);
    std::vector <double> dist;

    std::vector <Point> point(num_entries);

    for (i = 0; i < num_entries; i ++){
        a[i] = i;
        b[i] = i^2;
        c[i] = i/2.;
        point[i].x = i;
        point[i].y = i^2;
        point[i].z = i/2.;
        
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) dist = distance_to_point_vector(a,b,c );
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> vector_duration = end_time - start_time;

    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_repeat ; i ++) dist = distance_to_point_struct(point);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> struct_duration = end_time - start_time;

    std::cout << "Vector run time: " << vector_duration.count() / n_repeat << " seconds" << std::endl;
    std::cout << "Struct run time: " << struct_duration.count() / n_repeat << " seconds" << std::endl;
    std::cout << "Speedup: " <<  vector_duration.count()  /struct_duration.count() << std::endl;
    

}