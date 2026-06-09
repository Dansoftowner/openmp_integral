#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define M_PIl   3.141592653589793238462643383279502884L 
#define M_PI_2l 1.570796326794896619231321691639751442L 

long double gettime()
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (long double)ts.tv_sec + (long double)ts.tv_nsec*0.000000001;
}

long double gauss2d(double x, double y) {
    return expl(-x * x - y * y);
}

long double gauss2d_transformed(long double u, long double v) {
    long double tan_u = tanl(u);
    long double tan_v = tanl(v);
    long double cos_u = cosl(u);
    long double cos_v = cosl(v);

    return expl(-(tan_u*tan_u + tan_v*tan_v)) / (cos_u*cos_u * cos_v*cos_v);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Használat: %s <N>\n", argv[0]);
        return 1;
    }

    long long N_samples = atoll(argv[1]);

    const long double a = -M_PI_2l;
    const long double b = M_PI_2l;
    
    const long double area = (b - a) * (b - a);

    long double start_time = gettime();
    long double sum = 0.0;

    #pragma omp parallel reduction(+:sum)
    {
        unsigned int seed = time(NULL) + omp_get_thread_num();

        #pragma omp for
        for (long long i = 0; i < N_samples; i++) {

            long double u = a + (long double)rand_r(&seed) / RAND_MAX * (b - a);
            long double v = a + (long double)rand_r(&seed) / RAND_MAX * (b - a);
            
            sum += gauss2d_transformed(u, v);
        }
    }

    long double result = area * (sum / N_samples);

    long double end_time = gettime();
    long double seq_time = end_time - start_time;

    printf("Minták száma: %lld\n", N_samples);
    printf("Eredmény: %.37Lf\n", result);
    printf("Hiba: %.37Lf\n", fabsl(M_PIl - result));
    printf("Idő: %.9Lf mp\n", seq_time);

    return 0;
}