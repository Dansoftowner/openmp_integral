#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define M_PIl 3.141592653589793238462643383279502884L 

long double gettime()
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (long double)ts.tv_sec + (long double)ts.tv_nsec*0.000000001;
}

long double gauss2d(long double x, long double y) {
    return expl(-x * x - y * y);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        printf("Használat: %s <R_sugar> <N_r_felosztas> <N_phi_felosztas>\n", argv[0]);
        return 1;
    }

    long double R = atof(argv[1]);
    long long N_r = atoll(argv[2]);
    long long N_phi = atoll(argv[3]);

    long double dr = R / N_r;
    long double dphi = (2.0L * M_PIl) / N_phi;

    long double start_time = gettime();

    long double sum = 0.0;
    
    #pragma omp parallel for collapse(2) reduction(+:sum)
    for (long long i = 0; i < N_r; i++) {
        for (long long j = 0; j < N_phi; j++) {
            
            long double r_mid = (i + 0.5L) * dr;
            //long double phi_mid = (j + 0.5L) * dphi;

            long double f_val = expl(-r_mid * r_mid);

            /*
            long double x = r_mid * cosl(phi_mid);
            long double y = r_mid * sinl(phi_mid);
            long double f_val = gauss2d(x, y);
            */

            sum += f_val * r_mid;
        }
    }
    
    sum *= (dr * dphi);

    long double end_time = gettime();
    long double parallel_time = end_time - start_time;

    // pi * (1 - e^(-R^2))
    long double analytical = M_PIl * (1.0L - expl(-R * R));

    printf("Eredmény: %.37Lf\n", sum);
    printf("Analitikus: %.37Lf\n", analytical);
    printf("Hiba: %.37Lf\n", fabsl(analytical - sum));
    printf("Idő: %.9Lf másodperc\n", parallel_time);

    return 0;
}
