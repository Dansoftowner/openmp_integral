#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define	M_PIl	3.141592653589793238462643383279502884L 
#define	M_PI_2l	1.570796326794896619231321691639751442L 

long double gettime()
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (long double)ts.tv_sec + (long double)ts.tv_nsec*0.000000001;
}

long double gauss2d(double x, double y) {
    return expl(-x * x - y * y);
}

long double gauss2d_transformed(double u, double v) {
    long double tan_u = tanl(u);
    long double tan_v = tanl(v);
    long double cos_u = cosl(u);
    long double cos_v = cosl(v);

    return expl(-(tan_u*tan_u + tan_v*tan_v)) / (cos_u*cos_u * cos_v*cos_v);
}

int main(int argc, char** argv) {
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);

    const long double a = -M_PI_2l;
    const long double b = M_PI_2l;
    const long double c = -M_PI_2l;
    const long double d = M_PI_2l;

    long double dx = (b - a) / N;
    long double dy = (d - c) / M;

    long double start_time = gettime();

    long double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            long double x = a + (i + 0.5) * dx;
            long double y = c + (j + 0.5) * dy;
            sum += gauss2d_transformed(x, y);
        }
    }
    sum *= dx * dy;

    long double end_time = gettime();
    long double seq_time = end_time - start_time;

    printf("Eredmény: %.37Lf\n", sum);
    printf("Hiba: %.37Lf\n", fabsl(M_PIl - sum));
    printf("Idő: %.9Lf\n", seq_time);

    return 0;
}