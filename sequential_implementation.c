#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923

long double gettime()
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (long double)ts.tv_sec + (long double)ts.tv_nsec*0.000000001;
}

double f(double x, double y) {
    return exp(-x * x - y * y);
}

double f_transformed(double u, double v) {
    double tan_u = tan(u);
    double tan_v = tan(v);
    double cos_u = cos(u);
    double cos_v = cos(v);

    return exp(-(tan_u*tan_u + tan_v*tan_v)) / (cos_u*cos_u * cos_v*cos_v);
}

int main(int argc, char** argv) {
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);

    const double a = -M_PI_2;
    const double b = M_PI_2;
    const double c = -M_PI_2;
    const double d = M_PI_2;

    double dx = (b - a) / N;
    double dy = (d - c) / M;

    long double start_time = gettime();

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            double x = a + (i + 0.5) * dx;
            double y = c + (j + 0.5) * dy;
            sum += f_transformed(x, y);
        }
    }
    sum *= dx * dy;

    long double end_time = gettime();
    long double seq_time = end_time - start_time;

    printf("Eredmény: %.15f\n", sum);
    printf("Hiba: %.15f\n", fabs(M_PI - sum));
    printf("Idő: %.9Lf\n", seq_time);

    return 0;
}