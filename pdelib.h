#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <sys/time.h>

#define pi 3.14159265359
typedef double t_float;

// Linear System Pentadiagonal Type
typedef struct t_LS5Diag {
    double *main_diagonal,
           *bottom_diagonal,
           *upper_diagonal,
           *away_bottom_diagonal,
           *away_upper_diagonal;
    // Termos independentes
    double *b;
    // Número de equações do SL
    unsigned int n;
} t_LS5Diag;

t_float f(t_float xi, t_float yj);

void allocatePentadiagonalLinearSystem(t_LS5Diag *SL, int nx, int ny);

void gaussSeidel(t_LS5Diag *SL, t_float *x, t_float error, int max_iterations, double *iterations_timestamp, t_float *residues);

void show_help(char *name);

void generateOuputFile(t_float *x, int n, double averageTimeGS, char *filename, t_float *residues, int max_iterations);

double timestamp(void);

double averageTimeGaussSeidel(double *iterations_timestamp, int max_iterations);
