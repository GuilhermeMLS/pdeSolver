#include "pdelib.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Função f(xi, yi)
t_float f(t_float xi, t_float yj)
{
    return 4 * pow(pi, 2) * (
        sin(2 * pi * xi) *
        sinh(pi * yj) +
        sin(2 * pi * (pi - xi)) *
        sinh(pi * (pi - yj))
    );
}

// Retorna tempo em milisegundos
double timestamp(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

void allocatePentadiagonalLinearSystem(t_LS5Diag *SL, int nx, int ny)
{
    SL->away_upper_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 2));
    SL->upper_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 1));
    SL->main_diagonal = malloc(sizeof(t_float) * (nx * ny));
    SL->bottom_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 1));
    SL->away_bottom_diagonal = malloc(sizeof(t_float) * ((nx * ny) - 2));
    // Termos independentes
    SL->b = malloc(sizeof(t_float) * (nx * ny));
    // Número de equações do Sistema Linear
    SL->n = nx * ny;
}

void gaussSeidel(
    t_LS5Diag *SL,
    t_float *x,
    t_float error,
    int max_iterations,
    double *iterations_timestamp,
    t_float *residues,
    int nx
) {
    double norma, diff, xk, initial_time, final_time;
    int k = 1, i;

    do {
        initial_time = timestamp();
        // primeira equação fora do laço
        i = 0;
        xk = (SL->b[i] - (SL->upper_diagonal[i] * x[i+1]) -
            SL->away_upper_diagonal[i] * x[i+2]
        ) / SL->main_diagonal[i];
        norma = fabs(xk - x[0]);
        x[i] = xk;

        // Segunda equação fora do laço
        i++;
        xk = (SL->b[i] - (
            SL->bottom_diagonal[i-1] * x[i-1]) -
            SL->upper_diagonal[i] * x[i+1] -
            SL->away_upper_diagonal[i] * x[i+2]
        ) / SL->main_diagonal[i];
        norma = fabs(xk - x[1]);
        x[i] = xk;

        for (i = 2; i < (SL->n) - 2; i++) {
            xk = (SL->b[i] - SL->away_bottom_diagonal[i-2] * x[i-2] -
                SL->bottom_diagonal[i-1] * x[i-1] -
                SL->upper_diagonal[i] * x[i+1] -
                SL->away_upper_diagonal[i] * x[i+2]
            ) / SL->main_diagonal[i];
            diff = fabs(xk - x[i]);
            norma = (diff > norma) ? (diff) : (norma);
            x[i] = xk;
        }

        // Penúltima equação fora do laço
        i++;
        xk = (SL->b[i] - SL->away_bottom_diagonal[i-2]*x[i-2] -
            SL->bottom_diagonal[i-1]*x[i-1] -
            SL->upper_diagonal[i]*x[i+1]
        ) / SL->main_diagonal[i];
        diff = fabs(xk - x[i]);
        norma = (diff > norma) ? (diff) : (norma);
        x[i] = xk;

        // Última equação fora do laço
        i++;
        xk = (SL->b[i] - SL->away_bottom_diagonal[i-2]*x[i-2] -
            SL->bottom_diagonal[i-1] * x[i-1]
        ) / SL->main_diagonal[i];
        diff = fabs(xk - x[i]);
        norma = (diff > norma) ? (diff) : (norma);
        x[i] = xk;

        final_time = timestamp();
        iterations_timestamp[k] = final_time - initial_time;

        //TODO: verificar esse resíduo

        // t_float *R = calloc(SL->n, sizeof(t_float));
        residues[k-1] = calculateResidue(SL, x, nx);

        k++;
        //TODO: implementar critério de parada com base no erro
    //} while (norma > error && k < max_iterations) ;
    } while (k < max_iterations);
}

void show_help(char *name)
{
    fprintf(stderr, "\
        [uso] %s <opções>\n\
        --nx INTEGER     Número de pontos da malha na direção de x.\n\
        --ny INTEGER     Número de pontos da malha na direção de y.\n\
        --i INTEGER      Número máximo de iterações do método de Gauss-Seidel.\n\
        --o STRING       Nome do arquivo de saída.\n\
        ", name);
    exit(-1);
}

void generateOuputFile(
    t_float *x,
    int n,
    double averageTimeGS,
    char *filename,
    t_float *residues,
    int max_iterations,
    int nx,
    int ny,
    int hx,
    int hy
) {
    FILE *file_pointer;
    file_pointer = fopen(filename, "w");

    // Firulas iniciais
    fprintf(file_pointer, "###########\n");
    fprintf(file_pointer, "# Tempo Método GS: %lf\n", averageTimeGS);
    fprintf(file_pointer, "#\n");
    fprintf(file_pointer, "# Norma L2 do Residuo\n");
    for (int i = 0; i < max_iterations; ++i) {
        fprintf(file_pointer, "# i=%d: %lf\n", i, residues[i]);
    }
    fprintf(file_pointer, "###########\n\n");

    //valores
    fprintf(file_pointer, "| x | y | z |\n");
//    for (int i = 0; i < n; i++) {
//        fprintf(file_pointer, "%lf %lf %lf",x[i]);
//    }
    for (int i = 1; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            fprintf(file_pointer, "%d %d %lf\n", i*hx, j*hy, x[i]);
        }
    }

    fclose(file_pointer);
}

double averageTimeGaussSeidel(double *iterations_timestamp, int max_iterations)
{
    double average = 0;
    for(int i = 0; i < max_iterations; i++) {
        average += iterations_timestamp[i] / max_iterations;
    }

    return average;
}

t_float calculateResidue (t_LS5Diag *SL, t_float *x, int nx)
{
    t_float R[SL->n];
    for (int j = 0; j < SL->n; ++j) {
        R[j] = 0.0;
    }
    int i = 1;

    // primeira equação
    R[i] = SL->b[i] - (SL->main_diagonal[i] * x[i] + SL->upper_diagonal[i] * x[i+1] + SL->away_upper_diagonal[i] * x[i+nx]);

    // equações centrais
    for (int i = 2; i < SL->n; ++i) {
        if (i > nx) {
            R[i] = SL->b[i] - ( SL->main_diagonal[i] * x[i] +
                    SL->away_bottom_diagonal[i] * x[i - nx] +
                    SL->bottom_diagonal[i] * x[i-1] +
                    SL->upper_diagonal[i] * x[i+1] +
                    SL->away_upper_diagonal[i] * x[i+nx]);
        } else {
            R[i] = SL->b[i] - ( SL->main_diagonal[i] * x[i] +
                                SL->bottom_diagonal[i] * x[i-1] +
                                SL->upper_diagonal[i] * x[i+1] +
                                SL->away_upper_diagonal[i] * x[i+nx]);
        }
    }
    // última equação
    R[i] = SL->b[i] - (SL->main_diagonal[i] * x[i] +
            SL->away_bottom_diagonal[i] * x[i-nx] +
            SL->bottom_diagonal[i] * x[i-1]);

    t_float sum = 0.0;
    for (i = 1; i <= SL->n; i++) {
        sum += R[i] * R[i];
    }

    return sqrt(sum);
}