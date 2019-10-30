/**
 * \file main.c
 * \author Guilherme M. Lopes
 * \brief Função principal do programa, utiliza a biblioteca <code>pdelib</code> para resolver uma PDE.
 *
 * \mainpage A program to find an approximate solution to a Partial Differential Equation using the Gauss-Seidel iterative method.
 *
 * \section introSec About the program
 * Given a Partial Differential Equation, this program discretizes the domain using the Central Finite Differences method.\n
 * Then, it uses Gauss-Seidel method to find a solution of the generated linear system.
 *
 * **Matéria/Course:** \n
 * Introdução à Computação Científica <i>(Introduction to Scientific Computing)</i> \n
 * **Professor:** \n
 * Daniel Weingaertner \n
 * **Aluno/Student:** \n
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *
 * \date 24 out 2019
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pdelib.h"
//#define DEBUG

int main(int argc, char *argv[]) {
    int nx; // M
    int ny; // N
    int max_iterations;
    char *output_file;

    if (argc < 2) {
        show_help(argv[0]);
    }

    get_options(argc, argv, &nx, &ny, &max_iterations, &output_file);

    t_LS5Diag *SL; // Sistema Linear Pentadiagonal
    allocate_and_start_linear_system(&SL, nx, ny);

    t_float *u; // Solução
    allocate_and_start_solution(&u, SL);

    // Residuos
    t_float *residues = malloc(sizeof(t_float)*max_iterations);

    double initial_time;
    double gauss_seidel_total_time = 0;

    unsigned int k = 0;
    do {
        k++;
        initial_time = timestamp();
        gaussSeidel(&SL, &u);
        gauss_seidel_total_time += (timestamp() - initial_time);
        // Calcula Residuo da iteração
        residues[k-1] = calculate_residues(SL, u);
    } while(((k <= max_iterations) && (residues[k-1] >= MAXIMUM_ERROR)));

//    if (residues[k-1] >= MAXIMUM_ERROR) {
//        perror("Não convergiu\n");
//
//        return (-1);
//    }

    generateOutputFile(k, output_file, gauss_seidel_total_time, residues, SL, u);

    return 0;
}
