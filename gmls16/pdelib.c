/**
 * \file pdelib.c
 * \author Guilherme M. Lopes
 * \brief Implementação da das funções da biblioteca <code>pdelib</code>.
 */
#include "pdelib.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * \fn int calculate_index(int i, int j, int n)
 * \brief Função utilitária para acessar a posição (i, j) em matrizes alocadas como vetores.
 * \param i: Linha
 * \param j: Coluna
 * \param n: Número de linhas da matriz
 * \details Matrizes alocadas como vetores possuem todos os seus campos alocados de forma\n
 * contínua, como se fossem <code>arrays</code>. Para acessar a posição (i, j) de uma matriz\n
 * alocada desta forma, basta usar esta função.
 */
int calculate_index(int i, int j, int n) {
    return i * n + j;
}

/**
 * \fn double timestamp(void)
 * \brief Função utilitária para calcular o tempo de execução.
 */
double timestamp(void) {
    struct timeval time_pointer;
    gettimeofday(&time_pointer, NULL);

    return((double)(time_pointer.tv_sec * 1000.0 + time_pointer.tv_usec / 1000.0));
}

/**
 * \fn int get_options(int argc, char **argv, int *nx, int *ny, int *max_iterations, char **output_file)
 * \brief Lê e guarda os argumentos passados pela linha de comando.
 * \param argc: Quantidade de argumentos passados;
 * \param argv: Argumentos;
 * \param nx: Referência à variável nx, que representa a quantidade de pontos a serem representados no eixo X;
 * \param ny: Referência à variável ny, que representa a quantidade de pontos a serem representados no eixo Y;
 * \param max_iterations: Número máximo de iterações para o Método de Gauss-Seidel;
 * \param output_file: Nome do arquivo de saída que deve conter a solução;
 */
int get_options(int argc, char **argv, int *nx, int *ny, int *max_iterations, char **output_file) {
    // Argumentos do programa
    const struct option stopcoes[] = {
        {"nx", required_argument,  0, 'a'},
        {"ny", required_argument,  0, 'b'},
        {"i", required_argument,   0, 'i'},
        {"o", required_argument,   0, 'o'},
        {0, 0, 0, 0},
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "ha:b:i:o:", stopcoes, NULL)) > 0) {
        switch (opt) {
            case 'h': /* -h ou --help */
                show_help(argv[0]);
                break;
            case 'a': /* -a ou --nx */
                *nx = atoi(optarg);
                break;
            case 'b': /* -b ou --nx */
                *ny = atoi(optarg);
                break;
            case 'i': /* -i ou --i */
                *max_iterations = atoi(optarg);
                break;
            case 'o': /* -o ou --o */
                *output_file = optarg;
                break;
            default:
                fprintf(stderr, "Opção inválida ou argumento faltando: `%c'\n", optopt) ;
                return -1;
        }
    }

    return 0;
}

/**
 * \fn void allocate_and_start_linear_system(t_LS5Diag **SL, int nx, int ny)
 * \brief Aloca e inicializa o sistema linear.
 * \param SL: Ponteiro para a estrutura do sistema linear;
 * \param nx: Quantidade de pontos para discretização no eixo das abscissas;
 * \param ny: Quantidade de pontos para discretização no eixo das ordenadas.
 */
void allocate_and_start_linear_system(t_LS5Diag **SL, int nx, int ny) {
    *SL = malloc(sizeof(t_LS5Diag));
    (*SL)->b = malloc(sizeof(t_float) * (nx + 2 * ny + 2));
    (*SL)->n = nx * ny;
    const float_t hx = pi / (nx + 1);
    const float_t hy = pi / (ny + 1);

    // Vetor B
    t_float x = 0;
    t_float y = 0;
    for (int i = 0; i < (nx + 2); i++) {
        y = 0;
        for (int j = 0; j < (ny + 2); j++) {
            (*SL)->b[calculate_index(i, j, (ny + 2))] = 4 * pi * pi * (
                sin(2 * pi * x) *
                sinh(pi * y) +
                sin(2 * pi * (pi - x)) *
                sinh(pi * (pi - y))
            );
            y += hy;
        }
        x += hx;
    }
    (*SL)->ny = ny;
    (*SL)->nx = nx;

    // Matriz A
    (*SL)->away_upper_diagonal = - ((1 / (hx * hx)) + (1 / (2 * hx)));
    (*SL)->upper_diagonal = - ((1 / (hy * hy)) + (1 / (2 * hy)));
    (*SL)->main_diagonal = ((4 * pi * pi) + (2 / (hx * hx)) + (2 / (hy * hy)));
    (*SL)->bottom_diagonal = - ((1 / (hy * hy)) - (1 / (2 * hy)));
    (*SL)->away_upper_diagonal = - ((1 / (hx * hx)) - (1 / (2 * hx)));
}

/**
 * \fn void gaussSeidel(t_LS5Diag **SL, t_float **u )
 * \brief Método de Gauss-Seidel para encontrar uma solução para o sistema linear de forma iterativa.
 * \param SL: Estrutura representando um sistema linear;
 * \param u: vetor solução;
 */
void gaussSeidel(t_LS5Diag **SL, t_float **u ) {
    t_float Uij;
    for (int i = 1; i < (((*SL)->nx) + 1); ++i) {
        for (int j = 1; j < (((*SL)->ny)+1); ++j) {

            Uij = (*SL)->b[ calculate_index(i,j,(((*SL)->ny)+2)) ];

            Uij -= (((*SL)->away_bottom_diagonal) * (*u)[ calculate_index(i, j - 1, ((*SL)->ny) + 2) ]);

            Uij -= (((*SL)->bottom_diagonal) * (*u)[ calculate_index(i - 1, j, ((*SL)->ny) + 2) ]);

            Uij -= (((*SL)->away_upper_diagonal) * (*u)[ calculate_index(i, j + 1, ((*SL)->ny) + 2) ]);

            Uij -= (((*SL)->upper_diagonal) * (*u)[ calculate_index(i + 1, j, ((*SL)->ny) + 2) ]);

            Uij /= ((*SL)->main_diagonal);

            (*u)[ calculate_index(i, j, ((*SL)->ny) + 2) ] = Uij;

        }
    }
}

/**
 * \fn void show_help(char *name)
 * \brief Função que imprime a ajuda na saída padrão.
 * \param name: Nome do programa;
 */
void show_help(char *name) {
    fprintf(stderr, "\
        [uso] %s <opções>\n\
        --nx INTEGER     Número de pontos da malha na direção de x.\n\
        --ny INTEGER     Número de pontos da malha na direção de y.\n\
        --i INTEGER      Número máximo de iterações do método de Gauss-Seidel.\n\
        --o STRING       Nome do arquivo de saída.\n\
        ", name);

    exit(-1);
}

/**
 * \fn void generateOutputFile(
 *          unsigned int iterations,
 *          char* output_file,
 *          double gauss_seidel_total_time,
 *          t_float *residues,
 *          t_LS5Diag *SL,
 *          t_float *u
 *      )
 * \brief Gera o arquivo de saída.
 *  \param iterations: Número de iterações;
 *  \param output_file: Nome do arquivo de saída;
 *  \param gauss_seidel_total_time: Tempo de execução do método de Gauss-Seidel;
 *  \param residues: Vetor de resíduos;
 *  \param SL: Estrutura do Sistema Linear;
 *  \param u: Vetor de solução.
 *  \details Esta função gera um arquivo de saída com o nome <code>output_file</code> que serve de entrada\n
 *  para o programa Gnuplot. No começo do arquivo, esta função escreve um comentário com o tempo de execução do\n
 *  método de Gauss-Seidel e o resíduo para cada iteração.
 */
void generateOutputFile(
    unsigned int iterations,
    char* output_file,
    double gauss_seidel_total_time,
    t_float *residues,
    t_LS5Diag *SL,
    t_float *u
) {
    const t_float hx = pi/((SL->nx)+1);
    const t_float hy = pi/((SL->ny)+1);
    t_float x = 0;
    t_float y = 0;

    FILE *output = fopen(output_file, "w+");

    fprintf(output, "###########\n");
    fprintf(output, "# Tempo Gauss-Seidel: %f ms.\n", ( gauss_seidel_total_time / iterations ));
    fprintf(output, "#\n# Norma L2 do resíduo\n");

    for (unsigned int i = 1; i <= iterations; i++) {
        fprintf(output, "# i = %d: %lf\n", i, residues[i-1]);
    }
    fprintf(output, "###########\n\n");

    for (int i = 0; i < ((SL->nx) + 2); i++) {
        y = 0;
        for(int j = 0; j < ((SL->ny) + 2); j++) {
            fprintf(output, "%lf %lf %lf\n", x, y, u[ calculate_index(i, j, (SL->ny) + 2) ]);
            y += hy;
        }
        fprintf(output, "\n" );
        x += hx;
    }

    fclose(output);
}


/**
 * \fn void allocate_and_start_solution(t_float **u, t_LS5Diag *SL)
 * \brief Método de Gauss-Seidel para encontrar uma solução para o sistema linear de forma iterativa.
 * \param SL: Estrutura representando um sistema linear;
 * \param u: vetor solução;
 */
void allocate_and_start_solution(t_float **u, t_LS5Diag *SL) {
    *u = malloc(sizeof(t_float) * ((SL->nx) + 2) * ((SL->ny) + 2));
    const t_float hx = pi / ((SL->nx) + 1);
    const t_float hy = pi / ((SL->ny) + 1);
    t_float x = 0;
    t_float y = 0;

    for (int i = 0; i < ((SL->nx) + 2); i++) {
        y = 0;
        for (int j = 0; j < ((SL->ny) + 2); j++) {
            if (i == 0 || i == ((SL->nx) + 1)) {
                (*u)[calculate_index(i, j, (SL->ny) + 2)] = 0;
            } else if (j == 0) {
                (*u)[calculate_index(i, j, (SL->ny) + 2)] = sin(2 * pi * (pi - x)) * sinh(pi * pi);
            } else if (j == ((SL->ny) + 1)) {
                (*u)[calculate_index(i, j, (SL->ny) + 2)] = sin(2 * pi * x) * sinh(pi * pi);
            } else {
                (*u)[calculate_index(i, j, (SL->ny) + 2)] = 0;
            }
            y += hy;
        }
        x += hx;
    }
}

/**
 * \fn t_float calculate_residues(t_LS5Diag *SL, t_float *u)
 * \brief Calcula os resíduos de cada iteração e retorna a norma L2.
 * \param SL: Estrutura representando um sistema linear;
 * \param u: vetor solução;
 */
t_float calculate_residues(t_LS5Diag *SL, t_float *u) {
    t_float *residue = malloc(sizeof(t_float) * ((SL->nx) + 2) * ((SL->ny) + 2));
    t_float absolute_value = 0;

    for(int i = 1; i < ((SL->nx) + 1); i++) {
        for(int j = 1; j < ((SL->ny) + 1); j++) {
            residue[ calculate_index(i, j, (SL->ny) + 2) ] = (SL->b[ calculate_index(i, j, (SL->ny) + 2) ]);
            residue[ calculate_index(i, j, (SL->ny) + 2) ] -= ((SL->away_bottom_diagonal) * u[ calculate_index(i, j - 1, (SL->ny) + 2) ]);
            residue[ calculate_index(i, j, (SL->ny) + 2) ] -= ((SL->bottom_diagonal) * u[ calculate_index(i - 1, j, (SL->ny) + 2) ]);
            residue[ calculate_index(i, j, (SL->ny) + 2) ] -= ((SL->away_upper_diagonal) * u[ calculate_index(i, j + 1,( SL->ny) + 2) ]);
            residue[ calculate_index(i, j, (SL->ny) + 2) ] -= ((SL->upper_diagonal) * u[ calculate_index(i + 1, j, (SL->ny) + 2) ]);
            residue[ calculate_index(i, j, (SL->ny) + 2) ] -= (SL->main_diagonal) * u[ calculate_index(i, j, (SL->ny) + 2) ];

            absolute_value += residue[ calculate_index(i, j, (SL->ny) + 2) ] * residue[ calculate_index(i, j, (SL->ny) + 2) ];
        }
    }
    absolute_value = sqrt(absolute_value);

    return absolute_value;
}