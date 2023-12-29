#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>

#define MIN_RANDOM 0
#define MAX_RANDOM 100

double* vectorXMatrix(int column, int row, const double* vector, double** matrix, double* timeElapsed);

// Funzioni accessorie
double getRandomDoubleNumberInRange(int min, int max);
double** getMatrixOfRandomNumbersOfSize(int column, int row);
double* getVectorOfRandomNumbersOfSize(int size);
bool parseInt(char* str, int* val);

/**
 * Gli argomenti vengono passati nella forma: column row numberOfThreads
 * @param column numero di colonne della matrice
 * @param row numero di righe della matrice
 */
int main(int argc, char** argv) {
    srand(time(NULL));

    // Dichiarazione delle variabili
    int i, column, row;
    double* result = NULL;
    double timeElapsed = 0;

    // assegnazione di valori alle dimensioni
    if (!(parseInt(argv[1], &column)) || !(parseInt(argv[2], &row))) {
        fprintf(stderr, "ERROR - Cannot parse the number of column or rows with values provided.\n"
                        "column:%s\nrow:%s\n\n", argv[1], argv[2]);
        return 1;
    }

    if (column < 1 || row < 1) {
        fprintf(stderr, "ERROR - Column and row number less than 1 are not allowed.\ncolumn:%d\nrow:%d\n", column, row);
        return 1;
    }

    // assegnazione dei valori alla matrice A ed al vettore vector
    double** matrix = getMatrixOfRandomNumbersOfSize(column, row);
    double* vector = getVectorOfRandomNumbersOfSize(column);

    // Ottiene e stampa il risultato
    result = vectorXMatrix(column, row, vector, matrix, &timeElapsed);
    for (i = 0; i < row; ++i) {
        printf("%f\n", result[i]);
    }
    printf("\nTime elapsed: %e\n\n", timeElapsed);

    free(result);
    free(vector);
    for (i = 0; i < column; ++i) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}

/**
 * Effettua il calcolo matrice per vettore e ritorna il risultato
 * @param column il numero di colonne della matrice
 * @param row il numero di righe della matrice
 * @param vector il vettore di dimensione column
 * @param matrix la matrice di dimensione column*row
 * @param timeElapsed valore di ritorno che rappresenta il tempo impiegato
 * @return il vettore risultato di dimensione column
 */
double* vectorXMatrix(int column, int row, const double* vector, double** matrix, double* timeElapsed) {
    int i, j;
    double* result = malloc(sizeof(double) * column);
    struct timeval timeVal;
    double startTime, endTime;

    // Tempo iniziale
    gettimeofday(&timeVal, NULL);
    startTime = timeVal.tv_sec + (timeVal.tv_usec/1000000.0);

    // Effettua il calcolo matrice x vettore utilizzando openmp
    #pragma omp parallel for default(none) shared(column, row, vector, matrix, result, timeElapsed) private(i, j)
    for (i = 0; i < row; i++) {
        for (j = 0; j < column; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    // Tempo finale
    gettimeofday(&timeVal, NULL);
    endTime = timeVal.tv_sec + (timeVal.tv_usec/1000000.0);

    // Tempo totale ottenuto dalla differenza finale - iniziale
    *timeElapsed = endTime - startTime;

    return result;
}

/**
 * Converte una stringa in input in int
 * @param str la stringa da convertire
 * @param val dove viene salvato il risultato della conversione
 * @return true se la conversione termina con successo, falso altrimenti
 */
bool parseInt(char* str, int* val) {
    char *temp;
    bool result = true;
    errno = 0;
    long ret = strtol(str, &temp, 0);

    if (temp == str || *temp != '\0' || ((ret == LONG_MIN || ret == LONG_MAX) && errno == ERANGE)) {
        result = false;
    }

    *val = (int) ret;
    return result;
}

/**
 * Ritorna un numero casuale nel range definito
 * @param min il numero minimo, incluso
 * @param max il numero massimo, incluso
 * @return float il numero casuale
 */
double getRandomDoubleNumberInRange(int min, int max) {
    return (double) min + rand() / (double) RAND_MAX * max - min;
}

/**
 * Alloca una matrice di dimensione column*row con valori casuali
 * @param column il numero di colonne
 * @param row il numero di righe
 * @return la matrice allocata e riempita casualmente
 */
double** getMatrixOfRandomNumbersOfSize(int column, int row) {
    int i, j;
    double** matrix = malloc(sizeof(double*) * column);
    for (i = 0; i < column; ++i) {
        matrix[i] = malloc(sizeof(double*) * row);
        for (j = 0; j < row; ++j) {
            matrix[i][j] = getRandomDoubleNumberInRange(MIN_RANDOM, MAX_RANDOM);
        }
    }
    return matrix;
}

/**
 * Alloca un vettore di dimensione size con valori casuali
 * @param size la dimensione del vettore
 * @return il vettore allocato e riempito casualmente
 */
double* getVectorOfRandomNumbersOfSize(int size) {
    int i;
    double* vector = malloc(sizeof(double) * size);
    for (i = 0; i < size; ++i) {
        vector[i] = getRandomDoubleNumberInRange(MIN_RANDOM, MAX_RANDOM);
    }
    return vector;
}
