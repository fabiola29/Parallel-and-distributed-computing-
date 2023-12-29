#include "utils.h"

bool parseInt(char* str, int* val) {
    char *temp;
    bool result = true;
    errno = 0;
    long ret = strtol(str, &temp, 0);

    if (temp == str || *temp != '\0' || ((ret == LONG_MIN || ret == LONG_MAX) && errno == ERANGE))
        result = false;

    *val = (int) ret;
    return result;
}

double getRandomDoubleNumberInRange(int min, int max) {
    return (double) min + rand() / (double) RAND_MAX * max - min;
}

bool isPerfectSquare(int number) {
    int s = sqrt(number);
    return (s * s) == number;
}

int mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

double** getMatrixOfRandomNumbersOfSize(int column, int row, int min, int max) {
    int i, j;
    double** matrix = malloc(sizeof(double*) * column);
    for (i = 0; i < column; ++i) {
        matrix[i] = malloc(sizeof(double*) * row);
        for (j = 0; j < row; ++j) {
            matrix[i][j] = getRandomDoubleNumberInRange(min, max);
        }
    }
    return matrix;
}

void printSquareMatrix(double** matrix, int size) {
    int i, j;

    for(i = 0; i < size; ++i) {
        for(j = 0; j < size; ++j) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}


