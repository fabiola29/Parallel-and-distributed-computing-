#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>

#ifndef PDC_4_JANUARY2_UTILS_H
#define PDC_4_JANUARY2_UTILS_H

/* **************************************************************** *
 * Allocazioni di matrici e vettori casuali con generatore casuale  *
 * **************************************************************** */

/**
 * Ritorna un numero casuale nel range definito
 * @param min il numero minimo, incluso
 * @param max il numero massimo, incluso
 * @return float il numero casuale
 */
double getRandomDoubleNumberInRange(int min, int max);

/**
 * Alloca una matrice di dimensione column*row con valori casuali
 * @param column il numero di colonne
 * @param row il numero di righe
 * @return la matrice allocata e riempita casualmente
 */
double** getMatrixOfRandomNumbersOfSize(int column, int row, int min, int max);

/**
 * Stampa la matrice in ingresso in standard output
 * @param matrix la matrice da stampare
 */
void printSquareMatrix(double** matrix, int dimension);



/* **************************************************************** *
 * Funzioni per il parsing di argomenti ottenuti da riga di comando *
 * **************************************************************** */

/**
 * Converte una stringa in input in int
 * @param str la stringa da convertire
 * @param val dove viene salvato il risultato della conversione
 * @return true se la conversione termina con successo, falso altrimenti
 */
bool parseInt(char* arg, int* output);



/* **************************************************************** *
 * Funzioni matematiche                                             *
 * **************************************************************** */

/**
 * Ritorna se il numero in ingresso ha una radice quadrata intera
 * @param number il valore da verificare
 * @return true se ha radice intera, false altrimenti
 */
bool isPerfectSquare(int number);

/**
 * Ritorna il modulo b di a
 * @param a il valore da cui ottenere il modulo b
 * @param b il modulo da utilizzare
 * @return a modulo b
 */
int mod(int a, int b);

#endif //PDC_4_JANUARY2_UTILS_H
