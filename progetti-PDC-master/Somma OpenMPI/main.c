#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define MAX_RANDOM_NUMBER 1000
#define TAG_START 100

float distributeNumbersAndGetPartialSum(char** argv, double* startTime, int processorId, int numberOfProcessors);

// Strategie
float strategyOne(float sum, int processorId, int numberOfProcessors, int masterId);
float strategyTwo(float sum, int processorId, int numberOfProcessors, int masterId);
float strategyThree(float sum, int processorId, int numberOfProcessors);

// Funzioni accessorie
float getRandomFloatNumberInRange(int min, int max);
bool parseInt(char* str, int* val);
bool parseFloat(char* str, float* val);
bool isPowerOfTwo(unsigned long x);

/**
 * Gli argomenti vengono passati nella forma: nInput valori[nInput] strategia masterId
 * @param nInput numero di valori da sommare
 * @param numbers[nInput] se nInput<=20 sono i numeri da sommare, altrimenti è ignorato
 * @param strategy la strategia da utilizzare, valore compreso tra 1 e 3
 * @param masterId il processore che effettua i calcoli e stampa la somma, se strategy=3 e -1 stampano tutti
 */
int main(int argc, char** argv) {
    int processorId, numberOfProcessors;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &processorId);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcessors);

    /* I tempi dei singoli processori, ottenuti tramite differenza, mentre il tempo totale come il
     * massimo valore tra i processori */
    double startTime, endTime, processorTime, totalTime = 0.0;

    /* Prova ad convertire da argv il valore di strategy e masterId */
    int masterId, strategy;
    if(!parseInt(argv[argc - 2], &strategy) || !parseInt(argv[argc - 1], &masterId)) {
        fprintf(stderr, "ERROR - Cannot parse the masterId or strategy with values provided.\n"
                        "strategy:%s\nprinter:%s\n\n", argv[argc - 2], argv[argc - 1]);
        return 1;
    }

    if(strategy < 0 || strategy > 3) {
        fprintf(stderr, "ERROR - Strategy number not allowed, must be a value between 1 and 3.\nstrategy:%d\n", strategy);
        return 1;
    }

    if(masterId < -1 || masterId > numberOfProcessors) {
        printf("WARNING - inserted masterId is not valid, the processor 0 will print instead.\n");
        masterId = 0;
    } else if(masterId == -1 && (strategy == 1 || strategy == 2)) {
        printf("WARNING - Print of all processor with strategy 1 or 2 is not available, the processor 0 will print the result instead.\n");
        masterId = 0;
    }

    /* Invoca la distribuzione dei numeri, assegna lo startTime ed ottiene la propria somma parziale */
    float sum = distributeNumbersAndGetPartialSum(argv, &startTime, processorId, numberOfProcessors);

    /* Per le strategie 2 e 3 è un requisito essenziale che il numero di processori sia una potenza di due */
    if(strategy == 3 && isPowerOfTwo(numberOfProcessors)) {
        sum = strategyThree(sum, processorId, numberOfProcessors);
    } else if(strategy == 2 && isPowerOfTwo(numberOfProcessors)) {
        sum = strategyTwo(sum, processorId, numberOfProcessors, masterId);
    } else {
        sum = strategyOne(sum, processorId, numberOfProcessors, masterId);
    }

    /* Calcola il tempo finale di esecuzione dalla somma parziale */
    endTime = MPI_Wtime();
    processorTime = endTime - startTime;

    /* Passa al masterId il tempo maggiore impiegato */
    int sendTo = masterId == -1 ? 0 : masterId;
    MPI_Reduce(&processorTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, sendTo, MPI_COMM_WORLD);

    if(masterId == processorId || masterId == -1) {
        printf("[Completed] Processor ID %d, total sum: %f\n", processorId, sum);
        if(totalTime > 0.0) {
            printf("Time elapsed: %e seconds\n", totalTime);
        }
    }

    MPI_Finalize();
    return 0;
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

    if (temp == str || *temp != '\0' || ((ret == LONG_MIN || ret == LONG_MAX) && errno == ERANGE))
        result = false;

    *val = (int) ret;
    return result;
}

/**
 * Converte una stringa in input in float
 * @param str la stringa da convertire
 * @param val dove viene salvato il risultato della conversione
 * @return true se la conversione termina con successo, falso altrimenti
 */
bool parseFloat(char* str, float* val) {
    *val = atof(str);
    return true;
}

/**
 * Verifica se il numero fornito è una potenza di due
 * @param x il valore da verificare
 * @return true se è una potenza di 2, false altrimenti
 */
bool isPowerOfTwo(unsigned long x) {
    return (x & (x - 1)) == 0;
}

/**
 * Ritorna un numero casuale nel range definito
 * @param min il numero minimo, incluso
 * @param max il numero massimo, incluso
 * @return float il numero casuale
 */
float getRandomFloatNumberInRange(int min, int max) {
    return (float) min + rand() / (float) RAND_MAX * max - min;
}

/**
 * Distribuisce o riceve i numeri tra i vari processori e ritorna la propria somma parziale
 * @param argv gli argomenti del programma, da cui leggere il numero di input ed i valori
 * @param startTime valore di ritorno che rappresenta l'inizio delle operazioni del processore
 * @param processorId l'id del processore corrente, per discriminare se deve inviare o ricevere
 * @param numberOfProcessors il numero totale di processori
 */
float distributeNumbersAndGetPartialSum(char** argv, double* startTime, int processorId, int numberOfProcessors) {
    int tmp, sentNumbers, tag, i;
    float sum = 0;

    /* Prova a leggere da argv il numero di valori da sommare */
    int inputSize;
    if(!parseInt(argv[1], &inputSize) || inputSize < 0) {
        fprintf(stderr, "ERROR - Input size not allowed, must be a value between 0 and INT_MAX\ninputSize:%s\n", argv[1]);
        return 1;
    }

    // Numero dei valori da distribuire per ogni processore
    int amountOfNumbers = inputSize / numberOfProcessors;

    // Numeri non multipli del numero di processori da ripartire
    int rest = inputSize % numberOfProcessors;
    float* numbers = NULL;

    // Alcuni processori sommeranno più di un numero
    if(processorId < rest) {
        amountOfNumbers++;
    }

    // Il processore 0 distribuisce i valori
    if(processorId == 0) {
        numbers = malloc(sizeof(float) * inputSize);

        if(inputSize <= 20) {
            // Legge i valori da riga di comando
            int k = 2;
            for (i = 0; i < inputSize; i++) {
                if(!parseFloat(argv[k++], &numbers[i])) {
                    fprintf(stderr, "ERROR - Cannot parse a value to sum.\nvalue[%d]: %s\n", i, argv[k-1]);
                    return 1;
                }
            }
        } else {
            // Genera valori casuali
            srand(time(NULL));
            for (i = 0; i <= inputSize; i++) {
                numbers[i] = getRandomFloatNumberInRange(0, MAX_RANDOM_NUMBER);
            }
        }

        tmp = amountOfNumbers;
        sentNumbers = 0;

        // Invia i valori ai processori
        for (i = 1; i < numberOfProcessors; ++i) {
            sentNumbers += tmp;
            tag = TAG_START + i;

            if (i == rest) {
                tmp--;
            }

            MPI_Send(&numbers[sentNumbers], tmp, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
        }
    } else { // Riceve i valori
        numbers = malloc(sizeof(float) * amountOfNumbers);
        MPI_Status status;

        tag = TAG_START + processorId;
        MPI_Recv(numbers, amountOfNumbers, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
    }

    // Si sincronizza con gli altri processori
    MPI_Barrier(MPI_COMM_WORLD);
    *startTime = MPI_Wtime();

    // Effettua la propria somma parziale
    for (i = 0; i < amountOfNumbers; ++i) {
        sum += numbers[i];
    }

    free(numbers);
    return sum;
}

/**
 * Applica la strategia uno per la somma, il processore master riceve i valori dagli altri processori
 * @param sum la somma corrente del processore
 * @param processorId l'id del processore
 * @param numberOfProcessors il numero di processori, deve essere una potenza di due
 * @param masterId l'id del processo incaricato di fare i calcoli
 */
float strategyOne(float sum, int processorId, int numberOfProcessors, int masterId) {
    float partialSum;
    int tag, i;
    MPI_Status status;

    if(processorId == masterId) {
        for (i = 0; i < numberOfProcessors; i++) {
            if(masterId == i) {
                continue;
            }

            tag = TAG_START + i;
            MPI_Recv(&partialSum, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
            sum += partialSum;
        }
    } else {
        tag = TAG_START + processorId;
        MPI_Send(&sum, 1, MPI_FLOAT, masterId, tag, MPI_COMM_WORLD);
    }

    return sum;
}

/**
 * Applica la strategia due per la somma, questa utilizza un albero per la risoluzione, quindi un pre-requisito
 * del metodo è che il numero di processori sia una potenza di due
 * @param sum la somma corrente del processore
 * @param processorId l'id del processore
 * @param numberOfProcessors il numero di processori, deve essere una potenza di due
 * @param masterId l'id del processo incaricato di fare i calcoli
 */
float strategyTwo(float sum, int processorId, int numberOfProcessors, int masterId) {
    float partialSum = 0;
    int tag, i;
    MPI_Status status;

    int logicId = processorId + (numberOfProcessors - masterId);
    logicId = logicId % numberOfProcessors;

    for (i = 0; i < (int) log2(numberOfProcessors); ++i) {
        tag = TAG_START + i;
        if((logicId % (int) pow(2, i)) == 0) {
            if((logicId % (int) pow(2, i + 1)) == 0) {
                int senderId = ((int) (processorId + pow(2, i)) % numberOfProcessors);
                MPI_Recv(&partialSum, 1, MPI_FLOAT, senderId, tag, MPI_COMM_WORLD, &status);
                sum += partialSum;
            } else {
                int receiverId = processorId - (int) pow(2, i);
                if(receiverId < 0) {
                    receiverId += numberOfProcessors;
                }
                MPI_Send(&sum, 1, MPI_FLOAT, receiverId, tag, MPI_COMM_WORLD);
            }
        }
    }
    return sum;
}

/**
 * Applica la strategia tre per la somma, questa utilizza un albero per la risoluzione, quindi un pre-requisito
 * del metodo è che il numero di processori sia una potenza di due
 * @param sum la somma corrente del processore
 * @param processorId l'id del processore
 * @param numberOfProcessors il numero di processori, deve essere una potenza di due
 */
float strategyThree(float sum, int processorId, int numberOfProcessors) {
    float partialSum = 0;
    int tag, i;
    MPI_Status status;

    for (i = 0; i < (int) log2(numberOfProcessors); ++i) {
        tag = TAG_START + i;
        if((processorId % (int) pow(2, i + 1)) < (int) pow(2, i)) {
            int otherId = processorId + (int) pow(2, i);
            MPI_Send(&sum, 1, MPI_FLOAT, otherId, tag, MPI_COMM_WORLD);
            MPI_Recv(&partialSum, 1, MPI_FLOAT, otherId, tag, MPI_COMM_WORLD, &status);
        } else {
            int otherId = processorId - (int) pow(2, i);
            MPI_Send(&sum, 1, MPI_FLOAT, otherId, tag, MPI_COMM_WORLD);
            MPI_Recv(&partialSum, 1, MPI_FLOAT, otherId, tag, MPI_COMM_WORLD, &status);
        }
        sum += partialSum;
    }

    return sum;
}
