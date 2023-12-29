#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "utils.h"

#define MIN_RANDOM_NUMBER 0
#define MAX_RANDOM_NUMBER 15
#define MASTER_ID 0

#define TAG_MATRIX_FIRST 1
#define TAG_MATRIX_SECOND 1

void distributeMatrixes(int matrixSize, int gridSize, MPI_Comm mpi_comm);
double** receiveMatrix(int matrixSize, int tag, MPI_Comm mpi_comm);
double** getPartialMatrixResult(int numberOfValues, int gridSize, double* startTime, MPI_Comm mpi_comm_grid);
void printResult(int matrixSize, int gridSize, double** partialResult, MPI_Comm mpi_comm);
int validateInput(int numberOfProcessors, int *matrixSize, int argc, char** argv, int processId);

/**
 * Gli argomenti vengono passati nella forma: matrixSize
 * @param matrixSize la dimensione della matrice MxM
 */
int main(int argc, char** argv) {
    int matrixSize, processId, numberOfProcessors;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcessors);

    /* I tempi dei singoli processori, ottenuti tramite differenza, mentre il tempo totale come il
     * massimo valore tra i processori */
    double startTime, endTime, processorTime, totalTime = 0.0;

     // ********************************************
    // ****** Controllo input
    // ********************************************

    // Validazione dei parametri di input
    if (validateInput(numberOfProcessors, &matrixSize, argc, argv, processId)) {
        // Errore nei parametri di input, esce dal programma
        MPI_Finalize();
        return 1;
    }

    // ********************************************
    // ****** Creazione della griglia di processori
    // ********************************************

    int gridSize;

    MPI_Bcast(&gridSize, 1, MPI_DOUBLE, MASTER_ID, MPI_COMM_WORLD);

    gridSize = (int) sqrt(numberOfProcessors);
    int numberOfValues = matrixSize/gridSize;

    // Siccome è una matrice quadrata le due dimensioni sono uguali
    int* sizesOfGrid = calloc(2, sizeof(int));
    sizesOfGrid[0] = sizesOfGrid[1] = gridSize;

    int* gridPeriod = calloc(2, sizeof(int));
    gridPeriod[0] = gridPeriod[1] = 0;

    MPI_Comm mpi_comm_grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, sizesOfGrid, gridPeriod, 0, &mpi_comm_grid);
    MPI_Comm_rank(mpi_comm_grid, &processId);

    free(sizesOfGrid);
    free(gridPeriod);

    // ********************************************
    // ****** Calcolo prodotto matrice x matrice
    // ********************************************

    // Il processore con id MASTER_ID invia le matrici generate casualmente
    if(processId == MASTER_ID) {
        distributeMatrixes(matrixSize, gridSize, mpi_comm_grid);
    }

    // Effettuo il prodotto parziale della matrice
    double** partialResult = getPartialMatrixResult(numberOfValues, gridSize, &startTime, mpi_comm_grid);

    // Effettua la stampa del risultato finale
    printResult(matrixSize, gridSize, partialResult, mpi_comm_grid);

    // Salvo l'istante di tempo finale di esecuzione dalla somma parziale, faccio stampare al MASTER_ID
    endTime = MPI_Wtime();
    processorTime = endTime - startTime;

    // Passa al masterId il tempo maggiore impiegato
    MPI_Reduce(&processorTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, MASTER_ID, MPI_COMM_WORLD);

    // Il MASTER_ID fa la stampa del tempo impiegato
    if(processId == MASTER_ID) {
        printf("\nTime elapsed: %f seconds\n", totalTime);
    }

    // Libero la memoria dinamica utilizzata
    int i;
    for(i = 0; i < numberOfValues; ++i) {
        free(partialResult[i]);
    }
    free(partialResult);

    MPI_Finalize();
    return 0;
}

/**
 * Distribuisce le matrici generate casualmente tra tutti i processori, compreso il chiamante
 * @param matrixSize la dimensione della matrice da generare casualmente
 * @param gridSize la dimensione della matrice di processori
 * @param mpi_comm il Communicator, deve essere una griglia di dimensione gridSize
 */
void distributeMatrixes(int matrixSize, int gridSize, MPI_Comm mpi_comm) {
    int i, j, k, l;
    int numberOfValues = matrixSize/gridSize;

    // Inizializzazione del seed dei valori casuali
    srand(time(NULL));

    // Generazione e stampa delle due matrici casuali
    printf("First random matrix:\n");
    double** firstMatrix = getMatrixOfRandomNumbersOfSize(matrixSize, matrixSize, MIN_RANDOM_NUMBER, MAX_RANDOM_NUMBER);
    printSquareMatrix(firstMatrix, matrixSize);

    printf("\nSecond random matrix:\n");
    double** secondMatrix = getMatrixOfRandomNumbersOfSize(matrixSize, matrixSize, MIN_RANDOM_NUMBER, MAX_RANDOM_NUMBER);
    printSquareMatrix(secondMatrix, matrixSize);

    // Distribuzione di porzioni delle due matrici firstMatrix e secondMatrix a tutti i processori
    for(i = 0; i < gridSize; ++i) {
        for(j = 0; j < gridSize; ++j) {
            int receiverCoordinate[2] = { i, j };
            int receiverProcessor;
            MPI_Cart_rank(mpi_comm, receiverCoordinate, &receiverProcessor);

            for(k = (i * numberOfValues); k < (i * numberOfValues) + numberOfValues; ++k) {
                for(l = (j * numberOfValues); l < (j * numberOfValues) + numberOfValues; ++l) {
                    // Per l'invio vengono utilizzati due tag differenti per identificare i valori
                    MPI_Send(&firstMatrix[k][l], 1, MPI_DOUBLE, receiverProcessor, TAG_MATRIX_FIRST, mpi_comm);
                    MPI_Send(&secondMatrix[k][l], 1, MPI_DOUBLE, receiverProcessor, TAG_MATRIX_SECOND, mpi_comm);
                }
            }
        }
    }

    free(firstMatrix);
    free(secondMatrix);
}

/**
 * Riceve una matrice inviata con un determinato tag
 * @param matrixSize il numero M di righe e colonne della matrice MxM
 * @param tag il tag da cui ricevere
 * @param mpi_comm il Communicator
 * @return la matrice ricevuta
 */
double** receiveMatrix(int matrixSize, int tag, MPI_Comm mpi_comm) {
    int i, j;
    MPI_Status status;

    double** matrix = malloc(sizeof(double*) * matrixSize);
    for(i = 0; i < matrixSize; ++i) {
        matrix[i] = malloc(sizeof(double*) * matrixSize);

        for(j = 0; j < matrixSize; ++j) {
            MPI_Recv(&matrix[i][j], 1, MPI_DOUBLE, MASTER_ID, tag, mpi_comm, &status);
        }
    }

    return matrix;
}

/**
 * Effettua un prodotto parziale matrice x matrice
 * @param numberOfValues la dimensione parziale della matrice
 * @param gridSize la dimensione della griglia dei processori
 * @param startTime il valore in output dell'istante di inizio dei calcoli
 * @param mpi_comm_grid il Communicator, deve essere una griglia di dimensione gridSize
 * @return il prodotto parziale, una matrice di dimensione numberOfValues x numberOfValues
 */
double** getPartialMatrixResult(int numberOfValues, int gridSize, double* startTime, MPI_Comm mpi_comm_grid) {
    int i, j, k;
    MPI_Status status;

    int processId;
    MPI_Comm_rank(mpi_comm_grid, &processId);
    
    int* coordinates = calloc(2, sizeof(int));
    MPI_Cart_coords(mpi_comm_grid, processId, 2, coordinates);

    // Tutti i processori ricevono una porzione di matrice
    double** matrixOne = receiveMatrix(numberOfValues, TAG_MATRIX_FIRST, mpi_comm_grid);
    double** matrixTwo = receiveMatrix(numberOfValues, TAG_MATRIX_SECOND, mpi_comm_grid);

    // Alloco la matrice risultato parziale
    double** result = calloc(numberOfValues, sizeof(double*));
    for(i = 0; i < numberOfValues; ++i) {
        result[i] = calloc(numberOfValues, sizeof(double*));
    }

    // Mi sincronizzo con tutti i processori che hanno ricevuto le loro porzioni e salvo l'istante di tempo iniziale
    MPI_Barrier(MPI_COMM_WORLD);
    *startTime = MPI_Wtime();

    // Se è sulla diagonale invia alle varie righe i suoi valori, altrimenti riceve dalla diagonale
    int iteration;
    for(iteration = 0; iteration < gridSize; ++iteration) {

        double** receivedMatrix = NULL;
        if(coordinates[0] == mod(coordinates[1] - iteration, gridSize)) {

            int rowProcessor;
            for(rowProcessor = 0; rowProcessor < gridSize; ++rowProcessor) {
                int receiverCoordinate[2] = { coordinates[0], rowProcessor };
                int receiverId;
                MPI_Cart_rank(mpi_comm_grid, receiverCoordinate, &receiverId);

                // Non invia a se stesso
                if(processId != receiverId) {
                    //printf("Sono %d, devo inviare a %d che SULLA GRIGLIA sta a (%d, %d)\n\n", processId, receiverId, coordinates[0], rowProcessor);
                    for(i = 0; i < numberOfValues; ++i) {
                        for(j = 0; j < numberOfValues; ++j) {
                            MPI_Send(&matrixOne[i][j], 1, MPI_DOUBLE, receiverId, coordinates[0], mpi_comm_grid);
                        }
                    }
                }
            }
        } else {

            // Calcolo da chi devo ricevere
            int diagonalCoordinates[2] = {coordinates[0], mod(coordinates[0] + iteration, gridSize)};
            int diagonalProcessor;
            MPI_Cart_rank(mpi_comm_grid, diagonalCoordinates, &diagonalProcessor);

            receivedMatrix = malloc(sizeof(double*) * numberOfValues);

            for(i = 0; i < numberOfValues; ++i) {
                receivedMatrix[i] = malloc(sizeof(double*) * numberOfValues);

                for(j = 0; j < numberOfValues; ++j) {
                    MPI_Recv(&receivedMatrix[i][j], 1, MPI_DOUBLE, diagonalProcessor, coordinates[0], mpi_comm_grid, &status);
                }
            }
        }

        // Non è necessario eseguirlo la prima volta
        if(iteration > 0) {
            // Invio la mia matrice B al processore (i - 1, j), nella riga precedente
            int sendTo;
            int sendCoordinate[2] = { mod(coordinates[0] - 1, gridSize), coordinates[1] };
            MPI_Cart_rank(mpi_comm_grid, sendCoordinate, &sendTo);

            for(i = 0; i < numberOfValues; ++i) {
                for(j = 0; j < numberOfValues; ++j) {
                    MPI_Send(&matrixTwo[i][j], 1, MPI_DOUBLE, sendTo, coordinates[1], mpi_comm_grid);
                }
            }

            // Ricevo la matrice B dal processore (i + 1, j), nella riga successiva
            int receiveFrom;
            int receiveCoordinate[2] = { mod(coordinates[0] + 1, gridSize), coordinates[1] };
            MPI_Cart_rank(mpi_comm_grid, receiveCoordinate, &receiveFrom);

            for(i = 0; i < numberOfValues; ++i) {
                for(j = 0; j < numberOfValues; ++j) {
                    MPI_Recv(&matrixTwo[i][j], 1, MPI_DOUBLE, receiveFrom, coordinates[1], mpi_comm_grid, &status);
                }
            }
        }

        // Creo un puntatore alla matrice che voglio usare, per i processori diagonali è MatrixOne altrimenti quella ricevuta
        double** matrixToUse = (receivedMatrix == NULL) ? matrixOne : receivedMatrix;

        // Calcolo del prodotto parziale
        for(i = 0; i < numberOfValues; ++i) {
            for(j = 0; j < numberOfValues; ++j) {
                for(k = 0; k < numberOfValues; ++k) {
                    result[i][j] += (matrixToUse[i][k] * matrixTwo[k][j]);
                }
            }
        }

        // Pulisco le matrici che allo step successivo non servono più
        if(receivedMatrix != NULL) {
            for(i = 0; i < numberOfValues; ++i) {
                free(receivedMatrix[i]);
            }
            free(receivedMatrix);
        }
    }

    for(i = 0; i < numberOfValues; ++i) {
        free(matrixOne[i]);
        free(matrixTwo[i]);
    }
    free(matrixOne);
    free(matrixTwo);
    free(coordinates);

    return result;
}


/**
 * Verifica la correttezza dei parametri di input
 * @param numberOfProcessors il numero di processori
 * @param matrixSize la dimensione della matrice (output)
 * @param argc il numero di argomenti passati al programma
 * @param argv gli argomenti del programma
 * @param processId l'ID del processo
 * @return 0 se i parametri sono validi, 1 altrimenti
 */
int validateInput(int numberOfProcessors, int *matrixSize, int argc, char** argv, int processId) {
      // Check if there are enough arguments
    if (argc != 2) {
        if (processId == MASTER_ID) {
            fprintf(stdout, "Provide only one argument: the number of rows (and columns) of the matrices\n", argv[0]);
        }
        return 1;
    }


    // Tests whether the number of processors is a square root of integers
    if (!isPerfectSquare(numberOfProcessors)) {
        if (processId == MASTER_ID) {
            fprintf(stdout, "ERROR - Number of processors provided cannot be distribuited in a NxN grid.\n"
                            "number of processors: %d\n\n", numberOfProcessors);
        }
        return 1;
    }

    // Check the correctness of the size of the matrix passed as input
    if (!parseInt(argv[1], matrixSize)) {
        if (processId == MASTER_ID) {
            fprintf(stdout, "ERROR - Cannot parse the matrixSize with value provided.\n"
                            "matrixSize:%s\n\n", argv[1]);
        }
        return 1;
    }

    //Check that matrixSize is equal to numberOfProcessors or a multiple thereof
    if (*matrixSize < numberOfProcessors || *matrixSize % numberOfProcessors != 0) {
        if (processId == MASTER_ID) {
            fprintf(stdout, "ERROR - Matrix size must be >= of the grid and a multiple of number of processors.\n"
                            "matrixSize:%d\nnumberOfProcessors:%d\n", *matrixSize, numberOfProcessors);
        }
        return 1;
    }

    return 0;
}



/**
 * Stampa il risultato del prodotto matrice per matrice
 * @param matrixSize la dimensione originale della matrice
 * @param gridSize la dimensione della griglia di processori
 * @param partialResult il risultato parziale del chiamante
 * @param mpi_comm il Communicator, deve essere una griglia di dimensione gridSize
 */
void printResult(int matrixSize, int gridSize, double** partialResult, MPI_Comm mpi_comm_grid) {
    int i, j, k, l;
    MPI_Status status;
    int numberOfValues = matrixSize/gridSize;
    int processId;
    MPI_Comm_rank(mpi_comm_grid, &processId);

    // Tutti i processori inviano al MASTER_ID il proprio prodotto parziale
    for(i = 0; i < numberOfValues; i++) {
        for(j = 0; j < numberOfValues; j++) {
            MPI_Send(&partialResult[i][j], 1, MPI_DOUBLE, MASTER_ID, MASTER_ID, mpi_comm_grid);
        }
    }

    // Il MASTER_ID riceve tutta la matrice e la stampa ordinata
    if(processId == MASTER_ID) {

        double** resultMatrix = malloc(sizeof(double*) * matrixSize);
        for (i = 0; i < matrixSize; ++i) {
            resultMatrix[i] = malloc(sizeof(double*) * matrixSize);
        }

        for (i = 0; i < gridSize; ++i) {
            for (j = 0; j < gridSize; ++j) {
                int receiverCoordinate[2] = {i, j};
                int receiverProcessor;
                MPI_Cart_rank(mpi_comm_grid, receiverCoordinate, &receiverProcessor);

                for (k = (i * numberOfValues); k < (i * numberOfValues) + numberOfValues; ++k) {
                    for (l = (j * numberOfValues); l < (j * numberOfValues) + numberOfValues; ++l) {
                        MPI_Recv(&resultMatrix[k][l], 1, MPI_DOUBLE, receiverProcessor, MASTER_ID, mpi_comm_grid, &status);
                    }
                }
            }
        }

        printf("\nRESULT: \n");
        printSquareMatrix(resultMatrix, matrixSize);

        for (i = 0; i < matrixSize; ++i) {
            free(resultMatrix[i]);
        }
        free(resultMatrix);
    }

}
