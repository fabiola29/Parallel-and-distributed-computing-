#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


// Dichiarazione Funzioni
int* matxvet(int,int,int*,int**);
void check_input(int,char**);

int main(int argc, char** argv){
    int **A, *x, *ris;
    int n, m, i, j;

    check_input(argc,argv); //Chiamata a funzione che verifica se i dati passati in input rispettano le restrizioni indicate nell'intestazione
    
    n = atoi(argv[1]); //Prendo in input il numero di righe della matrice
    m = atoi(argv[2]); //Prendo in input il numero di colonne della matrice e la dimensione del vettore

    //------ INIZIO ALLOCAZIONI

    A = (int**) malloc (n*sizeof(int*)); //Alloco le righe della matrice
    for(i=0;i<n;i++)
        A[i] = (int*) malloc (m*sizeof(int)); //Alloco le colonne della matrice
    
    x = (int*) malloc (m*sizeof(int)); //Alloco il vettore

    //------ FINE ALLOCAZIONI


    //------ INIZIO INIZIALIZZAZIONE MATRICE E VETTORE

    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
            A[i][j] = (i*m) + j + 1; //Inizializzo la matrice da 1 a NxM

    for(i=0; i<m; i++)
        x[i] = i+1; //Inizializzo il vettore da 1 a M

    //------ FINE INIZIALIZZAZIONE MATRICE E VETTORE

    ris = matxvet(n,m,x,A); //Chiamo la funzione che effettua il prodotto matrice-vettore Ax

    for(i=0;i<n;i++)
        printf("x[%d]: %d\n",i,ris[i]); //Stampo il risultato

    return 0;
}


/**
 * Calcola il prodotto matrice-vettore in modo parallelo utilizzando OpenMP.
 *
 * @param n La dimensione delle righe della matrice.
 * @param m La dimensione delle colonne della matrice e del vettore.
 * @param x Il vettore di input.
 * @param A La matrice di input.
 * @return int* Il vettore risultato del prodotto matrice-vettore.
 *
 * Questa funzione prende in input la dimensione delle righe (n), delle colonne (m), un vettore (x)
 * e una matrice (A). Alloca dinamicamente un vettore risultato (ris) e inizializza i suoi elementi a 0.
 * Utilizza la libreria OpenMP per parallelizzare il calcolo del prodotto matrice-vettore, distribuendo
 * il lavoro tra i thread disponibili. Infine, restituisce il vettore risultato.
 *
 * La funzione stampa anche il tempo impiegato per l'esecuzione del calcolo in modalità parallela.
 */
int* matxvet(int n, int m, int *x, int **A){
    int i,j;
    int *ris;
    ris = (int*) calloc (n,sizeof(int)); //Alloco il vettore risultato e inizializzo i suoi elementi a 0
    struct timeval time;
    double inizio, fine;

    gettimeofday(&time, NULL);
    inizio=time.tv_sec+(time.tv_usec/1000000.0); //Prendo il tempo iniziale

    //------- INIZIO REGIONE PARALLELA

    #pragma omp parallel for default(none) shared(A,x,m,n,ris) private(i,j)
        for(i=0; i<n; i++)
            for(j=0; j<m; j++)
                ris[i] += A[i][j]*x[j]; //Effettuo il prodotto matrice-vettore in parallelo

    //------- FINE REGIONE PARALLELA

    gettimeofday(&time, NULL); 
    fine=time.tv_sec+(time.tv_usec/1000000.0); //Prendo il tempo finale
    printf("tempo impiegato: %f\n", fine-inizio); //Stampo il tempo di esecuzione del prodotto matrice-vettore

    return ris; //Ritorno il vettore risultato
}

/**
 * Verifica la correttezza degli argomenti forniti da riga di comando.
 * 
 * @param argc Il numero totale di argomenti passati da riga di comando.
 * @param argv Un array di stringhe contenente gli argomenti passati da riga di comando.
 * 
 * Questa funzione accetta il numero totale di argomenti (argc) e un array di stringhe (argv)
 * e controlla che siano presenti esattamente tre argomenti. Inoltre, verifica che i valori numerici
 * associati al numero di righe e colonne siano entrambi maggiori di zero. In caso di violazioni di
 * queste condizioni, stampa messaggi di errore appropriati e termina il programma con un codice di errore.
 *
 */
void check_input(int argc, char **argv){
    if(argc!=3){
        printf("Inserire i 2 argomenti richiesti: \n1. Numero di righe della matrice\n2. Numero di colonne della matrice (e quindi dimensione   vettore)\n");
        exit(1);
    }
    if(atoi(argv[1]) <= 0){
        printf("Il numero di righe della matrice deve essere un valore maggiore di 0");
        exit(2);
    }
    if(atoi(argv[2]) <= 0){
        printf("Il numero di colonne della matrice (e quindi dimensione vettore) deve essere un valore maggiore di 0");
        exit(3);
    }
}

