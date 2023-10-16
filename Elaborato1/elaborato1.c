#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <time.h>

/*Algoritmo per il calcolo della somma di N numeri reali, 
in ambiente di calcolo parallelo su architettura MIMD a memoria distribuita*/

/*----------------------------------------------------------
I parametri in input vanno inseriti in questo ordine:
   1. ID Processo che deve stampare il risultato
   2. Numero strategia da utilizzare
   3. Numero di valori totali in input N
   4. Se N<=20 inserire i valori da sommare, altrimenti inserire solo i 3 parametri precedenti
-----------------------------------------------------------*/


void check_input(int argc, char *argv[], int nproc);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int menum, nproc, n, nloc, tag, i, rest, strategia, id, start, tmp, id_tot_sum;
    double *x;
    double *xloc;
    double t0, t1, timep, timetot, sum_parz, sum = 0.0;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &menum); // Prendo il rank
    MPI_Comm_size(MPI_COMM_WORLD, &nproc); // Prendo il numero di processi

    // --- INIZIO DISTRIBUZIONE INPUT
    if (menum == 0) {

        check_input(argc, argv, nproc); // Verifico che i parametri rispettino le restrizioni indicate. In caso contrario MPI_Abort().

        n = atoi(argv[3]); // Assegno il numero di valori alla variabile n
        x = (double *)malloc(n * sizeof(double)); // Alloco l'array x

        // --- ASSEGNO VALORI IN ARRAY x A SECONDA DEL NUMERO DI DATI CONSIDERATI
        if (n > 20) {
            srand(time(NULL));
            for (i = 0; i < n; i++)
                x[i] = ((double)rand() * 10.0) / (double)RAND_MAX; // Genero un double random compreso tra 0 e 10
        } else {
            for (i = 0; i < n; i++)
                x[i] = atof(argv[i + 4]); // Assegno i valori passati come argomento ad x
        }

        id = atoi(argv[1]); // Assegno rank processo che dovrà stampare il risultato
        strategia = atoi(argv[2]); // Assegno la strategia da utilizzare

        if (strategia != 1 && (nproc & (nproc - 1)) != 0) // SE LA STRATEGIA E' DIVERSA DA 1 E "nproc" NON E' UNA POTENZA DI 2 ALLORA FORZO STRATEGIA A 1
            strategia = 1;
        printf("STRATEGIA UTILIZZATA PER LA SOMMA: %d\n", strategia); // STAMPO STRATEGIA UTILIZZATA (STAMPA SOLO PROCESSO RANK 0)

        // --- FINE ASSEGNAZIONE

        for (i = 1; i < nproc; i++) {
            tag = 10 + i;
            MPI_Send(&n, 1, MPI_INT, i, tag, MPI_COMM_WORLD); // INVIO VALORE n AI PROCESSI DIVERSI DA 0
            tag = 11 + i;
            MPI_Send(&strategia, 1, MPI_INT, i, tag, MPI_COMM_WORLD); // INVIO STRATEGIA AI PROCESSI DIVERSI DA 0
            tag = 12 + i;
            MPI_Send(&id, 1, MPI_INT, i, tag, MPI_COMM_WORLD); // INVIO ID PROCESSO CHE DEVE STAMPARE AI PROCESSI DIVERSI DA 0
        }
    } else {
        tag = 10 + menum;
        MPI_Recv(&n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status); // RICEZIONE VALORE n DA PROCESSO 0
        tag = 11 + menum;
        MPI_Recv(&strategia, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status); // RICEZIONE STRATEGIA DA PROCESSO 0
        tag = 12 + menum;
        MPI_Recv(&id, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status); // RICEZIONE ID PROCESSO STAMPA DA PROCESSO 0
    }
    // --- INIZIO DISTRIBUZIONE INPUT

    nloc = n / nproc; // Calcolo n locale
    rest = n % nproc; // calcolo eventuale resto divisione
    if (menum < rest) nloc = nloc + 1; // fix n locale
    xloc = (double *)malloc(nloc * sizeof(double)); // Alloco array xloc per quanti valori il processo i-esimo ne ha bisogno (nloc)

    // --- INIZIO DISTRIBUZIONE VALORI DA SOMMARE AI VARI PROCESSI
    if (menum == 0) {
        xloc = x;
        tmp = nloc;
        start = 0;
        for (i = 1; i < nproc; i++) {
            start = start + tmp;
            tag = 22 + i;
            if (i == rest) tmp = tmp - 1;
            MPI_Send(&x[start], tmp, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }
    } else {
        tag = 22 + menum;
        MPI_Recv(xloc, nloc, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }
    // --- FINE DISTRIBUZIONE VALORI DA SOMMARE

    MPI_Barrier(MPI_COMM_WORLD); // BARRIERA PER CALCOLARE POI I TEMPI
    t0 = MPI_Wtime(); // Prendo il primo tempo

    // --- INIZIO CALCOLO SOMMA LOCALE
    for (i = 0; i < nloc; i++) {
        sum = sum + xloc[i];
    }
    // --- FINE CALCOLO SOMMA LOCALE

    // --- INIZIO PRIMA STRATEGIA
    if (strategia == 1) {
        id_tot_sum = id;
        if (id_tot_sum == -1) id_tot_sum = 0; // Se devono stampare tutti, rank 0 avrà la somma totale
        if (menum == id_tot_sum) {
            for (i = 0; i < nproc; i++) {
                if (menum != i) { // ricevo solo dai processi diversi da me stesso
                    tag = 80 + i;
                    MPI_Recv(&sum_parz, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status); // Ricevo dal processo i
                    sum = sum + sum_parz; // Aggiungo la somma parziale inviatami alla mia somma parziale
                }
            }
        } else {
            tag = 80 + menum;
            MPI_Send(&sum, 1, MPI_DOUBLE, id_tot_sum, tag, MPI_COMM_WORLD); // Invio sum al processo che dovrà stampare alla fine (se -1 invio a 0)
        }
    }
    // --- FINE PRIMA STRATEGIA
    // --- INIZIO SECONDA STRATEGIA (STRATEGIA COSTRUITA IN MANIERA TALE CHE L'ALBERO RUOTI NEL CASO IN CUI DEBBANO STAMPARE
    //                              PROCESSI DIVERSI DA 0)
    else if (strategia == 2) {
        int recv, send;
        id_tot_sum = id;
        if (id_tot_sum == -1) id_tot_sum = 0; // Se devono stampare tutti, rank 0 avrà la somma totale
        for (i = 0; i < (log10(nproc) / log10(2)); i++) {
            if (fmod((menum - id_tot_sum + nproc) % nproc, (pow(2, i))) == 0) { // chi partecipa (rotazione data dal primo membro di fmod)
                if (fmod((menum - id_tot_sum + nproc) % nproc, (pow(2, i + 1))) == 0) { // decido chi deve ricevere e chi inviare
                    tag = menum;
                    recv = (int)(menum + pow(2, i)) % nproc; // da chi ricevo la somma parziale
                    MPI_Recv(&sum_parz, 1, MPI_DOUBLE, recv, tag, MPI_COMM_WORLD, &status);
                    sum = sum + sum_parz;
                } else {
                    send = (((int)(menum - pow(2, i)) % nproc) + nproc) % nproc; // a chi invio la somma parziale (formula modulo non negativo)
                    tag = send;
                    MPI_Send(&sum, 1, MPI_DOUBLE, send, tag, MPI_COMM_WORLD);
                }
            }
        }
    }
    // --- FINE SECONDA STRATEGIA
    // --- INIZIO TERZA STRATEGIA. Send e Recv sono invertite negli if per evitare Dead Lock.
    // I processi nell'else non fanno nessuna somma ma aggiornano solo il proprio "sum" col valore prelevato dalla Recv
    else if (strategia == 3) {
        for (i = 0; i < (log10(nproc) / log10(2)); i++) {
            if (fmod(menum, (pow(2, i + 1))) < pow(2, i)) {
                tag = menum;
                MPI_Recv(&sum_parz, 1, MPI_DOUBLE, menum + pow(2, i), tag, MPI_COMM_WORLD, &status);
                tag = menum + pow(2, i);
                MPI_Send(&sum, 1, MPI_DOUBLE, menum + pow(2, i), tag, MPI_COMM_WORLD);
            } else {
                tag = menum - pow(2, i);
                MPI_Send(&sum, 1, MPI_DOUBLE, menum - pow(2, i), tag, MPI_COMM_WORLD);
                tag = menum;
                MPI_Recv(&sum_parz, 1, MPI_DOUBLE, menum - pow(2, i), tag, MPI_COMM_WORLD, &status);
            }
            sum = sum + sum_parz; // Ogni processo fa la somma col valore sum_parz che ha ricevuto dall'apposito processo
        }
    }
    // --- FINE TERZA STRATEGIA

    t1 = MPI_Wtime(); // Prendo il secondo tempo
    timep = t1 - t0; // Calcolo la differenza per capire il tempo impiegato dal singolo processo
    MPI_Reduce(&timep, &timetot, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // Calcolo il tempo massimo e lo pongo in "timetot" al processo 0

    // --- INIZIO STAMPA RISULTATI
    if (id == -1) {
        if (strategia != 3)
            MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Invio in broadcast da 0 la somma a tutti (quando id = -1 la somma totale la ha solo rank 0 per convenzione)

        printf("Somma calcolata da %d è %lf\n", menum, sum); // stampano tutti i processi la somma totale
    } else if (menum == id) {
        printf("Somma calcolata da %d è %lf\n", menum, sum); // stampa solo il processo "id" che avrà la somma totale
    }
    // --- FINE STAMPA RISULTATI

    if (menum == 0) {
        printf("Tempo massimo %lf\n", timetot); // Stampo tempo massimo di somma + strategia adottata
    }

    MPI_Finalize();
    return 0;
}

void check_input(int argc, char *argv[], int nproc) {
    if (argc < 4) {
        printf("INSERIRE ALMENO 3 ARGOMENTI:\n1. ID processo stampa\n2. Strategia da utilizzare\n3. Numero di dati in input\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (atoi(argv[1]) >= nproc) {
        printf("NUMERO PROCESSO NON VALIDO!!! INSERIRE UN VALORE COMPRESO TRA -1 e %d\n", nproc - 1);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    if (atoi(argv[2]) < 1 || atoi(argv[2]) > 3) {
        printf("NUMERO STRATEGIA NON VALIDO!!! INSERIRE VALORE COMPRESO TRA 1 E 3\n");
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    if (atoi(argv[3]) <= 0) {
        printf("NUMERO DATI INPUT NON VALIDO!!! INSERIRE UN VALORE NON NEGATIVO\n");
        MPI_Abort(MPI_COMM_WORLD, 4);
    }
    if (atoi(argv[3]) > 20 && argc > 4) {
        printf("INSERITI DATI IN INPUT CON VALORE N MAGGIORE DI 20. NON INSERIRE NESSUN DATO IN INPUT\n");
        MPI_Abort(MPI_COMM_WORLD, 5);
    }
    if ((atoi(argv[3]) > 0 && atoi(argv[3]) <= 20) && argc - 4 != atoi(argv[3])) {
        printf("NUMERO DI DATI IN INPUT NON CORRETTO CON IL VALORE INDICATO N\n");
        MPI_Abort(MPI_COMM_WORLD, 6);
    }
}
