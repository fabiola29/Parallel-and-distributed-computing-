#include <stdio.h>        
#include <stdlib.h>       
#include <string.h>       
#include <math.h>           
#include "mpi.h"
#include <sys/time.h>

// Dichiarazione funzioni
int     Controllo_input(int input[],int argc,char *argv[],int nproc);   
float** Crea_mat_random(int righe,int colonne);
float*  Crea_vet_random(int dim);
void    Ripartizione_mat_vet(MPI_Comm griglia,float** A,int m,int n,float subA[],int mloc,int nloc,int r,int c,float X[],float subX[],int coordinate[]);
void    Stampa_vet (float *X,int n);
void    Stampa_mat(float **A,int m,int n);
void    Crea_griglia(MPI_Comm *griglia, MPI_Comm *grigliar,MPI_Comm *grigliac, int menum, int nproc,int riga, int col, int *coordinate);
void    Prodotto( int mloc, int nloc, float *subA, float *subX, float *subB);

int main(int argc, char *argv[]) { 
    int menum,    //id processore nel communicator MPI_COMM_WORLD
      nproc,  //numero processori
      input[4],  //risultato controllo input,righe della griglia di calcolo , righe matrice, colonne matrice
      mloc,   //righe parziali assegnate ad ogni processore
      nloc,   //colonne parziali assegnate ad ogni processore 
      mloc0,  //righe parziali di P0
      nloc0,  //colonne parziali di P0
      r,     //resto della divisione  tra M e righe della griglia per capire se le righe della matrice sono ripartibili tra i processori equamente
      c,    //resto della divisione tra N e colonne della griglia per capire se le colonne della matrice sono ripartibili tra i processori equamente
      colonne_griglia,  //colonne della griglia di calcolo
      dim,  //dimensione della griglia di calcolo
      *ndim, //puntatore al vettore dimensione che conterra' i valori di ogni dimensione della griglia
      reorder, //indica la possibilita' di riordinare la griglia
      *period,  //periodicita' della griglia
      coordinate[2], //coordinate dei processori nella griglia 
      i,j;  //indice utilizzato per i cicli for
  float  **A, //puntatore alla matrice 
         *X, //puntatore al vettore da moltiplicare alla matrice
         *B,  //puntatore al vettore risultante  
         *subA, //puntatore alla sottomatrice creata in P0 ma assegnata agli altri processori
         *subX, //puntatore al vettore parziale  
         *subB,  //vettore risultante parziale
         *subA0,  //sottomatrice assegnata a P0
         *subX0, //vettore parziale assegnato a P0
         *sumB;   //vettore d'apppoggio per eseguire le somme dei vettori parziali
  double tempo1,tempo2; //variabili per il calcolo del tempo impiegato
  MPI_Status status;
  MPI_Comm griglia, //communicator griglia completa
           grigliar, //communicator griglia divisa per righe
           grigliac;  //communicator griglia divisa per colonne
        
  MPI_Init(&argc,&argv);   //Inizio programma MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&menum);   //Acquisizione del proprio identificativo	
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);   //Acquisizione del numero di processori 
  
  if(menum == 0)  //FASE DI INIZIALIZZAZIONE
    input[0]=Controllo_input(input,argc,&(*argv),nproc); //controllo input il cui risultato va in input[0]. 0 input ok, 1 input errato             
  MPI_Bcast(input,4,MPI_INT,0,MPI_COMM_WORLD);    //invio in broadcast da 0 , degli input e risultato controllo_input
  if(input[0] == 0){                            //se l'input e' corretto proseguo l'elaborazione
    colonne_griglia=nproc/input[1];      //colonne griglia
    Crea_griglia (&griglia,&grigliar,&grigliac, menum,nproc,input[1],colonne_griglia,coordinate);  // crezione griglia
    MPI_Comm_rank(griglia,&menum);   //Acquisizione del proprio identificativo nella griglia	
    MPI_Cart_coords (griglia,menum,2,coordinate);  //prendo le coordinate del processore nella griglia 
    mloc=input[2]/input[1];       //mloc = M/ righe_griglia
    nloc=input[3]/colonne_griglia;  //nloc = N/ colonne_griglia
    r=input[2]%input[1];   //controllo se le righe della matrice sono divisibili per le righe della griglia
    c=input[3]%colonne_griglia;  //controllo se le colonne della matrice sono divisibili per le colonne della griglia
    if(r != 0){      //ogni processore calcola il suo mloc ed nloc cone le righe e colonne aggiunte in caso m ed n nn siano divisibili equamente
      if(coordinate[0] == 0) //in questo modo assegno le righe in piu' della matrice ai processori della prima parte della griglia
        mloc+=r;  
    }
    if(c != 0){
      if(coordinate[1] == 0) //in questo modo assegno le colonne in piu' della matrice ai processori della prima parte della griglia
        nloc+=c;     
    } 
    subA=(float*)malloc(mloc*nloc*(sizeof(float))); //allocazione matrice parziale che ogni processore ricevera' da P0
    subX=(float*)malloc(nloc*(sizeof(float))); //allocazione vettore parziale che ogni processore ricevera' da P0
    subB=(float*)malloc(mloc*(sizeof(float))); //allocazione vettore risultante parizale calcolato su ogni processore
    sumB=(float*)malloc(mloc*(sizeof(float)));  //allocazione del vettore di appoggio
    B=(float*)malloc(input[2]*(sizeof(float)));  //alloco il vettore risultante  
    for (i=0;i<mloc;i++){
       subB[i]=0;
       sumB[i]=0;
    }
    for (i=0;i<input[2];i++)
	  B[i]=0;      
    if( menum != 0){    //se non sono P0 mi metto in attesa di ricevere dati
      MPI_Recv(subA,(mloc*nloc),MPI_FLOAT,0,80,griglia,&status);  //ricevo sottomatrice da P0 sottoforma di vettore
      MPI_Recv(subX,(nloc),MPI_FLOAT,0,80,griglia,&status);    //ricevo vettore parziale da P0
    }
    else if (menum == 0){   //se sono in P0
      srand(time(NULL));    //inizializzo seme per la random       
      A=Crea_mat_random(input[2],input[3]);  //creo matrice con elementi random
      X=Crea_vet_random(input[3]);    //creo vettore con elementi random
      subA0=(float*)malloc(mloc*nloc*(sizeof(float)));  //alloco matrice parziale che P0 si assegnera' durante la ripartizione
      subX0=(float*)malloc(nloc*(sizeof(float)));   //vettore parziale ch P0 si assegnera' durante la ripartizione
      for(i=0; i<nproc; i++){        //ripeto il ciclo per ogni processore
        MPI_Cart_coords(griglia,i,2,coordinate);   //prendo le coordinate nella griglia del processore i
        if(i==0){   //se sono P0
          mloc0=mloc;  
          nloc0=nloc; 
          Ripartizione_mat_vet(griglia,A,input[2],input[3],subA0,mloc0,nloc0,r,c,X,subX0,coordinate); //ripartizione 
        }
        else{
          mloc=input[2]/input[1];       //mloc = M/ righe_griglia
          nloc=input[3]/colonne_griglia;  //nloc = N/ colonne_griglia
          if(r != 0){      //ogni processore calcola il suo mloc ed nloc cone le righe e colonne aggiunte in caso m ed n nn siano divisibili equamente
            if(coordinate[0] == 0) //in questo modo assegno le righe in piu' della matrice ai processori della prima parte della griglia
              mloc+=r;  
          }
          if(c != 0){
            if(coordinate[1] == 0) //in questo modo assegno le colonne in piu' della matrice ai processori della prima parte della griglia
              nloc+=c;     
          }   
          Ripartizione_mat_vet(griglia,A,input[2],input[3],subA,mloc,nloc,r,c,X,subX,coordinate);
          MPI_Cart_rank(griglia,coordinate,&menum); //funzione inversa a cart_coords cioe' tramite le coordinate nella griglia conosco il rank del processore 
          MPI_Send(subA,(mloc*nloc),MPI_FLOAT,menum,80,griglia); //invio sottomatrice  
          MPI_Send(subX,nloc,MPI_FLOAT,menum,80,griglia);   //invio vettore parziale
        }
      } 
      free(A);  //dealloco lo spazio occupato dalla matrice e dal vettore, iniziali, creati solamente in P0
      free(X); 
    }
    MPI_Barrier(griglia);
    MPI_Comm_rank(griglia,&menum);   //Acquisizione del proprio identificativo nella griglia perche P0 perde il suo menum durante l'elaborazione
    tempo1=MPI_Wtime();  //inizio registrazione tempi di calcolo
    if(menum == 0){
      Prodotto(mloc0,nloc0,subA0,subX0,subB); //Calcolo dei prodotti parziali    
      mloc=mloc0;
      nloc=nloc0;
    }
    else{ 
      Prodotto(mloc,nloc,subA,subX,subB); //Calcolo dei prodotti parziali    
    }
    MPI_Allreduce(subB,sumB,mloc,MPI_FLOAT,MPI_SUM,grigliar);   //somma dei vettori parziali secondo la riga dei processori 
    MPI_Allgather(sumB,mloc,MPI_FLOAT,B,mloc,MPI_FLOAT,grigliac);  //raccolta dei risultati
    tempo2=MPI_Wtime();  //fine registrazione tempi di calcolo
    if (input[2]<=50 && menum == 0){    
      printf ("\nRisultato = ",menum);
      Stampa_vet (B,input[2]);
    }
    else
      printf ("\nIl vettore risultante ha troppi elementi e per questo non verra' stampato\n");

    printf ("\nTempo impiegato = %lf sec\n",tempo2-tempo1);
    free(B);   //deallocazione risorse
    free(subA);
    free(subX);
    if(menum == 0 ){
      free(subA0);
      free(subX0);         
    }
  }     
  
  MPI_Finalize();
}

/**
 *
 * Questa funzion si occupa di effettuare i principali controllo sugli input inseriti tramite PBS dall'utente
 * Ritorna 1 in caso di errore e 0 in caso di successo, ed il valore di ritorno viene raccolto in input[0]
 *
 */
int Controllo_input(int input[], int argc, char *argv[], int nproc ){  
  int resto;  
  if(argc >= 4){    //controllo dei parametri necessario o le atoi causano un core dump
    input[1]=atoi(argv[1]); //righe griglia   
    input[2]=atoi(argv[2]); //numero righe matrice
    input[3]=atoi(argv[3]); //numero colonne matrice
    resto=nproc%input[1];
    if(nproc < 1 || nproc > 8){   //controllo numero processori utilizzati
      printf("ERRORE: il numero di processori da impiegare deve essere compreso tra 2 e 8!!!\n");
      return 1; 
    }
    else if(input[1] > nproc || input[1] < 1){  
      printf("ERRORE: il parametro 'righe griglia' deve essere compreso tra 1 e nproc!!!\n");
      return 1;}
    else if(resto != 0){  
      printf("ERRORE: il numero di processori deve essere uguale al numero di processori impiegati nella griglia\n");     
      return 1;}
    else if((input[2]*input[3]) > 1000000000 || (input[2]*input[3]) < nproc ){
      printf("ERRORE: il numero totale di elementi della matrice deve essere compreso tra nproc e 1.000.000.000!!!\n");     
      return 1;}
    else { return 0;}
  }
  else{
    printf("ERRORE: errore nei parametri!\nassicurarsi di aver inserito -righe griglia -righe matrice -colonne matrice.\n");
    return 1;}
}

/**
 *
 * Questa funzione si occupa di creare la griglia di calcolo restituendo come output griglia,cioe' la griglia di calcolo completa;
 * inoltre ritorna in output altri due communicator cioe' grigliar e grigliac che sarebbero la griglia di calcolo divisa per righe e per
 * colonne.Questi sottocommunicator verranno utilizzato dalle funzioni MPI Allreduce e Allgather
 *
*/
void Crea_griglia(MPI_Comm *griglia, MPI_Comm *grigliar, MPI_Comm *grigliac, int menum, int nproc, int riga, int col, int *coordinate){
  int dim,    //dimensione griglia
      *ndim,
      reorder,
      *period,
      vc[2];
  dim = 2;
  ndim = (int*) calloc (dim, sizeof(int));
  ndim[0] = riga;
  ndim[1] = col;
  period = (int*) calloc (dim, sizeof(int));
  period[0] = period [1] = 0;
  reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, dim, ndim,period, reorder, &(*griglia));
  vc[0] = 0;
  vc[1] = 1;
  MPI_Cart_sub(*griglia, vc,&(*grigliar));
  vc[0] = 1;
  vc[1] = 0;
  MPI_Cart_sub(*griglia, vc, &(*grigliac));
}

/**
 *
 * Questa funzione si occupa di creare una matrice di dimensione righexcolonne costituita da elementi del campo dei numeri reali
 * creati in maniera casuale e che come valore massimo hanno 100.0. Il valore di ritorno di questa funzione e' il puntatore alla
 * matrice
 *
*/
float** Crea_mat_random(int righe, int colonne){
  int i,j;  //indici
  float **A;  //puntatore alla matrice
  A = (float**)malloc(righe*(sizeof(float*)));  //inizio allocazione matrice
  for (i=0; i<righe; i++){
    A[i] = (float*)malloc(colonne*sizeof(float));  
  }	                      //fine allocazione matrice
  for(i=0;i < righe;i++){
    for(j=0;j < colonne;j++){       
      A[i][j]=(float)rand()/((float)RAND_MAX/100.0);  //creazione elementi casuali
  }}   
  return A;     //restituzione puntatore alla matrice
}

/**
 *
 * Questa funzione crea un vettore di dimensione dim costituito da elementi del campo dei numeri reali, creati casualmente.Il valore
 * restituito da questo funzione e' il puntatore al vettore creato
 *
*/
float* Crea_vet_random(int dim){
  int i;
  float *V;
  V = (float*)malloc(dim*(sizeof(float)));	
  for(i=0;i < dim;i++){       
    V[i]=(float)rand()/((float)RAND_MAX/100.0); 
  }   
  return V;     
}


void Stampa_mat(float **A,int m,int n){ 
  int i,j;   
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      printf("%f  ",A[i][j]);
    }
    printf("\n");
  }    
}

void Stampa_vet(float *X, int n) {
    int i;
    printf("[ ");
    for (i = 0; i < n; i++) {
        printf("%f ", X[i]);
    }
    printf("]\n");
}


/**
 *
 * Questa funzione si occupa di ripartire la matrice e il vettore tra i vari processori.Viene eseguita solo da P0 nella fase
 * di inizializzazione, e gli output restituiti da questa funzione, cioe' la matrice parziale e il vettore parziale, tramite un calcolo
 * viene inviato da P0 al processore designato. In particolare partendo dalla matrice A e dal vettore X, questa funzione calcola, tramite le coordinate
 * dei processori all'interno della griglia, quale porzione di A e di X leggere e le va a scrivere nella sottomatrice subA e nel vettore parziale subX.
 * La matrice e il vettore parziale in P0 vengono allocati una sola volta e sono sempre gli stessi quindi non c'e' spreco di memoria in quanto vengono sempre
 * sovrascritti, visto che ripartizione_mat_vet e' chiamate tante volte quanti sono i processori.
 *
*/
void Ripartizione_mat_vet(MPI_Comm griglia, float** A, int m, int n, float subA[], int mloc, int nloc, int r, int c, float X[], float subX[], int coordinate[]){
  int i,j,z,   //indici
      startR,  //segnala da quale riga iniziare
      startC,  //segnala da quale colonna iniziare
      menum;  //id processore
      
  startR=coordinate[0]*mloc; //serve per sapere da quale riga dobbiamo iniziare a leggere gli elementi di A 
  startC=coordinate[1]*nloc;
  if(r != 0){      //ogni processore calcola il suo mloc ed nloc cone le righe e colonne aggiunte in caso m ed n nn siano divisibili equamente
      if(coordinate[0] != 0) //in questo modo assegno le righe in piu' della matrice ai processori della prima parte della griglia
        startR+=r;    //cioe' a tutti i processori la cui coordinata riguardante le righe della griglia e' 0
    }
    if(c != 0){
      if(coordinate[1] != 0) //in questo modo assegno le colonne in piu' della matrice ai processori della prima parte della griglia
        startC+=c;      //stesso ragionamento di startR ma sulle colonne
    } 
  i=0;
  for(j=0;j<mloc;j++){ //ripartisco la matrice
    for(z=0;z<nloc;z++){ //j e z sono usati come offset da sommare a startR e startC. 
      subA[i]=A[startR+j][startC+z];
      i++;
    }                    
  }
 for(i=0;i<nloc;i++){  //ripartisco il vettore X
    subX[i]=X[startC+i];                      
 }
}


/**
 *
 * Questa funzione si occupa di svolgere il prodotto matrice vettore.come output restituisce un vettore contenente i prodotti.
 *
*/
void Prodotto(int mloc, int nloc, float *subA, float *subX, float *subB){ 	
   int i,j;
   for (i=0; i<mloc; i++){ 		
      for (j=0; j<nloc; j++){ 			
         subB[i] += subA[i*nloc+j] * subX[j]; 
      }
   }
}
