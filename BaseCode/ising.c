#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#define N 32
//#define Temperatura 10
//#define DEBUG

typedef struct{
  int **spins;
  float energy;
  float spin_sum;
  float J;
}Ising;

int ** allocate_spins(int N){
  int **matrix = malloc(N * sizeof(int *));
  for(int i=0; i<N; i++){
    matrix[i] = malloc(N * sizeof(int));
  }
  return matrix;
}

void free_spins(int **matrix, int N){
  for(int i=0; i<N; i++){
    free(matrix[i]);
  }
  free(matrix);
}


void create_spin_lattice(Ising *model, int proportion, int N){
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      float random = (rand() % 1000) + 1;
      if(random >= proportion)
        model->spins[i][j] = 1;
      else
        model->spins[i][j] = -1;
    }
  }
}

void calculate_spin_sum(Ising *model, int N){
 model->spin_sum = 0;
 for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    model->spin_sum+=model->spins[i][j];
}

void calculate_hamiltonian(Ising *model, int N){
  model->energy = 0;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      int current_spin = model->spins[i][j];
      if(i!=0)
        model->energy+= model->spins[i-1][j]*current_spin;
      if(i!=N-1)
        model->energy+= model->spins[i+1][j]*current_spin;
      if(j!=0)
        model->energy+= model->spins[i][j-1]*current_spin;
      if(j!=N-1)
        model->energy+= model->spins[i][j+1]*current_spin;
    }
  }
  model->energy*=-model->J;
}

void isingStep(Ising *X_t, float Temp, int N){
    float dE=0, dEy=0, dEx=0;
 
    //Rotate random spin
    int i=(rand() % N);
    int j=(rand() % N);
    #ifdef DEBUG
    printf("Rotate (%d,%d)\n",i,j);
    #endif
    int current_Y_spin = X_t->spins[i][j]*-1;

    //Calculate delta E
    if(i!=0){
      dEy+= X_t->spins[i-1][j]*current_Y_spin;
      dEx+= X_t->spins[i-1][j]*X_t->spins[i][j];
    }
    if(i!=N-1){
      dEy+= X_t->spins[i+1][j]*current_Y_spin;
      dEx+= X_t->spins[i+1][j]*X_t->spins[i][j];
    }
    if(j!=0){
      dEy+= X_t->spins[i][j-1]*current_Y_spin;
      dEx+= X_t->spins[i][j-1]*X_t->spins[i][j];
    }
    if(j!=N-1){
      dEy+= X_t->spins[i][j+1]*current_Y_spin;
      dEx+= X_t->spins[i][j+1]*X_t->spins[i][j];      
    }
    dE = -1*X_t->J*(dEy - dEx);
    
    //Metropolis step
    if(dE<0){
      #ifdef DEBUG
      printf("dE<0");
      #endif
      X_t->spins[i][j] = current_Y_spin;
      X_t->energy += dE;
      X_t->spin_sum += 2*current_Y_spin;
    }else{
      float U = ((rand() % 1000) + 1)*0.001;
      float boltz_dist = -dE/Temp;
      float probability = exp(boltz_dist);
      float alpha = fmin(1,probability);
      if(U<alpha){
        #ifdef DEBUG
        printf("dE>0 e U<alpha");
        #endif
        X_t->spins[i][j] = current_Y_spin;
        X_t->energy += dE;
        X_t->spin_sum += 2*current_Y_spin;
      }
    }
}

void print_model(Ising *X, int N){
  for(int j=0; j<N; j++){
    for(int i=0; i<N; i++){
      if(X->spins[i][j] == -1)
        printf("%d,",X->spins[i][j]);
      else
        printf(" %d,",X->spins[i][j]);
    }
    printf("\n");
  }
}

int main(int argc, char **argv) {

    int proportion = 750;
    unsigned int steps = 1000000;
    int N = 32;
    float Temperatura = 5;

    clock_t t,t_total;
    t_total = clock();
    srand(time(NULL));
    FILE *statistics;
    
    //Prepara o modelo
    Ising X;
    X.J = 1;
    
    for (int i = 0; i < argc; i++) {
      printf("Argumento [%d]: %s\n", i, argv[i]);
  
      if(strcmp(argv[i], "-P") == 0){
        proportion = atoi(argv[i+1]);
      }  

      if(strcmp(argv[i], "-steps") == 0){
        steps = atoi(argv[i+1]);
      }      
      
      if(strcmp(argv[i], "-N") == 0){
        N = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-J") == 0){
        X.J = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-T") == 0){
        Temperatura = atof(argv[i+1]);
      }
    }

    X.spins = allocate_spins(N);
    

    t = clock();
    create_spin_lattice(&X,proportion,N);
    t = clock() - t;
    printf("Tempo create_lattice:%f\n",(float)t/CLOCKS_PER_SEC);
    X.energy = 0;
    X.spin_sum = 0;

    calculate_hamiltonian(&X,N);
    calculate_spin_sum(&X,N);
    #ifdef DEBUG
    printf("X\nenergy: %f\n",X.energy);
    //print_model(&X,N);
    #endif
    
    //Executa o simulador
    statistics = fopen("statistics.txt","w");
    if(statistics == NULL){
      printf("Can't open file");
    }else{
      t = clock();
      for(int i=0; i<steps; i++){
        isingStep(&X,Temperatura,N);
        //Escreve no documento
        char mensagem[50];
        sprintf(mensagem, "%d,%.5f,%.3f\n", i, X.energy, X.spin_sum/(N*N));
        fputs(mensagem, statistics);
      }    
      t = clock() - t;
      printf("Tempo isingStep:%f\n",(float)t/CLOCKS_PER_SEC);
    }
    fclose(statistics);
    #ifdef DEBUG
    printf("X\nenergy: %f\n",X.energy);
    //print_model(&X,N);
    #endif

    free_spins(X.spins,N);

    t_total = clock() - t_total;
    printf("Tempo total:%f\n",(float)t_total/CLOCKS_PER_SEC);
    return 0;
}
