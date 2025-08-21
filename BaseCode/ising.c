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

int ** allocate_spins(int M,int N){
  int **matrix = malloc(M * sizeof(int *));
  for(int i=0; i<M; i++){
    matrix[i] = malloc(N * sizeof(int));
  }
  return matrix;
}

void free_spins(int **matrix, int M){
  for(int i=0; i<M; i++){
    free(matrix[i]);
  }
  free(matrix);
}


void create_spin_lattice(Ising *model, int proportion, int M, int N){
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      float random = (rand() % 1000) + 1;
      if(random >= proportion)
        model->spins[i][j] = 1;
      else
        model->spins[i][j] = -1;
    }
  }
}

void calculate_spin_sum(Ising *model, int M, int N){
 model->spin_sum = 0;
 for(int i=0; i<M; i++)
  for(int j=0; j<N; j++)
    model->spin_sum+=model->spins[i][j];
}

void calculate_hamiltonian_classic_boundary(Ising *model, int M, int N){

  model->energy = 0;
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      int current_spin = model->spins[i][j];
      if(i!=0)
        model->energy+= model->spins[i-1][j]*current_spin;
      if(i!=M-1)
        model->energy+= model->spins[i+1][j]*current_spin;
      if(j!=0)
        model->energy+= model->spins[i][j-1]*current_spin;
      if(j!=N-1)
        model->energy+= model->spins[i][j+1]*current_spin;
    }
  }
  model->energy*=-model->J;
}

float calculate_dE_classic_boundary(Ising *model ,int i, int j, int M, int N){
    
    float dE = 0;
    float dEy = 0;
    float dEx = 0;
    int new_spin = model->spins[i][j]*-1;
    int current_spin = model->spins[i][j];
    //Calculate delta E
    if(i!=0){
      dEy+= model->spins[i-1][j]*new_spin;
      dEx+= model->spins[i-1][j]*current_spin;
    }
    if(i!=M-1){
      dEy+= model->spins[i+1][j]*new_spin;
      dEx+= model->spins[i+1][j]*current_spin;
    }
    if(j!=0){
      dEy+= model->spins[i][j-1]*new_spin;
      dEx+= model->spins[i][j-1]*current_spin;
    }
    if(j!=N-1){
      dEy+= model->spins[i][j+1]*new_spin;
      dEx+= model->spins[i][j+1]*current_spin;      
    }

    dE = -1*model->J*(dEy - dEx);

    return dE;
}

void calculate_hamiltonian_torus_boundary(Ising *model, int M, int N){
  float energy = 0;
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      int current_spin = model->spins[i][j];
      energy += 
      model->spins[(i+M-1)%M][j]*current_spin +
      model->spins[(i+1)%M][j]*current_spin + 
      model->spins[i][(j+N-1)%N]*current_spin +
      model->spins[i][(j+1)%N]*current_spin;
    }
  }
  energy*=-model->J;
  model->energy = energy;
}

float calculate_dE_torus_boundary(Ising *model ,int i, int j, int M, int N){
    float dE = 0;
    //Calculate delta E
    double neighbor_interaction =
    model->spins[(i + M - 1) % M][j] + model->spins[(i + 1) % M][j]+
    model->spins[i][(j + N - 1) % N] + model->spins[i][(j + 1) % N];

    dE = 2.0 *  model->spins[i][j]* neighbor_interaction;

    return dE;
}

void isingStep(Ising *X_t, float Temp, int M, int N){
    float dE=0;
 
    //Rotate random spin
    int i=(rand() % M);
    int j=(rand() % N);
    #ifdef DEBUG
    printf("Rotate (%d,%d)\n",i,j);
    #endif
    int current_Y_spin = X_t->spins[i][j]*-1;

    //Calculate delta E
    /*if(i!=0){
      dEy+= X_t->spins[i-1][j]*current_Y_spin;
      dEx+= X_t->spins[i-1][j]*X_t->spins[i][j];
    }
    if(i!=M-1){
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
    }*/
    dE = calculate_dE_torus_boundary(X_t, i, j, M, N);//-1*X_t->J*(dEy - dEx);
    
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

void print_model(Ising *X, int M, int N){
  for(int j=0; j<M; j++){
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
    unsigned int mcs = 1000000; //flips por spin
    int M = 32, N = 32;
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

      if(strcmp(argv[i], "-mcs") == 0){
        mcs = atoi(argv[i+1]);
      }      

      if(strcmp(argv[i], "-M") == 0){
        M = atoi(argv[i+1]);
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

    int steps = M*N*mcs;
    //int initialization_steps = M*N*mcs/5;
    //int specific_heat_collecting_steps = M*N*mcs - initialization_steps;

    X.spins = allocate_spins(M,N);
    

    t = clock();
    create_spin_lattice(&X,proportion,M,N);
    t = clock() - t;
    printf("Tempo create_lattice:%f\n",(float)t/CLOCKS_PER_SEC);
    X.energy = 0;
    X.spin_sum = 0;

    calculate_hamiltonian_torus_boundary(&X,M,N);
    calculate_spin_sum(&X,M,N);
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
        isingStep(&X,Temperatura,M,N);
        //Escreve no documento
        char mensagem[50];
        sprintf(mensagem, "%d,%.5f,%.3f\n", i, X.energy, X.spin_sum/(M*N));
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
