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
  double energy;
  double spin_sum;
  double J;
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
      double random = (rand() % 1000) + 1;
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

void calculate_hamiltonian(Ising *model, int N, double a){
  model->energy = 0;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      int current_spin = model->spins[i][j];
      double J0 = model->J;
      double J1 = (1-a)*J0;
      double J2 = -1*J0;


      model->energy += 
      J0*model->spins[i][(j+N-1)%N]*current_spin +
      J0*model->spins[i][(j+1)%N]*current_spin +
      J1*model->spins[(i+N-1)%N][j]*current_spin +
      J1*model->spins[(i+1)%N][j]*current_spin +
      J2*model->spins[(i+N-2)%N][j]*current_spin +
      J2*model->spins[(i+2)%N][j]*current_spin;
    }
  }
  model->energy*=-0.5;
}

void isingStep(Ising *X_t, double Temp, int N, double a){
    double dE=0, dEy=0, dEx=0;
 
    //Rotate random spin
    int i=(rand() % N);
    int j=(rand() % N);
    #ifdef DEBUG
    printf("Rotate (%d,%d)\n",i,j);
    #endif
    int current_Y_spin = X_t->spins[i][j]*-1;

    double J0 = X_t->J;
    double J1 = (1-a)*J0;
    double J2 = -1*J0;
    //Calculate delta E
double neighbor_interaction =
    J0 * (X_t->spins[i][(j + N - 1) % N] + X_t->spins[i][(j + 1) % N]) +
    J1 * (X_t->spins[(i + N - 1) % N][j] + X_t->spins[(i + 1) % N][j]) +
    J2 * (X_t->spins[(i + N - 2) % N][j] + X_t->spins[(i + 2) % N][j]);

    dE = 2.0 *  X_t->spins[i][j]* neighbor_interaction;
    
    //Metropolis step
    if(dE<0){
      #ifdef DEBUG
      printf("dE<0");
      #endif
      X_t->spins[i][j] = current_Y_spin;
      X_t->energy += dE;
      X_t->spin_sum += 2*current_Y_spin;
    }else{
      double U = (rand() / (double)RAND_MAX);
      double boltz_dist = -dE/Temp;
      double probability = exp(boltz_dist);
      double alpha = fmin(1,probability);
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

    int proportion = 500;
    unsigned int mcs = 20000; // monte carlo per spin
    int N = 20;
    double Temperatura = 5;
    double a = 0;

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
      
      if(strcmp(argv[i], "-N") == 0){
        N = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-J") == 0){
        X.J = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-T") == 0){
        Temperatura = atof(argv[i+1]);
      }

      if(strcmp(argv[i], "-a") == 0){
        a = atof(argv[i+1]);
      }
    }

    X.spins = allocate_spins(N);
    
    unsigned int total_flips = mcs*N*N;
    unsigned int equilibrium_flips = total_flips/5;
    unsigned int oficial_flips = total_flips - equilibrium_flips;

    //Executa o simulador
    statistics = fopen("statistics.txt","w");
    if(statistics == NULL){
      printf("Can't open file");
    }else{
      Temperatura = 0.03;
      while(Temperatura<=3){

        t = clock();
        create_spin_lattice(&X,proportion,N);
        t = clock() - t;
        printf("Tempo create_lattice:%f\n",(double)t/CLOCKS_PER_SEC);
        X.energy = 0;
        X.spin_sum = 0;
        
        calculate_hamiltonian(&X,N,a);
        #ifdef DEBUG
        printf("X\nenergy: %f\n",X.energy);
        //print_model(&X,N);
        #endif
        
        for(int i=0; i<equilibrium_flips; i++){
          isingStep(&X,Temperatura,N,a);          
        }
        double mean_energy = 0;
        double mean_square_energy = 0;

        t = clock();
        for(int i=0; i<oficial_flips; i++){
          isingStep(&X,Temperatura,N,a);
          mean_energy += X.energy;
          mean_square_energy += X.energy*X.energy;
          //Escreve no documento
        }    
        t = clock() - t;
        printf("Tempo isingStep:%f\n",(double)t/CLOCKS_PER_SEC);
        mean_energy/=oficial_flips;
        mean_square_energy/=oficial_flips;
        double calor_especifico = (mean_square_energy - mean_energy*mean_energy)/(Temperatura*Temperatura*N*N);
        char mensagem[50];
        sprintf(mensagem, "%.3f,%.5f\n", Temperatura, calor_especifico);
        fputs(mensagem, statistics);
        printf("Temperatura:%f\n",Temperatura);
        Temperatura+=0.03;
      }
    }
    fclose(statistics);
    #ifdef DEBUG
    printf("X\nenergy: %f\n",X.energy);
    //print_model(&X,N);
    #endif

    free_spins(X.spins,N);

    t_total = clock() - t_total;
    printf("Tempo total:%f\n",(double)t_total/CLOCKS_PER_SEC);
    return 0;
}
