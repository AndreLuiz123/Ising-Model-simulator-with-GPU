#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
//#define N 32
//#define Temperatura 10
//#define DEBUG

typedef struct{
  int **spins;
  double energy;
  double spin_sum;
  double J;
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
      double random = (rand() % 1000) + 1;
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

void calculate_hamiltonian_torus_boundary(Ising *model, int M, int N){
  double energy = 0;
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

void calculate_hamiltonian(Ising *model, int M, int N, bool torus){
  if(torus){
    calculate_hamiltonian_torus_boundary(model, M, N);
  }else{
    calculate_hamiltonian_classic_boundary(model, M, N);
  }
}

double calculate_dE_classic_boundary(Ising *model ,int i, int j, int M, int N){
    
    double dE = 0;
    double dEy = 0;
    double dEx = 0;
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

double calculate_dE_torus_boundary(Ising *model ,int i, int j, int M, int N){
    double dE = 0;
    //Calculate delta E
    double neighbor_interaction =
    model->spins[(i + M - 1) % M][j] + model->spins[(i + 1) % M][j]+
    model->spins[i][(j + N - 1) % N] + model->spins[i][(j + 1) % N];

    dE = 2.0 *  model->spins[i][j]* neighbor_interaction;

    return dE;
}

double calculate_dE(Ising *model ,int i, int j, int M, int N, bool torus){
  double result = 0;
  if(torus){
    result = calculate_dE_torus_boundary(model, i, j, M, N);
  }else{
    result = calculate_dE_classic_boundary(model, i, j, M, N);
  }
  return result;
}

void isingStep(Ising *X_t, double Temp, int M, int N, bool torus){
    double dE=0;
 
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
    dE = calculate_dE(X_t, i, j, M, N, torus);//-1*X_t->J*(dEy - dEx);
    
    //Metropolis step
    if(dE<0){
      #ifdef DEBUG
      printf("dE<0");
      #endif
      X_t->spins[i][j] = current_Y_spin;
      X_t->energy += dE;
      X_t->spin_sum += 2*current_Y_spin;
    }else{
      double U = ((rand() % 1000) + 1)*0.001;
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

void monte_carlo_execution(Ising *model, int Temperatura, int M, int N,int steps, bool torus){
    FILE *statistics;
    statistics = fopen("statistics.txt","w");
    if(statistics == NULL){
      printf("Can't open file");
    }else{
      for(int i=0; i<steps; i++){
        isingStep(model,Temperatura,M,N, torus);
        //Escreve no documento
        char mensagem[50];
        sprintf(mensagem, "%d,%.5f,%.3f\n", i, model->energy, model->spin_sum/(M*N));
        fputs(mensagem, statistics);
      }    
    }
    fclose(statistics); 
}

void specific_heat_calculation(Ising *model, int M, int N, int steps, bool torus, int proportion, double temp0, double tempMax, double tempAdd){
    
    FILE *statistics;
    statistics = fopen("specific_heat.txt","w");
    if(statistics == NULL){
      printf("Can't open file");
    }else{
      int equilibrium_steps = steps/5;
      int specific_heat_collecting_steps = steps - equilibrium_steps;

      while(temp0 <= tempMax){

        double mean_energy = 0;
        double mean_square_energy = 0;
        printf("Temperatura:%f\n",temp0);  

        create_spin_lattice(model,proportion,M,N);
        model->energy = 0;
        model->spin_sum = 0;
        
        calculate_hamiltonian(model,M,N,torus);
        
        monte_carlo_execution(model, temp0, M, N, equilibrium_steps,torus);
        
        for(int i=0; i<specific_heat_collecting_steps; i++){
          isingStep(model,temp0,M,N,torus);
          mean_energy += model->energy;
          mean_square_energy += model->energy*model->energy;
        }    

        mean_energy/=specific_heat_collecting_steps;
        mean_square_energy/=specific_heat_collecting_steps;
        double calor_especifico = (mean_square_energy - mean_energy*mean_energy)/(temp0*temp0*M*N);

        fprintf(statistics, "%.3f,%.5f\n", temp0, calor_especifico);

        temp0+=tempAdd;
      }
      
    }
    fclose(statistics); 
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
    unsigned int mcs = 1000; //flips por spin
    int M = 32, N = 32;
    double Temperatura = 5;
    bool torus = true;

    clock_t t,t_total;
    t_total = clock();
    srand(time(NULL));
    
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

      if(strcmp(argv[i], "-B") == 0){
        torus = atoi(argv[i+1]);
      }
    }

    int steps = M*N*mcs;

    X.spins = allocate_spins(M,N);
    

    t = clock();
    create_spin_lattice(&X,proportion,M,N);
    t = clock() - t;
    printf("Tempo create_lattice:%f\n",(double)t/CLOCKS_PER_SEC);
    X.energy = 0;
    X.spin_sum = 0;

    calculate_hamiltonian(&X,M,N,torus);
    calculate_spin_sum(&X,M,N);
    #ifdef DEBUG
    printf("X\nenergy: %f\n",X.energy);
    //print_model(&X,N);
    #endif
    
    //Executa o simulador
    t = clock();
    monte_carlo_execution(&X,Temperatura,M,N,steps,torus);
    t = clock() - t;
    printf("Tempo isingStep:%f\n",(double)t/CLOCKS_PER_SEC);
    
    #ifdef DEBUG
    printf("X\nenergy: %f\n",X.energy);
    //print_model(&X,N);
    #endif

    free_spins(X.spins,N);

    t_total = clock() - t_total;
    printf("Tempo total:%f\n",(double)t_total/CLOCKS_PER_SEC);
    return 0;
}
