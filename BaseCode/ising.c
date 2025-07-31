#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define N 32
#define Temperatura 0.5
//#define DEBUG

typedef struct{
  int spins[N][N];
  float energy;
  float spin_sum;
  float J;
}Ising;

void create_spin_lattice(Ising *model, int proportion){
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

void calculate_spin_sum(Ising *model){
 model->spin_sum = 0;
 for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    model->spin_sum+=model->spins[i][j];
}

void calculate_hamiltonian(Ising *model){
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

void isingStep(Ising *X_t, float Temp){
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

void print_model(Ising *X){
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

int main() {
    clock_t t,t_total;
    t_total = clock();
    srand(time(NULL));
    unsigned int executions = 1000000;
    FILE *statistics;
    
    //Prepara o modelo
    Ising X;
    X.J = 1;
    
    t = clock();
    create_spin_lattice(&X,750);
    t = clock() - t;
    printf("Tempo create_lattice:%f\n",(float)t/CLOCKS_PER_SEC);
    X.energy = 0;
    X.spin_sum = 0;
    calculate_hamiltonian(&X);
    calculate_spin_sum(&X);
    #ifdef DEBUG
    printf("X\nenergy: %f\n",X.energy);
    //print_model(&X);
    #endif
    
    //Executa o simulador
    statistics = fopen("statistics.txt","w");
    if(statistics == NULL){
      printf("Can't open file");
    }else{
      t = clock();
      for(int i=0; i<executions; i++){
        isingStep(&X,Temperatura);
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
    //print_model(&X);
    #endif
    t_total = clock() - t_total;
    printf("Tempo total:%f\n",(float)t_total/CLOCKS_PER_SEC);
    return 0;
}
