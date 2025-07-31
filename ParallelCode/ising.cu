#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <cassert>
#include <iostream>
//#define N 32
//#define Temperatura 0.5
//#define DEBUG
#include <curand_kernel.h>

using namespace std;

__global__ void setup_kernel(curandState *state, int N,unsigned long long seed) {
    int tid = threadIdx.x + threadIdx.y * N;
    curand_init(seed, tid, 0, &state[tid]);
}

__global__ void generate_randoms(curandState *state, int N,float *randoms) {
    int tid = threadIdx.x + threadIdx.y * N;
    curandState localState = state[tid];
    randoms[tid] = curand_uniform(&localState); // entre 0.0 e 1.0
}


//Estratégias de update
__global__ void update_spin(int *lattice, int N){
  
  if(threadIdx.x < N && threadIdx.y < N){

      if(threadIdx.x%2 == 0 && threadIdx.y%2 == 0)
      {
        lattice[threadIdx.y*N + threadIdx.x] = 0;   
      }  
    
  }
}

__global__ void update_spin_checkboard(int *lattice, float *randoms, int J, float Temp,int N, bool white){
  
  if(threadIdx.x < N && threadIdx.y < N){
    int current_spin = lattice[threadIdx.y*N + threadIdx.x];
    int new_spin = lattice[threadIdx.y*N + threadIdx.x]*-1;
    int dE=0;
    int dE_current = 0; 
    int dE_new = 0;

    if(white)
    {
      if((threadIdx.x%2 == 0 && threadIdx.y%2 == 0) || (threadIdx.x%2 != 0 && threadIdx.y%2 != 0)){

        //Calculate delta E
        if(threadIdx.x != 0){
          dE_new += new_spin*(lattice[threadIdx.y*N + threadIdx.x - 1]);
          dE_current += current_spin*(lattice[threadIdx.y*N + threadIdx.x - 1]);
        }
        if(threadIdx.x != N-1){
          dE_new += new_spin*(lattice[threadIdx.y*N + threadIdx.x + 1]);
          dE_current += current_spin*(lattice[threadIdx.y*N + threadIdx.x + 1]);
        }
        if(threadIdx.y != 0){
          dE_new += new_spin*(lattice[(threadIdx.y-1)*N + threadIdx.x]);
          dE_current += current_spin*(lattice[(threadIdx.y-1)*N + threadIdx.x]);
        }
        if(threadIdx.y != N-1){
          dE_new += new_spin*(lattice[(threadIdx.y+1)*N + threadIdx.x]);
          dE_current += current_spin*(lattice[(threadIdx.y+1)*N + threadIdx.x - 1]);   
        }
        
        dE = -1*J*(dE_new - dE_current);

        if(dE <= 0){
          lattice[threadIdx.y*N + threadIdx.x] = new_spin; 
        }else{
          float boltz_dist = -dE/Temp; 
          float probability = expf(boltz_dist);
          float alpha = fminf(1,probability);
          if(randoms[threadIdx.y*N + threadIdx.x] < alpha){
            lattice[threadIdx.y*N + threadIdx.x] = new_spin; 
          }
        }
      }
    }else{
      if((threadIdx.x%2 != 0 && threadIdx.y%2 == 0) || (threadIdx.x%2 == 0 && threadIdx.y%2 != 0)){

        //Calculate delta E
        if(threadIdx.x != 0){
          dE_new += new_spin*(lattice[threadIdx.y*N + threadIdx.x - 1]);
          dE_current += current_spin*(lattice[threadIdx.y*N + threadIdx.x - 1]);
        }
        if(threadIdx.x != N-1){
          dE_new += new_spin*(lattice[threadIdx.y*N + threadIdx.x + 1]);
          dE_current += current_spin*(lattice[threadIdx.y*N + threadIdx.x + 1]);
        }
        if(threadIdx.y != 0){
          dE_new += new_spin*(lattice[(threadIdx.y-1)*N + threadIdx.x]);
          dE_current += current_spin*(lattice[(threadIdx.y-1)*N + threadIdx.x]);
        }
        if(threadIdx.y != N-1){
          dE_new += new_spin*(lattice[(threadIdx.y+1)*N + threadIdx.x]);
          dE_current += current_spin*(lattice[(threadIdx.y+1)*N + threadIdx.x - 1]);   
        }
        
        dE = -1*J*(dE_new - dE_current);

        if(dE <= 0){
          //printf("dE(%d) < 0\n", dE);
          lattice[threadIdx.y*N + threadIdx.x] = new_spin; 
        }else{
          float boltz_dist = -dE/Temp; 
          float probability = expf(boltz_dist);
          float alpha = fminf(1,probability);
          if(randoms[threadIdx.y*N + threadIdx.x] < alpha){
            //printf("random(%f) < alpha(%f) - (probability = %f; boltz = %f) - (dE = %d; Temp = %f)\n",randoms[threadIdx.y*N + threadIdx.x],alpha,probability,boltz_dist,dE, Temp);
            lattice[threadIdx.y*N + threadIdx.x] = new_spin; 
          }
        }
      }
    }
  }
}
//Calcula os valores de capa spin hamiltoniano
__global__ void calculate_spin_hamiltonian(int *lattice, int N, int *hamiltonian){
  
  if(threadIdx.x < N && threadIdx.y < N){

    int index = threadIdx.y*N + threadIdx.x;
    int current_spin = lattice[index];

      if(threadIdx.x!=0)
        hamiltonian[index]+= lattice[threadIdx.y*N + threadIdx.x-1]*current_spin;
      if(threadIdx.x!=N-1)
        hamiltonian[index]+= lattice[threadIdx.y*N + threadIdx.x+1]*current_spin;
      if(threadIdx.y!=0)
        hamiltonian[index]+= lattice[(threadIdx.y-1)*N + threadIdx.x]*current_spin;
      if(threadIdx.y!=N-1)
        hamiltonian[index]+= lattice[(threadIdx.y+1)*N + threadIdx.x]*current_spin;

  }
}

//Soma de elementos do vetor hamiltoniano - paralelizar depois
int calculate_hamiltonian(int *hamiltonian, int N, int J){
  
  int result=0;

  for(int i=0; i<N*N; i++){
    result += hamiltonian[i];
  }

  result*=-J;

  return result;

}

//Soma de elementos do modelo ising e tira a métia - paralelizar depois
float calculate_magnetization(int *lattice, int N){
  
  float result=0;

  for(int i=0; i<N*N; i++){
    result += lattice[i];
  }

  result = result/(N*N);

  return result;
}



void print_lattice(int *lattice, int N){
    for(int j=0; j<N; j++){
      for(int i=0; i<N; i++){
        if(lattice[j*N + i] == -1)
          printf("%d,",lattice[j*N + i]);
        else
          printf(" %d,",lattice[j*N + i]);
      }
      printf("\n");
    }  
}


void init_matrix(int *lattice, int N, int proportion=0){
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      float random = (rand() % 1000) + 1;
      if(random >= proportion)
        lattice[j*N + i] = 1;
      else
        lattice[j*N + i] = -1;
    }
  }
}

void init_white_black_matrix(int *lattice, int *white_matrix, int *black_matrix,int N){
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if((i%2==0 && j%2==0) || (i%2!=0 && j%2!=0)){
        white_matrix[j*N/2 + i] = lattice[j*N + i];
      }else{
        black_matrix[j*N/2 + i] = lattice[j*N + i];
      }
    }
  }
}

int calculate_hamiltonian2(int *lattice, int N, int J){
  int energy = 0;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      int current_spin = lattice[N*j + i];
      if(i!=0)
        energy+= lattice[N*j + i-1]*current_spin;
      if(i!=N-1)
        energy+= lattice[N*j + i+1]*current_spin;
      if(j!=0)
        energy+= lattice[N*(j-1) + i]*current_spin;
      if(j!=N-1)
        energy+= lattice[N*(j+1) + i]*current_spin;
    }
  }
  energy*=-1*J;
  return energy;
}

int main(int argc, char *argv[]){
  
  printf("Número de argumentos: %d\n", argc);
  //Parâmetros do sistema
  int steps = 1000;
  int J = 1;
  int N = 32;
  float Temp = 5;

  for (int i = 0; i < argc; i++) {
      printf("Argumento [%d]: %s\n", i, argv[i]);

      if(strcmp(argv[i], "-steps") == 0){
        steps = atoi(argv[i+1]);
      }      
      
      if(strcmp(argv[i], "-N") == 0){
        N = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-J") == 0){
        J = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-T") == 0){
        Temp = atof(argv[i+1]);
      }
  }

  printf("J:%d, Temp:%f, N:%d, steps:%d\n",J, Temp, N, steps);

  //Valores estatísticos
  int energy = 0;
  float magnetization;
  
  //Matrizes
  int *ising_matrix;
  int *hamiltonian_matrix;
  float *ising_rand;

  //Valores auxiliares
  bool white = false;
  curandState *ising_rand_states;

  //Valores para as Estruturas de Dados
  size_t bytes = N*N*sizeof(int);
  size_t bytes_float = N*N*sizeof(float);
  size_t bytes_curand = N*N*sizeof(curandState);
  
   // Verifique o retorno de cada chamada CUDA
  cudaError_t err;
  err = cudaMallocManaged(&ising_matrix,bytes);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }

  err = cudaMallocManaged(&ising_rand,bytes_float);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }

  err = cudaMallocManaged(&ising_rand_states,bytes_curand);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }

  /*err = cudaMallocManaged(&white_matrix,bytes_half);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }

  err = cudaMallocManaged(&black_matrix,bytes_half);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }*/

  err = cudaMallocManaged(&hamiltonian_matrix,bytes);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }

  cudaDeviceSynchronize();


  cout<<"Initializating ising matrix:"<<endl;
  init_matrix(ising_matrix,N,350);

  //init_white_black_matrix(ising_matrix, white_matrix, black_matrix, N);
  
  cudaDeviceSynchronize();

  dim3 blocos(1);
  dim3 threads_blocos(N,N);

  setup_kernel<<<blocos, threads_blocos>>>(ising_rand_states, N,time(NULL));
  generate_randoms<<<blocos, threads_blocos>>>(ising_rand_states, N, ising_rand);
   
  //Checkboard monte-carlo
  FILE *statistics;
  statistics = fopen("statistics.txt","w");
  for(int i=0; i<steps; i++){
    
    setup_kernel<<<blocos, threads_blocos>>>(ising_rand_states, N,time(NULL)+i+100);
    generate_randoms<<<blocos, threads_blocos>>>(ising_rand_states, N, ising_rand);
    //printf("first random: %f\n", ising_rand[0]);
    cudaDeviceSynchronize();
    update_spin_checkboard<<<blocos,threads_blocos>>>(ising_matrix,ising_rand, J, Temp, N, white);
    cudaDeviceSynchronize();
    update_spin_checkboard<<<blocos,threads_blocos>>>(ising_matrix,ising_rand, J, Temp, N, !white);
    cudaDeviceSynchronize();

    //Calculate hamiltonian
    //calculate_spin_hamiltonian<<<blocos,threads_blocos>>>(ising_matrix, N, hamiltonian_matrix);
    //cudaDeviceSynchronize();
    energy = calculate_hamiltonian2(ising_matrix,N,J);
    //energy = calculate_hamiltonian2(hami


    //Calculate magnetization
    magnetization = calculate_magnetization(ising_matrix,N);
    
    //Print results:
    //printf("---------------------------------\n");
    //print_lattice(ising_matrix,N);
    //printf("Energy:%d; Magnetization:%f\n",energy,magnetization);
    //printf("---------------------------------\n");
    
    char mensagem[50];
    sprintf(mensagem, "%d,%d,%.3f\n", i, energy, magnetization);
    fputs(mensagem, statistics);
    
  }
  fclose(statistics);

  cudaFree(ising_matrix);
  cudaFree(ising_rand);
  cudaFree(ising_rand_states);
  cudaFree(hamiltonian_matrix);
  
  cout<<"PROGRAM COMPLETED OK"<<endl;
}