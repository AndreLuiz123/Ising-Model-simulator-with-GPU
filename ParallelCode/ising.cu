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
  
  if(blockDim.x*blockIdx.x + threadIdx.x < N && blockDim.y*blockIdx.y + threadIdx.y < N){
    int current_index = (blockDim.y*blockIdx.y + threadIdx.y)* N + blockDim.x*blockIdx.x + threadIdx.x;
    int current_spin = lattice[current_index];
    int new_spin = lattice[current_index]*-1;
    int dE=0;
    int dE_current = 0; 
    int dE_new = 0;
    bool accept = false;

    if(white){
      if(((blockDim.x*blockIdx.x + threadIdx.x)%2 == 0 && (blockDim.y*blockIdx.y + threadIdx.y)%2 == 0) || ((blockDim.x*blockIdx.x + threadIdx.x)%2 != 0 && (blockDim.y*blockIdx.y + threadIdx.y)%2 != 0))
        accept = true;
    }else{
      if(((blockDim.x*blockIdx.x + threadIdx.x)%2 != 0 && (blockDim.y*blockIdx.y + threadIdx.y)%2 == 0) || ((blockDim.x*blockIdx.x + threadIdx.x)%2 == 0 && (blockDim.y*blockIdx.y + threadIdx.y)%2 != 0))
        accept = true;
    }

    if(accept){
      //Calculate delta E
      
        if(blockDim.x*blockIdx.x + threadIdx.x != 0){
          dE_new += new_spin*(lattice[current_index - 1]);
          dE_current += current_spin*(lattice[current_index - 1]);
        }
        if(blockDim.x*blockIdx.x + threadIdx.x  != N-1){
          dE_new += new_spin*(lattice[current_index + 1]);
          dE_current += current_spin*(lattice[current_index + 1]);
        }
        if(blockDim.y*blockIdx.y + threadIdx.y != 0){
          dE_new += new_spin*(lattice[current_index - N]);
          dE_current += current_spin*(lattice[current_index - N]);
        }
        if(blockDim.y*blockIdx.y + threadIdx.y  != N-1){
          dE_new += new_spin*(lattice[current_index + N]);
          dE_current += current_spin*(lattice[current_index + N]);   
        }
        
        dE = -1*J*(dE_new - dE_current);

        if(dE <= 0){
          lattice[current_index] = new_spin; 
        }else{
          float boltz_dist = -dE/Temp; 
          float probability = expf(boltz_dist);
          float alpha = fminf(1,probability);
          if(randoms[current_index] < alpha){
            lattice[current_index] = new_spin; 
          }
        }
        
        //lattice[current_index] = 0;
    }
  }
}

__global__ void update_spin_double_buffering(int *lattice, int *old_lattice, float *randoms, int J, float Temp,int N){
  if(threadIdx.x < N && threadIdx.y < N){
      int current_spin = old_lattice[threadIdx.y*N + threadIdx.x];
      int new_spin = old_lattice[threadIdx.y*N + threadIdx.x]*-1;
      int dE=0;
      int dE_current = 0; 
      int dE_new = 0;

      //Calculate delta E
      if(threadIdx.x != 0){
        dE_new += new_spin*(old_lattice[threadIdx.y*N + threadIdx.x - 1]);
        dE_current += current_spin*(old_lattice[threadIdx.y*N + threadIdx.x - 1]);
      }
      if(threadIdx.x != N-1){
        dE_new += new_spin*(old_lattice[threadIdx.y*N + threadIdx.x + 1]);
        dE_current += current_spin*(old_lattice[threadIdx.y*N + threadIdx.x + 1]);
      }
      if(threadIdx.y != 0){
        dE_new += new_spin*(old_lattice[(threadIdx.y-1)*N + threadIdx.x]);
        dE_current += current_spin*(old_lattice[(threadIdx.y-1)*N + threadIdx.x]);
      }
      if(threadIdx.y != N-1){
        dE_new += new_spin*(old_lattice[(threadIdx.y+1)*N + threadIdx.x]);
        dE_current += current_spin*(old_lattice[(threadIdx.y+1)*N + threadIdx.x]);   
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
}

__global__ void update_old_lattice(int *lattice, int *old_lattice, int N){

  if(threadIdx.x < N && threadIdx.y < N){
    old_lattice[threadIdx.y*N + threadIdx.x] = lattice[threadIdx.y*N + threadIdx.x];
  }

}

//Calcula os valores de cada spin hamiltoniano
__global__ void calculate_spin_hamiltonian(int *lattice, int N, int *hamiltonian){
  
  if(blockDim.x*blockIdx.x + threadIdx.x < N && blockDim.y*blockIdx.y + threadIdx.y < N){
    int index = (blockDim.y*blockIdx.y + threadIdx.y)* N + blockDim.x*blockIdx.x + threadIdx.x;
    int current_spin = lattice[index];
    hamiltonian[index] = 0;

      if(blockDim.x*blockIdx.x + threadIdx.x!=0)
        hamiltonian[index]+= lattice[index-1]*current_spin;
      if(blockDim.x*blockIdx.x + threadIdx.x!=N-1)
        hamiltonian[index]+= lattice[index+1]*current_spin;
      if(blockDim.y*blockIdx.y + threadIdx.y!=0)
        hamiltonian[index]+= lattice[index - N]*current_spin;
      if(blockDim.y*blockIdx.y + threadIdx.y!=N-1)
        hamiltonian[index]+= lattice[index + N]*current_spin;
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

//Soma de elementos do modelo ising e tira a média - paralelizar depois
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

int main(int argc, char *argv[]){
  
  printf("Número de argumentos: %d\n", argc);
  //Parâmetros do sistema
  int approach = 0;
  int proportion = 750;
  int steps = 1000;
  int J = 1;
  int N = 32;
  float Temp = 5;

  for (int i = 0; i < argc; i++) {
      printf("Argumento [%d]: %s\n", i, argv[i]);

      if(strcmp(argv[i], "-approach") == 0){
        approach = atoi(argv[i+1]);
      }  

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
        J = atoi(argv[i+1]);
      }

      if(strcmp(argv[i], "-T") == 0){
        Temp = atof(argv[i+1]);
      }
  }

  printf("approach:%d(0:checkboard; 1:doublebuffering), proportion: %d, J:%d, Temp:%f, N:%d, steps:%d\n",approach, proportion, J, Temp, N, steps);

  //Valores estatísticos
  int energy = 0;
  float magnetization;
  
  //Matrizes
  int *ising_matrix;
  int *ising_matrix_old;
  int *hamiltonian_matrix;
  float *ising_rand;

  //Valores auxiliares
  bool white = true;
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

  err = cudaMallocManaged(&ising_matrix_old,bytes);
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

  err = cudaMallocManaged(&hamiltonian_matrix,bytes);
  if (err != cudaSuccess) {
    cout << "Failed to allocate memory for ising_matrix: " << cudaGetErrorString(err) << endl;
    return -1;
  }

  cudaDeviceSynchronize();


  cout<<"Initializating ising matrix:"<<endl;
  init_matrix(ising_matrix,N,proportion);
  
  cudaDeviceSynchronize();

  dim3 blocos(2,2);
  dim3 threads_blocos(N/2,N/2);

  setup_kernel<<<blocos, threads_blocos>>>(ising_rand_states, N,time(NULL));
  generate_randoms<<<blocos, threads_blocos>>>(ising_rand_states, N, ising_rand);
   
  clock_t t,t_total;
  t_total = clock();
  //Checkboard monte-carlo
  FILE *statistics;
  statistics = fopen("statistics.txt","w");
  t= clock();
  if(approach == 0){
    //Checkboard
    print_lattice(ising_matrix,N);
    for(int i=0; i<steps; i++){
      
      //Ising step
      setup_kernel<<<blocos, threads_blocos>>>(ising_rand_states, N,time(NULL)+i+100);
      generate_randoms<<<blocos, threads_blocos>>>(ising_rand_states, N, ising_rand);
      
      cudaDeviceSynchronize();
      update_spin_checkboard<<<blocos,threads_blocos>>>(ising_matrix,ising_rand, J, Temp, N, white);
      
      cudaDeviceSynchronize();
      update_spin_checkboard<<<blocos,threads_blocos>>>(ising_matrix,ising_rand, J, Temp, N, !white);
      cudaDeviceSynchronize();
            

      //Calculate energy
      calculate_spin_hamiltonian<<<blocos,threads_blocos>>>(ising_matrix,N,hamiltonian_matrix);
      cudaDeviceSynchronize();
      energy = calculate_hamiltonian(hamiltonian_matrix,N,J);    
      
      //Calculate magnetization
      magnetization = calculate_magnetization(ising_matrix,N);
      
      char mensagem[50];
      sprintf(mensagem, "%d,%d,%.3f\n", i, energy, magnetization);
      fputs(mensagem, statistics);
      
    }
  }
  
  t = clock() - t;
  printf("Tempo isingStep:%f\n",(float)t/CLOCKS_PER_SEC);
  fclose(statistics);

  print_lattice(ising_matrix,N);

  cudaFree(ising_matrix);
  cudaFree(ising_matrix_old);
  cudaFree(ising_rand);
  cudaFree(ising_rand_states);
  cudaFree(hamiltonian_matrix);
  
  t_total = clock() - t_total;
  printf("Tempo total:%f\n",(float)t_total/CLOCKS_PER_SEC);
  cout<<"PROGRAM COMPLETED OK"<<endl;
}