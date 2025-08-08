# Nome do executável
TARGET = main

# Arquivo fonte
SRC = ising.cu

# Compilador CUDA
NVCC = nvcc

# Flags de compilação - caso seu cuda seja diferente do sm_86, trocar
NVCC_FLAGS = -arch=sm_86

# Regra padrão
all: $(TARGET)

# Como compilar
$(TARGET): $(SRC)
	$(NVCC) $(NVCC_FLAGS) -o $(TARGET) $(SRC)

# Limpar arquivos gerados
clean:
	rm -f $(TARGET)
