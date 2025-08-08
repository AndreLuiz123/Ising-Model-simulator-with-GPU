# Nome do executável
TARGET = main

# Arquivo fonte
SRC = ising.c

# Compilador C
CC = gcc

# Flags de compilação
CFLAGS = -Wall -O2

# Bibliotecas

LIBS = -lm

# Regra padrão
all: $(TARGET)


# Como compilar
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LIBS)

# Limpar arquivos gerados
clean:
	rm -f $(TARGET)
