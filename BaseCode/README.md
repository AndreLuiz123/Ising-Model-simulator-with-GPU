### Implementações não paralelas

### Como executar
Compile o binário com o comando `make`.

Rode `./main`. Será gerado um arquivo `./statistics.txt` com o valor da energia e da magnetização do sistema ao longo do tempo. Esses dados podem ser melhor visualizados ao executar o script `statistics.py`.

É possível adicionar parâmetros para o `./main`:
- `-mcs` : Determina quantos flips cada spin do sistema pode potencialmente realizar. Se o parâmetro for `./main -mcs 1000`, por exemplo, quer dizer que os passos Monte Carlo serão executados `<número de spins>*1000` vezes, o que permite que cada spin seja escolhido 1000 vezes para ser flipado durante a simulação;
- `-T` : Determina a temperatura do sistema;
- `-M` : Determina a dimensão horizontal do modelo Ising a ser gerado;
- `-N` : Determina a dimensão vertical do modelo Ising a ser gerado;
- `-proportion` : Número entre 0 e 1000 que determina a razão entre spins 1 e -1 do estado inicial da malha de spins. Quanto mais próximo de 1000, mais spins iguais a 1 aparecem;
- `-J` : Energia entre spins do sistema
- `-B` : Determina condição de borda. Para `-B 0` a condição de borda é a clássica, para `-B 1` a condição de borda é a de torus

### Observações
- O tratamento das entradas ainda precisa ser aprimorado. Sendo assim, até o momento, sempre adicione um número após a chamada dos parâmetros do ./main.
- Para executar o script `statistics.py`, sugiro criar uma venv python com o comando `python -m venv minhas_estatisticas`, acessar a venv com o comando `source minhas_estatisticas/bin/activate`, e instalar o matplot lib na venv com os comandos `pip install matplotlib` e `pip install PyQt6` (ou `pip install PyQt5`). Depois de acessar as estatísticas, basta usar o comando `deactivate` no seu terminal. Depois da primeira vez executando esses comandos, para acessar as estatísticas novamente basta usar o comando `source minhas_estatisticas/bin/activate` de novo.
