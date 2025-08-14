### Implementações paralelas

### Como executar
Compile o binário com o comando `make`.

Rode `./main`. Será gerado um arquivo `./statistics.txt` com o valor da energia e da magnetização do sistema ao longo do tempo. Esses dados podem ser melhor visualizados ao executar o script `./statistics.py`.

É possível adicionar parâmetros para o `./main`:
- `-approach` : Determina se o script irá executar o checkboard algorithm(`-approach 0`) ou o double-buffering(`-approach 0`);
- `-steps` : Determina em quantos ciclos o simulador irá atuar. Por exemplo, `-steps 1000` indica que o simulador demonstra a evolução de um modelo Ising ao longo de 1000 iterações;
- `-T` : Determina a temperatura do sistema;
- `-N` : Determina as dimensões do modelo Ising a ser gerado;
- `-proportion` : Número entre 0 e 1000 que determina a razão entre spins 1 e -1 do estado inicial da malha de spins. Quanto mais próximo de 1000, mais spins iguais a 1 aparecem;
- `-J` : Energia entre spins do sistema


### Observações
O tratamento das entradas ainda precisa ser aprimorado. Sendo assim, até o momento, sempre adicione um número após a chamada dos parâmetros do `./main`. Preferencialmente, utilize números múltiplos de 32 para o parâmetro `-N`.
