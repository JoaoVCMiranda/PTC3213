# PTC3213 - Eletromagnetismo

Aqui estão algumas das nossas ideias e desenvolvimento acerca do assunto de Eletromagnetismo.

Vou documentar os meus aprendizados de MATLAB aqui nesse arquivo para o caso de vocês precisarem de alguma referência rápida depois.

```octave

```
## MATLAB/octave

O que é o MATLAB/octave, é uma linguagem de descrição matemática(high-level) usada para simulações numéricas, assim como o LEAN e outros.

O octave é uma alternativa opensource para o MATLAB, que é uma linguagem de programação proprietária, e eles compartilham entre si muitas similaridades, e aparentemente são intercambiáveis para o que vamos fazer

Os códigos não aparentam ter necessidade de serem compilados para executar(ao menos para mim usuário).
E podem ser rodados iterativamente, ou como foi descrito na referência, _"REPL(Read Evaluate Print Loop)"_

## Dicas de MATLAB

Para abrir a GUI do octave via terminal
```sh
octave --gui 
```

```octave
% Comentários em matlab começam com "%" tudo o que vier depois será comentado
```

```octave
% O pkg é um gerenciador de pacotes para códigos em matlab
% assim como o pip está para o python 

% Um pacote que precisa ser instalado para executarmos as simulações é o matgeom
pkg install -local -forge matgeom;
```
> Obs: lembrar de colocar o ponto e vírgula

```octave
% O ponto e vírgula é opcional, mas ele tem uma função.
% Esse atua como um supressor de saída
% Com ";" não há print do que está armazenado na variável em si
% Sem ";" sim.
% Exemplo
>> A = [1,2,3]
A =

   1   2   3
>> A = [1,2,3];
>> 
```
