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

```octave
% Condicionais
% Os delimitadores no octave dos laços são os próprios operadores
% Ao invés de termos:
% if(condicao){ //entre delimitadores }

% temos
if(A == B) % condicao <- {A == B}
% delimitado pelos operadores
end
% O mesmo vale para os laços!
for % a forma de fazer laços será vista a seguir
% codigo
end
```

```octave
% Laços!!
for k=1:size()
% o loop for, está utilizando como variável de iteração k
% se a sintaxe for com é na linguagem c
% está fazendo k=1 e passar por todos os valores que estão dentro do retorno dessa função size.
end
```


```octave
% Para fazer uma lista em octave é simples!

1:x % faz uma lista que começa em 1 e vai até x

1:x - 1 % faz uma lista que começa em 1 e vai até x - 1 
(1:x) -1 % faz uma lista que vai de 0 até x-1, ou seja remove 1 de cada item da lista
```

```octave
% é uma função um pouco estranha...
% primeiro, se você executar ela com duas arrays quaisquer X e Y
% ela repete a array X, y vezes
meshgrid(X, Y);
% Porém se você chamar assim
[a,b] = meshgrid(X,Y)
% a variável "a" recebe X, y vezes
% e a variável "b" recebe Y, x vezes, mas na vertical
% exemplo:
[u,v] = meshgrid((1:2), (1:3))
%u =
%   1   2
%   1   2
%   1   2
%v =
%   1   1
%   2   2
%   3   3
% curioso não?
```

```octave
% o octave tem outra coisa interessante
% * as arrays são indexadas no 1. e não no 0 como é de costume
% se fizermos a seguinte array2d/matriz
u = [0 0;1 0;1 1;0 1;0 0];
% para acessar os seus items dessa forma
u(X)
% vou fazer uma função explicativa, u(x) = x
u = [1 6;2 7;3 8;4 9;5 10]
% agora se você quiser acessar a array com dois parâmetros
u(y,x)
% farei outra função explicativa tal que u(y,x) = "xy" x concatenado com y
u = [1 2;1 2;1 2;1 2;1 2]
```

## Física

### E agora, como discretizar na prática a eq de laplace?

Veja que já temos a função potencial phi...

Pela definição de derivada temos que

$$
\frac{\partial \Phi(x,y)}{\partial x} = lim_{h\rightarrow 0} \frac{\Phi (x+h,y) - \Phi (x,y)}{h}
$$
Como não podemos colocar esse limite na prática, o mais próximo que temos é o dx, então teremos

a diferença entre dois valores de potencial adjacentes, em relação a x dividido por dx...


