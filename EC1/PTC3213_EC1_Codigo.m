% ///////////////////////////////////////////////////////////////////////////
%			PTC3213 - EC1 - 2024 - DIFERENÇAS FINITAS

%    #14754306 Arthur Soresini Bedin
%    #14582927 João Victor Cavalcante Miranda
%    #NUSP3 Marcos Gabriel
%			Complete os campos com [Preencher Aqui]
%
%  /////////////////////////////////////////////////////////////////////////
clear;
clf;
warning ("off");
pkg install -local -forge  matgeom; 
pkg load matgeom;
clc;

%   Dados de entrada

NUSP = 14654306 ; % NUSP do 1o aluno (ordem alfab.)

a=   11; 
b=   5 ;
c=   4 ;
d=   b-3; % 2
g=   3 ;

h=(b-d)/2;

epsr=   2.5   ; % adimensional
sigma=    2.5   ;  % S/m
sigma_dual=   3.5   ; % S/m


eps0=  8.85418717*10^(-12) ;  % F/m

Vmin=  0 ;     % Volts, de acordo com o enunciando.pdf
Vmax=  100 ;   % Volts, idem

%            Definicao do dominio
% A variavel dx abaixo é a discretizacao utilizada. Valores diferentes
% daqueles sugeridos abaixo nao funcionarao. Diminua o dx para gerar a 
%versao final a ser entregue. Ou seja, aumentar a resolução da simulação

% Definições para o cálculo usando o MDF

dx=0.5;   %% Sugestão do Prof: Mude para dx=0.25 somente quando for gerar os resultados finais!!!
erro=0.0;
start=start_Dual= (Vmin + Vmax)/2; % um bom palpite é sempre a média, a partir daí melhoraremos iterativamente.
iter=0;
dy=dx;
lx=a;
ly=b;
% "Resolução" da tela, quantidade de quadradinhos que devemos ter em cada dimensão
Nx=round(lx/dx)+1;
Ny=round(ly/dy)+1; 
% Geometria do problema

% São retangulos, mas tem 5 pontos pois o último é igual ao primeiro para fechar o anel
% Os aneis vão em direções contrárias
% Os aneis são os "vertices dos retângulos um dentro do outro"
ring1= [0 0; lx 0; lx ly; 0 ly; 0 0];
ring2=[g h; g h+d; g+c h+d; g+c h; g h];

%novo anel, onde serão realizados os cálculos
ring3 = [0+1 0+1; lx-1 0+1; lx-1 ly-1; 0+1 ly-1; 0+1 0+1];

% array de arrays
polyg={ring1,ring2};
% É a união dos dois poligonos
% simplesmente as duas arrays ring1 e ring2 concatenadas
verts = polygonVertices(polyg);
% sintaxe para fazer uma grid
xgv=((1:Nx)-1)*dx;
% faz uma sequência de 0 até Nx-1 e multiplica por dx, para assim ter todos os pontos em x até
ygv=((1:Ny)-1)*dx;
% tipo um produto cartesiano
[x,y]=meshgrid(xgv,ygv);
% O x tem todos os valores de x para cada ponto no conjunto xgv cartesiano ygv

% aparentemente essas linhas são redundantes, pois ring1 e verts1 contém exatamente a mesma informação, talvez pudéssemos chamar ringN de vertsN sem maiores problemas
verts1 = polygonVertices(ring1);
verts2 = polygonVertices(ring2);


xv1=verts1(:,1);
yv1=verts1(:,2);
xv2=verts2(:,1);
yv2=verts2(:,2);


% vertices do polígono de calculo
verts3 = polygonVertices(ring3);
xv3 = verts3(:,1);
yv3 = verts3(:,2);
[in3,on3] = inpolygon(x,y,xv3,yv3);
% Só para deixar logo tudo junto


% Essa função gera dois valores, in1 e on1, que são arrays booleanas que seguem o seguinte padrão:
% in1 é 1 para todos os valores do espaço (x,y) dado que estiver dentro do poligono de vertices xv1 e yv1, e 0 caso contrário. Como uma área hachurada
% on1 é a borda do polígono. Como um perímetro.
[in1,on1] = inpolygon(x,y,xv1,yv1);

[in2,on2] = inpolygon(x,y,xv2,yv2);

% Atribui Condicoes de contorno

r=find(in1&~in2|on2); % tudo

p=find(in1&~on1&~in2); %so  nos internos = tudo que está na placa e não está na borda de fora, mas pega a borda de dentro

q=find(on1|on2); %so fronteira = só os dois perimetros, um retangulo dentro do outro

iVmax=find(on2);

% Dica do prof: pode ser usado ! no lugar de ~ como uma extensão do octave.
iFuro=find(in2&~on2);

% inicialização da função potencial
Phi_prev=zeros(size(x));
Phi_new=zeros(size(x));


% a borda exterior, definimos que Vmax=100, e isso não será mais alterado
Phi_new(iVmax)= Vmax;
Phi_new(iFuro)= NaN;

% Um jeito interessantíssimo de definir os valores dentro de uma variável...
% da mesma forma que é feito para contra domínio
% Para todo ponto no domínio p, atribui-se o valor start à imagem
Phi_new(p)= start;

% Contador de iteracoes
iter=0;

% Erro maximo entre Phi_new e Phi_prev
% Phi_new - Phi_prev é uma subtração de matrizes
% abs apenas tira o valor absoluto de cada um deles
% max vai pegar o máximo de cada coluna
% o segundo max vai pegar o máximo das colunas
% Ou seja, estamos calculando o máximo erro em qualquer que seja o ponto
erro=max(max(abs(Phi_new-Phi_prev)));
 
%            Laco iterativo - Metodo das Diferencas Finitas
% agora começa a brincadeira

while(erro > 1e-4 && iter < 1e4)% Executa ate convergir ou atingir o maximo de iteracoes
    iter=iter+1; % Incrementa iteracao

%	Atualiza o potencial dos nos internos pela media dos 4 vizinhos - Eq. Laplace - M.D.F.
%	Funciona especialmente bem para a Eq. de laplace, pois a ideia que é consequência dessa é
%	O cálculo do laplaciano
%	Que é o somatório da segunda derivada do potencial em cada dimensão
%	E segundas derivadas tem tudo a ver com médias(segundo o Richard P. Feynman)

%	Esse laço é apenas passando por todos os elementos em Phi
    for k=1:size(p,1);
        [i,j]=ind2sub(size(x),p(k));
            Phi_new(i,j)=(Phi_new(i-1,j)+Phi_new(i+1,j)+Phi_new(i,j-1)+Phi_new(i,j+1))/4;
    end
% Calcula maximo erro entre Phi_atual e Phi_prev de todo o dominio

    erro=max(max(abs(Phi_new-Phi_prev)));
    % coloca na conta do eps em quanto que está o erro agora.
    eps(iter)=erro;

%    Atualiza a matriz de potenciais
    Phi_prev=Phi_new;
end

niter1=iter;

if (niter1 == 1e4 && erro > 1e-4)
	disp([' Numero maximo de iteracoes atingido sem convergencia :', num2stg(niter1), '  iteracoes \? Erro: \n', num2str(erro), 'Os resultados podem nao ter significado!\n']);
end

% Para calcular as derivadas parciais... Tem um jeitinho...

% Define-se as funções dos potenciais...
campoEx = zeros(size(x));
campoEy = zeros(size(y));

campoEx(iFuro)= NaN;
campoEx(on1) = NaN;
campoEy(iFuro)= NaN;
campoEy(on1) = NaN;



% p são as "coordenadas" em uma dimensão de cada ponto válido da tela discretizada...

% para cada elemento em p, começando em 1
for k=1:size(p,1)
	% pegar as coordenadas, em uma dela do formado dado por size(x) na coordenada p(k)
	[i,j] = ind2sub(size(x), p(k));
	% parece que está trocado, mas é assim mesmo pois é indexado na coluna e depois na linha.
	campoEx(i,j) = ( (Phi_new(i,j)-Phi_new(i,j+1)) + (Phi_new(i+1,j)-Phi_new(i+1,j+1)) )/(2*dx);
	campoEy(i,j) = ( (Phi_new(i,j)-Phi_new(i+1,j)) + (Phi_new(i,j+1)-Phi_new(i+1,j+1)) )/(2*dy);
end


% agora que temos os campos
% Como pela geometria do problema os campos ficarão sempre na mesma direção da superfície podemos tirar o valor absoluto.
Sigma = 0;
for k=1:size(on3,1)
    [i,j]=ind2sub(size(x),find(on3)(k));
    if(i==ring3(1,1)|i==ring3(3,1))
      Sigma += abs(campoEy(i,j))*dy;
    end
    if(j==ring3(1,2)|j==ring3(3,2))
      Sigma += abs(campoEx(i,j))*dx;
    end
end
% Atenção às unidades
I = sigma*1*Sigma/1000; % [sigma] = (kohm*m)^-1, 1 m ,[Sigma] = Volts % [I] = A
R = (Vmax-Vmin)/I; % ohm
Cap =  (1e15*eps0*epsr)/(R*sigma); % [C] = Farad

% dx=.5 ->    I = 50.69
% dx=.25 ->   I =  0.76
% dx=.125 ->  I =  3.42
% dx=.0625 -> I =  0.00




% Problema Dual (Somente para tracado dos Quadrados Curvilineos!)
% Atribui Condicoes de Contorno

iyDual=find( (y(:,:) < ly/1.999) & (y(:,:) > ly/2.001) );
iVmaxdual=find( (x(iyDual) > (-0.01)) & (x(iyDual) < (1.0001*g)));
i0=find( (x(iyDual)> (0.9999*(g+c))) & (x(iyDual)< (1.0001*lx)) );
xfe=find(  x(iVmax)< 1.0001*min(x(iVmax)) ); xfd=find(  x(iVmax)> 0.9999*max(x(iVmax)) );
yfa=find(  y(iVmax)> 0.9999*max(y(iVmax)) ); yfb=find(  y(iVmax)< 1.0001*min(y(iVmax)) );
tol=1e-4;
for k=1:size(iVmax,1);
    if ( abs( x(iVmax(k))-min(x(iVmax)) )< tol && abs( y(iVmax(k))-min(y(iVmax)) )< tol)
             [ieb,jeb]=ind2sub(size(x), iVmax(k));
     elseif (abs( x(iVmax(k))-min(x(iVmax)) )< tol && abs( y(iVmax(k))-max(y(iVmax)) )< tol)
            [iea,jea]=ind2sub(size(x), iVmax(k));
     elseif ( abs( x(iVmax(k))-max(x(iVmax)) )< tol && abs( y(iVmax(k))-min(y(iVmax)) )< tol)
             [idb,jdb]=ind2sub(size(x), iVmax(k));
     elseif (abs( x(iVmax(k))-max(x(iVmax)) )< tol && abs( y(iVmax(k))-max(y(iVmax)) )< tol)
            [ida,jda]=ind2sub(size(x), iVmax(k));
    end
 end
  
Dual_prev=zeros(size(x));
Dual_new=Dual_prev;
Dual_new(r)= -1;
Dual_new(iFuro)= NaN;
Dual_new(iyDual(iVmaxdual))=Vmax;
Dual_new(iyDual(i0))=Vmin;
p2=find(Dual_new(p) < 0);
Dual_new(r)= start_Dual;
Dual_new(iFuro)= NaN;
Dual_new(iyDual(iVmaxdual))=Vmax;
Dual_new(iyDual(i0))=Vmin;

% Contador de iteracoes - dual
iter2=0;

% Erro maximo entre Phi_new e Phi_prev (Dual)
erro2=max(max(abs(Dual_new-Dual_prev)));

%       Laco iterativo (Problema Dual) - MDF(Não é a placa de madeira, pode significar Método das Diferenças Finitas)

while(erro2 > 1e-3 && iter2 < 1e4)% Executa ate convergir ou atingir o maximo de iteracoes
    iter2=iter2+1; % Incrementa iteracao
%	Atualiza o potencial das fronteiras
%	Da mesma forma que foi feito com a matrix de potenciais ?

    Dual_new(1,:)=Dual_prev(2,:);
    Dual_new(Ny,:)=Dual_prev(Ny-1,:);
    Dual_new(:,1)=Dual_prev(:,2);
    Dual_new(2:Ny-1,Nx)=Dual_prev(2:Ny-1,Nx-1);   
    for k=2:size(xfe,1)-1
        [ie,je]=ind2sub(size(Dual_new), iVmax(xfe(k)));
        Dual_new(ie,je)=Dual_new(ie,je-1);
    end
    for k=2:size(xfd,1)-1
        [id,jd]=ind2sub(size(Dual_new), iVmax(xfd(k)));
        Dual_new(id,jd)=Dual_new(id,jd+1)
    end
    for k=2:size(yfb,1)-1
        [ib,jb]=ind2sub(size(Dual_new), iVmax(yfb(k)));
        Dual_new(ib,jb)=Dual_new(ib-1,jb);
    end
    for k=2:size(yfa,1)-1
        [ia,ja]=ind2sub(size(Dual_new), iVmax(yfa(k)));
        Dual_new(ia,ja)=Dual_new(ia+1,ja);
    end
    Dual_new(iyDual(iVmaxdual))=Vmax;
    Dual_new(iyDual(i0))=Vmin;

% Atualiza o potencial dos nos internos pela media dos 4 vizinhos - Eq. Laplace - M.D.F.

    for k=1:size(p2,1);
        [i,j]=ind2sub(size(x),p(p2(k)));
        Dual_new(i,j)=(Dual_new(i-1,j)+Dual_new(i+1,j)+Dual_new(i,j-1)+Dual_new(i,j+1))/4;
    end

% Cantos não-gregorianos (efeito das pontas)

    Dual_new(ieb,jeb)=(Dual_new(ieb-1,jeb)+Dual_new(ieb+1,jeb)+Dual_new(ieb,jeb-1)+Dual_new(ieb,jeb+1))/4;
    Dual_new(iea,jea)=(Dual_new(iea-1,jea)+Dual_new(iea+1,jea)+Dual_new(iea,jea-1)+Dual_new(iea,jea+1))/4;
    Dual_new(idb,jdb)=(Dual_new(idb-1,jdb)+Dual_new(idb+1,jdb)+Dual_new(idb,jdb-1)+Dual_new(idb,jdb+1))/4;
    Dual_new(ida,jda)=(Dual_new(ida-1,jda)+Dual_new(ida+1,jda)+Dual_new(ida,jda-1)+Dual_new(ida,jda+1))/4;

% Calcula maximo erro entre Phi_atual e Phi_prev de todo o dominio
    erro2=max(max(abs(Dual_new-Dual_prev)));
    eps2(iter2)=erro2;
% Atualiza a matriz de potenciais
    Dual_prev=Dual_new;
end
niter2=iter2;
if (niter2 == 1e4 && erro2 > 1e-3)
	disp([' Numero maximo de iteracoes atingido sem convergencia :', num2stg(niter2), '  iteracoes \? Erro: \n', num2str(erro2), 'Interprete este resultado com ressalvas!\n']);
end

% Dados de Saída - São coisas que eu irei calcular


%      CORRENTE TOTAL (A)
Somat=sum(Phi_new(2,:))+sum(Phi_new(Ny-1,:))+sum(Phi_new(:,2))+sum(Phi_new(:,Nx-1));

I = I;

%       RESISTENCIA em ohms
R=  R;

%        CAPACITANCIA em pF
Cap=  Cap;

%     RESISTENCIA DUAL em ohms
%Rdual=  [Preencher Aqui]????????????????  ;

%    VETOR DESLOCAMENTO
Dn=[Phi_new(2,1:Nx-1),Phi_new(1:Ny-1,Nx-1)',Phi_new(Ny-1,1:Nx-1),Phi_new(1:Ny-1,2)']*epsr*eps0/dx*100;

%   Densidade de carga mínima em nC/m^2
%Rho_s_min=  ;

%  Numero de tubos de corrente
nsnp= 19;
% Correção
ntubos=10/nsnp;

%              IMPRESSAO DE RESULTADOS NO TERMINAL 
%                  ATENCAO para as unidades:
%          R e Rdual em ohms     Cap em pF    Rho_s  em nC/m^2
%fprintf('\n\n nUSP: %d\n R = %g ohms\n C = %g pF\n Rho_s_min = %g nC/m^2\n Rdual = %g ohms\n Tubos: %g\n', NUSP, R, Cap, Rho_s_min,Rdual,floor(ntubos) );
fprintf('\n\n nUSP: %d\n R = %g ohms\n C = %g pF\n', NUSP, R, Cap);
FIG=figure (1);

%           TRACADO DE EQUIPOTENCIAIS

V=0:10:Vmax;
colormap cool;
[C,H]=contour(x,y,Phi_new,V);
clabel(C,V);
axis('equal');
hold on

%   EQUIPOTENCIAIS PROBLEMA DUAL (para tracado dos quadrados curvilineos)

deltaV=   10;
V=0:deltaV:Vmax;
colormap jet;
contour(x,y,Dual_new , V );
axis('equal');
strusp=sprintf('%d',NUSP);
titulo=['Mapa de Quadrados Curvilineos (EC1 2024) - ', strusp, ' - ', date()];
title(titulo);
hold off

%      ARQUIVO DE SAIDA COM O MAPA DOS QUADRADOS CURVILINEOS
%(Grava na pasta exibida no Navegador de Arq. da interface gráfica do Octave)

arq=['EC1_2024_QC_',strusp,'.png'];
print(FIG,arq);

%% Foi divertido! A partir de dois "%" o meu interpretador de linguagem faz um outro tipo de highlight, por isso removi todos para o que parecia mais simples para mim. Obrigado pelo código, foi um desafio proveitoso
