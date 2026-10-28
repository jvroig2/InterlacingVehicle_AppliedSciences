function [Reg_sep]=regdiagsep(Reg,n)
% function [Reg_sep]=regdiagsep(Reg,n)
%Función para separar los estados del sistema diagonal en subespacio
%1º) Se obtiene la forma diagonal(modal)
%2º) Se ordenan los estados de más rápido a menos rápido
%3º) Se separan en los subespacio que indica n=[nr nl1 nl2 ...]
%Cada subespacio es independiente del resto
% tendrá como entradas el resto de subespacios
%Reg: Sistema discreto a separar (tipo ss)
%n: vector con la dimensión de cada subespacio
%Reg_sep. estructura con todos los subsistemas Red_sep.Sx (con x=1,2,...)
%NOTA: La matriz D irá al 1º subespacio

nT=size(Reg.A,1);%Número total de estados.

if sum(n) ~= nT
    disp('La dimensión de los subespacio no coincide con el número de estados del sistema')
    return;
end
m=size(Reg.B,2);%Número de entradas
sal=size(Reg.C,1);%Número de salidas

p=length(n);%Número de partes a separar
T=Reg.Ts;%Periodo de muestreo base

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RegD=canon(Reg,'modal'); %Forma modal diagonal
%Esto debería hacerse fuera elegir las dimensiones de los subespacios y no
%separar polos complejos
%Falta ordenar polos diagonales de más a menos rápido.
polos=pole(Reg);
[~, idx] = sort(abs(polos), 'ascend');
%Definir Tx para RegD que sea un cambio de orden
Tx=[];
for i=1:length(idx)
    Tx=[Tx; zeros(1,idx(i)-1) 1 zeros(1,length(idx)-idx(i))];
end
% Reordenar matriz diagonal de más rápido a más lento
RegDo=ss2ss(RegD,Tx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:p
    if k==1
        f=1:n(k);
    else
        f=sum(n(1:k-1))+1:sum(n(1:k));%Filas del bloque k
    end;
    A=RegDo.A(f,f);
    B=RegDo.B(f,:);
    C=RegDo.C(:,f);
    if k==1
        D=RegDo.D;%D a la parte rápida
    else
        D=zeros(sal,m);
    end;
    eval(['Reg_sep.S' num2str(k) '= ss(A,B,C,D,T);'])
end