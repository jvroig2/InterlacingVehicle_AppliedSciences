%Revista Applied Sciences (MDPI)
%Comparativa de diferentes criterios para ordenar las dinámicas
%Rapidez, Ganancia estática, Residuo, HVS, HVS en coordenadas modales
%Aplicación: Vehiculo
% salida yaw rate entrada steering
% El control está diseñado para considerar una velocidad lineal entre 4 y 20 m/s (circuito urbano)
%Regulador diseñado en continuo con H_inf más agresivo para modelo vehículo de Zhang
%El regulador discreto se separa en 9 Dinámicas: 4 complejos + 5 reales
%En la forma digonal son sistemas de S1 a S9
%Analisis de los modos del regulador

clear
%Perido de muestreo
T=10e-3;%10ms

load('regmu_cochezhang')

%K: regulador continuo de orden 13


Gr=c2d(K,T,'zoh');
Reg=ss(Gr);
SR_Diag=canon(Reg,'modal')
[pole(SR_Diag) abs(pole(SR_Diag))] %No Ordenados de rápido a lento

%Ordenar polos de mayor a menor rapidez (modulo pi)
[RegDiag_sep]=regdiagsep(SR_Diag, [1 2 2 2 1 1 1 2 1]);%separar todos los subespacios
S1=RegDiag_sep.S1;%real
p1=pole(S1)%p1=pole(SR_Diag.A(1,1))
S2=RegDiag_sep.S2;%complejo
p2p3=pole(S2)%=eig(SR_Diag.A(2:3,2:3))
S3=RegDiag_sep.S3;%complejo
p4p5=pole(S3)%p4p5=eig(SR_Diag.A(4:5,4:5))
S4=RegDiag_sep.S4;%complejo
p6p7=pole(S4)%p6p7=eig(SR_Diag.A(6:7,6:7))
S5=RegDiag_sep.S5;%real
p8=pole(S5)%p8=eig(SR_Diag.A(8,8))
S6=RegDiag_sep.S6;%complejo
p9=pole(S6)%p9p10=eig(SR_Diag.A(9:10,9:10))
S7=RegDiag_sep.S7;%real
p10=pole(S7)%p11=eig(SR_Diag.A(11,11))
S8=RegDiag_sep.S8;%real
p11p12=pole(S8)%p12=eig(SR_Diag.A(12,12))
S9=RegDiag_sep.S9;%real
p13=pole(S9)%p13=eig(SR_Diag.A(13,13))
%9 Dinámicas: 4 complejos + 5 reales
%En la forma digonal son sistemas de S1 a S9

ns=9,%nº de modos
%%%Implemetar código para obtener ns a partir de Reg

p=[];
%Construcción de vectores para análisis de modos
for i=1:ns
    nomSist = sprintf('S%d', i);   % 'S1', 'S2', ...
    Si = evalin('base', nomSist);  % Recupera la variable S1, S2, ...
    pols = pole(Si);
    p(i) = pols(1);%Polo del modo
    speed(i)=abs(p(i)); %Más rápido a menor valor
    gain(i)=norm(dcgain(Si));
    residuo(i)=norm(Si.C*Si.B);
    [Sibal sigma]=balreal(Si);
    HVS(i)=sigma(1); 
end

%%%Automatizar la construcción de Sd
%%Integrar Sd dentro de regdiagsep()
Ad=blkdiag(S1.A,S2.A,S3.A,S4.A,S5.A,S6.A,S7.A,S8.A,S9.A);
Bd=[S1.B;S2.B;S3.B;S4.B;S5.B;S6.B;S7.B;S8.B;S9.B];
Cd=[S1.C S2.C S3.C S4.C S5.C S6.C S7.C S8.C S9.C];
Sd=ss(Ad,Bd,Cd,Gr.D,T);
[sysbal,sigma,TL,TR] = balreal(Sd);
%TL matriz de transformación para pasar a sysbal 

% sigma: HSV globales, ordenados de mayor a menor
% TR: coordenadas balanceadas -> coordenadas modales
% blks: dimensiones de los bloques modales
blks = [1 2 2 2 1 1 1 2 1];%Dimensiones de los bloques modales

hbal = zeros(numel(blks),1);
ini  = 1;
for i = 1:numel(blks)
    idx = ini:(ini + blks(i) - 1);

    alpha = sum(abs(TR(idx,:)).^2,1) ./ sum(abs(TR).^2,1);

    hbal(i) = alpha * sigma(:);

    ini = ini + blks(i);
end

%Mostrar resultados
[p' speed' gain' residuo' HVS' hbal]



