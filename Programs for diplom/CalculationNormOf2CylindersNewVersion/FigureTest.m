%fs=500;
%t = -20:1/fs:20;

%x = -20:1/fs:20;
%w_0=1.1:4.6;
%M = zeros(1,length(w_0)); 
%i=1;
%for %w_0=1.1 :4.6
    %R=3+i;
    %M(i)=R;
    %i=i+1;
%end
%5plot(t,x);

%%%%% Эта программа позовляет вычислять норму
%%%%% п

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters

n = 2000;
a = -2000;  
b = 2000;
cr = -2000;
d = 2000;
h = (b-a) / n; 
x=a:h:b;
y=cr:h:d;
L = 3 * a_0;
Outside = 0;
N_max = 0;
N_cylinders = 2;
p = 4;
wR = 1.1*w_H:1.1*w_H:4* w_H;
M = zeros(1,length(wR)); 
wS=1*w_H:w_H:2* w_H;
i=1;
for wm=1.1*w_H:w_H:4* w_H
    BD = crystalCoefficients2D_eigenwave(N_max, wm, a_0, cylXY, k_0, EE, GG, HH, c, p);
    RateCalculationMatrix1 = RateCalculationMatrix(N_max, k_0, p, a_0, EE, GG, HH , BD, N_cylinders );
    %RateCalculationMatrixOutInt = RateCalculationMatrixOutInt(N_max, k_0, p, a_0, EE, GG, HH , BD, N_cylinders );
end
