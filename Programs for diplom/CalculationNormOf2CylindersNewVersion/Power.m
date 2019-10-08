%%%%% Эта программа считает мощность
%%%%% 

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters


p=4; 
N_max = 0;
[BD] = crystalCoefficients2D_eigenwave(N_max, w_0, a_0, cylXY, k_0, EE, GG, HH, c, p);
E_usually = ElectricField(N_max, BD, k_0, p, a_0, EE, GG, HH);% + ElectricFieldEphi(N_max, BD, k_0, p, a_0, EE, GG, HH);
%IntegralFull = -1271246357504.24 - 967098839412.943i; %p = 4; 
IntegralFull = 8702974.82643714 + 2133155041.75737i; %p = 4; 
%IntegralFull = 24178992888100.8 - 2411885172.30695i; %p=8.298
%IntegralFull = 3014528.84134231 - 812540412.909618i %p=8.298
%Задаем  J
d = 2;
I_0 = 4;
j = I_0 / d;

E_unussualy = ElectricField(N_max, BD, k_0, -p, a_0, EE, -GG, HH); %ElectricFieldEphi(N_max, BD, k_0, -p, a_0, EE, -GG, HH);




%2 формулы для P
a_p =  a_0 * 2 *pi *I_0*E_unussualy/(IntegralFull*1i*k_0*p)*(exp(1i*k_0*p*d)-exp(-1i*k_0*p*d));
P_1 = -1/4 * IntegralFull * a_p^2;%[796.877305450313 - 606.223264834755i]
P_2 = a_0 *pi *j*a_p*E_usually*(exp(1i*k_0*p*d)-exp(-1i*k_0*p*d))/(1i*k_0*p);%[-33.8572446978129 - 2006644.89641296i]
R = 2*P_2/(I_0^2);
