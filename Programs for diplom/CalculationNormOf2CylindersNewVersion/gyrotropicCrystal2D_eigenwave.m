%%%%% Эта программа позовляет вычислять норму
%%%%% п

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters

n = 1000;
a = -1000;  
b = 1000;
cr = -1000;
d = 1000;
h = (b-a) / n; 
x=a:h:b;
y=cr:h:d;
L = 3 * a_0;
Outside = 0;
N_max = 0;
%[-1.48525898195422e+25 - 1.04231343344304e+24i]
N_cylinders = 2;
p = 8.298;
wR = 1.1*w_H:1.1*w_H:4* w_H;
M = zeros(1,length(wR)); 
i=1;
wm = 1.1 * w_H;
cylXY = [(-L) 0;(0) (0.001)];

BD = crystalCoefficients2D_eigenwave(N_max, wm, a_0, cylXY, k_0, EE, GG, HH, c, p);
RateCalculationMatrix = RateCalculationMatrix(N_max, k_0, p, a_0, EE, GG, HH , BD, N_cylinders );
RateCalculationMatrixOutInt = RateCalculationMatrixOutInt(N_max, k_0, p, a_0, EE, GG, HH , BD, N_cylinders );

Integral = Integral(RateCalculationMatrix);
BeenFieldEX =  BeenFieldEX(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
BeenFieldEY =  BeenFieldEY(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
BeenFieldHX =  BeenFieldHX(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
BeenFieldHY =  BeenFieldHY(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
VectorProduct = VectorProduct( BeenFieldHX, BeenFieldHY, BeenFieldEX, BeenFieldEY );
Sum = 0;
Man = sum(BeenFieldEX)+sum(BeenFieldEY)+sum(BeenFieldHX)+sum(BeenFieldHY);
Man(1,1001) = 0;
%IntegralOut = IntegralOutside(Sum , x, y);
IntOut = IntegralOutside(Man, x, y);
IntegralOutIntVent = IntegralNew(RateCalculationMatrixOutInt);
IntegralFull = Integral + IntOut - IntegralOutIntVent;


%%%%%Рассчет R(омега)
%E_usually = ElectricField(N_max, BD, k_0, p, a_0, EE, GG, HH);
%d = 2;
%I_0 = 4;
%j = I_0 / d;
%E_unussualy = ElectricField(N_max, BD, k_0, -p, a_0, EE, -GG, HH);
%a_p =  a_0 * 2 *pi *I_0*E_unussualy/(IntegralFull*1i*k_0*p)*(exp(1i*k_0*p*d)-exp(-1i*k_0*p*d));
%P_1 = 1/4 * IntegralFull * a_p^2;
%R = 2*P_1/(I_0^2);
%M(i)=R;
%i=i+1;

  

    




   