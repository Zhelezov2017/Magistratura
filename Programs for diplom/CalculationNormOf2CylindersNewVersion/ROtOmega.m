%%%%% Эта программа позовляет вычислять норму
%%%%% п

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders p_nRe
GPC_systemParameters
findDispersionCurvesOfModes_ArrayCyls

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
N_cylinders = 2;
p = 8.298;
wR = 1.1*w_H:0.1*w_H:4* w_H;
M = zeros(1,length(wR)); 
i=1;
cylXY = [(-L) 0;(0) (0.001)];

for wm=1.1*w_H:0.1*w_H:4*w_H
    k_0 = wm / c;
  
    
    % %%%%% EE при учете потерь
    

    
    
    
    BD = crystalCoefficients2D_eigenwave(N_max, wm, a_0, cylXY, k_0, EE, GG, HH, c, p_nRe(2,i));
    rateCalculationMatrix = RateCalculationMatrix(N_max, k_0, p_nRe(2,i), a_0, EE, GG, HH , BD, N_cylinders );
    rateCalculationMatrixOutInt = RateCalculationMatrixOutInt(N_max, k_0, p_nRe(2,i), a_0, EE, GG, HH , BD, N_cylinders );
    integral = Integral(rateCalculationMatrix);
    beenFieldEX =  BeenFieldEX(N_cylinders, N_max, cylXY, BD,k_0, p_nRe(2,i), x, y);
    beenFieldEY =  BeenFieldEY(N_cylinders, N_max, cylXY, BD,k_0, p_nRe(2,i), x, y);
    beenFieldHX =  BeenFieldHX(N_cylinders, N_max, cylXY, BD,k_0, p_nRe(2,i), x, y);
    beenFieldHY =  BeenFieldHY(N_cylinders, N_max, cylXY, BD,k_0, p_nRe(2,i), x, y);
    vectorProduct = VectorProduct( beenFieldHX, beenFieldHY, beenFieldEX, beenFieldEY );
    Sum = 0;
    Man = sum(beenFieldEX)+sum(beenFieldEY)+sum(beenFieldHX)+sum(beenFieldHY);
    Man(1,501) = 0;
    %IntegralOut = IntegralOutside(Sum , x, y);
    IntOut = IntegralOutside(Man, x, y);
    IntegralOutIntVent = IntegralNew(rateCalculationMatrixOutInt);
    IntegralFull = integral + IntOut - IntegralOutIntVent;


    %%%%%Рассчет R(омега)
    E_usually = ElectricField(N_max, BD, k_0, p_nRe(2,i), a_0, EE, GG, HH);
    d = 2;
    I_0 = 4;
    j = I_0 / d;
    E_unussualy = ElectricField(N_max, BD, k_0, -p_nRe(2,i), a_0, EE, -GG, HH);
    a_p =  a_0 * 2 *pi *I_0*E_unussualy/(IntegralFull*1i*k_0*p_nRe(2,i))*(exp(1i*k_0*p_nRe(2,i)*d)-exp(-1i*k_0*p_nRe(2,i)*d));
    P_1 = -1/4 * IntegralFull * a_p^2;
    R = 2*P_1/(I_0^2);
    M(i)=R*898755199;
    i=i+1;
end
plot(wR/w_H,real(M));