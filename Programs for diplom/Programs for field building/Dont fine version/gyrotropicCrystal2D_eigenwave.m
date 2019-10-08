%%%%% Эта программа позовляет вычислять распределения поля при рассеяние
%%%%% плоской электомагнитной волны TE типа на нескольких цилиндрах

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p
GPC_systemParameters

p = 4.9;
p = 3.2;
p = 3.08;
p = 2.786684413716663;
% p = 1.881115534409880;
p = 3.156;
N_max = 4;


% p = 7.866648559789426; N_max = 0;%%% L = 3 * a_0;
p = 8.298; N_max = 3;%%% L = 3 * a_0;
%p = 12.2; N_max = 1;%%% L = 3 * a_0;

%%%% - здесь не верно - start
%  p = 28.941988095093570; N_max = 1;%%% L = 3 * a_0;
%  p = 32.190805143748310; N_max = 3;%%% L = 3 * a_0;
%  p = 33.207896533060810; N_max = 4;%%% L = 3 * a_0;
%  p = 33.461563322108020; N_max = 5;%%% L = 3 * a_0;
%  p = 33.516639015053980; N_max = 6;%%% L = 3 * a_0;
%  p = 33.527680512473594; N_max = 7;%%% L = 3 * a_0;
%  p = 33.529788115239610; N_max = 8;%%% L = 3 * a_0;
%  p = 33.530177110875470; N_max = 9;%%% L = 3 * a_0; 
%  p = 33.530247127319570; N_max = 10;%%% L = 3 * a_0;
%  p = 33.530261608010235; N_max = 12;%%% L = 3 * a_0;
%%%% - здесь не верно - end

%p = 5.871521662116040; N_max = 0; %%% L = 10 * a_0;
%p = 5.873890418824283; N_max = 1; %%% L = 10 * a_0;
% p = 5.286542845723133; N_max = 1; %%% L = 10 * a_0;
%p = 5.873926183388964; N_max = 3;  %%% L = 10 * a_0;
% % p = 5.286578399794019; N_max = 3;  %%% L = 10 * a_0;
% % ошибочно p = 5.935074304863453; N_max = 1;%%% L = 10 * a_0;
% % ошибочно p = 5.283865266442298; N_max = 1;%%% L = 10 * a_0;
% % ошибочно p = 5.937033970898658; N_max = 3;%%% L = 10 * a_0;
% % p = 5.937034355153082; N_max = 5;%%% L = 10 * a_0;
% % p = 5.937034355202124; N_max = 9;%%% L = 10 * a_0;
 
 
%  p = 8.065893043214102; N_max = 5; %%% L = 5 * a_0;
 
%  p = 5.588458097319342; N_max = 1; %%% L = 30 * a_0;
%  p = 5.588458297910931; N_max = 5; %%% L = 30 * a_0; 
 
% p = 5.596860278775192; N_max = 1; %%% L = 25 * a_0;

% p = 5.583360260119487; N_max = 5; %%% L = 300 * a_0;
% p = 5.583360260119487; N_max = 5; %%% L = 300 * a_0;

%  p = 1.207911561297986; N_max = 3; %%% L = 30 * a_0; a_0 = 1000;

% p = 30.3006; N_max = 5;%%% L = 3 * a_0; w_H = - 0.00001;
% p = 30.330261164313598; N_max = 12;%%% L = 3 * a_0; w_H = - 0.00001;
 
 

%%%% вычисляем коэффициены рассеяния (D), коэффициенты падения (C) и
%%%% коэффициенты рассеяния одиночного цилиндра (S)
    [BD] = crystalCoefficients2D_eigenwave(N_max, w_0, a_0, cylXY, k_0, EE, GG, HH, c, p);

    N_cylinders = size(cylXY,1)
%%%% задаём границы в которых будем строить поле
    bound_x = 200;
    bound_y = 200;
    
    L1 = 3 * a_0;
    bound_x = 2 * L1;
    bound_y = 2 * L1;
    
%     bound_x = a_0;
%     bound_y = a_0;
    
    step_x = 0.01 * bound_x;
    step_y = 0.01 * bound_y;
    [x y] = meshgrid([- bound_x :step_x: bound_x],[- bound_y :step_y: bound_y]);
    mMax = 4 * (2 * N_max + 1);
    m  = [-N_max:N_max];
    
%%%% строим "трафарет" для внутренних областей цилиндра
rho = sqrt(x.^2 + y.^2);
cylinderPlace = (sqrt((x + cylXY(1,1)).^2 + (y + cylXY(1,2)).^2) < a_0);
for ncyl = 1:N_cylinders-1
    rho1 = (sqrt((x + cylXY(ncyl + 1,1)).^2 + (y + cylXY(ncyl + 1,2)).^2));
    cylinderPlace = cylinderPlace +  (rho1 < a_0);
end

%%%% вычисляем поля внутри и вне цилидров. По отдельности для каждого 
%%%% цилиндра и затем суммируем
[H_z_in, H_z_out, dHz_dx_out, dHz_dx_in, dHz_dy_out, dHz_dy_in,...
 E_z_in, E_z_out, dEz_dx_out, dEz_dx_in, dEz_dy_out, dEz_dy_in] = fieldsGyrotropicCrystal2D_eigenwave(x + cylXY(1,1), y + cylXY(1,2), m, a_0, ee, k_0, p, w_0, BD(1:mMax));
for ncyl = 1:N_cylinders-1
    ncoeff = ncyl * mMax + 1:(ncyl + 1)* mMax;
    [H_z_in1, H_z_out1, dHz_dx_out, dHz_dx_in, dHz_dy_out, dHz_dy_in,...
     E_z_in1, E_z_out1, dEz_dx_out, dEz_dx_in, dEz_dy_out, dEz_dy_in] = fieldsGyrotropicCrystal2D_eigenwave(x + cylXY(ncyl + 1,1), y + cylXY(ncyl + 1,2), m, a_0, ee, k_0, p, w_0, BD(ncoeff));

    H_z_in  = H_z_in  + H_z_in1;
    H_z_out = H_z_out + H_z_out1;
    
    E_z_in  = E_z_in  + E_z_in1;
    E_z_out = E_z_out + E_z_out1;
end

%%%% накладиваем "трафарет" на рассеянное поле и суммируем с внутренним
%%%% полем цилиндров
Hz = H_z_in.* cylinderPlace + H_z_out.* (1 - cylinderPlace);
Ez = E_z_in.* cylinderPlace + E_z_out.* (1 - cylinderPlace);

%%%% построение реальной и мнимой частей магнитного поля, а также его модуля поля
pcolorReImModule(x / a_0, 'x / a_0', y / a_0, 'y / a_0', Hz, 'H_z')
%%%% отрисовка контуров цилиндров
drawnCountor(a_0, step_x, cylXY, N_cylinders)

%%%% построение реальной и мнимой частей электрического поля, а также его модуля поля
pcolorReImModule(x / a_0, 'x / a_0', y / a_0, 'y / a_0', Ez, 'E_z')
%%%% отрисовка контуров цилиндров
drawnCountor(a_0, step_x, cylXY, N_cylinders)


 
% %%%% можно посмотреть анимацию реальной части поля в зависимости от времени
% t = [0:0.02:2] * 2 * pi / w_0;
% videoReal(t, w_0, x / a_0, 'x / a_0', y / a_0, 'y / a_0', Ez, 'E_z')
% % 
% %%%% можно посмотреть анимацию реальной части поля в зависимости от времени
% t = [0:0.02:6] * 2 * pi / w_0;
% videoReal(t, w_0, x / a_0, 'x / a_0', y / a_0, 'y / a_0', Hz, 'H_z')

Bm1 = BD(([1:2*N_max + 1] - 1) * 4 + 1);
Bm2 = BD(([1:2*N_max + 1] - 1) * 4 + 2);
DE1 = BD(([1:2*N_max + 1] - 1) * 4 + 3);
DH1 = BD(([1:2*N_max + 1] - 1) * 4 + 4);

figure(11)
hold on
bar(m, abs(DE1))
hold off
figure(12)
hold on
bar(m, abs(DH1))
hold off
% 
% 
% 
% %%%% energies are carried by the separated waves
% q = sqrt(1 - p.^2);
% size_DE1 = size(DE1,1);
% P_n0 = zeros(size_DE1,1);
% for im = 1:size_DE1
%     P_n0(im) = wavesPower_of_discreteModes(q, p, 1, k_0, a_0, EE, GG, HH, m(im), c, Bm1(im), Bm2(im), DE1(im), DH1(im));
% end
% 
% figure(111)
% hold on
% bar(m, abs(P_n0))
% hold off







toc
