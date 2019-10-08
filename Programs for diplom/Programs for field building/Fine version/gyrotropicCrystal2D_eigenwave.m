%%%%% Ёта программа позовл€ет вычисл€ть распределени€ пол€ при рассе€ние
%%%%% плоской электомагнитной волны TE типа на нескольких цилиндрах

clear all
tic

%%%% задаЄм параметры системы
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
p = 8.247632676510769; N_max = 1;%%% L = 3 * a_0;
p1 = 7.7; N_max = 0;%%% L = 3 * a_0;
p2 = 7.8; N_max = 1;%%% L = 3 * a_0;

    [BD1] = crystalCoefficients2D_eigenwave(N_max, w_0, a_0, cylXY, k_0, EE, GG, HH, c, p1);
    [BD2] = crystalCoefficients2D_eigenwave(N_max, w_0, a_0, cylXY, k_0, EE, GG, HH, c, p2);
    N_cylinders = size(cylXY,1)
%%%% задаЄм границы в которых будем строить поле
    bound_x = 500;
    bound_y = 500;
    
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
    
%%%% строим "трафарет" дл€ внутренних областей цилиндра
rho = sqrt(x.^2 + y.^2);
cylinderPlace = (sqrt((x + cylXY(1,1)).^2 + (y + cylXY(1,2)).^2) < a_0);
for ncyl = 1:N_cylinders-1
    rho1 = (sqrt((x + cylXY(ncyl + 1,1)).^2 + (y + cylXY(ncyl + 1,2)).^2));
    cylinderPlace = cylinderPlace +  (rho1 < a_0);
end

%%%% вычисл€ем пол€ внутри и вне цилидров. ѕо отдельности дл€ каждого 
%%%% цилиндра и затем суммируем
[H_z_in, H_z_out, dHz_dx_out, dHz_dx_in, dHz_dy_out, dHz_dy_in,...
 E_z_in, E_z_out, dEz_dx_out, dEz_dx_in, dEz_dy_out, dEz_dy_in] = fieldsGyrotropicCrystal2D_eigenwave(x + cylXY(1,1), y + cylXY(1,2), m, a_0, ee, k_0, p1, w_0, BD1(1:mMax));
for ncyl = 1:N_cylinders-1
    ncoeff = ncyl * mMax + 1:(ncyl + 1)* mMax;
    [H_z_in1, H_z_out1, dHz_dx_out, dHz_dx_in, dHz_dy_out, dHz_dy_in,...
     E_z_in1, E_z_out1, dEz_dx_out, dEz_dx_in, dEz_dy_out, dEz_dy_in] = fieldsGyrotropicCrystal2D_eigenwave(x + cylXY(ncyl + 1,1), y + cylXY(ncyl + 1,2), m, a_0, ee, k_0, p1, w_0, BD1(ncoeff));

    H_z_in  = H_z_in  + H_z_in1;
    H_z_out = H_z_out + H_z_out1;
    
    E_z_in  = E_z_in  + E_z_in1;
    E_z_out = E_z_out + E_z_out1;
end
[H_z_in2, H_z_out2, dHz_dx_out, dHz_dx_in, dHz_dy_out, dHz_dy_in,...
 E_z_in2, E_z_out2, dEz_dx_out, dEz_dx_in, dEz_dy_out, dEz_dy_in] = fieldsGyrotropicCrystal2D_eigenwave(x + cylXY(1,1), y + cylXY(1,2), m, a_0, ee, k_0, p2, w_0, BD2(1:mMax));
for ncyl = 1:N_cylinders-1
    ncoeff = ncyl * mMax + 1:(ncyl + 1)* mMax;
    [H_z_in3, H_z_out3, dHz_dx_out, dHz_dx_in, dHz_dy_out, dHz_dy_in,...
     E_z_in3, E_z_out3, dEz_dx_out, dEz_dx_in, dEz_dy_out, dEz_dy_in] = fieldsGyrotropicCrystal2D_eigenwave(x + cylXY(ncyl + 1,1), y + cylXY(ncyl + 1,2), m, a_0, ee, k_0, p2, w_0, BD2(ncoeff));

    H_z_in2  = H_z_in2  + H_z_in3;
    H_z_out2 = H_z_out2 + H_z_out3;
    
    E_z_in2  = E_z_in2  + E_z_in3;
    E_z_out2 = E_z_out2 + E_z_out3;
end
%%%% накладиваем "трафарет" на рассе€нное поле и суммируем с внутренним
%%%% полем цилиндров
Hz1 = H_z_in.* cylinderPlace + H_z_out.* (1 - cylinderPlace);
Ez1 = E_z_in.* cylinderPlace + E_z_out.* (1 - cylinderPlace);
Hz2 = H_z_in2.* cylinderPlace + H_z_out2.* (1 - cylinderPlace);
Ez2 = E_z_in2.* cylinderPlace + E_z_out2.* (1 - cylinderPlace);

Hz = Hz1 + Hz2;
Ez = Ez1 + Ez2;
%%%% построение реальной и мнимой частей магнитного пол€, а также его модул€ пол€
%pcolorReImModule(x / a_0, 'x / a_0', y / a_0, 'y / a_0', Hz, 'H_z')
%%%% отрисовка контуров цилиндров
%drawnCountor(a_0, step_x, cylXY, N_cylinders)

%%%% построение реальной и мнимой частей электрического пол€, а также его модул€ пол€
%pcolorReImModule(x / a_0, 'x / a_0', y / a_0, 'y / a_0', Ez1, 'E_z')
%%%% отрисовка контуров цилиндров

%%%% построение реальной и мнимой частей электрического пол€, а также его модул€ пол€
%pcolorReImModule(x / a_0, 'x / a_0', y / a_0, 'y / a_0', Ez2, 'E_z')
%%%% отрисовка контуров цилиндров

%%%% построение реальной и мнимой частей электрического пол€, а также его модул€ пол€
pcolorReImModule(x / a_0, 'x / a_0', y / a_0, 'y / a_0', Ez1, 'E_z')
%%%% отрисовка контуров цилиндров
drawnCountor(a_0, step_x, cylXY, N_cylinders)
 
% %%%% можно посмотреть анимацию реальной части пол€ в зависимости от времени
% t = [0:0.02:2] * 2 * pi / w_0;
% videoReal(t, w_0, x / a_0, 'x / a_0', y / a_0, 'y / a_0', Ez, 'E_z')
% % 
% %%%% можно посмотреть анимацию реальной части пол€ в зависимости от времени
% t = [0:0.02:6] * 2 * pi / w_0;
% videoReal(t, w_0, x / a_0, 'x / a_0', y / a_0, 'y / a_0', Hz, 'H_z')




%Bm1 = BD(([1:2*N_max + 1] - 1) * 4 + 1);
%Bm2 = BD(([1:2*N_max + 1] - 1) * 4 + 2);
%DE1 = BD(([1:2*N_max + 1] - 1) * 4 + 3);
%DH1 = BD(([1:2*N_max + 1] - 1) * 4 + 4);

%figure(11)
%hold on
%bar(m, abs(DE1))
%hold off
%figure(12)
%hold on
%bar(m, abs(DH1))
%hold off
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

