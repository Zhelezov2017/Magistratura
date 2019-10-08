% /*===========================================================================
% 
% DESCRIPTION
%       Программа вычисляет зависимость "диспесионной" функции от
%       действительной и мнимой частей постоянной распространения при
%       фиксированной частоте.
% 
%   Copyright (c) 2005-2013 by Vasiliy Es'kin. All Rights Reserved.
% ===========================================================================*/
% 
%                       EDIT HISTORY FOR FILE
% 
%   This section contains comments describing changes made to the module.
%   Notice that changes are listed in reverse chronological order.
% 
% when       who              what, where, why
% --------   ---       ----------------------------------------------------------
% 11/09/2007 Vasiliy Es'kin   Create programma.
% ==========================================================================*/

clear all

tic
GPC_systemParameters
% p_n___dielectricArrayWaveguides;
% p_n___gyrotropicArrayOfWaveguides
p_n___dielectricArrayWaveguides;
m = M_max;


LL = [L: -0.01 * a_0: 2* a_0];
LL = [L: 0.01 * a_0: 40 * a_0];
% LL = [L: 0.001 * a_0: 30 * a_0];

I = 0;
x0 = p_n;
x1 = [4.22   0
      7.86   0];    
% % %%%%% EE при учете потерь
[EE, GG, HH, c] = channelparameters_sources(H0, typeOfmedia, w_0, 0, 0, []);

step_p = 0.01;
for L1 = LL
        
    %%% waveguide consists of two cylinders
    cylXY = [-L1 (0.001); (0) (0.001)];
    
%     %%% waveguide consists of six cylinders
%     LS = L1;
%     cylXY =  - [LS*cos(0*2*pi/6 + pi/6)     LS*sin(0*2*pi/6 + pi/6)+0.001*a_0;
%             LS*cos(1*2*pi/6 + pi/6)     LS*sin(1*2*pi/6 + pi/6)+0.001*a_0;
%             LS*cos(2*2*pi/6 + pi/6)     LS*sin(2*2*pi/6 + pi/6)+0.001*a_0;
%             LS*cos(3*2*pi/6 + pi/6)     LS*sin(3*2*pi/6 + pi/6)+0.001*a_0;
%             LS*cos(4*2*pi/6 + pi/6)     LS*sin(4*2*pi/6 + pi/6)+0.001*a_0;
%             LS*cos(5*2*pi/6 + pi/6)     LS*sin(5*2*pi/6 + pi/6)+0.001*a_0;
%             ];
    
    
    I=I+1;
    step_p = abs(x0(1,1) - x0(2,1))/300;
    if(step_p < 1e-4) step_p = 1e-4; end
    
    %%% ищем моды с помощью матлабовской функции fminserch
    for J=1:size(x0,1)
        
%         if(x0(J,1) >= 1.02)
%             x0(J,:) = fminsearch(@(x) (dispeq_gyrotropic_cylindersArrayFull(m, a_0, cylXY, k_0,  x(1) + 1i * x(2), EE, GG, HH)), [x0(J,:)], optimset('TolX',1e-15, 'TolFun',1e-15));
            x0(J,1) = fminbnd(@(x) log10(dispeq_gyrotropic_cylindersArrayFull(M_max, a_0, cylXY, k_0, x, EE, GG, HH)./...
                 ((x - Pcin).* x).^((2 * M_max + 1) * size(cylXY,1)/2)), x0(J,1)-step_p, x0(J,1)+step_p);
            
%             x = [x0(J,1)-step_p: step_p/20: x0(J,1)+step_p];
%             for ix=1:size(x,2)
%                 ff(ix) = log10(dispeq_gyrotropic_cylindersArrayFull(m, a_0, cylXY, k_0, x(ix), EE, GG, HH));
%             end
%             imin = min(ff);
%             x0(J,1) = x(ff==imin);
% %             figure(1)
% %             hold on
% %              plot(x, ff)
% %              hold off
            
%         else
%             x0(J,:) = x0(J,:) * 0;
%         end 

         x1(J,:) = fminbnd(@(x) (dispeq_gyrotropic_cylindersArrayTwoCyls_Approx(a_0, L1, k_0, x, EE, GG, HH)), (1-step_p) * x1(J,1), (1+step_p) * x1(J,1));
  end
    
    p_nRe(:,I)=x0(:,1);
    p_nIm(:,I)=x0(:,2);
    
        p_nReApp(:,I)=x1(:,1);
        p_nImApp(:,I)=x1(:,2);
end

figure(23)
hold on
plot(LL / a_0, real(p_nRe))
plot(LL / a_0, real(p_nReApp), 'b')
hold off

% figure(24)
% hold on
% plot(LL / a_0, imag(p_nIm))
% hold off

%%
%%%%%% вычисляем элементы матрицы рассеяния для внешней области цилиндра
% euler = double(eulergamma)
p = [4:0.01:8];
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k_0.* a_0 * q;
euler = 0.577215664901533;
    H2m  =  besselh(0, 2, Q);
    H2M =  besselh(1, 2, Q);
    Jm = besselj(0, Q);
    JM = besselj(1, Q);
    
  L1 = 2./(k_0 * q).* exp(- euler - pi/2 * 1i * (1 - (H2M - HH / 2 * Q.* H2m)./ (JM - HH / 2 * Q.* Jm)));
  L2 = 2./(k_0 * q).* exp(- euler - pi/2 * 1i * (1 + (H2M - HH / 2 * Q.* H2m)./ (JM - HH / 2 * Q.* Jm)));
  
figure(23)
hold on
plot(real(L1)./ a_0, (p), 'r* ')
plot(real(L2)./ a_0, (p), 'b* ')
hold off

 L1appr = 2./(k_0 * q).* exp(- euler - pi/2 * 1i * (1 - (1i * 2./Q/pi - HH / 2 * Q.* (1 - 1i * 2 / pi * (log(Q / 2) + euler)))./...
     (Q/2 - HH / 2 * Q)));
 L2appr = 2./(k_0 * q).* exp(- euler - pi/2 * 1i * (1 + (1i * 2./Q/pi - HH / 2 * Q.* (1 - 1i * 2 / pi * (log(Q / 2) + euler)))./...
     (Q/2 - HH / 2 * Q)));
  
figure(23)
hold on
plot(real(L1appr)./ a_0, (p), 'r- ')
plot(real(L2appr)./ a_0, (p), 'b- ')
hold off

toc












