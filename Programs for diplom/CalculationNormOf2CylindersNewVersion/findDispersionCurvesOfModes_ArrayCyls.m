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
p_n___dielectricArrayWaveguides;
m = M_max;

% k = [1:-0.00001:0.9] * k_0;
% k = [1:-0.001:0.1] * k_0;
I = 0;
x0 = p_n;
x1 = [4.22   0
      7.86   0];
% for D_Q = k

w0 = [3.5: -0.01 : 1.1] * w_H;
w0 = [1.1: 0.1 : 4] * w_H;
for D_Q = w0


    w_0 = D_Q;
    k_0 = w_0 / c;
    R = a_0 * k_0;
    
    % %%%%% EE при учете потерь
    [EE1, GG1, HH1, c] = channelparameters_sources(H0, typeOfmedia, w_0, 0, 1, []);
    EE = EE1;
    GG = GG1;
    HH = HH1;
    Pcin  = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));

    I=I+1;
%     step_p = 0.01;
    step_p = abs(x0(1,1) - x0(2,1))/30;
    if(step_p < 1e-4) step_p = 1e-4; end
    
%     k_0 = D_Q;
%     deltax = 0.1;
    
    Qpts(I) = w_0 / w_H;
    
%     Pin(I)  = (sqrt(EE-GG));
%     Pcin(I) = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
%          2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));

    %%% ищем моды с помощью матлабовской функции fminserch
    for J=1:size(x0,1)

        %%%% single cylinder
%         x0(J,:) = fminsearch(@(x) dispeq_gyrotropic_cylinder(0, x(1) + 1i * x(2), EE, GG, HH, k_0, a_0), [x0(J,:)], optimset('TolX',1e-8));
        
%         %%%% multiple cylinder waveguide
%         x0(J,:) = fminsearch(@(x) dispeq_gyrotropic_cylindersArrayFull(M_max, a_0, cylXY, k_0,  x(1) + 1i * x(2), EE, GG, HH), [x0(J,:)], optimset('TolX',1e-8));
        
        if(x0(J,1) >= 1.02)
            x0(J,:) = fminbnd(@(x) log10(dispeq_gyrotropic_cylindersArrayFull(M_max, a_0, cylXY, k_0, x, EE, GG, HH)./...
                                  ((x - Pcin).* x).^((2 * M_max + 1) * size(cylXY,1)/2)),...
                                  (1-step_p) * x0(J,1), (1+step_p) * x0(J,1));
        else
            x0(J,:) = x0(J,:) * 0;
        end
        p_nRe(:,I)=x0(:,1);
        p_nIm(:,I)=x0(:,2);
        
        x1(J,:) = fminbnd(@(x) (dispeq_gyrotropic_cylindersArrayTwoCyls_Approx(a_0, L, k_0, x, EE, GG, HH)), (1-step_p) * x1(J,1), (1+step_p) * x1(J,1));
        p_nReApp(:,I)=x1(:,1);
        p_nImApp(:,I)=x1(:,2);
    end

end
% 
% figure(23)
% hold on
% plot(k, real(p_nRe))
% hold off
% 
% figure(24)
% hold on
% plot(k, imag(p_nIm))
% hold off


figure(23)
hold on
plot(w0/w_H, real(p_nRe))
%plot(w0/w_H, real(p_nReApp))
hold off

% figure(24)
% hold on
% plot(w0/w_H, imag(p_nIm))
% plot(w0/w_H, imag(p_nImApp))
% hold off

toc



