% /*===========================================================================
% 
% DESCRIPTION
%       Программа вычисляет зависимость "диспесионной" функции от
%       действительной и мнимой частей постоянной распространения при
%       фиксированной частоте. В каждый момент времени строится это
%       распределение для конкретной частоты.
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

clear all;     %%% очищаем память перед новыми вычислениями
tubeparam;    %%% здесь задаются параметры плазменного столба
global EE GG HH m
EE = 10;
HH = EE;
HH = 1;
GG = 1;
m = 1;
p_0 = 1/sqrt(2);
q_0 = sqrt(1-p_0^2);
a_0 = (3.8317)./ (k_0 * q_0);
a_0 = 2.4048./ (k_0 * q_0);
R   = k * a_0;

I=0;

%%% задаём постоянные распространения мод, найденные для начального
%%% значения частоты. (найденные, например, визуально с использованием
%%% mesh_de_h)

eps = 0.0001
y0  = 0
x0  = 1
xmin =  x0-eps;
xmax =  x0+eps;
ymin =  y0-eps;
ymax =  y0+eps;
Npntx = 1000;
Npnty = 100;

px = [xmin : (xmax - xmin) / Npntx :xmax];

clear Qpts p_ p__ Po Pi

GG1=[0.6:0.001:1];
% GG = 1e-3;
for D=GG1
      
% HH = D;
GG = D;
    
m = 1;
    zz1 = (dispeq_gyrotropic_cylinder(px));
    zz2 = abs(zz1);
m = -1;
    zz1 = (dispeq_gyrotropic_cylinder(px));
    zz3 = abs(zz1); 
     
figure(1)
plot(px, log(zz2),px, log(zz3),'r')
% set(gca,'ZScale','log');
title(['g =' num2str(GG)],'FontSize',16)
% shading interp
% colorbar('location','eastoutside')
 getframe(gcf)
end;

