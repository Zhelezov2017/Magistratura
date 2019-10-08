function [EE, GG, HH, c] = channelparameters_sources(H0, typeOfmedia, w, d, EE_0, structParameters)

%%%%%%%% for ionosphere %%%%%%%%%%%%%%%%
% plasma chanel parameters:
if(strcmp(typeOfmedia, 'Plasma'))
n_e = 1e6;
% w_H = 1.7e7 * H0;
% w_p = 5.5e4*(n_e)^(1/2);

c = 3e10;   % velocity of light
e_0 = 4.8e-10;
m_e = 9.1e-28;

n_e = 1e13; % концентрация электронов

a_0 = 2.5;

n_e = 1e6; % концентрация электронов
H_0 = 0.5;
a_0 = 100;

w_H = (e_0 * H0) / (m_e * c);
w_p = sqrt(4 * pi * e_0^2 * n_e / m_e);

% e_0 = 4.8e-10;
% T = 2.1494;
% T_e = T * 1.6e-12;
% r_D = (T_e / (4 * pi * e_0^2 * n_e))^(1/2);
% m_e = 9.1e-28;
M_Ar=39.95*1.66e-24;
O_p = sqrt(4 * pi * e_0^2 * n_e / M_Ar);
O_H = w_H * m_e / M_Ar;
% Par = 1;
% c = 3e10;   % velocity of light

if(isfield(structParameters,'Nu_e'))
    Nu_e =  structParameters.Nu_e;
else
    Nu_e = 0.0*w_H;
end
Nu_i = 0;
w_LH = w_H * sqrt((O_p^2 + O_H^2) / (w_p^2 + w_H^2));
w_UH = sqrt(w_p^2 + w_H^2);

w_0 = w;
EE = (w_0.^2 - w_UH.^2).* (w_0.^2 - w_LH.^2).*...
            ((w_0.^2  - w_H^2).* (w_0.^2 - O_H^2)).^(-1) -...
            w_p.^2 * i*Nu_e./ ((w_H^2 - w_0.^2).* w_0);
        
%%%%% EE при учете потерь
EE = 1 + w_p.^2 * (w_0 - 1i*Nu_e)./ ((w_H.^2 - (w_0 - 1i*Nu_e).^2).* w_0) +...
         O_p.^2 * (w_0 - 1i*Nu_i)./ ((O_H.^2 - (w_0 - 1i*Nu_i).^2).* w_0);
        
GG =  -w_p.^2 * w_H./ ((w_H.^2 - (w_0 - 1i*Nu_e).^2).* w_0) +...
              O_p.^2 * O_H./ ((O_H.^2 - (w_0).^2).* w_0);
          
HH = 1 - w_p.^2./ ((w_0 - 1i*Nu_e).* w_0) - O_p.^2./ ((w_0).* w_0);


Pcin  = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));
    %prop.cnst for conic refr. waves
Pin=sqrt(EE-GG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif(strcmp(typeOfmedia, 'Graphen'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Расчёт диэлектрической проницаемсоти графена-диэлектрической
%%%%%%%%%%%% слоистой среды

%%%%% поверхностная проводимость графена. Взято из стать
% @ARTICLE{Gusynin2007,
%   author = {V P Gusynin and S G Sharapov and J P Carbotte},
%   title = {Magneto-optical conductivity in graphene},
%   journal = {Journal of Physics: Condensed Matter},
%   year = {2007},
%   volume = {19},
%   pages = {026222},
%   number = {2},
%   abstract = {Landau level quantization in graphene reflects the Dirac nature of
% 	its quasiparticles and has been found to exhibit an unusual integer
% 	quantum Hall effect. In particular, the lowest Landau level can be
% 	thought of as shared equally by electrons and holes, and this leads
% 	to characteristic behaviour of the magneto-optical conductivity as
% 	a function of frequency О© for various values of the chemical potential
% 	Ој. Particular attention is paid to the optical spectral weight under
% 	various absorption peaks and its redistribution as Ој is varied.
% 	We also provide results for magnetic field B as well as chemical
% 	potential sweeps at selected fixed frequencies, which can be particularly
% 	useful for possible measurements in graphene. Both diagonal and Hall
% 	conductivities are considered.},
%   url = {http://stacks.iop.org/0953-8984/19/i=2/a=026222}
% }
c        = 3e10;                         %%% скоротсть света в см/с
Delta    = 0;
vF       = 10^8;                         %%% в см/c
e_0      = 4.8e-10;                      %%% в СГСЭ

if(isempty(H0))
B_tesla  = 3;                            %%% внешнее магнитное поле в теслах
else
B_tesla  = H0;                           %%% внешнее магнитное поле в теслах
end

B_0      = B_tesla * 1e4;                %%% внешнее магнитное поле в эрстедах

mu_eV    = 58.67 * 1e-3;                 %%% химический потенциал в эВ
mu       = mu_eV * (1.6 * 1e-12);        %%% химический потенциал в эрг

% mu_K     = 680;                          %%% химический потенциал в кельвинах
% mu       = mu_K * 1.3805 * 1e-16;        %%% химический потенциал в эрг

T_K      = 10;                           %%% температура в кельвинах
T        = T_K * 1.3805 * 1e-16;         %%% температура в эрг
h        = 6.6748 * 1e-27;               %%% постоянная Планка в  эрг с
h_dirac  = 1.0545 * 1e-27;               %%% дираковская постоянная Планка в  эрг с

% h = 6.582119 * 1e-16;     %%% постоянная Планка  h = 6.582119 * 1e-16 эВ с

GaM_n_eV = 1.3 * 1e-3;                   %%% энергия рассения электрона в эВ
GaM_n    = GaM_n_eV * (1.6 * 1e-12);     %%% энергия рассения электрона в эрг
% GaM_n_K  = 15;                           %%% энергия рассения электрона в кельвинах
% GaM_n    = GaM_n_K  * 1.3805 * 1e-16;    %%% энергия рассения электрона в эрг
GaM_n1   = GaM_n;

% w0_eV    = 1e-2;                         %%% энергия падающего излучения в эВ
w0       = w * h_dirac;                    %%% энергия падающего излучения в эрг
% w0       = w0_eV * (1.6 * 1e-12);          %%% энергия падающего излучения в эрг

% w_H      = 1.7e7 * T_K * 1e4;            %%% циклотронная частот
% E_land   = h_dirac * w_H;                %%% первый уровень Ландау (полуклассический)

%%% функция вичисления энергии уровня Ландау
Mn     = @(n) sign(n) * sqrt(Delta.^2 + 2 *  abs(n) * h_dirac * vF.^2 * e_0 * B_0 / c);
%%% функция вичисления распределения Ферми
nF     = @(Mnn) 1./ (exp((Mnn - mu)/T) + 1);

% w_0 = w0;
% 
% w_0 = [0.1:0.1:30]*w0;
% 
% for iw = 1:size(w_0,2)
%     w0 = w_0(iw);

%     w0 = w_0;
n_max        = 25; %%% максимальный рассматриваемый уровень Ландау
sigmaDial    = 0;
sigmaNondiag = 0;
for n = 0:n_max
    
% M_n  = Mn(n) + 1e-30; %%% прибавляем эту малую величину для избавления деления на ноль
% M_n1 = Mn(n + 1);
%     
% sigmaDial = sigmaDial - e_0^2 * vF.^2 * e_0 * B_0 * h_dirac / (h * 1i * c ).*...
%       ((1-Delta^2/(1e-20+M_n.* M_n1)).* ((nF(M_n) - nF(M_n1)) + (nF(-M_n1)-nF(-M_n))).*...
%        (1./(M_n - M_n1 - w0 + 1i*(GaM_n + GaM_n1)) - 1./(M_n - M_n1 + w0 -1i*(GaM_n + GaM_n1)))./(M_n1 - M_n) +...
%        (1+Delta^2/(1e-20+M_n.* M_n1)).* ((nF(-M_n)- nF(M_n1)) + (nF(-M_n1)-nF(M_n))).*...
%        (1./(M_n + M_n1 - w0 + 1i*(GaM_n + GaM_n1)) - 1./(M_n + M_n1 + w0 -1i*(GaM_n + GaM_n1)))./(M_n1 + M_n));
%    
% sigmaNondiag = sigmaNondiag + e_0^2 * vF^2 * e_0* B_0  * h_dirac /(h * c).*...
%     (((nF(M_n)-nF(M_n1))-(nF(-M_n1)-nF(-M_n))).*...
%     ((1-Delta.^2/(M_n.* M_n1)).*(1./(M_n - M_n1 - w0 + 1i * (GaM_n + GaM_n1)) + 1./(M_n - M_n1 + w0 - 1i*(GaM_n+GaM_n1)))/(M_n1-M_n)-...
%      (1+Delta.^2/(M_n.* M_n1)).*(1./(M_n + M_n1 - w0 + 1i * (GaM_n + GaM_n1)) + 1./(M_n + M_n1 + w0 - 1i*(GaM_n+GaM_n1)))/(M_n1+M_n)));

%%% перепишем предыдущее для ускорения вычислений
M_n  =   sign(n) * sqrt(Delta.^2 + 2 * abs(n)   * h_dirac * vF.^2 * e_0 * B_0 / c) + 1e-30; %%% прибавляем эту малую величину для избавления деления на ноль
M_n1 = sign(n+1) * sqrt(Delta.^2 + 2 * abs(n+1) * h_dirac * vF.^2 * e_0 * B_0 / c);

nF_M_n    = 1./ (exp(( M_n  - mu)/T) + 1);%%nF(M_n)
nF_M_n1   = 1./ (exp(( M_n1 - mu)/T) + 1);%%nF(M_n1)
nF__M_n   = 1./ (exp((-M_n  - mu)/T) + 1);%%nF(-M_n)
nF__M_n1  = 1./ (exp((-M_n1 - mu)/T) + 1);%%nF(-M_n1)
    
sigmaDial = sigmaDial - e_0^2 * vF.^2 * e_0 * B_0 * h_dirac / (h * 1i * c ).*...
      ((1-Delta^2/(1e-20 + M_n.* M_n1)).* ((nF_M_n - nF_M_n1) + (nF__M_n1 - nF__M_n)).*...
       (1./(M_n - M_n1 - w0 + 1i*(GaM_n + GaM_n1)) - 1./(M_n - M_n1 + w0 -1i*(GaM_n + GaM_n1)))./(M_n1 - M_n) +...
       (1+Delta^2/(1e-20+M_n.* M_n1)).* ((nF__M_n- nF_M_n1) + (nF__M_n1-nF_M_n)).*...
       (1./(M_n + M_n1 - w0 + 1i*(GaM_n + GaM_n1)) - 1./(M_n + M_n1 + w0 -1i*(GaM_n + GaM_n1)))./(M_n1 + M_n));
   
sigmaNondiag = sigmaNondiag + e_0^2 * vF^2 * e_0* B_0  * h_dirac /(h * c).*...
    (((nF_M_n - nF_M_n1) - (nF__M_n1 - nF__M_n)).*...
    ((1-Delta.^2/(M_n.* M_n1)).*(1./(M_n - M_n1 - w0 + 1i * (GaM_n + GaM_n1)) + 1./(M_n - M_n1 + w0 - 1i*(GaM_n+GaM_n1)))/(M_n1-M_n)-...
     (1+Delta.^2/(M_n.* M_n1)).*(1./(M_n + M_n1 - w0 + 1i * (GaM_n + GaM_n1)) + 1./(M_n + M_n1 + w0 - 1i*(GaM_n+GaM_n1)))/(M_n1+M_n)));


end
% sigm(iw)        = sigmaDial;
% sigmaNon(iw)    = sigmaNondiag;
% end
sigm        = sigmaDial;
sigmaNon    = sigmaNondiag;


% norm_w     = sqrt(2 *  h_dirac * vF.^2 * e_0 * B_0 / c);
% % norm_w     = (1.6 * 1e-12);

% norm_sigma = e_0^2 / h;
% norm_sigma = e_0^2 / h / (GaM_n/sqrt(h_dirac * vF.^2 * e_0 * B_0 / c));
% % norm_sigma = c;
if(isempty(EE_0))
EE_0 = 4.5;              %%% диэлектрическа проницамость диэлектрика между графеновыми дисками
                         %%% В нашем случае SiO2
end

if(isempty(d))
d    = 1e-6;             %%% толщина диэлектрического слоя в см
end
w0   = w0 / h_dirac;     %%% частота падающего излучения в радиан/c
       
EE = EE_0 - 1i * 4 * pi * sigm./ (w0 * d);
GG = 4 * pi * sigmaNon./ (w0 * d);
HH = EE_0* ones(size(w0));

% EE = sigm / c;
% GG = sigmaNon / c;


end


P  = sqrt(EE - GG);
Pc  = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));
     
Pb  = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  -...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));


