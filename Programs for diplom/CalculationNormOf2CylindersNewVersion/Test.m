clear all
tic

%%%% ����� ��������� �������
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters

p=6;
N_max =0;
m = [-N_max: N_max];  
if(N_max == 0) m = 0; end
wm = 1.1 * w_H;
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k_0.* a_0 * q;
BD = crystalCoefficients2D_eigenwave(N_max, wm, a_0, cylXY, k_0, EE, GG, HH, c, p);    
    %������� ������ ��������
    valueSet = zeros(4,(2*N_max+1)*N_cylinders);
    for jj = 1: N_cylinders
        for jm = 1:2*N_max+1
            DC = zeros(1,4);
            for jh = 1:4
                for jn = 1:4
                    DC(1,jn)=BD(4*(2*N_max+1)*jj-4*jm-jn+5);
                end
                valueSet(jh,(2*N_max+1)*jj-jm+1)=DC(1,5-jh);
            end 
            
        end
    end  
   
    
    IntegralForm1 = 0;
    IntegralForm2 = 0;
    for jj = 1:2*N_max+1
        %R= values(M,{2});
        R(1:4)=valueSet(1:4,jj);
        M(1:4)=valueSet(1:4,jj+2*N_max+1);
        H2m  = besselh(m(jj), 2, Q);
        dH2m = (H2m.* m(jj))./ Q  - besselh(m(jj) + 1, 2, Q);
        FieldOutsideEPHI1 =-1i*R(3)*dH2m-1i*R(4)*p*m(jj)*H2m/Q;
        FieldOutsideERHO1 = p*R(4)*dH2m+R(3)*m(jj)*H2m/Q;
        FieldOutsideHPHI1 =-R(4)*dH2m+R(3)*p*m(jj)*H2m/Q;
        FieldOutsideHRHO1 = 1i*p*R(3)*dH2m+1i*R(4)*m(jj)*H2m/Q;
        
        FieldOutsideEPHI2 =-1i*M(3)*dH2m-1i*M(4)*p*m(jj)*H2m/Q;
        FieldOutsideERHO2 = p*M(4)*dH2m+M(3)*m(jj)*H2m/Q;
        FieldOutsideHPHI2 =-M(4)*dH2m+M(3)*p*m(jj)*H2m/Q;
        FieldOutsideHRHO2 = 1i*p*M(3)*dH2m+1i*M(4)*m(jj)*H2m/Q;
        
        MagnetOnElect1 = 2*pi*(FieldOutsideERHO1*FieldOutsideHPHI1-FieldOutsideEPHI1*FieldOutsideHRHO1);
        MagnetOnElect2 = 2*pi*(FieldOutsideERHO2*FieldOutsideHPHI2-FieldOutsideEPHI2*FieldOutsideHRHO2);
        InternalForm1 = IntegralForm1 + MagnetOnElect1;     
        InternalForm2 = IntegralForm2 + MagnetOnElect2;
        RateCalculationMatrixOutInt =+ InternalForm1 + InternalForm2;
    end