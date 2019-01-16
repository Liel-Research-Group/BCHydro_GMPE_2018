function [ Sa, sig ] = BCHydro2018( T, M, R, Ztor, Vs30, F )

% Coded by Robert Chase - University of Colorado Boulder
% Last Modified 1/16/2019

% Inputs
% T = Period [s]
% M = Magnitude [Mw]
% R = Rupture Distance [km]
% Ztor = Depth to Slab [km]
% Vs30 = Average shear wave velocity in top 30 m [m/s]
% F is a inslab flag; F = 1 => Inslab, F = 0 => Interface

% Outputs
% Sa = Spectral Acceleration in g
% sig = Sigma in ln space

Tvec = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];

tol = 0.00001;
PGA1000 = 0.00001;
dif = 100;
count = 0;
while dif > tol
 PGA1000_guess = PGA1000;
 PGA1000 = BCHydro2018_PGA1000( M, R, Ztor, F, PGA1000_guess );
 dif = abs(PGA1000 - PGA1000_guess);
 count = count + 1;
 if count == 10000
     break
 end
end

% Non period dependent coefficients
n = 1.18; c = 1.88; C4 = 10; a5 = 0; a9 = 0.4; a10 = 1.73; C1_slab = 7.2;  
%a3 = -0.1;%report
a3 = 0.1; %spreadsheet - correct input

% Period Dependent Coefficients
a1_vec = [2.340 2.360 2.384 2.446 2.751 3.019 3.349 3.284 3.211 3.145 2.997 2.839 2.658 2.346 1.851 1.216 0.649 0.082 -0.369 -1.034 -1.520 -1.810 -2.173 -2.712];
a2_vec = [-1.044 -1.044 -1.08 -1.11 -1.11 -1.11 -1.084 -1.027 -0.983 -0.947 -0.89 -0.845 -0.809 -0.76 -0.698 -0.612 -0.55 -0.501 -0.46 -0.455 -0.45 -0.45 -0.45 -0.45];
a4_vec = [0.59 0.59 0.59 0.59 0.59 0.59 0.59 0.62 0.64 0.66 0.68 0.68 0.68 0.68 0.68 0.68 0.68 0.68 0.68 0.68 0.73 0.78 0.84 0.93];
a6_vec = [-0.00705 -0.00707 -0.00710 -0.00725 -0.00758 -0.00788 -0.00820 -0.00835 -0.00835 -0.00828 -0.00797 -0.00770 -0.00740 -0.00698 -0.00645 -0.00570 -0.00510 -0.00465 -0.00430 -0.00390 -0.00370 -0.00357 -0.00340 -0.00327];
a11_vec = [0.0170 0.0170 0.0170 0.0180 0.0180 0.0180 0.0175 0.0170 0.0160 0.0152 0.0140 0.0130 0.0122 0.0113 0.0100 0.0082 0.0070 0.0060 0.0052 0.0040 0.0030 0.0022 0.0013 0.0000];
a12_vec = [0.818 0.857 0.921 1.007 1.225 1.457 1.849 2.082 2.240 2.341 2.415 2.359 2.227 1.949 1.402 0.329 -0.487 -0.770 -0.700 -0.607 -0.540 -0.479 -0.393 -0.350];
a13_vec = [-0.0135 -0.0135 -0.0135 -0.0138 -0.0142 -0.0145 -0.0153 -0.0162 -0.0172 -0.0183 -0.0206 -0.0231 -0.0256 -0.0296 -0.0363 -0.0493 -0.061 -0.0711 -0.0798 -0.0935 -0.098 -0.098 -0.098 -0.098];
a14_vec = [-0.223 -0.196 -0.128 -0.130 -0.130 -0.130 -0.156 -0.172 -0.184 -0.194 -0.210 -0.223 -0.233 -0.245 -0.261 -0.285 -0.301 -0.313 -0.323 -0.282 -0.250 -0.250 -0.250 -0.250];

Vlin_vec = [865.1 865.1 907.8 1053.5 1085.7 1032.5 877.6 748.2 654.3 587.1 503 456.6 430.3 410.5 400 400 400 400 400 400 400 400 400 400];
b_vec = [-1.186 -1.219 -1.273 -1.346 -1.471 -1.624 -1.931 -2.188 -2.381 -2.518 -2.657 -2.669 -2.599 -2.401 -1.955 -1.025 -0.299 0 0 0 0 0 0 0];
C1_inter_vec = [8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.2 8.15 8.1 8.05 8 7.95 7.9 7.85 7.8 7.8 7.8 7.8]; 
Adj_Inter_vec = [1.04398 1.05 1.23 1.34 1.32 1.32 1.21 1.14 1.05 0.95 0.79 0.66 0.54 0.36 0.24 -0.08 -0.21 -0.21 -0.22 -0.06 0.06 0.09 0.14 0.31];
Adj_Inslab_vec = [0.83429 0.79 0.71 0.98 0.99 1.00 0.92 0.88 0.81 0.75 0.62 0.54 0.46 0.33 0.23 -0.01 -0.09 -0.05 -0.02 0.06 0.06 0.10 0.13 0.24];

phi_vec = [0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62 0.62];
tau_vec = [0.58 0.58 0.58 0.58 0.58 0.58 0.56 0.54 0.52 0.505 0.48 0.46 0.45 0.45 0.45 0.45 0.45 0.45 0.45 0.45 0.45 0.45 0.45 0.45];

% Interpolation depending on period
a1 = interp1(Tvec,a1_vec,T,'linear','extrap');
a2 = interp1(Tvec,a2_vec,T,'linear','extrap');
a4 = interp1(Tvec,a4_vec,T,'linear','extrap');
a6 = interp1(Tvec,a6_vec,T,'linear','extrap');
a11 = interp1(Tvec,a11_vec,T,'linear','extrap');
a12 = interp1(Tvec,a12_vec,T,'linear','extrap');
a13 = interp1(Tvec,a13_vec,T,'linear','extrap');
a14 = interp1(Tvec,a14_vec,T,'linear','extrap');

Vlin = interp1(Tvec,Vlin_vec,T,'linear','extrap');
b = interp1(Tvec,b_vec,T,'linear','extrap');
C1_inter = interp1(Tvec,C1_inter_vec,T,'linear','extrap');
Adj_Inter = interp1(Tvec,Adj_Inter_vec,T,'linear','extrap');
Adj_Inslab = interp1(Tvec,Adj_Inslab_vec,T,'linear','extrap');

phi = interp1(Tvec,phi_vec,T,'linear','extrap');
tau = interp1(Tvec,tau_vec,T,'linear','extrap');

% fmag calculation
C1 = F*C1_slab+(1-F)*C1_inter;
SourceAdj = F*Adj_Inslab+(1-F)*Adj_Inter;

if M <= C1
    fmag = a4*(M-C1)+a13*(10-M)^2;
elseif M > C1
    fmag = a5*(M-C1)+a13*(10-M)^2;
end

% fsite calculation
if Vs30 > 1000
    Vs = 1000;
elseif Vs30 <= 1000
    Vs = Vs30;
end

if Vs30 < Vlin
    fsite = a12*log(Vs/Vlin)-b*log(PGA1000+c)+b*log(PGA1000+c*(Vs/Vlin)^n); %spreadsheet this is right
    %fsite = a12*log(Vs/Vlin)-b*log(PGA1000+c)*b*log(PGA1000+c*(Vs/Vlin)^n); %report
elseif Vs30 >= Vlin
    fsite = (a12+b*n)*log(Vs/Vlin);
end

% fztor calculation
if Ztor <= 100
    fztor = a11*(Ztor-60)*F;
elseif Ztor > 100
    fztor = a11*(100-60)*F;
end

lnPSA = a1+a4*(C1_slab-C1_inter)*F+(a2+a14*F+a3*(M-7.8))*log(R+C4*exp((M-6)*a9))+a6*R+a10*F+fmag+fztor+fsite;

%Outputs
Sa = exp(lnPSA)*exp(SourceAdj); % Spectral Acceleration in g
sig = sqrt(phi^2+tau^2); % Sigma in ln space

end

