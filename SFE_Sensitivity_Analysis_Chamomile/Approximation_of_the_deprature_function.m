clc, close all
clear all

addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters_cell         = table2cell(Parameters_table(:,3));

%% 1) Define function H = f(T,P).
%f = @(T,P) 1000 + T.*P - 0.02*(T.^2);  % Example

T           = MX.sym('T');
P           = MX.sym('P');

Z           = Compressibility( T, P,         Parameters_cell );
rho         = rhoPB_Comp(      T, P, Z,      Parameters_cell );
h           = SpecificEnthalpy(T, P, Z, rho, Parameters_cell );

g = Function('g',{[T, P]},{h});

%% 2) Generate a grid of (T,P) and compute H.
Tmin = 30+273;  Tmax = 50+273;
Pmin = 80;      Pmax = 300;
nT   = 100;     nP   = 200;

Tvec = linspace(Tmin, Tmax, nT);
Pvec = linspace(Pmin, Pmax, nP);

[TT, PP] = meshgrid(Tvec, Pvec);

HH       = nan(nP, nT);

TT       = TT(:);
PP       = PP(:);
HH       = HH(:);

for ii=1:length(HH)
    HH(ii) = full(g([TT(ii), PP(ii)]));
end

TT = reshape(TT(:),nP,nT);
PP = reshape(PP(:),nP,nT);
HH = reshape(HH(:),nP,nT);

%% 3) Fit a polynomial surface for T in terms of (H,P).
%   Flatten the data, then use 'fit' with a chosen polynomial type, e.g. 'poly22'.
[surfaceFit, gof] = fit([HH(:), PP(:)], TT(:), 'poly33');
gof

%% 4) Create side-by-side plots: 
figure('Name','Two Surfaces','Units','normalized','Position',[0.1 0.1 0.8 0.4]);

% --- (a) Left plot: original surface H=f(T,P)
subplot(1,3,1);
surf(HH, PP, TT-273, 'EdgeColor','none');
xlabel('H'); ylabel('P'); zlabel('T');
title('Surface 1: H = f(T,P)');
view(2);
colormap jet; colorbar
axis tight; grid off

% --- (b) Right plot: fitted surface T ~ g(H,P)
subplot(1,3,2);

% Evaluate the fitted polynomial
TTplot = surfaceFit(HH, PP);

surf(HH, PP, TTplot-273, 'EdgeColor','none');
view(2)
xlabel('H'); ylabel('P'); zlabel('T');
title('Surface 2: T ~ g(H,P) from polynomial fit');
colormap jet; colorbar
axis tight; grid off

% Show difference between dataset and predictions
subplot(1,3,3)
surf(HH, PP, TTplot - TT, 'EdgeColor','none');
view(2)
xlabel('H'); ylabel('P'); zlabel('T');
title('Surface 3: Difference f and g functions');
colormap jet; colorbar
axis tight; grid off

%% 5) (Optional) Check the fit or display results
disp('Fitted polynomial model (T as function of H,P):');
disp(surfaceFit);

%% introduce Casadi variables
%{
P = 80;

h=[];
for T=[linspace(40,80)]+273
    Z           = Compressibility( T, P,         Parameters_cell );
    rho         = rhoPB_Comp(      T, P, Z,      Parameters_cell );
    h           =[h, SpecificEnthalpy(T, P, Z, rho, Parameters_cell )];
end

T_s             = MX.sym('T_s',numel(h),1);

Z               = Compressibility( T_s, P,         Parameters_cell );
rho             = rhoPB_Comp(      T_s, P, Z,      Parameters_cell );
h_sym           = SpecificEnthalpy(T_s, P, Z, rho, Parameters_cell );

H               = h' - h_sym;

g = Function('g',{T_s},{H});
G = rootfinder('G','newton',g);

tic
( G(40+273)-273 - linspace(40,80)' )
toc

T_full(G(40+273)-273)
%}