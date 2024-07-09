startup;
delete(gcp('nocreate'));

addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%%
AA = xlsread('Regression.xlsx');
DI = [0.701472275, 1.331443055, 2.239307889, 2.711813187, 1.32629228, 1.485504345, 1.73827467, 2.59502961, 0.48656241, 1.363499511, 0.72227, 0.756214019];
GG = [4.274825704, 2.189390368, 2.552240039, 1.365163176, 2.830760407, 2.573487374, 1.642279591, 1.906200052, 4.287215235, 2.723682117, 3.82240, 3.35589348];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
TT = [313, 313, 313, 313, 303, 303, 303, 303, 303, 303, 313, 313];
Tr = TT ./ 304;
RHO= [630, 691, 793, 840, 772, 802, 856, 891, 772, 891, 691, 793];

%%
Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

LabResults              = xlsread('wpd_datasets.xlsx');

which_k                 = [44, 46];

%% Load paramters
m_total                 = 3.5;

% Bed geometry
before                  = 0.04;                                             % Precentage of length before which is empty
bed                     = 0.92;                                             % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 600;
timeStep                = 10;                                               % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

SAMPLE                  = LabResults(21:34,1);

% Check if the number of data points is the same for both the dataset and the simulation
N_Sample                = [];
for i = 1:numel(SAMPLE)
    N_Sample            = [N_Sample ; find(round(Time,3) == round(SAMPLE(i))) ];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end

%% Specify parameters to estimate
nstages                 = Parameters{1};

nstagesbefore           = 1:floor(before*nstages);
nstagesbed              = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter            = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

%% Number of variables
Nx                      = 3 * nstages+2;                                    % 3*Nstages(C_f, C_s, H) + P(t) + yield
Nu                      = 3 + numel( Parameters );                          % T_in, P, F + numel(Parameters)

%% Extractor geometry
r                       = Parameters{3};                                    % Radius of the extractor  [m]
epsi                    = Parameters{4};                                    % Fullness [-]
L                       = Parameters{6};                                    % Total length of the extractor [m]

L_nstages               = linspace(0,L,nstages);
V                       = L  * pi * r^2;                                    % Total volume of the extractor [m3]
A                       = pi *      r^2;                                    % Extractor cross-section

%--------------------------------------------------------------------
V_slice                 = (L/nstages) * pi * r^2;

V_before                = V_slice * numel(nstagesbefore);
V_after                 = V_slice * numel(nstagesafter);
V_bed                   = V_slice * numel(nstagesbed);                      % Volume of the fixed bed [m3]

V_before_solid          = repmat(V_before * 0          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_solid             = repmat(V_bed    * epsi       / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_solid           = repmat(V_after  * 0          / numel(nstagesbed)   , numel(nstagesafter) ,1);

V_solid                 = [V_before_solid; V_bed_solid; V_after_solid];

V_before_fluid          = repmat(V_before * 1          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_fluid             = repmat(V_bed    * (1 - epsi) / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_fluid           = repmat(V_after  * 1          / numel(nstagesafter) , numel(nstagesafter) ,1);

V_fluid                 = [V_before_fluid; V_bed_fluid; V_after_fluid];

L_bed_after_nstages     = L_nstages(nstagesbed(1):end);
L_bed_after_nstages     = L_bed_after_nstages - L_bed_after_nstages(1);
L_end                   = L_bed_after_nstages(end);

%% symbolic variables
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

%% Set Integrator
f                       = @(x, u) modelSFE(x, u, bed_mask, timeStep_in_sec);
xdot                    = modelSFE(x, u, bed_mask, timeStep_in_sec);

% Integrator
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

%% Set inital state and inital conditions
msol_max                = m_total;                                          % g of product in solid and fluid phase
mSol_ratio              = 1;

mSOL_s                  = msol_max*mSol_ratio;                              % g of product in biomass
mSOL_f                  = msol_max*(1-mSol_ratio);                          % g of biomass in fluid

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                % Solid phase kg / m^3
Parameters{2}           = C0solid;

G                       =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

m_fluid                 = G(L_bed_after_nstages)*( L_bed_after_nstages(2) ); % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
m_fluid                 = [zeros(1,numel(nstagesbefore)) m_fluid];
C0fluid                 = m_fluid * 1e-3 ./ V_fluid';

for jj = 9:12
    which_dataset           = jj;

    data_org                = LabResults(21:34,which_dataset+1)';
    data_diff               = diff(data_org);

    %% Set operating conditions
    T0homog                 = LabResults(1,which_dataset+1);                    % K
    feedPress               = LabResults(2,which_dataset+1) * 10;               % MPa -> bar
    Flow                    = LabResults(3,which_dataset+1) * 1e-5 ;            % kg/s

    Z                       = Compressibility( T0homog, feedPress,         Parameters );
    rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );

    enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

    feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;     % Kelvin

    feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars

    feedFlow                = Flow * ones(1,length(Time_in_sec));               % kg/s

    uu                      = [feedTemp', feedPress', feedFlow'];

    % Initial conditions
    x0                      = [ C0fluid'                         ;
                                C0solid         * bed_mask       ;
                                enthalpy_rho    * ones(nstages,1);
                                feedPress(1)                     ;
                                0                                ;
                                ];

    %KOUT = DATA_K_OUT(:,ii);
    %KOUT = [ 0.7, 3.8 ];
    KOUT = [ DI(jj) , GG(jj) ];
    Parameters_opt = Parameters;
    for i=1:numel(which_k)
        Parameters_opt{which_k(i)}  = KOUT(i);
    end

    Parameters_init_time   = [uu repmat(cell2mat(Parameters_opt),1,N_Time)'];
    [xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );

    hold on; 
    plot(Time,xx_0(end,:),'LineWidth',2, 'DisplayName',[num2str(round(T0homog-273)),'$[^\circ C]$, ',num2str(feedPress(1)),'[bar]']); 
    plot(SAMPLE, data_org,'ko','LineWidth',2,'HandleVisibility','off'); 
    hold off
    
    xlabel('t~[min]')
    ylabel('$y~[g]$')

    set(gca,'FontSize',12)
end

legend boxoff 
lgd = legend('Location','northoutside', 'Orientation','horizontal');
lgd.FontSize = 10;

%exportgraphics(figure(1), ['Fit_Di_Gamma_9_12.png'], "Resolution",300);
%close all