startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*
%%
Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

LabResults              = xlsread('wpd_datasets.xlsx');
which_dataset           = 4;

SAMPLE                  = LabResults(21:34,1);
data_org                = LabResults(21:34,which_dataset+1)';

%% Load paramters
m_total                 = 3.0;

% Bed geometry
before                  = 0.04;                                             % Precentage of length before which is empty
bed                     = 0.92;                                             % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 15000;
timeStep                = 1;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

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
dp                      = Parameters{5};                                    % Paritcle diameter
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

RES                     = [];

for PP = [100,125,150,175,200]
    %% Set operating conditions
    T0homog                 = 35+273;                    % K
    feedPress               = PP;                        % bar
    Flow                    = 5 * 1e-5 ;                 % kg/s
    
    Z                       = Compressibility( T0homog, feedPress,         Parameters );
    rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );
    
    enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );
    
    feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;     % Kelvin
    
    feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars
    
    feedFlow                = Flow * ones(1,length(Time_in_sec));               % kg/s
    
    uu                      = [feedTemp', feedPress', feedFlow'];
    
    MU                      =     Viscosity(T0homog,rho);
    VELOCITY                =     Velocity(Flow, rho, Parameters);
    RE                      =     dp .* rho .* VELOCITY ./ MU;
    
    % Initial conditions
    x0                      = [ C0fluid'                         ;
                                C0solid         * bed_mask       ;
                                enthalpy_rho    * ones(nstages,1);
                                feedPress(1)                     ;
                                0                                ;
                                ];
    
    %% Set the inital simulation and plot it against the corresponding dataset
    %Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
    %[xx_0]                 = simulateSystem(F, [], x0, Parameters_init_time );
    
    %hold on
    %plot(Time,xx_0(end,:))
    %plot(SAMPLE,data_org,'o')
    %hold off
    
    %% Set sensitivity analysis
    
    name_v = {'T_{in}', 'P_{in}', 'F'                                   };
    name_s = {'c_f'   , 'c_s'   , '(h\times\rho)' , 'P_{t-1}'     , 'y' };
    name_p = {'CF'    , 'CS'    , 'H'             , 'P_{in}'      , 'Y' };
    
    My_Font = 14;
    num_levels = 100;
    
    for ii = 2     
    
            %% Sensitivities calculations
            Parameters{8} = ii;
            Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];
            [S,p,Sdot]              = Sensitivity(x, xdot, u, ii );
    
            % Initial conditions
            x0_SA                   = [ x0; zeros(length(S)-length(xdot),1) ];
    
            f_SA = @(S, p) Sdot(S, p, bed_mask);
            Results = Integrator_SS(Time*60, x0_SA, S, p, Sdot, Parameters_init_time);
            Res = Results(Nx+1:end,:);

            RES = [RES; Res(end,:)];
    
            %% Sensitivities plot 
            %{
            for ind=0:2
                figure()
                %subplot(2,1,1)
                imagesc(Time, L_nstages, Res(ind*nstages+1:(ind+1)*nstages,:)); cb = colorbar; colormap turbo;
                %subplot(2,1,2)
                %contourf(Time, L_nstages, flip(Res(ind*nstages+1:(ind+1)*nstages,:)),'EdgeColor','none','Levels',200); cb = colorbar; colormap turbo;
    
                hold on
                yline([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))],'k--');
                plot([ExtractionTime-40, ExtractionTime-40], [L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))], 'k--');
                text(ExtractionTime-55,[mean([L_nstages(nstagesbed(end)) L_nstages(nstagesbed(1))])],'fixed bed', 'Interpreter', 'latex', 'Color', 'black', 'HorizontalAlignment','center','VerticalAlignment','middle', 'Rotation', 90);
                hold off
    
                title(cb, ['$\frac{d',name_s{ind+1},'}{d',name_v{ii},'}$'], 'Interpreter', 'latex'); cb.TickLabelInterpreter = 'latex'; 
                cb.Label.Rotation = 0; % to rotate the text
                xlabel('Time [min]'); ylabel('Length [m]'); 
                set(gca,'FontSize',My_Font)
                %xscale log
                
                %exportgraphics(figure(1),[name_p{ind+1},'_',name_v{ii},'.png'], "Resolution",300);
                %close all;
            end
    
            %% Sensitivities plot - P and y
            indx = 0;
            for ind = 4:-1:3
                figure()
                hold on
                plot(Time, Res(end - indx,:), LineWidth=2); 
                %yline(0, LineWidth=2)
                xlabel('Time [min]'); ylabel(['$\frac{d ',name_s{ind+1},'}{d',name_v{ii},'}$'])
                hold off
                set(gca,'FontSize',My_Font)
                
                %exportgraphics(figure(1),[name_p{ind+1},'_',name_v{ii},'.png'], "Resolution",300);
                %close all;
                indx = indx + 1;
            end     
            %}   
    end
end
%%
plot(Time, RES','LineWidth',2)
legend('100 bar', '125 bar', '150 bar', '175 bar', '200 bar', 'Location','best')
legend box off
xlabel('Time [min]'); 
ylabel(['$\frac{d ',name_s{5},'}{d',name_v{2},'}$'])
set(gca,'FontSize',My_Font)
%exportgraphics(figure(1),['Yield_multiple.png'], "Resolution",300);
%close all;

%}









































