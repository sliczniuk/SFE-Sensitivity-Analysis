startup;
delete(gcp('nocreate'));

%%
AA = xlsread('Regression.xlsx');
DI = [0.701472275, 1.331443055, 2.239307889, 2.711813187, 1.32629228, 1.485504345, 1.73827467, 2.59502961, 0.48656241, 1.363499511, 0.72227, 0.756214019];
GG = [4.274825704, 2.189390368, 2.552240039, 1.365163176, 2.830760407, 2.573487374, 1.642279591, 1.906200052, 4.287215235, 2.723682117, 3.82240, 3.35589348];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
TT = [313, 313, 313, 313, 303, 303, 303, 303, 303, 303, 313, 313];
Tr = TT ./ 304;
RHO= [630, 691, 793, 840, 772, 802, 856, 891, 772, 891, 691, 793];

delta = 2.79 .* sqrt(7.38) .* Tr.^(1/4) .* (RHO./470);

FF = [6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 3.33, 3.33, 3.33, 3.33];

%%
x1 = RE(1:8) ; y1 = DI(1:8);
x2 = RE(9:12); y2 = DI(9:12);

[mdl1,stuff1] = fit(x1',y1','poly1');
[mdl2,stuff2] = fit(x2',y2','poly1');

p1 = coeffvalues(mdl1);
p2 = coeffvalues(mdl2);

xtest  = linspace(0.1,0.5);
ytest1 = polyval(p1,xtest);
ytest2 = polyval(p2,xtest);

hold on
plt1 = plot(x1,y1,'ko', 'LineWidth',2);
plt2 = plot(x2,y2,'kd', 'LineWidth',2);
plot(xtest,ytest1, 'LineWidth',2);
plot(xtest,ytest2, 'LineWidth',2);
yline(0, 'k', 'LineWidth',2);
hold off

%text(0.3, 4.0, sprintf('~$~\\Upsilon = %.3f \\cdot + Re %.3f$ \n $R^2=%.3f$',[p1, stuff1.rsquare]))
%text(0.15, 6.0, sprintf('~$~ \\Upsilon = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p2, stuff2.rsquare]))

text(0.35, 2.0, sprintf('$~ D_i^R = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p1, stuff1.rsquare]))
text(0.15, -1.0, sprintf('$~ D_i^R = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p2, stuff2.rsquare]))

legend([plt1(1),plt2(1)],'$6.67\times 10^{-5}[kg/s]$','$3.33\times 10^{-5}[kg/s]$','Location','best');
legend box off

ylabel('$D_i^R \times 10^{-13} [m^2/s]$')
%ylabel('$\Upsilon [-]$')
xlabel('Re [-]')
set(gca,'FontSize',12)

%exportgraphics(figure(1), ['Correlation_Di_Re.png'], "Resolution",300);
%close all

%% gamma function plot
%{\
C0 = 1;
XX = C0:-0.01:0;

figure()
for ii=1:12

    %RE_temp = sort(RE, 'descend');
    %jj = find(RE_temp(ii)==RE);

    f = @(x) (DI(ii) * 1e-13) * exp( -GG(ii) * (1-(x/(C0))) );    
    hold on
    %plot(XX, f(XX) ,'LineWidth', 2, 'DisplayName',['Re = ', num2str(round(RE(ii),2)), ', $\rho$ =', num2str(RHO(ii)), '$[kg/m^3]$'])
    plot(XX, f(XX) ,'LineWidth', 5, 'DisplayName', ['$\delta_H$ = ', num2str(round(delta(ii),2))])
    hold off
    text(XX(1)+0.01, f(XX(1)), ['$\delta_H$ = ', num2str(round(delta(ii),2)),'$[MPa^{1/2}]$, Re = ', num2str(round(RE(ii),2)), '[-]'], 'FontSize',28 );
end

%legend('Location','northwest')
%legend box off

xlim([0 1.6])
ylabel('$D_i^R \gamma(\Upsilon, C_s) [m^2/s] $')
xlabel('Normalized $C_s [kg/m^3]$')
set(gca,'FontSize',42)
%exportgraphics(figure(1), ['Gamma_function.png'], "Resolution",300);
%close all
%}

%% Fit a surface in RE and F space

figure()
[ml3, sf3] = fit([RHO', FF'],DI','poly11');
h3 = plot(ml3, [RHO', FF'],DI');
set(h3,'linestyle','none');
alpha(.5)

p3 = coeffvalues(ml3);
text(-0.25, 2.7, 0.1, sprintf('$D_i^R = %.3f %.3f \\cdot Re + %.3f \\cdot F$ \n $R^2 = %.3f $',[p3, sf3.rsquare]))
xlabel('Re[-]',Position=[0.12, 0.5, 0.7])
ylabel('$F \times 10^{-5}$[kg/s]',Position=[0, 4.5, -1.7])
zlabel('$D_i^R \times 10^{-13} [m^2/s]$')

exportgraphics(figure(1), ['Di_Re_F.png'], "Resolution",300);
close all

figure()
[ml4, sf4] = fit([RHO', FF'],GG','poly11');
h4 = plot(ml4, [RHO', FF'],GG');
set(h4,'linestyle','none');
alpha(.5)

p4 = coeffvalues(ml4);
text(-0.18, 3.8, 1, sprintf('$\\Upsilon ~~= %.3f + %.3f \\cdot Re %.3f \\cdot F$ \n $R^2 = %.3f $',[p4, sf4.rsquare]))
xlabel('Re[-]', Position=[0.3, 3, -0.57])
ylabel('$F \times 10^{-5}$[kg/s]',Position=[0, 4.5, 0.35])
zlabel('$\Upsilon[-]$')

exportgraphics(figure(1), ['Gamma_Re_F.png'], "Resolution",300);
close all




































