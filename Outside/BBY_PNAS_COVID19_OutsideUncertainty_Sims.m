% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%Outside Uncertainty Simulations 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
close all
clear all
clc

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

set(0,'defaulttextinterpreter','tex')
set(0, 'defaultAxesTickLabelInterpreter','tex'); 
set(0, 'defaultLegendInterpreter','tex');

for jj = 1:16

    case_val = jj;

    direct = pwd;
    directory = [direct];

    if case_val == 1
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1111']); 
    elseif case_val ==  2
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2111']); 
    elseif case_val == 3
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1211']); 
    elseif case_val == 4
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1121']); 
    elseif case_val == 5
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1112']); 
    elseif case_val == 6
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2211']); 
    elseif case_val == 7
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2121']); 
    elseif case_val == 8
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2112']); 
    elseif case_val == 9
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1211']); 
    elseif case_val == 10
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1212']); 
    elseif case_val == 11
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1122']); 
    elseif case_val == 12
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2221']); 
    elseif case_val == 13
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2122']); 
    elseif case_val == 14
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1222']); 
    elseif case_val == 15
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2212']); 
    elseif case_val == 16
        filename2 =([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2222']); 
        
    end  


filename_keep = filename2;

load(filename2);
directory = pwd;

filename3 = [filename_keep,'_sims'];

tic;

%Simulation with policy function
 
Nyears = 2 ; % for simulations and figures

time= 365*Nyears ; 
period=[1:1:time];

dts = 1/30;

R0 = R_0;

Stime_Lck=zeros(time,1); 
Itime_Lck=zeros(time,1); 
Rtime_Lck=zeros(time,1); 
Ntime_Lck=zeros(time,1);  
Dtime_Lck=zeros(time,1); 
ztime_Lck=zeros(time,1); 

% initial conditions 
Itime_Lck(1)=0.02;
Rtime_Lck(1)=0.03;
Stime_Lck(1)=1-Itime_Lck(1)-Rtime_Lck(1); 
Ntime_Lck(1) = Stime_Lck(1) + Itime_Lck(1) + Rtime_Lck(1);

policy_L=L_optimal; 

pol_interp = griddedInterpolant(S_mat,I_mat,policy_L,'linear');  
% pol_interp = griddedInterpolant(S_mat,I_mat,policy_L,'cubic');


Lck(1) =  pol_interp(Stime_Lck(1),Itime_Lck(1));

for t=2:time 

    Stime_Lck(t) =  Stime_Lck(t-1) - beta*Stime_Lck(t-1)*Itime_Lck(t-1)*( 1-theta*Lck(t-1))^2*dts - (-Itime_Lck(t-1)*(delta+5.*delta.*Itime_Lck(t-1)))*Stime_Lck(t-1)*dts ;
    Itime_Lck(t) =  Itime_Lck(t-1) - ggamma*Itime_Lck(t-1)*dts - (-Itime_Lck(t-1)*(delta+5.*delta.*Itime_Lck(t-1)))*Itime_Lck(t-1)*dts + beta*Stime_Lck(t-1)*Itime_Lck(t-1)*( 1-theta*Lck(t-1))^2*dts ;    
    Dtime_Lck(t) =  Dtime_Lck(t-1) + Itime_Lck(t-1)*(delta+5.*delta.*Itime_Lck(t-1))*dts  ;
    Rtime_Lck(t) =  Rtime_Lck(t-1) + Itime_Lck(t-1)*(ggamma - (delta+5.*delta.*Itime_Lck(t-1)) )*dts ;
    ztime_Lck(t) =  ztime_Lck(t-1) - aa*Lck(t-1)*varsigma*dts - ztime_Lck(t-1)*varsigma*dts ;
    Ntime_Lck(t) =  Stime_Lck(t)   + Itime_Lck(t) + Rtime_Lck(t);
    
    Lck(t) = pol_interp(Stime_Lck(t),Itime_Lck(t)) ; 
    
    
end

toc;


save(filename3)


end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%Outside Uncertainty Simulation Plots
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

close all
clear all
clc
direct = pwd;
directory = [pwd];

set(0,'defaulttextinterpreter','tex')
set(0, 'defaultAxesTickLabelInterpreter','tex'); 
set(0, 'defaultLegendInterpreter','tex');

 
Aa=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1111','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights','time','Nyears','dts');
Bb=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1121','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Cc=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1112','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Dd=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1122','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');

Ee=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2111','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Ff=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2121','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Gg=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2112','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Hh=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2122','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');

Ii=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1211','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Jj=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1221','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Kk=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1212','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Ll=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1222','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');

Mm=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2211','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Nn=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2221','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Oo=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2212','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');
Pp=load([directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_2222','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','beta','delta','theta','aa','ggamma','weights');

color1 = [0    0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];
color4 = [0.4940    0.1840    0.5560];
color5 = [0.4660    0.6740    0.1880];
color6 = [0.3010    0.7450    0.9330];
color7 = [0.6350    0.0780    0.1840];

scale = 2;

time = Aa.time;
time_alt = Aa.time./scale;
Nyears = Aa.Nyears; 
Nyears_alt = Aa.Nyears./scale;
dt = Aa.dts;
time_vec = linspace(0,Nyears_alt,time_alt);
time_vec2 = linspace(0,52*Nyears_alt,time_alt);
Alt_end = round(length(Aa.Stime_Lck)./scale);

% stop



figure;
hold on;
plot(time_vec2,(Aa.Dtime_Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Aa.theta,Aa.aa))
plot(time_vec2,(Bb.Dtime_Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Bb.theta,Bb.aa))
plot(time_vec2,(Cc.Dtime_Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Cc.theta,Cc.aa))
plot(time_vec2,(Dd.Dtime_Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Dd.theta,Dd.aa))
title(sprintf('\\textbf{Dead:} $\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Aa.R0,Aa.CFR), 'Interpreter', 'latex')
ylabel('$$d_t$$', 'Interpreter', 'latex')
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.04])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northwest');
print('covid_RR_D1_SA','-dpng')

figure;
hold on;
plot(time_vec2,(Ee.Dtime_Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ee.theta,Ee.aa))
plot(time_vec2,(Ff.Dtime_Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ff.theta,Ff.aa))
plot(time_vec2,(Gg.Dtime_Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Gg.theta,Gg.aa))
plot(time_vec2,(Hh.Dtime_Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Hh.theta,Hh.aa))
title(sprintf('\\textbf{Dead:} $\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Ee.R0,Ee.CFR), 'Interpreter', 'latex')
ylabel('$$d_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.04])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northwest');
print('covid_RR_D2_SA','-dpng')

figure;
hold on;
plot(time_vec2,(Ii.Dtime_Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ii.theta,Ii.aa))
plot(time_vec2,(Jj.Dtime_Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Jj.theta,Jj.aa))
plot(time_vec2,(Kk.Dtime_Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Kk.theta,Kk.aa))
plot(time_vec2,(Ll.Dtime_Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ll.theta,Ll.aa))
title(sprintf('\\textbf{Dead:} $\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Ii.R0,Ii.CFR), 'Interpreter', 'latex')
ylabel('$$d_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.04])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northwest');
print('covid_RR_D3_SA','-dpng')

figure;
hold on;
plot(time_vec2,(Mm.Dtime_Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Mm.theta,Mm.aa))
plot(time_vec2,(Nn.Dtime_Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Nn.theta,Nn.aa))
plot(time_vec2,(Oo.Dtime_Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Oo.theta,Oo.aa))
plot(time_vec2,(Pp.Dtime_Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Pp.theta,Pp.aa))
title(sprintf('\\textbf{Dead:} $\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Mm.R0,Mm.CFR), 'Interpreter', 'latex')
ylabel('$$d_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.04])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'southeast');
print('covid_RR_D4_SA','-dpng')



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure;
hold on;
plot(time_vec2,(Aa.Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Aa.theta,Aa.aa))
plot(time_vec2,(Bb.Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Bb.theta,Bb.aa))
plot(time_vec2,(Cc.Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Cc.theta,Cc.aa))
plot(time_vec2,(Dd.Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Dd.theta,Dd.aa))
title(sprintf('\\textbf{Quarantine:} $$\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Aa.R0,Aa.CFR), 'Interpreter', 'latex')
ylabel('$$q_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.85])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Q1_SA','-dpng')

figure;
hold on;
plot(time_vec2,(Ee.Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ee.theta,Ee.aa))
plot(time_vec2,(Ff.Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ff.theta,Ff.aa))
plot(time_vec2,(Gg.Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Gg.theta,Gg.aa))
plot(time_vec2,(Hh.Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Hh.theta,Hh.aa))
title(sprintf('\\textbf{Quarantine:} $$\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Ee.R0,Ee.CFR), 'Interpreter', 'latex')
ylabel('$$q_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.85])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Q2_SA','-dpng')

figure;
hold on;
plot(time_vec2,(Ii.Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ii.theta,Ii.aa))
plot(time_vec2,(Jj.Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Jj.theta,Jj.aa))
plot(time_vec2,(Kk.Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Kk.theta,Kk.aa))
plot(time_vec2,(Ll.Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Ll.theta,Ll.aa))
title(sprintf('\\textbf{Quarantine:} $$\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Ii.R0,Ii.CFR), 'Interpreter', 'latex')
ylabel('$$q_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.85])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Q3_SA','-dpng')

figure;
hold on;
plot(time_vec2,(Mm.Lck(1:Alt_end)),'Color',color1,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Mm.theta,Mm.aa))
plot(time_vec2,(Nn.Lck(1:Alt_end)),'Color',color2,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Nn.theta,Nn.aa))
plot(time_vec2,(Oo.Lck(1:Alt_end)),'Color',color3,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Oo.theta,Oo.aa))
plot(time_vec2,(Pp.Lck(1:Alt_end)),'Color',color4,'LineWidth',3, 'DisplayName', sprintf('$\\zeta = %0.2f, \\alpha = %0.1f$',Pp.theta,Pp.aa))
title(sprintf('\\textbf{Quarantine:} $$\\mathcal{R}_0 =$ %0.1f, $CFR =$ %0.3f',Mm.R0,Mm.CFR), 'Interpreter', 'latex')
ylabel('$$q_t$$', 'Interpreter', 'latex') 
xlabel('Weeks', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
axis([0 52*Nyears_alt 0 0.85])
xticks([0 10 20 30 40 50])
ax = gca;
ax.YAxis.Exponent = 0;
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Q4_SA','-dpng')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


