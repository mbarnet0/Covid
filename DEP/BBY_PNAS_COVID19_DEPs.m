% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%DEP Simulations 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

close all
clear all
clc

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

set(0,'defaulttextinterpreter','tex')
set(0, 'defaultAxesTickLabelInterpreter','tex'); 
set(0, 'defaultLegendInterpreter','tex');

filename_loada = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_41']; % equal-weighted prior (4), uncertainty neutral (1)
filename_loadb = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_42']; % equal-weighted prior (4), uncertainty averse (2)

filename_keep = filename_loada;
filename_save= [filename_loada,'_deps'];

A=load(filename_loada);
B=load(filename_loadb);

n=0;
rng(seed);


%Number of shocks
nShocks = 3;
    
%Simulation with policy function
Nyears = 1 ; % for simulations and figures
Nmonths = 12*Nyears ; % for simulations and figures
periods = 30;
time= periods*Nmonths ; 

period=[1:1:time];

dts = 1/periods;
% dts = 1;

R0 = A.R_0;


policy_L_base=A.L_optimal; 
policy_L_tilde=B.L_optimal;

S_mat = B.S_mat;
I_mat = B.I_mat;
pi_tilde = B.pi_tilde;
weights = A.weights;
beta_vals = B.beta_vals;
delta_vals = B.delta_vals;
theta_vals = B.theta_vals;
aa_vals = B.aa_vals;
ggamma = B.ggamma;
varsigma = B.varsigma;

sigma_i = B.sigma_i;
sigma_d = B.sigma_d;
sigma_z = B.sigma_z;

pol_interp_base = griddedInterpolant(S_mat,I_mat,policy_L_base,'linear');  
pol_interp_tilde = griddedInterpolant(S_mat,I_mat,policy_L_tilde,'linear');  

pitlde1_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,1),'linear');  
pitlde2_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,2),'linear');  
pitlde3_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,3),'linear');  
pitlde4_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,4),'linear');  
pitlde5_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,5),'linear');  
pitlde6_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,6),'linear');  
pitlde7_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,7),'linear');  
pitlde8_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,8),'linear');  
pitlde9_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,9),'linear');  
pitlde10_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,10),'linear');  
pitlde11_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,11),'linear');  
pitlde12_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,12),'linear');  
pitlde13_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,13),'linear');  
pitlde14_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,14),'linear');  
pitlde15_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,15),'linear');  
pitlde16_interp = griddedInterpolant(S_mat,I_mat,pi_tilde(:,:,16),'linear'); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

base_beta =   beta_vals'*weights;                    
base_delta =  delta_vals'*weights;    
base_theta =  theta_vals'*weights;         
base_aa =     aa_vals'*weights; 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
              
alt_beta =   @(x1,x2)   beta_vals(1).*pitlde1_interp(x1,x2)+beta_vals(2).*pitlde2_interp(x1,x2)+beta_vals(3).*pitlde3_interp(x1,x2)...
             +beta_vals(4).*pitlde4_interp(x1,x2)+beta_vals(5).*pitlde5_interp(x1,x2)+beta_vals(6).*pitlde6_interp(x1,x2)...
             +beta_vals(7).*pitlde7_interp(x1,x2)+beta_vals(8).*pitlde8_interp(x1,x2)+beta_vals(9).*pitlde9_interp(x1,x2)...
             +beta_vals(10).*pitlde10_interp(x1,x2)...
             +beta_vals(11).*pitlde11_interp(x1,x2)+beta_vals(12).*pitlde12_interp(x1,x2)+beta_vals(13).*pitlde13_interp(x1,x2)...
             +beta_vals(14).*pitlde14_interp(x1,x2)+beta_vals(15).*pitlde15_interp(x1,x2)+beta_vals(16).*pitlde16_interp(x1,x2);

alt_delta =   @(x1,x2)  delta_vals(1).*pitlde1_interp(x1,x2)+delta_vals(2).*pitlde2_interp(x1,x2)+delta_vals(3).*pitlde3_interp(x1,x2)...
             +delta_vals(4).*pitlde4_interp(x1,x2)+delta_vals(5).*pitlde5_interp(x1,x2)+delta_vals(6).*pitlde6_interp(x1,x2)...
             +delta_vals(7).*pitlde7_interp(x1,x2)+delta_vals(8).*pitlde8_interp(x1,x2)+delta_vals(9).*pitlde9_interp(x1,x2) ...
             +delta_vals(10).*pitlde10_interp(x1,x2)...
             +delta_vals(11).*pitlde11_interp(x1,x2)+delta_vals(12).*pitlde12_interp(x1,x2)+delta_vals(13).*pitlde13_interp(x1,x2)...
             +delta_vals(14).*pitlde14_interp(x1,x2)+delta_vals(15).*pitlde15_interp(x1,x2)+delta_vals(16).*pitlde16_interp(x1,x2);

alt_theta =   @(x1,x2)  theta_vals(1).*pitlde1_interp(x1,x2)+theta_vals(2).*pitlde2_interp(x1,x2)+theta_vals(3).*pitlde3_interp(x1,x2)...
             +theta_vals(4).*pitlde4_interp(x1,x2)+theta_vals(5).*pitlde5_interp(x1,x2)+theta_vals(6).*pitlde6_interp(x1,x2)...
             +theta_vals(7).*pitlde7_interp(x1,x2)+theta_vals(8).*pitlde8_interp(x1,x2)+theta_vals(9).*pitlde9_interp(x1,x2)...
             +theta_vals(10).*pitlde10_interp(x1,x2)...
             +theta_vals(11).*pitlde11_interp(x1,x2)+theta_vals(12).*pitlde12_interp(x1,x2)+theta_vals(13).*pitlde13_interp(x1,x2)...
             +theta_vals(14).*pitlde14_interp(x1,x2)+theta_vals(15).*pitlde15_interp(x1,x2)+theta_vals(16).*pitlde16_interp(x1,x2);

alt_aa =    @(x1,x2)    aa_vals(1).*pitlde1_interp(x1,x2)+aa_vals(2).*pitlde2_interp(x1,x2)+aa_vals(3).*pitlde3_interp(x1,x2)...
             +aa_vals(4).*pitlde4_interp(x1,x2)+aa_vals(5).*pitlde5_interp(x1,x2)+aa_vals(6).*pitlde6_interp(x1,x2)...
             +aa_vals(7).*pitlde7_interp(x1,x2)+aa_vals(8).*pitlde8_interp(x1,x2)+aa_vals(9).*pitlde9_interp(x1,x2)...
             +aa_vals(10).*pitlde10_interp(x1,x2)...
             +aa_vals(11).*pitlde11_interp(x1,x2)+aa_vals(12).*pitlde12_interp(x1,x2)+aa_vals(13).*pitlde13_interp(x1,x2)...
             +aa_vals(14).*pitlde14_interp(x1,x2)+aa_vals(15).*pitlde15_interp(x1,x2)+aa_vals(16).*pitlde16_interp(x1,x2);
            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% stop
sims_iters = num_its;


S_initial = B.Stime_Lck(seed);
I_initial = B.Itime_Lck(seed);
N_initial = B.Ntime_Lck(seed);
z_initial = B.ztime_Lck(seed);

rA = zeros(1,sims_iters);
rB = zeros(1,sims_iters);

tic;

pc = parcluster('local');
pc.JobStorageLocation = strcat('/home/user/');
% pc.JobStorageLocation = strcat('/home/user/',getenv('SLURM_JOB_ID'));


parfor sims=1:sims_iters

Stime=zeros(time,1); 
Itime=zeros(time,1); 
Ntime=zeros(time,1);
ztime=zeros(time,1); 

% initial conditions 
Stime(1)= S_initial; 
Itime(1) = I_initial;
Ntime(1) = N_initial;
ztime(1) = z_initial;


Stime_tilde=zeros(time,1); 
Itime_tilde=zeros(time,1); 
Ntime_tilde=zeros(time,1);  
ztime_tilde=zeros(time,1);  

% Same initial conditions
Stime_tilde(1)=Stime(1); 
Itime_tilde(1)=Itime(1); 
Ntime_tilde(1)=Ntime(1); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Lck=zeros(time,1); 
Lck_p=zeros(time,1); 
Lck_tilde=zeros(time,1); 
Lck_tilde_p=zeros(time,1);

r_A=zeros(time,1); 
r_B=zeros(time,1); 

      
% tic;

for t=2:time 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % baseline process  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

shock = normrnd(0,1, nShocks, 1);
sigma_S = [-Stime(t-1).*Itime(t-1).*sigma_i, Stime(t-1).*Itime(t-1).*sigma_d, 0].*sqrt(dts);
sigma_I = [ Stime(t-1).*Itime(t-1).*sigma_i, - Itime(t-1).*(1-Itime(t-1)).*sigma_d, 0].*sqrt(dts);
sigma_N = [0, -Itime(t-1).*sigma_d.*Ntime(t-1), 0].*sqrt(dts);
sigma_Z = [0,0,sigma_z].*sqrt(dts);
sigma = [sigma_S;sigma_I;sigma_N;sigma_Z];
sigma_p = sigma';
sigma_p_sigma_inv = pinv(sigma_p*sigma);
sigma_p_sigma_inv_p = sigma_p_sigma_inv';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Lck(t-1) = pol_interp_base(Stime(t-1),Itime(t-1)); 

    Stime(t) =  Stime(t-1) - base_beta*Stime(t-1)*Itime(t-1)*( 1-base_theta*Lck(t-1))^2*dts ...
        - (-Itime(t-1)*(base_delta+5.*base_delta.*Itime(t-1)))*Stime(t-1)*dts ...
        + sigma_S*shock;
    Itime(t) =  Itime(t-1) + base_beta*Stime(t-1)*Itime(t-1)*( 1-base_theta*Lck(t-1))^2*dts  ...
        - (-Itime(t-1)*(base_delta+5.*base_delta.*Itime(t-1)))*Itime(t-1)*dts ...
        - ggamma*Itime(t-1)*dts + sigma_I*shock;    
    Ntime(t) =  Ntime(t-1) - Itime(t-1)*(base_delta+5.*base_delta.*Itime(t-1)).*Ntime(t-1)*dts  + sigma_N*shock;
    ztime(t) =  ztime(t-1) - base_aa*Lck(t-1)*varsigma*dts - ztime(t-1)*varsigma*dts  + sigma_Z*shock;   
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    Stimep = - base_beta*Stime(t-1)*Itime(t-1)*( 1-base_theta*Lck(t-1))^2*dts ...
             - (-Itime(t-1)*(base_delta+5.*base_delta.*Itime(t-1)))*Stime(t-1)*dts;
    Itimep =   base_beta*Stime(t-1)*Itime(t-1)*( 1-base_theta*Lck(t-1))^2*dts...
             - (-Itime(t-1)*(base_delta+5.*base_delta.*Itime(t-1)))*Itime(t-1)*dts ...
             - ggamma*Itime(t-1)*dts;    
    Ntimep = - Itime(t-1)*(base_delta+5.*base_delta.*Itime(t-1)).*Ntime(t-1)*dts ;
    ztimep = - base_aa*Lck(t-1)*varsigma*dts - ztime(t-1)*varsigma*dts;  
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Lck_p(t-1) = pol_interp_tilde(Stime(t-1),Itime(t-1)); 

beta_alt = round(alt_beta(Stime(t-1),Itime(t-1)),7);
delta_alt = round(alt_delta(Stime(t-1),Itime(t-1)),7);
theta_alt = round(alt_theta(Stime(t-1),Itime(t-1)),7);
aa_alt = round(alt_aa(Stime(t-1),Itime(t-1)),7);

    Stime_p = - beta_alt*Stime(t-1)*Itime(t-1)*( 1-theta_alt*Lck_p(t-1))^2*dts ...
              - (-Itime(t-1)*(delta_alt+5.*delta_alt.*Itime(t-1)))*Stime(t-1)*dts ;
    Itime_p =   beta_alt*Stime(t-1)*Itime(t-1)*( 1-theta_alt*Lck_p(t-1))^2*dts ...
              - (-Itime(t-1)*(delta_alt+5.*delta_alt.*Itime(t-1)))*Itime(t-1)*dts ...
              - ggamma*Itime(t-1)*dts ;    
    Ntime_p = - Itime(t-1)*(delta_alt+5.*delta_alt.*Itime(t-1)).*Ntime(t-1)*dts ;
    ztime_p = - aa_alt*Lck_p(t-1)*varsigma*dts - ztime(t-1)*varsigma*dts ;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % worst-case process
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

shock_tilde = normrnd(0,1, nShocks, 1);
sigma_S_tilde = [-Stime_tilde(t-1).*Itime_tilde(t-1).*sigma_i, Stime_tilde(t-1).*Itime_tilde(t-1).*sigma_d, 0].*sqrt(dts);
sigma_I_tilde = [ Stime_tilde(t-1).*Itime_tilde(t-1).*sigma_i, - Itime_tilde(t-1).*(1-Itime_tilde(t-1)).*sigma_d, 0].*sqrt(dts);
sigma_N_tilde = [0, -Itime_tilde(t-1).*sigma_d.*Ntime_tilde(t-1), 0].*sqrt(dts);
sigma_Z_tilde = [0,0,sigma_z].*sqrt(dts);
sigma_tilde = [sigma_S_tilde;sigma_I_tilde;sigma_N_tilde;sigma_Z_tilde];
sigma_tilde_p = sigma_tilde';
sigma_tilde_p_sigma_tilde_inv = pinv(sigma_tilde_p*sigma_tilde);
sigma_tilde_p_sigma_tilde_inv_p = sigma_tilde_p_sigma_tilde_inv';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

Lck_tilde(t-1) = pol_interp_tilde(Stime_tilde(t-1),Itime_tilde(t-1)); 
beta_alt = round(alt_beta(Stime(t-1),Itime(t-1)),7);
delta_alt = round(alt_delta(Stime(t-1),Itime(t-1)),7);
theta_alt = round(alt_theta(Stime(t-1),Itime(t-1)),7);
aa_alt = round(alt_aa(Stime(t-1),Itime(t-1)),7);

    Stime_tilde(t) =  Stime_tilde(t-1) - beta_alt*Stime_tilde(t-1)*Itime_tilde(t-1)*( 1-theta_alt*Lck_tilde(t-1))^2*dts ...
        - (-Itime_tilde(t-1)*(delta_alt+5.*delta_alt.*Itime_tilde(t-1)))*Stime_tilde(t-1)*dts ...
        + sigma_S_tilde*shock_tilde;
    Itime_tilde(t) =  Itime_tilde(t-1) + beta_alt*Stime_tilde(t-1)*Itime_tilde(t-1)*( 1-theta_alt*Lck_tilde(t-1))^2*dts ...
        - (-Itime_tilde(t-1)*(delta_alt+5.*delta_alt.*Itime_tilde(t-1)))*Itime_tilde(t-1)*dts ...
        - ggamma*Itime_tilde(t-1)*dts + sigma_I_tilde*shock_tilde;        
    Ntime_tilde(t) =  Ntime_tilde(t-1) - Itime_tilde(t-1)*(delta_alt+5.*delta_alt.*Itime_tilde(t-1)).*Ntime_tilde(t-1)*dts  + sigma_N_tilde*shock_tilde;
    ztime_tilde(t) =  ztime_tilde(t-1) - aa_alt*Lck_tilde(t-1)*varsigma*dts - ztime_tilde(t-1)*varsigma*dts + sigma_Z_tilde*shock_tilde;    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

    Stime_tildep =  - beta_alt*Stime_tilde(t-1)*Itime_tilde(t-1)*( 1-theta_alt*Lck_tilde(t-1))^2*dts ...
                    - (-Itime_tilde(t-1)*(delta_alt+5.*delta_alt.*Itime_tilde(t-1)))*Stime_tilde(t-1)*dts;
    Itime_tildep =    beta_alt*Stime_tilde(t-1)*Itime_tilde(t-1)*( 1-theta_alt*Lck_tilde(t-1))^2*dts ...
                    - (-Itime_tilde(t-1)*(delta_alt+5.*delta_alt.*Itime_tilde(t-1)))*Itime_tilde(t-1)*dts ...
                    - ggamma*Itime_tilde(t-1)*dts;
    Ntime_tildep =  - Itime_tilde(t-1)*(delta_alt+5.*delta_alt.*Itime_tilde(t-1)).*Ntime_tilde(t-1)*dts;
    ztime_tildep =  - aa_alt*Lck_tilde(t-1)*varsigma*dts - ztime_tilde(t-1)*varsigma*dts;       
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

Lck_tilde_p(t-1) = pol_interp_base(Stime_tilde(t-1),Itime_tilde(t-1)); 

    Stime_tilde_p =  - base_beta*Stime_tilde(t-1)*Itime_tilde(t-1)*( 1-base_theta*Lck_tilde_p(t-1))^2*dts ...
                     - (-Itime_tilde(t-1)*(base_delta+5.*base_delta.*Itime_tilde(t-1)))*Stime_tilde(t-1)*dts ;
    Itime_tilde_p =    base_beta*Stime_tilde(t-1)*Itime_tilde(t-1)*( 1-base_theta*Lck_tilde_p(t-1))^2*dts ...
                     - (-Itime_tilde(t-1)*(base_delta+5.*base_delta.*Itime_tilde(t-1)))*Itime_tilde(t-1)*dts ...
                     - ggamma*Itime_tilde(t-1)*dts;    
    Ntime_tilde_p =  - Itime_tilde(t-1)*(base_delta+5.*base_delta.*Itime_tilde(t-1)).*Ntime_tilde(t-1)*dts ;
    ztime_tilde_p =  - base_aa*Lck_tilde_p(t-1)*varsigma*dts - ztime_tilde(t-1)*varsigma*dts ;
  
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

delta_mu_A = [Stime_p-Stimep;Itime_p-Itimep;Ntime_p-Ntimep;ztime_p-ztimep];
h_val_A = sigma_p_sigma_inv*sigma_p*delta_mu_A.*(dts);
h_A_p_h_A = (h_val_A')*h_val_A;

ra = 0.5*h_A_p_h_A - h_val_A'*shock;
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

delta_mu_B = -[Stime_tilde_p-Stime_tildep;Itime_tilde_p-Itime_tildep;Ntime_tilde_p-Ntime_tildep;ztime_tilde_p-ztime_tildep];
h_val_B = sigma_tilde_p_sigma_tilde_inv*sigma_tilde_p*delta_mu_B.*(dts);
h_B_p_h_B = (h_val_B')*h_val_B;

rb = 0.5*h_B_p_h_B + h_val_B'*shock_tilde;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
        
    r_A(t) = r_A(t-1) + ra*(1./(time));
    r_B(t) = r_B(t-1) + rb*(1./(time));    

end


% toc;

rA(sims)=r_A(end);
rB(sims)=r_B(end);

% pause
sims

end

toc;

R_A = rA;
R_B = rB;


p_A = R_A<0;
p_B = R_B<0;

dep = mean([mean(p_A',1); mean(p_B',1)],1);

save(filename_save,'R_A','R_B','sims_iters','p_A','p_B','dep')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
