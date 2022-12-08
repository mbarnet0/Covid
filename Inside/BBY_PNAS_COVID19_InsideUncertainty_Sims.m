% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%Inside Uncertainty Simulations 
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


for ll = 1:4
    case_val = ll;
    
for jj = 1:2
    uncertain_val = jj;

direct = pwd;
directory = [direct];    


       if uncertain_val == 1
 

        if case_val == 1              
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_11'];   
        elseif case_val == 2
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_21'];  
        elseif case_val == 3             
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_31'];              
        elseif case_val == 4
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_41'];             
        end
        
    elseif uncertain_val == 2
        
        if case_val == 1              
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_12'];              
        elseif case_val == 2
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_22'];
        elseif case_val == 3             
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_32'];               
        elseif case_val == 4
            file_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_42'];               
        end
    end

filename_keep = file_load;

load(file_load);
directory = pwd;


filename_save = [filename_keep,'_sims'];

tic;
%Simulation with policy function

Nyears = 2 ; % for simulations and figures

time= 365*Nyears ; 
period=[1:1:time];

dts = 1/30;

R0 = R_0;


Stime=zeros(time,1); 
Itime=zeros(time,1); 
Rtime=zeros(time,1); 
Ntime=zeros(time,1);
Dtime=zeros(time,1); 
ztime=zeros(time,1); 



% initial conditions 
Itime(1) = 0.02;
Rtime(1) = 0.03;
Stime(1)= 1-Itime(1)-Rtime(1); 
Lck=[];

% Ntime(1) = Stime(1) + Itime(1) + Rtime(1);


Stime_Lck=zeros(time,1); 
Itime_Lck=zeros(time,1); 
Rtime_Lck=zeros(time,1); 
Ntime_Lck=zeros(time,1);  
Dtime_Lck=Dtime;
ztime_Lck=zeros(time,1);  


% Same initial conditions
Itime_Lck(1)=Itime(1); 
Stime_Lck(1)=Stime(1); 
Rtime_Lck(1)=Rtime(1); 

Ntime_Lck(1) = Stime_Lck(1) + Itime_Lck(1) + Rtime_Lck(1);

policy_L=L_optimal; 

pol_interp = griddedInterpolant(S_mat,I_mat,policy_L,'linear');  

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

         
alt_beta =    (beta_vals(1).*pi_tilde(:,:,1)+beta_vals(2).*pi_tilde(:,:,2)+beta_vals(3).*pi_tilde(:,:,3)...
             +beta_vals(4).*pi_tilde(:,:,4)+beta_vals(5).*pi_tilde(:,:,5)+beta_vals(6).*pi_tilde(:,:,6)...
             +beta_vals(7).*pi_tilde(:,:,7)+beta_vals(8).*pi_tilde(:,:,8)+beta_vals(9).*pi_tilde(:,:,9)...
             +beta_vals(10).*pi_tilde(:,:,10)...
             +beta_vals(11).*pi_tilde(:,:,11)+beta_vals(12).*pi_tilde(:,:,12)+beta_vals(13).*pi_tilde(:,:,13)...
             +beta_vals(14).*pi_tilde(:,:,14)+beta_vals(15).*pi_tilde(:,:,15)+beta_vals(16).*pi_tilde(:,:,16))./ggamma;
alt_R0 = griddedInterpolant(S_mat,I_mat,alt_beta,'linear');  

I_hat_interp = griddedInterpolant(S_mat,I_mat,I_hat,'linear'); 

Lck(1) =  pol_interp(Stime_Lck(1),Itime_Lck(1));

PiTilde1(1) =  pitlde1_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde2(1) =  pitlde2_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde3(1) =  pitlde3_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde4(1) =  pitlde4_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde5(1) =  pitlde5_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde6(1) =  pitlde6_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde7(1) =  pitlde7_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde8(1) =  pitlde8_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde9(1) =  pitlde9_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde10(1) =  pitlde10_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde11(1) =  pitlde11_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde12(1) =  pitlde12_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde13(1) =  pitlde13_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde14(1) =  pitlde14_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde15(1) =  pitlde15_interp(Stime_Lck(1),Itime_Lck(1));
PiTilde16(1) =  pitlde16_interp(Stime_Lck(1),Itime_Lck(1));

Alt_R0(1) =  alt_R0(Stime_Lck(1),Itime_Lck(1));

Ihat(1) = I_hat_interp(Stime_Lck(1),Itime_Lck(1)); 

beta = mean(beta_vals);
delta = mean(delta_vals);
theta = mean(theta_vals);
aa = mean(aa_vals);


for t=2:time 
    % benchmark uncontrolled process
    Stime(t) =  Stime(t-1) - beta*Stime(t-1)*Itime(t-1)*(1-0)*dts - (-Itime(t-1)*(delta+5.*delta.* Itime(t-1)))*Stime(t-1)*dts ;
    Itime(t) =  Itime(t-1) - ggamma*Itime(t-1)*dts - (-Itime(t-1)*(delta+5.*delta.* Itime(t-1)))*Itime(t-1)*dts + beta*Stime(t-1)*Itime(t-1)*(1-0)*dts ;
    Dtime(t) =  Dtime(t-1) + Itime(t-1)*(delta+5.*delta.* Itime(t-1))*dts ;
    Rtime(t) =  Rtime(t-1) + Itime(t-1)*(ggamma - (delta+5.*delta.* Itime(t-1)) )*dts ;
    ztime(t) =  ztime(t-1) - ztime(t-1)*varsigma*dts ;
    Ntime(t) =  Stime(t) + Itime(t) + Rtime(t); 
    
    % Controlled process
    % fatality ctrolled

    Stime_Lck(t) =  Stime_Lck(t-1) - beta*Stime_Lck(t-1)*Itime_Lck(t-1)*( 1-theta*Lck(t-1))^2*dts - (-Itime_Lck(t-1)*(delta+5.*delta.*Itime_Lck(t-1)))*Stime_Lck(t-1)*dts ;
    Itime_Lck(t) =  Itime_Lck(t-1) - ggamma*Itime_Lck(t-1)*dts - (-Itime_Lck(t-1)*(delta+5.*delta.*Itime_Lck(t-1)))*Itime_Lck(t-1)*dts + beta*Stime_Lck(t-1)*Itime_Lck(t-1)*( 1-theta*Lck(t-1))^2*dts ;    
    Dtime_Lck(t) =  Dtime_Lck(t-1) + Itime_Lck(t-1)*(delta+5.*delta.*Itime_Lck(t-1))*dts  ;
    Rtime_Lck(t) =  Rtime_Lck(t-1) + Itime_Lck(t-1)*(ggamma - (delta+5.*delta.*Itime_Lck(t-1)) )*dts ;
    ztime_Lck(t) =  ztime_Lck(t-1) - aa*Lck(t-1)*varsigma*dts - ztime_Lck(t-1)*varsigma*dts ;
    Ntime_Lck(t) =  Stime_Lck(t)   + Itime_Lck(t) + Rtime_Lck(t);
    

    Lck(t) =  pol_interp(Stime_Lck(t),Itime_Lck(t));
    
    PiTilde1(t) = pitlde1_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde2(t) = pitlde2_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde3(t) = pitlde3_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde4(t) = pitlde4_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde5(t) = pitlde5_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde6(t) = pitlde6_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde7(t) = pitlde7_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde8(t) = pitlde8_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde9(t) = pitlde9_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde10(t) = pitlde10_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde11(t) = pitlde11_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde12(t) = pitlde12_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde13(t) = pitlde13_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde14(t) = pitlde14_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde15(t) = pitlde15_interp(Stime_Lck(t),Itime_Lck(t)) ;  
    PiTilde16(t) = pitlde16_interp(Stime_Lck(t),Itime_Lck(t)) ;     
    Alt_R0(t) =  alt_R0(Stime_Lck(t),Itime_Lck(t));
    Ihat(t) = I_hat_interp(Stime_Lck(t),Itime_Lck(t)); 
    
    
sigma_S = [-Stime_Lck(t-1).*Itime_Lck(t-1).*sigma_i, Stime_Lck(t-1).*Itime_Lck(t-1).*sigma_d, 0];
sigma_I = [ Stime_Lck(t-1).*Itime_Lck(t-1).*sigma_i, - Itime_Lck(t-1).*(1-Itime_Lck(t-1)).*sigma_d, 0];
sigma_N = [0, -Itime_Lck(t-1).*sigma_d.*Ntime_Lck(t-1), 0];
sigma_Z = [0,0,sigma_z];
sigma = [sigma_S;sigma_I;sigma_N;sigma_Z];
sigma_p = sigma';
sigma_p_sigma_inv = pinv(sigma_p*sigma);
sigma_p_sigma_inv_p = sigma_p_sigma_inv';

mu_s_bar =  Stime_Lck(t) - Stime_Lck(t-1);
mu_i_bar =  Itime_Lck(t) - Itime_Lck(t-1);
mu_z_bar =  ztime_Lck(t) - ztime_Lck(t-1);
mu_N_bar =  Ntime_Lck(t) - Ntime_Lck(t-1);

alt_beta(t-1) =    beta_vals(1).*PiTilde1(t-1)+beta_vals(2).*PiTilde2(t-1)+beta_vals(3).*PiTilde3(t-1)...
             +beta_vals(4).*PiTilde4(t-1)+beta_vals(5).*PiTilde5(t-1)+beta_vals(6).*PiTilde6(t-1)...
             +beta_vals(7).*PiTilde7(t-1)+beta_vals(8).*PiTilde8(t-1)+beta_vals(9).*PiTilde9(t-1)...
             +beta_vals(10).*PiTilde10(t-1)...
             +beta_vals(11).*PiTilde11(t-1)+beta_vals(12).*PiTilde12(t-1)+beta_vals(13).*PiTilde13(t-1)...
             +beta_vals(14).*PiTilde14(t-1)+beta_vals(15).*PiTilde15(t-1)+beta_vals(16).*PiTilde16(t-1);
alt_delta(t-1) =    delta_vals(1).*PiTilde1(t-1)+delta_vals(2).*PiTilde2(t-1)+delta_vals(3).*PiTilde3(t-1)...
             +delta_vals(4).*PiTilde4(t-1)+delta_vals(5).*PiTilde5(t-1)+delta_vals(6).*PiTilde6(t-1)...
             +delta_vals(7).*PiTilde7(t-1)+delta_vals(8).*PiTilde8(t-1)+delta_vals(9).*PiTilde9(t-1)...
             +delta_vals(10).*PiTilde10(t-1)...
             +delta_vals(11).*PiTilde11(t-1)+delta_vals(12).*PiTilde12(t-1)+delta_vals(13).*PiTilde13(t-1)...
             +delta_vals(14).*PiTilde14(t-1)+delta_vals(15).*PiTilde15(t-1)+delta_vals(16).*PiTilde16(t-1);    
         
alt_theta(t-1) =    theta_vals(1).*PiTilde1(t-1)+theta_vals(2).*PiTilde2(t-1)+theta_vals(3).*PiTilde3(t-1)...
             +theta_vals(4).*PiTilde4(t-1)+theta_vals(5).*PiTilde5(t-1)+theta_vals(6).*PiTilde6(t-1)...
             +theta_vals(7).*PiTilde7(t-1)+theta_vals(8).*PiTilde8(t-1)+theta_vals(9).*PiTilde9(t-1)...
             +theta_vals(10).*PiTilde10(t-1)...
             +theta_vals(11).*PiTilde11(t-1)+theta_vals(12).*PiTilde12(t-1)+theta_vals(13).*PiTilde13(t-1)...
             +theta_vals(14).*PiTilde14(t-1)+theta_vals(15).*PiTilde15(t-1)+theta_vals(16).*PiTilde16(t-1);
         
alt_aa(t-1) =      aa_vals(1).*PiTilde1(t-1)+aa_vals(2).*PiTilde2(t-1)+aa_vals(3).*PiTilde3(t-1)...
             +aa_vals(4).*PiTilde4(t-1)+aa_vals(5).*PiTilde5(t-1)+aa_vals(6).*PiTilde6(t-1)...
             +aa_vals(7).*PiTilde7(t-1)+aa_vals(8).*PiTilde8(t-1)+aa_vals(9).*PiTilde9(t-1)...
             +aa_vals(10).*PiTilde10(t-1)...
             +aa_vals(11).*PiTilde11(t-1)+aa_vals(12).*PiTilde12(t-1)+aa_vals(13).*PiTilde13(t-1)...
             +aa_vals(14).*PiTilde14(t-1)+aa_vals(15).*PiTilde15(t-1)+aa_vals(16).*PiTilde16(t-1);
    
mu_s_tilde =  - alt_beta(t-1)*Stime_Lck(t-1)*Itime_Lck(t-1)*( 1-alt_theta(t-1)*Lck(t-1))^2*dts - (-Itime_Lck(t-1)*(alt_delta(t-1)+5.*alt_delta(t-1).*Itime_Lck(t-1)))*Stime_Lck(t-1)*dts ;
mu_i_tilde =  - ggamma*Itime_Lck(t-1)*dts - (-Itime_Lck(t-1)*(alt_delta(t-1)+5.*alt_delta(t-1).*Itime_Lck(t-1)))*Itime_Lck(t-1)*dts + alt_beta(t-1)*Stime_Lck(t-1)*Itime_Lck(t-1)*( 1-alt_theta(t-1)*Lck(t-1))^2*dts ;    
mu_r_tilde =    Itime_Lck(t-1)*(ggamma - (alt_delta(t-1)+5.*alt_delta(t-1).*Itime_Lck(t-1)) )*dts ;
mu_z_tilde =  - alt_aa(t-1)*Lck(t-1)*varsigma*dts - ztime_Lck(t-1)*varsigma*dts ;
mu_N_tilde =    mu_s_tilde + mu_i_tilde + mu_r_tilde;  


delta_mu = [mu_s_tilde-mu_s_bar;mu_i_tilde-mu_i_bar;mu_N_tilde-mu_N_bar;mu_z_tilde-mu_z_bar];
h_val = sigma_p_sigma_inv*sigma_p*delta_mu;
h_p_h = (h_val')*h_val;
prob2(t-1) = 0.5.*exp(-12.*h_p_h./8);   

    
end

Frac_Pop_Lckdwn= (Stime_Lck + Itime_Lck + Rtime_Lck).*Lck';

toc;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

close all
clear all
% clc

set(0,'defaulttextinterpreter','tex')
set(0, 'defaultAxesTickLabelInterpreter','tex'); 
set(0, 'defaultLegendInterpreter','tex');

directory = pwd;

for cases = 1:4

direct = pwd;
directory = [direct];    


    if cases == 1
A=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_11','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','time','Nyears','dts',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');

B=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_11','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');

    elseif cases == 2
A=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_21','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','time','Nyears','dts',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');

B=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_22','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');

    elseif cases == 3
A=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_31','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','time','Nyears','dts',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');

B=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_32','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');

    elseif cases == 4
A=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_41','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR','time','Nyears','dts',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa','N_S','k','Alt_R0');

B=load([directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_42','_sims'],'Stime_Lck','Itime_Lck','Rtime_Lck','Dtime_Lck','Lck','R0','CFR',...
    'PiTilde1','PiTilde2','PiTilde3','PiTilde4','PiTilde5','PiTilde6','PiTilde7','PiTilde8','PiTilde9','PiTilde10','PiTilde11','PiTilde12','PiTilde13','PiTilde14','PiTilde15','PiTilde16', ...
    'Ihat','beta_vals','delta_vals','aa_vals','theta_vals','ggamma','weights','prob2','alt_beta','alt_delta','alt_theta','alt_aa');
    end  

color1 = [0    0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];
color4 = [0.4940    0.1840    0.5560];
color5 = [0.4660    0.6740    0.1880];
color6 = [0.3010    0.7450    0.9330];
color7 = [0.6350    0.0780    0.1840];

scale = 2;

time = A.time;
time_alt = A.time./scale;
Nyears = A.Nyears; 
Nyears_alt = A.Nyears./scale;
dt = A.dts;
time_vec = linspace(0,Nyears_alt,time_alt);
time_vec2 = linspace(0,52*Nyears_alt,time_alt);
Alt_end = round(length(A.Stime_Lck)./scale);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure;
hold on;
plot(time_vec2,(A.Itime_Lck(1:Alt_end)),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
plot(time_vec2,(B.Itime_Lck(1:Alt_end)),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('\textbf{Infected}', 'Interpreter', 'latex')
ylabel('$$i_t$$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
axis([0 52*Nyears_alt 0 0.25])
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
print('covid_RR_I_SA_Correct','-dpng')
elseif cases == 2
print('covid_RR_I_SA_Over','-dpng')    
elseif cases == 3
print('covid_RR_I_SA_Under','-dpng')    
elseif cases == 4
print('covid_RR_I_SA_True','-dpng')    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure;
hold on;
plot(time_vec2,(A.Dtime_Lck(1:Alt_end)),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
plot(time_vec2,(B.Dtime_Lck(1:Alt_end)),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('\textbf{Dead}', 'Interpreter', 'latex')
ylabel('$$d_t$$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northwest');
axis([0 52*Nyears_alt 0 0.025])
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
print('covid_RR_D_SA_Correct','-dpng')
elseif cases == 2
print('covid_RR_D_SA_Over','-dpng')    
elseif cases == 3
print('covid_RR_D_SA_Under','-dpng')    
elseif cases == 4
print('covid_RR_D_SA_True','-dpng')   
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure;
hold on;
plot(time_vec2,(A.Lck(1:Alt_end)),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
plot(time_vec2,(B.Lck(1:Alt_end)),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('\textbf{Quarantine}', 'Interpreter', 'latex')
ylabel('$$q_t$$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
axis([0 52*Nyears_alt 0 0.75])
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])

if cases == 1
print('covid_RR_Q_SA_Correct','-dpng')
elseif cases == 2
print('covid_RR_Q_SA_Over','-dpng')    
elseif cases == 3
print('covid_RR_Q_SA_Under','-dpng')   
elseif cases == 4
print('covid_RR_Q_SA_True','-dpng')    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if cases == 2 || cases == 4
pi_tilde_iters = [B.PiTilde1',B.PiTilde2',B.PiTilde3',B.PiTilde5',B.PiTilde6',B.PiTilde7',B.PiTilde8', ...
                  B.PiTilde9',B.PiTilde10',B.PiTilde11',B.PiTilde13',B.PiTilde14',B.PiTilde15',B.PiTilde16']';

sorted = sort(pi_tilde_iters, 1);
middle_line = mean(sorted, 1);
bottom_line = sorted(end, :);
top_line = sorted(1,:);
xfill = [time_vec2, fliplr(time_vec2)];
inBetween = [bottom_line(1:Alt_end), fliplr(top_line(1:Alt_end))];
elseif cases == 1 || cases == 3
    pi_tilde_iters = [B.PiTilde1',B.PiTilde2',B.PiTilde3',B.PiTilde5',B.PiTilde6',B.PiTilde7',B.PiTilde8', ...
                  B.PiTilde9',B.PiTilde10',B.PiTilde11',B.PiTilde13',B.PiTilde14',B.PiTilde15',B.PiTilde16']';

sorted = sort(pi_tilde_iters, 1);
middle_line = mean(sorted, 1);
bottom_line = sorted(end, :);
top_line = sorted(1,:);
xfill = [time_vec2, fliplr(time_vec2)];
inBetween = [bottom_line(1:Alt_end), fliplr(top_line(1:Alt_end))];
end

figure
hold on
plot(time_vec2,(B.PiTilde12(1:Alt_end)),'LineWidth',3,'Color',color3,'LineStyle','--', 'DisplayName', sprintf('$$\\mathcal{R}_0 = %0.1f, CFR = %0.3f, \\zeta = %0.2f, \\alpha = %0.1f$$', ...
    B.beta_vals(12)./B.ggamma,B.delta_vals(12)./B.ggamma,B.theta_vals(12),B.aa_vals(12)))
plot(time_vec2,(B.PiTilde4(1:Alt_end)),'LineWidth',3,'Color',color2,'LineStyle','--', 'DisplayName', sprintf('$$\\mathcal{R}_0 = %0.1f, CFR = %0.3f, \\zeta = %0.2f, \\alpha = %0.1f$$', ...
    B.beta_vals(4)./B.ggamma,B.delta_vals(4)./B.ggamma,B.theta_vals(4),B.aa_vals(4)))
if cases == 1 || cases == 3
plot(time_vec2,(B.PiTilde5(1:Alt_end)),'LineWidth',3,'Color',color5,'LineStyle','--', 'DisplayName', sprintf('$$\\mathcal{R}_0 = %0.1f, CFR = %0.3f, \\zeta = %0.2f, \\alpha = %0.1f$$', ...
    B.beta_vals(5)./B.ggamma,B.delta_vals(5)./B.ggamma,B.theta_vals(5),B.aa_vals(5)))
end
fill(xfill, inBetween, color2, 'FaceAlpha', .25, 'DisplayName', 'Remaining Weights');
hl = legend('show','Interpreter', 'latex', 'Location', 'northeast');
title('\textbf{Distorted Probability Weights}', 'Interpreter', 'latex')
ylabel('$$\tilde{\pi}_t$$','Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
axis([0 52*Nyears_alt 0 1.0])
print('covid_RR_pis_SA_Correct','-dpng')
elseif cases == 2
axis([0 52*Nyears_alt 0 1.0])
    set(hl, 'Interpreter', 'latex', 'Location', 'southeast');
print('covid_RR_pis_SA_Over','-dpng')    
elseif cases == 3
axis([0 52*Nyears_alt 0 1.0])
print('covid_RR_pis_SA_Under','-dpng')   
elseif cases == 4
axis([0 52*Nyears_alt 0 1.0])
print('covid_RR_pis_SA_True','-dpng')    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        
         
alt_beta =    B.beta_vals(1).*B.PiTilde1+B.beta_vals(2).*B.PiTilde2+B.beta_vals(3).*B.PiTilde3...
             +B.beta_vals(4).*B.PiTilde4+B.beta_vals(5).*B.PiTilde5+B.beta_vals(6).*B.PiTilde6...
             +B.beta_vals(7).*B.PiTilde7+B.beta_vals(8).*B.PiTilde8+B.beta_vals(9).*B.PiTilde9...
             +B.beta_vals(10).*B.PiTilde10...
             +B.beta_vals(11).*B.PiTilde11+B.beta_vals(12).*B.PiTilde12+B.beta_vals(13).*B.PiTilde13...
             +B.beta_vals(14).*B.PiTilde14+B.beta_vals(15).*B.PiTilde15+B.beta_vals(16).*B.PiTilde16;
         
base_beta =   A.beta_vals(1).*A.weights(1)+A.beta_vals(2).*A.weights(2)+A.beta_vals(3).*A.weights(3)...
             +A.beta_vals(4).*A.weights(4)+A.beta_vals(5).*A.weights(5)+A.beta_vals(6).*A.weights(6)...
             +A.beta_vals(7).*A.weights(7)+A.beta_vals(8).*A.weights(8)+A.beta_vals(9).*A.weights(9)...
             +A.beta_vals(10).*A.weights(10)...
             +A.beta_vals(11).*A.weights(11)+A.beta_vals(12).*A.weights(12)+A.beta_vals(13).*A.weights(13)...
             +A.beta_vals(14).*A.weights(14)+A.beta_vals(15).*A.weights(15)+A.beta_vals(16).*A.weights(16);         

alt_R0 = alt_beta./B.ggamma;
base_R0 = base_beta./B.ggamma;
alt2_R0 = B.alt_beta./B.ggamma;

figure
plot(time_vec2,base_R0.*ones(size(A.Itime_Lck(1:Alt_end))),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
hold on
plot(time_vec2,alt2_R0(1:Alt_end),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('$\mathcal{R}_0$ \textbf{Comparison}', 'Interpreter', 'latex')
ylabel('$\mathcal{R}_0$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
axis([0 52*Nyears_alt min(B.beta_vals)./B.ggamma max(B.beta_vals)./B.ggamma]) 
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_R0_SA_Correct','-dpng')
elseif cases == 2
    set(hl, 'Interpreter', 'latex', 'Location', 'southeast');
print('covid_RR_R0_SA_Over','-dpng')    
elseif cases == 3
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_R0_SA_Under','-dpng')    
elseif cases == 4
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_R0_SA_True','-dpng')    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        
         
base_delta =  A.delta_vals(1).*A.weights(1)+A.delta_vals(2).*A.weights(2)+A.delta_vals(3).*A.weights(3)...
             +A.delta_vals(4).*A.weights(4)+A.delta_vals(5).*A.weights(5)+A.delta_vals(6).*A.weights(6)...
             +A.delta_vals(7).*A.weights(7)+A.delta_vals(8).*A.weights(8)+A.delta_vals(9).*A.weights(9) ...
             +A.delta_vals(10).*A.weights(10)...
             +A.delta_vals(11).*A.weights(11)+A.delta_vals(12).*A.weights(12)+A.delta_vals(13).*A.weights(13)...
             +A.delta_vals(14).*A.weights(14)+A.delta_vals(15).*A.weights(15)+A.delta_vals(16).*A.weights(16)...
         +5.*(A.delta_vals(1).*A.weights(1)+A.delta_vals(2).*A.weights(2)+A.delta_vals(3).*A.weights(3)...
             +A.delta_vals(4).*A.weights(4)+A.delta_vals(5).*A.weights(5)+A.delta_vals(6).*A.weights(6)...
             +A.delta_vals(7).*A.weights(7)+A.delta_vals(8).*A.weights(8)+A.delta_vals(9).*A.weights(9)...
             +A.delta_vals(10).*A.weights(10)...
             +A.delta_vals(11).*A.weights(11)+A.delta_vals(12).*A.weights(12)+A.delta_vals(13).*A.weights(13)...
             +A.delta_vals(14).*A.weights(14)+A.delta_vals(15).*A.weights(15)+A.delta_vals(16).*A.weights(16)).*A.Itime_Lck ;

base_delta2 =  B.delta_vals(1).*B.weights(1)+B.delta_vals(2).*B.weights(2)+B.delta_vals(3).*B.weights(3)...
             +B.delta_vals(4).*B.weights(4)+B.delta_vals(5).*B.weights(5)+B.delta_vals(6).*B.weights(6)...
             +B.delta_vals(7).*B.weights(7)+B.delta_vals(8).*B.weights(8)+B.delta_vals(9).*B.weights(9) ...
             +B.delta_vals(10).*B.weights(10)...
             +B.delta_vals(11).*B.weights(11)+B.delta_vals(12).*B.weights(12)+B.delta_vals(13).*B.weights(13)...
             +B.delta_vals(14).*B.weights(14)+B.delta_vals(15).*B.weights(15)+B.delta_vals(16).*B.weights(16)...
         +5.*(B.delta_vals(1).*B.weights(1)+B.delta_vals(2).*B.weights(2)+B.delta_vals(3).*B.weights(3)...
             +B.delta_vals(4).*B.weights(4)+B.delta_vals(5).*B.weights(5)+B.delta_vals(6).*B.weights(6)...
             +B.delta_vals(7).*B.weights(7)+B.delta_vals(8).*B.weights(8)+B.delta_vals(9).*B.weights(9)...
             +B.delta_vals(10).*B.weights(10)...
             +B.delta_vals(11).*B.weights(11)+B.delta_vals(12).*B.weights(12)+B.delta_vals(13).*B.weights(13)...
             +B.delta_vals(14).*B.weights(14)+B.delta_vals(15).*B.weights(15)+B.delta_vals(16).*B.weights(16)).*B.Itime_Lck ;

alt_delta_temp =   B.delta_vals(1).*B.PiTilde1+B.delta_vals(2).*B.PiTilde2+B.delta_vals(3).*B.PiTilde3...
             +B.delta_vals(4).*B.PiTilde4+B.delta_vals(5).*B.PiTilde5+B.delta_vals(6).*B.PiTilde6...
             +B.delta_vals(7).*B.PiTilde7+B.delta_vals(8).*B.PiTilde8+B.delta_vals(9).*B.PiTilde9 ...
             +B.delta_vals(10).*B.PiTilde10...
             +B.delta_vals(11).*B.PiTilde11+B.delta_vals(12).*B.PiTilde12+B.delta_vals(13).*B.PiTilde13...
             +B.delta_vals(14).*B.PiTilde14+B.delta_vals(15).*B.PiTilde15+B.delta_vals(16).*B.PiTilde16;
         
alt_delta =   alt_delta_temp + 5.*alt_delta_temp.*B.Itime_Lck;
alt2_delta =   B.alt_delta + 5.*B.alt_delta.*B.Itime_Lck;
     
alt_cfr = alt_delta./B.ggamma;
alt2_cfr = alt2_delta./B.ggamma;

base_cfr = base_delta./B.ggamma;
base_cfr2 = base_delta2./B.ggamma;

figure
plot(time_vec2,base_cfr(1:Alt_end),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
hold on
plot(time_vec2,alt_cfr(1:Alt_end),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('$CFR$ \textbf{Comparison}', 'Interpreter', 'latex')
ylabel('$CFR$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
ax = gca;
ax.YAxis.Exponent = 0;
axis([0 52*Nyears_alt 0.005 0.03])  
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_CFR_SA_Correct','-dpng')
elseif cases == 2
    set(hl, 'Interpreter', 'latex', 'Location', 'southeast');
print('covid_RR_CFR_SA_Over','-dpng')    
elseif cases == 3
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_CFR_SA_Under','-dpng')    
elseif cases == 4
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_CFR_SA_True','-dpng')  
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

alt_theta =    B.theta_vals(1).*B.PiTilde1+B.theta_vals(2).*B.PiTilde2+B.theta_vals(3).*B.PiTilde3...
             +B.theta_vals(4).*B.PiTilde4+B.theta_vals(5).*B.PiTilde5+B.theta_vals(6).*B.PiTilde6...
             +B.theta_vals(7).*B.PiTilde7+B.theta_vals(8).*B.PiTilde8+B.theta_vals(9).*B.PiTilde9...
             +B.theta_vals(10).*B.PiTilde10...
             +B.theta_vals(11).*B.PiTilde11+B.theta_vals(12).*B.PiTilde12+B.theta_vals(13).*B.PiTilde13...
             +B.theta_vals(14).*B.PiTilde14+B.theta_vals(15).*B.PiTilde15+B.theta_vals(16).*B.PiTilde16;
         
base_theta =   A.theta_vals(1).*A.weights(1)+A.theta_vals(2).*A.weights(2)+A.theta_vals(3).*A.weights(3)...
             +A.theta_vals(4).*A.weights(4)+A.theta_vals(5).*A.weights(5)+A.theta_vals(6).*A.weights(6)...
             +A.theta_vals(7).*A.weights(7)+A.theta_vals(8).*A.weights(8)+A.theta_vals(9).*A.weights(9)...
             +A.theta_vals(10).*A.weights(10)...
             +A.theta_vals(11).*A.weights(11)+A.theta_vals(12).*A.weights(12)+A.theta_vals(13).*A.weights(13)...
             +A.theta_vals(14).*A.weights(14)+A.theta_vals(15).*A.weights(15)+A.theta_vals(16).*A.weights(16);         


figure
plot(time_vec2,base_theta.*ones(size(A.Itime_Lck(1:Alt_end))),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
hold on
plot(time_vec2,alt_theta(1:Alt_end),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('$$\zeta$$ \textbf{Comparison}', 'Interpreter', 'latex')
ylabel('$$\zeta$$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
axis([0 52*Nyears_alt min(B.theta_vals) max(B.theta_vals)]) 
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Zeta_SA_Correct','-dpng')
elseif cases == 2
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Zeta_SA_Over','-dpng')    
elseif cases == 3
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Zeta_SA_Under','-dpng')  
elseif cases == 4
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Zeta_SA_True','-dpng')  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



alt_aa =      B.aa_vals(1).*B.PiTilde1+B.aa_vals(2).*B.PiTilde2+B.aa_vals(3).*B.PiTilde3...
             +B.aa_vals(4).*B.PiTilde4+B.aa_vals(5).*B.PiTilde5+B.aa_vals(6).*B.PiTilde6...
             +B.aa_vals(7).*B.PiTilde7+B.aa_vals(8).*B.PiTilde8+B.aa_vals(9).*B.PiTilde9...
             +B.aa_vals(10).*B.PiTilde10...
             +B.aa_vals(11).*B.PiTilde11+B.aa_vals(12).*B.PiTilde12+B.aa_vals(13).*B.PiTilde13...
             +B.aa_vals(14).*B.PiTilde14+B.aa_vals(15).*B.PiTilde15+B.aa_vals(16).*B.PiTilde16;
         
base_aa =     A.aa_vals(1).*A.weights(1)+A.aa_vals(2).*A.weights(2)+A.aa_vals(3).*A.weights(3)...
             +A.aa_vals(4).*A.weights(4)+A.aa_vals(5).*A.weights(5)+A.aa_vals(6).*A.weights(6)...
             +A.aa_vals(7).*A.weights(7)+A.aa_vals(8).*A.weights(8)+A.aa_vals(9).*A.weights(9)...
             +A.aa_vals(10).*A.weights(10)...
             +A.aa_vals(11).*A.weights(11)+A.aa_vals(12).*A.weights(12)+A.aa_vals(13).*A.weights(13)...
             +A.aa_vals(14).*A.weights(14)+A.aa_vals(15).*A.weights(15)+A.aa_vals(16).*A.weights(16);         


figure
plot(time_vec2,base_aa.*ones(size(A.Itime_Lck(1:Alt_end))),'-r','LineWidth',3, 'DisplayName', 'Uncertainty Neutral');
hold on
plot(time_vec2,alt_aa(1:Alt_end),'-b','LineWidth',3, 'DisplayName', 'Uncertainty Averse');
title('$$\alpha$$ \textbf{Comparison}', 'Interpreter', 'latex')
ylabel('$$\alpha$$', 'Interpreter', 'latex')
xlabel('\textbf{Weeks}', 'Interpreter', 'latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold')
hl = legend('show');
set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
axis([0 52*Nyears_alt min(B.aa_vals) max(B.aa_vals)]) 
xticks([0 10 20 30 40 50])
set(gcf,'Position',[50 50 650 450])
if cases == 1
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Alpha_SA_Correct','-dpng')
elseif cases == 2
    set(hl, 'Interpreter', 'latex', 'Location', 'southeast');
print('covid_RR_Alpha_SA_Over','-dpng')    
elseif cases == 3
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Alpha_SA_Under','-dpng')    
elseif cases == 4
    set(hl, 'Interpreter', 'latex', 'Location', 'northeast');
print('covid_RR_Alpha_SA_True','-dpng')    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

