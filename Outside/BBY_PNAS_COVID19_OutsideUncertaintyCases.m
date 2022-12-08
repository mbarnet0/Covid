% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% False Transient Method For Simple Climate Example
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all
close all
clc

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Filename set-up
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

directory = pwd;

% Setting for saving data and figures
savedata   = 0 ; 
initialize = 1 ; % == 0 if already have VF

filename_load = [directory, '/Results_BBY_PNAS_COVID_InsideUncertainty_41']; % equal-weighted prior (4), uncertainty neutral (1)
filename_save = [directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1111']; %beta case 1, delta case 1, zeta case 1, alpha case 1

smartguess = 0;
if smartguess == 1
SmartGuess = load(filename_load, 'V','S_mat','I_mat','beta','varphi','dV_dI','dV_dS');
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Parameters for the model
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% RATES - ALLL MEASURED AS PER-YEAR RATES 
r  = 0.03./12;    % annual interest rate
nu = 1./1.0./12;       % annual vaccine arrival rate

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Vectors for SA
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
ggamma  = 30/18;    %  avg 18 DAYS to get out of infection state (either recovered or dead)

% Uncomment for the case you want to solve
% Birth Number of Pandemic
% R_0 = 1.5;
% R_0 = 4.5;

% Case Fatality Rate of Pandemic
% CFR = 0.005;
% CFR = 0.02;

% Mitigation Effectiveness
% theta = 0.35;
% theta = 0.65;

% Mitigation Costs Parameter
% aa = 0.0;
% aa = 0.8;

beta = ggamma.*R_0;
delta = ggamma.*CFR;

upsilon = 5*delta;

% Welfare cost parameters
w     = 1   ;  % ANNUAL GDP per capita 
vsl   = 1  ;  % value of  stistical life as multiple of ANNUAL earnings (since w=1 ANNUAL OUTPUT)
omega = 0.0;

varphi = 1.0;
varsigma = 1;

xi = 100000;

% Volatility
sigma_s = 0.0;
sigma_i = 0.075;
sigma_r = 0.0;
sigma_d = 0.030;
sigma_z = 0.016./sqrt(12);

% policy vector
phi        =  0.4;   % productivity for infected
L_max      =  0.99;  % bench = 0.70;

% Initial conditions for graphs
I0 = 0.005; % S(1) = 1- I(1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Grids for state variables
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

N_S = 650;
k = 5;
N_I = (N_S-1)*k + 1 ;

Delta_S = 1/(N_S-1);
Delta_I = 1/(N_I-1);

hs = Delta_S;
hi = Delta_I;

ns = N_S;
ni = N_I;

S        = linspace(0,1,N_S); 
I        = linspace(0,1,N_I); 

[S_mat,I_mat] = ndgrid(S,I); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Kushner-Dupuis for comparison
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if initialize==1

    V_initial = r.*log(1-(1-phi)*I_mat)./(r+nu) -r.*log(1-S_mat.*I_mat)./(r+nu); % PV in daily output units w=1
%     V_initial = r.*log(1-(1-phi)*I_mat)./(r+nu); % PV in daily output units w=1

    if smartguess == 1
    V_initial_temp = SmartGuess.V;
    SmartGuess.check_val = SmartGuess.S_mat+SmartGuess.I_mat;
    SmartGuess.checkval = SmartGuess.check_val.*(SmartGuess.check_val<=1);
    SmartGuess.inds2=find(SmartGuess.checkval==0);
    SmartGuess.inds2a = SmartGuess.inds2(2:end);
    V_initial_temp(SmartGuess.inds2a)= 0;
    dV_dI_initial_temp(SmartGuess.inds2a)= 0;
    dV_dS_initial_temp(SmartGuess.inds2a)= 0;
    

      V_initial = interpn(SmartGuess.S_mat,SmartGuess.I_mat,V_initial_temp,S_mat,I_mat);
    end
    
    V = V_initial;
    V_guess = V; 
    V_next = V_guess;
    
 
   policy=zeros(ns,ni);  policy_indx=policy; policy_L=policy_indx; 
   iteration=0;
end




Diff_Der_I_m_S = zeros(ns,ni);

misses = zeros(ns,ni); 

check_val = S_mat+I_mat;
checkval = check_val.*(check_val<=1);
inds2=find(checkval==0);
inds2a = inds2(2:end);
Diff_Der_I_m_S(inds2a) = NaN;
V_initial(inds2a)= NaN;
V(inds2a) = NaN;
V_guess(inds2a) = NaN;
V_next(inds2a)=NaN;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Probabilites
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

delta_t = 1./((r+nu) ...
    + (beta.*S_mat.*I_mat+S_mat.*I_mat.*(delta+upsilon.*I_mat)).*(1./hs) + (beta.*S_mat.*I_mat+(ggamma-(delta+upsilon.*I_mat).*I_mat).*I_mat).*(1./hi) ...
    + ((S_mat.*I_mat).^2.*(sigma_d-sigma_i).^2).*(1./(hs.^2)) ...
    + ((sigma_i.*I_mat.*S_mat-(1-I_mat).*I_mat.*sigma_d).^2).*(1./(hi.^2)));

max_dt = min(min(min(delta_t)));

dt = 0.95 * max_dt;


L_optimal = zeros(ns,ni);
L_old     = L_optimal ; 
L_check     = L_optimal ; 
check_foc = zeros(ns,ni);
Convex_HBJ   = zeros(ns,ni); 

check_val = S_mat+I_mat;
checkval = check_val.*(check_val<=1);
inds=find(checkval~=0);
inds2=find(checkval==0);
inds2a = inds2(2:end);
V(inds2a)=NaN;
V_next(inds2a)=NaN;
L_check(inds2a) = NaN;
check_val(inds2a) = NaN;
delta_t(inds2a) = NaN;
Diff_Der_I_m_S(inds2a) = NaN;

% stop
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%__________________________________________________________________________
%
% Convergence parameters
%__________________________________________________________________________
% convergence tolerance level
tol = 1e-6;
tol2 = 1e-12;

diffmax=10; diffmax_policy= diffmax; 
%Value function iteration

diff_v = 0; iters = 0; iter_check = 0;


countchecknan = max(sum(sum(isnan(V_next))));
[val1,val2] = size(V_next);
countcheck = max(val1*val2-countchecknan);

dV_dS = Diff_Der_I_m_S;
dV_dI = Diff_Der_I_m_S;

beta_theta_theta_hat = beta.*(theta.^2).*ones(size(S_mat));


while  diffmax>tol2 || (diffmax_policy > tol) 

%     tic;

    iteration=iteration+1   ;
    

for s= 1:ns-1    

     for i= 2 : ni-1-(s-1)*k
         iter_check = iter_check+1;

  %% 
        if s == 1 
                                    
            L_optimal(s,i) = 0;
            L_val = L_optimal(s,i) ; 
            dV_dI(s,i) = (1./hi).*(V(s,i)-V(s,i-1));

                V_next(s,i)  =  ...
                    r.*log(1-(1-phi)*I(i)).*dt - 0.5.*(I(i).*sigma_d).^2.*dt ...
                    -(delta+5.*delta.*I(i)).*I(i).*dt ...
                    + ( (1-dt.*(r+nu)) - ( (ggamma-(delta+5.*delta.*I(i)).*I(i)).*I(i).*(dt/hi) + ((-sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2)) ) ).*V(s,i) ...
                    + ( 0.5.*((-sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2)) ).* V(s,i+1)   ...
                    + ( (ggamma-(delta+5.*delta.*I(i)).*I(i)).*I(i).*(dt./hi) + 0.5.*((-sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2))  ).*V(s,i-1);                 
                
                diff_v = diff_v + abs(V_next(s,i)-V(s,i)); if isnan(diff_v) == 1; stop; end; 
                iters = iters+1;
                
  %%             
        elseif s > 1 && s < ns

              Diff_Der_I_m_S(s,i)  = 1/(hi) * (V(s,i+1) - V(s,i)) - 1/(hs) * (V(s,i) - V(s-1,i)) ;  
              dV_dS(s,i) = 1/(hs) * (V(s,i) - V(s-1,i));
              dV_dI(s,i) = 1/(hi) * (V(s,i+1) - V(s,i));

              if isnan(Diff_Der_I_m_S(s,i)) == 1; stop; end;                      
              if S(s+1)+I(i+1) <= 1
                  Partial_Diff_I_S(s,i) = (V(s+1,i+1)-V(s+1,i-1)-V(s-1,i+1)+V(s-1,i-1))./(4*hs*hi);
              else
                  Partial_Diff_I_S(s,i) = (V(s,i+1)-V(s,i-1)-V(s-1,i+1)+V(s-1,i-1))./(2*hs*hi);
              end
                  
            if  (r.*varphi).*(1-L_old(s,i)).^(-2)./ ...
                (2.*beta_theta_theta_hat(s,i).*S(s).*I(i).*(Diff_Der_I_m_S(s,i) )) < 1            
            Aa = 1;
            Bb = -(beta.*theta+beta_theta_theta_hat(s,i))./beta_theta_theta_hat(s,i) - (r.*varsigma./(r+varsigma)).*aa./(2.*beta_theta_theta_hat(s,i).*S(s).*I(i).*Diff_Der_I_m_S(s,i));
            Cc = (r.*varphi + (r.*varsigma./(r+varsigma)).*aa)./(2.*beta_theta_theta_hat(s,i).*S(s).*I(i).*Diff_Der_I_m_S(s,i)) + beta.*theta./beta_theta_theta_hat(s,i);
            Lfoc	= (-Bb -sqrt(Bb.^2 -4.*Aa.*Cc)) ./(2.*Aa);             
            L_optimal(s,i) = min(L_max, max(Lfoc, 0));
            L_val = L_optimal(s,i) ;        
            else
            L_optimal(s,i) = 0;
            L_val = L_optimal(s,i) ;  
            end
            
            

    i_val = max(i-k,1);

                V_next(s,i)  =  ...
                    r.*log(1-(1-phi)*I(i)).*dt + (r.*varphi)*log(1-L_val).*dt - 0.5.*(I(i).*sigma_d).^2.*dt ...
                    - Partial_Diff_I_S(s,i).*(sigma_i.^2.*I(i).^2.*S(s).^2 + sigma_d.^2.*I(i).^2.*(1-I(i)).*S(s) ).*dt ...  
                    -(delta+5.*delta.*I(i)).*I(i).*dt - (r.*varsigma./(r+varsigma)).*aa.*L_val.*dt ...
                    + (1-dt.*(r+nu)).*V(s,i) ...
                    +  I(i).*S(s).*(beta - 2.*beta.*theta.*L_val + beta_theta_theta_hat(s,i).*L_val.^2).*(dt./hs).* V(s-1,i)  ...
                    -  I(i).*S(s).*(beta - 2.*beta.*theta.*L_val + beta_theta_theta_hat(s,i).*L_val.^2).*(dt./hs).*V(s,i) ...  
                    + S(s).*I(i).*(delta+5.*delta.*I(i)).*(dt./hs).* V(s+1,i_val) ...                    
                    - S(s).*I(i).*(delta+5.*delta.*I(i)).*(dt./hs).*V(s,i_val) ...
                    + I(i).*S(s).*(beta - 2.*beta.*theta.*L_val + beta_theta_theta_hat(s,i).*L_val.^2).*(dt./hi).* V(s,i+1)   ...                   
                    - I(i).*S(s).*(beta - 2.*beta.*theta.*L_val + beta_theta_theta_hat(s,i).*L_val.^2).*(dt/hi).*V(s,i) ...
                    + (ggamma-(delta+5.*delta.*I(i)).*I(i)).*I(i).*(dt./hi).*V(s,i-1)...                                        
                    - (ggamma-(delta+5.*delta.*I(i)).*I(i)).*I(i).*(dt/hi).*V(s,i) ...
                    - ((sigma_d.^2+sigma_i.^2).*I(i).^2.*S(s).^2).*(dt./(hs.^2)).*V(s,i_val) ...                    
                    + ( 0.5.*((sigma_d.^2+sigma_i.^2).*I(i).^2.*S(s).^2).*(dt./(hs.^2)) ).* V(s+1,i_val) ...
                    + ( 0.5.*((sigma_d.^2+sigma_i.^2).*I(i).^2.*S(s).^2).*(dt./(hs.^2)) ).* V(s-1,i_val)  ...                                        
                    - ((sigma_i.*I(i).*S(s)).^2 + (sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2)).*V(s,i) ...                    
                    + ( 0.5.*((sigma_i.*I(i).*S(s)).^2 + (sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2))  ).* V(s,i+1)   ...
                    + ( 0.5.*((sigma_i.*I(i).*S(s)).^2 + (sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2))  ).*V(s,i-1); 
                
                diff_v = diff_v + abs(V_next(s,i)-V(s,i)); if isnan(diff_v) == 1; stop; end; 
                iters = iters+1;
                                if isnan(L_optimal(s,i)) == 1; stop; end;      

        end
                     
                 
        
        
        
     end
         
        i = ni-(s-1)*k;
        if s==1 
                 
            L_optimal(s,i) = 0;
            L_val = L_optimal(s,i) ; 
            dV_dI(s,i) = (1./hi).*(V(s,i)-V(s,i-1));
  
                V_next(s,i)  =  ...
                   r.*log(1-(1-phi)).*dt - 0.5.*(sigma_d).^2.*dt...
                    -(delta+5.*delta.*I(i)).*I(i).*dt ...                
                    + ( (1-dt.*(r+nu)) - ((ggamma-(delta+5.*delta)).*(dt/hi)  ) ).*V(s,i) ...
                    + ((ggamma-(delta+5.*delta)).*(dt./hi) ).*V(s,i-1);   
                
  
                diff_v = diff_v + abs(V_next(s,i)-V(s,i)); if isnan(diff_v) == 1; stop; end; 
                iters = iters+1;
                
        elseif s > 1 && s < ns
            

             Diff_Der_I_m_S(s,i)  = 1/(hs) * (V(s-1,i+k) - V(s,i) ) ;   
             Partial_Diff_I_S(s,i) = (V(s,i)-V(s,i-1)-V(s-1,i)+V(s-1,i-1))./(hs*hi);
             dV_dS(s,i) = 1/(hs) * (V(s,i) - V(s-1,i));
             dV_dI(s,i) = 1/(hi) * (V(s-1,i+1) - V(s-1,i));  
             if isnan(Diff_Der_I_m_S(s,i)) == 1; stop; end; 

           
      
            if  (r.*varphi).*(1-L_old(s,i)).^(-2)./ ...
                (2.*beta_theta_theta_hat(s,i).*S(s).*I(i).*(Diff_Der_I_m_S(s,i) )) < 1            
            Aa = 1;
            Bb = -(beta.*theta+beta_theta_theta_hat(s,i))./beta_theta_theta_hat(s,i) - (r.*varsigma./(r+varsigma)).*aa./(2.*beta_theta_theta_hat(s,i).*S(s).*I(i).*Diff_Der_I_m_S(s,i));
            Cc = (r.*varphi + (r.*varsigma./(r+varsigma)).*aa)./(2.*beta_theta_theta_hat(s,i).*S(s).*I(i).*Diff_Der_I_m_S(s,i)) + beta.*theta./beta_theta_theta_hat(s,i);
            Lfoc	= (-Bb -sqrt(Bb.^2 -4.*Aa.*Cc)) ./(2.*Aa);             
            L_optimal(s,i) = min(L_max, max(Lfoc, 0));
            L_val = L_optimal(s,i) ;        
            else
            L_optimal(s,i) = 0;
            L_val = L_optimal(s,i) ;  
            end

                V_next(s,i)  =  ...
                    r.*log(1-(1-phi)*I(i)).*dt + (r.*varphi)*log(1-L_val).*dt - 0.5.*(I(i).*sigma_d).^2.*dt ...
                    - Partial_Diff_I_S(s,i).*(sigma_i.^2.*I(i).^2.*S(s).^2 + sigma_d.^2.*I(i).^2.*(1-I(i)).*S(s) ).*dt ...  
                    -(delta+5.*delta.*I(i)).*I(i).*dt - (r.*varsigma./(r+varsigma)).*aa.*L_val.*dt ...
                    + (1-dt.*(r+nu)).*V(s,i) ...
                    +  I(i).*S(s).*(beta - 2.*beta.*theta.*L_val + beta_theta_theta_hat(s,i).*L_val.^2).*(dt./hs).* V(s-1,i+k)   ...   ...
                    -  I(i).*S(s).*(beta - 2.*beta.*theta.*L_val + beta_theta_theta_hat(s,i).*L_val.^2).*(dt./hs).*V(s,i) ...  
                    - S(s).*I(i).*(delta+5.*delta.*I(i)).*(dt./hs).* V(s+1,i-k) ...                    
                    + S(s).*I(i).*(delta+5.*delta.*I(i)).*(dt./hs).*V(s,i-k) ...
                    + (ggamma-(delta+5.*delta.*I(i)).*I(i)).*I(i).*(dt./hi).*V(s,i-1)...                                        
                    - (ggamma-(delta+5.*delta.*I(i)).*I(i)).*I(i).*(dt/hi).*V(s,i) ...
                    - ((sigma_d.^2+sigma_i.^2).*I(i).^2.*S(s).^2).*(dt./(hs.^2)).*V(s,i-k) ...                    
                    + ( 0.5.*((sigma_d.^2+sigma_i.^2).*I(i).^2.*S(s).^2).*(dt./(hs.^2)) ).* V(s+1,i-k) ...
                    + ( 0.5.*((sigma_d.^2+sigma_i.^2).*I(i).^2.*S(s).^2).*(dt./(hs.^2)) ).* V(s-1,i-k)  ...                                        
                    - ((sigma_i.*I(i).*S(s)).^2 + (sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2)).*V(s-1,i) ...                    
                    + ( 0.5.*((sigma_i.*I(i).*S(s)).^2 + (sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2))  ).* V(s-1,i+1)   ...
                    + ( 0.5.*((sigma_i.*I(i).*S(s)).^2 + (sigma_d.*I(i).*(1-I(i))).^2).*(dt./(hi.^2))  ).*V(s-1,i-1);  

             
                diff_v = diff_v + abs(V_next(s,i)-V(s,i)); if isnan(diff_v) == 1; stop; end; 
                iters = iters+1;
                                if isnan(L_optimal(s,i)) == 1; stop; end;          

        end
      
      
     end


    diffmax_policy = norm(L_optimal(:)-L_old(:))/norm(L_old(:)) ;
    diffmaxpolicy(iteration) = diffmax_policy;

    misses(inds2) = NaN;
    inds_fx = find(misses == 0);
    
    diffmax = r*diff_v / iters ;

    diffmax_iters(iteration) = diffmax;
    
    V = V_next ;
    
    L_old = L_optimal ;
    L_check = L_old;

    diff_v = 0 ;
    iters = 0;
    iter_check = 0;
    
    diffmax
    iteration
    diffmax_policy

%     toc;
    
end


save(filename_save);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
