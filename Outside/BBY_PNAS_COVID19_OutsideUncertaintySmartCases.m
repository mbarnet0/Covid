% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Markov Approx. Method For Covid Policy with SA Aversion
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

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
savedata   = 1 ; 
initialize = 1 ; % == 0 if already have VF

filename_load = [directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1111']; % beta case 1, delta case 1, zeta case 1, alpha case 1
filename_new = [directory, '/Results_BBY_PNAS_COVID_OutsideUncertainty_1111_update']; % beta case 1, delta case 1, zeta case 1, alpha case 1


load(filename_load);


R_0
CFR
theta
aa

%__________________________________________________________________________
%
% Convergence parameters
%__________________________________________________________________________
% convergence tolerance level
tol = 5e-8;
tol2 = 1e-12;

diffmax = 100;
diffmax_policy = 100;

while  diffmax>tol2 || (diffmax_policy > tol) 


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

% stop
save(filename_new);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
