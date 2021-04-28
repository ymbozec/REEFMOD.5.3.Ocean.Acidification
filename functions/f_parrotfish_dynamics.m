%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created May 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [densities, biomasses, GI, size_distri_simulated, mean_adult_size, annual_Z_old, Catches, Harvest_rate] = ...
            f_parrotfish_dynamics(time_max, time_start_fishing, time_stop_fishing, F, L_min_select, L_min_exploit, do_stochastic)

% Assign parameter values
% -----------------------
% Define this first!
time_step = 4; % in month, time step of the model
nb_steps = time_max/time_step ;  % total number of time steps
nb_steps_before_fishing = time_start_fishing/time_step ;
nb_steps_stop_fishing = time_stop_fishing/time_step ;

P1 = 0.05 ; % to generate values within truncated Gaussian distribution 1-2*P % confidence interval 
% P1 = 0.4; %NEW: to reduce output variability which is huge with 0.05!
P2 = 1 - P1 ;

if do_stochastic == 1
    L_star = norminv(P1+(P2-P1)*rand(nb_steps,1),28.4,0.729);
    mu_1 = norminv(P1+(P2-P1)*rand(nb_steps,1),0.0081,0.002);
    mu_2 = norminv(P1+(P2-P1)*rand(nb_steps,1),-0.325,0.106);
    r = norminv(P1+(P2-P1)*rand(nb_steps,1),23.50,5.775);
else
    L_star = 28.4*ones(nb_steps,1);
    mu_1 = 0.0081*ones(nb_steps,1);
    mu_2 = -0.325*ones(nb_steps,1);
    r = 23.50*ones(nb_steps,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With Bermuda parametrisation
% Linf = 36 ;
% rand_L_star = rand_L_star - 2;
% rand_mu_2 = rand_mu_2 - 0.05 ;
% rand_r = rand_r/10 ;
 
% With Bonaire parametrisation
Linf = 39 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 1.99698 - 0.03665 * Linf ; % empirical K/Linf relationship (see R script 'GROWTH_ANALYSIS_NEW.R")
L = 1:1:Linf ; % mid_point values of the length classes

% Parameters of the L-W relationship
a = 0.025 ; % Bohnsack and Harper (1988) -> see FishBase
b = 2.921 ; % Bohnsack and Harper (1988) -> see FishBase
W = a*L'.^b ; % vector of individual weight (g) at length (cm FL)

% Grazing calculations; correspondence between sizes and phase are derived from Bruggemann et al. (1996) (Table 1):
% 73% of the IP individuals are in the range 15-29 cm and 74% of the TP are >=30cm FL. 
coeff_phase = zeros(Linf,1);
coeff_phase(1:14,1) = 0.84 ; % for juveniles
coeff_phase(15:29,1) = 1 ; % for intermediate phase
coeff_phase(30:Linf,1) = 0.80 ; % for terminal phase (but minimum size of TP can be 25cm)
BR = coeff_phase.*((1088 - 17.12.*L')-56) ; % Bite rate in hour-1 from Mumby (2006). Note formulas 2a and 2b
% in Mumby (2006) are actually for Sparisoma and Scarus, respectively (and
% not the opposite as published). Based on Bruggemann et al. (1994 II) Fig.4B
% -> feeding rate (bites.h-1 = 1088.84 - 17.12*FL
% -> correction factor (deep reef): Fig 4C
areabitten = 5.839*0.0001*L.^2 ; % area bitten per bite, in cm2 from Mumby (2006) based on Bruggemann et al. (1994).
% Note that Bruggeman et al. (1994, 1996) used coeff 0.0001 when Mumby used 0.001
% Note that factor is 5.839 (SD=0.174, n=414) in Bruggeman et al. (1994 I) but 5.257 (SD=0.174, n=414) in Bruggeman et al. (1996)
IGR = BR.*areabitten'/10000 ; % individual grazing rate = area bitten in m2 per hour for a fish of a given length


% Build the growth matrix (note the GM accounts for unsertainty already!)
% -------------------------------------------------------------------------
% (code derived from the function growtrans of the R package "fishmethods")
L_min = 1; % Mid-point of starting size class.
L_max = Linf + 1; % Mid-point of end size class. This should be one increment larger than Linf.
Linc = 1 ; % size increments

K = K*time_step/12 ; % growth rate scaled to the model time step
SEK = 0.07605494 *time_step/12 ; % standard error from the empirical K magnified by 'mag' (and scaled to time step)
SELinf = 1.754034 ; % standard error from the empirical Linf magnified by 'mag'
rhoLinfK = -0.8452311 ; % correlation coefficient K/Linf (not affected by time scale)

rhoLinfK = abs(rhoLinfK) ;
Ln = L_min:Linc:L_max+1 ;
COV = rhoLinfK * SEK * SELinf ;
DL = (Linf - Ln) * (1 - exp(-K)) ;  % equation #6 of Chen et al. (2003)
DL = DL(1:length(DL(DL>=0))+1) ;
Ln = Ln(1:length(DL(DL>=0))+1) ;
DL(DL<0)=0 ;

VL = SELinf^2 * (1 - exp(-K))^2 + ((Linf - Ln).^2) * SEK^2 * ...
    exp(-2 * K) - 2 * COV * (1 - exp(-K)) * (Linf - Ln) * exp(-K) ;%  correct equation (8)

growmat = zeros(length(Ln)-1, length(Ln)-1) ;

for m = 1:(length(Ln)-1) % semi-vectorized version of the original loop
    temp = normcdf(Ln(m+1), Ln+DL, sqrt(VL(m)))' - normcdf(Ln(m), Ln+DL, sqrt(VL(m)))' ;
    growmat(1:m,m) = temp(1:m,1) ;
end

rowsums = sum(growmat,2) ;
GM_bermuda = growmat./ rowsums(:,ones(1,size(growmat,1))) ;

% Set up recruitment
% -------------------
add_recruit = zeros(size(L,2), nb_steps);
add_recruit(1,2:nb_steps) = r(2:nb_steps) ;  % defined for every time step

% Set up instantaneous rate of natural mortality
% -----------------------------------------------
M = zeros(length(L),nb_steps);
for k = 1: length(L)
    M(k,:) = mu_1.*((L_star-L(k)).^2) + mu_2 ; % instantaneous natural mortality rate (time scale of the model)
end
M(M<0)=0 ;

% Set up instantaneous rate of fishing mortality
% -----------------------------------------------
selectivity = zeros(length(L),nb_steps);
selectivity(L_min_exploit:Linf,nb_steps_before_fishing:nb_steps_stop_fishing)= 1 ; % select harvested body lengths during fishing scenario

F1 = selectivity.*F(ones(1,length(L)),:) ; % instantaneous fishing mortality rate for every body length selected by the fishery

% Deduce survivorship 
% -----------------------------------------------
Z = M + F1 ; % Total instantaneous mortality rate at every time step
S = exp(-Z) ; % survivorship for every fish length (proportion of fish that survive during the time step)

% Run simulations
% ---------------
densities = zeros(size(L,2), nb_steps);
densities(:,1) = add_recruit(:,1); % initialise population with the recruitment of one cohort
total_survival_old=zeros(nb_steps,1); % to track actual survivorship of fish older than 2yr (>=15cm)
total_survival_young=zeros(nb_steps,1); % to track actual survivorship of fish younger than 2yr (<15cm)
L_survival = find(L==L_min_select) ;
mean_adult_size=zeros(nb_steps,1); % to track average size of fish older than 2 yr
Catches = zeros(size(L,2), nb_steps);
Harvest_rate = zeros(nb_steps,1);

for t = 2:nb_steps

    % Calculate survival of the population (individuals larger than L_survival)
    POP_before = densities(:,t-1) + add_recruit(:,t);
    
    % Grow first
    POP_after_growth = (POP_before'*GM_bermuda)' ;
    
    % Then survive (total mortality)
    POP_after_death1 = POP_after_growth.*S(:,t) ;
    
    % Record population at the end of the time step
    densities(:,t) = POP_after_death1 ;
    total_survival_old(t) = sum(POP_after_death1(L_survival:length(L)))/sum(POP_after_growth(L_survival:length(L)));
    total_survival_young(t) = sum(POP_after_death1(1:(L_survival-1)))/sum(POP_after_growth(1:(L_survival-1)));

    % Estimate catches using Baranov catch equation
    Catches(:,t) = (F1(:,t)./(M(:,t)+F1(:,t))).*POP_after_growth.*(1-S(:,t)).*W ; % numbers-at-length in the catch times individual weights  
    Catches(isnan(Catches(:,t))==1,t)=0; % remove NaNs generated when F1=0
    
    % Estimate harvest rate = proportion of the exploitable (> L_min_exploit) biomass harvested
    Harvest_rate(t) = sum(Catches(:,t),1)/sum(selectivity(:,t).*POP_after_growth.*W) ; % 
    
    % Estimate mean body legnth of the surviving stock 
    mean_adult_size(t) = sum(densities(L_survival:Linf,t).*[L_survival:1:Linf]')/sum(densities(L_survival:Linf,t));
    
end

% Arrange model outputs
% ---------------------

% Group size classes according to the size decomposition of Bonaire
size_distri_simulated = [ sum(densities(5:9,:)) ; sum(densities(10:14,:)) ; sum(densities(15:19,:)) ; sum(densities(20:24,:)) ; ...
            sum(densities(25:29,:)) ; sum(densities(30:34,:)) ; sum(densities(35:Linf,:))]' ;

% Calculate total annual survival and mortality (equivalent to a 4mo time step)
annual_survival_old = total_survival_old.^(12/time_step); % yr-1
% annual_survival_young = total_survival_young.^(12/time_step); % yr-1

% Calculate total instantaneous mortality rate (yr-1)
annual_Z_old = -log(annual_survival_old) ; % instantaneous total mortality yr-1 for each fish length
% annual_Z_young = -log(annual_survival_young) ; % instantaneous total mortality yr-1 for each fish length

% Estimate biomass for every fish length at every time step
biomasses = densities.*W(:,ones(1,nb_steps));

% Estimate grazing intensity in per cent area grazed per hour for each size class
% this accounts for fish density in each size class - we divide by the
% reference area that was used to calculate density (100 m2)
% multiplied by 100%
GI = 100*densities.*IGR(:,ones(1,nb_steps))/100; 

%%%%%%%%%%%%%%%%%
%%%% Et voila!...
%%%%%%%%%%%%%%%%%