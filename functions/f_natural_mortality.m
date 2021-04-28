%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified in August 2013 after changes in parameterisation as shown in Mumby et al. 2013 appendix

function [coral_cm2, algal_cm2]=f_natural_mortality(coral_cm2, algal_cm2, CORAL, species_ID, META)

%%%% This is new stuff (Aug 2013) for implementing species-specific natural mortalities 
[i,j]=size(coral_cm2) ;
sensitivity_whole_natural = zeros(i,j);
sensitivity_partial_natural = sensitivity_whole_natural;
extent_partial_natural = sensitivity_whole_natural ;

col_start = 1;
col_stop = 0;

for s = 1:META.nb_coral_types
    
    col_stop = col_stop + species_ID(s) ; 
    sensitivity_whole_natural(:,col_start:col_stop)= CORAL.sensitivity_whole_natural(s);
    sensitivity_partial_natural(:,col_start:col_stop)= CORAL.sensitivity_partial_natural(s);
    extent_partial_natural(:,col_start:col_stop)= CORAL.extent_partial_natural(s);
    col_start = col_start + species_ID(s) ;
    
end

%----------------------------------
% WHOLE COLONY MORTALITY
%----------------------------------
id1 = ones(size(coral_cm2)) ; % assigns 1 to every colony (faster than spones, here)

id1(coral_cm2 < META.adol_size) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)
id_adult = id1;
id_adult(coral_cm2 < META.adult_size) = 0;
id_adol = id1 - id_adult;

rand_mort = rand(size(coral_cm2));
rand_mort_adol = rand_mort.*id_adol ;
rand_mort_adult = rand_mort.*id_adult ;

% switch to 0 the adolescents escaping mortality due to chance
id_adol(rand_mort_adol > CORAL.adol_whole_mortality_rate*sensitivity_whole_natural) = 0 ;

% switch to 0 the adults escaping mortality by chance
id_adult(rand_mort_adult > CORAL.adult_whole_mortality_rate*sensitivity_whole_natural) = 0 ;

% dead corals are those remaining (id=1)
id_dead = id_adol + id_adult;

% Update turf
% algal_cm2(:,1) = algal_cm2(:,1) + sum(coral_cm2.*id_dead,2) ;

% Update coral cover (turn into negative size the dead colonies)
coral_cm2(id_dead==1) = - coral_cm2(id_dead==1) ;

%----------------------------------
% PARTIAL COLONY MORTALITY
%----------------------------------
rand_mort2 = rand(size(coral_cm2));
rand_mort2(coral_cm2<0) = 0 ; % turn mortality to 0 for those already dead (negative sizes)

% allocate space
area_lost = 0*rand_mort2 ;
log_calc = area_lost;

log_calc(coral_cm2>0) = log(coral_cm2(coral_cm2>0));

proba = 1-(CORAL.partial_mortality_inci_int + ...
    CORAL.partial_mortality_inci_gra * log_calc)/100 ;

% Adjust the mortality for each species
proba = proba.*sensitivity_partial_natural ;

lost_factor(rand_mort2 < proba) = CORAL.partial_mortality_area_int + ...
    CORAL.partial_mortality_area_gra * log_calc(rand_mort2 < proba) ;

lost_factor(lost_factor<0)=0;  % because the formula above generates negatives

% Since Mumby et al. (2013) the area lost for every affected colony is:
area_lost(rand_mort2 < proba) = (exp(lost_factor(rand_mort2 < proba))-1)/100 ;

area_lost = round(extent_partial_natural.*area_lost);

% cannot loose more than the colony area
area_lost(area_lost > coral_cm2 & coral_cm2>0) = coral_cm2(area_lost > coral_cm2 & coral_cm2>0) ;

%%%% Update the grid
coral_cm2 = coral_cm2 - area_lost;  % remove lost tissues
% algal_cm2(:,1) = algal_cm2(:,1) + sum(area_lost,2) ; % update turf with coral loss
