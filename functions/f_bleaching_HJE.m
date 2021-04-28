%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: Aug 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coral, algal] = f_bleaching_HJE (coral, algal, bleaching_mortality, t, CORAL, META)
% MODIFIED FOR USE IN RESILIENCE MAPPING - INPUT DHW INCLUDES ALL TIME STEPS
% Simulates bleaching of corals at elevated sea surface temperatures.  
% Simulates partial and total mortality, and does not model the recovery of initially affected corals.
% Brooders and spawners may be affected differently by raised temperatures.

% What do we mean by survival now that there are bleaching events of differing severity?
% Previously survival meant that you had been present (and not died) during an event in
% which DHMs exceeded say 3. Now all individuals are present during minor bleaching events. 
% If never bleached before have higher risk of dying-> history = 1
% If don't die then history becomes -1
% If previously exposed have a lower risk. Stay -1 or die with lower prob.

% Parameters of bleaching mortalities are drawn from a distribution (in 'runmodel.m'),
% having had the number of DHWs passed to it in 'main.m'.
% -> bleaching_mortality.whole_colony(t) 
% -> bleaching_mortality.partial_on_brooders(t)
% -> bleaching_mortality.partial_on_spawners(t)


% Extract data from the structures (need to be filled again when leaving)
algal_cm2 = [algal.cover_cm2] ;
[coral_cm2, surface_cm2, volume_cm3, clade, colony_ID, species_ID] = f_struct_deploy (coral);

% Locate brooders and spawners
id_brooders = zeros(size(coral_cm2)) ; % initialize brooders id
id_spawners = zeros(size(coral_cm2)) ; % initialize spawners id

%%%% This is new stuff (August 2013) for implementing species-specific bleaching mortalities 
[i,j]=size(coral_cm2) ;
sensitivity_bleaching = zeros(i,j);
extent_bleaching = sensitivity_bleaching ;
proba_switching = sensitivity_bleaching ;

col_start = 1;
col_stop = 0;

for s = 1:META.nb_coral_types
    
    col_stop = col_stop + species_ID(s) ; 
    sensitivity_bleaching(:,col_start:col_stop)= CORAL.sensitivity_bleaching(s);
    extent_bleaching(:,col_start:col_stop)= CORAL.bleaching_partial_reduction(s);
    proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
    
    if CORAL.is_brooder(s) == 1
        id_brooders(:,col_start:col_stop) = 1 ;
    else
        id_spawners(:,col_start:col_stop) = 1 ;
    end
    
    col_start = col_start + species_ID(s) ;
    
end

% Clade-induced tolerance to thermal stress
sensitivity_bleaching(clade==2) = sensitivity_bleaching(clade==2) * CORAL.bleaching_tolerance_clade ;

%________________________________
%
% Total colony mortality
%________________________________

% Do this first as there's no point in making a colony suffer partial mortality before it dies

% Generate random mortalities for adol + adults
id1 = ones(size(coral_cm2)) ; % assigns 1 to every colony
id1(coral_cm2 < META.adol_size) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)
rand_mort1 = rand(size(id1)) ; % Generates random probability for adol + adults

prob_whole_mortality = id1.*sensitivity_bleaching * bleaching_mortality.whole_colony(t);

id_dead = id1 ;
id_dead(rand_mort1 > prob_whole_mortality) = 0 ; % exclude the survivors
algal_cm2(:,1) = algal_cm2(:,1) + sum((coral_cm2.*id_dead), 2) ;
coral_cm2(id_dead==1) = - coral_cm2(id_dead==1);  % now dead (negatives)

%________________________________
%
% Partial mortality
%________________________________

% These are used to update the status of each coral. If a coral has never been bleached
% and then survives a bleaching event it's status changes from 0 to 1.
% NOTE from JH: the record of previous bleaching doesn't  have any effect on likelihood
% of partial mortality whereas it does on total mortality, doesn't make sense.
% NOTE from YM: we now (01/2015) apply the same reduction to the probability of partial mortality 

id2 = id1 - id_dead ;
rand_mort2 = rand(size(id2)) ; % Generates random probability for adol + adults

id_brooders = id_brooders.*id2 ;
id_spawners = id_spawners.*id2 ;

prob_partial_mortality = id2 ;
prob_partial_mortality(id_brooders==1) = sensitivity_bleaching(id_brooders==1) * bleaching_mortality.partial_on_brooders(t) ;
prob_partial_mortality(id_spawners==1) = sensitivity_bleaching(id_spawners==1) * bleaching_mortality.partial_on_spawners(t) ;

id_part = id2 ;
id_part(rand_mort2 > prob_partial_mortality) = 0 ; 
bleach_extent = floor(extent_bleaching.* coral_cm2 .* id_part) ; 
algal_cm2(:,1) = algal_cm2(:,1) + sum(bleach_extent, 2) ;
coral_cm2 = coral_cm2 - bleach_extent ;

%________________________________
%
% Clade switching
%________________________________

rand_switch = rand(size(id2)) ; % Generates random probability of switching to the thermally-tolerant clade (clade 2)
clade(rand_switch < proba_switching & coral_cm2>0) = 2 ; % NOTE THIS INDEPENDENT OF BLEACHING MORTALITY


%%%%%%%% Before leaving, store the new covers into 'coral' and 'algal'%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, META.nb_coral_types, META.doing_3D);

for a=1:size(algal_cm2, 2) 
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end
