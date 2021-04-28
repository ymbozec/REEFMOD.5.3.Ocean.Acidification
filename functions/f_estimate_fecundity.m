function [fecundity_adol, fecundity_adult] = f_estimate_fecundity (coral_cm2, META)

% Calculate fecundity output for a reef in terms of total number of larvae produced over 
% the grid for a coral species (will be transformed to larval input for other reefs using
% the connectivity matrix).
% 
% calculate fecundity for adol
id1 = zeros(size(coral_cm2)) ;
id2 = id1 ;

id1(coral_cm2 >= META.adol_size & coral_cm2 < META.adult_size) = 1 ;
fecundity_adol = sum(sum(0.25*2*pi*((sqrt(coral_cm2.*id1/pi)).^2)*216)) ;
% 
% % calculate fecundity for adults'
id2(coral_cm2 >= META.adult_size) = 1 ;
fecundity_adult = sum(sum(2*pi*((sqrt(coral_cm2.*id2/pi)).^2)*216)) ;


