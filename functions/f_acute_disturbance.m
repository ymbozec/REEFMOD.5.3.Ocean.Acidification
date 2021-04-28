% Created for Moorea simulations
function [coral, algal]=f_acute_disturbance(coral, algal, mortalities, META)

% Extract data from the structures (need to be filled again when leaving)
algal_cm2 = [algal.cover_cm2] ;
[coral_cm2, surface_cm2, volume_cm3, clade, colony_ID, species_ID] = f_struct_deploy (coral);

%%%% This is new stuff (Aug 2013) for implementing species-specific natural mortalities 
[i,j]=size(coral_cm2) ;
mortality = zeros(i,j);

col_start = 1;
col_stop = 0;

for s = 1:META.nb_coral_types
    
    col_stop = col_stop + species_ID(s) ; 
    mortality(:,col_start:col_stop)= mortalities(s);
    col_start = col_start + species_ID(s) ;
    
end

id = ones(size(coral_cm2)) ; % assigns 1 to every colony (faster than spones, here)
rand_mort = rand(size(coral_cm2));
id(rand_mort > mortality) = 0 ;
coral_cm2(id==1) = - coral_cm2(id==1) ;
algal_cm2(:,1) = algal_cm2(:,1) + sum((coral_cm2.*id), 2) ;

%%%%%%%% Before leaving, store the new covers into 'coral' and 'algal'%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, META.nb_coral_types, META.doing_3D);

for a=1:size(algal_cm2, 2) 
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end