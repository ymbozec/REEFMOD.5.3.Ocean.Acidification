%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2015.
% Last modified: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TEMPORARY -> will be replaced by Karlo's code of COTs dynamical consumption
% These are the COTS feeding preferences for the six coral types.
% Tabular=0.393;
% CorymboseAcro=0.258;
% BranchingThicket=0.12;
% CorymboseOther=0.12;
% SmallMassive=0.084;
% LargeMassive=0.025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To implement here the generation of food preferences with Karlo's code
% DOUBLE-CHECK THE STRUCTURE OF THE PREFERENCE MATRIX REGARDING THE ORDER
% OF SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eaten_coral_cm2] = f_COTs_dynamics(reef_coral_cm2,Reef_areas_cm2,feeding_rate,COTs_density)

% Reef_areas_cm2 is the total area of every reef
% reef_coral_cm2 gives the cover area (cm2) of each coral species (columns) in every reef (rows)

nb_reefs = size(reef_coral_cm2,1) ;
nb_coral_species = size(reef_coral_cm2,2) ; %morphological groups used in Ortiz et al. 2014 Nature Climate change

reef_coral_cover = reef_coral_cm2./Reef_areas_cm2(:,ones(1,nb_coral_species));
reef_coral_relative = reef_coral_cover/sum(reef_coral_cover);

COTs_food_preference = [0.12 ; 0.393 ; 0.258 ; 0.12 ; 0.084 ; 0.025];

food_pref_matrix = COTs_food_preference(:,ones(1,nb_reefs))' ;

cots_density_per_reef = Reef_areas_cm2(:,ones(1,nb_coral_species))*COTs_density/(1e6*10000) ; %for the modelled reef
eaten_coral_cm2 = round(feeding_rate*cots_density_per_reef.*food_pref_matrix.*reef_coral_relative/sum(food_pref_matrix.*reef_coral_relative));
