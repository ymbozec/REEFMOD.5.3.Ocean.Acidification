%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 31/05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    SETTINGS OF CORAL COMPETITION TO REDUCE CORAL GROWTH RATE IN A CELL   
%__________________________________________________________________________

% This estimates a vector of percent reduction to apply to growth rates for
% a given value of coral cover. Based on empirical relationships derived
% from observations and aquarium experiments in Moorea

% Coral_cover = 1:100 ;
% Coral_contact = 0*Coral_cover;
% % 1) Empirical relationship giving the average percentage contact of a
% % coral colony with any other coral (regardless of the species) as a
% % function of coral cover (from Nick's quadrats in Moorea).
% Coral_contact(Coral_cover<5) = 2.71 * Coral_cover(Coral_cover<5) ; % (Colony size was found to not affect this relationship)
% Coral_contact(Coral_cover>=5) = 10.87 + 0.46 * Coral_cover(Coral_cover>=5) ; % (Colony size was found to not affect this relationship)
% Coral_contact(Coral_contact>100)=100 ; % force to 100% just in case
Coral_contact = 0:100 ; %101 values (now includes 0)

% 2) Empirical relationship between percent growth reduction and percent
% contact (from Nick's experiment in Moorea). Needs to differentiate btw
% conspecifics (Pocillopora) and heterospecifics (Acropora) - will be
% weighted according to the relative cover of the two species
CORAL.growth_reduction_cspec = zeros(length(Coral_contact),META.nb_coral_types);
CORAL.growth_reduction_hspec = zeros(length(Coral_contact),META.nb_coral_types);

% % ambient relationships %%%%%%%%%%%%%%%%%%%%%%%%%%%
% proportional change (decrease) in lateral extension as the percent contact increases - scaled to 0-1
% For POCILLOPORA
CORAL.growth_reduction_cspec(:,1) = 1 - 0.005 * Coral_contact ; % with conspecifics
CORAL.growth_reduction_hspec(:,1) = 1 * exp(-0.018 * Coral_contact) ; % with heterospecifics

% For ACROPORA
CORAL.growth_reduction_cspec(:,2) = 1 - 0.005 * Coral_contact ; % with conspecifics
CORAL.growth_reduction_hspec(:,2) = 1 * exp(-0.01 * Coral_contact) ; % with heterospecifics

% For MONTIPORA
CORAL.growth_reduction_cspec(:,3) = 1 - 0.005 * Coral_contact ; % with conspecifics
CORAL.growth_reduction_hspec(:,3) = 1 * exp(-0.022 * Coral_contact) ; % with heterospecifics

% For PORITES
CORAL.growth_reduction_cspec(:,4) = 1 - 0.005 * Coral_contact ; % with conspecifics
CORAL.growth_reduction_hspec(:,4) = 1 * exp(-0.026 * Coral_contact) ; % with heterospecifics


% % % OA relationships %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % For POCILLOPORA
% CORAL.growth_reduction_cspec(:,1) = 1 * exp(-0.0097 * Coral_contact) ; % with conspecifics
% CORAL.growth_reduction_hspec(:,1) = 1 * exp(-0.063 * Coral_contact) ; % with heterospecifics
% 
% % For ACROPORA
% CORAL.growth_reduction_cspec(:,2) = 1 * exp(-0.0097 * Coral_contact) ; % with conspecifics
% CORAL.growth_reduction_hspec(:,2) = 1 * exp(-0.035 * Coral_contact) ; % with heterospecifics
% 
% % For MONTIPORA
% CORAL.growth_reduction_cspec(:,3) = 1 * exp(-0.0097 * Coral_contact) ; % with conspecifics
% CORAL.growth_reduction_hspec(:,3) = 1 * exp(-0.077 * Coral_contact) ; % with heterospecifics
% 
% % For PORITES
% CORAL.growth_reduction_cspec(:,4) = 1 * exp(-0.0097 * Coral_contact) ; % with conspecifics
% CORAL.growth_reduction_hspec(:,4) = 1 * exp(-0.091 * Coral_contact) ; % with heterospecifics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%