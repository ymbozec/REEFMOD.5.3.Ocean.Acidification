%_________________________________________________________________________________________
%
%       DISEASE SETTINGS
%_________________________________________________________________________________________

% THE EFFECT OF DISEASE IS NOT IMPLEMENTED YET
% THIS IS JUST FOR FUTURE IMPLEMENTATION:

META.disease_names = { 'BBD','YBD','WPII','WBD' } ;
META.nb_disease_types = length(META.disease_names) ;

% INITIAL disease lesion size in cm2
META.disease_start_lesion_size = 1 ; % NEED TO PARAMETERIZE

%%%% BELOW ARE CORAL SPECIES SPECIFIC PARAMETERS
% -> needs as many values as there are coral species

% Min size (cm2) at which coral vulnerable get disease
META.disease_coral_min_size = [ 0 ; 0 ; 20 ; 100 ; 100 ];

% NOTE: CORAL.is_disease_vulnerable = [0, 0, 1, 1, 1]; to be added ????????????

% Currently parameterised by mean of Ernestos data from 1999-2005; also shown are data 
% from his review article with ranges.

% Prevalence of black band disease vector (0.002-0.06)
META.disease_BBD_prev = [ 0 ; 0 ; 0 ; 0.0029 ; 0 ];

% Prevalence of yellow blotch disease vector (0.01-0.45)
META.disease_YBD_prev = [ 0 ; 0 ; 0 ; 0.012 ; 0 ];

% Prevalence of white plague II vector - from Ernesto over all sites;
%WPII0.01 WPI in FLorida 0.036 - NEEDS PARAMETERISATION
META.disease_WPII_prev = [ 0 ; 0 ; 0.024 ; 0.018 ; 0 ];

% Prevalence of white band disease vector - Ernesto 0.0001-0.0011
META.disease_WBD_prev = [ 0 ; 0 ; 0 ; 0 ; 0.046 ];

% Linear extension rate in cm per 6 month iteration
META.disease_BBD_prog = 54 ; % 0.3-1 mm per d which is 54 cm per 6 mo (30*6)
META.disease_YBD_prog = 2.8 ; % 2.8 from bruno; ernesto review 0.01-0.12 mm per d which is the same as 1.8 to 21 cm/6 mo
META.disease_WPII_prog = 54 ; % 0.3-3 which is 54 cm to 540 cm per 6 months
META.disease_WBD_prog = 0.3 ; % unknown
