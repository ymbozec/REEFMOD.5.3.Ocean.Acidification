% N.R. Evensen, MSEL, created May 2016.
% Last modified: 31/05/2016
%
% GENERATES DEFAULT PARAMETERS FOR REEFMOD
%
% This is default parameters. It is recommended to keep those values and to declare specific parametrisation
% in a separate file
%_________________________________________________________________________________________
%_________________________________________________________________________________________
%
%       SET UP META-LEVEL PARAMETERS
%_________________________________________________________________________________________

META.nb_reefs = 1 ; % Number of populations (= reefs)

META.nb_time_steps = 10 ; % Number of time steps (= 6 mo) to run the model for (default is 20 years)

META.nb_simul = 100 ; % Number of simulations for each population

META.max_colonies = 15 ; % Maximum number of colonies of a given species in a cell

META.recruitment_type = 0 ; % set to 0 for full, 1 for random (when only 1 reef, no connectivity)

META.space_limited_settlement = 1 ; % set to 0 for no effect of available space on the probaility of recruitment

META.use_connectivity_matrix = 0 ; % set to 0/1 for using/not using a connectivity matrix

META.do_parrotfish_fishing = 0 ; % set to 1 to allow grazing impact changing as a response to fishing
% NOTE this requires running f_parrotfish_dynamics to estimate grazing over time

META.doing_COTs = 0 ; % set to 1 for simulating coral predation by the crown-of-thorn starfish (Pacific)

META.do_size_frequency = 0 ; % set to 1 to estimate the number of colonies in different size classes

META.grid_x_count = 20; % number of grid cells along the x-edge
META.grid_y_count = 20; % number of grid cells along the y-edge

META.cell_x_size = 100 ; % size of a single cell in centimeters along its x-edge
META.cell_y_size = 100 ; % size of a single cell in centimeters along its y-edge

META.track_populations = 0 ;% set to 1 for tracking the size of every single coral colony (for John Hedley)
META.doing_coral_competition = 1 ;

%_________________________________________________________________________________________
%
%       SET UP REEF PARAMETERS
%_________________________________________________________________________________________

REEF.exposure = 0 ; % (replaces the previous 'is_leeward' that was affecting Dictyota)
% proportion of exposure between 0 (leeward glovers) and 1 (windward glovers).
% exposure used to slow down rate of growth of Dictyota to maximum of 43% of windward value (when exposure =0);

REEF.dictyota_declines_seasonally = 1 ; % set to 0/1 to switch off/on Dictyota dying back in winter

REEF.nongrazable_substratum = 0.05 ; % non-grazeable cover

REEF.herbivory = 1 ; % total amount of algae that can be herbivorised by fish per time step
%expressed as a proportion of the reef area - now acts as a constraint to max amount of seabed

REEF.diadema = 0 ; % total amount of algae that can be herbivorised by Diadema per time step
%expressed as a proportion of the reef area.

%_________________________________________________________________________________________
%
%       SET UP CORALS
%_________________________________________________________________________________________

%%%%% 1) SPECIES-SPECIFIC PARAMETERS (need as many values there are coral types)
% Parameters are in the given order:
% (1) POCILLOPORA , (2) ACROPORA, (3) MONTIPORA, (4) PORITES [MASSIVE]

% Code indicating coral types that are brooder (1) or spawner (0)
CORAL.is_brooder = [0 ; 0 ; 0 ; 0];

% Initial coral covers may vary from 2.5 to 25 percent - keep initial coral cover the same for each species - WHY?
CORAL.initial_cover = [0.00 ; 0.00 ; 0.00 ; 0.2]; 

CORAL.growth_rate = [1.6 ; 2.5; 1.5 ; 0.65]; % Horizontal growth rate (radius extension) in cm per time interval
% 1- POCILLOPORA SPP. (based on P. verrucosa - in my Marine Bio. paper)
% 2- ACROPORA SPP. (based on young A. hyacinthus and A. nasuta - from Coral Trait database)
% 3- MONTIPORA (based on Browne 2012 (reported in my MEPS paper))
% 4- Massive Porites (based on JC's parameterization for massive corals)

% CORAL.growth_rate = [1.12 ; 1.575; 0.855 ; 0.65]; % UNDER OA
% 1- POCILLOPORA SPP. - 30% decrease compared to ambient growth (from my Marine Bio. paper)
% 2- ACROPORA SPP. - 37% decrease compared to ambient growth (from where???)
% 3- MONTIPORA - 43% decrease compared to ambient growth (from my MEPS paper)
% 4- Massive Porites - 0% decrease compared to ambient growth (from Fabricius et al. 2011)

%%%%%%% NEW (08/2015): Maximum density of settlers after recruitment %%%%%%%%%%%%%%%%%%%%%%%
% This replaces the previous parameter of probability of recruitment
%( CORAL.recruit_prob =[ 0.025 ; 0.85 ; 0.72 ; 1 ; 0.48 ; 0.06 ] )
META.max_density_settlers = [1 ; 0.25; 0.5; 1]; 
% This is the maximum number of larvae able to settle in a cell if space is 100% available.
% The actual number of settlers in a cell is (stochastically) determined by the proportion of turf in that cell.
% Obtained by multiplying the previous vector of CORAL.recruit_prob by a given factor (to be adjusted).
% Note this also depends on the total number of colonies allowed to populate a cell (cannot recruit more than META.max_colonies)
% From now, corals recruit BEFORE any other process (growth, mortality etc.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum diameter of coral colonies to estimate max size (cm2)
% Note this is necessary for large cells (> 100 cm)
CORAL.max_diameter = [50 ; 100; 50; 250]; % NOTE this will be further limited by the size of a cell

% Define the initial size distribution of coral colonies
% 1 for a random selection of corals of different sizes - i.e. equal probability per size
% 2 for lognormal distribution (based on Meesters et al. 2001 in curacao)
CORAL.size_model = [2 ; 2 ; 2 ; 2]; % USE MODEL 2 FROM NOW, OTHERWISE THE NEW PARTIAL MORTALITY
% HAS A DISPROPORTIONATE EFFECT ON LARGE COLONIES (colonies with a 100x100 cm cell size are much bigger
% compared to a 50x50 cm cell size model

% Below are parameters for generating coral colonies at initial steps based
% on a lognormal distribution (using Meesters et al 2001 estimates).
% Note that generated sizes will be further limited by CORAL.max_diameter and META.cell_area_cm2
% (the size of a colony cannot exceed the size of a cell). 
% Note that CORAL.size_mean / CORAL.size_var can be used to parameterize other kind of distribution
% given a new model is specified in f_derivecoralcover

% Mean of colony size (cm2) - based on Nick's quadrat measurements,
% therefore on the smaller end of the spectrum as it's early recovery
CORAL.size_mean = [log(65.07) ; log(73.73); log(50); log(50)]; % BASED ON CURACAO - Meesters (surface area)


% Variance of colony size - based on Nick's quadrat measurements,
% therefore on the smaller end of the spectrum as it's early recovery
CORAL.size_var = [log(1.74^2) ; log(1.75^2) ; log(1.75^2) ; log(1.75^2)]

% Rate at which coral recedes when in contact with Lobophora;
% 4 from Agaricia / Porites; % 0.75 from Meandrina;
% 0.75 from Mycetophyllia (but use value from Agaricia if assume few Mycetophyllia);
% 0.5 from Montastraea
CORAL.lobophora_reduce_rate = [0.75 ; 0.75 ; 0.75 ; 0.75]/7; 
% CORAL.lobophora_reduce_rate = [ 8 ; 0.75 ; 0.75 ; 0.5]/7; % JC's parameterisation

% Rate at which coral recedes when in contact with Lobophora;
% = 0.25*6 months from Lirman 2001 (ranges from 0.25 cm per month to 0.43) (replaces DICT_OVER_SPAWNER?)
CORAL.dictyota_reduce_rate = [1.5 ; 1.5 ; 1.5 ; 1.5]; % for spawners only
% CORAL.dictyota_reduce_rate = [ 0 ; 2.58 ; 0 ; 2.58 ]; % JC's parameterisation


%%%%% 2) THE FOLLOWING CORAL PARAMETERS ARE NOT SPECIES-SPECIFIC %%%%%%%%

% Percent of colony growth that is actively overgrow other corals
% set arbitrarily so that up to 20% of the colony's projected growth can overgrow smaller corals
CORAL.fractional_overgrowth = 0.2 ; 

% size of full fecundity
CORAL.adult_size = 50 ; % allow for greater sizes?

% size of partial fecundity
CORAL.adol_size = 30 ;

% Define the bins of the size classes for (1) juvenile, (2) adolescent and (3) adult colonies
CORAL.size_bins = [5 ; 10 ; 100];

% Intercept for incidence relationship from Meesters - need to find this
%  CORAL.partial_mortality_inci_int = 60; %60 in C++ ; Using this requires changes in the calculation of probability
CORAL.partial_mortality_inci_int = 88.9;% Meesters et al. 1997 (see Mumby et al. 2013)

% Gradient for incidence relationship from Meesters - need to find this
%  CORAL.partial_mortality_inci_gra = -12; %-12 in C++ ; Using this requires changes in the calculation of probability
CORAL.partial_mortality_inci_gra = -11.2; % % Meesters et al. 1997 (see Mumby et al. 2013)

% Intercept for area relationship from Meesters - need to find this
%  CORAL.partial_mortality_area_int = -0.5 ; % old way
CORAL.partial_mortality_area_int = -2.9 ; % Mumby et al. 2013

% Gradient for area relationship from Meesters - need to find this
%  CORAL.partial_mortality_area_gra = 1.1 ; % old way
CORAL.partial_mortality_area_gra = 1.59 ; % Mumby et al. 2013

CORAL.adult_whole_mortality_rate = 0.01 ; % usually 0.01

CORAL.adol_whole_mortality_rate = 0.02 ; % usually 0.02

% probability of coral recruits removed
CORAL.parrotfish_predation = 0.05 ;

% max coral size in cm2 in which parrotfish predation occurs
CORAL.threshold_predation_size = 6 ;

% New stuff (August 2013) -> species-specific sensitivity to natural mortality
CORAL.sensitivity_whole_natural = [4 ; 4.5 ; 3 ; 1] ;
CORAL.sensitivity_partial_natural = [1 ; 1 ; 1 ; 1] ;
CORAL.extent_partial_natural = [1 ; 2 ; 1.5 ; 1] ; % multiplicator affects the extent of area lost. This affects (magnitude) the output calculation of size-dependent area lost


%__________________________________________________________________________________________
%
%       SET UP ALGAE
%__________________________________________________________________________________________
% Initial algal covers in the following order
% (1) TURF ; (2) DICTYOTA ; (3) LOBOPHORA ; (4) HALIMEDA, but HALIMEDA processes are not implemented yet.
ALGAL.initial_cover = [0 ; 0 ; 0 ; 0]; % Turf MUST be set to 0 here

% Relative proportion of consumption of the different algal types by fish herbivory
% [0.6 0.36 .03 0.01]  from separate analysis of fish community and grazing pattern
% - shifted hal to lob from glovers west fish community)
ALGAL.herbivory_props = [0.60 ; 0.36 ; 0.04 ; 0];

% Relative proportion of consumption of the different algal types by Diadema
ALGAL.diadema_props = [0.8 ; 0.2 ; 0 ; 0];

% Proportional reduction in algal growth rate from contact with coral (actually used as 1-0.25 in the code)
ALGAL.coral_reduce_macrogrowth = 0.25 ;

% Local neighborhoud of algal that smothers adolescent corals
ALGAL.critical_algal_contact = 0.40 ;

% use 0.8 but could drop to 50% based on Foster et al
ALGAL.vcritical_algal_contact = 0.80 ; %0.5 now
% ALGAL.vcritical_algal_contact = 0.50 ; % Juan's parameterisation

% Effect of macroalgae on growth rate of coral recruits (based on Box) == PROB_DICTYOTA in parameters2.m ???
ALGAL.macroalgal_coral_recruit_growth_rate = 0.3 ;

% Effect of macroalgae on growth rate of adol and adult corals (based on Lirman for Dicty vs Ag and Pa)
ALGAL.macroalgal_coral_growth_rate = 0.1 ;

% Min percent cover in winter (ATTENTION: 0.02??)
ALGAL.dict_winter_min = 2 ; % Not used in the original code

%__________________________________________________________________________________________
%
%       SET UP HURRICANES
%__________________________________________________________________________________________

META.doing_hurricanes = 0 ; % set to 0/1 to switch off/on hurricanes

META.randomize_hurricane_chronology = 0 ; % set to 1 to randomize hurricane chronology for each simulation

META.random_hurricanes = 1 ; % if no prescribed scenario, just random hurricanes
REEF.hurr_strike_proba = 0.1; % Frequency at which hurricanes occur in summer, if no prescribed hurricane scenario

META.randomize_hurricanes = 0 ; % set to 1 to randomize hurricane for each simulation

META.hurricane_effect_on_macroalgae = 0.9 ; % proportion of macroalgae removed by a hurricane

META.hurricane_effect_on_recruits = 0.8 ; % prob that recruits die if less than 10 cm diameter

CORAL.sensitivity_hurricane = [6.36 ; 15; 6.66 ; 1]; % relative sensitivity of coral species to hurricanes - based on Adjeroud et al. 2009 data, effect of 1991 cyclone on each coral genus
CORAL.extent_hurricane = [1 ; 1 ; 1 ; 1]; % multiplicator affecting the extent of partial mortality due to hurricanes
% A multiplicator of the hurricane mortality probabilities (whole and partial colonies)
% -> dampens hurricane damages for robust species

META.hurri_whole_mort_dead = 0; % Multiplicator of CORAL.hurr_mortality_rate for the dead colonies
% 0 produces no effect of hurricane on dead colonies (they are not removed)
% 1 applies the same relative effect than specified for living colonies
% 2 doubles the relative effect for every species
% etc.

%__________________________________________________________________________________________
%
%       SET UP CORAL BLEACHING
%__________________________________________________________________________________________
META.doing_bleaching = 0 ; % set to 0/1 to switch off/on bleaching

CORAL.sensitivity_bleaching = [10 ; 10 ; 5 ; 10 ]; % relative sensitivity of coral species to bleaching
CORAL.bleaching_partial_reduction = [0.5 ; 0.3 ; 0.1 ; 0.5]; % extent of area lost (default is 70% of area lost = same as before)
% A multiplicator of the bleaching mortality probabilities (whole and partial colonies).
% Reduction of the size of a colony affected by partial mortality due to bleaching
% If 0.7, must be read as "the colony looses 70% of its size". Note there is a confusion here with Edwards et al. (2011) which says:
%"Partial mortality corresponded to colonies being reduced in area by 30%." whereas in the orginal code of bleaching_effect_HJE
% this is implemented as a loss of 70% of the pre-bleaching size.

CORAL.bleaching_mortality_old_corals = 0.6 ; % Using Porites and Favia from Van Woesik data
% Below are parameters not currently in use - probably left-over from previous versions (before HJE changes on bleaching?)
% CORAL.bleaching_partial_on_new_brooders = 0.0825 ;  % Based on Agaricia and Porites in McField spreadsheet
% CORAL.bleaching_partial_on_new_spawners = 0.0738 ;  % Based on Montastraea and Siderastrea in McField
REEF.past_bleached_prop = 0.5;


CORAL.clade_prop = 0.95 ; % proportion of clade thermally sensitive over tolerant at initial step
CORAL.clade_reduced_growth = 0.5 ; % proportional reduction of coral growth rate when clade 2

CORAL.proba_switching = [0.5 ; 0.5 ; 0.5 ; 0.5]; % probability of a switching clade from sensitive to tolerant after a bleaching event
% NOTE THIS IS INDEPENDENT OF BLEACHING MORTALITY

CORAL.bleaching_tolerance_clade = 0.29 ; % Reduces mortality due to bleaching if clade 2
% NOTE this now replaces the reduction in mortality due to past bleaching experience


%__________________________________________________________________________________________
%
%       SET UP 3D
%_________________________________________________________________________________________

META.doing_3D = 0 ; % set to 0/1 to switch off/on coral 'calcification'

REEF.initial_rugosity = 1.5 ; % Rugosity of the non-living substrate for initialisation
% Actual rugosity will be higher after adding living coral cover, so that several trials are
% be necessary to obtain the desired reef rugosity. 
% Reef dominated by Montastraea: R~1.2 for coral cover=5% (Alvarez-Filip et al. 2011),
% then add 0.2 units for every 5 additional percent of coral cover

REEF.min_height = 3; % minimum height of dead coral colonies after erosion (when below, colony is removed)
REEF.min_diameter = 2 ; % minimum diameter of dead coral colonies (when below, colony is ignored, 
% and exposes the underlying substrate)

%__________________________________________________________________________________________
%
%       SET UP RANDOMIZATION
%_________________________________________________________________________________________

META.doing_randomization = 0 ; % set to 1 to randomize input parameters (see f_randomize_inputs.m)
% If ==1, the standard deviations velow must be valued

% Standard deviation of initial coral cover for randomization
CORAL.initial_cover_SD2 = [0 ; 0 ; 0 ; 0]; % null because no randomization by default -> needs parametrisation

REEF.nongrazable_substratum_SD2 = 0 ; % standard deviation of the non-grazable cover

CORAL.growth_rate_SD2 = [0 ; 0 ; 0 ; 0]; % Standard deviation of coral growth rate

REEF.herbivory_SD2 = 0 ; % Standard deviation of hebivory

REEF.fish_bioerosion_SD2 = 0; % for randomization; null by default

%__________________________________________________________________________________________
%
%       SET UP CALCIFICATION * (OLD WAY FOR THE CARIBBEAN)
%_________________________________________________________________________________________
% * For the moment only estimates changes in coral calcification due to SST rises
% Parameters for carbonate budget may be implemented later

META.doing_calcification = 0 ; % set to 0/1 to switch off/on coral 'calcification'

%__________________________________________________________________________________________
%
%       SET UP BIOEROSION
%_________________________________________________________________________________________
% Bioeorsion on dead substratum (mixing past and recently dead)
REEF.fish_bioerosion = 0.1989*9.33*365/(2*10000); % in cm3/cm2/6 month
% 0.1989 cm3/m2/hour obtained from  extractdataforgrazers(SWC_unfished_dat,SWC_unfished_txt)
% 9.33 is the number of daylight hours (see Pete_bioerosion2.m) and /10000 to convert in cm2
% -> gives 0.0339 cm3/cm2/6mo
% Mallela & Perry (2007): 20 g CaCO3/m2/yr ->  20/(1.5*10000) = 0.0013 cm3/cm2/year
 
REEF.sponge_bioerosion = 0.507*1000/(10000*2*1.5);
% 0.507 kg/m2/year which is erosion rate on dead colonies (see Pete_bioerosion2.m)
% 1.5 is coral density (g/cm3) (average from Pete_bioerosion2.m). Use the species specific densities?
% 1000 converts kg to g
% /10000 to convert in cm2
%NOTE this produces bioerosion rate for fish much higher than for sponges(3 times higher)

%__________________________________________________________________________________________
%
%       SET UP SEDIMENTATION - implementation to be checked (desactivated until further notice)
%__________________________________________________________________________________________

REEF.doing_sedimentation = 0 ; % set to 0/1 to switch off/on CHRONIC sedimentation effects

REEF.is_sedimentation_seasonal = 0 ; % set to 0/1 to switch off/on sedimentation seasonally in summer

REEF.proba_sedimentation = 0 ; % use for a stochastic probability of a sedimentation event

% Probability of recruitment for each coral under sedimentation (0.2*CORAL.recruitment_proba)
recruit_prob_sed = [ 0.01 ; 0.01 ; 0.01 ; 0.01];

% Coral growth rates under sedimentation
growth_rate_sed = [ 0.50 ; 0.50 ; 0.30 ; 0.20];

% Intercept for incidence relationship for sedimentation
partial_mortality_inci_int_sed = 60 ; % NEEDS A NEW PARAMETERISATION!

partial_mortality_inci_gra_sed = -12 ;

partial_mortality_area_int_sed = -0.5 ;

partial_mortality_area_gra_sed = 1.1 ;

%__________________________________________________________________________________________
%
%       SET UP CORAL DISEASES - not implemented yet
%__________________________________________________________________________________________

META.doing_disease = 0 ; % set to 0/1 to switch off/on diseases

% specific parametrisation has to be done in settings_DISEASES.m

%_________________________________________________________________________________________
%
%       LOADING FILES
%_________________________________________________________________________________________

% Load algal growth parameters
% NOTE: declared as 16-bit unsigned integers for light storage (range from 0 to 65535)
load('data/dict_and_lob_for_model.mat')
META.lobophora_dynamics = lob_dynamics; 
META.dictyota_dynamics = dict_dynamics;
clear lob_dynamics dict_dynamics

% Set up algal growth under sedimentation stress (NOT GENERATED YET)
% load('data/s_dict_and_lob_for_model.mat'); % STILL NEED TO CREATE THESE DATASETS
