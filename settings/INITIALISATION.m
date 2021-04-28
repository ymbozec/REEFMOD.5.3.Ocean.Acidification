% Y.-M. Bozec, MSEL, created Jan 2012.
% Last modified: 29/08/2012
%
%
% INITIALISATION OF PARAMETERS.
% 
% This performs automatic generation of parameters, 
% such as loading files, basics calculation.
%_________________________________________________________________________________________

%_________________________________________________________________________________________
%
%       MEMORY OPTIMIZATION
%_________________________________________________________________________________________

REEF.dictyota_declines_seasonally=uint8(REEF.dictyota_declines_seasonally);
REEF.doing_sedimentation=uint8(REEF.doing_sedimentation);
REEF.is_sedimentation_seasonal=uint8(REEF.is_sedimentation_seasonal);

CORAL.is_brooder = uint8(CORAL.is_brooder) ;
CORAL.size_model = uint8(CORAL.size_model) ;

%_________________________________________________________________________________________
%
%       DISPLAYING THE SELECTED OPTIONS
%_________________________________________________________________________________________

disp('------------------------------------------')

if META.nb_reefs == 1
    disp('---> SINGLE REEF SIMULATIONS');
    if REEF.dictyota_declines_seasonally==0
        disp('---> OFF seasonality');
    else
        disp('---> ON  seasonality');
    end
    
    if REEF.doing_sedimentation==0
        disp('---> OFF sedimentation');
    else
        disp('---> ON  sedimentation');
        
        if REEF.is_sedimentation_seasonal==0
            disp('---> sedimentation IS NOT seasonal');
        else
            disp('---> sedimentation IS seasonal');
        end
    end
    
else
    disp('---> MULTIPLE REEF SIMULATIONS');
    disp('---> ON  seasonality')
    disp('---> OFF sedimentation');
end


if META.track_populations == 1
    
    disp('---> TRACKING COLONY SIZES');
    
    if META.nb_simul>1
        
        META.nb_simul = 1 ; % forces to only one simulation
        disp(' ########################################################')
        disp(' ## PROGRAM FORCED TO 1 SIMULATION BECAUSE OF TRACKING ##')
        disp(' ########################################################')
        
    end
end

%_________________________________________________________________________________________
%
%       AUTOMATIC GENERATIONS
%_________________________________________________________________________________________

% Number of coral species in the model
META.nb_coral_types = length(CORAL.initial_cover);

% Number of algal species in the model
META.nb_algal_types = length(ALGAL.initial_cover);

% Planar area of a cell in cm2
META.cell_area_cm2 = META.cell_x_size * META.cell_y_size ;

% Planar area of the reef (grid) in cm2
META.total_area_cm2 = META.grid_x_count * META.grid_y_count * META.cell_area_cm2 ; 

% Maximum age of coral colonies in time steps of 6 months so 30 years but constrained by
% the area of the cell - this is now OLD STUFF
% CORAL.max_age = floor(((META.cell_area_cm2/pi)^(1/2))./CORAL.growth_rate);
% Replaced by maximum size (necessary for initialing populations (f_derivecoralcover)
CORAL.max_size = floor(pi*(CORAL.max_diameter/2).^2) ;
CORAL.max_size(CORAL.max_size>pi*(META.cell_x_size/2)^2)=floor(pi*(META.cell_x_size/2)^2); % diameter can't be bigger than a cell

% Generate a scenario of constant herbivory over time if no parrotfish fishing
if META.do_parrotfish_fishing == 0
    REEF.herbivory = REEF.herbivory(1,ones(META.nb_time_steps,1));
end


% Set up sedimenation parameters
if REEF.doing_sedimentation == 1
    % Collate the new rates -> column 1 = no sedim, column 2 = sedim
    CORAL.growth_rate_sed = [CORAL.growth_rate  growth_rate_sed ];
    CORAL.partial_mortality_inci_int = [CORAL.partial_mortality_inci_int  partial_mortality_inci_int_sed];
    CORAL.partial_mortality_inci_gra = [CORAL.partial_mortality_inci_gra  partial_mortality_inci_gra_sed];
    CORAL.partial_mortality_area_int = [CORAL.partial_mortality_area_int  partial_mortality_area_int_sed];
    CORAL.partial_mortality_area_gra = [CORAL.partial_mortality_area_gra  partial_mortality_area_gra_sed];
end
  

for n = 1:META.nb_reefs % This must be done for every reef
    
    REEF(n)=REEF(1);
    CORAL(n)=CORAL(1);
    ALGAL(n)=ALGAL(1);
        % Check consistency in parametrisation
    if META.nb_coral_types ~= length(CORAL(n).initial_cover) || META.nb_coral_types ~= length(CORAL(n).growth_rate) ||...
            META.nb_coral_types ~= length(CORAL(n).max_diameter) ||...
            META.nb_coral_types ~=length(CORAL(n).size_model) || META.nb_coral_types ~= length(CORAL(n).size_mean) ||...
            META.nb_coral_types ~=length(CORAL(n).size_var) || META.nb_coral_types ~= length(CORAL(n).lobophora_reduce_rate) ||...
            META.nb_coral_types ~=length(CORAL(n).dictyota_reduce_rate) || META.nb_coral_types ~=length(CORAL(n).sensitivity_hurricane) ||...
            META.nb_coral_types ~=length(CORAL(n).initial_cover_SD2) || META.nb_coral_types ~=length(CORAL(n).sensitivity_bleaching) ||...
            META.nb_coral_types ~=length(CORAL(n).growth_rate_SD2)  

        error('Error ReefMod: inconsistent number of coral types !!!!');
    end
    
    if sum(CORAL(n).initial_cover)+sum(ALGAL(n).initial_cover)+REEF(n).nongrazable_substratum > 1
        error('Error ReefMod: sum of initial covers is > 1 !!!!')
    end
    
    if ((REEF(n).doing_sedimentation==1) && (REEF(n).is_sedimentation_seasonal==1)) ||...
            ((REEF(n).doing_sedimentation==1) && (REEF(n).proba_sedimentation==1))
        error('too much sedimentation you idiot'); % prevents use of incompatible rules
    end
end

% Create the neighbouring "environment" for every cell (TORUS)
% This identifies the cells which are around every cell, based on the type defined in the function
% f_create_environment (currently 4-cell neighboring for the cellular automaton - but could be
% changed in the future (8-cell neighboring?)
META.environ = f_create_environ (META.grid_x_count, META.grid_x_count) ;

% clear all the temporary variables
clear partial_mortality_inci_int_sed partial_mortality_inci_gra_sed
clear partial_mortality_area_int_sed partial_mortality_area_gra_sed
clear recruit_prob_sed growth_rate_sed

%_________________________________________________________________________________________
%
%       SETTING UP BLEACHING, HURRICANES, CONNECTIVITY, DISEASE
%_________________________________________________________________________________________


%%%% Set up connectivity parameters %%%%%%%%%%%%%%%%%
if META.use_connectivity_matrix==0
    disp('---> OFF connectivity');
    META.connectivity = 0 ;
    
    if META.recruitment_type==0
        disp('---> FULL recruitment');

    else
        disp('---> RANDOM recruitment');
    end
    
else
    settings_CONNECTIVITY;
    
    if META.nb_reefs == 1
        error('REEFMOD error: Define the number of reefs for connectivity (in META.nb_reefs)')
    end
    
    for s=1:META.nb_coral_types
        
        if size(META.connectivity(s).matrix,1) ~= META.nb_reefs
            error('REEFMOD error: connectivity matrices wrong size (must equal the number of reefs!)');
        end
    end
    
    disp('---> ON  connectivity');
end

%%%% Set up bleaching parameters %%%%%%%%%%%%%%%%%%%%
if META.doing_bleaching == 0
    disp('---> OFF  bleaching'); 
else
    settings_BLEACHING;
    disp('---> ON  bleaching');
end

%%%% Set up hurricane parameters %%%%%%%%%%%%%%%%%%%%
if META.doing_hurricanes==0
    disp('---> OFF hurricanes');
else
    settings_HURRICANES;
    disp('---> ON  hurricanes');
end

%%%% Set up 3D parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
if META.doing_3D==0
    disp('---> OFF 3D reef');
else
    settings_3D
    disp('---> ON  3D reef');
end

%%%% Set up disease parameters %%%%%%%%%%%%%%%%%%%%%%
if META.doing_disease == 0
    disp('---> OFF diseases');    
else
    settings_DISEASES
    disp('---> ON  diseases');
end

%%%% Message concenrng calcification
if META.doing_calcification==0
    disp('---> OFF calcification');
else
    disp('---> ON  calcification');
end

%%%% Set up coral-coral competition and acidification
if META.doing_coral_competition==0
    disp('---> OFF coral/coral competition');
else
    settings_CORAL_COMPETITION;
    disp('---> ON  coral/coral competition');
    
end

%%%% Set up randomization of input parameters %%%%%%%%%%%%%%%%%%%%%%
if META.doing_randomization==0
    disp('---> OFF randomize input parameters');
    RANDOM_INPUTS = 0 ;

else
    disp('---> ON  randomize input parameters');

    RANDOM_INPUTS = f_randomize_inputs(META,REEF,CORAL) ;

end

% Display the size of the cell
str_x_cell = num2str(META.cell_x_size);
str_y_cell = num2str(META.cell_y_size);
str_tot_cell = num2str(META.cell_area_cm2/10000);
disp( ['---> CELL SIZE: '  str_x_cell 'cm * '  str_y_cell 'cm' ' = ' str_tot_cell ' m2'] );

% Display the size of the grid
str_x_grid = num2str(META.cell_x_size * META.grid_x_count/100);
str_y_grid = num2str(META.cell_y_size * META.grid_y_count/100);
str_tot_grid = num2str(META.total_area_cm2/10000);
disp( ['---> GRID SIZE: '  str_x_grid 'm * '  str_y_grid 'm' ' = ' str_tot_grid ' m2'] );

clear str_tot_cell str_tot_grid str_x_cell str_x_grid str_y_cell str_y_grid
