%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 17/12/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REEFMOD RUN FILE
%
% This function runs 1 simulation of ReefMod, either for a single reef or multiple reefs
% The number of reefs is specified in META.nb_reefs (default is 1).
%
% Coral and algal covers for reef #n are stored respectively in the struct arrays:
% - metapop(n).coral
% - metapop(n).algal
% where ~.coral is itself a struct array that gives for every coral species 's':
%       - ~.coral(s).cover_cm2 = the 2D planimetric area (coral size) in cm2 of every colony in every cell
%               -> positive values for live corals
%               -> negative values for dead corals
%       - ~.coral(s).surface_cm2 = the 3D surface area in cm2 (paraboloid) of every colony in every cell
%       - ~.coral(s).volume_cm3 = the volume in cm3 (paraboloid) of every colony in every cell
%       - ~.coral(s).clade = clade 1 (thermally sensitive) or 2 (tolerant)
% and where ~.algal is itself a struct array that gives for every algal type 'a':
%       - ~.algal(a).cover_cm2 = the size (planimetric area in cm2) of their patch in every cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RESULT, RECORD] = f_runmodel(META, REEF, CORAL, ALGAL)

%___________________________________________________________________________________________________
%
%       SPACE ALLOCATION
%___________________________________________________________________________________________________

metapop(META.nb_reefs).coral = {} ;
metapop(META.nb_reefs).algal = {} ;

% matrices of bleaching mortality (reef specific because linked to DHW)
bleaching_mortality(META.nb_reefs).whole_colony ={} ;
bleaching_mortality(META.nb_reefs).partial_on_brooders ={} ;
bleaching_mortality(META.nb_reefs).partial_on_spawners ={} ;

% This records coral and algal cover (%) over the grid for each species at every time steps for
% every reef. Thus colony cover is not tracked over time, just population cover.
RESULT.coral_pct2D = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.coral_pct3D = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.algal_pct = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_algal_types);
RESULT.clade_prop = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.nb_settlers = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.total_fecundity = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.input_larvae = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);

if META.doing_3D==1
    RESULT.live_CaCO3 = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    RESULT.dead_CaCO3 = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    RESULT.substrate_cm2 = zeros(META.nb_reefs, META.nb_time_steps+1);
end

if META.do_size_frequency==1
    RESULT.juv_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(0:META.size_bins(1):META.juvenile_max_size));
    RESULT.adol_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(0:META.size_bins(2):META.adult_size));
    RESULT.adult_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(0:META.size_bins(3):floor(pi*(META.cell_x_size/2)^2)));
end

% Recording mortality events
RECORD.bleaching_events = zeros(META.nb_reefs, META.nb_time_steps) ; 
RECORD.hurricane_events = zeros(META.nb_reefs, META.nb_time_steps) ;
% NOTE could be done for predation, hurricane etc. -> see past versions for attempt of implementation

SEDIMENTATION = zeros(META.nb_reefs, META.nb_time_steps);

% Memory optimization (just extract connectivity matrices from META and store it temporally)
if META.use_connectivity_matrix==1
    connectivity.matrix =[];
    % For connectivity stuff
    FECUNDITY_ADOL = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    FECUNDITY_ADULT = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    INPUT_LARVAE = zeros(META.nb_reefs,META.nb_coral_types);
    for s=1:META.nb_coral_types
        connectivity(s).matrix = META.connectivity(s).matrix; %just rename to ease calculations
    end
end
META = rmfield(META, 'connectivity');

ID_colony_tracking = zeros(1,META.nb_coral_types);
colony_list = zeros(1,6);
environ_list = zeros(1,5);

%______________________________________________________________________________________________________________________________________________
%
% INITIALISATION
%______________________________________________________________________________________________________________________________________________
% disp('Initialising reef states...')

nb_cells = META.grid_x_count * META.grid_y_count ;
t = 0;

for n = 1:META.nb_reefs % This must be done for every reef before time simulations
    
    %%%%%% RE-INITIALISE GROWTH RATE IN CASE IT PREVIOUSLY CHANGED %%%%%%
    if META.doing_calcification == 1
        % Re-initialise coral growth rates because they may have changed
        CORAL(n).growth_rate = init_coral_growth_rate ;
    end
    
    %%%%%% DEFINE UNGRAZABLE CELLS FIRST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    list_cell = uint32(1:nb_cells) ;
    REEF(n).grazable_cell = ones(nb_cells,1) ; % create a list of grazable cells (grazable==1)
    id = randsample(list_cell, uint32(REEF(n).nongrazable_substratum * nb_cells)) ; % select randomly cell indices
    REEF(n).grazable_cell(id) = 0 ; % switch to 0 the randomly selected cells (nongrazable==0)
    
    % Define the surface area (SA) of the reef substrate of each grazable cell
    % By default, every cell is flat. If doing 3D, the substrate SA will change over time,
    % otherwise will stay flat. Non-grazable cells (sand) are flat (SA = cell area)
    REEF(n).substrate_SA_cm2 = ones(nb_cells,1)*META.cell_area_cm2;
    REEF(n).floor_SA_cm2 = ones(nb_cells,1)*META.cell_area_cm2;
  
    %%%%%% INITIALISE THE GRID WITH CORALS AND ALGAE %%%%%%%%%%%%%%%%%%%%%%
    [metapop(n).coral, metapop(n).algal] = ...
        f_initialise_population(META, REEF(n), CORAL(n), ALGAL(n).initial_cover) ;
    
    [metapop(n).coral] = f_struct_arrange(metapop(n).coral, META.nb_coral_types, META.doing_3D);
        
    %%%%%% SET UP 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if META.doing_3D == 1
        % Estimate surface, volume of living and dead (after erosion) colonies.
        % Produces also different surface areas for each cell in REEF that will vary over time:
        % - surface area of the substrate underneath live corals (REEF.substrate_SA_cm2)
        % - surface area of the reef floor, ie including living colonies (REEF.floor_SA_cm2)
        [metapop(n).coral, metapop(n).algal, REEF] = f_initialise_rugosity(META, REEF(n), CORAL(n), metapop(n).coral, metapop(n).algal, ALGAL(n).initial_cover);
        
        RESULT.volume_eroded(n,1) = 0 ;
        RESULT.substrate_cm2(n,1) = sum(REEF(n).substrate_SA_cm2); % total 3D area of the substrate
        RESULT.floor_cm2(n,1) = sum(REEF(n).floor_SA_cm2); % total 3D area of the sea floor (including coral 3D surface)
    end
    
    %%%%%% SET UP HURRICANES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If randomized with the chosen regime (makes the chronology starting at different dates)
    if META.doing_hurricanes ==1 && META.randomize_hurricane_chronology == 1
        
        [c]=size(REEF(n).hurricane_chronology,2); % c is the number of time steps (seasons) of the chronology
        date = randsample(c,1);% Pick up a date randomly for starting the new chronology
        % Cut the chronology in two pieces at the chosen date, and collate the left piece after the right piece
        REEF(n).hurricane_chronology = [REEF(n).hurricane_chronology(:,date:c) REEF(n).hurricane_chronology(:,1:c-1)] ;
        
    end
    
    %%%%%% SET UP BLEACHING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if META.doing_bleaching == 1
        
        % Generate mortality probability due to bleaching:
        bleaching_mortality(n).whole_colony = Inf;

        while sum(bleaching_mortality(n).whole_colony)==Inf % Because the function below produces Inf values sometimes            
            bleaching_mortality(n).whole_colony = f_generate_bleaching_mortalities(REEF(n).predicted_DHWs) ;
        end
        
        bleaching_mortality(n).partial_on_brooders = 0.0027*REEF(n).predicted_DHWs + 0.0717 ;
        bleaching_mortality(n).partial_on_spawners = 0.0027*REEF(n).predicted_DHWs + 0.0630 ; 
        
    end

    %%%%%% STORE INITIAL COVERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) Algal cover: total cover (%) of each algal type over the grid at initial step
    RESULT.algal_pct(n,1,1:META.nb_algal_types) = 100*full(sum([metapop(n).algal(1:META.nb_algal_types).cover_cm2]))./sum(REEF(n).substrate_SA_cm2) ;
 
%     [environ_list] = f_record_environ(environ_list,metapop(n).algal,t,META,REEF(n).grazable_cell);

    % 2) Coral cover: total cover (%) of each coral type at initial step
    for s=1:META.nb_coral_types
        
        living_planar_coral_cover_cm2 = metapop(n).coral(s).cover_cm2(metapop(n).coral(s).cover_cm2>0) ;
        RESULT.coral_pct2D(n,1,s) = 100*sum(living_planar_coral_cover_cm2)./sum(REEF(n).substrate_SA_cm2) ;

        living_coral_clades = metapop(n).coral(s).clade(metapop(n).coral(s).cover_cm2>0) ;
        RESULT.clade_prop(n,1,s) = nnz(living_coral_clades==1)/nnz(living_coral_clades);
        
        % Estimaton of living and dead carbonate volumes (per species)
        if META.doing_3D == 1
            living_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2>0) ;
            dead_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2<0) ;
            RESULT.live_CaCO3(n,1,s)= sum(living_coral_volume_cm3);
            RESULT.dead_CaCO3(n,1,s)= sum(dead_coral_volume_cm3);
        end
        
        % Coral size structure at initial step
        if META.do_size_frequency == 1
            [count, class] = f_count_sizefreq(living_planar_coral_cover_cm2, META) ;
            RESULT.juv_count(n,1,s,:)=count.juv;
            RESULT.adol_count(n,1,s,:)=count.adol;
            RESULT.adult_count(n,1,s,:)=count.adult;
        end
        
        % Record the ID max
        if sum(sum(metapop(n).coral(s).cover_cm2))~=0            
            ID_colony_tracking(1,s) = max(max(metapop(n).coral(s).colony_ID));
        end

    end
    
    % 3) Estimate reef rugosity    
    if META.doing_3D == 1      
        coral_sp_pct2D = squeeze(squeeze(RESULT.coral_pct2D(n,1,:)));
        [rugosity, SI_reef] = f_estimate_rugosity(REEF(n), META, CORAL(n), coral_sp_pct2D);
        RESULT.rugosity(n,1)=rugosity;
        RESULT.SI_reef(n,1)=SI_reef;
    end
    
    if META.track_populations == 1
        [environ_list,colony_list] = f_track_populations(environ_list,colony_list,metapop(n).algal,metapop(n).coral,t,META,REEF(n).grazable_cell);
    end
    
    %%%% If doing connectivity, then estimate coral fecundity of the reef
    %%%% --------------------------------------------------------------
    if META.use_connectivity_matrix == 1
        
        for s=1:META.nb_coral_types
            [FECUNDITY_ADOL(n, 1, s), FECUNDITY_ADULT(n, 1, s)] = f_estimate_fecundity (metapop(n).coral(s).cover_cm2, META);
            % = total number of larvae produced per species for each reef at every time step
        end
        
        RESULT.total_fecundity(n,1,:) = FECUNDITY_ADOL(n, 1, :) + FECUNDITY_ADULT(n, 1, :) ;
    end
    
end

if META.use_connectivity_matrix == 1
    MAX_LARVAL_POOL = zeros(1, META.nb_coral_types);
    
    for s=1:META.nb_coral_types
        MAX_LARVAL_POOL(1,s) = median((connectivity(s).matrix * RESULT.total_fecundity(:,1,s)));
    end
end
%______________________________________________________________________________________________________________________________________________
%
% RUN TIME SIMULATIONS
%______________________________________________________________________________________________________________________________________________
% disp('Running time simulations...')

arrange_step = zeros(1,META.nb_time_steps);
arrange_step(1:3:META.nb_time_steps)=1;
arrange_step(1)=0;
eaten_coral_cm2 = zeros(META.nb_reefs,META.nb_coral_types);

for t = 1:META.nb_time_steps
    
    %%%%%%%%% ADD SEASONALITY FOR DICTYOTA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    season = iseven(t); % 1 is winter, 0 is summer (first time step is summer)

    %%%%%%%%% ESTIMATE LARVAL POOL ARRIVING AT EVERY REEF %%%%%%%%%%%%%%%%%
    if META.nb_reefs > 1 && META.use_connectivity_matrix == 1
        
        for s=1:META.nb_coral_types %
            % First weigh output larvae by proportional suitable habitat for Montastraea
            output_larvae = RESULT.total_fecundity(:,t,s).*META.connectivity_area_habitat;
            % Now larvae are distributed over the connected reefs
            RESULT.input_larvae(:,t,s) = transp(output_larvae)*connectivity(s).matrix ;
            % Matrix multiplication1 (Connectivity matrices: rows = sink reefs, columns = source)
        end
    end
    
    %%%%%%%%% PROCESS EVERY REEF ONE AFTER THE OTHER  %%%%%%%%%%%%%%%%%%%%%
    for n = 1:META.nb_reefs
               
        %%%% --------------------------------------------------------------
        %%%% Recruitment 
        %%%% --------------------------------------------------------------
        
        if META.nb_reefs > 1 && META.use_connectivity_matrix == 1 % --> CONNECTIVITY-based recruitment
            % Max density of settlers is weighted by loca larval input, proportional to maximum larval input
            max_density_settlers = META.max_density_settlers.*transp(squeeze(RESULT.input_larvae(n,t,:)))./MAX_LARVAL_POOL ;
            [metapop(n).coral, total_settled, ID_colony_tracking, metapop(n).algal] = f_apply_recruitment_new(metapop(n).coral, metapop(n).algal,...
                META,REEF(n).grazable_cell, REEF(n).floor_SA_cm2, max_density_settlers, CORAL(n).clade_prop,ID_colony_tracking);
            
        else   % --> FIXED probability of recruitment = MOST COMMON SCENARIO FOR A SINGLE REEF [SEE OLDER SCRIPTS FOR OTHERS]          
            max_density_settlers = META.max_density_settlers(:,t); %
            [metapop(n).coral, total_settled, ID_colony_tracking, metapop(n).algal] = f_apply_recruitment_new(metapop(n).coral, metapop(n).algal,...
                META, REEF(n).grazable_cell, REEF(n).floor_SA_cm2, max_density_settlers, CORAL(n).clade_prop,ID_colony_tracking);

        end
        
        RESULT.nb_settlers(n,t,:) = total_settled ; % record the number of settlers before they are processed

        %%%% --------------------------------------------------------------
        %%%% POPULATION DYNAMICS 
        %%%% --------------------------------------------------------------
        
        %%%%%%%%% EFFECTS OF SST ON CORAL GROWTH RATE %%%%%%%%%%%%%%%%%%%%%
        % This has to occur before processing population because will modify coral growth
        % Estimated growth rate at t is directly stored in META but it is OK since initial growth rates have been recorded earlier
        % NOTE: use the pre-specified scenario of growth rates as in Bozec & Mumby (2015)
        if META.doing_calcification == 1          
            % This is the predicted SST for the current time step, in the exact same way Pete_bioerosion2.m does
            current_sst = REEF(n).SST(1,t);  % this is the SST of the current time step
            % sd_sst = 1.47 % useful ?? 
            % Then calculate the relative reduction of calcification rate for each coral species
            Rel_change = exp( - 0.5 *( (current_sst - CORAL(n).SST_OPT) ./ CORAL(n).SD_relative_calci) .^2 );
            % Update the coral growth rates
            CORAL(n).growth_rate = init_coral_growth_rate .* CORAL(n).relative_calci .* Rel_change;   
        end
        

        %%%%%%%%% THEN PROCESS CORAL POPULATION OF REEF N %%%%%%%%%%%%%%%%%   
        [metapop(n).coral, metapop(n).algal, last_surface_area_grazed] = ...
            f_process_population (metapop(n).coral, metapop(n).algal, season, REEF(n).herbivory(1,t), META, REEF(n), CORAL(n), ALGAL(n));
                % first, re-arrange the array of corals to speed-up the code
        
        %%%%%%%%% Re-arrange the coral matrix for optimization (at the chosen time step)
        if arrange_step(t)==1 % only occurs at the selected time steps)
            [metapop(n).coral] = f_struct_arrange(metapop(n).coral, META.nb_coral_types, META.doing_3D);
        end
        

        %%%% --------------------------------------------------------------
        %%%% HURRICANES
        %%%% --------------------------------------------------------------   
        if META.doing_hurricanes == 1 && season == 0
            
            if META.random_hurricanes == 1 && rand(1)< REEF.hurr_strike_proba               
                RECORD.hurricane_events(n,t) = randi(5);  % random category hurricane
%                 RECORD.hurricane_events(n,t)  = 4 ; %TEMP -> limit to category 1 (but very high frequency)
                
            else % then use the hurricane chronology
                RECORD.hurricane_events(n,t) = REEF(n).hurricane_chronology(n,t);
            end
            
            % Estimate hurricane impact with the specified category (includes category 0 = no hurricane -> no effect
            [metapop(n).coral, metapop(n).algal] = f_hurricane_effect...
                (metapop(n).coral, metapop(n).algal, RECORD.hurricane_events(n,t), META, CORAL(n)) ;
        end
        
        %%%% --------------------------------------------------------------
        %%%% BLEACHING
        %%%% --------------------------------------------------------------        
        if META.doing_bleaching == 1 && season == 0
            
            if REEF(n).predicted_DHWs(t)>4 && RECORD.hurricane_events(n,t) == 0 % only bleach if no hurricane
                
                RECORD.bleaching_events(n,t) = 1 ; % holds a 0 (no bleaching) or a 1 (bleaching) 

                [metapop(n).coral, metapop(n).algal] = ...
                    f_bleaching_HJE (metapop(n).coral, metapop(n).algal, bleaching_mortality(n), t, CORAL(n), META);
            end
        end  

        %%%% --------------------------------------------------------------
        %%%% CROWN-OF-THORN STARFISH PREDATION
        %%%% --------------------------------------------------------------        
        % to be applied first because was calculated on the final coral cover at the previous time step
        if META.doing_COTs == 1 
            
             COTs_density = META.total_area_cm2*META.COTs_density ; % to decide later whether should be reef specific
             COTs_feeding_rate = META.COTs_feeding_rate(season+1); % pick up the seasonal feeding rate
            [metapop(n).coral, metapop(n).algal] = ...
                f_COTs_predation(metapop(n).coral, metapop(n).algal, COTs_density, COTs_feeding_rate, META.COTs_feeding_prefs) ;       
        end
        
        %%%% --------------------------------------------------------------
        %%%% Indeterminate acute disturbance for Moorea (COTS/cyclone)
        %%%% --------------------------------------------------------------        
        % Bleaching has been negligible in the past
        
        if META.doing_acute_disturbances == 1 && rand(1)< REEF.disturbance_event_proba
            
            [metapop(n).coral, metapop(n).algal] = ...
                f_acute_disturbance(metapop(n).coral, metapop(n).algal, REEF.mortality_proba(:,1+iseven(randi(2)))', META);           
            
        end
        
        %%%% --------------------------------------------------------------
        %%%% FINAL ESTIMATES FOR REEF N
        %%%% --------------------------------------------------------------      
        
        %%%% If doing connectivity, then estimate coral fecundity of the reef
        %%%% --------------------------------------------------------------
        if META.use_connectivity_matrix == 1
            
            for s=1:META.nb_coral_types % this is to be used for the next time step, so t+1
                [FECUNDITY_ADOL(n, t+1, s), FECUNDITY_ADULT(n, t+1, s)] = f_estimate_fecundity (metapop(n).coral(s).cover_cm2, META);
                % = total number of larvae produced per species for each reef at every time step
            end
            
            RESULT.total_fecundity(n,t+1,:) = FECUNDITY_ADOL(n, t+1, :) + FECUNDITY_ADULT(n, t+1, :) ;
        end         
        
        %%%% Estimate 3D areas and RUGOSITY
        %%%% --------------------------------------------------------------
        if META.doing_3D == 1
            do_erosion=1;
            [metapop(n).coral,volume_dead_eroded_cm3] = f_estimate_3D_colony(metapop(n).coral, do_erosion, REEF(n), META, last_surface_area_grazed);
            [metapop(n).coral, metapop(n).algal, REEF(n)]= f_estimate_3D_reef(metapop(n).coral, metapop(n).algal, META.cell_area_cm2, REEF(n),META);
            
            RESULT.volume_eroded(n,t+1) = sum(volume_dead_eroded_cm3) ;
            RESULT.substrate_cm2(n,t+1) = sum(REEF(n).substrate_SA_cm2);
        end
        
                
        %%%% Report current reef status
        %%%% --------------------------------------------------------------
        
        % 1) Algal cover: total cover (%) of each algal type over the grid at time step t
        RESULT.algal_pct(n,t+1,:) = 100*full(sum([metapop(n).algal(:).cover_cm2]))./sum(REEF(n).substrate_SA_cm2)  ;

        % 2) Coral cover: size structure of coral populations at t
        for s=1:META.nb_coral_types
            
            living_planar_coral_cover_cm2 = metapop(n).coral(s).cover_cm2(metapop(n).coral(s).cover_cm2>0) ;
            RESULT.coral_pct2D(n,t+1,s) = 100*sum(living_planar_coral_cover_cm2)./sum(REEF(n).substrate_SA_cm2)  ;
            
            living_coral_clades = metapop(n).coral(s).clade(metapop(n).coral(s).cover_cm2>0) ;
            RESULT.clade_prop(n,t+1,s) = nnz(living_coral_clades==1)/nnz(living_coral_clades);
          
            % Estimaton of living and dead carbonate volumes (per species)
            if META.doing_3D == 1
                living_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2>0) ;
                dead_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2<0) ;
                RESULT.live_CaCO3(n,t+1,s)= sum(living_coral_volume_cm3);
                RESULT.dead_CaCO3(n,t+1,s)= sum(dead_coral_volume_cm3);              
            else % if not doing 3D we remove the dead colonies to speed up the code
                metapop(n).coral(s).cover_cm2(metapop(n).coral(s).cover_cm2<0)=0;
            end
            
            % Coral size structure at time step t
            if META.do_size_frequency == 1
                [count, class] = f_count_sizefreq(living_planar_coral_cover_cm2, META) ;
                RESULT.juv_count(n,t+1,s,:)=count.juv;
                RESULT.adol_count(n,t+1,s,:)=count.adol;
                RESULT.adult_count(n,t+1,s,:)=count.adult;
            end
                     
        end        
     
        % 3) Rugosity at time step t
        if META.doing_3D == 1
            coral_sp_pct2D = squeeze(squeeze(RESULT.coral_pct2D(n,t+1,:)));
            [rugosity, SI_reef] = f_estimate_rugosity(REEF(n), META, CORAL(n), coral_sp_pct2D);
            RESULT.rugosity(n,t+1)=rugosity;
            RESULT.SI_reef(n,t+1)=SI_reef;
        end

    end % end of the loop over the n reefs
    
    if META.track_populations == 1
        [environ_list,colony_list] = f_track_populations(environ_list,colony_list,metapop(n).algal,metapop(n).coral,t,META,REEF(n).grazable_cell);      
    end
               
end % end of the time steps loop for a single simulation


if META.track_populations == 1     
    f_generate_track_files(META, REEF, CORAL, ALGAL, RECORD, colony_list, environ_list)     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%