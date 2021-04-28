%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 20/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Populate the reef grid randomly with coral colonies of random sizes
% Number of colonies within a cell for each species is random and bounded by META.max_colonies
% Now integrates clade C/D, which replaces the bleaching history of each colony

function [coral, algal] = f_initialise_population(META, REEF, CORAL, algal_initial_cover)    

% Dimensions of the grid
m = META.grid_x_count ;
n = META.grid_y_count ;

% Initialise metrics of coral colonies
coral(META.nb_coral_types).cover_cm2 = sparse(zeros(m*n, 1)) ; % planar area in cm2 (~size)

% Initialise the cover of each algal type
algal(META.nb_algal_types).cover_cm2 = sparse(zeros(m*n, 1)) ;

% Initialise temporary variables
coral_cm2 = zeros(m*n, META.nb_coral_types,META.max_colonies) ;
coral_ID = coral_cm2 ;
algal_cm2 = zeros(m*n, META.nb_algal_types) ;

% List of ID of grazable cells (exclude sand in loops)
list_cell = uint16(1:(m*n)) ;
list_cell = list_cell(REEF.grazable_cell==1) ;

%%%%% ADD CORALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If initial coral cover is non-null, lets fill the grid with coral colonies
if sum([CORAL.initial_cover])~=0
    
    all_colony_sizes=[];
    all_colony_types=[];
    all_colony_ID=[];
    
    for type = 1:META.nb_coral_types
        
        % PULL IN FREQUENCY FOR EACH SIZE CLASS
        colony_size = f_derivecoralcover(type, REEF, CORAL); % generate colonies of different sizes (cm2)

        % STORE THE CORRESPONDING CORAL TYPE
        colony_type = type*ones(size(colony_size));
        
        % CREATE an ID for each colony
        colony_ID = 1:1:length(colony_type);

        % CONCATENATE
        all_colony_sizes = [all_colony_sizes colony_size]; % list of sizes of all colonies, all species combined
        all_colony_types = [all_colony_types colony_type]; % corresponding coral type
        all_colony_ID = [all_colony_ID colony_ID]; %corresponding ID (to be combined later with colony type for full identification)

    end

    [all_colony_sizes,rank] = sort(all_colony_sizes,'descend'); % sort colony size by decreasing value (largest colonies get settled first)
    all_colony_types = all_colony_types(rank); % sort colony type accordingly
    all_colony_ID = all_colony_ID(rank); % sort colony type accordingly
      
    % THEN FILL THE GRAZABLE CELLS WITH CORAL COLONIES (largest get settled first)  
    i = 0; % iterator for visiting the cells
    r = randperm(length(list_cell)) ;
    list_cell2 = list_cell(r) ;
       
    while length(list_cell2)<100000       
        list_cell2 = [list_cell2, list_cell2, list_cell2, list_cell2, list_cell2]; % duplicate the random list of cells        
    end

    current_total_cover = zeros(m*n,1); % initialise the counter of total coral cover per cell
    current_total_colonies = zeros(m*n,META.nb_coral_types); % counter of number of colonies per cell
    
    for c = 1:length(all_colony_sizes) % pick up every colony, one by one
    
        accept = 1;  % decide whether the colony has been processed or not

        while accept == 1           

            i = i+1; % iterator of cell visiting -> goes to the next cell until accept = 2           
            cell = list_cell2(i); 
            colony_count = current_total_colonies(cell,all_colony_types(c))+1;
                           
            if (colony_count <= META.max_colonies) ...
                    && ((all_colony_sizes(c)+current_total_cover(cell,1))<=REEF.substrate_SA_cm2(cell))

                coral_cm2(cell, all_colony_types(c), colony_count) = all_colony_sizes(c) ; % allocate it to the cell
                coral_ID(cell, all_colony_types(c), colony_count) = all_colony_ID(c) ;
                current_total_cover(cell,1) = current_total_cover(cell,1) + all_colony_sizes(c) ; %update total coral cover in that cell
                current_total_colonies(cell,all_colony_types(c)) = colony_count;
                accept = 2; % process the next colony
   
            end
        end
    end
end


%%%%% ADD MACROALGAE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note this works in such a way that after initialisation all cells can only contain one type of algae
% and where macroaglae is present, together with the coral cover the cell is 100% covered
% only the last cell allocated macroalgae for each macroalgal type will leave any space for turf

% total amount of corals to restrict allocation of macroalgae
total_coral_cm2 = round(sum(sum(coral_cm2,3),2));

for a = 2:META.nb_algal_types % process Dictyota, then Lobophora and Halimeda (Turf is processed at the end)
    
    % total amount of macroalgae to allocate
     total_macroalgal_cm2 = algal_initial_cover(a)*sum(REEF.substrate_SA_cm2);

    i = 0 ;  
    r = randperm(length(list_cell)) ;
    list_cell2 = list_cell(r) ;

    while total_macroalgal_cm2 > 0
        
        i = i+1;   
        cell = list_cell2(i) ;

        % Estimate the current available space that can be filled with that macroalgae
        free_space = REEF.substrate_SA_cm2(cell) - total_coral_cm2(cell) - sum(algal_cm2(cell,:)) ;
        % NOTE: summing over all algal types is OK because turf is not defined yet
        
        if  free_space > 0 % if there is some space for the current macroalgae
            
            if total_macroalgal_cm2 - free_space > 0
                % fill that space with macroalgae
                algal_cm2(cell,a) = free_space ;
                % update the total of macroalgae to allocate
                total_macroalgal_cm2 = total_macroalgal_cm2 - free_space ;
                
            else % more space available than macroalgae left
                algal_cm2(cell,a) = total_macroalgal_cm2;
                break
                
            end
        end
    end
end

%%%%% FINAL ALLOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:META.nb_coral_types
    
    coral(s).cover_cm2 = squeeze(coral_cm2(:,s,:)) ;
    coral(s).colony_ID = squeeze(coral_ID(:,s,:)) ;
    
    % Assign a clade (affects mortality to thermal stress and growth rate)
    coral(s).clade = spones(coral(s).cover_cm2); %default is 1 (thermally sensitive, clade C)
    rand_clade = sprand(coral(s).cover_cm2) ;
    coral(s).clade(rand_clade > CORAL.clade_prop) = 2 ; % alternative is 2 (thermally tolerant, clade D)
    
    if META.doing_3D==1    %Initialise 3D variables
        
        coral(s).surface_cm2 = zeros(size(coral(s).cover_cm2)) ;
        coral(s).volume_cm3 = coral(s).surface_cm2 ;
        
    else
        coral(s).surface_cm2 = 0 ;
        coral(s).volume_cm3 = 0 ;
        
    end
end

for a = 2:META.nb_algal_types
    
    tmp_a = zeros(m*n,1) ; 
    tmp_a(:) = round(algal_cm2(:,a)) ;
    algal(a).cover_cm2 = sparse(tmp_a) ;
    
end

% Fill the rest of each cell with turf
algal(1).cover_cm2 = REEF.substrate_SA_cm2 ...
	- sum([algal(2:META. nb_algal_types).cover_cm2],2) - sum([coral.cover_cm2],2) ;
algal(1).cover_cm2(REEF.grazable_cell==0) = 0 ;
algal(1).cover_cm2 = sparse(algal(1).cover_cm2);

% Final check
total = sum(cat(2,algal.cover_cm2),2) + sum(cat(2,coral.cover_cm2),2);
total(REEF.grazable_cell==0) = REEF.substrate_SA_cm2(REEF.grazable_cell==0);

if sum(total) ~= sum(REEF.substrate_SA_cm2)
    error('inconsistent filling across the grid')
end