%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coral, algal, REEF] = f_initialise_rugosity(META, REEF, CORAL, coral, algal, algal_initial_cover)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First set up the desired surface area of the reef substrate (dead corals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate colony 3D colony features (height, surface, volume) before turning the existing corals into dead ones
do_erosion = 0;
last_surface_area_grazed = zeros(size(REEF.grazable_cell));
[coral,volume_dead_eroded_cm3] = f_estimate_3D_colony(coral, do_erosion, REEF, META, last_surface_area_grazed);
temp_coral = coral ; %store the live coral cover for temporary use

% Kill everything
for s=1:META.nb_coral_types
    coral(s).cover_cm2 = - coral(s).cover_cm2 ;
end

 % Then estimate 3D reef features (substrate area, floor area)
[coral, algal,  REEF]= f_estimate_3D_reef(coral, algal, META.cell_area_cm2, REEF,META);
coral_sp_pct2D = zeros(META.nb_coral_types,1);
[rugosity, SI_reef] = f_estimate_rugosity(REEF, META, CORAL, coral_sp_pct2D) ;

% Repeat for increasing substrate SA until it reaches the minimum desired
list_cell = uint32(1:(META.grid_y_count*META.grid_x_count)) ;
list_cell = list_cell(REEF.grazable_cell==1) ;

i = 0; % counter to quit the while loop 

while rugosity < REEF.initial_rugosity

    i = i+1;
    
    if i==50 % stop after 50 iterations
        % another way is to control i in the while command so we give up if we can't reach initial
        % rugosity: while rugosity_1D < REEF.initial_rugosity && i<50
        % but then there is no error message
        rugosity
        error('Initial rugosity cannot be reached -> consider decreasing it or increase coral cover')    
    end
    
    for s=1:META.nb_coral_types
        r = uint16(randperm(length(list_cell))) ; % randomize the order of cell to visit
        coral(s).cover_cm2(list_cell,:) = coral(s).cover_cm2(list_cell(r),:) ;
        coral(s).surface_cm2(list_cell,:) = coral(s).surface_cm2(list_cell(r),:) ;
        coral(s).volume_cm3(list_cell,:) = coral(s).volume_cm3(list_cell(r),:) ;
        coral(s).clade(list_cell,:) = coral(s).clade(list_cell(r),:) ;
        
        coral(s).cover_cm2 = [-abs(coral(s).cover_cm2) -temp_coral(s).cover_cm2] ;
        coral(s).surface_cm2 = [coral(s).surface_cm2 temp_coral(s).surface_cm2] ;
        coral(s).volume_cm3 = [coral(s).volume_cm3 temp_coral(s).volume_cm3] ;
        coral(s).clade = [coral(s).clade temp_coral(s).clade] ;
    end

    % then estimate the reef surface area
    [coral, algal,  REEF]= f_estimate_3D_reef(coral, algal, META.cell_area_cm2, REEF,META);
    [rugosity, SI_reef] = f_estimate_rugosity(REEF, META, CORAL, coral_sp_pct2D);
%     rugosity_2D = sum(REEF.floor_SA_cm2)/META.total_area_cm2;
%     rugosity_1D = 0.88*rugosity_2D;
end

% Because dead colonies where added in blocks (sum of all deads = initial cover)
% we may have overtaken the desired rugosity (defined by REEF.initial_rugosity)
% so that it is necessary to erase some dead colonies to get down rugosity.
% First thing is to define a size threshold for removing dead colonies.
% This assumes that smaller dead colonies are less represented in the substrate
% Actuallty we don't know if the size distribution of the deads is similar to that of the living
% corals, but it is reasonable to assume there are less small colonies (because eroded rapidly)

cover_max = 3000;

 while rugosity > REEF.initial_rugosity   % as long as we are above the defined rugosity

    cover_max = cover_max - 50 ;
    for s=1:META.nb_coral_types

        id_large = find(abs(coral(s).cover_cm2)> cover_max);
        coral(s).cover_cm2(id_large) = 0 ;
        coral(s).surface_cm2(id_large) = 0 ;
        coral(s).volume_cm3(id_large) = 0 ;
        coral(s).clade(id_large) = 0 ;
    end
    
    % then estimate the reef surface area
    [coral, algal,  REEF]= f_estimate_3D_reef(coral, algal, META.cell_area_cm2, REEF,META);
    
    [rugosity, SI_reef] = f_estimate_rugosity(REEF, META, CORAL, coral_sp_pct2D);
%     rugosity_2D = sum(REEF.floor_SA_cm2)/META.total_area_cm2;
%     rugosity_1D = 0.88*rugosity_2D;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then perform the final initialization of live corals on the new substrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[new_coral, algal] = f_initialise_population(META, REEF, CORAL, algal_initial_cover) ;

% Estimate surface area for the new live corals
do_erosion=0;
[new_coral,volume_dead_eroded_cm3] = f_estimate_3D_colony(new_coral, do_erosion, REEF, META, last_surface_area_grazed);

for s=1:META.nb_coral_types

    coral(s).cover_cm2 = [-abs(coral(s).cover_cm2) new_coral(s).cover_cm2] ;
    coral(s).surface_cm2 = [coral(s).surface_cm2 new_coral(s).surface_cm2] ;
    coral(s).volume_cm3 = [coral(s).volume_cm3 new_coral(s).volume_cm3] ;
    coral(s).clade = [coral(s).clade new_coral(s).clade] ;
end

% then estimate reef surface area
[coral, algal,  REEF]= f_estimate_3D_reef(coral, algal, META.cell_area_cm2, REEF,META);
