%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Oct 2012.
% Last modified: Aug 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   [coral,volume_dead_eroded_cm3]= f_estimate_3D_colony(coral, do_erosion, REEF, META, last_surface_area_grazed)

% Performs erosion of dead corals, then estimate various 3D features for both living and dead colonies.

nb_coral_types = size(coral,2) ;
coral_cm2 = [coral.cover_cm2];
volume_cm3 = [coral.volume_cm3];
id_dead_tot = zeros(size(volume_cm3));
id_dead_tot(coral_cm2<0)=1;

prop = zeros(size(volume_cm3,1),1);

if do_erosion == 1 && sum(last_surface_area_grazed)>0
    
    % Estimate the current reef rugosity to scale up/down parrotfish erosion
%     rugosity_2D = sum(REEF.floor_SA_cm2)/META.total_area_cm2;
%     rugosity_1D = 0.88*rugosity_2D;

    % Then estimate relative erosion similarly to relative grazing
    % This accounts for variations in abundance x grazing rate as a response to changing rugosity
%     a=1.55; b=1.37;
%     relative_erosion = (1/0.6)*a*(rugosity_1D-1)/(1+a*b*(rugosity_1D-1));
    
    % First estimate the total volume of dead CaCO3 to remove from the reef, based on its planimetric area
%     reef_volume_to_erode = relative_erosion*META.total_area_cm2*REEF.fish_bioerosion ; %(1,1)
%     reef_volume_to_erode = META.total_area_cm2*REEF.fish_bioerosion ; %(1,1)
    
%     total_surface_area_grazed = sum(last_surface_area_grazed,1);
%     volume_removed_per_cm2 = reef_volume_to_erode/total_surface_area_grazed; % spread erosion rate over the area grazed

    % Now derive the volume to remove in every cell (only grazable cells)
%     volume_removed_per_cell = volume_removed_per_cm2*last_surface_area_grazed; % then converted to a volume per cell
    volume_removed_per_cell = REEF.fish_bioerosion*last_surface_area_grazed;
    
    volume_dead_tot_cm3 = sum(volume_cm3.*id_dead_tot,2) ;
    id_cells = ones(size(volume_dead_tot_cm3)) ;
    id_cells(volume_dead_tot_cm3==0) = 0 ;

    % this is the proportion of the volume of every dead in every cell that needs to be eroded
    prop(id_cells==1) = volume_removed_per_cell(id_cells==1)./volume_dead_tot_cm3(id_cells==1) ;
    
else
    volume_dead_eroded_cm3 = 0; % Otherwise there is no erosion this time
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proceed with one coral species at a time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Total_V_removed = 0;

for s=1:nb_coral_types 

    % First check if cover_cm2 and volume_cm3 are still concordant
    % (some living colonies may have disappeared so that corresponding volume must be erased)
    id_disappear = find(coral(s).cover_cm2==0 & coral(s).volume_cm3>0);
    coral(s).volume_cm3(id_disappear)=0;
    coral(s).surface_cm2(id_disappear)=0;
  
    % Estimate diameter (cm) for live and dead altogether
    D = floor(2*sqrt(abs(coral(s).cover_cm2)/pi)); 
    % Note that when cover=0cm2 (no colony), D=0 as well; when cover=1cm2 (recruits), D=1cm (approximation)

    % then separate live/dead colonies
    id_live=find(coral(s).cover_cm2>0);
    id_dead=find(coral(s).cover_cm2<0);
    
    %%%%%%%% HEIGHT OF LIVING COLONIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate the height of living corals after they have grown/shrunk in f_process_population

    H=0*D;
    LOOK = META.H_from_D(:,s); % Use species-specific lookup table to estimate H from D
    H(id_live)=LOOK(D(id_live)); % For this reason D must be positive integer

    %%%%%%%% HEIGHT OF DEAD COLONIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bioerosion reduces the volume of dead colonies, then their height assuming a paraboloid

    if do_erosion == 0
        % This is just for passing this step when there is no bioerosion (e.g., at initial step).
        H(id_dead)=LOOK(D(id_dead)); % Not sure it is still necessary though
        
    else % if non-null erosion, then reduces H
        
        % IMPORTANT: some H values can be disproportionately high compared to D.
        % This concerns corals that first shrunk (partial mortality) then died completely (whole mortality). 
        % Because volume is not recalculated for deads, those recently deads kept the volume they had before shrinking,
        % and therefore this volume does not agree with the current diameter. A simple assumption is to consider 
        % that H also shrunk with D to update volumes before erosion
        
        % First estimate H from volume and D
        H(id_dead) = ceil((8/pi)*coral(s).volume_cm3(id_dead)./((D(id_dead).^2))) ;
        
        % Then look for excessive H
        id_H=0*D;
        id_H(D<H)=1; % normally D is always > H, so D<H inidcates those corals that first shrunk then died completely
%         id_H(id_live)=0;
        
        if sum(sum(id_H))>0
            
            H(id_H==1)=ceil(3*D(id_H==1)/4); % Fix H to 3/4 of D (true on average for living colonies)
            % Using ceil ensures there will be no 0 values for H when D=1cm, but it can produce H=D for small colonies
            
            % Finally update volume using the lookup table
            idx1=0*D;
            idx1(id_H==1)=sub2ind(size(META.V_from_H_D), D(id_H==1), H(id_H==1));
            coral(s).volume_cm3(id_H==1)=META.V_from_H_D(idx1(id_H==1));
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now, bioerosion reduces the volume of dead corals

        prop_array = prop(:,ones(size(coral(s).volume_cm3,2),1)); % makes the vector 'prop' an array with same dimensions as volume_cm3
        V_removed = ceil(coral(s).volume_cm3.*prop_array);  % estimates the volume to remove in every dead colony
        V_removed(id_live)=0;
        coral(s).volume_cm3(id_dead) = coral(s).volume_cm3(id_dead) - V_removed(id_dead); % note this may lead to negative volumes for some dead corals
Total_V_removed = Total_V_removed+sum(sum(V_removed));
        % check complete removal of dead volumes
        id0 = find(coral(s).cover_cm2 < 0 & coral(s).volume_cm3 < 1);
        D(id0) = 0 ;
        coral(s).cover_cm2(id0) = 0 ;
        coral(s).volume_cm3(id0) = 0 ; % also allows removing negatives
        
        % update the id of dead colonies
        id_live=find(coral(s).cover_cm2>0);
        id_dead=find(coral(s).cover_cm2<0);
        
        % estimate H from new volume
        H(id_dead) = ceil((8/pi)*coral(s).volume_cm3(id_dead)./((D(id_dead).^2))) ;
        
        % Get rid of dead colonies with very low height to speed-up calculations for faster runs (the threshold is setup in REEF.height_min)        
        id_low = 0*coral(s).cover_cm2 ;
        id_low(coral(s).cover_cm2<0 & H <= REEF.min_height) = 1; % dead colonies with low height have a 1
        id_low(coral(s).cover_cm2<0 & D <= REEF.min_diameter) = 1; % dead colonies with small height have a 1        
        H(id_low==1)=0 ;
        D(id_low==1)=0 ;
        coral(s).cover_cm2(id_low==1) = 0 ;
        coral(s).volume_cm3(id_low==1) = 0 ;
        
    end
      
    %%%%%%%% NOW ESTIMATE SURDACE AREA FOR ALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Surface area is calculated for dead and live colonies at the same time
    id_all = find(coral(s).cover_cm2~=0);
    idx2=0*D; 

    idx2(id_all) = sub2ind(size(META.SA_from_H_D), D(id_all), H(id_all));
    coral(s).surface_cm2(id_all)=META.SA_from_H_D(idx2(id_all));

    %%%%%%%% NOW ESTIMATE VOLUME FOR LIVE CORALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % (has been done for dead colonies already)
    idx3=0*D;
    idx3(id_live) = sub2ind(size(META.V_from_H_D), D(id_live), H(id_live));
    coral(s).volume_cm3(id_live)=META.V_from_H_D(idx3(id_live)); 

end

volume_dead_eroded_cm3 = Total_V_removed;