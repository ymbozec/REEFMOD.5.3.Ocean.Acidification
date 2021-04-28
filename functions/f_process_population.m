%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: Aug 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coral, algal, last_surface_area_grazed] = f_process_population (coral, algal, season, grazing, META, REEF, CORAL, ALGAL)

% This function processes all neighbourhood type interactions and succession for every cell
% in the grid once and once only. All the cells are processed simultaneously (vectorisation).
% Also keeps records of mortality events and calculates the fecundity total for the population
% (based on the cells after they have transformed).

% Dimensions of the grid
m = META.grid_x_count ;
n = META.grid_y_count ;

% List of grazable cells (exclude sand for cell loops)
list_cell = 1:(m*n) ;
list_cell = list_cell(REEF.grazable_cell==1) ;

% algal(1).cover_cm2(1,1) %%%%%%TEMP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extracting corals and algae records for internal use
algal_cm2 = full([algal.cover_cm2]);  % NOT a sparse matrix otherwise slow down the code
[coral_cm2, surface_cm2, volume_cm3, clade, colony_ID,species_ID] = f_struct_deploy(coral);
clear coral algal


% Create an id matrix (1 for every colony)
id1 = zeros(size(coral_cm2));
id1(coral_cm2 > 0) = 1; 

% total coral cover in every cell
total_coral_cm2 = sum(coral_cm2.*id1,2) ;

%_________________________________________________________________________________________
%
% SET UP THE SURROUNDING ENVIRONMENT AT t-1
%_________________________________________________________________________________________

% Coral and algal environments are set before processing - they are fully determined by covers at t-1
% NOTE: what about dictlob_cm2 in the 4 surrounding cells?
environ = META.environ; % changing name for tractability

% Total coral cover over 5 cells (cm2)
coral_env_cm2 = total_coral_cm2(environ(:,1)) + ...
    total_coral_cm2(environ(:,2)) + total_coral_cm2(environ(:,3)) + ...
    total_coral_cm2(environ(:,4)) + total_coral_cm2(environ(:,5)) ;

coral_env_cm2(REEF.grazable_cell==0)=0;

% Proportion of coral cover over 5 cells
% From now we work out the effect of actual surface area
substrate_env_cm2=REEF.substrate_SA_cm2(environ(:,1)) + ...
    REEF.substrate_SA_cm2(environ(:,2)) + REEF.substrate_SA_cm2(environ(:,3)) + ...
    REEF.substrate_SA_cm2(environ(:,4)) + REEF.substrate_SA_cm2(environ(:,5)) ;
coral_env_prop = coral_env_cm2./substrate_env_cm2 ;

% Neighborhood coral cover limit macroalgal growth 
coral_reduce_macrogrowth = coral_env_prop*ALGAL.coral_reduce_macrogrowth ;

% Algal covers over 5 cells (cm2)
algal_env_cm2 = algal_cm2(environ(:,1),:) + ...
    algal_cm2(environ(:,2),:) + algal_cm2(environ(:,3),:) + ...
    algal_cm2(environ(:,4),:) + algal_cm2(environ(:,5),:) ;
algal_env_cm2(REEF.grazable_cell==0)=0;

% Proportion (prop) of lob cover over 5 cells
lob_env_prop = algal_env_cm2(:,3)./substrate_env_cm2;

%_________________________________________________________________________________________
%
% PROCESS ALGAE
%_________________________________________________________________________________________

% First estimate the balance amount of grazing needed (previously ALGALREMOVAL)
% Converts proportion of total foods to actual amount removed
algal_removal = (REEF.diadema*ALGAL.diadema_props) + (grazing*ALGAL.herbivory_props);
% algal_removal contains the total amount (in %) of each alga that should be removed,
% based on grazing preferences and grazing intensity of fish/urchin (=potential effect of grazing).
actual_algal_consumpt_pct = f_algal_removal(algal_cm2, algal_removal, sum(REEF.substrate_SA_cm2)) ;

% actual_algal_consumpt_pct gives the total amount (in %) of each alga to removed based on their availability
% Note that REEF.substrate_SA_cm2 is the surface area of the substrate underneath live corals
% -> we use now the actual surface area of the reef instead of the planar area (META.total_area_cm2)

% Converts in total amount of cm2 algae that can be consumed
max_fish_consump = actual_algal_consumpt_pct*sum(REEF.substrate_SA_cm2);

% Set up a id matrix of algae to consume in every cell
eaten_algal = zeros(m*n, META.nb_algal_types) ;
% e.g., -> eaten_algal(:,3) return the cells where lobophora is grazed

%%%%%%%%%%%%%%%%%% Set up grazing for every cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = uint16(randperm(length(list_cell)));  % randomize the order of cell to visit
cumsum_algal = max_fish_consump(ones(size(list_cell(r),2),1),:) - cumsum(algal_cm2(list_cell(r),:), 1) ;
eaten_algal(list_cell(r),:) = sign(cumsum_algal); 
eaten_algal(algal_cm2==0)=0;
eaten_algal(eaten_algal<0)=0; % just remove the negatives% GRAZABLE=REEF.grazable_SA_cm2(1:3,1)
% NOTE in the last selected cell the alga is completely removed, not just the required cover.

% Record total algal area to calculate the surface area grazed in every cell
last_surface_area_grazed = sum(algal_cm2.*eaten_algal,2);

%%%%%%%%%%%%%%%%%% Remove LOBOPHORA due to grazing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amount of dict sitting above lob already
% initial_colocation_cm2 = sum(algal_cm2,2) + total_coral_cm2 - cell_areas_cm2 ;
initial_colocation_cm2 = max(sum(algal_cm2,2) + total_coral_cm2 - REEF.substrate_SA_cm2,0);
% note max(x, 0) gives a 0 when x is negative

% where lobophora is eaten so is the dictyota that is with it
algal_cm2(eaten_algal(:,3)==1,2) = max(algal_cm2(eaten_algal(:,3)==1,2) - initial_colocation_cm2(eaten_algal(:,3)==1),0) ;
% all lobophora eaten so it turns to turf
algal_cm2(eaten_algal(:,3)==1,1) = algal_cm2(eaten_algal(:,3)==1,1) + algal_cm2(eaten_algal(:,3)==1,3) ;
% no more lobophora
algal_cm2(eaten_algal(:,3)==1,3)= 0 ;


%%%%%%%%%%%%%%%%%% Now process LOBOPHORA growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lob can grow only in cells that were not grazed, for any kind of algae
id_nongrazed = zeros(m*n,1) ;
id_nongrazed(sum(eaten_algal,2)==0) = 1 ; % id of cells that were not grazed
id_nongrazed = sparse(id_nongrazed.*REEF.grazable_cell);  % exclude non-grazable cells

a = 3 ; % id of Lobophora

[algal_cm2, lob_over_dict_cm2] = f_algal_growth(algal_cm2, id_nongrazed, a, ...
   META.lobophora_dynamics, lob_env_prop, coral_reduce_macrogrowth, initial_colocation_cm2, REEF.substrate_SA_cm2, REEF.exposure);

%%%%%%%%%%%%%%%%%% Now process DICTYOTA grazing and growth%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the amount of dictyota which co-located with lobophora already
% existing_colocation_cm2 = sum(algal_cm2,2) + total_coral_cm2 - cell_areas_cm2;
existing_colocation_cm2 = max(sum(algal_cm2,2) + total_coral_cm2 - REEF.substrate_SA_cm2,0);
% note max(x, 0) gives a 0 when x is negative


% check existing_colocation_cm2 with lob_over_dict_cm2 t see how big the difference is

if (REEF.dictyota_declines_seasonally == 1 && season == 1)  %season = 1|0 -> winter|summer
   
   % if winter, all dictyota killed so it turns to turf
   algal_cm2(:,1) = algal_cm2(:,1) + algal_cm2(:,2) - existing_colocation_cm2 ;
   % no more dictyota, any that was with lobophora is just lobophora now% GRAZABLE=REEF.grazable_SA_cm2(1:3,1)

   algal_cm2(:,2) = 0 ;
   dict_over_lob_cm2 = algal_cm2(:,2) ; % null

else % otherwise all will be eaten "including dislodging any lob underneath"
   
   % all dictyota killed so it turns to turf
   algal_cm2(eaten_algal(:,2)==1,1) = algal_cm2(eaten_algal(:,2)==1,1) + algal_cm2(eaten_algal(:,2)==1,2) ;
   % also reduce any co-located lobophora
   algal_cm2(eaten_algal(:,2)==1,3) = max(algal_cm2(eaten_algal(:,2)==1,3) - existing_colocation_cm2(eaten_algal(:,2)==1),0) ;
   % no more dictyota
   algal_cm2(eaten_algal(:,2)==1,2) = 0 ;
   
   % Dictyota can grow in cells that were not grazed (whatever the algae)
	a = 2; % id of Dictyota
   dict_env_prop = 0 ; % note: no need of dict_env_prop because growth of Dictyota is driven by wave exposure (see code inside)
   
   [algal_cm2, dict_over_lob_cm2] = f_algal_growth(algal_cm2, id_nongrazed, a, ...
       META.dictyota_dynamics, dict_env_prop, coral_reduce_macrogrowth, existing_colocation_cm2, REEF.substrate_SA_cm2, REEF.exposure);

end

% Final check to feed correct amount of dictlob/lobdict to reconciliation - indispensable
% dictlob_cm2 = sum(algal_cm2,2) + total_coral_cm2 - cell_areas_cm2;
% existing_colocation_cm2=sparse(existing_colocation_cm2)
% dict_over_lob_cm2
% dictlob_cm2 = dict_over_lob_cm2;
dictlob_cm2 = max(sum(algal_cm2,2) + total_coral_cm2 - REEF.substrate_SA_cm2,0);


% remove used objects 
clear r eaten_algal coral_env_prop max_fish_consump coral_reduce_macrogrowth  id_nongrazed algal_env_cm2 coral_env_cm2
%_________________________________________________________________________________________
%
% PROCESS CORALS
%_________________________________________________________________________________________

%%%%%%% 1) PARROTFISH PREDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_pred = id1; % id1 is id matrix of living colonies (see above)
% assigns a random proba of mortality to every colony
rand_mort = rand(size(id_pred));
% switch to 0 colonies escaping predation due to large size
id_pred(coral_cm2 > CORAL.threshold_predation_size) = 0 ;
% % switch to 0 colonies escaping predation by chance
id_pred(rand_mort > CORAL.parrotfish_predation) = 0 ;

% update the grid (id_pred == 1 identifies the colonies eaten)
algal_cm2(:,1) = algal_cm2(:,1) + sum(coral_cm2.*id_pred,2) ; % update turf with coral loss
coral_cm2(id_pred==1) = 0 ; % Eaten colonies are not kept - flat substratum (complete removal of small colonies)

%%%%%%% 2) MACROALGAL OVERGROWTH OF CORALS %%%%%%%%%%%%%%%%%%
% Update the proportions of macroalgae over the 5 cells after grazing and algal growth
% local macroalgal abundance for competitive interactions

% New algal covers over 5 cells (cm2)
algal_env_cm2 = algal_cm2(environ(:,1),:) + ...
    algal_cm2(environ(:,2),:) + algal_cm2(environ(:,3),:) + ...
    algal_cm2(environ(:,4),:) + algal_cm2(environ(:,5),:) ;

algal_env_cm2(REEF.grazable_cell==0)=0; %just to be sure there is no algae on sand

% Calculate proportion of macroalgal cover over 5 cells
lob_env_prop = algal_env_cm2(:,3)./substrate_env_cm2; % Proportion (prop) of lob cover over 5 cells
dict_env_prop = algal_env_cm2(:,2)./substrate_env_cm2; % Proportion (prop) of dict cover over 5 cells
macroalgal_env_prop = lob_env_prop + dict_env_prop ; % Proportion (prop) of all macroalgae

% Lobophora overgrowth
[coral_cm2, algal_cm2] = f_macroalgae_overgrowth_corals(coral_cm2, algal_cm2, species_ID, 3, lob_env_prop, CORAL.lobophora_reduce_rate, META.nb_coral_types); % /7 is unclear (mult in the reference code)

% Dictyota overgrowth
[coral_cm2, algal_cm2] = f_macroalgae_overgrowth_corals(coral_cm2, algal_cm2, species_ID, 2, dict_env_prop, CORAL.dictyota_reduce_rate, META.nb_coral_types);

% Note that overgrowth leads to shrinking so that if a colony dies it becomes 0, not negative

%%%%%%% 3) CORAL GROWTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE that the reference code compares potential coral growth with the total cell area,
% while available space should account for macroalgal cover
% Here available space is the free space for growth (should it be turf?)
total_coral_cm2 = sum(coral_cm2.*id1,2) ; % update total coral cover in every cell
avail_cell_areas_cm2 = REEF.substrate_SA_cm2 - sum(algal_cm2(:,2:4),2) + dictlob_cm2 - total_coral_cm2;
avail_cell_areas_cm2(REEF.grazable_cell==0)=0; % cannot grow on sand

[coral_cm2] = f_coral_growth(coral_cm2, species_ID, clade, macroalgal_env_prop, total_coral_cm2, avail_cell_areas_cm2,REEF.substrate_SA_cm2, META, CORAL, ALGAL) ;
% NOTE: algal turf not updated here but at the end (f_adjust_total_cover)


%%%%%%% 4) NATURAL MORTALITY   %%%%%%%%%%%%%%%%%%%%%%%%%
% Includes now (27/08/13) partial and whole colony mortalities 
[coral_cm2, algal_cm2]=f_natural_mortality(coral_cm2, algal_cm2, CORAL, species_ID, META);

%_________________________________________________________________________________________
%
%       FINALLY ADJUST THE TOTALS FOR CORAL GROWTH AND MORTALITY
%_________________________________________________________________________________________

% %%%%%%% Calculate overshoot for adjustment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1(coral_cm2<=0)=0;  % update the identity matrix of living colonies

overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - dictlob_cm2 - REEF.substrate_SA_cm2 ;
overshoot(REEF.grazable_cell==0)= 0;

% IF OVERSHOOT IS NEGATIVE -> fill the empty space with turf
% IF OVERSHOOT IS POSITIVE -> reduce turf covers
algal_cm2(:,1) = algal_cm2(:,1) - overshoot ;

% This may have generated negative turf size if the sum of all other covers
% is still too large for the cell. Then force turf to be null...
algal_cm2(algal_cm2(:,1)<0,1) = 0 ;
% ... and calculate the new overshoot
overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - dictlob_cm2 - REEF.substrate_SA_cm2 ;
overshoot(REEF.grazable_cell==0)= 0;

if sum(overshoot)>0
    
    % FOR CELLS WITH POSITIVE OVERSHOOT -> runs through hierarchical sequence
    % adjusting macroalgae for cell size
    list_cell2 = list_cell(overshoot(REEF.grazable_cell==1) > 0) ;
    
    for cell_id=list_cell2
        [algal_cm2(cell_id,:), dictlob_cm2(cell_id,1)] = ...
            f_adjust_total_cover (algal_cm2(cell_id,:), dictlob_cm2(cell_id,1), overshoot(cell_id,1)) ;
    end
    
end


% Then store the new algal cover
for a=1:META.nb_algal_types
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end


%%%%%%% Record the new coral cover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, META.nb_coral_types, META.doing_3D);
