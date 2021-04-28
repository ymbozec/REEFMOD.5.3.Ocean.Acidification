%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Oct 2012.
% Last modified: Aug 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   [coral, algal, REEF]= f_estimate_3D_reef(coral, algal, cell_area_cm2, REEF, META)

[coral_cm2, surface_cm2, volume_cm3, clade, colony_ID, species_ID] = f_struct_deploy(coral);
algal_cm2 = [algal.cover_cm2] ;

% test_a = sign(algal_cm2);
% it=find(test_a<0);
% if it~=0
% algal_cm2(it, :)
% disp('before ajusting')
% stop
% end
algal_cm2(algal_cm2<0)=0 ; % sometimes turf takes a -1 value for some reason, so we arbitrary fix the pb

initial_substrate_SA_cm2 = REEF.substrate_SA_cm2;

coral_PA_cm2 = coral_cm2;  % individual planar areas in every cell
coral_SA_cm2 = surface_cm2;  % individual surface areas in every cell

id_live = zeros(size(coral_PA_cm2));
id_dead = id_live ;

id_live(coral_PA_cm2>0) = 1 ;
id_dead(coral_PA_cm2<0) = 1 ;

live_PA_cm2 = sum(coral_PA_cm2.*id_live,2);  % total planar area of living colonies in every cell
live_SA_cm2 = sum(coral_SA_cm2.*id_live,2);  % total surface area of living colonies in every cell

dead_PA_cm2 = -sum(coral_PA_cm2.*id_dead,2); % total planar area of dead colonies in every cell
dead_SA_cm2 = sum(coral_SA_cm2.*id_dead,2); % total surface area of dead colonies in every cell
        
% Estimate the new area of every cell - basically this is the actual
% surface area of the substrate underneath live corals - this is for space occupation
REEF.substrate_SA_cm2 = dead_SA_cm2 + cell_area_cm2 - dead_PA_cm2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then adjust the benthic cover to the new substrate area %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_algal_types = size(algal_cm2,2) ;
diff_substrate_SA_cm2 = REEF.substrate_SA_cm2 - initial_substrate_SA_cm2;
% This is positive in cells where substratum has increased -> must be colonized by turf
% This is negative in cells where substratum has shrunk -> benthic cover must be reduced
id_diff = sign(diff_substrate_SA_cm2); % looking at increase (+1) or decrease (-1) in surface

%%%%% 2) Reduction of substrate lead to a reduction of algae
% Note after this, there won't be any overtopping macroalgae in the selected cells
% check = sum(algal_cm2,2) + live_PA_cm2 - REEF.substrate_SA_cm2;
% check = - check .*id_diff ; % ensure we select only the cells that have shrink
% id_check = sign(check); % must have only positive 1 = id of cells for which benthic 
% first_check=sum(id_check==-1)
% while sum(id_check)>0
id_diff(id_diff>0)=0; % we process only cells that have shrunk

for a = 1:nb_algal_types
    
    new_algal_cm2 = 0*id_diff ;
    % This is how much algal type 'a' will be left after reduction
    new_algal_cm2(id_diff==-1,1) = algal_cm2(id_diff==-1,a) + diff_substrate_SA_cm2(id_diff==-1,1);
    % Can't remove more than the current alga
    new_algal_cm2(new_algal_cm2<0,1)=0;
    % update the differential to overshoot
    diff_substrate_SA_cm2(id_diff==-1,1)= diff_substrate_SA_cm2(id_diff==-1,1)+ ...
        algal_cm2(id_diff==-1,a) - new_algal_cm2(id_diff==-1,1);
    
    % update the algal type
    algal_cm2(id_diff==-1,a) = new_algal_cm2(id_diff==-1,1);
    % update the check
    id_diff(diff_substrate_SA_cm2==0)=0;
end


if sum(id_diff)<0

    % if benthic cover is still too high, then reduce corals (minor shrinkage)
    id_reduce = id_diff(:,ones(size(coral_cm2,2),1)).*id_live; % selects the corals that must be reduced
    diff_matrix = diff_substrate_SA_cm2(:,ones(size(coral_cm2,2),1));
    
    % reduce each coral type relative to their cover in the cell
    prop_coral = sparse(zeros(size(coral_cm2)));
    new_coral_cm2 = prop_coral ;

    total_living_coral_cm2 = live_PA_cm2(:,ones(size(coral_cm2,2),1));
    prop_coral(id_reduce==-1) = coral_cm2(id_reduce==-1)./total_living_coral_cm2(id_reduce==-1);

    new_coral_cm2(id_reduce==-1) = coral_cm2(id_reduce==-1) + ...
        floor(prop_coral(id_reduce==-1).*diff_matrix(id_reduce==-1));
    % NOTE: using floor ensures that we won't remove less than what is
    % necessary (floor of a negative). A little bit more of coral than required might be
    % removed but this concerns few cm2 for each
    
    new_coral_cm2(id_reduce==-1)=max(new_coral_cm2(id_reduce==-1),0);
    % Using max(,0) ensures no colony will have negative cover and
    % id_reduce only identfies living corals
    
    % update corals
    coral_cm2(id_reduce==-1) = new_coral_cm2(id_reduce==-1);
   
    idx1 = find(id_reduce==-1 & coral_cm2==0); % look for colonies that disappeared
    surface_cm2(idx1)=0;
    volume_cm3(idx1)=0;
    id_reduce(idx1)=0;
    
    % update the total planar area of living corals per cell
    live_PA_cm2 = sum(coral_cm2.*id_live,2);
    
    % Update H -> this is an approximation: H = 3*D/4 (true on on average)
    % Note this is temporary anyway since H will be re-estimated with
    % f_estimate_3D_colony properly (with the species-specific relationships)
    % Note this only affects the live corals ; we do that pragmatic way because
    % estimating H from D with lookup tables would have to be done species per species
    D = zeros(size(new_coral_cm2)) ;
    H = D;
    D(id_reduce==-1) = ceil(2*sqrt(new_coral_cm2(id_reduce==-1)/pi));
    H(id_reduce==-1) = ceil(3*D(id_reduce==-1)/4);

    % Update surface area of living corals
    idx2=0*new_coral_cm2;
    idx2(id_reduce==-1) = sub2ind(size(META.SA_from_H_D), D(id_reduce==-1), H(id_reduce==-1));
    surface_cm2(id_reduce==-1) = META.SA_from_H_D(idx2(id_reduce==-1));
    
    % update the total surface area of living corals per cell
    live_SA_cm2 = sum(surface_cm2.*id_live,2);
    
    % Update volume of living corals
    idx3=0*new_coral_cm2;
    idx3(id_reduce==-1) = sub2ind(size(META.V_from_H_D), D(id_reduce==-1), H(id_reduce==-1));
    volume_cm3(id_reduce==-1) = META.V_from_H_D(idx3(id_reduce==-1));

end  

% Estimate the total SA (including the SA of live corals) - this is for rugosity calculation
REEF.floor_SA_cm2 = REEF.substrate_SA_cm2 - live_PA_cm2  + live_SA_cm2 ;

%%%%% 2) Addition of grazable surface is immediately covered by turf
algal_cm2(:,1) = REEF.substrate_SA_cm2 - (live_PA_cm2 + sum(algal_cm2(:,2:nb_algal_types),2));
algal_cm2(REEF.grazable_cell==0)=0;
algal_cm2(algal_cm2<0)=0; % because produces sometimes negative cm2 for turf
%Note processing this here allow to overshoot the shrunk cells for which too much coral has been removed

%%%%% Then rebuild the coral/ algal matrices
[coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, META.nb_coral_types, META.doing_3D);

for a=1:nb_algal_types
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end

% test_a = sign(algal_cm2);
% it=find(test_a<0);
% if it~=0
% algal_cm2(it, :)
% disp('after ajusting')
% stop
% end