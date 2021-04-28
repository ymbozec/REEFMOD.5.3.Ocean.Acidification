%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2015.
% Last modified: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now multiple recruits per cell

function [coral, total_settled, ID_colony_tracking, algal] = f_apply_recruitment_new(coral, algal, META, grazable_cell, floor_SA_cm2, max_density_settlers, clade_prop,ID_colony_tracking)

nb_coral_types = size(coral,2) ;
total_settled = zeros(1,nb_coral_types) ; 
random_ct = randperm(nb_coral_types) ; % list coral types in random order

turf_cm2 = algal(1).cover_cm2 ;


for s = random_ct
  
    [l,c] = size(coral(s).cover_cm2);
    id1 = zeros(l,c);
    id1(coral(s).cover_cm2>0)=1;
    
    colony_count = sum(id1,2) ; % number of colonies of a given species for each cell

    lambda = full(max_density_settlers(s) * turf_cm2./ floor_SA_cm2); % mean parameter of the poisson distribution of number of recruits
    % Here, lambda is adjusted to available space (prop of turf) and larval input  
    new_settler_count = poissrnd(lambda, l, 1).*grazable_cell ;  %cannot settle on sand
    new_settler_count(colony_count + new_settler_count > META.max_colonies) = META.max_colonies - colony_count(colony_count + new_settler_count > META.max_colonies);
    MAX = max(new_settler_count); % actual max nb of candidates for settlement in a cell
    
    if MAX==0
        continue % No recruitment for species s, so jump to the next
    end
    
    add_settler = zeros(l,MAX) ;
    
    for j = 1:MAX     
        add_settler(:,j)=1; % add a settler in all cell
        add_settler(sum(add_settler,2)>new_settler_count,j)=0; % delete settlers in excess
    end

    coral(s).cover_cm2(:,c+(1:MAX))= add_settler ;

    total_settled(1,s) = sum(new_settler_count); % record the total number of settlers for this species
    
    settler_ID = zeros(size(add_settler));
    settler_ID(add_settler==1)=[1:1:total_settled(1,s)]+ID_colony_tracking(s);
    coral(s).colony_ID(:,c+(1:MAX)) = settler_ID;% assigns ID to the new colonies, starting from ID max   
    ID_colony_tracking(1,s) = ID_colony_tracking(1,s)+total_settled(1,s);% Record the new ID max
    
    rand_clade = rand(size(add_settler)).*add_settler ; % by default, all recruits have the sensitive clade (clade = 1)
    add_settler(rand_clade > clade_prop) = 2 ; % recruits have the same proportion of clade C over D than at initial step
    coral(s).clade(:,c+(1:MAX)) = add_settler ;
    
    turf_cm2 = turf_cm2 - new_settler_count; % Update turf (only one recruit per species per cell)
  
    if META.doing_3D == 1
        
        coral(s).surface_cm2(:,c+(1:MAX))= add_settler;
        coral(s).volume_cm3(:,c+(1:MAX))= add_settler;
        
    end    
end

algal(1).cover_cm2 = turf_cm2;
