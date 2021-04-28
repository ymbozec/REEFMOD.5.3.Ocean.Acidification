%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: Jan 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only 1 recruit per cell per species

function [coral, algal] = f_apply_recruitment(coral, algal, max_colonies, grazable_cell, floor_SA_cm2, prob_recruit, clade_prop,space_limited_settlement)

% This function applies recruitment to turf patches based on the larval
% Can have up to 2 brooders and 1 of others have a go. But should
% really randomise the order in which corals are processed to give equal
% chance. So each get one go.

[i,j] = size(grazable_cell) ;

nb_coral_types = size(coral,2) ;
random_ct = randperm(nb_coral_types) ; % list coral types in random order

turf_cm2 = algal(1).cover_cm2 ;

for s = random_ct
    
    if space_limited_settlement==1
        
        P_settle = turf_cm2./ floor_SA_cm2 ;
        % the probability that a competent larva settles in a cell is the proportion of available space
        % Must be recalculated each time a species has been processed (updated turf also)
    else
        P_settle = ones(size(floor_SA_cm2)) ;
    end
    
    coral_cm2 = coral(s).cover_cm2 ;
    surface_cm2 = coral(s).surface_cm2 ;
    volume_cm3 = coral(s).volume_cm3 ;
    clade = coral(s).clade ;

    
    [l,c] = size(coral_cm2);
    coral_cm2(:,c+1) = zeros(l,1); % adds a new column for potential recruits

    rand_recruit = rand(i,1).*grazable_cell ; % assigns a random number between 0-1 to every cells
    proba = prob_recruit(s)*P_settle.*grazable_cell ; % note that in case of full recruitment proba=1 for every cell

    % Recruits settle in the new colum
    test_col1 = sign(coral_cm2);
    test_col1(test_col1==-1)=0;
    test_col2 = sum(test_col1,2);
%     test_neg = find(test_col2> max_colonies); % cannot have more colonies than the max per species
test_neg = find(test_col2> max_colonies); % cannot have more colonies than the max per species
    
coral_cm2(proba>rand_recruit,c+1) = 1;  % gives a 1 cm2 cover to a recruit
    coral_cm2(test_neg,c+1) = 0;
new_recruits = coral_cm2(:,c+1);

    surface_cm2(:,c+1) = new_recruits; % gives a surface of 1 cm2 to a recruit
    volume_cm3(:,c+1) = new_recruits; % gives a volume of 1 cm3 to a recruit
    
    rand_clade = sprand(new_recruits) ; % by default, all recruits have the sensitive clade (clade = 1)
    new_recruits(rand_clade > clade_prop) = 2 ; % recruits have the same proportion of clade C over D than at initial step
    clade(:,c+1) = new_recruits ;

    
    % Update turf (only one recruit per species per cell)
    turf_cm2 = turf_cm2 - coral_cm2(:,c+1);

    % Store the new (optimized) coral cover
    coral(s).cover_cm2 = sparse(coral_cm2) ;
    coral(s).surface_cm2 = sparse(surface_cm2) ;
    coral(s).volume_cm3 = sparse(volume_cm3) ;
    coral(s).clade = sparse(clade) ;
    
end

algal(1).cover_cm2 = turf_cm2;
