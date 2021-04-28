%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: Aug 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [coral, algal] = f_hurricane_effect(coral, algal, hurricaneCat, META, CORAL)

if hurricaneCat > 0
    
    % IF YOU NEED TO ELEVATE THE PROBABILITY DISTRIBUTION, JUST ADD AND EXTRA 0.2 OR SOMETHING
    % SO PROB = ANSWER+0.2
    % THIS ONE IS MODIFIED BY PETE IN APRIL 07 TO REDUCE IMPACT. LATER MODIFIED BY HELEN (HJE)
    % TO SIMULATE MORTALITY BEING DEPENDENT UPON THE STRENGTH (CATEGORY) OF HURRICANE
    % Re-written by YMB: 15/03/2012
    
    %%%%%%%%%%%%%%%%%%%%%%%%   HJE MODIFICATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%
    % Original hurricane probability equation (function of colony size):
    %       hurricane_prob=(-3e-007*pop(x,y,c)^2)+(0.0007*pop(x,y,c))+0.0551;
    %
    % Rearranging into the more familiar format:
    %       hurricane_prob = a(pop(x,y,c)-h)^2 + k;
    % where the x and y coordinates of the vertex are the values of h and k
    % respectively.
    %
    % Multiplying this out we get ax^2-2ahx+ah^2+k = -3e-007*x^2+0.0007*x+0.0551
    % Therefore
    %           a = -3e-007
    %           h = 0.0007/(-2*a) = 1166.7
    %           k = 0.0551-ah^2 = 0.4138
    
    % To modify mortality based on hurricane category we modify a and k
    % First make the function symmetrical (so that we can lower the curve
    % without lowering the end points):
    h = 1250; %Value was previously ~1167
    a = -3e-007;
    k_cat5 = 0.0551-a*((0.0007/(-2*a))^2); %(Uses old value for h to achieve same maximum as before)
    % The above values are appropriate for category 5 hurricanes. For the other
    % categories we reduce a and k based on the relative predicted impacts
    % (using Madin, Zimmerman etc).
    relHurrMort = [0.0461 0.1183 0.2504 0.5677 1]; %e.g. a cat 4 has 57% of the impact as a cat 5
    a = a*relHurrMort(hurricaneCat);
    k = k_cat5*relHurrMort(hurricaneCat);
    
    %%%%%%%%%%%%%%%%%%%%%%  END OF HJE MODIFICATIONS   %%%%%%%%%%%%%%%%%%%%%
    
    % TODO: this should be an output of the function
    % hu_info = zeros(META.nb_coral_types,1); % records the number of adults being pushed into adol size class
    
    %_______________________________________________________________________
    %
    % Hurricane effects on macroalgae
    %_______________________________________________________________________
    
    algal_cm2 = [algal.cover_cm2] ;
    [m,n] = size(algal_cm2) ;
    
    id1 = zeros(m,n-1) ; % exclude turf
    rand_effect = rand(m,n-1) ;
    id1(rand_effect < META.hurricane_effect_on_macroalgae) = 1 ;
    macroalgal_loss_cm2 = id1 .* algal_cm2(:,2:n) ;
    
    % remove fixed proportion of macroalgal and set to turf
    algal_cm2(:,1) = algal_cm2(:,1) + sum(macroalgal_loss_cm2,2) ;
    % update macroalgal cover for losses
    algal_cm2(:,2:n) = algal_cm2(:,2:n) - macroalgal_loss_cm2 ;
    
    %_______________________________________________________________________
    %
    % Hurricane effects on corals
    %_______________________________________________________________________
    
    [coral_cm2, surface_cm2, volume_cm3, clade,colony_ID, species_ID] = f_struct_deploy (coral);
    
    %%%% This is new stuff (March 2013) for implementing species-specific hurricane mortalities
    [i,j]=size(coral_cm2) ;
    sensitivity_hurricane = zeros(i,j);
    extent_partmort = sensitivity_hurricane ;
    
    col_start = 1;
    col_stop = 0;
    
    for s = 1:META.nb_coral_types
        
        col_stop = col_stop + species_ID(s) ;
        sensitivity_hurricane(:,col_start:col_stop)= CORAL.sensitivity_hurricane(s);
        extent_partmort(:,col_start:col_stop) = CORAL.extent_hurricane(s);
        col_start = col_start + species_ID(s) ;
        
    end
    
    % Initialisation of colony localization on the grid
    id0 = sparse(zeros(size(coral_cm2)));
    
    % WHOLE-COLONY MORTALITY OF LIVE CORALS (=REMOVAL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    id01 = id0 ;
    id02 = id0 ;
    
    rand_prop = sprand(coral_cm2) ; % generate a random probability for every colony
    
    % Whole mortality of recruits: 80% due to sand scouring
    id01(coral_cm2 > 0 & coral_cm2 < META.adol_size) = 1 ; % flag the recruits ("<=" in the reference code)
    id01(rand_prop >= META.hurricane_effect_on_recruits) = 0 ; % exclude the survivors
    
    % Whole mortality of adults
    id02(coral_cm2 >= META.adol_size) = 1 ; % flag the others
    hurricane_prob = (a*((coral_cm2.*id02) - h).^2) + k ;  %% MODIFIED - BASED ON STRENGTH OF HURRICANE
    hurricane_prob = sensitivity_hurricane.*hurricane_prob; % species-specific rates FOR COZUMEL
    id02(rand_prop >= hurricane_prob) = 0 ; % exclude the survivors
    
    coral_loss = coral_cm2.*(id01 + id02);
    algal_cm2(:,1) = algal_cm2(:,1) + sum(coral_loss,2); % turn into turf the lost colonies
    coral_cm2(coral_loss>0) = 0 ; % according to Pete, colonies are removed so do not turn into dead
    
    % WHOLE-COLONY MORTALITY OF DEAD CORALS (=REMOVAL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    id04=id0;
    id04(coral_cm2 < 0) = 1 ; % flag the deads
    hurricane_prob_deads = (a*(abs(coral_cm2.*id04) - h).^2) + k  ;
    hurricane_prob_deads = META.hurri_whole_mort_dead* sensitivity_hurricane.*hurricane_prob_deads ;
    id04(rand_prop >= hurricane_prob_deads) = 0 ; % exclude the deads not affected by chance
    algal_cm2(:,1) = algal_cm2(:,1) - sum(coral_cm2.*id04,2); % minus because coral_cm2 is negative for deads
    coral_cm2(id04==1)=0; % unlikely deads are removed from the grid
    
    %%% NOTE: surface_cm2 and volume_cm3 are updated at the end (in f_struct_rebuild)
    
    % PARTIAL MORTALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    id03 = id0;
    id03(coral_cm2 >= META.adult_size) = 1 ; % flag the adults (">" in the reference code)
    
    hu_partmort_prob = 0.30*relHurrMort(hurricaneCat) + ...
        0.2*relHurrMort(hurricaneCat).*rand_prop ; %% MODIFIED - BASED ON STRENGTH OF HURRICANE
    
    hu_partmort_prob = extent_partmort.*hu_partmort_prob ; % Cozumel: higher susceptibility for agaricia and porites
    hu_partmort_prob(hu_partmort_prob > 1) = 1 ;
    hu_mort = round(coral_cm2 .* id03 .* hu_partmort_prob);
    
    algal_cm2(:,1) = algal_cm2(:,1) + sum(hu_mort,2) ; % dead tissue becomes turf
    coral_cm2 = coral_cm2 - hu_mort ; % update colonies with losses
    
    % 	hu_info(s) = hu_info(s) + nnz(hu_mort) ;
    % NOTE: the hu_info now record the adults pushed back to adolescent size for each coral type as in the C++ version
    
    
    %%%%%%%% Before leaving, store the new covers into 'coral' and 'algal'%%%%%%%%%%%%%%%%%%%%
    
    [coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, META.nb_coral_types, META.doing_3D);
    
    for a=1:size(algal_cm2, 2)
        algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
    end
    
end
