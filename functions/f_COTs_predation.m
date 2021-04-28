%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2015.
% Last modified: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coral,algal] = f_COTs_predation(coral, algal, COTs_density, COTs_feeding_rate, COTs_feeding_prefs)

% eaten_coral_cm2 is the total amount of cm2 coral that can be consumed
% over the grid for each coral species and during 6 months

% First need to estimate relative prop of each coral species
coral_relative_prop = zeros(1, size(coral,2)) ;
% total coral cover on the reef
coral_cm2 = [coral.cover_cm2] ;
total_coral_cm2 = sum(sum(coral_cm2(coral_cm2>0))) ;

for s = 1:size(coral,2)
    coral_relative_prop(1,s) = sum(sum(coral(s).cover_cm2(coral(s).cover_cm2>0)))/total_coral_cm2;
end

D = COTs_feeding_prefs.*coral_relative_prop' ;

% Determine total amount (cm2) of each coral species eaten by COTs
eaten_coral_cm2 = round( COTs_density*COTs_feeding_rate*D/sum(D) ) ;
% total_eaten = sum(eaten_coral_cm2)
% actual_consumption = zeros(size(eaten_coral_cm2));

for s = 1:size(coral,2)
    
    if eaten_coral_cm2(s)==0
        continue
    else
        
        cover_cm2 = coral(s).cover_cm2 ; % temporary storage to speed up code
        I = find(cover_cm2>0) ; % Spot the living colonies
        which_eaten = zeros(size(I)) ;
        r = uint32(randperm(length(I))) ;
        
        % Calculate the differential from total consumption and the cumulative sum of all colony areas (randomly picked-up)
        cumsum_coral = eaten_coral_cm2(s)*ones(length(I),1) - cumsum(cover_cm2(I(r)), 1) ;
        % Positive values in cumsum_coral indicate coral colonies that will
        % be removed to reach the total area of coral consumed
        
        last_eaten = find(cumsum_coral<0,1,'first'); % the last eaten colony is the first to produce a negative differential
        % (the size of this colony overtakes the required total of areas consummed by COTs).
        % Here we assume this colony is entirely eaten, which makes the total effectively consumed a little bit higher than
        % expected (conservative assumption). As a result, there is no partial mortality due to COTs (they
        % finish the job before moving to the next colony (or reef)
        
        cumsum_coral(last_eaten)= 1;
        which_eaten(I(r),1) = sign(cumsum_coral);
        which_eaten(which_eaten<0,1)=0; % excludes negatives (coral cm2 in excess from consumption)

        temp = zeros(size(cover_cm2));
        temp(which_eaten==1) = cover_cm2(which_eaten==1) ;
        algal(1).cover_cm2 = algal(1).cover_cm2 + sum(temp,2) ;
         
        coral(s).cover_cm2(which_eaten==1) = - cover_cm2(which_eaten==1) ;
%         actual_consumption(s)= sum(sum(temp));
    end

end

% total_consummed = sum(actual_consumption)