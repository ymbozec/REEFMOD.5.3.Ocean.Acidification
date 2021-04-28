function [count,class] = f_count_sizefreq(cover_cm2, META)
% CALCULATES THE NUMBER OF CORALS OF A GIVEN SPECIES IN A SIZE CLASS

colony_sizes = cover_cm2(:) ;

juveniles = colony_sizes(colony_sizes < META.adol_size & colony_sizes > 1) ;
adolescents = colony_sizes(colony_sizes < META.adult_size & colony_sizes >= META.adol_size) ;
adults = colony_sizes(colony_sizes >= META.adult_size) ;

[count.juv, class.juv] = hist(juveniles,0:META.size_bins(1):META.adol_size);
[count.adol, class.adol] = hist(adolescents,0:META.size_bins(2):META.adult_size);
[count.adult, class.adult] = hist(adults,0:META.size_bins(3):floor(pi*(META.cell_x_size/2)^2));
