%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implements the growth of coral colonies according to macroalgae and coral-coral competition

function [coral_cm2] = f_coral_growth(coral_cm2, species_ID, clade, macroalgal_env_prop, total_coral_cm2,avail_cell_areas_cm2, cell_area_cm2, META, CORAL, ALGAL)

id1=spones(coral_cm2);
id1(coral_cm2<0)=0 ; % exclude dead colonies (negatives)

[i,j]=size(coral_cm2) ; % size indicators of the coral matrix

% turn into matrices:
macroalgal_env_prop = id1.*macroalgal_env_prop(:,ones(1,j)); % macroalgal environment of every colony (same env if same cell)
% total_coral_cm2 = id1.*total_coral_cm2(:,ones(1,j)); % total coral cover in every cell for coral-coral interactions

if META.doing_coral_competition == 1
    load('LM_simulated_competition.mat')
    total_coral_pct = full(100 * total_coral_cm2./cell_area_cm2);
    Coral_contact = random(LinearModelCompetition,total_coral_pct);
    Coral_contact(Coral_contact<0)=0;
    Coral_contact(Coral_contact>100)=100;
    Coral_contact = uint8(Coral_contact+1); % needs to add 1 to get the right index
    % (0 contact = first value of growth reduction)
    
else
    
    Coral_contact = uint8(1); % if no contact then picks first value of growth reduction (ie, for 0% contact)
    
end


%%

algal_contact = id1;
algal_contact(macroalgal_env_prop >= ALGAL.critical_algal_contact)=2;
algal_contact(macroalgal_env_prop >= ALGAL.vcritical_algal_contact)=3;

deltagrowth = id1 ;% preallocating space for the proportional reduction in coral growth due to macroalgae
deltagrowth(coral_cm2 > 0 & coral_cm2 <= META.adol_size & algal_contact==2) = ALGAL.macroalgal_coral_recruit_growth_rate ;
deltagrowth(coral_cm2 > 0 & coral_cm2 <= META.adol_size & algal_contact==3) = 0 ;
deltagrowth(coral_cm2 > META.adol_size & algal_contact==2) = ALGAL.macroalgal_coral_growth_rate ;

growth_rates = zeros(i,j);
max_size = zeros(i,j);
over_cell = zeros(i,1);

col_start = 1;
col_stop = 0;


for s = 1:META.nb_coral_types
    
    col_stop = col_stop + species_ID(s) ;
    growth_rates(:,col_start:col_stop)= CORAL.growth_rate(s);
    max_size(:,col_start:col_stop)= CORAL.max_size(s);

    % further adjust for coral-coral competition
    if META.doing_coral_competition==1
        
        % First determine the relative proportion of coral conspecifics
        conspec_cm2 = coral_cm2(:,col_start:col_stop).*id1(:,col_start:col_stop);
        prop_conspec = sum(conspec_cm2,2)./total_coral_cm2 ;
        growth_reduction = prop_conspec.*CORAL.growth_reduction_cspec(Coral_contact,s) + ...
            (1 - prop_conspec).*CORAL.growth_reduction_hspec(Coral_contact,s);
        growth_reduction(isnan(growth_reduction)==1)=0;
        growth_reduction = growth_reduction(:,ones(1,col_stop-col_start+1));  % duplicate columns
        growth_rates(:,col_start:col_stop) = growth_reduction.*growth_rates(:,col_start:col_stop);
    end

    col_start = col_start + species_ID(s) ;
end
  
% adjust growth rate following coral clade
id_clade2 = (clade.*id1) - id1 ; % eliminates id of clade 1 ("1")
id_clade1 = id1 - id_clade2 ; % eliminates id of clade 2 ("2")
growth_rates = (growth_rates .* id_clade1) + (growth_rates .* id_clade2)*CORAL.clade_reduced_growth ;
   
% increases overall size by one growth increment -> this gives a first potential growth for a hypothetical free space
newsize = pi*(growth_rates + sqrt(coral_cm2 .* id1/pi)).^2 ; % note final sizes will be rounded at the end

% just check if new sizes do not exceed maximum size of each species (cannot grow bigger than the max size)
check_size = max_size - newsize ;
newsize(check_size<0) = max_size(check_size<0) ; % if too big just assign the maximum size

% Now, reduce the potential growth due to algae. Note deltagrowth takes values 1 (growth is not reduced), 0.1 or 0.3
growth_cm2 = (newsize - (coral_cm2 .* id1)).*deltagrowth ;

% identify cells where corals will overload the cell area
total_growth = sum(growth_cm2,2); % sum of all potential growth per cell
over_cell(avail_cell_areas_cm2 - total_growth < 0) = 1 ;

id_over = over_cell(:,ones(1,j)) .* id1 ; % gives 1 for cells where total_growth overtake available space
id_ok = id1 - id_over ; % identify with 1 the cells where total_growth can take the available space

% Corals fill the available space proportionately to their growth. If no available space, corals don't grow 
total_growth = total_growth(:,ones(1,j)).*id_over; % tranform into a matrix with same dimensions as coral colonies
avail_cell_areas_cm2 = avail_cell_areas_cm2(:,ones(1,j)).*id_over ; %turn it into a matrix form
coral_cm2(id_over==1) = coral_cm2(id_over==1) + avail_cell_areas_cm2(id_over==1).*growth_cm2(id_over==1)./total_growth(id_over==1) ;

% for the OK cells, all corals just grow freely
coral_cm2 = coral_cm2 + growth_cm2.*id_ok ;

% Finally round all the sizes after growth
% Flooring let a couple of cm2 free of corals so future recruits may settle again but won't ever grow
coral_cm2 = floor(coral_cm2);
