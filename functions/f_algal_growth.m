%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [algal_cm2, lobdict_cm2] = f_algal_growth(algal_cm2, id_nongrazed, a, ...
    algal_dynamics, algal_env_prop, coral_reduce_macrogrowth, colocation_cm2, cell_area_cm2, exposure)

% Determine the potential growth for macroalgal of type a in the nongrazed cells

% First thing to do is to determine the type of macroalgal
% This implements only lobophora (a=3) and Dictyota (a=2) growths)
% get the current macroalgal cover as an integer from 0 to 100;

algal_pct = uint8(full(100 * algal_cm2(:,a).*id_nongrazed./ cell_area_cm2));

% Need to unsparse matrix to account fo 0 pct?

algal_pct(algal_pct>100)=uint8(100);

lobdict_cm2 = 0*id_nongrazed;
% this apparently 'adjusts using method in productivity paper'
% this is what the macroalgal cover wants to become
switch a
    
    case 3 % the algae to be processed is Lobophora
        % algal_dynamics here corresponds to 'lob_dynamics' (101x3) uploaded
        % from the data file 'dict_and_lob_for_model.mat', where growth is
        % fitted to the 4-cell neighbourhood coverage of Lobophora
        % based on de Ruyter van Steveninck and Breeman (1987)
        
        % col 1 is ignored (id of current %cover of the algae)
        % col 2 gives the new areal cover of Lobophora in a 2500cm2 cell
        % for a sparse lobophora vegetation - this could be based on Henk's data
        % but unclear what the relationship is (this is not 0.1666*x^2+0.081959*x+0.0665)
        % col 3 is the same but for a dense lobophora vegetation, based on de Ruyter VS
        % Lobophora cover (%) is given by = 32.027*log(x) - 0.3847 , where x is months from 0 to 5

        b = 2; % the antagonistic algae is then Dictyota

        algal_new_cm2 = round(algal_dynamics(algal_pct+1, 2) + ...
            algal_env_prop.*(algal_dynamics(algal_pct+1, 3) - ...
            algal_dynamics(algal_pct+1, 2))) ;
        % Note YM: to what I know, this hasn't been published. In earlier
        % version of the code this is referred to the "productivity paper"
        % which to what I know was the first implementation of this.
        % Estimates from sparse and dense vegetation are derived from logistic growth
        % with unknown parameters. Max cover is 100%.
        % Note YM2: wave exposure does not affect Lobophora but dictyota.
        
        
        % Now scale the predicted algal cover to the current cell area
        % This is because algal covers in dict_and_lob_for_model.mat are
        % given for a 2500 cm2 cell
        algal_new_cm2 = round(cell_area_cm2.*algal_new_cm2/2500);
        
    case 2 % the algae to be processed is Dictyota
        % algal_dynamics here corresponds to 'dic_dynamics' (101x3) uploaded
        % from the data file 'dict_and_lob_for_model.mat'.
        % This lookup table was built from Henk's caged data from Belize
        % assuming the cover (%) of Dictyota (exposed site) is given by:
        % 14.712*log(x))+33.286 ,  where x is months from 0 to 5
        % For non-exposed site growth rate is 43% that of the exposed site
        
        % col 1 is ignored (id of current %cover of the algae)
        % col 2 gives the new areal cover of Dictyota in a 2500cm2 cell
        % for a windward (productive) reef
        % col 3 is the same but for a leeward (less-productive)reef
        
        b = 3; % the antagonistic algae is then Lobophora
 
        algal_new_cm2 = round(algal_dynamics(algal_pct+1, 2) - ...
            (1-exposure)* (algal_dynamics(algal_pct+1, 2) - ...
            algal_dynamics(algal_pct+1, 3)) );   
        
        % Note YM: this is the new formulation scaled to wave exposure for
        % the 'mapping resilience' paper in Cons. Lett. 2013
        % OLD stuff: alternative method based on Henks data
        %         if is_leeward == 1
        %             algal_new_cm2 = round(algal_dynamics(algal_pct+1, 3));
        %         else
        %             algal_new_cm2 = round(algal_dynamics(algal_pct+1, 2));
        %         end
        
        % Now scale the predicted algal cover to the current cell area
        % This is because algal covers in dict_and_lob_for_model.mat are
        % given for a 2500 cm2 cell
        algal_new_cm2 = round(cell_area_cm2.*algal_new_cm2/2500);
       

end

algal_new_cm2 = algal_new_cm2.*id_nongrazed ;

% this is the potential macroalgal growth if it really gets to become algal_new_cm2
pot_algal_growth_cm2 = (algal_new_cm2 - algal_cm2(:,a)).*id_nongrazed;

pot_algal_growth_cm2(pot_algal_growth_cm2<0)=0; % required because for Dictyota can be negative - not clear why...

% %%% TEST FOR NEGATIVE GROWTH OF DICTYOTA
% test = find(pot_algal_growth_cm2<0);
% if sum(test)~=0
%     a
%     current_area = cell_area_cm2(test,1)
%     current_cover = algal_cm2(test,a)
%     current_cover_pct = algal_pct(test,1)
%     cover_after_growth = algal_new_cm2(test,1)
%     potential = pot_algal_growth_cm2(test,1)
%    stop
% end      
%%%

% coral competition reduces potential growth
pot_algal_growth_cm2 = round(pot_algal_growth_cm2.*(1 - coral_reduce_macrogrowth)) ;
% this is how much the macroalgae will become after coral competition is accounted for
algal_new_cm2 = sparse((algal_cm2(:,a) + pot_algal_growth_cm2).*id_nongrazed) ;

% can't be bigger than the cell area 
algal_new_cm2(algal_new_cm2 > cell_area_cm2) = cell_area_cm2(algal_new_cm2 > cell_area_cm2) ;
% algal_new_cm2(algal_new_cm2 > cell_area_cm2-total_coral_cm2) = cell_area_cm2-total_coral_cm2 ;
% if it grew into the turf area this is how much turf space would be left
turf_space_avail_cm2 = sparse(algal_cm2(:,1).*id_nongrazed - pot_algal_growth_cm2) ; %% JUSQUE LA TOUT VA BIEN!!!!!!!!!!!!!!!

% Now we need to identify the cells where there is enough space to grow over turf
id_space = id_nongrazed ; % id of cells not previously grazed and containing turf
id_space(pot_algal_growth_cm2 <= 0) = 0;
id_nospace = id_space ;
% disp('%%% cells with enough space for growing:')
id_space(turf_space_avail_cm2 < 0) = 0 ; % exclude cells with no available space
% disp('%%% cells with too much algae so that colocation will occur:')
id_nospace(turf_space_avail_cm2 >= 0) = 0 ; % exclude cells with available space

% 1) First process the cells where there is enough turf for growing
% Reduce the cover of turf
algal_cm2(id_space==1,1) = algal_cm2(id_space==1,1) - pot_algal_growth_cm2(id_space==1);
% write the new cover of macroalgae "a" 
algal_cm2(id_space==1,a) = algal_new_cm2(id_space==1) ;

% 2) Then process the cells where "a" will overtop "b" because not enough turf available
% Amount left after overgrowing available turf
lobdict_cm2(id_nospace==1) = pot_algal_growth_cm2(id_nospace==1) - algal_cm2(id_nospace==1,1);
% macroalgal "a" overgrows the turf that is there
algal_cm2(id_nospace==1,a) = algal_cm2(id_nospace==1,a) + algal_cm2(id_nospace==1,1);
% no more turf
algal_cm2(id_nospace==1,1) = 0;

% Two cases arise: a does overtop b completely or partially 
% 1)if overgrows b but b remains
test = (algal_cm2(:,b) - colocation_cm2).*id_nospace;
algal_cm2(lobdict_cm2 <= test,a) = algal_cm2(lobdict_cm2 <= test,a) + lobdict_cm2(lobdict_cm2 <= test);

% 2)otherwise limit further expansion of a to area of b
algal_cm2(lobdict_cm2 > test,a) = algal_cm2(lobdict_cm2 > test,a) + algal_cm2(lobdict_cm2 > test,b) ...
    - colocation_cm2(lobdict_cm2 > test);


% algal_cm2(id2==1,a) = algal_cm2(id2==1,a) + algal_cm2(id2==1,b) ...
%     - colocation_cm2(id2==1);
% coloc_b2 = colocation_cm2.*id2
% ensure that lobdict_cm2 now represents the amount of underlying b
% lobdict_cm2(id2==1) = algal_cm2(id2==1,b) - colocation_cm2(id2==1);
% lobdict_cm2 = lobdict_cm2 + (algal_cm2(:,b) - colocation_cm2).*id2;
% lobdict_cm2 = lobdict_cm2 + (algal_cm2(:,b).*id2 - colocation_cm2.*id2);
% lobdict_cm2(id2==1) = lobdict_cm2(id2==1) + colocation_cm2(id2==1);
