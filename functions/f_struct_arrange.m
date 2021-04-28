%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created May 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coral] = f_struct_arrange(coral, nb_coral_types, doing_3D)

% Re-arrange the matrix coral_cm2 for optimization
% We try to limit the number of columns by moving colonies back to first colums 
% if space is available, then delete the last column if all colonies have been displaced
% Note a displaced colony stay in its cell (same row)

% Probably not worth doing this at every time step

for c = 1:nb_coral_types
    
    [rows cols pages] = size(coral(c).cover_cm2);
    a = [1:rows]' ;
    R = a(:,ones(cols,1)) ;
    [temp_coral_cm2, I] = sort(coral(c).cover_cm2,2,'descend') ;
    nIdx = R + (I-1)*rows ;
    % then re-organize the matrix of ID and clades accordingly
    temp_colony_ID = coral(c).colony_ID(nIdx);
    temp_clade = coral(c).clade(nIdx);
    % check if the last column(s) has zeros
    id_col = spones(temp_coral_cm2);
    id_sum=sum(id_col,1);
    % then delete last column(s)
    temp_coral_cm2(:,id_sum==0)=[] ;
    temp_colony_ID(:,id_sum==0)=[] ;
    temp_clade(:,id_sum==0)=[] ;
    
    % Store the new (optimized) coral cover
    coral(c).cover_cm2 = sparse(temp_coral_cm2) ;
    coral(c).colony_ID = sparse(temp_colony_ID) ;
    coral(c).clade = sparse(temp_clade) ;
    
    if doing_3D == 1 % then re-arrange the other matrices accordingly
        
        temp_surface_cm2 = coral(c).surface_cm2(nIdx);
        temp_volume_cm3 = coral(c).volume_cm3(nIdx);
        
        temp_surface_cm2(:,id_sum==0)=[] ;
        temp_volume_cm3(:,id_sum==0)=[] ;
        
        coral(c).surface_cm2 = sparse(temp_surface_cm2) ;
        coral(c).volume_cm3 = sparse(temp_volume_cm3) ;
        
    end
    
end
