%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mortalities = f_generate_bleaching_mortalities(predicted_DHWs)

%%%% YMB: can be optimized (slow simulations when multiple reefs)

%  Function to generate whole-colony mortality values
%  for each year from 2010 to 2100 for the current
%  location. Whole-colony mortality is dependent upon
%  the predicted number of degree heating weeks (DHWs). 
%  
%  If predicted DHW is less than 17, use window based
%  on 1 DHW either side of predicted value to generate 
%  relative cumulative frequency distribution. 
%  Mortality is drawn from this distribution.
% 
%  For higher DHWs (for which we have no data), 
%  mortality drawn from normal distribution with mean 
%  from regression line and SD calculated from average 
%  of SDs in windows.

% For each year we have a predicted DHW value. Use a window associated
% with (centred on??) this value.
num_windows = length(predicted_DHWs);

%For the x-axis we have min, mid and max dhw for each window
min_dhw = NaN(num_windows,1);
mid_dhw = predicted_DHWs;
max_dhw = NaN(num_windows,1);

for i = 1:num_windows
    %Make sure that we don't define a window that has negative DHWs
    if mid_dhw(i) < 1
        min_dhw(i) = 0;
    else
        min_dhw(i) = mid_dhw(i) - 1;
    end
    max_dhw(i) = mid_dhw(i) + 1;         
end 


%%%%%% YMB
% mid_dhw = predicted_DHWs;
% min_dhw = mid_dhw - 1;
% max_dhw = mid_dhw + 1;
% %Make sure that we don't define a window that has negative DHWs
% min_dhw(mid_dhw < 1) = 0;
%%%%%%

% Next, load the DHW and corresponding mortality data
% col 3 = total stress (dhw),   col 2 = prop coral cover mortality
load('data/DHW_dist_window_data.mat');
MORTDATA = DHW_dist_window_data(:,2);
DHWDATA = DHW_dist_window_data(:,3);

% Create an output mortalities for each window
mortalities = NaN(num_windows,1);

for win = 1:num_windows
    if mid_dhw(win) <=16    
        %CASE 1: use the actual data
        
        % get all data for the current window
        I = find( DHWDATA >= min_dhw(win) & DHWDATA <= max_dhw(win) );        
        numElementsInWindow = length(I);
        % there are numElementsInWindow DHW and corresponding MORT records
        
        % sort the mortality data into ascending order
        asc = sort(MORTDATA(I));
        diff = max(MORTDATA(I)) - min(MORTDATA(I));
        % we will sort into 10 bins for the current window
        freq = NaN(10,1);     % frequency for each of the 10 bins
        relfreq = NaN(10,1);  % relative freq (freq / numElementsInWindow)
        min_mort = NaN(10,1);  % min mort for bin i
        max_mort = NaN(10,1);  % max mort for bin i
        mid_mort = NaN(10,1);  % mid mort for bin i
        counter = 0;
        for i=0:0.1:0.9 % create ten bins for the current window
            counter = counter + 1;
            min_mort(counter) = (asc(1)+(diff*i));
            max_mort(counter) = (asc(1)+(diff*(i+0.1)));
            mid_mort(counter) = min_mort(counter) + ((max_mort(counter) - min_mort(counter)) / 2);
            if i==0.9
                freq(counter) = length(find( (MORTDATA(I) >= min_mort(counter)) & (MORTDATA(I) <= max_mort(counter)) ));
            else
                freq(counter) = length(find( (MORTDATA(I) >= min_mort(counter)) & (MORTDATA(I) < max_mort(counter)) ));
            end
        end
        %Find the relative frequency of all the values
        relfreq =  freq ./ numElementsInWindow;
        
        % sum relfreq = 1
        % [sum(freq) numElementsInWindow] % equals the same

        % freq and relative freq contain the count (and relative count) for
        % each of the 10 bins on the y-axis for the current window
        % Next we need the cumulative relative frequency
        crf = cumsum(relfreq);
        
        % NOW - WE WILL CHOOSE A RANDOM NUMBER rnd BETWEEN 0 and 1
        % from our cumulative distribution (relfreq sums to 1) so the
        % final cumulative frequency (where all are summed) equals 1
        %                       _._
        %  1 |                 |   |
        %    |        _._      |   |    (Imagine the dots are connected)
        %    |    _._|   |     |   |
        %    |_._|   |   |     |   |
        %  0 .___|___|___|_____|___|
        %    | 1 | 2 | 3 | ... | 10|
        %
        %  If rnd equals one of the crf values then we can return the
        %  mid-point (mid_mort) value of the corresponding bin.
        %
        %  crf(end) = 1, so if rnd also equals 1 then the mid point of the
        %  final bin will be returned.
        %
        %  Otherwise we must interpolate by constructing a line between
        %  the two mid-points to get at the required mortality value
        %
        %  Note that rnd < crf(1) is a special case since there is no
        %  preceeding midpoint.  In this case we use the min_mort(1) point
        
        rnd = rand;
        I = find(crf==rnd);
        if ~isempty(I)
            % the random number falls exactly on a midpoint, therefore
            % return the mort for this bin as the answer for the current
            % window
            mortalities(win) = mid_mort(I(1));
        elseif rnd < crf(1)
            % special case: interpolate between min_mort(1) and mid_mort(1)
            % if rnd == 0 then min_mort(1) is returned
            % if rnd == crf(1) then mid_mort(1) is returned
            % Interpolate:
            % since the line between the two points is linear then if we
            % calc the proportion of the rnd distance (y component) of the
            % line then we can get the mort by calculating the same
            % proportion of the mort distance (x component) of the line
            %
            % What is rnd as a proportion of crf(1)?
            y_proportion = rnd/crf(1);
            
            x_component = mid_mort(1) - min_mort(1);
            delta_x = y_proportion * x_component;
            mortalities(win) = min_mort(1) + delta_x;            
        else
            % find the crf points between which the rnd number lies
            % interpolate between these points.            
            for i=1:10
                if crf(i) < rnd
                    lower_point = i;
                end
            end
            if lower_point==10
                upper_point = 10;
            else
                upper_point = lower_point+1;
            end
            
            if lower_point < 1 | upper_point > 10
                error('Unexpected error, points should not be <=1 or >10');
            end
            % Interpolate between lower_point and upper_point
            % Use the same method as above
            %
            % What is rnd - crf(lower_point) as a proportion of
            % crf(upper_point)-crf(lower_point)?
            y_proportion = (rnd - crf(lower_point)) / (crf(upper_point)-crf(lower_point));
            x_component = mid_mort(upper_point) - min_mort(lower_point);
            delta_x = y_proportion * x_component;
            mortalities(win) = min_mort(lower_point) + delta_x;            
        end
    else
        %CASE 2: use linear regression line for DHW > 16 
        %mean of stdevs used to fit a distribution around the
        %regression line randn(mean=regressionline, std=meanstds)
   
        %total_mortality_GIS_ISO_cut_joined4km_output_ANALYSIS.xls (in this
        %folder) and the 
        %stdev calculated previously in DHW_dist_from_window
        %(meanstds = 0.0249)

        meandhw = 0.0027*mid_dhw(win) - 0.0033;
        if meandhw < 0
            warning('problem here');
        end
        stdevdhw = 0.0249;
        mortalities(win) = (stdevdhw*randn)+meandhw;        
        if mortalities(win)<0
            mortalities(win) = 0;
        end
    end
end