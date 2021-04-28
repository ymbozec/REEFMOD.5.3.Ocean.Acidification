% Y.-M. Bozec, MSEL, created Feb 2012.
% Last modified: 07/08/2012
%
% SETTING DATA FOR HURRICANES
%_________________________________________________________________________________________
%_________________________________________________________________________________________


%%%%%%%% Multiple-reefs scenarios %%%%%%%%%%%%%%%%%%%%%
% NEEDS DOUBLE CHECK
% load('Hurr_Hist.mat') % This is for the analysis of connectivity (Nick)
% % Includes Hurr_Max = 3156 (reefs) x 100 (years) which gives the Max hurricane category
% % (1-5) experienced by each reef for years 1909-2008. 0 means no hurricane.
% % Uses Ian Elliot Keim et al 2007 method

% for n=1:META.nb_reefs
%     REEF(n).hurricane_chronology = Hurr_Max(n,:);
% end

%%%%%%%% Single reefs scenarios %%%%%%%%%%%%%%%%%%%%%

% load('Hurr_Hist_336.mat') % This is for the analysis of resilience in Belize (Nick)
% Includes Hurr_Max = 336 (reefs) x 100 (years) which gives the Max hurricane category
% (1-5) experienced by each reef for years 1909-2008. 0 means no hurricane.
% Uses Ian Elliot Keim et al 2007 method

REEF.hurricane_chronology = zeros(1, META.nb_time_steps);

if META.random_hurricanes == 0   % if not random, pickup a hurricane regime (Caribbean) 
    
    load('Hurr_Hist.mat')
    nreef = 3100 ; % select a reef (3100 is Cozumel in Nick's reef dataset)

    REEF(1).hurricane_chronology(1,1:2:META.nb_time_steps) = Hurr_Max(nreef,1:META.nb_time_steps/2); % feed the chronology to the model
    
    clear nreef Hurr_Max
end



    
        
