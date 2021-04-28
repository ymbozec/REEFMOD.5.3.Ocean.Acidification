% Y.-M. Bozec, MSEL, created Feb 2012.
% Last modified: 07/08/2012
%
% BLEACHING AND CLIMATE CHANGE SETTINGS
%_________________________________________________________________________________________
%_________________________________________________________________________________________
%
% One reef scenarios
%_________________________________________________________________________________________

%%% CARIBBEAN AVERAGE (from GCM.m processed by Edwards, based on Halloran data)
load('data/RCP_use') ; % Gives SSTi, out26 and out85 for different RCP (Caribbean average)
% -> SSTi is the optimal SST value for coral calcification based on 1985-1991 mean
% -> out26: SST, DHM, Aragonite for every 6 monthly period, from 2006 to 2099 based on RCP2.6
% -> out85: SST, DHM, Aragonite for every 6 monthly period, from 2006 to 2099 based on RCP8.5
% 
% % Optimal SST for coral calcification
% CORAL.SST_OPT(1:4,1)=SSTi ; % use SSTi as the optimal temperature for every coral species
% CORAL.SST_OPT = [27.2 ; 24 ; 24 ; 23 ] ; % alternative inputs per species (see Pete_bioerosion2.m) 
% CORAL.relative_calci = [1 ; 1 ; 1 ; 1] ; % Relative calcification for each species
% CORAL.SD_relative_calci = [5.6 ; 3.35 ; 5.6 ; 3.35] ; %StDev (i.e. width of the calcification curve for each species)
% CORAL.SD_relative_calci = [3.35 ; 3.35 ; 3.35 ; 3.35] ;

% % Trajectory of future SST values; currently Caribbean average 2010-2100
% REEF.SST = out26(1,:); % Note you have to change out26 to out85  depending on greenhouse emissions
% 
% % Predicted number of degree heating weeks, for each year from 2010 to 2100.
% the predicted chronology starts in winter -> we exclude the first point
% because simulations starts in summer
REEF.predicted_DHWs = 4*out85(2,2:META.nb_time_steps+1); clear out26 % converts to dhws.
% REEF.predicted_DHWs = 4*out26(2,2:META.nb_time_steps+1); clear out85 % converts to dhws.


%_________________________________________________________________________________________
%
% Multiple-reefs scenarios
%_________________________________________________________________________________________

% load('data/Nick_connectivity/input2p6.mat')
% DHW = DHW2p6 ;
% % includes DHW2p6 = 3156 x 180 double matrix (Reefs x years 2010:0.5:2099.5 (every 6 months)).
% % Degree Heating Weeks for RCP2.6, aggressive greenhouse gas mitigation.
% REEF.predicted_DHWs2 = DHW2p6(1,:);

% load('data/Nick_connectivity/input8p5.mat');
% DHW = DHW8p5 ;
% % includes DHW8p5 = 3156 x 180 double matrix (Reefs x years 2010:0.5:2099.5 (every 6 months)).
% % Degree Heating Weeks for RCP8.5, business as usual greenhouse gas mitigation.
% % Year_DHW = 1 x 180 double array (years 2010:0.5:2099.5 (every 6 months))
% 
% for n=1:META.nb_reefs
%     REEF(n).predicted_DHWs = DHW(n,:) ; % select 2.6 or 8.5 above
% end

%%%%%%%% One-reef scenarios but from multiple reef data %%%%%%%%%%%%%%%%%%%%%
% load('input8p5_336.mat') % generates DHW8p5 (reef cells x time steps) and year (1 x time)
% clear year ; %(we dont need year)
% %  Predicted number of degree heating weeks, for each year from 2010 to 2100.
%  REEF.predicted_DHWs = DHW8p5(1,:);
