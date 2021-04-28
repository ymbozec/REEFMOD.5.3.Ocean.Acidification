%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 31/05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%__________________________________________________________________________
%
% REEFMOD RUN FILE
%__________________________________________________________________________
%__________________________________________________________________________

clear
PARAMETERS_DEFAULT
PARAMETERS_MOOREA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine here any parameter as needed
META.nb_simul = 10; % matches output of 38 replicates per year from field data
META.nb_time_steps = 16 ; % 1 step = 6 months
META.doing_coral_competition = 1 ;

META.max_colonies = 15 ; % use 1 to compare with Ortiz et al 2014 (= only 1 colony per 1m2 cell per species)

% Annual recruit densities observed on tiles in Moorea from 2010 to 2014
Recruit_density_Pocillo = [ 3.74 3.77 3.15 3.13 0.87 1.47 ];
Recruit_density_Acropora = [ 0.1 0 0.5 0.064 0.12 0.18 ];
Recruit_density_Montipora = [ 0.1 0.133 1 0.266 0.615 0.723 ];
Recruit_density_Porites = [ 0.20 0.35 0.64 0.20 0.34 0.65 ];
% scaling factor to calibrate recruitment rate in the model
scale_recruit_Pocillo = 1.15 ;
scale_recruit_Acropora = 0.6 ;
scale_recruit_Montipora = 1.65 ;
scale_recruit_Porites = 3.1 ;

% probability of coral recruits removed
CORAL.parrotfish_predation = 0.05 ;

% max coral size in cm2 in which parrotfish predation occurs
CORAL.threshold_predation_size = 5 ;

CORAL.initial_cover=0*CORAL.initial_cover; % no intitial coral cover
REEF.herbivory=1; % full herbivory
ALGAL.initial_cover=0*ALGAL.initial_cover; % no initial algal cover

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.doing_acute_disturbances = 0 ; % set to 0/1 to switch off/on
REEF.disturbance_event_proba = 1/2; % Frequency at which hurricanes occur in summer, if no prescribed hurricane scenario

mortality_cyclone = 1 - [11/25 ; 3/12 ; 3/7 ; 1 ] ; % observed mortalities after category 3 cyclone Wasa in 1991-1992 from Lamy et al. 2016
mortality_COTS = 1 - [12/19 ; 2/15 ; 0.2/0.3 ; 6.5/10 ] ; % observed mortalities due to COTS between 2006-2008 from LTER data
% Need assumption for COTS:
% 1- we simulate over 6 month the outcome of a full COTS outbreak
% 2- we assume the observed losses are relative (i.e. from 19% to 12% means
% a 37% loss)
% 3- while for Montipora the picture might be wrong because cover before COTS was negligible, the resulting loss is comparable to Pocillo and Porites
REEF.mortality_proba = [mortality_cyclone  mortality_COTS] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INITIALISATION ;

for simul = 1:META.nb_simul
    
    META.max_density_settlers = zeros(2,META.nb_time_steps+1);

    META.max_density_settlers(1,:) = scale_recruit_Pocillo*Recruit_density_Pocillo(randi(6,1,META.nb_time_steps+1));
    META.max_density_settlers(2,:) = scale_recruit_Acropora*Recruit_density_Acropora(randi(6,1,META.nb_time_steps+1));
    META.max_density_settlers(3,:) = scale_recruit_Montipora*Recruit_density_Montipora(randi(6,1,META.nb_time_steps+1));
    META.max_density_settlers(4,:) = scale_recruit_Porites*Recruit_density_Porites(randi(6,1,META.nb_time_steps+1));

    
    [RESULT(simul), RECORD] = f_runmodel(META, REEF, CORAL, ALGAL) ;
    
    META.nb_simul-simul;
    
end
coral_cover_tot = sum(cat(1,RESULT.coral_pct2D),3) ; % total coral cover for each simulation
coral_cover_per_taxa = cat(1,RESULT.coral_pct2D) ; % species coral cover for each simulation

M_coral_cover_tot = mean(coral_cover_tot, 1) ;
M_coral_cover_per_taxa = mean(coral_cover_per_taxa, 1) ;

algal_cover = cat(1,RESULT.algal_pct) ;
years = (0:META.nb_time_steps)/2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot coral cover

figure
whitebg([0 .5 .6])
subplot(1,2,1)
hold on
plot(years,  M_coral_cover_per_taxa(:,:,1), 'blue', 'LineWidth',3) %Pocillopora
plot(years,  M_coral_cover_per_taxa(:,:,2), 'red', 'LineWidth',3) %Acropora
plot(years,  M_coral_cover_per_taxa(:,:,3), 'green', 'LineWidth',3) %Montipora
plot(years,  M_coral_cover_per_taxa(:,:,4), 'black', 'LineWidth',3) %Porites
legend('Pocillopora','Acropora','Montipora','Porites')

axis([years(1) years(META.nb_time_steps) 0 90])
xlabel('years','FontSize',12)
ylabel('Species coral cover (%)','FontSize',12)    

subplot(1,2,2)
plot(years, coral_cover_tot, 'white')
axis([years(1) years(META.nb_time_steps) 0 90])
hold on
plot(years, M_coral_cover_tot, 'blue', 'LineWidth',3)
xlabel('years','FontSize',12)
ylabel('Total coral cover (%)','FontSize',12)    
