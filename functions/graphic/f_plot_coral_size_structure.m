function [M]=f_plot_coral_size_structure(RESULT, META, nb_time_steps, nb_coral_types)

years = 0:nb_time_steps;

figure
whitebg([0 .5 .6])

% RECRUIT=cat(1,RESULT.recruit_count)
JUV=cat(1,RESULT.juv_count);
ADOL=cat(1,RESULT.adol_count);
ADULT=cat(1,RESULT.adult_count);

size_juv = 0:META.size_bins(1):META.adol_size ;
size_adol = 0:META.size_bins(2):META.adult_size ;
size_adult = 0:META.size_bins(3):META.max_size ;

% M_recruit=zeros(META.nb_time_steps+1,META.nb_coral_types)
M_juv=zeros(nb_time_steps+1,nb_coral_types,length(size_juv));
M_adol=zeros(nb_time_steps+1,nb_coral_types,length(size_adol));
M_adult=zeros(nb_time_steps+1,nb_coral_types,length(size_adult));

for c=1:nb_coral_types

% 	M_recruit(:,c) = squeeze(mean(RECRUIT(:,:,c),1))
	M_juv(:,c,:) = squeeze(mean(JUV(:,:,c,:),1))
	M_adol(:,c,:) = squeeze(mean(ADOL(:,:,c,:),1))
	M_adult(:,c,:) = squeeze(mean(ADULT(:,:,c,:),1))
 subplot(1,4,c)
 plot(size_adol,squeeze(M_adol(60,c,:)))

end




% %%%% 1) Plot all the trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(years, coral, 'white')
% axis([years(1) years(nb_time_steps) 0 100])
% xlabel('years','FontSize',12)
% ylabel('Total coral cover (%)','FontSize',12)
% title('ReefMod 3.1','FontSize',16)
% 
% %%%% 2) Determine mean trajectory (coral cover in %) %%%%%%%%%%%%%%
% M = mean(coral,1);
% hold on
% plot(years, M,'blue','LineWidth',3)


