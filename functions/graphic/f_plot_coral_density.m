function [M]=f_plot_coral_density(density, nb_time_steps)

years = 0:nb_time_steps ;

figure
whitebg([0 .5 .6])

% set(gcf,'PaperUnits','centimeters')
% % This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
% xSize = 8; ySize = 12;

for c = 1:META.nb_coral_types
    
    recruit = cat(1, RESULT.size_struct)
    mean_recruit(c,1:(nb_time_steps+1) = mean(
    coral = sum(cat(1,RESULT.coral_pct),3)

%%%% 1) Plot density of recruits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(years, coral, 'white')
axis([years(1) years(nb_time_steps) 0 100])
xlabel('years','FontSize',12)
ylabel('Total coral cover (%)','FontSize',12)
title('ReefMod 3.1','FontSize',16)

%%%% 2) Determine mean trajectory (coral cover in %) %%%%%%%%%%%%%%
M = mean(coral,1);
hold on
plot(years, M,'blue','LineWidth',3)

%%%% 3) Compute standard deviation for plotting the envelope %%%%%%
% SD = std(coral,1,1);
% L = 1.96*SD;
% figure
% confplot(0:META.nb_time_steps, M,L,L)
% axis([0 META.nb_time_steps 0 100])
% grid on
% title('REEFMOD.2.5.2')
% xlabel('Time steps (6 mo)')
% ylabel('Total coral cover (%)')

