function [M]=f_plot_mean_coral_cover(coral, nb_time_steps, herbivory)

years = (0:nb_time_steps)/2 ;

whitebg([0 .5 .6])

% set(gcf,'PaperUnits','centimeters')
% % This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
% xSize = 8; ySize = 12;

%%%% 1) Plot all the trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(years, coral, 'white')
axis([years(1) years(nb_time_steps) 0 50])
xlabel('years','FontSize',12)
ylabel('Total coral cover (%)','FontSize',12)
% title('ReefMod 3.2.1.1','FontSize',16)

%%%% 2) Determine mean trajectory (coral cover in %) %%%%%%%%%%%%%%
M = mean(coral,1);
hold on
% plot(years, M,'blue','LineWidth',3)

% text(5,90,['herbivory = ', num2str(herbivory)],'FontSize',24,'FontWeight','bold','Color','black')
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

