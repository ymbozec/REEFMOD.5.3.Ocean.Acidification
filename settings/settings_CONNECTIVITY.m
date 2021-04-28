%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Mar 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%_________________________________________________________________________________________
%
%       CONNECTIVITY SETTINGS
%_________________________________________________________________________________________

% Make sure you have a matrix with equal number of columns to number of reefs
% See previous REEFMOD version for the use of other connectivity matrices

%%%% MATRIX 2x2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('data/Connectivity_Matrix.mat'); % produces CONNECTIVITY.spawner [2x2 double] and CONNECTIVITY.brooder [2x2 double]

%%%% MATRIX 3202x3202 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% load('data/MA_2006_Norm.mat'); % connectivity matrix 3202x3202 for Mont.
% annularis & P. astreoides from Holstein & Paris (see Readme Connectivity Matrices.doc).

%%%% MATRIX 3156x3156 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% load('data/Connectivity_3156.mat'); 
% % connectivity matrix (3156 reefs - excluding Gulf of Mexico) for Montastraea annularis
% % and Porites astreoides from Holstein & Paris (see Readme Connectivity Matrices.doc).
% n= META.nb_reefs;
% conn_mtrx_PA = conn_mtrx_PA(1:n,1:n) ; % TEMP: limit the size of the matrix
% conn_mtrx_MA = conn_mtrx_MA(1:n,1:n) ; % TEMP: limit the size of the matrix
% 
% META.connectivity(1).matrix = sparse(conn_mtrx_PA) ; % using Porites astr. connectivity for BROODER MASS RES
% META.connectivity(2).matrix = sparse(conn_mtrx_MA) ; % using Mont. ann. connectivity for SPAWNER MASS RES
% META.connectivity(3).matrix = sparse(conn_mtrx_PA) ; % using Porites astr. connectivity for BROODER MASS VUL
% META.connectivity(4).matrix = sparse(conn_mtrx_MA) ; % using Mont. ann. connectivity for SPAWNER MASS VUL
% 
% clear conn_mtrx_PA conn_mtrx_MA

%%%% SMOOTHED MATRIX 2441 x 2441 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
load('data/YM_connectivity.mat'); 
n = META.nb_reefs;
META.connectivity(1).matrix = sparse(Connect_mtrx_PA(1:n,1:n)) ; % using Porites astr. connectivity for BROODER MASS RES
META.connectivity(2).matrix = sparse(Connect_mtrx_MA(1:n,1:n)) ; % using Mont. ann. connectivity for SPAWNER MASS RES
META.connectivity_area_reef = Connect_habitat(1:n,2) ; % select proportion of polygon area covered by reefs
META.connectivity_area_habitat = Connect_habitat(1:n,3) ; % select proportion of polygon area covered by suitable habitat (wave exposure)
META.connectivity_polygons = Connect_polygons(1:n,:) ; % coordinates and region of the 2441 reef polygons
META.connectivity_region = Connect_regions(1:n,1);
clear Connect_habitat Connect_polygons Connect_mtrx_PA Connect_mtrx_MA Connect_regions
