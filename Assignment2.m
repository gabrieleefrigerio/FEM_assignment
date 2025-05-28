clc
close all
clear
addpath('/Users/francescogiglio/Desktop/Dynamics/Assignment /Assignment2/Functions')

%% Assignment 2
%% Data
data = init_data();

%% Check length
check_all_element_lengths(data);

%% Input file
[file_name, nodes, elements] = inputA2(data);

%% Total mass 
total_mass = TotalMassComputation(nodes, elements, data);

%% Load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure_i(file_name);

%% Draw structure
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);

%% Assemble mass and stiffness matrices
[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);

%% Partition the matrices to get MFF and KFF
[MFF, MCF, MFC, MCC, KFF, KCF, KFC, KCC] = extract_matrices(M, K, ndof);

%% Compute natural frequencies and mode shapes
[freq_below_fmax, freq0, modes, omega] = compute_natural_frequencies(MFF, KFF, data.f_max);

%% Mode shape plot
plot_mode_shapes(freq0, modes, data.f_max, incidenze, l, gamma, posiz, idb, xy);

%% Define the damping matrix
[C, CFF, CCF] = damping_matrix(M, K, data.alpha, data.beta, ndof);

%% Frequency Response Function
[X, f, om] = compute_frequency_response(MFF, CFF, KFF, ndof, idb, data);

% Plot FRF
plot_frequency_response(X, f, idb);

%% FRF comparison
frf_comparison(MFF, CFF, KFF, modes, om, f, idb, ndof)

%% Considering the weight of the structure
[V_def_A, max_vertical_disp] = compute_static_deflection_due_to_gravity(M, KFF, idb, ndof, incidenze, l, gamma, posiz, xy)

%% Reducing maximum deflection of A (new structure)
new_structure(data, total_mass)