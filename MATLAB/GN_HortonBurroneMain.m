%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script loads and analysis data from inhib. and excit. synapses: 
% - distance along branch; 
% - intensity;
% - branch length;
%
% Authors: Gustavo Caetano & Paulo Aguiar
% Email: pauloaguiar@ineb.up.pt
% Neuroengineering and Computational Neuroscience Lab
% https://www.i3s.up.pt/neuroengineering-and-computational-neuroscience
% i3S/INEB
%
%   Imports Data into a DataStructure used by the Main_RhoBin_Data function
% The data is stored in .csv files in \ExampleData.

clear;
close all;
clc;

% Load Data
[dataStruct, ~] = LoadData;

%%%% Run Bootstrap Analysis for cumulative Fluorescence
%%%%%%%%% Decide the number of analysis repeats
%%% More Repeats is slower but less accurate
N_repeats=1000;

Density=0;
%%% If density analysis is required use 
% Density=1;

[BootstrapSummary,DistributionBootstrap]=RhoWholeBranch_Bootstrap(dataStruct,Density,N_repeats);
BootsrapFigure=figure;
x = [ones(1,N_repeats)];
swarmchart(x,DistributionBootstrap);

%%%% Run Binning Analysis for cumulative Fluorescence

[BinningSummary,DistributionBins]=Main_RhoBin_Data(dataStruct,Density,N_repeats);


BinningFigure=figure;
plot(BinningSummary(:,1),BinningSummary(:,2),'*r');
xlim([0,30]);

%%%% Run Binning Analysis with simulated branches
age=21;
[BinningSummary_Rand,DistributionBins_Rand]=Main_RhoBin_SimBr(age,Density,N_repeats);

hold on;
plot(BinningSummary_Rand(:,1),BinningSummary_Rand(:,2),'ok');
