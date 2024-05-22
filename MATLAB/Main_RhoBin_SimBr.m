
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Guilherme Neves & Paulo Aguiar   2024
%% guilherme.neves@kcl.ac.uk
% Juan Burrone Lab
% MRC Center for NeuroDevelopmental Disorders
% IoPPN, King's College London
%% pauloaguiar@i3s.up.pt
% Neuroengineering and Computational Neuroscience Lab
% i3S, University of Porto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Usage [Summary,Distribution]=Main_RhoBin_SimBr(7,0);
% age represents the version of the simulation you want to run
% Matched to the characteristics of 7,10,14,21 and EM (use age=22) datasets
% Density decides if you want to run the analysis using density (1) or 
% Cumulative Size (0) as a synapse metric
%%% Runs Spearmann rank correlation on randomly picked regions of size
%%% (BinSize) from randomly generated simulated branches
%%%% Summary Returns a matrix with columns: 1-Bin Size (microns)
%   2 - Mean Rho, 3 - Standard Deviation Rho 4 - Median Rho (Named Summary)
% And a separate matrix (Named Distribution) with 1 comlumn for each Bin Size defined in the Bin
% Size array (default 5). Each Row represents a separate Repeat with a randomly
% generated branch
% Requires Auxiliary functions RhoBinCalc_SimBr, GenerateSimBranches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Summary,Distribution]=Main_RhoBin_SimBr(age,Density)
%%%% Usage [Summary,Distribution]=Main_RhoBin_SimBr(7,0);
% age represents the version of the simulation you want to run
% Matched to the characteristics of 7,10,14,21 and EM (use age=22) datasets
% Density decides if you want to run the analysis using density (1) or 
% Cumulative Size (0) as a synapse metric
N_Repeats=1000;
N_bins=20;
BinSize=[5 10 15 20 25];
N_Sizes=size(BinSize,2);
Summary=zeros(N_Sizes,4);
Distribution=zeros(N_Repeats,N_Sizes);
for s=1:N_Sizes
    Size_Bins=BinSize(s);
    Summary(s,1) = Size_Bins;
    [Summary(s,2:4),Distribution(:,s)]  = RhoBinCalc_SimBr(age,Size_Bins,N_bins,Density);
end

end
