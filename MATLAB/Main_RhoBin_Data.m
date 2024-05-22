%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Guilherme Neves  2024
%% guilherme.neves@kcl.ac.uk
% Juan Burrone Lab
% MRC Center for NeuroDevelopmental Disorders
% IoPPN, King's College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Usage [Summary,Distribution]=Main_RhoBin_Data(dataStruct,0)
% dataStruct represents the data Structure loaded using LoadData Function
% Density decides if you want to run the analysis using density (1) or 
% Cumulative Size (0) as a synapse metric
%%% Runs Spearmann rank correlation on randomly picked regions of size
%%% (BinSize) from dendritic branches
%%%% Summary Returns a matrix with columns: 1-Bin Size (microns)
%   2 - Mean Rho, 3 - Standard Deviation Rho 4 - Median Rho (Named Summary)
% And a separate matrix (Named Distribution) with 1 column for each Bin Size defined in the Bin
% Size array (default 5). Each Row represents a separate Repeat randomly
% generated start bin regions throuought the dendrite
% Requires Auxiliary functions RhoBinning_Data, LoadData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Summary,Distribution]=Main_RhoBin_Data(dataStruct,Density)


N_repeats=1000;
BinSize=[5 10 15 20 25];
N_Sizes=size(BinSize,2);
Summary=zeros(N_Sizes,4);
Distribution=zeros(N_repeats,N_Sizes);
for s=1:N_Sizes
    Size_Bins=BinSize(s);
    Summary(s,1) = Size_Bins;
    [Summary(s,2:4),Distribution(:,s)]  = RhoBinning_Data(dataStruct,Size_Bins,10,0,Density);
end

end
