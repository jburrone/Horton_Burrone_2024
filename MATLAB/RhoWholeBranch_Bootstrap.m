%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Guilherme Neves   2024
%% guilherme.neves@kcl.ac.uk
% Juan Burrone Lab
% MRC Center for NeuroDevelopmental Disorders
% IoPPN, King's College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Summary,Distribution]=RhoWholeBranch_Bootstrap(dataStruct,Density,N_repeats)
%%%% Usage [Summary,Distribution]=RhoWholeBranch_Bootstrap(dataStruct,0,1000);
% DataStruct is obtained using the LoadData function
% Density decides if you want to run the analysis using density (1) or 
% Cumulative Size (0) as a synapse metric
%%%% Summary Returns a vector with fields 1 - Mean Rho, 2 - Median Rho
%   3 - Standard Deviation Rho 
% And a separate matrix (Named Distribution) where each Row represents a separate
% random re-sampling of the branches (always keeping the total branch
% number constant. The value reported is the Rho value of the Spearmann
% Rank Correlation between inhibitory and excitatory synapse metrics
% Control the number of resampling using the N_repeats variable
%% Slice dataStruct regarding type of synapse
dataStruct_spine = dataStruct(strcmp({dataStruct.Synapse}, 'spines'));
dataStruct_shaft = dataStruct(strcmp({dataStruct.Synapse}, 'shaft'));

Nbranches = size(dataStruct_shaft, 1); % This assumes there are equal number of branches in exci and inhi data
Inhib = struct('Distances',[], 'Sizes', [], 'Totdensity', [], 'TotSum', []);
    
Spines = struct('Distances',[], 'Sizes',[], 'Totdensity', [],'TotSum', []);
    
%N_repeats=10000;
Distribution=zeros(N_repeats,1);
Summary=zeros(1,3);
for n = 1:N_repeats
    bt=ceil(Nbranches*(rand(1,Nbranches)));
%% Populate structs
    pool_Exc=[];
    pool_Inh=[];
    for i = 1:Nbranches
        b=bt(1,i);
        Size_Branch=dataStruct_shaft(b).Data(1,3);
        Inhib(b).Distances = dataStruct_shaft(b).Data(:,1);
        Inhib(b).Sizes = dataStruct_shaft(b).Data(:,2);
        N_Shafts = size(Inhib(b).Distances,1);    
        Inhib(b).Totdensity=N_Shafts/Size_Branch;
        Inhib(b).TotSum=sum(Inhib(b).Sizes)/Size_Branch;
        Spines(b).Distances = dataStruct_spine(b).Data(:,1);
        Spines(b).Sizes = dataStruct_spine(b).Data(:,2);
        N_Spines = size(Spines(b).Distances,1);
        Spines(b).Totdensity=N_Spines/Size_Branch;
        Spines(b).TotSum=sum(Spines(b).Sizes)/Size_Branch;
        if (Density==1)
            Spine_DensSum=Spines(b).Totdensity(:,1).';
            pool_Exc=[pool_Exc, Spine_DensSum];
            Shaft_DensSum=Inhib(b).Totdensity(:,1).';
            pool_Inh=[pool_Inh,Shaft_DensSum];
        elseif (Density==0)
            Spine_Sum=Spines(b).TotSum(:,1).';
            pool_Exc=[pool_Exc, Spine_Sum];
            Shaft_Sum=Inhib(b).TotSum(:,1).';
            pool_Inh=[pool_Inh,Shaft_Sum];
        else
            print "Choose 1 for density or 0 for sum"
            return
        end    
    end        
    [Distribution(n,1), ~] = corr( pool_Inh.', pool_Exc.', 'Type', 'Spearman');
end

Summary(1,1)=mean(Distribution(:,1));
Summary(1,2)=median(Distribution(:,1));
Summary(1,3)=std(Distribution(:,1));

end

