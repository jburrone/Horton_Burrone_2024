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
%%%% Usage
%%%% [Results,RHO]=RhoBinCalc_SimBr(age,Size_Bins,N_bins,Density);
%% Designed to be called from Main_RhoBin_SimBr fuction
%%%%%%%%% Requires GenerateSimBranches function
function [Results,RHO]=RhoBinCalc_SimBr(age,Size_Bins,N_bins,Density)

N_Repeats_1=1000;
N_Repeats_2=50;
Results=zeros(1,3);
RHO = zeros(N_Repeats_1,1);
RHOTemp = zeros(N_Repeats_2,1);
for i=1:N_Repeats_1
    [BranchSynth]=GenerateSimBranches(age);
    N_branches=size(BranchSynth,2);
    for r=1:N_Repeats_2
        pool_Exc=[];
        pool_Inh=[];
        for b=1:N_branches
            % Get branch size
            Lenght=BranchSynth(b).lenght_um;
        % Get distances
            inhib_distances = BranchSynth(b).loc_Inh_um(:,1);
            spine_distances = BranchSynth(b).loc_Exc_um(:,1);
        % Get sizes
            inhib_sizes = BranchSynth(b).size_Inh(:,1);
            spine_sizes = BranchSynth(b).size_Exc(:,1);
            if (N_bins>2*Lenght/Size_Bins)
                N_bins=ceil(2*Lenght/Size_Bins);
            end 
            for bin = 1:N_bins
            %Generate a random start position ensuring the bin is
            %completely enclosed within the measured branch
                LowerBound=floor((Lenght-Size_Bins)*rand);
                UpperBound=LowerBound+Size_Bins;
                Spineindx = find( spine_distances <= UpperBound & spine_distances > LowerBound);
                BranchSynth(b).SpineDens_Bins(bin,1) = size((Spineindx),1) /Size_Bins;
                Spine_sizes = spine_sizes(Spineindx);
                BranchSynth(b).SumExc_Bin(bin,1)=sum(Spine_sizes)/Size_Bins;
                Shaftindx = find( inhib_distances <= UpperBound & inhib_distances > LowerBound);
                BranchSynth(b).InhibDens_Bins(bin,1) = size((Shaftindx),1) /Size_Bins;
                Shaft_sizes = inhib_sizes(Shaftindx);
                BranchSynth(b).SumInh_Bin(bin,1)=sum(Shaft_sizes)/Size_Bins;
            end
            Spine_BinDens=BranchSynth(b).SpineDens_Bins(:,1).';           
            Shaft_BinDens=BranchSynth(b).InhibDens_Bins(:,1).';
            Spine_BinSum=BranchSynth(b).SumExc_Bin(:,1).';           
            Shaft_BinSum=BranchSynth(b).SumInh_Bin(:,1).';
            if (Density==1)
                pool_Exc=[pool_Exc, Spine_BinDens];
                pool_Inh=[pool_Inh,Shaft_BinDens];
            elseif (Density==0)
                pool_Exc=[pool_Exc, Spine_BinSum];
                pool_Inh=[pool_Inh,Shaft_BinSum];
            else
                print "Choose 1 for Density and 0 for Sum"
            end
        end
        [RHOTemp(r,1),~] = corr( pool_Inh.', pool_Exc.', 'Type', 'Spearman');
    end
    RHO(i,1)=mean(RHOTemp(:,1));
end
Results(1,1)=mean(RHO);
Results(1,2)=std(RHO);
Results(1,3)=median(RHO);
end