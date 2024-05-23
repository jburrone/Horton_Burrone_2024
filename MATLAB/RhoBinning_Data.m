function [RHO, RHO_temp]=RhoBinning_Data(dataStructure,Size_Bins,N_bins,Size_Ign,Density,N_repeats)

%% Slice dataStruct regarding type of synapse
dataStruct_spine = dataStructure(strcmp({dataStructure.Synapse}, 'spines'));
dataStruct_shaft = dataStructure(strcmp({dataStructure.Synapse}, 'shaft'));

Nbranches = size(dataStruct_shaft, 1); % This assumes there are equal number of branches in exci and inhi data
Inhib = struct('Distances',cell(Nbranches,1), 'Sizes', [], 'BranchSize', [], 'Dens_Bins', [], ...
    'Sum_Bin', [], 'Lower_Bound', [], 'Upper_Bound', []);
Spines = struct('Distances',cell(Nbranches,1), 'Sizes',[], 'Dens_Bins', [], ...
 'Sum_Bin', []);
%% Populate structs
%N_repeats=1000;
RHO_temp=zeros(N_repeats,1);
RHO=zeros(1,3);
for r=1:N_repeats
    pool_Exc=[];
    pool_Inh=[];
    for b = 1:Nbranches
        Size_Branch=dataStruct_shaft(b).Data(1,3);
        if (Size_Branch-2*Size_Ign>Size_Bins)
            Inhib(b).Distances = dataStruct_shaft(b).Data(:,1);
            Inhib(b).Sizes = dataStruct_shaft(b).Data(:,2);
            Inhib(b).BranchSize = dataStruct_shaft(b).Data(:,3);
            Spines(b).Distances = dataStruct_spine(b).Data(:,1);
            Spines(b).Sizes = dataStruct_spine(b).Data(:,2);
            Size_Branch_corr=Size_Branch-2*Size_Ign;
            if (N_bins>2*Size_Branch_corr/Size_Bins)
                N_bins=ceil(2*Size_Branch_corr/Size_Bins);
            end    
            for bin = 1:N_bins
                LowerBound=floor((Size_Branch_corr-Size_Bins)*rand)+Size_Ign;
                UpperBound=LowerBound+Size_Bins;
                Inhib(b).Lower_Bound(bin,1)=LowerBound;
                Inhib(b).Upper_Bound(bin,1)=UpperBound;
                Spineindx = find( Spines(b).Distances <= UpperBound & Spines(b).Distances > LowerBound);
                Spines(b).Dens_Bins(bin,1) = size((Spineindx),1) /Size_Bins;
                Spine_sizes = Spines(b).Sizes(Spineindx);
                Spines(b).Sum_Bin(bin,1)=sum(Spine_sizes)/Size_Bins;
                Shaftindx = find( Inhib(b).Distances <= UpperBound & Inhib(b).Distances > LowerBound);
                Inhib(b).Dens_Bins(bin,1) = size((Shaftindx),1) /Size_Bins;
                Inhib_sizes = Inhib(b).Sizes(Shaftindx);
                Inhib(b).Sum_Bin(bin,1)=sum(Inhib_sizes)/Size_Bins;
            end

            if (Density==1)
                Spine_DensSum=Spines(b).Dens_Bins(:,1).';
                pool_Exc=[pool_Exc, Spine_DensSum];
                Shaft_DensSum=Inhib(b).Dens_Bins(:,1).';
                pool_Inh=[pool_Inh,Shaft_DensSum];
            elseif (Density==0)
                Spine_BinSum=Spines(b).Sum_Bin(:,1).';
                pool_Exc=[pool_Exc, Spine_BinSum];
                Shaft_BinSum=Inhib(b).Sum_Bin(:,1).';
                pool_Inh=[pool_Inh,Shaft_BinSum];
            else
                print "Choose 1 for density or 0 for sum"
                return
            end
        end
    end
    X_forCorrel=pool_Inh.';
    Y_forCorrel=pool_Exc.';
    [RHO_temp(r,1),~] = corr( X_forCorrel, Y_forCorrel, 'Type', 'Spearman');
end

RHO(1,1)=mean(RHO_temp);
RHO(1,2)=std(RHO_temp);
RHO(1,3)=median(RHO_temp);
end

