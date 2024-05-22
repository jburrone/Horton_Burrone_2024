%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthetic Synapse (Inh/Exc) Spatial Distributions
% Paulo Aguiar
% pauloaguiar@i3s.up.pt
% Neuroengineering and Computational Neuroscience Lab
% i3S, University of Porto
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [branch]=GenerateSimBranches(age)
%%%% Usage
%%%% [branch]=GenerateSimBranches(age)
%% Designed to be called from Main_RhoBin_SimBr fuction

%% PARAMETERS for Branches matching the different datasets
%% P21vEM Use 22
if (age==22)
    lambda_Inh_um = 2.88;
    lambda_Exc_um = 0.51;

    branches_length_mean_um = 51.9;
    branches_length_std_um  = 7.5;

    sizes_Exc_mean = 0.25;
    sizes_Exc_std = 0.17;
    sizes_Inh_mean = 0.31;
    sizes_Inh_std = 0.21;

    N_branches=5;

%% P21
elseif (age==21)
    lambda_Inh_um = 2.88;
    lambda_Exc_um = 0.51;

    branches_length_mean_um = 131.7;
    branches_length_std_um  = 31;

    sizes_Exc_mean = 1.1;
    sizes_Exc_std = 0.5;
    sizes_Inh_mean = 1.03;
    sizes_Inh_std = 0.5;

    N_branches=21;
elseif (age==7)
%% P7
    lambda_Inh_um = 5.55;
    lambda_Exc_um = 2.4;
    
    branches_length_mean_um = 92.0;
    branches_length_std_um  = 33.1;

    sizes_Exc_mean = 1;
    sizes_Exc_std = 0.3;
    sizes_Inh_mean = 1.1;
    sizes_Inh_std = 0.6;

    N_branches=18;
elseif (age==10)
%% P10
    lambda_Inh_um = 4.78;
    lambda_Exc_um = 0.97;
    
    branches_length_mean_um = 72.8;
    branches_length_std_um  = 29.4;

    sizes_Exc_mean = 1.2;
    sizes_Exc_std = 0.6;
    sizes_Inh_mean = 1.2;
    sizes_Inh_std = 0.6;

    N_branches=13;
elseif (age==14)
%% P14
    lambda_Inh_um = 3.9;
    lambda_Exc_um = 1.13;
    
    branches_length_mean_um = 84;
    branches_length_std_um  = 34.5;

    sizes_Exc_mean = 1.5;
    sizes_Exc_std = 1.4;
    sizes_Inh_mean = 1.2;
    sizes_Inh_std = 0.9;

    N_branches=13;
else
    print "Please choose age 7 10 14 21 or 22 for EM data";
end
branch=struct('lenght_um',[],'loc_Inh_um',[], 'loc_Exc_um',[], 'N_Inh',[], 'N_Exc',[],...
                'D_Inh', [], 'D_Exc', [],'size_Inh',[],'size_Exc',[]);
%% GENERATE SYNTHETIC BRANCHES

b=1;
while b <=N_branches
    % toss branch lenght
    branch_L_um = branches_length_mean_um + branches_length_std_um * randn;
    % remove branches smaller than 40 microns
    if (branch_L_um > 40)
    % place INH
        n = round( branch_L_um / lambda_Inh_um );
    % toss synapse size
        size_Inh = sizes_Inh_mean + sizes_Inh_std * randn(1, 5*n);
        loc_Inh_um = exprnd( lambda_Inh_um, 1, 5*n ) + lambda_Inh_um * ones(1, 5*n);
        loc_Inh_um = cumsum( loc_Inh_um );
    % place EXC
        n = round( branch_L_um / lambda_Exc_um );
        size_Exc = sizes_Exc_mean + sizes_Exc_std * randn(1, 5*n);
        loc_Exc_um = exprnd( lambda_Exc_um, 1, 5*n ) + lambda_Exc_um * ones(1, 5*n);
        loc_Exc_um = cumsum( loc_Exc_um );
    % store everything
        branch(b).lenght_um  = branch_L_um;
        branch(b).loc_Inh_um = loc_Inh_um.';
        branch(b).size_Inh = size_Inh.';
        branch(b).size_Inh(size_Inh < 0.1)=0.1;
        branch(b).loc_Inh_um(loc_Inh_um > branch_L_um)=[];
        branch(b).size_Inh(loc_Inh_um > branch_L_um)=[];
        branch(b).loc_Exc_um = loc_Exc_um.';
        branch(b).size_Exc = size_Exc.';
        branch(b).size_Exc(size_Exc < 0.1)=0.1;
        branch(b).loc_Exc_um(loc_Exc_um > branch_L_um)=[];
        branch(b).size_Exc(loc_Exc_um > branch_L_um)=[];
        branch(b).N_Inh      = numel( branch(b).loc_Inh_um );
        branch(b).N_Exc      = numel( branch(b).loc_Exc_um );
        branch(b).D_Inh      = branch(b).N_Inh/branch_L_um;
        branch(b).D_Exc      = branch(b).N_Exc/branch_L_um;
        b=b+1;
    end

end
end


