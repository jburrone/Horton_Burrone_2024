%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LoadData
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
%%  usage [dataStruct, ~] = LoadData;
function [dataStruct, error_flag] = LoadData()
%% Set assumptions on file format
data_delimiter = ',';
data_Nheaderlines = 1;
fname_types = 'Basal|Oblique|spines|shaft';
error_flag = false;
max_numbers = 4;
max_types = 2;

% File name rules:
% - 4 ordered numbers (max.): age, brain, cell and branch
% - 2 ordered key-words: 'Basal' or 'Oblique' AND 'spines' or 'shaft'

% Valid examples:
% Pn14_Br01_Ce030_Basal_br001_shaft.csv
% PostNatal_14_Brain_01_Cell_030_oblique_branch_001_shaft.csv
% 14_01_030_001_Oblique_shaft.csv

%% Locate files
data_folder = uigetdir();


%% List all files found
all_flist = dir([data_folder, '\*.csv']);
all_Nfiles = numel(all_flist);

% Error?
if all_Nfiles == 0
    error_flag = true;
    disp(['Error: No .csv files were found in ', data_folder])
    return
end
%% Get file name descriptors
number_matrix = zeros(all_Nfiles,4);
type_cell = cell(all_Nfiles, 2);

for f = 1:all_Nfiles
    % Get file name
    fname = all_flist(f).name;
    fnumbers = regexp(fname,'\d*','match');
    ftypes = regexp(fname,fname_types,'match');
    
    % Error?
%     if size(fnumbers,2)~=4 || size(ftypes,2)~=2
%         error_flag = true;
%         disp('Error: File name format is unexpected.')
%         return
%     end

    for n = 1:size(fnumbers,2)
        number_matrix(f,max_numbers-size(fnumbers,2)+n) = str2double(fnumbers{n});  % e.g.: if 2 numbers, cell and branch;
    end                                                                             % if 3, brain, cell and branch; etc.
    
    for n = 1:size(ftypes,2)
        type_cell{f,max_types-size(ftypes,2)+n} = ftypes{n};
    end
    
end

age_numbers = number_matrix(:,1);
brain_numbers = number_matrix(:,2);
cell_numbers = number_matrix(:,3);
branch_numbers = number_matrix(:,4);
region_types = type_cell(:,1);
synapse_types = type_cell(:,2);

%% Request files to analyse based on descriptors
%  Age
chosen_age_idxs = ones(all_Nfiles,1);

% Brain
chosen_brain_idxs = ones(all_Nfiles,1);

% Cell
chosen_cell_idxs = ones(all_Nfiles,1);

% Region
chosen_region_idxs = ones(all_Nfiles,1);

% Branch
chosen_branch_idxs = ones(all_Nfiles,1);

chosen_idxs = chosen_age_idxs & chosen_brain_idxs & chosen_cell_idxs & chosen_region_idxs & chosen_branch_idxs;
Nfiles = sum(chosen_idxs(:) == 1);


%% Create struct array, with each row for each file
%  Note: cell(Nfiles,1) defines how many rows the arrays will have.
dataStruct = struct('Age', cell(Nfiles,1), 'Brain', [],'Cell', [],'Region',[],'Branch',[],'Synapse',[],'Data', []);

%% Import data from files
chosen_ages = age_numbers(chosen_idxs);
chosen_brains = brain_numbers(chosen_idxs);
chosen_cells = cell_numbers(chosen_idxs);
chosen_regions = region_types(chosen_idxs);
chosen_branches = branch_numbers(chosen_idxs);
chosen_synapses = synapse_types(chosen_idxs);

chosen_flist = all_flist(chosen_idxs);

for f = 1:Nfiles
    % Populate struct
    dataStruct(f).Age = chosen_ages(f);
    dataStruct(f).Brain = chosen_brains(f);
    dataStruct(f).Cell = chosen_cells(f);
    dataStruct(f).Region = chosen_regions{f};
    dataStruct(f).Branch = chosen_branches(f);
    dataStruct(f).Synapse = chosen_synapses{f};
    
    fname = chosen_flist(f).name;
    temp = importdata([data_folder,'\',fname],data_delimiter,data_Nheaderlines);
    data = sortrows(temp.data);   % Sort data by ascending distances
    dataStruct(f).Data = data;
end

%% Correct missing branch lengths
for f = 1:Nfiles
    % Check if branch length is missing
    data = dataStruct(f).Data;
    
    if size(data,2) == 2    % and not 3
        % Select this row and row with complementary synapse type
        age = dataStruct(f).Age;
        brain = dataStruct(f).Brain;
        cel = dataStruct(f).Cell;
        region = dataStruct(f).Region;
        branch = dataStruct(f).Branch;
        synapse = dataStruct(f).Synapse;
        
        dataStruct_2syn = dataStruct([dataStruct.Age]==age &  [dataStruct.Brain]==brain & ...
            [dataStruct.Cell]==cel & strcmp({dataStruct.Region}, region) & ...
            [dataStruct.Branch]==branch);
        
        this_synapse_logical = strcmp({dataStruct_2syn.Synapse}, synapse);
        other_synapse_logical = ~this_synapse_logical;
        
        % Error?
        if size(dataStruct_2syn(other_synapse_logical).Data,2)~=3
            error_flag = true;
            disp('Error: No branch length was found.')
            return
        end
        
        % Get branch length from row of complementary synapse type
        branch_length = dataStruct_2syn(other_synapse_logical).Data(1,3);
        
        % Insert branch length in this row
        new_data = [data, branch_length.*ones(size(data,1),1)];
        dataStruct(f).Data = new_data;
    end
end
    
%% Garantee struct is sorted by Cell, then Region, then Branch, then type of Synapse
[~, ind] = sortrows([{dataStruct.Age}', {dataStruct.Brain}',{dataStruct.Cell}', {dataStruct.Region}', {dataStruct.Branch}', {dataStruct.Synapse}']);
dataStruct = dataStruct(ind);

% Error?
if isempty(dataStruct)
    error_flag = true;
    disp('Error: No data was selected.')
    return
end

end