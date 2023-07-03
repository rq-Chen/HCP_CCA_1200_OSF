%% summarize_movement - Get head movement summary
%
% Ruiqi Chen, 07/03/2023
%
% Output is a CSV file with the following columns: Subject, Scan, MovementRelativeRMSmean

HCP_dir = '/net/10.27.136.121/hcpdb/packages/unzip/HCP_1200';  % containing subject folders
outfile = '../data/HCP1200-Motion-Summaries-rfMRI.csv';

scans = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};

% Get subject list
sublist = dir(HCP_dir);
sublist = sublist(3:end);  % remove . and ..
sublist = {sublist.name};

% Get movement summary
Subject = {};
Scan = {};
MovementRelativeRMSmean = [];
for i = 1:length(sublist)
    sub = sublist{i};
    for j = 1:length(scans)
        scan = scans{j};
        try
            % Get movement summary
            movement = readmatrix(fullfile(HCP_dir, sub, 'MNINonLinear', 'Results', scan, 'Movement_RelativeRMS_mean.txt'));
            % Append to output
            Subject{end+1} = sub;
            Scan{end+1} = scan;
            MovementRelativeRMSmean(end+1) = movement;
        catch
            warning('Failed to read movement summary for subject %s scan %s', sub, scan);
        end
    end
end

% Write to CSV
T = table(Subject', Scan', MovementRelativeRMSmean');
T.Properties.VariableNames = {'Subject', 'Scan', 'MovementRelativeRMSmean'};
writetable(T, outfile);