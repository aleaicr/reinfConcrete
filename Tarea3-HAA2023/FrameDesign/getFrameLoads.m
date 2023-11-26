function [data, frameTable] = getFrameLoads(fileDir, frameIDs, loadComb)
% This function reads a SAP2000 excel file and returns the a struct with 
% the array of all internal loads (N, Vy, Vx, My, Mx) for each frame in 
% the frameIDs vector
%
% INPUTS
% fileName: Name of the file with the SAP2000 results (at least should have
%           the 'Element Froces - Frames' sheet
% frameIDs: ID number of each frame element 
% loadComb: String with 'LRFD' or 'ASD' (It works also for the exact name)
%
% OUTPUTS
% data: Struct para cada frame ID que 
%
%
%

% Read table
frameTable = readtable(fileDir, 'Sheet', 'Element Forces - Frames');
frameTable = frameTable(contains(frameTable.OutputCase, loadComb), :);

% Transform vector to cell array
frameIDs_cell = arrayfun(@(x) {num2str(x)}, frameIDs);

% Init variables
data = struct();

% get frame loads struct for each frame
for i = 1:length(frameIDs)
    filteredTable = frameTable(ismember(frameTable.Frame, frameIDs_cell{i}), :);
    data(i).frameID = frameIDs(i);
    data(i).frameTable = filteredTable;
end
end
