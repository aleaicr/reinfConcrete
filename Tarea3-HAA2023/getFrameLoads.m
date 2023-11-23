function data = getFrameLoads(fileDir, frameIDs)
% This function reads a SAP2000 excel file and returns the array of the
% all the internal loads (N, Vy, Vx, My, Mx) for each frame in the frameIDs
% vector
%
% INPUTS
% fileName: Name of the file with the SAP2000 results (at least should have
%           the 'Element Froces - Frames' sheet
% frameIDs: ID number of each frame element 
%
% OUTPUTS
% data: Struct para cada frame ID que 
%
%
%

% Read table
frameTable = readtable(fileDir, 'Sheet', 'Element Forces - Frames');
nElems = length(frameIDs);

data = struct();

for i = 1: nElems
    data(i).frameID = frameIDs(i);
    data(i).frameLoads = frameTable(contains(frameTable.Text, frameIDs(i)));
end

end

