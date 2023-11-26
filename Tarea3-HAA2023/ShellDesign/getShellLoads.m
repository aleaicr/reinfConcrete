function [M11,M22,V13,V23] = getShellLoads(fileDir, loadComb, cNodes)
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
% Notes:
% *this function does not have frameIDs because we find maximum of all only
%

% Read table
ShellTable = readtable(fileDir, 'Sheet', 'Element Forces - Area Shells');
ShellTable_original = ShellTable;
ShellTable = ShellTable(contains(ShellTable.OutputCase, loadComb), :);      % {'ASD1'} -> contains              
ShellTable = ShellTable(strcmp(ShellTable.StepType, 'Max'), :);             % 'max' -> strcmp

% M11
M11 = struct();
M11.all = ShellTable.M11;
M11.max = max(M11.all);
M11.min = min(M11.all);
M11.max_OutputCase = ShellTable.OutputCase(ShellTable.M11 == max(ShellTable.M11));
M11.min_OutputCase = ShellTable.OutputCase(ShellTable.M11 == min(ShellTable.M11));

% M11
M22 = struct();
M22.all = ShellTable.M22;
M22.max = max(M22.all);
M22.min = min(M22.all);
M22.max_OutputCase = ShellTable.OutputCase(ShellTable.M22 == max(ShellTable.M22));
M22.min_OutputCase = ShellTable.OutputCase(ShellTable.M22 == min(ShellTable.M22));

% Filter again (not under columns)
validIndex = ~ismember(str2double(ShellTable.Joint), cNodes);
ShellTable_ = ShellTable;
ShellTable = ShellTable(validIndex, :);                                     % Quit frame-shell connection node.

% V13
V13 = struct();
V13.all = ShellTable.V13;
V13.max = max(abs(V13.all));
V13.max_OutputCase = ShellTable.OutputCase(ShellTable.V13 == max(abs(V13.all)));
if isequal(V13.max_OutputCase, cell(0,1))
    V13.max_OutputCase = ShellTable.OutputCase(ShellTable.V13 == -max(abs(V13.all)));
end
% V23
V23 = struct();
V23.all = ShellTable.V23;
V23.max = max(abs(V23.all));
V23.max_OutputCase = ShellTable.OutputCase(ShellTable.V23 == max(abs(V23.all)));
if isequal(V23.max_OutputCase, cell(0,1))
    V23.max_OutputCase = ShellTable.OutputCase(ShellTable.V23 == -max(abs(V23.all)));
end
end
