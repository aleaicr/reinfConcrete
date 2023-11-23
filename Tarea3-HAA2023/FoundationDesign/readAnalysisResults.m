function [story, loads, pier, loads_selected] = readAnalysisResults(file_name, Piers, stories)
% Tarea 2 - Hormigón Armado Avanzado
% Departamento de Obras Civiles - Universidad Técnica Federico Santa María
% Alexis Contreras R. - Gabriel Ramos V.
%%
% Inputs
% file_Name: char or string with the name of the file (the direction)
% Piers: cell with strings of the piers to be selected
% stories: vector of the 
% Outputs
% 
%
% Notes
%
%
% Read data from the Excel file with tab delimiters
data = readtable(file_name);

% Extract the "PIER" column
pier = data.Pier;

% Create a "Loads" matrix consisting of "P, V2, V3, T, M2, M3"
loads = [data.P, data.V2, data.V3, data.T, data.M2, data.M3];

% Get the "Story" column
story = data.Story;

% Get selected loads
loads_selected = loads((ismember(pier, Piers) & ismember(story, stories)) , :);

end



