function p = project_paths()
% PROJECT_PATHS  Absolute paths of the project's data, results and figures folders.
%   p = project_paths() returns a struct with fields
%     p.data    - input data (HMD deaths and exposures)
%     p.results - saved simulation outputs (*.mat), read by the *_plots scripts
%     p.figures - where the *_plots scripts export figures
root = fileparts(fileparts(mfilename('fullpath')));
p.data    = fullfile(root, 'data');
p.results = fullfile(root, 'results');
p.figures = fullfile(root, 'figures');
end
