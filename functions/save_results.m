function save_results(file, res)
% SAVE_RESULTS  Save a results struct, with the MATLAB version and date, to a .mat file.
res.matlabVersion = version;
res.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
save(file, '-struct', 'res', '-v7');
fprintf('Saved %s\n', file);
end
