function v = get_option(opts, name, default)
% GET_OPTION  Field NAME of struct OPTS if present, otherwise DEFAULT.
if isstruct(opts) && isfield(opts, name)
    v = opts.(name);
else
    v = default;
end
end
