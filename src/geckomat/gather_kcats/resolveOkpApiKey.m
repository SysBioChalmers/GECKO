function key = resolveOkpApiKey(argKey, basePath)
% resolveOkpApiKey  Resolve the OpenKineticsPredictor API key.
%
% Internal helper shared by submitOpenKineticsPredictor and
% fetchOpenKineticsPredictor; not intended to be called directly.
%
% Checked in order: the argKey input, the OKP_API_KEY environment variable,
% then a plain-text file data/okpApiKey.txt under basePath (this name is
% git-ignored by GECKO). The key is a secret and must never be read from the
% model adapter, which is shared/committed.
%
% Parameters
% ----------
% argKey : char
%     an explicitly provided key, or empty to fall through to the
%     environment variable / file.
% basePath : char
%     the model adapter's params.path, under which data/okpApiKey.txt is
%     looked up.
%
% Returns
% -------
% key : char
%     the resolved API key.
%
% See also
% --------
% submitOpenKineticsPredictor, fetchOpenKineticsPredictor

if ~isempty(argKey)
    key = strtrim(char(argKey));
    return
end
envKey = getenv('OKP_API_KEY');
if ~isempty(envKey)
    key = strtrim(envKey);
    return
end
keyFile = fullfile(basePath,'data','okpApiKey.txt');
if exist(keyFile,'file')
    key = strtrim(fileread(keyFile));
    return
end
error(['No OpenKineticsPredictor API key found. Provide it as the apiKey ' ...
       'argument, set the OKP_API_KEY environment variable, or place it in ' ...
       '%s. Generate a key at https://predictor.openkinetics.org/.'], keyFile);
end
