function varargout = struct2array(structure,variables)
% This function fills a gap in Matlab's built-in set of functions. It
% allocates a set of values from a structure to a set of variables in an
% array of the same size. The inputs are the name of the structure followed
% by a string array containing the fieldnames. The output array must be in
% the same order as the variables array! Warning: if a fieldname is not
% recognised, the output is an empty cell.

if length(variables)~=nargout, error('Number of input and output variables must be the same.'); end

varargout = cell(1,nargout);
for i = 1:length(variables)
    if isfield(structure,variables{i})
        varargout{i} = structure.(variables{i});
    else % if not a valid fieldname, output empty cell
        varargout{i} = [];
    end
end

end
