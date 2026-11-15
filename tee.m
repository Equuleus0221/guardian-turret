function tee(fid, varargin)
%TEE  fprintf to the command window and to the file fid.
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
