function [data] = fileload(varargin)
%fileload Summary of this function goes here
%   Detailed explanation goes here

% check inputs
if nargin == 0; error('fileload(filename) or fileload(filename,format)'); end
filename = varargin{1};
if numel(varargin) > 1
    format = varargin{2};
    if numel(varargin) > 2
        delimiter = varargin{3};
    end
else
    format = '%f';
end

% check if file exists
if ~exist(filename,'file')
    error('File does not exist.')
end

% open file
fid = fopen(filename);

% read data
if strcmp(format,'bin')
    data = fread(fid);
elseif exist('delimiter','var')==1
    data = textscan(fid,format,'delimiter',delimiter);
    data = data{1};
else
    data = textscan(fid,format);
    data = data{1};
end

% close file
fclose(fid);


end

