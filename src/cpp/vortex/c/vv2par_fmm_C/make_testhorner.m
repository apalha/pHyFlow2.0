function make_testhorner(varargin)

defs = {'-DUNROLLHORNER' ...   % (1)
        '-DUNROLLLEVEL=1' ...  % (2)
        '-DSSE2INTRINSIC' ...  % (3)
        '-DUNROLLHORNEROLD' ...% (4)
        '-DCROSSORDER' ...     % (5)
        '-DFULLCROSSORDER'};   % (6)
% path, include and link path, source
s = pwd;
include = ['-I' s '/source'];
link = ['-L' s '/source'];
source = {[s '/source/testhornershift.cpp'] ...
          [s '/source/hornershift.cpp']};

define=defs([1,2]);
if nargin > 0, define = varargin; end
fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
if strcmp(mexext,'mexw32')
  mex('COMPFLAGS=-c -DMATLAB_MEX_FILE -msse2 -mstackrealign',...
    '-largeArrayDims',define{:},include,link,'-outdir',s,source{:});
else
  mex('-largeArrayDims',define{:},include,link,'-outdir',s,source{:});
end
