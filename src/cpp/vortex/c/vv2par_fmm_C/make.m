function make(pustr,varargin)
%MAKE Makefile for FMM2D (CPU/GPU-versions).
%   MAKE('CPU' | 'GPU', ARG1,...) compiles FMM2D for CPU/GPU with the
%   compiler arguments ARG1, ... using the MEX-command. With no
%   arguments, default platform-dependent define arguments are used
%   (passed to MEX through the -D switch).
%
%   Examples:
%     make cpu % make for CPU using all default defines
%     make gpu -DCUDASUPPORT % "minimal" make for GPU
%
%   See also MEX.

% S. Engblom 2011-07-07 (CPU/GPU-syntax added)
% S. Engblom 2011-04-01

% determine platform
if nargin == 0
  fprintf(1,'Defaulting to GPU-make.\n');
  cpu = false;
elseif strcmpi(pustr,'CPU')
  cpu = true;
elseif strcmpi(pustr,'GPU')
  cpu = false;
else
  error('First argument must be a valid processing unit (CPU/GPU).');
end

% Caveats:
%
% -making for GPU/CUDA requires NVIDIAs nvcc to be properly installed
% and available on the path,
%
% -the extension '.cu' should be handled by mex as C++-code (for this you
% need to manually edit the 'mex' bash file, typically residing in
% [matlabroot '/bin']).
%
% Compiler-defines:
%
% Debug:
% ------
% -DEBUGPANELSORT     debug information on panels
%
% Platforms:
% ----------
% -DmwIndex=int       for platforms where mex is not automatically
% -DmwSize=int        linked against the library defining the
%                     mwIndex and mwSize types
%
% FMM-algorithm:
% --------------
% -DINLINECOMPLEX         use inline complex arithmetic (for CPU-performance)
% -DINLINECOMPLEX_DIRECT  use inline complex arithmetic also in direct
%                         evaluations (for CPU-performance)
% -DTHETACUTOFFCHECK      respect cutoff in theta criterion on CPU (always
%                         checked on the GPU)
% -DFMM_THETA=<value>     (default 0.5) sets the constant of
%                         well-separatedness in the theta-criterion
% -DFMM_MESHRATIO=<value> (default 2.0) tunes the criterion for
%                         splitting long edges in the mesh
% -DCHANNELPOT            support for channel potential
%
% Panels:
% -------
% -DPANELSORT                 uses panel positions in the sort algorithm,
%                             necessary when there are only panels and no
%                             point potentials
% -DPANELDUMMYFACTOR=<value>  weight factor (default 2) for panels versus 
%                             potentials
% -DPANELSHRINKBOX            decreases box size after each split
%
% Experimental:
% -------------
% -DRADIALSHRINK              use radial distance to the vortex furthest 
%                             away in a box instead of the diagonal 
%                             (requires PANELSHRINKBOX or CUDASORT)
% -DRADIALSHRINKLIMIT=<value> radial shrink will only be performed for 
%                             boxes containing less than this amount of
%                             points
% -DCUDATIMESORT              print timings for the sorting phase on the
%                             GPU (requires CUDASUPPORT and REORDERCUDA)
% -DASYMMETRIC_CPU            store connection matrices in an asymmetric 
%                             fashion
% -DUNROLLHORNER              use manual unrolling of horner shifts
% -DUNROLLLEVEL=<value>       how long the unroll should be, can be between
%                             1-7 for normal UNROLL or between 1-6 for old
%                             version
% -DSSE2INTRINSIC             use SSE2 instructions for horner shifts
% -DSSE2M2PSSCALING           use SSE2 in m2ps scaling (small time
%                             difference for a64)
% -DSSE2DIRECT                use SSE2 in the direct evaluation
% -DSSEEVALUNROLL             unroll in direct multipole evaluation
% -DSSE2INIT                  use SSE2 in initialization
% -DCROSSORDER                alternates the additions in shift, only
%                             applies if DSSE2INTRINSIC is false
%                             wksp[(j)*2]   += wksp[(j-1)*2];
%                             wksp[(j)*2+1] += wksp[(j-1)*2+1];
% -DFULLCROSSORDER            alternates additions in m2ps shift according
%                             to
%                             wksp1[(j)*2] += wksp1[(j-1)*2];
%                             wksp2[(j)*2] += wksp2[(j-1)*2];
% -DUNROLLHORNEROLD           older more primitive version of unrolled m2ps
%                             shifts. Less efficient, but works better for
%                             old compilers which does not align stack to
%                             16 bytes in 32 bit mode
% -DM2PSINBOX                 enables m2ps interaction between boxes with
%                             same parents (CPU only, always enabled on
%                             cuda)
% -DCHECKNANINPUT             checks for NaN in input vectors and aborts,
%                             CPU only, always true on GPU
% -DOSEENSSEALTERNATIVE       alternative implementation of Oseen with SSE2
%                             as gcc 4.4 compiles incorrect code with
%                             original implementation
%
% CUDA:
% -----
% -DCUDATIME           use timing at the GPU (this is output 1--10 in the
%                      'out'='time'-syntax, see also the 'printtime'
%                      argument)
% -DCUDASUPPORT        compile for CUDA (direct evaluation ported to
%                      the GPU)
% -DMULTIPOLEEVALCUDA  multipole evaluation on the GPU
%                      (mpexp_eval and mpexp_evalmp, requires CUDASUPPORT)
% -DMULTIPOLEINITCUDA  initialization and upward pass on the GPU
%                      (requires CUDASUPPORT)
% -DMULTIPOLESHIFTCUDA downward pass on the GPU (shiftm2ps and shiftp2p,
%                      requires CUDASUPPORT and MULTIPOLEINITCUDA)
% -DCUDASORT           sorting code on the GPU (requires CUDASUPPORT and 
%                      REORDERCUDA)
% -DSORTLIMIT          activates support for sort limits
% -DSORTLIMITVALUES
%    =<string>         default values, see fmm2d.m

mx = mexext;
ver = version;

% different categories of defines
debug_defs = {'-DEBUGPANELSORT' ...       % (1)
             };
platform_defs = {'-DmwIndex=int' ...      % (1)
                 '-DmwSize=int' ...       % (2)
                };
fmm_defs = {'-DINLINECOMPLEX' ...         % (1)
            '-DINLINECOMPLEX_DIRECT' ...  % (2)
            '-DTHETACUTOFFCHECK' ...      % (3)
            '-DFMM_THETA=0.5' ...         % (4)
            '-DFMM_MESHRATIO=2.0' ...     % (5)
            '-DCHANNELPOT' ...            % (6)
           };
panel_defs = {'-DPANELSORT' ...           % (1)
              '-DPANELDUMMYFACTOR=2' ...  % (2)
              '-DPANELSHRINKBOX' ...      % (3)
             };
exp_defs = {'-DRADIALSHRINK' ...          % (1)
            '-DRADIALSHRINKLIMIT=512' ... % (2)
            '-DCUDATIMESORT' ...          % (3)
            '-DASYMMETRIC_CPU' ...        % (4)
            '-DUNROLLHORNER' ...          % (5)
            '-DUNROLLLEVEL=6' ...         % (6)
            '-DSSE2INTRINSIC' ...         % (7)
            '-DSSE2M2PSSCALING' ...       % (8)
            '-DSSE2DIRECT' ...            % (9)
            '-DSSEEVALUNROLL' ...         % (10)
            '-DSSE2INIT' ...              % (11)
            '-DCROSSORDER' ...            % (12)
            '-DFULLCROSSORDER' ...        % (13)
            '-DUNROLLHORNEROLD' ...       % (14)
            '-DM2PSINBOX' ...             % (15)
            '-DCHECKNANINPUT' ...         % (16)
            '-DOSEENSSEALTERNATIVE'};     % (17)
if ~cpu
  cuda_defs = {'-DCUDATIME' ...           % (1)
               '-DCUDASUPPORT' ...        % (2)
               '-DMULTIPOLEEVALCUDA' ...  % (3)
               '-DMULTIPOLEINITCUDA' ...  % (4)
               '-DMULTIPOLESHIFTCUDA' ... % (5)
               '-DCUDASORT' ...           % (6)
               '-DSORTLIMIT' ...          % (7)
               '-DSORTLIMITVALUES={1,0.5,-0.5,1,-1,0.5}' ... % (8)
              };
else
  cuda_defs = {}; % not meaningful without GPU
end

% path, include and link path, source
s = pwd;
include = ['-I' s '/source'];
link = ['-L' s '/source'];
source = {[s '/source/fmm2d.cpp'] ...
          [s '/source/fmm.cpp'] ...
          [s '/source/fmmsort.cpp'] ...
          [s '/source/fmmshift.cpp'] ...
          [s '/source/hornershift.cpp'] ...
          [s '/source/direct.cpp'] ...
          [s '/source/panel.cpp'] ...
          [s '/source/mexpanel.cpp'] ...
          [s '/source/channelpot.cpp'] ...
          [s '/source/directchannelpot.cpp'] ...
          [s '/source/channelpotpanel.cpp'] ...
          [s '/source/cudasort.cu'] ...
          [s '/source/cudaeval.cu'] ...
          [s '/source/cudatiming.cu']};
if cpu
  % not meaningful without GPU
  source = source(cellfun('isempty',strfind(source,'.cu')));
end

if strcmp(mx,'mexa64')
  if ~strncmp(version,'7.11',4) && ~strncmp(version,'7.13',4)
    warning(['Extension .' mexext ' tested with Matlab version(s) ' ...
	     '7.11 and 7.13 only.']);
  end
  cudaroot = '/usr/local/cuda'; % path to cuda
  if cpu
    cxxflags = ['CXXFLAGS=-O3 -msse3 -fPIC -fno-omit-frame-pointer -pthread'];
    cxxlibs = ['CXXLIBS=-Wl,-rpath-link,' matlabroot '/bin/glnxa64 ' ...
               '-L' matlabroot '/bin/glnxa64 ' ...
               '-lmx -lmex -lmat -lm -lstdc++ -lmwblas'];
    cxx='CXX=gcc';
  else
    cxxflags = ['CXXFLAGS=-arch=sm_13 ' ...
                '-Xcompiler=''-O3 -msse3 -fPIC -fno-omit-frame-pointer -pthread'''];
    cxxlibs = ['CXXLIBS=-Wl,-rpath-link,' matlabroot '/bin/glnxa64 ' ...
               '-L' matlabroot '/bin/glnxa64 -L' cudaroot '/lib64 ' ...
               '-lmx -lmex -lmat -lm -lstdc++ -lmwblas -lcudart'];
    cxx='CXX=nvcc';
  end
  ldflags = ['LDFLAGS=-pthread -shared -m64 -Wl,--version-script,' ...
             matlabroot '/extern/lib/glnxa64/mexFunction.map ' ...
             '-Wl,--no-undefined'];
  include = ['-I' cudaroot '/include,' include]; % must include cuda first
  link = [link ' -L' matlabroot '/sys/os/glnxa64' ' -L' cudaroot '/lib64'];

  if ~cpu
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2]) cuda_defs([1:6])];
  else
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2 5:7 8 9:11 13 17])];
  end
  if nargin > 1, define = varargin; end

  fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
  mex(cxx,cxxflags,cxxlibs, ...
      'LD=gcc',ldflags, ...
      include,link, ...
      '-largeArrayDims', ...
      define{:}, ...
      '-outdir',s, ...
      source{:});
elseif strcmp(mx,'mexglx')
  if ~strncmp(version,'7.11',4)
    warning(['Extension .' mexext ' tested with Matlab version(s) 7.11 only.']);
  end
  cudaroot = '/usr/local/cuda'; % path to cuda
  cxxflags = ['CXXFLAGS=-O3 -msse3 -mstackrealign'];
  cxxlibs = ['CXXLIBS=-Wl,-rpath-link,' matlabroot '/bin/glnxa64 ' ...
             '-L' matlabroot '/bin/glnxa64 -L' cudaroot '/lib64 ' ...
             '-lmx -lmex -lmat -lm -lstdc++ -lmwblas -lcudart'];
  ldflags = ['LDFLAGS=-pthread -shared -m64 -Wl,--version-script,' ...
             matlabroot '/extern/lib/glnxa64/mexFunction.map ' ...
             '-Wl,--no-undefined'];
  include = ['-I' cudaroot '/include,' include]; % must include cuda first
  link = [link ' -L' matlabroot '/sys/os/glnxa64' ' -L' cudaroot '/lib64'];

  if ~cpu
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2]) cuda_defs(1:6)];
  else
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2 5:7 8 9:11])];
  end
  if nargin > 1, define = varargin; end

  fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
  mex( '-largeArrayDims',cxxflags, ...
      define{:}, ...
      '-outdir',s, ...
      source{:});
elseif strcmp(mx,'mexmaci')
  if ~strncmp(version,'7.8',3)
    warning(['Extension .' mexext ...
             ' tested with Matlab version(s) 7.8 only.']);
  end
  if ~cpu
    error('No GPU-make currently available for this platform.');
  end
  link = [link ' -lmwblas'];

  define = [debug_defs([]) platform_defs([]) fmm_defs ...
            panel_defs exp_defs([])];
  if nargin > 1, define = varargin; end
  
  fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
  mex(define{:},include,link, ...
      '-outdir',s, ...
      source{:});
elseif strcmp(mx,'mexmaci64')
  if ~strncmp(version,'7.10',4) && ~strncmp(version,'7.11',4)
    warning(['Extension .' mexext ...
             ' tested with Matlab version(s) 7.10 and 7.11 only.']);
  end
  if ~cpu
    error('No GPU-make currently available for this platform.');
  end
  link = [link ' -lmwblas'];

  define = [debug_defs([]) platform_defs([]) fmm_defs ...
            panel_defs exp_defs([])];
  if nargin > 1, define = varargin; end

  fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
  mex(define{:},include,link, ...
      '-largeArrayDims', ...
      '-outdir',s, ...
      source{:});
elseif strcmp(mx,'mexw32')
  if ~strncmp(version,'7.8',4)
    warning(['Extension .' mexext ...
             ' tested with Matlab version(s) 7.8']);
  end
  
  if ~cpu
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2]) cuda_defs(1:6)];
  else
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2 5:7 8 9:11])];
  end
  if nargin > 1, define = varargin; end

  fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
  mex('COMPFLAGS=-c -DMATLAB_MEX_FILE',define{:},include,link, ...
      '-outdir',s, ...
      source{:});
elseif strcmp(mx,'mexw64')
  if ~strncmp(version,'7.11',4)
    warning(['Extension .' mexext ' tested with Matlab version(s) 7.11 only.']);
  end

  if ~cpu
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2]) cuda_defs(1:6)];
  else
    define = [debug_defs([]) platform_defs([]) fmm_defs ...
              panel_defs exp_defs([1:2 5:7 8 9:11 13 16])];
  end
  if nargin > 1, define = varargin; end

  fprintf(1,['MAKE ' sprintf('%s ',define{:}) '\n']);
  mex(define{:}, ...
      '-largeArrayDims', ...
      '-outdir',s, ...
      source{:});
else
  error('No make currently available for this platform.');
end
