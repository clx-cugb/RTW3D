function ResIntep=SeisIntep1D_RBF(nodesLoc,nodesValue,inteX,rbf,const,smooth)
%
% Function SeisIntep1D_RBF
% 
% Input:
%    nodesLoc: Location of interpolation nodes, Dim*N matrix where Dim is dimenstion, N is number of nodes
%        NOTE: for 1D interpolation, Dim=1
%    nodesValue: Value of interpolation nodes, 1*N vector
%    inteX: Grid in x direction, i.e. Xi vector, 1*Nx vector
%    rbf: Name of radial basis function, e.g gaussian, multiquadric,cubic,linear and thinplate
%    const: Constant for RBF of Gaussian and multiquadric, usually is average variance of 'Distance', default value 2
%    smooth: Smoothing parameter, default value is 0. Might be helpful for noisy data
%
% Output:
%    ZI: Interpolation results, 1*N vector 
%    ResIntep: Interpolation results, Nx*Ny Matrix
%=============================================================================

% PARAMETERS
% Options
if nargin<5
    fprintf('Error: Missing parameters\n');
    return;
elseif nargin<6
    const=2; % Default value
elseif nargin<7
    smooth=0; % Default value
end
rbf=lower(rbf);
rbf=deblank(rbf);  % rbf method

% Dimension and Number of Nodes
[Dim N]=size(nodesLoc);

%=====-------------------------------------------------------------------=====
% RBF INTERPOLATION
op=rbfcreate(nodesLoc(1:Dim,:),nodesValue,'RBFFunction',rbf,'RBFConstant',const,'RBFSmooth',smooth);
%op=rbfcreate(nodesLoc(1:Dim,:),nodesValue,'RBFFunction',rbf,'RBFSmooth',smooth);  % Use default const

% Checking, i.e. Show difference between nodes (before and after interpolation) --> Validity of RBF
% rbfcheck(op);

XI=inteX;
ResIntep=rbfinterp(XI(:)',op);
