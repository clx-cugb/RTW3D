function [VI,ResIntep]=SeisIntep3D_RBF(nodesLoc,nodesValue,inteT,inteX,inteY,rbf,const,smooth)
%
% Function SeisIntep3D_RBF
%    3D Radial Basis Function interpolation for seismic data 
%
%*****05/23/2020
%=====3D interpolation for wave field
%
% Input:
%    nodesLoc: Location of interpolation nodes, Dim*N matrix where Dim is dimenstion, N is number of nodes
%        NOTE: for 3D interpolation, Dim=3, where in  x, y and t direction
%    nodesValue: Value of interpolation nodes, 1*N vector
%    inteT: Grid in x direction, i.e. Xi vector, 1*Nx vector
%    inteX: Grid in x direction, i.e. Xi vector, 1*Nx vector
%    inteY: Grid in y dierction, i.e. Yi vector, 1*Ny vector
%    rbf: Name of radial basis function, e.g. gaussian, multiquadric,cubic,linear and thinplate
%    const: Constant for RBF of Gaussian and multiquadric, usually is average variance of 'Distance', default value 2
%    smooth: Smoothing parameter, default value is 0. Might be helpful for noisy data
%
% Output:
%    VI: Interpolation results, 1*N vector 
%    ResIntep: Interpolation results, Nx*Ny Matrix
%=============================================================================

% PARAMETERS
% Options
if nargin<6
    fprintf('Error: Missing parameters\n');
    return;
elseif nargin<7
    const=2; % Default value
elseif nargin<8
    smooth=0; % Default value
end
rbf=lower(rbf);
rbf=deblank(rbf);

%Dimension and Number of Nodes
[Dim N]=size(nodesLoc);

%=====-------------------------------------------------------------------=====
%RBF INTERPOLATION
op=rbfcreate(nodesLoc(1:Dim,:),nodesValue,'RBFFunction',rbf,'RBFConstant',const,'RBFSmooth',smooth);
%op=rbfcreate(nodesLoc(1:Dim,:),nodesValue,'RBFFunction',rbf,'RBFSmooth',smooth);  % Use default const

%Checking, i.e. Show difference of nodes (before and after interpolation) --> Validity of RBF
%rbfcheck(op);

[TI,XI,YI]=meshgrid(inteT,inteX,inteY);
VI=rbfinterp([TI(:)';XI(:)';YI(:)'],op);
ResIntep=reshape(VI,size(TI));
