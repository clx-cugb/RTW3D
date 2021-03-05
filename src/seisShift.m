function seisShift(fname,shift)
%
% Function seisShift: shifting of seismic record 
%
% Input: 
%   fname, filename of seismic data, in sac format
%   shift, shifting value, unit: s, shift>0 (right, -->)
%                                   shift<0 (left, <--)
% Output:
%   seismic data file named as ['shf' fname], in sac format
%
% Note:
%   Data after shifting is in seis_1
%
%=====------------------------------------------------------------=====
% Read input data (fname)
[seis(:),head(:)]=readsac(fname,0,'b');
srate=head.DELTA;

% Shift
nt=length(seis(:));
f=abs(round(shift/srate));

seis_1=zeros(nt);
if shift>0
    seis_1(f:nt)=seis(1:nt+1-f);
elseif shift<0
    seis_1(1:nt+1-f)=seis(f:nt);
else
    seis_1=seis;
end

% Output
outfname=['shf' fname];
writesac(seis_1(:),head(:),outfname);
