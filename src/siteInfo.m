function siteInfo()
%
% Function siteInfo
%   Read and output head informations (name, longitude, latitude, etc) from the seismic data(in sac format)
%
% Input:
%   sac format files in current directiory
%
% Output:
%   *.loc file
%
% Note(1): seismic data is in 'seis_z'; and head information is in 'head', users could obtain more informations from these two arrays with easy coding
%
% Note(2): for data reading, choose appropriate format
%   'b' for data saved on Solaris read into Linux
%   'l' for data saved on Linux read into Lin
%---------------------------------------------------------------------------
%
% Read seismic data (in the current directory, *.bhz, *.bhe, *.bhn)
% Here only *.bhz format as example
Mydir=pwd;
AllF=dir([Mydir,'/*.bhz']);

num_F=length(AllF); 
if num_F<1
    error(sprintf('No .bhz file'));
end

Loc=zeros(3,num_F); % x, y and z

fprintf('Reading *.bhz files\n');
filename=AllF(num_F).name;
[seis_z(:,num_F),head(:,num_F)]=readsac(filename,0,'b'); % Initial, for processing speed 
for i=1:num_F
    filename=AllF(i).name;
    [seis_z(:,i),head(:,i)]=readsac(filename,0,'b');
    %[seis_z(:,i),head(:,i)]=readsac(filename,0,'l');
end

% Data normalization
%maxA=max(abs(seis_z(norm0:norm1,1:num_F))); % local norm
%  for i=1:num_F
%    seis_z(:,i)=seis_z(:,i)/maxA(i);
%end

%=====-----------------------------------------------------------------=====
% Location information
for i=1:num_F
    Loc(1,i)=head(i).STLO; % x value, usually is longtitude
    Loc(2,i)=head(i).STLA; % y value, usually is latitude
    Loc(3,i)=head(i).STEL; % z value, usually is elevation
end

%Output the location data
outF=[num2str(num_F) 'files.loc'];
outID=fopen(outF,'w');
for i=1:num_F
    fprintf(outID,'%s %8.4f %8.4f %7.2f\n',deblank(tempF(i).name),Loc(1,i),Loc(2,i),Loc(3,i));
end
fclose(outID);
