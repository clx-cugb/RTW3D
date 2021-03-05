function [longi,lati,win,dir,lstFN,outDir,dimType,rbf]=readPar(inF)
% Function readPar: read system parameters 

% Input: 
%   inF, Parameters file name
%
% Output
%   longi[1,3], mininum and maximum value of longitude in longi[1,1] and longi[1,2], and interval in logi[1,3]
%   lati[1,3], min and max value of latitude, and also interval
%   win[1,2] min and max value of normalization window
%   dir, directory of data files
%   lstFN{}, list files' name
%   outDir, directory of output files
%   dim, dimension of interpolation method
%   rbf, radial basis function, e.g. multiquare, gaussian, ...

if ~exist(inF,'file')
    error(sprintf('Missing parameter file: %s\n',inF));
end

inID=fopen(inF,'r');
while (~feof(inID))
    strtemp=fgetl(inID);
    if(strtemp(1)=='#')
        continue;  % The comments
    end
    if(strtemp(1)=='>')
        if ~isempty(strfind(strtemp,'NORM'))
            longi=fscanf(inID,'%f\n',[1,3]);
            lati=fscanf(inID,'%f\n',[1,3]);
            win=fscanf(inID,'%d\n',[1,2]);
            checkNorm(win(1,1),win(1,2));
            continue;
        end
        if ~isempty(strfind(strtemp,'DATAF'))
            dir=fgetl(inID);
            nF=fscanf(inID,'%d\n');
            lstFN=cell(nF,1);
            for i=1:nF
                lstFN{i}=fgetl(inID);
            end
            outDir=fgetl(inID);
            continue;
        end
        if ~isempty(strfind(strtemp,'DIM-TYPE'))
           dimType=fgetl(inID); 
           if ~isempty(strfind(dimType,'RBF'))
               rbf=fgetl(inID);
           else
               rbf='';
           end
           continue;
        end
    end
end
fclose(inID);

%=====-------------------------------------------------------------------------=====
function checkNorm(norm0,norm1)
% Function checkNorm
% Check whether the definition of normalization window  is suitable (Integer)

if norm0~=fix(abs(norm0))
    input('WRONG:begining index should be a positive integer','s'); 
end
if norm1~=fix(abs(norm1))
    input('WRONG:ending index should be a positive integer','s'); 
end
if norm0>=norm1
    input('WRONG:ending index should be larger than begining index','s'); 
end
