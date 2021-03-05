%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     %
% RTW3D (Regularization of 3D teleseismic wavefield)  %
%                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*****=====================================================================================*****
% 
% This code is to regularize the randomly distributed seismic stations onto
% a uniform grid. This can greatly reduce the artifacts caused by uneven
% stations intervals, and can avoid the step effect due to insufficient
% station numbers. Both the theory and application examples can be found in our papers.
% 
% This code was firstly finished by Jinhai Zhang, Aug, 2013, and modified by Chengliang Xie, June 2020.
% This code can be freely used for personal and non-commercial purpose. 
%
% Related papers:
%   Zhang, J., T. Zheng, 2015. Receiver Function Imaging with Reconstructed Wavefields from Sparsely Scattered Stations. Seismological Research Letters 86, 165-172.
%   Xie,C., Y. Fang, and J. Zhang, 2021, Regularizing 3D teleseismic wavefield for receiver function imaging using radial basis function. Journal of Geophysical Research: Solid Earth (submitted)
%
% Users can send an email to authors for any comment or error.
%
% -------------------------------------
% -------------------------------------
%  Dr. Chengliang Xie
%  China University of Geosciences (Beijing)
%  No. 29, Xueyuan Road, Haidian District
%  Beijing, 100083
%  CHINA
%  Phone: 86-15201018722
%  Email: clx@cugb.edu.cn
%
% -------------------------------------
%  Dr. JinHai Zhang
%  Institute of Geology and Geophysics,
%  Chinese Academy of Sciences.
%  P. O. Box 9825
%  Beitucheng Western Road, 19#
%  Beijing, 100029
%  CHINA
%  Phone: 86-13552750211
%  Email: zjh@mail.iggcas.ac.cn, geophysics.zhang@gmail.com  
%  https://www.researchgate.net/profile/Jin-Hai-Zhang
%
% -------------------------------------
%
% V1.0  Aug 2013
%   2D-CUBIC, 2D interpolation using cubic method in MATLAB (Zhang and Zheng, 2015)
%
% V2.0  June 2020
%   Rdial basis function (RBF) was supplied as an option for interpolation methods.
%   Several interpolation options are available as follows
%   1D-RBF,   1D interpolation using RBF for profile data
%   2D-CUBIC, 2D interpolation using cubic method in MATLAB for array data
%   2D-RBF,   2D interpolation using RBF for array data
%   3D-RBF,   3D interpolation using RBF for array data
% 
% Acknowledgments:  
%   RBF function is from Matlab Central by Alex Chirokov (alex.chirokov@gmail.com). 
% 
%   We also thank  C.D. Saragiotis (csar@auth.gr), Y.M. Altman (altmany@gmail.com) and F.J. Simons (fjsimons@alum.mit.edu) 
%   for their contributions on codes associated with the sac format files. 
% 
%*****=====================================================================================*****

clear all;
clc;
% Display general informations of this code
InfoDisp();

% Read configure information
inF='../sp/test.cfg'  % System parameter file
[ndlongi,ndlati,win,dir,lstFN,outDir,dimType,rbf]=readPar(inF); % Read system parameters

%=====-------------------------------------------------------------------------------------=====
minLong=ndlongi(1,1); % Minimum longitude for interpolation area
maxLong=ndlongi(1,2); % Maximum longitude for interpolation area
intLong=ndlongi(1,3); % Interval of interpolation in longitude direction

minLat=ndlati(1,1); % Minimum latitude for interpolation area
maxLat=ndlati(1,2); % Maximum latitude for interpolation area
intLat=ndlati(1,3); % Interval of interpolation in latitude direction

norm0 =win(1,1); % Minimum index of normalization window
norm1 =win(1,2); % Maximum index of normalization window

%=====-------------------------------------------------------------------------------------=====
% Read records and conduct interpolation
[nLst,temp]=size(lstFN); % nLst: number of list files
for i=1:nLst
    % For each *.lst file
    lstName=char(lstFN{i});  % filename of each list file
    lstID=fopen(lstName,'r'); % Open the list file
    numF=0; % number of data files
    while(~feof(lstID))
        % Read all data files'name (e.g *.BHR/BHZ/BHT) in current list file
        strtemp=fgetl(lstID);
        numF=numF+1;
        fname{numF}=[dir strtemp]; % Get the data filename with directory
    end
    fclose(lstID);

    %------------------------------------------------------------------------------------------
    % Read Data
    % The following two lines are for initialization of seis and head, that would improve efficience
    filename=char(fname{numF});
    % Note: the option 'b' or 'l' depends on the type of sac file
    % [seis(:,numF),head(:,numF)]=readsac(filename,0,'b');
    [seis(:,numF),head(:,numF)]=readsac(filename,0,'l'); % Read sac data

    % Check whether maximum normalization window(norm1) is less than total length of data
    if norm1>size(seis,1)
        error(sprintf('Ending index should be smaller than the maximum number of samples\n'));
    end

    % Determine the data type (*.bhe, *.bhn or *.bhz)
    inx=strfind(filename,'.');
    inx2=inx(length(inx));
    addr=filename(inx2+1:length(filename)); % addr: file address
    fprintf('Running on %s files...\n',addr);

    % Main part for data reading
    for j=1:numF-1
        filename=char(fname{j}); % filename of sac file
        %[seis(:,j),head(:,j)]=readsac(filename,0,'b'); % Read data
        [seis(:,j),head(:,j)]=readsac(filename,0,'l');
    end

    %------------------------------------------------------------------------------------------
    % DATA normalization
    maxA=max(abs(seis(norm0:norm1,1:numF))); % local norm
    %for i=1:numF
    %    fprintf('maxA= %f file=%s\n',maxA(i),char(fname{i}));
    %end
    for i=1:numF
        seis(:,i)=seis(:,i)/maxA(i); % Data normalization
    end
    
    % Location of seismic stations
    for i=1:numF
        y(i)=head(i).STLO; % Longitude
        x(i)=head(i).STLA; % Latitude
    end

    %------------------------------------------------------------------------------------------
    % Grid generating (Grids)
    lati = minLat :intLat :maxLat ; % Latitude of grid
    long = minLong:intLong:maxLong; % Longitude of grid
    [XI,YI] = meshgrid(lati,long); % Generate the grid

    %------------------------------------------------------------------------------------------
    %*** Interpolation ***
    t0=1; % Start index for time slices
    t1=size(seis,1); % End index for time slices

    % Data reconstruction according to the user definition
    if ~isempty(strfind(dimType,'RBF'))
        % Using RBF method

        % Calculating the average distance of the stations as the shape parameter
        [sP,nearInx]=AverDis(x,y); % SP, shape parameter
        smooth=0; % Smoothing parameter

        % *** 1D ***, orinted in x direction
        if ~isempty(strfind(dimType,'1D'))
            XI=XI'; % Transpose
            YI=YI';
            z=zeros(1,numF);
            Loc=zeros(1,numF);
            Loc(1,:)=x(:); % Get the location in x direction
            dat=zeros(t1,length(lati),length(long)); % Initialize dataset

            for it=t0:t1
                z=seis(it,:); % seismic data
                ZI=SeisIntep1D_RBF(Loc,z,lati,rbf,sP,smooth); % Interpolation
                ZI=ZI';
                dat(it,:,:)=ZI; % Get the result into 'dat'
                fprintf(' %dth layer\n',it);
            end
        end

        % *** 2D ***
        if ~isempty(strfind(dimType,'2D'))
            z=zeros(1,numF);
            Loc=zeros(2,numF);
            Loc(1,:)=x(:); % Location in x direction (Latitude)
            Loc(2,:)=y(:); % Location in y direction (Longitude)
            dat=zeros(t1,length(long),length(lati));

            for it=t0:t1
                z=seis(it,:); % Seismic data
                [temp,ZI]=SeisIntep2D_RBF(Loc,z,lati,long,rbf,sP,smooth); % Interpolation
                dat(it,:,:)=ZI; % Results
                fprintf(' %dth layer\n',it);
            end
        end

        % *** 3D ***
        if ~isempty(strfind(dimType,'3D'))
            winLen=t1;
            % winLen, window length, i.e. number of time slices in each window. This could be user defined
            % Note: 
            % (1) Generally, datset has to be divided into several windows (in t direction) because of the memory limitation
            % (2) If set winLen=t1, means no division for datset. Acturally this could be easily computed
            %     by the following commented code ( User could search 'No division')
            % Tip: set winLen to let mod(nt,winLen)=0, might reduce troubles...

            nt=t1;  % Total number of slices
            nWin=floor(nt/winLen);  % Number of windows
            Last=mod(nt,winLen); % Number of slices in the last window

            % --------------------------------------------------------------------------------
            % *** Construct 3D data set ***
            % --> Location in Loc
            N=length(x); % Number of stations
            inods=1; % index for each node 
            timeDet==head(1).DELTA; % Sampling interval, here we suppose all stations have same 'timeDet'
            for it=1:nt % For the ith time slice
                for jst=1:N
                    % For the jth station
                    Loc(1,inods)=timeDet*it; % Location in t direction (time) 
                    Loc(2,inods)=x(jst); % Location in x direction (Latitude)
                    Loc(3,inods)=y(jst); % Location in x direction (Longitude)
                    inods=inods+1;
                end
            end

            % --> Value in Node
            inods=1;  % index for node 
            for it=1:nt
                for jst=1:N
                    Node(inods)=seis(it,jst); % Seismic data in each interpolated node
                    inods=inods+1;
                end
            end
            fprintf('Data preparation is done\n')

            % --------------------------------------------------------------------------------
            % *****INTERPOLATION***** 
            % Grid and resconstructed data set
            tT=1*timeDet:timeDet:nt*timeDet;
            lenT=length(tT);
            lenX=length(long); % Number of interpolation data (in x/latitude direction)
            lenY=length(lati); % Number of interpolation data (in y/longitude direction)
            dat=zeros(lenT,lenX,lenY); %Data after interpolation, in 3D matrix format
            ngrd=lenT*lenX*lenY; %Number of grids
            VI=zeros(1,ngrd); %Data after iterpolation, in 1D vector format

            % --------------------------------------------------------------------------------
            fprintf(' 3D interpolation is running, might take very long time...\n');
            % ***** Interpolation in each window

            % Initialization for the following arrays
            loctemp=zeros(3,winLen*N);
            nodetemp=zeros(1,winLen*N);
            vStemp=zeros(winLen,lenX,lenY);
            vItemp=zeros(1,winLen*lenX*lenY);

            sxSlice=1; % Start index of slice
            edSlice=winLen; % End index of slice

            sIx=1; % Start index of node
            eIx=winLen*N; % End index of node

            sxGrid=1; % Start index of grid
            edGrid=winLen*lenX*lenY; % End index of grid

            for iw=1:nWin
                loctemp(1:3,:)=Loc(1:3,sIx:eIx); % Location data for interpolation
                nodetemp(:)=Node(sIx:eIx); % Node value for interpolation

                inteTtemp=sxSlice:1:edSlice; % number of slices for each interpolation
    
                [vItemp,vStemp]=SeisIntep3D_RBF(loctemp,nodetemp,inteTtemp,long,lati,rbf,sP,smooth); % Interpolation

                % Save results
                VI(sxGrid:edGrid)=vItemp;
                %vConst(sxSlice:edSlice,:,:)=vStemp;
                fprintf('%dth window (length=%d) is done\n',iw,winLen);

                % Update index
                sxSlice=edSlice+1;  % For slice index
                edSlice=edSlice+winLen;

                sIx=eIx+1; % Node index
                eIx=eIx+winLen*N;

                sxGrid=edGrid+1; % Grid index
                edGrid=edGrid+winLen*lenX*lenY;
            end
            % ------
            %The last window that length is less than winLen
            if Last
                edSlice=nt;
                eIx=nt*N;
                edGrid=nt*lenX*lenY;

                loctemp=zeros(3,Last*N);  % Local node
                nodetemp=zeros(1,Last*N);
                vStemp=zeros(Last,lenX,lenY); % Local grid
                vItemp=zeros(1,Last*lenX*lenY);

                loctemp(1:3,:)=Loc(1:3,sIx:eIx);  % For node
                nodetemp(:)=Node(sIx:eIx);
                inteTtemp=sxSlice:1:edSlice;
    
                [vItemp,vStemp]=SeisIntep3D_RBF(loctemp,nodetemp,inteTtemp,long,lati,rbf,sP,smooth);

                % Save results
                VI(sxGrid:edGrid)=vItemp;
                %vConst(sxSlice:edSlice,:,:)=vStemp;

                fprintf('Last window (length=%d) is done\n',Last);
            end

            % --------------------------------------------------------------------------------
            % Obtain the final Matrix
            % Array arrangement as the output necessary
            igrid=1;
            istart=1;
            for iw=1:nWin
                iend=iw*winLen;
                for i=1:lenX
                    for j=istart:iend
                        for k=1:lenY
                            dat(j,i,k)=VI(igrid);
                            igrid=igrid+1;
                        end
                    end
                end
                istart=iend+1;
            end
            % -----
            if Last
                iend=nt;
                for i=1:lenX
                    for j=istart:iend
                        for k=1:lenY
                            dat(j,i,k)=VI(igrid);
                            igrid=igrid+1;
                        end
                    end
                end
            end

            % --------------------------------------------------------------------------------
            % No division ( if memory is large enough for data set) 
            % [VI,temp]=SeisIntep3D_RBF(Loc,Node,inteT,inteX,inteY,rbf,const,smooth);
            % dat(:,:,:)=temp;

            % --------------------------------------------------------------------------------
            fprintf('Successful!\n');
            ZI=ones(length(long),length(lati));
        end
    elseif ~isempty(strfind(dimType,'CUBIC'))
        % Using Matlab-cubic method
        % The following three lines: Initial array for saving computational cost
        z=seis(1,:);
        ZI=griddata(x,y,z,XI,YI,'cubic'); % Generate grid
        dat(t1,:,:)=ZI(:,:);  % data set

        for it = t0:t1
            z=seis(it,:); % data for each time slice
            ZI=griddata(x,y,z,XI,YI,'cubic'); % Interpolation
            dat(it,:,:)=ZI(:,:); % Save the results
        end
    end

    %------------------------------------------------------------------------------------------
    % *** Head information for griddata***
    % Note: some information might be missed in present version, thank you for your reminding and improvements
    hd(1:size(ZI,1),1:size(ZI,2))=head(1);

    % AZ 
    for i1=1:numF
        zz(i1)=head(i1).AZ;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4'); 
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).AZ=ZZ(ix,iy);
        end
    end

    % BAZ
    for i1=1:numF
        zz(i1)=head(i1).BAZ;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).BAZ=ZZ(ix,iy);
        end
    end

    % DIST
    for i1=1:numF
        zz(i1)=head(i1).DIST;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).DIST=ZZ(ix,iy);
        end
    end

    % EVLA
    for i1=1:numF
        zz(i1)=head(i1).EVLA;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).EVLA=ZZ(ix,iy);
        end
    end

    %EVLO
    for i1=1:numF
        zz(i1)=head(i1).EVLO;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).EVLO=ZZ(ix,iy);
        end
    end

    % GCARC
    for i1=1:numF
        zz(i1)=head(i1).GCARC;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).GCARC=ZZ(ix,iy);
        end
    end

    % STLA
    for i1=1:numF
        zz(i1)=head(i1).STLA;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).STLA=ZZ(ix,iy);
        end
    end

    % STLO
    for i1=1:numF
        zz(i1)=head(i1).STLO;
    end
    ZZ=griddata(x,y,zz,XI,YI,'v4');
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            hd(ix,iy).STLO=ZZ(ix,iy);
        end
    end

    %------------------------------------------------------------------------------------------
    % Output result files
    % Checking wheter there is folder './Res/' (Note: for windows, '\' might be necessary)
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    % Write the results into output 'sac files'
    for ix=1:size(ZI,1)
        for iy=1:size(ZI,2)
            if ~isnan(hd(ix,iy).AZ) && ~isnan(dat(1,ix,iy))
                outfname=[num2str(ix,'%04d'),'_',num2str(iy,'%04d'),'.',addr]; % Output filename
                writesac(dat(:,ix,iy),hd(ix,iy),[outDir outfname]); % Write into the file
            end
        end
    end

    %------------------------------------------------------------------------------------------
end

%=====-------------------------------------------------------------------------------------=====
disp('Done.')
disp(' ')
%input('Press any key to quit.','s'); 
% exit % For running background
