function preEval(conf)
% Function preEval: 
%   Pre-evaluation for seismic stations, to determine whether the specific station is suitable for interpolation 
%   Two major steps:
%   --> Low pass wave filtering
%   --> Evaluation based on correlation degree
%
%----------------------------------------------------------------------------------------
%   Input: 
%       conf, filename of configure file (*.cfg)
%   Output:
%       reports of evaluation in *.out file
%
% Note: 
%   (1) For low pass filter, default frequency is 0.1 Hz, user revision could be made by searching 'Wp='
%   (2) The evaluation standard is important, and user definition is necessary by searching 'STANDARD' 
%=====******************************************************************************=====
%
% Read configure file
[ndlongi,ndlati,win,dir,lstFN,outDir,dimType,rbf]=readPar(conf);

% Window index for normalization
norm0 =win(1,1);
norm1 =win(1,2);

%----------------------------------------------------------------------------------------
% Read seismic data 
[nLst,temp]=size(lstFN);
for i=1:nLst
    % For each *.lst file
    lstName=char(lstFN{i});
    lstID=fopen(lstName,'r');
    numF=0; %number of data files
    while(~feof(lstID))
        % Read all data files' name in current list file
        strtemp=fgetl(lstID);
        numF=numF+1;
        fname{numF}=[dir strtemp]; % File name (include directory)
        stname{numF}=strtemp; % Station name
    end
    fclose(lstID);

    %-----------------------------------
    outName=strrep(lstName,'.lst','.out'); % Output file's name
    outID=fopen(outName,'w');

    %-------------------------------------------------------------------------------------
    % The following two lines are for initial of seis and head, that would improve efficience
    filename=char(fname{numF});
    [seis(:,numF),head(:,numF)]=readsac(filename,0,'b');
    %[seis(:,numF),head(:,numF)]=readsac(filename,0,'1');

    % Check whether maximum normalization window(norm1) is less than total length of data
    % 检查最大的正则窗是否大于数据的总长？
    if norm1>size(seis,1)
        error(sprintf('Ending index should be smaller than the maximum number of samples\n'));
    end

    % Main part for data reading
    for j=1:numF-1
        filename=char(fname{j});
        [seis(:,j),head(:,j)]=readsac(filename,0,'b');
        %[seis(:,j),head(:,j)]=readsac(filename,0,'1');
    end
    fprintf('Data reading is done\n');
    fprintf('Data evaluation...\n');

    %-------------------------------------------------------------------------------------
    % Normalization
    maxA=max(abs(seis(norm0:norm1,1:numF)));
    for i=1:numF
        seis(:,i)=seis(:,i)/maxA(i); % Data normalization
        longi(i)=head(i).STLO; % Longitude
        lati(i)=head(i).STLA; % Latitude
        srate(i)=head(i).DELTA; % Sampling interval (unit: s)
    end

    %=====-------------------------------------------------------------------------=====
    % Low pass filter
    fs=1./srate(1); % Sampling rate
    Wp=0.1; % Cut off frequency of the low pass filter
    [b,a]=butter(4,Wp*2/fs,'low'); % Butterworth IIR filter 
    for i=1:numF
        seis_filt(:,i)=filter(b,a,seis(:,i)); % Filtering for each station
    end

    %=====-------------------------------------------------------------------------=====
    % Parameters setting for evaluation
    % Calculate distance between two stations 
    dist=zeros(numF,numF);
    for i=1:numF
        for j=i+1:numF
            dist(i,j)=sqrt((longi(i)-longi(j))^2+(lati(i)-lati(j))^2); % Calculte the distance between ith an jth station
            dist(j,i)=dist(i,j); % Assignment for symmetrical element
        end
    end

    % Set evaluation standard: search radius and correlation threshold etc
    % *** STANDARD *** 
    srd=0.5;  % srd, search radius, evaluate stations belong in srd
    cor=0.6;  % cor, correlation threshold, larger than cor means 'good'
    thN=1./3; % thN, proportion threshold, satisfy this means 'good'

    % Data band (times) [StBand, EdBand] to be evaluated
    StBand=20;  % Start time, unit: sec.
    EdBand=40;  % End Time, unit: sec.
    StN=round(StBand/srate(1)); % Index associated with StBand
    EdN=round(EdBand/srate(1)); % Index associated with EdBand
    % ***

    

    %=====-------------------------------------------------------------------------=====
    % Evaluate each station and give suggestion 
    for i=1:numF
        % Determine station index that satisfied the search radius
        indx=find(dist(i,:)<=srd);
        ns=length(indx); % number of stations in the search range
        ncount=0;
        isgood=0;

        %------------------------------------------------------------------------------
        for j=1:ns
            if i==indx(j)
                continue;  % Station itself
            end
            RR=corr2(seis_filt(:,i),seis_filt(:,indx(j))); % correlation degree
            RR=abs(RR);
            if RR>=cor
                ncount=ncount+1;  % Count stations with statisfied correlation degree that has RR>=cor 
            end
            if ncount>=thN*ns
                %'Good station'
                isgood=1;
                break;
            end
        end

        %------------------------------------------------------------------------------
        % For the possible 'bad' stations
        if ~isgood
            if ns==1
                fprintf(outID,'=== Caution: No station near station %s\n',char(stname{i}));
                %No stations in research range (for ith station), only itself
                continue;
            end
            %--------------------------------------------------------------------------
            % Give suggestions that make 'bad station' to be 'good'
            % --> Firstly, find a 'best' station that has most 'good relationship' 
            ntemp=zeros(ns,1);
            for j=1:ns
                for k=j+1:ns
                    rtemp=corr2(seis_filt(:,indx(j)),seis_filt(:,indx(k))); % Calculate each correlation degree
                    rtemp=abs(rtemp);
                    if rtemp>=cor
                        ntemp(j)=ntemp(j)+1; % Count for the jth station while it has rtemp>=cor
                        ntemp(k)=ntemp(k)+1; % Count for the kth station
                    end
                end
            end

            maxtemp=max(ntemp);
            inxbest=find(ntemp==maxtemp); % index of the best station
            ix=inxbest(1);

            % --> Secondly, determine displacement to make 'ith station better' (reference: inxbest station)
            %Find the maximum aplitude
            [maxA,inA]=max(seis(StN:EdN,indx(ix)));
            [maxB,inB]=max(seis(StN:EdN,i));

            %*****Version 1: shift according to correlation degree
            %sA=inA(1)-shiftN;
            %eA=inA(1)+shiftN;
            %sB=inB(1)-shiftN;
            %eB=inB(1)+shiftN;

            %[COR,lag]=xcorr(seis(sA:eA,indx(ix)),seis(sB:eB,i),'coeff');
            %maxCor=max(COR);
            %shiftValue=lag(find(COR==maxCor));
            %*****Version 1 End

            %*****Version 2: shift directly from the inx of maxA and maxB
            shiftValue=inA(1)-inB(1);
            %*****Version 2 End

            % --> Thirdly, Checking
            nt=length(seis(:,i));
            newData=zeros(nt,1);
            tempData=zeros(nt,1);
            f=abs(shiftValue);
            if shiftValue>0
                %Shift data of ith station to right (-->)
                newData(f:nt)=seis(1:nt+1-f,i);
                tempData=filter(b,a,newData);
                corTemp=abs(corr2(seis_filt(f:nt,indx(ix)),tempData(f:nt))); % Corr. degree after shifting
            elseif shiftValue<0
                %Shift data  to left (<--)
                newData(1:nt+1-f)=seis(f:nt,i);
                tempData=filter(b,a,newData);
                corTemp=abs(corr2(seis_filt(1:nt+1-f,indx(ix)),tempData(1:nt+1-f))); % Corr. degree after shifting
            else
                fprintf(outID,'=== Caution: too few stations (Num=%d) around %s\n',ns-1,char(stname{i}));
                continue;  %The most situation is i==indx(ix);
            end

            % Give shifting suggestion
            fprintf(outID,'*** Suggestion: Shift %s with %f s, COR=%f with station %s\n',char(stname{i}),srate(i)*shiftValue,corTemp,char(stname{indx(ix)}));
        end
    end
fclose(outID);
end

%----------------------------------------------------------------------------------------
fprintf('=====Done=====\n');
