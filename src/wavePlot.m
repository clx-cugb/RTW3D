function [Rms,Cor]=wavePlot(sName,cflag,cfgF)

% Function wavePlot
% Read seismic data and give simple illustration of wave curves
%
%*****=================================================================*****
%
% 07/20/2020 Modification:
%   (1) Add output parameters: Rms, Cor. The two value is for evaluating quality of interpolation
%   (2) Write the rms and correlation information in the pdf file
%   Note: Rms and Cor only makes sence while length(sName)=2 and cflag~=0
%
%*****=================================================================*****
%
% 07/13/2020
% Input:
%   sName: cell array, input seismic files
%   cflag (optional): interger, if it is non-zero, plot all curves in one panel
%   cfgF: filename of system parameters (for normalization if necessary)
%
% Output:
%   Rms, root-mean-square of defference between original and interpolated seismic stations
%   Cor, correlation between original and interpolated stations
%   *.pdf

% Note(1): here we assume numF=length(sName).
%   Especially for the comparision between observed and interpolated 
%   records, we consider the first file is observed one and second file is
%   interpolated one
%
% Note(2): for data reading, appropriate option is necessary 
%   'b' for data saved on Solaris read into Linux
%   'l' for data saved on Linux read into Linux
%*****=================================================================*****

% Defalut setting, if number numF<=7 (numF=length(sName)), use the following  default colors, otherwise use rand() to generate a random color
dClr=['k','r','b','g','c','m','y'];  % length(sName)<=7

numF=length(sName);

if nargin==1
    cflag=0;
end

% Rms: rms misfit; Cor: correlation degree
Rms=0;
Cor=0;

%=====-----------------------------------------------------------------=====
% Special case: comparision between observed and interpolated stations
if numF==2 && cflag~=0
    % Read information of normalization
    inF=cfgF;
    [ndlongi,ndlati,win,dir,lstFN,outDir,dimType,rbf]=readPar(inF);
    norm0 =win(1,1);
    norm1 =win(1,2);

    % ----------------------------------
    % Read data
    filename=char(sName{1});
    [seis_z(:,1),head(:,1)]=readsac(filename,0,'b');  % Observed data
    filename=char(sName{2});
    [seis_z(:,2),head(:,2)]=readsac(filename,0,'l');  % Interpolted data

    % ----------------------------------
    % DelayT: relative time
    [nt,temp]=size(seis_z(:,1));
    srate=head(1).DELTA; %Sample rate, Unit: s
    DelayT=zeros(1,nt);

    for i=2:nt
        DelayT(i)=DelayT(i-1)+srate;
    end

    % ----------------------------------
    % Normalize real station: seis_z(:,1)
    maxA=max(abs(seis_z(norm0:norm1,1)));%local norm
    %fprintf('maxA=%f file=%s\n',maxA,char(sName{1}));
    seis_z(:,1)=seis_z(:,1)./maxA;

    % ----------------------------------
    % Calculate rms
    Rms=0;
    for i=1:nt
        Rms=Rms+(seis_z(i,1)-seis_z(i,2))^2;
    end
    Rms=sqrt(Rms/nt);

    % ----------------------------------
    % Calculating correlation
    % Note: (1) take a low pass filter; (2) get correlation
    %
    % Version 1: Use corr2 function
    %fs=1./srate;
    %Wp=10; 
    %[b,a]=butter(4,Wp*2/fs,'low');
    %Sig1=filter(b,a,seis_z(:,1));
    %Sig2=filter(b,a,seis_z(:,2));

    %Cor1=corr2(seis_z(:,1),seis_z(:,2))
    %Cor=corr2(Sig1,Sig2);

    % Version 2: Obtain grey correlation
    greyCor=grc(seis_z(:,1),seis_z(:,2))

    % ----------------------------------
    % Plotting
    figure;
    figure('visible','off');

    plot(DelayT,seis_z(:,1),['-',dClr(1)],'LineWidth',1);
    hold on;

    plot(DelayT,seis_z(:,2),['-',dClr(2)],'LineWidth',1);
    xlabel('Delay time (s)');
    title(char(sName{1}));
    legend('Real Station','Synthetic Station');

    str1=sprintf('rms=%f',Rms);
    str2=sprintf('corr=%f',greyCor);
    xxt=DelayT(end)*3/4;
    yyt=max(max(seis_z(:,1)))*3/4;
    text(xxt,yyt,{str1,str2});

    % Save pdf file
	temp=char(sName{2});
	indx=strfind(temp,'.');
    saveas(gcf,[temp(1:indx-1) temp(indx+1:end)],'pdf');
    close(gcf);
    return;
end

%=====-----------------------------------------------------------------=====
% For general cases 
for i=1:numF
    filename=char(sName{i});
    [seis_z(:,i),head(:,i)]=readsac(filename,0,'b');
    %[seis_z(:,i),head(:,i)]=readsac(filename,0,'l');
    if (i<=7) 
        clr=dClr(i);
    else
        clr=rand(1,3);
    end
    [nt,temp]=size(seis_z(:,i));
    srate=head(i).DELTA; % Sample rate, Unit: s
    DelayT=zeros(1,nt);
    for j=2:nt
        DelayT(j)=DelayT(j-1)+srate;
    end

    % ----------------------------------
    if ~cflag
        % Plot curves in separate pdf
        figure;
        figure('visible','off');
        h=plot(DelayT,seis_z(:,i));
        set(h,'color',clr,'LineWidth',1);
        xlabel('Delay time (s)');
        title(filename);

        % Save pdf file
	    indx=strfind(filename,'.');
        saveas(gcf,filename(1:indx-1),'pdf');
        close(gcf);
    else
        % Plot curves in a same pdf
        h=plot(DelayT,seis_z(:,i));
        set(h,'color',clr,'LineWidth',1);
        hold on;
    end
end

if cflag
    xlabel('Delay time (s)');
    % Save pdf file
    tempstr=[num2str(numF),'files'];
    saveas(gcf,tempstr,'pdf');
    close(gcf);
end
