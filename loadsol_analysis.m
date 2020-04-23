% Loadsol Analysis Program
% Written by Alex Peebles (apeebles@vt.edu) at Virginia Tech

% If you use this program for research, please cite this paper: 
%       Luftglass AL, Peebles AT, Miller TK, and Queen RM, "The impact of
%       standardized footwear on load and load symmetry: implications for
%       clinical assessment." Manuscript in preparation as of 4.9.2020

% Please refer to the user manual and tutorial guide for information
% related to how to use this program

function varargout = loadsol_analysis(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loadsol_analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @loadsol_analysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function loadsol_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.fpath = [pwd];
guidata(hObject, handles);

function varargout = loadsol_analysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadData_Callback(hObject, eventdata, handles)
    disp(['Select loadsol Data']);
    [handles.fname,handles.fpath] = uigetfile([handles.fpath,'*.tsv;*.*'],'MultiSelect', 'on');
    clc
    handles.closedtrials = [];
    IC_Thresh = str2num(get(handles.IC_Thresh,'String'));
    TO_Thresh = str2num(get(handles.TO_Thresh,'String'));
    window = str2num(get(handles.EventWindow,'String'));
    handles.nTrials = length(handles.fname);
    set(handles.BW_Enter,'String','Enter weight');
    handles.BW_norm = 0; 
    
    TrialNames = [];
    for t = 1:handles.nTrials
    TrialNames = [TrialNames;'Trial ' num2str(t) ];
    
    % Load data and gap fill
    tempstr = strcat(handles.fpath,string(handles.fname(t)));
    open_LS = importdata(tempstr,' ',4);
    data = [];
    data(:,1:2) = open_LS.data(:,1:2);
    data(:,3) = open_LS.data(:,4);
    [data_out] = GapFill(data);
    
    eval(['handles.T' num2str(t) '_data = data_out;']);
    eval(['handles.T' num2str(t) '_data_r = data_out;']);
    
    % Identify events
    [EventData] = GaitEvents(data_out,IC_Thresh,TO_Thresh,window);    
    eval(['handles.T' num2str(t) '_IC_2 = EventData.IC_events_2;'])
    eval(['handles.T' num2str(t) '_TO_2 = EventData.TO_events_2;'])
    eval(['handles.T' num2str(t) '_IC_3 = EventData.IC_events_3;'])
    eval(['handles.T' num2str(t) '_TO_3 = EventData.TO_events_3;'])
    eval(['handles.T' num2str(t) '_EventFlag = 0;'])
    
    % Load names for each trial
    temp = char(open_LS.textdata(3,1));
    temp2 = find(isspace(temp)==0);
    temp3 = find(diff(temp2)>1);
    eval(['handles.T' num2str(t) '_names(1,:) = temp(temp2(1:temp3));']);
    eval(['handles.T' num2str(t) '_names(2,:) = temp(temp2(temp3+1:end));']);
    
    % Pre-allocate outcomes
    eval(['handles.T' num2str(t) '_PIF_2 = [];'])
    eval(['handles.T' num2str(t) '_PIF_3 = [];'])
    eval(['handles.T' num2str(t) '_PIF_LSI = [];'])
    
    eval(['handles.T' num2str(t) '_IMP_2 = [];'])
    eval(['handles.T' num2str(t) '_IMP_3 = [];'])
    eval(['handles.T' num2str(t) '_IMP_LSI = [];'])
    
    end
    
    LSI_Method_text = {'Difference over average';...
        [handles.T1_names(1,:) ' / ' handles.T1_names(2,:)];...
        [handles.T1_names(2,:) ' / ' handles.T1_names(1,:)]};
    
    handles.name2 = handles.T1_names(1,:);
    handles.name3 = handles.T1_names(2,:); 
    
    handles.SamplingFreq = 1/(mean(diff(handles.T1_data(:,1))));
    set(handles.SampFreqReport,'string',handles.SamplingFreq);
    set(handles.TrialNumber,'string',TrialNames);
    set(handles.LSI_Method,'string',LSI_Method_text);
    handles.comp_flag = 0;
    guidata(hObject,handles);
    plot_data(hObject, eventdata, handles)
        
function [data_out] = GapFill(data)

data2 = data(:,1:2);
data3 = [data(:,1),data(:,3)];

vect2 = find(data2(:,2)==-1);
vect3 = find(data3(:,2)==-1);

if ~isempty(vect2)
vect2t = data2(vect2,1);    
data2(vect2,:) = [];
X2 = data2(:,1);
V2 = data2(:,2);
vq2 = interp1(X2,V2,vect2t,'linear');
data2 = [data2(1:vect2(1)-1,2);vq2;data2(vect2(1):end,2)];
data(:,2) = data2;
end

if ~isempty(vect3)
vect3t = data3(vect3,1);
data3(vect3,:) = [];
X3 = data3(:,1);
V3 = data3(:,2);
vq3 = interp1(X3,V3,vect3t,'linear');
data3 = [data3(1:vect3(1)-1,2);vq3;data3(vect3(1):end,2)];
data(:,3) = data3;
end

data_out = data;

function BW_Enter_Callback(hObject, eventdata, handles)
BW = str2num(get(handles.BW_Enter,'String'));

    for t = 1:handles.nTrials
    eval(['data = handles.T' num2str(t) '_data_r;']);
    data(:,2:3) = data(:,2:3)/BW;
    eval(['handles.T' num2str(t) '_data = data;']);
    end
    handles.BW_norm = 1;
    guidata(hObject,handles);
    plot_data(hObject, eventdata, handles)    
    
function [EventData] = GaitEvents(data,IC_Thresh,TO_Thresh,window)
EventData.IC_events_2 = [];
EventData.TO_events_2 = [];
EventData.IC_events_3 = [];
EventData.TO_events_3 = [];

for i=2:length(data(:,1))-window
%     if i>
%     apple='banana'
    if data(i-1,2) < IC_Thresh & data(i:i+window,2) >= IC_Thresh
        EventData.IC_events_2 = [EventData.IC_events_2;i-1];
    elseif data(i-1,2) > TO_Thresh & data(i:i+window,2) <= TO_Thresh
        EventData.TO_events_2 = [EventData.TO_events_2;i];
    end

    if data(i-1,3) < IC_Thresh & data(i:i+window,3) >= IC_Thresh
        EventData.IC_events_3 = [EventData.IC_events_3;i-1];
    elseif data(i-1,3) > TO_Thresh & data(i:i+window,3) <= TO_Thresh
        EventData.TO_events_3 = [EventData.TO_events_3;i];
    end
    
end

function plot_data(hObject, eventdata, handles)
axes(handles.axes1)
hold off

IC_Thresh = str2num(get(handles.IC_Thresh,'String'));
TO_Thresh = str2num(get(handles.TO_Thresh,'String'));
firstrow = 50;
secondrow = 100;

if handles.BW_norm == 1
    BW = str2num(get(handles.BW_Enter,'String'));
    IC_Thresh = IC_Thresh/BW;
    TO_Thresh = TO_Thresh/BW;
    firstrow = firstrow/BW;
    secondrow = secondrow/BW;
end

contents = cellstr(get(handles.TrialNumber,'String'));
temp = contents{get(handles.TrialNumber,'Value')};

eval(['data = handles.T' temp(end) '_data;'])
eval(['name = handles.T' temp(end) '_names;'])
eval(['IC_events_2 = handles.T' temp(end) '_IC_2;'])
eval(['TO_events_2 = handles.T' temp(end) '_TO_2;'])
eval(['IC_events_3 = handles.T' temp(end) '_IC_3;'])
eval(['TO_events_3 = handles.T' temp(end) '_TO_3;'])

plot(data(:,1),data(:,2),'-r.',data(:,1),data(:,3),'-b.','MarkerSize',10)
hold on
plot(data(:,1),ones(length(data(:,1)),1)*IC_Thresh,'c',...
    data(:,1),ones(length(data(:,1)),1)*TO_Thresh,'y')



if ~isempty(IC_events_2)
plot(data(IC_events_2,1),-firstrow*ones(length(IC_events_2)),'r ^','MarkerSize',10,'LineWidth',2)
end
if ~isempty(TO_events_2)
plot(data(TO_events_2,1),-firstrow*ones(length(TO_events_2)),'r v','MarkerSize',10,'LineWidth',2)
end

if ~isempty(IC_events_3)
plot(data(IC_events_3,1),-secondrow*ones(length(IC_events_3)),'b ^','MarkerSize',10,'LineWidth',2)
end
if ~isempty(TO_events_3)
plot(data(TO_events_3,1),-secondrow*ones(length(TO_events_3)),'b v','MarkerSize',10,'LineWidth',2)
end

if handles.comp_flag == 1
    eval(['PIF2 = handles.plotpeaks_2_' temp(end) ';'])    
    eval(['PIF3 = handles.plotpeaks_3_' temp(end) ';']) 
   
    
    if ~isempty(PIF2)
    plot(data(PIF2(:,1),1),PIF2(:,2),'k o','MarkerSize',10,'LineWidth',2)
    end
    if ~isempty(PIF3)
    plot(data(PIF3(:,1),1),PIF3(:,2),'k o','MarkerSize',10,'LineWidth',2)
    end
end

legend(name(1,:),name(2,:),'IC threshold','TO threshold',[name(1,:) ' IC events'],...
    [name(1,:) ' TO events'],[name(2,:) ' IC events'],[name(2,:) ' TO events'])
   
function IC_Thresh_Callback(hObject, eventdata, handles)
    IC_Thresh = str2num(get(handles.IC_Thresh,'String'));
    TO_Thresh = str2num(get(handles.TO_Thresh,'String'));
    window = str2num(get(handles.EventWindow,'String'));
    if handles.BW_norm == 1
        BW = str2num(get(handles.BW_Enter,'String'));
        IC_Thresh = IC_Thresh/BW;
        TO_Thresh = TO_Thresh/BW;
    end
    for t = 1:handles.nTrials
    eval(['data = handles.T' num2str(t) '_data;']);
    % Identify events
    [EventData] = GaitEvents(data,IC_Thresh,TO_Thresh,window);    
    eval(['handles.T' num2str(t) '_IC_2 = EventData.IC_events_2;'])
    eval(['handles.T' num2str(t) '_TO_2 = EventData.TO_events_2;'])
    eval(['handles.T' num2str(t) '_IC_3 = EventData.IC_events_3;'])
    eval(['handles.T' num2str(t) '_TO_3 = EventData.TO_events_3;'])
    
    end
    guidata(hObject,handles);
    plot_data(hObject, eventdata, handles)

function TO_Thresh_Callback(hObject, eventdata, handles)
    IC_Thresh = str2num(get(handles.IC_Thresh,'String'));
    TO_Thresh = str2num(get(handles.TO_Thresh,'String'));
    window = str2num(get(handles.EventWindow,'String'));
    if handles.BW_norm == 1
        BW = str2num(get(handles.BW_Enter,'String'));
        IC_Thresh = IC_Thresh/BW;
        TO_Thresh = TO_Thresh/BW;
    end
    for t = 1:handles.nTrials
    eval(['data = handles.T' num2str(t) '_data;']);
    % Identify events
    [EventData] = GaitEvents(data,IC_Thresh,TO_Thresh,window);    
    eval(['handles.T' num2str(t) '_IC_2 = EventData.IC_events_2;'])
    eval(['handles.T' num2str(t) '_TO_2 = EventData.TO_events_2;'])
    eval(['handles.T' num2str(t) '_IC_3 = EventData.IC_events_3;'])
    eval(['handles.T' num2str(t) '_TO_3 = EventData.TO_events_3;'])
    
    end
    guidata(hObject,handles);
    plot_data(hObject, eventdata, handles)

function EventWindow_Callback(hObject, eventdata, handles)
    IC_Thresh = str2num(get(handles.IC_Thresh,'String'));
    TO_Thresh = str2num(get(handles.TO_Thresh,'String'));
    window = str2num(get(handles.EventWindow,'String'));
    if handles.BW_norm == 1
        BW = str2num(get(handles.BW_Enter,'String'));
        IC_Thresh = IC_Thresh/BW;
        TO_Thresh = TO_Thresh/BW;
    end
    for t = 1:handles.nTrials
    eval(['data = handles.T' num2str(t) '_data;']);
    % Identify events
    [EventData] = GaitEvents(data,IC_Thresh,TO_Thresh,window);    
    eval(['handles.T' num2str(t) '_IC_2 = EventData.IC_events_2;'])
    eval(['handles.T' num2str(t) '_TO_2 = EventData.TO_events_2;'])
    eval(['handles.T' num2str(t) '_IC_3 = EventData.IC_events_3;'])
    eval(['handles.T' num2str(t) '_TO_3 = EventData.TO_events_3;'])
    
    end
    guidata(hObject,handles);
    plot_data(hObject, eventdata, handles)
        
function Del_IC_Callback(hObject, eventdata, handles)
contents = cellstr(get(handles.TrialNumber,'String'));
temp = contents{get(handles.TrialNumber,'Value')};

eval(['IC_events_2 = handles.T' temp(end) '_IC_2;'])
eval(['IC_events_3 = handles.T' temp(end) '_IC_3;'])

clicks = ginput();
clicks = clicks*handles.SamplingFreq;

for i=1:size(clicks,1)
[val2,ind2] = nanmin(abs(IC_events_2-clicks(i)));
[val3,ind3] = nanmin(abs(IC_events_3-clicks(i)));

if isnan(val2) 
    IC_events_3(ind3) = nan;  
elseif isempty(val2) 
    IC_events_3(ind3) = nan;
elseif isnan(val3) 
    IC_events_2(ind2) = nan; 
elseif isempty(val3)
    IC_events_2(ind2) = nan;  
elseif val2 < val3
    IC_events_2(ind2) = nan;
else
    IC_events_3(ind3) = nan;    
end

end
IC_events_2(isnan(IC_events_2)) = [];
IC_events_3(isnan(IC_events_3)) = [];

eval(['handles.T' temp(end) '_IC_2 = IC_events_2;'])
eval(['handles.T' temp(end) '_IC_3 = IC_events_3;'])
guidata(hObject,handles);     
plot_data(hObject, eventdata, handles)
    
function Del_TO_Callback(hObject, eventdata, handles)
contents = cellstr(get(handles.TrialNumber,'String'));
temp = contents{get(handles.TrialNumber,'Value')};
eval(['TO_events_2 = handles.T' temp(end) '_TO_2;'])
eval(['TO_events_3 = handles.T' temp(end) '_TO_3;'])

clicks = ginput();
clicks = clicks*handles.SamplingFreq;

for i=1:size(clicks,1)
[val2,ind2] = nanmin(abs(TO_events_2-clicks(i)));
[val3,ind3] = nanmin(abs(TO_events_3-clicks(i)));

if isnan(val2) 
    TO_events_3(ind3) = nan;  
elseif isempty(val2) 
    TO_events_3(ind3) = nan;
elseif isnan(val3) 
    TO_events_2(ind2) = nan; 
elseif isempty(val3)
    TO_events_2(ind2) = nan; 
elseif val2 < val3
    TO_events_2(ind2) = nan;
else
    TO_events_3(ind3) = nan;    
end

end
TO_events_2(isnan(TO_events_2)) = [];
TO_events_3(isnan(TO_events_3)) = [];

eval(['handles.T' temp(end) '_TO_2 = TO_events_2;'])
eval(['handles.T' temp(end) '_TO_3 = TO_events_3;'])

guidata(hObject,handles);    
plot_data(hObject, eventdata, handles) 

function Compute_Callback(hObject, eventdata, handles)
MovementType = get(handles.Movement,'Value');

    Trial = [1:handles.nTrials]';
    Trial(handles.closedtrials)=[];
    Trial = num2cell(Trial);
    Trial{end+1} = 'Mean';
    Trial{end+1} = 'StDev';

if MovementType == 1 %Bilateral landing
    
LSI_method = get(handles.LSI_Method,'Value');
    StanceCut = 30; % Eventually read from slider
    TerminalCut = 200; 
    
    c = 1;
    for t = 1:handles.nTrials

    if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)    
    eval(['data = handles.T' num2str(t) '_data;']);
    eval(['IC_events_2 = handles.T' num2str(t) '_IC_2;'])
    eval(['IC_events_3 = handles.T' num2str(t) '_IC_3;'])
    eval(['TO_events_2 = handles.T' num2str(t) '_TO_2;'])
    eval(['TO_events_3 = handles.T' num2str(t) '_TO_3;'])
     
% First landing
complength =  floor((TO_events_2(1) - IC_events_2(1))*.01*StanceCut);
[peaks,locs] = findpeaks(data(IC_events_2(1):IC_events_2(1)+complength,2));
if isempty(locs) %Then just take max
    First_PIF_2(c,1) = IC_events_2(1)+complength
    First_PIF_2(c,2) = data(IC_events_2(1)+complength,2);
else
    [V,I] = sort(peaks,'descend');
    First_PIF_2(c,1) = IC_events_2(1)+locs(I(1));
    First_PIF_2(c,2) = V(1);
end

First_IMP30_2(c) = trapz(data(IC_events_2(1):IC_events_2(1)+complength,1),...
    data(IC_events_2(1):IC_events_2(1)+complength,2));
First_IMP_2(c) = trapz(data(IC_events_2(1):TO_events_2(1),1),...
    data(IC_events_2(1):TO_events_2(1),2));


complength =  floor((TO_events_3(1) - IC_events_3(1))*.01*StanceCut);
[peaks,locs] = findpeaks(data(IC_events_3(1):IC_events_3(1)+complength,3));
if isempty(locs) %Then just take max
    First_PIF_3(c,1) = IC_events_3(1)+complength;
    First_PIF_3(c,2) = data(IC_events_3(1)+complength,3);
else
    [V,I] = sort(peaks,'descend');
    First_PIF_3(c,1) = IC_events_3(1)+locs(I(1));
    First_PIF_3(c,2) = V(1);
end

First_IMP30_3(c) = trapz(data(IC_events_3(1):IC_events_3(1)+complength,1),...
    data(IC_events_3(1):IC_events_3(1)+complength,3));
First_IMP_3(c) = trapz(data(IC_events_3(1):TO_events_3(1),1),...
    data(IC_events_3(1):TO_events_3(1),3));


%Second landing
    complength =  floor(TerminalCut*.001*handles.SamplingFreq);
    [peaks,locs] = findpeaks(data(IC_events_2(2):IC_events_2(2)+complength,2));
    if isempty(locs) %Then just take max
        Second_PIF_2(c,1) = IC_events_2(2)+complength;
        Second_PIF_2(c,2) = data(IC_events_2(2)+complength,2);
    else
        [V,I] = sort(peaks,'descend');
        Second_PIF_2(c,1) = IC_events_2(2)+locs(I(1));
        Second_PIF_2(c,2) = V(1);
    end

    Second_IMP_2(c) = trapz(data(IC_events_2(2):IC_events_2(2)+complength,1),...
        data(IC_events_2(2):IC_events_2(2)+complength,2));


    complength =  floor(TerminalCut*.001*handles.SamplingFreq);
    [peaks,locs] = findpeaks(data(IC_events_3(2):IC_events_3(2)+complength,3));
    if isempty(locs) %Then just take max
        Second_PIF_3(c,1) = IC_events_3(2)+complength;
        Second_PIF_3(c,2) = data(IC_events_3(2)+complength,3);
    else
        [V,I] = sort(peaks,'descend');
        Second_PIF_3(c,1) = IC_events_3(2)+locs(I(1));
        Second_PIF_3(c,2) = V(1);
    end

    Second_IMP_3(c) = trapz(data(IC_events_3(2):IC_events_3(2)+complength,1),...
        data(IC_events_3(2):IC_events_3(2)+complength,3));
    
  eval(['handles.plotpeaks_2_' num2str(t)  ' = [First_PIF_2(c,:);Second_PIF_2(c,:)];'])
  eval(['handles.plotpeaks_3_' num2str(t)  ' = [First_PIF_3(c,:);Second_PIF_3(c,:)];'])
  c=c+1;
    end


    
    end
    

    
    First_PIF_2 = First_PIF_2(:,2);
    Second_PIF_2 = Second_PIF_2(:,2);
    First_PIF_3 = First_PIF_3(:,2);
    Second_PIF_3 = Second_PIF_3(:,2);    
    
    if LSI_method == 1
        First_PIF_LSI = (abs(First_PIF_2 - First_PIF_3)./((First_PIF_2 + First_PIF_3)/2))*100;
        Second_PIF_LSI = (abs(Second_PIF_2 - Second_PIF_3)./((Second_PIF_2 + Second_PIF_3)/2))*100;

        First_IMP30_LSI = (abs(First_IMP30_2 - First_IMP30_3)./((First_IMP30_2 + First_IMP30_3)/2))*100;
        First_IMP_LSI = (abs(First_IMP_2 - First_IMP_3)./((First_IMP_2 + First_IMP_3)/2))*100;
        Second_IMP_LSI = (abs(Second_IMP_2 - Second_IMP_3)./((Second_IMP_2 + Second_IMP_3)/2))*100;

    elseif LSI_method == 2
        First_PIF_LSI = (First_PIF_2./First_PIF_3)*100;
        Second_PIF_LSI = (Second_PIF_2./Second_PIF_3)*100;

        First_IMP30_LSI = (First_IMP30_2./First_IMP30_3)*100;
        First_IMP_LSI = (First_IMP_2./First_IMP_3)*100;
        Second_IMP_LSI = (Second_IMP_2./Second_IMP_3)*100;
    else
        First_PIF_LSI = (First_PIF_3./First_PIF_2)*100;
        Second_PIF_LSI = (Second_PIF_3./Second_PIF_2)*100;

        First_IMP30_LSI = (First_IMP30_3./First_IMP30_2)*100;
        First_IMP_LSI = (First_IMP_3./First_IMP_2)*100;
        Second_IMP_LSI = (Second_IMP_3./Second_IMP_2)*100;
    end
   
    
    First_PIF_2 = [First_PIF_2;mean(First_PIF_2);std(First_PIF_2)];
    First_PIF_3 = [First_PIF_3;mean(First_PIF_3);std(First_PIF_3)];
    First_PIF_LSI = [First_PIF_LSI;mean(First_PIF_LSI);std(First_PIF_LSI)];
    
    First_IMP30_2 = [First_IMP30_2';mean(First_IMP30_2);std(First_IMP30_2)];
    First_IMP30_3 = [First_IMP30_3';mean(First_IMP30_3);std(First_IMP30_3)];
    First_IMP30_LSI = [First_IMP30_LSI';mean(First_IMP30_LSI);std(First_IMP30_LSI)];
    
    First_IMP_2 = [First_IMP_2';mean(First_IMP_2);std(First_IMP_2)];
    First_IMP_3 = [First_IMP_3';mean(First_IMP_3);std(First_IMP_3)];
    First_IMP_LSI = [First_IMP_LSI';mean(First_IMP_LSI);std(First_IMP_LSI)];
    
    Second_PIF_2 = [Second_PIF_2;mean(Second_PIF_2);std(Second_PIF_2)];
    Second_PIF_3 = [Second_PIF_3;mean(Second_PIF_3);std(Second_PIF_3)];
    Second_PIF_LSI = [Second_PIF_LSI;mean(Second_PIF_LSI);std(Second_PIF_LSI)];
    
    Second_IMP_2 = [Second_IMP_2';mean(Second_IMP_2);std(Second_IMP_2)];
    Second_IMP_3 = [Second_IMP_3';mean(Second_IMP_3);std(Second_IMP_3)];
    Second_IMP_LSI = [Second_IMP_LSI';mean(Second_IMP_LSI);std(Second_IMP_LSI)];
    
    eval(['First_PIF_' handles.name2(end) '=First_PIF_2;'])
    eval(['First_PIF_' handles.name3(end) '=First_PIF_3;'])
    eval(['First_IMP30_' handles.name2(end) '=First_IMP30_2;'])
    eval(['First_IMP30_' handles.name3(end) '=First_IMP30_3;'])
    eval(['First_IMP_' handles.name2(end) '=First_IMP_2;'])
    eval(['First_IMP_' handles.name3(end) '=First_IMP_3;'])
    
    eval(['Second_PIF_' handles.name2(end) '=Second_PIF_2;'])
    eval(['Second_PIF_' handles.name3(end) '=Second_PIF_3;'])
    eval(['Second_IMP_' handles.name2(end) '=Second_IMP_2;'])
    eval(['Second_IMP_' handles.name3(end) '=Second_IMP_3;'])   
                                            

    eval(['T = table(Trial,First_PIF_' handles.name2(end) ',First_PIF_' handles.name3(end) ',First_PIF_LSI,'...
        'First_IMP30_' handles.name2(end) ',First_IMP30_' handles.name3(end) ',First_IMP30_LSI,'...
        'First_IMP_' handles.name2(end) ',First_IMP_' handles.name3(end) ',First_IMP_LSI,'...
        'Second_PIF_' handles.name2(end) ',Second_PIF_' handles.name3(end) ',Second_PIF_LSI,'...
        'Second_IMP_' handles.name2(end) ',Second_IMP_' handles.name3(end) ',Second_IMP_LSI)'])

    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));
    
    writetable(T,filename,'Sheet',trialname)
    
    
elseif MovementType == 2 %Unilateral landing
        StanceCut = 30; % Eventually read from slider
    TerminalCut = 200; 
    
    c = 1;
    for t = 1:handles.nTrials

    if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)    
    eval(['data = handles.T' num2str(t) '_data;']);
    eval(['IC_events_2 = handles.T' num2str(t) '_IC_2;'])
    eval(['IC_events_3 = handles.T' num2str(t) '_IC_3;'])
    eval(['TO_events_2 = handles.T' num2str(t) '_TO_2;'])
    eval(['TO_events_3 = handles.T' num2str(t) '_TO_3;'])
    
    if length(IC_events_2)>length(IC_events_3) % Analyze 2
     
    complength =  floor(TerminalCut*.001*handles.SamplingFreq);
    [peaks,locs] = findpeaks(data(IC_events_2(1):IC_events_2(1)+complength,2));
    if isempty(locs) %Then just take max
        PIF_2(c,1) = IC_events_2(1)+complength;
        PIF_2(c,2) = data(IC_events_2(1)+complength,2);
    else
        [V,I] = sort(peaks,'descend');
        PIF_2(c,1) = IC_events_2(1)+locs(I(1));
        PIF_2(c,2) = V(1);
    end

    IMP_2(c) = trapz(data(IC_events_2(1):IC_events_2(1)+complength,1),...
        data(IC_events_2(1):IC_events_2(1)+complength,2));
    
  eval(['handles.plotpeaks_2_' num2str(t)  ' = [PIF_2(c,:)];'])
  eval(['handles.plotpeaks_3_' num2str(t)  ' = [];'])
 c=c+1;
    else
        
    complength =  floor(TerminalCut*.001*handles.SamplingFreq);
    [peaks,locs] = findpeaks(data(IC_events_3(1):IC_events_3(1)+complength,3));
    if isempty(locs) %Then just take max
        PIF_3(c,1) = IC_events_3(1)+complength;
        PIF_3(c,2) = data(IC_events_3(1)+complength,3);
    else
        [V,I] = sort(peaks,'descend');
        PIF_3(c,1) = IC_events_3(1)+locs(I(1));
        PIF_3(c,2) = V(1);
    end

    IMP_3(c) = trapz(data(IC_events_3(1):IC_events_3(1)+complength,1),...
        data(IC_events_3(1):IC_events_3(1)+complength,3));
    
  eval(['handles.plotpeaks_2_' num2str(t)  ' = [];'])
  eval(['handles.plotpeaks_3_' num2str(t)  ' = [PIF_3(c,:)];'])
  c=c+1;
    end
end

end
    



    if length(IC_events_2)>length(IC_events_3)
    PIF_2 = [PIF_2(:,2);mean(PIF_2(:,2));std(PIF_2(:,2))];
    IMP_2 = [IMP_2';mean(IMP_2);std(IMP_2)];
    eval(['PIF_' handles.name2(end) '=PIF_2;'])
    eval(['IMP_' handles.name2(end) '=IMP_2;'])
    eval(['T = table(Trial,PIF_' handles.name2(end) ',IMP_' handles.name2(end) ')']) 
    
    else
        
    PIF_3 = [PIF_3(:,2);mean(PIF_3(:,2));std(PIF_3(:,2))];
    IMP_3 = [IMP_3';mean(IMP_3);std(IMP_3)];
    eval(['PIF_' handles.name2(end) '=PIF_3;'])
    eval(['IMP_' handles.name2(end) '=IMP_3;'])
    eval(['T = table(Trial,PIF_' handles.name2(end) ',IMP_' handles.name2(end) ')'])   
    end                                     



    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));
    
    writetable(T,filename,'Sheet',trialname)
    
elseif MovementType == 3 %Running
    
end

handles.comp_flag = 1;
guidata(hObject,handles);
plot_data(hObject, eventdata, handles)
 
function TrialNumber_Callback(hObject, eventdata, handles)
    plot_data(hObject, eventdata, handles)

function Zoom_Callback(hObject, eventdata, handles)
zoom on

function Pan_Callback(hObject, eventdata, handles)
pan on

function Cursor_Callback(hObject, eventdata, handles)
datacursormode on

function CloseTrial_Callback(hObject, eventdata, handles)
contents = cellstr(get(handles.TrialNumber,'String'));
temp = contents{get(handles.TrialNumber,'Value')};
trial_close = temp(end);
contents = cellstr(get(handles.TrialNumber,'string'));

contents(str2num(trial_close))=[];
% contents(7)=[];

set(handles.TrialNumber,'value',1);
set(handles.TrialNumber,'string',contents);
handles.closedtrials = [handles.closedtrials;str2num(temp(end))];
guidata(hObject,handles);
plot_data(hObject, eventdata, handles)

function Filter_Callback(hObject, eventdata, handles)
contents = cellstr(get(handles.Filter_type,'String'));
temp = contents{get(handles.Filter_type,'Value')};

if temp(1) == 'L'
    filttype = 'low'
else
    filttype = 'High'
end
    

fc = str2num(get(handles.Filter_freq,'String'));
fs = round(handles.SamplingFreq); 

[b,a] = butter(4,fc/(fs/2),filttype);

for t = 1:handles.nTrials
eval(['data = handles.T' num2str(t) '_data;']);
eval(['handles.T' num2str(t) '_Stored_data = data;']);

data(:,2:3) = round(filtfilt(b,a,data(:,2:3)));

eval(['handles.T' num2str(t) '_data = data;']);
end


guidata(hObject,handles);
plot_data(hObject, eventdata, handles)

function Unfilter_Callback(hObject, eventdata, handles)
for t = 1:handles.nTrials
eval(['handles.T' num2str(t) '_data = handles.T' num2str(t) '_Stored_data;']);
end

guidata(hObject,handles);
plot_data(hObject, eventdata, handles)




% Future need maybe? 
function Save_Callback(hObject, eventdata, handles)

function SubMin_Callback(hObject, eventdata, handles)

function Mod_TO_Callback(hObject, eventdata, handles)

function Mod_IC_Callback(hObject, eventdata, handles)

function Mod_IP_Callback(hObject, eventdata, handles)




% Not needed
function IC_Thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IC_Thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
function TO_Thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TO_Thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TrialNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrialNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Filter_type_Callback(hObject, eventdata, handles)

function Filter_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filter_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Filter_freq_Callback(hObject, eventdata, handles)

function Filter_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filter_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EventWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EventWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LSI_Method_Callback(hObject, eventdata, handles)

function LSI_Method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LSI_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BW_Enter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BW_Enter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Movement.
function Movement_Callback(hObject, eventdata, handles)
% hObject    handle to Movement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Movement contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Movement


% --- Executes during object creation, after setting all properties.
function Movement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Movement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
