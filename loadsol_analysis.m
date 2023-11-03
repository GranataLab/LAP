
% Load Analysis Program
% Written by Alex Peebles (apeebles@vt.edu) at Virginia Tech

% If you use this program for research, please cite this paper: 
%       Luftglass AL, Peebles AT, Miller TK, and Queen RM, "The impact of
%       standardized footwear on load and load symmetry: implications for
%       clinical assessment." Manuscript in preparation as of 4.9.2020

% Please refer to the user manual and tutorial guide for information
% related to how to use this program
%% Code Adjustments
% Expanded on by Nicole Stark (nestark@vt.edu) at Virginia Tech 

% 5/24/21
% The code was expanded to include GRF curves, 1st peak (p-GRF), 2nd
% peak (WA-GRF), and Stance time pre a trial 
% 6/1/21 
% This code was expanded on again to output GRF curves into a secondary
% excel

%Theresa - change on line 78 to accommodate for new loadsol files
%% 
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
    if iscell(handles.fname)
        handles.nTrials = size(handles.fname,2);
    else
        handles.nTrials = 1;
    end
     
    set(handles.BW_Enter,'String','Enter weight');
    handles.BW_norm = 0; 
    
    TrialNames = [];
    for t = 1:handles.nTrials
    TrialNames = [TrialNames;'Trial ' num2str(t) ];
    
    % Load data and gap fill
    if iscell(handles.fname)
        tempstr = strcat(handles.fpath,string(handles.fname(:,t)));
    else
        tempstr = strcat(handles.fpath,string(handles.fname));
    end
    
    open_LS = importdata(tempstr); %,' ',4 -- taken out to accommodate for new loadsol files
    data = [];
    data(:,1:2) = open_LS.data(:,1:2);
    %data(:,3) = open_LS.data(:,4);
    %comment the line above and use the line below when processing control subjects (switch for ACLR) 
    data(:,3) = open_LS.data(:,6);
    [data_out] = GapFill(data);
    
    eval(['handles.T' num2str(t) '_data = data_out;']);
    eval(['handles.T' num2str(t) '_data_r = data_out;']);
    eval(['handles.T' num2str(t) '_Stored_data = data_out;']);
    
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
    eval(['handles.plotpeaks_2_' num2str(t) ' = [];'])
    eval(['handles.plotpeaks_3_' num2str(t) ' = [];'])
    end
    
    DominantLimb = {handles.T1_names(1,:);handles.T1_names(2,:)};

    
    handles.name2 = handles.T1_names(1,:);
    handles.name3 = handles.T1_names(2,:); 
    
    handles.SamplingFreq = 1/(mean(diff(handles.T1_data(:,1))));
    set(handles.SampFreqReport,'string',handles.SamplingFreq);
    set(handles.TrialNumber,'string',TrialNames);
    set(handles.DominantLimb,'string',DominantLimb);
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

if eval(['~isempty(handles.plotpeaks_2_' temp(end) ');'])
   eval(['PIF2 = handles.plotpeaks_2_' temp(end) ';']) 
   plot(data(PIF2(:,1),1),PIF2(:,2),'k o','MarkerSize',10,'LineWidth',2)
end

if eval(['~isempty(handles.plotpeaks_3_' temp(end) ');'])
   eval(['PIF3 = handles.plotpeaks_3_' temp(end) ';']) 
   plot(data(PIF3(:,1),1),PIF3(:,2),'k o','MarkerSize',10,'LineWidth',2)
end

legend(name(1,:),name(2,:))
   
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

StanceCut_IMP = round(get(handles.PercentGCP_IMP, 'Value')*100);
StanceCut_PIF = round(get(handles.PercentGCP_PIF, 'Value')*100);
TimeCut_IMP = round(get(handles.TimeFC_IMP, 'Value')*100);
TimeCut_PIF = round(get(handles.TimeFC_PIF, 'Value')*100);


if MovementType == 1 %Bilateral landing
    
Trial = [1:handles.nTrials]';
Trial(handles.closedtrials)=[];
Trial = num2cell(Trial);
Trial{end+1} = 'Mean';
Trial{end+1} = 'StDev';

    c = 1;
    for t = 1:handles.nTrials

        if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)
            eval(['data = handles.T' num2str(t) '_data;']);
            eval(['IC_events_2 = handles.T' num2str(t) '_IC_2;'])
            eval(['IC_events_3 = handles.T' num2str(t) '_IC_3;'])
            eval(['TO_events_2 = handles.T' num2str(t) '_TO_2;'])
            eval(['TO_events_3 = handles.T' num2str(t) '_TO_3;'])

            % First landing
            
            if isempty(TO_events_2) %No toe of for this step, use time
                
            PIFLength2 =  floor(TimeCut_PIF*.001*handles.SamplingFreq);
            IMPLength2 =  floor(TimeCut_IMP*.001*handles.SamplingFreq);
            PIFLength3 =  floor(TimeCut_PIF*.001*handles.SamplingFreq);
            IMPLength3 =  floor(TimeCut_IMP*.001*handles.SamplingFreq);
                
            else %There is a TO, use %stance
                
            PIFLength2 =  floor((TO_events_2(1) - IC_events_2(1))*.01*StanceCut_PIF);
            IMPLength2 =  floor((TO_events_2(1) - IC_events_2(1))*.01*StanceCut_IMP);
            PIFLength3 =  floor((TO_events_3(1) - IC_events_3(1))*.01*StanceCut_PIF);
            IMPLength3 =  floor((TO_events_3(1) - IC_events_3(1))*.01*StanceCut_IMP);
            
            end
            
            
            [peaks,locs] = findpeaks(data(IC_events_2(1):IC_events_2(1)+PIFLength2,2));
            if isempty(locs) %Then just take max
                First_PIF_2(c,1) = IC_events_2(1)+PIFLength2
                First_PIF_2(c,2) = data(IC_events_2(1)+PIFLength2,2);
                First_LR_2(c) = First_PIF_2(c,2)/PIFLength2*handles.SamplingFreq;
            else
                [V,I] = sort(peaks,'descend');
                First_PIF_2(c,1) = IC_events_2(1)+locs(I(1));
                First_PIF_2(c,2) = V(1);
                First_LR_2(c) = First_PIF_2(c,2)/locs(I(1))*handles.SamplingFreq;
            end

            First_IMP_2(c) = trapz(data(IC_events_2(1):IC_events_2(1)+IMPLength2,1),...
                data(IC_events_2(1):IC_events_2(1)+IMPLength2,2));

            [peaks,locs] = findpeaks(data(IC_events_3(1):IC_events_3(1)+PIFLength3,3));
            if isempty(locs) %Then just take max
                First_PIF_3(c,1) = IC_events_3(1)+PIFLength3;
                First_PIF_3(c,2) = data(IC_events_3(1)+PIFLength3,3);
                First_LR_3(c) = First_PIF_3(c,2)/PIFLength3*handles.SamplingFreq;
            else
                [V,I] = sort(peaks,'descend');
                First_PIF_3(c,1) = IC_events_3(1)+locs(I(1));
                First_PIF_3(c,2) = V(1);
                First_LR_3(c) = First_PIF_3(c,2)/locs(I(1))*handles.SamplingFreq;
            end

            First_IMP_3(c) = trapz(data(IC_events_3(1):IC_events_3(1)+IMPLength3,1),...
                data(IC_events_3(1):IC_events_3(1)+IMPLength3,3));
            
            %Second landing
            if length(IC_events_2)>1 && length(IC_events_3)>1
            
            PIFLength =  floor(TimeCut_PIF*.001*handles.SamplingFreq);
            IMPLength =  floor(TimeCut_IMP*.001*handles.SamplingFreq);
            
            [peaks,locs] = findpeaks(data(IC_events_2(2):IC_events_2(2)+PIFLength,2));
            if isempty(locs) %Then just take max
                Second_PIF_2(c,1) = IC_events_2(2)+PIFLength;
                Second_PIF_2(c,2) = data(IC_events_2(2)+PIFLength,2);
                Second_LR_2(c) = Second_PIF_2(c,2)/PIFLength*handles.SamplingFreq;
            else
                [V,I] = sort(peaks,'descend');
                Second_PIF_2(c,1) = IC_events_2(2)+locs(I(1));
                Second_PIF_2(c,2) = V(1);
                Second_LR_2(c) = Second_PIF_2(c,2)/locs(I(1))*handles.SamplingFreq;
            end
            
            
            Second_IMP_2(c) = trapz(data(IC_events_2(2):IC_events_2(2)+IMPLength,1),...
                data(IC_events_2(2):IC_events_2(2)+IMPLength,2));

            
            [peaks,locs] = findpeaks(data(IC_events_3(2):IC_events_3(2)+PIFLength,3));
            if isempty(locs) %Then just take max
                Second_PIF_3(c,1) = IC_events_3(2)+PIFLength;
                Second_PIF_3(c,2) = data(IC_events_3(2)+PIFLength,3);
                Second_LR_3(c) = Second_PIF_3(c,2)/PIFLength*handles.SamplingFreq;
            else
                [V,I] = sort(peaks,'descend');
                Second_PIF_3(c,1) = IC_events_3(2)+locs(I(1));
                Second_PIF_3(c,2) = V(1);
                Second_LR_3(c) = Second_PIF_3(c,2)/locs(I(1))*handles.SamplingFreq;
            end

            Second_IMP_3(c) = trapz(data(IC_events_3(2):IC_events_3(2)+IMPLength,1),...
                data(IC_events_3(2):IC_events_3(2)+IMPLength,3));

            eval(['handles.plotpeaks_2_' num2str(t)  ' = [First_PIF_2(c,:);Second_PIF_2(c,:)];'])
            eval(['handles.plotpeaks_3_' num2str(t)  ' = [First_PIF_3(c,:);Second_PIF_3(c,:)];'])
            else
            eval(['handles.plotpeaks_2_' num2str(t)  ' = [First_PIF_2(c,:)];'])
            eval(['handles.plotpeaks_3_' num2str(t)  ' = [First_PIF_3(c,:)];'])                
            end
         
            c=c+1;
        end
    end



    First_PIF_2 = First_PIF_2(:,2)';
    First_PIF_3 = First_PIF_3(:,2)';
    
    if length(IC_events_2)>1 && length(IC_events_3)>1
     Second_PIF_2 = Second_PIF_2(:,2)';   
     Second_PIF_3 = Second_PIF_3(:,2)';
    end
    
    LSI_method = get(handles.LSI_Method,'Value');
    if LSI_method == 1%LSI
        
        if get(handles.DominantLimb,'Value') == 1
        First_PIF_LSI = (First_PIF_3./First_PIF_2)*100;
        First_LR_LSI = (First_LR_3./First_LR_2)*100;
        First_IMP_LSI = (First_IMP_3./First_IMP_2)*100;
        if length(IC_events_2)>1 && length(IC_events_3)>1
        Second_PIF_LSI = (Second_PIF_3./Second_PIF_2)*100;
        Second_LR_LSI = (Second_LR_3./Second_LR_2)*100;
        Second_IMP_LSI = (Second_IMP_3./Second_IMP_2)*100;
        end
        else
        First_PIF_LSI = (First_PIF_2./First_PIF_3)*100;
        First_LR_LSI = (First_LR_2./First_LR_3)*100;
        First_IMP_LSI = (First_IMP_2./First_IMP_3)*100;
        if length(IC_events_2)>1 && length(IC_events_3)>1
        Second_PIF_LSI = (Second_PIF_2./Second_PIF_3)*100;
        Second_LR_LSI = (Second_LR_2./Second_LR_3)*100;
        Second_IMP_LSI = (Second_IMP_2./Second_IMP_3)*100;  
        end
        end
        
    elseif LSI_method == 2 %ASI
        
        First_PIF_LSI = (abs(First_PIF_2 - First_PIF_3)./((First_PIF_2 + First_PIF_3)/2))*100;
        First_LR_LSI = (abs(First_LR_2 - First_LR_3)./((First_LR_2 + First_LR_3)/2))*100;
        First_IMP_LSI = (abs(First_IMP_2 - First_IMP_3)./((First_IMP_2 + First_IMP_3)/2))*100;
        if length(IC_events_2)>1 && length(IC_events_3)>1
        Second_PIF_LSI = (abs(Second_PIF_2 - Second_PIF_3)./((Second_PIF_2 + Second_PIF_3)/2))*100;
        Second_LR_LSI = (abs(Second_LR_2 - Second_LR_3)./((Second_LR_2 + Second_LR_3)/2))*100;
        Second_IMP_LSI = (abs(Second_IMP_2 - Second_IMP_3)./((Second_IMP_2 + Second_IMP_3)/2))*100;
        end
        
    else %NSI
        
        for j=1:length(First_IMP_2)
            if get(handles.DominantLimb,'Value') == 1 %2 is dom
                First_PIF_LSI(j) = 100*(First_PIF_2(j) - First_PIF_3(j))/max([First_PIF_2,First_PIF_3]);
                First_LR_LSI(j) = 100*(First_LR_2(j) - First_LR_3(j))/max([First_LR_2,First_LR_3]);
                First_IMP_LSI(j) = 100*(First_IMP_2(j) - First_IMP_3(j))/max([First_IMP_2,First_IMP_3]);
                if length(IC_events_2)>1 && length(IC_events_3)>1
                Second_PIF_LSI(j) = 100*(Second_PIF_2(j) - Second_PIF_3(j))/max([Second_PIF_2,Second_PIF_3]);
                Second_LR_LSI(j) = 100*(Second_LR_2(j) - Second_LR_3(j))/max([Second_LR_2,Second_LR_3]);
                Second_IMP_LSI(j) = 100*(Second_IMP_2(j) - Second_IMP_3(j))/max([Second_IMP_2,Second_IMP_3]);
                end
            else
                First_PIF_LSI(j) = 100*(First_PIF_3(j) - First_PIF_2(j))/max([First_PIF_2,First_PIF_3]);
                First_LR_LSI(j) = 100*(First_LR_3(j) - First_LR_2(j))/max([First_LR_2,First_LR_3]);
                First_IMP_LSI(j) = 100*(First_IMP_3(j) - First_IMP_2(j))/max([First_IMP_2,First_IMP_3]);
                if length(IC_events_2)>1 && length(IC_events_3)>1
                Second_PIF_LSI(j) = 100*(Second_PIF_3(j) - Second_PIF_2(j))/max([Second_PIF_2,Second_PIF_3]);
                Second_LR_LSI(j) = 100*(Second_LR_3(j) - Second_LR_2(j))/max([Second_LR_2,Second_LR_3]);
                Second_IMP_LSI(j) = 100*(Second_IMP_3(j) - Second_IMP_2(j))/max([Second_IMP_2,Second_IMP_3]);
                end
            end
        end

    end


    First_PIF_2 = [First_PIF_2';mean(First_PIF_2);std(First_PIF_2)];
    First_PIF_3 = [First_PIF_3';mean(First_PIF_3);std(First_PIF_3)];
    First_PIF_LSI = [First_PIF_LSI';mean(First_PIF_LSI);std(First_PIF_LSI)];
    
    First_LR_2 = [First_LR_2';mean(First_LR_2);std(First_LR_2)];
    First_LR_3 = [First_LR_3';mean(First_LR_3);std(First_LR_3)];
    First_LR_LSI = [First_LR_LSI';mean(First_LR_LSI);std(First_LR_LSI)];
    
    First_IMP_2 = [First_IMP_2';mean(First_IMP_2);std(First_IMP_2)];
    First_IMP_3 = [First_IMP_3';mean(First_IMP_3);std(First_IMP_3)];
    First_IMP_LSI = [First_IMP_LSI';mean(First_IMP_LSI);std(First_IMP_LSI)];
    
    eval(['First_PIF_' handles.name2(end) '=First_PIF_2;'])
    eval(['First_PIF_' handles.name3(end) '=First_PIF_3;'])
    eval(['First_LR_' handles.name2(end) '=First_LR_2;'])
    eval(['First_LR_' handles.name3(end) '=First_LR_3;'])
    eval(['First_IMP_' handles.name2(end) '=First_IMP_2;'])
    eval(['First_IMP_' handles.name3(end) '=First_IMP_3;'])
    
    
    if length(IC_events_2)>1 && length(IC_events_3)>1
    Second_PIF_2 = [Second_PIF_2';mean(Second_PIF_2);std(Second_PIF_2)];
    Second_PIF_3 = [Second_PIF_3';mean(Second_PIF_3);std(Second_PIF_3)];
    Second_PIF_LSI = [Second_PIF_LSI';mean(Second_PIF_LSI);std(Second_PIF_LSI)];

    Second_LR_2 = [Second_LR_2';mean(Second_LR_2);std(Second_LR_2)];
    Second_LR_3 = [Second_LR_3';mean(Second_LR_3);std(Second_LR_3)];
    Second_LR_LSI = [Second_LR_LSI';mean(Second_LR_LSI);std(Second_LR_LSI)];
    
    Second_IMP_2 = [Second_IMP_2';mean(Second_IMP_2);std(Second_IMP_2)];
    Second_IMP_3 = [Second_IMP_3';mean(Second_IMP_3);std(Second_IMP_3)];
    Second_IMP_LSI = [Second_IMP_LSI';mean(Second_IMP_LSI);std(Second_IMP_LSI)];

    eval(['Second_PIF_' handles.name2(end) '=Second_PIF_2;'])
    eval(['Second_PIF_' handles.name3(end) '=Second_PIF_3;'])
    eval(['Second_LR_' handles.name2(end) '=Second_LR_2;'])
    eval(['Second_LR_' handles.name3(end) '=Second_LR_3;'])
    eval(['Second_IMP_' handles.name2(end) '=Second_IMP_2;'])
    eval(['Second_IMP_' handles.name3(end) '=Second_IMP_3;'])
    
    eval(['T = table(Trial,First_PIF_' handles.name2(end) ',First_PIF_' handles.name3(end) ',First_PIF_LSI,'...
        'First_LR_' handles.name2(end) ',First_LR_' handles.name3(end) ',First_LR_LSI,'...
        'First_IMP_' handles.name2(end) ',First_IMP_' handles.name3(end) ',First_IMP_LSI,'...
        'Second_PIF_' handles.name2(end) ',Second_PIF_' handles.name3(end) ',Second_PIF_LSI,'...
        'Second_LR_' handles.name2(end) ',Second_LR_' handles.name3(end) ',Second_LR_LSI,'...
        'Second_IMP_' handles.name2(end) ',Second_IMP_' handles.name3(end) ',Second_IMP_LSI)'])
    else
    eval(['T = table(Trial,First_PIF_' handles.name2(end) ',First_PIF_' handles.name3(end) ',First_PIF_LSI,'...
        'First_LR_' handles.name2(end) ',First_LR_' handles.name3(end) ',First_LR_LSI,'...
        'First_IMP_' handles.name2(end) ',First_IMP_' handles.name3(end) ',First_IMP_LSI)'])
    end

    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));

    writetable(T,filename,'Sheet',trialname)


elseif MovementType == 2 %Unilateral landing
PIFLength =  floor(TimeCut_PIF*.001*handles.SamplingFreq);
IMPLength =  floor(TimeCut_IMP*.001*handles.SamplingFreq);    
c2 = 1;
c3 = 1;
for t = 1:handles.nTrials

if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)
    eval(['data = handles.T' num2str(t) '_data;']);
    eval(['IC_events_2 = handles.T' num2str(t) '_IC_2;'])
    eval(['IC_events_3 = handles.T' num2str(t) '_IC_3;'])
    eval(['TO_events_2 = handles.T' num2str(t) '_TO_2;'])
    eval(['TO_events_3 = handles.T' num2str(t) '_TO_3;'])
    
                if length(IC_events_2)>length(IC_events_3) % Analyze 2

                [peaks,locs] = findpeaks(data(IC_events_2(1):IC_events_2(1)+PIFLength,2));
                if isempty(locs) %Then just take max
                    PIF_2(c2,1) = IC_events_2(1)+PIFLength;
                    PIF_2(c2,2) = data(IC_events_2(1)+PIFLength,2);
                    LR_2(c2) = PIF_2(c2,2)/PIFLength*handles.SamplingFreq;
                else
                    [V,I] = sort(peaks,'descend');
                    PIF_2(c2,1) = IC_events_2(1)+locs(I(1));
                    PIF_2(c2,2) = V(1);
                    LR_2(c2) = PIF_2(c2,2)/locs(I(1))*handles.SamplingFreq;
                end
                
                time_TO_2(c2)=data(TO_events_2,1);
                time_IC_2(c2)=data(IC_events_2,1);
                %time_PIF_2(c2)=data(PIF_2(1),1);
                                
                IMP_2(c2) = trapz(data(IC_events_2(1):IC_events_2(1)+IMPLength,1),...
                    data(IC_events_2(1):IC_events_2(1)+IMPLength,2));
                
                Curve2(:,c2) = data(IC_events_2(1):IC_events_2(1)+PIFLength,2);
                
                eval(['handles.plotpeaks_2_' num2str(t)  '=PIF_2(c2,:);']);
                c2=c2+1;
                
                else
                    
                [peaks,locs] = findpeaks(data(IC_events_3(1):IC_events_3(1)+PIFLength,3));
                if isempty(locs) %Then just take max
                    PIF_3(c3,1) = IC_events_3(1)+PIFLength;
                    PIF_3(c3,2) = data(IC_events_3(1)+PIFLength,3);
                    LR_3(c3) = PIF_3(c3,2)/PIFLength*handles.SamplingFreq;
                else
                    [V,I] = sort(peaks,'descend');
                    PIF_3(c3,1) = IC_events_3(1)+locs(I(1));
                    PIF_3(c3,2) = V(1);
                    LR_3(c3) = PIF_3(c3,2)/locs(I(1))*handles.SamplingFreq;
                end
                
                time_TO_3(c3)=data(TO_events_3,1);
                time_IC_3(c3)=data(IC_events_3,1);
                %time_PIF_3(c3)=data(PIF_3(1),1);
                
                IMP_3(c3) = trapz(data(IC_events_3(1):IC_events_3(1)+IMPLength,1),...
                    data(IC_events_3(1):IC_events_3(1)+IMPLength,3));

                Curve3(:,c3) = data(IC_events_3(1):IC_events_3(1)+PIFLength,3);

                eval(['handles.plotpeaks_3_' num2str(t)  '=PIF_3(c3,:);']);
                c3=c3+1;
                end
                
                %outputting PPF and time points for TO and IC
                
                
end

end

    PIF_2 = PIF_2(:,2);
    PIF_3 = PIF_3(:,2);
    
    if length(PIF_2)>length(PIF_3)
        trim2 = length(PIF_2)-length(PIF_3);
        trim3 = 0;
    else
        trim2 = length(PIF_2)-length(PIF_3);
        trim3 = 0;        
    end
    
    
    LSI_method = get(handles.LSI_Method,'Value');
    if LSI_method == 1%LSI
        
        if get(handles.DominantLimb,'Value') == 1
        PIF_Symm = (mean(PIF_3(1:end-trim3))./mean(PIF_2(1:end-trim2)))*100;
        LR_Symm = (mean(LR_3(1:end-trim3))./mean(LR_2(1:end-trim2)))*100;
        IMP_Symm = (mean(IMP_3(1:end-trim3))./mean(IMP_2(1:end-trim2)))*100;
        else
        PIF_Symm = (mean(PIF_2(1:end-trim2))./mean(PIF_3(1:end-trim3)))*100;
        LR_Symm = (mean(LR_2(1:end-trim2))./mean(LR_3(1:end-trim3)))*100;
        IMP_Symm = (mean(IMP_2(1:end-trim2))./mean(IMP_3(1:end-trim3)))*100;
        end
        
    elseif LSI_method == 2 %ASI
        
        PIF_Symm = (abs(mean(PIF_2(1:end-trim2)) - mean(PIF_3(1:end-trim3)))./((mean(PIF_2(1:end-trim2)) + mean(PIF_3(1:end-trim3)))/2))*100;
        LR_Symm = (abs(mean(LR_2(1:end-trim2)) - mean(LR_3(1:end-trim3)))./((mean(LR_2(1:end-trim2)) + mean(LR_3(1:end-trim3)))/2))*100;
        IMP_Symm = (abs(mean(IMP_2(1:end-trim2)) - mean(IMP_3(1:end-trim3)))./((mean(IMP_2(1:end-trim2)) + mean(IMP_3(1:end-trim3)))/2))*100;

        
    else %NSI
        
        for j=1:min([length(PIF_2),length(PIF_3)])
            if get(handles.DominantLimb,'Value') == 1 %2 is dom
                PIF_Symm(j) = 100*(PIF_2(j) - PIF_3(j))/max([PIF_2;PIF_3]);
                ILR_Symm(j) = 100*(LR_2(j) - LR_3(j))/max([LR_2;LR_3]);
                IMP_Symm(j) = 100*(IMP_2(j) - IMP_3(j))/max([IMP_2;IMP_3]);
            else
                PIF_Symm(j) = 100*(PIF_3(j) - PIF_2(j))/max([PIF_2;PIF_3]);
                ILR_Symm(j) = 100*(LR_3(j) - LR_2(j))/max([LR_2;LR_3]);
                IMP_Symm(j) = 100*(IMP_3(j) - IMP_2(j))/max([IMP_2;IMP_3]);
            end
        end

    end
    
    Step = [1:min([length(PIF_2),length(PIF_3)])]';
    Step = num2cell(Step);
    Step{end+1} = 'Mean';
    Step{end+1} = 'StDev';
    
    if LSI_method == 3 %NSI
    PIF_Symm = [PIF_Symm;mean(PIF_Symm);std(PIF_Symm)];
    LR_Symm = [LR_Symm';mean(LR_Symm);std(LR_Symm)];
    IMP_Symm = [IMP_Symm';mean(IMP_Symm);std(IMP_Symm)];
    else
    PIF_Symm = [nan(min([length(PIF_2),length(PIF_3)]),1);PIF_Symm;nan];
    LR_Symm = [nan(min([length(PIF_2),length(PIF_3)]),1);LR_Symm;nan];
    IMP_Symm = [nan(min([length(PIF_2),length(PIF_3)]),1);IMP_Symm;nan];
    end

    PIF_2 = [PIF_2(1:end-trim2);mean(PIF_2(1:end-trim2));std(PIF_2(1:end-trim2))];
    PIF_3 = [PIF_3(1:end-trim3);mean(PIF_3(1:end-trim3));std(PIF_3(1:end-trim3))];

    LR_2 = [LR_2(1:end-trim2)';mean(LR_2(1:end-trim2));std(LR_2(1:end-trim2))];
    LR_3 = [LR_3(1:end-trim3)';mean(LR_3(1:end-trim3));std(LR_3(1:end-trim3))];
    
    IMP_2 = [IMP_2(1:end-trim2)';mean(IMP_2(1:end-trim2));std(IMP_2(1:end-trim2))];
    IMP_3 = [IMP_3(1:end-trim3)';mean(IMP_3(1:end-trim3));std(IMP_3(1:end-trim3))];
   
    time_TO_2=[time_TO_2'; mean(time_TO_2); std(time_TO_2)];
    time_TO_3=[time_TO_3'; mean(time_TO_3); std(time_TO_3)];
    time_IC_2=[time_IC_2'; mean(time_IC_2); std(time_IC_2)];
    time_IC_3=[time_IC_3'; mean(time_IC_3); std(time_IC_3)];
    
    eval(['PIF_' handles.name2(end) '=PIF_2;'])
    eval(['PIF_' handles.name3(end) '=PIF_3;'])
    eval(['LR_' handles.name2(end) '=LR_2;'])
    eval(['LR_' handles.name3(end) '=LR_3;'])
    eval(['IMP_' handles.name2(end) '=IMP_2;'])
    eval(['IMP_' handles.name3(end) '=IMP_3;'])
    eval(['Curve_' handles.name2(end) '=Curve2;'])
    eval(['Curve_' handles.name3(end) '=Curve3;'])
    eval(['time_TO_' handles.name2(end) '=time_TO_2;'])
    eval(['time_TO_' handles.name3(end) '=time_TO_3;'])
    eval(['time_IC_' handles.name2(end) '=time_IC_2;'])
    eval(['time_IC_' handles.name3(end) '=time_IC_3;'])
    
    
   % add symmetry tab 
    
     if LSI_method == 1
        LSI = 1;
        ASI = 0;
        NSI = 0;
     elseif LSI_method == 2
        LSI = 0;
        ASI = 1;
        NSI = 0;
     else
        LSI = 0;
        ASI = 0;
        NSI = 1;
     end

%     eval(['T = table(Step,PIF_' handles.name2(end) ',PIF_' handles.name3(end) ',PIF_Symm,'...
%         'LR_' handles.name2(end) ',LR_' handles.name3(end) ',LR_Symm,'...
%         'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',IMP_Symm)'])
  
    
    eval(['T = table(Step,PIF_' handles.name2(end) ',PIF_' handles.name3(end) ',PIF_Symm,'...
        'LR_' handles.name2(end) ',LR_' handles.name3(end) ',LR_Symm,'...
        'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',time_TO_' handles.name2(end) ',time_TO_' handles.name3(end) ',time_IC_' handles.name2(end) ',time_IC_' handles.name3(end) ')'])
    
%     eval(['T = table(Step,PIF_' handles.name2(end) ',PIF_' handles.name3(end)...
%         'LR_' handles.name2(end) ',LR_' handles.name3(end)...
%         'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',time_TO_' handles.name2(end) ',time_TO_' handles.name3(end) ',IMP_Symm)' ])
    
%     eval(['T = table(Step,PIF_' handles.name2(end) ',PIF_' handles.name3(end)...
%         'LR_' handles.name2(end) ',LR_' handles.name3(end)...
%         'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',IMP_Symm)' ])
    
    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));

    writetable(T,filename,'Sheet',trialname)

    eval(['T1 = table(Curve_' handles.name2(end) ',Curve_' handles.name3(end) ')']);
    
    answer1 = inputdlg({'input file name','input condition name'});
    filename1 = [selpath '/' char(answer1(1)) '.xlsx'];
    trialname1 = char(answer1(2));
    writetable(T1,filename1,'Sheet',trialname1);

%     eval(['TS = table(LSI' ',ASI'  ',NSI' ')'])
%     answerS = inputdlg({'input file name','input condition name'});
%     filenameS = [selpath '/' char(answerS(1)) '.xlsx'];
%     trialnameS = char(answerS(2));
% 
%     writetable(TS,filenameS,'Sheet',trialnameS)

elseif MovementType == 3 %Gait
    
    c2 = 1;
    c3 = 1;
    for t = 1:handles.nTrials
        
        if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)
            eval(['data = handles.T' num2str(t) '_data;']);
            eval(['IC_events_2 = handles.T' num2str(t) '_IC_2;']);
            eval(['IC_events_3 = handles.T' num2str(t) '_IC_3;']);
            eval(['TO_events_2 = handles.T' num2str(t) '_TO_2;']);
            eval(['TO_events_3 = handles.T' num2str(t) '_TO_3;']);
            
            for Step2 = 1:length(IC_events_2)
                if Step2 == length(IC_events_2)
                    ind = find(TO_events_2>IC_events_2(Step2));
                else
                    ind = find(TO_events_2>IC_events_2(Step2) & TO_events_2<IC_events_2(Step2+1));
                end
                
                if isempty(ind) %No toe of for this step, use time
                    
                    PIFLength =  floor(TimeCut_PIF*.001*handles.SamplingFreq);
                    IMPLength =  floor(TimeCut_IMP*.001*handles.SamplingFreq);
                    PPFLength =  floor((100-TimeCut_PIF)*.001*handles.SamplingFreq);
                    
                else %There is a TO, use %stance
                    
                    PIFLength =  floor((TO_events_2(Step2) - IC_events_2(Step2))*.01*StanceCut_PIF);
                    IMPLength =  floor((TO_events_2(Step2) - IC_events_2(Step2))*.01*StanceCut_IMP);
                    PPFLength =  floor((TO_events_2(Step2) - IC_events_2(Step2))- PIFLength);
                    
                end
                
                [peaks,locs] = findpeaks(data(IC_events_2(Step2):IC_events_2(Step2)+PIFLength,2));
                frame1 = IC_events_2(Step2)+PPFLength;
                frame2 = TO_events_2(Step2);
                [peaks2,locs2] = findpeaks(data((IC_events_2(Step2)+PPFLength):TO_events_2(Step2),2));
               
                
                if isempty(locs) %Then just take max
                    PIF_2(c2,1) = IC_events_2(Step2)+PIFLength
                    PIF_2(c2,2) = data(IC_events_2(Step2)+PIFLength,2);
                    PPF_2(c2,1) = TO_events_2(Step2)-PIFLength
                    PPF_2(c2,2) = data(TO_events_2(Step2)-PIFLength,2);
                    pklocs = PIFLength;
                else
                    [V,I] = sort(peaks,'descend');
                    PIF_2(c2,1) = IC_events_2(Step2)+locs(I(1));
                    PIF_2(c2,2) = V(1);
                    pklocs = locs(I(1));
                    [V_P,I_P] = sort(peaks2,'descend');
                    PPF_2(c2,1) = IC_events_2(Step2)+locs2(I_P(1));
                    PPF_2(c2,2) = V_P(1);
                    pklocs2 = locs2(I_P(1))+PIFLength;
                end
                
                IMP_2(c2) = trapz(data(IC_events_2(Step2):IC_events_2(Step2)+IMPLength,1),...
                    data(IC_events_2(Step2):IC_events_2(Step2)+IMPLength,2));
                
                % Inst Loading Rate
                ILR_2(c2) = max(diff(data(IC_events_2(Step2):IC_events_2(Step2)+pklocs,2)))/(1/handles.SamplingFreq);
                
                % Aver loading rate (20-80%)
                twentypercent = floor(.2*pklocs)+IC_events_2(Step2)-1;
                eightypercent = floor(.8*pklocs)+IC_events_2(Step2)-1;
                ALR_2(c2) = (data(eightypercent,2)-data(twentypercent,2))/((eightypercent-twentypercent)/handles.SamplingFreq);
                
                %GRF curve norm to 101 steps 
                normstep = ((((IC_events_2(Step2):1:TO_events_2(Step2))-IC_events_2(Step2))/(TO_events_2(Step2)-IC_events_2(Step2)))*100);
                normlength = [0:1:100];
                Curve2{:,c2} = data(IC_events_2(Step2):TO_events_2(Step2),2);  
                
                %Stance time 
                st_start_time =  data(IC_events_2(Step2),1);
                st_end_time = data(TO_events_2(Step2),1);
                st_2(c2,1) = st_end_time - st_start_time;
                
            
                c2=c2+1;
                
                
            end
            
             gg2 = size(Curve2);
            for cc2 = 1:gg2(2)
                sizes2(:,cc2) = length(Curve2{cc2});
            end
            [curve2_leng, curve2_num] = max(sizes2);

            Curve2new = NaN(curve2_leng,gg2(2));
            for cc2 = 1:gg2(2)
                if length(Curve2{cc2}) < curve2_leng
                    Curve2new(:,cc2) = [Curve2{cc2};NaN(curve2_leng-length(Curve2{cc2}),1)];
                else
                    Curve2new(:,cc2) = Curve2{cc2};
                end
            end
            
            for Step3 = 1:length(IC_events_3)
                if Step3 == length(IC_events_3)
                    ind = find(TO_events_3>IC_events_3(Step3));
                else
                    ind = find(TO_events_3>IC_events_3(Step3) & TO_events_3<IC_events_3(Step3+1));
                end
                
                if isempty(ind) %No toe of for this step, use time
                    
                    PIFLength =  floor(TimeCut_PIF*.001*handles.SamplingFreq);
                    IMPLength =  floor(TimeCut_IMP*.001*handles.SamplingFreq);
                    
                    
                else %There is a TO, use %stance
                    
                    PIFLength =  floor((TO_events_3(Step3) - IC_events_3(Step3))*.01*StanceCut_PIF);
                    IMPLength =  floor((TO_events_3(Step3) - IC_events_3(Step3))*.01*StanceCut_IMP);
                    PPFLength =  floor((TO_events_3(Step3) - IC_events_3(Step3))- PIFLength);
                end
                
                [peaks,locs] = findpeaks(data(IC_events_3(Step3):IC_events_3(Step3)+PIFLength ,3));
                [peaks2,locs2] = findpeaks(data((IC_events_3(Step3)+PIFLength):TO_events_3(Step3),3));
                
                if isempty(locs) %Then just take max
                    PIF_3(c3,1) = IC_events_3(Step3)+PIFLength;
                    PIF_3(c3,2) = data(IC_events_3(Step3)+PIFLength,3);
                    pklocs = PIFLength;
                    PPF_3(c3,1) = IC_events_3(Step3)+PIFLength;
                    PPF_3(c3,2) = data(IC_events_3(Step3)+PIFLength,3);
                  
                else
                    [V,I] = sort(peaks,'descend');
                    PIF_3(c3,1) = IC_events_3(Step3)+locs(I(1));
                    PIF_3(c3,2) = V(1);
                    pklocs = locs(I(1));
                    [V_P,I_P] = sort(peaks2,'descend');
                    PPF_3(c3,1) = IC_events_3(Step3)+locs2(I_P(1));
                    PPF_3(c3,2) = V_P(1);
                    pklocs2 = locs2(I_P(1))+PIFLength;
                end
                
                IMP_3(c3) = trapz(data(IC_events_3(Step3):IC_events_3(Step3)+IMPLength,1),...
                    data(IC_events_3(Step3):IC_events_3(Step3)+IMPLength,3));
                
                % Inst Loading Rate
                ILR_3(c3) = max(diff(data(IC_events_3(Step3):IC_events_3(Step3)+pklocs,3)))/(1/handles.SamplingFreq);
                
                % Aver loading rate (20-80%)
                twentypercent = floor(.2*pklocs)+IC_events_3(Step3)-1;
                eightypercent = floor(.8*pklocs)+IC_events_3(Step3)-1;
                
                ALR_3(c3) = (data(eightypercent,3)-data(twentypercent,3))/((eightypercent-twentypercent)/handles.SamplingFreq);
                
                normstep = ((((IC_events_3(Step3):1:TO_events_3(Step3))-IC_events_3(Step3))/(TO_events_3(Step3)-IC_events_3(Step3)))*100);
                normlength = [0:1:100];
                
                Curve3{:,c3} = data(IC_events_3(Step3):TO_events_3(Step3),3); 
                
                %Stance time 
                st_start_time =  data(IC_events_3(Step3),1);
                st_end_time = data(TO_events_3(Step3),1);
                st_3(c3,1) =  st_end_time - st_start_time;
                
                
                c3=c3+1;
            end
            

            gg3 = size(Curve3);
            for cc3 = 1:gg3(2)
                sizes3(:,cc3) = length(Curve3{cc3});
            end
            [curve3_leng, curve3_num] = max(sizes3);

            Curve3new = NaN(curve3_leng,gg3(2));
            for cc3 = 1:gg3(2)
                if length(Curve3{cc3}) < curve3_leng
                    Curve3new(:,cc3) = [Curve3{cc3};NaN(curve3_leng-length(Curve3{cc3}),1)];
                else
                    Curve3new(:,cc3) = Curve3{cc3};
                end
            end

            if length(Curve2new) > length(Curve3new)
                Curve3new(length(Curve3new)+1:length(Curve2new),1:5) = NaN; 
            elseif length(Curve3new) > length(Curve2new)
                Curve2new(length(Curve2new)+1:length(Curve3new),1:5) = NaN;
            end

            eval(['handles.plotpeaks_3_' num2str(t)  ' = [PIF_3];'])
            eval(['handles.plotpeaks_2_' num2str(t)  ' = [PIF_2];'])

        end
    end

    
    PIF_2 = PIF_2(:,2)';
    PIF_3 = PIF_3(:,2)';
    PPF_2 = PPF_2(:,2)';
    PPF_3 = PPF_3(:,2)';
    
    if length(PIF_2)>length(PIF_3)
        trim2 = length(PIF_2)-length(PIF_3);
        trim3 = 0;
    else
        trim2 = length(PIF_2)-length(PIF_3);
        trim3 = 0;
    end
    if length(PPF_2)>length(PPF_3)
        trim2 = length(PPF_2)-length(PPF_3);
        trim3 = 0;
    else
        trim2 = length(PPF_2)-length(PPF_3);
        trim3 = 0;
    end
    
    
    LSI_method = get(handles.LSI_Method,'Value');
    if LSI_method == 1%LSI
        
        if get(handles.DominantLimb,'Value') == 1
            PIF_Symm = (PIF_3(1:end-trim3)./PIF_2(1:end-trim2))*100; % R/L*100 %Surg/Nonsurg
            PPF_Symm = (PPF_3(1:end-trim3)./PPF_2(1:end-trim2))*100;
            ILR_Symm = (ILR_3(1:end-trim3)./ILR_2(1:end-trim2))*100;
            ALR_Symm = (ALR_3(1:end-trim3)./ALR_2(1:end-trim2))*100;
            IMP_Symm = (IMP_3(1:end-trim3)./IMP_2(1:end-trim2))*100;
        else
            PIF_Symm = (PIF_2(1:end-trim2)./PIF_3(1:end-trim3))*100; %L/R*100 %
            PPF_Symm = (PPF_2(1:end-trim3)./PPF_3(1:end-trim2))*100;
            ILR_Symm = (ILR_2(1:end-trim2)./ILR_3(1:end-trim3))*100;
            ALR_Symm = (ALR_2(1:end-trim2)./ALR_3(1:end-trim3))*100;
            IMP_Symm = (IMP_2(1:end-trim2)./IMP_3(1:end-trim3))*100;
        end
        
    elseif LSI_method == 2 %ASI
        
        PIF_Symm = (abs(PIF_2(1:end-trim2) - PIF_3(1:end-trim3))./((PIF_2(1:end-trim2) + PIF_3(1:end-trim3))/2))*100;
        PPF_Symm = (abs(PPF_2(1:end-trim2) - PPF_3(1:end-trim3))./((PPF_2(1:end-trim2) + PPF_3(1:end-trim3))/2))*100;
        ILR_Symm = (abs(ILR_2(1:end-trim2) - ILR_3(1:end-trim3))./((ILR_2(1:end-trim2) + ILR_3(1:end-trim3))/2))*100;
        ALR_Symm = (abs(ALR_2(1:end-trim2) - ALR_3(1:end-trim3))./((ALR_2(1:end-trim2) + ALR_3(1:end-trim3))/2))*100;
        IMP_Symm = (abs(IMP_2(1:end-trim2) - IMP_3(1:end-trim3))./((IMP_2(1:end-trim2) + IMP_3(1:end-trim3))/2))*100;
        
        
    else %NSI
        
        for j=1:min([length(PIF_2),length(PIF_3)])
            if get(handles.DominantLimb,'Value') == 1 %2 is dom
                PIF_Symm(j) = 100*(PIF_2(j) - PIF_3(j))/max([PIF_2,PIF_3]);
                PPF_Symm(j) = 100*(PPF_2(j) - PPF_3(j))/max([PPF_2,PPF_3]);
                ILR_Symm(j) = 100*(ILR_2(j) - ILR_3(j))/max([ILR_2,ILR_3]);
                ALR_Symm(j) = 100*(ALR_2(j) - ALR_3(j))/max([ALR_2,ALR_3]);
                IMP_Symm(j) = 100*(IMP_2(j) - IMP_3(j))/max([IMP_2,IMP_3]);
            else
                PIF_Symm(j) = 100*(PIF_3(j) - PIF_2(j))/max([PIF_2,PIF_3]);
                PPF_Symm(j) = 100*(PPF_3(j) - PPF_2(j))/max([PPF_2,PPF_3]);
                ILR_Symm(j) = 100*(ILR_3(j) - ILR_2(j))/max([ILR_2,ILR_3]);
                ALR_Symm(j) = 100*(ALR_3(j) - ALR_2(j))/max([ALR_2,ALR_3]);
                IMP_Symm(j) = 100*(IMP_3(j) - IMP_2(j))/max([IMP_2,IMP_3]);
            end
        end
        
    end
    
    Step = [1:min([length(PIF_2),length(PIF_3)])]';
    Step = num2cell(Step);
    Step{end+1} = 'Mean';
    Step{end+1} = 'StDev';
    
    PIF_2 = [PIF_2(1:end-trim2)';mean(PIF_2(1:end-trim2));std(PIF_2(1:end-trim2))];
    PIF_3 = [PIF_3(1:end-trim3)';mean(PIF_3(1:end-trim3));std(PIF_3(1:end-trim3))];
    PIF_Symm = [PIF_Symm';mean(PIF_Symm);std(PIF_Symm)];
    
    PPF_2 = [PPF_2(1:end-trim2)';mean(PPF_2(1:end-trim2));std(PPF_2(1:end-trim2))];
    PPF_3 = [PPF_3(1:end-trim3)';mean(PPF_3(1:end-trim3));std(PPF_3(1:end-trim3))];
    PPF_Symm = [PPF_Symm';mean(PPF_Symm);std(PPF_Symm)];

    ILR_2 = [ILR_2(1:end-trim2)';mean(ILR_2(1:end-trim2));std(ILR_2(1:end-trim2))];
    ILR_3 = [ILR_3(1:end-trim3)';mean(ILR_3(1:end-trim3));std(ILR_3(1:end-trim3))];
    ILR_Symm = [ILR_Symm';mean(ILR_Symm);std(ILR_Symm)];
    
    ALR_2 = [ALR_2(1:end-trim2)';mean(ALR_2(1:end-trim2));std(ALR_2(1:end-trim2))];
    ALR_3 = [ALR_3(1:end-trim3)';mean(ALR_3(1:end-trim3));std(ALR_3(1:end-trim3))];
    ALR_Symm = [ALR_Symm';mean(ALR_Symm);std(ALR_Symm)];
    
    IMP_2 = [IMP_2(1:end-trim2)';mean(IMP_2(1:end-trim2));std(IMP_2(1:end-trim2))];
    IMP_3 = [IMP_3(1:end-trim3)';mean(IMP_3(1:end-trim3));std(IMP_3(1:end-trim3))];
    IMP_Symm = [IMP_Symm';mean(IMP_Symm);std(IMP_Symm)];
    
    
    eval(['PIF_' handles.name2(end) '=PIF_2;'])
    eval(['PIF_' handles.name3(end) '=PIF_3;'])
    eval(['PPF_' handles.name2(end) '=PPF_2;'])
    eval(['PPF_' handles.name3(end) '=PPF_3;'])
    eval(['ILR_' handles.name2(end) '=ILR_2;'])
    eval(['ILR_' handles.name3(end) '=ILR_3;'])
    eval(['ALR_' handles.name2(end) '=ALR_2;'])
    eval(['ALR_' handles.name3(end) '=ALR_3;'])
    eval(['IMP_' handles.name2(end) '=IMP_2;'])
    eval(['IMP_' handles.name3(end) '=IMP_3;'])
    eval(['st_' handles.name2(end) '=st_2;'])
    eval(['st_' handles.name3(end) '=st_3;'])
    eval(['Curve_' handles.name2(end) '=Curve2new;'])
    eval(['Curve_' handles.name3(end) '=Curve3new;'])
    
        % add symmetry tab 
    
     if LSI_method == 1
        LSI = 1;
        ASI = 0;
        NSI = 0;
     elseif LSI_method == 2
        LSI = 0;
        ASI = 1;
        NSI = 0;
     else
        LSI = 0;
        ASI = 0;
        NSI = 1;
     end

    eval(['T = table(Step,PIF_' handles.name2(end) ',PIF_' handles.name3(end) ...
        ',PPF_' handles.name2(end) ',PPF_' handles.name3(end) ',PIF_LSI,'...
        'ILR_' handles.name2(end) ',ILR_' handles.name3(end) ',ILR_LSI,'...
        'ALR_' handles.name2(end) ',ALR_' handles.name3(end) ',ALR_LSI,'...
        'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',IMP_LSI)']);
    
   selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));
    
    writetable(T,filename,'Sheet',trialname);
   
    
    eval(['T1 = table(Curve_' handles.name2(end) ',Curve_' handles.name3(end) ')']);
    
    answer1 = inputdlg({'input file name','input condition name'});
    filename1 = [selpath '/' char(answer1(1)) '.xlsx'];
    trialname1 = char(answer1(2));
    writetable(T1,filename1,'Sheet',trialname1);

    eval(['T2 = table(st_' handles.name2(end) ',st_' handles.name3(end) ')']);
    answer2 = inputdlg({'input file name','input condition name'});
    filename2 = [selpath '/' char(answer2(1)) '.xlsx'];
    trialname2 = char(answer2(2));
    writetable(T2,filename2,'Sheet',trialname2);

    eval(['TS = table(LSI' ',ASI'  ',NSI' ')'])
    answerS = inputdlg({'input file name','input condition name'});
    filenameS = [selpath '/' char(answerS(1)) '.xlsx'];
    trialnameS = char(answerS(2));

    writetable(TS,filenameS,'Sheet',trialnameS)
    
elseif MovementType == 4 %Squatting
    
    % Compute peak force and impulse between peak forces
    for t = 1:handles.nTrials

    if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)
        eval(['data = handles.T' num2str(t) '_data;']);
        
        [pks2,lks2] = findpeaks(data(:,2),'MinPeakHeight',0.75*max(data(:,2)),'MinPeakDistance',40);
        [pks3,lks3] = findpeaks(data(:,3),'MinPeakHeight',0.75*max(data(:,3)),'MinPeakDistance',40);

        plot(data(lks2,1),pks2,'k o',data(lks3,1),pks3,'ko','MarkerSize',15)
        
        set(handles.UserMessageBox,'String',['Click on any missidentified peaks for ' handles.name2 ' or click after the end of the trial if no misidentified peaks'])
        clicks = ginput(10);
        

        if clicks(1,1)<data(end,1)
        clicks=clicks*handles.SamplingFreq;
        for i=1:size(clicks,1)
        [~,ind] = nanmin(abs(lks2-clicks(i,1)));
        lks2(ind) = nan;  
        pks2(ind) = nan;
        end
        pks2(isnan(pks2)) = [];
        lks2(isnan(lks2)) = [];
        end
        
        set(handles.UserMessageBox,'String',['Click on any missidentified peaks for ' handles.name3 ' or click after the end of the trial if no misidentified peaks'])
        clicks = ginput(10);
        

        if clicks(1,1)<data(end,1)
        clicks=clicks*handles.SamplingFreq;
        for i=1:size(clicks,1)
        [~,ind] = nanmin(abs(lks3-clicks(i,1)));
        lks3(ind) = nan;  
        pks3(ind) = nan;
        end
        pks3(isnan(pks3)) = [];
        lks3(isnan(lks3)) = [];
        end
        
        plot_data(hObject, eventdata, handles)
        plot(data(lks2,1),pks2,'k o',data(lks3,1),pks3,'ko','MarkerSize',15)
        
        %%%%%%%%%%%%%%%%%%%%%%
        
        for i = 1:length(lks2)-1
           PF_2(i) = pks2(i+1);
           PF_3(i) = pks3(i+1);
           IMP_2(i) = trapz(data(lks2(i):lks2(i+1),1),data(lks2(i):lks2(i+1),2));
           IMP_3(i) = trapz(data(lks3(i):lks3(i+1),1),data(lks3(i):lks3(i+1),3));
        end
        
            eval(['handles.plotpeaks_3_' num2str(t)  ' = [lks2,pks2];'])
            eval(['handles.plotpeaks_2_' num2str(t)  ' = [lks3,pks3];'])
            
            
    LSI_method = get(handles.LSI_Method,'Value');
    if LSI_method == 1%LSI
        
        if get(handles.DominantLimb,'Value') == 1
        PF_LSI = (PF_3./PF_2)*100;
        IMP_Symm = (IMP_3./IMP_2)*100;
        else
        PF_LSI = (PF_2./PF_3)*100;
        IMP_Symm = (IMP_2./IMP_3)*100;
        end
        
    elseif LSI_method == 2 %ASI
        
        PF_LSI = (abs(PF_2 - PF_3)./((PF_2 + PF_3)/2))*100;
        IMP_Symm = (abs(IMP_2 - IMP_3)./((IMP_2 + IMP_3)/2))*100;
        
    else %NSI
        
        for j=1:length(IMP_2)
            if get(handles.DominantLimb,'Value') == 1 %2 is dom
                PF_LSI(j) = 100*(PF_2(j) - PF_3(j))/max([PF_2,PF_3]);
                IMP_Symm(j) = 100*(IMP_2(j) - IMP_3(j))/max([IMP_2,IMP_3]);
            else
                PF_LSI(j) = 100*(PF_3(j) - PF_2(j))/max([PF_2,PF_3]);
                IMP_Symm(j) = 100*(IMP_3(j) - IMP_2(j))/max([IMP_2,IMP_3]);
            end
        end

    end
    
    Step = [1:min([length(PF_2),length(PF_3)])]';
    Step = num2cell(Step);
    Step{end+1} = 'Mean';
    Step{end+1} = 'StDev';

    PF_2 = [PF_2';mean(PF_2);std(PF_2)];
    PF_3 = [PF_3';mean(PF_3);std(PF_3)];
    PF_LSI = [PF_LSI';mean(PF_LSI);std(PF_LSI)];
    
    IMP_2 = [IMP_2';mean(IMP_2);std(IMP_2)];
    IMP_3 = [IMP_3';mean(IMP_3);std(IMP_3)];
    IMP_Symm = [IMP_Symm';mean(IMP_Symm);std(IMP_Symm)];
    
    eval(['PF_' handles.name2(end) '=PF_2;'])
    eval(['PF_' handles.name3(end) '=PF_3;'])
    eval(['IMP_' handles.name2(end) '=IMP_2;'])
    eval(['IMP_' handles.name3(end) '=IMP_3;'])
    


    eval(['T = table(Step,PF_' handles.name2(end) ',PF_' handles.name3(end) ',PF_LSI,'...
        'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',IMP_LSI)'])

    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));

    writetable(T,filename,'Sheet',trialname)
    
    
        %%%%%%%%%%%%%%%%%%%%%%

    end
    end
    elseif MovementType == 4 %Squatting
    
    % Compute peak force and impulse between peak forces
    for t = 1:handles.nTrials

    if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)
        eval(['data = handles.T' num2str(t) '_data;']);
        
        [pks2,lks2] = findpeaks(data(:,2),'MinPeakHeight',0.75*max(data(:,2)),'MinPeakDistance',40);
        [pks3,lks3] = findpeaks(data(:,3),'MinPeakHeight',0.75*max(data(:,3)),'MinPeakDistance',40);

        plot(data(lks2,1),pks2,'k o',data(lks3,1),pks3,'ko','MarkerSize',15)
        
        set(handles.UserMessageBox,'String',['Click on any missidentified peaks for ' handles.name2 ' or click after the end of the trial if no misidentified peaks'])
        clicks = ginput(10);
        

        if clicks(1,1)<data(end,1)
        clicks=clicks*handles.SamplingFreq;
        for i=1:size(clicks,1)
        [~,ind] = nanmin(abs(lks2-clicks(i,1)));
        lks2(ind) = nan;  
        pks2(ind) = nan;
        end
        pks2(isnan(pks2)) = [];
        lks2(isnan(lks2)) = [];
        end
        
        set(handles.UserMessageBox,'String',['Click on any missidentified peaks for ' handles.name3 ' or click after the end of the trial if no misidentified peaks'])
        clicks = ginput(10);
        

        if clicks(1,1)<data(end,1)
        clicks=clicks*handles.SamplingFreq;
        for i=1:size(clicks,1)
        [~,ind] = nanmin(abs(lks3-clicks(i,1)));
        lks3(ind) = nan;  
        pks3(ind) = nan;
        end
        pks3(isnan(pks3)) = [];
        lks3(isnan(lks3)) = [];
        end
        
        plot_data(hObject, eventdata, handles)
        plot(data(lks2,1),pks2,'k o',data(lks3,1),pks3,'ko','MarkerSize',15)
        
        %%%%%%%%%%%%%%%%%%%%%%
        
        for i = 1:length(lks2)-1
           PF_2(i) = pks2(i+1);
           PF_3(i) = pks3(i+1);
           IMP_2(i) = trapz(data(lks2(i):lks2(i+1),1),data(lks2(i):lks2(i+1),2));
           IMP_3(i) = trapz(data(lks3(i):lks3(i+1),1),data(lks3(i):lks3(i+1),3));
        end
        
            eval(['handles.plotpeaks_3_' num2str(t)  ' = [lks2,pks2];'])
            eval(['handles.plotpeaks_2_' num2str(t)  ' = [lks3,pks3];'])
            
            
    LSI_method = get(handles.LSI_Method,'Value');
    if LSI_method == 1%LSI
        
        if get(handles.DominantLimb,'Value') == 1
        PF_LSI = (PF_3./PF_2)*100;
        IMP_Symm = (IMP_3./IMP_2)*100;
        else
        PF_LSI = (PF_2./PF_3)*100;
        IMP_Symm = (IMP_2./IMP_3)*100;
        end
        
    elseif LSI_method == 2 %ASI
        
        PF_LSI = (abs(PF_2 - PF_3)./((PF_2 + PF_3)/2))*100;
        IMP_Symm = (abs(IMP_2 - IMP_3)./((IMP_2 + IMP_3)/2))*100;
        
    else %NSI
        
        for j=1:length(IMP_2)
            if get(handles.DominantLimb,'Value') == 1 %2 is dom
                PF_LSI(j) = 100*(PF_2(j) - PF_3(j))/max([PF_2,PF_3]);
                IMP_Symm(j) = 100*(IMP_2(j) - IMP_3(j))/max([IMP_2,IMP_3]);
            else
                PF_LSI(j) = 100*(PF_3(j) - PF_2(j))/max([PF_2,PF_3]);
                IMP_Symm(j) = 100*(IMP_3(j) - IMP_2(j))/max([IMP_2,IMP_3]);
            end
        end

    end
    
    Step = [1:min([length(PF_2),length(PF_3)])]';
    Step = num2cell(Step);
    Step{end+1} = 'Mean';
    Step{end+1} = 'StDev';

    PF_2 = [PF_2';mean(PF_2);std(PF_2)];
    PF_3 = [PF_3';mean(PF_3);std(PF_3)];
    PF_LSI = [PF_LSI';mean(PF_LSI);std(PF_LSI)];
    
    IMP_2 = [IMP_2';mean(IMP_2);std(IMP_2)];
    IMP_3 = [IMP_3';mean(IMP_3);std(IMP_3)];
    IMP_Symm = [IMP_Symm';mean(IMP_Symm);std(IMP_Symm)];
    
    eval(['PF_' handles.name2(end) '=PF_2;'])
    eval(['PF_' handles.name3(end) '=PF_3;'])
    eval(['IMP_' handles.name2(end) '=IMP_2;'])
    eval(['IMP_' handles.name3(end) '=IMP_3;'])
    


    eval(['T = table(Step,PF_' handles.name2(end) ',PF_' handles.name3(end) ',PF_LSI,'...
        'IMP_' handles.name2(end) ',IMP_' handles.name3(end) ',IMP_LSI)'])

    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));

    writetable(T,filename,'Sheet',trialname)
    
    
        %%%%%%%%%%%%%%%%%%%%%%

    end
    end
    
elseif MovementType == 5 %Sit to Stand
    
    c2 = 1;
    c3 = 1;
    for t = 1:handles.nTrials
        
        if isempty(handles.closedtrials) || ~any(handles.closedtrials == t)
            eval(['data = handles.T' num2str(t) '_data;']);
            eval(['IC_events_2 = handles.T' num2str(t) '_IC_2;'])
            eval(['IC_events_3 = handles.T' num2str(t) '_IC_3;'])
            eval(['TO_events_2 = handles.T' num2str(t) '_TO_2;'])
            eval(['TO_events_3 = handles.T' num2str(t) '_TO_3;'])
            
            % _2 is Left
            for Step2 = 1:length(IC_events_2)
                if Step2 == length(IC_events_2)
                    ind = find(TO_events_2>IC_events_2(Step2));
                else
                    ind = find(TO_events_2>IC_events_2(Step2) & TO_events_2<IC_events_2(Step2+1));
                end
                
                    %initial peak
                    PIFLength =  floor((TO_events_2(Step2) - IC_events_2(Step2))*.01*StanceCut_PIF); %settings should be set to 15%
                    %total peak
                    PPFLength =  floor((TO_events_2(Step2) - IC_events_2(Step2))); %total peak

                %peaks for initial peak
                [peaks,locs] = findpeaks(data(IC_events_2(Step2):IC_events_2(Step2)+PIFLength,2));
                frame1 = IC_events_2(Step2)
                frame2 = TO_events_2(Step2)
                %peaks for total peak
                [peaks2,locs2] = findpeaks(data((IC_events_2(Step2):TO_events_2(Step2)),2));
                
                if isempty(locs) %Then just take max
                    PIF_2(c2,1) = IC_events_2(Step2)+PIFLength
                    PIF_2(c2,2) = data(IC_events_2(Step2)+PIFLength,2);
                    PPF_2(c2,1) = TO_events_2(Step2)-PIFLength
                    PPF_2(c2,2) = data(IC_events_2(Step2):TO_events_2(Step2),2);
                    pklocs = PIFLength;
                else
                    [V,I] = sort(peaks,'descend');
                    %initial peak
                    PIF_2(c2,1) = IC_events_2(Step2)+locs(I(1));
                    PIF_2(c2,2) = V(1);
                    pklocs = locs(I(1)); %index of the peak relative to the IC event set to 1
                    %total peak
                    [V_P,I_P] = sort(peaks2,'descend');
                    PPF_2(c2,1) = IC_events_2(Step2)+locs2(I_P(1));
                    PPF_2(c2,2) = V_P(1);
                    pklocs2 = locs2(I_P(1));
                end

                
                Curve2{:,c2} = data(IC_events_2(Step2):TO_events_2(Step2),2);  
                
                %Stance time 
                st_start_time =  data(IC_events_2(Step2),1);
                st_end_time = data(TO_events_2(Step2),1);
                st_2(c2,1) = st_end_time - st_start_time;

                %Time to initial peak
                st_start_time =  data(IC_events_2(Step2),1);
                st_IP_time = data(PIF_2(c2,1),1);
                it_2(c2,1) = st_IP_time - st_start_time;

                %Time to total peak
                st_start_time =  data(IC_events_2(Step2),1);
                st_TP_time = data(PPF_2(c2,1),1);
                pt_2(c2,1) = st_TP_time - st_start_time;
                
                c2=c2+1;
            end
            
            gg2 = size(Curve2);
            for cc2 = 1:gg2(2)
                sizes2(:,cc2) = length(Curve2{cc2});
            end
            [curve2_leng, curve2_num] = max(sizes2);

            Curve2new = NaN(curve2_leng,gg2(2));
            for cc2 = 1:gg2(2)
                if length(Curve2{cc2}) < curve2_leng
                    Curve2new(:,cc2) = [Curve2{cc2};NaN(curve2_leng-length(Curve2{cc2}),1)];
                else
                    Curve2new(:,cc2) = Curve2{cc2};
                end
            end
            
            % _3 is Right
            for Step3 = 1:length(IC_events_3)
                if Step3 == length(IC_events_3)
                    ind = find(TO_events_3>IC_events_3(Step3));
                else
                    ind = find(TO_events_3>IC_events_3(Step3) & TO_events_3<IC_events_3(Step3+1));
                end
                
                 %initial peak
                    PIFLength =  floor((TO_events_3(Step3) - IC_events_3(Step3))*.01*StanceCut_PIF); %settings should be set to 15%
                    %total peak
                    PPFLength =  floor((TO_events_3(Step2) - IC_events_3(Step3))); %total peak
                
                
                %peaks for initial peak
                [peaks,locs] = findpeaks(data(IC_events_3(Step3):IC_events_3(Step3)+PIFLength,3));
                %peaks for total peak
                [peaks2,locs2] = findpeaks(data((IC_events_3(Step3):TO_events_3(Step3)),3));
                
                if isempty(locs) %Then just take max
                    PIF_3(c3,1) = IC_events_3(Step3)+PIFLength;
                    PIF_3(c3,2) = data(IC_events_3(Step3)+PIFLength,3);
                    pklocs = PIFLength;
                    PPF_3(c3,1) = IC_events_3(Step3)+PIFLength;
                    PPF_3(c3,2) = data(IC_events_3(Step3)+PIFLength,3);
                  
                else
                    [V,I] = sort(peaks,'descend');
                    PIF_3(c3,1) = IC_events_3(Step3)+locs(I(1));
                    PIF_3(c3,2) = V(1);
                    pklocs = locs(I(1));
                    [V_P,I_P] = sort(peaks2,'descend');
                    PPF_3(c3,1) = IC_events_3(Step3)+locs2(I_P(1));
                    PPF_3(c3,2) = V_P(1);
                    pklocs2 = locs2(I_P(1));
                end

                Curve3{:,c3} = data(IC_events_3(Step3):TO_events_3(Step3),3); 
                
                %Stance time 
                st_start_time =  data(IC_events_3(Step3),1);
                st_end_time = data(TO_events_3(Step3),1);
                st_3(c3,1) =  st_end_time - st_start_time;

                %Time to initial peak
                st_start_time =  data(IC_events_3(Step3),1);
                st_IP_time = data(PIF_3(c3,1),1);
                it_3(c3,1) = st_IP_time - st_start_time;

                %Time to total peak
                st_start_time =  data(IC_events_3(Step3),1);
                st_TP_time = data(PPF_3(c3,1),1);
                pt_3(c3,1) = st_TP_time - st_start_time;
                
                c3=c3+1;
            end
            
            %added to write different length curves
            gg3 = size(Curve3);
            for cc3 = 1:gg3(2)
                sizes3(:,cc3) = length(Curve3{cc3});
            end
            [curve3_leng, curve3_num] = max(sizes3);

            Curve3new = NaN(curve3_leng,gg3(2));
            for cc3 = 1:gg3(2)
                if length(Curve3{cc3}) < curve3_leng
                    Curve3new(:,cc3) = [Curve3{cc3};NaN(curve3_leng-length(Curve3{cc3}),1)];
                else
                    Curve3new(:,cc3) = Curve3{cc3};
                end
            end

            if length(Curve2new) > length(Curve3new)
                Curve3new(length(Curve3new)+1:length(Curve2new),1:5) = NaN; 
            elseif length(Curve3new) > length(Curve2new)
                Curve2new(length(Curve2new)+1:length(Curve3new),1:5) = NaN;
            end

            

            eval(['handles.plotpeaks_3_' num2str(t)  ' = [PIF_3];'])
            eval(['handles.plotpeaks_2_' num2str(t)  ' = [PIF_2];'])

        end
    end
    
    %inital peak
    PIF_2 = PIF_2(:,2)';
    PIF_3 = PIF_3(:,2)';
    %total peak
    PTF_2 = PPF_2(:,2)';
    PTF_3 = PPF_3(:,2)';
    
    if length(PIF_2)>length(PIF_3)
        trim2 = length(PIF_2)-length(PIF_3);
        trim3 = 0;
    else
        trim2 = length(PIF_2)-length(PIF_3);
        trim3 = 0;
    end
    if length(PTF_2)>length(PTF_3)
        trim2 = length(PTF_2)-length(PTF_3);
        trim3 = 0;
    else
        trim2 = length(PTF_2)-length(PTF_3);
        trim3 = 0;
    end
    
    
    LSI_method = get(handles.LSI_Method,'Value');
    if LSI_method == 1%LSI
        
        if get(handles.DominantLimb,'Value') == 1
            PIF_Symm = (PIF_3(1:end-trim3)./PIF_2(1:end-trim2))*100; % R/L*100 %Surg/Nonsurg
            PTF_Symm = (PTF_3(1:end-trim3)./PTF_2(1:end-trim2))*100;
        else
            PIF_Symm = (PIF_2(1:end-trim2)./PIF_3(1:end-trim3))*100; %L/R*100 %
            PTF_Symm = (PTF_2(1:end-trim3)./PTF_3(1:end-trim2))*100;
        end
        
    elseif LSI_method == 2 %ASI
        
        PIF_Symm = (abs(PIF_2(1:end-trim2) - PIF_3(1:end-trim3))./((PIF_2(1:end-trim2) + PIF_3(1:end-trim3))/2))*100;
        PTF_Symm = (abs(PTF_2(1:end-trim2) - PTF_3(1:end-trim3))./((PTF_2(1:end-trim2) + PTF_3(1:end-trim3))/2))*100;

    else %NSI
        
        for j=1:min([length(PIF_2),length(PIF_3)])
            if get(handles.DominantLimb,'Value') == 1 %2 is dom
                PIF_Symm(j) = 100*(PIF_2(j) - PIF_3(j))/max([PIF_2,PIF_3]);
                PTF_Symm(j) = 100*(PTF_2(j) - PTF_3(j))/max([PTF_2,PTF_3]);
            else
                PIF_Symm(j) = 100*(PIF_3(j) - PIF_2(j))/max([PIF_2,PIF_3]);
                PTF_Symm(j) = 100*(PTF_3(j) - PTF_2(j))/max([PTF_2,PTF_3]);
            end
        end
        
    end
    
    Step = [1:min([length(PIF_2),length(PIF_3)])]';
    Step = num2cell(Step);
    Step{end+1} = 'Mean';
    Step{end+1} = 'StDev';
    
    PIF_2 = [PIF_2(1:end-trim2)';mean(PIF_2(1:end-trim2));std(PIF_2(1:end-trim2))];
    PIF_3 = [PIF_3(1:end-trim3)';mean(PIF_3(1:end-trim3));std(PIF_3(1:end-trim3))];
    PIF_Symm = [PIF_Symm';mean(PIF_Symm);std(PIF_Symm)];
    
    PTF_2 = [PTF_2(1:end-trim2)';mean(PTF_2(1:end-trim2));std(PTF_2(1:end-trim2))];
    PTF_3 = [PTF_3(1:end-trim3)';mean(PTF_3(1:end-trim3));std(PTF_3(1:end-trim3))];
    PTF_Symm = [PTF_Symm';mean(PTF_Symm);std(PTF_Symm)];

    
    eval(['PIF_' handles.name2(end) '=PIF_2;'])
    eval(['PIF_' handles.name3(end) '=PIF_3;'])
    eval(['PTF_' handles.name2(end) '=PTF_2;'])
    eval(['PTF_' handles.name3(end) '=PTF_3;'])
    eval(['st_' handles.name2(end) '=st_2;'])
    eval(['st_' handles.name3(end) '=st_3;'])
    eval(['it_' handles.name2(end) '=it_2;'])
    eval(['it_' handles.name3(end) '=it_3;'])
    eval(['pt_' handles.name2(end) '=pt_2;'])
    eval(['pt_' handles.name3(end) '=pt_3;'])
    eval(['Curve_' handles.name2(end) '=Curve2new;'])
    eval(['Curve_' handles.name3(end) '=Curve3new;'])
    
    % add symmetry tab 
    
     if LSI_method == 1
        LSI = 1;
        ASI = 0;
        NSI = 0;
     elseif LSI_method == 2
        LSI = 0;
        ASI = 1;
        NSI = 0;
     else
        LSI = 0;
        ASI = 0;
        NSI = 1;
     end

    %writing symmetry measures out
    eval(['T = table(Step,PIF_' handles.name2(end) ',PIF_' handles.name3(end) ',PIF_Symm,' ...
        'PTF_' handles.name2(end) ',PTF_' handles.name3(end) ',PTF_Symm)'])
    
    selpath = uigetdir(handles.fpath);
    answer = inputdlg({'input file name','input condition name'});
    filename = [selpath '/' char(answer(1)) '.xlsx'];
    trialname = char(answer(2));
    
    writetable(T,filename,'Sheet',trialname)
    
    %writing table of GRF curves
    eval(['T1 = table(Curve_' handles.name2(end) ',Curve_' handles.name3(end) ')'])
    
    answer1 = inputdlg({'input file name','input condition name'});
    filename1 = [selpath '/' char(answer1(1)) '.xlsx'];
    trialname1 = char(answer1(2));
    writetable(T1,filename1,'Sheet',trialname1)

    %writing table of t
    eval(['T2 = table(st_' handles.name2(end) ',st_' handles.name3(end)  ',it_' handles.name2(end)...
        ',it_' handles.name3(end) ',pt_' handles.name2(end) ',pt_' handles.name3(end) ')'])
    
    answer2 = inputdlg({'input file name','input condition name'});
    filename2 = [selpath '/' char(answer2(1)) '.xlsx'];
    trialname2 = char(answer2(2));
    writetable(T2,filename2,'Sheet',trialname2)

    eval(['TS = table(LSI' ',ASI'  ',NSI' ')'])
    answerS = inputdlg({'input file name','input condition name'});
    filenameS = [selpath '/' char(answerS(1)) '.xlsx'];
    trialnameS = char(answerS(2));

    writetable(TS,filenameS,'Sheet',trialnameS)
  

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
    filttype = 'Low'
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

function SubMin_Callback(hObject, eventdata, handles)
    for t = 1:handles.nTrials
    eval(['data = handles.T' num2str(t) '_data;']);
    data(:,2) = data(:,2)-min(data(:,2));
    data(:,3) = data(:,3)-min(data(:,3));
    eval(['handles.T' num2str(t) '_data = data;']);
    
    eval(['data = handles.T' num2str(t) '_data_r;']);
    data(:,2) = data(:,2)-min(data(:,2));
    data(:,3) = data(:,3)-min(data(:,3));
    eval(['handles.T' num2str(t) '_data_r = data;']);
    
    eval(['data = handles.T' num2str(t) '_Stored_data;']);
    data(:,2) = data(:,2)-min(data(:,2));
    data(:,3) = data(:,3)-min(data(:,3));
    eval(['handles.T' num2str(t) '_Stored_data = data;']);
    end
    
    guidata(hObject,handles);
    IC_Thresh_Callback(hObject, eventdata, handles);
    guidata(hObject,handles);
    plot_data(hObject, eventdata, handles) 

function Movement_Callback(hObject, eventdata, handles)
MovementType = get(handles.Movement,'Value');

if MovementType == 1 %Bilateral landing
    set(handles.PercentGCP_IMP, 'Value',1)
    set(handles.PercentGCP_PIF, 'Value',.3)
    set(handles.TimeFC_IMP, 'Value',2)
    set(handles.TimeFC_PIF, 'Value',2)
    
    set(handles.Impulse_PercentShow, 'String', '100 / 100%')
    set(handles.ImpactPeak_PercentShow, 'String', '30 / 100%')
    set(handles.Impulse_TimeShow, 'String', '200 / 500 ms')
    set(handles.ImpactPeak_TimeShow, 'String', '200 / 500 ms')
elseif MovementType == 2 %Unilateral landing
    set(handles.PercentGCP_IMP, 'Value',1)
    set(handles.PercentGCP_PIF, 'Value',.3)
    set(handles.TimeFC_IMP, 'Value',2)
    set(handles.TimeFC_PIF, 'Value',2)
    
    set(handles.Impulse_PercentShow, 'String', '100 / 100%')
    set(handles.ImpactPeak_PercentShow, 'String', '30 / 100%')
    set(handles.Impulse_TimeShow, 'String', '200 / 500 ms')
    set(handles.ImpactPeak_TimeShow, 'String', '200 / 500 ms')
elseif MovementType == 3 %Walking 
    set(handles.PercentGCP_IMP, 'Value',1)
    set(handles.PercentGCP_PIF, 'Value',.6)
    set(handles.TimeFC_IMP, 'Value',2)
    set(handles.TimeFC_PIF, 'Value',2)   
    
    set(handles.Impulse_PercentShow, 'String', '100 / 100%')
    set(handles.ImpactPeak_PercentShow, 'String', '60 / 100%')
    set(handles.Impulse_TimeShow, 'String', '200 / 500 ms')
    set(handles.ImpactPeak_TimeShow, 'String', '200 / 500 ms')
elseif MovementType == 4 %squatting

elseif MovementType == 5 %Sit to Stand
    set(handles.PercentGCP_IMP, 'Value',1)
    set(handles.PercentGCP_PIF, 'Value',.6)
    set(handles.TimeFC_IMP, 'Value',2)
    set(handles.TimeFC_PIF, 'Value',2)

    set(handles.Impulse_PercentShow, 'String', '100 / 100%')
    set(handles.ImpactPeak_PercentShow, 'String', '60 / 100%')
    set(handles.Impulse_TimeShow, 'String', '200 / 500 ms')
    set(handles.ImpactPeak_TimeShow, 'String', '200 / 500 ms')
    
end

function TimeFC_PIF_Callback(hObject, eventdata, handles)
set(handles.ImpactPeak_TimeShow, 'String', [num2str(round(get(handles.TimeFC_PIF, 'Value')*100)) ' / 500 ms'])

function TimeFC_IMP_Callback(hObject, eventdata, handles)
set(handles.Impulse_TimeShow, 'String', [num2str(round(get(handles.TimeFC_IMP, 'Value')*100)) ' / 500 ms'])

function PercentGCP_PIF_Callback(hObject, eventdata, handles)
set(handles.ImpactPeak_PercentShow, 'String', [num2str(round(get(handles.PercentGCP_PIF, 'Value')*100)) ' / 100%'])

function PercentGCP_IMP_Callback(hObject, eventdata, handles)
set(handles.Impulse_PercentShow, 'String', [num2str(round(get(handles.PercentGCP_IMP, 'Value')*100)) ' / 100%'])




% Not needed
function DominantLimb_Callback(hObject, eventdata, handles)

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

function Movement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Movement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PercentGCP_IMP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PercentGCP_IMP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function PercentGCP_PIF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PercentGCP_PIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function TimeFC_IMP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeFC_IMP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function TimeFC_PIF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeFC_PIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function DominantLimb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DominantLimb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
