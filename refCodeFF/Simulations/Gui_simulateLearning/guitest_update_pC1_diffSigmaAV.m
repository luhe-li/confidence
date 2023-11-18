function varargout = guitest_update_pC1_diffSigmaAV(varargin)
% GUITEST_UPDATE_PC1_DIFFSIGMAAV MATLAB code for guitest_update_pC1_diffSigmaAV.fig
%      GUITEST_UPDATE_PC1_DIFFSIGMAAV, by itself, creates a new GUITEST_UPDATE_PC1_DIFFSIGMAAV or raises the existing
%      singleton*.
%
%      H = GUITEST_UPDATE_PC1_DIFFSIGMAAV returns the handle to a new GUITEST_UPDATE_PC1_DIFFSIGMAAV or the handle to
%      the existing singleton*.
%
%      GUITEST_UPDATE_PC1_DIFFSIGMAAV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITEST_UPDATE_PC1_DIFFSIGMAAV.M with the given input arguments.
%
%      GUITEST_UPDATE_PC1_DIFFSIGMAAV('Property','Value',...) creates a new GUITEST_UPDATE_PC1_DIFFSIGMAAV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guitest_update_pC1_diffSigmaAV_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guitest_update_pC1_diffSigmaAV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guitest_update_pC1_diffSigmaAV

% Last Modified by GUIDE v2.5 26-Feb-2021 00:25:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guitest_update_pC1_diffSigmaAV_OpeningFcn, ...
                   'gui_OutputFcn',  @guitest_update_pC1_diffSigmaAV_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes just before guitest_update_pC1_diffSigmaAV is made visible.
function guitest_update_pC1_diffSigmaAV_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guitest_update_pC1_diffSigmaAV (see VARARGIN)

% Choose default command line output for guitest_update_pC1_diffSigmaAV
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guitest_update_pC1_diffSigmaAV wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guitest_update_pC1_diffSigmaAV_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
showLaTexSymbols(handles)

% --- Executes on button press in button_start.
function button_start_Callback(hObject, eventdata, handles)
%clear the progress bar and reset the plots
clearProgressBar(handles); 
for i = 1:4; resetPlots(handles, i);end

%--------------------------------------------------------------------------
%                           Get subject info
%--------------------------------------------------------------------------
% savedP      = NaN(7,7);
% savedGrid   = cell(7,7);
% savedD      = cell(7,2);
% cell_savedP = {savedP, savedGrid, savedD};
% save('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
global inputN param_name set_name d d_bool
dict_subjN       = [3,4,5,6,8,9,11,12,13,15,16,17,18];
dict_subjI       = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX','ZL'};
dict_cond_order  = [2,1;1,2;2,1;1,2;2,1;1,2;2,1;2,1;1,2;2,1;1,2;1,2;2,1];%1: cong; 2: incong
subjI            = handles.popupmenu_subjI.String{handles.popupmenu_subjI.Value}; 
inputN           = find(contains(dict_subjI, subjI));
subjN            = dict_subjN(inputN); 
sesNum           = dict_cond_order(inputN,:); %incongruent condition
Cond             = {'congruent','incongruent'};
numCond          = length(Cond);
Modality         = {'A','V'};

%--------------------------------------------------------------------------
%                     Load the best-fitting M and Theta
%--------------------------------------------------------------------------
%best model:
%'PW':full-MA-MAP,       'SW':samePC1cond-MA-posterior, 'HL':samePC1cond-MS-posterior,
%'YZ':full-MA-MAP,       'NH':samePC1pre-MA-MAP,        'ZZ':full-MA-MAP,
%'BB':full-MS-posterior, 'ZY':full-MA-MAP,              'MR':samePC1pre-MA-MAP,
%'AD':full-MA-MAP,       'SM':full-MA-MAP,              'SX':samePC1pre-MA-m
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                'ModelFitting/Fits']));
Mcmp        = load('Results_modelComparison_overall.mat', 'ModelComparison');
bestM       = Mcmp.ModelComparison{end}{inputN};
%e.g., 'SW-samePC1cond-strategyMAP_MA_strategyUnity_posteriorC1'
dash_idx    = find(bestM == '-');
model_pC1   = bestM((dash_idx(1)+1):(dash_idx(2)-1)); %full model, samePC1pre, samePC1cond
seg_idx     = find(bestM == '_'); 
ds_locResp  = bestM((seg_idx(1)+1):(seg_idx(2)-1)); %'MS' or 'MA'
ds_unityJdg = bestM((seg_idx(3)+1):end); %measurements, MAP, posteriorC1

%load files for data from the bimodal spatial discrimination task
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                        'ModelFitting/Fits/',subjI,'/Fits_locResp_',...
                        'conditioned_unityJdg_jointly']));
C                     = load(['ModelFitting_updatePrior_',model_pC1,...
                        '_strategyMAP_', ds_locResp,'_strategyUnity_',...
                        ds_unityJdg,'_sub', num2str(subjN), '.mat'],...
                        'ModelFitting');
%get the following parameter estimates from model fitting
estimatedP            = C.ModelFitting{end}.P_f;
param.sigmaP_spatial  = 100;
param.muP_spatial     = 0;
param.criterion       = 0.5;
param.a_A             = estimatedP(1);
param.b_A             = estimatedP(2);
c_unityJdg            = estimatedP(7);
%given the best model, get the estimated pCommon
switch model_pC1
    case 'fullModel'
        emp_pCommon_pre  = [estimatedP(8), estimatedP(10)];%1st: cong; 2nd: incong
        emp_pCommon_post = [estimatedP(9), estimatedP(11)];
    case 'samePC1pre'
        emp_pCommon_pre  = [estimatedP(8), estimatedP(8)];
        emp_pCommon_post = [estimatedP(9), estimatedP(10)];   
    case 'samePC1cond'
        emp_pCommon_pre  = [estimatedP(8), estimatedP(9)];
        emp_pCommon_post = [estimatedP(8), estimatedP(9)];       
end

%--------------------------------------------------------------------------
%                 Load the data from the learning phase
%--------------------------------------------------------------------------
empReportingC1 = NaN(1,numCond);
%get the data we need
for i = 1:numCond 
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Adaptation v2/Data/',subjI]));
    if i == 1
        F = load(['Adaptation_',Cond{i}, '_sub', num2str(subjN),'_session',...
            num2str(sesNum(i)),'.mat'], 'Adaptation_data');
        D.AVpairs(i,1,:)  = F.Adaptation_data{4}.arrangedLocs_deg; %A loc
        D.AVpairs(i,3,:)  = F.Adaptation_data{3}.arrangedLocs_deg; %V loc
        D.AVpairs(i,2,:)  = zeros(1,160); D.AVpairs(i,4,:) = zeros(1,160);
        empReportingC1(i) = sum(F.Adaptation_data{end}.unity==1)/...
            sum(~isnan(F.Adaptation_data{end}.unity));
    else
        F = load(['Adaptation_',Cond{i}, '_sub', num2str(subjN),'_session',...
            num2str(sesNum(i)),'.mat'], ['Adaptation_',Cond{i},'_data']);
        D.AVpairs(i,1,:)  = F.Adaptation_incongruent_data{4}.arrangedLocs_deg; %A loc
        D.AVpairs(i,3,:)  = F.Adaptation_incongruent_data{3}.arrangedLocs_deg; %V loc
        D.AVpairs(i,2,:)  = -F.Adaptation_incongruent_data{3}.timing_relative./2; %A time relative
        D.AVpairs(i,4,:)  = F.Adaptation_incongruent_data{3}.timing_relative./2; %V time relative
        empReportingC1(i) = sum(F.Adaptation_incongruent_data{end}.unity==1)/...
            sum(~isnan(F.Adaptation_incongruent_data{end}.unity));
    end
end
d_emp         = vertcat(emp_pCommon_post, empReportingC1);
D.totalTrials = size(D.AVpairs,3);
D.numSims     = 1e2;

%--------------------------------------------------------------------------
%                     Read off the free parameters
%--------------------------------------------------------------------------
param_name = {'sigma_deltaT','muP_deltaT_C1','sigmaP_deltaT_C1',...
    'muP_deltaT_C2','sigmaP_deltaT_C2','alpha_pC1_congruent',...
    'alpha_pC1_incongruent'};
roundup_digit = [0,0,0,0,0,3,3];

if get(handles.radiobutton_sliders,'Value') == 1 %read off from the sliders
    for i = 1:length(param_name)
        eval(['param.',param_name{i},'=handles.slider_',param_name{i},'.Value;']);
        %change the values in the text boxes accordingly
        eval(['handles.text_', param_name{i}, '.String = num2str(round(param.',...
            param_name{i},',',num2str(roundup_digit(i)),'));']);
    end
elseif get(handles.radiobutton_texts,'Value') == 1 %read off from the text boxes
    for i = 1:length(param_name)
        eval(['param.',param_name{i},'=str2double(handles.text_', param_name{i},...
            '.String);']);
        %change the positions of the sliders accordingly
        eval(['handles.slider_', param_name{i},'.Value=param.',param_name{i},';']);
    end
elseif get(handles.radiobutton_saved,'Value') == 1
    %if we have already found the combinations of parameters that can
    %predict our results, load the file
    %------------------This saved file has 3 main cells--------------------
    %1. saved parameters, as specified in the cell "param_name"
    %2. saved setting information, as specificed in the cell "set_name"
    %3. saved data
    %   (1). simulated p_{C=1,post} and p(reporting C=1)
    %   (2). boolean (whether simulated value is close to the empirical data)
    %----------------------------------------------------------------------
    C = load('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
    p = C.cell_savedP{1}(inputN,:);
    if isnan(p(1)); errordlg('This file does not exist!'); return; 
    else
        for i = 1:length(param_name)
            eval(['[param.', param_name{i},',handles.slider_',param_name{i},...
                '.Value, handles.text_', param_name{i}, '.String] = deal(p(',...
                num2str(i),'));']);
        end
    end
end

%--------------------------------------------------------------------------
%                     Read off the grid information
%--------------------------------------------------------------------------
set_name = {'sigma_AV_A_lb_ub_congruent','sigma_AV_V_lb_ub_congruent',...
    'sigma_AV_A_lb_ub_incongruent','sigma_AV_V_lb_ub_incongruent','tol_pC1',...
    'tol_reportingC1','numBins'};
if get(handles.radiobutton_saved,'Value') == 1 %load saved setting information
    C = load('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
    setting = C.cell_savedP{2}(inputN,:);
    for i = 1:length(set_name)
        eval([set_name{i}, '= setting{i};']);
        eval(['handles.text_', set_name{i}, '.String = ',set_name{i},';']);
    end
else
    for i = 1:length(set_name)
        eval([set_name{i},'= handles.text_', set_name{i},'.String;']);
    end
end
for i = length(set_name):-1:1
    if i < 5
        eval(['[lb,ub] = getBds(',set_name{i},');']);
        eval([strcat([set_name{i}(1:11), set_name{i}(18:end)]),...
            ' = linspace(lb,ub,numBins);']);
    else
        eval([set_name{i},'= str2double(', set_name{i},');']);
    end
end

%progress bar
progress_mark   = 0.1:0.1:1;
idx_progressBar = arrayfun(@(idx) find(abs((1:numBins)./numBins-idx)<1e-3,1),...
    progress_mark);

%--------------------------------------------------------------------------
%                        Main code (run simulation)
%--------------------------------------------------------------------------
if get(handles.radiobutton_saved,'Value') == 1
    C      = load('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
    d      = C.cell_savedP{3}{inputN,1};
    d_bool = C.cell_savedP{3}{inputN,2};
    for i = 1:length(progress_mark); addProgressBar(handles, i);end
else
    [propC1, pCommon_end] = deal(NaN(numCond, numBins,numBins,length(D.numSims)));
    for i = 1:numBins %A
        %load the progress bar
        if ismember(i, idx_progressBar); addProgressBar(handles, find(i==idx_progressBar));end
        for j = 1:numBins %V
            for k = 1:D.numSims
                for l = 1:numCond
                    eval(['param.sigma_spatial_AV_A = sigma_AV_A_', Cond{l},'(i);']);
                    eval(['param.sigma_spatial_AV_V = sigma_AV_V_', Cond{l},'(j);']);
                    eval(['param.alpha = param.alpha_pC1_', Cond{l}, ';']);
                    param.pCommon = emp_pCommon_pre(l); 
                    [propC1(l,i,j,k), pCommon_end(l,i,j,k)] = ...
                        simulateLearningPhase_update_pC1(...
                        param, D, l, ds_unityJdg, c_unityJdg, ds_locResp);
                end
            end
        end
    end

    %compare simulations with the empirical data
    sim_pCommon_post = mean(pCommon_end,4); 
    bool_withinTol_pCommon = abs(sim_pCommon_post - repmat(emp_pCommon_post',...
        [1, numBins, numBins])) < tol_pC1;

    simReportingC1 = mean(propC1,4); 
    bool_withinTol_ReportingC1 = abs(simReportingC1 - repmat(empReportingC1',...
        [1, numBins, numBins])) < tol_reportingC1;
    
    d = {sim_pCommon_post, simReportingC1};
    d_bool = {bool_withinTol_pCommon, bool_withinTol_ReportingC1};
end

%--------------------------------------------------------------------------
%                                Plotting
%--------------------------------------------------------------------------
%before plotting, create custom colormap
c_low = [36, 116, 185]./255; c_high = [253, 111, 63]./255;
for i = 1:numCond %2 conditions
    for j = 1:length(d)
        [cticks, cMap] = customize_cmap(min(min(d{j}(i,:,:))), ...
            max(max(d{j}(i,:,:))),d_emp(j,i), c_low, c_high);
        eval(['axes(handles.axes',num2str((i-1)*numCond+j),')']); 
        imagesc(eval(['sigma_AV_V_',Cond{i}]), eval(['sigma_AV_A_',Cond{i}]),...
            squeeze(d{j}(i,:,:))); 
        c = colorbar; set(c,'XTick',cticks); caxis([cticks(1),cticks(end)]); 
        eval(['colormap(handles.axes',num2str((i-1)*numCond+j),',cMap)']); hold on;
        drawLines(eval(['sigma_AV_V_',Cond{i}]), eval(['sigma_AV_A_',Cond{i}]), ...
            squeeze(d_bool{j}(i,:,:))); hold off;
        xlabel('\sigma_{AV,V}'); ylabel('\sigma_{AV,A}');
        eval(['yticks(round(sigma_AV_V_', Cond{i},'(1:3:end),2));']);
        eval(['yticks(round(sigma_AV_A_', Cond{i},'(1:3:end),2));']);
        if j == 1
            title(['p_{C=1, pre} = ', num2str(round(emp_pCommon_pre(i),3)),...
                ', p_{C=1, post} = ', num2str(round(d_emp(j,i),3))]);
        else
            title(['p(reporting C=1) = ',num2str(round(d_emp(j,i),3))]);
        end
        set(gca,'FontSize',15);
    end
end

%--------------------------------------------------------------------------
%                       Showing Latex symbols
%--------------------------------------------------------------------------
function showLaTexSymbols(handles)
symbs1 = {'\sigma_{\Deltat}', '\mu_{P_{\Deltat},C=1}','\sigma_{\Deltat,C=1}',...
    '\mu_{P_{\Deltat},C=2}','\sigma_{\Deltat,C=2}',...
    '\alpha_{p_{C=1}} (cong)','\alpha_{p_{C=1}} (incong)'};
for i = 1:length(symbs1)
    eval(['axes(handles.axes',num2str(20 + i),')']);
    text(0,0.5,symbs1{i},'fontSize',16);
end
symbs2 = {'\sigma_{AV,A,cong} (lb,ub):', '\sigma_{AV,V,cong} (lb,ub):',...
    '\sigma_{AV,A,incong} (lb,ub):','\sigma_{AV,V,incong} (lb,ub):',...
    'tol_1:', 'tol_2:','#grids:'};
for j = 1:length(symbs2)
    eval(['axes(handles.axes',num2str(30 + j),')']);
    text(0,0.5,symbs2{j},'fontSize',16);    
end

function [lb,ub] = getBds(s)
comma_idx = find(s==',');
lb = str2num(s(2:(comma_idx-1))); ub = str2num(s((comma_idx+1):(end-1)));


%--------------------------------------------------------------------------
%                     Helping functions for plotting
%--------------------------------------------------------------------------
function addProgressBar(handles, idx)
eval(['axes(handles.axes',num2str(10 + idx),')']);
patch([0, 1, 1, 0], [0, 0, 1, 1], [0.47,0.67,0.19],...
    'EdgeColor',[0.47,0.67,0.19],'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);

function clearProgressBar(handles)
for idx = 1:10
    eval(['axes(handles.axes',num2str(10 + idx),')']);
    patch([0, 1, 1, 0], [0, 0, 1, 1], [0.94,0.94,0.94],...
        'EdgeColor',[0.94,0.94,0.94]);
end

function resetPlots(handles,plot_idx)
eval(['axes(handles.axes',num2str(plot_idx),')']); 
cla reset; xticks([]); yticks([]); 
eval(['handles.axes', num2str(plot_idx), '.XColor = [1,1,1];']);
eval(['handles.axes', num2str(plot_idx), '.YColor = [1,1,1];']);
text(0.35,0.5,'Loading...','fontSize',20,'Color',[0.47,0.67,0.19]);

function [c_ticks, c_Mat] = customize_cmap(v_min, v_max, v_emp, c_low, c_high)
%first check if the empirical value is within the range of the minimum and
%the maximum simulated values
n = 256;
p1 = [linspace(c_low(1),1,n)', linspace(c_low(2),1,n)',linspace(c_low(3),1,n)'];
p2 = [linspace(1,c_high(1),n)', linspace(1,c_high(2),n)',linspace(1,c_high(3),n)'];
if v_emp > v_min && v_emp < v_max
    %calculate the position (0 -> v_min -> v_emp -> v_max -> 0.99)
    d2 = v_max - v_emp; d1 = v_emp - v_min;
    if d2 > d1; c_low = p1(floor(n-d1/d2*n),:); %if v_emp is closer to v_min
    else; c_high = p2(ceil(d2/d1*n),:); end
    steps_low   = round(d1/(d1+d2)*n); steps_high = n - steps_low;
    c_Mat_part1 = [linspace(c_low(1),1,steps_low)',...
        linspace(c_low(2),1,steps_low)',linspace(c_low(3),1,steps_low)'];
    c_Mat_part2 = [linspace(1,c_high(1),steps_high)',...
        linspace(1,c_high(2),steps_high)',linspace(1,c_high(2),steps_high)'];
    c_Mat = vertcat(c_Mat_part1, c_Mat_part2);
    c_ticks = round([v_min, v_emp, v_max],2);
else
    if v_emp <= v_min %(0 -> v_emp -> v_min -> v_max -> 0.99)
        d2 = v_max - v_min; d1 = v_min - v_emp;
        c_low = p2(ceil(n-d2/(d1 + d2)*n),:);
        c_Mat = [linspace(c_low(1),c_high(1),n)',linspace(c_low(2),c_high(2),n)',...
                 linspace(c_low(3),c_high(3),n)'];
        c_ticks = round([v_min, v_max],2);
    else %(0 -> v_min -> v_max -> v_emp -> 0.99)
        d2 = v_emp - v_max; d1 = v_max - v_min;
        c_high = p1(floor(d1/(d1 + d2)*n),:);
        c_Mat = [linspace(c_low(1),c_high(1),n)',linspace(c_low(2),c_high(2),n)',...
                 linspace(c_low(3),c_high(3),n)'];
        c_ticks = round([v_min, v_max],2);
    end
end

function [c_ticks, c_Mat] = customize_cmap_new(v_min, v_max, v_emp, c_low, c_high)
%first check if the empirical value is within the range of the minimum and
%the maximum simulated values
n = 256;
p1 = [linspace(c_low(1),1,n)', linspace(c_low(2),1,n)',linspace(c_low(3),1,n)'];
p2 = [linspace(1,c_high(1),n)', linspace(1,c_high(2),n)',linspace(1,c_high(3),n)'];
lin_vmin = linspace(v_emp-1, v_emp,n);
lin_vmax = linspace(v_emp, v_emp+1,n);
[~,idx_lb] = min(abs(lin_vmin - v_min)); %p1_lb = p1(idx_lb,:);
[~,idx_ub] = min(abs(lin_vmax - v_max)); %p1_ub = p2(idx_ub,:);
c_Mat = vertcat(p1(idx_lb:end,:), p2(idx_ub:end,:));
c_ticks = round([v_min, v_emp, v_max],2);

    
function drawLines(xTicks, yTicks, bool)
binsize_x = diff(xTicks(1:2))/2;
binsize_y = diff(yTicks(1:2))/2;
for i = 1:size(bool,1)
    %left
    idx_first_1_left = find(bool(i,:)==1,1);
    if ~isempty(idx_first_1_left)
        plot([xTicks(idx_first_1_left), xTicks(idx_first_1_left)]-binsize_x,...
            [yTicks(i)- binsize_y, yTicks(i)+ binsize_y] , 'Color',...
            [0.47,0.67,0.19],'lineWidth',5); hold on;
    end
    idx_strpattern_left = strfind(bool(i,:), [1,0]);
    if ~isempty(idx_strpattern_left)
        for ii = 1:length(idx_strpattern_left)
            plot([xTicks(idx_strpattern_left(ii)), xTicks(idx_strpattern_left(ii))]+...
                binsize_x, [yTicks(i)- binsize_y, yTicks(i)+ binsize_y],...
                'Color',[0.47,0.67,0.19],'lineWidth',5); hold on;
        end
    end
    %right
    idx_first_1_right = size(bool,1)- find(fliplr(bool(i,:))==1,1) + 1;
    if ~isempty(idx_first_1_right)
        plot([xTicks(idx_first_1_right), xTicks(idx_first_1_right)]+binsize_x,...
            [yTicks(i)- binsize_y, yTicks(i)+ binsize_y] , 'Color',...
            [0.47,0.67,0.19],'lineWidth',5); hold on;
    end 
    idx_strpattern_right =strfind(bool(i,:), [0,1]);
    if ~isempty(idx_strpattern_right)
        for ii = 1:length(idx_strpattern_right)
            plot([xTicks(idx_strpattern_right(ii)), ...
                xTicks(idx_strpattern_right(ii))]+binsize_x, ...
                [yTicks(i)- binsize_y, yTicks(i)+ binsize_y],...
                'Color',[0.47,0.67,0.19],'lineWidth',5); hold on;
        end
    end
end
for j = 1:size(bool,2)
    %up
    idx_first_1_up = find(bool(:,j)==1,1);
    if ~isempty(idx_first_1_up)
        plot([xTicks(j)-binsize_x, xTicks(j)+binsize_x], [yTicks(idx_first_1_up),...
            yTicks(idx_first_1_up)] - binsize_y,'Color',[0.47,0.67,0.19],...
            'lineWidth',5); hold on;
    end 
    idx_strpattern_up = strfind(bool(:,j)', [1,0]);
    if ~isempty(idx_strpattern_up)
        for jj = 1:length(idx_strpattern_up)
            plot([xTicks(j)-binsize_x, xTicks(j)+binsize_x], ...
                [yTicks(idx_strpattern_up(jj)),yTicks(idx_strpattern_up(jj))] +...
                binsize_y, 'Color',[0.47,0.67,0.19],'lineWidth',5); hold on;
        end
    end
    %down
    idx_first_1_down = size(bool,2)- find(flipud(bool(:,j))==1,1) + 1;
    if ~isempty(idx_first_1_down)
        plot([xTicks(j)-binsize_x, xTicks(j)+binsize_x], [yTicks(idx_first_1_down),...
            yTicks(idx_first_1_down)] + binsize_y, 'Color',[0.47,0.67,0.19],...
            'lineWidth',5); hold on;
    end
    idx_strpattern_down = strfind(bool(:,j)', [0,1]);
    if ~isempty(idx_strpattern_down)
        for jj = 1:length(idx_strpattern_down)
            plot([xTicks(j)-binsize_x, xTicks(j)+binsize_x], ...
                [yTicks(idx_strpattern_down(jj)),yTicks(idx_strpattern_down(jj))] +...
                binsize_y, 'Color',[0.47,0.67,0.19],'lineWidth',5); hold on;
        end
    end
end


% --- Executes on slider movement.
function slider_sigmaP_deltaT_C1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sigmaP_deltaT_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_sigmaP_deltaT_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sigmaP_deltaT_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_sigmaP_deltaT_C2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider_sigmaP_deltaT_C2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu_subjI.
function popupmenu_subjI_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_subjI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),... 
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_sigma_deltaT_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigma_deltaT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text_sigmaP_deltaT_C1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigmaP_deltaT_C1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text_sigmaP_deltaT_C2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigmaP_deltaT_C2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_sliders.
function radiobutton_sliders_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.radiobutton_texts.Value = 0;
    handles.radiobutton_saved.Value = 0;
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in radiobutton_texts.
function radiobutton_texts_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.radiobutton_sliders.Value = 0;
    handles.radiobutton_saved.Value = 0;
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on button press in radiobutton_saved.
function radiobutton_saved_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.radiobutton_sliders.Value = 0;
    handles.radiobutton_texts.Value = 0;
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on slider movement.
function slider_muP_deltaT_C1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider_muP_deltaT_C1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_muP_deltaT_C2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider_muP_deltaT_C2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'),... 
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function text_muP_deltaT_C2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_muP_deltaT_C2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text_muP_deltaT_C1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_muP_deltaT_C1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_sigma_deltaT_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider_sigma_deltaT_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
global inputN param_name set_name d d_bool
if get(hObject,'Value')
    C           = load('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
    cell_savedP = C.cell_savedP;
    %------------------This saved file has 3 main cells--------------------
    %1. saved parameters, as specified in the cell "param_name"
    for i = 1:length(param_name)
        cell_savedP{1}(inputN, i) = eval(['handles.slider_',param_name{i},'.Value']);
    end
    %2. saved setting information, as specificed in the cell "set_name"
    for j = 1:length(set_name)
        cell_savedP{2}{inputN, j} = eval(['handles.text_', set_name{j},'.String']);
    end
    %3. saved data
    %   (1). simulated p_{C=1,post} and p(reporting C=1)
    %   (2). boolean (whether simulated value is close to the empirical data)   
    cell_savedP{3}{inputN,1} = d; cell_savedP{3}{inputN,2} = d_bool;  
    save('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
    %---------------------------------------------------------------------- 
    handles.pushbutton_save.Value = 0;
end


function text_tol_pC1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_tol_pC1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_tol_reportingC1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_tol_reportingC1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_alpha_pC1_incongruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider_alpha_pC1_incongruent_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function text_alpha_pC1_incongruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_alpha_pC1_incongruent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider_alpha_pC1_congruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider_alpha_pC1_congruent_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function text_alpha_pC1_congruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_alpha_pC1_congruent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text_sigma_AV_A_lb_ub_congruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigma_AV_A_lb_ub_congruent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text_sigma_AV_V_lb_ub_congruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigma_AV_V_lb_ub_congruent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_sigma_AV_A_lb_ub_incongruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigma_AV_A_lb_ub_incongruent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),... 
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_sigma_AV_V_lb_ub_incongruent_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_sigma_AV_V_lb_ub_incongruent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_numBins_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_numBins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
