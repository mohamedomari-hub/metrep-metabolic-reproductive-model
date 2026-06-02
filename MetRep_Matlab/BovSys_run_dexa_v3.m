function BovSys_run_dexa_v3()
% Matlab code: MetRep 
%
%% Introduction:
%   This is an interactive code that allows user to simulate a mechanistic ODE model of a dairy cow that couples:
%	  1) an estrous-cycle model (hypothalamus–pituitary–ovary feedbacks that generate ~21-day cyclicity)
%	  2) a glucose–insulin model with compartments that describe glucose flow (blood, liver, storage, fat, plus milk/maintenance usage)
%     3) and a dexamethasone PK/PD module to reproduce observed metabolic responses after an IM injection.
%
%% Goal:
%   To build a robust mechanistic model structure that enables qualitative comparison with experimental observations and exploration of
%   system-level responses to nutritional and pharmacological perturbations, rather than precise quantitative prediction.
%
%% Dependencies (must be on MATLAB path):
%     1) BovSys_para_dexa_v3()         : returns parameter vector 'par'
%     2) BovSys_Equa_dexa_v3(t,y,...)  : ODE right-hand side (the model)
%     3) DM.mat containing DMdataset   : weekly DMI data
%     4) ML.mat containing MLKdataset  : weekly milk yield data
%     5) Data_P4_Hol.mat / Data_IGF_Hol.mat : experimental hormone data (comparison plots)
%
%% State vectors (length 25):
%   1..16  endocrine/reproductive states
%   17..22 metabolic states
%   23..25 PK/PD states: A_dep, A_cent, C_e
%
%% Simulation scenarios:
%   The ODE solver (ode15s) allows simulation of:
%
%       1) Non-lactating: 
%                - standard diet (11.7kg DMI/day), 60 days, PK/PD off
%                - standard diet, 3 days, PK/PD on (single dose administered at day 0)
%                - acute diet (DMI reduced to 33% for 15 days), 90 days, PK/PD off
%                - chronic diet (DMI reduced to 58% for 30 weeks), 330 days, PK/PD off
%
%       2) Lactating:
%                – Two values of glucose content in the DMI (c0 = 20% and 35%), PK/PD off  
%                - One value of glucose content in the DMI (c0 = 20%) PK/PD on (single dose administered at day 50)
%
%% Start
clc; close all;      % clear command window + close all existing figures
tic                  % start timer to measure runtime

par = BovSys_para_dexa_v3();   % load full parameter vector used in the ODE

%% ===================== INITIAL CONDITIONS =====================
% y0 = initial state vector (25x1)
% IMPORTANT: the order must match BovSys_Equa_dexa_v3 exactly.
% Any mismatch here will change the biological meaning of each entry.
y0 = [ ...
    0.667; 0.551; 0.316; 0.395; 1.000; 0.642; 1.000; 0.00506; 0.0;   0.004; ...
    0.89;  0.826; 0;     0.0183;0.35;  0.48;  15.5;  0.48;   15e4;  110; ...
    535;   105.;  0;     0;     0];

%% Names used for plot titles / labeling only (no effect on simulation)
% Indices:
%   10 = P4, 11 = E2, 16 = IGF-I
%   17 = Insulin, 18 = Blood glucose
%   23 = A_dep, 24 = A_cent, 25 = C_e
Ynames = { ...
    'GnRhH','GnRhP','FSHP','Blood-FSH','LHP','LH','Follicle','PGF', ...
    'CL','P4','E2','INH','ENZ','OXT','IOF','IGF-I', ...
    'Insulin','Glucose_{blood}','Fat','Glucose_{liver}', ...
    'Glu_{store}','Glucagon','A_{dep}','A_{cent}','C_e'};

%% ===================== DATA =====================
% Load weekly datasets:
% DMdataset(:,1) expected = dry matter intake (DMI)
% MLKdataset(:,1) expected = time (weeks)
% MLKdataset(:,2) expected = milk yield
load('DM.mat','DMdataset');
load('ML.mat','MLKdataset');

dm_week = DMdataset(:,1);      % weekly dry matter intake
mlk_week = MLKdataset(:,2);    % weekly milk yield
wk_time = MLKdataset(:,1);     % weekly time axis (weeks)
at_week_days = wk_time * 7;    % convert weeks -> days for interpolation

%% ===================== EXPERIMENTAL DATA =====================
% These are only used in the NON-LACTATING STANDARD case when PK/PD is OFF
% They allow comparing model P4 (state 10) and IGF (state 16) with data points
load('Data_P4_Hol.mat','Data_P4_Hol');
load('Data_IGF_Hol.mat','Data_IGF_Hol');

%% ===================== TIME GRID =====================
% User chooses simulation horizon in days
day = input('Simulate how many days? ');

% Solver step (fixed output grid):
% 1/48 day = 30 minutes
STEP = 1/48;

% Time grid where solution will be reported
% Note: ode15s internally adapts step size, but outputs at these times.
tspan = 0:STEP:day;

% 'at' is passed to ODE as an alias of tspan (used for indexing/interp inside ODE)
at = tspan;

% Interpolate weekly DMI + milk yield onto daily/half-hour grid
% 'extrap' allows values beyond provided weekly data range
dm_l  = interp1(at_week_days, dm_week,  tspan,'linear','extrap');
mlk_l = interp1(at_week_days, mlk_week, tspan,'linear','extrap');

%% ===================== DOSE =====================
% BW is body weight in kg (used for dose scaling)
BW = 600;

% Dose_ng is dexamethasone dose in nanograms:
% 0.02 (mg/kg) * BW (kg) * 1e6 (ng/mg) = total ng dose
Dose_ng = 0.02 * BW * 1e6;

%% ===================== MODE =====================
% User chooses lactating vs non-lactating.
% This affects:
%   - the forcing data (milk yield)
%   - the dosing day (0 vs 50)
%   - initial condition adjustment (y0(14)=2.5 for lactating)
str = input('Mode lactating (l) or non-lactating (n)? l/n: ','s');

switch str

%% =======================================================================
case 'n'
    % NON-LACTATING branch
    % c0 fixed at 0.08 here (baseline condition parameter)
    c0 = 0.08;

    % No milk in non-lactating mode
    mlk_n = zeros(size(tspan));

    % User selects feeding scenario
    mode = input('Non-lactating: standard (s) / acute (a) / chronic (c): ','s');

    if mode=='s'
        % STANDARD feeding: constant DMI profile
        dm_n = 11700 * ones(size(tspan));

        % Ask whether PK/PD is linked
        % If yes: run a PK/PD ON simulation (T1,Y1)
        % Always: run a PK/PD OFF simulation (T0,Y0)
        use_pkpd = input('Link PK/PD? y/n: ','s')=='y';

        % Standard case uses dose at day 0 (immediate)
        dose_day = 0;

        % Run PK/PD OFF (always done, used as baseline and for plots)
        [T0,Y0] = solve_with_dose(tspan,y0,par,'n',at,dm_n,mlk_n,c0,false,dose_day,Dose_ng,STEP);

        % If PK/PD requested: run PK/PD ON, otherwise copy OFF results
        if use_pkpd
            [T1,Y1] = solve_with_dose(tspan,y0,par,'n',at,dm_n,mlk_n,c0,true,dose_day,Dose_ng,STEP);
        else
            T1=T0; Y1=Y0;
        end

        % Plot endocrine hormone dynamics 1..16
        % Overlays OFF vs ON (if ON exists; otherwise OFF vs OFF identical)
        plot_all_requested('Non-lactating | standard',T0,Y0,T1,Y1,Ynames,par,use_pkpd);

        % ===== EXTRA PhD DATA COMPARISON =====
        % Only performed when PK/PD is OFF (because the PhD comparison figure
        % was built for the original baseline condition).
        if ~use_pkpd
            % time shift (sft) aligns data sampling days with model time axis
            sft = 18;
            figure('Name','P4 & IGF-1 vs data'); clf;

            % P4 is state 10
            subplot(2,2,1);
            plot(T0,Y0(:,10),'LineWidth',2); hold on;
            % Data scaling: 0.3*Data_P4_Hol(:,2) used in your PhD figure convention
            plot(Data_P4_Hol(1:8,1)+sft,0.3*Data_P4_Hol(1:8,2),'ro','MarkerFaceColor','r');
            ylabel('P4'); grid off;

            % IGF-1 is state 16
            subplot(2,2,2);
            plot(T0,Y0(:,16),'LineWidth',2); hold on;
            % Data scaling: /300 used in your PhD figure convention
            plot(Data_IGF_Hol(1:8,1)+sft,Data_IGF_Hol(1:8,2)/300,'ro','MarkerFaceColor','r');
            ylabel('IGF-1'); grid off;
        end

        % ===== METABOLIC (PKPD OFF) =====
        % Plots metabolic states 17..22 only for OFF case
        if ~use_pkpd
            idx = [17 18 19 20 21 22];
            figure('Name','Non-lactating | Metabolic'); clf;
            for k=1:6
                subplot(2,3,k); plot(T0,Y0(:,idx(k)),'LineWidth',1.8);
                title(Ynames{idx(k)}); grid off;
            end
        end

        % ===== METABOLIC + PKPD (ON) =====
        % If PK/PD is ON:
        %   - plot metabolic states from ON simulation
        %   - plot PK readouts: plasma concentration and effect compartment
        if use_pkpd
            idx = [17 18 19 20 21 22];
            figure('Name','Non-lactating | Metabolic + PKPD'); clf;
            for k=1:6
                subplot(2,3,k); plot(T1,Y1(:,idx(k)),'LineWidth',1.8);
                title(Ynames{idx(k)}); grid off;
            end

            % Vd used to convert A_cent (ng) -> C_plasma (ng/mL)
            % Vd = 1.105 * BW * 1e3 (mL)
            Vd = 1.105*BW*1e3;

            figure('Name','Non-lactating | PK'); clf;
            % State 24 = A_cent. Plasma concentration = A_cent / Vd
            subplot(2,1,1); plot(T1,Y1(:,24)/Vd,'LineWidth',2); title('C_{plasma}');
            % State 25 = C_e (effect compartment concentration)
            subplot(2,1,2); plot(T1,Y1(:,25),'LineWidth',2); title('C_e');
        end

    elseif mode=='a'
        % ACUTE feeding scenario:
        % build_dmi_acute constructs a time-varying DMI profile
        dm_n = build_dmi_acute(tspan,STEP,11700,0.33,11700,46,15);

        % In acute case here, PK/PD is OFF (use_dexa=false)
        [T,Y] = solve_with_dose(tspan,y0,par,'n',at,dm_n,mlk_n,c0,false,0,Dose_ng,STEP);

        % Plot endocrine variables (1..16)
        plot_single('Non-lactating | acute',T,Y,Ynames);

        % Plot metabolic variables (17..22)
        idx=[17 18 19 20 21 22];
        figure('Name','Acute | Metabolic'); clf;
        for k=1:6
            subplot(2,3,k); plot(T,Y(:,idx(k)),'LineWidth',1.8);
            title(Ynames{idx(k)}); grid off;
        end

    elseif mode=='c'
        % CHRONIC feeding scenario:
        % build_dmi_chronic constructs a longer restriction/refeeding profile
        dm_n = build_dmi_chronic(tspan,STEP,11700,0.58,11700*1.6,45,210,70);

        % PK/PD OFF in chronic case here
        [T,Y] = solve_with_dose(tspan,y0,par,'n',at,dm_n,mlk_n,c0,false,0,Dose_ng,STEP);

        % Plot endocrine variables (1..16)
        plot_single('Non-lactating | chronic',T,Y,Ynames);

        % Plot metabolic variables (17..22)
        idx=[17 18 19 20 21 22];
        figure('Name','Chronic | Metabolic'); clf;
        for k=1:6
            subplot(2,3,k); plot(T,Y(:,idx(k)),'LineWidth',1.8);
            title(Ynames{idx(k)}); grid off;
        end
    end

%% =======================================================================
case 'l'
    % LACTATING branch
    % Modify initial condition for state 14 (OXT) to lactating baseline
    y0(14)=2.5;

    % Ask whether PK/PD is linked
    % If yes: solve_with_dose will apply dose at day 50 via split integration
    use_pkpd = input('Link PK/PD? y/n: ','s')=='y';

    % Number of c0 values: allows overlay of multiple initial conditions
    nvals = input('How many c0 values? ');

    % Preallocate c0 array
    c0s = zeros(1,nvals);

    % Read each c0 as percentage and convert to fraction
    for i=1:nvals, c0s(i)=input(sprintf('c0 #%d (%%): ',i))/100; end

    % Lactating dose time: day 50
    dose_day = 50;

    % Store outputs for each c0 in cell arrays
    Tall=cell(nvals,1); 
    Yall=cell(nvals,1);

    % Run one simulation per c0 value
    for i=1:nvals
        [T,Y] = solve_with_dose(tspan,y0,par,'l',at,dm_l,mlk_l,c0s(i),use_pkpd,dose_day,Dose_ng,STEP);
        Tall{i}=T; 
        Yall{i}=Y;
    end

    % Plot overlay results for lactating case
    plot_overlay_lactating(Tall,Yall,c0s,Ynames,use_pkpd,dose_day,at_week_days,dm_week,mlk_week);
end

toc    % stop timer and display runtime
end

%% ======================= SOLVER =======================
function [T,Y]=solve_with_dose(tspan,y0,par,str,at,dm,mlk,c0,use_dexa,dose_day,Dose_ng,STEP)
% This helper function wraps ode15s and handles dosing consistently.
%
% Key behavior:
%   - If use_dexa = false: model is simulated with PK/PD disabled (no dosing)
%   - If dose_day <= 0: dose is applied immediately by adding to y0(23)
%   - If dose_day > 0: two-phase solve:
%         phase 1: integrate [0, dose_day]
%         add dose to A_dep (state 23)
%         phase 2: integrate [dose_day, end]
%
% This ensures correct continuity of state variables at dosing time.

opt=odeset('nonnegative',1:25);    % enforce all states to remain nonnegative

% --- PK/PD OFF: integrate full system with use_dexa flag = false ---
if ~use_dexa
    [T,Y]=ode15s(@BovSys_Equa_dexa_v3,tspan,y0,opt,par,str,at,dm,mlk,c0,false); return;
end

% --- Dose at day 0: simply add to A_dep and integrate once ---
if dose_day<=0
    y0(23)=y0(23)+Dose_ng;   % state 23 = A_dep
    [T,Y]=ode15s(@BovSys_Equa_dexa_v3,tspan,y0,opt,par,str,at,dm,mlk,c0,true); return;
end

% --- Dose after time 0: split integration at dose_day ---
t1=tspan(tspan<=dose_day); 
t2=tspan(tspan>=dose_day);

% Phase 1: integrate until dose time
[T1,Y1]=ode15s(@BovSys_Equa_dexa_v3,t1,y0,opt,par,str,at,dm,mlk,c0,true);

% Take final state at dose time and add dose to A_dep
yd=Y1(end,:)'; 
yd(23)=yd(23)+Dose_ng;

% If t2 begins exactly at dose_day, remove that duplicate point (prevents repeated time)
if abs(t2(1)-dose_day)<STEP/2, t2=t2(2:end); end

% Phase 2: integrate from dose time to end
[T2,Y2]=ode15s(@BovSys_Equa_dexa_v3,[dose_day t2],yd,opt,par,str,at,dm,mlk,c0,true);

% Merge outputs into one continuous trajectory
T=[T1;T2(2:end)]; 
Y=[Y1;Y2(2:end,:)];
end

%% ======================= PLOTS =======================
function plot_overlay_lactating(T,Y,c0s,Ynames,show_pk,dose_day,at_week_days,dm_week,mlk_week)
% Overlays multiple lactating simulations (each one corresponds to a different c0)
% Each color corresponds to one c0 value.

cols=lines(numel(c0s));   % MATLAB default colormap for distinct lines

% --- Plot endocrine/reproductive hormones/states 1..16 ---
figure('Name','Lactating | Hormones'); clf;
for k=1:16
    subplot(4,4,k); hold on;
    for i=1:numel(c0s)
        plot(T{i},Y{i}(:,k),'Color',cols(i,:),'LineWidth',1.8);
    end
    title(Ynames{k}); grid off;
end

% Legend shows which curve corresponds to which c0
legend(arrayfun(@(x)sprintf('c0=%.1f%%',100*x),c0s,'Uni',0),'Location','bestoutside');

% ============================================================
% EXTRA PLOTS WHEN PK/PD IS OFF (show_pk == false)
%   1) Insulin + Blood glucose (single figure)
%   2) P4, E2, IGF-I (single figure with 3 subplots)
% ============================================================
if ~show_pk

    
    % ----- (1) Insulin, Glucose, DMI & Milk  -----
    figure('Name','Lactating | Glucose - Insulin - Intake & Milk (PK/PD OFF)'); clf;
    
    % === Panel 1: Insulin ===
    subplot(1,3,1); hold on;
    for i = 1:numel(c0s)
        plot(T{i}, Y{i}(:,17), 'LineWidth', 2);   % Insulin
    end
    xlabel('Days'); ylabel('mU/L');
    legend(arrayfun(@(x)sprintf('c0=%.1f%%',100*x),c0s,'Uni',0),'Location','best');
    title('Insulin (PK/PD OFF)');
    grid off;
    
    % === Panel 2: Glucose ===
    subplot(1,3,2); hold on;
    for i = 1:numel(c0s)
        plot(T{i}, Y{i}(:,18), 'LineWidth', 2);   % Blood glucose
    end
    xlabel('Days'); ylabel('g/L');
    legend(arrayfun(@(x)sprintf('c0=%.1f%%',100*x),c0s,'Uni',0),'Location','best');
    title('Glucose (PK/PD OFF)');
    grid off;

   % === Panel 3: DMI & Milk ===
    % === Panel 3: DMI & Milk ===
    subplot(1,3,3); hold on;
    
    yyaxis left
    plot(at_week_days, dm_week/1000, 'LineWidth', 2, 'Color', [0 0.45 0.74]);
    ylabel('DMI (kg/day)');
    
    yyaxis right
    plot(at_week_days, mlk_week, 'LineWidth', 2, ...
         'Color', [0.85 0.33 0.10], 'LineStyle', '--');
    ylabel('Milk yield (L/day)');    
    
    xlabel('Days');
    title('Intake & Milk yield');
    legend({'DMI (kg/day)','Milk (L/day)'}, 'Location','best');
    grid off;
    % ----- (2) P4 + E2 + IGF-I -----
    vars = [10 11 16];  % P4, E2, IGF-I indices
    figure('Name','Lactating | P4 E2 IGF (PK/PD OFF)'); clf;

    for k=1:3
        subplot(3,1,k); hold on;
        for i=1:numel(c0s)
            plot(T{i}, Y{i}(:,vars(k)), 'LineWidth', 2);
        end
        title(Ynames{vars(k)}); 
        xlabel('Days'); ylabel('Relative Level');grid off;
        legend(arrayfun(@(x)sprintf('c0=%.1f%%',100*x),c0s,'Uni',0),'Location','bestoutside');

    end

end

% Additional plots only when PK/PD is ON
if show_pk
    % --- Plot P4, E2, IGF-I with dose time line ---
    vars=[10 11 16];    % indices: P4=10, E2=11, IGF-I=16
    figure('Name','Lactating with 0.02 mg/kg IM Dose | P4 E2 IGF'); clf;
    for k=1:3
        subplot(3,1,k); hold on;
        for i=1:numel(c0s)
            plot(T{i},Y{i}(:,vars(k)),'LineWidth',2);
        end
        xline(dose_day,'k--','Dose');   % vertical marker at dosing time
        xlabel('Days'); ylabel('Relative Level');
        title(Ynames{vars(k)}); grid off;
    end

        % ============================================================
    % Metabolic variables (PK/PD ON) – 17..22
    % ============================================================
    idx = [17 18 19 20 21 22];
    figure('Name','Lactating case with 0.02 mg/kg IM Dose | Metabolic (PK/PD ON)'); clf;

    for k = 1:6
        subplot(2,3,k); hold on;
        for i = 1:numel(c0s)
            plot(T{i}, Y{i}(:,idx(k)), 'LineWidth', 2);
        end
        xline(dose_day,'k--','Dose');
        title(Ynames{idx(k)});
        xlabel('Days'); 
        grid off;
    end 
end
end

function plot_all_requested(tag,T0,Y0,T1,Y1,Ynames,~,~)
% Standard endocrine plot overlay for PK/PD OFF vs ON.
% Plots only states 1..16.
% NOTE: 'par' and 'show_pk' are passed but not used here (kept for interface consistency).
figure('Name',[tag ' | Hormones']); clf;
for k=1:16
    subplot(4,4,k); plot(T0,Y0(:,k),T1,Y1(:,k),'LineWidth',1.6);
    title(Ynames{k}); grid off;
end
end

function plot_single(tag,T,Y,Ynames)
% Plot only endocrine/reproductive states 1..16 for a single simulation.
figure('Name',[tag ' | Hormones']); clf;
for k=1:16
    subplot(4,4,k); plot(T,Y(:,k),'LineWidth',1.6);
    title(Ynames{k}); grid off;
end
end

function dm=build_dmi_acute(t,STEP,dm1,f,dm2,s,l)
% Acute DMI profile generator:
%   dm1 = baseline intake
%   f   = reduction factor during restriction window
%   dm2 = intake after restriction
%   s   = start time (days)
%   l   = duration (days)
%
% The output dm is a column vector aligned with time vector t.
dm=dm1*ones(numel(t),1); 
i1=round(s/STEP)+1; 
i2=round((s+l)/STEP)+1;
dm(i1:i2)=f*dm1; 
dm(i2:end)=dm2;
end

function dm=build_dmi_chronic(t,STEP,dm1,f,dm2,s,l,r)
% Chronic DMI profile generator:
%   dm1 = baseline intake
%   f   = restriction factor during restriction window
%   dm2 = intake during refeeding window
%   s   = start time (days)
%   l   = restriction duration (days)
%   r   = refeeding duration (days)
%
% Timeline:
%   [0..s] baseline dm1
%   [s..s+l] restriction f*dm1
%   [s+l..s+l+r] refeeding dm2
%   [after] baseline dm1
dm=dm1*ones(numel(t),1);
i1=round(s/STEP)+1; 
i2=round((s+l)/STEP)+1; 
i3=round((s+l+r)/STEP)+1;
dm(i1:i2)=f*dm1; 
dm(i2:i3)=dm2; 
dm(i3:end)=dm1;
end
