function [T,Y,OvT]=BovSys_run_v3_baseline(par, frac2, day, daystr, RefP4, str, STEP)
% BovSys_run_v3_baseline
%
% Non-interactive v3 baseline runner for Bayesian experimental design.
%
% This wrapper lets the BED workflow use the published v3 model equations
% while keeping Dexa PK/PD switched off. It provides the non-interactive
% BED-facing interface:
%
%     [T,Y,OvT] = runner(par, frac2, day, daystr, RefP4, str, STEP)
%
% Inputs:
%   par     : v3 parameter vector from BovSys_para_dexa_v3()
%   frac2   : baseline c0/feed-glucose scaling used by the model
%   day     : simulation horizon in days
%   daystr  : number of final days used for P4 pattern comparison
%   RefP4   : reference P4 segment used to screen cyclic simulations
%   str     : 'l' lactating or 'n' non-lactating
%   STEP    : output time step in days
%
% Output:
%   T       : solver time vector
%   Y       : state trajectories, using the 25-state v3 layout
%   OvT     : inferred ovulation time, or NaN if the simulated cycle fails
%             the same P4-pattern screen used by the original BED script
%
% Note:
%   Dexa is disabled here. This is a baseline-model BED runner, not a Dexa
%   perturbation experiment.

this_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(this_dir,'..','..','..');
model_dir = fullfile(repo_root,'MetRep_Matlab');
addpath(model_dir);

y0 = [ ...
    0.667; 0.551; 0.316; 0.395; 1.000; 0.642; 1.000; 0.00506; 0.0;   0.004; ...
    0.89;  0.826; 0;     0.0183;0.35;  0.48;  15.5;  0.48;   15e4;  110; ...
    535;   105.;  0;     0;     0];

tspan = 0:STEP:day;
at = tspan;
c0 = frac2;
use_dexa = false;

if str=='l'
    y0(14)=2.5;
    load(fullfile(model_dir,'DM.mat'),'DMdataset');
    load(fullfile(model_dir,'ML.mat'),'MLKdataset');
    dm_week = DMdataset(:,1);
    mlk_week = MLKdataset(:,2);
    wk_time = MLKdataset(:,1);
    at_week_days = wk_time * 7;
    dm = interp1(at_week_days, dm_week, tspan,'linear','extrap');
    mlk = interp1(at_week_days, mlk_week, tspan,'linear','extrap');
else
    dm = 11700 * ones(size(tspan));
    mlk = zeros(size(tspan));
end

odeopt = odeset('nonnegative',1:25);
[T,Y] = ode15s(@BovSys_Equa_dexa_v3, tspan, y0, odeopt, par, str, at, dm, mlk, c0, use_dexa);

OvT = nan;

if numel(RefP4)<=1
    return;
end

start_idx = max(1, floor((day-daystr)*(1/STEP)));
end_idx = min(size(Y,1), floor(day*(1/STEP)));
g = Y(start_idx:end_idx,10);
f = RefP4(:);

if isempty(g) || isempty(f) || norm(f)==0 || norm(g)==0
    return;
end

if str=='l'
    LAG=90;
else
    LAG=0;
end

[acor,~] = xcorr(f,g,LAG);
croc = acor/(norm(f)*norm(g));

NAD = abs(sum(f) - sum(g))/sum(f);
NSD = abs(norm(f)^2 - norm(g)^2)/norm(f)^2;

if str=='l'
    CNL=(any(croc > 0.5) && NAD < 0.16 && NSD < 0.3);
else
    CNL=(any(croc > 0.8) && NAD < 0.16 && NSD < 0.3);
end

if CNL
    [pks,locs] = findpeaks(Y(:,10));
    TPeak=T(locs(find(pks>0.9)));
    if ~isempty(TPeak)
        OvT = TPeak(1) - 14;
    end
end

end
