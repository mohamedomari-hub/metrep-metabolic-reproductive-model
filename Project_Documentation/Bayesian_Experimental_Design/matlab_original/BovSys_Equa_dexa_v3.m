function f = BovSys_Equa_dexa_v3(t, y, par, str, at, dm, mlk, c0, use_dexa)
% BovSys_Equa_dexa_v3
% =========================================================================
% This function defines the RIGHT-HAND SIDE (RHS) of the ODE system:
%     dy/dt = f(t, y; parameters, inputs)
%
% It combines:
%   (1) Estrous-cycle endocrine model (states 1..16)
%   (2) Metabolic model (states 17..22)
%   (3) Dexamethasone PK (IM depot + central) and PD effect-site (23..25)
%
% The solver (ode15s) repeatedly calls this function to compute derivatives
% at many time points. The output f is a 25x1 vector of time derivatives.
% =========================================================================
%
% Inputs:
%   t         : current simulation time (days)
%   y         : current state vector (25x1)
%   par       : parameter vector (all physiological + PK/PD constants)
%   str       : mode flag (e.g., 'l' lactating, 'n' non-lactating)
%   at        : time grid for external inputs dm/mlk (days)
%   dm        : dry matter intake values over time grid at (same length as at)
%   mlk       : milk production values over time grid at
%   c0        : scaling factor for glucose pool derived from feed intake
%   use_dexa  : boolean flag (true/false) toggling Dexa PK/PD ON/OFF
%
% Output:
%   f         : derivatives dy/dt (25x1)
%
% State layout (n=25):
%   1..16  Estrous-cycle hormones/mediators
%   17..22 Metabolic states
%   23     Dexa depot amount (ng)      [IM injection site]
%   24     Dexa central amount (ng)    [systemic circulation]
%   25     Dexa effect-site conc (ng/mL) [PD link compartment]
% =========================================================================

% Pre-allocate derivative vector 
% This avoids dynamic resizing and is standard practice for ODE RHS functions.
f = zeros(25,1);

%% =========================================================================
%% Unpack states
%% =========================================================================
% For readability, each element of y is assigned a descriptive variable name.
% --- Estrous-cycle endocrine states (1..16) ---
y_gnrhh  = y(1);   % GnRH "hypothalamic store" (or hypothalamic pool)
y_gnrh   = y(2);   % GnRH in portal blood / pituitary exposure
y_fshp   = y(3);   % pituitary FSH pool (releasable)
y_fsh    = y(4);   % circulating (blood) FSH
y_lhp    = y(5);   % pituitary LH pool (releasable)
y_lh     = y(6);   % circulating (blood) LH
y_fol    = y(7);   % follicle functional size/activity variable
y_pg     = y(8);   % PGF2α proxy variable
y_cl     = y(9);   % corpus luteum functional mass/activity
y_p4     = y(10);  % progesterone (P4)
y_e2     = y(11);  % estradiol (E2)
y_ih     = y(12);  % inhibin (INH)
y_enz    = y(13);  % enzyme proxy controlling PGF production
y_ot1    = y(14);  % oxytocin (OXT)
y_iof    = y(15);  % IOF (Intra-ovarian factors)
y_igf_b  = y(16);  % IGF-1 in blood (or systemic IGF proxy)

% --- Metabolic states (17..22) ---
y_ins    = y(17);  % insulin
y_glu    = y(18);  % blood glucose
y_fat    = y(19);  % fat reserve pool (energy store)
y_lv     = y(20);  % liver glucose pool / hepatic glucose availability
y_glu_st = y(21);  % glucose storage pool (glycogen-like)
y_gluca  = y(22);  % glucagon

% --- Dexa PK/PD states (23..25) ---
Adep     = y(23);  % IM depot amount of Dexa (ng)
Acent    = y(24);  % central (systemic) amount of Dexa (ng)
Ce       = y(25);  % effect-site concentration of Dexa (ng/mL)

%% =========================================================================
%% Inputs (time-varying exogenous drivers)
%% =========================================================================
% DMI and Milk are defined as time series on grid "at".
% We interpolate them to the exact solver time "t" using linear interpolation.
% 'extrap' ensures values are available even if solver queries outside range.
DMI  = interp1(at, dm,  t, 'linear', 'extrap');   % Dry matter intake at time t
Milk = interp1(at, mlk, t, 'linear', 'extrap');   % Milk yield at time t

%% ========================================================================
% Estrous Cycle Model (UNCHANGED except tiny consistency edits)
% ========================================================================
% This block defines dynamics of GnRH, FSH, LH, follicle, CL, P4, E2, INH,
% and mediators (ENZ, OXT, IOF, PGF, IGF).
%
% It uses Hill-type activation/inhibition terms:
%   - "hp_*" : Hill-positive (activation)
%   - "hm_*" : Hill-minus (inhibition)
% and mass-balance style equations:
%   synthesis - release/clearance, etc.

% -------------------------------------------------------------------------
% GnRH module
% -------------------------------------------------------------------------
% hm_e2_1: inhibitory modulation by E2 (form: K^2/(K^2 + E2^2))
% hm_p4_1: inhibitory modulation by P4
% hm_p4_2: an additional P4-related component (scaled) contributing to release
% hp_e2_1: strong positive E2 effect (power 5) controlling GnRH "pulse" transmission
hm_e2_1 = par(4)^2 / (par(4)^2 + y_e2^2);
hm_p4_1 = par(5)^2 / (par(5)^2 + y_p4^2);
hm_p4_2 = par(6) * par(7)^2 / (par(7)^2 + y_p4^2);
hp_e2_1 = par(8) * y_e2^5 / (par(9)^5 + y_e2^5);

% gnrh_syn: replenishment of hypothalamic GnRH store towards capacity par(1)
%           (1 - y_gnrhh/par(1)) gives saturation (no growth at capacity)
gnrh_syn = par(2) * (1.0 - y_gnrhh / par(1));

% gnrh_rel: effective release rate from hypothalamic store
% The term (hm_p4_1 + hm_e2_1 - hm_e2_1*hm_p4_1) acts like OR-combination of inhibitors
% plus an additional P4-dependent component hm_p4_2.
gnrh_rel = par(3) * (hm_p4_1 + hm_e2_1 - hm_e2_1 * hm_p4_1) + hm_p4_2;

% f(1): d/dt of hypothalamic store = synthesis - release*store
f(1) = gnrh_syn - gnrh_rel * y_gnrhh;

% f(2): d/dt of circulating/portal GnRH = release*store*E2-dependent gating - clearance
f(2) = gnrh_rel * y_gnrhh * hp_e2_1 - par(10) * y_gnrh;

% -------------------------------------------------------------------------
% FSH module + insulin coupling
% -------------------------------------------------------------------------
% efsh: Hill exponent controlling insulin -> FSH stimulation steepness
efsh = 2;

% hp_ins_FSH: insulin-dependent positive modulation of FSH synthesis
% Form: max_gain * Ins^e / (Ins^e + K^e)
hp_ins_FSH = par(71) * (y_ins^efsh / (y_ins^efsh + par(72)^efsh));

% Additional endocrine modulators:
% hp_p4_1   : P4 stimulates (?) release-related term
% hm_e2_2   : E2 inhibition (or modulation) term
% hp_gnrh_1 : GnRH activation term (Michaelis-Menten form)
hp_p4_1   = par(13) * y_p4^2    / (par(14)^2 + y_p4^2);
hm_e2_2   = par(15) * par(16)^2 / (par(16)^2 + y_e2^2);
hp_gnrh_1 = par(17) * y_gnrh    / (par(18) + y_gnrh);

% fsh_syn: baseline pituitary synthesis term inhibited by inhibin (power 5)
%          (par(12)^5)/(par(12)^5 + INH^5) decreases synthesis when INH high.
fsh_syn = par(11) * par(12)^5 / (par(12)^5 + y_ih^5);

% fsh_rel: effective release flux from pituitary pool to blood (depends on modulators)
fsh_rel = (par(20) + hp_p4_1 + hm_e2_2 + hp_gnrh_1) * y_fshp;

% f(3): pituitary pool dynamics = synthesis*(insulin modulation) - release
f(3) = fsh_syn * hp_ins_FSH - fsh_rel;

% f(4): blood FSH dynamics = release - clearance
f(4) = fsh_rel - par(19) * y_fsh;

% -------------------------------------------------------------------------
% LH module + insulin coupling
% -------------------------------------------------------------------------
% elh: Hill exponent controlling insulin -> LH stimulation
elh = 3;

% hp_ins_LH: insulin-dependent positive modulation of LH synthesis
hp_ins_LH = par(73) * (y_ins^elh / (y_ins^elh + par(74)^elh));

% hp_e2_2   : E2-dependent stimulation
% hm_p4_3   : P4-dependent inhibition/modulation
% hp_gnrh_2 : GnRH-dependent activation (very steep power 5)
hp_e2_2   = par(21) * y_e2^2    / (par(22)^2 + y_e2^2);
hm_p4_3   = par(23) * par(24)^2 / (par(24)^2 + y_p4^2);
hp_gnrh_2 = par(25) * y_gnrh^5  / (par(26)^5 + y_gnrh^5);

% lh_syn: combined E2 and P4 modulation (additive)
lh_syn = (hp_e2_2 + hm_p4_3);

% lh_rel: release from pituitary pool, boosted by GnRH
lh_rel = (par(27) + hp_gnrh_2) * y_lhp;

% f(5): pituitary LH pool dynamics = synthesis*(insulin modulation) - release
f(5) = lh_syn * hp_ins_LH - lh_rel;

% f(6): blood LH = release - clearance
f(6) = lh_rel - par(28) * y_lh;

% -------------------------------------------------------------------------
% Follicle dynamics + IGF coupling
% -------------------------------------------------------------------------
% hm_IGF: inhibitory term driven by IGF (note: form K^2/(K^2 + IGF^2))
%         This modifies how LH acts on follicle (through hp_lh)
hm_IGF = par(69) * (par(70)^2) / (par(70)^2 + y_igf_b^2);

% hp_lh: LH effect on follicle, reduced by hm_IGF (IGF-dependent scaling of LH sensitivity)
hp_lh  = par(34) * y_lh^2 / (hm_IGF^2 + y_lh^2);

% hp_p4_2: P4-related inhibitory (or decay) effect on follicle (power 5)
hp_p4_2 = par(32) * y_p4^5 / (par(33)^5 + y_p4^5);

% hm_foll: self-limiting term by follicle size (K^2/(K^2+fol^2))
hm_foll = (par(31)^2) / (par(31)^2 + y_fol^2);

% hp_fsh_mod: FSH-driven follicle growth, modulated by hm_foll (FSH sensitivity reduced at large fol)
hp_fsh_mod = par(29) * y_fsh^2 / ((par(30) * hm_foll)^2 + y_fsh^2);

% f(7): follicle growth - (loss terms)*follicle
f(7) = hp_fsh_mod - (hp_p4_2 + hp_lh) * y_fol;

% -------------------------------------------------------------------------
% PGF2α proxy dynamics
% -------------------------------------------------------------------------
% hp_enz: ENZ activates PGF production
% hp_ot : OXT gate (very steep, power 10) produces a pulse-like response
hp_enz = par(36) * y_enz / (par(37) + y_enz);
hp_ot  = y_ot1^10 / (par(38)^10 + y_ot1^10);

% f(8): production - clearance
f(8) = hp_enz * hp_ot - par(39) * y_pg;

% -------------------------------------------------------------------------
% Corpus luteum (CL) dynamics
% -------------------------------------------------------------------------
% hp_cl_1: CL auto-support term (very steep power 30) ~ switch-like
% hp_iof : IOF-driven luteolysis
hp_cl_1 = par(41) * y_cl^30 / (y_cl^30 + par(42)^30);
hp_iof  = par(43) * y_iof^5 / (par(44)^5 + y_iof^5);

% f(9): CL gains from LH-driven follicle-to-CL formation + self term - IOF-induced loss
f(9) = par(40) * hp_lh * y_fol + hp_cl_1 - hp_iof * y_cl;

% -------------------------------------------------------------------------
% Progesterone (P4) dynamics
% -------------------------------------------------------------------------
% Basal production par(47) + CL-dependent production - clearance
f(10) = par(47) + par(45) * y_cl^2 - par(46) * y_p4;

% -------------------------------------------------------------------------
% Estradiol (E2) dynamics
% -------------------------------------------------------------------------
% Follicle produces E2, and E2 is cleared linearly
f(11) = par(48) * y_fol^2 - par(49) * y_e2;

% -------------------------------------------------------------------------
% Inhibin (INH) dynamics
% -------------------------------------------------------------------------
% Follicle produces inhibin, and inhibin is cleared linearly
f(12) = par(50) * y_fol^2 - par(51) * y_ih;

% -------------------------------------------------------------------------
% ENZ dynamics
% -------------------------------------------------------------------------
% ENZ is stimulated by P4 and cleared linearly
hp_p4_3 = par(52) * y_p4 / (par(53) + y_p4);
f(13) = hp_p4_3 - par(54) * y_enz;

% -------------------------------------------------------------------------
% Oxytocin (OXT) dynamics
% -------------------------------------------------------------------------
% hp_e2_3: E2-dependent OXT stimulation
hp_e2_3 = par(55) * y_e2^2 / (par(56)^2 + y_e2^2);

% input_oxt_L: optional lactation input pulse (only for lactating animals)
% This adds an external transient source term when str == 'l'
if str == 'l'
    input_oxt_L = par(58) * exp(-par(60) * (t)^2);
else
    input_oxt_L = 0;
end

% f(14): OXT dynamics = external pulse + E2/CL-driven production - clearance
f(14) = input_oxt_L + hp_e2_3 * y_cl^2 - par(57) * y_ot1;

% -------------------------------------------------------------------------
% IOF dynamics (intermediate luteolysis factor)
% -------------------------------------------------------------------------
% hp_pg: PGF drives IOF strongly (power 10)
% hp_cl_3: CL-dependent gate for IOF formation
hp_pg   = par(61) * y_pg^10 / (y_pg^10 + par(62)^10);
hp_cl_3 = y_cl / (y_cl + par(63));

% f(15): IOF production - clearance
f(15) = hp_pg * hp_cl_3 - par(64) * y_iof;

% -------------------------------------------------------------------------
% IGF dynamics with insulin coupling
% -------------------------------------------------------------------------
% hp_ins_IGF: insulin-dependent stimulation of IGF
hp_ins_IGF = par(75) * (y_ins^1 / (y_ins^1 + par(76)^1));

% hm_p4_3: P4-dependent modulation term (power 4), reducing IGF drive
hm_p4_3 = par(65) * par(66)^4 / (par(66)^4 + y_p4^4);

% f(16): IGF = basal + (P4-modulated)*(insulin drive) - clearance
f(16) = par(68) + hm_p4_3 * hp_ins_IGF - par(67) * y_igf_b;

%% ========================================================================
% Dexa PK/PD block with TRUE ON/OFF behavior
% ========================================================================
% This block toggles dexamethasone pharmacokinetics (PK) and pharmacodynamics (PD).
%
% If use_dexa == false:
%   - Dexa compartments (Adep, Acent, Ce) are "frozen" (no dynamics)
%   - PD multipliers are set to 1 (no effect on metabolism)
%
% If use_dexa == true:
%   - Adep -> Acent via first-order absorption (ka)
%   - Acent eliminated first-order (ke)
%   - Ce follows central concentration C via effect-site equilibration (keo)
%   - PD:
%       Effect_gluca > 1  increases glucagon secretion
%       Effect_bt    < 1  reduces glucose utilization/uptake terms
if ~use_dexa
    % Freeze Dexa states and neutralize PD
    f(23) = 0;     % dAdep/dt = 0 (no depot change)
    f(24) = 0;     % dAcent/dt = 0 (no central change)
    f(25) = 0;     % dCe/dt = 0 (no effect-site change)

    % PD multipliers set to unity => metabolism behaves like "no Dexa"
    Effect_gluca = 1;
    Effect_bt    = 1;
else
    % ---------------------------------------------------------------------
    % Convert body weight to distribution volume in mL
    % Vd_mL is used to convert amount (ng) to concentration (ng/mL)
    % ---------------------------------------------------------------------
    BW_kg = 600;
    Vd_mL = 1.105 * BW_kg * 1e3;

    % PK parameters:
    % ka  : absorption rate constant from depot to central
    % ke  : elimination rate constant from central
    % F   : bioavailability fraction from depot to central
    % keo : effect-site equilibration rate constant
    ka  = 13.4352;
    ke  = 2.7086;
    F   = 0.72;
    keo = 0.7;

    % Central concentration derived from amount and volume
    C = Acent / Vd_mL;

    % Dexa depot amount decays by absorption
    f(23) = -ka * Adep;

    % Central amount increases from absorbed depot (scaled by F), decreases by elimination
    f(24) =  F * ka * Adep - ke * Acent;

    % Effect-site concentration tracks central concentration (Ce -> C) with keo
    f(25) =  keo * (C - Ce);

    % ---------------------------------------------------------------------
    % PD (Emax-style) functions:
    % Effect_gluca: increases glucagon secretion (Emax up-regulation)
    % Effect_bt   : decreases "uptake/utilization" terms (inhibitory effect)
    % ---------------------------------------------------------------------
    Emax_PD = 3;   % maximum fractional stimulation for glucagon secretion (above baseline)
    Ca      = 1.8; % potency parameter for glucagon PD curve (Ce producing half effect, in Hill form)
    Cb      = 1.8; % potency parameter for "bt" inhibition curve

    % Glucagon stimulation multiplier: 1 + Emax * Ce^10/(Ce^10 + Ca^10)
    Effect_gluca = 1 + Emax_PD * (Ce^10) / (Ce^10 + Ca^10);

    % "bt" inhibition multiplier: 1 - Ce^7/(Ce^7 + Cb^7)
    % This reduces certain metabolic fluxes when Ce increases.
    Effect_bt    = 1 - (Ce^7) / (Ce^7 + Cb^7);
end

%% ========================================================================
% Metabolic model (apply PD consistently)
% ========================================================================
% This block models glucose/insulin/glucagon/fat/liver/storage dynamics.
%
% Key concepts:
%   - Feed intake (DMI) creates a "glucose pool" that is split into:
%       * glu_feed_gng : flow into liver / gluconeogenic/absorptive supply
%       * glu_feed_bl  : direct contribution to blood glucose
%   - Insulin secretion increases with blood glucose; insulin degrades linearly.
%   - Glucagon secretion increases when glucose is low, and is multiplied by Dexa PD.
%   - Liver produces glucose depending on liver state and glucagon.
%   - Several exchange fluxes connect blood, liver, storage, and fat.
%   - Dexa PD multiplier Effect_bt inhibits selected insulin-dependent uptake/utilization terms.

% -------------------------------------------------------------------------
% Feed-derived glucose pool
% -------------------------------------------------------------------------
% glu_pool: glucose availability derived from intake (scaled by c0)
glu_pool = c0 * DMI;

% par(77) splits glucose pool between blood and liver-directed pathways
glu_feed_gng = (1 - par(77)) * glu_pool;  % portion routed to liver supply
glu_feed_bl  = par(77) * glu_pool;        % portion routed directly to blood

% -------------------------------------------------------------------------
% Insulin dynamics (secretion and degradation)
% -------------------------------------------------------------------------
% ins_sec: glucose-stimulated insulin secretion (very steep power 10)
ins_sec = par(96) * (y_glu^10 / (y_glu^10 + par(91)^10));

% ins_deg: first-order insulin degradation
ins_deg = par(97) * y_ins;

% -------------------------------------------------------------------------
% Glucagon dynamics (secretion and degradation) + Dexa stimulation
% -------------------------------------------------------------------------
% gluca_sec: glucose-inhibited secretion (higher when glucose is low),
% multiplied by Effect_gluca (Dexa-stimulated when Dexa ON)
gluca_sec = par(84) * (par(86)^2 / (y_glu^2 + par(86)^2)) * Effect_gluca;

% gluca_deg: first-order glucagon degradation
gluca_deg = par(85) * y_gluca;

% -------------------------------------------------------------------------
% Hepatic glucose production
% -------------------------------------------------------------------------
% glu_prod: glucagon-driven liver glucose production (scaled by liver state y_lv)
glu_prod = par(87) * y_lv * y_gluca;

% -------------------------------------------------------------------------
% Blood-to-liver flux (insulin-dependent)
% -------------------------------------------------------------------------
% glu_bl_lv: insulin-dependent transport/uptake from blood to liver
% NOTE: comment in code: currently NOT inhibited by Dexa (no Effect_bt here).
glu_bl_lv = par(88) * (y_glu^10/(y_glu^10 + par(83)^10)) * y_ins;

% -------------------------------------------------------------------------
% Liver-to-storage flux (glycogen-like storage formation)
% -------------------------------------------------------------------------
% glu_lv_st: storage formation depends on:
%   - remaining storage capacity (1 - y_glu_st/1000)
%   - milk dependence gate (par(81)^2/(par(81)^2 + Milk^2))
%   - liver state y_lv and insulin y_ins
%   - inhibited by Dexa via Effect_bt
glu_lv_st = par(82) * (1 - (y_glu_st/1000)) ...
          * (par(81)^2/(par(81)^2 + Milk^2)) ...
          * y_lv * y_ins ...
          * Effect_bt;

% -------------------------------------------------------------------------
% Storage-to-liver mobilization (glycogen breakdown proxy)
% -------------------------------------------------------------------------
% glu_st_lv: glucagon-driven mobilization from storage to liver
glu_st_lv = par(80) * y_gluca * (y_glu_st^10 / (y_glu_st^10 + par(79)^10));

% -------------------------------------------------------------------------
% Fat-related fluxes
% -------------------------------------------------------------------------
% glu_fat_lv: fat mobilization contributing to liver glucose availability,
% driven by glucagon and modulated by storage and fat saturation terms
glu_fat_lv = par(78) ...
    * (par(59)^10 / (y_glu_st^10 + par(59)^10)) ...
    * (y_fat / (y_fat + par(98))) ...
    * y_gluca;

% glu_lv_fat: conversion of liver glucose toward fat (insulin-dependent),
% inhibited by Dexa via Effect_bt
glu_lv_fat = par(35) * (par(89)^1/(par(89)^1 + Milk^1)) ...
           * (y_glu_st^10/(y_glu_st^10 + par(90)^10)) ...
           * y_lv * y_ins ...
           * Effect_bt;

% -------------------------------------------------------------------------
% Blood glucose usage (peripheral utilization) + milk component
% -------------------------------------------------------------------------
% The code intends Effect_bt to inhibit the whole utilization term.
glu_bl_usage =  par(93) * (y_glu^10 / (y_glu^10 + par(92)^10)) + par(95) * Milk  * Effect_bt;

% Liver usage (baseline hepatic utilization)
glu_lv_usage = par(94) * y_lv;

% -------------------------------------------------------------------------
% Unit conversion between fat and glucose equivalents
% -------------------------------------------------------------------------
% fat_to_gl: converts liver->fat flux to fat pool change (energy equivalence)
% gl_to_fat: inverse conversion for fat->glucose related flux
fat_to_gl = 9.3 / 4.1;
gl_to_fat = 1 / fat_to_gl;

% -------------------------------------------------------------------------
% Blood volume to convert net blood glucose flux into concentration change
% -------------------------------------------------------------------------
Vblood_L = 22.8;

% -------------------------------------------------------------------------
% Final metabolic state derivatives (17..22)
% -------------------------------------------------------------------------
% f(17): insulin = secretion - degradation
f(17) = ins_sec - ins_deg;

% f(18): blood glucose = (inputs + production - transfers - usage) / blood volume
f(18) = (glu_feed_bl + glu_prod - glu_bl_lv - glu_bl_usage) / Vblood_L;

% f(19): fat pool = (liver->fat converted) - (fat->liver mobilization)
f(19) = fat_to_gl * glu_lv_fat - glu_fat_lv;

% f(20): liver glucose pool balance:
%   + feed to liver
%   - hepatic glucose production (export)
%   + storage mobilization - storage formation
%   + fat conversion contributions
%   + blood->liver transfer
%   - liver usage
f(20) = glu_feed_gng - glu_prod + glu_st_lv - glu_lv_st ...
        + gl_to_fat * glu_fat_lv - glu_lv_fat ...
        + glu_bl_lv - glu_lv_usage;

% f(21): storage pool = formation - mobilization
f(21) = glu_lv_st - glu_st_lv;

% f(22): glucagon = secretion - degradation
f(22) = gluca_sec - gluca_deg;

end