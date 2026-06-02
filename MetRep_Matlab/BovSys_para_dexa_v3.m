function par= BovSys_para_dexa_v3() 

par=zeros(98,1);

%% Estrous cycle model parameters

  %-------------GnRH------------
    par( 1) = 16;      % gnrh_hypo_max 
    par( 2) = 2.75;    % gnrh_syn
    par( 3) = 2.05;    % h(hm_e2_1+hm_p4_1-hm_e2_1*hm_p4_1) 
    par( 4) = 0.0972;   % T(hm_e2_1)  
    par( 5) = 0.35;     % T(hm_p4_1)
    par( 6) = 1.91;    % h(hm_p4_2)   
    par( 7) = 0.252;    % T(hm_p4_2) 
    par( 8) = 0.99;     % h(hp_e2_1) 
    par( 9) = 0.648;    % T(hp_e2_1) 
    par(10) = 1.63;    % gnrh_clear 
    
    
    %  ----------- fsh ----------- 
    par(11) = 4.21;     % h(hm_ih_2)
    par(12) = 0.118;    % T(hm_ih_2)
    par(13) = 0.293;    % h(hp_p4_2)
    par(14) = 0.152;    % T(hp_p4_2)
    par(15) = 0.396;    % h(hm_e2_2)
    par(16) = 0.312;      % T(hm_e2_2)
    par(17) = 1.23;     % h(hp_gnrh_2)
    par(18) = 0.0708;    % T(hp_gnrh_2)
    par(19) = 2.73;      % clear(fsh_blood) 
    par(20) = 0.948;    % basal FSH release
    % ------------- lh -----------
    par(21) = 0.376;   % h(hp_e2)
    par(22) = 0.243;   % T(hp_e2)
    par(23) = 2.71;    % h(hn_p4)
    par(24) = 0.0269;  % T(hn_p4)
    par(25) = 2.22;    % h(hp_gnrh)
    par(26) = 0.69;    % T(hp_gnrh)  
    par(27) = 0.0141;  % basal LH rel
    par(28) = 2.0;    % clearance rate of LH blood
    
    
    % ------------ follicles --------------    
    par(29) = 0.562;   % h(hp_fsh_mod)
    par(30) = 1.497;    % T(hp_fsh_mod)
    par(31) = 0.322;    % T(hm_foll)                     This change number of the waves
    
    
    par(32) = 1.1;     % h(hp_p4)
    par(33) = 0.126;   % T(hp_p4)
    par(34) = 3.49;    % h(hp_lh)
    par(35) = 0.1;   % T(hp_lh)
    
    
    % -------------------- prostaglandin ------------------
    par(36) = 53.91;   % h enz-pg
    par(37) = 1.43;    % T
    par(38) = 1.087;   % T otr
    par(39) = 1.23;    % cl pg   
    
    
   % -------------------- corpus luteum ----------------   
    par(40) = 0.4;    % SF
    par(41) = 0.0335; % h(hp_cl)
    par(42) = 0.2807;    % T(hp_cl)  
    par(43) = 41.39;  % h_iof
    par(44) = 1.82;   % T_iof
    
    
   % ---------------------- progesterone  ---------------   
    par(45) = 2.25;   % linear factor cl to p4
    par(46) = 1.41;   % decay rate of progesterone  
    par(47) = 0.1;   % baseline of progesterone  

    
    % -------------------- estradiol -------------------------
    par(48) = 2.19;  % growth rate e2 dep. on foll
    par(49) = 1.23;  % decay rate e2
    
    
    % -------------------- inhibin ------------------
    par(50) = 1.41;  % rate of growth of ih(mult. with y_foll_tau)
    par(51) = 0.475; % decay rate  
    
    
    % ------------------ enzymes ------------------------  
    par(52) = 3.58;  % h enz
    par(53) = 0.77;  % T 
    par(54) = 2.98;  %cl

   
    % ------------------- OTR ------------------
    par(55) = 1.59;   % h(hp_e2&&cl)
    par(56) = 0.143;  % T(hp_e2_4)   
    par(57) = 0.644;  % clear(OTR)
    par(58) = 1.5;
    par(59) = 10;
    par(60) = 0.0007;
    
    % ----------- inter-ovarian factors ---------------
    par(61) = 39.68;  % h(hp_pg2)
    par(62) = 1.22;   % T(hp_pg2)
    par(63) = 0.6;    % T(hp_cl_3a) self destroy
    par(64) = 0.298;  % cl_iofz

    %-------------Insulin-like growth factor I------------
    par(65) = 2.5;  % h(hp_lh_2)*h(hm_p4_3)
    par(66) = 0.3;  % T(hm_p4_3)
    par(67) = 1.7;   % clearance rate
    par(68) = 0.4;   % basal IGF_Bld release  
    par(69) = 1;     % IGF
    par(70) = 0.5;     % IGF threshold 


    %-------------- Effect of Insulin on FSH & LH & IGF------------------
    par(71) = 3.;    % effect of insulin on FSH   % 2.8
    par(72) = 15.;   % threshold 
    par(73) = 1.05;  % effect of insulin on LH
    par(74) = 16.;   % threshold 
    par(75)  = 0.4;  % effect of insulin on IGF
    par(76)  = 15.;  % threshold


   
%% Metabolic model parameters
% ---------------------------

    par(77)= 0.08;         % fraction of glucose goes directly to the blood
    par(78)= 8;       % glucose transfer from blood to storage (maximal rate)
    par(79)= 10;          % Threshold
    par(80)= 1350;         % glucose transfer from storage to blood (maximal rate)
    par(81)= 10;          % Threshold
    par(82)= 180;         % Threshold
    par(83)= 0.45;          % Threshold
    
    par(84)=70182;
    par(85)=350.87;
    par(86)=0.5;
    
    par(87)= 0.0684;           % glucose release from the liver (constant rate)
    par(88)= 50;          % Threshold
    
    par(89)=10;           
    par(90)=1000;         % Storage threshold
    par(91)=0.5;
    
    par(92)= 0.5;     % scaling factor
    par(93)=1000;
    par(94)= 5;          % metabolism in the liver (constant rate)
    par(95)= 72;           % Glucose in Milk: 72 g of glucose needed to produce 1 L Milk
    par(96)= 84211;        % insulin synthesis rate
    par(97)= 2105;          % Insulin clearance rate
    par(98)= 15e4;       % transfer of body fat to glucose in the blood (maximal rate)


    









end

   


