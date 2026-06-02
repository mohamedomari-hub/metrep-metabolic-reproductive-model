
% BED_1M_ALL
%
% Clean GitHub-facing Bayesian experimental design script for the BovSys/MetRep
% model using the published v3 baseline model interface.
%
% Important:
%   - this script uses the published v3 model equations through
%     BovSys_run_v3_baseline();
%   - Dexa PK/PD is switched off in the runner, so this is a baseline-model
%     BED workflow rather than a Dexa perturbation workflow;
%   - the original historical n2 BED script was moved out of this folder to
%     avoid confusing two model versions in the GitHub repository.
%
% Readability/port notes for GitHub curation:
%   - comments mark the biological setup, Monte Carlo simulation, mutual
%     information calculations, posterior calculations, and figure outputs;
%   - the first active workflow evaluates lactating/ovulation-time designs;
%   - the second active workflow evaluates the non-lactating KDE design for
%     a selected parameter.
%
% Main dependencies: BovSys_para_dexa_v3, BovSys_Equa_dexa_v3,
% BovSys_run_v3_baseline, Parallel Computing Toolbox for parfor, and
% Statistics and Machine Learning Toolbox functions such as unifrnd, normrnd,
% ksdensity, mvksdensity, mvnpdf, and datasample.

function BED_1M_ALL()
tic

datetime('now')   
clear all
close all

this_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(this_dir,'..','..','..');
addpath(this_dir);
addpath(fullfile(repo_root,'MetRep_Matlab'));

display('----Sampling-----')

% -------------------------------------------------------------------------
% Workflow 1: lactating/ovulation-time design
% -------------------------------------------------------------------------
% Goal: use simulated model ensembles to identify which measurement day and
% which species carry the most information about ovulation timing W.

%% State variable of the model

Ynames =  { 'GnRhH'; 'GnRhP';'FSHP'; 'FSH'; 'LHP';  'LH'; 'Foll'; 'PGF';   'CL';    'P4';  'E2'; .....
           'INH'; 'ENZ'; 'OXT';  'IOF'; 'IGF';'Ins'; 'Glu'; 'Fat';   'Glucose in Liver'; 'Stored Glucose'};

%% Defining Species to include in simulation

        Species= 1:21; % Species to simulate in Metabolism (Species 1 and 2) and Reproduction model (the rest) 
        
        % Species retained for the BED calculation:
        % FSH, PGF, P4, E2, INH, IGF, Ins, and Glu.
        Species= Species([4,8,10,11,12,16,17,18]);

        
        [~, spel]=ismember([4,8,10,11,12,16,17,18],Species);
     
     
        %% Time of sampling
        
                    % Candidate measurement days. Each value becomes a
                    % possible design option for the mutual information
                    % ranking.
                    tpoint = 7:7:63;                                                  
    
                    Comb=reshape(tpoint,length(tpoint),1);
                    
                    numComb=size(Comb,1);
                    
                    
        
        %% All possible Combination of pair-Species mutual information
        
                    % Pairs are used to test whether two measured species
                    % provide redundant or complementary information.
                    CombZ = nchoosek(1:length(spel),2);

                    numCombZ = size(CombZ,1);
                    
                    str='l';
                    
                    display('----  Calcule de information Mutuel-----')
    
        
        %% Number of samples

        % Monte Carlo ensemble size for parameter sampling and simulated
        % model outputs.
        N0= 10000;

        %% Argument to pass to the solver

        frac2=0.25; 
        STEP=1/6;
        day=100;
        daystr = 95;

                

        %% Reading Paramaters

        par=BovSys_para_dexa_v3();

        
        %% Sampling theta from Uniform dis...

        % Defining the boundary of parameters interval
        
        % Parameter uncertainty is represented by a uniform interval around
        % the reference parameter vector.
        devS = 0.1;
        
        %Sens= [12    50    92    40    48    49    29    34   101    58 ...
         %   57    38    60    81    51    79    72   100    76    88    ...
         %   74    46    47    55    41    87    98    62    35    36    33    56];
    
    
        aRep=nan(length(par),1);
        bRep=nan(length(par),1);

        %NotSens = setxor(1:length(par),Sens);

        aRep = par - par*devS;
        bRep = par + par*devS;
        
        %devNotS = 0.1;
        
        %aRep(NotSens)=par(NotSens) - par(NotSens)*devNotS;
        %bRep(NotSens)=par(NotSens) + par(NotSens)*devNotS;      
        

        %% Check if parameters are negative

        if any(aRep<0)
            aRep(find(aRep<0))=0;
        end

        %% Sampling from parameters from uniform distribution ...
        
        repa=repmat(aRep,1,N0);
        repb=repmat(bRep,1,N0);
    
        thetaRep0 = unifrnd(repa,repb);

        %% Run the standard model to calculate the error ...
        
                     % The reference simulation supplies the output scale
                     % used to define synthetic measurement error.
                     RefP4=1;
    
                     [T,YRepS,~]=BovSys_run_v3_baseline(par,frac2,day,daystr,RefP4,str,STEP);
                                                           
                     RefP4=YRepS(((day-daystr)*(1/STEP)):(day*(1/STEP)),10);

                     Out = YRepS;
          
                     % Measurement error is approximated as 10% of the mean
                     % simulated value for each selected species.
                     Abs_Error = 0.1*mean(Out(:,Species),1);
                                          
                     save('Abs_Error.mat','Abs_Error')
        
                                                figure(111)

                                                for i =[1 4]
                                                    subplot(2,2,i)
                                                    plot(T,YRepS(:,i+6),'k')
                                                    hold on
                                                end
       


        %% Run the model for different theta (smapled parameters)
        
                    % Each sampled parameter vector is propagated through
                    % the model. Outputs at candidate days are retained for
                    % the BED calculations.
                    %% First define matrix for input simulation

                    OvS0  = zeros(length(tpoint),(length(Species))) ;
                    OVmat = nan(N0,1);



tic
        
display('----Generate models -----')


       
       %colorlist={'k','r','b','g'};
      
        parfor (h=1:N0,16)
        %for h=1:N0
            
                    %% Run the model for different theta (smapled parameters)

                     [T,Y0,OvT0]=BovSys_run_v3_baseline(thetaRep0(:,h),frac2,day,daystr, RefP4,str,STEP);
                     
                     Out = Y0;
                     
                     %% check if the model return nan

                     if isnan(OvT0)

                             OvS0(:,:,h)=nan(length((tpoint)),(length(Species))) ;                       

                     else

                             OvS0(:,:,h) = Out((tpoint)*(1/STEP),Species);
                             OVmat(h,1)  = OvT0;
                   
                     end                 
                        
       end
       
 % 

       %% Select the successfully simulated samples and reject the nan 
       
                    % Failed simulations are removed from the ensemble
                    % before density and information estimates are computed.
                    [~, ~, dim]= ind2sub(size(OvS0),find(isnan(OvS0)));

                    if ~isempty(dim)
                        
                        OvS0(:,:,unique(dim))=[];
                      
                        save('OvS0.mat','OvS0')

                        display('............')
     
                        OVmat=OVmat(~isnan(OVmat));

                        save('OVmat.mat','OVmat')
                        
                        %% update the number of samples

                        N0 = N0 - length(unique(dim));

                    else
                        
                        display('...llll....')

                        save('OvS0.mat','OvS0');
                        size(OvS0)
                        
                        %save('OVmat.mat','OVmat')

        
                    end

                    AcceptedOutput=sprintf('Selected samples N0=%d ... Number of Rejected samples are: %d', N0, length(unique(dim)))

                    
toc
datetime('now')            
                    
                    
                    
      
                    
        %% defining matrices for input

        Info_OvT = nan(1,length(tpoint));
        Info_TS = nan(length(tpoint),length(spel));


        
        EntropyZ=nan(numComb,length(spel));
        EntropyW=nan(numComb,1);
        Info_ZZ = nan(numComb,numCombZ);     
        dataMI= nan(numComb,N0);

        
        
        NT=N0;
        
for Ni=[0]
                    
sprintf('-------- Ni=%d ------',Ni)
    
    numSam=[];
    % Loop over candidate sampling days. For each day, exclude simulations
    % where ovulation already occurred before the hypothetical measurement.
    for p=1:9
        
        p

         
             %% Number of samples

              %N0 = 31470;   

              %N0=5000;


             %% loading models

             data1=load('OvS0.mat');
             OvS0=data1.OvS0(:,spel,:);

             data2=load('OVmat.mat');
             OVmat=data2.OVmat;
             data3=load('Abs_Error.mat');
             Abs_Error=data3.Abs_Error(spel);


             %% For sampling measurement: Filering ovulation time before measurements time

             OVmatZ=OVmat;
             indFilter=OVmatZ<=tpoint(p);
             OVmatZ(indFilter)=[];

             OvS0Z=OvS0;
             OvS0Z(:,:,indFilter)=[];


             %% For mutual information calculations: Filering ovulation time before measurements time
             OVmatS=OVmat(1:N0);
             indFilter=OVmatS<=tpoint(p);
             OVmatS(indFilter)=[];

             OvS0=OvS0(:,:,(1+Ni):(N0+Ni));
             OvS0(:,:,indFilter)=[];


             %% Updating number of samples

             
             N0=N0-length(find(indFilter==1))
             
             numSam=[numSam N0];


                        indTT=p;

                        Info_ns=0;
                        Info_ns_ind=0;
                        EntZ=0;
                        Entw=0;


                        subComb = Comb(indTT,:);
                        [~,a]=ismember(subComb,tpoint);
                        indT=a(a~=0);


                        mvML=nan(numCombZ,1);

                        ns_ZZ=zeros(1,numCombZ);

                                % Monte Carlo estimate of the density terms
                                % needed for mutual information:
                                %   W = ovulation time;
                                %   Z = simulated noisy measurements;
                                %   I(W;Z) = log p(W,Z) - log p(W) - log p(Z).

                                for i=1:N0
                                    
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %% Calculate Joint mutual information: I(w,Z) where Z represents more than species
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                                            %% Evaluate Samples at f(w) dist using KDE    

                                            % p(W): kernel density estimate
                                            % for the sampled ovulation time.
                                            [~,~,bw] = ksdensity(OVmatS);    

                                            mv=nan(1,N0);

                                            for ii = 1:N0

                                                mv(1,ii) =  normpdf(OVmatS(i), OVmatS(ii), bw);

                                            end

                                            out_w = (1/N0)*sum(mv);

                                            %% Entropy of w

                                            Entw  =  Entw  +  log(OVmatS(i));

                                            %% Evaluate Sample at f(z) dist using KDE for all z

                                            %  z follows N(y,sigma2)
                                            % A synthetic observation is
                                            % generated by adding Gaussian
                                            % measurement noise to model
                                            % outputs at the candidate day.

                                            Measurements =  normrnd(OvS0(indT,spel,i),Abs_Error);

                                            %Cov = repmat(Abs_Error(:), length(indT), 1);

                                            if length(Measurements)==1

                                                mv=nan(1,N0);

                                                for ii = 1:N0

                                                    mv(1,ii) =  mvnpdf( Measurements(:)', reshape(OvS0(indT,spel,ii),1,length(spel)*length(indT)), diag((Abs_Error) ));

                                                end

                                            else

                                                 mv=nan(1,N0);

                                                for ii = 1:N0

                                                    mv(1,ii) =  mvnpdf( Measurements(:)', reshape(OvS0(indT,spel,ii),1,length(spel)*length(indT)), diag((Abs_Error).^2 ));

                                                end
                                            end

                                            out_z = (1/N0)*sum(mv);                                

                                            %% Evaluate Sample at Joint dist f(w,z) for all z (species)

                                            % p(W,Z): joint density of
                                            % ovulation time and all selected
                                            % measurements.
                                            wz=[OVmatS(i) Measurements(:)'];

                                            mv=nan(1,N0);

                                            for ii = 1:N0

                                                mv(1,ii) =  mvnpdf(wz, [OVmatS(ii) reshape(OvS0(indT,spel,ii),1,length(spel)*length(indT))],  diag(([bw Abs_Error]).^2));

                                            end

                                            out_wz = (1/N0)*sum(mv);   
                                            
                                            %% Mutual information of I(W,Z)     

                                            Info_ns = Info_ns + (log(out_wz) - log(out_w) - log(out_z));
                                            dataMI(p,i) = (log(out_wz) - log(out_w) - log(out_z));


                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %% Calculate individual mutual information: I(W,Z) where Z represents only one species
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                            % Repeat the information
                                            % calculation one species at a
                                            % time to rank measured outputs.
                                            out_wz_ind = nan(1,length(spel));
                                            out_wz_ind2 = nan(1,length(spel));


                                            for iii=1:length(spel)

                                                                %% Evaluate Sample at Joint dist f(w,z) for every z

                                                                Measure_oneS =  Measurements(:,iii);

                                                                wz=[OVmatS(i) Measure_oneS(:)'];

                                                                %Cov = repmat(Abs_Error(iii), length(indT), 1);

                                                mvWZ=nan(1,N0);

                                                for ii = 1:N0

                                                                mvWZ(1,ii)=mvnpdf(wz, [OVmatS(ii) reshape(OvS0(indT,iii,ii),1,length(indT))], diag(([bw Abs_Error(iii)]).^2) );
                                                end
                                                                out_wz_ind(1,iii) = (1/N0)*sum(mvWZ);

                                                                %% Evaluate Sample at f(z) dist for every z  

                                                mvZ=nan(1,N0);

                                                for ii = 1:N0

                                                                mvZ(1,ii)=mvnpdf(Measure_oneS(:)', reshape(OvS0(indT,iii,ii),1,length(indT)), diag(Abs_Error(iii).^2) );
                                                end
                                                                out_wz_ind2(1,iii) = (1/N0)*sum(mvZ);


                                            end

                                            %% Mutual information I(W,Z) for only one species
                                            
                                            Info_ns_ind = Info_ns_ind + (log(out_wz_ind(1,:)) - log(out_w) - log(out_wz_ind2(1,:)));                                         

                                            %% Entropy of every Z

                                             EntZ  =  EntZ  +  log(out_wz_ind2(1,:));
                                             

                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                            %% Mutual information I(W,Z) between W and all possible Combination of couple i.e. Z represents only 2 species.
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                            % Pairwise measurement terms are
                                            % used to identify redundant or
                                            % complementary species pairs.
                                            for iii=1:numCombZ

                                                ZZ=reshape(Measurements(:,CombZ(iii,:)),1,numel(Measurements(:,CombZ(iii,:))));

                                                iCov=repmat(Abs_Error(CombZ(iii,:)), length(indT), 1);

                                                %CovZZ=diag((reshape(iCov,1,numel(iCov))).^2);

                                                ZZcom=nan(1,N0);

                                                for ii = 1:N0

                                                    ZZcom(1,ii) =  mvnpdf( ZZ, reshape(OvS0(indT,CombZ(iii,:),ii),1, numel(OvS0(indT,CombZ(iii,:),ii))), diag((reshape(iCov,1,numel(iCov))).^2));
                                                end

                                                mvML(iii,1)=(1/N0)*sum(ZZcom);

                                                %% Mutual information of pair-Species Z-Z

                                                ns_ZZ(1,iii) = ns_ZZ(1,iii) + (log(mvML(iii,1)) - log(out_wz_ind2(1,CombZ(iii,1))) - log(out_wz_ind2(1,CombZ(iii,2))) );



                                            end


                                end


                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %% Normalizing the results
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                Info_TS(indTT,:)=(1/N0) * Info_ns_ind;

                                Info_OvT(1,indTT)=(1/N0) * Info_ns;

                                EntropyZ(indTT,:) = - (1/N0)*EntZ;
                                EntropyW(indTT,:) = - (1/N0)*Entw;

                                Info_ZZ(indTT,:) = (1/N0)*ns_ZZ;



     display('----plotting-----')


    %---------------------------------------------------------------------

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot joint distribution of every species and time of OV at time where mutual information is high      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % Rearranging the matrix OvS0 of 3D to 2D i.e. the dimension of OvS0(n,m,d) ---> becomes  OvS0(n*d,m) 
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                                A=[];

                                for j = 1:N0

                                    M_mean=cat(1,A,OvS0(p,:,j));

                                    A=M_mean;

                                end


                            %% Plotting

                            f100=figure(3);
                            
                            f4 = figure(tpoint(p) + 200);   


                            %% Defining color
                            
                            colorlist={[1 0 0], [1 0.6 0],[1 0.9 0],  [0.3 0.4 0.1], [0.7 0.8 0], [0 1 0], [0 0.1 0.7], [0 0.5 1], [0 1 1], [0.6 0.3 0.9], [1 0 1], [0 0 0] ...
                                [0.8 0.6 1], [0.1 0.4 0.5], [0.7 0.2 0.3]};



                            for indS = 1:length(spel)

                                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                           % Bivariate distribution density of W and Z    
                                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                           
                                           if str=='l'

                                                   set(0, 'CurrentFigure', f4)
                                                   subplot(3,3,indS)

                                                   ksdensity([OVmatS(1:size(A,1)) A(:,indS)], 'PlotFcn','contour');       
                                                   hold on
                                                   plot( OVmatS(1:size(A,1)), A(:,indS),'k.','MarkerSize',2);

                                                   xlabel('Ovulation time','fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                                   ylabel(Ynames(Species(spel(indS))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
                                                   title(sprintf('Day %d pp',tpoint(p)),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');

            
                                                    RHO = corr2(OVmatS(1:size(A,1)),A(:,indS));

                                                    Com = Info_TS(p,indS);

                                                    Corr = sign(Com)*sqrt(1 - exp(-2*abs(Com)));

                                                    xl=get(gca,'xlim');
                                                    yl=get(gca,'ylim');
                                                    text(xl(1),yl(2),sprintf('I=%.3f, p=%.3f',Corr,RHO),'fontweight','bold')
                                                   
                                           end
                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                   %% Distrubition density of Z 
                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                   set(0, 'CurrentFigure', f100)
                                                   subplot(3,3,indS)
                                                   hold on


                                                   [f,xi] = ksdensity(A(:, indS));

                                                   plot(xi,f,'color',colorlist{p}, 'LineWidth',2);

                                                   zlabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                                   ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                                   xlabel(Ynames(Species(spel(indS))) ,'fontsize',22,'Fontname','IPAPMincho');




                            end

                             legend(sprintf('D%d pp', Comb(1)),sprintf('D%d pp', Comb(2)),sprintf('D%d pp',Comb(3)),sprintf('D%d pp', tpoint(4))...
                                               ,sprintf('D%d pp', tpoint(5)),sprintf('D%d pp', tpoint(6)),sprintf('D%d pp', tpoint(7)),sprintf('D%d pp', tpoint(8)))

                                            %% Distrubition density of W 
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                            subplot(4,3,12)
                                            hold on

                                            [~,xi] = ksdensity(OVmat(1:NT));

                                            rang= 1:0.2:max(xi);

                                            [~, indf] = min(abs(rang-0));

                                            [f,xi] = ksdensity(OVmat(1:NT),rang);
                                            plot(xi(1:end),f(1:end),'k:', 'LineWidth',2);
                                            ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                            xlabel('Ovulation time','fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');

                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %% Distrubition density of W defined on the support interval [x1 xn], where x1 is changing
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            
                                            hold on   

                                            fS=f*(NT/N0);                                    
                                            
                                            [~, indfS] = min(abs(rang-tpoint(p)));
                                            
                                            fS(1:indfS+1)=0;

                                            plot(xi,fS,'color',colorlist{p}, 'LineWidth',2);
                                            ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                            

                                            


                                            
    
                                            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Selecting high mutual information and estimate the posterior density
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
   
        
        %% Randomly Selecting observed data 
        
                numMeas_to_plot=30;
                numx=36;
                numz=36;

                if N0 < numz          
                    numz=N0;
                    numx=N0;
                    numMeas_to_plot=N0;
                end

                xq = datasample(1:N0,numz,'Replace',false);           
 
        %% Ploting random mutual information using the random observed data xq
        
                figure(600)

                for i = 1:numx

                        subplot(numx,1,i)
                        bar(tpoint(p),dataMI(p,xq(i)),'FaceColor',[0.5 0.9 0.13]);
                        hold on
                end

                set(gca,'xticklabel',tpoint,'fontsize',8);
                ylabel(('I(\Theta,Y)'),'fontsize',22,'Fontname','IPAPMincho');
                xlabel('Period' ,'fontsize',22,'Fontname','IPAPMincho');
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Estimating posterior density for every selected data: P(W/Z*), where Z* is the observed data %%
        %%
        %%                            P(W,Z*)
        %% The idea is: P(W/Z*) =   ------------
        %%                          S P(w,z*)dw 
        %%                          "
        %%                          "
        %%                          "
        %%                          V
        %%                          S means integral
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %% Matrix representing the joint density at the selected time p of W and Z: Density_WZ = [W Z]  
        %% i.e. matrix of all possible value of W (ovulation time: represented by the vector OVmatS) 
        %% and Z (measurement: represented by the matrix A) at time of sampling p
        
                    Density_WZ = [OVmatS A];

        %% Discretization: increasing the number of points at which we want to estimate our density
        
                    range_w = Comb(p):100;   
                    gridx1 = unique([OVmatS; range_w']);

        %% Ranking the species containing high information at the time of sampling p
        
                     rank=p;
                     [~, maxspp]=sort(Info_TS(rank,:),'descend');

        
        %% We start estimating posterior p(w/z*), with different measurements z*.  
        %% numMeas_to_plot is the number of measurements
        % This block converts the design result into an interpretable
        % posterior plot: how the prior uncertainty about ovulation time
        % changes after a hypothetical measurement z* is collected.
        
                    figure(tpoint(p)+100)                                       
                    clf 

        for sam=1:numMeas_to_plot
            
            %% Generating data z*= Zval

                    Zval =  normrnd(OvS0(p,:,xq(sam)),repmat(Abs_Error, size(OvS0(p,:,xq(sam)),1),1));

                    Zval(Zval<0)=0;
            
            
            %% Discretization:            
            %% Matrix xii define point at which we want to estimate our density
            %% first column is for time ovulation, the other column are for measurement point 
            
                    xii = [gridx1, ones(length(gridx1),length(spel))];
            
            %% update xii matrix that have points at which you want to evaluate your posterior 
        
                    xii(:,2:end)= bsxfun(@times,Zval,xii(:,2:end));
            
            
            %% Estimate the joint density at w and the selected z* : p(w,z*) using the function 'mvksdensity'
            %% Dimension of z* is high: it includes all species
            %% f1: estimation of p(w,z*) where z* includes all species
            %% f2: estimation of p(w,z*) where z* includes all species except one which has lower rank 
            %% f3: ........................................................... the 2 last ones .......
            %% f4 .................................................................3..................
            %% .......................................................................................
            %% f8 .................................................................7..................
            
            
                    Cov = repmat(Abs_Error(:), length(p), 1);            

                    [f1,~] =  mvksdensity(Density_WZ,xii,'bandwidth',([bw Cov(:)']));
                    [f2,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-1)+1]),xii(:,[1 maxspp(1:(end)-1)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-1))']);
                    [f3,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-2)+1]),xii(:,[1 maxspp(1:(end)-2)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-2))']);
                    [f4,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-3)+1]),xii(:,[1 maxspp(1:(end)-3)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-3))']);
                    [f5,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-4)+1]),xii(:,[1 maxspp(1:(end)-4)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-4))']);
                    [f6,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-5)+1]),xii(:,[1 maxspp(1:(end)-5)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-5))']);
                    [f7,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-6)+1]),xii(:,[1 maxspp(1:(end)-6)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-6))']);
                    [f8,~]  = mvksdensity(Density_WZ(:,[1 maxspp(1:(end)-7)+1]),xii(:,[1 maxspp(1:(end)-7)+1]),'bandwidth',[bw Cov(maxspp(1:(end)-7))']);


            %% Calculating the normalization constant: S P(w,z*)dw       
                   
                    mv=nan(8,N0);

                    for ii = 1:N0

                        mv(1,ii) =  mvnpdf( Zval, reshape(OvS0(p,:,ii),1,length(spel)*length(p)), diag((Cov(:)').^2 ));

                        mv(2,ii) =  mvnpdf( Zval(maxspp(1:(end)-1)), OvS0(p,maxspp(1:(end)-1),ii), diag(Cov(maxspp(1:(end)-1))).^2);
                        mv(3,ii) =  mvnpdf( Zval(maxspp(1:(end)-2)), OvS0(p,maxspp(1:(end)-2),ii), diag(Cov(maxspp(1:(end)-2))).^2);
                        mv(4,ii) =  mvnpdf( Zval(maxspp(1:(end)-3)), OvS0(p,maxspp(1:(end)-3),ii), diag(Cov(maxspp(1:(end)-3))).^2);
                        mv(5,ii) =  mvnpdf( Zval(maxspp(1:(end)-4)), OvS0(p,maxspp(1:(end)-4),ii), diag(Cov(maxspp(1:(end)-4))).^2);
                        mv(6,ii) =  mvnpdf( Zval(maxspp(1:(end)-5)), OvS0(p,maxspp(1:(end)-5),ii), diag(Cov(maxspp(1:(end)-5))).^2);
                        mv(7,ii) =  mvnpdf( Zval(maxspp(1:(end)-6)), OvS0(p,maxspp(1:(end)-6),ii), diag(Cov(maxspp(1:(end)-6))).^2);
                        mv(8,ii) =  mvnpdf( Zval(maxspp(1:(end)-7)), OvS0(p,maxspp(1:(end)-7),ii), diag(Cov(maxspp(1:(end)-7))).^2);
                    end

                    I1 = (1/N0)*sum(mv(1,:));  
                    I2 = (1/N0)*sum(mv(2,:));
                    I3 = (1/N0)*sum(mv(3,:));
                    I4 = (1/N0)*sum(mv(4,:));
                    I5 = (1/N0)*sum(mv(5,:));  
                    I6 = (1/N0)*sum(mv(6,:));
                    I7 = (1/N0)*sum(mv(7,:));
                    I8 = (1/N0)*sum(mv(8,:));

            

            %% Estimating posterior density for every selected data: P(W/Z*), where Z* is the observed data %%
            %%
            %%                            P(W,Z*)
            %% The idea is: P(W/Z*) =   ------------
            %%                          S P(w,z*)dw 
            
                    Post1 = f1/I1;
                    Post2 = f2/I2;
                    Post3 = f3/I3;
                    Post4 = f4/I4;
                    Post5 = f5/I5;
                    Post6 = f6/I6;
                    Post7 = f7/I7;
                    Post8 = f8/I8;

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %% Potting Prior and Posterior Distrubition of W after collecting data Z
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


           subplot(6,6,sam)
           hold on

                    
                    %% Prior density of W 
                    
                    %plot(xi(indfS:end),fS(indfS:end),'-.b', 'LineWidth',2);
                    plot(xi,fS,'-.b', 'LineWidth',2);
                    Ar0=trapz(xi,fS);
                    
                    ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    xlim([Comb(p)-10 100])
                    box on
                    hold on


                    % Posterior density of W
       
                    plot(gridx1,Post1,'r', 'LineWidth',2);
                    Ar1=trapz(gridx1,Post1);
                    hold on
                    plot(gridx1,Post2,'k', 'LineWidth',2);
                    Ar2=trapz(gridx1,Post2);
                    hold on
                    plot(gridx1,Post3,'g', 'LineWidth',2);
                    Ar3=trapz(gridx1,Post3);
                    hold on
                    plot(gridx1,Post4,'m', 'LineWidth',2);
                    Ar4=trapz(gridx1,Post4);
                    hold on
                    plot(gridx1,Post5,'--k', 'LineWidth',2);
                    Ar5=trapz(gridx1,Post5);
                    hold on
                    plot(gridx1,Post6,'--g', 'LineWidth',2);
                    Ar6=trapz(gridx1,Post6);
                    hold on
                    plot(gridx1,Post7,'--m', 'LineWidth',2);
                    Ar7=trapz(gridx1,Post7);
                    hold on
                    plot(gridx1,Post8,'--r', 'LineWidth',2);
                    Ar8=trapz(gridx1,Post8);
                    hold on
                    
                    legend(sprintf('A0= %d', Ar0),sprintf('A1= %d', Ar1),sprintf('A2= %d', Ar2),sprintf('A3= %d', Ar3),sprintf('A4= %d', Ar4)...
                    ,sprintf('A5= %d', Ar5),sprintf('A6= %d', Ar6),sprintf('A7= %d', Ar7),sprintf('A8= %d', Ar8))

                    
                    ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    xlabel('Ovulation time','fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');



                    %% Defining time windows through which the posterior predicts well the ovulation
                    
                    if abs(OVmatZ(xq(sam))-xii(find(Post1==max(Post1)))) <= 5

                        set(gca,'Color',[0.85 0.85 0.85])

                    end

                    hold on
                    line([tpoint(p) tpoint(p)],[0 max([fS,Post1'])],'color','m','LineStyle','--')
                    line([OVmatZ(xq(sam)) OVmatZ(xq(sam))],[0 max([fS,Post1'])],'color','m','LineStyle','--')
                    
                    %suptitle(sprintf('measurement at days %d',tpoint(p)) ,'fontsize',22,'Fontname','IPAPMincho');

                   annotation('textbox', [0.5 0.9 0.1 0.1],'String', sprintf('Posterior P(W/z*), z* includes all species at day %d',tpoint(p)),'fontsize',20,'Fontname','IPAPMincho', 'HorizontalAlignment', 'center')

                    

        end
        
                     idnames1=Ynames(Species(maxspp(1:(end))));
                     idnames2=Ynames(Species(maxspp(1:(end)-1)));
                     idnames3=Ynames(Species(maxspp(1:(end)-2)));
                     idnames4=Ynames(Species(maxspp(1:(end)-3)));
                     idnames5=Ynames(Species(maxspp(1:(end)-4)));
                     idnames6=Ynames(Species(maxspp(1:(end)-5)));
                     idnames7=Ynames(Species(maxspp(1:(end)-6)));
                     idnames8=Ynames(Species(maxspp(1:(end)-7)));



                     h_legend=legend('Prior',sprintf(' %s',idnames1{:}),sprintf(' %s',idnames2{:}),sprintf(' %s',idnames3{:}),...
                         sprintf(' %s',idnames4{:}),sprintf(' %s',idnames5{:}),sprintf(' %s',idnames6{:})...
                         ,sprintf(' %s',idnames7{:}),sprintf(' %s',idnames8{:}));
                    set(h_legend,'FontSize',6)   


                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                    
                    
                    
                    
                    
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Selecting high mutual information recorded for a species and estimate its the posterior density
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
         rank=p;
         [~, maxsp]=sort(Info_TS(rank,:),'descend');

         figure(tpoint(p))                                       
         clf

         %% loading models
             
         Density_WZmax = [OVmatS A(:,maxsp(1))];    
         Density_WZmed = [OVmatS A(:,maxsp(2))];    
         Density_WZmed3 = [OVmatS A(:,maxsp(3))];    
         Density_WZmin = [OVmatS A(:,maxsp(end))];    


         %xq=datasample(1:size(OvS0,3),numMeas,'Replace',false);


         for sam=1:numMeas_to_plot
             
             
             
                    %% Collect data Z at random 

                    xiimax = [gridx1, ones(length(gridx1),1)];
                    xiimed = [gridx1, ones(length(gridx1),1)];
                    xiimed3 = [gridx1, ones(length(gridx1),1)];
                    xiimin = [gridx1, ones(length(gridx1),1)];

                    Zmax =  normrnd(OvS0(p,maxsp(1),xq(sam)),repmat(Abs_Error(maxsp(1)), size(OvS0(p,maxsp(1),xq(sam)),1),1));
                    Zmax(Zmax<0)=0;
                    
                    Zmed =  normrnd(OvS0(p,maxsp(2),xq(sam)),repmat(Abs_Error(maxsp(2)), size(OvS0(p,maxsp(2),xq(sam)),1),1));
                    Zmed(Zmed<0)=0;
                    
                    Zmed3 =  normrnd(OvS0(p,maxsp(3),xq(sam)),repmat(Abs_Error(maxsp(3)), size(OvS0(p,maxsp(3),xq(sam)),1),1));
                    Zmed3(Zmed3<0)=0;
                    
                    Zmin =  normrnd(OvS0(p,maxsp(end),xq(sam)),repmat(Abs_Error(maxsp(end)), size(OvS0(p,maxsp(end),xq(sam)),1),1));
                    Zmin(Zmin<0)=0;

%                     Zmax=OvS0(rank,maxsp(1),xq(sam));  
%                     Zmed=OvS0(rank,maxsp(2),xq(sam)) ;   
%                     Zmed3=OvS0(rank,maxsp(3),xq(sam)) ;                                       
%                     Zmin=OvS0(rank,maxsp(end),xq(sam)) ;
                     

                    %% update your matrix that have ponit at which you want to evaluate your posterior 

                    xiimax(:,2)= bsxfun(@times,Zmax,xiimax(:,2));
                    xiimed(:,2)= bsxfun(@times,Zmed,xiimed(:,2));
                    xiimed3(:,2)= bsxfun(@times,Zmed3,xiimed3(:,2));
                    xiimin(:,2)= bsxfun(@times,Zmin,xiimin(:,2));


                    %% Normalizing the joint density f 

                    mv=nan(4,N0);

                    Covmax = Abs_Error(maxsp(1));
                    Covmed = Abs_Error(maxsp(2));    
                    Covmed3 = Abs_Error(maxsp(3));    
                    Covmin = Abs_Error(maxsp(end));    

                    for ii = 1:N0

                        mv(1,ii) =  normpdf( Zmax, OvS0(p,maxsp(1),ii), Covmax);
                        mv(2,ii) =  normpdf( Zmed, OvS0(p,maxsp(2),ii), Covmed);
                        mv(3,ii) =  normpdf( Zmed3, OvS0(p,maxsp(3),ii), Covmed3);
                        mv(4,ii) =  normpdf( Zmin, OvS0(p,maxsp(end),ii), Covmin);

                    end

                    Imax = (1/N0)*sum(mv(1,:));  
                    Imed = (1/N0)*sum(mv(2,:)); 
                    Imed3 = (1/N0)*sum(mv(3,:)); 
                    Imin = (1/N0)*sum(mv(4,:)); 

                    %% Calculate joint density at W and the selected Z 
                    
                    [fmax,xiimax] = ksdensity(Density_WZmax,xiimax,'bandwidth',([bw Covmax]));
                    [fmed,xiimed] = ksdensity(Density_WZmed,xiimed,'bandwidth',([bw Covmed]));
                    [fmed3,xiimed3] = ksdensity(Density_WZmed3,xiimed3,'bandwidth',([bw Covmed3]));
                    [fmin,xiimin] = ksdensity(Density_WZmin,xiimin,'bandwidth',([bw Covmin]));
                    
                    %f = mvksdensity(Density_WZ,xii,'bandwidth',([bw Cov(:)']));

                    %% Posterior is 

                     Postmax = fmax/Imax;
                     Postmed = fmed/Imed;
                     Postmed3 = fmed3/Imed3;
                     Postmin = fmin/Imin;
                     
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %% Potting Prior and Posterior Distrubition of W after collecting data Z
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           subplot(6,5,sam)
           hold on

                    
                    %% Prior density of W 
                    
              
                    
                    plot(xi(indfS:end),fS(indfS:end),'g+', 'LineWidth',2);
                    ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    xlim([Comb(p)-10 100])
                    box on
                    hold on


                    % Posterior density of W
       
                    plot(gridx1,Postmax,'r', 'LineWidth',2);

                    hold on
                    
                    if abs(OVmatZ(xq(sam))-xii(find(Postmax==max(Postmax)))) <= 5

                        set(gca,'Color',[0.85 0.85 0.85])

                    end
                    
                    % Posterior density of W
       
                    plot(gridx1,Postmed,'b', 'LineWidth',2);
                    
                    hold on            
                    
                    % Posterior density of W
       
                    plot(gridx1,Postmed3,'color',[0 0.6 0.3], 'LineWidth',2);
                    
                    hold on
                
                    % Posterior density of W
       
                    plot(gridx1,Postmin,'k', 'LineWidth',2);
                    
                    hold on



                    ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    xlabel('Ovulation time','fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');

                    hold on
                    line([tpoint(p) tpoint(p)],[0 max([fS,Postmax',Postmed',Postmed3',Postmin'])],'color',[0.2 0 0],'LineStyle','--')
                    line([OVmatZ(xq(sam)) OVmatZ(xq(sam))],[0 max([fS,Postmax',Postmed',Postmed3',Postmin'])],'color',[0.2 0 0],'LineStyle','--')



                    subplot(6,5,26)
                    hold on       
                    bar(sam,Zmin,'k')
                                    ylabel(Ynames(Species(spel(maxsp(end)))),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                    xlabel(sprintf('measurement at days %d',tpoint(rank)),'fontweight','bold','fontsize',16,'Fontname','IPAPMincho');
                                    box on


                    subplot(6,5,27)
                    hold on       
                    bar(sam,Zmed3,'FaceColor',[0 0.6 0.3])
                                    ylabel(Ynames(Species(spel(maxsp(3)))),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                    xlabel(sprintf('measurement at days %d',tpoint(rank)),'fontweight','bold','fontsize',16,'Fontname','IPAPMincho');
                                    box on

                    subplot(6,5,28)
                    hold on       
                    bar(sam,Zmed,'b')
                                    ylabel(Ynames(Species(spel(maxsp(2)))),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                    xlabel(sprintf('measurement at days %d',tpoint(rank)),'fontweight','bold','fontsize',16,'Fontname','IPAPMincho');
                                    box on

                    subplot(6,5,29)
                    hold on 

                    if abs(OVmatZ(xq(sam))-xi(find(Postmax==max(Postmax)))) <= 5

                        bar(sam,Zmax,'FaceColor', [0.75 0.75 0.75])
                        box on
                    else

                        bar(sam,Zmax,'r')
                        box on
                    end

                                    ylabel(Ynames(Species(spel(maxsp(1)))),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                                    xlabel(sprintf('measurement at days %d',tpoint(rank)),'fontweight','bold','fontsize',16,'Fontname','IPAPMincho');

                   annotation('textbox', [0.5 0.9 0.1 0.1],'String', sprintf('Posterior P(W/z*), z* includes only one species at day %d',tpoint(p)),'fontsize',20,'Fontname','IPAPMincho', 'HorizontalAlignment', 'center')


         end
         
         
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                    
         
        
        
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Selecting samples that the prediction of the model is close to the the prediction of the posterior 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         
         ConSp=[];
         for sam=1:numz
             
                    mv=nan(1,N0);

                    Covmax = Abs_Error(maxsp(1));
                    
                    Zmax =  normrnd(OvS0(p,maxsp(1),xq(sam)),repmat(Abs_Error(maxsp(1)), size(OvS0(p,maxsp(1),xq(sam)),1),1));
                    Zmax(Zmax<0)=0;
                    
                    for ii = 1:N0

                        mv(1,ii) =  normpdf( Zmax, OvS0(p,maxsp(1),ii), Covmax);

                    end

                    Imax = (1/N0)*sum(mv(1,:));  

                    %% Calculate joint density at W and the selected Z 
                    
                    [fmax,xiimax] = ksdensity(Density_WZmax,xiimax);
                 
     
                    %f = mvksdensity(Density_WZ,xii,'bandwidth',([bw Cov(:)']));

                    %% Posterior is 

                     Postmax = fmax/Imax;
    

                    if abs(OVmatZ(xq(sam))-xii(find(Postmax==max(Postmax)))) <= 2

                        ConSp=[ConSp Zmax];

                    end


         end


         subplot(6,5,30)
         hold on

         [~,xi] = ksdensity(ConSp);
         rang= min(xi):0.002:max(xi);
         [f,xi] = ksdensity(ConSp,rang);
         plot(xi,f,'r', 'LineWidth',2);
         xlabel(strcat(Ynames(Species(spel(maxsp(1)))),' at OVtime'),'fontweight','bold','fontsize',16,'Fontname','IPAPMincho');
         ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
         box on







        %% Plotting Mutual information between species W

        figure(4)

            subplot(3,4,p)
            bar(Info_ZZ(p,:))
            title(sprintf('D%d',Comb(p,:)),'fontsize',22,'Fontname','IPAPMincho');


        suptitle('Mutual information between species')


        figure(5)

            subplot(3,4,p)
            bar(EntropyZ(p,:))
            title(sprintf('D%d',Comb(p,:)),'fontsize',22,'Fontname','IPAPMincho');


        suptitle('Entropy of every species')


       figure(800+p)
       clf
              for i=1:size(CombZ,1)
           
                    subplot(5,6,i)
           
                    ksdensity([A(:,CombZ(i,1)) A(:,CombZ(i,2))], 'PlotFcn','contour'); 
                    hold on
                    plot(A(:,CombZ(i,1)),  A(:,CombZ(i,2)),'k.','MarkerSize',2);
                    
                    xlabel(Ynames(Species(CombZ(i,1))),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    ylabel(Ynames(Species(CombZ(i,2))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
                    
                    RHO = corr2(A(:,CombZ(i,1)),A(:,CombZ(i,2)));
                     
                    Com = Info_ZZ(p,i);
                    
                    Corr = sign(Com)*sqrt(1 - exp(-2*abs(Com)));
                    
                    xl=get(gca,'xlim');
                    yl=get(gca,'ylim');
                    text(xl(1),yl(2),sprintf('I=%.3f, p=%.3f',Corr,RHO),'fontweight','bold')

                    
                    
                  
              end
       
              suptitle('Dependency from high mutual information');

    datetime('now')
    
    end
    
    


    
    

%% Plot mutual information at different ts for all species       
       
%% Figure 1

figure(1)
clf


    hb=bar(Info_OvT);
   
    hold on

    [~, inMI]=max(Info_OvT);
    hbr = bar(inMI,Info_OvT(1,inMI),'r');

    set(gca,'xticklabel',tpoint,'fontsize',8);
    
    ylabel(('I(Tov,Y)'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
    xlabel('Day' ,'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
    title('Mutual Information','fontsize',22,'Fontname','IPAPMincho');

%% Plot mutual information at different ts for one species       

    figure(2)

    for p=1:length(tpoint)
        
        subplot(3,4,p) 
        bar(Info_TS(p,:));
        set(gca,'xticklabel',Ynames(Species(spel)),'fontsize',8);

        ylabel(('I(Tov,Y_i)'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
        %xlabel('Model components' ,'fontsize',22,'Fontname','IPAPMincho');
        title(sprintf('Day %d',Comb(p,:)),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
        %title(sprintf('D%d',length(tpoint),'fontsize',22,'Fontname','IPAPMincho');
    end

    hold on 
    
    subplot(3,4,inMI) 

    hbr = bar(Info_TS(inMI,:),'r');
    
    
    figure(6)

    bar(EntropyW(:,1))
    set(gca,'xticklabel',tpoint,'fontsize',8);
    suptitle('Entropy of time of Ovulation')
    


    


    %AA=cat(2,(1:numComb)', Comb)

    h = get(0, 'children');

  if Ni==0 
      
    for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_L/Figures_L1/figure' num2str(i)], 'fig');
    end
    
  elseif Ni==3000

          for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_L/Figures_L2/figure' num2str(i)], 'fig');
          end
    
  elseif Ni==6000
          for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_L/Figures_L3/figure' num2str(i)], 'fig');
          end
          
  elseif Ni==9000
          for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_L/Figures_L4/figure' num2str(i)], 'fig');
          end

  end
  
  numSam
  
end

    
    
    
    
% display('=========================================== Non_Lactating paper approach ===================================================')
% 
% 
%  clear all
% close all
% 
% 
% display(' -------------------------  Lactating ------------------------')
% 
% 
% display('----Sampling-----')
% 
% %% State variable of the model
% 
% Ynames =  { 'GnRhH'; 'GnRhP';'FSHP'; 'FSH'; 'LHP';  'LH'; 'Foll'; 'PGF';   'CL';    'P4';  'E2'; .....
%            'INH'; 'ENZ'; 'OXT';  'IOF'; 'IGF';'Ins'; 'Glu'; 'Fat';   'Glucose in Liver'; 'Stored Glucose'};
%        
% %% Defining Species to include in simulation
% 
%         Species= 1:21; % Species to simulate in Metabolism (Species 1 and 2) and Reproduction model (the rest) 
%         
%         Species= Species([4,8,10,11,12,16,17,18]);
%         
%  for Ni=[1]
%        
%      sprintf('-------- Ni=%d ------',Ni)
% 
%      
%         %% Number of samples
% 
%         N0= 20;
%         
%         %% Argument to pass to the solver
% 
%         frac2=0.25; 
%         STEP=1/6;
%         day=100;
%         daysVOR=17;
%                 
% 
%         %% Reading Paramaters
% 
%         par=BovSys_para_1M();
%         
%         
%         %% Sampling theta from Uniform dis...
% 
%         % Defining the boundary of parameters interval
%         
%         devS = 0.01;
%         
%         Sens= [12    50    92    40    48    49    29    34   101    58 ...
%             57    38    60    81    51    79    72   100    76    88    ...
%             74    46    47    55    41    87    98    62    35    36    33    56];
%     
%     
%         aRep=nan(length(par),1);
%         bRep=nan(length(par),1);
% 
%         NotSens = setxor(1:length(par),Sens);
% 
%         aRep(Sens)=par(Sens) - par(Sens)*devS;
%         bRep(Sens)=par(Sens) + par(Sens)*devS;
%         
%         devNotS = 0.05;
%         
%         aRep(NotSens)=par(NotSens) - par(NotSens)*devNotS;
%         bRep(NotSens)=par(NotSens) + par(NotSens)*devNotS;        
%         
% 
%         %% Check if parameters are negative
% 
%         if any(aRep<0)
%             aRep(find(aRep<0))=0;
%         end
% 
%         %% Sampling from parameters from uniform distribution ...
%         
%         repa=repmat(aRep,1,N0);
%         repb=repmat(bRep,1,N0);
%     
%         thetaRep0 = unifrnd(repa,repb);
% 
%         %% Run the standard model to calculate the error ...
%         
%                      str='n';   % we switch to non-lactating to calculate standard P4 model
%     
%                      RefP4=0;
%     
%                      [T,YRepS,OvTS]=BovSys_run_1M(par,frac2,day,RefP4,str,STEP);
%                      
%                      day=day-daysVOR;
%                                         
%                      RefP4=sum(YRepS(((day-daystr)*(1/STEP)):(day*(1/STEP)),10));
% 
%                      Out = YRepS;
%           
%                      Abs_Error_NL = 0.1*mean(Out(:,Species),1);
%                      
%                      
%                      %str= input('Lactating, l or nonlactating, n: ','s');   % switch back to lactating case!!!
%                
%                                           
% %                      figure(55)
% %                         
% %                      plot(T,YRepS(:,10),'k')
% %                      hold on
% %                      plot(T,0.4*YRepS(:,10),'r')
%        
%                     str='n';
%                     
%                     %% Nonlactating case
%                     
%                     tpoint = 19:1:78;
% 
%                     Comb=reshape(tpoint,length(tpoint),1);
% 
%                     numComb=size(Comb,1);
%         
%                     
%         
%         %% All possible Combination of pair-Species mutual information
%         
% 
%                     CombZ = nchoosek(1:length(Species),2);
% 
%                     numCombZ = size(CombZ,1);
% 
%                     mvML=nan(numCombZ,1);
%     
% 
%     
% 
%         %% Run the model for different theta (smapled parameters)
%         
%                     %% First define matrix for input simulation
% 
%                     OvS0_NL  = zeros(length(tpoint),(length(Species))) ;
%                     PAR  = nan(length(par),N0) ;
% 
% 
%                     %OVmat_NL = nan(N0,1);
% 
% 
% 
% tic
%         
% display('----Generate models -----')
% 
% 
%        day=day+daysVOR;
%        
%        %colorlist={'k','r','b','g'};
%        
%        parfor (h=1:N0,16)
%             
%                     %% Run the model for different theta (smapled parameters)
% 
%                      
%                      
%                      [T,Y0,OvT0]=BovSys_run_1M(thetaRep0(:,h),frac2,day, RefP4,str,STEP);
%                      
%                      OT=(T-daysVOR);
%                      T=OT(find(OT==0):end);
%                      Y0=Y0(find(OT==0):end,:);
%                      
%                              %figure(1000)
%         
%                              %plot(T,Y0(:,10),'color',colorlist{h})
%                              %hold on
%                              
%                      Out = Y0;
%                      
%                      %% check if the model return nan
% 
%                      if isnan(OvT0)
% 
%                              OvS0_NL(:,:,h)=nan(length((tpoint)),(length(Species))) ;                       
% 
%                      else
% 
%                              OvS0_NL(:,:,h) = Out((tpoint)*(1/STEP),Species);
%                              %OVmat_NL(h,1)  = OvT0;
%                              PAR(:,h) =thetaRep0(:,h);
%                    
%                      end                 
%                         
%        end
%        
%        
%        PAR=reshape(PAR(~isnan(PAR)),length(par),length(PAR(~isnan(PAR)))/length(par));
%        
%         
%        day=day-daysVOR;
% 
%        %% Select the successfully simulated samples and reject the nan 
%        
%                     [~, ~, dim]= ind2sub(size(OvS0_NL),find(isnan(OvS0_NL)));
% 
%                     if ~isempty(dim)
%                         
%                         OvS0_NL(:,:,unique(dim))=[];
%                         
%                         
%                           if Ni==0  
%                             
%                               %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL/Data_1/Abs_Error_NL.mat','Abs_Error_NL')
%                               %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL/Data_1/OvS0_NL.mat','OvS0_NL')
% 
% 
%                           elseif Ni==1
% 
%                               %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL/Data_2/Abs_Error_NL.mat','Abs_Error_NL')
%                               %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL/Data_2/OvS0_NL.mat','OvS0_NL')
% 
% 
%                           elseif Ni==2daysVOR
%                               
%                               %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL/Data_3/Abs_Error_NL.mat','Abs_Error_NL')
%                               %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL/Data_3/OvS0_NL.mat','OvS0_NL')
% 
%                           end
%                         
% 
%                         display('............')
% 
%                         
%                         %OVmat_NL=OVmat_NL(~isnan(OVmat_NL));
% 
%                         %save('OVmat_NL.mat','OVmat_NL')
%                         
%                         %% update the number of samples
% 
%                         N0 = N0 - length(unique(dim));
% 
%                     else
%                         
%                         display('...llll....')
% 
%                         save('OvS0_NL.mat','OvS0_NL');
% 
%                         
%                         %save('OVmat_NL.mat','OVmat_NL')
% 
%         
%                     end
% 
%                     AcceptedOutput=sprintf('Selected samples N0=%d ... Number of Rejected samples are: %d', N0, length(unique(dim)))
% 
% toc
% datetime('now')
% % 
% %         %% here we start calculating pdfs:  ---------------------------------------------------------------------------------
% 
% 
%        
% %---------------------------------------------------- 1 -------------------------------------------------------------------  
% 
% 
% 
% 
%             N00=N0-round(N0/2)-1;
%             N1=N0-(N00+1);
%         
%             Info_Dtime = nan(1,numComb);
%             Info_TS = nan(numComb,length(Species));
%             Info_ZZ = nan(numComb,numCombZ);
%             dataMI= nan(numComb,N0);
% 
% 
%             
%             
% tic
% 
% 
%             for ind=1:numComb
%                 
%                 subComb = Comb(ind,:);
%                 [~,a]=ismember(subComb,tpoint);
%                 indT=a(a~=0);
%                 
%                 Info_ns=0;
%                 Info_ind=zeros(1,length(Species));
%                 ns_ZZ=zeros(1,numCombZ);
% 
%                 for i=1:N0
% 
%                     %% Select time point and its corresponding model output at which the experiment is carried out 
% 
%                                 %% Perturbate the selected model output   
% 
%                                 
%                                 Measurements =  normrnd(OvS0_NL(indT,:,i),repmat(Abs_Error_NL, size(OvS0_NL(indT,:,i),1),1));
%                                 
%                                 Measurements(Measurements<0)=0;
% 
%                                 Cov = repmat(Abs_Error_NL, length(indT), 1);
%                                 
%                                 %% Calculate the pdf
%                                         
%                                 if length(Measurements)==1
% 
%                                    Pdf0=  normpdf( Measurements(:)', reshape(OvS0_NL(indT,:,i),1,length(Species)*length(indT)), diag(Cov(:)'));
%                                             
%                                 else
% 
%                                    Pdf0 = mvnpdf( Measurements(:)', reshape(OvS0_NL(indT,:,i),1,length(Species)*length(indT)), diag((Cov(:)').^2 ));
%                                    
%                                   
% 
%                                 end
%                                         
% 
%                                 Pdf1 = 0;
% 
%                                 for j=1:N0
% 
% 
%                                     if length(Measurements)==1
% 
%                                         Pdf1 = Pdf1 +  normpdf( Measurements(:)', reshape(OvS0_NL(indT,:,j),1,length(Species)*length(indT)), diag(Cov(:)'));
% 
%                                     else
% 
%                                         Pdf1 = Pdf1 + mvnpdf( Measurements(:)', reshape(OvS0_NL(indT,:,j),1,length(Species)*length(indT)), diag((Cov(:)').^2 ));
% 
%                                     end
%                                 
% 
%                                 end
%                                 
%                                 dataMI(ind,i)=(log(Pdf0) - log(Pdf1/N1));
% 
%                                 Info_ns = Info_ns + (log(Pdf0) - log(Pdf1/N1));
%                                 
%                                 %% Mutual information for every species
%                                 
%                                 
%                                 for iii=1:length(Species)
%                                                             
%                                             Measure_oneS =  Measurements(:,iii);                                                           
%                                                             
%                                             Cov = repmat(Abs_Error_NL(iii), length(indT), 1);
% 
%                                             %% Calculate the pdf
% 
%                                             if length(Measurements)==1
% 
%                                                Pdf0=  normpdf( Measure_oneS(:)', reshape(OvS0_NL(indT,iii,i),1,length(indT)), diag(Cov(:)'));
% 
%                                             else
% 
%                                                Pdf0 = mvnpdf( Measure_oneS(:)', reshape(OvS0_NL(indT,iii,i),1,length(indT)), diag((Cov(:)').^2 ));
% 
%                                             end
% 
% 
%                                             Pdf1 = 0;
% 
%                                             for j=1:N0
% 
% 
%                                                 if length(Measurements)==1
% 
%                                                     Pdf1 = Pdf1 +  normpdf( Measure_oneS(:)', reshape(OvS0_NL(indT,iii,j),1,length(indT)), diag(Cov(:)'));
% 
%                                                 else
% 
%                                                     Pdf1 = Pdf1 + mvnpdf( Measure_oneS(:)', reshape(OvS0_NL(indT,iii,j),1,length(indT)), diag((Cov(:)').^2 ));
% 
%                                                 end
% 
% 
%                                             end
% 
%                                             Info_ind(1,iii) = Info_ind(1,iii) + (log(Pdf0) - log(Pdf1/N1));
% 
%                                 end
%                                 
%                                         %% individual mutual information
%                                         
%                                         out_wz_ind2 = nan(1,length(Species));
%                                                             
% 
%                                         for iii=1:length(Species)
% 
%                                                             
%                                             %% Evaluate Sample at f(z) dist for every z  
%                                                             
%                                             mvZ=nan(1,N0);
%    
%                                             for ii = 1:N0
% 
%                                                             mvZ(1,ii)=mvnpdf(Measure_oneS(:)', reshape(OvS0_NL(indT,iii,ii),1,length(indT)), diag((Cov(:)').^2) );
%                                             end
%                                                             out_wz_ind2(1,iii) = (1/N0)*sum(mvZ);
%                             
%                                                             
%                                         end
%                                 
%                                 
%                                         %% All possible Combination of couple mutual information
%                                         
%                                                                                 
%                                         for iii=1:numCombZ
% 
%                                                 ZZ=reshape(Measurements(:,CombZ(iii,:)),1,numel(Measurements(:,CombZ(iii,:))));
% 
%                                                 iCov=repmat(Abs_Error_NL(CombZ(iii,:)), length(indT), 1);
% 
%                                                 CovZZ=diag((reshape(iCov,1,numel(iCov))).^2);
% 
%                                                 ZZcom=nan(1,N0);
% 
%                                                 for ii = 1:N0
% 
%                                                     ZZcom(1,ii) =  mvnpdf( ZZ, reshape(OvS0_NL(indT,CombZ(iii,:),ii),1, numel(OvS0_NL(indT,CombZ(iii,:),ii))), CovZZ);
%                                                 end
% 
%                                                 mvML(iii,1)=(1/N0)*sum(ZZcom);
% 
%                                                 %% Mutual information of pair-Species Z-Z
% 
%                                                 ns_ZZ(1,iii) = ns_ZZ(1,iii) + (log(mvML(iii,1)) - log(out_wz_ind2(1,CombZ(iii,1))) - log(out_wz_ind2(1,CombZ(iii,2))) );                                         
% 
%                                         end                             
%                                 
%                                 
%                                 
% 
%                 end
% 
%                 Info_Dtime(1,ind)=(1/N0) * Info_ns;
% 
%                 Info_TS(ind,:)=(1/N0) * Info_ind;
%                 
%                 Info_ZZ(ind,:) = (1/N0)*ns_ZZ;
% 
% 
% 
%             end
%             
%    
% 
% 
% %% Pigure 1
% 
% 
% figure(200)
% clf;
% 
% subplot(2,1,1)
% 
%     bar(Info_Dtime,'FaceColor',[0.93 0.69 0.13]);
%    
%     hold on
%     %set(gca,'xticklabel',1:length(tpoint),'fontsize',8);
%     
%     ylabel(('I(\Theta,Y)'),'fontsize',22,'Fontname','IPAPMincho');
%     xlabel('Period' ,'fontsize',22,'Fontname','IPAPMincho');
% 
% figure(300)
% %clf;
% 
% for np=1:16%numComb
% 
%     subplot(4,4,np)
%     
%             bar(Info_TS(np,:),'r');
% 
%     set(gca,'xticklabel',Ynames(Species),'fontsize',8);
%     
%     ylabel(('I(Tov,Y)'),'fontsize',18,'Fontname','IPAPMincho');
%     %xlabel('Model components' ,'fontsize',18,'Fontname','IPAPMincho');
%     title(sprintf('D%d',Comb(np,:)),'fontsize',18,'Fontname','IPAPMincho');
%     %title(sprintf('D%d',length(tpoint),'fontsize',18,'Fontname','IPAPMincho');
% 
% end
% 
% 
%     figure(800)
%     
%     for j=1:16%numComb
% 
%         subplot(4,4,j)
%         bar(Info_ZZ(j,:))
%         title(sprintf('D%d',Comb(j,:)),'fontsize',18,'Fontname','IPAPMincho');
% 
%     
%     end
% 
% 
% 
% 
%        
%     %%-----------------------------------------------------------------------------------------------------
%        
% 
%     %% Selecting high mutual information
% 
%      rank=find(Info_Dtime==max(Info_Dtime));
%      rankn=find(Info_Dtime==min(Info_Dtime));
%      
%      rank=rank(1);
%      rankn=rankn(1);
%      
%      [~, maxsp]=sort(Info_TS(rank,:),'descend');
%      [~, minsp]=sort(Info_TS(rankn,:),'descend');
% 
%                 
%                 subComb = Comb(rank,:);
%                 [~,a]=ismember(subComb,tpoint);
%                 indmax=a(a~=0);
%                 
% 
%                 subComb = Comb(rankn,:);
%                 [~,a]=ismember(subComb,tpoint);
%                 indmin=a(a~=0);
%      
%    
%        
% 
%        
% 
%                             Amax=[];
%                             Amin=[];
%                             Amaxmin=[];
%                             Aminmax=[];
%                             Amaxall=[];
%                             Aminall=[];
% 
%                             for j = 1:N0
%                                 
%                                 M_meanmax=cat(1,Amax,OvS0_NL(indmax,maxsp(1),j)'); 
%                                 Amax=M_meanmax;
%                                 
%                                 M_meanmaxmin=cat(1,Amaxmin,OvS0_NL(indmax,maxsp(end),j)'); 
%                                 Amaxmin=M_meanmaxmin;
%                                 
%                                 M_meanmin=cat(1,Amin,OvS0_NL(indmin,minsp(end),j)'); 
%                                 Amin=M_meanmin;
%                                 
%                                 M_meanminmax=cat(1,Aminmax,OvS0_NL(indmin,minsp(1),j)'); 
%                                 Aminmax=M_meanminmax;
%                                 
%                                 maxall=cat(1,Amaxall,OvS0_NL(indmax,:,j)); 
%                                 Amaxall=maxall;
%                                 
%                                 minall=cat(1,Aminall,OvS0_NL(indmin,:,j)); 
%                                 Aminall=minall;
% 
%                             end
%            
%                             
% 
%                                     %% loading models
%                                     
%         figure(1000)
%         
%         for i=1:16
%              
%              Density_Amax = [PAR(Sens(i),:)' Amax];    
%              Density_Amaxmin = [PAR(Sens(i),:)' Amaxmin];    
%              Density_Amin = [PAR(Sens(i),:)' Amin];    
%              Density_Aminmax = [PAR(Sens(i),:)' Aminmax];    
%              
%              Density_Amaxall = [PAR(Sens(i),:)' Amaxall];    
%              Density_Aminall = [PAR(Sens(i),:)' Aminall]; 
%              
%         
%              gridx1 = unique(PAR(Sens(i),:))';
%              
%              gridx1 = linspace((min(gridx1) - min(gridx1)/60),(max(gridx1)+max(gridx1)/60),500)';
% 
% 
%   
%              
%                     xqmax=datasample(1:N0,1,'Replace',false);
%                     %xqmin=datasample(1:N0,1,'Replace',false);
%              
%                     %% Collect data Z at random 
% 
%                     xiiAmax = [gridx1, ones(length(gridx1),1)];
%                     xiiAmaxmin = [gridx1, ones(length(gridx1),1)];
%                     xiiAmin = [gridx1, ones(length(gridx1),1)];
%                     xiiAminmax = [gridx1, ones(length(gridx1),1)];
%                     xiiAmaxall = [gridx1, ones(length(gridx1),length(Species))];
%                     xiiAminall = [gridx1, ones(length(gridx1),length(Species))];
%                     
% 
%                        Covmax = repmat(Abs_Error_NL(maxsp(1)), length(indmax), 1);    
%                        Zmax=OvS0_NL(indmax,maxsp(1),xqmax); 
%                     
%                        Covmaxmin = repmat(Abs_Error_NL(maxsp(end)), length(indmax), 1);    
%                        Zmaxmin=OvS0_NL(indmax,maxsp(end),xqmax);    
% 
% 
%                        Covmin = repmat(Abs_Error_NL(minsp(end)), length(indmin), 1); 
%                        Zmin=OvS0_NL(indmin,minsp(end),xqmax) ;       
% 
% 
%                        Covminmax = repmat(Abs_Error_NL(minsp(1)), length(indmin), 1); 
%                        Zminmax=OvS0_NL(indmin,minsp(1),xqmax) ;    
%                        
%                        Covmaxminall = Abs_Error_NL(:)'; 
%                        Zmaxall=OvS0_NL(indmax,:,xqmax) ;       
%                        Zminall=OvS0_NL(indmin,:,xqmax) ;   
% 
%                      
% 
%                     %% update your matrix that have ponit at which you want to evaluate your posterior 
% 
%                     xiiAmax(:,2)= Zmax.*xiiAmax(:,2);
%                     xiiAmaxmin(:,2)= Zmaxmin.*xiiAmaxmin(:,2);
%                     xiiAmin(:,2)= Zmin.*xiiAmin(:,2);
%                     xiiAminmax(:,2)= Zminmax.*xiiAminmax(:,2);
%                     xiiAmaxall(:,2:end)= Zmaxall.*xiiAmaxall(:,2:end);
%                     xiiAminall(:,2:end)= Zminall.*xiiAminall(:,2:end);
% 
% 
%                     %% Normalizing the joint density f 
% 
%                     mv=nan(6,N0);
% 
% %                 
% 
%                     for ii = 1:N0
% 
%                         mv(1,ii) =  normpdf( Zmax, OvS0_NL(indmax,maxsp(1),ii), Covmax);
%                         mv(2,ii) =  normpdf( Zmaxmin, OvS0_NL(indmax,maxsp(end),ii), Covmaxmin);
%                         mv(3,ii) =  normpdf( Zmin, OvS0_NL(indmin,minsp(end),ii), Covmin);
%                         mv(4,ii) =  normpdf( Zminmax, OvS0_NL(indmin,minsp(1),ii), Covminmax);
%                         mv(5,ii) =  mvnpdf( Zmaxall, OvS0_NL(indmax,:,ii), diag(Covmaxminall).^2);
%                         mv(6,ii) =  mvnpdf( Zminall, OvS0_NL(indmin,:,ii), diag(Covmaxminall).^2);
% 
%                     end
% 
%                     Imax = (1/N0)*sum(mv(1,:));  
%                     Imaxmin = (1/N0)*sum(mv(2,:)); 
%                     Imin = (1/N0)*sum(mv(3,:)); 
%                     Iminmax = (1/N0)*sum(mv(4,:)); 
%                     Imaxall = (1/N0)*sum(mv(5,:)); 
%                     Iminall = (1/N0)*sum(mv(6,:)); 
% 
%                     %% Calculate joint density at W and the selected Z 
%                     
%                     [~,~,bwp] = ksdensity(PAR(Sens(i),:)');
%                     
%                     [fmax,xiiAmax] = ksdensity(Density_Amax,xiiAmax,'bandwidth',([bwp Covmax]));
%                     [fmaxmin,xiiAmaxmin] = ksdensity(Density_Amaxmin,xiiAmaxmin,'bandwidth',([bwp Covmaxmin]));
%                     [fmin,xiiAmin] = ksdensity(Density_Amin,xiiAmin,'bandwidth',([bwp Covmin]));
%                     [fminmax,xiiAminmax] = ksdensity(Density_Aminmax,xiiAminmax,'bandwidth',([bwp Covminmax]));
%                     
%                     
%                     [fmaxall,xiiAmaxall] = mvksdensity(Density_Amaxall,xiiAmaxall,'bandwidth',[bwp Covmaxminall]);
%                     [fminall,xiiAminall] = mvksdensity(Density_Aminall,xiiAminall,'bandwidth',[bwp Covmaxminall]);
%                     %f = mvksdensity(Density_WZ,xii,'bandwidth',([bwp Cov(:)']));
% 
%                     %% Posterior is 
% 
%                      Postmax = fmax/Imax;
%                      Postmaxmin = fmaxmin/Imaxmin;
%                      Postmin = fmin/Imin;
%                      Postminmax = fminmax/Iminmax;
%                      Postmaxall = fmaxall/Imaxall;
%                      Postminall = fminall/Iminall;
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            %% Potting Prior and Posterior Distrubition of W after collecting data Z
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%            subplot(4,4,i)
%            hold on
% 
%                     
%                     %% Prior density of W 
%                     
%                     
%                     [f1,gridx1] = ksdensity(PAR(Sens(i),:)',gridx1);
%                     plot(gridx1,f1,'k+', 'LineWidth',2,'MarkerSize',8);
%                     %xlim([Comb(p)-10 100])
%              
%                     hold on
% 
% 
%                     % Posterior density of W
%        
%                     plot(gridx1,Postmax,'r', 'LineWidth',5);
% 
%                     hold on
% 
%                     
%                     % Posterior density of W
%        
%                     plot(gridx1,Postmaxmin,'r+', 'LineWidth',2,'MarkerSize',3);
%                     
%                     hold on            
%                     
%                     
%                     % Posterior density of W
%        
%                     plot(gridx1,Postminmax,'b+', 'LineWidth',0.5);
%                     
%                     hold on
%                     
%                     % Posterior density of W
%        
%                     plot(gridx1,Postmin,'b', 'LineWidth',1);
%                     
%                     hold on       
%                     
% 
%                     
%                     plot(gridx1,Postmaxall,'color',[0 0.6 0.3], 'LineWidth',3);
% 
%                     hold on
%                     
%                     plot(gridx1,Postminall,'y+', 'LineWidth',2,'MarkerSize',1);
%                     
%                     hold on
% 
% 
% 
% 
%                    xlabel('par value','fontweight','bold','fontsize',16,'Fontname','IPAPMincho'); 
%                    ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
% 
% 
%            
%                    annotation('textbox', [0.5 0.9 0.1 0.1],'String', 'Posterior P(\theta/z*)','fontsize',20,'Fontname','IPAPMincho', 'HorizontalAlignment', 'center')
%             
%         end                
%                             
%                  
%                    h_legend=legend('Prior','Post max from high info','Post min from high info','Post max from low info','Post min from low info','Post maxall from high info','Post maxall from low info');
%                    set(h_legend,'FontSize',6)
%                             
%                             
%                             
%                             
%                             
%                             
%                             
%                             
%              
% 
%        figure(600)
%        clf
%               for i=1:16
%            
%                     subplot(4,4,i)
%            
%                     ksdensity([PAR(Sens(i),:)' M_meanmax(:,1)], 'PlotFcn','contour'); 
%                     hold on
%                     plot(PAR(Sens(i),:)',  M_meanmax(:,1),'k.','MarkerSize',2);
%                     
%                     xlabel(sprintf('para%d',i),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
%                     ylabel(Ynames(Species(maxsp(1))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
% 
%               end
%               suptitle('high info from High mutual information');
%               
%        figure(601)
%        clf
%               for i=1:16
%            
%                     subplot(4,4,i)
%            
%                     ksdensity([PAR(Sens(i),:)' M_meanmaxmin(:,1)], 'PlotFcn','contour'); 
%                     hold on
%                     plot(PAR(Sens(i),:)',  M_meanmaxmin(:,1),'k.','MarkerSize',2);
%                     
%                     xlabel(sprintf('para%d',i),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
%                     ylabel(Ynames(Species(maxsp(end))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
% 
%               end
%               suptitle('low info from High mutual information');
% 
%               
%        figure(700)
%        clf
%               for i=1:16
%            
%                     subplot(4,4,i)
%            
%                     ksdensity([PAR(Sens(i),:)' M_meanmin(:,1)], 'PlotFcn','contour'); 
%                     hold on
%                     plot(PAR(Sens(i),:)',  M_meanmin(:,1),'k.','MarkerSize',2);
%                     
%                     xlabel(sprintf('para%d',i),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
%                     ylabel(Ynames(Species(minsp(end))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
% 
%                     
%                     
%               end
%        
%               suptitle('Low info from Low mutual information');
%               
%               
%        figure(701)
%        clf
%               for i=1:16
%            
%                     subplot(4,4,i)
%            
%                     ksdensity([PAR(Sens(i),:)' M_meanminmax(:,1)], 'PlotFcn','contour'); 
%                     hold on
%                     plot(PAR(Sens(i),:)',  M_meanminmax(:,1),'k.','MarkerSize',2);
%                     
%                     xlabel(sprintf('para%d',i),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
%                     ylabel(Ynames(Species(minsp(1))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
% 
%                     
%                     
%               end
%        
%               suptitle('High info from Low mutual information');
% 
%        
%        
% figure(201)
% clf;
% 
% subplot(3,1,1)
% 
%     bar(dataMI(:,datasample(1:N0,1,'Replace',false)),'FaceColor',[0.5 0.9 0.13]);
%    
%     hold on
%     
% subplot(3,1,2)
% 
%     bar(dataMI(:,datasample(1:N0,1,'Replace',false)),'FaceColor',[0.5 0.9 0.13]);
%    
%     hold on
%    
% subplot(3,1,3)
% 
%     bar(dataMI(:,datasample(1:N0,1,'Replace',false)),'FaceColor',[0.5 0.9 0.13]);
%     
%     ylabel(('I(\Theta,Y)'),'fontsize',22,'Fontname','IPAPMincho');
%     xlabel('Period' ,'fontsize',22,'Fontname','IPAPMincho');
%        
%        
%        
%        
%        
%        
%        
% 
% 
% 
% 
% 
% h = get(0, 'children');
% 
% 
%   if Ni==0  
%     for i=1:length(h)
%         saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_NL/Figures_L1/figure' num2str(i)], 'fig');
%     end
%     
%   elseif Ni==1
% 
%           for i=1:length(h)
%         saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_NL/Figures_L2/figure' num2str(i)], 'fig');
%           end
%     
%   elseif Ni==2
%           for i=1:length(h)
%         saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_NL/Figures_L3/figure' num2str(i)], 'fig');
%           end
% 
%   end
% toc
% datetime('now')
%  end
%  
%  
  display('=========================================== Non_Lactating KDE approach ===================================================')

% -------------------------------------------------------------------------
% Workflow 2: non-lactating KDE design
% -------------------------------------------------------------------------
% Goal: evaluate which measurements are most informative about a selected
% parameter in the non-lactating scenario. This section uses KDE-based
% density estimates and posterior plots for the parameter-focused BED result.

clear all

close all


display('----Sampling-----')

%% State variable of the model

Ynames =  { 'GnRhH'; 'GnRhP';'FSHP'; 'FSH'; 'LHP';  'LH'; 'Foll'; 'PGF';   'CL';    'P4';  'E2'; .....
           'INH'; 'ENZ'; 'OXT';  'IOF'; 'IGF';'Ins'; 'Glu'; 'Fat';   'Glucose in Liver'; 'Stored Glucose'};
       
%% Defining Species to include in simulation

        % str='n' selects the non-lactating simulation branch in the model.
        str='n'; 

        Species= 1:21; % Species to simulate in Metabolism (Species 1 and 2) and Reproduction model (the rest) 
        
        %Species= Species([4,8,10]);
        Species= Species([4,8,10,11,12,16,17,18]);
        
        %[~, spel]=ismember([4,8,10,11,12,16,17,18],Species);


        
        
for Ni=[1 ]
    
    sprintf('-------- Ni=%d ------',Ni)

        
        %% Number of samples

        % Smaller ensemble than workflow 1 because this KDE/posterior block
        % is computationally expensive.
        N0= 2000;

        %% Argument to pass to the solver

        frac2=0.08; 
        STEP=1/6;
        day=100;
        daystr = 90;

                

        %% Reading Paramaters

        par=BovSys_para_dexa_v3();
        % Parameter selected as the BED target in this workflow.
        parnum=48;

        
        %% Sampling theta from Uniform dis...

        % Defining the boundary of parameters interval
        
        % Non-lactating parameter uncertainty interval used for the KDE
        % design calculation.
        devS = 0.05;
        
        %Sens= [12    50    92    40    48    49    29    34   101    58 ...
         %   57    38    60    81    51    79    72   100    76    88    ...
         %   74    46    47    55    41    87    98    62    35    36    33    56];
    
    
        aRep=nan(length(par),1);
        bRep=nan(length(par),1);

        %NotSens = setxor(1:length(par),Sens);

        aRep = par - par*devS;
        bRep = par + par*devS;
        
        %devNotS = 0.1;
        
        %aRep(NotSens)=par(NotSens) - par(NotSens)*devNotS;
        %bRep(NotSens)=par(NotSens) + par(NotSens)*devNotS;      
        

        %% Check if parameters are negative

        if any(aRep<0)
            aRep(find(aRep<0))=0;
        end

        %% Sampling from parameters from uniform distribution ...
        
        repa=repmat(aRep,1,N0);
        repb=repmat(bRep,1,N0);
    
        thetaRep0 = unifrnd(repa,repb);

        %% Run the standard model to calculate the error ...
        
                     RefP4=1;
    
                     [T,YRepS,~]=BovSys_run_v3_baseline(par,frac2,day,daystr,RefP4,str,STEP);
                                                           
                     RefP4=YRepS(((day-daystr)*(1/STEP)):(day*(1/STEP)),10);

                     Out = YRepS;
          
                     Abs_Error_NL = 0.1*mean(Out(:,Species),1);
                                          
                     save('Abs_Error_NL.mat','Abs_Error_NL')
                     
                     %str= input('Lactating, l or nonlactating, n: ','s');   % switch back to lactating case!!!
               
                                          
                                                figure(111)

                                                for i =[1 4]
                                                    subplot(2,2,i)
                                                    plot(T,YRepS(:,i+6),'k')
                                                    hold on
                                                end
%                      hold on
%                      plot(T,0.4*YRepS(:,10),'r')
       
                    
                    %% Nonlactating case
                    
                    % Candidate sampling days for the non-lactating design.
                    tpoint = 54:1:89;  

                    Comb=reshape(tpoint,length(tpoint),1);

                    numComb=size(Comb,1);
        
                    
        
        %% All possible Combination of pair-Species mutual information
        

                    CombZ = nchoosek(1:length(Species),2);

                    numCombZ = size(CombZ,1);

                    mvML=nan(numCombZ,1);
    

    

        %% Run the model for different theta (smapled parameters)
        
                    %% First define matrix for input simulation

                    OvS0_NL  = zeros(length(tpoint),(length(Species))) ;
                    PAR  = nan(length(par),N0) ;


                    %OVmat_NL = nan(N0,1);



tic
        
display('----Generate models -----')


       
       %colorlist={'k','r','b','g'};
      
       % Simulate the non-lactating ensemble and retain both model outputs
       % and accepted parameter vectors. Failed simulations are filtered out
       % before information calculations.
       parfor (h=1:N0,16)
            
                    %% Run the model for different theta (smapled parameters)

                     [T,Y0,OvT0]=BovSys_run_v3_baseline(thetaRep0(:,h),frac2,day,daystr, RefP4,str,STEP);
                     
                             %figure(1000)
        
                             %plot(T,Y0(:,10),'color',colorlist{h})
                             %hold on
                             
                     Out = Y0;
                     
                     %% check if the model return nan

                     if isnan(OvT0)

                             OvS0_NL(:,:,h)=nan(length((tpoint)),(length(Species))) ;                       

                     else

                             OvS0_NL(:,:,h) = Out((tpoint)*(1/STEP),Species);
                             %OVmat_NL(h,1)  = OvT0;
                             PAR(:,h) =thetaRep0(:,h);
                   
                     end                 
                        
       end
       
       
       PAR=reshape(PAR(~isnan(PAR)),length(par),length(PAR(~isnan(PAR)))/length(par));


       %% Select the successfully simulated samples and reject the nan 
       
                    [~, ~, dim]= ind2sub(size(OvS0_NL),find(isnan(OvS0_NL)));

                    if ~isempty(dim)
                        
                        OvS0_NL(:,:,unique(dim))=[];
                        
                        if Ni==0  
                            
                              %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL_KDE/Data_1/Abs_Error_NL.mat','Abs_Error_NL')
                              %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL_KDE/Data_1/OvS0_NL.mat','OvS0_NL')


                          elseif Ni==1

                              %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL_KDE/Data_2/Abs_Error_NL.mat','Abs_Error_NL')
                              %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL_KDE/Data_2/OvS0_NL.mat','OvS0_NL')


                          elseif Ni==2
                              
                              %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL_KDE/Data_3/Abs_Error_NL.mat','Abs_Error_NL')
                              %save('/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Data_NL_KDE/Data_3/OvS0_NL.mat','OvS0_NL')

                         end
                        

                        display('............')
     
                        %OVmat_NL=OVmat_NL(~isnan(OVmat_NL));

                        %save('OVmat_NL.mat','OVmat_NL')
                        
                        %% update the number of samples

                        N0 = N0 - length(unique(dim));

                    else
                        
                        display('...llll....')

                        save('OvS0_NL.mat','OvS0_NL');
                        size(OvS0_NL)
                        
                        %save('OVmat_NL.mat','OVmat_NL')

        
                    end

                    AcceptedOutput=sprintf('Selected samples N0=%d ... Number of Rejected samples are: %d', N0, length(unique(dim)))

toc
datetime('now')
% 
%         %% here we start calculating pdfs:  ---------------------------------------------------------------------------------


       
%---------------------------------------------------- 1 -------------------------------------------------------------------  


        
            Info_Dtime = nan(1,numComb);
            Info_TS = nan(numComb,length(Species));
            Info_WZ_parnum_Ent = nan(numComb,length(Species));
            Info_ZZ = nan(numComb,numCombZ);
            Info_ZZ_Ent = nan(numComb,numCombZ);
            Info_Z_Ent= nan(numComb,length(Species));
            dataMI= nan(numComb,N0);

            
            Entropypara=nan(size(PAR,1),1);
            
            
            Entw=zeros(size(PAR,1),1);
            
            % Prior entropy proxy for each sampled parameter. The selected
            % target parameter is later compared with measurement-informed
            % posterior densities.
            
            for i=1:N0
            
                Entw  =  Entw  +  log(PAR(:,i));
                
            end

            Entropypara(:,1) = - (1/N0)*Entw;


            
            
tic


            for ind=1:numComb
                
                datetime('now')
                sprintf('Index ind=%d', ind)
                sprintf('N0=%d', N0)
                
                subComb = Comb(ind,:);
                [~,a]=ismember(subComb,tpoint);
                indT=a(a~=0);
                
                Info_ns=0;
                Info_ns_ind=zeros(1,length(Species));
                Ent_ns_WZ_parnum=zeros(1,length(Species));
                ns_ZZ=zeros(1,numCombZ);
                ns_Z=zeros(1,length(Species));
                Ent_ns_ZZ=zeros(1,numCombZ);
                Ent_ns_Z=zeros(1,length(Species));
                
                % For each candidate sampling day, estimate information in
                % all species together, each species separately, and every
                % species pair. The same ensemble is also used to evaluate
                % information about the chosen parameter parnum.

                for i=1:N0

                             %% Evaluate Sample at f(w) dist using KDE    

                                        
                                            %% Evaluate Sample at f(w) dist using KDE    
                                            BW=[];
                                            
                                            for u=1:size(PAR,1)

                                                [~,~,bw] = ksdensity(PAR(u,:)');   

                                                BW= [BW bw];
                                            
                                            end

                                            % p(theta): multivariate KDE
                                            % approximation for the sampled
                                            % parameter vector.
                                            mv=nan(1,N0);
 
                                        for ii = 1:N0
                                        
                                            mv(1,ii) =  mvnpdf(PAR(:,i)', PAR(:,ii)', diag(BW).^2);
                                            
                                        end
                                        
                                        out_w = (1/N0)*sum(mv);


                                %% Perturbate the selected model output   

                                % Synthetic measurement Z for this design:
                                % model output plus Gaussian measurement
                                % noise, truncated at zero.
                                
                                Measurements =  normrnd(OvS0_NL(indT,:,i),repmat(Abs_Error_NL, size(OvS0_NL(indT,:,i),1),1));

                                Measurements(Measurements<0)=0;

                                Cov = repmat(Abs_Error_NL, length(indT), 1);
                                
                                %% Calculate the pdf
                                        
                                if length(Measurements)==1

                                   mv=nan(1,N0);
                                   
                                   for ii = 1:N0

                                        mv(1,ii)=  normpdf( Measurements(:)', reshape(OvS0_NL(indT,:,ii),1,length(Species)*length(indT)), diag(Cov(:)'));
                                   
                                   end
                                            
                                else
                                    
                                   mv=nan(1,N0);
                                    
                                   for ii = 1:N0

                                        mv(1,ii) = mvnpdf( Measurements(:)', reshape(OvS0_NL(indT,:,ii),1,length(Species)*length(indT)), diag((Cov(:)').^2 ));
                                   end

                                end
                                        
                                out_z = (1/N0)*sum(mv);                                


                                
                                        %% Evaluate Sample at Joint dist f(w,z) for all z (species)
                                        
                                        % Joint density of all sampled
                                        % parameters and all selected
                                        % measurements.
                                        wz=[PAR(:,i)' Measurements(:)'];
                                        
                                        mv=nan(1,N0);
   
                                        for ii = 1:N0
                                       
                                            mv(1,ii) =  mvnpdf(wz, [PAR(:,ii)' reshape(OvS0_NL(indT,:,ii),1,length(Species)*length(indT))],  diag(([BW Cov(:)']).^2));
                                        
                                        end
                                        
                                        out_wz = (1/N0)*sum(mv);   
                                        
                                        %% individual mutual information
                                        
                                        out_wz_ind = nan(1,length(Species));
                                        out_wz_parnum_ind = nan(1,length(Species));
                                        out_wz_ind2 = nan(1,length(Species));
                                                            

                                        for iii=1:length(Species)

                                                                %% Evaluate Sample at Joint dist f(w,z) for every z

                                                                Measure_oneS =  Measurements(:,iii);

                                                                wz=[PAR(:,i)' Measure_oneS(:)'];
                                                                
                                                                wz_parnum=[PAR(parnum,i)' Measure_oneS(:)'];

                                                                Cov = repmat(Abs_Error_NL(iii), length(indT), 1);

                                                % Compare the full
                                                % parameter-vector joint
                                                % density with the selected
                                                % parameter-only joint
                                                % density.
                                                mvWZ=nan(1,N0);
                                                mvWZ_parnum=nan(1,N0);
                                                
                                                for ii = 1:N0

                                                                mvWZ(1,ii)=mvnpdf(wz, [PAR(:,ii)' reshape(OvS0_NL(indT,iii,ii),1,length(indT))], diag(([BW Cov(:)']).^2) );
                                                                
                                                                mvWZ_parnum(1,ii)=mvnpdf(wz_parnum, [PAR(parnum,ii)' reshape(OvS0_NL(indT,iii,ii),1,length(indT))], diag(([BW(parnum) Cov(:)']).^2) );

                                                end
                                                                out_wz_ind(1,iii) = (1/N0)*sum(mvWZ);
                                                                
                                                                out_wz_parnum_ind(1,iii) = (1/N0)*sum(mvWZ_parnum);

                                                                %% Evaluate Sample at f(z) dist for every z  

                                                mvZ=nan(1,N0);

                                                for ii = 1:N0

                                                                mvZ(1,ii)=mvnpdf(Measure_oneS(:)', reshape(OvS0_NL(indT,iii,ii),1,length(indT)), diag((Cov(:)').^2) );
                                                end
                                                                out_wz_ind2(1,iii) = (1/N0)*sum(mvZ);
                                                                
                                                                Ent_ns_Z(1,iii) = ns_Z(1,iii) + log(out_wz_ind2(1,iii));

                            
                                                            
                                        end

                                        %% Mutual information of W - Z    
                                        
                                        Ent_ns_WZ_parnum(1,:) = Ent_ns_WZ_parnum(1,:) + log(out_wz_parnum_ind);
                                        dataMI(ind,i) = (log(out_wz) - log(out_w) - log(out_z));
                                        Info_ns = Info_ns + (log(out_wz) - log(out_w) - log(out_z));
                                        Info_ns_ind(1,:) = Info_ns_ind(1,:) + (log(out_wz_ind(1,:)) - log(out_w) - log(out_wz_ind2(1,:)));  

                                
                                        %% All possible Combination of couple mutual information
                                        
                                                                                
                                        for iii=1:numCombZ

                                                ZZ=reshape(Measurements(:,CombZ(iii,:)),1,numel(Measurements(:,CombZ(iii,:))));

                                                iCov=repmat(Abs_Error_NL(CombZ(iii,:)), length(indT), 1);

                                                CovZZ=diag((reshape(iCov,1,numel(iCov))).^2);

                                                ZZcom=nan(1,N0);

                                                for ii = 1:N0

                                                    ZZcom(1,ii) =  mvnpdf( ZZ, reshape(OvS0_NL(indT,CombZ(iii,:),ii),1, numel(OvS0_NL(indT,CombZ(iii,:),ii))), CovZZ);
                                                end

                                                mvML(iii,1)=(1/N0)*sum(ZZcom);

                                                %% Mutual information of pair-OVmat_NL Z-Z

                                                ns_ZZ(1,iii) = ns_ZZ(1,iii) + (log(mvML(iii,1)) - log(out_wz_ind2(1,CombZ(iii,1))) - log(out_wz_ind2(1,CombZ(iii,2))) );  
                                                
                                                Ent_ns_ZZ(1,iii) = ns_ZZ(1,iii) + log(mvML(iii,1));

                                        end                             
                                
                                
                                

                end

                Info_Dtime(1,ind)=(1/N0) * Info_ns;

                Info_TS(ind,:)=(1/N0) * Info_ns_ind;
                
                Info_ZZ(ind,:) = (1/N0)*ns_ZZ;
                
                Info_ZZ_Ent(ind,:) = - (1/N0)*Ent_ns_ZZ;
                
                Info_Z_Ent(ind,:) = - (1/N0)*Ent_ns_Z;
                
                Info_WZ_parnum_Ent(ind,:) = - (1/N0)*Ent_ns_WZ_parnum;
 
                
            end
            
            
            
            
            
   


%% Pigure 1


figure(200)
clf;

subplot(2,1,1)

    bar(tpoint, Info_Dtime,'FaceColor',[0.93 0.69 0.13]);

   
    hold on
    %set(gca,'xticklabel',tpoint,'fontsize',8);
    
    ylabel(('I(\Theta,Y)'),'fontsize',22,'Fontname','IPAPMincho');
    xlabel('Period' ,'fontsize',22,'Fontname','IPAPMincho');
    
    
    
figure(201)
clf;

numz=16;
numx=16;

if N0<numz   
    numz=N0;
    numx=N0;
end

samv = datasample(1:N0,numz,'Replace',false);



for i = 1:numx
    
        subplot(numx,1,i)

        bar(tpoint,dataMI(:,samv(i)),'FaceColor',[0.5 0.9 0.13]);
        %set(gca,'xticklabel',tpoint,'fontsize',8);
        hold on

end

    ylabel(('I(\Theta,Y)'),'fontsize',22,'Fontname','IPAPMincho');
    xlabel('Period' ,'fontsize',22,'Fontname','IPAPMincho');
    

figure(300)
%clf;

for np=1:numComb

    subplot(6,7,np)
    
            bar(Info_TS(np,:),'r');

    set(gca,'xticklabel',Ynames(Species),'fontsize',8);
    
    ylabel(('I(Tov,Y)'),'fontsize',18,'Fontname','IPAPMincho');
    %xlabel('Model components' ,'fontsize',18,'Fontname','IPAPMincho');
    title(sprintf('D%d',Comb(np,:)),'fontsize',18,'Fontname','IPAPMincho');
    %title(sprintf('D%d',length(tpoint),'fontsize',18,'Fontname','IPAPMincho');

end


    figure(400)
    
    for j=1:numComb

        subplot(6,7,j)
        bar(Info_ZZ(j,:))
        title(sprintf('D%d',Comb(j,:)),'fontsize',18,'Fontname','IPAPMincho');

    
    end



       
    %%-----------------------------------------------------------------------------------------------------
       
    % Selecting high mutual information

     rank=find(Info_Dtime==max(Info_Dtime));
     rankn=find(Info_Dtime==min(Info_Dtime));
     
     rank=rank(1);
     rankn=rankn(1);
     
     [~, maxsp]=sort(Info_TS(rank,:),'descend');
     [~, minsp]=sort(Info_TS(rankn,:),'descend');

                
                subComb = Comb(rank,:);
                [~,a]=ismember(subComb,tpoint);
                indmax=a(a~=0);
                

                subComb = Comb(rankn,:);
                [~,a]=ismember(subComb,tpoint);
                indmin=a(a~=0);
                
                


                            Amax=[];
                            Amin=[];
                            Amaxmin=[];
                            Aminmax=[];
                            Amaxall=[];
                            Aminall=[];
                            %Aicol=[];

                            

                            for j = 1:N0
                                
                                M_meanmax=cat(1,Amax,OvS0_NL(indmax,maxsp(1),j)'); 
                                Amax=M_meanmax;
                                
                                M_meanmaxmin=cat(1,Amaxmin,OvS0_NL(indmax,maxsp(end),j)'); 
                                Amaxmin=M_meanmaxmin;
                                
                                M_meanmin=cat(1,Amin,OvS0_NL(indmin,minsp(end),j)'); 
                                Amin=M_meanmin;
                                
                                M_meanminmax=cat(1,Aminmax,OvS0_NL(indmin,minsp(1),j)'); 
                                Aminmax=M_meanminmax;
                                
                                maxall=cat(1,Amaxall,OvS0_NL(indmax,:,j)); 
                                Amaxall=maxall;
                                
                                minall=cat(1,Aminall,OvS0_NL(indmin,:,j)); 
                                Aminall=minall;
                                
                                %sels=cat(1,Aicol,OvS0_NL(indmax,icol,j)); 
                                %Aicol=sels;
                            end
           
                            

                                    %% loading models
                                    
                                    

        
        xqmax=samv;
 
        figure(1000)
        clf;
        
            for j=1:numz



                 Density_Amax = [PAR(parnum,:)' Amax];    
                 Density_Amaxmin = [PAR(parnum,:)' Amaxmin];    
                 Density_Amin = [PAR(parnum,:)' Amin];    
                 Density_Aminmax = [PAR(parnum,:)' Aminmax];    

                 Density_Amaxall = [PAR(parnum,:)' Amaxall];    
                 Density_Aminall = [PAR(parnum,:)' Aminall]; 

                 %Density_Aicol = [PAR(parnum,:)' Aicol]; 


                 gridx1 = unique(PAR(parnum,:))';

                 gridx1 = linspace((min(gridx1) - min(gridx1)/60),(max(gridx1)+max(gridx1)/60),500)';

                 % Posterior grid for the selected parameter. The posterior
                 % is evaluated after conditioning on different hypothetical
                 % informative measurements.


                        %% Collect data Z at random 

                        xiiAmax = [gridx1, ones(length(gridx1),1)];
                        xiiAmaxmin = [gridx1, ones(length(gridx1),1)];
                        xiiAmin = [gridx1, ones(length(gridx1),1)];
                        xiiAminmax = [gridx1, ones(length(gridx1),1)];
                        xiiAmaxall = [gridx1, ones(length(gridx1),length(Species))];
                        xiiAminall = [gridx1, ones(length(gridx1),length(Species))];
                        %xiiAicol = [gridx1, ones(length(gridx1),length(icol))];



                           Covmax = repmat(Abs_Error_NL(maxsp(1)), length(indmax), 1);    
                           Zmax=OvS0_NL(indmax,maxsp(1),xqmax(j)); 

                           Covmaxmin = repmat(Abs_Error_NL(maxsp(end)), length(indmax), 1);    
                           Zmaxmin=OvS0_NL(indmax,maxsp(end),xqmax(j));    


                           Covmin = repmat(Abs_Error_NL(minsp(end)), length(indmin), 1); 
                           Zmin=OvS0_NL(indmin,minsp(end),xqmax(j)) ;       


                           Covminmax = repmat(Abs_Error_NL(minsp(1)), length(indmin), 1); 
                           Zminmax=OvS0_NL(indmin,minsp(1),xqmax(j)) ;    

                           Covmaxminall = Abs_Error_NL(:)'; 
                           Zmaxall=OvS0_NL(indmax,:,xqmax(j)) ;       
                           Zminall=OvS0_NL(indmin,:,xqmax(j)) ;   

                           %Covicol = Abs_Error_NL(icol); 
                           %Zicol=OvS0_NL(indmax,icol,xqmax(j)) ;       



                        %% update your matrix that have ponit at which you want to evaluate your posterior 

                        xiiAmax(:,2)= Zmax.*xiiAmax(:,2);
                        xiiAmaxmin(:,2)= Zmaxmin.*xiiAmaxmin(:,2);
                        xiiAmin(:,2)= Zmin.*xiiAmin(:,2);
                        xiiAminmax(:,2)= Zminmax.*xiiAminmax(:,2);
                        xiiAmaxall(:,2:end)= Zmaxall.*xiiAmaxall(:,2:end);
                        xiiAminall(:,2:end)= Zminall.*xiiAminall(:,2:end);
                        %xiiAicol(:,2:end)= Zicol.*xiiAicol(:,2:end);



                        %% Normalizing the joint density f 




                        mv=nan(9,N0);

    %                 

                        for ii = 1:N0

                            mv(1,ii) =  normpdf( Zmax, OvS0_NL(indmax,maxsp(1),ii), Covmax);
                            mv(2,ii) =  normpdf( Zmaxmin, OvS0_NL(indmax,maxsp(end),ii), Covmaxmin);
                            mv(3,ii) =  normpdf( Zmin, OvS0_NL(indmin,minsp(end),ii), Covmin);
                            mv(4,ii) =  normpdf( Zminmax, OvS0_NL(indmin,minsp(1),ii), Covminmax);
                            mv(5,ii) =  mvnpdf( Zmaxall, OvS0_NL(indmax,:,ii), diag(Covmaxminall).^2);
                            mv(6,ii) =  mvnpdf( Zminall, OvS0_NL(indmin,:,ii), diag(Covmaxminall).^2);

                            mv(7,ii) =  mvnpdf( Zmaxall(maxsp(1:(end-3)-1)), OvS0_NL(indmax,maxsp(1:(end-3)-1),ii), diag(Covmaxminall(maxsp(1:(end-3)-1))).^2);
                            mv(8,ii) =  mvnpdf( Zmaxall(maxsp(1:(end-3)-2)), OvS0_NL(indmax,maxsp(1:(end-3)-2),ii), diag(Covmaxminall(maxsp(1:(end-3)-2))).^2);
                            mv(9,ii) =  mvnpdf( Zmaxall(maxsp(1:(end-3)-3)), OvS0_NL(indmax,maxsp(1:(end-3)-3),ii), diag(Covmaxminall(maxsp(1:(end-3)-3))).^2);


    %                         mv(7,ii) =  mvnpdf( Zicol, OvS0_NL(indmax,icol,ii), diag(Covicol).^2);               
    %                         mv(8,ii) =  mvnpdf( Zicol(1:end-1), OvS0_NL(indmax,icol(1:end-1),ii), diag(Covicol(1:end-1)).^2);
    %                         mv(9,ii) =  mvnpdf( Zicol(1:end-2), OvS0_NL(indmax,icol(1:end-2),ii), diag(Covicol(1:end-2)).^2);

                        end

                        Imax = (1/N0)*sum(mv(1,:));  
                        Imaxmin = (1/N0)*sum(mv(2,:)); 
                        Imin = (1/N0)*sum(mv(3,:)); 
                        Iminmax = (1/N0)*sum(mv(4,:)); 
                        Imaxall = (1/N0)*sum(mv(5,:)); 
                        Iminall = (1/N0)*sum(mv(6,:)); 

                        Iicol = (1/N0)*sum(mv(7,:));
                        Iicol1 = (1/N0)*sum(mv(8,:));
                        Iicol2 = (1/N0)*sum(mv(9,:));

                        %% Calculate joint density at W and the selected Z 

                        % Estimate p(theta,z*) for each candidate
                        % measurement subset. These values are normalized
                        % below to give p(theta | z*).
                        [~,~,bwp] = ksdensity(PAR(parnum,:)');

                        [fmax,xiiAmax] = ksdensity(Density_Amax,xiiAmax,'bandwidth',([bwp Covmax]));
                        [fmaxmin,xiiAmaxmin] = ksdensity(Density_Amaxmin,xiiAmaxmin,'bandwidth',([bwp Covmaxmin]));
                        [fmin,xiiAmin] = ksdensity(Density_Amin,xiiAmin,'bandwidth',([bwp Covmin]));
                        [fminmax,xiiAminmax] = ksdensity(Density_Aminmax,xiiAminmax,'bandwidth',([bwp Covminmax]));


                        [fmaxall,xiiAmaxall] = mvksdensity(Density_Amaxall,xiiAmaxall,'bandwidth',[bwp Covmaxminall]);
                        [fminall,xiiAminall] = mvksdensity(Density_Aminall,xiiAminall,'bandwidth',[bwp Covmaxminall]);

                        [ficol,xiiAicol]   = mvksdensity(Density_Amaxall(:,[1 maxsp(1:(end-3)-1)+1]),xiiAmaxall(:,[1 maxsp(1:(end-3)-1)+1]),'bandwidth',[bwp Covmaxminall(maxsp(1:(end-3)-1))]);
                        [ficol1,xiiAicol1] = mvksdensity(Density_Amaxall(:,[1 maxsp(1:(end-3)-2)+1]),xiiAmaxall(:,[1 maxsp(1:(end-3)-2)+1]),'bandwidth',[bwp Covmaxminall(maxsp(1:(end-3)-2))]);
                        [ficol2,xiiAicol2] = mvksdensity(Density_Amaxall(:,[1 maxsp(1:(end-3)-3)+1]),xiiAmaxall(:,[1 maxsp(1:(end-3)-3)+1]),'bandwidth',[bwp Covmaxminall(maxsp(1:(end-3)-3))]);


    %                     [ficol,xiiAicol] = mvksdensity(Density_Aicol,xiiAicol,'bandwidth',[bwp Covicol]);                                     
    %                     [ficol1,xiiAicol1] = mvksdensity(Density_Aicol(:,1:end-1),xiiAicol(:,1:end-1),'bandwidth',[bwp Covicol(1:end-1)]);
    %                     [ficol2,xiiAicol2] = mvksdensity(Density_Aicol(:,1:end-2),xiiAicol(:,1:end-2),'bandwidth',[bwp Covicol(1:end-2)]);


                        %% Posterior is 

                         % Posterior density for the selected parameter
                         % after observing each synthetic measurement set.
                         Postmax = fmax/Imax;
                         Postmaxmin = fmaxmin/Imaxmin;
                         Postmin = fmin/Imin;
                         Postminmax = fminmax/Iminmax;
                         Postmaxall = fmaxall/Imaxall;
                         Postminall = fminall/Iminall;
                         Posticol = ficol/Iicol;
                         Posticol1 = ficol1/Iicol1;
                         Posticol2 = ficol2/Iicol2;

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %% Potting Prior and Posterior Distrubition of W after collecting data Z
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               subplot(4,4,j)
               hold on


                        %% Prior density of W 


                        [f1,gridx1] = ksdensity(PAR(parnum,:)',gridx1);
                        plot(gridx1,f1,'k+', 'LineWidth',2,'MarkerSize',8);
                        %xlim([Comb(p)-10 100])

                        hold on


                        % Posterior density of W

                        plot(gridx1,Postmax,'r', 'LineWidth',5);

                        hold on


                        % Posterior density of W

                        plot(gridx1,Postmaxmin,'r+', 'LineWidth',2,'MarkerSize',3);

                        hold on            


                        % Posterior density of W

                        plot(gridx1,Postminmax,'b+', 'LineWidth',0.5);

                        hold on

                        % Posterior density of W

                        plot(gridx1,Postmin,'b', 'LineWidth',1);

                        hold on       



                        plot(gridx1,Postmaxall,'color',[0 0.6 0.3], 'LineWidth',3);

                        hold on

                        plot(gridx1,Postminall,'y+', 'LineWidth',2,'MarkerSize',1);

                        hold on

                        plot(gridx1,Posticol,'color',[0.8 0.5 0], 'LineWidth',3);

                        hold on

                        plot(gridx1,Posticol1,'color',[0.3 0.9 0.5], 'LineWidth',2);

                        hold on

                       plot(gridx1,Posticol2,'color',[0.1 0.2 0.6], 'LineWidth',1);

                        hold on



                       xlabel('par value','fontweight','bold','fontsize',16,'Fontname','IPAPMincho'); 
                       ylabel(('pdf'),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');



                       annotation('textbox', [0.5 0.9 0.1 0.1],'String', 'Posterior P(\theta/z*)','fontsize',20,'Fontname','IPAPMincho', 'HorizontalAlignment', 'center')

            end                


                        idnames=Ynames(Species(maxsp(1:(end-3)-1)));
                        idnames1=Ynames(Species(maxsp(1:(end-3)-2)));
                        idnames2=Ynames(Species(maxsp(1:(end-3)-3)));


                       h_legend=legend('Prior','Post maxsp from high info','Post minsp from high info','Post maxsp from low info','Post minsp from low info','Post allsp from high info','Post allsp from low info',sprintf(' from max %s',idnames{:}) ...
                           ,sprintf(' from max %s',idnames1{:}),sprintf(' from max %s',idnames2{:}));
                       set(h_legend,'FontSize',6)     

                               

       figure(600)
       clf
              for i=1:length(Species)
           
                    subplot(2,4,i)
           
                    ksdensity([PAR(parnum,:)' maxall(:,i)], 'PlotFcn','contour'); 
                    hold on
                    plot(PAR(parnum,:)',  maxall(:,i),'k.','MarkerSize',2);
                    
                    xlabel(sprintf('para%d',parnum),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    ylabel(Ynames(Species(i)),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
                    
                    RHO = corr2(PAR(parnum,:)',maxall(:,i));

                    Com = Info_Z_Ent(indmax,i) + Entropypara(parnum,1) - Info_WZ_parnum_Ent(indmax,i);
                    
                    Corr = sign(Com)*sqrt(1 - exp(-2*abs(Com)));
                    
                    xl=get(gca,'xlim');
                    yl=get(gca,'ylim');
                    text(xl(1),yl(2),sprintf('I=%.3f, p=%.3f',Corr,RHO),'fontweight','bold')

              end
              
              suptitle('From high mutual information');
              
       figure(601)
       clf
              for i=1:length(Species)
           
                    subplot(2,4,i)
           
                    ksdensity([PAR(parnum,:)' minall(:,i)], 'PlotFcn','contour'); 
                    hold on
                    plot(PAR(parnum,:)',  minall(:,i),'k.','MarkerSize',2);
                    
                    xlabel(sprintf('para%d',parnum),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    ylabel(Ynames(Species(i)),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
                    
                    RHO = corr2(PAR(parnum,:)',minall(:,i));
                    
                    Com = Info_Z_Ent(indmin,i) + Entropypara(parnum,1) -  Info_WZ_parnum_Ent(indmin,i);
                    
                    Corr = sign(Com)*sqrt(1 - exp(-2*abs(Com)));
                    
                    xl=get(gca,'xlim');
                    yl=get(gca,'ylim');
                    text(xl(1),yl(2),sprintf('I=%.3f, p=%.3f',Corr,RHO),'fontweight','bold')

              end
              suptitle('From low mutual information');

              
%        figure(700)
%        clf
%               for i=1:16
%            
%                     subplot(4,4,i)
%            
%                     ksdensity([PAR(parnum,:)' M_meanmin(:,1)], 'PlotFcn','contour'); 
%                     hold on
%                     plot(PAR(parnum,:)',  M_meanmin(:,1),'k.','MarkerSize',2);
%                     
%                     xlabel(sprintf('para%d',i),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
%                     ylabel(Ynames(Species(minsp(end))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
% 
%                     
%                     
%               end
%        
%               suptitle('Low info from Low mutual information');
%               
%               
%        figure(701)
%        clf
%               for i=1:16
%            
%                     subplot(4,4,i)
%            
%                     ksdensity([PAR(parnum,:)' M_meanminmax(:,1)], 'PlotFcn','contour'); 
%                     hold on
%                     plot(PAR(parnum,:)',  M_meanminmax(:,1),'k.','MarkerSize',2);
%                     
%                     xlabel(sprintf('para%d',i),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
%                     ylabel(Ynames(Species(minsp(1))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
% 
%                     
%                     
%               end
%        
%               suptitle('High info from Low mutual information');
%               
              
              
       figure(800)
       clf
              for i=1:size(CombZ,1)
           
                    subplot(6,7,i)
           
                    ksdensity([Amaxall(:,CombZ(i,1)) Amaxall(:,CombZ(i,2))], 'PlotFcn','contour'); 
                    hold on
                    plot(Amaxall(:,CombZ(i,1)),  Amaxall(:,CombZ(i,2)),'k.','MarkerSize',2);
                    
                    xlabel(Ynames(Species(CombZ(i,1))),'fontweight','bold','fontsize',22,'Fontname','IPAPMincho');
                    ylabel(Ynames(Species(CombZ(i,2))),'fontweight','bold' ,'fontsize',22,'Fontname','IPAPMincho');
                    
                    RHO = corr2(Amaxall(:,CombZ(i,1)),Amaxall(:,CombZ(i,2)));
                     
                    Com = sum (Info_Z_Ent(indmax, CombZ(i,:))) -  Info_ZZ_Ent(indmax,i);
                    
                    Corr = sign(Com)*sqrt(1 - exp(-2*abs(Com)));
                    
                    xl=get(gca,'xlim');
                    yl=get(gca,'ylim');
                    text(xl(1),yl(2),sprintf('I=%.3f, p=%.3f',Corr,RHO),'fontweight','bold')

                    
                    
                  
              end
       
              suptitle('Dependency from high mutual information');


       
       
       
       
       figure(3000)
       
       bar(Entropypara)
       
       suptitle('Parameters Entropy');

       
       
       
       

    
    
    CombZi=ones(size(CombZ,1),3);
    
    CombZi(:,1)=1:size(CombZ,1);
    CombZi(:,2:3)=CombZ;
    
    figure(202)
    clf;
   
    text(0.5,0.5,num2str(CombZi,2))
    
    

       
       
       





h = get(0, 'children');


  if Ni==0  
    for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_NL/KDE/Figures_L1/figure' num2str(i)], 'fig');
    end
    
  elseif Ni==1

          for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_NL/KDE/Figures_L2/figure' num2str(i)], 'fig');
          end
    
  elseif Ni==2
          for i=1:length(h)
        saveas(h(i), ['/nfs/datanumerik/bzfomari/PhD/Project3_BED/BED_1_M_models_last version/Figures_NL/KDE/Figures_L3/figure' num2str(i)], 'fig');
          end

  end

  toc
  datetime('now')
end
 
    




toc
end
