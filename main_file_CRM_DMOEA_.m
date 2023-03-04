%in the name of god 
%Copyright (c) 2023, github (www.github.com)
% All rights reserved. Please read the "license.txt" for license terms.
% Project Title: A Novel Combinational Response Mechanism for Dynamic Multi-Objective Optimization
% Publisher: github (www.github.com)
% Developers: Zahra Aliniya and  Seyed Hossien Khsteh
% Contact Info:Z.aliniya@student.alzahra.ac.ir, Zahra.Aliniya.M@gmail.com

function CRM_DMOEA_main()
% load ('sp.mat')
%  for  tot=1:2
 clear all
close all
global ideal_point  t nPop  itrCounter TestProblem step window  CostFunction nVar VarMin VarMax numOfObj ns nt   Cr mp change_conter 

%  rng(0);
 run_t1=0;
 for  run_t=1:2

step =5;%nT  intensity
window =5;%Tt  ferequency  (Tt,nT)(5,10) (10,10)
itrCounter = 1;
TestProblem = 31; %FDA1
nPop = 100;
ns=round(nPop/2);%sample from source domain
nt=round(nPop);%sample from target domain
pool = round(nPop/2);
tour = 2;
mu = 7;% dc:10  distribution index for crossover 
mum = 10;%dm: 20  distrib0ution index for mutation
Cr=0.7; %0.9  probability perform crossove
mp=0.3;%probability perform mutation
move_plot=0;
value_move=0;
% figure('Color',[1,1,1]);
% maxIt = round(step*window);
%  maxIt = round(window*100);
[numOfObj, nVar, VarMin, VarMax] = objective_description_function();
%  VarMin(1,:)=0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ICA
%% Problem Definition
% CostFunction=@(x) Sphere(x);        % Cost Function
% nVar=5;             % Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
% VarMin=-10;         % Lower Bound of Variables
% VarMax= 10;         % Upper Bound of Variables
%% ICA Parameters
MaxIt=10;         % Maximum Number of Iterations
flag6=true;
% nPop=100;            % Population Size
nEmp=3;            % Number of Empires/Imperialists
alpha=1;            % Selection Pressure
beta=1.5;           % Assimilation Coefficient
pRevolution=0.1;   % Revolution Probability
mu=0.1;             % Revolution Rate
zeta=0.2;           % Colonies Mean Cost Coefficient
%% Globalization of Parameters and Settings
global ProblemSettings;
ProblemSettings.CostFunction=CostFunction;
ProblemSettings.nVar=nVar;
ProblemSettings.VarSize=VarSize;
ProblemSettings.VarMin=VarMin;
ProblemSettings.VarMax=VarMax;
ProblemSettings.numOfObj=numOfObj;
global ICASettings;
ICASettings.MaxIt=MaxIt;
ICASettings.nPop=nPop;
ICASettings.nEmp=nEmp;
ICASettings.alpha=alpha;
ICASettings.beta=beta;
ICASettings.pRevolution=pRevolution;
ICASettings.mu=mu;
ICASettings.zeta=zeta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************
ideal_point = zeros(numOfObj, 1)';
%% Initialization
vx=5;
t=vx;
% Initialize Empires
emp=CreateInitialEmpires();
%Archive
Archive=emp;
% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,2);
value_move=0;
%^^^^^^^^^^^^^^^^^^^^^^^^^^^main^^^^^^^^^^^^^^^^^^^^
change_conter=0;
memory=[];
pool=70;
tour=2;
flagg=false;
change_conter=1;
chromosome = [];
chromosome = initialize_variables1(chromosome, 1);
chromosome = non_domination_sort_mod1(chromosome, numOfObj, nVar);

     for itrCounter=1:round(step*window)
         itrCounter;
%     for itrCounter=1:30
     t =floor(itrCounter/window)/step;
  
    %change Detection
    flag_change=ChangeDetection(chromosome);
    if flag_change==true 
%         t;
     t2 =floor((itrCounter-window)/window)/step;
%             if      t>=0.5 && t<=1%%%JY3

%       if      t>=0.5 && t<=1%%%JY9
%                 if      t>=0.5 && t<=1%%%JY7
%                  if  t<=0.5%%%JY6

       % plot front solution
%  [MS(change_conter),SP1(change_conter),IGD1(change_conter),GD1(change_conter)]=PlotFrontSolution(PF_front1,
%  value_move,t)%jy6
 [MS(change_conter),SP1(change_conter),IGD1(change_conter),GD1(change_conter)]=PlotFrontSolution(PF_front1, value_move,t2);
GD1;
         change_conter=change_conter+1;
     
 value_move=value_move+0.2;
pause(0.1)    
        % Response machanism
        % Archive is all population from t-1 invironment
        % memory is the pareto front solution until now
        emp=ResponseMechanism(emp,Archive,memory);
        if change_conter>1
        %memory contain best solution until now
        memory=Update_Memory(memory,PF_front1,round(nPop/3));
        end
    
    
 %%%%%%%%%%%%%%%covert emp to chromosome
 
     nEmp=numel(emp);
         pop=[emp(1).Imp];
         for j=1:emp(1).nCol
         pop=[pop;emp(1).Col(j)];
         end
        for k=2:nEmp
         pop=[pop;emp(k).Imp];
         for j=1:emp(k).nCol
         pop=[pop;emp(k).Col(j)];
         end
        end
        
   %%%%%%%%%%%
   for i=1: numel(pop)
   % [[pop(i).Position],[pop(i).Cost],[pop(i).Front],[pop(i).crowded_dis]]
x(i,1:nVar+numOfObj+3)=[[pop(i).Position],[pop(i).Cost],[pop(i).Front],[pop(i).crowded_dis],[pop(i).violation]];
   end
%%%%%%%%%%%%%%%%
    chromosome=x;
    chromosome = non_domination_sort_mod1(chromosome, numOfObj, nVar);

   %  end%JY9
%         end%JY3
%                  end %jy6
%                 end
    end
        flagg=true;
 %%%%%%%%% %% ICA Main Loop  
% %  if MaxIt==20
% %      if flag6==true
% %      MaxIt=20;
% %      flag6=false;
% %      end
% %  else
% %      MaxIt=20;
% %  end
   for itr=1:MaxIt
%        tic;
    parent_chromosome = tournament_selection(chromosome, pool, tour);
    offspring_chromosome = ...
        genetic_operator(parent_chromosome, ...
        numOfObj, nVar, mu, mum, VarMin, VarMax, itrCounter);
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    clear temp
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : numOfObj+nVar) = ...
        offspring_chromosome;
    intermediate_chromosome = ...
        non_domination_sort_mod1(intermediate_chromosome, numOfObj, nVar);
    chromosome = replace_chromosome(intermediate_chromosome, numOfObj, nVar, nPop);
%        time_run1=toc

   end
  
   
   %%%%%%%%%convert    chromosome to emp
   
   emp=convert_chromosome_emp(chromosome);
   
%%%%%%%%%%%%%%%%%%   
             [Archive , PF_front1]=finalPFpreviouseChenge(emp);
             
   [o1,o2,o3,total_GD1(itrCounter)]=PlotFrontSolution1(PF_front1, value_move,t);

% t=t+5;
     end
      
GD1(find(GD1==0))=[];
IGD1(find(IGD1==0))=[];
SP1(find(SP1==0))=[];
MS(find(MS==0))=[];
% figure
% save  dmop3_GD_10_30   total_GD1
for  i=1:size(total_GD1,2)
    P_GD(1:2,i)=[total_GD1(1,i);i];
end
meanGD=mean(GD1);
StdGD=std(GD1);

meanIGD=mean(IGD1);
StdIGD=std(IGD1);
SP1=SP1;
meanSP=mean(SP1);
StdSP=std(SP1);
meanMS=mean(MS);
StdMS=std(MS);
run_t=run_t
disp(['GD_Metric=   '  num2str(   GD1)   ,    '  meanGD  =   ',num2str(meanGD),'   StdGD =   ',num2str(StdGD)]);
disp(['IGD_Metric=   '  num2str(   IGD1)   ,    '  meanIGD  =   ',num2str(meanIGD),'   StdIGD =   ',num2str(StdIGD)]);
disp(['SP1 Metric  =   ',num2str(SP1),'   meanSP = ',num2str(meanSP),'   StdSP  = ',num2str(StdSP)]);
disp(['MS Metric =  ',num2str(MS),'   meanMS  =  ',num2str(meanMS),'  StdMS   =  ',num2str(StdMS)]);
total_GD(run_t)=meanGD;
total_IGD(run_t)=meanIGD;
total_sp(run_t)=meanSP;
total_ms(run_t)=meanMS;
F10_P_GD=P_GD;
% save F10_P_GD
% plot(P_GD(2,4:end),P_GD(1,4:end),'k'); 

% if numel(IGD1)==step
%     run_t1=run_t1+1;
% T_IGD1(run_t1,1:step)=IGD1;
% T_SP1(run_t1,1:step)=SP1;
% T_MS(run_t1,1:step)=MS;
% meanT_IGD1= (T_IGD1)
% meanT_SP1= mean (T_SP1)
% meanT_MS=mean (T_MS)
% end

% MS=[];
% SP1=[];
% IGD1=[];
% run_t
end
% total_IGD_mean=-mean(total_IGD)
% total_sp_mean=mean(total_sp)
% total_ms_mean=mean(total_ms)
% total_GD_mean=mean(total_GD);
% total_GD_std=std(total_GD);
total_IGD_mean=mean(total_IGD);
total_IGD_std=std(total_IGD);
total_sp_mean=mean(total_sp);
total_sp_std=std(total_sp);
total_ms_mean=mean(total_ms);
total_ms_std=std(total_ms);
% disp([   '  meanGD  =   ',num2str(total_GD_mean),'   StdGD =   ',num2str(total_GD_std)]);
disp([   '  meanIGD  =   ',num2str(total_IGD_mean),'   StdIGD =   ',num2str(total_IGD_std)]);
disp(['   meanSP = ',num2str(total_sp_mean),'   StdSP  = ',num2str(total_sp_std)]);
disp(['   meanMS  =  ',num2str(total_ms_mean),'  StdMS   =  ',num2str(total_ms_std)]);

            % Update Total Cost of Empires
%     emp=UpdateTotalCost(emp);
    
    % Inter-Empire Competition
%     emp=InterEmpireCompetition(emp);
    
    % Update Best Solution Ever Found


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    


% %     imp=[emp.Imp];
% %     best_imp = non_domination_sort_mod(imp, numOfObj, nVar);
% % %     [~, BestImpIndex]=min([imp.Cost]);
% % %     BestSol=imp(BestImpIndex);
% % %     
% %      BestSol=best_imp(1);
% % 
% %     % Update Best Cost
% %     
% %     BestCost(it,:)=BestSol.Cost;
% %     
% %     % Show Iteration Information
% %     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it,:))]);
% %     
% %  
