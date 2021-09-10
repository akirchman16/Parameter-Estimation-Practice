clearvars;
close all;

% This code will use the Aguilera method of estimating the following
% reaction rate constants for a simple reaction: A <--> B.
%                       k_f = 0.2
%                       k_r = 0.12

%%
disp('Generating Artificial Data');
% Generate Artificial Data - SSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_f_Real = 1.41; %real value for k_f
k_r_Real = 0.92;    %real value for k_r

Datasets = 100;    %number of data sets generated
MaxTime = 5;   %time value that the simulations run until

% Memory Allocation
t_Artificial = zeros(Datasets,1);   %time
A_Artificial = zeros(Datasets,1);   %A population
B_Artificial = zeros(Datasets,1);   %B population
dt_Artificial = zeros(Datasets,1);  %time steps
j_Artificial = zeros(Datasets,1);   %reaction identifier

Set_Count = 0;  %counter for which data set is being run
while Set_Count < Datasets
    Set_Count = Set_Count+1;    %advances dataset counter
    
    t_Artificial(Set_Count,1) = 0;    %initial time
    A_Artificial(Set_Count,1) = 100;    %initial population of A
    B_Artificial(Set_Count,1) = 0;    %initial population of B
    
    Event = 0;  %counter for how many events have occured in this dataset
    while max(t_Artificial(Set_Count,:)) < MaxTime
        Event = Event+1;    %advances event counter
        
        p_Artificial = [k_f_Real*A_Artificial(Set_Count,Event), k_r_Real*B_Artificial(Set_Count,Event)];    %propensity functions
        Randoms_Artificial = [rand, rand]; %random numbers for SSA
        dt_Artificial(Set_Count,Event) = (1/sum(p_Artificial))*log(1/Randoms_Artificial(1)); %time step before reaction occurs
        
        if p_Artificial(1) >  Randoms_Artificial(2)*sum(p_Artificial)
            j_Artificial(Set_Count,Event) = 1; %forward reaction occurs
            A_Artificial(Set_Count,Event+1) = A_Artificial(Set_Count,Event)-1;  %decreases A population
            B_Artificial(Set_Count,Event+1) = B_Artificial(Set_Count,Event)+1;  %increases B population
        else
            j_Artificial(Set_Count,Event) = 2; %reverse reaction occurs
            A_Artificial(Set_Count,Event+1) = A_Artificial(Set_Count,Event)+1;  %increases A population
            B_Artificial(Set_Count,Event+1) = B_Artificial(Set_Count,Event)-1;  %decreases B population
        end
        
        t_Artificial(Set_Count,Event+1) = t_Artificial(Set_Count,Event)+dt_Artificial(Set_Count,Event);    %advances timer
    end
    
    figure(1);  %plot trajectory of populations from SSA
    subplot(2,2,1);
    scatter(t_Artificial(Set_Count,:),A_Artificial(Set_Count,:),1,'r','filled');
    hold on;
    scatter(t_Artificial(Set_Count,:),B_Artificial(Set_Count,:),1,'b','filled');
    
end

figure(1);
subplot(2,2,1)
xlabel('Time, t');
xlim([0 MaxTime]);
ylabel('Population');
ylim([0 inf]);
title('A \leftrightarrow B (Artificial Data)');
legend('A','B');

A_FinalState_Artificial = zeros(1,Datasets);    %memory allocation for final measured state of A
B_FinalState_Artificial = zeros(1,Datasets);    %memory allocation for final measured state of B
for i = 1:Datasets
    A_FinalState_Artificial(1,i) = A_Artificial(i,find(t_Artificial(i,:) <= MaxTime & t_Artificial(i,:) > 0,1,'last'));    %finds last measurement value of population A
    B_FinalState_Artificial(1,i) = B_Artificial(i,find(t_Artificial(i,:) <= MaxTime & t_Artificial(i,:) > 0,1,'last'));    %finds last measurement value of population B
end

figure(1); %histograms of final state
subplot(2,2,3); %A data
hA_Artificial = histfit(A_FinalState_Artificial,round(sqrt(Datasets)),'kernel');    %kernel fitting of histogram
hA_Artificial(1).FaceColor = [0.95,0.40,0.40]; %new color of histogram
hA_Artificial(2).Color = [0 0 0];
hold on;
xlabel('Population');
ylabel('Count');
title('Final A');
subplot(2,2,4); %B data
hB_Artificial = histfit(B_FinalState_Artificial,round(sqrt(Datasets)),'kernel');    %kernel fitting of histogram
hB_Artificial(1).FaceColor = [0.44,0.44,0.88]; %color of histogram
hB_Artificial(2).Color = [0 0 0];
hold on;
xlabel('Population');
ylabel('Count');
title('Final B');

% Center of Bins at Peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_kernelPeak = hA_Artificial(2).XData(hA_Artificial(2).YData == max(hA_Artificial(2).YData));    %maximum peak of the kernel fit for artificial A data
[A_val,A_idx] = min(abs(hA_Artificial(1).XData - A_kernelPeak));    %value of minimum difference between kernel peak and centers of bins
A_Max_BinCenter = hA_Artificial(1).XData(A_idx);    %center of the bin which is closest to where the peak of the kernel fit is
B_kernelPeak = hB_Artificial(2).XData(hB_Artificial(2).YData == max(hB_Artificial(2).YData));    %maximum peak of the kernel fit for artificial B data
[B_val,B_idx] = min(abs(hB_Artificial(1).XData - B_kernelPeak));    %value of minimum difference between kernel peak and centers of bins
B_Max_BinCenter = hB_Artificial(1).XData(B_idx);    %center of the bin which is closest to where the peak of the kernel fit is

disp('Done');
%% 
disp('Selecting Random Parameters');
% Parameter Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameter_Sets = 1000;   %number of parameter sets to generate
range_f = [0 3];    %range of values for k_f
range_r = [0 3];    %range of values for k_r

figure(1);  %scatter plot of all parameter values
subplot(2,2,2);
yline(k_r_Real,'r');
hold on;
xline(k_f_Real,'r');
xlabel('k_f');
ylabel('k_r');
xlim(range_f);
ylim(range_r);
title('Parameter Values');
box on;

k_f_Test = [];  %initial empty array of k_f parameters that will be tested
k_r_Test = [];  %initial empty array of k_r parameters that will be tested
for j = 1:Parameter_Sets
    k_f_Test = [k_f_Test, range_f(1)+(range_f(2)-range_f(1))*rand];
    k_r_Test = [k_r_Test, range_r(1)+(range_r(2)-range_r(1))*rand];
end
figure(1);  %scatter plot of all parameter values
subplot(2,2,2);
scatter(k_f_Test,k_r_Test,3,'k','filled');
hold on;

disp('Done');
%% 
disp('Running ODE Models and Testing Deterministic Precondition');
% ODE Model to Test Deterministic Precondition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_High = 1.05;   %upper limit multiplier for deterministic precondition
beta_Low = 0.95;    %lower limit multiplier for deterministic precondition
Good_Parameters = [];   %initializes list of parameter sets which satisfy precondition

dt = 0.001; %time step between Euler steps
Steps = MaxTime/dt; %number steps taken in each numerical solution
for k = 1:Parameter_Sets
    k_f = k_f_Test(k);  %forward rate constant to be examined
    k_r = k_r_Test(k);  %reverse rate constant to be examined
    Parameters = [k_f;k_r]; %vector of parameters
    
    X = zeros(2,Steps+1);   %allocation of vectors to represent populations [A;B]
    t = zeros(1,Steps+1);    %allocation of time array
    
    X(:,1) = [100;0];   %initial conditions
    
    for l = 1:Steps
        X(:,l+1) = [1-k_f*dt, k_r*dt; k_f*dt, 1-k_r*dt]*X(:,l);  %Euler step
        t(l+1) = t(l)+dt;   %advances time
    end
    
    ODE_FinalState = round(X(:,end));
    if (beta_Low*ODE_FinalState(1)) <= A_Max_BinCenter & A_Max_BinCenter <= (beta_High*ODE_FinalState(1)) %tests deterministic precondition for A
        if (beta_Low*ODE_FinalState(2)) <= B_Max_BinCenter & B_Max_BinCenter <= (beta_High*ODE_FinalState(2))   %tests deterministic precondition for B
            Good_Parameters = [Good_Parameters,Parameters]; %records that these were good parameters
            
            figure(1);
            subplot(2,2,2);     %changes color of good parameters to red
            scatter(Good_Parameters(1,:),Good_Parameters(2,:),3,'r','filled');
            hold on;
        end
    end
end
disp('Done');
%%
disp('Running SSA of "Reasonable" Parameters');

AveragesAmt = Datasets;  %number of simulations run at each parameter set (same as measurements of aritficial data)

% figure(2);
% subplot(2,2,2);
% scatter(Good_Parameters(1,:),Good_Parameters(2,:),3,'k','filled');
% hold on;
% yline(k_r_Real,'r');
% xline(k_f_Real,'r');
% title('Reasonable Parameters');
% xlabel('k_f');
% ylabel('k_r');
% xlim(range_f);
% ylim(range_r);
% box on;

% Stochastic Simulation of each Parameter that passed Precondition %%%%%%%%
Good_Sets = numel(Good_Parameters)/2;   %number of stochastic datasets to generate

% Memory Allocation
t_Precondition = zeros(Good_Sets,1);
A_Precondition = zeros(Good_Sets,1);
B_Precondition = zeros(Good_Sets,1);
dt_Precondition = zeros(Good_Sets,1);
j_Precondition = zeros(Good_Sets,1);
A_FinalState_Precondition = zeros(1,AveragesAmt);
B_FinalState_Precondition = zeros(1,AveragesAmt);
ObjFunction = zeros(1,Good_Sets);

for m = 1:Good_Sets
    k_f_Precondition = Good_Parameters(1,m);    %selects parameters for this set of stochastic data
    k_r_Precondition = Good_Parameters(2,m);
    hA_Precondition = [];
    t_Precondition = [];
    A_Precondition = [];
    B_Precondition = [];
    for Avg_Sets = 1:AveragesAmt    %runs multiple simulations of each acceptable parameter set
        t_Precondition(Avg_Sets,1) = 0;  %initial time
        A_Precondition(Avg_Sets,1) = 100;    %initial population of A
        B_Precondition(Avg_Sets,1) = 0;  %initial population of B

        Event_Precondition = 0;
        while max(t_Precondition(Avg_Sets,:)) < MaxTime
            Event_Precondition = Event_Precondition+1;  %advances the event counter

            p_Precondition = [k_f_Precondition*A_Precondition(Avg_Sets,Event_Precondition),k_r_Precondition*B_Precondition(Avg_Sets,Event_Precondition)]; %propensity function(s)
            Randoms_Precondition = [rand, rand];    %random numbers for SSA
            dt_Precondition(Avg_Sets,Event_Precondition) = (1/sum(p_Precondition))*log(1/Randoms_Precondition(1)); %determining the time step for event

            if p_Precondition(1) > Randoms_Precondition(2)*sum(p_Precondition)
                j_Precondition(Avg_Sets,Event_Precondition) = 1; %forward reaction occurs
                A_Precondition(Avg_Sets,Event_Precondition+1) = A_Precondition(Avg_Sets,Event_Precondition)-1;  %decreases A population
                B_Precondition(Avg_Sets,Event_Precondition+1) = B_Precondition(Avg_Sets,Event_Precondition)+1;  %increases B population
            else
                j_Precondition(Avg_Sets,Event_Precondition) = 2; %reverse reaction occurs
                A_Precondition(Avg_Sets,Event_Precondition+1) = A_Precondition(Avg_Sets,Event_Precondition)+1;  %increases A population
                B_Precondition(Avg_Sets,Event_Precondition+1) = B_Precondition(Avg_Sets,Event_Precondition)-1;  %decreases B population
            end

            t_Precondition(Avg_Sets,Event_Precondition+1) = t_Precondition(Avg_Sets,Event_Precondition)+dt_Precondition(m,Event_Precondition);    %advances the timer
        end
        
        A_FinalState_Precondition(Avg_Sets) = A_Precondition(Avg_Sets,find(t_Precondition(Avg_Sets,:) <= MaxTime & t_Precondition(Avg_Sets,:) > 0,1,'last')); %finds the last measured value for each simulation of A and B
        B_FinalState_Precondition(Avg_Sets) = B_Precondition(Avg_Sets,find(t_Precondition(Avg_Sets,:) <= MaxTime & t_Precondition(Avg_Sets,:) > 0,1,'last')); %based on data from parameters which passed deterministic precondition
    end
    hA_Precondition = histfit(A_FinalState_Precondition,round(sqrt(AveragesAmt)),'kernel');   %kernel fitting of this parameter set
    ObjFunction(m) = sum((hA_Artificial(1).YData-hA_Precondition(1).YData).^2);    %Objective Function
end

disp('Done');