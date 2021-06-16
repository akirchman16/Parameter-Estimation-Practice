clearvars;
close all;

% This code will use the Aguilera method of estimating the following
% reaction rate constants for a simple reaction: A <--> B.
%                       k_f = 0.2
%                       k_r = 0.12

% Generate Artificial Data - SSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_f_Real = 0.2; %real value for k_f
k_r_Real = 0.12;    %real value for k_r

Datasets = 100;    %number of data sets generated
MaxTime = 20;   %time value that the simulations run until

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
        dt_Artificial(Set_Count,Event) = (1/sum(p_Artificial))*log(1/rand(1)); %time step before reaction occurs
        
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
    scatter(t_Artificial(Set_Count,:),A_Artificial(Set_Count,:),1,'r','filled');
    hold on;
    scatter(t_Artificial(Set_Count,:),B_Artificial(Set_Count,:),1,'b','filled');
    
end

figure(1);
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

figure(2); %histograms of final state
subplot(1,2,1); %A data
histogram(A_FinalState_Artificial,round(sqrt(Datasets)),'FaceColor','r','Normalization','probability');
hold on;
xlabel('Final A Population');
ylabel('Probability');
title('A');
subplot(1,2,2); %B data
histogram(B_FinalState_Artificial,round(sqrt(Datasets)),'FaceColor','b','Normalization','probability');
hold on;
xlabel('Final B Population');
title('B');

