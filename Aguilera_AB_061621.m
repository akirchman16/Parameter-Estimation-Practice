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
hA_Artificial = histfit(A_FinalState_Artificial,round(sqrt(Datasets)),'kernel');    %kernel fitting of histogram
hA_Artificial(1).FaceColor = [0.95,0.40,0.40]; %new color of histogram
hA_Artificial(2).Color = [0 0 0];
hold on;
xlabel('Final A Population');
ylabel('Probability');
title('A');
subplot(1,2,2); %B data
hB_Artificial = histfit(B_FinalState_Artificial,round(sqrt(Datasets)),'kernel');    %kernel fitting of histogram
hB_Artificial(1).FaceColor = [0.44,0.44,0.88]; %color of histogram
hB_Artificial(2).Color = [0 0 0];
hold on;
xlabel('Final B Population');
title('B');

% Center of Bins at Peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_kernelPeak = hA_Artificial(2).XData(hA_Artificial(2).YData == max(hA_Artificial(2).YData));    %maximum peak of the kernel fit for artificial A data
[A_val,A_idx] = min(abs(hA_Artificial(1).XData - A_kernelPeak));    %value of minimum difference between kernel peak and centers of bins
A_Max_BinCenter = hA_Artificial(1).XData(A_idx);    %center of the bin which is closest to where the peak of the kernel fit is
B_kernelPeak = hB_Artificial(2).XData(hB_Artificial(2).YData == max(hB_Artificial(2).YData));    %maximum peak of the kernel fit for artificial B data
[B_val,B_idx] = min(abs(hB_Artificial(1).XData - B_kernelPeak));    %value of minimum difference between kernel peak and centers of bins
B_Max_BinCenter = hB_Artificial(1).XData(B_idx);    %center of the bin which is closest to where the peak of the kernel fit is

% Parameter Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameter_Sets = 1000;   %number of parameter sets to generate
range_f = [0 5];    %range of values for k_f
range_r = [0 5];    %range of values for k_r

k_f_Test = [];  %initial empty array of k_f parameters that will be tested
k_r_Test = [];  %initial empty array of k_r parameters that will be tested
for j = 1:Parameter_Sets
    k_f_Test = [k_f_Test, range_f(1)+(range_f(2)-range_f(1))*rand];
    k_r_Test = [k_r_Test, range_r(1)+(range_r(2)-range_r(1))*rand];
end
figure(3);  %scatter plot of all parameter values
scatter(k_f_Test,k_r_Test,3,'k','filled');
xlabel('k_f');
ylabel('k_r');
xlim(range_f);
ylim(range_r);
title('Parameter Values');