clearvars;
close all;

% This will generate artificial data using the Gillespie Algorithm (SSA - 
% First Reaction Method) for a chemical reaction of A+B <--> C. This reaction
% has two parameters (rate constants): k_f (A+B-->C) and k_r (C-->A+B). It will
% track the total number of molecules in a sample at various times. Future
% codes will then use various parameter estimation methods to determine the
% given parameters:
%                       k_f = 0.2
%                       k_r = 0.12

k_f = 0.2;  %kinetic rate constant for forward reaction
k_r = 0.12;  %kinetic rate constant for reverse reaction

Iterations = 1000;    %number of experimental trials are done for data
MaxTime = 20;        %time value each iteration will run until

% Allocation of Variables
t = zeros(Iterations,1);
A = zeros(Iterations,1);
B = zeros(Iterations,1);

I_Count = 0;    %I_Count tracks which iteration the code is currently on
while I_Count < Iterations
    I_Count = I_Count+1;    %increases iteration counter
    
    t(I_Count,1) = 0;   %initial time
    A(I_Count,1) = 100; %initial number of A molecules
    B(I_Count,1) = 0; %initial number of B molecules
    
    Events = 0;     %Counts which event the current iteration is on
    while max(t(I_Count,:)) <= MaxTime
        Events = Events+1;  %increases events counter
        
        a = [k_f*A(I_Count,Events),k_r*B(I_Count,Events)];    %propensity functions
        R = [rand, rand];   %random numbers to be used for time and reaction selection
        TimeToNext = [(1/a(1))*log(1/R(1)),(1/a(2))*log(1/R(2))];   %calculates the time to next reaction for each reaction
        tau = min(TimeToNext);  %time until next reaction (first reaction method)
        j = find(TimeToNext == min(TimeToNext));    %selection of reaction
        
        if j == 1   %forward reaction occurs
            A(I_Count,Events+1) = A(I_Count,Events)-1;  %reaction decreases A count by 1
            B(I_Count,Events+1) = B(I_Count,Events)+1;  %reaction decreases B count by 1
        elseif j == 2   %reverse reaction occurs
            A(I_Count,Events+1) = A(I_Count,Events)+1;  %reaction increases A count by 1
            B(I_Count,Events+1) = B(I_Count,Events)-1;  %reaction increases B count by 1
        end
        
        t(I_Count,Events+1) = t(I_Count,Events)+tau;    %increases time by chosen interval
    end
    
    figure(1);  %plots populations of species over time
    scatter(t(I_Count,:),A(I_Count,:),1,'r','filled');
    hold on;
    scatter(t(I_Count,:),B(I_Count,:),1,'b','filled');
    xlabel('Time, t');
    xlim([0 MaxTime]);
    ylabel('Population');
    ylim([0 inf]);
    legend('A','B');
    box on;
end

figure(1);
title('A \leftrightarrow B');

NumMeasurements = 20;   %number of measurements taken in the experiment
Measurement_Step = MaxTime/NumMeasurements; %time steps between measurements

% Memory Allocation
Measurement_State_A = zeros(Iterations,NumMeasurements);
Measurement_State_B = zeros(Iterations,NumMeasurements);

t_Measurement = Measurement_Step*[1:1:NumMeasurements];    %records time for each measurement


for m = 1:Iterations
    for n = 1:NumMeasurements
        Measurement_State_A(m,n) = A(m,find(t(m,:) <= n*Measurement_Step & t(m,:) > 0,1,'last'));   %records the state of the system (A population) at each measurement
        Measurement_State_B(m,n) = B(m,find(t(m,:) <= n*Measurement_Step & t(m,:) > 0,1,'last'));   %records the state of the system (B population) at each measurement
    end
%     figure(2);  %plot of measurement data
%     scatter(t_Measurement,Measurement_State_A(m,:),3,'r','filled');
%     hold on;
%     scatter(t_Measurement,Measurement_State_B(m,:),3,'b','filled');
%     xlabel('Time');
%     ylabel('State(s)');
%     title('Measurement Data');
%     legend('A','B');
%     box on;
end

Final_State_A = reshape(Measurement_State_A(:,NumMeasurements),[1,Iterations]); %array of final state for A with each iteration
Final_State_B = reshape(Measurement_State_B(:,NumMeasurements),[1,Iterations]); %array of final state for B with each iteration

figure(3);
subplot(2,1,1); %histogram of A final state
histogram(Final_State_A,ceil(sqrt(Iterations)),'FaceColor','r');
xlabel('Final State (Population)');
ylabel('Count');
title('A');
subplot(2,1,2); %histogram of B final state
histogram(Final_State_B,ceil(sqrt(Iterations)),'FaceColor','b');
xlabel('Final State (Population)');
ylabel('Count');
title('B');