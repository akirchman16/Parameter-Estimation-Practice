clearvars;
close all;

% This will provide an ODE model for the reaction between A<-->B. This will
% be used as the basis for the "Deterministic Precondition" as described by
% Aguilera, et al. An euler method will be used to continuously solve for
% the state of the system.

A(1) = 100; %initial population for A
B(1) = 0;   %initial population for B

k_f = 0.2;  %kinetic rate constant for a forward reaction
k_r = 0.12; %kinetic rate constant for a reverse reaction

MaxTime = 20;    %maximum time the simulation will run to
dt = 0.001; %time step between Euler steps

Steps = MaxTime/dt; %number of steps taken

X = zeros(2,Steps+1); %creation of zero vector for state vectors (memory allocation)
t = zeros(1,Steps+1); %memory allocation for time
X(:,1) = [A(1);B(1)]; %State Vector

for i = 1:Steps
    X(:,i+1) = [1-k_f*dt, k_r*dt; k_f*dt, 1-k_r*dt]*X(:,i); %Euler Method for reaction
    t(i+1) = t(i)+dt;   %time advancement
end

figure(1);
scatter(t,X(1,:),3,'r','filled');
hold on;
scatter(t,X(2,:),3,'b','filled');
xlabel('Time, t');
xlim([0 MaxTime]);
ylabel('Populations');
ylim([0 inf]);
legend('A','B');
title('A \leftrightarrow B');

FinalState = round(X(:,end));
disp(['A_f = ', num2str(FinalState(1))]);
disp(['B_f = ', num2str(FinalState(2))]);