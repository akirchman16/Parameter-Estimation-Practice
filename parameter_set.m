function [k_1,k_2] = parameter_set(X)
    k_1 = [];
    k_2 = [];
    for i = 1:X
        k_1 = [k_1, 5*rand];    %selects random values for parameter set
        k_2 = [k_2, 5*rand];    %selects random values for parameter set
    end
    
    figure(1);
    scatter(k_1,k_2,3,'k','filled');    %scatterplot of parameters that will be tested
    xlabel('k_1');
    ylabel('k_2');
    box on;
end