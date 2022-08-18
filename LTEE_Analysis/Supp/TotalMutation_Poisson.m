%% Load data and process
load Total_mutation_data.mat

NumMutations_s = transpose(0:10);
NumGenesWithMutations_s = PerStrain(1:11,1);
NumGenesWithMutations_SD_s = PerStrain(1:11,2);
PercentGenesWithMutations_s = PerStrain(1:11,3);
PercentGenesWithMutations_SD_s = PerStrain(1:11,4);

%% Initial mutation frequency plot
figure(1)
clf
hold on

% Plot data
bar(NumMutations_s,NumGenesWithMutations_s)

% Error
plot_error = errorbar(NumMutations_s,NumGenesWithMutations_s,NumGenesWithMutations_SD_s);
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';

% Aesthetic changes
title('Mutation frequency at generation 50,000')
xlabel('Number of muations')
ylabel('Number of genes')
xlim([-1,10])
ylim([0,60000])

hold off

%% Initial mutation frequency plot as percent of all genes
figure(2)
clf
hold on

TotalGenes = 12*4712;

% Plot data
bar(NumMutations_s,PercentGenesWithMutations_s)

% Error
plot_error = errorbar(NumMutations_s,PercentGenesWithMutations_s,PercentGenesWithMutations_SD_s);
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';

% Aesthetic changes
title('Mutation frequency at generation 50,000')
xlabel('Number of muations')
ylabel('Percent of all genes')
xlim([-1,10])
ylim([0,60000/TotalGenes])

hold off

%% Lump the last set of mutations together
figure(3)
clf
hold on

% Process data
NumMutations_lumped = [NumMutations_s(1:end-5); NumMutations_s(end)];
PercentGenesWithMutations_lumped = [PercentGenesWithMutations_s(1:end-5); sum(PercentGenesWithMutations_s(end-5:end))];
PercentGenesWithMutations_SD_lumped = PercentGenesWithMutations_SD_s(1:end-4);

% Plot data
bar(NumMutations_lumped,PercentGenesWithMutations_lumped)

% Error
plot_error = errorbar(NumMutations_lumped,PercentGenesWithMutations_lumped,PercentGenesWithMutations_SD_lumped);
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';

% Aesthetic changes
title('Mutation frequency at generation 50,000 lumped for >13')
xlabel('Number of muations')
ylabel('Percent of all genes')
xlim([-1,10])
ylim([0,60000/TotalGenes])

hold off

%% Fmincon fitting

lambda0 = 1;

F=@(x)poissfit(x,[NumMutations_s,PercentGenesWithMutations_s],PercentGenesWithMutations_SD_s.^2);

lambda_fmin = fmincon(F, lambda0);

%% Plot the Poisson fit versus measured data
figure(4)
clf
hold on
% Process data
Percent_Poisson = poisspdf(NumMutations_s,lambda_fmin);

% Plot data
plot(NumMutations_s,PercentGenesWithMutations_s,'ko',NumMutations_s, Percent_Poisson,'b-')

% Error
plot_error = errorbar(NumMutations_s,PercentGenesWithMutations_s,PercentGenesWithMutations_SD_s);
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';

% Aesthetic changes
title('Mutation frequency at generation 50,000 fit to Poisson distribution')
xlabel('Number of muations')
ylabel('Percent of all genes')
legend('Measured frequency', 'Calculated frequency','Location','northeast')
xlim([0,10])
ylim([0,60000/TotalGenes])


hold off
%% Plot as bar

figure(5)
clf
hold on

% Plot data
bar(NumMutations_s,[PercentGenesWithMutations_s, Percent_Poisson])

% Error
plot_error = errorbar(NumMutations_s-0.14,PercentGenesWithMutations_s,PercentGenesWithMutations_SD_s);
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';


% Aesthetic changes
title('Mutation frequency at generation 50,000')
xlabel('Number of muations')
ylabel('Percent of all genes')
legend('Measured frequency', 'Calculated frequency','Location','northwest')
xlim([-0.5,10])
ylim([0,60000/TotalGenes])

hold off

% Plot inset
NumMutations_inset = NumMutations_s(5:end);
PercentGenesWithMutations_inset = PercentGenesWithMutations_s(5:end);
Percent_Poisson_inset = Percent_Poisson(5:end);

axes('Position',[0.3, 0.3, 0.5, 0.5])
hold on

bar(NumMutations_inset,[PercentGenesWithMutations_inset,Percent_Poisson_inset])

plot_error = errorbar(NumMutations_inset-0.14,PercentGenesWithMutations_inset,PercentGenesWithMutations_SD_s(5:end));
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';

hold off
%print -painters -depsc output.eps

%% Kolmogorov-Smirnov test

s1 = PercentGenesWithMutations_s;
s2 = Percent_Poisson;
[h,p] = kstest2(s1,s2);

%% Plot average number of mutations of any gene
figure(6)
clf
hold on
bar(Avg_Mutations_DataCopy(:,1),Avg_Mutations_DataCopy(:,2));
poiss = poisspdf(1:12,3);
%bar([1:12]/12,poiss)

hold off

inset = Avg_Mutations_DataCopy(21:end,1);
Percent_inset = Avg_Mutations_DataCopy(21:end,2);

axes('Position',[0.3, 0.3, 0.5, 0.5])
hold on

bar(inset,Percent_inset)

hold off



%% Score function
function score = poissfit(lambda, dataIn, variance)
    NumMutations = dataIn(:,1);
    fitted_data = poisspdf(NumMutations,lambda);
    score = (dataIn(:,2) - fitted_data).^2;
    score(variance == 0) = [];
    variance(variance == 0) = [];
    score = sum(score./variance);
end


