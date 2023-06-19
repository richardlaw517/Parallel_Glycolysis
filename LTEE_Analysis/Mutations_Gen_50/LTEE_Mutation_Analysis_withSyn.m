%% Import data
tic
MutationsA = cell(2500,12);
MutationsB = cell(2500,12);
Mutations_Order_A = cell(1,12);
Mutations_Order_B = cell(1,12);

for StrainNum = 1:6
    Mutations_Order_A(StrainNum) = {['ara-',num2str(StrainNum),'_a']};
    Mutations_Order_A(StrainNum+6) = {['ara+',num2str(StrainNum),'_a']};
    Mutations_Order_B(StrainNum) = {['ara-',num2str(StrainNum),'_b']};
    Mutations_Order_B(StrainNum+6) = {['ara+',num2str(StrainNum),'_b']};
    
    MutationsA_entry_minus = readcell(['ara-',num2str(StrainNum),'a','50ks.csv']);
    MutationsA_entry_minus = cat(1,MutationsA_entry_minus,cell(2500 - length(MutationsA_entry_minus),1));

    MutationsB_entry_minus = readcell(['ara-',num2str(StrainNum),'b','50ks.csv']);
    MutationsB_entry_minus = cat(1,MutationsB_entry_minus,cell(2500 - length(MutationsB_entry_minus),1));

    MutationsA_entry_plus = readcell(['ara+',num2str(StrainNum),'a','50ks.csv']);
    MutationsA_entry_plus = cat(1,MutationsA_entry_plus,cell(2500 - length(MutationsA_entry_plus),1));

    MutationsB_entry_plus = readcell(['ara+',num2str(StrainNum),'b','50ks.csv']);
    MutationsB_entry_plus = cat(1,MutationsB_entry_plus,cell(2500 - length(MutationsB_entry_plus),1));
    
    MutationsA(:,StrainNum) = MutationsA_entry_minus;
    MutationsB(:,StrainNum) = MutationsB_entry_minus;
    MutationsA(:,StrainNum+6) = MutationsA_entry_plus;
    MutationsB(:,StrainNum+6) = MutationsB_entry_plus;
end
disp('Import done')
toc

%% Count number of mutations in each gene for each strain
tic
MutationsA_thenB = [MutationsA, MutationsB];
MutationsA_thenB(cellfun(@isempty,MutationsA_thenB))={'000'};
Uniques = unique(MutationsA_thenB);
Uniques = Uniques(2:end);
Unique_Counts = zeros(size(Uniques,1)-1,size(MutationsA_thenB,2));

for UniqueNum = 1:size(Uniques,1)
    Unique_Counts(UniqueNum,:) = sum(strcmp(MutationsA_thenB,Uniques(UniqueNum)));
end

Unique_Counts_A = Unique_Counts(:,1:12);
Unique_Counts_B = Unique_Counts(:,13:24);

disp('Mutations in genes counted')
toc

%% Add up the number of mutations in a gene across all 12 strains

Total_Mutations_Each_Gene = sum(Unique_Counts,2);
numMutations = [unique(Total_Mutations_Each_Gene(:));inf];
geneCounts = histcounts(Total_Mutations_Each_Gene(:),numMutations);
geneCounts = geneCounts.';
numMutations = [0; numMutations(1:end-1)];
geneCounts = [4386 - sum(geneCounts);geneCounts];

%% Plot the frequency of mutations of any gene
figure(1)
clf
hold on
bar(numMutations,geneCounts);
xlim([-1,51])
hold off

inset = numMutations(21:end);
Count_inset = geneCounts(21:end);

axes('Position',[0.3, 0.3, 0.5, 0.5])
hold on
xlim([19,51])
ylim([0,8])

bar(inset,Count_inset)

hold off

%print -painters -depsc FrequencyOfMutations.eps

%% Generate all combinations of the mutations in twelve strains (considering each gene from each strain as unique)
tic
a = [0,1]; % 0 = A, 1 = B
Combo_maker = combvec(a,a,a,a,a,a,a,a,a,a,a,a); % Binary list of combinations
NumCombos = size(Combo_maker,2);
Mutations_Combos = zeros(12*(size(Uniques,1)),NumCombos);
for ComboNum = 1:NumCombos % Makes all 2^12 = 4096 combinations
    Mutations_combo_entry = [];
    for StrainNum = 1:12 % Adds the number of mutations of 12 strains into single cell vector
        if Combo_maker(StrainNum,ComboNum) == 0
            Mutations_combo_entry = cat(1,Mutations_combo_entry,Unique_Counts_A(:,StrainNum));
        else
            Mutations_combo_entry = cat(1,Mutations_combo_entry,Unique_Counts_B(:,StrainNum));
        end
    end
    Mutations_Combos(:,ComboNum) = Mutations_combo_entry;
end

disp('Combinations generated')
toc

%% Count the number of genes with X mutations in the combinations
tic
MaxNumMutations = max(Mutations_Combos,[],'all');

NumGenesWithMutations = zeros(MaxNumMutations,NumCombos);

for MutationsNum = 1:MaxNumMutations
    NumGenesWithMutations(MutationsNum,:) = sum(Mutations_Combos==MutationsNum,1);
end

% Count the number of genes with zero mutations
TotalGenes = 12*4386;
NumZero = TotalGenes - sum(NumGenesWithMutations,1);
NumGenesWithMutations = [NumZero;NumGenesWithMutations];

disp('Number of mutations counted')
toc

%% Take the average and standard deviation across all combos
tic
Average_Num = mean(NumGenesWithMutations,2);
Average_Percent = Average_Num/(TotalGenes);
Std_Num = std(NumGenesWithMutations,0,2);
Std_Percent = Std_Num/(TotalGenes);

disp('Average and SD calculated')
toc

%% Fmincon fitting

num_Mutations = transpose(0:MaxNumMutations);

lambda0 = 1;

F=@(x)poissfit(x,[num_Mutations,Average_Percent],Std_Percent.^2);

lambda_fmin = fmincon(F, lambda0);

%% Plot fitted versus observed data

figure(2)
clf
hold on

PoissonFit = poisspdf(num_Mutations,lambda_fmin);

% Plot data
bar(num_Mutations,[Average_Percent, PoissonFit])

% Error
plot_error = errorbar(num_Mutations-0.14,Average_Percent,Std_Percent);
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
inset_Start = 5;

axes('Position',[0.3, 0.3, 0.5, 0.5])
hold on
xlim([3.5,10.5])
ylim([0,9E-4])

bar(num_Mutations(inset_Start:end),[Average_Percent(inset_Start:end),PoissonFit(inset_Start:end)])

plot_error = errorbar(num_Mutations(inset_Start:end)-0.14,Average_Percent(inset_Start:end),Std_Percent(inset_Start:end));
plot_error.Color = [0,0,0];
plot_error.LineStyle = 'none';

hold off
%print -painters -depsc output.eps

%% Kolmogorov-Smirnov test

s1 = Average_Percent;
s2 = PoissonFit;
[h,p] = kstest2(s1,s2);


%% Score function
function score = poissfit(lambda, dataIn, variance)
    NumMutations = dataIn(:,1);
    fitted_data = poisspdf(NumMutations,lambda);
    score = (dataIn(:,2) - fitted_data).^2;
    score(variance == 0) = [];
    variance(variance == 0) = [];
    score = sum(score./variance);
end