%% Code to analyze amino acid sequence data, downloaded from RCSB Protein Data Bank

fileID = fopen('dna_seq_linkers.txt'); % Open text file with DNA/amino acid sequence (spaces between letters) HERE: Eco RpoD 
C = textscan(fileID,'%s'); % scan all letters
fclose(fileID);
DNA_raw=C{1,1}; % array of letters
%DNA=[DNA_raw{2:end,1}];

Linker_DNA = [DNA_raw{2,1},DNA_raw{5,1},DNA_raw{6,1},DNA_raw{7,1},DNA_raw{8,1},DNA_raw{9,1},DNA_raw{12,1},DNA_raw{15,1},DNA_raw{18,1},DNA_raw{21,1},DNA_raw{24,1}, DNA_raw{27,1},DNA_raw{30,1},DNA_raw{33,1},DNA_raw{36,1},DNA_raw{39,1},DNA_raw{42,1},DNA_raw{45,1},DNA_raw{48,1},DNA_raw{51,1}, DNA_raw{54,1}, DNA_raw{57,1}, DNA_raw{61,1}, DNA_raw{64,1}, DNA_raw{67,1}, DNA_raw{70,1}, DNA_raw{73,1}, DNA_raw{76,1}, DNA_raw{79,1}, DNA_raw{85,1}, DNA_raw{88,1}, DNA_raw{91,1}, DNA_raw{94,1} ];

expression = '\s';
split = regexp(Linker_DNA, expression, 'split');

splitdata_linker_DNA = regexp(split, '([A-Z])', 'tokens'); % using regular expressions to separate all letters into single cells

figure 
f = aacount(Linker_DNA, 'chart','bar');
title('Amino Acids Present in Disordered Linkers');
xlabel('Amino Acids');
saveas(gcf,'DNA_sequence_analysis_linkers_aas.pdf','pdf') 

% save for gnuplot
aa = aacount(Linker_DNA);

writetable(struct2table(aa), 'DNA_sequence_linker_aa_analysis.txt')
%% Make matrix where every column is a different linker
a = splitdata_linker_DNA;
m = numel(a);
k = cellfun(@numel,a);
Z = zeros(max(k),m);
DNA_sequences_Linkers_ALL = string(Z);
for i = 1:m
    DNA_sequences_Linkers_ALL(1:k(i),i) = a{i}(:);
end


%% Analysis of amino acids
%Index_R = find(contains(DNA_sequences_Linkers_ALL,'R'));
%idx = all(ismember(DNA_sequences_Linkers_ALL,'R','rows'))

% length
length_linkers = k;

% positive aa
[LiaR, LocbR] = ismember(DNA_sequences_Linkers_ALL,'R');

[LiaK, LocbK] = ismember(DNA_sequences_Linkers_ALL,'K');

[LiaH, LocbH] = ismember(DNA_sequences_Linkers_ALL,'H');

amount_positive = sum([LocbR; LocbK; LocbH]);

% negative aa
[LiaD, LocbD] = ismember(DNA_sequences_Linkers_ALL,'D');

[LiaE, LocbE] = ismember(DNA_sequences_Linkers_ALL,'E');


amount_negative = sum([LocbD; LocbE]);

% order promoting aa vs disorder promoting
[LiaI, LocbI] = ismember(DNA_sequences_Linkers_ALL,'I');

[LiaL, LocbL] = ismember(DNA_sequences_Linkers_ALL,'L');

[LiaV, LocbV] = ismember(DNA_sequences_Linkers_ALL,'V');

[LiaW, LocbW] = ismember(DNA_sequences_Linkers_ALL,'W');

[LiaY, LocbY] = ismember(DNA_sequences_Linkers_ALL,'Y');

[LiaF, LocbF] = ismember(DNA_sequences_Linkers_ALL,'F');

[LiaC, LocbC] = ismember(DNA_sequences_Linkers_ALL,'C');

[LiaN, LocbN] = ismember(DNA_sequences_Linkers_ALL,'N');

%disorder promoting
[LiaR, LocbR] = ismember(DNA_sequences_Linkers_ALL,'R');

[LiaG, LocbG] = ismember(DNA_sequences_Linkers_ALL,'G');

[LiaQ, LocbQ] = ismember(DNA_sequences_Linkers_ALL,'Q');

[LiaS, LocbS] = ismember(DNA_sequences_Linkers_ALL,'S');

[LiaP, LocbP] = ismember(DNA_sequences_Linkers_ALL,'P');

[LiaE, LocbE] = ismember(DNA_sequences_Linkers_ALL,'E');

[LiaK, LocbK] = ismember(DNA_sequences_Linkers_ALL,'K');

amount_dis_order_promoting = [sum([LocbI; LocbL; LocbV; LocbW; LocbY; LocbF; LocbC; LocbN]); sum([LocbR; LocbG; LocbQ; LocbS; LocbP; LocbE; LocbK])];

%% Plotting the data 

% positive aas

figure
H1=histfit(amount_positive,21,'kernel');
%set(H1(2),'Color',rand(1,3)); 
%set(H1(1),'FaceColor',rand(1,3)); 
set(H1(1),'FaceAlpha',.25);

% negative aas

figure
histfit(amount_negative, 21,'kernel')


%positive vs negative

figure
HP = histfit(amount_positive, 21,'kernel');
hold on 
HN =histfit(amount_negative, 21, 'kernel');
xlabel('Amount of Charged Amino Acids per Linker')
legend('Positively charged amino acids', '', 'Negatively charged amino acids', '')
ylim([0 10])
set(HP(1),'FaceAlpha',.4);
set(HN(2),'Color','red');
set(HN(1),'FaceAlpha',.4);
set(HP(2),'Color','blue');
saveas(gcf, 'DNA_sequence_analysis_linkers_charged_aminoacids.pdf','pdf')

%% length
figure
histfit(length_linkers,21,'kernel')
xlabel('Number of amino acids in linker')
saveas(gcf, 'DNA_sequence_analysis_linkers_length.pdf','pdf')

%% Order promoting/disorder promoting aas
figure
O = histfit(amount_dis_order_promoting(1,:),21,'kernel');
hold on
D = histfit(amount_dis_order_promoting(2,:),21,'kernel');
xlabel('Amount of dis-/order promoting amino acids per linker')
legend('Order promoting amino acids', '', 'Disorder promoting amino acids', '')
ylim([0 7])
xlim([-2 40])
set(O(1),'FaceAlpha',.4);
set(O(2),'Color','blue');
set(D(1),'FaceAlpha',.4);
set(D(2),'Color','red');
saveas(gcf, 'DNA_sequence_analysis_linkers_dis_order_aa.pdf','pdf')

% Order-Disorder
ratio_disorder = amount_dis_order_promoting(2,:).\amount_dis_order_promoting(1,:);
figure
R = histogram(ratio_disorder, 'BinWidth', 0.05);

%% order-disorder ratio vs length
figure
scatter(length_linkers, ratio_disorder,'*')
hold on
l = lsline;
l.LineWidth = 2;
r = corrcoef(ratio_disorder(1,2:end), length_linkers(1,2:end));
disp(r(1,2));
str=['r= ',num2str(r(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box on
grid on
xlabel('length of linkers')
ylabel('ratio of disorder vs order promoting amino acids')
saveas(gcf, 'DNA_sequence_analysis_linkers_corr_length_order.pdf','pdf')
%% Structure breaking
structure_breaking= sum([LocbG;LocbP]);

figure
scatter(length_linkers, structure_breaking,'*')
hold on
l = lsline;
l.LineWidth = 2;
r = corrcoef(structure_breaking(1,2:end), length_linkers(1,2:end));
disp(r(1,2));
str=['r= ',num2str(r(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box on
grid on
xlabel('length of linkers')
ylabel('structure breaking amino acids')
saveas(gcf, 'DNA_sequence_analysis_linkers_corr_length_strbreak.pdf','pdf')
%% positive charges vs length
figure
scatter(length_linkers, amount_positive,'*')
hold on
l = lsline;
l.LineWidth = 2;
r = corrcoef(length_linkers(1,2:end), amount_positive(1,2:end));
disp(r(1,2));
str=['r= ',num2str(r(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box on
grid on
xlabel('length of linkers')
ylabel('positively charged amino acids');
saveas(gcf, 'DNA_sequence_analysis_linkers_corr_length_positive.pdf','pdf')

%% negative charges vs length
figure
scatter(length_linkers, amount_negative,'*')
hold on
l = lsline;
l.LineWidth = 2;
r = corrcoef(length_linkers(1,2:end).', amount_negative(1,2:end));
disp(r(1,2));
str=['r= ',num2str(r(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box on
grid on
xlabel('length of linkers')
ylabel('negatively charged amino acids')
saveas(gcf, 'DNA_sequence_analysis_linkers_corr_length_negative.pdf','pdf')


%% negative charges vs length
figure
scatter(amount_positive, amount_negative,'*')
hold on
l = lsline;
l.LineWidth = 2;
r = corrcoef(amount_positive(1,2:end).', amount_negative(1,2:end));
disp(r(1,2));
str=['r= ',num2str(r(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
box on
grid on
xlabel('positively charged amino acids')
ylabel('negatively charged amino acids')
%saveas(gcf, 'DNA_sequence_analysis_linkers_corr_length_negative.pdf','pdf')



%% save matlab data to text file (for plotting with gnuplot)- Negative charges vs. length
fileID = fopen('DNA_sequence_ana_linkers_corr_length_negative.txt','w');
fprintf(fileID,'%6s %12s\n','length','negative');
fprintf(fileID,'%6.2f %12.8f\n',[length_linkers; amount_negative]);
fclose(fileID);

% Positive charges vs.length
fileID = fopen('DNA_sequence_ana_linkers_corr_length_positive.txt','w');
fprintf(fileID,'%6s %12s\n','length','positive');
fprintf(fileID,'%6.2f %12.8f\n',[length_linkers; amount_positive]);
fclose(fileID);

% Disorder promoting vs. length
fileID = fopen('DNA_sequence_ana_linkers_corr_length_order.txt','w');
fprintf(fileID,'%6s %12s\n','length','order');
fprintf(fileID,'%6.2f %12.8f\n',[length_linkers; ratio_disorder]);
fclose(fileID);

% Structure breaking vs. length
fileID = fopen('DNA_sequence_ana_linkers_corr_length_strbreak.txt','w');
fprintf(fileID,'%6s %12s\n','length','structure_breaking');
fprintf(fileID,'%6.2f %12.8f\n',[length_linkers; structure_breaking]);
fclose(fileID);

%% plot the average position of positively charged amino acids in different domains
%position_pos =sum(LiaR,2)+sum(LiaH,2)+sum(LiaK,2);
pos_2_pos = LiaR+LiaH+LiaK;
position_pos = pos_2_pos(:,2:end);


sub_domain_length = length_linkers(1,2:end)/5;
rounded = round(sub_domain_length);
region_1 = zeros(1,22);

for i = 1:22
    region_1(:,i) = sum(position_pos(1:rounded(1,i),i));
end

region_2 = zeros(1,22);

for i = 1:22
    region_2(:,i) = sum(position_pos(rounded(1,i):2*rounded(1,i),i));
end

region_3 = zeros(1,22);

for i = 1:22
    region_3(:,i) = sum(position_pos(2*rounded(1,i):3*rounded(1,i),i));
end

region_4 = zeros(1,22);

for i = 1:22
    region_4(:,i) = sum(position_pos(3*rounded(1,i):4*rounded(1,i),i));
end


region_5 = zeros(1,22);

for i = 1:22
    region_5(:,i) = sum(position_pos(4*rounded(1,i):length_linkers(1,i+1),i));
end


% averages 
region_1_mean = mean(region_1);
region_2_mean = mean(region_2);
region_3_mean = mean(region_3);
region_4_mean = mean(region_4);
region_5_mean = mean(region_5);
    
means = [region_1_mean region_2_mean region_3_mean region_4_mean region_5_mean];

x = linspace(1,5,5);

figure
plot(x,means)
%% negatively charged amino acids positioning
pos_2_neg = LiaD+LiaE;
position_neg = pos_2_neg(:,2:end);

region_1_neg = zeros(1,22);

for i = 1:22
    region_1_neg(:,i) = sum(position_neg(1:rounded(1,i),i));
end

region_2_neg = zeros(1,22);

for i = 1:22
    region_2_neg(:,i) = sum(position_neg(rounded(1,i):2*rounded(1,i),i));
end

region_3_neg = zeros(1,22);

for i = 1:22
    region_3_neg(:,i) = sum(position_neg(2*rounded(1,i):3*rounded(1,i),i));
end

region_4_neg = zeros(1,22);

for i = 1:22
    region_4_neg(:,i) = sum(position_neg(3*rounded(1,i):4*rounded(1,i),i));
end


region_5_neg = zeros(1,22);

for i = 1:22
    region_5_neg(:,i) = sum(position_neg(4*rounded(1,i):length_linkers(1,i+1),i));
end


% averages 
region_1_mean_neg = mean(region_1_neg);
region_2_mean_neg = mean(region_2_neg);
region_3_mean_neg = mean(region_3_neg);
region_4_mean_neg = mean(region_4_neg);
region_5_mean_neg = mean(region_5_neg);
    
means_neg = [region_1_mean_neg region_2_mean_neg region_3_mean_neg region_4_mean_neg region_5_mean_neg];

figure 
plot(x, means,'--*')
hold on 
plot(x, means_neg,'--o')


% save for gnuplot
fileID = fopen('DNA_sequence_linker_position_pos.txt','w');
fprintf(fileID,'%6s %12s\n','x','amount');
fprintf(fileID,'%6.2f %12.8f\n',[x; means]);
fclose(fileID);

fileID = fopen('DNA_sequence_linker_position_neg.txt','w');
fprintf(fileID,'%6s %12s\n','x','amount');
fprintf(fileID,'%6.2f %12.8f\n',[x; means_neg]);
fclose(fileID);

%% taking out the shortest linkers 


position_pos_short = pos_2_pos(:,[2 9:end]);


sub_domain_length_short = length_linkers(1,2:end)/8;
rounded_short = round(sub_domain_length_short);
region_1_short = zeros(1,16);

for i = 1:16
    region_1_short(:,i) = sum(position_pos_short(1:rounded_short(1,i),i));
end

region_2_short = zeros(1,16);

for i = 1:16
    region_2_short(:,i) = sum(position_pos_short(rounded_short(1,i):2*rounded_short(1,i),i));
end

region_3_short = zeros(1,16);

for i = 1:16
    region_3_short(:,i) = sum(position_pos_short(2*rounded_short(1,i):3*rounded_short(1,i),i));
end

region_4_short = zeros(1,16);

for i = 1:16
    region_4_short(:,i) = sum(position_pos_short(3*rounded_short(1,i):4*rounded_short(1,i),i));
end


region_5_short = zeros(1,16);

for i = 1:16
    region_5_short(:,i) = sum(position_pos_short(4*rounded_short(1,i):5*rounded_short(1,i),i));
end

% if want for more domains:
region_6_short = zeros(1,16);

for i = 1:16
    region_6_short(:,i) = sum(position_pos_short(5*rounded_short(1,i):6*rounded_short(1,i),i));
end


region_7_short = zeros(1,16);

for i = 1:16
    region_7_short(:,i) = sum(position_pos_short(6*rounded_short(1,i):7*rounded_short(1,i),i));
end


region_8_short = zeros(1,16);

for i = 1:16
    region_8_short(:,i) = sum(position_pos_short(7*rounded_short(1,i):length_linkers(1,i+1),i));
end

% averages 
region_1_mean_short = mean(region_1_short);
region_2_mean_short = mean(region_2_short);
region_3_mean_short = mean(region_3_short);
region_4_mean_short = mean(region_4_short);
region_5_mean_short = mean(region_5_short);
region_6_mean_short = mean(region_6_short);
region_7_mean_short = mean(region_7_short);
region_8_mean_short = mean(region_8_short);
    
means_short = [region_1_mean_short region_2_mean_short region_3_mean_short region_4_mean_short region_5_mean_short region_6_mean_short region_7_mean_short region_8_mean_short];

x = linspace(1,8,8);

figure
plot(x,means_short)



position_neg_short = pos_2_neg(:,[2 9:end]);


sub_domain_length_short = length_linkers(1,2:end)/8;
rounded_short = round(sub_domain_length_short);
region_1_short_neg = zeros(1,16);

for i = 1:16
    region_1_short_neg(:,i) = sum(position_neg_short(1:rounded_short(1,i),i));
end

region_2_short_neg = zeros(1,16);

for i = 1:16
    region_2_short_neg(:,i) = sum(position_neg_short(rounded_short(1,i):2*rounded_short(1,i),i));
end

region_3_short_neg = zeros(1,16);

for i = 1:16
    region_3_short_neg(:,i) = sum(position_neg_short(2*rounded_short(1,i):3*rounded_short(1,i),i));
end

region_4_short_neg = zeros(1,16);

for i = 1:16
    region_4_short_neg(:,i) = sum(position_neg_short(3*rounded_short(1,i):4*rounded_short(1,i),i));
end


region_5_short_neg = zeros(1,16);

for i = 1:16
    region_5_short_neg(:,i) = sum(position_neg_short(4*rounded_short(1,i):5*rounded_short(1,i),i));
end

% if want for more domains:
region_6_short_neg = zeros(1,16);

for i = 1:16
    region_6_short_neg(:,i) = sum(position_neg_short(5*rounded_short(1,i):6*rounded_short(1,i),i));
end


region_7_short_neg = zeros(1,16);

for i = 1:16
    region_7_short_neg(:,i) = sum(position_neg_short(6*rounded_short(1,i):7*rounded_short(1,i),i));
end


region_8_short_neg = zeros(1,16);

for i = 1:16
    region_8_short_neg(:,i) = sum(position_neg_short(7*rounded_short(1,i):length_linkers(1,i+1),i));
end

% averages 
region_1_mean_short_neg = mean(region_1_short_neg);
region_2_mean_short_neg = mean(region_2_short_neg);
region_3_mean_short_neg = mean(region_3_short_neg);
region_4_mean_short_neg = mean(region_4_short_neg);
region_5_mean_short_neg = mean(region_5_short_neg);
region_6_mean_short_neg = mean(region_6_short_neg);
region_7_mean_short_neg = mean(region_7_short_neg);
region_8_mean_short_neg = mean(region_8_short_neg);
    
means_short_neg = [region_1_mean_short_neg region_2_mean_short_neg region_3_mean_short_neg region_4_mean_short_neg region_5_mean_short_neg region_6_mean_short_neg region_7_mean_short_neg region_8_mean_short_neg];

x = linspace(1,8,8);

figure
plot(x,means_short)
hold on 
plot(x, means_short_neg)

%probability density function (here for positively charged amino acids with
%errorbars)
regions_all_short = [region_1_short' region_2_short' region_3_short' region_4_short' region_5_short' region_6_short' region_7_short' region_8_short'];
y = normpdf(regions_all_short, means_short);
y_mean = mean(y);
err_y_pos_short =[max(y(:,1))-y_mean(1,1) max(y(:,2))-y_mean(1,2) max(y(:,3))-y_mean(1,3) max(y(:,4))-y_mean(1,4) max(y(:,5))-y_mean(1,5) max(y(:,6))-y_mean(1,6) max(y(:,7))-y_mean(1,7) max(y(:,8))-y_mean(1,8)]; 

err_y_neg_short =[y_mean(1,1)-min(y(:,1)) y_mean(1,2)-min(y(:,2)) y_mean(1,3)-min(y(:,3)) y_mean(1,4)-min(y(:,4)) y_mean(1,5)-min(y(:,5)) y_mean(1,6)-min(y(:,6)) y_mean(1,7)-min(y(:,7)) y_mean(1,8)-min(y(:,8))];


figure 
errorbar(x,y_mean,err_y_neg_short, err_y_pos_short, '--o')


% save for gnuplot
fileID = fopen('DNA_sequence_linker_position_pos_short.txt','w');
fprintf(fileID,'%6s %12s\n','x','amount');
fprintf(fileID,'%6.2f %12.8f\n',[x; means_short]);
fclose(fileID);

fileID = fopen('DNA_sequence_linker_position_neg_short.txt','w');
fprintf(fileID,'%6s %12s\n','x','amount');
fprintf(fileID,'%6.2f %12.8f\n',[x; means_short_neg]);
fclose(fileID);
%%

%probability density function (here for positively charged amino acids with
%errorbars)
regions_all = [region_1' region_2' region_3' region_4' region_5'];
y = normpdf(regions_all, means);
y_mean = mean(y);
err_y_pos =[max(y(:,1))-y_mean(1,1) max(y(:,2))-y_mean(1,2) max(y(:,3))-y_mean(1,3) max(y(:,4))-y_mean(1,4) max(y(:,5))-y_mean(1,5)];

err_y_neg =[y_mean(1,1)-min(y(:,1)) y_mean(1,2)-min(y(:,2)) y_mean(1,3)-min(y(:,3)) y_mean(1,4)-min(y(:,4)) y_mean(1,5)-min(y(:,5))];

x_2 = [1 2 3 4 5];
figure 
errorbar(x_2,y_mean,err_y_neg, err_y_pos, '--o')


%probability density function (here for negatively charged amino acids with
%errorbars)
regions_all_neg = [region_1_neg' region_2_neg' region_3_neg' region_4_neg' region_5_neg'];
y_neg = normpdf(regions_all_neg, means_neg);
y_mean_neg = mean(y_neg);
err_y_pos_neg =[max(y_neg(:,1))-y_mean_neg(1,1) max(y_neg(:,2))-y_mean_neg(1,2) max(y_neg(:,3))-y_mean_neg(1,3) max(y_neg(:,4))-y_mean_neg(1,4) max(y_neg(:,5))-y_mean_neg(1,5)];

err_y_neg_neg =[y_mean_neg(1,1)-min(y_neg(:,1)) y_mean_neg(1,2)-min(y_neg(:,2)) y_mean_neg(1,3)-min(y_neg(:,3)) y_mean_neg(1,4)-min(y_neg(:,4)) y_mean_neg(1,5)-min(y_neg(:,5))];

x_2 = [1 2 3 4 5];
figure 
errorbar(x_2,y_mean_neg,err_y_neg_neg, err_y_pos_neg, '--o')
hold on 
errorbar(x_2, y_mean, err_y_neg, err_y_pos , '--*')




fileID = fopen('DNA_sequence_linker_prob_dens.txt','w');
fprintf(fileID,'%s %s %s %s %s %s %s\n ','x','probability_positive','err_pos_pos', 'err_pos_neg', 'probability_negative', 'err_neg_pos', 'err_neg_neg');
fprintf(fileID,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %12.8f\n',[x_2; y_mean; err_y_pos; err_y_neg; y_mean_neg; err_y_pos_neg; err_y_neg_neg]);
fclose(fileID);


%% Statistical analysis based on Fleming et al. 

% Counts of length of protein chain separating successive appearances of
% the same amino acid for each protein. Start with trial example. 

trial = [2 1 1 2 1 2 1 1 1 1 1 1 2 2]; % trial vector where the amino acid in question is denoted by "2" and all other amino acids are "1"

trial_incices=find(trial==2); % find the indices of amino acid in question

trial_distance = trial_incices(2:end) - trial_incices(1:end-1)-1; % calculate how many other amino acids are in between. 


% can't find better way: make separate matrix for each protein (since they all have different lengths etc. Can't add too many zeros as it would influence the distribution). Distance from 1st amino acid to 1st wanted aa is also not taken into account. Start with
% Arginine

indices_R_P1 = [find(LiaR(:,2)==1); length_linkers(:,2)];
indices_R_P2 = [find(LiaR(:,3)==1); length_linkers(:,3)];
indices_R_P3 = [find(LiaR(:,4)==1); length_linkers(:,4)];
indices_R_P4 = [find(LiaR(:,5)==1); length_linkers(:,5)];
indices_R_P5 = [find(LiaR(:,6)==1); length_linkers(:,6)];
indices_R_P6 = [find(LiaR(:,7)==1); length_linkers(:,7)];
indices_R_P7 = [find(LiaR(:,8)==1); length_linkers(:,8)];
indices_R_P8 = [find(LiaR(:,9)==1); length_linkers(:,9)];
indices_R_P9 = [find(LiaR(:,10)==1); length_linkers(:,10)];
indices_R_P10 = [find(LiaR(:,11)==1); length_linkers(:,11)];
indices_R_P11 = [find(LiaR(:,12)==1); length_linkers(:,12)];
indices_R_P12 = [find(LiaR(:,13)==1); length_linkers(:,13)];
indices_R_P13 = [find(LiaR(:,14)==1); length_linkers(:,14)];
indices_R_P14 = [find(LiaR(:,15)==1); length_linkers(:,15)];
indices_R_P15 = [find(LiaR(:,16)==1); length_linkers(:,16)];
indices_R_P16 = [find(LiaR(:,17)==1); length_linkers(:,17)];
indices_R_P17 = [find(LiaR(:,18)==1); length_linkers(:,18)];
indices_R_P18 = [find(LiaR(:,19)==1); length_linkers(:,19)];
indices_R_P19 = [find(LiaR(:,20)==1); length_linkers(:,20)];
indices_R_P20 = [find(LiaR(:,21)==1); length_linkers(:,21)];
indices_R_P21 = [find(LiaR(:,22)==1); length_linkers(:,22)];
indices_R_P22 = [find(LiaR(:,23)==1); length_linkers(:,23)];


distance_R_P1 =[indices_R_P1(1,1)-1; indices_R_P1(2:end)-indices_R_P1(1:end-1)-1];
distance_R_P2 = [indices_R_P2(2:end) -indices_R_P2(1:end-1)-1];
distance_R_P3 = indices_R_P3(2:end) -indices_R_P3(1:end-1)-1;
distance_R_P4 = indices_R_P4(2:end) -indices_R_P4(1:end-1)-1;
distance_R_P5 = indices_R_P5(2:end) -indices_R_P5(1:end-1)-1;
distance_R_P6 = indices_R_P6(2:end) -indices_R_P6(1:end-1)-1;
distance_R_P7 = [indices_R_P7(1,1)-1; indices_R_P7(2:end)-indices_R_P7(1:end-1)-1];
distance_R_P8 = [indices_R_P1(1,1)-1; indices_R_P8(2:end)-indices_R_P8(1:end-1)-1];
distance_R_P9 = indices_R_P9(2:end) -indices_R_P9(1:end-1)-1;
distance_R_P10 = [indices_R_P10(1,1)-1;indices_R_P10(2:end)-indices_R_P10(1:end-1)-1];
distance_R_P11 = indices_R_P11(2:end) -indices_R_P11(1:end-1)-1;
distance_R_P12 = indices_R_P12(2:end) -indices_R_P12(1:end-1)-1;
distance_R_P13 = indices_R_P13(2:end) -indices_R_P13(1:end-1)-1;
distance_R_P14 = indices_R_P14(2:end) -indices_R_P14(1:end-1)-1;
distance_R_P15 = indices_R_P15(2:end) -indices_R_P15(1:end-1)-1;
distance_R_P16 = indices_R_P16(2:end) -indices_R_P16(1:end-1)-1;
distance_R_P17 = indices_R_P17(2:end) -indices_R_P17(1:end-1)-1;
distance_R_P18 = [indices_R_P18(1,1)-1;indices_R_P18(2:end)-indices_R_P18(1:end-1)-1];
distance_R_P19 = indices_R_P19(2:end) -indices_R_P19(1:end-1)-1;
distance_R_P20 = indices_R_P20(2:end) -indices_R_P20(1:end-1)-1;
distance_R_P21 = [indices_R_P21(1,1)-1;indices_R_P21(2:end)-indices_R_P21(1:end-1)-1];
distance_R_P22 = [indices_R_P22(1,1)-1;indices_R_P22(2:end)-indices_R_P22(1:end-1)-1];


[h_R, p_R] = chi2gof([distance_R_P1; distance_R_P10; distance_R_P18; distance_R_P8; distance_R_P7; distance_R_P21; distance_R_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.


R_position_distribution = [distance_R_P1; distance_R_P10; distance_R_P18; distance_R_P8; distance_R_P7; distance_R_P21; distance_R_P22];

figure
histogram(R_position_distribution)

%% Histidine

indices_H_P1 = [find(LiaH(:,2)==1); length_linkers(:,2)];
indices_H_P2 = [find(LiaH(:,3)==1); length_linkers(:,3)];
indices_H_P3 = [find(LiaH(:,4)==1); length_linkers(:,4)];
indices_H_P4 = [find(LiaH(:,5)==1); length_linkers(:,5)];
indices_H_P5 = [find(LiaH(:,6)==1); length_linkers(:,6)];
indices_H_P6 = [find(LiaH(:,7)==1); length_linkers(:,7)];
indices_H_P7 = [find(LiaH(:,8)==1); length_linkers(:,8)];
indices_H_P8 = [find(LiaH(:,9)==1); length_linkers(:,9)];
indices_H_P9 = [find(LiaH(:,10)==1); length_linkers(:,10)];
indices_H_P10 = [find(LiaH(:,11)==1); length_linkers(:,11)];
indices_H_P11 = [find(LiaH(:,12)==1); length_linkers(:,12)];
indices_H_P12 = [find(LiaH(:,13)==1); length_linkers(:,13)];
indices_H_P13 = [find(LiaH(:,14)==1); length_linkers(:,14)];
indices_H_P14 = [find(LiaH(:,15)==1); length_linkers(:,15)];
indices_H_P15 = [find(LiaH(:,16)==1); length_linkers(:,16)];
indices_H_P16 = [find(LiaH(:,17)==1); length_linkers(:,17)];
indices_H_P17 = [find(LiaH(:,18)==1); length_linkers(:,18)];
indices_H_P18 = [find(LiaH(:,19)==1); length_linkers(:,19)];
indices_H_P19 = [find(LiaH(:,20)==1); length_linkers(:,20)];
indices_H_P20 = [find(LiaH(:,21)==1); length_linkers(:,21)];
indices_H_P21 = [find(LiaH(:,22)==1); length_linkers(:,22)];
indices_H_P22 = [find(LiaH(:,23)==1); length_linkers(:,23)];


distance_H_P1 = [indices_H_P1(1,1)-1;indices_H_P1(2:end)-indices_H_P1(1:end-1)-1];
distance_H_P2 = indices_H_P2(2:end) -indices_H_P2(1:end-1)-1;
distance_H_P3 = indices_H_P3(2:end) -indices_H_P3(1:end-1)-1;
distance_H_P4 = indices_H_P4(2:end) -indices_H_P4(1:end-1)-1;
distance_H_P5 = indices_H_P5(2:end) -indices_H_P5(1:end-1)-1;
distance_H_P6 = indices_H_P6(2:end) -indices_H_P6(1:end-1)-1;
distance_H_P7 = indices_H_P7(2:end) -indices_H_P7(1:end-1)-1;
distance_H_P8 = indices_H_P8(2:end) -indices_H_P8(1:end-1)-1;
distance_H_P9 = indices_H_P9(2:end) -indices_H_P9(1:end-1)-1;
distance_H_P10 = indices_H_P10(2:end) -indices_H_P10(1:end-1)-1;
distance_H_P11 = [indices_H_P11(1,1)-1;indices_H_P11(2:end)-indices_H_P11(1:end-1)-1];
distance_H_P12 = indices_H_P12(2:end) -indices_H_P12(1:end-1)-1;
distance_H_P13 = indices_H_P13(2:end) -indices_H_P13(1:end-1)-1;
distance_H_P14 = indices_H_P14(2:end) -indices_H_P14(1:end-1)-1;
distance_H_P15 = indices_H_P15(2:end) -indices_H_P15(1:end-1)-1;
distance_H_P16 = [indices_H_P16(1,1)-1;indices_H_P16(2:end)-indices_H_P16(1:end-1)-1];
distance_H_P17 = indices_H_P17(2:end) -indices_H_P17(1:end-1)-1;
distance_H_P18 = [indices_H_P18(1,1)-1;indices_H_P18(2:end)-indices_H_P18(1:end-1)-1];
distance_H_P19 = indices_H_P19(2:end) -indices_H_P19(1:end-1)-1;
distance_H_P20 = indices_H_P20(2:end) -indices_H_P20(1:end-1)-1;
distance_H_P21 = [indices_H_P21(1,1)-1; indices_H_P21(2:end)-indices_H_P21(1:end-1)-1];
distance_H_P22 = [indices_H_P22(1,1)-1;indices_H_P22(2:end)-indices_H_P22(1:end-1)-1];


[h_H, p_H] = chi2gof([distance_H_P1; distance_H_P2; distance_H_P3; distance_H_P4; distance_H_P5; distance_H_P6; distance_H_P8; distance_H_P7;distance_H_P9; distance_H_P10; distance_H_P11; distance_H_P12; distance_H_P13; distance_H_P14; distance_H_P15; distance_H_P16; distance_H_P17; distance_H_P18; distance_H_P19; distance_H_P20; distance_H_P21; distance_H_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.


%% Lysine 

indices_K_P1 = [find(LiaK(:,2)==1); length_linkers(:,2)];
indices_K_P2 = [find(LiaK(:,3)==1); length_linkers(:,3)];
indices_K_P3 = [find(LiaK(:,4)==1); length_linkers(:,4)];
indices_K_P4 = [find(LiaK(:,5)==1); length_linkers(:,5)];
indices_K_P5 = [find(LiaK(:,6)==1); length_linkers(:,6)];
indices_K_P6 = [find(LiaK(:,7)==1); length_linkers(:,7)];
indices_K_P7 = [find(LiaK(:,8)==1); length_linkers(:,8)];
indices_K_P8 = [find(LiaK(:,9)==1); length_linkers(:,9)];
indices_K_P9 = [find(LiaK(:,10)==1); length_linkers(:,10)];
indices_K_P10 = [find(LiaK(:,11)==1); length_linkers(:,11)];
indices_K_P11 = [find(LiaK(:,12)==1); length_linkers(:,12)];
indices_K_P12 = [find(LiaK(:,13)==1); length_linkers(:,13)];
indices_K_P13 = [find(LiaK(:,14)==1); length_linkers(:,14)];
indices_K_P14 = [find(LiaK(:,15)==1); length_linkers(:,15)];
indices_K_P15 = [find(LiaK(:,16)==1); length_linkers(:,16)];
indices_K_P16 = [find(LiaK(:,17)==1); length_linkers(:,17)];
indices_K_P17 = [find(LiaK(:,18)==1); length_linkers(:,18)];
indices_K_P18 = [find(LiaK(:,19)==1); length_linkers(:,19)];
indices_K_P19 = [find(LiaK(:,20)==1); length_linkers(:,20)];
indices_K_P20 = [find(LiaK(:,21)==1); length_linkers(:,21)];
indices_K_P21 = [find(LiaK(:,22)==1); length_linkers(:,22)];
indices_K_P22 = [find(LiaK(:,23)==1); length_linkers(:,23)];


distance_K_P1 = [indices_K_P1(1,1)-1;indices_K_P1(2:end)-indices_K_P1(1:end-1)-1];
distance_K_P2 = [indices_K_P2(1,1)-1;indices_K_P2(2:end)-indices_K_P2(1:end-1)-1];
distance_K_P3 = [indices_K_P3(1,1)-1;indices_K_P3(2:end)-indices_K_P3(1:end-1)-1];
distance_K_P4 = [indices_K_P4(1,1)-1;indices_K_P4(2:end)-indices_K_P4(1:end-1)-1];
distance_K_P5 = [indices_K_P5(1,1)-1;indices_K_P5(2:end)-indices_K_P5(1:end-1)-1];
distance_K_P6 = [indices_K_P6(1,1)-1;indices_K_P6(2:end)-indices_K_P6(1:end-1)-1];
distance_K_P7 = [indices_K_P7(1,1)-1;indices_K_P7(2:end)-indices_K_P7(1:end-1)-1];
distance_K_P8 = [indices_K_P8(1,1)-1;indices_K_P8(2:end)-indices_K_P8(1:end-1)-1];
distance_K_P9 = [indices_K_P9(1,1)-1;indices_K_P9(2:end)-indices_K_P9(1:end-1)-1];
distance_K_P10 = indices_K_P10(2:end) -indices_K_P10(1:end-1)-1;
distance_K_P11 = indices_K_P11(2:end) -indices_K_P11(1:end-1)-1;
distance_K_P12 = indices_K_P12(2:end) -indices_K_P12(1:end-1)-1;
distance_K_P13 = indices_K_P13(2:end) -indices_K_P13(1:end-1)-1;
distance_K_P14 = indices_K_P14(2:end) -indices_K_P14(1:end-1)-1;
distance_K_P15 = indices_K_P15(2:end) -indices_K_P15(1:end-1)-1;
distance_K_P16 = indices_K_P16(2:end) -indices_K_P16(1:end-1)-1;
distance_K_P17 = indices_K_P17(2:end) -indices_K_P17(1:end-1)-1;
distance_K_P18 = indices_K_P18(2:end) -indices_K_P18(1:end-1)-1;
distance_K_P19 = indices_K_P19(2:end) -indices_K_P19(1:end-1)-1;
distance_K_P20 = indices_K_P20(2:end) -indices_K_P20(1:end-1)-1;
distance_K_P21 = [indices_K_P21(1,1)-1;indices_K_P21(2:end)-indices_K_P21(1:end-1)-1];
distance_K_P22 = [indices_K_P22(1,1)-1;indices_K_P22(2:end)-indices_K_P22(1:end-1)-1];


[h_K, p_K] = chi2gof([distance_K_P1; distance_K_P2; distance_K_P3; distance_K_P4; distance_K_P5; distance_K_P6; distance_K_P8; distance_K_P7;distance_K_P9; distance_K_P10; distance_K_P11; distance_K_P12; distance_K_P13; distance_K_P14; distance_K_P15; distance_K_P16; distance_K_P17; distance_K_P18; distance_K_P19; distance_K_P20; distance_K_P21; distance_K_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.


%% Apartic Acid D
indices_D_P1 = [find(LiaD(:,2)==1); length_linkers(:,2)];
indices_D_P2 = [find(LiaD(:,3)==1); length_linkers(:,3)];
indices_D_P3 = [find(LiaD(:,4)==1); length_linkers(:,4)];
indices_D_P4 = [find(LiaD(:,5)==1); length_linkers(:,5)];
indices_D_P5 = [find(LiaD(:,6)==1); length_linkers(:,6)];
indices_D_P6 = [find(LiaD(:,7)==1); length_linkers(:,7)];
indices_D_P7 = [find(LiaD(:,8)==1); length_linkers(:,8)];
indices_D_P8 = [find(LiaD(:,9)==1); length_linkers(:,9)];
indices_D_P9 = [find(LiaD(:,10)==1); length_linkers(:,10)];
indices_D_P10 = [find(LiaD(:,11)==1); length_linkers(:,11)];
indices_D_P11 = [find(LiaD(:,12)==1); length_linkers(:,12)];
indices_D_P12 = [find(LiaD(:,13)==1); length_linkers(:,13)];
indices_D_P13 = [find(LiaD(:,14)==1); length_linkers(:,14)];
indices_D_P14 = [find(LiaD(:,15)==1); length_linkers(:,15)];
indices_D_P15 = [find(LiaD(:,16)==1); length_linkers(:,16)];
indices_D_P16 = [find(LiaD(:,17)==1); length_linkers(:,17)];
indices_D_P17 = [find(LiaD(:,18)==1); length_linkers(:,18)];
indices_D_P18 = [find(LiaD(:,19)==1); length_linkers(:,19)];
indices_D_P19 = [find(LiaD(:,20)==1); length_linkers(:,20)];
indices_D_P20 = [find(LiaD(:,21)==1); length_linkers(:,21)];
indices_D_P21 = [find(LiaD(:,22)==1); length_linkers(:,22)];
indices_D_P22 = [find(LiaD(:,23)==1); length_linkers(:,23)];


distance_D_P1 = [indices_D_P1(1,1)-1;indices_D_P1(2:end)-indices_D_P1(1:end-1)-1];
distance_D_P2 = indices_D_P2(2:end) -indices_D_P2(1:end-1)-1;
distance_D_P3 = indices_D_P3(2:end) -indices_D_P3(1:end-1)-1;
distance_D_P4 = indices_D_P4(2:end) -indices_D_P4(1:end-1)-1;
distance_D_P5 = indices_D_P5(2:end) -indices_D_P5(1:end-1)-1;
distance_D_P6 = [indices_D_P6(1,1)-1;indices_D_P6(2:end)-indices_D_P6(1:end-1)-1];
distance_D_P7 = indices_D_P7(2:end) -indices_D_P7(1:end-1)-1;
distance_D_P8 = indices_D_P8(2:end) -indices_D_P8(1:end-1)-1;
distance_D_P9 = indices_D_P9(2:end) -indices_D_P9(1:end-1)-1;
distance_D_P10 = [indices_D_P10(1,1)-1;indices_D_P10(2:end)-indices_D_P10(1:end-1)-1];
distance_D_P11 = indices_D_P11(2:end) -indices_D_P11(1:end-1)-1;
distance_D_P12 = [indices_D_P12(1,1)-1;indices_D_P12(2:end)-indices_D_P12(1:end-1)-1];
distance_D_P13 = [indices_D_P13(1,1)-1;indices_D_P13(2:end)-indices_D_P13(1:end-1)-1];
distance_D_P14 = [indices_D_P14(1,1)-1;indices_D_P14(2:end)-indices_D_P14(1:end-1)-1];
distance_D_P15 = [indices_D_P15(1,1)-1;indices_D_P15(2:end)-indices_D_P15(1:end-1)-1];
distance_D_P16 = [indices_D_P16(1,1)-1;indices_D_P16(2:end)-indices_D_P16(1:end-1)-1];
distance_D_P17 = [indices_D_P17(1,1)-1;indices_D_P17(2:end)-indices_D_P17(1:end-1)-1];
distance_D_P18 = [indices_D_P18(1,1)-1;indices_D_P18(2:end)-indices_D_P18(1:end-1)-1];
distance_D_P19 = [indices_D_P19(1,1)-1;indices_D_P19(2:end)-indices_D_P19(1:end-1)-1];
distance_D_P20 = [indices_D_P20(1,1)-1;indices_D_P20(2:end)-indices_D_P20(1:end-1)-1];
distance_D_P21 = indices_D_P21(2:end) -indices_D_P21(1:end-1)-1;
distance_D_P22 = [indices_D_P22(1,1)-1;indices_D_P22(2:end)-indices_D_P22(1:end-1)-1];


[h_D, p_D] = chi2gof([distance_D_P1; distance_D_P2; distance_D_P3; distance_D_P4; distance_D_P5; distance_D_P6; distance_D_P8; distance_D_P7;distance_D_P9; distance_D_P10; distance_D_P11; distance_D_P12; distance_D_P13; distance_D_P14; distance_D_P15; distance_D_P16; distance_D_P17; distance_D_P18; distance_D_P19; distance_D_P20; distance_D_P21; distance_D_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.

%% Glutamic acid E

indices_E_P1 = [find(LiaE(:,2)==1); length_linkers(:,2)];
indices_E_P2 = [find(LiaE(:,3)==1); length_linkers(:,3)];
indices_E_P3 = [find(LiaE(:,4)==1); length_linkers(:,4)];
indices_E_P4 = [find(LiaE(:,5)==1); length_linkers(:,5)];
indices_E_P5 = [find(LiaE(:,6)==1); length_linkers(:,6)];
indices_E_P6 = [find(LiaE(:,7)==1); length_linkers(:,7)];
indices_E_P7 = [find(LiaE(:,8)==1); length_linkers(:,8)];
indices_E_P8 = [find(LiaE(:,9)==1); length_linkers(:,9)];
indices_E_P9 = [find(LiaE(:,10)==1); length_linkers(:,10)];
indices_E_P10 = [find(LiaE(:,11)==1); length_linkers(:,11)];
indices_E_P11 = [find(LiaE(:,12)==1); length_linkers(:,12)];
indices_E_P12 = [find(LiaE(:,13)==1); length_linkers(:,13)];
indices_E_P13 = [find(LiaE(:,14)==1); length_linkers(:,14)];
indices_E_P14 = [find(LiaE(:,15)==1); length_linkers(:,15)];
indices_E_P15 = [find(LiaE(:,16)==1); length_linkers(:,16)];
indices_E_P16 = [find(LiaE(:,17)==1); length_linkers(:,17)];
indices_E_P17 = [find(LiaE(:,18)==1); length_linkers(:,18)];
indices_E_P18 = [find(LiaE(:,19)==1); length_linkers(:,19)];
indices_E_P19 = [find(LiaE(:,20)==1); length_linkers(:,20)];
indices_E_P20 = [find(LiaE(:,21)==1); length_linkers(:,21)];
indices_E_P21 = [find(LiaE(:,22)==1); length_linkers(:,22)];
indices_E_P22 = [find(LiaE(:,23)==1); length_linkers(:,23)];


distance_E_P1 = [indices_E_P1(1,1)-1;indices_E_P1(2:end)-indices_E_P1(1:end-1)-1];
distance_E_P2 = [indices_E_P2(1,1)-1;indices_E_P2(2:end)-indices_E_P2(1:end-1)-1];
distance_E_P3 = [indices_E_P3(1,1)-1;indices_E_P3(2:end)-indices_E_P3(1:end-1)-1];
distance_E_P4 = indices_E_P4(2:end) -indices_E_P4(1:end-1)-1;
distance_E_P5 = [indices_E_P5(1,1)-1;indices_E_P5(2:end)-indices_E_P5(1:end-1)-1];
distance_E_P6 = indices_E_P6(2:end) -indices_E_P6(1:end-1)-1;
distance_E_P7 = indices_E_P7(2:end) -indices_E_P7(1:end-1)-1;
distance_E_P8 = indices_E_P8(2:end) -indices_E_P8(1:end-1)-1;
distance_E_P9 = [indices_E_P9(1,1)-1;indices_E_P9(2:end)-indices_E_P9(1:end-1)-1];
distance_E_P10 = [indices_E_P10(1,1)-1;indices_E_P10(2:end)-indices_E_P10(1:end-1)-1];
distance_E_P11 = [indices_E_P11(1,1)-1;indices_E_P11(2:end)-indices_E_P11(1:end-1)-1];
distance_E_P12 = [indices_E_P12(1,1)-1;indices_E_P12(2:end)-indices_E_P12(1:end-1)-1];
distance_E_P13 = [indices_E_P13(1,1)-1;indices_E_P13(2:end)-indices_E_P13(1:end-1)-1];
distance_E_P14 = [indices_E_P14(1,1)-1;indices_E_P14(2:end)-indices_E_P14(1:end-1)-1];
distance_E_P15 = [indices_E_P15(1,1)-1;indices_E_P15(2:end)-indices_E_P15(1:end-1)-1];
distance_E_P16 = [indices_E_P16(1,1)-1;indices_E_P16(2:end)-indices_E_P16(1:end-1)-1];
distance_E_P17 = indices_E_P17(2:end) -indices_E_P17(1:end-1)-1;
distance_E_P18 = [indices_E_P18(1,1)-1;indices_E_P18(2:end)-indices_E_P18(1:end-1)-1];
distance_E_P19 =[indices_E_P19(1,1)-1; indices_E_P19(2:end)-indices_E_P19(1:end-1)-1];
distance_E_P20 = [indices_E_P20(1,1)-1;indices_E_P20(2:end)-indices_E_P20(1:end-1)-1];
distance_E_P21 = [indices_E_P21(1,1)-1;indices_E_P21(2:end)-indices_E_P21(1:end-1)-1];
distance_E_P22 = [indices_E_P22(1,1)-1;indices_E_P22(2:end)-indices_E_P22(1:end-1)-1];


[h_E,p_E] = chi2gof([distance_E_P1; distance_E_P2; distance_E_P3; distance_E_P4; distance_E_P5; distance_E_P6; distance_E_P8; distance_E_P7;distance_E_P9; distance_E_P10; distance_E_P11; distance_E_P12; distance_E_P13; distance_E_P14; distance_E_P15; distance_E_P16; distance_E_P17; distance_E_P18; distance_E_P19; distance_E_P20; distance_E_P21; distance_E_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.
% The returned value h = 1 indicates that chi2gof rejects the null hypothesis at the default 5% significance level.

%% Serine S 


indices_S_P1 = [find(LiaS(:,2)==1); length_linkers(:,2)];
indices_S_P2 = [find(LiaS(:,3)==1); length_linkers(:,3)];
indices_S_P3 = [find(LiaS(:,4)==1); length_linkers(:,4)];
indices_S_P4 = [find(LiaS(:,5)==1); length_linkers(:,5)];
indices_S_P5 = [find(LiaS(:,6)==1); length_linkers(:,6)];
indices_S_P6 = [find(LiaS(:,7)==1); length_linkers(:,7)];
indices_S_P7 = [find(LiaS(:,8)==1); length_linkers(:,8)];
indices_S_P8 = [find(LiaS(:,9)==1); length_linkers(:,9)];
indices_S_P9 = [find(LiaS(:,10)==1); length_linkers(:,10)];
indices_S_P10 = [find(LiaS(:,11)==1); length_linkers(:,11)];
indices_S_P11 = [find(LiaS(:,12)==1); length_linkers(:,12)];
indices_S_P12 = [find(LiaS(:,13)==1); length_linkers(:,13)];
indices_S_P13 = [find(LiaS(:,14)==1); length_linkers(:,14)];
indices_S_P14 = [find(LiaS(:,15)==1); length_linkers(:,15)];
indices_S_P15 = [find(LiaS(:,16)==1); length_linkers(:,16)];
indices_S_P16 = [find(LiaS(:,17)==1); length_linkers(:,17)];
indices_S_P17 = [find(LiaS(:,18)==1); length_linkers(:,18)];
indices_S_P18 = [find(LiaS(:,19)==1); length_linkers(:,19)];
indices_S_P19 = [find(LiaS(:,20)==1); length_linkers(:,20)];
indices_S_P20 = [find(LiaS(:,21)==1); length_linkers(:,21)];
indices_S_P21 = [find(LiaS(:,22)==1); length_linkers(:,22)];
indices_S_P22 = [find(LiaS(:,23)==1); length_linkers(:,23)];


distance_S_P1 = [indices_S_P1(1,1)-1;indices_S_P1(2:end)-indices_S_P1(1:end-1)-1];
distance_S_P2 = indices_S_P2(2:end) -indices_S_P2(1:end-1)-1;
distance_S_P3 = indices_S_P3(2:end) -indices_S_P3(1:end-1)-1;
distance_S_P4 = indices_S_P4(2:end) -indices_S_P4(1:end-1)-1;
distance_S_P5 = [indices_S_P5(1,1)-1;indices_S_P5(2:end)-indices_S_P5(1:end-1)-1];
distance_S_P6 = indices_S_P6(2:end) -indices_S_P6(1:end-1)-1;
distance_S_P7 = [indices_S_P7(1,1)-1;indices_S_P7(2:end)-indices_S_P7(1:end-1)-1];
distance_S_P8 = [indices_S_P8(1,1)-1;indices_S_P8(2:end)-indices_S_P8(1:end-1)-1];
distance_S_P9 = [indices_S_P9(1,1)-1;indices_S_P9(2:end)-indices_S_P9(1:end-1)-1];
distance_S_P10 = [indices_S_P10(1,1)-1;indices_S_P10(2:end)-indices_S_P10(1:end-1)-1];
distance_S_P11 = [indices_S_P11(1,1)-1;indices_S_P11(2:end)-indices_S_P11(1:end-1)-1];
distance_S_P12 = [indices_S_P12(1,1)-1;indices_S_P12(2:end)-indices_S_P12(1:end-1)-1];
distance_S_P13 = [indices_R_P13(1,1)-1;indices_S_P13(2:end)-indices_S_P13(1:end-1)-1];
distance_S_P14 = [indices_S_P14(1,1)-1;indices_S_P14(2:end)-indices_S_P14(1:end-1)-1];
distance_S_P15 = [indices_S_P15(1,1)-1;indices_S_P15(2:end)-indices_S_P15(1:end-1)-1];
distance_S_P16 = [indices_S_P16(1,1)-1;indices_S_P16(2:end)-indices_S_P16(1:end-1)-1];
distance_S_P17 = [indices_S_P17(1,1)-1;indices_S_P17(2:end)-indices_S_P17(1:end-1)-1];
distance_S_P18 = [indices_S_P18(1,1)-1;indices_S_P18(2:end)-indices_S_P18(1:end-1)-1];
distance_S_P19 = [indices_S_P19(1,1)-1;indices_S_P19(2:end)-indices_S_P19(1:end-1)-1];
distance_S_P20 = [indices_S_P20(1,1)-1;indices_S_P20(2:end)-indices_S_P20(1:end-1)-1];
distance_S_P21 = [indices_S_P21(1,1)-1;indices_S_P21(2:end)-indices_S_P21(1:end-1)-1];
distance_S_P22 = [indices_S_P22(1,1)-1;indices_S_P22(2:end)-indices_S_P22(1:end-1)-1];


[h_S, p_S] = chi2gof([distance_S_P1; distance_S_P2; distance_S_P3; distance_S_P4; distance_S_P5; distance_S_P6; distance_S_P8; distance_S_P7;distance_S_P9; distance_S_P10; distance_S_P11; distance_S_P12; distance_S_P13; distance_S_P14; distance_S_P15; distance_S_P16; distance_S_P17; distance_S_P18; distance_S_P19; distance_S_P20; distance_S_P21; distance_S_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.


%% Glutamine (as this will be the amino acid to be incorporated in the linkers as filler 


indices_Q_P1 = [find(LiaQ(:,2)==1); length_linkers(:,2)];
indices_Q_P2 = [find(LiaQ(:,3)==1); length_linkers(:,3)];
indices_Q_P3 = [find(LiaQ(:,4)==1); length_linkers(:,4)];
indices_Q_P4 = [find(LiaQ(:,5)==1); length_linkers(:,5)];
indices_Q_P5 = [find(LiaQ(:,6)==1); length_linkers(:,6)];
indices_Q_P6 = [find(LiaQ(:,7)==1); length_linkers(:,7)];
indices_Q_P7 = [find(LiaQ(:,8)==1); length_linkers(:,8)];
indices_Q_P8 = [find(LiaQ(:,9)==1); length_linkers(:,9)];
indices_Q_P9 = [find(LiaQ(:,10)==1); length_linkers(:,10)];
indices_Q_P10 = [find(LiaQ(:,11)==1); length_linkers(:,11)];
indices_Q_P11 = [find(LiaQ(:,12)==1); length_linkers(:,12)];
indices_Q_P12 = [find(LiaQ(:,13)==1); length_linkers(:,13)];
indices_Q_P13 = [find(LiaQ(:,14)==1); length_linkers(:,14)];
indices_Q_P14 = [find(LiaQ(:,15)==1); length_linkers(:,15)];
indices_Q_P15 = [find(LiaQ(:,16)==1); length_linkers(:,16)];
indices_Q_P16 = [find(LiaQ(:,17)==1); length_linkers(:,17)];
indices_Q_P17 = [find(LiaQ(:,18)==1); length_linkers(:,18)];
indices_Q_P18 = [find(LiaQ(:,19)==1); length_linkers(:,19)];
indices_Q_P19 = [find(LiaQ(:,20)==1); length_linkers(:,20)];
indices_Q_P20 = [find(LiaQ(:,21)==1); length_linkers(:,21)];
indices_Q_P21 = [find(LiaQ(:,22)==1); length_linkers(:,22)];
indices_Q_P22 = [find(LiaQ(:,23)==1); length_linkers(:,23)];


distance_Q_P1 = indices_Q_P1(2:end)-indices_Q_P1(1:end-1)-1;
distance_Q_P2 = indices_Q_P2(2:end) -indices_Q_P2(1:end-1)-1;
distance_Q_P3 = indices_Q_P3(2:end) -indices_Q_P3(1:end-1)-1;
distance_Q_P4 = [indices_Q_P4(1,1)-1;indices_Q_P4(2:end)-indices_Q_P4(1:end-1)-1];
distance_Q_P5 = indices_Q_P5(2:end)-indices_Q_P5(1:end-1)-1;
distance_Q_P6 = indices_Q_P6(2:end) -indices_Q_P6(1:end-1)-1;
distance_Q_P7 = indices_Q_P7(2:end)-indices_Q_P7(1:end-1)-1;
distance_Q_P8 = [indices_Q_P8(1,1)-1;indices_Q_P8(2:end)-indices_Q_P8(1:end-1)-1];
distance_Q_P9 = indices_Q_P9(2:end)-indices_Q_P9(1:end-1)-1;
distance_Q_P10 = [indices_Q_P10(1,1)-1;indices_Q_P10(2:end)-indices_Q_P10(1:end-1)-1];
distance_Q_P11 = indices_Q_P11(2:end)-indices_Q_P11(1:end-1)-1;
distance_Q_P12 = indices_Q_P12(2:end)-indices_Q_P12(1:end-1)-1;
distance_Q_P13 = [indices_Q_P13(1,1)-1;indices_Q_P13(2:end)-indices_Q_P13(1:end-1)-1];
distance_Q_P14 = indices_Q_P14(2:end)-indices_Q_P14(1:end-1)-1;
distance_Q_P15 = indices_Q_P15(2:end)-indices_Q_P15(1:end-1)-1;
distance_Q_P16 = [indices_Q_P16(1,1)-1;indices_Q_P16(2:end)-indices_Q_P16(1:end-1)-1];
distance_Q_P17 = [indices_Q_P17(1,1)-1;indices_Q_P17(2:end)-indices_Q_P17(1:end-1)-1];
distance_Q_P18 = [indices_Q_P18(1,1)-1;indices_Q_P18(2:end)-indices_Q_P18(1:end-1)-1];
distance_Q_P19 = indices_Q_P19(2:end)-indices_Q_P19(1:end-1)-1;
distance_Q_P20 = [indices_Q_P20(1,1)-1;indices_Q_P20(2:end)-indices_Q_P20(1:end-1)-1];
distance_Q_P21 = [indices_Q_P21(1,1)-1;indices_Q_P21(2:end)-indices_Q_P21(1:end-1)-1];
distance_Q_P22 = [indices_Q_P22(1,1)-1;indices_Q_P22(2:end)-indices_Q_P22(1:end-1)-1];


[h_Q, p_Q] = chi2gof([distance_Q_P1; distance_Q_P2; distance_Q_P3; distance_Q_P4; distance_Q_P5; distance_Q_P6; distance_Q_P8; distance_Q_P7;distance_Q_P9; distance_Q_P10; distance_Q_P11; distance_Q_P12; distance_Q_P13; distance_Q_P14; distance_Q_P15; distance_Q_P16; distance_Q_P17; distance_Q_P18; distance_Q_P19; distance_Q_P20; distance_Q_P21; distance_Q_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.



%% Glycine 


indices_G_P1 = [find(LiaG(:,2)==1); length_linkers(:,2)];
indices_G_P2 = [find(LiaG(:,3)==1); length_linkers(:,3)];
indices_G_P3 = [find(LiaG(:,4)==1); length_linkers(:,4)];
indices_G_P4 = [find(LiaG(:,5)==1); length_linkers(:,5)];
indices_G_P5 = [find(LiaG(:,6)==1); length_linkers(:,6)];
indices_G_P6 = [find(LiaG(:,7)==1); length_linkers(:,7)];
indices_G_P7 = [find(LiaG(:,8)==1); length_linkers(:,8)];
indices_G_P8 = [find(LiaG(:,9)==1); length_linkers(:,9)];
indices_G_P9 = [find(LiaG(:,10)==1); length_linkers(:,10)];
indices_G_P10 = [find(LiaG(:,11)==1); length_linkers(:,11)];
indices_G_P11 = [find(LiaG(:,12)==1); length_linkers(:,12)];
indices_G_P12 = [find(LiaG(:,13)==1); length_linkers(:,13)];
indices_G_P13 = [find(LiaG(:,14)==1); length_linkers(:,14)];
indices_G_P14 = [find(LiaG(:,15)==1); length_linkers(:,15)];
indices_G_P15 = [find(LiaG(:,16)==1); length_linkers(:,16)];
indices_G_P16 = [find(LiaG(:,17)==1); length_linkers(:,17)];
indices_G_P17 = [find(LiaG(:,18)==1); length_linkers(:,18)];
indices_G_P18 = [find(LiaG(:,19)==1); length_linkers(:,19)];
indices_G_P19 = [find(LiaG(:,20)==1); length_linkers(:,20)];
indices_G_P20 = [find(LiaG(:,21)==1); length_linkers(:,21)];
indices_G_P21 = [find(LiaG(:,22)==1); length_linkers(:,22)];
indices_G_P22 = [find(LiaG(:,23)==1); length_linkers(:,23)];


distance_G_P1 = [indices_G_P1(1,1)-1;indices_G_P1(2:end)-indices_G_P1(1:end-1)-1];
distance_G_P2 = [indices_G_P2(1,1)-1;indices_G_P2(2:end)-indices_G_P2(1:end-1)-1];
distance_G_P3 = [indices_G_P3(1,1)-1;indices_G_P3(2:end)-indices_G_P3(1:end-1)-1];
distance_G_P4 = [indices_G_P4(1,1)-1;indices_G_P4(2:end)-indices_G_P4(1:end-1)-1];
distance_G_P5 = indices_G_P5(2:end)-indices_G_P5(1:end-1)-1;
distance_G_P6 = [indices_G_P6(1,1)-1;indices_G_P6(2:end)-indices_G_P6(1:end-1)-1];
distance_G_P7 = indices_G_P7(2:end)-indices_G_P7(1:end-1)-1;
distance_G_P8 = [indices_G_P8(1,1)-1;indices_G_P8(2:end)-indices_G_P8(1:end-1)-1];
distance_G_P9 = [indices_G_P9(1,1)-1;indices_G_P9(2:end)-indices_G_P9(1:end-1)-1];
distance_G_P10 = [indices_G_P10(1,1)-1;indices_G_P10(2:end)-indices_G_P10(1:end-1)-1];
distance_G_P11 = [indices_G_P11(1,1)-1;indices_G_P11(2:end)-indices_G_P11(1:end-1)-1];
distance_G_P12 = [indices_G_P12(1,1)-1;indices_G_P12(2:end)-indices_G_P12(1:end-1)-1];
distance_G_P13 = [indices_G_P13(1,1)-1;indices_G_P13(2:end)-indices_G_P13(1:end-1)-1];
distance_G_P14 = [indices_G_P14(1,1)-1;indices_G_P14(2:end)-indices_G_P14(1:end-1)-1];
distance_G_P15 = [indices_G_P15(1,1)-1;indices_G_P15(2:end)-indices_G_P15(1:end-1)-1];
distance_G_P16 = [indices_G_P16(1,1)-1;indices_G_P16(2:end)-indices_G_P16(1:end-1)-1];
distance_G_P17 = [indices_G_P17(1,1)-1;indices_G_P17(2:end)-indices_G_P17(1:end-1)-1];
distance_G_P18 = indices_G_P18(2:end)-indices_G_P18(1:end-1)-1;
distance_G_P19 = [indices_G_P19(1,1)-1;indices_G_P19(2:end)-indices_G_P19(1:end-1)-1];
distance_G_P20 = [indices_G_P20(1,1)-1;indices_G_P20(2:end)-indices_G_P20(1:end-1)-1];
distance_G_P21 = [indices_G_P21(1,1)-1;indices_G_P21(2:end)-indices_G_P21(1:end-1)-1];
distance_G_P22 = [indices_G_P22(1,1)-1;indices_G_P22(2:end)-indices_G_P22(1:end-1)-1];


[h_G, p_G] = chi2gof([distance_G_P1; distance_G_P2; distance_G_P3; distance_G_P4; distance_G_P5; distance_G_P6; distance_G_P8; distance_G_P7;distance_G_P9; distance_G_P10; distance_G_P11; distance_G_P12; distance_G_P13; distance_G_P14; distance_G_P15; distance_G_P16; distance_G_P17; distance_G_P18; distance_G_P19; distance_G_P20; distance_G_P21; distance_G_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.

%% Proline

indices_P_P1 = [find(LiaP(:,2)==1); length_linkers(:,2)];
indices_P_P2 = [find(LiaP(:,3)==1); length_linkers(:,3)];
indices_P_P3 = [find(LiaP(:,4)==1); length_linkers(:,4)];
indices_P_P4 = [find(LiaP(:,5)==1); length_linkers(:,5)];
indices_P_P5 = [find(LiaP(:,6)==1); length_linkers(:,6)];
indices_P_P6 = [find(LiaP(:,7)==1); length_linkers(:,7)];
indices_P_P7 = [find(LiaP(:,8)==1); length_linkers(:,8)];
indices_P_P8 = [find(LiaP(:,9)==1); length_linkers(:,9)];
indices_P_P9 = [find(LiaP(:,10)==1); length_linkers(:,10)];
indices_P_P10 = [find(LiaP(:,11)==1); length_linkers(:,11)];
indices_P_P11 = [find(LiaP(:,12)==1); length_linkers(:,12)];
indices_P_P12 = [find(LiaP(:,13)==1); length_linkers(:,13)];
indices_P_P13 = [find(LiaP(:,14)==1); length_linkers(:,14)];
indices_P_P14 = [find(LiaP(:,15)==1); length_linkers(:,15)];
indices_P_P15 = [find(LiaP(:,16)==1); length_linkers(:,16)];
indices_P_P16 = [find(LiaP(:,17)==1); length_linkers(:,17)];
indices_P_P17 = [find(LiaP(:,18)==1); length_linkers(:,18)];
indices_P_P18 = [find(LiaP(:,19)==1); length_linkers(:,19)];
indices_P_P19 = [find(LiaP(:,20)==1); length_linkers(:,20)];
indices_P_P20 = [find(LiaP(:,21)==1); length_linkers(:,21)];
indices_P_P21 = [find(LiaP(:,22)==1); length_linkers(:,22)];
indices_P_P22 = [find(LiaP(:,23)==1); length_linkers(:,23)];


distance_P_P1 = [indices_P_P1(1,1)-1;indices_P_P1(2:end)-indices_P_P1(1:end-1)-1];
distance_P_P2 = [indices_P_P2(1,1)-1;indices_P_P2(2:end)-indices_P_P2(1:end-1)-1];
distance_P_P3 = indices_P_P3(2:end)-indices_P_P3(1:end-1)-1;
distance_P_P4 = [indices_P_P4(1,1)-1;indices_P_P4(2:end)-indices_P_P4(1:end-1)-1];
distance_P_P5 = [indices_P_P5(1,1)-1;indices_P_P5(2:end)-indices_P_P5(1:end-1)-1];
distance_P_P6 = [indices_P_P6(1,1)-1;indices_P_P6(2:end)-indices_P_P6(1:end-1)-1];
distance_P_P7 = [indices_P_P7(1,1)-1;indices_P_P7(2:end)-indices_P_P7(1:end-1)-1];
distance_P_P8 = [indices_P_P8(1,1)-1;indices_P_P8(2:end)-indices_P_P8(1:end-1)-1];
distance_P_P9 = [indices_P_P9(1,1)-1;indices_P_P9(2:end)-indices_P_P9(1:end-1)-1];
distance_P_P10 = [indices_P_P10(1,1)-1;indices_P_P10(2:end)-indices_P_P10(1:end-1)-1];
distance_P_P11 = [indices_P_P11(1,1)-1;indices_P_P11(2:end)-indices_P_P11(1:end-1)-1];
distance_P_P12 = [indices_P_P12(1,1)-1;indices_P_P12(2:end)-indices_P_P12(1:end-1)-1];
distance_P_P13 = [indices_P_P13(1,1)-1;indices_P_P13(2:end)-indices_P_P13(1:end-1)-1];
distance_P_P14 = [indices_P_P14(1,1)-1;indices_P_P14(2:end)-indices_P_P14(1:end-1)-1];
distance_P_P15 = [indices_P_P15(1,1)-1;indices_P_P15(2:end)-indices_P_P15(1:end-1)-1];
distance_P_P16 = [indices_P_P16(1,1)-1;indices_P_P16(2:end)-indices_P_P16(1:end-1)-1];
distance_P_P17 = [indices_P_P17(1,1)-1;indices_P_P17(2:end)-indices_P_P17(1:end-1)-1];
distance_P_P18 = [indices_P_P18(1,1)-1;indices_P_P18(2:end)-indices_P_P18(1:end-1)-1];
distance_P_P19 = [indices_P_P19(1,1)-1;indices_P_P19(2:end)-indices_P_P19(1:end-1)-1];
distance_P_P20 = [indices_P_P20(1,1)-1;indices_P_P20(2:end)-indices_P_P20(1:end-1)-1];
distance_P_P21 = [indices_P_P21(1,1)-1;indices_P_P21(2:end)-indices_P_P21(1:end-1)-1];
distance_P_P22 = [indices_P_P22(1,1)-1;indices_P_P22(2:end)-indices_P_P22(1:end-1)-1];


[h_P, p_P] = chi2gof([distance_P_P1; distance_P_P2; distance_P_P3; distance_P_P4; distance_P_P5; distance_P_P6; distance_P_P8; distance_P_P7;distance_P_P9; distance_P_P10; distance_P_P11; distance_P_P12; distance_P_P13; distance_P_P14; distance_P_P15; distance_P_P16; distance_P_P17; distance_P_P18; distance_P_P19; distance_P_P20; distance_P_P21; distance_P_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.

%% Leucine 


indices_L_P1 = [find(LiaL(:,2)==1); length_linkers(:,2)];
indices_L_P2 = [find(LiaL(:,3)==1); length_linkers(:,3)];
indices_L_P3 = [find(LiaL(:,4)==1); length_linkers(:,4)];
indices_L_P4 = [find(LiaL(:,5)==1); length_linkers(:,5)];
indices_L_P5 = [find(LiaL(:,6)==1); length_linkers(:,6)];
indices_L_P6 = [find(LiaL(:,7)==1); length_linkers(:,7)];
indices_L_P7 = [find(LiaL(:,8)==1); length_linkers(:,8)];
indices_L_P8 = [find(LiaL(:,9)==1); length_linkers(:,9)];
indices_L_P9 = [find(LiaL(:,10)==1); length_linkers(:,10)];
indices_L_P10 = [find(LiaL(:,11)==1); length_linkers(:,11)];
indices_L_P11 = [find(LiaL(:,12)==1); length_linkers(:,12)];
indices_L_P12 = [find(LiaL(:,13)==1); length_linkers(:,13)];
indices_L_P13 = [find(LiaL(:,14)==1); length_linkers(:,14)];
indices_L_P14 = [find(LiaL(:,15)==1); length_linkers(:,15)];
indices_L_P15 = [find(LiaL(:,16)==1); length_linkers(:,16)];
indices_L_P16 = [find(LiaL(:,17)==1); length_linkers(:,17)];
indices_L_P17 = [find(LiaL(:,18)==1); length_linkers(:,18)];
indices_L_P18 = [find(LiaL(:,19)==1); length_linkers(:,19)];
indices_L_P19 = [find(LiaL(:,20)==1); length_linkers(:,20)];
indices_L_P20 = [find(LiaL(:,21)==1); length_linkers(:,21)];
indices_L_P21 = [find(LiaL(:,22)==1); length_linkers(:,22)];
indices_L_P22 = [find(LiaL(:,23)==1); length_linkers(:,23)];


distance_L_P1 = [indices_L_P1(1,1)-1;indices_L_P1(2:end)-indices_L_P1(1:end-1)-1];
distance_L_P2 = indices_L_P2(2:end)-indices_L_P2(1:end-1)-1;
distance_L_P3 = indices_L_P3(2:end)-indices_L_P3(1:end-1)-1;
distance_L_P4 = indices_L_P4(2:end)-indices_L_P4(1:end-1)-1;
distance_L_P5 = indices_L_P5(2:end)-indices_L_P5(1:end-1)-1;
distance_L_P6 = indices_L_P6(2:end)-indices_L_P6(1:end-1)-1;
distance_L_P7 = [indices_L_P7(1,1)-1;indices_L_P7(2:end)-indices_L_P7(1:end-1)-1];
distance_L_P8 = indices_L_P8(2:end)-indices_L_P8(1:end-1)-1;
distance_L_P9 = [indices_L_P9(1,1)-1;indices_L_P9(2:end)-indices_L_P9(1:end-1)-1];
distance_L_P10 = [indices_L_P10(1,1)-1;indices_L_P10(2:end)-indices_L_P10(1:end-1)-1];
distance_L_P11 = [indices_L_P11(1,1)-1;indices_L_P11(2:end)-indices_L_P11(1:end-1)-1];
distance_L_P12 = [indices_L_P12(1,1)-1;indices_L_P12(2:end)-indices_L_P12(1:end-1)-1];
distance_L_P13 = [indices_L_P13(1,1)-1;indices_L_P13(2:end)-indices_L_P13(1:end-1)-1];
distance_L_P14 = [indices_L_P14(1,1)-1;indices_L_P14(2:end)-indices_L_P14(1:end-1)-1];
distance_L_P15 = [indices_L_P15(1,1)-1;indices_L_P15(2:end)-indices_L_P15(1:end-1)-1];
distance_L_P16 = [indices_L_P16(1,1)-1;indices_L_P16(2:end)-indices_L_P16(1:end-1)-1];
distance_L_P17 = [indices_L_P17(1,1)-1;indices_L_P17(2:end)-indices_L_P17(1:end-1)-1];
distance_L_P18 = indices_L_P18(2:end)-indices_L_P18(1:end-1)-1;
distance_L_P19 = [indices_L_P19(1,1)-1;indices_L_P19(2:end)-indices_L_P19(1:end-1)-1];
distance_L_P20 = [indices_L_P20(1,1)-1;indices_L_P20(2:end)-indices_L_P20(1:end-1)-1];
distance_L_P21 = [indices_L_P21(1,1)-1;indices_L_P21(2:end)-indices_L_P21(1:end-1)-1];
distance_L_P22 = [indices_L_P22(1,1)-1;indices_L_P22(2:end)-indices_L_P22(1:end-1)-1];


[h_L, p_L] = chi2gof([distance_L_P1; distance_L_P2; distance_L_P3; distance_L_P4; distance_L_P5; distance_L_P6; distance_L_P8; distance_L_P7;distance_L_P9; distance_L_P10; distance_L_P11; distance_L_P12; distance_L_P13; distance_L_P14; distance_L_P15; distance_L_P16; distance_L_P17; distance_L_P18; distance_L_P19; distance_L_P20; distance_L_P21; distance_L_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.

%% Isoleucine 


indices_I_P1 = [find(LiaI(:,2)==1); length_linkers(:,2)];
indices_I_P2 = [find(LiaI(:,3)==1); length_linkers(:,3)];
indices_I_P3 = [find(LiaI(:,4)==1); length_linkers(:,4)];
indices_I_P4 = [find(LiaI(:,5)==1); length_linkers(:,5)];
indices_I_P5 = [find(LiaI(:,6)==1); length_linkers(:,6)];
indices_I_P6 = [find(LiaI(:,7)==1); length_linkers(:,7)];
indices_I_P7 = [find(LiaI(:,8)==1); length_linkers(:,8)];
indices_I_P8 = [find(LiaI(:,9)==1); length_linkers(:,9)];
indices_I_P9 = [find(LiaI(:,10)==1); length_linkers(:,10)];
indices_I_P10 = [find(LiaI(:,11)==1); length_linkers(:,11)];
indices_I_P11 = [find(LiaI(:,12)==1); length_linkers(:,12)];
indices_I_P12 = [find(LiaI(:,13)==1); length_linkers(:,13)];
indices_I_P13 = [find(LiaI(:,14)==1); length_linkers(:,14)];
indices_I_P14 = [find(LiaI(:,15)==1); length_linkers(:,15)];
indices_I_P15 = [find(LiaI(:,16)==1); length_linkers(:,16)];
indices_I_P16 = [find(LiaI(:,17)==1); length_linkers(:,17)];
indices_I_P17 = [find(LiaI(:,18)==1); length_linkers(:,18)];
indices_I_P18 = [find(LiaI(:,19)==1); length_linkers(:,19)];
indices_I_P19 = [find(LiaI(:,20)==1); length_linkers(:,20)];
indices_I_P20 = [find(LiaI(:,21)==1); length_linkers(:,21)];
indices_I_P21 = [find(LiaI(:,22)==1); length_linkers(:,22)];
indices_I_P22 = [find(LiaI(:,23)==1); length_linkers(:,23)];


distance_I_P1 = [indices_I_P1(1,1)-1;indices_I_P1(2:end)-indices_I_P1(1:end-1)-1];
distance_I_P2 = indices_I_P2(2:end)-indices_I_P2(1:end-1)-1;
distance_I_P3 = indices_I_P3(2:end)-indices_I_P3(1:end-1)-1;
distance_I_P4 = indices_I_P4(2:end)-indices_I_P4(1:end-1)-1;
distance_I_P5 = indices_I_P5(2:end)-indices_I_P5(1:end-1)-1;
distance_I_P6 = indices_I_P6(2:end)-indices_I_P6(1:end-1)-1;
distance_I_P7 = indices_I_P7(2:end)-indices_I_P7(1:end-1)-1;
distance_I_P8 = indices_I_P8(2:end)-indices_I_P8(1:end-1)-1;
distance_I_P9 = [indices_I_P9(1,1)-1;indices_I_P9(2:end)-indices_I_P9(1:end-1)-1];
distance_I_P10 = [indices_I_P10(1,1)-1;indices_I_P10(2:end)-indices_I_P10(1:end-1)-1];
distance_I_P11 = [indices_I_P11(1,1)-1;indices_I_P11(2:end)-indices_I_P11(1:end-1)-1];
distance_I_P12 = [indices_I_P12(1,1)-1;indices_I_P12(2:end)-indices_I_P12(1:end-1)-1];
distance_I_P13 = [indices_I_P13(1,1)-1;indices_I_P13(2:end)-indices_I_P13(1:end-1)-1];
distance_I_P14 = indices_I_P14(2:end)-indices_I_P14(1:end-1)-1;
distance_I_P15 = indices_I_P15(2:end)-indices_I_P15(1:end-1)-1;
distance_I_P16 = [indices_I_P16(1,1)-1;indices_I_P16(2:end)-indices_I_P16(1:end-1)-1];
distance_I_P17 = indices_I_P17(2:end)-indices_I_P17(1:end-1)-1;
distance_I_P18 = indices_I_P18(2:end)-indices_I_P18(1:end-1)-1;
distance_I_P19 = indices_I_P19(2:end)-indices_I_P19(1:end-1)-1;
distance_I_P20 = indices_I_P20(2:end)-indices_I_P20(1:end-1)-1;
distance_I_P21 = indices_I_P21(2:end)-indices_I_P21(1:end-1)-1;
distance_I_P22 = indices_I_P22(2:end)-indices_I_P22(1:end-1)-1;


[h_I, p_I] = chi2gof([distance_I_P1; distance_I_P2; distance_I_P3; distance_I_P4; distance_I_P5; distance_I_P6; distance_I_P8; distance_I_P7;distance_I_P9; distance_I_P10; distance_I_P11; distance_I_P12; distance_I_P13; distance_I_P14; distance_I_P15; distance_I_P16; distance_I_P17; distance_I_P18; distance_I_P19; distance_I_P20; distance_I_P21; distance_I_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.

%% Asparagine  


indices_N_P1 = [find(LiaN(:,2)==1); length_linkers(:,2)];
indices_N_P2 = [find(LiaN(:,3)==1); length_linkers(:,3)];
indices_N_P3 = [find(LiaN(:,4)==1); length_linkers(:,4)];
indices_N_P4 = [find(LiaN(:,5)==1); length_linkers(:,5)];
indices_N_P5 = [find(LiaN(:,6)==1); length_linkers(:,6)];
indices_N_P6 = [find(LiaN(:,7)==1); length_linkers(:,7)];
indices_N_P7 = [find(LiaN(:,8)==1); length_linkers(:,8)];
indices_N_P8 = [find(LiaN(:,9)==1); length_linkers(:,9)];
indices_N_P9 = [find(LiaN(:,10)==1); length_linkers(:,10)];
indices_N_P10 = [find(LiaN(:,11)==1); length_linkers(:,11)];
indices_N_P11 = [find(LiaN(:,12)==1); length_linkers(:,12)];
indices_N_P12 = [find(LiaN(:,13)==1); length_linkers(:,13)];
indices_N_P13 = [find(LiaN(:,14)==1); length_linkers(:,14)];
indices_N_P14 = [find(LiaN(:,15)==1); length_linkers(:,15)];
indices_N_P15 = [find(LiaN(:,16)==1); length_linkers(:,16)];
indices_N_P16 = [find(LiaN(:,17)==1); length_linkers(:,17)];
indices_N_P17 = [find(LiaN(:,18)==1); length_linkers(:,18)];
indices_N_P18 = [find(LiaN(:,19)==1); length_linkers(:,19)];
indices_N_P19 = [find(LiaN(:,20)==1); length_linkers(:,20)];
indices_N_P20 = [find(LiaN(:,21)==1); length_linkers(:,21)];
indices_N_P21 = [find(LiaN(:,22)==1); length_linkers(:,22)];
indices_N_P22 = [find(LiaN(:,23)==1); length_linkers(:,23)];


distance_N_P1 = indices_N_P1(2:end)-indices_N_P1(1:end-1)-1;
distance_N_P2 = indices_N_P2(2:end)-indices_N_P2(1:end-1)-1;
distance_N_P3 = [indices_N_P3(1,1)-1;indices_N_P3(2:end)-indices_N_P3(1:end-1)-1];
distance_N_P4 = [indices_N_P5(1,1)-1;indices_N_P4(2:end)-indices_N_P4(1:end-1)-1];
distance_N_P5 = indices_N_P5(2:end)-indices_N_P5(1:end-1)-1;
distance_N_P6 = indices_N_P6(2:end)-indices_N_P6(1:end-1)-1;
distance_N_P7 = indices_N_P7(2:end)-indices_N_P7(1:end-1)-1;
distance_N_P8 = indices_N_P8(2:end)-indices_N_P8(1:end-1)-1;
distance_N_P9 = [indices_N_P9(1,1)-1;indices_N_P9(2:end)-indices_N_P9(1:end-1)-1];
distance_N_P10 = indices_N_P10(2:end)-indices_N_P10(1:end-1)-1;
distance_N_P11 = [indices_N_P11(1,1)-1;indices_N_P11(2:end)-indices_N_P11(1:end-1)-1];
distance_N_P12 = [indices_N_P12(1,1)-1;indices_N_P12(2:end)-indices_N_P12(1:end-1)-1];
distance_N_P13 = [indices_N_P13(1,1)-1;indices_N_P13(2:end)-indices_N_P13(1:end-1)-1];
distance_N_P14 = [indices_N_P14(1,1)-1;indices_N_P14(2:end)-indices_N_P14(1:end-1)-1];
distance_N_P15 = [indices_N_P15(1,1)-1;indices_N_P15(2:end)-indices_N_P15(1:end-1)-1];
distance_N_P16 = [indices_N_P16(1,1)-1;indices_N_P16(2:end)-indices_N_P16(1:end-1)-1];
distance_N_P17 = [indices_N_P17(1,1)-1;indices_N_P17(2:end)-indices_N_P17(1:end-1)-1];
distance_N_P18 = indices_N_P18(2:end)-indices_N_P18(1:end-1)-1;
distance_N_P19 = indices_N_P19(2:end)-indices_N_P19(1:end-1)-1;
distance_N_P20 = [indices_N_P20(1,1)-1;indices_N_P20(2:end)-indices_N_P20(1:end-1)-1];
distance_N_P21 = [indices_N_P21(1,1)-1;indices_N_P21(2:end)-indices_N_P21(1:end-1)-1];
distance_N_P22 = indices_N_P22(2:end)-indices_N_P22(1:end-1)-1;


[h_N, p_N] = chi2gof([distance_N_P1; distance_N_P2; distance_N_P3; distance_N_P4; distance_N_P5; distance_N_P6; distance_N_P8; distance_N_P7;distance_N_P9; distance_N_P10; distance_N_P11; distance_N_P12; distance_N_P13; distance_N_P14; distance_N_P15; distance_N_P16; distance_N_P17; distance_N_P18; distance_N_P19; distance_N_P20; distance_N_P21; distance_N_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.


%% Valine


indices_V_P1 = [find(LiaV(:,2)==1); length_linkers(:,2)];
indices_V_P2 = [find(LiaV(:,3)==1); length_linkers(:,3)];
indices_V_P3 = [find(LiaV(:,4)==1); length_linkers(:,4)];
indices_V_P4 = [find(LiaV(:,5)==1); length_linkers(:,5)];
indices_V_P5 = [find(LiaV(:,6)==1); length_linkers(:,6)];
indices_V_P6 = [find(LiaV(:,7)==1); length_linkers(:,7)];
indices_V_P7 = [find(LiaV(:,8)==1); length_linkers(:,8)];
indices_V_P8 = [find(LiaV(:,9)==1); length_linkers(:,9)];
indices_V_P9 = [find(LiaV(:,10)==1); length_linkers(:,10)];
indices_V_P10 = [find(LiaV(:,11)==1); length_linkers(:,11)];
indices_V_P11 = [find(LiaV(:,12)==1); length_linkers(:,12)];
indices_V_P12 = [find(LiaV(:,13)==1); length_linkers(:,13)];
indices_V_P13 = [find(LiaV(:,14)==1); length_linkers(:,14)];
indices_V_P14 = [find(LiaV(:,15)==1); length_linkers(:,15)];
indices_V_P15 = [find(LiaV(:,16)==1); length_linkers(:,16)];
indices_V_P16 = [find(LiaV(:,17)==1); length_linkers(:,17)];
indices_V_P17 = [find(LiaV(:,18)==1); length_linkers(:,18)];
indices_V_P18 = [find(LiaV(:,19)==1); length_linkers(:,19)];
indices_V_P19 = [find(LiaV(:,20)==1); length_linkers(:,20)];
indices_V_P20 = [find(LiaV(:,21)==1); length_linkers(:,21)];
indices_V_P21 = [find(LiaV(:,22)==1); length_linkers(:,22)];
indices_V_P22 = [find(LiaV(:,23)==1); length_linkers(:,23)];


distance_V_P1 = [indices_V_P1(1,1)-1;indices_V_P1(2:end)-indices_V_P1(1:end-1)-1];
distance_V_P2 = indices_V_P2(2:end)-indices_V_P2(1:end-1)-1;
distance_V_P3 = indices_V_P3(2:end)-indices_V_P3(1:end-1)-1;
distance_V_P4 = indices_V_P4(2:end)-indices_V_P4(1:end-1)-1;
distance_V_P5 = indices_V_P5(2:end)-indices_V_P5(1:end-1)-1;
distance_V_P6 = indices_V_P6(2:end)-indices_V_P6(1:end-1)-1;
distance_V_P7 = indices_V_P7(2:end)-indices_V_P7(1:end-1)-1;
distance_V_P8 = [indices_V_P8(1,1)-1;indices_V_P8(2:end)-indices_V_P8(1:end-1)-1];
distance_V_P9 = [indices_V_P9(1,1)-1;indices_V_P9(2:end)-indices_V_P9(1:end-1)-1];
distance_V_P10 = [indices_V_P10(1,1)-1;indices_V_P10(2:end)-indices_V_P10(1:end-1)-1];
distance_V_P11 = [indices_V_P11(1,1)-1;indices_V_P11(2:end)-indices_V_P11(1:end-1)-1];
distance_V_P12 = indices_V_P12(2:end)-indices_V_P12(1:end-1)-1;
distance_V_P13 = [indices_V_P13(1,1)-1;indices_V_P13(2:end)-indices_V_P13(1:end-1)-1];
distance_V_P14 = [indices_V_P14(1,1)-1;indices_V_P14(2:end)-indices_V_P14(1:end-1)-1];
distance_V_P15 = indices_V_P15(2:end)-indices_V_P15(1:end-1)-1;
distance_V_P16 = [indices_V_P16(1,1)-1;indices_V_P16(2:end)-indices_V_P16(1:end-1)-1];
distance_V_P17 = [indices_V_P17(1,1)-1;indices_V_P17(2:end)-indices_V_P17(1:end-1)-1];
distance_V_P18 = [indices_V_P18(1,1)-1;indices_V_P18(2:end)-indices_V_P18(1:end-1)-1];
distance_V_P19 = [indices_V_P19(1,1)-1;indices_V_P19(2:end)-indices_V_P19(1:end-1)-1];
distance_V_P20 = indices_V_P20(2:end)-indices_V_P20(1:end-1)-1;
distance_V_P21 = indices_V_P21(2:end)-indices_V_P21(1:end-1)-1;
distance_V_P22 = [indices_V_P22(1,1)-1;indices_V_P22(2:end)-indices_V_P22(1:end-1)-1];


[h_V, p_V] = chi2gof([distance_V_P1; distance_V_P2; distance_V_P3; distance_V_P4; distance_V_P5; distance_V_P6; distance_V_P8; distance_V_P7;distance_V_P9; distance_V_P10; distance_V_P11; distance_V_P12; distance_V_P13; distance_V_P14; distance_V_P15; distance_V_P16; distance_V_P17; distance_V_P18; distance_V_P19; distance_V_P20; distance_V_P21; distance_V_P22]); % The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level.


%% Save statistical data about amino acid sequence in txt file: 

amino_acids = ["D"; "E"; "G"; "I"; "K"; "L"; "N"; "P"; "Q"; "R"; "S"; "V"];
h = [h_D; h_E; h_G; h_I; h_K; h_L; h_N; h_P; h_Q; h_R; h_S; h_V];
p = [p_D; p_E; p_G; p_I; p_K; p_L; p_N; p_P; p_Q; p_R; p_S; p_V];

% save in table
fileID = fopen('DNA_sequence_linker_stat_ana.txt','w');


amino_acid_statistics={amino_acids,h,p};
for k1 = 1:size(amino_acids,1);
    fprintf(fileID,'%.4s \t %.4f \t %f \n',char(amino_acid_statistics{1}(k1)), amino_acid_statistics{2}(k1,:),amino_acid_statistics{3}(k1,:))
end
fclose(fileID);
