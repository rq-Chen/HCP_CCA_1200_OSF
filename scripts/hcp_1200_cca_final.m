% Canonical correlation analysis between HCP 1200 subject measures and MINDy parameters
%
% Adapted by Ruiqi Chen, 2023/07/12
%
% Instead of functional connectivity, we use the nonzero entries of MINDy W matrix, decay, and curvature.
%
% HCP 1200 Computational Replication
% Original Code by Stephen Smith, FMRIB Analysis Group, Oxford (https://www.fmrib.ox.ac.uk/datasets/HCP-CCA/)
% Adapted by Nikhil Goyal, National Instite of Mental Health, 2019-2020
% Note, "grot" is just a temp variable used in the code.

% Additional matlab toolboxes required (these are packaged in the 'dependencies/' folder included in the repo)
% FSLNets     http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets
% PALM        http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
% nearestSPD  http://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

% PLEASE NOTE THAT YOU MUST RUN THIS CODE FROM THE /HCP1200/scripts folder!

% Input and Output
nParcels = 200;
indir = '/data/nil-external/ccp/chenr/MINDy_Analysis/data/';
outdir = fullfile('..', 'matlab_outputs', ['MINDy' num2str(nParcels)]);
if ~exist(outdir, 'dir')
  mkdir(outdir)
end

%% Add the dependencies folders to the PATH, and read in necessary data
addpath(genpath('./dependencies/'));

% Read in data, set some variables, create confounds matrix
VARS=readmatrix('../processed_data/VARS.txt');  % Subjects X SMs text file
VARS(:,sum(~isnan(VARS))<130)=NaN;            % Pre-delete any variables in VARS that have lots of missing data (fewer than 130 subjects have measurements)

% Get MINDy parameters
mindy = load(fullfile(indir, ['HCP_Mdl' num2str(nParcels) '.mat']), 'allMdl', 'sublist');
mindy.subs = int32(str2double(extractBetween(mindy.sublist, 'sub', 'Y')));
if nParcels == 100
  load(fullfile(indir, 'atlas', 'Wmask_RC.mat'), 'Wmask');  % Mask for MINDy connectomes
else
  Wmask = true(nParcels);
end
NET = cellfun(@(x) [x.Param{5}(Wmask(:)); x.Param{6}; x.Param{2}]', mindy.allMdl, 'UniformOutput', false);
NET = cell2mat(NET);  % Models from the same subject are concatenated
NET = normalize(NET); % Normalize the parameters since they have very different scales

% Get subjects with both VARS and MINDy
[~, ia, ib] = intersect(VARS(:,1), mindy.subs);
VARS = VARS(ia, :);
NET = NET(ib, :);

clear mindy ia ib

% Number of PCA and CCA components
Nkeep=100;
% Number of permutations
Nperm=10000;

% Set up confounds matrix (this is based on variables selected by Smith et al. in the original study). Confounds are demeaned, any missing data is set to 0
conf=palm_inormal([ VARS(:,[7 14 15 22 23 25]) VARS(:,[265 266]).^(1/3) ]);   % Gaussianise
conf(isnan(conf)|isinf(conf))=0;                % impute missing data as zeros
conf=nets_normalise([conf conf(:,2:end).^2]);   % add on squared terms and renormalise
conf(isnan(conf)|isinf(conf))=0;                % again convert NaN/inf to 0 (above line makes them reappear for some reason)

%% Generate permutation scheme using PALM - for more details see:
% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/ExchangeabilityBlocks#Generating_the_set_of_permutations_only
% in the paper Smith et al. used 100000, but their code says 10,000 should be enough
EB=hcp2blocks('../processed_data/r1200_m.csv', [ ], false, VARS(:,1));  % Input is the raw restricted file downloaded from Connectome DB
PAPset=palm_quickperms([ ], EB, Nperm);                                 % the final matrix of permuations

% Since subjects are dropped by hcp2blocks, we need to drop them from the other matrices (VARS, NET, varsQconf) to avoid errors
subs = EB(:,5);             % Pull list of subject IDs in EB from column 5 (EB is what returned from hcp2blocks), these are the subjects we KEEP
LIA = ismember(VARS,subs);  % LIA = Logical Index Array, an array with logical true (1) where data in VARS is found in subs
rows_keep = LIA(:,1);       % Row #'s to keep (1=keep, 0=drop)

% Now drop all but (rows) of subjects we want to keep
S1 = VARS(rows_keep,:);
N0 = NET(rows_keep,:);
conf = conf(rows_keep,:);

%% Prepare the netmat matrix for the CCA (N1, N2, N3, N4, N5), following steps outlined in Smith et al.
fprintf("Calculating netmat matrices N1 through N5\n")
% N1, formed by 1. demean, 2. globally variance normalize
N1=nets_demean(N0);   % 1. Demean
N1=N1/std(N1(:));     % 2. variance normalize
% We remove the N2 part since our MINDy parameters have to be and have been normalized column-wise
N2 = [];
% N3, formed by horizontally concat N1 and N2
N3=[N1 N2]; % Concat horizontally
% N4, formed by regressing the confounds matrix out of N3
N4=nets_demean(N3-conf*(pinv(conf)*N3));
% N5
[N5,ss1,vv1]=nets_svds(N4,Nkeep); % 100-dim PCA of netmat via SVD reduction

%% Prepare the SM matrix - apply quantitative exclusion criteria
fprintf("Calculating SM matrices S2 through S5\n")

% Remove "bad" SMs, which are defined as:
  % 1. large outlier
  % 2. too many missing (more than 250 subjects missing measurement for an SM)
  % 3. or not enough distinct values (defined as >95% of subjects having the same SM value)
badvars=[];
for i=1:size(S1,2)                          % Iterate over the SMs of S1 (i.e. the columns)
  Xs=S1(:,i);                               % Vector of the values for this SM
  measure_present=~isnan(Xs);                 % How many are elements present? >250 needed (this is a vector where 1=present for a subject)
  Ys=(Xs(measure_present)-median(Xs(measure_present))).^2;  % Of the values present, calculate vector Ys = (Xs -median(Xs))^2, extreme outlier if max(Ys) > 100*mean(Ys) (or max(Ys/mean(Ys)) > 100 is extreme)
  ratio=max(Ys/mean(Ys));
  if (sum(measure_present)>250) && (std(Xs(measure_present))>0) && (max(sum(nets_class_vectomat(Xs(measure_present))))/length(Xs(measure_present))<0.95) && (ratio<100)
      % First criteria: >250 values?
      % Second criteria: std dev of the values > 0?
      % Third criteria: is the size of largest equal-values-group too large? (i.e. > 95% of subjects)
      % Fourth criteria: is there an extreme value?
      % if (1 & 2 & 3 & 4)=True, then keep the SM
    continue; % do nothing
  else
    % A criteria for drop is met, so add the SM to badvars (i = index of the column for an SM)
    badvars=[badvars i];
  end
end

% Get list of the SMs we want to feed into the CCA.
% Found by comparing a list of 1,2,3...478 w/ the indices of the SMs to drop (using setdiff()) to get the indices of SMs to keep
varskeep=setdiff(1:size(S1,2),[1 6 267:457 ...                              % SMs we generally ignore (ID, race, FreeSurfer)
 2 7 14 15 22 23 25 265 266  ...                                              % confound SMs
 11 12 13 17 19 27 29 31 34 40 204 205 212:223 229:233 236 238 242 477 ...    % some more SMs to ignore for the CCA
 3 4 5 8 9 10 16 18 20 21 24 26 28 30 32 33 35:39 458 459 460 463 464 ...     % some more SMs to ignore for the CCA
 badvars]);                                                                   % the "bad" SMs auto-detected above

% Now, prepare the final SM matrix and run the PCA
% S2, formed by gaussianizing the SMs we keep
S2=palm_inormal(S1(:,varskeep)); % Gaussianise

% Now generate S3 (aka varsd), formed by deconfounding the 17 confounds out of S2
S3=S2;
for i=1:size(S3,2) % deconfound ignoring missing data
  grot=(~isnan(S3(:,i)));
  grotconf=nets_demean(conf(grot,:));
  S3(grot,i)=normalize(S3(grot,i)-grotconf*(pinv(grotconf)*S3(grot,i)));
end

% Next, we need to generate S4 (1003x149 matrix)
% First, estimate the SubjectsXSubjects covariance matrix (where for any two subjects, SMs missing for either subject are ignored)
% The approximate covariance matrix (varsdCOV) is then projected onto the nearest valid covariance matrix using nearestSPD toolbox.
% This method avoids imputation, and S4 can be fed into PCA.
varsdCOV=zeros(size(S3,1));
for i=1:size(S3,1) % estimate "pairwise" covariance, ignoring missing data
  for j=1:size(S3,1)
    grot=S3([i j],:);
    grot=cov(grot(:,sum(isnan(grot))==0)');
    varsdCOV(i,j)=grot(1,2);
  end
end
S4=nearestSPD(varsdCOV); % project onto the nearest valid covariance matrix. This method avoids imputation (we can't have any missing values before running the PCA)

% Generate S5, the top 100 eigenvectors for SMs, to avoid overfitting and reduce dimensionality
[uu,dd]=eigs(S4,Nkeep);       % SVD (eigs actually)
S5=uu-conf*(pinv(conf)*uu);   % deconfound again just to be safe

%% CCA
fprintf("Running CCA on matrices S5 and N5\n")
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(N5,S5);

%% Permutation Testing
fprintf("Permutation testing - this may take a while, please be patient\n")
% Use a temporary version of S1 for the null testing
grotvars=palm_inormal(S1);
grotvars(:,std(grotvars)<1e-10)=[];
grotvars(:,sum(~isnan(grotvars))<20)=[];

% permutation testing
grotRp=zeros(Nperm,Nkeep+1);
clear grotRpval;
nullNETr = nan(size(N0, 2), Nperm);
nullSMr = nan(size(grotvars, 2), Nperm);
nullNETv = nan(size(grotU, 2), Nperm);
nullSMv = nan(size(grotV, 2), Nperm);
tic
for j=1:Nperm
  fprintf('Permutation %d\n', j');
  [grotAr,grotBr,grotRp(j,1:end-1),grotUr,grotVr,grotstatsr]=canoncorr(N5,S5(PAPset(:,j),:));
  grotRp(j,end)=mean(grotRp(j,1:end-1));
  nullNETr(:, j) = corr(grotUr(:,1),N0)';
  nullSMr(:, j) = corr(grotVr(:,1),grotvars(PAPset(:,j),:),'rows','pairwise')';
  nullNETv(:, j) = sum(corr(grotUr,N0).^2,2);
  nullSMv(:, j) = sum(corr(grotVr,grotvars(PAPset(:,j),:),'rows','pairwise').^2,2);
end
toc
% figure; plot([mean(grotRp)' prctile(grotRp,95)' prctile(grotRp,99)' prctile(grotRp,99.9)' grotR'])   % final point is mean of all Nkeep R values
grotRpval=zeros(1,Nkeep);
for i=1:Nkeep  % show corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of significant CCA components
grotRpval(1:Ncca)

% Based on null distributions, what are the thresholds for significance in % variance explained?
prctile( max(abs(nullSMr)) ,95)
prctile( max(abs(nullNETr)) ,95)

%% Calculate all CCA Mode weights for netmats and SMs
% Netmat weights for CCA modes
grotAA = corr(grotU,N0)';
 % or
grotAAd = corr(grotU,N4(:,1:size(N0,2)))'; % weights after deconfounding

%%% SM weights for CCA modes
grotBB = corr(grotV,palm_inormal(S1),'rows','pairwise')';
 % or 
varsgrot=palm_inormal(S1);
for i=1:size(varsgrot,2)
  grot=(~isnan(varsgrot(:,i))); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grotBBd = corr(grotV,varsgrot,'rows','pairwise')'; % weights after deconfounding

%% Plot SMs vs. Connectome edges (canonical variables)
% Color plot using fluid intelligence values for plot
SMs=importdata('../data/column_headers.txt');  % Load in the column headers file (list of 478 SMs)
index=find(contains(SMs,'PMAT24_A_CR'));    % Get column index of 'PMAT24_A_CR' (fluid intelligence), so we can find the correct column in S1
g_f=S1(:,index);                            % g_f vector: fluid intelligence values for 1003 subjects for the plot of CCA weights

% Plot of CCA Weights (SMs) vs. CCA Weights (connectomes)
figure;
g_f_n = palm_inormal(g_f);
scatter(grotU(:,1),grotV(:,1),15, g_f_n,'filled') % NOTE: negated vectors are plotted
cmap=colormap(jet(20));               % 20 shades from the 'jet' colormap option
cb=colorbar();
cb.Ticks = [min(g_f_n),max(g_f_n)];   % min and max ticks for the colorbar taken from the fluid intelligence values (min 4, max 24)
cb.TickLabels = {'4','24'};
set(cb,'FontSize',17);
ylabel(cb,'Fluid Intelligence Score');
th = title('Scatter plot of SM weights versus connectome weights');
titlePos = get(th,'position');
titlePos1 = titlePos + [0 0.15 0];
set(th, 'position', titlePos1);
xlabel('CCA Weights (Connectomes)');
ylabel('CCA Weights (Subject measures)');
set(gca,'FontSize',16)

% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
width=6;
height=5;
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

print('-dpdf',fullfile(outdir,'SMs_vs_edges.pdf'));
saveas(gcf, fullfile(outdir,"SMs_vs_edges.png"));

figure
scatter(grotU(:,1),grotV(:,1))
xlabel('CCA scores (connectomes)');
ylabel('CCA scores (subject measures)');
lsline();
title(sprintf('Correlation = %.2f', corr(grotU(:,1),grotV(:,1))));
saveas(gcf, fullfile(outdir, "regression.png"));

%% Variance analyses - how much of the total percent variance do the first 20 CCA modes explain?
grot_NET=[ sum(corr(grotU,N0).^2,2) prctile(nullNETv,5,2) mean(nullNETv,2) prctile(nullNETv,95,2) sum(corr(N5,N0).^2,2) ] * 100 / size(N0,2);
grot_VARS=[ sum(corr(grotV,grotvars,'rows','pairwise').^2,2) prctile(nullSMv,5,2) mean(nullSMv,2) prctile(nullSMv,95,2) ...
  sum(corr(S5,grotvars,'rows','pairwise').^2,2)  ] * 100 / size(grotvars,2);
I=1:20;

% Connectomes variance
figure;
subplot(2,1,1); 
hold on;
% Draw the rectangles for null distributions per mode
for i=1:length(I)
  rectangle('Position',[i-0.5 grot_NET(i,2) 1 grot_NET(i,4)-grot_NET(i,2)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
plot(grot_NET(I,3),'k'); plot(grot_NET(I,1),'b'); plot(grot_NET(I,1),'b.'); % plot(grot_NET(I,5),'g');  % turned off showing the PCA equivalent plots
% th1 = title({'Connectome total % variance explained by CCA modes';''});
ylabel('%% variance (connectomes)')
xlabel('CCA Mode')
xlim([1 20])
% ylim([0.3 0.9])
% yticks([0.3 0.35 0.4 0.45 0.5 0.55])
set(gca,'FontSize',15)

% Subject measures variance
subplot(2,1,2); 
hold on;
for i=1:length(I)
  rectangle('Position',[i-0.5 grot_VARS(i,2) 1 grot_VARS(i,4)-grot_VARS(i,2)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
plot(grot_VARS(I,3),'k'); plot(grot_VARS(I,1),'b'); plot(grot_VARS(I,1),'b.'); % plot(grot_VARS(I,5),'g');
% th2 = title({'';'Subject measures total % variance explained by CCA modes';''});
ylabel('%% variance (SMs)')
xlabel('CCA Mode')
xlim([1 20])
% ylim([0 2])
% yticks([0 0.5 1 1.5 2])
set(gca,'FontSize',15)

% Here we preserve the size of the image and save it (pdf and png)
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
width=6;
height=5;
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('-dpdf',fullfile(outdir, 'CCAvarianceexplained.pdf'));
saveas(gcf, fullfile(outdir, "CCAvarianceexplained.png"));


%% Output important results to file and terminal
% Calculate Z-scores for the CCA modes
z_scores_NET=sort((sum(corr(grotU,N0).^2,2)-  mean(nullNETv,2)) ./ std(nullNETv,[],2))';
z_scores_SM=sort((sum(corr(grotV,grotvars,'rows','pairwise').^2,2) - mean(nullSMv,2)) ./ std(nullSMv,[],2))';

% Calculate the r value for CCA SM weights vs. connectome weights
reg_coef = corrcoef(grotU(:,1),grotV(:,1));
reg_coef = norm(reg_coef(1,2));

% Print to file
id=fopen(fullfile(outdir, 'results.txt'),'a'); % append if file exists
str=string(datetime('now','Format','yyyy-MM-dd''_T''HH-mm-ss')); %datetime format ex. "2020-02-11_T11-21-16"
fprintf(id,'------ HCP1200 CCA Ran at %s ------\n',str)
fprintf(id,'Number of permutations: %d\n', Nperm)
fprintf(id,"Pearson's r coefficient: %.4f\n", reg_coef)
fprintf(id,'Number of FWE-significant CCA components: %d\n', Ncca)
fprintf(id,'Number subjects fed into CCA: %d\n',size(N0,1))
fprintf(id,'Number of subject measures fed into CCA: %d\n', size(varskeep,2))
fprintf(id,'%% connectome variance explained by mode 1: %.4f\n',grot_NET(1,1))
fprintf(id,'%% SM variance explained by mode 1: %.4f\n', grot_VARS(1,1))
fprintf(id,'%% connectome variance threshold: %.4f\n', prctile( max(abs(nullNETr)) ,95))
fprintf(id,'%% SM variance threshold: %.4f\n', prctile( max(abs(nullSMr)) ,95))
fprintf(id,'Z-Scores for connectome variance:\n')
fprintf(id,'\thighest: %.4f\n \tsecond highest:%.4f\n \tratio:%.4f\n', max(z_scores_NET),max(z_scores_NET(1:99)),max(z_scores_NET)/max(z_scores_NET(1:99)))
fprintf(id,'Z-Scores for SM variance:\n')
fprintf(id,'\thighest: %.4f\n \tsecond highest:%.4f\n \traio:%.4f\n', max(z_scores_SM),max(z_scores_SM(1:99)), max(z_scores_SM)/max(z_scores_SM(1:99)))
fclose(id);

%% Save workspace
str=string(datetime('now','Format','yyyy-MM-dd''_T''HH-mm-ss')); %datetime format ex. "2020-02-11_T11-21-16"
fname = fullfile(outdir, "hcp_1200_cca_MINDy_"+str);
save(fname, '-v7.3')