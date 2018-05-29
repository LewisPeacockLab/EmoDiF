AddPsychJavaPath;

%{

Experiment: IMDiF
Type      : Item-Method Directed Forgetting
Date      : February 26th, 2015
Version   : 4
Author    : Katerina Placek; Tracy Wang amended
Contact   : tracy.wang@utexas.edu

** ONE DAY STUDY **
* Timing
  - Localizer - 1-Back Task
    - Run 1-3
      - 5 mins per run = 15 mins
	- Study - (about 1 our)
			- Run 1-8
				- ~6.3 min per run = ~50.4s
    - Test Phase
      - self guided

- Experiment Details -
Localizer 
   Faces, Scenes and Objects presented in miniblocks
   - 90 Faces, 90 Scenes, 90 Object and 90 Rest Trials. - 30 of each
   subcategory Male/Female(Face), Indoor/Outdoor(Scene), and
   Artificial/Manmade(Object). 

Study 
    Faces and Scenes presented with Remember or Forget Cue - See Gdrive for
    document of counterbalancing.
        336 Stimuli 
            168 Scenes - 84 Indoor/ 84 Outdoor
            168 Faces  - 84 Male/ 84 Female
Test
    420 Stimuli
        336 Repeated from Study
        168 New 
            84 Scenes - 42 Indoor/ 42 Outdoor
            84 Faces  - 42 Male/ 42 Female

Conditions of Interest:
* 208 trials; 104 Scene, 104 Face - Remember(R)
* 128 trials; 64 Scene, 64 Face - Forget(F); 
    - 104 trials with Post F opportunity; 52 Scene, 52 Face

- Experimental Design -

1) Localizer
Trial Timing
        2s per trial - 1.5s STIM, .5s ITI
        9 trials per miniblock, 4 miniblocks (1 per category) per block
        4 BLOCKS per run (6s rest at the start and end of each block.
           *****MINIBLOCKS ARE TREATED LIKE TRIALS ********

2) Study Phase (BEHAVIORAL/SCANNED)
	- 8 runs
   	- Timing
        9s per trial- 3s STIM, 6s R/F CUE / ITI
   	- Run Setup
      - 42tr per run, 336 totaltr

3) Test Phase (BEHAVIORAL)
	- Timing
        - image displayed on screen; self-guided response

TODO: Must structure trials to demonstrate 1-2R for each block depending on
dependencies.
  %}

exp_dir = imdif_addpaths;

%where the experiment stim directory is
exp_stim_dir = '~/Dropbox (LewPeaLab)/DATA/imdif/experiment';

create_corresp = @(max_val, nConds) Shuffle(mod(1:max_val, nConds)+1)';

% Get subject info

header = struct('experiment', 'imdif', 'date', datestr(now));
header.version = 2;
header.setuptime = fix(clock);
header.exp_dir = exp_dir;
header.data_dir = '~/Dropbox (LewPeaLab)/DATA/imdif/experiment';

fprintf('\nSubject Information\n-------------------\n');
header.subnum    = input('Subject Number:  ');
header.subinit   = input('Subject Initials:  ', 's');
header.subgender = input('Is this male (1) or female (2)?  ');
header.datestr   = datestr(now,'yymmdd');
par.subscreen    = input('Screen Number (0:main; 1:external):  ');
par.fmri_mode    = input('fMRI mode? (0:no; 1:yes):  ');

fprintf('\n');
disp('---------------------------------------------');
disp('Please make sure this information is correct.');
disp('---------------------------------------------');
disp(['Subject Nr = ', num2str(header.subnum)]);
disp(['Subject Initials = ',header.subinit]);
disp(['Subject DateStr = ',header.datestr]);

if header.subgender == 1; 
    disp('This subject is male');
else
    disp('This subject is female');
end
disp(['Screen Nr = ', num2str(par.subscreen)]);  
disp(['fMRI Mode = ', num2str(par.fmri_mode)]);  
disp('---------------'); disp(''); disp('');

yn = input('Is this correct?  (y,n):  ', 's');

if(~isequal(upper(yn(1)), 'Y'))
    error('Please RE-RUN imdif_header!')
end

% Initialize rand seed - computes a random number (got rid of legacy
% rand('seed', sd))

rng(mod(header.subnum*1e6*pi,2^32)); % make sure the randomization is not the same every time.


% Set psychtoolbox parameters
par.fixcolor   = 0; %white
par.respcolor  = 1;
par.backcolor  = 128; %gray
par.txtcolor   = 0;
par.txtsize    = 100;
par.txtsizesm  = 50;
par.txtsizex   = 75;
par.labelsize  = 35;
par.respkey1   = '1!';
par.respkey2   = '2@'; 
par.respkey3   = '3#';
par.respkey4   = '4$';
par.confkeys   = {'1!', '2@', '3#', '4$'}; % keys used to indicate choice in test %

% Set General parameters (across all phases)
par.nptrials        = 18; % number of practice trials
par.nlptrials       = 36; % number of localizer practice trials
par.nlcats          = 4; % face, scene, object and rest
par.nlsubcats       = 8; % stimulus subcategories = natural/artificial, male/female, indoor/outdoor, rest/rest
par.nltrials        = 432; % number of localizer trials
par.nltrials_norpt  = 360; % number of unique localizer trials
par.nlmb            = 48;  % number of localizer mini-blocks
par.nlreptype       = 2; % number of repetition types (i.e. 1 and 2)
par.nlruns          = 3; % number of localizer runs
par.nltripersubcat  = par.nltrials/par.nlsubcats; %number of stims for subcategory
par.nmbtrials       = 9; % number of trials per miniblock - not the number of unique stims
par.reptype1        = 8; % number of unique stims for repetition of 1 item per miniblock
par.reptype2        = 7; % number of unique stims for reptition of 2 items per miniblock

par.sstimnum       = 504;
par.ncats          = 2; % face scene
par.nsubcats       = 2; % male/female, indoor/outdoor - CONFUSING because see above '6' subcategory/ shouldn't this be 4?
par.ntrials        = 504; %total study AND test trials
par.ntrifor        = 96; %number of forget trials
par.ntriforpost    = 78; %number of forget trials with post R opportunities
par.ntrirem        = 156; %number of remember trials
par.ntripersubcat  = par.ntrials/(par.ncats*par.nsubcats); %126(imgspersubcat)

% Set Practice parameters (YET TO DO)!!!!!!!!! OR possibly just use
% localizer paramters

% Set Specific Localizer parameters 
% Set Trial-type TO miniblocks - called miniblock trials
par.local.type = 'N-back (1-back) test';
par.local.nruns    = 3 ;%number of runs 
par.local.mb= 16; %number of miniblocks per run
par.local.mbtrial = 4; % number of miniblocks per block
par.local.block = 4; % number of blocks per run
par.local.trtime = 2; 
par.local.objects = 1;
par.local.faces = 2;
par.local.scenes = 3;
par.local.rest = 4;
% for contingency matrix
par.local.npres    = par.nltrials/par.nlruns; %number of trials PER run
par.local.stim     = 1.5; %number of TRs or seconds that the stimulus is on the screen
par.local.iti     = .5; %the ITI or response time
% par.local.rest = 6; %rest time at the beginning and at the end of a trial
% already used WaitSecs to do this
par.local.mbtrialtime = par.local.mbtrial*(par.local.stim + par.local.iti) ; 
par.local.runtime = par.local.rest*2 + (par.local.mbtrialtime*par.local.mb);
% Set Study Phase parameters
par.study.nruns  = 6; 
par.study.npres  = (par.ntrifor+par.ntrirem)/par.study.nruns;
par.study.stim   = 3; %stim presentation TR number
par.study.rf     = 6; %remember/forget cue TR number

par.study.totaltime = par.study.stim + par.study.rf; 

% Set Test Phase parameters
par.test.type    = 'Self-paced four-alternative forced-choice';
par.test.ntrials = par.ntrials;
par.test.newitems  = par.ntrials/2;



header.par = par;

% %%%%%%%%
%%%%%%%%
%%%%%%%%
%LOCALIZER HEADER%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%%% starting with BLOCKS

% LOCALIZER PARAMETERS
%   header.lmaster info
%   category: 1= Obj, 2= Face, 3= Scene
%   subcategory: 1= NO (Natural Object) , 2= AO (Artifical Object)
%   3= Female Face, 4= Male Face, 5= Inside Scene, 6= Outside Scene
%   Stim Rpt - 1 - 1st veiw, 2 - repeated view.
%      
%  

% header.master info
%   category: 1 = Face, 2 = Scene
%   subcategory: 1 = Female, Indoor, 2 = Male, Outdoor
%   remember/forget: 1 = remember, 2 = forget 


% CREATE Localizer Master Stimlist (HEADER.LMASTER) creates master stimlist, grabs locations
local_mastermat = local_master(par.nltrials_norpt, par.nlsubcats);
% CONVERTS master stimlist to TABLE format. 
header.lmaster = cell2table(local_mastermat, 'VariableNames',{'Orig_trial_Number' 'Category' 'Subcategory' 'namelocation' 'subcategorynum'});

% Localizer Run Setup
%ASSIGNS stims BY subcategory to the different runs (default number of runs
%= 3)
local_runmat = lseparate_run(header.lmaster, par.local.nruns, par.nltrials_norpt, par.nlsubcats);
header.lruns = local_runmat;

%assign repetition and expand miniblock - in ABOVE function

header.lruns = prep_repeats(header.lruns,par.nltrials_norpt, par.nlmb, par.nlsubcats, par.reptype1, par.reptype2);

% items here are now organized by category and prepared by miniblock. index
% by miniblock, use contingency matrix to prepare blocks and runs. add a
% stimlist to each run

%adds a stimlist field to header
% tmp=cell(size(header.lruns)); [header.lruns(:).stimlist]=deal(tmp{:});
% clearvars tmp
%1 = object, 2 = faces, 3 = scenes and 4 = rest.  this randomizes the
%categories and randomly permutates within each miniblock.  It iterates
%until there are no two category miniblocks of the same type that follow
%each other. 

local_stimlist = prep_lstimlist(header.lruns, par.nlcats, par.nlruns);
header.lruns = local_stimlist;
ltrialnum = (1:144)';
ltrialnum = table(ltrialnum);
for i = 1:3
    header.lruns(i).stimlist = horzcat(ltrialnum, header.lruns(i).stimlist);
end

%FINISHED WITH LOCALIZER STIMS

%%%%%%%%
%%%%%%%%
%Study/Test MASTER stilist) 
% CREATE Master Study/Test Stimlist (HEADER.MASTER)  
% creates matrix with cat, subcat, R/F/N assignments
mastermat = create_master(par.sstimnum, par.nsubcats);
header.master = mastermat;

%for each sequence type - randomly matched with a sequence type that is
%within 'set' - to offset unequal numbers of faces and scenes for each
%sequence. 

header.sequence = sequence_assign;

%fill in each run with items from mastermat 


st_stimlist = prep_run_mat(header);

%assigns output of prep_run_mat to the structure sruns (as opposed to
%localizer_runs (or lruns)

header.sruns = st_stimlist;

%clears temporary output of sequence assignation. 

clear header.sequence;

%clears some variables out of study_stimlists (randomization numbers)

for i = 1:6
    
header.sruns(i).stimlist(:,2) = [];
header.sruns(i).stimlist(:,6) = [];
header.sruns(i).stimlist(:,6) = [];
header.sruns(i).stimlist = cell2table(header.sruns(i).stimlist, 'VariableNames',{'trialnumber' 'Category' 'Subcategory' 'RememberForget' 'item_name'});
end




% % Test Phase
% % Create the test parameters
rand_trials   = randperm(par.test.ntrials);
testtrials    = header.master;
header.trun   = testtrials(rand_trials', :);

% header.trun_temp = header.master;
% header.trun_temp(:,6) = []
% 
% if mod(header.subnum,2) %tests if the subjNum is even or odd
%     %if odd this will run
%     ttemp1 = header.trun_temp(85:115,:); % this gets 31 male faces
%     ttemp2 = header.trun_temp(211:242,:); % this gets 32 female faces
%     ttemp3 = header.trun_temp(337:367,:); % this gets 31 indoor scenes
%     ttemp4 = header.trun_temp(463:494,:); % this gets 32 outdoor scenes
% else
%     %if even this will run
%     ttemp1 = header.trun_temp(85:116,:); % this gets 32 male faces
%     ttemp2 = header.trun_temp(211:241,:); % this gets 31 female faces
%     ttemp3 = header.trun_temp(337:368,:); % this gets 32 indoor scenes
%     ttemp4 = header.trun_temp(463:493,:); % this gets 31 outdoor scenes
% end
% 
% header.trun = vertcat(ttemp1, ttemp2, ttemp3, ttemp4);
% header.trun = vertcat(header.sruns(1).stimlist, header.sruns(2).stimlist, header.sruns(3).stimlist,header.sruns(4).stimlist,header.sruns(5).stimlist,header.sruns(6).stimlist, header.trun); 
% 
% clear header.trun(:,1)

%giving a trial idx before random permutation
trialidx = (1:par.test.ntrials)';
trialnidx = num2cell(trialidx);
header.trun(:,1)=trialnidx;

%random permutation
repeatedcat = 5;
repeatedcue = 5;

while ~(repeatedcat < 5 && repeatedcue < 5);
    randnum = randperm(504)';
    randnum = num2cell(randnum);
    header.trun(:,6)=randnum;
    header.trun = sortrows(header.trun,6);
    
    %pseudorandomize test trials - look for sequences fo 3
    repeatcat = header.trun(:,3);
    repeatcue = header.trun(:,4);%
    
%     repeatcat = table2cell(repeatcat);
%     repeatcue = table2cell(repeatcue);
    repeatcat = cell2mat(repeatcat);
    repeatcue = cell2mat(repeatcue);
    repeatedcue = max(diff(find([1,diff(repeatcue'),1]))); 
    repeatedcat = max(diff(find([1,diff(repeatcat'),1]))); % looks for max amount of repeats
    % finds where repeats are >3, indexes where they are.  
    % randomizes the following 15 places
    % repeats check
    
        if repeatedcat > 4
            where = ([1,diff(repeatcat'),1]);
            index = find(where);
            rptnum = diff(index);
            rptidx = find(rptnum>4);
            m = index(rptidx);
            for i = 1:numel(m);
                if m<490
                    randsub = randperm(15)';
                    randsub = num2cell(randsub);
                    header.trun(m(i):(m(i)+14),6) = randsub;
                    temp = header.trun(m(i):(m(i)+14),:);
                    temp = sortrows(temp,6);
                    header.trun(m(i):(m(i)+14),:) = temp;
                else
                    randsub = randperm(15)';
                    randsub = num2cell(randsub);
                    header.trun(490:504,6) = randsub;
                    temp = header.trun(490:504,:);
                    temp = sortrows(temp,6);
                    header.trun(490:504,:) = temp;
                end
                    
            end
        else
            where = ([1,diff(repeatcue'),1]);
            index = find(where);
            rptnum = diff(index);
            rptidx = find(rptnum>4);%
            m = index(rptidx);
            for i = 1:numel(m);
                if m <490
                    randsub = randperm(15)';
                    randsub = num2cell(randsub);
                    header.trun(m(i):(m(i)+14),6) = randsub;
                    temp = header.trun(m(i):(m(i)+14),:);
                    temp = sortrows(temp,6);
                    header.trun(m(i):(m(i)+14),:) = temp;
                else
                    randsub = randperm(15)';
                    randsub = num2cell(randsub);
                    header.trun(490:504,6) = randsub;
                    temp = header.trun(490:504,:);
                    temp = sortrows(temp,6);
                    header.trun(490:504,:) = temp;
                end
            end 
            end
            
    repeatcat = header.trun(:,3);
    repeatcue = header.trun(:,4);%
    
%     repeatcat = table2cell(repeatcat);
%     repeatcue = table2cell(repeatcue);
    repeatcat = cell2mat(repeatcat);
    repeatcue = cell2mat(repeatcue);                
    repeatedcue = max(diff(find([1,diff(repeatcue'),1]))); 
    repeatedcat = max(diff(find([1,diff(repeatcat'),1])));
    
    %Check to make sure no 3 of the same subcattype appear in a row (scene or face) -
end


%giving the list a trial number

trialnum = (1:par.test.ntrials)';
trialnum = num2cell(trialnum);
header.trun(:,1)= trialnum;


header.trun(:,6)= [];
header.trun = cell2table(header.trun, 'VariableNames',{'trialnumber' 'Category' 'Subcategory' 'RememberForgetNew' 'item_name'});
% header.trun.Properties.VariableNames{'RememberForget'} = 'RememberForgetNew';

%prep into stimlists
% % 
% % %CREATES ALL THREE PRACTICE SCRIPTS -
%practice localizer
%practice study 
%practice test

% this selects the practice trials and randomly sorts
prac_mastermat = prac_create_master(par.nlsubcats);
header.pmaster = prac_mastermat;



%creates localizer practice
prac_localmaster = local_practlist(prac_mastermat);
header.prac.lrun = prac_localmaster;
header.prac.lrun = cell2table(header.prac.lrun, 'VariableNames',{'Category' 'Subcategory' 'Item' 'randrep_assign' 'repetitionassign'});

%creates studytest master

prac_studytest_mastermat = create_studytestlist(header.pmaster);

header.prac.srun = prac_studytest_mastermat(1:12,:); 

header.prac.srun = cell2table(header.prac.srun, 'VariableNames',{'Category' 'Subcategory' 'Item' 'RF'});

prac_testlist = prac_studytest_mastermat;
randtest = randperm(length(prac_studytest_mastermat))'; 
randtest = num2cell(randtest);
prac_testlist = horzcat(randtest, prac_studytest_mastermat);
prac_testlist = sortrows(prac_testlist, 1);
header.prac.trun = prac_testlist;
header.prac.trun = cell2table(header.prac.trun, 'VariableNames',{ 'TrialNum' 'Category' 'Subcategory' 'Item' 'RFN' });

% % rand_trials   = randperm(18);
% % testtrials    = header.pmaster;
% % header.ptest  = testtrials(rand_trials', :);
% % 
% % %practice localizer - this selects the practice trials and randomly sorts
% % local_prac_master = local_prac_master(par.nlptrials,par.nlsubcats); 
% % 
% % 
% % rand_trials = randperm(12);
% % lptrials = local_prac_master;
% % header.plocal = lptrials(rand_trials', :);



% Save header
header.subpath = sprintf('%s/data/sub_%02d_%s', header.data_dir, header.subnum, header.datestr);

if ~exist(header.subpath, 'dir')
  mkdir(header.subpath)
else
  warning('Subject folder already exists!')
end

savefname = sprintf('%s/imdif_header_%02d_%s', header.subpath, header.subnum, header.datestr);
save(savefname,'header');

addpath(genpath(header.subpath))

clearvars -except header;    
