function beh_analysis(subs)

subs = [1:10];
n = length(subs);

% Dirs
d.base = '~/Documents/Classwork/Epstein_Rotation/LocalGeometry';
d.data = [d.base '/Behavior/Data'];

% Concatentate subject data
data = [];
for i = 1:length(subs)
    % Filename
    d.sub{i} =  sprintf('%s/%02d',d.data,subs(i));
    fn = sprintf('sub_%02d*.dat',subs(i));
    files = dir([d.sub{i} '/' fn]);
    
    % Find largest file
    if length(files) > 1
        disp(sprintf(['WARNING: multiple files for subject %02d... ' ...
                      'using largest one...'],subs(i)));
        fn = [files(find([files.bytes] == ...
                         max([files.bytes]))).name];
    else
        fn = files.name;
    end
    
    datafile = [d.sub{i} '/' fn];
    dat = importdata(datafile,' ');
    % tmp = cellfun(@str2num,dat.textdata(2:end,1:3));
    
    % data = [data; [tmp dat.data]];
    % objects = dat.textdata(:,4);
    labels = dat.textdata{1};
    data = [data; dat.data];
end

%% ERROR RATES
% DP
dp_idx = data(:,8) > 2;
dp_mean = grpstats(data(dp_idx,18),data(dp_idx,1));

ndp_idx = data(:,8) == 1;
ndp_mean = grpstats(data(ndp_idx,18),data(ndp_idx,1));

% Attend
att_idx = data(:,10) == 1;
att_mean = grpstats(data(att_idx,18),data(att_idx,1));

nat_idx = data(:,10) == 0 & data(:,5) == 0;
nat_mean = grpstats(data(nat_idx,18),data(nat_idx,1));

% Maze
m1_idx = data(:,7) == 1;
m1_mean = grpstats(data(m1_idx,19),data(m1_idx,1));

m2_idx = data(:,7) == 2;
m2_mean = grpstats(data(m2_idx,19),data(m2_idx,1));


%% RT
% DP
dp_RT = grpstats(data(dp_idx,15),data(dp_idx,1));
ndp_RT = grpstats(data(ndp_idx,15),data(ndp_idx,1));
[h, p] = ttest2(dp_RT,ndp_RT);

% Attend
att_RT = grpstats(data(att_idx,15),data(att_idx,1));
nat_RT = grpstats(data(nat_idx,15),data(nat_idx,1));

% Maze
m1_RT = grpstats(data(m1_idx,16),data(m1_idx,1));
m2_RT = grpstats(data(m2_idx,16),data(m2_idx,1));

dp_data = [grpstats(data(att_idx & dp_idx,15), data(att_idx & dp_idx,1)); ...
           grpstats(data(nat_idx & dp_idx,15), data(nat_idx & dp_idx,1))];
ndp_data = [grpstats(data(att_idx & ndp_idx,15), data(att_idx & ndp_idx,1)); ...
           grpstats(data(nat_idx & ndp_idx,15), data(nat_idx & ndp_idx,1))];
[p,tbl,stats] = anova2([dp_data ndp_data],n);



































% Initial Self Analysis
if 1 == 2
    %% ERROR RATES
    % DP
    dp_mean  = mean(data(data(:,7) > 2,17));
    dp_std   = std(data(data(:,7) > 2,17));
    % nDP
    ndp_mean = mean(data(data(:,7) == 1,17));
    ndp_std  = std(data(data(:,7) == 1,17));
    % t-test
    [h,p_dp,ci_dp,stats] = ttest(data(data(:,7) == 1,17),data(data(:,7) > 2,17));

    % Attended
    att_mean = mean(data(data(:,9) == 1,17));
    att_std  = std(data(data(:,9) == 1,17));
    % Unattended
    uat_mean = mean(data(data(:,9) == 0 & data(:,4) == 0,17));
    uat_std  = std(data(data(:,9) == 0 & data(:,4) == 0,17));
    [h,p_att,ci_att] = ttest(data(data(:,9) == 0 & data(:,4) == 0,17), ...
                             data(data(:,9) == 1,17));

    % Mazes
    [h,p_maze,ci_maze] = ttest(data(data(:,9) == 0


    %% REACTION TIMES
    dp_data = [data(data(:,7) > 2 & data(:,9) == 1,17); 
               data(data(:,7) > 2 & data(:,9) == 0 & data(:,4) == 0, ...
                    17)];
    ndp_data = [data(data(:,7) == 1 & data(:,9) == 1,17); 
                data(data(:,7) == 1 & data(:,9) == 0 & data(:,4) == 0, ...
                     17)];
    [p,tbl,stats] = anova2([dp_data ndp_data],192/2);
end





    


