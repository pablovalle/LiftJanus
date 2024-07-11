clc; clear; close all;
initialPath="..\Experiments_withDD\ExperimentsScenario1";
addpath(initialPath+"\guided\");
addpath(initialPath+"\unguided\");

numOfExecs=9;
numOfIterations=99;


%% Get the data from the excel files

for i=2:numOfExecs
    path_guided_exec=initialPath+"\guided\"+num2str(i)+"\report"+num2str(i)+".csv";
    path_unguided_exec=initialPath+"\unguided\"+num2str(i)+"\report"+num2str(i)+".csv";
    path_guided_archive=initialPath+"\guided\"+num2str(i)+"\results"+num2str(i)+".csv";
    path_unguided_archive=initialPath+"\unguided\"+num2str(i)+"\results"+num2str(i)+".csv";
    
    report_guided(i-1,:,:)=table2array(importReport(path_guided_exec));
    report_unguided(i-1,:,:)=table2array(importReport(path_unguided_exec));

    archive_guided(i-1,:,:)=table2array(importArchiveResults(path_guided_archive));
    archive_unguided(i-1,:,:)=table2array(importArchiveResults(path_unguided_archive));
end

%% Comparision of guided vs unguided for the on-running data

% Figure Min LWT
min_LWT_guided=mean(report_guided(:,:,7));
min_LWT_unguided=mean(report_unguided(:,:,7));

figure; plot(min_LWT_unguided,"r--"); hold on ;title("Min Longest Waiting Time Difference between members");xlabel("# of iterations");ylabel("LWT (s)"); plot(min_LWT_guided,"--b");legend('unguided','guided');

% Figure Max LWT
max_LWT_guided=mean(report_guided(:,:,8));
max_LWT_unguided=mean(report_unguided(:,:,8));

figure; plot(max_LWT_unguided,"r--"); hold on;title("Max Longest Waiting Time Difference between members");xlabel("# of iterations");ylabel("LWT (s)"); plot(max_LWT_guided,"--b");legend('unguided','guided');

% Figure AVG LWT
avg_LWT_guided=mean(report_guided(:,:,9));
avg_LWT_unguided=mean(report_unguided(:,:,9));

figure; plot(avg_LWT_unguided,"r--"); hold on;title("Average Longest Waiting Time Difference between members");xlabel("# of iterations");ylabel("LWT (s)");  plot(avg_LWT_guided,"--b");legend('unguided','guided');

% Figure Min dist
min_dist_guided=mean(report_guided(:,:,11));
min_dist_unguided=mean(report_unguided(:,:,11));

figure; plot(min_dist_unguided,"r--"); hold on;title("Min Distance between members");xlabel("# of iterations");ylabel("Distance");  plot(min_dist_guided,"--b");legend('unguided','guided');

% Figure Max dist
max_dist_guided=mean(report_guided(:,:,12));
max_dist_unguided=mean(report_unguided(:,:,12));

figure; plot(max_dist_unguided,"r--"); hold on;title("Max Distance between members");xlabel("# of iterations");ylabel("Distance");  plot(max_dist_guided,"--b");legend('unguided','guided');

% Figure AVG dist
avg_dist_guided=mean(report_guided(:,:,13));
avg_dist_unguided=mean(report_unguided(:,:,13));

figure; plot(avg_dist_unguided,"r--"); hold on;title("Average Distance between members");xlabel("# of iterations");ylabel("Distance");plot(avg_dist_guided,"--b");legend('unguided','guided');

%% Comparision of guided vs unguided for the archives

% LWT
lwt_guided=getLWT(archive_guided);
lwt_unguided=getLWT(archive_unguided);

[H, pValue, SWstatistic]=swtest(lwt_guided,0.05);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(lwt_unguided,0.05);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
xx=[lwt_guided(:,1),lwt_unguided(:,1) ];
p_lwt=anova1(xx);

A12_lwt = Atest(lwt_guided(:,1),lwt_unguided(:,1));
d_lwt=computeCohen_d(lwt_guided(:,1),lwt_unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

figure
boxplot([lwt_guided(:,1),lwt_unguided(:,1)],'Notch','on','Labels',{'Guided','Unguided'})
title('LWT difference guided vs unguided')

%% Dist to member
dist_guided=reshape(archive_guided(:,:,5),[],1);
dist_unguided=reshape(archive_unguided(:,:,5),[],1);


[H, pValue, SWstatistic]=swtest(dist_guided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(dist_unguided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
xx=[dist_guided(:,1),dist_unguided(:,1) ];
p_dist=anova1(xx);

A12_dist = Atest(dist_guided(:,1),dist_unguided(:,1));
d_dist=computeCohen_d(dist_guided(:,1),dist_unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(1-A12_dist)+" with a significance value of "+num2str(p_dist)+" and cohen values of: "+num2str(d_dist))

figure
boxplot([dist_guided(:,1),dist_unguided(:,1)],'Notch','on','Labels',{'Guided','Unguided'})

title('Distance to member guided vs unguided')


%% ScatterPlot of # of different parameters vs LWT
[lwt_guided, changes_guided]=getDataArchive(archive_guided);
[lwt_unguided, changes_unguided]=getDataArchive(archive_unguided);

figure;
hold on
plot(changes_guided, lwt_guided,'sb','MarkerFaceColor','blue', 'MarkerSize', 7)
plot(changes_unguided,lwt_unguided, '.r', 'MarkerSize', 20)
title("LWT difference and hamming distance between members")
xlabel('Hamming distance', 'FontWeight','bold') 
ylabel('LWT difference between pairs (s)', 'FontWeight','bold') 
legend({'guided','unguided'});
%xlim([0,25]);
xlim([0,0.005]);
ylim([0,350]);

%% # of different parameters
changes_guided=reshape(archive_guided(:,:,6),[],1);
changes_unguided=reshape(archive_unguided(:,:,6),[],1);

[H, pValue, SWstatistic]=swtest(changes_guided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(changes_unguided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
xx=[changes_guided(:,1),changes_unguided(:,1) ];
[p_changes,tbl,stats] = kruskalwallis(xx,[],'off');

A12_changes = Atest(changes_guided(:,1),changes_unguided(:,1));
d_changes=computeCohen_d(changes_guided(:,1),changes_unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(1-A12_changes)+" with a significance value of "+num2str(p_changes)+" and cohen values of: "+num2str(d_changes))
figure
boxplot([changes_guided(:,1),changes_unguided(:,1)],'Notch','on','Labels',{'Guided','Unguided'})
ylim([0 30]);
title('# of changes guided vs unguided')
%%
% AUC (area under the curve)
changes_guided=[10477.1375, 10115.474999999999, 9852.285, 9580.210000000001, 10350.642499999998, 9562.9075, 11332.915, 10642.245000000003]';
changes_unguided=[5539.451249999999, 7879.374374999999, 3996.2875000000013, 5515.171250000001, 5952.593125, 6559.913124999999, 5789.652500000001, 6325.2112499999985]';

[H, pValue, SWstatistic]=swtest(changes_guided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(changes_unguided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
xx=[changes_guided(:,1);changes_unguided(:,1) ];
[p_changes,tbl,stats] = kruskalwallis(xx,[],'off');

A12_changes = Atest(changes_guided(:,1),changes_unguided(:,1));
d_changes=computeCohen_d(changes_guided(:,1),changes_unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_changes)+" with a significance value of "+num2str(p_changes)+" and cohen values of: "+num2str(d_changes))
figure
boxplot([changes_guided(:,1),changes_unguided(:,1)],'Notch','on','Labels',{'Guided','Unguided'})
title('Area Under the Curve')

%% Functions

function lwt=getLWT(archive)

    for i=1:size(archive,1)
        for j=1:size(archive,2)/2
            index1=j*2-1;
            index2=j*2;
            lwt_matrices(i,j)=abs(archive(i,index1,3)-archive(i,index2,3));
            
        end
        lwt_matrices2(i)=mean(lwt_matrices(i,:));
    end

    lwt=reshape(lwt_matrices2,[],1);
end

function [lwt,diff_params]=getDataArchive(archive)

    for i=1:size(archive,1)
        for j=1:size(archive,2)/2
            index1=j*2-1;
            index2=j*2;
            lwt_matrices(i,j)=abs(archive(i,index1,3)-archive(i,index2,3));
            diff_params(i,j)=archive(i,index1,5);
        end
    end

    lwt=reshape(lwt_matrices,[],1);
    diff_params=reshape(diff_params,[],1);
end

function [cohen_d,d]=getCohenValue(A12)


    d=2*abs(A12-0.5);
    if d<0.147
        cohen_d="negligible";
    elseif d<0.33 
        cohen_d="small";
    elseif d<0.474
        cohen_d="medium";
    elseif d>=0.474
        cohen_d="large";
    end
    

end
