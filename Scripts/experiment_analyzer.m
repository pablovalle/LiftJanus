clc; clear; close all;
initialPath="..\Experiments_withDD\ExperimentsInstallation2";
addpath(initialPath+"\guided\");
addpath(initialPath+"\unguided\");

numOfExecs=9;
numOfIterations=99;
baselineRuns=7;

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


for i=0:baselineRuns
    baselinePath="..\Baseline\ExperimentsInstallation2\Run_"+num2str(i)+"\";
    archive_baseline_mass{i+1}=table2array(importBaseline(baselinePath+"Pareto_mass.csv"));
    archive_baseline_capacity{i+1}=table2array(importBaseline(baselinePath+"Pareto_capacity.csv"));
    archive_baseline_mass_capacity{i+1}=table2array(importBaseline(baselinePath+"Pareto_mass_capacity.csv"));

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
lwt_baseline_mass=getbaselineLWT(archive_baseline_mass);
lwt_baseline_capacity=getbaselineLWT(archive_baseline_capacity);
lwt_baseline_mass_capacity=getbaselineLWT(archive_baseline_mass_capacity);

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
[H, pValue, SWstatistic]=swtest(lwt_baseline_mass,0.05);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(lwt_baseline_capacity,0.05);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(lwt_baseline_mass_capacity,0.05);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
disp("Approach VS Unguided")
xx=[lwt_guided(:,1),lwt_unguided(:,1)];
p_lwt=anova1(xx);
A12_lwt = Atest(lwt_guided(:,1),lwt_unguided(:,1));
d_lwt=computeCohen_d(lwt_guided(:,1),lwt_unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

disp("Approach VS Baseline C1")
xx=[lwt_guided(:,1);lwt_baseline_mass(:,1)];
p_lwt = anova1(xx);
A12_lwt = Atest(lwt_guided(:,1),lwt_baseline_mass(:,1));
d_lwt=computeCohen_d(lwt_guided(:,1),lwt_baseline_mass(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

disp("Approach VS Baseline C2")
xx=[lwt_guided(:,1),lwt_baseline_capacity(:,1)];
p_lwt=anova1(xx);
A12_lwt = Atest(lwt_guided(:,1),lwt_baseline_capacity(:,1));
d_lwt=computeCohen_d(lwt_guided(:,1),lwt_baseline_capacity(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

disp("Approach VS Baseline C3")
xx=[lwt_guided(:,1),lwt_baseline_mass_capacity(:,1)];
p_lwt=anova1(xx);
A12_lwt = Atest(lwt_guided(:,1),lwt_baseline_mass_capacity(:,1));
d_lwt=computeCohen_d(lwt_guided(:,1),lwt_baseline_mass_capacity(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

%% Dist to member
dist_guided=getHD(archive_guided);
dist_unguided=getHD(archive_unguided);
dist_baseline_mass=getbaselineHD(archive_baseline_mass);
dist_baseline_capacity=getbaselineHD(archive_baseline_capacity);
dist_baseline_mass_capacity=getbaselineHD(archive_baseline_mass_capacity);

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
[H, pValue, SWstatistic]=swtest(dist_baseline_mass);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(dist_baseline_capacity);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(dist_baseline_mass_capacity);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end

disp("Approach VS Unguided")
xx=[dist_guided(:,1),dist_unguided(:,1)];
p_hd = anova1(xx, [],'off');
A12_hd = Atest(dist_guided(:,1),dist_unguided(:,1));
d_hd=computeCohen_d(dist_guided(:,1),dist_unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_hd)+" with a significance value of "+num2str(p_hd)+" and cohen values of: "+num2str(d_hd))

disp("Approach VS Baseline C1")
xx=[dist_guided(:,1),dist_baseline_mass(:,1)];
p_hd = anova1(xx, [],'off');
A12_hd = Atest(dist_guided(:,1),dist_baseline_mass(:,1));
d_hd=computeCohen_d(dist_guided(:,1),dist_baseline_mass(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_hd)+" with a significance value of "+num2str(p_hd)+" and cohen values of: "+num2str(d_hd))

disp("Approach VS Baseline C2")
xx=[dist_guided(:,1),dist_baseline_capacity(:,1)];
p_hd = anova1(xx, [],'off');
A12_hd = Atest(dist_guided(:,1),dist_baseline_capacity(:,1));
d_hd=computeCohen_d(dist_guided(:,1),dist_baseline_capacity(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_hd)+" with a significance value of "+num2str(p_hd)+" and cohen values of: "+num2str(d_hd))

disp("Approach VS Baseline C3")
xx=[dist_guided(:,1),dist_baseline_mass_capacity(:,1)];
p_hd = anova1(xx, [],'off');
A12_hd = Atest(dist_guided(:,1),dist_baseline_mass_capacity(:,1));
d_hd=computeCohen_d(dist_guided(:,1),dist_baseline_mass_capacity(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_hd)+" with a significance value of "+num2str(p_hd)+" and cohen values of: "+num2str(d_hd))


%% ScatterPlot of # of different parameters vs LWT
[lwt_guided, changes_guided]=getDataArchive(archive_guided);
[lwt_unguided, changes_unguided]=getDataArchive(archive_unguided);
[lwt_baseline_mass, changes_baseline_mass]=getDataArchiveBaseline(archive_baseline_mass);
[lwt_baseline_capacity, changes_baseline_capacity]=getDataArchiveBaseline(archive_baseline_capacity);
[lwt_baseline_mass_capacity, changes_baseline_mass_capacity]=getDataArchiveBaseline(archive_baseline_mass_capacity);


figure;
hold on
plot(changes_guided, lwt_guided,'sb','MarkerFaceColor','blue', 'MarkerSize', 7)
plot(changes_unguided,lwt_unguided, '.r', 'MarkerSize', 20)
plot(changes_baseline_mass,lwt_baseline_mass, '*', 'MarkerSize', 7)
plot(changes_baseline_capacity,lwt_baseline_capacity, '+', 'MarkerSize', 7)
plot(changes_baseline_mass_capacity,lwt_baseline_mass_capacity, '^', 'MarkerSize', 7, 'MarkerFaceColor',[0.4660 0.6740 0.1880])

title("LWT difference and hamming distance between members")
xlabel('Hamming distance', 'FontWeight','bold') 
ylabel('LWT difference between pairs (s)', 'FontWeight','bold') 
legend({'LiftJanus','BaseJnaus', 'UncerRobua - C1', 'UncerRobua - C2', 'UncerRobua - C3'}, "FontSize",7); 

%xlim([0,25]);
%xlim([0,0.005]);
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
%Installation 1
guided=[10477.1375, 10115.474999999999, 9852.285, 9580.210000000001, 10350.642499999998, 9562.9075, 11332.915, 10642.245000000003]';
unguided=[5539.451249999999, 7879.374374999999, 3996.2875000000013, 5515.171250000001, 5952.593125, 6559.913124999999, 5789.652500000001, 6325.2112499999985]';
mass=[0 0 0 0 0 0 0 0 ]';
capacity=[7959.42 7863.00125 8807.29875 8650.31625 7936.8125 9437.98375 9170.22875 8592.37875]';
mass_capacity=[10010.0675 7944.41375 8583.58 9117.70875 8397.07875 8280.65875 7870.01125 8883.9075 ]';
%Installation 2
% guided=[10379.2375  9187.78   11437.765  10469.6625  9303.1825 11028.435 10929.0675 11606.6775]';
% unguided=[6004.894375 7096.6475   6833.973125 5918.1525   3239.0175   6124.79 5870.27125  7163.315625]';
% mass=[105.08125 0. 329.86 0. 173.075 0. 49.9875 27.4125]';
% capacity= [11374.325 13316.8575 14015.97 12543.3425 12869.87375 14056.73125 12657.2075 13655.725]';
% mass_capacity= [13954.6475 13654.1525 13247.2525 13039.255 12727.63 12421.01625 12411.79875 13647.60375]';

[H, pValue, SWstatistic]=swtest(guided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(unguided);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(mass);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(capacity);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end
[H, pValue, SWstatistic]=swtest(mass_capacity);
if H==1 
    disp("Not normally distributed")
else
    disp("Normally distributed")
end

disp("Approach VS Unguided")
xx=[guided(:,1),unguided(:,1)];
p_lwt=anova1(xx);
A12_lwt = Atest(guided(:,1),unguided(:,1));
d_lwt=computeCohen_d(guided(:,1),unguided(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

disp("Approach VS Baseline C1")
xx=[guided(:,1),mass(:,1)];
groupLabels = [ones(size(guided(:,1))); 2*ones(size(mass(:,1)))];
p_lwt = kruskalwallis(xx, [],'off');
A12_lwt = Atest(guided(:,1),mass(:,1));
d_lwt=computeCohen_d(guided(:,1),mass(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

disp("Approach VS Baseline C2")
xx=[guided(:,1),capacity(:,1)];
groupLabels = [ones(size(guided(:,1))); 2*ones(size(mass(:,1)))];
p_lwt = kruskalwallis(xx, [],'off');
A12_lwt = Atest(guided(:,1),capacity(:,1));
d_lwt=computeCohen_d(guided(:,1),capacity(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))

disp("Approach VS Baseline C3")
xx=[guided(:,1),mass_capacity(:,1)];
p_lwt=anova1(xx);
A12_lwt = Atest(guided(:,1),mass_capacity(:,1));
d_lwt=computeCohen_d(guided(:,1),mass_capacity(:,1));
disp("The guided version dominates with an A12 value of "+num2str(A12_lwt)+" with a significance value of "+num2str(p_lwt)+" and cohen values of: "+num2str(d_lwt))



% xx=[changes_guided(:,1);changes_unguided(:,1) ];
% [p_changes,tbl,stats] = kruskalwallis(xx,[],'off');
% 
% A12_changes = Atest(changes_guided(:,1),changes_unguided(:,1));
% d_changes=computeCohen_d(changes_guided(:,1),changes_unguided(:,1));
% disp("The guided version dominates with an A12 value of "+num2str(A12_changes)+" with a significance value of "+num2str(p_changes)+" and cohen values of: "+num2str(d_changes))
% figure
% boxplot([changes_guided(:,1),changes_unguided(:,1)],'Notch','on','Labels',{'Guided','Unguided'})
% title('Area Under the Curve')

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
function lwt=getbaselineLWT(archive)

    for i=1:size(archive,2)
        for j=1:size(archive{i},1)
            lwt_matrices(i,j)=archive{i}(j,3);
        end
        lwt_matrices2(i)=mean(lwt_matrices(i,:));
    end
 
    lwt=reshape(lwt_matrices2,[],1);
end

function hamming_dist=getbaselineHD(archive)

    for i=1:size(archive,2)
        for j=1:size(archive{i},1)
            hd_matrices(i,j)=archive{i}(j,4);
        end
        hd_matrices2(i)=mean(hd_matrices(i,:));
    end
 
    hamming_dist=reshape(hd_matrices2,[],1);
end

function hamming_dist=getHD(archive)

    for i=1:size(archive,1)
        for j=1:size(archive,2)
            hd_matrices(i,j)=archive(i,j,5);
        end
        hd_matrices2(i)=mean(hd_matrices(i,:));
    end
 
    hamming_dist=reshape(hd_matrices2,[],1);
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

function [lwt,diff_params]=getDataArchiveBaseline(archive)
    count=1;
    for i=1:size(archive,2)
        for j=1:size(archive{i},1)
            lwt_matrices(count)=archive{i}(j,3);
            diff_params(count)=archive{i}(j,4);
            count=count+1;
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
