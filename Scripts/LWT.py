import pandas as pd
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
import numpy as np
import matplotlib.pyplot as plt

runs=[2,3,4,5,6,7,8,9]
runs_withoutDD=[2,3,4,5]
guided=list()
unguided=list()
guided_withoutDD=list()
LWT_unguided=list()
LWT_guided=list()
LWT_guided_withoutDD=list()
mass=list()
capacity=list()
mass_capacity=list()
installation="ExperimentsInstallation2"

for run in runs:
    file_guided=pd.read_csv(os.path.join("../Experiments_withDD",installation,"guided",str(run),"report"+str(run)+".csv"))
    guided.append(file_guided)
    LWT_guided.append(np.asarray(file_guided['LWT_avg']))


    file_unguided=pd.read_csv(os.path.join("../Experiments_withDD",installation,"unguided",str(run),"report"+str(run)+".csv"))
    unguided.append(file_unguided)
    LWT_unguided.append(np.asarray(file_unguided['LWT_avg']))

for run in runs_withoutDD:
    file_guided_withoutDD=pd.read_csv(os.path.join("../Experiments_withoutDD",installation,"guided",str(run),"report"+str(run)+".csv"))
    guided_withoutDD.append(file_guided_withoutDD)
    LWT_guided_withoutDD.append(np.asarray(file_guided_withoutDD['LWT_avg']))

for i in range(0, len(runs)):
    mass_file=pd.read_csv(os.path.join("../Baseline", installation, "Run_"+str(i), "Pareto_iteration_mass.csv"))
    mass.append(np.asarray(mass_file['LWT']))
    capacity_file=pd.read_csv(os.path.join("../Baseline", installation, "Run_"+str(i), "Pareto_iteration_capacity.csv"))
    capacity.append(np.asarray(capacity_file['LWT']))
    mass_capacity_file=pd.read_csv(os.path.join("../Baseline", installation, "Run_"+str(i), "Pareto_iteration_mass_capacity.csv"))
    mass_capacity.append(np.asarray(mass_capacity_file['LWT']))

mean_guided=list()
std_guided=list()

mean_mass=list()
std_mass=list()

mean_capacity=list()
std_capacity=list()

mean_mass_capacity=list()
std_mass_capacity=list()

mean_unguided=list()
std_unguided=list()

mean_unguided_withoutDD=list()
std_unguided_withoutDD=list()


auc_guided=list()
auc_unguided=list()
auc_mass=list()
auc_capacity=list()
auc_mass_capacity=list()


time_intervals=[i for i in range(1,100)]

for i in range(0,len(LWT_guided[0])):
    guided_values=list()
    unguided_values=list()
    unguided_values_withoutDD=list()
    mass_values=list()
    capacity_values=list()
    mass_capacity_values=list()
    for j in range(0,len(LWT_guided)):
        guided_values.append(LWT_guided[j][i])
        unguided_values.append(LWT_unguided[j][i])
        mass_values.append(mass[j][i])
        capacity_values.append(capacity[j][i])
        mass_capacity_values.append(mass_capacity[j][i])
    for j in range(0,len(LWT_guided_withoutDD)):
        unguided_values_withoutDD.append(LWT_guided_withoutDD[j][i])
    mean_guided.append(np.mean(guided_values))
    mean_unguided.append(np.mean(unguided_values))
    mean_mass.append(np.mean(mass_values))
    mean_capacity.append(np.mean(capacity_values))
    mean_mass_capacity.append(np.mean(mass_capacity_values))
    mean_unguided_withoutDD.append(np.mean(unguided_values_withoutDD))
    std_guided.append(np.std(guided_values))
    std_unguided.append(np.std(unguided_values))
    std_unguided_withoutDD.append(np.std(unguided_values_withoutDD))
    std_mass.append(np.std(mass_values))
    std_capacity.append(np.std(capacity_values))
    std_mass_capacity.append(np.std(mass_capacity_values))
for i in range(0,len(LWT_guided)):
    auc_guided.append(np.trapezoid(y=LWT_guided[i], x=time_intervals))
    auc_unguided.append(np.trapezoid(y=LWT_unguided[i], x=time_intervals))
    auc_mass.append(np.trapezoid(y=mass[i][:-1], x=time_intervals))
    auc_capacity.append(np.trapezoid(y=capacity[i][:-1], x=time_intervals))
    auc_mass_capacity.append(np.trapezoid(y=mass_capacity[i][:-1], x=time_intervals))

print(np.asarray(auc_guided))
print(np.asarray(auc_unguided))
print(np.asarray(auc_mass))
print(np.asarray(auc_capacity))
print(np.asarray(auc_mass_capacity))
print(mean_guided[-1])
print(mean_unguided[-1])
print(mean_unguided_withoutDD[-1])

#Lipings baseline Installation 1
#mass=[0, 0, 0, 0, 0, 46.1000, 0, 0]
#capacity=[123.0200, 104.7000, 48.6714, 114.7250, 78.1250, 68.2625, 116.5500, 65.5500]
#mass_capacity=[166.3250, 126.6000, 91.9500, 132.4400, 124.6333, 122.1333, 78.6167, 105.5700]

#Lipings baseline Installation 2
#mass=[19.9750, 11.8000, 38.1250, 11.9750, 11.8000, 13.5500, 11.8000, 14.6000]
#capacity=[164.2286, 89.8000, 79.9143, 128.5857, 144.4875, 150.4444, 101.5889, 75.7444]
#mass_capacity=[168.3300, 82.1100, 53.8100, 133.5300, 58.4100, 79.3600, 152.3100, 63.9300]


plt.plot(time_intervals, mean_guided, label='LiftJanus', color='#1f77b4')
plt.fill_between(time_intervals, np.asanyarray(mean_guided)-np.asanyarray(std_guided), np.asanyarray(mean_guided)+np.asanyarray(std_guided), color='#1f77b4', alpha=0.3)

plt.plot(time_intervals, mean_unguided,'--', label='BaseJanus', color='red')
plt.fill_between(time_intervals, np.asanyarray(mean_unguided)-np.asanyarray(std_unguided), np.asanyarray(mean_unguided)+np.asanyarray(std_unguided), color='red', alpha=0.3)


plt.ylabel("Average LWT difference (s)",weight='bold')
plt.xlabel("# of iteration",weight='bold')
plt.title("Average Longest Waiting Time difference between members",weight='bold')
plt.xlim([0,100])
plt.ylim([0,200])
plt.plot(time_intervals, mean_unguided_withoutDD, '.',label='LiftJanus non-minimized', color='green')
plt.fill_between(time_intervals, np.asanyarray(mean_unguided_withoutDD)-np.asanyarray(std_unguided_withoutDD), np.asanyarray(mean_unguided_withoutDD)+np.asanyarray(std_unguided_withoutDD), color='green', alpha=0.3)

plt.plot(time_intervals, mean_mass, 'v',label='UncerRobua C1' , markersize=4, color='#e377c2')
plt.fill_between(time_intervals, np.asanyarray(mean_mass)-np.asanyarray(std_mass), np.asanyarray(mean_mass)+np.asanyarray(std_mass), color='#e377c2', alpha=0.1)

plt.plot(time_intervals, mean_capacity, 's',label='UncerRobua C2', markersize=3, color='purple')
plt.fill_between(time_intervals, np.asanyarray(mean_capacity)-np.asanyarray(std_capacity), np.asanyarray(mean_capacity)+np.asanyarray(std_capacity),color='purple', alpha=0.1)

plt.plot(time_intervals, mean_mass_capacity, '*',label='UncerRobua C3', markersize=4,color='orange')
plt.fill_between(time_intervals, np.asanyarray(mean_mass_capacity)-np.asanyarray(std_mass_capacity), np.asanyarray(mean_mass_capacity)+np.asanyarray(std_mass_capacity),color='orange', alpha=0.1)

plt.legend(fontsize="10",loc ="best")
plt.savefig("RQ1_Scenario2_new.pdf", format="pdf", bbox_inches="tight")
plt.show()


