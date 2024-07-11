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

for run in runs:
    file_guided=pd.read_csv(os.path.join("..\Experiments_withDD","ExperimentsScenario2","guided",str(run),"report"+str(run)+".csv"))
    guided.append(file_guided)
    LWT_guided.append(np.asarray(file_guided['LWT_avg']))


    file_unguided=pd.read_csv(os.path.join("..\Experiments_withDD","ExperimentsScenario2","unguided",str(run),"report"+str(run)+".csv"))
    unguided.append(file_unguided)
    LWT_unguided.append(np.asarray(file_unguided['LWT_avg']))

for run in runs_withoutDD:
    file_guided_withoutDD=pd.read_csv(os.path.join("..\Experiments_withoutDD","ExperimentsScenario2","guided",str(run),"report"+str(run)+".csv"))
    guided_withoutDD.append(file_guided_withoutDD)
    LWT_guided_withoutDD.append(np.asarray(file_guided_withoutDD['LWT_avg']))


mean_guided=list()
std_guided=list()

mean_unguided=list()
std_unguided=list()

mean_unguided_withoutDD=list()
std_unguided_withoutDD=list()


auc_guided=list()
auc_unguided=list()


time_intervals=[i for i in range(1,100)]

for i in range(0,len(LWT_guided[0])):
    guided_values=list()
    unguided_values=list()
    unguided_values_withoutDD=list()
    for j in range(0,len(LWT_guided)):
        guided_values.append(LWT_guided[j][i])
        unguided_values.append(LWT_unguided[j][i])
    for j in range(0,len(LWT_guided_withoutDD)):
        unguided_values_withoutDD.append(LWT_guided_withoutDD[j][i])
    mean_guided.append(np.mean(guided_values))
    mean_unguided.append(np.mean(unguided_values))
    mean_unguided_withoutDD.append(np.mean(unguided_values_withoutDD))
    std_guided.append(np.std(guided_values))
    std_unguided.append(np.std(unguided_values))
    std_unguided_withoutDD.append(np.std(unguided_values_withoutDD))
for i in range(0,len(LWT_guided)):
    auc_guided.append(np.trapz(y=LWT_guided[i], x=time_intervals))
    auc_unguided.append(np.trapz(y=LWT_unguided[i], x=time_intervals))

print(auc_guided)
print(auc_unguided)
print(mean_guided[-1])
print(mean_unguided[-1])
print(mean_unguided_withoutDD[-1])





plt.plot(time_intervals, mean_guided, label='Guided')
plt.fill_between(time_intervals, np.asanyarray(mean_guided)-np.asanyarray(std_guided), np.asanyarray(mean_guided)+np.asanyarray(std_guided), alpha=0.3)
plt.ylabel("Average LWT difference (s)",weight='bold')
plt.xlabel("# of iteration",weight='bold')
plt.title("Average Longest Waiting Time difference between members",weight='bold')
plt.ylim([0, None])
plt.xlim([0,100])

plt.plot(time_intervals, mean_unguided,'--', label='Unguided')
plt.fill_between(time_intervals, np.asanyarray(mean_unguided)-np.asanyarray(std_unguided), np.asanyarray(mean_unguided)+np.asanyarray(std_unguided), alpha=0.3)
plt.legend()
plt.plot(time_intervals, mean_unguided_withoutDD, '.',label='Guided non-minimized')
plt.fill_between(time_intervals, np.asanyarray(mean_unguided_withoutDD)-np.asanyarray(std_unguided_withoutDD), np.asanyarray(mean_unguided_withoutDD)+np.asanyarray(std_unguided_withoutDD), alpha=0.3)
plt.legend()
plt.savefig("RQ1_Scenario2_new.pdf", format="pdf", bbox_inches="tight")
plt.show()


