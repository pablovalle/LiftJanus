import pandas as pd
import numpy as np

def time_to_hours(time_str):
    # Split the time string into hours, minutes, and seconds
    hours, minutes, seconds = map(float, time_str.split(':'))
    
    # Calculate the fraction of hours for minutes and seconds
    minutes_fraction = minutes / 60
    seconds_fraction = seconds / 3600
    
    # Calculate the total time in hours
    total_hours = hours + minutes_fraction + seconds_fraction
    
    return total_hours

def time_to_hours_withoutDD(time_str):
    # Split the time string into hours, minutes, and seconds
    days, hours, minutes, seconds = map(float, time_str.split(':'))
    
    # Calculate the fraction of hours for minutes and seconds
    days_fraction=days * 24
    minutes_fraction = minutes / 60
    seconds_fraction = seconds / 3600
    
    # Calculate the total time in hours
    total_hours = days_fraction + hours + minutes_fraction + seconds_fraction
    
    return total_hours


path="..\Experiments_withoutDD\ExperimentsScenario1\guided\\"
times=[]
for i in range(2,6):
    to_read=path+str(i)+"\\report"+str(i)+".csv"
    data=pd.read_csv(to_read)
    time_before=data.timestamp[len(data)-1]
    #Uncomment this for analyzing without DD ---------------
    time_before=time_before.replace(" days, ",":")
    time_after=time_to_hours_withoutDD(time_before)
    #--------------------------------------------------------
    #Uncomment this for analyzing with DD ---------------
    #time_after=time_to_hours(time_before)
    #--------------------------------------------------------
    times.append(time_after)

print(np.mean(times))
print(np.std(times))