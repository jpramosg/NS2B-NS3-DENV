import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
import seaborn as sns 
 
## indicate path 

path = r"/home/jpramosg/Desktop/MD/dengue/data/APO/results_analysis/"

## Load RMSD data from APO and HOLO state of NS3 monomer 
######### NS3 ################
#APO 
filename = r"RMSD_NS3_r1.dat" 
df_NS3_r1 = pd.read_csv(path+filename, sep="\t")
filename = r"RMSD_NS3_r2.dat" 
df_NS3_r2 = pd.read_csv(path+filename, sep="\t")
filename = r"RMSD_NS3_r3.dat" 
df_NS3_r3 = pd.read_csv(path+filename, sep="\t")

#HOLO 

path = r"/home/jpramosg/Desktop/MD/dengue/data/HOLO/results_analysis/"

filename = r"RMSDHo_NS3_r1.dat" 
df_NS3Ho_r1 = pd.read_csv(path+filename, sep="\t")
filename = r"RMSDHo_NS3_r2.dat" 
df_NS3Ho_r2 = pd.read_csv(path+filename, sep="\t")
filename = r"RMSDHo_NS3_r3.dat" 
df_NS3Ho_r3 = pd.read_csv(path+filename, sep="\t")

# get data RMSD and time of replicas 

x_NS3 =  df_NS3_r1["time"].dropna() 
y1_NS3 = df_NS3_r1["rmsd"].dropna()
y2_NS3 = df_NS3_r2["rmsd"].dropna()
y3_NS3 = df_NS3_r3["rmsd"].dropna()

x_NS3Ho =  df_NS3Ho_r1["time"].dropna()
y1_NS3Ho = df_NS3Ho_r1["rmsd"].dropna()
y2_NS3Ho = df_NS3Ho_r2["rmsd"].dropna()
y3_NS3Ho = df_NS3Ho_r3["rmsd"].dropna()

# get average from all replicas 
## Ensure all y arrays have the same length

num_replica = 3 

min_length = min(len(y1_NS3), len(y2_NS3), len(y3_NS3))

y1_NS3 = y1_NS3[:min_length]
y2_NS3 = y2_NS3[:min_length]
y3_NS3 = y3_NS3[:min_length]
x_NS3 = x_NS3[:min_length]

y_avg_NS3 = (y1_NS3 + y2_NS3 + y3_NS3) / num_replica

min_length = min(len(y1_NS3Ho), len(y2_NS3Ho), len(y3_NS3Ho))

y1_NS3Ho = y1_NS3Ho[:min_length]
y2_NS3Ho = y2_NS3Ho[:min_length]
y3_NS3Ho = y3_NS3Ho[:min_length]
x_NS3Ho = x_NS3Ho[:min_length]

y_avg_NS3Ho = (y1_NS3Ho + y2_NS3Ho + y3_NS3Ho) / num_replica

# bodacius color 
## 4 cause i have only three data set 
colors = sns.color_palette("rocket", 3)
seshadri = ["#c3121e", "#0348a1", "#ffb01c", "#027608", "#0193b0"]


fig = plt.figure(1, figsize=(10, 5))
plt.plot(x_NS3, y_avg_NS3, linestyle ="-", label = "NS3 Apo ", color = colors[0], mfc = "w", linewidth=1 )
plt.plot(x_NS3Ho, y_avg_NS3Ho, linestyle ="-", label = "NS3 Holo ", color = colors[1], mfc = "w", linewidth=1 )
#plt.plot(x1, y3, linestyle ="-", label = "NS2B R3", color = colors[2], mfc = "w", markersize=8 )

#plot params
plt.xlim([2,300]) ## range of values in x-axis 
plt.ylim([0,4]) ## range of values in y-axis 
plt.minorticks_on() ## between ranges 0 to 5, show little lines 
plt.tick_params(labelsize=14) ## size of numbers in axis 
xticks = np.arange(0,300.1,100) ## create a matrix with values 

plt.tick_params(direction="in", which="minor", length=4.5) ## size of lines minor in axis 
plt.tick_params(direction="in", which="major", length=9) ## size of lines major in axis 
plt.xticks(xticks)

plt.xlabel(r"Time [ns]", fontsize = 14)
plt.ylabel(r"RMSD [$\mathrm{\AA}$]", fontsize = 14)

plt.legend(fontsize=14)


plt.savefig("/home/jpramosg/Desktop/MD/dengue/plots/RMSD_NS3.png", dpi = 300, bbox_inches="tight")
