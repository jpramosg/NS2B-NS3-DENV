import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
import seaborn as sns 

### NS3 APO ### 

path = r"/home/jpramosg/Desktop/MD/dengue/data/APO/results_analysis/"

filename = r"RMSF_NS3_r1.dat"
RMSF_NS3_r1 = pd.read_csv(path+filename, sep="\t")
## load data of rmsd from replica 2 
filename = r"RMSF_NS3_r2.dat"
RMSF_NS3_r2 = pd.read_csv(path+filename, sep="\t")
## load data of rmsd from replica 3 
filename = r"RMSF_NS3_r3.dat"
RMSF_NS3_r3 = pd.read_csv(path+filename, sep="\t")

path = r"/home/jpramosg/Desktop/MD/dengue/data/HOLO/results_analysis/"

filename = r"RMSFHo_NS3_r1.dat"
RMSF_Ho_NS3_r1 = pd.read_csv(path+filename, sep="\t")
## load data of rmsd from replica 2 
filename = r"RMSFHo_NS3_r2.dat"
RMSF_Ho_NS3_r2 = pd.read_csv(path+filename, sep="\t")
## load data of rmsd from replica 3 
filename = r"RMSFHo_NS3_r3.dat"
RMSF_Ho_NS3_r3 = pd.read_csv(path+filename, sep="\t")

## Create date and getting average 
x1_RMSF_NS3 = RMSF_NS3_r1["resid"].dropna()
y1_RMSF_NS3 = RMSF_NS3_r1["rmsf"].dropna()
y2_RMSF_NS3 = RMSF_NS3_r2["rmsf"].dropna()
y3_RMSF_NS3 = RMSF_NS3_r3["rmsf"].dropna()


min_length = min(len(y1_RMSF_NS3), len(y2_RMSF_NS3), len(y3_RMSF_NS3))

x1_RMSF_NS3 = x1_RMSF_NS3[:min_length]
y1_RMSF_NS3 = y1_RMSF_NS3[:min_length]
y2_RMSF_NS3 = y2_RMSF_NS3[:min_length]
y3_RMSF_NS3 = y3_RMSF_NS3[:min_length]

y_avg_RMSF_NS3 = (y1_RMSF_NS3 + y2_RMSF_NS3 + y3_RMSF_NS3) / num_replica

###########

x1_RMSF_Ho_NS3 = RMSF_Ho_NS3_r1["resid"].dropna()
y1_RMSF_Ho_NS3 = RMSF_Ho_NS3_r1["rmsf"].dropna()
y2_RMSF_Ho_NS3 = RMSF_Ho_NS3_r2["rmsf"].dropna()
y3_RMSF_Ho_NS3 = RMSF_Ho_NS3_r3["rmsf"].dropna()


min_length = min(len(y1_RMSF_Ho_NS3), len(y2_RMSF_Ho_NS3), len(y3_RMSF_Ho_NS3))

x1_RMSF_Ho_NS3 = x1_RMSF_Ho_NS3[:min_length]
y1_RMSF_Ho_NS3 = y1_RMSF_Ho_NS3[:min_length]
y2_RMSF_Ho_NS3 = y2_RMSF_Ho_NS3[:min_length]
y3_RMSF_Ho_NS3 = y3_RMSF_Ho_NS3[:min_length]

y_avg_RMSF_Ho_NS3 = (y1_RMSF_Ho_NS3 + y2_RMSF_Ho_NS3 + y3_RMSF_Ho_NS3) / num_replica

# bodacius color 
## 4 cause i have only three data set 
colors = sns.color_palette("rocket", 3)
seshadri = ["#c3121e", "#0348a1", "#ffb01c", "#027608", "#0193b0"]

fig = plt.figure(1, figsize=(20, 5))
plt.plot(x1_RMSF_NS3, y_avg_RMSF_NS3, linestyle ="-", label = "NS3 Apo", color = colors[0], mfc = "w", markersize=8 )
plt.plot(x1_RMSF_Ho_NS3, y_avg_RMSF_Ho_NS3, linestyle ="-", label = "NS3 Holo", color = colors[1], mfc = "w", markersize=8 )

#plot params
plt.xlim([1,154]) ## range of values in x-axis 
#plt.ylim([0,7]) ## range of values in y-axis 
plt.minorticks_on() ## between ranges 0 to 5, show little lines 
plt.tick_params(labelsize=14) ## size of numbers in axis 
xticks = np.arange(1,160.1,40) ## create a matrix with values 

plt.tick_params(direction="in", which="minor", length=4.5) ## size of lines minor in axis 
plt.tick_params(direction="in", which="major", length=9) ## size of lines major in axis 
plt.xticks(xticks)

plt.xlabel(r"Residue index of protein", fontsize = 14)
plt.ylabel(r"RMSF [$\mathrm{\AA}$]", fontsize = 14)

plt.legend(fontsize=14)

plt.savefig("/home/jpramosg/Desktop/MD/dengue/plots/RMSF_NS3.png", dpi = 300, bbox_inches="tight")
