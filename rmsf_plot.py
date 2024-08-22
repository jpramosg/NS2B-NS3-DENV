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

num_replica = 3

y_avg_RMSF_Ho_NS3 = (y1_RMSF_Ho_NS3 + y2_RMSF_Ho_NS3 + y3_RMSF_Ho_NS3) / 3
###########################################################
##### getting of the max value in each region of RMSF #####
###########################################################

df = pd.DataFrame({
    'X': x1_RMSF_NS3,
    'Y': y_avg_RMSF_NS3
})

regions = {
    "R1": [12, 13, 14], 
    "R2": [43, 44, 45], 
    "R3": [69, 70, 71], 
    "R4": [100, 101, 102, 103, 104], 
    "R5": [137, 138, 139, 140, 141, 142, 143]    
}

max_values = []
# Extract rows where X is in values_of_interest

for region_name, R_i in regions.items():
    filt_R_i = df[df['X'].isin(R_i)]
    max_value = filt_R_i['Y'].max()
    max_values.append(max_value)

#################################
#################################

# Set the color for all regions
#region_color = "#FF9999"  # You can change this to any color you prefer

# bodacius color 
## 4 cause i have only three data set 
colors = sns.color_palette("rocket", 3)
seshadri = ["#c3121e", "#0348a1", "#ffb01c", "#027608", "#0193b0"]


fig = plt.figure(1, figsize=(10, 5))
plt.plot(x1_RMSF_NS3, y_avg_RMSF_NS3, linestyle ="-", label = "NS3 Apo", color = colors[0], mfc = "w", markersize=8 )
plt.plot(x1_RMSF_Ho_NS3, y_avg_RMSF_Ho_NS3, linestyle ="-", label = "NS3 Holo", color = colors[2], mfc = "w", markersize=8 )

legend = plt.legend(fontsize=14, loc="upper left", handlelength=1, frameon=False)
for legend_handle in legend.get_lines():
    legend_handle.set_linewidth(9.0)

for i, (region_name, residues) in enumerate(regions.items()):
    plt.axvspan(min(residues), max(residues), color="skyblue", alpha=0.3)
    
    # Get the max value for the current region
    max_value = max_values[i]
    
    # Position the label at the top of the peak
    plt.text((min(residues) + max(residues)) / 2, 
             max_value + 0.4,  # Position the label at the maximum value
             region_name, 
             color="black", 
             fontsize=12, 
             ha='center', 
             va='bottom')

#plot params
plt.xlim([1,154]) ## range of values in x-axis 
#plt.ylim([0,7]) ## range of values in y-axis 
plt.minorticks_on() ## between ranges 0 to 5, show little lines 
plt.tick_params(labelsize=14) ## size of numbers in axis 
xticks = np.arange(0,160.1,20) ## create a matrix with values 

plt.tick_params(direction="in", which="minor", length=4.5) ## size of lines minor in axis 
plt.tick_params(direction="in", which="major", length=9) ## size of lines major in axis 
plt.xticks(xticks)

plt.xlabel(r"Residue index of protein", fontsize = 14, weight='bold')
plt.ylabel(r"RMSF [$\mathrm{\AA}$]", fontsize = 14, weight='bold')

#plt.legend(fontsize=14)

plt.savefig("/home/jpramosg/Desktop/MD/dengue/plots/RMSF_NS3.png", dpi = 300, bbox_inches="tight")
