import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
import seaborn as sns 
from matplotlib.gridspec import GridSpec
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


fig = plt.figure(figsize=(10, 5))
gs = GridSpec(4, 4, hspace=0.4, wspace=0.0)  # Reduced wspace to bring the density plot closer

ax_main = fig.add_subplot(gs[:, :-1])
ax_density = fig.add_subplot(gs[:, -1], sharey=ax_main)

# Main plot
ax_main.plot(x_NS3, y_avg_NS3, linestyle ="-", label = "NS3 Apo", color = colors[0], linewidth=1)
ax_main.plot(x_NS3Ho, y_avg_NS3Ho, linestyle ="-", label = "NS3 Holo", color = colors[2], linewidth=1)

ax_main.set_xlim([2, 300])
ax_main.set_ylim([0, 4])
ax_main.minorticks_on()
ax_main.tick_params(labelsize=14)
xticks = np.arange(0, 300.1, 100)
yticks = np.arange(0, 6.1, 1)
ax_main.set_xticks(xticks)
ax_main.set_yticks(yticks)


ax_main.tick_params(direction="in", which="minor", length=3.0)
ax_main.tick_params(direction="in", which="major", length=8)

ax_main.set_xlabel(r"Time [ns]", fontsize=14,  weight='bold') 
ax_main.set_ylabel(r"RMSD [$\mathrm{\AA}$]", fontsize=14,  weight='bold')

legend = ax_main.legend(fontsize=14, loc="upper left", handlelength=1, frameon=False )

for legend_handle in legend.get_lines():
    legend_handle.set_linewidth(9.0)

sns.kdeplot(y=y_avg_NS3, ax=ax_density, color = colors[0], fill=True)
sns.kdeplot(y=y_avg_NS3Ho, ax=ax_density, color = colors[2], fill=True)

#sns.kdeplot(y_avg_NS3, ax=ax_density, color="red", vertical=True, fill=True)
#sns.kdeplot(y_avg_NS3Ho, ax=ax_density, color="blue", vertical=True, fill=True)

# Remove ticks and borders from the density plot
ax_density.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
ax_density.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_density.set_xlabel('')
ax_density.set_ylabel('')
#ax_density.set_yticks([])


# Remove borders from the density plot
ax_density.spines['top'].set_visible(False)
ax_density.spines['right'].set_visible(False)
ax_density.spines['bottom'].set_visible(False)
ax_density.spines['left'].set_visible(True)  # Ensure no y-ticks on the density plot

plt.savefig("/home/jpramosg/Desktop/MD/dengue/plots/RMSD_NS3.png", dpi = 300, bbox_inches="tight")
