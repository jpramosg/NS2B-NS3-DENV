import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors 

path = r"/home/jpramosg/Desktop/MD/dengue/data/HOLO/results_analysis/"

filename = r"FES_r1.dat"
FEL_r1 = pd.read_csv(path+filename, sep='\t', header=None)

filename = r"FES_r2.dat"
FEL_r2 = pd.read_csv(path+filename, sep='\t', header=None)

filename = r"FES_r3.dat"
FEL_r3 = pd.read_csv(path+filename, sep='\t', header=None)


## truncate the range of color with the truncate function 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
    
cmap = plt.get_cmap('viridis')
new_cmap = truncate_colormap(cmap, 0.01, 0.90)

# Asignar nombres a las columnas
FEL_r1.columns = ['RMSD', 'Rgyr', 'FreeEnergy']
FEL_r2.columns = ['RMSD', 'Rgyr', 'FreeEnergy']
FEL_r3.columns = ['RMSD', 'Rgyr', 'FreeEnergy']

# Calcular el promedio de las tres réplicas
average_FEL = (FEL_r1 + FEL_r2 + FEL_r3) / 3


heatmap_data = average_FEL.pivot_table(index='Rgyr', columns='RMSD', values='FreeEnergy') 

#rmsd = FEL['RMSD']
#rgyr = FEL['Rgyr']
#free_energy = FEL['FreeEnergy']

X, Y = np.meshgrid(heatmap_data.columns, heatmap_data.index)
Z = heatmap_data.values 

plt.figure(figsize=(8,6) )

### create a matrix 10 values from min to max 
contour_levels = np.linspace(Z.min(), Z.max(), 100)

### choose size of line depends values in Z
linewidths = [0.0 if level < 2.5 else 0.00 if level < 4.0 else 0.0 for level in contour_levels]

### add the contour lines with differents widths and lines styles 
plt.contour(X, Y, Z, levels=contour_levels, colors='black', linewidths=linewidths, linestyles='-')

plt.contourf(X, Y, Z, levels=contour_levels, cmap=new_cmap)

cbar = plt.colorbar() 

cbar.set_label(r"$\Delta G$ [kcal/mol]", fontsize=14, family='serif')

#plt.xlim([-0.1,3.5]) ## range of values in x-axis 
#plt.ylim([0,10]) ## range of values in y-axis 
plt.minorticks_on() ## between ranges 0 to 5, show little lines 
plt.tick_params(labelsize=14) ## size of numbers in axis 
xticks = np.arange(0.5,3.5,0.5) ## create a matrix with values 

plt.tick_params(direction="in", which="minor", length=4.5) ## size of lines minor in axis 
plt.tick_params(direction="in", which="major", length=9) ## size of lines major in axis 
plt.xticks(xticks)

plt.xlabel(r"RMSD [$\mathrm{\AA}$]", fontsize = 14)
plt.ylabel(r"Rgyr [$\mathrm{\AA}$]", fontsize = 14)
#plt.title('Gráfico de Curvas de Contorno')
print (linewidths)
# Mostrar el gráfico
#plt.show()

plt.savefig("/home/jpramosg/Desktop/MD/dengue/plots/FES_LIG_dengue.png", dpi = 300, bbox_inches="tight") 
