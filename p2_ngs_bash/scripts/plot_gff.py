
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

frecuencias = [8355, 485, 1, 8568, 8569, 410, 13, 1, 58, 133, 8402, 400, 22, 1]

categorias = ['CDS', 'mat_peptide', 'tmRNA', 'gene', 'exon', 'pseudogene', 
              'misc_feature', 'rep_origin', 'misc_RNA', 'tRNA', 'mRNA', 
              'repeat_region', 'rRNA', 'chromosome']
              

y_pos = np.arange(len(categorias))

# Gráfico de barras
plt.bar(y_pos, frecuencias)

# Nombres en el eje-x
plt.xticks(y_pos, categorias, rotation=90)

# Guardar la grafica
plt.plot()
plt.savefig('figures/barplot.png', bbox_inches='tight')


