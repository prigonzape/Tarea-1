#!/usr/bin/env python
# coding: utf-8

# # Primera tarea
# Célula objetivo: células CHO
# 
# Modelo: Mus musculus iMM1415

# # 1. FBA modelo iMM1415 predeterminado

# In[73]:


#Importar modelo a escala genómica
import cobra
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from cobra.flux_analysis import production_envelope
from numpy import zeros
model = cobra.io.read_sbml_model("iMM1415.xml.gz")


# In[74]:


print "Distribución de flujos mediante FBA"
print "===================================="
solution=model.optimize()
model.summary()


# In[75]:


print "Características del modelo" #(cantidad de reacciones, metabolitos y genes)
print "==========================="
print "               "

print "Reacciones"
print "=============="
print len(model.reactions)
print "               "

print "Metabolitos"
print "=============="
print len(model.metabolites)
print "               "

print "Genes involucrados"
print "====================="
print len(model.genes)
print "               "

print "Función obetivo del modelo"
print "====================="
print model.objective
print "               "

print "Valor de la solución objetivo"
print "================================"
print solution.objective_value 
print "               "

print "Reacción de biomasa"
print "================================"
biomass = model.reactions.get_by_id("BIOMASS_mm_1_no_glygln")
print biomass
print "=================================================================="


# In[76]:


print "Compartimentalización de las reacciones"
import re
def buscarMetabolito(nombreMetabolite):
    metabolites=[]
    for metabolite in model.metabolites:
        if re.match(nombreMetabolite,metabolite.id,re.IGNORECASE):    
            metabolites.append(metabolite) 
    if len(metabolites)==0:
        print "Not found metabolite"
        return None
    else: 
        for metabolite in metabolites:
            print "=========================================="
            print metabolite.name, metabolite.id
            print "=================================="
            for reaction in model.metabolites.get_by_id(metabolite.id).reactions:
                print reaction

result=buscarMetabolito("gln__L")
result=buscarMetabolito("glc__D")
result=buscarMetabolito("xoltri25_e")
result=buscarMetabolito("lac")
result=buscarMetabolito("nh4")


# In[77]:


print 
get_ipython().magic(u'matplotlib inline')
prod_env = production_envelope(model, ["EX_gln__L_e", "EX_glc__D_e"]) 
prod_env.head(20)


# In[86]:


x=np.unique(-1*prod_env["EX_glc__D_e"])
y=np.unique(-1*prod_env["EX_gln__L_e"])
z=prod_env["flux_maximum"]

MATRIZ_SOLUCION=np.zeros((20,20))
for i in range(datos):
    print prod_env["flux_maximum"][i*20:(i+1)*20]
    MATRIZ_SOLUCION[i,0:20]=z[i*20:(i+1)*20]
print "Matriz solución"
print "====================="
print MATRIZ_SOLUCION


# In[79]:


fig  = plt.figure()
ax   = Axes3D(fig)
x,y = np.meshgrid(x, y)
surf = ax.plot_surface(x, y, MATRIZ_SOLUCION, rstride=1, cstride=1, linewidth=1, antialiased=True, cmap=plt.cm.CMRmap)
ax.view_init(elev =None, azim =60)
ax.set_xlabel("Glucosa mmol/10^6 cel/h")
ax.set_ylabel("Glutamina mmol/10^6 cel/h")
ax.set_zlabel("Velocidad de crecimiento [1/h]")
plt.title("Crecimiento celular")
plt.show()


# In[61]:


xlabels=np.unique(prod_env["EX_glc__D_e"])
ylabels=np.unique(prod_env["EX_gln__L_e"])

x_int=list(np.int_(xlabels))
y_int=list(np.int_(ylabels))

MATRIZ=np.zeros((20,20))
for i in range(20):
    MATRIZ[i,0:20]=prod_env["flux_maximum"][i*20:(i+1)*20]
print MATRIZ


# In[80]:


fig,ax=plt.subplots()
ax.set_xticks(np.arange(len(x_int)))
ax.set_yticks(np.arange(len(y_int)))

ax.set_xticklabels(x_int)
ax.set_yticklabels(y_int)   

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
ax.set_xlabel("Glucosa")
ax.set_ylabel("Glutamina")
ax.set_title("Crecimiento celular")
MATRIZ[np.isnan(MATRIZ)]=0
ax.imshow(MATRIZ)


# # 2. Simulación del modelo iMM1415 con datos experimentales

# In[81]:


#Glutamina (EX_gln_D_e)
value = -0.019 #(mmol/10^6cel/h)
model.reactions.get_by_id("EX_gln__L_e").upper_bound=value+0.1*value 
model.reactions.get_by_id("EX_gln__L_e").lower_bound=value-0.1*value

#Glucosa (EX_glc_D_e), se consume por eso es signo negativo
value = -0.186 #(mmol/10^6cel/h)
model.reactions.get_by_id("EX_glc__D_e").upper_bound=value+0.1*value 
model.reactions.get_by_id("EX_glc__D_e").lower_bound=value-0.1*value

solution=model.optimize()
print "Distribución de flujos con datos experimentales mediante FBA"
print "============================================================="
model.summary()


# In[87]:


v_p = np.array([solution.fluxes["EX_glc__D_e"], solution.fluxes["EX_gln__L_e"]])
v_e= np.array ([0.186,0.019])
#Trasponer matrices
d = (v_e - v_p)
Norma_Euclideana = np.dot(d,d)
print "Error de ajuste de los datos experimentales al modelo iMM1415"
print "====================================================="
print Norma_Euclideana


# In[ ]:




