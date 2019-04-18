#!/usr/bin/env python
# coding: utf-8

# TAREA: Buscar la secuencia seq dentro de una gran secuencia de DNA. DNA='ACGTCTTTCCGCATCATTCTAACTCCGCATTGTCCTTACGATACCGATTACGATTCACGATACCGCATTTACGATACGTCAGATACTCCGACATTTACGATACTCC' seq='CATCATTCTAAC'

# In[1]:


DNA="ACGTCTTTCCGCATCATTCTAACTCCGCATTGTCCTTACGATACCGATTACGATTCACGATACCGCATTTACGATACGTCAGATACTCCGACATTTACGATACTCC"
sec="CATCATTCTAAC"

#Saber cual es la longitud de DNA y de sec
print len(DNA)
print len(sec)

#Saber si sec se encuentra dentro de DNA
print sec in DNA

#Dividir DNA (string) en la posición de la secuencia 
x= (DNA.split(sec))
print(x)

#Buscar la posición de la secuencia
len(x)
len (x[0])

#La secuencia se encuentra en la posición 12 y tiene una longitud de 12 bases nitrogenadas


# In[ ]:




