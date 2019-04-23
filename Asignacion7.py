import numpy as np
import cobra


def miFuncion():
	model = cobra.io.read_sbml_model("iMM904.xml")
	solution=model.optimize()	
	v_p = np.array([solution.fluxes["EX_glyc_e"], solution.fluxes["EX_ac_e"], solution.fluxes["EX_acald_e"]])
	v_e= np.array ([0.211,0.0078,0.0048])

	d = (v_e - v_p)
	Euclidean_norm= np.dot(d,d)
	print Euclidean_norm

    	return (1)

if __name__ == "__main__":
    miFuncion()





