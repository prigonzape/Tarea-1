import cobra

import cobra.test


def miFuncion():
 
    model = cobra.test.create_test_model("ecoli")
  

    print "Running FVA=============="

    solution = model.optimize()

    model.summary(fva=1)


    print "Running Geometric FBA ============"

    geometric_fba_sol = cobra.flux_analysis.geometric_fba(model)

    geometric_fba_sol

    return 1


if __name__ == "__main__":

    miFuncion()