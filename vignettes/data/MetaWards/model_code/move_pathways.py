from metawards.movers import go_stage
from metawards.utils import Console

def move_pathways(network, **kwargs):

    # extract user defined parameters
    params = network.params
    
    ## extract number of age-classes
    nage = params.user_params["nage"]

    ## moves out of E class
    pE = []
    pEA = []
    for j in range(nage):
        pE.append(params.user_params[f'pE_{j}'])
        pEA.append(params.user_params[f'pEA_{j}'])
    
    ## moves out of A class
    pA = []
    for j in range(nage):
        pA.append(params.user_params[f'pA_{j}'])
    
    ## moves out of I1 class
    pI1 = []
    pI1H = []
    pI1I2 = []
    for j in range(nage):
        pI1.append(params.user_params[f'pI1_{j}'])
        pI1H.append(params.user_params[f'pI1H_{j}'])
        pI1I2.append(params.user_params[f'pI1I2_{j}'])
    
    ## moves out of I2 class
    pI2 = []
    for j in range(nage):
        pI2.append(params.user_params[f'pI2_{j}'])
    
    ## moves out of H class
    pH = []
    pHR = []
    for j in range(nage):
        pH.append(params.user_params[f'pH_{j}'])
        pHR.append(params.user_params[f'pHR_{j}'])
        
    func = []
    
    ## moves in reverse order through the stages to 
    ## ensure correct mapping
    
    ## NOTE: looping over the age ranges don't work unless you use a 
    ## default parameter k = j in the lambda - solution derived from:
    ## https://stackoverflow.com/questions/10452770/python-lambdas-binding-to-local-values
                                      
    #######################################################
    #########              H MOVES                #########
    #######################################################
                                      
    ## move H hospital to R hospital
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHR = []
    for j in range(nage):
        tpHR.append(pH[j] * pHR[j])
        tpHR[j] = 1.0 if tpHR[j] > 1.0 else tpHR[j]
        tpHR[j] = 0.0 if tpHR[j] < 0.0 else tpHR[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'hospital{k + 1}',
                                      go_to=f'hospital{k + 1}',
                                      from_stage="IH",
                                      to_stage="RH",
                                      fraction=tpHR[j],
                                      **kwargs))
                                      
    ## move H hospital to D hospital
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHD = []
    for j in range(nage):
        tpHD.append(pH[j] * (1.0 - pHR[j]) / (1.0 - tpHR[j]))
        tpHD[j] = 0.0 if tpHR[j] == 1.0 else tpHD[j]
        tpHD[j] = 1.0 if tpHD[j] > 1.0 else tpHD[j]
        tpHD[j] = 0.0 if tpHD[j] < 0.0 else tpHD[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'hospital{k + 1}',
                                      go_to=f'hospital{k + 1}',
                                      from_stage="IH",
                                      to_stage="DH",
                                      fraction=tpHD[j],
                                      **kwargs))
                                      
    #######################################################
    #########              I2 MOVES               #########
    #######################################################

    ## move I2 genpop to R genpop
    tpI2R = []
    for j in range(nage):
        tpI2R.append(pI2[j])
        tpI2R[j] = 1.0 if tpI2R[j] > 1.0 else tpI2R[j]
        tpI2R[j] = 0.0 if tpI2R[j] < 0.0 else tpI2R[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="I2",
                                      to_stage="R",
                                      fraction=tpI2R[j],
                                      **kwargs))
    
    #######################################################
    #########              I1 MOVES               #########
    #######################################################

    ## move I1 genpop to H hospital
    tpI1H = []
    for j in range(nage):
        tpI1H.append(pI1[j] * pI1H[j])
        tpI1H[j] = 1.0 if tpI1H[j] > 1.0 else tpI1H[j]
        tpI1H[j] = 0.0 if tpI1H[j] < 0.0 else tpI1H[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'hospital{k + 1}',
                                      from_stage="I1",
                                      to_stage="IH",
                                      fraction=tpI1H[j],
                                      **kwargs))

    ## move I1 genpop to I2 genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpI1I2 = []
    for j in range(nage):
        tpI1I2.append(pI1[j] * pI1I2[j] / (1.0 - tpI1H[j]))
        tpI1I2[j] = 0.0 if tpI1H[j] == 1.0 else tpI1I2[j]
        tpI1I2[j] = 1.0 if tpI1I2[j] > 1.0 else tpI1I2[j]
        tpI1I2[j] = 0.0 if tpI1I2[j] < 0.0 else tpI1I2[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="I1",
                                      to_stage="I2",
                                      fraction=tpI1I2[j],
                                      **kwargs))

    ## move I1 genpop to D genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpI1D = []
    for j in range(nage):
        tpI1D.append(pI1[j] * (1 - pI1H[j] - pI1I2[j]) / (1.0 - pI1[j] * (pI1H[j] + pI1I2[j])))
        tpI1D[j] = 0.0 if (pI1[j] * (pI1H[j] + pI1I2[j])) == 1.0 else tpI1D[j]
        tpI1D[j] = 1.0 if tpI1D[j] > 1.0 else tpI1D[j]
        tpI1D[j] = 0.0 if tpI1D[j] < 0.0 else tpI1D[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="I1",
                                      to_stage="D",
                                      fraction=tpI1D[j],
                                      **kwargs))
    
    #######################################################
    #########              A MOVES                #########
    #######################################################

    ## move A asymp to R asymp
    tpAR = []
    for j in range(nage):
        tpAR.append(pA[j])
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'asymp{k + 1}',
                                      go_to=f'asymp{k + 1}',
                                      from_stage="IA",
                                      to_stage="RA",
                                      fraction=tpAR[j],
                                      **kwargs))
    
    #######################################################
    #########              E MOVES                #########
    #######################################################

    ## move E genpop to A asymp
    tpEA = []
    for j in range(nage):
        tpEA.append(pE[j] * pEA[j])
        tpEA[j] = 1.0 if tpEA[j] > 1.0 else tpEA[j]
        tpEA[j] = 0.0 if tpEA[j] < 0.0 else tpEA[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'asymp{k + 1}',
                                      from_stage="E",
                                      to_stage="IA",
                                      fraction=tpEA[j],
                                      **kwargs))

    ## move E genpop to I1 genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpEI1 = []
    for j in range(nage):
        tpEI1.append(pE[j] * (1.0 - pEA[j]) / (1.0 - tpEA[j]))
        tpEI1[j] = 0.0 if tpEA[j] == 1.0 else tpEI1[j]
        tpEI1[j] = 1.0 if tpEI1[j] > 1.0 else tpEI1[j]
        tpEI1[j] = 0.0 if tpEI1[j] < 0.0 else tpEI1[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="E",
                                      to_stage="I1",
                                      fraction=tpEI1[j],
                                      **kwargs))

    return func
