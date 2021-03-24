from metawards.movers import go_stage
from metawards.utils import Console

def move_pathways(network, **kwargs):

    # extract user defined parameters
    params = network.params
    
    ## extract number of age-classes
    nage = params.user_params["nage"]

    ## moves out of E class
    pE = []
    pEP = []
    for j in range(nage):
        pE.append(params.user_params[f'pE_{j + 1}'])
        pEP.append(params.user_params[f'pEP_{j + 1}'])
    
    ## moves out of A class
    pA = []
    for j in range(nage):
        pA.append(params.user_params[f'pA_{j + 1}'])
    
    ## moves out of P class
    pP = []
    for j in range(nage):
        pP.append(params.user_params[f'pP_{j + 1}'])
    
    ## moves out of I1 class
    pI1 = []
    pI1H = []
    pI1D = []
    for j in range(nage):
        pI1.append(params.user_params[f'pI1_{j + 1}'])
        pI1H.append(params.user_params[f'pI1H_{j + 1}'])
        pI1D.append(params.user_params[f'pI1D_{j + 1}'])
    
    ## moves out of I2 class
    pI2 = []
    for j in range(nage):
        pI2.append(params.user_params[f'pI2_{j + 1}'])
    
    ## moves out of H class
    pH = []
    pHD = []
    for j in range(nage):
        pH.append(params.user_params[f'pH_{j + 1}'])
        pHD.append(params.user_params[f'pHD_{j + 1}'])
        
    func = []
    
    ## moves in reverse order through the stages to 
    ## ensure correct mapping
    
    ## NOTE: looping over the age ranges don't work unless you use a 
    ## default parameter k = j in the lambda - solution derived from:
    ## https://stackoverflow.com/questions/10452770/python-lambdas-binding-to-local-values
                                      
    #######################################################
    #########              H MOVES                #########
    #######################################################
                                      
    ## move H genpop to DH genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHD = []
    for j in range(nage):
        tpHD.append(pH[j] * pHD[j])
        tpHD[j] = 1.0 if tpHD[j] > 1.0 else tpHD[j]
        tpHD[j] = 0.0 if tpHD[j] < 0.0 else tpHD[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="H",
                                      to_stage="DH",
                                      fraction=tpHD[k],
                                      **kwargs))
                                      
    ## move H genpop to R genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHR = []
    for j in range(nage):
        tpHR.append(pH[j] * (1.0 - pHD[j]) / (1.0 - tpHD[j]))
        tpHR[j] = 0.0 if tpHD[j] == 1.0 else tpHR[j]
        tpHR[j] = 1.0 if tpHR[j] > 1.0 else tpHR[j]
        tpHR[j] = 0.0 if tpHR[j] < 0.0 else tpHR[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="H",
                                      to_stage="RH",
                                      fraction=tpHR[k],
                                      **kwargs))
                                      
    #######################################################
    #########              I2 MOVES               #########
    #######################################################

    ## move I2 genpop to RI genpop
    tpI2R = []
    for j in range(nage):
        tpI2R.append(pI2[j])
        tpI2R[j] = 1.0 if tpI2R[j] > 1.0 else tpI2R[j]
        tpI2R[j] = 0.0 if tpI2R[j] < 0.0 else tpI2R[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="I2",
                                      to_stage="RI",
                                      fraction=tpI2R[k],
                                      **kwargs))
    
    #######################################################
    #########              I1 MOVES               #########
    #######################################################

    ## move I1 genpop to H genpop
    tpI1H = []
    for j in range(nage):
        tpI1H.append(pI1[j] * pI1H[j])
        tpI1H[j] = 1.0 if tpI1H[j] > 1.0 else tpI1H[j]
        tpI1H[j] = 0.0 if tpI1H[j] < 0.0 else tpI1H[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="I1",
                                      to_stage="H",
                                      fraction=tpI1H[k],
                                      **kwargs))

    ## move I1 genpop to DI genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpI1D = []
    for j in range(nage):
        tpI1D.append(pI1[j] * pI1D[j] / (1.0 - tpI1H[j]))
        tpI1D[j] = 0.0 if tpI1H[j] == 1.0 else tpI1D[j]
        tpI1D[j] = 1.0 if tpI1D[j] > 1.0 else tpI1D[j]
        tpI1D[j] = 0.0 if tpI1D[j] < 0.0 else tpI1D[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="I1",
                                      to_stage="DI",
                                      fraction=tpI1D[k],
                                      **kwargs))

    ## move I1 genpop to I2 genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpI1I2 = []
    for j in range(nage):
        tpI1I2.append(pI1[j] * (1.0 - pI1H[j] - pI1D[j]) / (1.0 - pI1[j] * (pI1H[j] + pI1D[j])))
        tpI1I2[j] = 0.0 if (pI1[j] * (pI1H[j] + pI1D[j])) == 1.0 else tpI1I2[j]
        tpI1I2[j] = 1.0 if tpI1I2[j] > 1.0 else tpI1I2[j]
        tpI1I2[j] = 0.0 if tpI1I2[j] < 0.0 else tpI1I2[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="I1",
                                      to_stage="I2",
                                      fraction=tpI1I2[k],
                                      **kwargs))
    
    #######################################################
    #########              P MOVES                #########
    #######################################################

    ## move P genpop to I1 genpop
    tpPI1 = []
    for j in range(nage):
        tpPI1.append(pP[j])
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="P",
                                      to_stage="I1",
                                      fraction=tpPI1[k],
                                      **kwargs))
    
    #######################################################
    #########              A MOVES                #########
    #######################################################

    ## move A genpop to R genpop
    tpAR = []
    for j in range(nage):
        tpAR.append(pA[j])
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="A",
                                      to_stage="RA",
                                      fraction=tpAR[k],
                                      **kwargs))
    
    #######################################################
    #########              E MOVES                #########
    #######################################################

    ## move E genpop to A genpop
    tpEA = []
    for j in range(nage):
        tpEA.append(pE[j] * (1.0 - pEP[j]))
        tpEA[j] = 1.0 if tpEA[j] > 1.0 else tpEA[j]
        tpEA[j] = 0.0 if tpEA[j] < 0.0 else tpEA[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="E",
                                      to_stage="A",
                                      fraction=tpEA[k],
                                      **kwargs))

    ## move E genpop to P genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpEP = []
    for j in range(nage):
        tpEP.append(pE[j] * pEP[j] / (1.0 - tpEA[j]))
        tpEP[j] = 0.0 if tpEA[j] == 1.0 else tpEP[j]
        tpEP[j] = 1.0 if tpEP[j] > 1.0 else tpEP[j]
        tpEP[j] = 0.0 if tpEP[j] < 0.0 else tpEP[j]
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'age{k + 1}',
                                      go_to=f'age{k + 1}',
                                      from_stage="E",
                                      to_stage="P",
                                      fraction=tpEP[k],
                                      **kwargs))

    return func
