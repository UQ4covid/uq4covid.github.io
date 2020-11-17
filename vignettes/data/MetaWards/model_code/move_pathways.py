from metawards.movers import go_stage
from metawards.utils import Console

def move_pathways(network, **kwargs):

    # extract user defined parameters
    params = network.params
    
    ## extract number of age-classes
    nage = params.user_params["nage"]

    ## moves out of E class
    pE = params.user_params["pE"]
    pEA = params.user_params["pEA"]
    
    ## moves out of A class
    pA = params.user_params["pA"]
    
    ## moves out of I class
    pI = params.user_params["pI"]
    pIH = params.user_params["pIH"]
    pIR = params.user_params["pIR"]
    
    ## moves out of H class
    pH = params.user_params["pH"]
    pHC = params.user_params["pHC"]
    pHR = params.user_params["pHR"]
    
    ## moves out of C class
    pC = params.user_params["pC"]
    pCR = params.user_params["pCR"]
        
    func = []
    
    ## movers in reverse order through the stages to 
    ## ensure correct mapping
    
    ## NOTE: looping over the age ranges don't work unless you use a 
    ## default parameter k = j in the lambda - solution derived from:
    ## https://stackoverflow.com/questions/10452770/python-lambdas-binding-to-local-values
    
    #######################################################
    #########              C MOVES                #########
    #######################################################
                                      
    ## move C critical to R critical
    tpCR = pC * pCR
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'critical{k + 1}',
                                      go_to=f'critical{k + 1}',
                                      from_stage="IC",
                                      to_stage="RC",
                                      fraction=tpCR,
                                      **kwargs))
                                      
    ## move C critical to D critical
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpCD = pC * (1.0 - pCR) / (1.0 - tpCR)
    tpCD = 1.0 if tpCD > 1.0 else tpCD
    tpCD = 0.0 if tpCD < 0.0 else tpCD
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'critical{k + 1}',
                                      go_to=f'critical{k + 1}',
                                      from_stage="IC",
                                      to_stage="DC",
                                      fraction=tpCD,
                                      **kwargs))
                                      
    #######################################################
    #########              H MOVES                #########
    #######################################################
    
    ## move H hospital to C critical
    tpHC = pH * pHC
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'hospital{k + 1}',
                                      go_to=f'critical{k + 1}',
                                      from_stage="IH",
                                      to_stage="IC",
                                      fraction=tpHC,
                                      **kwargs))
                                      
    ## move H hospital to R hospital
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHR = pH * pHR / (1.0 - tpHC)
    tpHR = 1.0 if tpHR > 1.0 else tpHR
    tpHR = 0.0 if tpHR < 0.0 else tpHR
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'hospital{k + 1}',
                                      go_to=f'hospital{k + 1}',
                                      from_stage="IH",
                                      to_stage="RH",
                                      fraction=tpHR,
                                      **kwargs))
                                      
    ## move H hospital to D hospital
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHD = pH * (1.0 - pHC - pHR) / (1.0 - pH * (pHC + pHR))
    tpHD = 1.0 if tpHD > 1.0 else tpHD
    tpHD = 0.0 if tpHD < 0.0 else tpHD
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'hospital{k + 1}',
                                      go_to=f'hospital{k + 1}',
                                      from_stage="IH",
                                      to_stage="DH",
                                      fraction=tpHD,
                                      **kwargs))
    
    #######################################################
    #########              I MOVES                #########
    #######################################################

    ## move I genpop to H hospital
    tpIH = pI * pIH
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'hospital{k + 1}',
                                      from_stage="I",
                                      to_stage="IH",
                                      fraction=tpIH,
                                      **kwargs))

    ## move I genpop to R genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpIR = pI * pIR / (1.0 - tpIH)
    tpIR = 1.0 if tpIR > 1.0 else tpIR
    tpIR = 0.0 if tpIR < 0.0 else tpIR
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="I",
                                      to_stage="R",
                                      fraction=tpIR,
                                      **kwargs))

    ## move I genpop to D genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpID = pI * (1 - pIH - pIR) / (1.0 - pI * (pIH + pIR))
    tpID = 1.0 if tpID > 1.0 else tpID
    tpID = 0.0 if tpID < 0.0 else tpID
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="I",
                                      to_stage="D",
                                      fraction=tpID,
                                      **kwargs))
    
    #######################################################
    #########              A MOVES                #########
    #######################################################

    ## move A asymp to R asymp
    tpAR = pA
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'asymp{k + 1}',
                                      go_to=f'asymp{k + 1}',
                                      from_stage="IA",
                                      to_stage="RA",
                                      fraction=tpAR,
                                      **kwargs))
    
    #######################################################
    #########              E MOVES                #########
    #######################################################

    ## move E genpop to A asymp
    tpEA = pE * pEA
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'asymp{k + 1}',
                                      from_stage="E",
                                      to_stage="IA",
                                      fraction=tpEA,
                                      **kwargs))

    ## move E genpop to I genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpEI = pE * (1.0 - pEA) / (1.0 - tpEA)
    tpEI = 1.0 if tpEI > 1.0 else tpEI
    tpEI = 0.0 if tpEI < 0.0 else tpEI
    for j in range(nage):
        func.append(lambda k = j, **kwargs: go_stage(go_from=f'genpop{k + 1}',
                                      go_to=f'genpop{k + 1}',
                                      from_stage="E",
                                      to_stage="I",
                                      fraction=tpEI,
                                      **kwargs))

    return func
