from metawards.movers import go_stage

def move_pathways(network, **kwargs):

    # extract user defined parameters
    params = network.params

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
    
    #######################################################
    #########              C MOVES                #########
    #######################################################
                                      
    ## move C critical to R critical
    tpCR = pC * pCR
    func.append(lambda **kwargs: go_stage(go_from="critical",
                                      go_to="critical",
                                      from_stage=1,
                                      to_stage=2,
                                      fraction=tpCR,
                                      **kwargs))
                                      
    ## move C critical to D critical
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpCD = pC * (1.0 - pCR) / (1.0 - tpCR)
    tpCD = 1.0 if tpCD > 1.0 else tpCD
    tpCD = 0.0 if tpCD < 0.0 else tpCD
    func.append(lambda **kwargs: go_stage(go_from="critical",
                                      go_to="critical",
                                      from_stage=1,
                                      to_stage=3,
                                      fraction=tpCD,
                                      **kwargs))
                                      
    #######################################################
    #########              H MOVES                #########
    #######################################################
    
    ## move H hospital to C critical
    tpHC = pH * pHC
    func.append(lambda **kwargs: go_stage(go_from="hospital",
                                      go_to="critical",
                                      from_stage=1,
                                      to_stage=1,
                                      fraction=tpHC,
                                      **kwargs))
                                      
    ## move H hospital to R hospital
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHR = pH * pHR / (1.0 - tpHC)
    tpHR = 1.0 if tpHR > 1.0 else tpHR
    tpHR = 0.0 if tpHR < 0.0 else tpHR
    func.append(lambda **kwargs: go_stage(go_from="hospital",
                                      go_to="hospital",
                                      from_stage=1,
                                      to_stage=2,
                                      fraction=tpHR,
                                      **kwargs))
                                      
    ## move H hospital to D hospital
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpHD = pH * (1.0 - pHC - pHR) / (1.0 - pH * (pHC + pHR))
    tpHD = 1.0 if tpHD > 1.0 else tpHD
    tpHD = 0.0 if tpHD < 0.0 else tpHD
    func.append(lambda **kwargs: go_stage(go_from="hospital",
                                      go_to="hospital",
                                      from_stage=1,
                                      to_stage=3,
                                      fraction=tpHD,
                                      **kwargs))
    
    #######################################################
    #########              I MOVES                #########
    #######################################################

    ## move I genpop to H hospital
    tpIH = pI * pIH
    func.append(lambda **kwargs: go_stage(go_from="genpop",
                                      go_to="hospital",
                                      from_stage=1,
                                      to_stage=1,
                                      fraction=tpIH,
                                      **kwargs))

    ## move I genpop to R genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpIR = pI * pIR / (1.0 - tpIH)
    tpIR = 1.0 if tpIR > 1.0 else tpIR
    tpIR = 0.0 if tpIR < 0.0 else tpIR
    func.append(lambda **kwargs: go_stage(go_from="genpop",
                                      go_to="genpop",
                                      from_stage=1,
                                      to_stage=2,
                                      fraction=tpIR,
                                      **kwargs))

    ## move I genpop to D genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpID = pI * (1 - pIH - pIR) / (1.0 - pI * (pIH + pIR))
    tpID = 1.0 if tpID > 1.0 else tpID
    tpID = 0.0 if tpID < 0.0 else tpID
    func.append(lambda **kwargs: go_stage(go_from="genpop",
                                      go_to="genpop",
                                      from_stage=1,
                                      to_stage=3,
                                      fraction=tpID,
                                      **kwargs))
    
    #######################################################
    #########              A MOVES                #########
    #######################################################

    ## move A asymp to R asymp
    tpAR = pA
    func.append(lambda **kwargs: go_stage(go_from="asymp",
                                      go_to="asymp",
                                      from_stage=1,
                                      to_stage=2,
                                      fraction=tpAR,
                                      **kwargs))
    
    #######################################################
    #########              E MOVES                #########
    #######################################################

    ## move E genpop to A asymp
    tpEA = pE * pEA
    func.append(lambda **kwargs: go_stage(go_from="genpop",
                                      go_to="asymp",
                                      from_stage=0,
                                      to_stage=1,
                                      fraction=tpEA,
                                      **kwargs))

    ## move E genpop to I genpop
    ## (denominator adjustment is due to operating on remainder
    ## as described in the vignette, also includes correction
    ## in case of rounding error)
    tpEI = pE * (1.0 - pEA) / (1.0 - tpEA)
    tpEI = 1.0 if tpEI > 1.0 else tpEI
    tpEI = 0.0 if tpEI < 0.0 else tpEI
    func.append(lambda **kwargs: go_stage(go_from="genpop",
                                      go_to="genpop",
                                      from_stage=0,
                                      to_stage=1,
                                      fraction=tpEI,
                                      **kwargs))

    return func
