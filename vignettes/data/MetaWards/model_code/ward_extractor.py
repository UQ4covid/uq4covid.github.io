from metawards.extractors import extract_default
from metawards import Networks, OutputFiles, Population, Workspace
from metawards.utils import Console
from copy import deepcopy

_zero_crossings = {}

def create_tables(network: Networks):
    ## return a function that creates the tables
    ## for the specified number of disease classes
    ## Primary key needs to be composite here
    
    ## set specific column names
    col_names = ["Einc", "Pinc", "I1inc", "I2inc", "RIinc", "DIinc", "Ainc", "RAinc", "Hinc", "RHinc", "DHinc"]
    col_names = [f"{i} int" for i in col_names]
    col_str = ','.join(col_names)
    
    ## set specific column names
    col_names_ini = ["E", "P", "I1", "I2", "RI", "DI", "A", "RA", "H", "RH", "DH"]
    col_names_ini = [f"{i} int" for i in col_names_ini]
    col_str_ini = ','.join(col_names_ini)

    def initialise(conn):
        c = conn.cursor()
        c.execute(f"create table compact(day int not null, ward int not null, {col_str}, "
                  f"primary key (day, ward))")
        c.execute(f"create table compact_ini(day int not null, ward int not null, {col_str_ini}, "
                  f"primary key (day, ward))")
        conn.commit()

    return initialise

def output_db(population: Population, network: Networks,
              workspace: Workspace, output_dir: OutputFiles, **kwargs):
    Console.print(f"Calling output_db for a {network.__class__} object")
    
    ## extract number of age-classes
    nage = network.params.user_params["nage"]
        
    ## setup marker for previous day
    for i, subnet in enumerate(network.subnets):
        ## if first day, then create a copy of the ward data
        ## this should produce incidences of zero on day 1
        if not hasattr(workspace.subspaces[i], "output_previous"):
            workspace.subspaces[i].output_previous = deepcopy(workspace.subspaces[i].ward_inf_tot)
    
    if population.day == 1:
        ## open connection to databases
        conn = []
        c = []
        for i in range(nage):
            conn.append(output_dir.open_db(f"age{i + 1}.db", initialise = create_tables(network)))
            c.append(conn[i].cursor())
            
            ## get data
            ward_inf_tot = workspace.subspaces[i].ward_inf_tot
            
            ## extract initial counts
            E = ward_inf_tot[0]
            P = ward_inf_tot[1]
            I1 = ward_inf_tot[2]
            I2 = ward_inf_tot[3]
            RI = ward_inf_tot[4]
            DI = ward_inf_tot[5]
            A = ward_inf_tot[6]
            RA = ward_inf_tot[7]
            H = ward_inf_tot[8]
            RH = ward_inf_tot[9]
            DH = ward_inf_tot[10]
                
            ## loop over wards and write to file
            wards = range(0, workspace.subspaces[0].nnodes + 1)
            day = [population.day] * len(wards)
            
            ## set column names
            col_names = ["day", "ward", "E", "P", "I1", "I2", "RI", "DI", "A", "RA", "H", "RH", "DH"]
            col_str = ','.join(col_names)
            
            ## write to file
            for day, ward, E, P, I1, I2, RI, DI, A, RA, H, RH, DH in\
            zip(day, wards, E, P, I1, I2, RI, DI, A, RA,\
            H, RH, DH):
                if ward not in _zero_crossings:
                    _zero_crossings[ward] = False
                    
                ## set up list of changes
                val = [day, ward, E, P, I1, I2, RI, DI, A, RA, H, RH, DH]
                keeps_str = ",".join([str(v) for v in val])
                qstring = f"insert into compact_ini ({col_str}) values ({keeps_str}) "
                    
                ## try to fudge a marker for first infections
                if any([ v > 0 for v in val[2:] ]) and _zero_crossings[ward] is False and ward != 0:
                    _zero_crossings[ward] = True
                    Console.print(f"Got first infection in ward {ward}")
                
                ## check for any changes in ward    
                if _zero_crossings[ward] is True and any([ v > 0 for v in val[2:] ]):
                    c[i].execute(qstring)
                    
            conn[i].commit()
        
        ## save today's data so that it can be used tomorrow
        for i, subnet in enumerate(network.subnets):
            workspace.subspaces[i].output_previous = deepcopy(workspace.subspaces[i].ward_inf_tot)
    else:
        ## open connection to databases
        conn = []
        c = []
        for i in range(nage):
            conn.append(output_dir.open_db(f"age{i + 1}.db", initialise = create_tables(network)))
            c.append(conn[i].cursor())
        
        ## NEED TO DO FOLLOWING CALCULATIONS IN ORDER
        
        ## loop over age classes
        for j in range(nage):
            
            ## get data
            ward_inf_previous = workspace.subspaces[j].output_previous
            ward_inf_tot = workspace.subspaces[j].ward_inf_tot
            
            ## ASYMPTOMATICS
            ## calculate incidence
            RAprime = []
            for old, new in zip(ward_inf_previous[7], ward_inf_tot[7]):
                RAprime.append(new - old)
            Aprime = []
            for old, new, Rinc in zip(ward_inf_previous[6], ward_inf_tot[6], RAprime):
                Aprime.append(new - old + Rinc)
            
            ## HOSPITAL
            ## calculate incidence
            RHprime = []
            for old, new in zip(ward_inf_previous[9], ward_inf_tot[9]):
                RHprime.append(new - old)
            DHprime = []
            for old, new in zip(ward_inf_previous[10], ward_inf_tot[10]):
                DHprime.append(new - old)
            Hprime = []
            for old, new, Rinc, Dinc in zip(ward_inf_previous[8], ward_inf_tot[8], RHprime, DHprime):
                Hprime.append(new - old + Rinc + Dinc)
            
            ## GENPOP
            ## calculate incidence
            RIprime = []
            for old, new in zip(ward_inf_previous[4], ward_inf_tot[4]):
                RIprime.append(new - old)
            DIprime = []
            for old, new in zip(ward_inf_previous[5], ward_inf_tot[5]):
                DIprime.append(new - old)
            I2prime = []
            for old, new, Rinc in zip(ward_inf_previous[3], ward_inf_tot[3], RIprime):
                I2prime.append(new - old + Rinc)
            I1prime = []
            for old, new, I2inc, Dinc, Hinc in zip(ward_inf_previous[2], ward_inf_tot[2], I2prime, DIprime, Hprime):
                I1prime.append(new - old + I2inc + Dinc + Hinc)
            Pprime = []
            for old, new, I1inc in zip(ward_inf_previous[1], ward_inf_tot[1], I1prime):
                Pprime.append(new - old + I1inc)
            Eprime = []
            for old, new, Pinc, Ainc in zip(ward_inf_previous[0], ward_inf_tot[0], Pprime, Aprime):
                Eprime.append(new - old + Pinc + Ainc)
                
            ## loop over wards and write to file
            wards = range(0, workspace.subspaces[0].nnodes + 1)
            day = [population.day] * len(wards)
            
            ## set column names
            col_names = ["day", "ward", "Einc", "Pinc", "I1inc", "I2inc", "RIinc", "DIinc", "Ainc", "RAinc", "Hinc", "RHinc", "DHinc"]
            col_str = ','.join(col_names)
            
            ## write to file
            for day, ward, Einc, Pinc, I1inc, I2inc, RIinc, DIinc, Ainc, RAinc, Hinc, RHinc, DHinc in\
            zip(day, wards, Eprime, Pprime, I1prime, I2prime, RIprime, DIprime, Aprime, RAprime,\
            Hprime, RHprime, DHprime):
                if ward not in _zero_crossings:
                    _zero_crossings[ward] = False
                    
                ## try to fudge a marker for first infections
                if Einc != 0 and _zero_crossings[ward] is False and ward != 0:
                    _zero_crossings[ward] = True
                    Console.print(f"Got first infection in ward {ward}")
                    
                ## set up list of changes
                val = [day, ward, Einc, Pinc, I1inc, I2inc, RIinc, DIinc, Ainc, RAinc, Hinc, RHinc, DHinc]
                keeps_str = ",".join([str(v) for v in val])
                qstring = f"insert into compact ({col_str}) values ({keeps_str}) "
                
                ## check for any changes in ward    
                if _zero_crossings[ward] is True and any([ v > 0 for v in val[2:] ]):
                    c[j].execute(qstring)
                    
            conn[j].commit()
        
        ## save today's data so that it can be used tomorrow
        for i, subnet in enumerate(network.subnets):
            workspace.subspaces[i].output_previous = deepcopy(workspace.subspaces[i].ward_inf_tot)

def extract_db(**kwargs):
    funcs = []
    funcs.append(output_db)
    return funcs
