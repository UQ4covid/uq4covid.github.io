from metawards.extractors import extract_default
from metawards import Networks, OutputFiles, Population, Workspace
from metawards.utils import Console
from copy import deepcopy

_zero_crossings = {}

def create_tables(network: Networks):
    # return a function that creates the tables
    # for the specified number of disease classes
    # Primary key needs to be composite here
    
    ## set specific column names
    col_names = ["Einc", "E", "Iinc", "I", "RI", "DI", "Ainc", "A", "RA", "Hinc", "H",\
        "RH", "DH", "Cinc", "C", "RC", "DC"]
    col_names = [f"{i} int" for i in col_names]
    col_str = ','.join(col_names)

    def initialise(conn):
        c = conn.cursor()
        table_name: str = "compact"
        values: List[str] = []
        c.execute(f"create table {table_name}(day int not null, ward int not null, {col_str}, "
                  f"primary key (day, ward))")
        conn.commit()

    return initialise


def output_db(population: Population, network: Networks,
              workspace: Workspace, output_dir: OutputFiles, **kwargs):
    Console.print(f"Calling output_db for a {network.__class__} object")
    conn3 = output_dir.open_db("stages.db", initialise=create_tables(network))
    c3 = conn3.cursor()
    
    ## setup marker for previous daya
    for i, subnet in enumerate(network.subnets):
        ## if first day, then create a copy of the ward data
        ## these should all be zero at day 0 and so not affect incidence
        if not hasattr(workspace.subspaces[i], "output_previous"):
            workspace.subspaces[i].output_previous = deepcopy(workspace.subspaces[i].ward_inf_tot)
            for j, sub in enumerate(workspace.subspaces[i].output_previous):
                for k, sub1 in enumerate(workspace.subspaces[i].output_previous[j]):
                   workspace.subspaces[i].output_previous[j][k] = 0
    
    ## extract Rprime and Dprime in each demographic
    Rprime = [[] for _ in range(4)]
    Dprime = [[] for _ in range(4)]
    for i, subnet in enumerate(network.subnets):

        ## get yesterday's data
        ward_inf_previous = workspace.subspaces[i].output_previous
        ## get today's data
        ward_inf_tot = workspace.subspaces[i].ward_inf_tot
        
        ## new deaths across wards
        if subnet.name != "asymp":
            for old, new in zip(ward_inf_previous[3], ward_inf_tot[3]):
                Dprime[i].append(new - old)
        
        ## new removals across wards
        for old, new in zip(ward_inf_previous[2], ward_inf_tot[2]):
            Rprime[i].append(new - old)
    
    ## calculate Iprime in each demographic
    Iprime = [[] for _ in range(4)]
    
    ## NEED TO DO FOLLOWING CALCULATIONS IN ORDER
    
    ## extract subnets names
    sub_names = [network.subnets[i].name for i in range(4)]
    
    ## ASYMPTOMATICS
    ia = sub_names.index("asymp")
    ## get data
    ward_inf_previous = workspace.subspaces[ia].output_previous
    ward_inf_tot = workspace.subspaces[ia].ward_inf_tot
    ## calculate incidence
    for old, new, Rinc in zip(ward_inf_previous[1], ward_inf_tot[1], Rprime[ia]):
        Iprime[ia].append(new - old + Rinc)
        
    ## CRITICAL
    ic = sub_names.index("critical")
    ## get data
    ward_inf_previous = workspace.subspaces[ic].output_previous
    ward_inf_tot = workspace.subspaces[ic].ward_inf_tot
    ## calculate incidence
    for old, new, Rinc, Dinc in zip(ward_inf_previous[1], ward_inf_tot[1], Rprime[ic], Dprime[ic]):
        Iprime[ic].append(new - old + Rinc + Dinc)
        
    ## HOSPITAL
    ih = sub_names.index("hospital")
    ## get data
    ward_inf_previous = workspace.subspaces[ih].output_previous
    ward_inf_tot = workspace.subspaces[ih].ward_inf_tot
    ## calculate incidence
    for old, new, Rinc, Dinc, Cinc in zip(ward_inf_previous[1], ward_inf_tot[1], Rprime[ih], Dprime[ih], Iprime[ic]):
        Iprime[ih].append(new - old + Rinc + Dinc + Cinc)
        
    ## GENPOP
    ig = sub_names.index("genpop")
    ## get data
    ward_inf_previous = workspace.subspaces[ig].output_previous
    ward_inf_tot = workspace.subspaces[ig].ward_inf_tot
    ## calculate incidence
    for old, new, Rinc, Dinc, Hinc in zip(ward_inf_previous[1], ward_inf_tot[1], Rprime[ig], Dprime[ig], Iprime[ih]):
        Iprime[ig].append(new - old + Rinc + Dinc + Hinc)
        
    ## calculate Eprime in GENPOP demographic
    Eprime = []
    for old, new, Iinc, Ainc in zip(ward_inf_previous[0], ward_inf_tot[0], Iprime[ig], Iprime[ia]):
        Eprime.append(new - old + Iinc + Ainc)
        
    ## loop over wards and write to file
    wards = range(0, workspace.subspaces[0].nnodes + 1)
    day = [population.day] * len(wards)
    # print(workspace.subspaces[0].nnodes)
    # print(len(day))
    # print(len(wards))
    # print(len(Eprime))
    
    ## set column names
    col_names = ["day", "ward", "Einc", "E", "Iinc", "I", "RI", "DI", "Ainc", "A", "RA", "Hinc", "H",\
        "RH", "DH", "Cinc", "C", "RC", "DC"]
    col_str = ','.join(col_names)
    
    ## extract demographics
    asymp_ward = workspace.subspaces[ia].ward_inf_tot
    genpop_ward = workspace.subspaces[ig].ward_inf_tot
    hospital_ward = workspace.subspaces[ih].ward_inf_tot
    critical_ward = workspace.subspaces[ic].ward_inf_tot
    
    ## write to file
    for day, ward, Einc, E, Iinc, I, RI, DI, Ainc, A, RA, Hinc, H, RH, RD, Cinc, C, RC, DC in\
    zip(day, wards, Eprime, genpop_ward[0], Iprime[ig], genpop_ward[1], genpop_ward[2], genpop_ward[3],\
    Iprime[ia], asymp_ward[1], asymp_ward[2], Iprime[ih], hospital_ward[1], hospital_ward[2],\
    hospital_ward[3], Iprime[ic], critical_ward[1], critical_ward[2], critical_ward[3]):
        if ward not in _zero_crossings:
            _zero_crossings[ward] = False
            
        ## try to fudge a marker for first infections
        if Einc != 0 and _zero_crossings[ward] is False and ward != 0:
            _zero_crossings[ward] = True
            Console.print(f"Got first infection in ward {ward}")

        val = [day, ward, Einc, E, Iinc, I, RI, DI, Ainc, A, RA, Hinc, H, RH, RD, Cinc, C, RC, DC]
        keeps_str = ",".join([str(v) for v in val])
        qstring = f"insert into compact ({col_str}) values ({keeps_str}) "
        if _zero_crossings[ward] is True:
            c3.execute(qstring)
            
    conn3.commit()
    
    ## save today's data so that it can be used tomorrow
    for i, subnet in enumerate(network.subnets):
        workspace.subspaces[i].output_previous = deepcopy(workspace.subspaces[i].ward_inf_tot)

def extract_db(**kwargs):
    funcs = []
    funcs.append(output_db)
    return funcs
