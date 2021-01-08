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
    col_names = ["Einc", "Iinc", "Rinc", "Dinc", "IAinc", "RAinc", "IHinc", "RHinc", "DHinc"]
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
    
    ## extract number of age-classes
    nage = network.params.user_params["nage"]
    
    ## open connection to databases
    conn = []
    c = []
    for i in range(nage):
        conn.append(output_dir.open_db(f"stages{i + 1}.db", initialise = create_tables(network)))
        c.append(conn[i].cursor())
    
    ## setup marker for previous day
    for i, subnet in enumerate(network.subnets):
        ## if first day, then create a copy of the ward data
        ## these should all be zero at day 0 and so not affect incidence
        if not hasattr(workspace.subspaces[i], "output_previous"):
            workspace.subspaces[i].output_previous = deepcopy(workspace.subspaces[i].ward_inf_tot)
            for j, sub in enumerate(workspace.subspaces[i].output_previous):
                for k, sub1 in enumerate(workspace.subspaces[i].output_previous[j]):
                   workspace.subspaces[i].output_previous[j][k] = 0
    
    ## extract Rprime and Dprime in each demographic
    Rprime = [[] for _ in range(len(network.subnets))]
    Dprime = [[] for _ in range(len(network.subnets))]
    for i, subnet in enumerate(network.subnets):

        ## get yesterday's data
        ward_inf_previous = workspace.subspaces[i].output_previous
        ## get today's data
        ward_inf_tot = workspace.subspaces[i].ward_inf_tot
        
        ## remove digits from demographic name
        subName = ''.join([j for j in subnet.name if not j.isdigit()])
        
        if subName == "genpop":
            ## new deaths across wards
            for old, new in zip(ward_inf_previous[3], ward_inf_tot[3]):
                Dprime[i].append(new - old)
            ## new removals across wards
            for old, new in zip(ward_inf_previous[2], ward_inf_tot[2]):
                Rprime[i].append(new - old)
        
        if subName == "asymp":
            ## new removals across wards
            for old, new in zip(ward_inf_previous[1], ward_inf_tot[1]):
                Rprime[i].append(new - old)
        
        if subName == "hospital":
            ## new deaths across wards
            for old, new in zip(ward_inf_previous[2], ward_inf_tot[2]):
                Dprime[i].append(new - old)
            ## new removals across wards
            for old, new in zip(ward_inf_previous[1], ward_inf_tot[1]):
                Rprime[i].append(new - old)
    
    ## calculate Iprime in each demographic
    Iprime = [[] for _ in range(len(network.subnets))]
    
    ## extract subnet names
    sub_names = [network.subnets[i].name for i in range(len(network.subnets))]
    
    ## NEED TO DO FOLLOWING CALCULATIONS IN ORDER
    
    ## loop over age classes
    for j in range(nage):
        
        ## ASYMPTOMATICS
        ia = sub_names.index(f'asymp{j + 1}')
        ## get data
        ward_inf_previous = workspace.subspaces[ia].output_previous
        ward_inf_tot = workspace.subspaces[ia].ward_inf_tot
        ## calculate incidence
        for old, new, Rinc in zip(ward_inf_previous[0], ward_inf_tot[0], Rprime[ia]):
            Iprime[ia].append(new - old + Rinc)
        
        ## HOSPITAL
        ih = sub_names.index(f'hospital{j + 1}')
        ## get data
        ward_inf_previous = workspace.subspaces[ih].output_previous
        ward_inf_tot = workspace.subspaces[ih].ward_inf_tot
        ## calculate incidence
        for old, new, Rinc, Dinc in zip(ward_inf_previous[0], ward_inf_tot[0], Rprime[ih], Dprime[ih]):
            Iprime[ih].append(new - old + Rinc + Dinc)
        
        ## GENPOP
        ig = sub_names.index(f'genpop{j + 1}')
        ## get data
        ward_inf_previous = workspace.subspaces[ig].output_previous
        ward_inf_tot = workspace.subspaces[ig].ward_inf_tot
        ## calculate incidence
        for old, new, Rinc, Dinc, Hinc in zip(ward_inf_previous[1], ward_inf_tot[1], Rprime[ig], Dprime[ig], Iprime[ih]):
            Iprime[ig].append(new - old + Rinc + Dinc + Hinc)
        
        # calculate Eprime in GENPOP demographic
        Eprime = []
        for old, new, Iinc, Ainc in zip(ward_inf_previous[0], ward_inf_tot[0], Iprime[ig], Iprime[ia]):
            Eprime.append(new - old + Iinc + Ainc)
            
        ## loop over wards and write to file
        wards = range(0, workspace.subspaces[0].nnodes + 1)
        day = [population.day] * len(wards)
        
        ## set column names
        col_names = ["day", "ward", "Einc", "Iinc", "Rinc", "Dinc", "IAinc", "RAinc", "IHinc", "RHinc", "DHinc"]
        col_str = ','.join(col_names)
        
        ## extract demographics
        asymp_ward = workspace.subspaces[ia].ward_inf_tot
        genpop_ward = workspace.subspaces[ig].ward_inf_tot
        hospital_ward = workspace.subspaces[ih].ward_inf_tot
        
        ## write to file
        for day, ward, Einc, Iinc, Rinc, Dinc, IAinc, RAinc, IHinc, RHinc, DHinc in\
        zip(day, wards, Eprime, Iprime[ig], Rprime[ig], Dprime[ig], Iprime[ia], Rprime[ia],\
        Iprime[ih], Rprime[ih], Dprime[ih]):
            if ward not in _zero_crossings:
                _zero_crossings[ward] = False
                
            ## try to fudge a marker for first infections
            if Einc != 0 and _zero_crossings[ward] is False and ward != 0:
                _zero_crossings[ward] = True
                Console.print(f"Got first infection in ward {ward}")
                
            ## set up list of changes
            val = [day, ward, Einc, Iinc, Rinc, Dinc, IAinc, RAinc, IHinc, RHinc, DHinc]
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
