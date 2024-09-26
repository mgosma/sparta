input_template = """
################################################################################
# 2d axi flow around a capsule ()
# initial simulation run
# Note:
#  - The "comm/sort” option to the “global” command is used to match MPI runs.
#  - The “twopass” option is used to match Kokkos runs.
# The "comm/sort" and "twopass" options should not be used for production runs.
################################################################################

###################################
# Trajectory inputs
###################################
variable            V equal {vstream}
variable            nden equal {nrho}
variable            temp equal {temp}

###################################
# Simulation initialization standards
###################################
variable            ppc equal 1
variable            cpmfp equal 10

###################################
# Parameter calculations
###################################
variable            mfp equal 1/(sqrt(2)*3.14159*(4.287e-10)^2*${nden}*(273/${temp})^(0.216))

variable            xmin equal -3.0
variable            xmax equal 5.0
variable            ymin equal 0.0
variable            ymax equal 5.0
variable            zmin equal -0.5
variable            zmax equal 0.5

variable            xncells equal floor((${xmax}-${xmin})/${mfp}*${cpmfp})
variable            yncells equal floor((${ymax}-${ymin})/${mfp}*${cpmfp})
#variable            zncells equal floor((${zmax}-${zmin})/${mfp}*${cpmfp})

variable            Fnum equal  ${nden}*(${xmax}-${xmin})*(${ymax}-${ymin})/${ppc}/${xncells}/${yncells}

variable            tstep equal (-${xmin}+${xmax})/${V}/${xncells}/30

###################################
# Print variable values to log file
###################################
print               " Velocity  = ${V}"
print               " Density  = ${nden}"
print               " Temp  = ${temp}"
print               " mean free path  = ${mfp}"
print               " cells per free stream mean free path = ${cpmfp}"
print               " x-min = ${xmin}"
print               " x-max = ${xmax}"
print               " y-min = ${ymin}"
print               " y-max = ${ymax}"
#print               " z-min = ${zmin}"
#print               " z-max = ${zmax}"
print               " x-cells = ${xncells}"
print               " y-cells = ${yncells}"
#print               " z-cells = ${zncells}"
print               " Simulation Ratio = ${Fnum}"
print               " Timestep = ${tstep}"

###################################
# Simulation parameters
###################################
seed	    	    847384
dimension   	    2
global		    nrho ${nden} fnum ${Fnum} gridcut 0.01
global		    comm/sort yes
timestep            ${tstep}

###################################
# Grid generation
###################################
boundary	    o ao p
create_box          ${xmin} ${xmax} ${ymin} ${ymax} ${zmin} ${zmax}
create_grid 	    ${xncells} ${yncells} 1 block * * * 
balance_grid        rcb cell
global          weight cell radius

#####################################
# Gas/Collision Model Specification #
#####################################
species             Titan.species N2 CH4 CH3 CH2 CH C2 H2 CN NH HCN N C H Ar e N2+ CN+ N+ C+ vibfile Titan.vib rotfile Titan.rot elecfile Titan.elec
mixture             all vstream ${V} 0 0 temp ${temp}
mixture             all N2 frac 0.978
mixture             all CH4 frac 0.022

collide             vss all Titan.vss
collide_modify      vremax 1000 yes nearcp yes 50 vibrate discrete
fix                 vib vibmode

compute tcell thermal/grid all all temp press
react        tce Titan.tce
react_modify rboost 1.0
react_modify partial_energy no c_tcell[1] 
react_modify partition yes

#####################################################
# Surface generation and collision specification
#####################################################
read_surf	    data.dragonfly group 1

compute             1 surf all all etot
fix                 1 ave/surf all 1 250 250 c_1[1] ave one
fix                 2 surf/temp all 250 f_1 100 0.9 temperature
surf_collide        1 diffuse s_temperature 1.0
surf_modify         1 collide 1

###################################
# Boundary conditions
###################################
fix                 in emit/face all xlo yhi twopass

###################################
# Initialize simulation
###################################
create_particles    all n 0 twopass
balance_grid        rcb part

###################################
# Unsteady Output
###################################
compute             2 grid all all nrho temp usq vsq wsq
fix                 3 ave/grid all 1 100 100 c_2[*]

compute             2b lambda/grid f_3[1] f_3[2] N2 kall
compute             MFP reduce min c_2b[1]
compute             Kn reduce min c_2b[2]

compute             tstep dt/grid all 1.0 1.0 &
                    c_2b[1] f_3[2] f_3[3] f_3[4] f_3[5]

fix                 DT dt/reset 1000 c_tstep 0.1 0

fix                 10 adapt 1000 all refine coarsen value c_2b[2] 0.8 1.6 &
                    combine min thresh less more cells 2 2 1
fix                 load balance 1000 1.1 rcb part

stats_style         step cpu np ncoll nreact nscoll c_MFP c_Kn ngrid maxlevel dt f_DT[*]
stats_modify        flush yes
stats               100

restart             10000 restart.sparta

run                 50000 pre no post no every 1000 "react_modify partial_energy no c_tcell[1]"
...
"""
import os

# Define test parameters
test_cases = [
    {
        "vstream": 7384.53,
        "temp": 165.93,
        "nrho": 2.58e18,
        "chemistry_model": "neutral"
    },
    {
        "vstream": 7392.23,
        "temp": 165.93,
        "nrho": 1.24e19,
        "chemistry_model": "ambipolar"
    },
    # Add more cases as needed
]

# Base directory for simulations
base_dir = "Dragonfly_Simulations"

def create_simulation_case(case, case_id, base_dir):
    # Create a unique directory for this case
    case_dir = os.path.join(base_dir, f"case_{case_id}")
    os.makedirs(case_dir, exist_ok=True)
    
    # Modify the input template with specific case values
    input_content = input_template.format(
        vstream=case["vstream"],
        temp=case["temp"],
        nrho=case["nrho"],
        chemistry_model=case["chemistry_model"]
    )
    
    # Write the input file to the directory
    input_file_path = os.path.join(case_dir, "input_file.txt")
    with open(input_file_path, 'w') as f:
        f.write(input_content)
    
    print(f"Created simulation case {case_id} in {case_dir}")

# Generate cases
for i, case in enumerate(test_cases):
    create_simulation_case(case, i, base_dir)

