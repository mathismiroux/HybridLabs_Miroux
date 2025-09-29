analysis_mode = "SOL103"  # Change to "SOL103" for modal analysis
from pathlib import Path
import math

beam_input = Path("IEA-15-240-RWT_BeamDyn_blade.dat")
out_dat = Path(f"beamdyn_{analysis_mode}_conm1.dat")



def floats_in_line(s):
    toks = s.replace(',', ' ').split()
    out = []
    for t in toks:
        try:
            out.append(float(t))
        except:
            pass
    return out

# Parse stations
with open(beam_input, "r", errors="ignore") as f:
    raw_lines = [l.rstrip() for l in f.readlines()]

stations = []
i = 0
nlines = len(raw_lines)
while i < nlines:
    line = raw_lines[i].strip()
    try:
        val = float(line)
        if 0.0 <= val <= 1.0:
            eta = val
            i += 1
            stiff = []
            while i < nlines and len(stiff) < 36:
                stiff += floats_in_line(raw_lines[i])
                i += 1
            while i < nlines and raw_lines[i].strip() == "":
                i += 1
            massm = []
            while i < nlines and len(massm) < 36:
                massm += floats_in_line(raw_lines[i])
                i += 1
            C = stiff[:36]
            M = massm[:36]
            # Check diagonal values for negatives
            diag_indices = [0, 7, 14, 21, 28, 35]
            neg_diag_C = [C[i] for i in diag_indices if i < len(C) and C[i] < 0]
            neg_diag_M = [M[i] for i in diag_indices if i < len(M) and M[i] < 0]
            if neg_diag_C:
                print(f"Negative diagonal in stiffness matrix at eta={eta:.6f}: {neg_diag_C}\nC={C}")
            if neg_diag_M:
                print(f"Negative diagonal in mass matrix at eta={eta:.6f}: {neg_diag_M}\nM={M}")
            stations.append({'eta': eta, 'C': C, 'M': M})
            continue
    except:
        pass
    i += 1


if len(stations) < 2:
    raise RuntimeError("Not enough stations parsed.")
blade_length = 117.0  # default; change if you want real blade length

# Calculate and print blade mass using original station data

station_mass = 0.0
for i in range(len(stations)-1):
    m11 = (stations[i]['M'][0]+stations[i+1]['M'][0])/2 if len(stations[i]['M']) > 0 else 0.0
    seg_length = (stations[i+1]['eta'] - stations[i]['eta']) * blade_length
    station_mass += m11 * seg_length
    print(f"Element {i+1}: eta={stations[i]['eta']:.6f}, m11={m11:.6f}, seg_length={seg_length:.3f}, seg_mass={m11*seg_length:.3f}")

print(f"Blade mass from original station data: {station_mass:.3f} kg")
print("Element eta and m11 (from station data):")
for i in range(len(stations)):
    eta = stations[i]['eta']
    m11 = stations[i]['M'][0] if len(stations[i]['M']) > 0 else 0.0
    print(f"Element {i+1}: eta={eta:.6f}, m11={m11:.6f}")

# Build nodes
nodes = []
for idx, s in enumerate(stations):
    nodes.append({'id': idx+1, 'x': s['eta'] * blade_length, 'station': s})

elements = []
for idx in range(len(nodes)-1):
    elements.append({'eid': idx+1, 'n1': nodes[idx]['id'], 'n2': nodes[idx+1]['id']})

# # Interpolate additional stations
# def interpolate_station(s1, s2, eta):
#     # Linear interpolation for eta, C, M
#     frac = (eta - s1['eta']) / (s2['eta'] - s1['eta'])
#     C = [(1-frac)*c1 + frac*c2 for c1, c2 in zip(s1['C'], s2['C'])]
#     M = [(1-frac)*m1 + frac*m2 for m1, m2 in zip(s1['M'], s2['M'])]
#     return {'eta': eta, 'C': C, 'M': M}

# # Number of additional stations between each pair
# n_interp = 2  # Change this for finer mesh
# refined_stations = []
# for i in range(len(stations)-1):
#     s1 = stations[i]
#     s2 = stations[i+1]
#     refined_stations.append(s1)
#     for k in range(1, n_interp+1):
#         eta_new = s1['eta'] + (s2['eta'] - s1['eta']) * k/(n_interp+1)
#         refined_stations.append(interpolate_station(s1, s2, eta_new))
# refined_stations.append(stations[-1])

def extract_section_props(station):
    C = station['C']
    M = station['M']
    EA = C[2*6+2]
    EI_edge = C[3*6+3]
    EI_flap = C[4*6+4]
    GJ = C[5*6+5]
    # Coupled terms
    edge_shear_torsional = C[1*6+5]
    flap_shear_torsional = C[0*6+5]
    axial_edge_bending = C[4*6+3]
    axial_flap_bending = C[4*6+2]
    flap_edge_shear = C[0*6+1]
    flap_edge_bending = C[4*6+3]
    m = M[0] if len(M)>0 else 0.0
    Xcm = Ycm = 0.0
    if abs(m) > 1e-12 and len(M) >= 12:
        Xcm = M[1*6+5] / m
        Ycm = - M[0*6+5] / m
    return {
        'EA':EA,
        'EI_edge':EI_edge,
        'EI_flap':EI_flap,
        'GJ':GJ,
        'edge_shear_torsional': edge_shear_torsional,
        'flap_shear_torsional': flap_shear_torsional,
        'axial_edge_bending': axial_edge_bending,
        'axial_flap_bending': axial_flap_bending,
        'flap_edge_shear': flap_edge_shear,
        'flap_edge_bending': flap_edge_bending,
        'm':m,
        'Xcm':Xcm,
        'Ycm':Ycm,
        'rawM':M
    }

elem_props = []
for idx, e in enumerate(elements):
    p1 = extract_section_props(nodes[idx]['station'])
    p2 = extract_section_props(nodes[idx+1]['station'])
    prop = {k: 0.5*(p1[k]+p2[k]) for k in p1.keys() if k != 'rawM'}
    # Manual input/adjustment for coupled terms (example: overwrite with user values)
    # prop['edge_shear_torsional'] = YOUR_VALUE_HERE
    # prop['flap_shear_torsional'] = YOUR_VALUE_HERE
    # prop['axial_edge_bending'] = YOUR_VALUE_HERE
    # prop['axial_flap_bending'] = YOUR_VALUE_HERE
    # prop['flap_edge_shear'] = YOUR_VALUE_HERE
    # prop['flap_edge_bending'] = YOUR_VALUE_HERE
    avgM = []
    M1 = p1['rawM']
    M2 = p2['rawM']
    span_length = nodes[idx+1]['x'] - nodes[idx]['x']
    if len(M1) == 36 and len(M2) == 36:
        avgM = [((a+b)*0.5) * span_length for a,b in zip(M1, M2)]
    prop['avgM'] = avgM
    # print(f"Element {idx+1} coupled terms:")
    # print(f"  edge_shear_torsional: {prop['edge_shear_torsional']}")
    # print(f"  flap_shear_torsional: {prop['flap_shear_torsional']}")
    # print(f"  axial_edge_bending: {prop['axial_edge_bending']}")
    # print(f"  axial_flap_bending: {prop['axial_flap_bending']}")
    # print(f"  flap_edge_shear: {prop['flap_edge_shear']}")
    # print(f"  flap_edge_bending: {prop['flap_edge_bending']}")
    elem_props.append(prop)

# utility: extract 21 lower-triangular terms from 6x6
def lower_triangle_21(matrix6):
    out = []
    for i in range(6):
        for j in range(i+1):
            out.append(matrix6[i*6 + j])
    return out

# user parameters
tip_force_value = -100000.0
tip_force_component = (0.0, 0.0, tip_force_value)
tip_node = nodes[-1]['id']

# Case control and bulk
lines = []
if analysis_mode == "SOL101":
    lines.append("SOL 101")
    lines.append("CEND")
    lines.append("TITLE = BeamDyn -> Nastran SOL 101 static displacement example with CONM1")
    lines.append("SUBCASE 1")
    lines.append("    LOAD = 100")
    lines.append("    SPC = 1")
    lines.append("    DISPLACEMENT(PLOT,PRINT) = ALL")
    lines.append("    STRESS(PLOT) = ALL")
    lines.append("    STRAIN(PLOT) = ALL")
    lines.append("PARAM,PRTMAXIM,YES")  
elif analysis_mode == "SOL103":
    lines.append("SOL 103")
    lines.append("CEND")
    lines.append("TITLE = BeamDyn -> Nastran SOL 103 modal analysis example with CONM1")
    lines.append("SUBCASE 1")
    lines.append("    SPC = 1")
    lines.append("    DISPLACEMENT(PLOT,PRINT) = ALL")
    lines.append("    METHOD = 1")
    lines.append("    ELSUM = ALL")
    # lines.append("    ELFORCE(PLOT,PRINT,HDF5) = ALL")
    # lines.append("    STRESS(PLOT,PRINT,HDF5) = ALL")

lines.append("")
lines.append("BEGIN BULK")
lines.append("PARAM,POST,-2")
lines.append("MDLPRM, HDF5,1")
lines.append("PARAM,PRTMAX,YES")
lines.append("MAT1,1,1.0,0.3")
if analysis_mode == "SOL103":
    # Add EIGRL card for modal analysis (extract 10 modes from 0 to 200 Hz)
    lines.append("EIGRL,1,0.0,200.0")
lines.append("")

# GRID
for n in nodes:
    lines.append(f"GRID,{n['id']},, {n['x']:.8E}, 0.0, 0.0")
lines.append("")
print()
# PBAR, CBEAM, PBEAM, CONM1 (mass distributed)
pid_start = 1
conm_eid = 1000
for idx, e in enumerate(elements):
    pid = pid_start + idx
    prop = elem_props[idx]
    lines.append(f"$ Element {e['eid']} between GRID {e['n1']} and {e['n2']}")
    # orientation vector: check if (0,0,1) aligns PBEAM I1/I2 with edge/flap. If modes look swapped, change orientation or swap I1/I2.
    lines.append(f"CBEAM,{e['eid']},{pid},{e['n1']},{e['n2']},0.0,1.0,0.0")
    E = 1.0
    G = 1.0
    A = prop['EA'] / E if E != 0 else 0.0
    I1 = prop['EI_edge'] / E if E != 0 else 0.0
    I2 = prop['EI_flap'] / E if E != 0 else 0.0
    # Coupled terms
    I12 = prop.get('flap_edge_bending', 0.0) / E if E != 0 else 0.0  # Example mapping
    # I12 = 0.0  # Example mapping
    print(f'{I12:.3E}')
    # You can add more coupled terms in extended fields or comments
    J = prop['GJ'] / G if G != 0 else 0.0
    E = 1.0
    G = 1.0
    A = prop['EA'] / E if E != 0 else 0.0
    I1 = prop['EI_edge'] / E if E != 0 else 0.0
    I2 = prop['EI_flap'] / E if E != 0 else 0.0
    # Coupled terms
    I12 = prop.get('flap_edge_bending', 0.0) / E if E != 0 else 0.0  # Example mapping
    J = prop['GJ'] / G if G != 0 else 0.0
    NSM = 0.0
    # Check for invalid values
    if A <= 0:
        print(f"WARNING: Element {e['eid']} has invalid area A={A:.3E}")
    if I1 <= 0:
        print(f"WARNING: Element {e['eid']} has invalid I1={I1:.3E}")
    if I2 <= 0:
        print(f"WARNING: Element {e['eid']} has invalid I2={I2:.3E}")
    if J <= 0:
        print(f"WARNING: Element {e['eid']} has invalid J={J:.3E}")
    if I1 * I2 - I12**2 <= 0:
        print(f"WARNING: Element {e['eid']} fails inertia matrix check: I1*I2 - I12^2 = {I1*I2 - I12**2:.3E}")
    # Add comments for coupled terms
    lines.append(f"$ Coupled terms: edge_shear_torsional={prop.get('edge_shear_torsional', 0.0):.3E}, flap_shear_torsional={prop.get('flap_shear_torsional', 0.0):.3E}, axial_edge_bending={prop.get('axial_edge_bending', 0.0):.3E}, axial_flap_bending={prop.get('axial_flap_bending', 0.0):.3E}, flap_edge_shear={prop.get('flap_edge_shear', 0.0):.3E}, flap_edge_bending={prop.get('flap_edge_bending', 0.0):.3E}")
    lines.append(f"PBEAM,{pid},1,{A:.6E},{I1:.6E},{I2:.6E},{I12:.3E},{J:.6E},{NSM:.1E}")
    lines.append("")

    # Mass: use full 6x6 from station, average between the two ends
    # M1 = nodes[e['n1']-1]['station']['M']
    # M2 = nodes[e['n2']-1]['station']['M']
    # span_length = nodes[idx+1]['x'] - nodes[idx]['x']

    # # avgM = [(a+b)*0.5 for a,b in zip(M1,M2)]
    # avgM = [(a+b)*0.5*span_length for a,b in zip(M1,M2)]
    m21 = lower_triangle_21(prop["avgM"])
    # print(m21[0])

    eid_conm1 = 1000 + idx
    g = e['n1']  # attach to first node (could also distribute)
    cid = 0
    # Write CONM1 card in Nastran format: first line 9 fields, then 8 fields per line
    # Prepare all fields
    conm1_fields = [f"CONM1", str(eid_conm1), str(g), str(cid)] + [f"{v:.3E}" for v in m21]
    # First line: 9 fields
    first_line = ",".join(conm1_fields[:9])
    lines.append(first_line)
    # Remaining lines: 8 fields per line
    idx = 9
    while idx < len(conm1_fields):
        next_line = ",".join(conm1_fields[idx:idx+8])
        lines.append(f',{next_line}')
        idx += 8
    lines.append("")

# SPC
lines.append("SPC1,1,123456,1")
# lines.append("SPC,1,1")
lines.append("")

# Load (only for SOL101)
if analysis_mode == "SOL101":
    # Apply tip force in Y direction using direction vector (0,1,0)
    lines.append(f"FORCE,100,{tip_node},0,{abs(tip_force_value):.6e},0.,0.,1.")
    lines.append("")

lines.append("ENDDATA")

with open(out_dat, "w") as f:
    f.write("\n".join(lines))

# # Calculate and print blade mass
# blade_mass = sum([prop['avgM'][0] for prop in elem_props if len(prop['avgM']) > 0])
# print(f"Total blade mass: {blade_mass:.3f} kg")

if analysis_mode == "SOL101":
    print(f"Wrote SOL101 static deck with CONM1 to: {out_dat}")
    print("Tip node id:", tip_node)
    print(f"Applied tip force at GRID {tip_node}: F = {tip_force_value} in Z direction (0,0,1)")
elif analysis_mode == "SOL103":
    print(f"Wrote SOL103 modal deck with CONM1 to: {out_dat}")