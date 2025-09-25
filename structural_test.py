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
            stations.append({'eta': eta, 'C': stiff[:36], 'M': massm[:36]})
            continue
    except:
        pass
    i += 1

if len(stations) < 2:
    raise RuntimeError("Not enough stations parsed.")

# Build nodes
blade_length = 117.0  # default; change if you want real blade length
nodes = []
for idx, s in enumerate(stations):
    nodes.append({'id': idx+1, 'x': s['eta'] * blade_length, 'station': s})

elements = []
for idx in range(len(nodes)-1):
    elements.append({'eid': idx+1, 'n1': nodes[idx]['id'], 'n2': nodes[idx+1]['id']})


def extract_section_props(station):
    C = station['C']
    M = station['M']
    EA = C[2*6+2]
    EI_edge = C[3*6+3]
    EI_flap = C[4*6+4]
    GJ = C[5*6+5]
    m = M[0] if len(M)>0 else 0.0
    Xcm = Ycm = 0.0
    if abs(m) > 1e-12 and len(M) >= 12:
        Xcm = M[1*6+5] / m
        Ycm = - M[0*6+5] / m
    return {'EA':EA, 'EI_edge':EI_edge, 'EI_flap':EI_flap, 'GJ':GJ, 'm':m, 'Xcm':Xcm, 'Ycm':Ycm, 'rawM':M}

elem_props = []
for idx, e in enumerate(elements):
    p1 = extract_section_props(nodes[idx]['station'])
    p2 = extract_section_props(nodes[idx+1]['station'])
    prop = {k: 0.5*(p1[k]+p2[k]) for k in p1.keys() if k != 'rawM'}
    # keep averaged raw mass matrix for distribution below (average element endpoints)
    avgM = []
    M1 = p1['rawM']
    M2 = p2['rawM']
    if len(M1) == 36 and len(M2) == 36:
        avgM = [(a+b)*0.5 for a,b in zip(M1, M2)]
    prop['avgM'] = avgM
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
lines.append("")
lines.append("BEGIN BULK")
lines.append("MAT1,1,1.0,1.0,0.3")
if analysis_mode == "SOL103":
    # Add EIGRL card for modal analysis (extract 10 modes from 0 to 200 Hz)
    lines.append("EIGRL,1,0.0,200.0")
lines.append("")

# GRID
for n in nodes:
    lines.append(f"GRID,{n['id']},, {n['x']:.8E}, 0.0, 0.0")
lines.append("")

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
    I12 = 0.0
    J = prop['GJ'] / G if G != 0 else 0.0
    NSM = 0.0
    lines.append(f"PBEAM,{pid},1,{A:.8E},{I1:.8E},{I2:.8E},{I12:.8E},{J:.8E},{NSM:.8E}")
    lines.append("")

    # Mass: use full 6x6 from station, average between the two ends
    M1 = nodes[e['n1']-1]['station']['M']
    M2 = nodes[e['n2']-1]['station']['M']
    avgM = [(a+b)*0.5 for a,b in zip(M1,M2)]
    m21 = lower_triangle_21(avgM)

    eid_conm1 = 1000 + idx
    g = e['n1']  # attach to first node (could also distribute)
    cid = 0
    mscale = 1.0
    # lines.append(f"CONM1,{eid_conm1},{g},{cid},{mscale}")
    # # Write continuation lines (up to 8 values per line)
    # for k in range(0, len(m21), 8):
    #     chunk = m21[k:k+8]
    #     vals = ",".join(f"{v:.6E}" for v in chunk)
    #     lines.append("+" + vals)
    # lines.append("")
    # First line (header)
    header = f"CONM1,{eid_conm1},{g},{cid},{mscale}"
    # Flatten the rest into one big comma-separated string
    body = ",".join(f"{v:.6E}" for v in m21)
    full = header + "," + body

    # Now wrap at 80 characters max per line
    while len(full) > 80:
        cut = full.rfind(",", 0, 80)  # split at last comma before 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        if cut == -1:
            cut = 80
        lines.append(full[:cut])
        full = full[cut+1:]  # skip the comma
        full = "," + full    # re-add leading comma for continuation

    lines.append(full)
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

if analysis_mode == "SOL101":
    print(f"Wrote SOL101 static deck with CONM1 to: {out_dat}")
    print("Tip node id:", tip_node)
    print(f"Applied tip force at GRID {tip_node}: F = {tip_force_value} in Z direction (0,0,1)")
elif analysis_mode == "SOL103":
    print(f"Wrote SOL103 modal deck with CONM1 to: {out_dat}")