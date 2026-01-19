import pandapower as pp
import pandapower.networks as nw
import pandapower.shortcircuit as sc
from pandapower.converter.matpower import from_mpc

sys_id = 69
fault_bus_idx = 10

print(f"\n=== PANDAPOWER VALIDATION (IEEE {sys_id}-Bus) ===")
print(f"Fault Location: Bus {fault_bus_idx}")

# 1. Load the Network
if sys_id == 33:
    net = nw.case33bw()
elif sys_id == 69:
    net = from_mpc('case69.m', f_hz=50) # Zhang et al. usually assume 50Hz or generic
else:
    raise ValueError("Only 33 or 69 bus systems supported.")


# 2. Configure External Grid for Short Circuit Analysis
# Pandapower needs to know the short-circuit capacity of the main grid.
# To simulate an "Infinite Bus" (ideal voltage source), we set s_sc_max_mva very high.
# This matches the theoretical assumption that the slack bus voltage doesn't drop.
net.ext_grid['s_sc_max_mva'] = 100  # 10 GVA (Very stiff grid)
net.ext_grid['rx_max'] = 0.1          # R/X ratio (typical value)

# zero sequence network parameters are used for single phase fault only
net.ext_grid['x0x_max'] = 1.0         # X0/X ratio
net.ext_grid['r0x0_max'] = 0.1        # R0/X0 ratio

# 3. Add Zero-Sequence Parameters for Lines (needed for 1ph and 2ph faults)
# Typical assumption: R0 = 3*R1, X0 = 3*X1
net.line['r0_ohm_per_km'] = net.line['r_ohm_per_km'] * 3
net.line['x0_ohm_per_km'] = net.line['x_ohm_per_km'] * 3
net.line['c0_nf_per_km'] = net.line['c_nf_per_km']   # C0 is often similar or smaller, but let's keep it simple


r_f = 0 # seems to have a bug if there is fault impedance. leave them as zero
        # the bug is that for 2ph (line to line), it accidentally takes double
        # the fault impedance.
x_f = 0

# Note: Pandapower uses 0-based indexing for buses.
# If you want Bus 18 (1-based), that is index 17 in pandapower.
pp_bus_idx = fault_bus_idx - 1


# ---------------------------------------------------------
# SCENARIO A: Single Line-to-Ground Fault (1ph)
# ---------------------------------------------------------
# 'fault'='1ph' -> Single phase to ground
# 'ip' -> true calculates peak current
# 'ith' -> true calculates thermal current
sc.calc_sc(net, bus=pp_bus_idx, fault="1ph", case="max", ip=True, r_fault_ohm=r_f, x_fault_ohm=x_f)

# Extract Results
# ikss_ka = Initial Symmetrical Short-Circuit Current (RMS)
ikss_1ph = net.res_bus_sc.at[pp_bus_idx, 'ikss_ka']

# Convert kA to p.u. for comparison with manual code
# Base Current I_base   = S_base / (sqrt(3) * V_base)
base_mva = net.sn_mva  # Usually 10 or 100 MVA
base_kv = net.bus.at[pp_bus_idx, 'vn_kv'] # Nominal voltage at that bus
I_base_kA = base_mva / (1.732 * base_kv)

ikss_1ph_pu = ikss_1ph / I_base_kA

print(f"\n[Scenario 1] Single Line-to-Ground (SLG)")
print(f"  > IEC 60909 Ik'' (RMS): {ikss_1ph:.4f} kA")
print(f"  > Converted to p.u.   : {ikss_1ph_pu:.4f} p.u.")


# ---------------------------------------------------------
# SCENARIO B: Line-to-Line Fault (2ph)
# ---------------------------------------------------------

sc.calc_sc(net, bus=pp_bus_idx, fault="2ph", case="max", ip=True, r_fault_ohm=r_f, x_fault_ohm=x_f)
ikss_2ph = net.res_bus_sc.at[pp_bus_idx, 'ikss_ka']
ikss_2ph_pu = ikss_2ph / I_base_kA

print(f"\n[Scenario 2] Line-to-Line (LL)")
print(f"  > IEC 60909 Ik'' (RMS): {ikss_2ph:.4f} kA")
print(f"  > Converted to p.u.   : {ikss_2ph_pu:.4f} p.u.")
