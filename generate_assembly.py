"""
Femoral Implant Multi-Block Assembly Generator
Patient-Specific FEA Post-Processing for Stress Shielding Analysis

Author: Jeffrey Husch, Senior Biomedical Engineer
Application: ParaView Orthopedic FEA Workflow
"""
import math
import os

# ============================================================================
# GRID PARAMETERS - High Resolution for FEA Post-Processing
# ============================================================================
nx, ny, nz = 45, 45, 100
spacing = 1.0  # mm
origin = (-22, -22, 0)

# Femur geometry parameters (proximal femur)
femur_length = 100  # mm (proximal portion)
femur_shaft_radius = 15  # mm
femur_head_radius = 24  # mm
femur_neck_angle = 125  # degrees (inclination)
cortical_thickness = 4  # mm

# Implant geometry (cemented stem)
stem_length = 85  # mm
stem_proximal_radius = 12  # mm
stem_distal_radius = 6  # mm
stem_taper = (stem_proximal_radius - stem_distal_radius) / stem_length

# Material properties (MPa)
CORTICAL_E = 17000  # Cortical bone Young's modulus
CANCELLOUS_E = 500  # Cancellous bone Young's modulus
TITANIUM_E = 110000  # Ti-6Al-4V implant
COBALT_CHROME_E = 210000  # CoCr femoral head

# Loading conditions
BODY_WEIGHT = 750  # N (75 kg patient)
HIP_LOAD_FACTOR = 2.5  # Peak load during gait
APPLIED_LOAD = BODY_WEIGHT * HIP_LOAD_FACTOR  # 1875 N

def distance_3d(x, y, z, cx, cy, cz):
    return math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)

def femur_geometry(x, y, z):
    """
    Define femoral geometry with anatomical features.
    Returns: (in_cortical, in_cancellous, in_marrow, region_label)
    """
    # Shaft region (z < 60)
    if z < 60:
        r = math.sqrt(x**2 + y**2)
        if r < femur_shaft_radius - cortical_thickness:
            if r < femur_shaft_radius - cortical_thickness - 3:
                return False, False, True, 3  # Marrow cavity
            return False, True, False, 2  # Cancellous
        elif r < femur_shaft_radius:
            return True, False, False, 1  # Cortical
        return False, False, False, 0  # Outside
    
    # Neck region (60 <= z < 85)
    elif z < 85:
        # Neck axis offset and angle
        neck_progress = (z - 60) / 25
        neck_offset_x = neck_progress * 15  # Lateral offset
        neck_offset_z = neck_progress * 10
        
        neck_center_x = neck_offset_x
        neck_center_y = 0
        
        r_neck = math.sqrt((x - neck_center_x)**2 + y**2)
        neck_radius = 14 - neck_progress * 2  # Tapers
        
        if r_neck < neck_radius - cortical_thickness:
            return False, True, False, 2
        elif r_neck < neck_radius:
            return True, False, False, 1
        return False, False, False, 0
    
    # Head region (z >= 85)
    else:
        head_center = (18, 0, 90)
        r_head = distance_3d(x, y, z, *head_center)
        
        if r_head < femur_head_radius - 3:
            return False, True, False, 2  # Cancellous head
        elif r_head < femur_head_radius:
            return True, False, False, 1  # Subchondral bone
        return False, False, False, 0

def implant_geometry(x, y, z):
    """
    Define hip prosthesis geometry (cemented femoral stem).
    Returns: (in_stem, in_head, in_cement, component_label)
    """
    # Stem region (0 < z < 85)
    if 5 < z < 85:
        stem_progress = z / stem_length
        current_radius = stem_proximal_radius - stem_taper * z
        
        r = math.sqrt(x**2 + y**2)
        
        # Stem body
        if r < current_radius:
            return True, False, False, 1  # Stem
        # Cement mantle (2mm around stem)
        elif r < current_radius + 2:
            return False, False, True, 3  # Cement
        return False, False, False, 0
    
    # Neck and head region
    elif z >= 85:
        neck_center_x = 15
        head_center = (18, 0, 95)
        
        # Femoral head (CoCr)
        r_head = distance_3d(x, y, z, *head_center)
        if r_head < 14:  # Prosthetic head radius
            return False, True, False, 2  # Head
        
        # Neck taper
        if z < 95:
            neck_r = math.sqrt((x - neck_center_x)**2 + y**2)
            if neck_r < 6:
                return True, False, False, 1  # Neck (part of stem)
        
        return False, False, False, 0
    
    return False, False, False, 0

def calculate_stress_field(x, y, z, in_bone, in_implant, has_implant=True):
    """
    Calculate von Mises stress distribution based on FEA principles.
    Implements analytical approximations for demonstration.
    """
    if not in_bone and not in_implant:
        return 0, 0, 0, 0, 0, 0
    
    r = math.sqrt(x**2 + y**2)
    
    # Axial stress from body weight
    if z < 60:  # Shaft
        area_bone = math.pi * (femur_shaft_radius**2 - (femur_shaft_radius - cortical_thickness)**2)
        if has_implant:
            area_implant = math.pi * (stem_proximal_radius - stem_taper * z)**2
            # Load sharing between implant and bone
            stiffness_ratio = TITANIUM_E / CORTICAL_E
            implant_load_fraction = (area_implant * stiffness_ratio) / (area_implant * stiffness_ratio + area_bone)
            bone_load_fraction = 1 - implant_load_fraction
            
            if in_implant:
                axial_stress = (APPLIED_LOAD * implant_load_fraction) / area_implant
            else:
                axial_stress = (APPLIED_LOAD * bone_load_fraction) / area_bone
        else:
            axial_stress = APPLIED_LOAD / area_bone
    else:
        # Neck/head region - bending stresses
        moment_arm = abs(x - 18)  # Distance from load axis
        bending_moment = APPLIED_LOAD * moment_arm
        I = math.pi * femur_shaft_radius**4 / 4
        axial_stress = bending_moment * r / I * 1000  # Scale factor
    
    # Add bending component
    bending_stress = abs(x) * APPLIED_LOAD / 5000
    
    # Hoop stress (from press-fit)
    if has_implant and in_bone and r < femur_shaft_radius:
        hoop_stress = 15 * (1 - r / femur_shaft_radius)  # Contact pressure effect
    else:
        hoop_stress = 0
    
    # Principal stresses (simplified)
    sigma_1 = axial_stress + bending_stress
    sigma_2 = hoop_stress
    sigma_3 = -0.3 * axial_stress  # Poisson effect
    
    # von Mises stress
    von_mises = math.sqrt(0.5 * ((sigma_1 - sigma_2)**2 + (sigma_2 - sigma_3)**2 + (sigma_3 - sigma_1)**2))
    
    # Hydrostatic stress
    hydrostatic = (sigma_1 + sigma_2 + sigma_3) / 3
    
    # Contact pressure (at interface)
    if has_implant and in_bone:
        stem_r = stem_proximal_radius - stem_taper * z
        dist_to_stem = abs(r - stem_r)
        if dist_to_stem < 3:
            contact_pressure = 8 * math.exp(-dist_to_stem)
        else:
            contact_pressure = 0
    else:
        contact_pressure = 0
    
    return von_mises, sigma_1, sigma_2, sigma_3, hydrostatic, contact_pressure

def calculate_strain_energy_density(von_mises, E_modulus):
    """
    Calculate Strain Energy Density (SED) for bone remodeling prediction.
    SED = sigma^2 / (2 * E)
    """
    if E_modulus > 0:
        sed = (von_mises ** 2) / (2 * E_modulus)
    else:
        sed = 0
    return sed

def write_vtk_structured_grid(filename, data_dict, title):
    """Write VTK structured grid file with multiple scalar fields."""
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        f.write(f"ORIGIN {origin[0]} {origin[1]} {origin[2]}\n")
        f.write(f"SPACING {spacing} {spacing} {spacing}\n")
        f.write(f"POINT_DATA {nx*ny*nz}\n")
        
        for name, values in data_dict.items():
            if name == "Displacement":
                f.write(f"VECTORS {name} float\n")
                for v in values:
                    f.write(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            else:
                f.write(f"SCALARS {name} float\nLOOKUP_TABLE default\n")
                for v in values:
                    f.write(f"{v:.6f}\n")

# ============================================================================
# GENERATE BONE MODEL (WITHOUT IMPLANT - INTACT FEMUR)
# ============================================================================
print("Generating intact femur model (pre-operative)...")

bone_intact = {"tissue_type": [], "von_mises_stress": [], "principal_1": [],
               "principal_3": [], "hydrostatic_stress": [], "SED": [],
               "E_modulus": [], "Displacement": []}

for k in range(nz):
    z = origin[2] + k * spacing
    for j in range(ny):
        y = origin[1] + j * spacing
        for i in range(nx):
            x = origin[0] + i * spacing
            
            in_cort, in_canc, in_marrow, region = femur_geometry(x, y, z)
            
            if in_cort:
                tissue = 1
                E = CORTICAL_E
            elif in_canc:
                tissue = 2
                E = CANCELLOUS_E
            elif in_marrow:
                tissue = 3
                E = 10  # Marrow
            else:
                tissue = 0
                E = 0
            
            bone_intact["tissue_type"].append(tissue)
            bone_intact["E_modulus"].append(E)
            
            if tissue > 0:
                vm, s1, s2, s3, hyd, cp = calculate_stress_field(x, y, z, True, False, has_implant=False)
                sed = calculate_strain_energy_density(vm, E) if E > 0 else 0
                
                # Displacement (scaled)
                disp_z = -vm / E * 0.1 if E > 0 else 0
                disp_x = -0.3 * disp_z * x / 20
                disp_y = -0.3 * disp_z * y / 20
            else:
                vm, s1, s3, hyd, sed = 0, 0, 0, 0, 0
                disp_x, disp_y, disp_z = 0, 0, 0
            
            bone_intact["von_mises_stress"].append(vm)
            bone_intact["principal_1"].append(s1)
            bone_intact["principal_3"].append(s3)
            bone_intact["hydrostatic_stress"].append(hyd)
            bone_intact["SED"].append(sed)
            bone_intact["Displacement"].append((disp_x, disp_y, disp_z))

write_vtk_structured_grid("femur_intact_stress.vtk", bone_intact, 
                          "Intact Femur Stress Analysis - Pre-Operative Baseline")
print("  Created: femur_intact_stress.vtk")

# ============================================================================
# GENERATE IMPLANTED MODEL (CORTICAL BONE BLOCK)
# ============================================================================
print("Generating implanted model - Cortical Bone block...")

cortical_data = {"tissue_type": [], "von_mises_stress": [], "stress_ratio": [],
                 "SED": [], "stress_shielding_index": [], "contact_pressure": [],
                 "bone_remodeling_stimulus": [], "Displacement": []}

# Reference SED from intact bone (for stress shielding calculation)
intact_sed_values = bone_intact["SED"]

idx = 0
for k in range(nz):
    z = origin[2] + k * spacing
    for j in range(ny):
        y = origin[1] + j * spacing
        for i in range(nx):
            x = origin[0] + i * spacing
            
            in_cort, in_canc, in_marrow, region = femur_geometry(x, y, z)
            in_stem, in_head, in_cement, comp = implant_geometry(x, y, z)
            
            # Bone with cavity for implant
            in_bone = (in_cort or in_canc) and not (in_stem or in_head or in_cement)
            
            if in_bone:
                if in_cort:
                    tissue = 1
                    E = CORTICAL_E
                else:
                    tissue = 2
                    E = CANCELLOUS_E
                
                vm, s1, s2, s3, hyd, cp = calculate_stress_field(x, y, z, True, False, has_implant=True)
                sed = calculate_strain_energy_density(vm, E)
                
                # Stress shielding index (ratio to intact)
                intact_sed = intact_sed_values[idx] if intact_sed_values[idx] > 0 else 0.001
                stress_ratio = sed / intact_sed if intact_sed > 0 else 1.0
                stress_shielding = max(0, 1 - stress_ratio)  # 0 = no shielding, 1 = complete shielding
                
                # Bone remodeling stimulus (Huiskes' lazy zone model)
                # Values < 0.75 * homeostatic -> resorption
                # Values > 1.25 * homeostatic -> formation
                sed_homeostatic = 0.02  # Reference SED (MPa)
                if sed < 0.75 * sed_homeostatic:
                    remodeling = -1  # Resorption risk
                elif sed > 1.25 * sed_homeostatic:
                    remodeling = 1   # Formation expected
                else:
                    remodeling = 0   # Equilibrium (lazy zone)
                
                disp_z = -vm / E * 0.15
                disp_x = -0.3 * disp_z * x / 20
                disp_y = -0.3 * disp_z * y / 20
                
            else:
                tissue, vm, sed, stress_ratio, stress_shielding = 0, 0, 0, 0, 0
                cp, remodeling = 0, 0
                disp_x, disp_y, disp_z = 0, 0, 0
            
            cortical_data["tissue_type"].append(tissue)
            cortical_data["von_mises_stress"].append(vm)
            cortical_data["stress_ratio"].append(stress_ratio)
            cortical_data["SED"].append(sed)
            cortical_data["stress_shielding_index"].append(stress_shielding)
            cortical_data["contact_pressure"].append(cp)
            cortical_data["bone_remodeling_stimulus"].append(remodeling)
            cortical_data["Displacement"].append((disp_x, disp_y, disp_z))
            
            idx += 1

write_vtk_structured_grid("cortical_bone_implanted.vtk", cortical_data,
                          "Cortical Bone with Femoral Implant - Stress Shielding Analysis")
print("  Created: cortical_bone_implanted.vtk")

# ============================================================================
# GENERATE PROSTHESIS BLOCK
# ============================================================================
print("Generating Prosthesis block...")

prosthesis_data = {"component_type": [], "von_mises_stress": [], "safety_factor": [],
                   "interface_stress": [], "E_modulus": [], "Displacement": []}

for k in range(nz):
    z = origin[2] + k * spacing
    for j in range(ny):
        y = origin[1] + j * spacing
        for i in range(nx):
            x = origin[0] + i * spacing
            
            in_stem, in_head, in_cement, comp = implant_geometry(x, y, z)
            
            if in_stem:
                component = 1  # Ti-6Al-4V stem
                E = TITANIUM_E
                yield_strength = 900  # MPa
            elif in_head:
                component = 2  # CoCr head
                E = COBALT_CHROME_E
                yield_strength = 500
            elif in_cement:
                component = 3  # PMMA cement
                E = 2500
                yield_strength = 80
            else:
                component = 0
                E = 0
                yield_strength = 1
            
            if component > 0:
                vm, s1, s2, s3, hyd, cp = calculate_stress_field(x, y, z, False, True, has_implant=True)
                safety_factor = yield_strength / max(vm, 0.1)
                
                # Interface stress (at cement-bone or stem-cement interface)
                r = math.sqrt(x**2 + y**2)
                stem_r = stem_proximal_radius - stem_taper * z
                interface_dist = abs(r - stem_r)
                if interface_dist < 2:
                    interface_stress = vm * 0.8 * (1 - interface_dist / 2)
                else:
                    interface_stress = 0
                
                disp_z = -vm / E * 0.12
                disp_x = -0.3 * disp_z * x / 20
                disp_y = -0.3 * disp_z * y / 20
            else:
                vm, safety_factor, interface_stress = 0, 0, 0
                disp_x, disp_y, disp_z = 0, 0, 0
            
            prosthesis_data["component_type"].append(component)
            prosthesis_data["von_mises_stress"].append(vm)
            prosthesis_data["safety_factor"].append(min(safety_factor, 10))
            prosthesis_data["interface_stress"].append(interface_stress)
            prosthesis_data["E_modulus"].append(E)
            prosthesis_data["Displacement"].append((disp_x, disp_y, disp_z))

write_vtk_structured_grid("prosthesis.vtk", prosthesis_data,
                          "Femoral Prosthesis - Stress and Safety Factor Analysis")
print("  Created: prosthesis.vtk")

# ============================================================================
# GENERATE MULTI-BLOCK VTM FILE
# ============================================================================
print("Generating Multi-Block VTM assembly...")

vtm_content = '''<?xml version="1.0"?>
<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">
  <vtkMultiBlockDataSet>
    <Block index="0" name="Cortical_Bone">
      <DataSet index="0" file="cortical_bone_implanted.vtk"/>
    </Block>
    <Block index="1" name="Prosthesis">
      <DataSet index="0" file="prosthesis.vtk"/>
    </Block>
    <Block index="2" name="Intact_Femur_Reference">
      <DataSet index="0" file="femur_intact_stress.vtk"/>
    </Block>
  </vtkMultiBlockDataSet>
</VTKFile>
'''

with open("femoral_implant_assembly.vtm", "w") as f:
    f.write(vtm_content)
print("  Created: femoral_implant_assembly.vtm")

# ============================================================================
# STATISTICS
# ============================================================================
print("\n" + "="*60)
print("GENERATION COMPLETE")
print("="*60)
print(f"Grid dimensions: {nx} x {ny} x {nz} = {nx*ny*nz:,} points per block")
print(f"Total data points: {nx*ny*nz*3:,} (3 blocks)")
print(f"\nFiles generated:")
print("  - femoral_implant_assembly.vtm (Multi-Block assembly)")
print("  - cortical_bone_implanted.vtk (Bone with stress shielding)")
print("  - prosthesis.vtk (Implant components)")
print("  - femur_intact_stress.vtk (Pre-operative baseline)")
print(f"\nLoading conditions: {APPLIED_LOAD:.0f} N ({HIP_LOAD_FACTOR}x body weight)")
