"""
Automated Stress Shielding Volume Calculator
ParaView Python Script for Orthopedic FEA Post-Processing

Author: Jeffrey Husch, Senior Biomedical Engineer
Purpose: Calculate total volume of bone subjected to stress shielding

Usage in ParaView:
  Tools > Python Shell > Run Script
  OR
  pvpython stress_shielding_volume_calculator.py
"""

# ============================================================================
# PARAVIEW SCRIPT - STRESS SHIELDING VOLUME CALCULATION
# ============================================================================

try:
    from paraview.simple import *
except ImportError:
    print("This script must be run within ParaView or pvpython")
    print("Standalone execution will demonstrate the algorithm logic")

import math

# Configuration Parameters
STRESS_SHIELDING_THRESHOLD = 0.5    # SED ratio below which bone is "shielded"
SEVERE_SHIELDING_THRESHOLD = 0.25  # Severe stress shielding threshold
SED_RESORPTION_THRESHOLD = 0.015   # SED (MPa) below which resorption occurs

def calculate_stress_shielding_volume_paraview():
    """
    ParaView implementation for stress shielding volume calculation.
    Requires active source with 'stress_shielding_index' scalar.
    """
    
    print("="*60)
    print("STRESS SHIELDING VOLUME ANALYSIS")
    print("="*60)
    
    # Get active source (should be cortical bone block)
    source = GetActiveSource()
    if source is None:
        print("ERROR: No active source. Load the cortical_bone_implanted.vtk first.")
        return
    
    # Update pipeline
    UpdatePipeline()
    
    # Get data information
    data_info = source.GetDataInformation()
    bounds = data_info.GetBounds()
    
    print(f"\nModel Bounds:")
    print(f"  X: {bounds[0]:.1f} to {bounds[1]:.1f} mm")
    print(f"  Y: {bounds[2]:.1f} to {bounds[3]:.1f} mm")
    print(f"  Z: {bounds[4]:.1f} to {bounds[5]:.1f} mm")
    
    # =========================================================================
    # METHOD 1: Threshold Filter Approach
    # =========================================================================
    print("\n--- Method 1: Threshold Filter Analysis ---")
    
    # Threshold for stress-shielded regions
    threshold_shielded = Threshold(Input=source)
    threshold_shielded.Scalars = ['POINTS', 'stress_shielding_index']
    threshold_shielded.LowerThreshold = STRESS_SHIELDING_THRESHOLD
    threshold_shielded.UpperThreshold = 1.0
    threshold_shielded.ThresholdMethod = 'Between'
    
    UpdatePipeline(proxy=threshold_shielded)
    
    # Integrate variables to get volume
    integrate_shielded = IntegrateVariables(Input=threshold_shielded)
    UpdatePipeline(proxy=integrate_shielded)
    
    # Get integrated volume
    integrated_data = servermanager.Fetch(integrate_shielded)
    shielded_volume = integrated_data.GetCellData().GetArray("Volume").GetValue(0) if integrated_data.GetNumberOfCells() > 0 else 0
    
    print(f"  Stress Shielding Index > {STRESS_SHIELDING_THRESHOLD}")
    print(f"  Volume affected: {shielded_volume:.2f} mm³")
    print(f"  Volume affected: {shielded_volume/1000:.2f} cm³")
    
    # Severe stress shielding
    threshold_severe = Threshold(Input=source)
    threshold_severe.Scalars = ['POINTS', 'stress_shielding_index']
    threshold_severe.LowerThreshold = 1.0 - SEVERE_SHIELDING_THRESHOLD
    threshold_severe.UpperThreshold = 1.0
    threshold_severe.ThresholdMethod = 'Between'
    
    UpdatePipeline(proxy=threshold_severe)
    integrate_severe = IntegrateVariables(Input=threshold_severe)
    UpdatePipeline(proxy=integrate_severe)
    
    severe_data = servermanager.Fetch(integrate_severe)
    severe_volume = severe_data.GetCellData().GetArray("Volume").GetValue(0) if severe_data.GetNumberOfCells() > 0 else 0
    
    print(f"\n  Severe Shielding (Index > {1-SEVERE_SHIELDING_THRESHOLD:.2f})")
    print(f"  Volume at risk: {severe_volume:.2f} mm³")
    
    # =========================================================================
    # METHOD 2: SED-Based Analysis
    # =========================================================================
    print("\n--- Method 2: Strain Energy Density Analysis ---")
    
    threshold_sed = Threshold(Input=source)
    threshold_sed.Scalars = ['POINTS', 'SED']
    threshold_sed.LowerThreshold = 0.0
    threshold_sed.UpperThreshold = SED_RESORPTION_THRESHOLD
    threshold_sed.ThresholdMethod = 'Between'
    
    UpdatePipeline(proxy=threshold_sed)
    integrate_sed = IntegrateVariables(Input=threshold_sed)
    UpdatePipeline(proxy=integrate_sed)
    
    sed_data = servermanager.Fetch(integrate_sed)
    resorption_volume = sed_data.GetCellData().GetArray("Volume").GetValue(0) if sed_data.GetNumberOfCells() > 0 else 0
    
    print(f"  SED < {SED_RESORPTION_THRESHOLD} MPa (Resorption threshold)")
    print(f"  Volume at resorption risk: {resorption_volume:.2f} mm³")
    
    # =========================================================================
    # REGIONAL ANALYSIS
    # =========================================================================
    print("\n--- Regional Analysis by Gruen Zone ---")
    
    # Define Gruen zones (simplified)
    gruen_zones = {
        "Zone 1 (Lateral Proximal)": (10, 25, 60, 85),
        "Zone 2 (Lateral Mid)": (10, 25, 30, 60),
        "Zone 3 (Lateral Distal)": (10, 25, 0, 30),
        "Zone 4 (Distal Tip)": (-10, 10, 0, 15),
        "Zone 5 (Medial Distal)": (-25, -10, 0, 30),
        "Zone 6 (Medial Mid)": (-25, -10, 30, 60),
        "Zone 7 (Medial Proximal)": (-25, -10, 60, 85),
    }
    
    for zone_name, (x_min, x_max, z_min, z_max) in gruen_zones.items():
        # Clip to zone
        clip_x_min = Clip(Input=source)
        clip_x_min.ClipType = 'Plane'
        clip_x_min.ClipType.Origin = [x_min, 0, 0]
        clip_x_min.ClipType.Normal = [1, 0, 0]
        
        clip_x_max = Clip(Input=clip_x_min)
        clip_x_max.ClipType = 'Plane'
        clip_x_max.ClipType.Origin = [x_max, 0, 0]
        clip_x_max.ClipType.Normal = [-1, 0, 0]
        
        clip_z_min = Clip(Input=clip_x_max)
        clip_z_min.ClipType = 'Plane'
        clip_z_min.ClipType.Origin = [0, 0, z_min]
        clip_z_min.ClipType.Normal = [0, 0, 1]
        
        clip_z_max = Clip(Input=clip_z_min)
        clip_z_max.ClipType = 'Plane'
        clip_z_max.ClipType.Origin = [0, 0, z_max]
        clip_z_max.ClipType.Normal = [0, 0, -1]
        
        # Threshold for shielding in this zone
        zone_threshold = Threshold(Input=clip_z_max)
        zone_threshold.Scalars = ['POINTS', 'stress_shielding_index']
        zone_threshold.LowerThreshold = STRESS_SHIELDING_THRESHOLD
        zone_threshold.UpperThreshold = 1.0
        
        UpdatePipeline(proxy=zone_threshold)
        
        zone_data = servermanager.Fetch(zone_threshold)
        zone_vol = 0
        if zone_data.GetNumberOfPoints() > 0:
            # Estimate volume from point count (approximate)
            zone_vol = zone_data.GetNumberOfPoints() * 1.0  # 1mm³ per point (1mm spacing)
        
        print(f"  {zone_name}: {zone_vol:.0f} mm³ affected")
        
        # Clean up
        Delete(clip_x_min)
        Delete(clip_x_max)
        Delete(clip_z_min)
        Delete(clip_z_max)
        Delete(zone_threshold)
    
    # =========================================================================
    # SUMMARY REPORT
    # =========================================================================
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total stress-shielded volume: {shielded_volume:.1f} mm³ ({shielded_volume/1000:.2f} cm³)")
    print(f"Severe shielding volume: {severe_volume:.1f} mm³ ({severe_volume/1000:.2f} cm³)")
    print(f"Bone resorption risk volume: {resorption_volume:.1f} mm³")
    print("\nClinical Interpretation:")
    if shielded_volume > 5000:
        print("  ⚠ HIGH RISK: Significant stress shielding detected")
        print("  Consider stem design modification or surface treatment")
    elif shielded_volume > 2000:
        print("  ⚡ MODERATE RISK: Notable stress shielding present")
        print("  Monitor with follow-up DXA scans")
    else:
        print("  ✓ LOW RISK: Minimal stress shielding")
        print("  Expected bone adaptation should be acceptable")
    
    # Clean up
    Delete(threshold_shielded)
    Delete(threshold_severe)
    Delete(threshold_sed)
    Delete(integrate_shielded)
    Delete(integrate_severe)
    Delete(integrate_sed)
    
    return shielded_volume, severe_volume, resorption_volume


def standalone_calculation():
    """
    Standalone calculation without ParaView for demonstration.
    Reads VTK file directly and calculates volumes.
    """
    print("="*60)
    print("STANDALONE STRESS SHIELDING CALCULATION")
    print("(For use outside ParaView environment)")
    print("="*60)
    
    import os
    
    vtk_file = "cortical_bone_implanted.vtk"
    if not os.path.exists(vtk_file):
        vtk_file = "../femoral-implant-assembly/cortical_bone_implanted.vtk"
    
    if not os.path.exists(vtk_file):
        print(f"VTK file not found: {vtk_file}")
        return
    
    # Parse VTK file
    shielding_values = []
    sed_values = []
    tissue_values = []
    
    reading_shielding = False
    reading_sed = False
    reading_tissue = False
    
    with open(vtk_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if "SCALARS stress_shielding_index" in line:
                reading_shielding = True
                reading_sed = False
                reading_tissue = False
                continue
            elif "SCALARS SED" in line:
                reading_sed = True
                reading_shielding = False
                reading_tissue = False
                continue
            elif "SCALARS tissue_type" in line:
                reading_tissue = True
                reading_shielding = False
                reading_sed = False
                continue
            elif line.startswith("SCALARS") or line.startswith("VECTORS"):
                reading_shielding = False
                reading_sed = False
                reading_tissue = False
                continue
            elif line.startswith("LOOKUP_TABLE"):
                continue
            
            try:
                val = float(line)
                if reading_shielding:
                    shielding_values.append(val)
                elif reading_sed:
                    sed_values.append(val)
                elif reading_tissue:
                    tissue_values.append(val)
            except ValueError:
                pass
    
    # Calculate volumes (1mm³ per voxel)
    voxel_volume = 1.0  # mm³
    
    # Count stress-shielded voxels
    moderate_shielding = sum(1 for s, t in zip(shielding_values, tissue_values) 
                             if s >= STRESS_SHIELDING_THRESHOLD and t > 0)
    severe_shielding = sum(1 for s, t in zip(shielding_values, tissue_values) 
                           if s >= (1 - SEVERE_SHIELDING_THRESHOLD) and t > 0)
    low_sed = sum(1 for s, t in zip(sed_values, tissue_values) 
                  if s < SED_RESORPTION_THRESHOLD and s > 0 and t > 0)
    total_bone = sum(1 for t in tissue_values if t > 0)
    
    shielded_vol = moderate_shielding * voxel_volume
    severe_vol = severe_shielding * voxel_volume
    resorption_vol = low_sed * voxel_volume
    total_vol = total_bone * voxel_volume
    
    print(f"\nTotal bone volume: {total_vol:.0f} mm³ ({total_vol/1000:.1f} cm³)")
    print(f"\nStress Shielding Analysis:")
    print(f"  Moderate shielding (>{STRESS_SHIELDING_THRESHOLD}): {shielded_vol:.0f} mm³ ({shielded_vol/total_vol*100:.1f}%)")
    print(f"  Severe shielding (>{1-SEVERE_SHIELDING_THRESHOLD:.2f}): {severe_vol:.0f} mm³ ({severe_vol/total_vol*100:.1f}%)")
    print(f"  Low SED (resorption risk): {resorption_vol:.0f} mm³ ({resorption_vol/total_vol*100:.1f}%)")
    
    return shielded_vol, severe_vol, resorption_vol


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    try:
        # Try ParaView execution
        from paraview.simple import GetActiveSource
        if GetActiveSource() is not None:
            calculate_stress_shielding_volume_paraview()
        else:
            print("No active source in ParaView. Running standalone calculation.")
            standalone_calculation()
    except ImportError:
        # Run standalone version
        standalone_calculation()
