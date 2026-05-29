import trimesh
import numpy as np
import os

def generate_dense_benchmarks(num_categories=10, pairs_per_cat=20, min_v=500, max_v=250000, seed=42):
    """
    Generates meshes with granular vertex counts using UV spheres 
    instead of icosphere subdivisions to avoid the 'power of 4' jump.
    """
    # Fix the random seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    root_path = "/home/yury/Projects/PythonMesh/Data/BigDataV2"
    
    # Target vertex counts (geometric progression)
    v_targets = np.geomspace(min_v, max_v, num=num_categories).astype(int)

    for cat_idx, target in enumerate(v_targets):
        cat_folder = os.path.join(root_path, str(cat_idx))
        os.makedirs(cat_folder, exist_ok=True)
        
        # For a UV sphere, vertices approx = (count * count)
        # So we take the square root of the target to find the count
        sphere_res = int(np.sqrt(target)) 

        print(f"Category {cat_idx}: Target {target} | Creating spheres with res {sphere_res}...")

        for i in range(pairs_per_cat):
            # 1. Create UV sphere (much more granular than icosphere subdivision)
            mesh_a = trimesh.creation.uv_sphere(radius=1.0, count=[sphere_res, sphere_res])
            
            # 2. Add noise to break symmetry (helps CGAL benchmark 'real world' cases)
            noise_amp = 0.03
            mesh_a.vertices += np.random.normal(0, noise_amp, mesh_a.vertices.shape)
            
            # 3. Create mesh B with slightly different noise
            mesh_b = mesh_a.copy()
            mesh_b.vertices += np.random.normal(0, 0.01, mesh_b.vertices.shape)

            # 4. Center both
            mesh_a.apply_translation(-mesh_a.centroid)
            mesh_b.apply_translation(-mesh_b.centroid)
            
            # 5. Overlap offset (35% of total Z-height)
            # This ensures a consistent intersection area across different scales
            z_offset = (mesh_a.extents[2] + mesh_b.extents[2]) * 0.35
            mesh_b.apply_translation([0, 0, z_offset])

            # 6. Export
            mesh_a.export(os.path.join(cat_folder, f"A_{i}.off"))
            mesh_b.export(os.path.join(cat_folder, f"B_{i}.off"))

    print("\nGeneration complete. Check your folders for a smooth gradient of file sizes!")

if __name__ == "__main__":
    # max_v=50000 is a good "medium" range for CGAL exact kernels
    # to run in reasonable time while still being "heavy".
    # Added seed parameter to ensure consistent outputs.
    generate_dense_benchmarks(num_categories=10, pairs_per_cat=20, min_v=500, max_v=250000, seed=42)
