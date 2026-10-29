import numpy as np
import os
import meshio
import pyvista as pv

def verify_samples_file_exists(meshes_folder, samples_folder, mesh_file, d):
    mesh_name = os.path.splitext(os.path.basename(mesh_file))[0]
    samples_filepath = os.path.join(samples_folder, f"{mesh_name}_d{d}.off")
    print(f"Checking if samples file exists: {samples_filepath}")
    if os.path.exists(samples_filepath):
        return True

    # Generate the file if it does not exist.
    input_mesh_filepath = os.path.join(meshes_folder, f"{mesh_name}.off")
    print(f"Generating samples for mesh: {input_mesh_filepath}")
    mesh_data = meshio.read(input_mesh_filepath)
    mesh_data_pv = pv.from_meshio(mesh_data)
    points, normals = poisson_sample_surface(mesh_data_pv, d)
    write_samples_off(points, normals, samples_filepath)


def poisson_sample_surface(mesh, d: float, seed=42):
    """
    Approximately Poisson-disk sample a triangle mesh surface.

    Samples lie on triangle faces, with minimum spacing approximately d.
    Returns sampled points and interpolated vertex normals.
    """
    if d <= 0:
        raise ValueError("d must be positive")

    mesh_min_x = np.min(np.asarray(mesh.points)[:, 0])
    mesh_min_y = np.min(np.asarray(mesh.points)[:, 1])
    mesh_min_z = np.min(np.asarray(mesh.points)[:, 2])

    mesh_max_x = np.max(np.asarray(mesh.points)[:, 0])
    mesh_max_y = np.max(np.asarray(mesh.points)[:, 1])
    mesh_max_z = np.max(np.asarray(mesh.points)[:, 2])

    max_diff = max(
        mesh_max_x - mesh_min_x,
        mesh_max_y - mesh_min_y,
        mesh_max_z - mesh_min_z
    )

    d = d * max_diff

    mesh = mesh.extract_surface().triangulate()
    mesh = mesh.compute_normals(
        point_normals=True,
        cell_normals=False,
        consistent_normals=True,
        auto_orient_normals=True,
    )

    vertices = np.asarray(mesh.points)
    vertex_normals = np.asarray(mesh.point_data["Normals"])

    # Extract triangle vertex indices.
    faces = np.asarray(mesh.faces).reshape(-1, 4)
    if not np.all(faces[:, 0] == 3):
        raise ValueError("Expected a triangulated mesh")
    triangles = faces[:, 1:]

    tri = vertices[triangles]
    cross = np.cross(tri[:, 1] - tri[:, 0],
                     tri[:, 2] - tri[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)

    valid = areas > 1e-15
    triangles = triangles[valid]
    areas = areas[valid]

    if len(areas) == 0:
        raise ValueError("Mesh contains no non-degenerate triangles")

    # Triangle selection probabilities proportional to area.
    probabilities = areas / areas.sum()
    rng = np.random.default_rng(seed)

    # Spatial hash: cells have side length d / sqrt(3).
    # Neighboring cells within the minimum-distance radius are checked.
    cell_size = d / np.sqrt(3.0)
    neighbor_range = int(np.ceil(d / cell_size))
    grid = {}
    samples = []
    sample_normals = []

    def cell_key(p):
        return tuple(np.floor(p / cell_size).astype(int))

    def accept(p):
        key = cell_key(p)

        for dx in range(-neighbor_range, neighbor_range + 1):
            for dy in range(-neighbor_range, neighbor_range + 1):
                for dz in range(-neighbor_range, neighbor_range + 1):
                    bucket = grid.get(
                        (key[0] + dx, key[1] + dy, key[2] + dz), []
                    )
                    for idx in bucket:
                        if np.linalg.norm(samples[idx] - p) < d:
                            return False

        idx = len(samples)
        samples.append(p)
        grid.setdefault(key, []).append(idx)
        return True

    # Generate candidates uniformly over triangle areas.
    # Stop after a long run of rejected candidates.
    max_attempts = max(10000, len(vertices) * 100)
    failures = 0
    max_failures = max(5000, len(vertices) * 50)

    while failures < max_failures and len(samples) < max_attempts:
        ti = rng.choice(len(triangles), p=probabilities)
        ids = triangles[ti]
        a, b, c = vertices[ids]

        # Uniform random point in a triangle.
        u, v = rng.random(2)
        su = np.sqrt(u)
        bary = np.array([1 - su, su * (1 - v), su * v])
        p = bary @ np.array([a, b, c])

        if accept(p):
            # Interpolate the vertex normals using barycentric coordinates.
            n = bary @ vertex_normals[ids]
            norm = np.linalg.norm(n)

            if norm > 1e-15:
                sample_normals.append(n / norm)
            else:
                sample_normals.append(
                    cross[valid][ti] / np.linalg.norm(cross[valid][ti])
                )
            failures = 0
        else:
            failures += 1

    return np.asarray(samples), np.asarray(sample_normals)

def write_samples_off(points, normals, filename="samples.off"):
    points = np.asarray(points)
    normals = np.asarray(normals)

    assert points.shape == normals.shape
    assert points.ndim == 2 and points.shape[1] == 3

    # make sure that the folder exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, "w") as f:
        f.write("OFF\n")
        f.write(f"{len(points)} 0 0\n")

        for p, n in zip(points, normals):
            f.write(
                f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g} "
                f"{n[0]:.17g} {n[1]:.17g} {n[2]:.17g}\n"
            )