import numpy as np


def create_discretized_sphere(radius=1.0, k=1):
    """
    Creates a discretized sphere based on the provided recursive method.
    This implementation ensures that the initial icosahedron is aligned with the
    Z-axis, including (0,0,±R) as vertices.

    Args:
        radius (float): Geometric radius of the sphere (R_g).
        k (int): Discretization factor, must be an integer >= 1.
                 k=1 corresponds to an icosahedron.

    Returns:
        numpy.ndarray: A numpy.ndarray with shape (N, 3), containing the 3D
                       coordinates of all discrete points.
                       Where N = 10 * 4**(k-1) + 2.
    """
    if not isinstance(k, int) or k < 1:
        raise ValueError("Discretization factor k must be an integer greater than or equal to 1.")
    if not isinstance(radius, (int, float)) or radius <= 0:
        raise ValueError("Radius must be a positive number.")

    # --- Step 1: Create a Z-axis-aligned icosahedron (k=1) ---
    # These vertices are generated programmatically to ensure precision and
    # correct Z-axis alignment. It includes the north pole (0,0,1) and
    # the south pole (0,0,-1).

    # 1. Define the initial positions of the 12 vertices (on a unit sphere)
    vertices_list = []
    # North pole
    vertices_list.append([0.0, 0.0, 1.0])

    # Calculate the height (z) and radius (r_xy) of the upper 5 vertices.
    # This angle can be derived from the geometry of the icosahedron.
    angle_from_z = np.arctan(2.0)
    z_coord = np.cos(angle_from_z)
    r_xy = np.sin(angle_from_z)

    # Add the 5 upper vertices
    for i in range(5):
        phi = (2.0 * np.pi * i) / 5.0
        vertices_list.append([r_xy * np.cos(phi), r_xy * np.sin(phi), z_coord])

    # Add the 5 lower vertices (rotated by pi/5 relative to the upper ones)
    for i in range(5):
        phi = (2.0 * np.pi * i + np.pi) / 5.0
        vertices_list.append([r_xy * np.cos(phi), r_xy * np.sin(phi), -z_coord])

    # South pole
    vertices_list.append([0.0, 0.0, -1.0])

    vertices = np.array(vertices_list, dtype=float)

    # 2. Define the 20 triangular faces (based on the new vertex order)
    # Vertex indices: 0=North Pole, 1-5=Upper Layer, 6-10=Lower Layer, 11=South Pole
    faces = [
        # Top 5 faces (connecting to the North Pole)
        [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1],
        # Middle 10 faces
        [1, 6, 2], [2, 7, 3], [3, 8, 4], [4, 9, 5], [5, 10, 1],
        [2, 6, 7], [3, 7, 8], [4, 8, 9], [5, 9, 10], [1, 10, 6],
        # Bottom 5 faces (connecting to the South Pole)
        [11, 7, 6], [11, 8, 7], [11, 9, 8], [11, 10, 9], [11, 6, 10]
    ]

    # --- Step 2: Recursively subdivide k-1 times ---
    for _ in range(k - 1):
        midpoint_cache = {}
        new_faces = []
        for face_indices in faces:
            v1_idx, v2_idx, v3_idx = face_indices

            # Define the three edges of the triangle
            edges = [(v1_idx, v2_idx), (v2_idx, v3_idx), (v3_idx, v1_idx)]
            midpoint_indices = []

            for p1_idx, p2_idx in edges:
                edge_key = tuple(sorted((p1_idx, p2_idx)))

                if edge_key not in midpoint_cache:
                    p1 = vertices[p1_idx]
                    p2 = vertices[p2_idx]

                    # Calculate the midpoint and project it back onto the unit sphere
                    mid_point = (p1 + p2) / 2.0
                    mid_point /= np.linalg.norm(mid_point)

                    # Add the new vertex to the vertices array
                    vertices = np.vstack([vertices, mid_point])
                    midpoint_cache[edge_key] = len(vertices) - 1

                midpoint_indices.append(midpoint_cache[edge_key])

            m1_idx, m2_idx, m3_idx = midpoint_indices

            # Replace the original large triangle with 4 new smaller ones
            new_faces.extend([
                [v1_idx, m1_idx, m3_idx],
                [v2_idx, m2_idx, m1_idx],
                [v3_idx, m3_idx, m2_idx],
                [m1_idx, m2_idx, m3_idx]
            ])

        faces = new_faces

    # --- Step 3: Scale by the specified radius and return the final vertices ---
    final_vertices = vertices * radius

    # Verify that the final number of vertices matches the formula
    expected_num_vertices = 10 * (4 ** (k - 1)) + 2
    if final_vertices.shape[0] != expected_num_vertices:
        print(f"Warning: Vertex count mismatch. Expected: {expected_num_vertices}, "
              f"Generated: {final_vertices.shape[0]}")

    return final_vertices


def main():
    """
    Main function to get user input and run the sphere discretization.
    """
    try:
        radius_str = input("Please enter the sphere radius (R_g), e.g., 1.0: ")
        radius = float(radius_str)

        k_str = input("Please enter the discretization factor (k, integer >= 1): ")
        k = int(k_str)

    except ValueError:
        print("\nError: Invalid input. Please enter valid numbers.")
        return
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        return

    print(f"\nGenerating sphere mesh with radius {radius} and factor k={k}...")
    try:
        discretized_points = create_discretized_sphere(radius, k)
    except ValueError as e:
        print(f"\nError: {e}")
        return

    # --- Use scientific notation and remove header ---
    output_filename = f"sphere_N_{discretized_points.shape[0]}_R_{radius_str}.txt"
    # Use '%.15e' format for scientific notation, with space as delimiter.
    # No header is written, to match the example format.
    np.savetxt(output_filename, discretized_points, fmt='%.15e', delimiter=' ')

    print("\n✅ Success!")
    print(f"Total vertices generated: {discretized_points.shape[0]}")
    print(f"Mesh data saved to file: {output_filename}")
    print(f"File format: One 3D vertex (x y z) per line, space-separated in scientific notation.")


if __name__ == "__main__":
    main()
