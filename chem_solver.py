import numpy as np
from numba import njit, prange
from numpy.polynomial.laguerre import laggauss

# =============================================================================
# Global Configurations for Gauss-Laguerre Quadrature (Local Part)
# Recommend QUAD_DEGREE = 3~5 for high precision numerical integration
# =============================================================================
QUAD_DEGREE = 3
GL_X_ARRAY, GL_W_ARRAY = laggauss(QUAD_DEGREE)

# Multiply the weight by e^(x_k) to perfectly cancel out the intrinsic e^{-x}
# assumption of the Gauss-Laguerre quadrature algorithm.
GL_W_EXP_X_ARRAY = GL_W_ARRAY * np.exp(GL_X_ARRAY)


# =============================================================================
# Helper Functions for Coordinate Transformations
# =============================================================================
@njit(fastmath=True)
def rotation_matrix_2d(theta):
    """
    Return the rotation matrix representing rotation
    by given an angle of rotation theta,
    which represents a rotation clockwise about the vector phi of magnitude phi.
    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])


@njit(fastmath=True)
def quaternion_to_rotation_matrix(q):
    """
    Converts a quaternion [x, y, z, w] to a 3x3 rotation matrix.
    """
    norm_q = np.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
    if norm_q > 1e-12:
        x = q[0] / norm_q
        y = q[1] / norm_q
        z = q[2] / norm_q
        w = q[3] / norm_q
    else:
        x = 0.0
        y = 0.0
        z = 0.0
        w = 1.0

    xx = x * x
    yy = y * y
    zz = z * z
    xy = x * y
    xz = x * z
    yz = y * z
    wx = w * x
    wy = w * y
    wz = w * z

    matrix = np.empty((3, 3))
    matrix[0, 0] = 1.0 - 2.0 * (yy + zz)
    matrix[0, 1] = 2.0 * (xy - wz)
    matrix[0, 2] = 2.0 * (xz + wy)

    matrix[1, 0] = 2.0 * (xy + wz)
    matrix[1, 1] = 1.0 - 2.0 * (xx + zz)
    matrix[1, 2] = 2.0 * (yz - wx)

    matrix[2, 0] = 2.0 * (xz - wy)
    matrix[2, 1] = 2.0 * (yz + wx)
    matrix[2, 2] = 1.0 - 2.0 * (xx + yy)
    return matrix


@njit(parallel=True, fastmath=True)
def _precompute_rotated_nodes_2d(n_history_steps, n_nodes, location_history, orientation_history, nodes_body_frame):
    """
    Precomputes the lab-frame positions of a 2D Janus particle's surface nodes across all history steps.
    """
    rotated_nodes_hist = np.zeros((n_history_steps, n_nodes, 2))
    for i in prange(n_history_steps):
        cos_t = orientation_history[i, 0]
        sin_t = orientation_history[i, 1]
        rot_matrix = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
        rotated_nodes_hist[i] = location_history[i] + np.dot(nodes_body_frame, rot_matrix.T)
    return rotated_nodes_hist


@njit(parallel=True, fastmath=True)
def _precompute_rotated_nodes_3d(n_history_steps, n_nodes, location_history, orientation_history, nodes_body_frame):
    """
    Precomputes the lab-frame positions of a 3D Janus particle's surface nodes across all history steps.
    """
    rotated_nodes_hist = np.zeros((n_history_steps, n_nodes, 3))
    for i in prange(n_history_steps):
        rot_matrix = quaternion_to_rotation_matrix(orientation_history[i])
        rotated_nodes_hist[i] = location_history[i] + np.dot(nodes_body_frame, rot_matrix.T)
    return rotated_nodes_hist


# =============================================================================
# 2D Kernel Evaluations (Numba Optimized)
# =============================================================================
@njit(parallel=True, fastmath=True)
def _history_part_point_2d_numba(target_nodes, n_history_steps, source_location_history, t_now, dt, D):
    n_targets = len(target_nodes)
    total_grad = np.zeros((n_targets, 2))
    for j in prange(n_targets):
        target_pos = target_nodes[j]
        grad_sum = np.zeros(2)
        for step in range(n_history_steps):
            tau = t_now - step * dt
            source_pos = source_location_history[step]
            r_vec = target_pos - source_pos
            r_sq = r_vec[0] ** 2 + r_vec[1] ** 2

            if r_sq > 1e-12:
                const = 4.0 * D * tau
                G = (1.0 / (np.pi * const)) * np.exp(-r_sq / const)
                grad_G = (-r_vec / (2.0 * D * tau)) * G

                weight = dt
                if step == 0 or step == n_history_steps - 1:
                    weight = 0.5 * dt
                grad_sum += grad_G * weight
        total_grad[j] = grad_sum
    return total_grad


@njit(parallel=True, fastmath=True)
def _history_part_dist_2d_numba_optimized(target_nodes, rotated_source_nodes_hist, sigma_dist, t_now, dt, D):
    n_targets = len(target_nodes)
    n_history_steps = rotated_source_nodes_hist.shape[0]
    n_sources = rotated_source_nodes_hist.shape[1]
    total_grad = np.zeros((n_targets, 2))
    sigma_dist_scaled = sigma_dist.reshape(-1, 1)

    # Calculate regularization parameter (epsilon squared) based on mesh spacing.
    # This desingularizes the point sources into continuous blobs.
    epsilon_sq = (2.0 * np.pi / n_sources) ** 2

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        grad_sum = np.zeros(2)
        for step in range(n_history_steps):
            tau = t_now - step * dt
            source_nodes_lab = rotated_source_nodes_hist[step]

            r_vecs = target_pos - source_nodes_lab

            # Add epsilon_sq to the squared distance to prevent logarithmic singularity
            # from self-interaction as dt approaches 0.
            r_sqs = r_vecs[:, 0] ** 2 + r_vecs[:, 1] ** 2 + epsilon_sq

            const = 4.0 * D * tau
            valid_indices = r_sqs > 1e-12

            G = np.zeros(n_sources)
            G[valid_indices] = (1.0 / (np.pi * const)) * np.exp(-r_sqs[valid_indices] / const)

            grad_G = (-r_vecs / (2.0 * D * tau)) * G.reshape(-1, 1)

            weight = dt
            if step == 0 or step == n_history_steps - 1:
                weight = 0.5 * dt
            grad_sum += np.sum(grad_G * sigma_dist_scaled, axis=0) * weight
        total_grad[j] = grad_sum
    return total_grad


@njit(parallel=True, fastmath=True)
def _local_part_dist_2d_numba_test(target_nodes, source_pos_now, source_orient_now,
                                   source_pos_before, source_nodes_body_frame, sigma_dist, dt, D):
    n_targets = len(target_nodes)
    n_sources = len(source_nodes_body_frame)
    total_grad_C_local = np.zeros((n_targets, 2))

    v_k = (source_pos_now - source_pos_before) / dt
    cos_t = source_orient_now[0]
    sin_t = source_orient_now[1]
    rot_matrix_now = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
    source_nodes_lab_now = source_pos_now + np.dot(source_nodes_body_frame, rot_matrix_now.T)
    lower_limit = 1.0 / dt

    # Calculate regularization parameter
    epsilon_sq = (2.0 * np.pi / n_sources) ** 2

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        grad_sum = np.zeros(2)

        for i in range(n_sources):
            if i == j:
                continue

            sigma_val = sigma_dist[i]
            if sigma_val == 0.0:
                continue

            source_pos_i = source_nodes_lab_now[i]

            R0_x = target_pos[0] - source_pos_i[0]
            R0_y = target_pos[1] - source_pos_i[1]

            # Add epsilon_sq to regularize the distance scaling for adjacent nodes
            R0_sq = R0_x ** 2 + R0_y ** 2 + epsilon_sq

            A_ij = R0_sq / (4.0 * D)

            pair_grad_sum_x = 0.0
            pair_grad_sum_y = 0.0

            for q in range(len(GL_X_ARRAY)):
                x_k = GL_X_ARRAY[q]
                w_exp_x = GL_W_EXP_X_ARRAY[q]

                y_k = x_k / A_ij
                u_val = lower_limit + y_k
                tau = 1.0 / u_val

                R_tau_x = R0_x + v_k[0] * tau
                R_tau_y = R0_y + v_k[1] * tau

                # Add epsilon_sq to regularize the dynamic distance
                R_tau_sq = R_tau_x ** 2 + R_tau_y ** 2 + epsilon_sq

                const = 4.0 * D * tau
                exponent = -R_tau_sq / const

                if exponent < -100.0:
                    continue

                G = (1.0 / (np.pi * const)) * np.exp(exponent)
                jacobian_scale = (tau ** 2) / A_ij

                grad_G_x = (-R_tau_x / (2.0 * D * tau)) * G * jacobian_scale
                grad_G_y = (-R_tau_y / (2.0 * D * tau)) * G * jacobian_scale

                pair_grad_sum_x += grad_G_x * w_exp_x
                pair_grad_sum_y += grad_G_y * w_exp_x

            grad_sum[0] += pair_grad_sum_x * sigma_val
            grad_sum[1] += pair_grad_sum_y * sigma_val

        total_grad_C_local[j] = grad_sum

    return total_grad_C_local


@njit(parallel=True, fastmath=True)
def _local_part_point_2d_numba(target_nodes, source_pos_now, source_pos_before, dt, D):
    n_targets = len(target_nodes)
    total_grad_C_local = np.zeros((n_targets, 2))
    v_k = (source_pos_now - source_pos_before) / dt
    lower_limit = 1.0 / dt

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        weighted_grad_sum = np.zeros(2)

        for q in range(len(GL_X_ARRAY)):
            x_k = GL_X_ARRAY[q]
            w_exp_x = GL_W_EXP_X_ARRAY[q]

            u_val = x_k + lower_limit
            tau = 1.0 / u_val

            if tau <= 1e-12:
                continue

            source_pos_approx = source_pos_now - v_k * tau
            r_vec = target_pos - source_pos_approx
            r_sq = r_vec[0] ** 2 + r_vec[1] ** 2

            if r_sq > 1e-12:
                const = 4.0 * D * tau
                G = (1.0 / (np.pi * const)) * np.exp(-r_sq / const)
                # Multiply by Jacobian (tau^2) due to u = 1/tau substitution
                grad_G = (-r_vec / (2.0 * D * tau)) * G * (tau ** 2)
                weighted_grad_sum += grad_G * w_exp_x

        total_grad_C_local[j] = weighted_grad_sum
    return total_grad_C_local


@njit(parallel=True, fastmath=True)
def _local_part_dist_2d_numba_optimized(target_nodes, source_pos_now, source_orient_now,
                                        source_pos_before, source_nodes_body_frame, sigma_dist, dt, D):
    n_targets = len(target_nodes)
    n_sources = len(source_nodes_body_frame)
    total_grad_C_local = np.zeros((n_targets, 2))

    v_k = (source_pos_now - source_pos_before) / dt
    cos_t = source_orient_now[0]
    sin_t = source_orient_now[1]
    rot_matrix_now = np.array([[cos_t, -sin_t], [sin_t, cos_t]])

    sigma_dist_scaled = sigma_dist.reshape(-1, 1)
    lower_limit = 1.0 / dt

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        weighted_grad_sum = np.zeros(2)

        for q in range(len(GL_X_ARRAY)):
            x_k = GL_X_ARRAY[q]
            w_exp_x = GL_W_EXP_X_ARRAY[q]

            u_val = x_k + lower_limit
            tau = 1.0 / u_val

            if tau <= 1e-12:
                continue

            source_pos_approx = source_pos_now - v_k * tau
            source_nodes_lab_approx = source_pos_approx + np.dot(source_nodes_body_frame, rot_matrix_now.T)

            r_vecs = target_pos - source_nodes_lab_approx
            r_sqs = r_vecs[:, 0] ** 2 + r_vecs[:, 1] ** 2

            const = 4.0 * D * tau
            valid_indices = r_sqs > 1e-12

            G = np.zeros(n_sources)
            G[valid_indices] = (1.0 / (np.pi * const)) * np.exp(-r_sqs[valid_indices] / const)

            grad_G = (-r_vecs / (2.0 * D * tau)) * G.reshape(-1, 1) * (tau ** 2)
            weighted_grad_sum += np.sum(grad_G * sigma_dist_scaled, axis=0) * w_exp_x

        total_grad_C_local[j] = weighted_grad_sum
    return total_grad_C_local


# =============================================================================
# 3D Kernel Evaluations (Numba Optimized)
# =============================================================================
@njit(parallel=True, fastmath=True)
def _history_part_point_3d_numba(target_nodes, n_history_steps, source_location_history, t_now, dt, D):
    n_targets = len(target_nodes)
    total_grad = np.zeros((n_targets, 3))
    for j in prange(n_targets):
        target_pos = target_nodes[j]
        grad_sum = np.zeros(3)
        for step in range(n_history_steps):
            tau = t_now - step * dt
            source_pos = source_location_history[step]
            r_vec = target_pos - source_pos
            r_sq = r_vec[0] ** 2 + r_vec[1] ** 2 + r_vec[2] ** 2

            if r_sq > 1e-12:
                const = 4.0 * np.pi * D * tau
                G = np.power(const, -1.5) * np.exp(-r_sq / const)
                grad_G = (-r_vec / (2.0 * D * tau)) * G

                weight = dt
                if step == 0 or step == n_history_steps - 1:
                    weight = 0.5 * dt
                grad_sum += grad_G * weight
        total_grad[j] = grad_sum
    return total_grad


@njit(parallel=True, fastmath=True)
def _history_part_dist_3d_numba_optimized(target_nodes, rotated_source_nodes_hist, sigma_dist, t_now, dt, D):
    n_targets = len(target_nodes)
    n_history_steps = rotated_source_nodes_hist.shape[0]
    n_sources = rotated_source_nodes_hist.shape[1]
    total_grad_C = np.zeros((n_targets, 3))
    sigma_dist_scaled = sigma_dist.reshape(-1, 1)

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        grad_sum = np.zeros(3)
        for step in range(n_history_steps):
            tau = t_now - step * dt
            source_nodes_lab = rotated_source_nodes_hist[step]

            r_vecs = target_pos - source_nodes_lab
            r_sqs = r_vecs[:, 0] ** 2 + r_vecs[:, 1] ** 2 + r_vecs[:, 2] ** 2

            const = 4.0 * np.pi * D * tau
            valid_indices = r_sqs > 1e-12

            G = np.zeros(n_sources)
            G[valid_indices] = np.power(const, -1.5) * np.exp(-r_sqs[valid_indices] / const)

            grad_G = (-r_vecs / (2.0 * D * tau)) * G.reshape(-1, 1)

            weight = dt
            if step == 0 or step == n_history_steps - 1:
                weight = 0.5 * dt
            grad_sum += np.sum(grad_G * sigma_dist_scaled, axis=0) * weight

        total_grad_C[j] = grad_sum
    return total_grad_C


@njit(parallel=True, fastmath=True)
def _local_part_point_3d_numba(target_nodes, source_pos_now, source_pos_before, dt, D):
    n_targets = len(target_nodes)
    total_grad_C_local = np.zeros((n_targets, 3))
    v_k = (source_pos_now - source_pos_before) / dt
    lower_limit = 1.0 / dt

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        weighted_grad_sum = np.zeros(3)

        for q in range(len(GL_X_ARRAY)):
            x_k = GL_X_ARRAY[q]
            w_exp_x = GL_W_EXP_X_ARRAY[q]

            u_val = x_k + lower_limit
            tau = 1.0 / u_val

            if tau <= 1e-12:
                continue

            source_pos_approx = source_pos_now - v_k * tau
            r_vec = target_pos - source_pos_approx
            r_sq = np.sum(r_vec ** 2)

            if r_sq > 1e-12:
                const = 4.0 * np.pi * D * tau
                G = np.power(const, -1.5) * np.exp(-r_sq / const)
                grad_G = (-r_vec / (2.0 * D * tau)) * G * (tau ** 2)
                weighted_grad_sum += grad_G * w_exp_x

        total_grad_C_local[j] = weighted_grad_sum
    return total_grad_C_local


@njit(parallel=True, fastmath=True)
def _local_part_dist_3d_numba_optimized(target_nodes, source_pos_now, source_orient_now,
                                        source_pos_before, source_nodes_body_frame, sigma_dist, dt, D):
    n_targets = len(target_nodes)
    n_sources = len(source_nodes_body_frame)
    total_grad_C_local = np.zeros((n_targets, 3))

    v_k = (source_pos_now - source_pos_before) / dt

    # Freeze rotation at the current time t_now to satisfy the O(dt^2)
    # local truncation error theory derived in the paper.
    rot_matrix_now = quaternion_to_rotation_matrix(source_orient_now)

    # Precompute the lab-frame positions of the source nodes at t_now
    source_nodes_lab_now = source_pos_now + np.dot(source_nodes_body_frame, rot_matrix_now.T)
    lower_limit = 1.0 / dt

    for j in prange(n_targets):
        target_pos = target_nodes[j]
        grad_sum = np.zeros(3)

        for i in range(n_sources):
            # Skip self-interaction (the distance R_0 = 0 would cause a division by zero
            # when A_ij is used as a denominator in the scaling step)
            if i == j:
                continue

            sigma_val = sigma_dist[i]
            if sigma_val == 0.0:
                continue

            source_pos_i = source_nodes_lab_now[i]

            # Calculate the initial distance vector R_0 at t_now (with frozen rotation)
            R0_x = target_pos[0] - source_pos_i[0]
            R0_y = target_pos[1] - source_pos_i[1]
            R0_z = target_pos[2] - source_pos_i[2]
            R0_sq = R0_x ** 2 + R0_y ** 2 + R0_z ** 2

            # The CRUCIAL scaling factor A_ij
            A_ij = R0_sq / (4.0 * D)

            pair_grad_sum_x = 0.0
            pair_grad_sum_y = 0.0
            pair_grad_sum_z = 0.0

            for q in range(len(GL_X_ARRAY)):
                x_k = GL_X_ARRAY[q]
                w_exp_x = GL_W_EXP_X_ARRAY[q]

                # Affine transformation mapping with appropriate scaling applied
                y_k = x_k / A_ij
                u_val = lower_limit + y_k
                tau = 1.0 / u_val

                R_tau_x = R0_x + v_k[0] * tau
                R_tau_y = R0_y + v_k[1] * tau
                R_tau_z = R0_z + v_k[2] * tau
                R_tau_sq = R_tau_x ** 2 + R_tau_y ** 2 + R_tau_z ** 2

                const = 4.0 * np.pi * D * tau
                exponent = -R_tau_sq / (4.0 * D * tau)

                # Prevent floating-point underflow for extremely small exponential values
                if exponent < -100.0:
                    continue

                # The coefficient for the 3D Gaussian kernel is (4 * pi * D * tau)^(-1.5)
                G = np.power(const, -1.5) * np.exp(exponent)

                # Multiply by the Jacobian (tau**2) from the u = 1/tau substitution
                # AND the affine transformation derivative (1.0 / A_ij)
                jacobian_scale = (tau ** 2) / A_ij

                grad_G_x = (-R_tau_x / (2.0 * D * tau)) * G * jacobian_scale
                grad_G_y = (-R_tau_y / (2.0 * D * tau)) * G * jacobian_scale
                grad_G_z = (-R_tau_z / (2.0 * D * tau)) * G * jacobian_scale

                pair_grad_sum_x += grad_G_x * w_exp_x
                pair_grad_sum_y += grad_G_y * w_exp_x
                pair_grad_sum_z += grad_G_z * w_exp_x

            grad_sum[0] += pair_grad_sum_x * sigma_val
            grad_sum[1] += pair_grad_sum_y * sigma_val
            grad_sum[2] += pair_grad_sum_z * sigma_val

        total_grad_C_local[j] = grad_sum

    return total_grad_C_local


# =============================================================================
# Main Dispatchers for 2D and 3D Systems
# =============================================================================
def history_local_compose_2d_multi_body(target_body, all_bodies, dt, *args, **kwargs):
    """
    Calculates the 2D chemical force and torque on a target body.
    Fully optimized to precompute rotating sources and process target arrays natively.
    """
    step = kwargs.get('step')
    if step == 0:
        return np.zeros(2), 0.0

    t_now = step * dt
    target_surface_nodes = target_body.get_surface_nodes()
    total_grad_on_target_nodes = np.zeros_like(target_surface_nodes)

    for source_body in all_bodies:
        D = 1.0 / source_body.peclet_number
        grad_C_history = np.zeros_like(target_surface_nodes)
        grad_C_local = np.zeros_like(target_surface_nodes)
        is_source_janus = source_body.is_janus

        # --- History Part ---
        if step > 1:
            n_history_steps = step
            if is_source_janus:
                rotated_source_nodes_hist = _precompute_rotated_nodes_2d(
                    n_history_steps, source_body.n_nodes,
                    np.array(source_body.location_history[:n_history_steps]),
                    np.array(source_body.orientation_history[:n_history_steps]),
                    source_body.nodes_body_frame
                )
                grad_C_history = _history_part_dist_2d_numba_optimized(
                    target_surface_nodes, rotated_source_nodes_hist,
                    source_body.sigma_distribution, t_now, dt, D
                )
            else:
                grad_C_history = _history_part_point_2d_numba(
                    target_surface_nodes, n_history_steps,
                    np.array(source_body.location_history[:n_history_steps]),
                    t_now, dt, D
                )

        # --- Local Part ---
        if step > 0:
            if is_source_janus:
                grad_C_local = _local_part_dist_2d_numba_test(
                    target_surface_nodes, source_body.location, source_body.orientation,
                    source_body.location_history[step - 1], source_body.nodes_body_frame,
                    source_body.sigma_distribution, dt, D
                )
            else:
                grad_C_local = _local_part_point_2d_numba(
                    target_surface_nodes, source_body.location,
                    source_body.location_history[step - 1], dt, D
                )

        total_grad_on_target_nodes += grad_C_history + grad_C_local

    # --- Net Force & Torque Calculation (Lab Frame) ---
    r_vectors_from_center = target_surface_nodes - target_body.location
    norm_r = np.linalg.norm(r_vectors_from_center, axis=1, keepdims=True)
    norm_r[norm_r < 1e-12] = 1.0
    normals = r_vectors_from_center / norm_r

    grad_dot_norm = np.sum(total_grad_on_target_nodes * normals, axis=1, keepdims=True)
    grad_normal_component = grad_dot_norm * normals
    tangential_grad = total_grad_on_target_nodes - grad_normal_component

    # Extract local mobility distribution of the target particle surface (reshape for broadcasting)
    local_alpha = target_body.mobility_distribution.reshape(-1, 1)

    # Calculate actual slip velocity: v_s = alpha_local * grad_s(C)
    weighted_tangential_grad = tangential_grad * local_alpha

    # Surface integral approximations
    final_chem_force = np.mean(weighted_tangential_grad, axis=0)

    # 2D Torque = r_x * F_y - r_y * F_x
    # torque_z = r_vectors_from_center[:, 0] * tangential_grad[:, 1] - r_vectors_from_center[:, 1] * tangential_grad[:, 0]

    # 2D Torque = r_x * F_y - r_y * F_x (Using weighted slip velocity for torque)
    torque_z = r_vectors_from_center[:, 0] * weighted_tangential_grad[:, 1] - r_vectors_from_center[:, 1] * weighted_tangential_grad[:, 0]
    final_chem_torque = np.mean(torque_z)

    return final_chem_force, final_chem_torque


def history_local_compose_3d_multi_body(target_body, all_bodies, dt, *args, **kwargs):
    """
    Calculates the 3D chemical force and torque on a target body.
    Correctly computes the cross-product for torque using Lab-Frame radial vectors.
    """
    step = kwargs.get('step')
    if step == 0:
        return np.zeros(3), np.zeros(3)

    t_now = step * dt
    target_surface_nodes = target_body.get_surface_nodes()
    total_grad_on_target_nodes = np.zeros_like(target_surface_nodes)

    for source_body in all_bodies:
        D = 1.0 / source_body.peclet_number
        grad_C_history = np.zeros_like(target_surface_nodes)
        grad_C_local = np.zeros_like(target_surface_nodes)
        is_source_janus = source_body.is_janus

        # --- History Part ---
        if step > 1:
            n_history_steps = step
            if is_source_janus:
                rotated_source_nodes_hist = _precompute_rotated_nodes_3d(
                    n_history_steps, source_body.n_nodes,
                    np.array(source_body.location_history[:n_history_steps]),
                    np.array(source_body.orientation_history[:n_history_steps]),
                    source_body.nodes_body_frame
                )
                grad_C_history = _history_part_dist_3d_numba_optimized(
                    target_surface_nodes, rotated_source_nodes_hist,
                    source_body.sigma_distribution, t_now, dt, D
                )
            else:
                grad_C_history = _history_part_point_3d_numba(
                    target_surface_nodes, n_history_steps,
                    np.array(source_body.location_history[:n_history_steps]),
                    t_now, dt, D
                )

        # --- Local Part ---
        if step > 0:
            if is_source_janus:
                grad_C_local = _local_part_dist_3d_numba_optimized(
                    target_surface_nodes, source_body.location, source_body.orientation,
                    source_body.location_history[step - 1], source_body.nodes_body_frame,
                    source_body.sigma_distribution, dt, D
                )
            else:
                grad_C_local = _local_part_point_3d_numba(
                    target_surface_nodes, source_body.location,
                    source_body.location_history[step - 1], dt, D
                )

        total_grad_on_target_nodes += grad_C_history + grad_C_local

    # --- Net Force & Torque Calculation (Lab Frame) ---
    r_vectors_from_center = target_surface_nodes - target_body.location
    norm_r = np.linalg.norm(r_vectors_from_center, axis=1, keepdims=True)
    norm_r[norm_r < 1e-12] = 1.0
    normals = r_vectors_from_center / norm_r

    grad_dot_norm = np.sum(total_grad_on_target_nodes * normals, axis=1, keepdims=True)
    grad_normal_component = grad_dot_norm * normals
    tangential_grad = total_grad_on_target_nodes - grad_normal_component

    # Surface integral approximations
    mean_tangential_grad_force = np.mean(tangential_grad, axis=0)

    # 3D Torque = r_lab x F_lab
    torque_vectors = np.cross(r_vectors_from_center, tangential_grad)
    mean_torque = np.mean(torque_vectors, axis=0)

    return mean_tangential_grad_force, mean_torque
