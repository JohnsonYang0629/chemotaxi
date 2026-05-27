import numpy as np
import os
import matplotlib.pyplot as plt


def smooth_step_1d(x, x_min, x_max, k):
    """
    Generate a 1D smooth step function using tanh
    x_min, x_max: region boundaries (must satisfy x_min < x_max)
    k: smoothness coefficient (larger = steeper transition)
    """
    return 0.5 * (np.tanh(k * (x - x_min)) + np.tanh(k * (x_max - x)))


def generate_2d_mesh(N_b, regions_theta, k, output_prefix="mesh_2d"):
    """
    Generate 2D circular mesh and chemical field (0 to 1)
    """
    d_theta = 2 * np.pi / N_b
    theta = d_theta * (np.arange(N_b) + 0.5) - np.pi

    nodes = np.zeros((N_b, 2))
    nodes[:, 0] = np.sin(theta)  # X = sin(theta), 0 degree at top
    nodes[:, 1] = np.cos(theta)  # Y = cos(theta)

    sigma = np.zeros(N_b)
    for (t_min, t_max) in regions_theta:
        t_mid = (t_min + t_max) / 2.0
        theta_eval = (theta - t_mid + np.pi) % (2 * np.pi) - np.pi + t_mid
        sigma += smooth_step_1d(theta_eval, t_min, t_max, k)

    sigma = np.clip(sigma, 0.0, 1.0)

    with open(f"{output_prefix}.vertex", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{nodes[i, 0]:.15f} {nodes[i, 1]:.15f}\n")

    with open(f"{output_prefix}.chem_dist.dat", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{sigma[i]:.15f}\n")

    print(f"2D Mesh generated: {output_prefix} with {N_b} nodes.")


def generate_2d_mesh_signed(N_b, regions_theta, k, output_prefix="mesh_2d_signed"):
    """
    Generate 2D circular mesh and chemical field supporting -1, 0, 1
    """
    d_theta = 2 * np.pi / N_b
    theta = d_theta * (np.arange(N_b) + 0.5) - np.pi

    nodes = np.zeros((N_b, 2))
    nodes[:, 0] = np.sin(theta)
    nodes[:, 1] = np.cos(theta)

    sigma = np.zeros(N_b)
    for (t_min, t_max, val) in regions_theta:
        t_mid = (t_min + t_max) / 2.0
        theta_eval = (theta - t_mid + np.pi) % (2 * np.pi) - np.pi + t_mid
        mask = smooth_step_1d(theta_eval, t_min, t_max, k)
        sigma += mask * val

    sigma = np.clip(sigma, -1.0, 1.0)

    with open(f"{output_prefix}.vertex", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{nodes[i, 0]:.15f} {nodes[i, 1]:.15f}\n")

    with open(f"{output_prefix}.chem_dist.dat", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{sigma[i]:.15f}\n")

    print(f"2D Signed Mesh generated: {output_prefix} with {N_b} nodes.")


# =====================================================================
# 新增: 生成迁移率 (Mobility) 分布的模块
# =====================================================================
def generate_2d_mobility_dist(N_b, regions_theta, k, default_alpha=1.0, output_prefix="mesh_2d"):
    """
    Generate 2D spatial mobility (alpha) distribution file.
    Args:
        N_b: Number of boundary nodes.
        regions_theta: List of tuples (theta_min, theta_max, target_alpha).
        k: Smoothness coefficient.
        default_alpha: The background mobility for areas not specified in regions_theta.
    """
    d_theta = 2 * np.pi / N_b
    theta = d_theta * (np.arange(N_b) + 0.5) - np.pi

    # Initialize the array with the default background mobility
    alpha_array = np.full(N_b, default_alpha, dtype=float)

    for (t_min, t_max, target_alpha) in regions_theta:
        t_mid = (t_min + t_max) / 2.0
        theta_eval = (theta - t_mid + np.pi) % (2 * np.pi) - np.pi + t_mid

        # Get the smooth mask (0 to 1) for this specific region
        mask = smooth_step_1d(theta_eval, t_min, t_max, k)

        # Smooth interpolation:
        # When mask = 0 (outside), it keeps its previous value.
        # When mask = 1 (inside), it becomes the target_alpha.
        alpha_array = alpha_array * (1.0 - mask) + target_alpha * mask

    # Write the mobility distribution to a .mob_dist.dat file
    with open(f"{output_prefix}.mob_dist.dat", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{alpha_array[i]:.15f}\n")

    print(f"2D Mobility Distribution generated: {output_prefix}.mob_dist.dat")


# =====================================================================


def visualize_2d_mesh_from_files(output_prefix):
    """
    Visualize the generated mesh, chemical field, and mobility field as separate figures.
    """
    vertex_file = f"{output_prefix}.vertex"
    chem_file = f"{output_prefix}.chem_dist.dat"
    mob_file = f"{output_prefix}.mob_dist.dat"

    if not os.path.exists(vertex_file) or not os.path.exists(chem_file):
        print(f"Error: Cannot find {vertex_file} or {chem_file}")
        return

    # 1. 读取节点坐标
    with open(vertex_file, "r") as f:
        lines = f.readlines()
        N_b = int(lines[0].strip())
        nodes = np.zeros((N_b, 2))
        for i in range(N_b):
            parts = lines[i + 1].split()
            nodes[i, 0] = float(parts[0])
            nodes[i, 1] = float(parts[1])

    # 2. 读取化学场分布 (Sigma)
    with open(chem_file, "r") as f:
        lines = f.readlines()
        sigma = np.zeros(N_b)
        for i in range(N_b):
            sigma[i] = float(lines[i + 1].strip())

    # 3. 读取迁移率分布 (Alpha)
    has_mob = os.path.exists(mob_file)
    if has_mob:
        with open(mob_file, "r") as f:
            lines = f.readlines()
            alpha = np.zeros(N_b)
            for i in range(N_b):
                alpha[i] = float(lines[i + 1].strip())

    theta_ref = np.linspace(-np.pi, np.pi, 100)

    # ================== 图 1: 化学场 (Sigma) 独立绘图 ==================
    max_abs_sigma = np.max(np.abs(sigma))
    if max_abs_sigma == 0:
        max_abs_sigma = 1.0  # 防止全零数据导致 vmin=vmax=0 报错

    plt.figure(figsize=(8, 6))
    plt.plot(np.cos(theta_ref), np.sin(theta_ref), 'k--', alpha=0.3, zorder=1)

    scatter1 = plt.scatter(
        nodes[:, 0], nodes[:, 1],
        c=sigma, cmap='coolwarm',
        vmin=-max_abs_sigma, vmax=max_abs_sigma,  # 自动提取上下限且 0 居中
        s=100, edgecolor='black', zorder=5
    )
    cbar1 = plt.colorbar(scatter1)
    cbar1.set_label('Release / Absorption Rate ($\sigma$)', fontsize=12)

    plt.axis('equal')
    plt.title(f'Chemical Source Dist ($\sigma$)\n({output_prefix})', fontsize=14)
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # 单独保存 chem_dist 图
    chem_vis_filename = f"{output_prefix}_chem_dist_visualization.png"
    plt.savefig(chem_vis_filename, dpi=300, bbox_inches='tight')
    print(f"Chemical distribution visualization saved to: {chem_vis_filename}")
    # plt.show()

    # ================== 图 2: 迁移率 (Alpha) 独立绘图 ==================
    if has_mob:
        max_abs_alpha = np.max(np.abs(alpha))
        if max_abs_alpha == 0:
            max_abs_alpha = 1.0

        plt.figure(figsize=(8, 6))
        plt.plot(np.cos(theta_ref), np.sin(theta_ref), 'k--', alpha=0.3, zorder=1)

        # 使用 PiYG (粉-绿) 区分物理场，或者改为 coolwarm
        scatter2 = plt.scatter(
            nodes[:, 0], nodes[:, 1],
            c=alpha, cmap='PiYG',
            vmin=-max_abs_alpha, vmax=max_abs_alpha,  # 自动提取上下限且 0 居中
            s=100, edgecolor='black', zorder=5
        )
        cbar2 = plt.colorbar(scatter2)
        cbar2.set_label('Mobility / Slip Coefficient ($\\alpha$)', fontsize=12)

        plt.axis('equal')
        plt.title(f'Mobility Dist ($\\alpha$)\n({output_prefix})', fontsize=14)
        plt.xlabel('X', fontsize=12)
        plt.ylabel('Y', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()

        # 单独保存 mob_dist 图
        mob_vis_filename = f"{output_prefix}_mob_dist_visualization.png"
        plt.savefig(mob_vis_filename, dpi=300, bbox_inches='tight')
        print(f"Mobility distribution visualization saved to: {mob_vis_filename}")
        # plt.show()


if __name__ == "__main__":
    # ===== 配置参数 =====
    N_2D = 40
    k_2D = N_2D / (2 * np.pi) * 1.5
    #prefix = "circle_R_1_N40_chiral_rotor"

    # 1. 定义表面化学释放分布 (Sigma)
    # 例如：标准的上下对称 Janus 粒子，上半圆释放，下半圆吸收
    #signed_regions_2d = [
    #    (-np.pi / 2, np.pi / 2, 1.0),  # 上半圆 (Y > 0) 为释放源 +1
    #    (np.pi / 2, 3 * np.pi / 2, -1.0)  # 下半圆 (Y < 0) 为吸收汇 -1
    #]

    # 2. 生成网格与化学场文件 (.vertex 和 .chem_dist.dat)
    #generate_2d_mesh_signed(
    #    N_b=N_2D,
    #    regions_theta=signed_regions_2d,
    #    k=k_2D,
    #    output_prefix=prefix
    #)

    # 3. 定义表面物理迁移率分布 (Alpha)
    # 核心魔法：左右不对称的表面材质！
    # 右半边材料 alpha = 1.0 (避化性，往推力反方向游)
    # 左半边材料 alpha = -1.0 (趋化性，往推力同方向游)
    #mobility_regions_2d = [
    #    (0.0, np.pi, 1.0),  # 右半圆 (X > 0) alpha = 1.0
    #    (-np.pi, 0.0, -1.0)  # 左半圆 (X < 0) alpha = -1.0
    #]

    # 4. 生成迁移率分布文件 (.mob_dist.dat)
    # 设定 default_alpha=0.0，如果在区域外则没有滑移。由于我们覆盖了全圆，default 实际上不会出现。
    #generate_2d_mobility_dist(
    #    N_b=N_2D,
    #    regions_theta=mobility_regions_2d,
    #    k=k_2D,
    #    default_alpha=0.0,
    #    output_prefix=prefix
    #)

    prefix = "circle_R_1_N40_slingshot_3_3"

    # 化学场 (Sigma)：头部释放，但左前强 (+1.0)，右前弱 (+0.6) -> 打破对称性
    signed_regions_2d = [
        (-np.pi / 2, 0.0, 4.0),  # 左前 1/4 圆
        (0.0, np.pi / 2, 0.2)  # 右前 1/4 圆
    ]

    # 物理场 (Alpha)：头部趋化 (-1.0，喜欢产物)，尾部避化 (保持默认 +1.5，讨厌产物)
    mobility_regions_2d = [
        (-np.pi / 2, np.pi / 2, -4.0)  # 前半圆趋化
    ]

    generate_2d_mesh_signed(N_2D, signed_regions_2d, k_2D, output_prefix=prefix)
    generate_2d_mobility_dist(N_2D, mobility_regions_2d, k_2D, default_alpha=6.0, output_prefix=prefix)
    # 5. 自动可视化化学场
    visualize_2d_mesh_from_files(prefix)