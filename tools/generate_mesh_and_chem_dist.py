import numpy as np
import os
import matplotlib.pyplot as plt


def smooth_step_1d(x, x_min, x_max, k):
    """
    使用 tanh 生成 1D 平滑阶跃函数
    x_min, x_max: 区域边界 (必须满足 x_min < x_max)
    k: 平滑系数 (越大越陡峭，趋近于绝对阶跃)
    """
    # 当 x 处于 (x_min, x_max) 内部时，两项均为 +1，均值为 1
    # 当 x 在外部时，一正一负，均值为 0
    return 0.5 * (np.tanh(k * (x - x_min)) + np.tanh(k * (x_max - x)))


def generate_2d_mesh(N_b, regions_theta, k, output_prefix="mesh_2d"):
    """
    生成 2D 圆形网格和化学场
    【注意：已修改为 theta=0 对应正上方(0,1)，顺时针方向增加角度】
    """
    d_theta = 2 * np.pi / N_b
    theta = d_theta * (np.arange(N_b) + 0.5) - np.pi

    nodes = np.zeros((N_b, 2))
    # ====== 修改坐标映射 ======
    nodes[:, 0] = np.sin(theta)  # X = sin(theta)
    nodes[:, 1] = np.cos(theta)  # Y = cos(theta)
    # ==========================

    sigma = np.zeros(N_b)
    for (t_min, t_max) in regions_theta:
        # 使用周期映射确保完美对称
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
    生成 2D 圆形网格和化学场 (支持 -1, 0, 1)
    【注意：已修改为 theta=0 对应正上方(0,1)，顺时针方向增加角度】
    """
    d_theta = 2 * np.pi / N_b
    theta = d_theta * (np.arange(N_b) + 0.5) - np.pi

    nodes = np.zeros((N_b, 2))
    # ====== 修改坐标映射 ======
    nodes[:, 0] = np.sin(theta)  # X = sin(theta)
    nodes[:, 1] = np.cos(theta)  # Y = cos(theta)
    # ==========================

    sigma = np.zeros(N_b)
    for (t_min, t_max, val) in regions_theta:
        t_mid = (t_min + t_max) / 2.0
        # 完美的周期映射，防止上下/左右边界截断
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


def generate_3d_mesh(N_b, regions_3d, k, output_prefix="mesh_3d"):
    """
    生成 3D 球面网格和化学场
    regions_3d: 列表，包含多个 ((theta_min, theta_max), (phi_min, phi_max)) 元组
    theta (方位角) in [-pi, pi], phi (极角) in [0, pi]
    """
    # 1. 使用斐波那契晶格 (Fibonacci Lattice) 生成极其均匀的球面点云
    nodes = np.zeros((N_b, 3))
    theta = np.zeros(N_b)
    phi = np.zeros(N_b)

    golden_ratio = (1 + 5 ** 0.5) / 2
    for i in range(N_b):
        theta[i] = 2 * np.pi * i / golden_ratio
        # 将 theta 映射到 [-pi, pi] 方便对齐习惯
        theta[i] = (theta[i] + np.pi) % (2 * np.pi) - np.pi

        # phi 极角 [0, pi]，采用 (i+0.5) 避免正好落在南北极点
        phi[i] = np.arccos(1 - 2 * (i + 0.5) / N_b)

        nodes[i, 0] = np.sin(phi[i]) * np.cos(theta[i])
        nodes[i, 1] = np.sin(phi[i]) * np.sin(theta[i])
        nodes[i, 2] = np.cos(phi[i])

    # 2. 计算平滑化学场
    sigma = np.zeros(N_b)
    for (t_range, p_range) in regions_3d:
        t_min, t_max = t_range
        p_min, p_max = p_range

        theta_eval = np.where(theta < t_min - np.pi, theta + 2 * np.pi, theta)

        # 3D 表面掩码 = 方位角掩码 * 极角掩码
        mask_theta = smooth_step_1d(theta_eval, t_min, t_max, k)
        mask_phi = smooth_step_1d(phi, p_min, p_max, k)
        sigma += mask_theta * mask_phi

    sigma = np.clip(sigma, 0.0, 1.0)

    # 3. 写入文件
    with open(f"{output_prefix}.vertex", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{nodes[i, 0]:.15f} {nodes[i, 1]:.15f} {nodes[i, 2]:.15f}\n")

    with open(f"{output_prefix}.chem_dist.dat", "w") as f:
        f.write(f"{N_b}\n")
        for i in range(N_b):
            f.write(f"{sigma[i]:.15f}\n")

    print(f"3D Mesh generated: {output_prefix} with {N_b} nodes.")


def visualize_2d_mesh_from_files(output_prefix):
    """
    读取生成的网格和化学场文件，并进行可视化绘制。
    用圆点代表离散网格点，颜色代表释放/吸收速率 (-1 到 1)。
    """
    vertex_file = f"{output_prefix}.vertex"
    chem_file = f"{output_prefix}.chem_dist.dat"

    # 1. 检查文件是否存在
    if not os.path.exists(vertex_file) or not os.path.exists(chem_file):
        print(f"Error: 找不到文件 {vertex_file} 或 {chem_file}")
        return

    # 2. 读取节点坐标数据
    with open(vertex_file, "r") as f:
        lines = f.readlines()
        N_b = int(lines[0].strip())
        nodes = np.zeros((N_b, 2))
        for i in range(N_b):
            parts = lines[i + 1].split()
            nodes[i, 0] = float(parts[0])
            nodes[i, 1] = float(parts[1])

    # 3. 读取化学场释放率数据
    with open(chem_file, "r") as f:
        lines = f.readlines()
        sigma = np.zeros(N_b)
        for i in range(N_b):
            sigma[i] = float(lines[i + 1].strip())

    # 4. 开始绘图
    plt.figure(figsize=(8, 6))

    # 画一个浅色的虚线圆作为边界参考线
    theta_ref = np.linspace(-np.pi, np.pi, 100)
    plt.plot(np.cos(theta_ref), np.sin(theta_ref), 'k--', alpha=0.3, zorder=1)

    # 绘制散点图：c=sigma 指定颜色数据，cmap='coolwarm' 使用冷暖渐变，vmin和vmax固定颜色范围
    scatter = plt.scatter(
        nodes[:, 0], nodes[:, 1],
        c=sigma, cmap='coolwarm', vmin=-1.0, vmax=1.0,
        s=100, edgecolor='black', zorder=5
    )

    # 添加颜色条
    cbar = plt.colorbar(scatter)
    cbar.set_label('Release / Absorption Rate ($\sigma$)', fontsize=12)

    # 设置图表格式
    plt.axis('equal')  # 确保X轴和Y轴比例一致，圆不会被压扁
    plt.title(f'2D Mesh Point Source Distribution\n({output_prefix})', fontsize=14)
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # 保存并显示图片
    vis_filename = f"{output_prefix}_visualization.png"
    plt.savefig(vis_filename, dpi=300)
    print(f"Visualization saved to: {vis_filename}")
    plt.show()


if __name__ == "__main__":
    # ===== 配置参数 =====
    N_2D = 40  # 提高节点数以获得更好的空间解析度
    #N_3D = 200

    # 自动计算平滑陡峭度 k (经验公式：使过渡带大约占据 3-4 个网格的宽度)
    k_2D = N_2D / (2 * np.pi) * 1.5
    #k_3D = np.sqrt(N_3D) / np.pi * 1.5

    # ===== 2D 示例 =====

    # 示例 1: 东半圆 (East Semicircle, X > 0)
    #generate_2d_mesh(N_2D, regions_theta=[(-np.pi / 2, np.pi / 2)], k=k_2D, output_prefix="circle_R_1_N20_east")

    # 示例 2: 多个离散区域叠加
    #multi_regions_2d = [
    #    (-np.pi, -np.pi / 2),
    #    (np.pi / 2, 2*np.pi/3)
    #]
    #generate_2d_mesh(N_2D, regions_theta=multi_regions_2d, k=k_2D, output_prefix="circle_R_1_N40_asym_NW_90_SW_30")

    # 定义带有幅值的离散区域 (theta_min, theta_max, value)
    # 例如：北半球部分区域为 1，南半球部分区域为 -1，其余区域默认为 0
    signed_regions_2d = [
        (-np.pi / 6, np.pi / 6, 1.0),  # 头部 (正上方) 30度释放源
        (np.pi / 2, 5 * np.pi / 6, -1.0),  # 右侧下方吸收汇
        (7 * np.pi / 6, 3 * np.pi / 2, -1.0)  # 左侧下方吸收汇
    ]
    prefix = "circle_R_1_N40_asym_type_5"
    generate_2d_mesh_signed(
        N_b=N_2D,
        regions_theta=signed_regions_2d,
        k=k_2D,
        output_prefix=prefix
    )
    visualize_2d_mesh_from_files(prefix)
    # ===== 3D 示例 =====

    # 示例 3: 3D 赤道带 (Equatorial band) 且只在东半球释放
    # 方位角 theta 在东半球 (-pi/2 到 pi/2)
    # 极角 phi 靠近赤道 (pi/4 到 3pi/4)
    multi_regions_3d = [
        ((-np.pi / 2, np.pi / 2), (np.pi / 4, 3 * np.pi / 4))
    ]
    #generate_3d_mesh(N_3D, regions_3d=multi_regions_3d, k=k_3D, output_prefix="sphere_equator_N200")