import numpy as np
import math as m
import os
import random
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt
#from mayavi.mlab import *
#from tool_storage_yz import *
#import sympy
#from sympy.abc import x,y
import shutil
from shutil import copyfile


def file_reader_2array_1free(filename: str, filepath: str):
    # erase first line
    # read file to store in array
    os.chdir(filepath)
    data_raw = open(filename, mode='r')
    line = data_raw.readline()
    data_list = []
    while line:
        num = list(map(float, line.split()))
        data_list.append(num)
        line = data_raw.readline()
    data_raw.close()

    data_array = np.zeros((len(data_list) - 1, len(data_list[1])))
    for i in range(0, len(data_list) - 1):
        data_array[i, :] = np.array(data_list[i+1])

    return data_array


def file_reader_2array(filename: str, filepath: str):
    # full pack
    # read file into array
    os.chdir(filepath)

    data_raw = open(filename, mode='r')
    line = data_raw.readline()
    data_list = []
    while line:
        num = list(map(float, line.split()))
        data_list.append(num)
        line = data_raw.readline()
    data_raw.close()

    data_array = np.zeros((len(data_list), len(data_list[0])))
    for i in range(0, len(data_list)):
        data_array[i, :] = np.array(data_list[i])

    return data_array


def file_writer_locate(filename: str, filepath: str, input_value: str, location: int):
    # input_value: string type
    # location: int type
    os.chdir(filepath)
    file = open(filename)
    lines = file.readlines()
    file.close()
    lines[location] = input_value
    file = open(filename, mode='w')
    file.writelines(lines)
    file.close()


def file_writer_4array(filename: str, filepath: str, input_value, number: int, blob_radius=0):
    # input_value: numpy array type
    os.chdir(filepath)
    length = np.shape(input_value)[0]
    if number == 0:
        file = open(filename, mode='w')
        file.write(str(1) + '\n')
        str_value = ' '.join(str(j) for j in input_value)
        line = str_value + '\n'
        file.write(line)
        file.close()
    elif number == -1:
        file = open(filename, mode='w')
        file.write(str(1) + '\n')
        file.close()
        file = open(filename, mode='a')
        str_value = ' '.join(str(j) for j in input_value)
        line = str_value + '\n'
        file.write(line)
        file.close()
    elif number == -2:
        file = open(filename, mode='w')
        for i in range(0, length):
            str_value = ' '.join(str(j) for j in input_value[i, :])
            line = str_value + '\n'
            file.write(line)
        file.close()
    elif number == -3:
        # flush the designated file
        file = open(filename, mode='w')
        file.write(str(''))
        file.close()

    else:
        file = open(filename, mode='w')
        file.write(str(number) + '\t' + str(blob_radius) + '\n')
        file.close()
        for i in range(0, length):
            file = open(filename, mode='a')
            str_value = ' '.join(str(j) for j in input_value[i, :])
            line = str_value + '\n'
            file.write(line)
            file.close()


def geo_dist_min(geo_array: list):
    # Euclidean metric
    list_length = np.shape(geo_array)[0]
    vec_dist = np.zeros(list_length-1)
    vec1 = geo_array[m.ceil((list_length-1)/2), :]
    for i in range(0, list_length-1):
        vec_dist[i] = np.linalg.norm(vec1-geo_array[i+1, :])
    vec_dist = np.sort(vec_dist)
    min_dist = vec_dist[1]
    return min_dist


def geo_dist_calculator(filename: str, filepath: str, array: list):
    # calculate the blob-blod distance with ratio r1, r2, r3, r4
    geo_data = file_reader_2array_1free(filename, filepath)
    dist = geo_dist_min(geo_data)
    dist_list = []
    for r in array:
        dist_list.append(r*dist)
    return dist_list


def text_reader(filename, filepath):
    os.chdir(filepath)
    str_file = []
    with open(filename, "r") as file:
        for line in file.readlines():
            str_file.append(line)

    return str_file


def text_reader_2array(filename: str, filepath: str):
    os.chdir(filepath)
    data_list = []
    with open(filename, "r") as file:
        for line in file.readlines():
            num = list(map(float, line.split()))
            data_list.append(num)

    return np.array(data_list)


def text_finder(filename, filepath, cue):
    # locate with given cue
    # output: cue-free string

    os.chdir(filepath)
    str_file = []
    count = 0
    with open(filename, "r") as file:
        for line in file.readlines():
            line = line.strip('\n')
            line = line.replace(' ', '')
            str_file.append(line)
            if cue in line:
                cue_pos = count
                break
            else:
                count = count + 1
    strip_cue_word = str_file[cue_pos].replace(cue, '')
    response = [count, strip_cue_word]
    return response


def text_replacer(filename, filepath, cue, input_value):
    # replace value after given cue

    os.chdir(filepath)
    pos = text_finder(filename, filepath, cue)[0]
    if cue == "structure1":
        cue = "structure"
    elif cue == "#SBATCH_job":
        cue = "#SBATCH -J"
    elif cue == "#SBATCH_core":
        cue = "#SBATCH -c"
    file_line = text_reader(filename, filepath)
    if type(input_value) == str:
        replacement = cue + '\t\t\t\t\t\t\t\t\t\t' + input_value + '\n'
    elif type(input_value) == float or type(input_value) == int:
        replacement = cue + '\t\t\t\t\t\t\t\t\t\t' + str(input_value) + '\n'
    else:
        replacement = cue + '\t\t\t\t\t\t\t\t\t\t' + 'error' + '\n'
        print('Error input value type')

    file_line[pos] = replacement
    with open(filename, "w") as file:
        for line in file_line:
            file.write(line)
    # print("input file rewritten successfully")


def file_switcher(scheme, filename, filepath):
    # switch input/output file names for different scheme

    os.chdir(filepath)
    file_name = text_finder(filename, filepath, 'output_name')[1]
    if scheme == "resistance":
        input_data_file_name = 'velocity_input.dat'
        output_data_file_name = file_name + '.force.dat'
    elif scheme == 'mobility':
        input_data_file_name = 'force_input.dat'
        output_data_file_name = file_name + '.velocity.dat'
    else:
        print('Wrong scheme')

    geo_file_name_raw = text_finder(filename, filepath, 'structure')[1]
    geo_file_name = geo_file_name_raw[1:-15]
    clone_file_name = 'squirmer.clones'
    save_file_name = file_name + '.result.dat'
    file_names = [input_data_file_name, output_data_file_name, geo_file_name, clone_file_name, save_file_name]
    return file_names


def generate_ellipsoid_spiral(filepath, num_pts):
    indices = np.arange(0, num_pts, dtype=float) + 0.5

    # 球坐标
    phi = np.arccos(1 - 2 * indices / num_pts)
    theta = np.pi * (1 + 5 ** 0.5) * indices

    # 参数化
    z, y, x = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), 2.0 * np.cos(phi)
    geo_array = np.concatenate((x, y, z), axis=0).reshape((3, num_pts)).transpose()

    geo_file_name = "ellipsoid_spiral_N_" + str(num_pts)+ "_abc_211.vertex"
    file_writer_4array(geo_file_name, filepath, geo_array, num_pts)


def read_input_file(filepath: str, filename: str):
    os.chdir(filepath)
    f = open(filename, 'r', errors='ignore')
    lines = f.readlines()
    dicts = {'blob_ratio': 0.0,
             'height_above_wall': 0.0,
             'test_input': 0.0,
             'clone_distance': 0.0,
             'clone_distribution': 0.0,
             'oscillating_frequency': 0.0,
             'oscillating_amplitude_A': 0.0,
             'magnetic_gamma': 0.0,
             'file_path': 0.0,
             'input_file': 0.0
             }

    dict_flag = {'blob_ratio': 1,
                 'height_above_wall': 1,
                 'test_input': 1,
                 'clone_distance': 1,
                 'clone_distribution': 2,
                 'oscillating_frequency': 1,
                 'oscillating_amplitude_A': 1,
                 'magnetic_gamma': 1,
                 'file_path': 0,
                 'input_file': 0
                 }

    for key in dicts:
        for line in lines:
            if key in line:
                line = line.replace(key, '')
                line = line.replace(' ', '')
                line = line.replace('\n', '')
                line = line.replace('\t', '')
                if dict_flag[key] == 1:
                    line_list = line.split(',')
                    line_num_list = list(map(float, line_list))
                    dicts[key] = line_num_list
                elif dict_flag[key] == 0:
                    dicts[key] = line
                elif dict_flag[key] == 2:
                    line_list = line.split(',')
                    line_num_list = list(map(int, line_list))
                    dicts[key] = line_num_list
                break

    return dicts


def create_clone_array(clone_distance: list, clone_distribution: list):

    col_dist, row_dist = clone_distance[0], clone_distance[1]
    col_num, row_num = clone_distribution[0], clone_distribution[1]
    clone_array = np.zeros((col_num*row_num, 7))

    counter = 0
    for idx_col in range(col_num):
        for idx_row in range(row_num):
            clone_array[counter, :4] = idx_col*np.array([col_dist, 0, 0, 0]) + \
                                       idx_row*np.array([0, row_dist, 0, 0]) + \
                                       np.array([0, 0, 1, 1])
            counter = counter+1
    return clone_array


def create_clone_quaternion(pos: list, quaternion: list):
    clone_pos_list = pos
    clone_q_list = quaternion
    clone_quaternion_list = clone_pos_list + clone_q_list
    clone_q_array = np.array(clone_quaternion_list)

    return clone_q_array


def div_green(r: float, theta: float, phi: float, theta0: float, phi0: float):
    """
    div_green: derivation from Green_Banff
    """

    div1 = 0
    div2 = (-m.cos(theta0) * m.sin(theta) + m.cos(phi - phi0) * m.cos(theta) * m.sin(theta0)) /\
           (4 * m.pi * r * (-1 + m.cos(theta) * m.cos(theta0) + m.cos(phi - phi0) * m.sin(theta) * m.sin(theta0)))
    div3 = (-m.sin(phi - phi0) * m.sin(theta0)) /\
           (4 * m.pi * r * (-1 + m.cos(theta) * m.cos(theta0) + m.cos(phi - phi0) * m.sin(theta) * m.sin(theta0)))

    div = np.array([div1, div2, div3])
    if len(div) < 3:
        div = np.array([0, 0, 0])
    else:
        div = div
    return div


def sph2cart_field(div: np.ndarray, r: float, theta: float, phi: float):
    """
    Convert spherical field matrix to cartesian coordinate
    """

    transformmatrix = np.array([[m.sin(theta) * m.cos(phi), m.cos(theta) * m.cos(phi), -m.sin(phi)],
                               [m.sin(theta) * m.cos(phi), m.cos(theta) * m.sin(phi), m.cos(phi)],
                               [m.cos(theta), -m.sin(theta), 0]])
    cart_field = r * np.dot(transformmatrix, div)

    cart_field = np.nan_to_num(cart_field, nan=0)
    return cart_field


def sph2cart_phy(r: float, theta: float, phi: float):
    """
    Convert spherical vector to cartesian vector using conventional physical notations
    """

    x = r * m.sin(theta) * m.cos(phi)
    y = r * m.sin(theta) * m.sin(phi)
    z = r * m.cos(theta)
    cart_vector = np.array([x, y, z])

    return cart_vector


def cart2sph_phy(x: float, y: float, z: float):
    """
    Convert cartesian vector to spherical vector using conventional physical notations
    """
    r = m.sqrt(m.pow(x, 2) + m.pow(y, 2) + m.pow(z, 2))
    theta = m.acos(z / r)
    phi = m.atan2(y, x)
    sph_vector = np.array([r, theta, phi])

    return sph_vector



def create_multidipoles(loc_vector: np.ndarray, pole_delta: float,
                        dipole_orientation: np.ndarray, blob_sph_loc: np.ndarray):
    """
    this function is used to manipulate multiple dipoles locating and placing.
    Required function: div_green
    ******
    loc_vector: [theta1,phi1;theta2,phi2;...] in nX2 list form.
    pole_delta: float, distance between two monopoles.
    dipole_orientation: normalized direction vector, determine the orientation of
                        the dipole postive direction, [1,0;1,0] for example.
    blob_sph_loc: recommend using function cart2sph_phy to preprocess the blob locations in Cartesian coordinate.
    """

    dipole_num = np.shape(loc_vector)[0]
    blob_slip_velocity = np.zeros([dipole_num, 3])

    for i in range(0, dipole_num):

        r = blob_sph_loc[0]
        theta = blob_sph_loc[1]
        phi = blob_sph_loc[2]

        theta_P0 = loc_vector[i, 0] + dipole_orientation[i, 0] * pole_delta     # source
        phi_P0 = loc_vector[i, 1] + dipole_orientation[i, 1] * pole_delta       # source
        theta_N0 = loc_vector[i, 0] - dipole_orientation[i, 0] * pole_delta     # sink
        phi_N0 = loc_vector[i, 1] - dipole_orientation[i, 1] * pole_delta       # sink

        derivative_P = div_green(r, theta, phi, theta_P0, phi_P0)               # source
        derivative_N = -div_green(r, theta, phi, theta_N0, phi_N0)              # sink

        blob_slip_velocity[i, :] = derivative_P + derivative_N

    if dipole_num == 1:
        total_blob_slip = blob_slip_velocity
    else:
        total_blob_slip = np.sum(blob_slip_velocity, axis=0)

    return total_blob_slip


def regularize_dipole_center(loc_vector: np.ndarray, blob_cart_loc: np.ndarray, control_length: float, enhancement=1):
    """
    manipulate the numerical divergence at the dipole center and surrounding slip boundaries
    ******
    Required function: sph2cart_phy
    ******
    loc_vector: the location of the dipoles centers (numpy array)
    control_length: control radius (float, in arc length)
    blob_cart_loc: the blob position in Cartesian coordinate (numpy array)
    """

    dipole_num = np.shape(loc_vector)[0]
    flag_list = np.zeros((dipole_num, 1))

    for i in range(0, dipole_num):
        radius = np.linalg.norm(blob_cart_loc)
        dipole_cart_loc = sph2cart_phy(radius, loc_vector[i, 0], loc_vector[i, 1])

        p1 = dipole_cart_loc
        p2 = blob_cart_loc
        arc_dist = radius * m.atan2(np.linalg.norm(np.cross(p1, p2)), np.dot(p1, p2))

        if arc_dist <= control_length:
            flag_list[i, 0] = enhancement * (arc_dist / control_length)
        else:
            flag_list[i, 0] = 1

    flag = np.min(flag_list)

    return flag


def get_quaternion_from_euler(roll, pitch, yaw):
    """
    Convert an Euler angle to a quaternion.

    Input
      :param roll: The roll (rotation around x-axis) angle in degrees.
      :param pitch: The pitch (rotation around y-axis) angle in degrees.
      :param yaw: The yaw (rotation around z-axis) angle in degrees.

    Output
      :return qx, qy, qz, qw: The orientation in quaternion [w,x,y,z] format
    """
    roll = roll * m.pi / 180
    pitch = pitch * m.pi / 180
    yaw = yaw * m.pi / 180

    qx = np.sin(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) - np.cos(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)
    qy = np.cos(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2)
    qz = np.cos(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2) - np.sin(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2)
    qw = np.cos(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)

    return [qw, qx, qy, qz]


def create_directory(directory_upper_path, directory_name, overwrite: bool):
    os.chdir(directory_upper_path)
    try:
        os.mkdir(directory_name)
    except FileExistsError:
        if overwrite:
            print("Previous version of [" + directory_name + "] has been overwritten.")
            shutil.rmtree(directory_name)
            os.mkdir(directory_name)
        else:
            print("A version of [" + directory_name + "] already existed.")
    except:
        print(" Directory [" + directory_name + "] creation failed.")


def create_zip_file(directory_upper_path, output_filename, overwrite: bool):
    os.chdir(directory_upper_path)
    try:
        shutil.make_archive(output_filename, 'zip', output_filename)
    except FileExistsError:
        if overwrite:
            print("Previous version of [" + output_filename + "] has been overwritten.")
            shutil.rmtree(output_filename)
            shutil.make_archive(output_filename, 'zip', output_filename)
        else:
            print("A version of [" + output_filename + "] already existed.")
    except:
        print(" Zip file [" + output_filename + "] creation failed.")


def create_parameter_tag(parameter: np.array):
    parameter_str = np.array(np.zeros(np.shape(parameter)), dtype=str)
    for i in range(0, np.size(parameter)):
        value = np.abs(parameter[i])
        int_str = str(int(m.floor(value)))
        decimal_str = str(int(round((value - m.floor(value)), 2) * 100))

        if decimal_str == "0":
            decimal_str = "00"
        if len(int_str) < 2:
            int_str = '00' + int_str
        elif len(int_str) == 2:
            int_str = '0' + int_str
        if len(decimal_str) < 2:
            decimal_str = '0' + decimal_str

        if parameter[i] >= 0:
            parameter_str[i] = int_str + 'd' + decimal_str
        else:
            parameter_str[i] = 'N' + int_str + 'd' + decimal_str
    return parameter_str


# test zone (not for retrieving tools)
if __name__ == '__main__':

    file_path = "/Users/johnson/Desktop/RigidMultiblobsWall-master/inputs"
    os.chdir(file_path)

    src_file_name = "inputfile_spring_rotor.dat"

    chi_list = np.linspace(5, 85, 9, True)

    #h_list = np.logspace(-1, 1, 10, True)
    h_list = np.linspace(1.5, 1.9, 5, True)
    h_str_list = []
    for i in range(0, 5):
        int_str = str(m.floor(h_list[i]))
        decimal_str = str(int(round((h_list[i] - m.floor(h_list[i])), 2) * 100))
        if decimal_str == "0":
            decimal_str = "00"
        h_str_list.append(int_str + 'd' + decimal_str)

    K_list = np.logspace(-1, 2, 10)
    K_str_list = []
    for i in range(0, 10):
        int_str = str(int(m.floor(K_list[i])))
        decimal_str = str(int(round((K_list[i] - m.floor(K_list[i])), 2) * 100))
        if decimal_str == "0":
            decimal_str = "00"
        K_str_list.append(int_str + 'd' + decimal_str)

    #K_list = [0.1, 0.46415888, 2.15443469, 10.0, 46.41588834]
    K_list = [10.00]
    #h_list = [1.16681005, 1.46415888, 2.29154967, 4.59381366, 11.0]
    #K_str_list = ["0d10", "0d46", "2d15", "10d00", "46d42"]
    K_str_list = ["10d00"]
    #h_str_list = ["1d17", "1d46", "2d29", "4d59", "11d00"]

    for i in range(0, 1):
        for j in range(0, 5):
            for k in range(0, 9):
                new_file_name = "K_" + K_str_list[i] + "_h_" + h_str_list[j] + "_chi_" + str(round(chi_list[k])) + ".dat"
                copyfile(src_file_name, new_file_name)

                if i == 4:
                    text_replacer(new_file_name, file_path, "dt", "0.25")
                    text_replacer(new_file_name, file_path, "n_steps", "200000")


                file_ID = "K_" + K_str_list[i] + "_h_" + h_str_list[j] + "_chi_" + str(round(chi_list[k]))
                text_replacer(new_file_name, file_path, "output_name", "shell_N1_SD_1P_S50_" + file_ID)

                initial_cg = [0, 0, h_list[j]]
                str_initial_cg = ' '.join(str(p) for p in initial_cg)
                text_replacer(new_file_name, file_path, "initial_cg", str_initial_cg)
                str_K = str(K_list[i])
                text_replacer(new_file_name, file_path, "lambda_spring", str_K)

                clone_name = file_ID + ".clones"

                str_clone = "../../Structures/shell_N_218_Rg_1_s_0.2095_cube2sphere.vertex  " + clone_name + \
                            "  ../../Slips/dipole_P_1_S_50_at_shell_N_218_Rg_1_cube2sphere.slip"
                text_replacer(new_file_name, file_path, "structure", str_clone)

                q = get_quaternion_from_euler(chi_list[k], 0, 0)
                cq = create_clone_quaternion(h_list[j], q)
                file_writer_4array(clone_name, file_path, cq, -1, 0.05)


    #file_writer_4array(clone_name, file_path, cq, 2, 0.03)
