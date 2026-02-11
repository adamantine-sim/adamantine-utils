import numpy as np
import os
import glob
import pyvista
import time
import cv2

def get_iteration_count(path_prefix):
    filename_pattern = path_prefix + '*.pvtu'
    files = glob.glob(filename_pattern)

    iteration_numbers = []
    for f in files:
        temp = f.replace('.pvtu', '')
        post_id_index = temp.rfind('.')
        iteration_string = temp[post_id_index+1:]
        iteration_count = int(iteration_string)
        #max_iteration_count = max(max_iteration_count, iteration_count)
        iteration_numbers.append(iteration_count)

    iteration_numbers.sort()

    return iteration_numbers

def get_melt_pool_time_series(path_to_adamantine_files, adamantine_filename):
    this_file_path = os.path.dirname(os.path.realpath(__file__))

    pl = pyvista.Plotter()

    # Auto-detect the number of MPI domains (should be generalized across workflow-scripts)
    filename_pattern = path_to_adamantine_files + adamantine_filename + '.*.*.vtu'  
    list_of_mpi_ranks = glob.glob(filename_pattern)
    raw_ranks = [int(os.path.basename(f).split('.')[2]) for f in list_of_mpi_ranks]
    mpi_domains = max(raw_ranks) + 1

    # Get the list of time steps and sort them (should be generalized across workflow-scripts)             
    sim_filename = adamantine_filename
    iteration_numbers = get_iteration_count(path_to_adamantine_files + sim_filename)

    # Cycle through the time steps
    cycle_index = []
    depth = []
    length = []
    width = []
    time = []

    for iteration_number in iteration_numbers:
        file_to_plot = path_to_adamantine_files + sim_filename + '.' + str(iteration_number) + '.pvtu'

        dataset = pyvista.read(file_to_plot)
        dataset.set_active_scalars("temperature")

        isovalue = 1670.0

        clipped_volume = dataset.clip_scalar(scalars="temperature", value=isovalue, invert=False)

        if len(clipped_volume.points) > 0:

            total_volume = clipped_volume.volume

            cv = pl.add_mesh(clipped_volume, show_edges=True, cmap='plasma')

            bounds = clipped_volume.bounds

            melt_pool_depth = bounds[4]-bounds[5]

            points = clipped_volume.points

            # Strip out z position from the points
            points_2d = np.zeros([len(points), 2])
            i = 0
            for pt in points:
                points_2d[i][0] = pt[0]
                points_2d[i][1] = pt[1]
                i = i+1

            points_2d = points_2d.astype(np.float32).reshape(-1, 1, 2)

            ellipse = cv2.fitEllipse(points_2d)
            (center_x, center_y), (major_axis, minor_axis), angle = ellipse

            length_val = np.max([major_axis, minor_axis])
            width_val = np.min([major_axis, minor_axis])

            cycle_index.append(iteration_number)
            depth.append(np.abs(melt_pool_depth))
            width.append(width_val)
            length.append(length_val)

        else:
            #print("zero points")

            cycle_index.append(iteration_number)
            depth.append(0.0)
            width.append(0.0)
            length.append(0.0)

        # Get the time elapsed from the VTK file
        f_path = path_to_adamantine_files + adamantine_filename + '.' + str(iteration_number) + '.0.vtu'
        with open(f_path, 'r') as file:
            for i, line in enumerate(file):
                if i == 8:  
                    start_tag_end = line.find('>') + 1  # Position right after the opening tag
                    end_tag_start = line.find('</', start_tag_end)  # Position of the closing tag
                    time_entry = float(line[start_tag_end:end_tag_start].strip())
                    break
        time.append(time_entry)

    return (depth, width, length, time)


def melt_pool_statistics(path_to_adamantine_files, adamantine_filename):
    (depth, width, length, time) = get_melt_pool_time_series(path_to_adamantine_files, adamantine_filename)

    filtered_depth = [value for value in depth if value > 1e-8]
    filtered_width = [value for value in width if value > 1e-8]
    filtered_length = [value for value in length if value > 1e-8]

    melt_pool_stats = {}
    melt_pool_stats['average_depth'] = np.mean(filtered_depth)
    melt_pool_stats['average_width'] = np.mean(filtered_width)
    melt_pool_stats['average_length'] = np.mean(filtered_length)
    melt_pool_stats['std_depth'] = np.std(filtered_depth)
    melt_pool_stats['std_width'] = np.std(filtered_width)
    melt_pool_stats['std_length'] = np.std(filtered_length)

    return melt_pool_stats


def melt_pool_analysis(path_to_adamantine_files, adamantine_filename, output_directory):
    import matplotlib
    import matplotlib.pyplot as plt

    print("Performing melt pool analysis")

    # Some hard-coded variables that we may want to open up to users
    time_series_plot_filename = output_directory + 'melt_pool_ts_plot.png'
    this_file_path = os.path.dirname(os.path.realpath(__file__))

    # ----------------------------------------------------------
    # Plotting parameters
    # ----------------------------------------------------------

    SMALL_SIZE = 12
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fWidth = 12
    fHeight = 7.5
    
    (depth, width, length, time) = get_melt_pool_time_series(path_to_adamantine_files, adamantine_filename)
    
    fig, ax = plt.subplots()
    plt.plot(time, depth, 'b-',label='depth')
    plt.plot(time, width, 'g-', label='width')
    plt.plot(time, length, 'r-', label='length')

    plt.xlabel("Time (s)")
    plt.ylabel("Dimensions (m)")
    plt.legend()
    plt.show()






    