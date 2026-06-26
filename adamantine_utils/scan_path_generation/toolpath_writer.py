import copy
import numpy as np
import re
import os
import shutil
import time
import json
import random
import pyvista as pv
from typing import List, Tuple
import logging
logging.basicConfig(level=logging.INFO)

# ------------------------------------------------------------
# FUNCTIONS TO READ ADAMANTINE INPUTS
# ------------------------------------------------------------
def get_adamantine_power(filename):
    f_original = open(filename, 'r')
    power = np.nan
    for line in f_original:
        split_line = line.split()
        if len(split_line) > 0:
            if split_line[0] == "max_power":
                power = float(split_line[1])
    return power

def get_adamantine_time_step(filename):
    f_original = open(filename, 'r')
    dt = np.nan
    for line in f_original:
        split_line = line.split()
        if len(split_line) > 0:
            if split_line[0] == "time_step":
                dt = float(split_line[1])
    return dt

def get_is_ensemble(filename):
    f_original = open(filename, 'r')
    is_ensemble = False
    for line in f_original:
        split_line = line.split()
        if len(split_line) > 0:
            if split_line[0] == "ensemble_simulation":
                is_ensemble = bool(split_line[1])
    return is_ensemble

def get_ensemble_count(filename):
    f_original = open(filename, 'r')
    count = None

    for line in f_original:
        split_line = line.split()
        if len(split_line) > 0:
            if split_line[0] == "ensemble_size":
                count = int(split_line[1])

    if count == None:
        logging.error("ERROR: Simulation template input file does not give the ensemble size.")
        sys.exit()

    return count

def get_adamantine_file_prefix(filename):
    f_original = open(filename, 'r')
    file_prefix = ""
    for line in f_original:
        split_line = line.split()
        if len(split_line) > 0:
            if split_line[0] == "filename_prefix":
                file_prefix = split_line[1]
    return file_prefix

def check_for_triggered_dwell(simulation_result_directory, output_file_prefix, expected_layer_end, trigger_type, trigger_temperature, sleep_time_inc = 2.0):
    # Search for all files that match the pattern for pvtu output
    last_adamantine_file = simulation_result_directory + "/" + output_file_prefix + "." + expected_layer_end + ".pvtu"
    file_check = True
    while (file_check):
        try:
            # Read the vtu files and get the max temperature from this time step
            adamantine_mesh = pv.read(last_adamantine_file)
            file_check = False
            if (trigger_type == "max"):
                # If max temperature is below the trigger, dwell is not triggered, otherwise dwell is triggered
                max_temperature = np.max(adamantine_mesh.point_data['temperature'])
                print("Max temperature for adamantine time step",expected_layer_end,"is",max_temperature,"K")
                if (max_temperature < trigger_temperature):
                    dwell = False
                else:
                    dwell = True
            elif (trigger_type == "mean"):
                # If mean temperature (excluding voxels that are inactive and have temperatures of zero) is below the trigger, dwell is not triggered, otherwise dwell is triggered
                clipped_mesh = adamantine_mesh.clip_scalar(scalars='temperature', value=0.0001, invert=False)
                mean_temperature = np.mean(clipped_mesh.point_data['temperature'])
                print("Mean temperature for adamantine time step",expected_layer_end,"is",mean_temperature,"K")
                if (mean_temperature < trigger_temperature):
                    dwell = False
                else:
                    dwell = True
        except:
            # File is not yet available to read - force code to wait and try again
            time.sleep(sleep_time_inc)
    return dwell

def get_time_position_power_inp(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    time_position_power = []

    for line in lines:
        if line != '\n':
            split_line = line.split(',')
            time = float(split_line[0])
            x = float(split_line[1])
            y = float(split_line[2])
            z = float(split_line[3])
            power = float(split_line[4])
            time_position_power.append((time, (x,y,z), power))

    return time_position_power

def extract_single_layer(time_position_power, target_layer):
    current_layer_height = None
    current_layer_index = 0
    current_layer_start_time = 0.0
    extracted_layer = []
    
    for entry in time_position_power:
        if entry[2] > 1e-3:
            z = entry[1][2]

            if current_layer_height == None:
                current_layer_height = z

            if (z - current_layer_height > 1e-3):
                current_layer_height = z
                current_layer_index = current_layer_index + 1
                current_layer_start_time = entry[0]

        if current_layer_index == target_layer:
            shifted_time = np.max([entry[0] - current_layer_start_time, 1e-12])
            shifted_entry = (shifted_time, (entry[1][0], entry[1][1], 0.0), entry[2])
            extracted_layer.append(shifted_entry)

    return extracted_layer

# ------------------------------------------------------------
# FUNCTIONS TO MODIFY ADAMANTINE INPUTS
# ------------------------------------------------------------
def set_adamantine_sim_end_time(filename, end_time):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("duration")
        if (starting_index > -1):
            ending_index = starting_index + len("duration")
            new_line = line[0:ending_index] + " " + str(end_time) + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    os.remove(filename)
    shutil.move("temp.info", filename)

def set_adamantine_sim_restart(filename, restart_file_prefix):
    f = open(filename, 'a')
    f.write("\n")
    f.write("restart\n")
    f.write("{\n")
    f.write("filename_prefix " + restart_file_prefix + "\n")
    f.write("}\n")

def set_adamantine_sim_checkpointing(filename, steps_between_checkpoints):
    f = open(filename, 'a')
    f.write("\n")
    f.write("checkpoint\n")
    f.write("{\n")
    f.write("time_steps_between_checkpoint " + str((int(steps_between_checkpoints))) + "\n")
    f.write("filename_prefix " + "cp" + "\n")
    f.write("overwrite_files " + "false" + "\n")

    f.write("}\n")

def set_adamantine_sim_power(filename, power):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("max_power")
        if (starting_index > -1):
            ending_index = starting_index + len("max_power")
            new_line = line[0:ending_index] + " " + str(power) + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    os.remove(filename)
    shutil.move("temp.info", filename)

def set_adamantine_no_ensemble(filename):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("ensemble_simulation")
        if (starting_index > -1):
            ending_index = starting_index + len("ensemble_simulation")
            new_line = line[0:ending_index] + " false" + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    os.remove(filename)
    shutil.move("temp.info", filename)

def set_adamantine_no_experiment(filename):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("read_in_experimental_data")
        if (starting_index > -1):
            ending_index = starting_index + len("read_in_experimental_data")
            new_line = line[0:ending_index] + " false" + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    os.remove(filename)
    shutil.move("temp.info", filename)

def set_adamantine_no_data_assimilation(filename):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("assimilate_data")
        if (starting_index > -1):
            ending_index = starting_index + len("assimilate_data")
            new_line = line[0:ending_index] + " false" + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    os.remove(filename)
    shutil.move("temp.info", filename)

def set_adamantine_no_scan_path_duration(filename):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("scan_path_for_duration")
        if (starting_index > -1):
            ending_index = starting_index + len("scan_path_for_duration")
            new_line = line[0:ending_index] + " false" + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    os.remove(filename)
    shutil.move("temp.info", filename)

def set_adamantine_sim_output_freq(filename, freq):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    for line in f_original:
        starting_index = line.find("time_steps_between_output")
        if (starting_index > -1):
            ending_index = starting_index + len("time_steps_between_output")
            new_line = line[0:ending_index] + " " + str(freq) + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)
    f_original.close()
    f_new.close()
    shutil.move("temp.info", filename)

# ------------------------------------------------------------
# FUNCTIONS TO MODIFY ADAMANTINE SCAN PATHS
# ------------------------------------------------------------
def time_position_power_dwell(start_time, position, duration):
    dwell_segment = [(start_time+duration, position, 0.0)]
    return dwell_segment

def fix_time_and_height(layer, start_time, start_height):
    updated_layer = ""
    for line in layer:
        split_line = line.split(',')
        updated_line = str(float(split_line[0]) + start_time) + ',' + split_line[1] + ',' + split_line[2] + ',' + str(float(split_line[3]) + start_height) + ',' + split_line[4]
        updated_layer = updated_layer + updated_line #+ "\n"

    return updated_layer

def fix_time_and_height_time_position_power(layer, start_time, start_height):
    current_start_time = layer[0][0][0]
    logging.info(f"{start_time}, {current_start_time}")
    time_shift = start_time - current_start_time + 1e-12

    current_start_height = layer[0][0][1][2]
    height_shift = start_height - current_start_height

    new_layer = []
    for entry in layer[0]:
        new_layer.append((entry[0]+time_shift, (entry[1][0], entry[1][1], entry[1][2]+height_shift), entry[2]))

    return new_layer

def fix_time_time_position_power(layer, start_time):
    current_start_time = layer[0][0][0]
    logging.info(f"{start_time}, {current_start_time}")
    time_shift = start_time - current_start_time + 1e-12

    new_layer = []
    for entry in layer[0]:
        new_layer.append((entry[0]+time_shift, (entry[1][0], entry[1][1], entry[1][2]), entry[2]))

    return new_layer

def get_layer_and_scan_labels(sorted_data, layer_limit):
    '''
    List the names of each scan path directory and associated layer number
    '''
    layer_list = []
    scan_list = []
    for layer in range(layer_limit):
        scan_names = sorted_data[layer]["components"]
        num_scans = len(scan_names)
        for scan in range(num_scans):
            layer_list.append(layer)
            scan_list.append(scan)
    return layer_list, scan_list

def drip_feed_toolpath_layers(base_filename, toolpath_info, toolpath_module, layer_index_to_add, chunk_size, layer_limit, sorted_data, simulation_template_dir, override_dwells=False, end_sim=False, extra_time_at_end=False):
    '''
    Adds a new layer or chunk of layers to the toolpath.
    '''
    get_time_position_power_inp = getattr(toolpath_module, 'get_time_position_power_inp')
    create_toolpath_layers = getattr(toolpath_module, 'create_toolpath_layers')
    write_event_series = getattr(toolpath_module, 'write_event_series')
    time_position_power_dwell = getattr(toolpath_module, 'time_position_power_dwell')

    if (layer_index_to_add == 0):
        old_scan_path = []
        start_time = 0.0
    else:
        old_scan_path = get_time_position_power_inp(base_filename)
        start_time = old_scan_path[-1][0]

    # Get the new layer from create_toolpath
    last_layer_to_add = min(layer_index_to_add + chunk_size, layer_limit)
    toolpath_info['selected_layers'] = [layer_index_to_add, last_layer_to_add]
    new_layer, actual_layer_variables = create_toolpath_layers(sorted_data, toolpath_info, simulation_template_dir, override_dwells)

    #shifted_new_layer = fix_time_and_height_time_position_power(new_layer, start_time, start_height)
    shifted_new_layer = fix_time_time_position_power([new_layer], start_time)

    new_toolpath = old_scan_path + shifted_new_layer

    # Get the final time (this is the real end of the layer, not the padded toolpath)
    end_time = new_toolpath[-1][0]

    # Pad the new layer with a short dwell to get around roundoff error for adamantine time steps
    if extra_time_at_end:
        extra = 1.0
        t_end = new_toolpath[-1][0]
        pos_end  = new_toolpath[-1][1]
        new_toolpath += time_position_power_dwell(t_end, pos_end, extra)

    # Add some buffer time if this is the end of the simulation (usually for a forward simulation)
    # NOTE: Adding this buffer time messes with parsing the end of the scan path for the time the last layer ends
    #if end_sim:
    #    buffer = 5.0
    #    new_layer_string = new_layer_string + str(end_time + buffer) + ",0.0,0.0,0.0,0.0\n"

    temp_filename = base_filename + '.tmp'

    write_event_series(new_toolpath, temp_filename, end_sim)

    os.rename(temp_filename, base_filename)
    print(end_time)
    return end_time, actual_layer_variables[0]

def drip_feed_toolpath_single_scan(base_filename, toolpath_info, toolpath_module, layer_list, scan_list, layer_scan_index, sorted_data, simulation_template_dir, override_dwells=False, end_sim=False, extra_time_at_end=False):
    '''
    Adds scan data from one file to the toolpath.
    '''
    get_time_position_power_inp = getattr(toolpath_module, 'get_time_position_power_inp')
    create_toolpath_single_scan = getattr(toolpath_module, 'create_toolpath_single_scan')
    write_event_series = getattr(toolpath_module, 'write_event_series')
    time_position_power_dwell = getattr(toolpath_module, 'time_position_power_dwell')
    num_layers = toolpath_info.get('num_layers')

    layer_idx = layer_list[layer_scan_index]
    scan_idx = scan_list[layer_scan_index]
    if (layer_scan_index == 0) or (layer_list[layer_scan_index-1] != layer_idx):
        first_scan_of_layer = True
    else:
        first_scan_of_layer = False

    if (layer_scan_index == len(layer_list)-1) or (layer_list[layer_scan_index+1] != layer_idx):
        last_scan_of_layer = True
    else:
        last_scan_of_layer = False

    # First scan of simulation - init new scan event list, otherwise reread existing events
    if (layer_scan_index == 0):
        old_scan_path = []
        scan_start_time = 0.0
    else:
        old_scan_path = get_time_position_power_inp(base_filename)
        scan_start_time = old_scan_path[-1][0]

    # Get print plan data for this layer and this scan path
    comp_name = sorted_data[layer_idx]["components"][scan_idx]
    plan_data = sorted_data[layer_idx][comp_name]
    
    # Unpack inputs
    layer_end_time_discretization = toolpath_info.get('layer_end_time_discretization', None)
    set_dwell_every_n_layers      = toolpath_info.get('set_dwell_every_n_layers')

    # Build the new time–position–power series
    new_tpp = []

    # quarters: 1/4, 2/4, 3/4 of the layer count
    q = num_layers // 4  # integer floor, stable
    quarter_layers_1based = {q*1, q*2, q*3}
            
    filename = simulation_template_dir + comp_name + '/layer_' + str(layer_idx) + '_scan_path.txt'

    # Scan path from one file (discretized to integer multiple of adamantine time step)
    new_scan_path, scan_end_time, additional = create_toolpath_single_scan(filename, plan_data, layer_idx, num_layers, scan_start_time, first_scan_of_layer, last_scan_of_layer, new_tpp, override_dwells, set_dwell_every_n_layers, layer_end_time_discretization, True)

    #shifted_new_layer = fix_time_and_height_time_position_power(new_layer, start_time, start_height)
    shifted_new_path = fix_time_time_position_power([new_scan_path], scan_start_time)
    new_toolpath = old_scan_path + shifted_new_path

    # Get the final time (this is the real end of the section, not the padded toolpath)
    scan_end_time = new_toolpath[-1][0]

    # Pad the new layer with a short dwell to get around roundoff error for adamantine time steps
    if extra_time_at_end:
        extra = 1.0
        t_end = new_toolpath[-1][0]
        pos_end  = new_toolpath[-1][1]
        new_toolpath += time_position_power_dwell(t_end, pos_end, extra)

    # Add some buffer time if this is the end of the simulation (usually for a forward simulation)
    # NOTE: Adding this buffer time messes with parsing the end of the scan path for the time the last layer ends
    #if end_sim:
    #    buffer = 5.0
    #    new_layer_string = new_layer_string + str(end_time + buffer) + ",0.0,0.0,0.0,0.0\n"

    temp_filename = base_filename + '.tmp'

    write_event_series(new_toolpath, temp_filename, end_sim)

    os.rename(temp_filename, base_filename)
    print(scan_end_time)

    return scan_end_time

def get_path_drip_first_segment(simulation_template_dir, target_filename, sorted_data, layer_list, scan_list, layer_scan_index, segment_start_time):

    # Layer and scan number, corresponding data for this scan path
    layer_idx = layer_list[layer_scan_index]
    scan_idx = scan_list[layer_scan_index]
    comp_name = sorted_data[0]["components"][scan_idx]
    plan_data = sorted_data[layer_idx][comp_name]

    # Scan path from this file
    scan_path_filename_read = simulation_template_dir + comp_name + '/layer_' + str(layer_idx) + '_scan_path.txt'
    print("Reading scan path file",scan_path_filename_read)

    # Read scan path data
    scan_path_data = get_time_position_power_inp(scan_path_filename_read)
    # start-of-slice power toggle for the first component
    eps = 1e-10
    scan_path_data = add_power_off_entry(scan_path_data, eps)

    # If power is given, update power for scan entries
    if "power" in plan_data:
        updated_power_val = plan_data["power"]
        scan_path_data = update_power(scan_path_data, updated_power_val)

    # Shift time
    scan_path_data = shift_time(scan_path_data, segment_start_time)

    # Number of segments in this file
    num_segments = len(scan_path_data)

    # Segment start times relative to prior segment start time
    relative_segment_start_times = np.zeros(num_segments, dtype=float)
    relative_segment_start_times[0] = 0.0
    for i in range(1, num_segments):
        relative_segment_start_times[i] = scan_path_data[i][0] - scan_path_data[i-1][0]

    # Drip power off entry into file
    dwell_string = str(scan_path_data[0][0]) + ',' + str(scan_path_data[0][1][0]) + ',' + str(scan_path_data[0][1][1]) + ',' + str(scan_path_data[0][1][2]) + ',' + str(scan_path_data[0][2]) + '\n'
    # Drip in first scan path segment
    scan_string = str(scan_path_data[1][0]) + ',' + str(scan_path_data[1][1][0]) + ',' + str(scan_path_data[1][1][1]) + ',' + str(scan_path_data[1][1][2]) + ',' + str(scan_path_data[1][2]) + '\n'

    temp_filename = target_filename + '.tmp'
    if (layer_scan_index == 0):
        fmode = 'w'
    else:
        fmode = 'a'
        shutil.copy(target_filename, temp_filename)

    f = open(temp_filename, fmode)
    f.write(dwell_string)
    f.write(scan_string)
    f.close()
    segment_start_time = segment_start_time + eps
    os.rename(temp_filename, target_filename)

    return num_segments, scan_path_data, segment_start_time, relative_segment_start_times

def drip_feed_segment(segment_num, old_scan_path_data, target_filename, scan_start_time, dwell_time_discretization, relative_segment_start_times):
    '''
    Adds a segment to a copy of the toolpath, which is then moved to the target toolpath
    '''
    temp_filename = target_filename + '.tmp'
    shutil.copy(target_filename, temp_filename)

    # Add scan time up through this point (including added dwells) to segment start time
    segment_end_time = scan_start_time + relative_segment_start_times[segment_num]

    if segment_end_time > 1e-4:
        time_string = str(segment_end_time)
    else:
        time_string = '{:.0e}'.format(segment_end_time)
    segment_end_loc = old_scan_path_data[segment_num][1]
    segment_power = old_scan_path_data[segment_num][2]
    if (segment_power == 0.0):
        currently_dwelling = True
    else:
        currently_dwelling = False
        # If next segment is a dwell, do not pad this segment
        if (old_scan_path_data[segment_num+1][2] == 0.0):
            currently_dwelling = True
    # Append to temporary file
    scan_path_string = time_string + ',' + str(segment_end_loc[0]) + ',' + str(segment_end_loc[1]) + ',' + str(segment_end_loc[2]) + ',' + str(segment_power) + '\n'
    f = open(temp_filename, 'a')
    f.write(scan_path_string)
    #print("Added line ",scan_path_string)
    # If previous segment was not a dwell and the current segment is not a dwell, pad scan with dwell so it ends on an integer multiple of the adamantine output
    if (not currently_dwelling):
        #print("segment end time before padding", segment_end_time)
        remainder = segment_end_time % dwell_time_discretization
        # add only the needed remainder; 0 if already on-grid
        additional = (dwell_time_discretization - remainder) % dwell_time_discretization
        if additional > 0.0:
            segment_end_time = segment_end_time + additional
            if segment_end_time > 1e-4:
                time_string = str(segment_end_time)
            else:
                time_string = '{:.0e}'.format(segment_end_time)
            if (segment_power == 0.0):
                # Use position of current segment for dwell (as this may be the last segment in the file, since the last segment is always a dwell) (?)
                dwell_string = time_string + ',' + str(old_scan_path_data[segment_num][1][0]) + ',' + str(old_scan_path_data[segment_num][1][1]) + ',' + str(old_scan_path_data[segment_num][1][2]) + ',' + str(0.0) + '\n'
            else:
                # Use position of next segment for dwell
                dwell_string = time_string + ',' + str(old_scan_path_data[segment_num+1][1][0]) + ',' + str(old_scan_path_data[segment_num+1][1][1]) + ',' + str(old_scan_path_data[segment_num+1][1][2]) + ',' + str(0.0) + '\n'
            f.write(dwell_string)
            #print("Added dwell pad line ",dwell_string)
            #print("segment end time after padding", segment_end_time, additional)
    f.close()

    # Rename temporary file to target file
    os.rename(temp_filename, target_filename)

    return segment_end_time, currently_dwelling

def drip_feed_dwell_copy(target_filename, old_scan_path_data, dwell_start_time, dwell_time_discretization):
    '''
    Adds a dwell to a copy of the toolpath, which is then moved to the target toolpath
    '''
    temp_filename = target_filename + '.tmp'
    shutil.copy(target_filename, temp_filename)

    # Add scan time up through this point (including added dwells) to dwell start time
    # Dwell increment is a multiple of the print increment
    dwell_end_time = dwell_start_time + dwell_time_discretization
    if dwell_end_time > 1e-4:
        time_string = str(dwell_end_time)
    else:
        time_string = '{:.0e}'.format(dwell_end_time)
    dwell_loc = old_scan_path_data[-1][1]

    # Append to temporary file
    dwell_string = time_string + ',' + str(dwell_loc[0]) + ',' + str(dwell_loc[1]) + ',' + str(dwell_loc[2]) + ',' + str(0.0) + '\n'
    f = open(temp_filename, 'a')
    f.write(dwell_string)
    f.close()

    # Rename temporary file to target file
    os.rename(temp_filename, target_filename)

    return dwell_end_time
    
def drip_feed_dwell(base_filename, toolpath_info, toolpath_module, dwell_time_discretization):
    '''
    Adds a dwell to the toolpath, which is written in its entirety to a new target scan path file
    '''
    get_time_position_power_inp = getattr(toolpath_module, 'get_time_position_power_inp')
    write_event_series = getattr(toolpath_module, 'write_event_series')

    old_scan_path = get_time_position_power_inp(base_filename)

    start_time = old_scan_path[-1][0]

    # Dwell increment is a multiple of the print increment
    end_time = start_time + dwell_time_discretization
    dwell_loc = old_scan_path[-1][1]
    new_dwell = ((end_time, (dwell_loc[0], dwell_loc[1], dwell_loc[2]), 0.0))
    new_toolpath = old_scan_path + [new_dwell]
    temp_filename = base_filename + '.tmp'

    # Don't write SCAN_PATH_END as more dwells may be necessary
    write_event_series(new_toolpath, temp_filename, False)

    os.rename(temp_filename, base_filename)

    return end_time

def drip_feed_end_toolpath(base_filename, toolpath_module):
    '''
    Adds scan path end to the toolpath
    '''
    get_time_position_power_inp = getattr(toolpath_module, 'get_time_position_power_inp')
    write_event_series = getattr(toolpath_module, 'write_event_series')

    old_scan_path = get_time_position_power_inp(base_filename)

    temp_filename = base_filename + '.tmp'
    # Old scan path + SCAN_PATH_END
    write_event_series(old_scan_path, temp_filename, True)
    os.rename(temp_filename, base_filename)
    return

def set_adamantine_deposition_time(filename, dwell_time_discretization):
    f_original = open(filename, 'r')
    f_new = open("temp.info", 'w')

    # Set adamantine deposition time to equal the layer discretization time to ensure accurate calculation of the final time step from deposited layers
    for line in f_original:
        starting_index = line.find("deposition_time")
        if (starting_index > -1):
            ending_index = starting_index + len("deposition_time")
            new_line = line[0:ending_index] + " " + str(dwell_time_discretization) + "\n"
            f_new.write(new_line)
        else:
            f_new.write(line)

    f_original.close()
    f_new.close()
    shutil.move("temp.info", filename)
    
def write_event_series(time_position_power, filename, include_end_message):
    scan_path_string = ""
    
    for entry in time_position_power:
        time = entry[0]
        if time > 1e-4:
            time_string = str(time)
        else:
            time_string = '{:.0e}'.format(time)

        line_string = time_string + ',' + str(entry[1][0]) + ',' + str(entry[1][1]) + ',' + str(entry[1][2]) + ',' + str(entry[2]) + '\n'

        scan_path_string = scan_path_string + line_string

    if include_end_message:
        scan_path_string = scan_path_string + "SCAN_PATH_END"
    else:
        scan_path_string = scan_path_string + "\n"
    
    f = open(filename, 'w')
    f.write(scan_path_string)
    f.close()

def split_into_layers(time_position_power, epsilon=0.001):
    # Assumes that the only z changes are between layers, although robot motions could have z moves too
    split_layers = []
    z_layer_start = time_position_power[0][1][2]

    layer = []
    for entry in time_position_power:
        z = entry[1][2]
        if z > z_layer_start+epsilon:
            z_layer_start = z
            split_layers.append(layer)
            layer = [entry]
        else:
            layer.append(entry)

    return split_layers

def update_power(time_position_power, power):
    # if max power is 0, assume we need to update all power
    if np.max([tpp[2] for tpp in time_position_power]) == 0:
        out = [(tpp[0], tpp[1], power) for tpp in time_position_power]
        return out

    # else update only non-zero powers
    out = []
    for entry in time_position_power:
        if entry[2] > 0:
            new_entry = (entry[0], entry[1], power)
            out.append(new_entry)
        else:
            out.append(entry)
    return out

def shift_time(layer, new_start_time, eps=1e-12):
    current_start_time = layer[0][0]
    shift = new_start_time - current_start_time + 1e-12
    new_layer = []
    for entry in layer:
        new_layer.append((entry[0]+shift, entry[1], entry[2]))
    return new_layer

def add_power_off_entry(layer, eps=1e-10):
    # New layer initial dwell - slightly perturb x position to ensure unique value
    new_layer_position = (layer[0][1][0] + eps, layer[0][1][1], layer[0][1][2])
    new_layer = [(layer[0][0], new_layer_position, 0.0)]
    new_layer = new_layer + [(layer[0][0]+ eps, layer[0][1], layer[0][2])]
    new_layer = new_layer + layer[1:]

    return new_layer

def flatten_layer_z(layer, reference_index=0):
    z_value = layer[reference_index][1][2]

    new_layer = []
    for entry in layer:
        new_layer.append((entry[0], (entry[1][0], entry[1][1], z_value), entry[2]))

    return new_layer

def strip_duplicate_locations(time_position_power, eps=1e-10):
    out = []
    for entry in reversed(time_position_power):
        if len(out) > 0:
            # Check whether this entry is at the same position as the one after it
            if entry[1] != out[-1][1]:
                # This entry is at a new position - append this event
                out.append(entry)
            else:
                # This entry is at the same position as the one after it
                # If both entries are dwells, they can be combined
                # Otherwise, combining entries will change the scan path itself
                # To avoid this, slightly perturb the position rather than combining entries
                # This is necessary as adamantine cannot handle consecutive entries with the same position
                if entry[2] != 0 or out[-1][2] != 0:
                    # Slightly perturb x position
                    modified_position = (entry[1][0] + eps, entry[1][1], entry[1][2])
                    modified_entry = (entry[0], modified_position, entry[2])
                    # Append this modified event
                    out.append(modified_entry)
        else:
            out.append(entry)

    return list(reversed(out))

# ------------------------------------------------------------
# FUNCTIONS TO EXTRACT PATH DATA
# ------------------------------------------------------------
def get_sorted_layer_files(directory, file_pattern):
    pattern = re.compile(r'layer_(\d+)_scan_path\.txt')
    matched_files = []

    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            layer_number = int(match.group(1))
            matched_files.append((layer_number, filename))

    # Sort by the extracted number
    matched_files.sort(key=lambda x: x[0])

    # Return only the filenames in order
    return [filename for _, filename in matched_files]

def get_chunked_value(vals, location, chunk_locations):
    val = vals[0]
    for i in range(0,len(chunk_locations)-1):
        if location > chunk_locations[i] and location < chunk_locations[i+1]:
            val = vals[i]
    return val

def pick_chunk(values: List[float], layer_idx: int, num_layers: int) -> float:
    if isinstance(values, (int, float)):
        return float(values)  # single value for all layers
    n_chunks = len(values)
    width = num_layers / n_chunks
    cidx = min(int(layer_idx // width), n_chunks - 1)
    return values[cidx]
  
def get_toolpath_info(simulation_template_dir, sorted_data, layer_end_time_discretization=None, layer_limit=None):
    toolpath_info = {}
    toolpath_info['scan_path_out'] = "scan_path.inp"
    toolpath_info['lump_size'] = 2
    toolpath_info['includes_end_message'] = True
    toolpath_info['layer_end_time_discretization'] = layer_end_time_discretization
    toolpath_info['num_layers'] = len(sorted_data)
    
    # Use print layers to get the layer height - don't store scan path data yet
    temp_base_layers_print = []
    for fn in range(0, layer_limit):
        # Assume that "print_layers" is given for each layer, otherwise, take the first component name given to search for a scan path
        if ("print_layers" in sorted_data[0]["components"]):
            comp_name = "print_layers"
        else:
            comp_name = sorted_data[0]["components"][0]
        filename = simulation_template_dir + comp_name + '/layer_' + str(fn) + '_scan_path.txt'
        temp_base_layers_print.append(get_time_position_power_inp(filename))

    # Currently we assume that all layers are the same height
    toolpath_info['layer_height'] = (temp_base_layers_print[toolpath_info['lump_size']][1][1][2] - temp_base_layers_print[0][1][1][2])/toolpath_info['lump_size']

    # By default, select all layers
    toolpath_info['selected_layers'] = (0, toolpath_info['num_layers'])
    return toolpath_info

# FIXME: Update to new print plan file format without hardcoded dwell_0, reheat_power, etc
def generate_print_plan_file(toolpath_info, plan_filename, actual_layer_variables=None):
    dwell_0       = toolpath_info['dwell_0']
    reheat_power  = toolpath_info['reheat_power']
    dwell_1       = toolpath_info['dwell_1']
    num_layers = toolpath_info.get('num_layers')

    data = []
    for layer_idx in range(num_layers):
        d0 = pick_chunk(dwell_0,      layer_idx+2, num_layers)
        rp = pick_chunk(reheat_power, layer_idx+2, num_layers)
        d1 = pick_chunk(dwell_1,      layer_idx+2, num_layers)

        if actual_layer_variables is not None:
            d1 = actual_layer_variables[layer_idx]['dwell_1']
            d0 = actual_layer_variables[layer_idx]['dwell_0']

        entry = {
            "layer": layer_idx,
            "power": rp,
            "dwell_0": d0,
            "dwell_1": d1
        }
        data.append(entry)

    with open(plan_filename, "w") as f:
        json.dump(data, f, indent=2)

# Called by create_toolpath_layers to create a toolpath from an existing file (a toolpath may consist of multiple files, in which case this function will be called multiple times)
def create_toolpath_single_scan(filename, plan_data, layer_idx, num_layers, section_start_time, first_comp, last_comp, new_tpp, override_dwells, set_dwell_every_n_layers, layer_end_time_discretization, discretize_scan=False):

    # Scan path from this file
    temp_store_layer_section = get_time_position_power_inp(filename)
        
    # If power is given, update power for scan entries
    if "power" in plan_data:
        updated_power_val = plan_data["power"]
        temp_store_layer_section = update_power(temp_store_layer_section, updated_power_val)

    # Shift time
    temp_store_layer_section = shift_time(temp_store_layer_section, section_start_time)
    
    # start-of-slice power toggle for the first component
    if (first_comp):
        temp_store_layer_section = add_power_off_entry(temp_store_layer_section)
    
    # ensure the very last print-point has power=0.0 # XXX is this done twice with add_power_off_entry?
    if temp_store_layer_section:
        t_last, pos_last, _ = temp_store_layer_section[-1]
        temp_store_layer_section[-1] = (t_last, pos_last, 0.0)
    new_tpp += temp_store_layer_section

    # If dwell is given, add the dwell (unless overriding). Dwells assumed to be zero for a scan if not given in print plan
    dwell_time = 0.0
    if ("dwell" in plan_data) and (not override_dwells):
        dwell_time = plan_data["dwell"]
    dwell_offset_time = 0.0
    if ("dwell_offset" in plan_data) and (not override_dwells):
        dwell_offset_time = plan_data["dwell_offset"]
    d0 = pick_chunk(dwell_time, layer_idx+2, num_layers) + dwell_offset_time

    # If we happen to be at the start of one of the sections and at the end of the current layer, add a 20 min dwell instead of the optimized dwell.
    if last_comp and set_dwell_every_n_layers and ((layer_idx + 2) in quarter_layers_1based):
        print(f'adjusting dwell to 20min for layer_idx {layer_idx}')
        d0 = 1200.0 + dwell_section_offset + dwell_offset_time # measured offset...

    additional = 0.0
    if last_comp or discretize_scan:
        layer_end_time = new_tpp[-1][0]
        pos = new_tpp[-1][1]
        if d0 > 0.0:
            new_tpp += time_position_power_dwell(layer_end_time, pos, d0)
            layer_end_time = new_tpp[-1][0]
        # Snap to discretization (optional)
        print("layer end time before padding", layer_end_time)
        if layer_end_time_discretization is not None:
            remainder = layer_end_time % layer_end_time_discretization
            # add only the needed remainder; 0 if already on-grid
            additional = (layer_end_time_discretization - remainder) % layer_end_time_discretization
            if additional > 0.0:
                new_tpp += time_position_power_dwell(layer_end_time, pos, additional)
                layer_end_time = new_tpp[-1][0]
            print("layer end time after padding", layer_end_time, additional)
        section_end_time = layer_end_time
    else:
        if d0 > 0.0:
            new_tpp += time_position_power_dwell(new_tpp[-1][0], new_tpp[-1][1], d0)
        section_end_time = new_tpp[-1][0]
    return new_tpp, section_end_time, additional

def create_toolpath_layers(sorted_data, toolpath_info, simulation_template_dir, override_dwells = False):
    """
    Builds a new scan_path.inp by:
      - loading each peeled layer from specified directories for each layer (via multiple calls to create_toolpath_single_scan)
      - inserting dwells as necessary, and using any specified powers to override power from files
      - shifting times accordingly to create a scan path from multiple files of data, each read via a separate call to create_toolpath_single_scan
    """
    # Unpack inputs
    dwell_section_offset          = toolpath_info.get('dwell_section_offset',0.0)
    includes_end_message          = toolpath_info.get('includes_end_message', True)
    layer_end_time_discretization = toolpath_info.get('layer_end_time_discretization', None)
    set_dwell_every_n_layers      = toolpath_info.get('set_dwell_every_n_layers')
    base_split_layers_print       = toolpath_info.get('base_split_layers_print')
    base_split_layers_reheat      = toolpath_info.get('base_split_layers_reheat')
    num_layers                    = toolpath_info.get('num_layers')
    selected_layers               = toolpath_info.get('selected_layers')  # tuple (start, end)

    # Infer the number of layers and selected layers if not set in parameters
    if num_layers is None:
        num_layers = len(base_split_layers_print)
        toolpath_info['num_layers'] = num_layers

    if selected_layers is None:
        selected_layers = (0, num_layers)
        toolpath_info['selected_layers'] = selected_layers

    # Build the new time–position–power series
    new_tpp = []
    section_start_time = 1e-10
    # quarters: 1/4, 2/4, 3/4 of the layer count
    q = num_layers // 4  # integer floor, stable
    quarter_layers_1based = {q*1, q*2, q*3}
    actual_layer_variables = []

    for layer_idx in range(*toolpath_info['selected_layers']):
        print(f'Writing toolpath for layer {layer_idx}')
        temp_store_layer_section = []
        comp_num = 0
        num_comps = len(sorted_data[layer_idx]["components"])
        for comp_name in sorted_data[layer_idx]["components"]:
            
            filename = simulation_template_dir + comp_name + '/layer_' + str(layer_idx) + '_scan_path.txt'
            plan_data = sorted_data[layer_idx][comp_name]
            first_comp = (comp_num == 0)
            last_comp = (comp_num == num_comps-1)

            # Scan path from this file
            new_tpp, section_start_time, additional = create_toolpath_single_scan(filename, plan_data, layer_idx, num_layers, section_start_time, first_comp, last_comp, new_tpp, override_dwells, set_dwell_every_n_layers, layer_end_time_discretization)
            
            # Last scan path file for this layer
            if (comp_num == num_comps-1):
                # Actual layer variables: variable names hardcoded for compatibility with opt-driver
                first_comp_name = sorted_data[layer_idx]["components"][0]
                if ("dwell" in sorted_data[layer_idx][first_comp_name]):
                    d0_actual_layer_var = sorted_data[layer_idx][first_comp_name]["dwell"]
                else:
                    d0_actual_layer_var = 0.0
                if (len(sorted_data[layer_idx]["components"]) >= 2):
                    second_comp_name = sorted_data[layer_idx]["components"][1]
                    if ("dwell" in sorted_data[layer_idx][second_comp_name]):
                        d1_actual_layer_var = sorted_data[layer_idx][second_comp_name]["dwell"]
                    else:
                        d1_actual_layer_var = 0.0
                    if ("power" in sorted_data[layer_idx][second_comp_name]):
                        reheat_power_actual_layer_var = sorted_data[layer_idx][second_comp_name]["power"]
                    else:
                        reheat_power_actual_layer_var = 0.0
                else:
                    d1_actual_layer_var = 0.0
                    reheat_power_actual_layer_var = 0.0
                actual_layer_variables.append({"dwell_0": [d0_actual_layer_var], "dwell_1": [d1_actual_layer_var+additional], "reheat_power": reheat_power_actual_layer_var})
                
            # Read next component
            comp_num = comp_num + 1

    # Clean up and return
    tpp_clean = strip_duplicate_locations(new_tpp)
    return tpp_clean, actual_layer_variables

# FIXME: Update based on new create_toolpath function
def write_toolpath(toolpath_info):
    
    tpp_clean, layer_end_times, layer_padding = create_toolpath_layers(toolpath_info)
    write_event_series(tpp_clean, toolpath_info['scan_path_out'], toolpath_info['includes_end_message'])

# FIXME: Remove hardcoded dwell_0, dwell_1, reheat_power variables
def generate_control_options(layer_time_discretization, num_forward_sims, sampling_strategy, toolpath_info, planned_controls):
    
    def get_uniform_random_discrete_value(bounds, discrete_step):
        x = random.uniform(bounds[0], bounds[1])
        out = round(x / discrete_step) * discrete_step
        return out


    dwell_0_bounds = (layer_time_discretization, 12*layer_time_discretization)
    dwell_1_bounds = (layer_time_discretization, 12*layer_time_discretization)
    reheat_power_bounds = (0.0, 1350.0)

    # The first entry should be the planned control
    modified_toolpath_info_list = []
    modified_toolpath_info_list.append(copy.deepcopy(toolpath_info))
    modified_toolpath_info_list[-1]['dwell_0'] = [planned_controls['dwell_0']]
    modified_toolpath_info_list[-1]['dwell_1'] = [planned_controls['dwell_1']]
    modified_toolpath_info_list[-1]['reheat_power'] = [planned_controls['power']]

    num_variables = 3
    num_adjacent_options = 2 * num_variables + 1


    match sampling_strategy:
        case 'uniform_random':
            for i in range(num_forward_sims-1):   
                
                unique_parameter_set_found = False
                counter = 0

                while not unique_parameter_set_found:
                    # Draw new values
                    dwell_0_val = get_uniform_random_discrete_value(dwell_0_bounds, layer_time_discretization)
                    dwell_1_val = get_uniform_random_discrete_value(dwell_1_bounds, layer_time_discretization)
                    reheat_val = get_uniform_random_discrete_value(reheat_power_bounds, 1.0)

                    # Compare to existing values
                    for entry in modified_toolpath_info_list:
                        if np.isclose(dwell_0_val, entry['dwell_0']) and np.isclose(dwell_1_val, entry['dwell_1']) and np.isclose(reheat_val, entry['reheat_power']):
                            counter = counter + 1
                            if counter > 1000:
                                print('Error: Unable to find enough unique random samples for the control simulations')
                                sys.exit()
                            continue

                    print("Control options:", dwell_0_val, dwell_1_val, reheat_val)

                    modified_toolpath_info_list.append(copy.deepcopy(toolpath_info))
                    modified_toolpath_info_list[-1]['dwell_0'] = [dwell_0_val]
                    modified_toolpath_info_list[-1]['dwell_1'] = [dwell_1_val]
                    modified_toolpath_info_list[-1]['reheat_power'] = [reheat_val]
                    unique_parameter_set_found = True

            return modified_toolpath_info_list
        
        case 'normal_random':
            # TODO
            return modified_toolpath_info_list

        case 'individual_perturbations':
            neighbors_to_sample_from = np.ceil(num_forward_sims/num_adjacent_options)
            # TODO
            return modified_toolpath_info_list
        
        case _:
            print('Error: Invalid sampling strategy chosen:', sampling_strategy)
            sys.exit()
    
    return modified_toolpath_info_list
    
# FIXME: main function is hardcoded and should be modified or removed
if __name__ == "__main__":
    cwd = os.getcwd()
    toolpath_info = {
        'print_path'    : os.path.join(cwd, "print_layers"),
        'reheat_path'   : os.path.join(cwd, "reheat_layers"),
        
        # U6
        'dwell_0'       : [5, 5, 5, 5],
        'reheat_power'  : [900, 900, 900, 900],
        'dwell_1'       : [5],

        # Preopt / O1
    	# 'dwell_0'       : [111.95093, 46.80592, 98.10942, 48.97092],
        # 'reheat_power'  : [990.72487, 687.35519, 531.07650, 1340.78436],
        # 'dwell_1'       : [55.20886],

        # experimentally tuned to minimize lags
        'dwell_0_offset'       : 0.586,
        'dwell_1_offset'       : 0.7375,
        'dwell_section_offset' : -1.94,

        'scan_path_out' : "scan_path_test.inp",
        'includes_end_message': True,
        'set_dwell_every_n_layers': True
    }    
    print("Running write_toolpath test with:")
    print(f"  dwell_0       = {toolpath_info['dwell_0']}")
    print(f"  reheat_power  = {toolpath_info['reheat_power']}")
    print(f"  dwell_1       = {toolpath_info['dwell_1']}")
    print(f"  writing -> {toolpath_info['scan_path_out']}")    
    # call the core function    
    write_toolpath(toolpath_info)    
    if os.path.exists(toolpath_info['scan_path_out']):
        print(":) scan_path_test.inp successfully written")
    else:
        print(":( failed to write scan_path_test.inp")
    
    
    '''
    preopt
    	'dwell_0'       : [111.95093, 46.80592, 98.10942, 48.97092],
        'reheat_power'  : [990.72487, 687.35519, 531.07650, 1340.78436],
        'dwell_1'       : [55.20886],
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
