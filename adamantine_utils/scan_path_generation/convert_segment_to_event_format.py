import numpy as np
import sys
from pathlib import Path

# Input argument 1 (required): scan path filename in the segment-based format
# Input argument 2 (optional): location of new scan path file in the event-based format
# If the optional input argument is given, the segment-based and event-based scan paths will have the same name, but will be located in different file system locations
# If the optional input argument is not given, the event-based scan path will be in the same directory as the segment-based scan path, but with "_events" appended to the end of the filename
if __name__ == "__main__":

    # Directory and filename for scan path in the event-based format
    input_dir_filename = str(Path(sys.argv[1]))
    parts = input_dir_filename.rsplit('/', 1)
    input_scan_path_directory = parts[0]
    input_scan_path_filename = parts[1]
    # Optional output scan path directory from command line (otherwise place new scan path in input scan path directory with new name)
    if (len(sys.argv) > 2):
        output_scan_path_directory = str(Path(sys.argv[2]))
        output_scan_path_filename = input_scan_path_filename
    else:
        output_scan_path_directory = input_scan_path_directory
        parts = input_scan_path_filename.rsplit('.', 1)
        output_scan_path_filename = parts[0] + "_events." + parts[1]
    output_dir_filename = output_scan_path_directory + "/" + output_scan_path_filename

    # Parsed scan path data
    with open(input_dir_filename, "r") as old_scan_path_datastream:
        old_scan_path_data = old_scan_path_datastream.readlines()
    total_path_time = 0.0
    # Second line gives the total number of path segments - substract 1 for SCAN_PATH_END
    tot_num_scans = int(old_scan_path_data[1]) - 1
    new_scan_path_data = [None] * tot_num_scans
    counter = 0
    # Line by line parsing of old scan path data, filling new_scan_path_data with event-based format data
    for line_num in range(3, len(old_scan_path_data)-1):
        line_parts = old_scan_path_data[line_num].split()
        segment_type = int(line_parts[0])
        segment_x = line_parts[1]
        segment_y = line_parts[2]
        segment_z = line_parts[3]
        nominal_power = line_parts[4]
        if (segment_type == 1):
            # Dwell
            dwell_time = float(line_parts[5])
            segment_endtime = total_path_time + dwell_time
        else:
            # Scan
            scan_speed = float(line_parts[5])
            last_line_parts = old_scan_path_data[line_num-1].split()
            last_segment_x_float = float(last_line_parts[1])
            last_segment_y_float = float(last_line_parts[2])
            last_segment_z_float = float(last_line_parts[3])
            segment_x_float = float(segment_x)
            segment_y_float = float(segment_y)
            segment_z_float = float(segment_z)
            scan_distance = np.sqrt((segment_x_float-last_segment_x_float)**2 + (segment_y_float-last_segment_y_float)**2 + (segment_z_float-last_segment_z_float)**2)
            scan_time = scan_distance / scan_speed
            segment_endtime = total_path_time + scan_time
        event_data = str(segment_endtime) + ' ' + segment_x + ' ' + segment_y + ' ' + segment_z + ' ' + nominal_power + '\n'
        new_scan_path_data[counter] = event_data
        counter = counter+1
        total_path_time = segment_endtime

    # Write new events-based scan path to file
    with open(output_dir_filename, "w") as file:
        file.writelines(new_scan_path_data)
        file.write("SCAN_PATH_END\n")
