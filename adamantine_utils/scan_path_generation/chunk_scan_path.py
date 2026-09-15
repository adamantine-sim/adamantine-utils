import os
import sys
import json
from pathlib import Path

# Input argument 1: a directory containing driver_config.json, which in turn contains "master_scan_path_filename" which is a multilayer scan path
# Input argument 2: either "segment" or "event", governing the format of the file specified by master_scan_path_filename
# Outputs: An updated driver_config.json file with the correct number of layers based on the parsed scan path; a new directory "print_layers" within the directory specified by input argument 1, containing a scan path file for each layer of the master scan path (parsed by Z coordinate); a new file "print_plan.json" within the directory specified by input argument 1
if __name__ == "__main__":

    # Directory containing driver_config.json, mesh, input file, and multilayer scan path file
    template_directory = str(Path(sys.argv[1]))
    # Segment or event based scan path format
    scan_format = str(sys.argv[2])
    
    # Find and parse config file
    config_filename = template_directory + "/driver_config.json"
    if os.path.exists(config_filename):
        with open(config_filename) as json_data:
            config = json.load(json_data)
            json_data.close()
    else:
        raise ValueError("Error: No driver_config.json file found in the template directory")

    # New directory inside of template directory for individual layer scan paths
    print_layers_dir = template_directory + "/" + "print_layers"
    Path(print_layers_dir).mkdir(parents=True, exist_ok=True)

    # Multilayer scan path filename from config file
    path_filename = template_directory + "/" + config["master_scan_path_filename"]
    with open(path_filename, "r") as path_datastream:
        lines = path_datastream.readlines()
    # Z coordinate of first line - segment-based format has 3 lines of header, event-based format does not
    if (scan_format == "event"):
        first_line_num = 0
    else:
        first_line_num = 3
    line_parts = lines[first_line_num].split()
    last_scan_z = float(line_parts[3])
    # Counter variables
    layer_num = 0
    # Write individual layer scan path files
    for line_num in range(first_line_num, len(lines)-1):
        if (scan_format == "segment"):
            line_parts = lines[line_num].split()
        else:
            line_parts = lines[line_num].split(",")
        current_scan_z = float(line_parts[3])
        if (current_scan_z != last_scan_z):
            lines_this_layer = lines[first_line_num:line_num]
            num_lines_this_layer = line_num - first_line_num
            output_filename = print_layers_dir + "/layer_" + str(layer_num) + "_scan_path.txt"
            with open(output_filename, "w") as file:
                if (scan_format == "segment"):
                    file.write(lines[0])
                    num_lines_this_layer_str = str(num_lines_this_layer) + "\n"
                    file.write(num_lines_this_layer_str)
                    file.write(lines[2])
                file.writelines(lines_this_layer)
                layer_num = layer_num + 1
                first_line_num = line_num
                last_scan_z = current_scan_z

    # Update number of layers in config file and write new config file
    config["layer_limit"] = layer_num
    with open(config_filename) as file:
        config_data = file.readlines()
    num_config_lines = len(config_data)
    for line_num in range(num_config_lines):
        if "layer_limit" in config_data[line_num]:
            config_data[line_num] = "    \"layer_limit\": " + str(layer_num) + ",\n"
    with open(config_filename, "w") as file:
        file.writelines(config_data)
    
    # Write print plan file
    print_plan_filename = template_directory + "/print_plan.json"
    with open(print_plan_filename, "w") as file:
        file.write("[\n")
        for print_layer in range(layer_num):
            file.write("    {\n")
            layer_string_write = "        \"layer\": " + str(print_layer) + ",\n"
            file.write(layer_string_write)
            file.write("        \"components\": [\"print_layers\"]\n")
            file.write("    }")
            if (print_layer < layer_num-1):
                file.write(",")
            file.write("\n")
        file.write("]\n")

