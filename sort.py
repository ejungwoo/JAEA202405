import csv
import sys
import os

def sort_ts_and_save(path_to_input,path_to_output,file_name):
    name_in = f"{path_to_input}/{file_name}"
    iRUN = name_in.find("RUN")
    #name_out = f"{path_to_output}/SORTED_{name_in[iRUN:iRUN+6]}_{name_in[-7:-4]}.dat"
    name_out = f"{path_to_output}/{name_in[iRUN:iRUN+6]}_{name_in[-7:-4]}.dat"
    print(f"  {name_in} -> {name_out}")
    with open(name_in,'r') as f1:
        reader = list(csv.reader(f1))

    data2 = []
    for line in reader:
        if len(line)==5:
            data2.append([line[0], line[1], line[2], int(line[3]), line[4]])

    data_sorted = sorted(data2, key=lambda row: row[3])

    with open(name_out,'w') as f2:
        for line in data_sorted:
            print(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]}",file=f2)

    print("finished!")

    import os

def get_files_by_run_number(path_to_input, run_number):
    matching_files = []
    prefix = f"RUN{run_number:03}_"

    # List all files in the directory
    for file_name in os.listdir(path_to_input):
        # Check if the file name starts with the specified run number prefix
        if file_name.startswith(prefix) and file_name.endswith(".dat"):
            matching_files.append(file_name)

    if matching_files:
        print(f"RUN {run_number}:")
        for file in matching_files:
            print("  "+file)
    else:
        print(f"No files found for RUN{run_number}.")

    return matching_files

# Example usage
if __name__ == "__main__":
    path_to_input = "data_input"
    path_to_output = "data_sorted"
    run_list = [10]
    #run_list = [88]
    #run_list = list(range(0,28))
    if len(sys.argv)>=2:
        run_list.clear()
        run_list.append(int(sys.argv[1]))
    for run_number in run_list:
        files = get_files_by_run_number(path_to_input, run_number)
        for file in files:
            sort_ts_and_save(path_to_input,path_to_output,file)
