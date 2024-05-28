import csv
import sys

def sort_ts_and_save(nameIn):
    #files = [f for f in os.listdir("log") if re.match(rf"RUN{runNo:03}.beam_rate.\d+\.log",f)]
    #exit()
    iRUN = nameIn.find("RUN")
    newRUN = int(nameIn[iRUN+3])+5
    nameOut = nameIn[:iRUN+3] + f"{newRUN}" + nameIn[iRUN+4:] 
    print(f"{nameIn} -> {nameOut}")

    with open(nameIn,'r') as f1:
        reader = list(csv.reader(f1))

    data2 = []
    for line in reader:
        if len(line)==5:
            data2.append([line[0], line[1], line[2], int(line[3]), line[4]])

    data_sorted = sorted(data2, key=lambda row: row[3])

    with open(nameOut,'w') as f2:
        for line in data_sorted:
            print(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]}",file=f2)

    print("finished!")

if __name__ == "__main__":
    sort_ts_and_save("../RUN012_202405191441_list_000.dat")
    #sort_ts_and_save(sys.argv[1])
    #sort_ts_and_save("../RUN088_202405212142_list_000.dat")
    #sort_ts_and_save("../RUN000_202405161147_list_000.dat")
    #sort_ts_and_save("../RUN001_202405161333_list_000.dat")
    #sort_ts_and_save("../RUN002_202405161147_list_000.dat")
    #sort_ts_and_save("../RUN003_202405161410_list_000.dat")
    #sort_ts_and_save("../RUN004_202405161604_list_000.dat")
    #sort_ts_and_save("../RUN005_202405161608_list_000.dat")
    #sort_ts_and_save("../RUN006_202405161611_list_000.dat")
    #sort_ts_and_save("../RUN007_202405161612_list_000.dat")
    #sort_ts_and_save("../RUN008_202405161637_list_000.dat")
    #sort_ts_and_save("../RUN009_202405161650_list_000.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_000.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_001.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_002.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_003.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_004.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_100.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_200.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_300.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_400.dat")
    #sort_ts_and_save("../RUN010_202405161747_list_408.dat")
    #sort_ts_and_save("../RUN011_202405170946_list_000.dat")
    #sort_ts_and_save("../RUN012_202405191441_list_000.dat")
    #sort_ts_and_save("../RUN013_202405191511_list_000.dat")
    #sort_ts_and_save("../RUN014_202405191606_list_000.dat")
    #sort_ts_and_save("../RUN015_202405200959_list_000.dat")
    #sort_ts_and_save("../RUN016_202405211917_list_000.dat")
    #sort_ts_and_save("../RUN017_202405211926_list_000.dat")
    #sort_ts_and_save("../RUN018_202405211938_list_000.dat")
    #sort_ts_and_save("../RUN019_202405212142_list_000.dat")
    #sort_ts_and_save("../RUN020_202405212343_list_000.dat")
    #sort_ts_and_save("../RUN021_202405212353_list_000.dat")
    #sort_ts_and_save("../RUN022_202405220052_list_000.dat")
    #sort_ts_and_save("../RUN023_202405220128_list_000.dat")
    #sort_ts_and_save("../RUN024_202405220158_list_000.dat")
    #sort_ts_and_save("../RUN025_202405220217_list_000.dat")
    #sort_ts_and_save("../RUN026_202405220235_list_000.dat")
    #sort_ts_and_save("../RUN027_202405220257_list_000.dat")
    #sort_ts_and_save("../RUN028_202405220335_list_000.dat")
