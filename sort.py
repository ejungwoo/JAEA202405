import csv
import sys

def soort(nameIn):
  nameOut = nameIn[:6] + "9" + nameIn[7:] 

  with open(nameIn,'r') as f1:

    reader = csv.reader(f1)
    data = list(reader)

  data2 = []
  for line in data:
    if len(line)==5:
      data2.append([line[0], line[1], line[2], int(line[3]), line[4]])

  datasorted = sorted(data2, key=lambda row: row[3])

  with open(nameOut,'w') as f2:
    for line in datasorted:
      print(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]}",file=f2)

  print(nameOut)

if __name__ == "__main__":
  soort(sys.argv[1])
