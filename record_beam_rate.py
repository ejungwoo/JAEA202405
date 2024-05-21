#!/usr/bin/env python3
import sys
import logging
import datetime
import os
import re

def get_file_creation_time(filepath):
  pass
  #run_files = [f for f in os.listdir(filepath) if re.match(rf"RUN{runNo:03}_\w{11}_\list_000.log",f)]
  #print(run_files)
  #try:
  #  if os.name=='nt':
  #    creation_time = os.path.getctime(filepath)
  #  else:
  #    stat = os.stat(filepath)
  #    creation_time = stat.st_birthtime if hasattr(stat, 'st_birthtimella') else stat.st_mtime
  #  return datetime.fromtimestamp(creation_time)
  #except Exeption as e:
  #  print(f"{e}")
  #  return None

def get_file_index(runNo):
  log_files = [f for f in os.listdir("log") if re.match(rf"RUN{runNo:03}.beam_rate.\d+\.log",f)]
  if not log_files:
    return 0

  indices = []
  for f in log_files:
    match=re.search(rf"RUN{runNo:03}.beam_rate.(\d+)\.log",f)
    if match:
      indices.append(int(match.group(1)))
  return max(indices) + 1 if indices else 0

def setup_logger(runNo):
  index = get_file_index(runNo)
  #ctime = get_file_creation_time(runNo)
  logname = f"log/RUN{runNo:03}.beam_rate.{index}.log"
  print()
  print(f"Starting {logname}")
  logging.basicConfig(filename=logname,
                      level=logging.INFO,
                      format='%(asctime)s  %(message)s',
                      datefmt='%y%m%d   %H:%M:%S')
  logger = logging.getLogger()
  return logger, logname

def main(runNo):
  logger, logname = setup_logger(runNo)
  logger.info("xx  0")
  print()
  print("Enter beam rate values. Type 'q' to stop")
  
  while True:
    beam_rate_unit = input("== Beam type (pc/si): ")
    if beam_rate_unit.lower() in ["end","q",".q","x","exit"]:
      print(f"\nEnd beam rate logger {logname}")
      break
    elif beam_rate_unit.lower() in ["pc","si"]:
      beam_rate = input("== Beam rate: ")
      if beam_rate.lower() in ["end","q",".q","x","exit"]:
        print(f"\nEnd beam rate logger {logname}")
        break
      else:
        logger.info(f"{beam_rate_unit}  {beam_rate}")
    else:
      print("Invalid input")
      continue

if __name__ == "__main__":
  print("=================================================================")
  if len(sys.argv)<2:
    print("Usage:")
    print("./record_beam_rate [RunNo]")
  else:
    main(int(sys.argv[1]))
  print("=================================================================")
