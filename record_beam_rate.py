#!/usr/bin/env python3
import sys
import logging
import datetime
import os
import re

def get_file_index(runNo):
  log_files = [f for f in os.listdir("log") if re.match(rf"beam_rate_RUN\({runNo:03}\)_\d+\.log",f)]
  if not log_files:
    return 0

  indices = []
  for f in log_files:
    match=re.search(r"beam_rate_RUN\({runNo:03}\)_(\d+)\.log",f)
    if match:
      indices.append(int(match.group(1)))
  return max(indices) + 1 if indices else 0

def setup_logger(runNo):
  index = get_file_index(runNo);
  logname = f"log/beam_rate_RUN{runNo:03}_{index}.log"
  print()
  print(f"Starting {logname}")
  logging.basicConfig(filename=logname,
                      level=logging.INFO,
                      format='%(asctime)s   %(message)s',
                      datefmt='%y%m%d   %H:%M:%S')
  logger = logging.getLogger()
  return logger, logname

def main(runNo):
  logger, logname = setup_logger(runNo)
  print()
  print("Enter beam rate values. Type 'q' to stop")
  
  while True:
    beam_rate_unit = input("== Beam type (pc/si): ")
    if beam_rate_unit.lower() in ["end","q",".q","x","exit"]:
      print()
      print(f"End beam rate logger {logname}")
      break
    elif beam_rate_unit.lower() in ["pc","si"]:
      beam_rate = input("== Beam rate: ")
      if beam_rate.lower() in ["end","q",".q","x","exit"]:
        print()
        print(f"End beam rate logger {logname}")
        break
      else:
        logger.info(f"{beam_rate_unit}  {beam_rate}")
        #try:
        #  beam_rate = float(beam_rate)
        #  logger.info(f"{beam_rate}")
        #  #print(datetime.datetime.now().strftime(f"++ %y%m%d   %H:%M:%S   {user_input}"))
        #except ValueError:
        #  print("!  Invalid input")
    else:
      print("Invalid input")
      continue

if __name__ == "__main__":
  if len(sys.argv)<2:
    print("ex) ./record_beam_rate 11")
  else:
    main(sys.argv[1])
