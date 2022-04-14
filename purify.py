#!/usr/bin/env python3
import os
import os.path as path
import argparse

def getargs():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  # plan to only go through one target at a time
  parser.add_argument('folder', type=str)
  return parser.parse_args()

def main():
  # The badzombie command runs over a list of files
  cmd = 'badzombie'
  args = getargs()
  truenames = []
  for dirpath, dirnames, filenames in os.walk(args.folder):
    for fname in filenames:
      if fname.endswith('.root'):
        truenames.append( path.join(dirpath, fname) )
  # Chunk size
  chnk = 100
  truenames = [truenames[i * chnk:(i+1) * chnk] for i in range((len(truenames) + chnk - 1) // chnk)]
  for fchunk in truenames:
    full_cmd = f"{cmd} {' '.join(fchunk)}"
    os.system(full_cmd)

if __name__ == '__main__':
  main()
