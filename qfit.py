#!/usr/bin/env python3
from shutil import copyfile
import os
import argparse

def getargs():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  # plan to only go through one target at a time
  parser.add_argument('infile', type=str)
  parser.add_argument('-o', '--outfile', type=str)
  parser.add_argument('--pdf', action='store_true')
  return parser.parse_args()

def main():
  # The badzombie command runs over a list of files
  args = getargs()
  cmd = 'qfit'
  if args.pdf:
    cmd = 'qfit --pdf'

  target = args.infile
  if args.outfile is not None:
    copyfile(args.infile, args.outfile)
    target = args.outfile

  full_cmd = f"{cmd} \"{target}\""
  print(f"{full_cmd}")
  os.system(full_cmd)

if __name__ == '__main__':
  main()
