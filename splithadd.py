#!/usr/bin/env python3
import os
import argparse
import numpy as np

def getargs():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-f', '--files', type=str, nargs='+', help='List of root files to hadd')
  parser.add_argument('-o', '--output', type=str, help='Output basename', default='output.root')
  parser.add_argument('-s', '--size', type=float, default='50.0',
      help='Size in gigabytes of each output')
  return parser.parse_args()

def main():
  args = getargs()
  flist = np.array(args.files)
  blist = np.array([os.path.getsize(f)/1e9 for f in flist])
  total_size = sum(blist)
  print(f'Running hadd on {len(flist)} files, totalling {total_size} GB')
  splits = int(np.ceil(total_size/args.size))
  print(f'-- Splitting into {splits} files, each <= {args.size}')
  counter = 0
  if splits > 1:
    while len(flist) > 0:
      cdf = np.cumsum(blist)
      mask = (cdf < args.size)
      outlist = flist[mask]
      outfile = (args.output).replace('.root', f'_{counter}.root')
      print(outfile)
      cmd = f'hadd -k {outfile} {" ".join(outlist)}'
      # update / remove
      flist = flist[~mask]
      blist = blist[~mask]
      counter += 1
      os.system(cmd)
  else:
    cmd = f'hadd -k {args.output} {" ".join(flist)}'
    os.system(cmd)

if __name__ == '__main__':
  main()
