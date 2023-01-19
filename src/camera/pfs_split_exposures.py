#!/usr/bin/python3
# Split the HDR exposure files into directories based on the capture time
# work in progress. Currently works only with Sony's ARW files

import os
import sys
import subprocess
import re
import argparse

parser = argparse.ArgumentParser(description='Move images files that belong to the same HDR exposure sequence into separate directories.')
parser.add_argument('dir', metavar='dir', nargs='?', help='the directory with images files', default='.')
parser.add_argument('--threshold', '-t', type=float, dest='threshold', default=5,
                    help='Exposures taken less than "threshold" seconds apart will belong to the same exposure')
parser.add_argument('--simulate', '-s', dest='simulate', action='store_true', 
                    help='Do not perform any action but display how files would be split')
parser.add_argument('--minexp', '-m', dest='minexp', default=3, 
                    help='The minimum number of exposures in each exposure stack. If there are fewer exposures than that number, files will not be moved.')
parser.add_argument('--basename', '-b', dest='basename', default='hdr', 
                    help='The base name for each created directory with an exposure stack. By default the directories will be hdr_000, hdr_001, ...')

args = parser.parse_args()

print("Searching for HDR exposure stacks in the directory '{}'".format(args.dir))
if args.simulate:
    print( '  Simulation mode - no action will be performed' )

def move_exp_to_dir( es, stack_no ):
    dir_name='{0}_{1:03d}'.format(args.basename, stack_no)
    if not args.simulate: 
        os.makedirs(dir_name,exist_ok=True)
    for file in es:
        print( 'Moving {} to {}'.format(file, os.path.join(dir_name,file)) )
        if not args.simulate: 
            os.rename(file, os.path.join(dir_name,file))

prev_capture_time=None
prev_shutter=None
prev_filename=None
exp_stack=list()
stack_no=1
for filename in os.listdir(args.dir):
    if filename[-4:] == '.ARW' or filename[-4:] == '.CR2':
        cp = subprocess.run( ['dcraw', '-v', '-i', '{}'.format(filename)], stdout=subprocess.PIPE, encoding="utf-8" )
        lines = cp.stdout.splitlines()
        capture_time=None
        for ll in lines:
            if ll[0:10] == "Timestamp:":
                ts = re.split('[ ]+', ll)
                th = ts[4].split(':')
                capture_time=float(th[0])*3600 + float(th[1])*60 + float(th[2])
            if ll[0:8] == "Shutter:":
                ts = ll.split(' ')
                th = ts[1].split('/')
                if len(th)>1:
                    shutter = 1/float(th[1])
                else:
                    shutter = float(ts[1])
                
        if capture_time!=None and prev_capture_time!=None:
            capt_time_diff = capture_time-prev_capture_time-shutter
            print( '  Files {0} and {1} were captured {2} seconds apart'.format(prev_filename, filename, capt_time_diff))
            if( capt_time_diff > args.threshold ):
                print( '===============')
                if len(exp_stack) < args.minexp:
                    print( 'Too few exposures to create an exposure stack')
                else:  
                    move_exp_to_dir(exp_stack, stack_no)
                    stack_no = stack_no+1
                exp_stack=list()
        sys.stdout.flush()
        prev_capture_time = capture_time
        prev_shutter=shutter
        prev_filename=filename
        exp_stack.append(filename)
        #print( lines )