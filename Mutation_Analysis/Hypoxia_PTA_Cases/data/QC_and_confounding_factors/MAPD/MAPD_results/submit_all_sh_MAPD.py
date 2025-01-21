#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 22:38:09 2024

@author: zhemingan
"""

import os
import sys

if __name__ == "__main__":
    command_folder = sys.argv[1]

# Directory containing the .sh files
directory = command_folder

# List all .sh files in the directory
sh_files = [file for file in os.listdir(directory) if file.endswith(".sh")]

# Submit each .sh file as a SLURM job
for file in sh_files:
    # Get the full path of the .sh file
    file_path = os.path.join(directory, file)
    
    # Generate the SLURM submission command
    submission_command = f"sbatch {file_path}"
    
    # Execute the submission command
    os.system(submission_command)
