#!/bin/bash

directory=`cat directory.txt`

sbatch -p cpu ${directory}/assembly/de_novo.sh
