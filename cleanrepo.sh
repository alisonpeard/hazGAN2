#!/bin/bash
# Clean up messy temp files
rm sbatch_dump/*.err
rm sbatch_dump/*.out
rm logs/*.log
find . -name "._*" -delete
find . -name ".DS_Store" -delete