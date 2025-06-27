#!/bin/bash
# Clean up messy temp files
find . -name "._*" -delete
find . -name ".DS_Store" -delete
find . -name "*.err" -delete
find . -name "*.out" -delete
find . -name "*.log" -delete