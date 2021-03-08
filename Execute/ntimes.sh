#!/bin/bash
# My example bash script
echo "Ejecucion multiple"
count=20
for i in $(seq $count); do
	source venv/bin/activate
    ./HybridH1 
    python insert.py
    deactivate
done