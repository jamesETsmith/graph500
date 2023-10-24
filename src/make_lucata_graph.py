import numpy as np
import subprocess
import os
import sys

EXE = "/home/james/apps/graph500/src/graph500_reference_bfs"

scale = 5
if len(sys.argv) == 1:
    print("Using default scale of 5")
elif len(sys.argv) == 2:
    scale = int(sys.argv[1])

edge_factor = 16
n_vertices = 2**scale
n_edges = edge_factor * n_vertices
intermediate_file = f"graph500_scale{scale}.bin.tmp"
final_file = f"graph500_scale{scale}.bin"

env = os.environ.copy()
# for i, s in enumerate(seeds):
#     env[f"SEED{i}"] = f"{s}"
env["TMPFILE"] = intermediate_file
env["REUSEFILE"] = "ON"
env["SKIP_BFS"] = "ON"

cli_commands = ["mpirun", "-np", "6", EXE, f"{scale}"]
subprocess.run(cli_commands, env=env)

header = f"--format el64 --num_edges {n_edges} --num_vertices {n_vertices} --is_undirected --seed0 2 --seed1 3\n"

with open(intermediate_file, "rb") as f:
    data = f.read()

with open(final_file, "w") as f2:
    f2.write(header)
with open(final_file, "ab") as f2:
    f2.write(data)