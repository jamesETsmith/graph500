"""
Hacky workaround to use the original graph500_reference_bfs to generate graph
that can be read by beedrill.
"""
import subprocess
import os
import sys
import argparse
import time

def main():
    #
    # CLI
    #
    parser = argparse.ArgumentParser(
        prog="make_lucata_graphs.py",
        description="Generates graphs useing the graph500 generator prepended with a header",
    )
    parser.add_argument("-e", "--exe", type=str, default="./graph500_generator/src/graph500_reference_bfs", help="Path to graph500_reference_bfs")
    parser.add_argument("-s", "--scale", type=int, default=5, help="Graph scale")
    parser.add_argument(
        "-o", "--output_dir", default="./", help="Output directory to write file."
    )
    parser.add_argument("-n", "--nprocs", type=int, default=6, help="Number of MPI processes to use")
    args = parser.parse_args()



    EXE = args.exe

    scale = args.scale
    final_file = f"{args.output_dir}/graph500_scale{scale}.bin"

    #
    # Derived or hardcoded settings
    #
    edge_factor = 16
    n_vertices = 2**scale
    n_edges = edge_factor * n_vertices
    intermediate_file = f"{final_file}.tmp"

    #
    # Use env variables to change graph generation settings
    #
    env = os.environ.copy()
    env["TMPFILE"] = intermediate_file
    env["REUSEFILE"] = "ON"
    env["SKIP_BFS"] = "ON"

    t_0 = time.time()
    cli_commands = ["mpirun", "-np", f"{args.nprocs}", EXE, f"{scale}"]
    subprocess.run(cli_commands, env=env)
    t_gen = time.time() - t_0


    #
    # Prepend existing file with header
    #

    t_0 = time.time()
    header = f"--format el64 --num_edges {n_edges} --num_vertices {n_vertices} --is_undirected --seed0 2 --seed1 3\n"
    with open(intermediate_file, "rb") as f:
        data = f.read()
    with open(final_file, "w") as f2:
        f2.write(header)
    with open(final_file, "ab") as f2:
        f2.write(data)
    os.remove(intermediate_file)
    t_prepend = time.time() - t_0

    # Time logging
    print(f"Time for graph generation {t_gen:.3e} (s)")
    print(f"Time for graph prepending {t_prepend:.3e} (s)")


if __name__ == "__main__":
    main()