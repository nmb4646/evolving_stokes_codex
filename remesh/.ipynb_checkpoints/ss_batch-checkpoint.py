import concurrent.futures
import os
import subprocess
import sys
import threading


def _filter_stderr(pipe):
    noisy_patterns = [
        "MathWorksServiceHost",
        "GLIBCXX_3.4.31",
        "GLIBCXX_3.4.32",
        "CXXABI_1.3.15",
        "libmwflnetwork.so",
        "libmwflcrypto.so",
        "libmwflcryptoutils.so",
        "libmwflcryptoopenssl.so",
        "libmwflurlmgrfactory.so",
    ]

    for line in pipe:
        if not any(pattern in line for pattern in noisy_patterns):
            sys.stderr.write(line)
            sys.stderr.flush()


def run(cmd):
    process = subprocess.Popen(
        cmd,
        shell=True,
        text=True,
        stdout=None,
        stderr=subprocess.PIPE,
    )

    thread = threading.Thread(target=_filter_stderr, args=(process.stderr,))
    thread.daemon = True
    thread.start()

    returncode = process.wait()
    thread.join(timeout=1)

    if returncode != 0:
        print(f"Command exited with code {returncode}", file=sys.stderr)


def last_frame(directory):
    frames = [
        int(file.split("geo")[1].split(".")[0])
        for file in os.listdir(directory)
        if file.startswith("geo") and file.endswith(".mat")
    ]
    return max(frames) if frames else 0


here = "/home/nickbroussinos/Code/evolving_stokes_codex/remesh/"

verbose = "true"
resume = False
supress_outputs = 0

subdivisions = 4
roughness = 0.05
remesh_size = 0
initial_remesh = 0
dt = 0.05
T = 200
k = 10

# Set to a .mat file containing P and M to start from a weird shape.
initial_geometry = ""

param = [roughness]

commands = []
for roughness_value in param:
    run_tag = f"rough_{roughness_value:.4g}_dt_{dt:.4g}_k_{k:.4g}".replace(".", "p")
    dir_path = here + f"data/ss_batch_data/{run_tag}"
    os.makedirs(dir_path, exist_ok=True)
    start = last_frame(dir_path) if resume else 0
    if resume:
        print(f"start: {start}")

    cmd = f"""export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 && \
matlab -nodisplay -nosplash -nodesktop -r "cd '{here}'; verbose = {verbose}; \
supress_outputs = {supress_outputs}; dir = '{dir_path}/'; \
p = struct(); p.start = {start}; p.subdivisions = {subdivisions}; \
p.roughness = {roughness_value}; p.dt = {dt}; p.T = {T}; p.k = {k}; \
p.remesh_size = {remesh_size}; p.initial_remesh = {initial_remesh}; \
p.initial_geometry = '{initial_geometry}'; \
ss_multi; exit" """
    commands.append(cmd)


max_workers = 1
with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(run, cmd) for cmd in commands]
    concurrent.futures.wait(futures)
