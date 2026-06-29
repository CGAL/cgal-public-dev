import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict
from functools import partial
import subprocess
import json
import multiprocessing
from itertools import cycle, permutations

root = "results"

def lines(file):
    proc = subprocess.run(["wc", "-l", file], encoding="utf-8", stdout=subprocess.PIPE)
    return int(proc.stdout.split()[0])

def no_ext(file):
    return os.path.splitext(os.path.basename(file))[0]

def expe_path(dirname, *args):
    os.path.join(root, dirname, *args)

def write_cache(data, path):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(path, "w") as f:
        f.write(json.dumps(data, indent=4))
    
def read_cache(path):
    with open(path, "r") as f:
        return json.load(f)
    
def expe_counts(dirname):
    data_dir = os.path.join(root, dirname, "data")
    cache = os.path.join(root, dirname, "cache", "counts")
    if os.path.exists(cache):
        return read_cache(cache)
    else:
        files = map(partial(os.path.join, root, dirname, "data"), os.listdir(data_dir))
        data = filter(lambda file: os.path.splitext(os.path.basename(file))[1] != ".off", files)
        output = {no_ext(file) + "s": lines(file) for file in data}
        write_cache(output, cache)
        output = defaultdict(int)
        return output

def expe_inputs(dirname):
    split = dirname.split("_")
    split2 = split[-1].split("-")
    return {"mesh" : split[0],
            "alpha" : int(split[1]),
            "offset" : int(split2[0]),
            "omega" : float(split2[1])}

def expe_data(dirname):
    return expe_inputs(dirname) | expe_counts(dirname)

def is_expe_path(dirname):
    is_dir = os.path.isdir(os.path.join(root, dirname))
    split = dirname.split("_")
    split2 = split[-1].split("-")
    split3 = split2[-1].split(".")
    has_format = split[1].isdigit() and split2[0].isdigit() and split3[0].isdigit() and split3[1].isdigit()
    return is_dir and has_format

def expe_stats(expe):
    rays = expe["segments"]
    relaxed_gain = expe.get("relaxed_steps", 0) * (expe["omega"] - 1)
    relaxed_loss = expe["oversteps"]
    relaxed_net = relaxed_gain - relaxed_loss
    return {
        "rays" : rays,
        "relaxed_gain" : relaxed_gain,
        "relaxed_gain/ray" : relaxed_gain / rays,
        "relaxed_loss" : relaxed_loss,
        "relaxed_loss/ray" : relaxed_loss / rays,
        "relaxed_net" : relaxed_net,
        "relaxed_net/ray" : relaxed_net / rays
    }

pool = multiprocessing.Pool()
expe = list(map(lambda e: e | expe_stats(e),
                pool.map(expe_data, filter(is_expe_path, os.listdir(root)))))

dd = defaultdict(list)

for ex in expe:
    mesh = ex["mesh"]
    dd[mesh].append((ex, ex.pop("mesh"))[0])

# print(dd)

def plot(constant_name, variable, value):
    constants = list(map(lambda e : e[constant_name], filter(lambda e : e[variable] == 300, dd["triceratops.off"])))

    fig, axs = plt.subplots(len(constants), 1)
    # A4 format
    fig.set_size_inches(8.27,11.69)

    colors = cycle(plt.get_cmap('tab20').colors)

    for i, constant in enumerate(constants):
        for mesh in dd:
            data = list(filter(lambda e : e[constant_name] == constant, dd[mesh]))
            x = np.array([e[variable] for e in data])
            y = np.array([e[value] for e in data])
            axs[i].plot(x, y, label=mesh, color=next(colors))
        axs[i].set_title(f'{constant_name}={constant}')

    for ax in axs.flat:
        ax.set(xlabel=variable, ylabel=value)
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        ax.label_outer()

    fig.tight_layout()

    # Make space for the legend
    for ax in axs.flat:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    # Add the legend
    handles, labels = axs.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1, 0.5))

    plt.savefig(f"{value.replace('/', '')}-given-{variable}.png")

for value in ["relaxed_gain/ray", "relaxed_loss/ray", "relaxed_net/ray"]:
    for variable, constant in permutations(["offset", "alpha"]):
        plot(constant, variable, value)
