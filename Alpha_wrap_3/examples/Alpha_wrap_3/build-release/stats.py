import csv
from collections import defaultdict, namedtuple
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import pickle
import os
from itertools import chain

experiment_keys = ("model", "alpha", "offset", "rays")
ray_keys = ("i", "source", "target", "intersection", "algos")
algo_keys = ("name", "parameter", "time", "steps")
step_keys = ("type", "point")

Experiment = namedtuple("Experiment", experiment_keys)
Ray = namedtuple("Ray", ray_keys)
Algo = namedtuple("Algo", algo_keys)
Step = namedtuple("Step", step_keys)

raw_step_keys = ("ray", "algo", "type", "x", "y", "z")

def convert(k,v):
    match k:
        case "ray":
            return int(v)
        case "x" | "y" | "z":
            return float(v)
        case _:
            return v

def row_point(row_data):
    return np.array((row_data["x"], row_data["y"], row_data["z"]))

def row_step(row_data):
    return {"type": row_data["type"],
            "point": row_point(row_data)}

def algo_with_parameter(algo):
    match algo["name"].split("-"):
        case ("relaxed", omega):
            algo["name"] = "relaxed"
            algo["parameter"] = float(omega)
        case _:
            algo["parameter"] = ""
    return algo

def make_ray(i):
    return {"i" : i, "algos" : []}

def make_algo(name):
    return {"name" : name, "steps" : []}

def read_steps(file_path):
    rays = []
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader)
        ray = make_ray(0)
        algo = make_algo("")
        for row in reader:
            row_data = {k:convert(k,v) for k,v in zip(raw_step_keys, row)}
            # new ray
            if row_data["ray"] != ray["i"]:
                rays.append(Ray(**ray))
                ray = make_ray(row_data["ray"])
            # new algo
            if row_data["algo"] != algo["name"]:
                if algo["name"] and algo["name"] != "tracer":
                    algo = Algo(**algo_with_parameter(algo))
                    ray["algos"].append(algo)
                algo = make_algo(row_data["algo"])
            match row_data:
                # ray data
                # source
                case {"algo": "tracer", "type": "source", **rest}:
                    ray["source"] = row_point(row_data)
                # target
                case {"algo": "tracer", "type": "target", **rest}:
                    ray["target"] = row_point(row_data)
                # intersection
                case {"algo": "tracer", "type": "intersection", **rest}:
                    ray["intersection"] = row_point(row_data)
                # algo data
                # time
                case {"type": "time", **rest}:
                    algo["time"] = row_data["x"]
                # step data
                case _:
                    algo["steps"].append(Step(**row_step(row_data)))
    return rays
                    
def number_of_steps(algo):
    return len(algo.steps)

def time_per_step(algo):
    return algo.time / number_of_steps(algo)

def ray_performence(ray):
    return [(algo.name, algo.parameter, number_of_steps(algo), time_per_step(algo))
            for algo in ray.algos]
    
def relaxed_trace(algo, ray):
    source = ray.source
    relaxed_steps = []
    normal_steps = []
    overshot = False
    for step in algo.steps:
        if step.type == "overstep":
            overshot = step.point
        if overshot is False:
            relaxed_steps.append(step.point)
        else:
            normal_steps.append(step.point)
    return (source, relaxed_steps, overshot, normal_steps)

def ray_segment(ray):
    return ray.target - ray.source

def ray_dir(ray):
    segment = ray_segment(ray)
    return segment / np.linalg.norm(segment)

def distance(v,w):
    return float(np.linalg.norm(v - w))

def get_algo(ray, name, param=""):
    def fltr(algo):
        return algo.name == name and algo.parameter == param
    return next(filter(fltr, ray.algos))

def relaxed_coherence(algo, ray):
    source, relaxed, overshot, normal = relaxed_trace(algo, ray)
    omega = algo.parameter
    rdir = ray_dir(ray)

    d = np.linalg.norm(overshot - source) / omega
    vanilla_computed = source + d * rdir
    overshot_computed = source + d * omega * rdir
    overshot_data = overshot
    vanilla_data = normal[0]
    vanilla_recovered_computed = overshot + d * (1 - omega) * rdir

    vanilla_rm = get_algo(ray, "sphere-marching").steps[0].point

    def compare(a,b):
        print(a, b, distance(a,b))
    
    compare(vanilla_computed, vanilla_data)
    compare(vanilla_computed, vanilla_rm)
    compare(vanilla_computed, vanilla_recovered_computed)
    compare(overshot_computed, overshot_data)


def output_steps(steps, filename):
    with open(filename, "w") as fl:
        for step in steps:
            point = list(map(str, step.point))
            print(" ".join(point), file=fl)

def step_param(step_p, source, unit_dir):
    return float(np.dot(unit_dir, (step_p - source)))
            
def steps_params(steps, ray):
    source = ray.source
    unit_dir = ray_dir(ray)
    return [step_param(step.point, source, unit_dir) for step in steps]

def algos_steps_plot(algos, ray, filename):
    fig = plt.figure()
    ax = fig.add_subplot()
    xmax = 0
    linestyle = "-"
    for algo in algos:
        params = [0] + steps_params(algo.steps, ray)
        if algo.name == "sphere-marching":
            linestyle = "--"
        plt.plot(range(len(params)), params, label=algo.name, linestyle=linestyle)
        xmax = max(xmax, len(params) - 1)
    target = step_param(ray.target, ray.source, ray_dir(ray))
    intersection = step_param(ray.intersection, ray.source, ray_dir(ray))
    # plt.hlines(0, label="source",
    #            xmin=0, xmax=xmax, color="red", linestyles="-")
    # plt.hlines(target, label="target",
    #            xmin=0, xmax=xmax, color="red", linestyles="--")
    # plt.hlines(intersection, label="intersection",
    #            xmin=0, xmax=xmax, color="black", linestyles=":")
    ax.set_xticks(range(xmax + 1))
    yticks = ax.set_yticks((0, target, intersection),
                           labels=("source", "target", "intersection"))
    def style_tick(tick, color, style):
        tick.label1.set_color(color)
        tick.tick1line.set_markeredgecolor(color)
        tick.gridline.set(visible=True, linestyle=style, linewidth=1.0, color=color)
    style_tick(yticks[0], "grey", "-")
    style_tick(yticks[1], "grey", "-")
    style_tick(yticks[2], "black", ":")
    plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(filename)

# difference between intersection points
# difference between steps by old and new sphere marching
# step length
# number of steps

def ray_algo(ray, name, param=None):
    for algo in ray.algos:
        if param is None:
            if algo.name == name:
                return algo
        else:
            if (algo.name, algo.parameter) == (name, param):
                return algo

def algo_steps_number(algo):
    return len(algo.steps)

def algo_time_per_step(algo):
    return algo.time / algo_steps_number(algo)

def algo_distribution(rays, algo, feature, param=None):
    return np.array([feature(ray_algo(r, algo, param)) for r in rays])    

def experiment_cache_file(e):
    return experiment_folder(e) + "steps.pkl"

def experiment_folder(e):
    return f"results/{e.model}_{e.alpha}_{e.offset}/"

def read_experiment_cache(e):
    with open(experiment_cache_file(e), "rb") as f:
        return pickle.load(f)

def write_experiment_cache(e, rays):
    with open(experiment_cache_file(e), "wb") as f:
        pickle.dump(rays, f)

def read_experiment(e):
    print(f"Reading {e}")
    if os.path.exists(experiment_cache_file(e)):
        rays = read_experiment_cache(e)
    else:
        rays = read_steps(experiment_folder(e) + "steps.csv")
        write_experiment_cache(e, rays)
    return Experiment(e.model, e.alpha, e.offset, rays)

def clip(a):
    "Remove extreme outliers."
    q = np.quantile(a, 0.99)
    return np.clip(a,a_min=0,a_max=q)

def plot_distributions(rays, algos, feature_getter, title, feature_name, output_file):
    distributions = np.array([clip(algo_distribution(rays,
                                                     algo.name,
                                                     feature_getter,
                                                     algo.parameter))
                              for algo in algos])
    x_rng = (distributions.min(), distributions.max())
    fig = plt.figure()
    ax = fig.add_subplot()
    colors = list(mpl.colormaps["tab10"].colors)[:len(algos)]
    bins = int(np.ceil(x_rng[1]))
    histograms = [np.histogram(distribution, range=x_rng, bins=bins)
                  for distribution in distributions]
    y_rng = (0, np.array([h[0] for h in histograms]).max() * 1.1)
    for distribution, histogram, algo, color in zip(distributions, histograms, algos, colors):
        counts, bins = histogram
        plt.stairs(counts, bins, label=algo.name, color=color)
        ax.fill_between(bins[:-1], counts, alpha=0.15, step="post", color=color)
        # median = np.median(distribution)
        # plt.vlines(median, 0, 1700, color=color, linestyle = ":")
        mean = np.mean(distribution)
        plt.vlines(mean, *y_rng, color=color, linestyle = "--")
    plt.legend()
    plt.title(title)
    plt.xlabel(feature_name)
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(output_file)

experiment = read_experiment(Experiment("bones", 80, 80, []))
algos = [algo for algo in experiment.rays[0].algos if algo.name != "relaxed" or algo.parameter == 1.1]

def experiment_from_folder(folder):
    split = folder.split("_")
    model = split[0]
    alpha = int(split[1])
    offset = int(split[2])
    expe = Experiment(model, alpha, offset, [])
    return read_experiment(expe)

experiments = [experiment_from_folder(f) for f in os.listdir("results") if f.split("_")[0] == "bones"]
rays = list(chain(*[e.rays for e in experiments]))

plot_distributions(rays, algos, algo_steps_number, "Ray steps distribution", "steps", "ray_steps.png")

# plot_distributions(rays, algos, algo_steps_number, "Ray steps distribution", "steps", "ray_steps.png")

# r = rays[0]
# print(ray_ranking(r))
# for a in r.algos:
#     if a.name == "relaxed":
#         print(relaxed_coherence(a, r))
        # output_steps(a.steps, f"relaxed-{a.parameter}.xyz")
    # print(steps_params(a.steps, r))
# algos_steps_plot([a for a in r.algos if not a.name == "relaxed"], r, "algo_steps-no-relaxed.png")
# algos_steps_plot(r.algos, r, "algo_steps.png")

