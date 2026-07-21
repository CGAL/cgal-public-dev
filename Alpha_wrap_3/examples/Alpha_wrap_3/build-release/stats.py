import csv
from collections import defaultdict
import numpy as np

experiment_keys = ("model", "alpha", "offset", "rays")
ray_keys = ("i", "source", "target", "intersection", "algos")
algo_keys = ("name", "parameter", "time", "steps")
step_keys = ("type", "point")

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
    return {"i" : 0, "algos" : []}

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
                rays.append(ray)
                ray = make_ray(row_data["ray"])
            # new algo
            if row_data["algo"] != algo["name"]:
                if algo["name"] and algo["name"] != "tracer":
                    ray["algos"].append(algo_with_parameter(algo))
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
                    algo["steps"].append(row_step(row_data))
    return rays
                    
rays = read_steps("results/bones_80_80/steps.csv")

def number_of_steps(algo):
    return len(algo["steps"])

def time_per_step(algo):
    return algo["time"] / number_of_steps(algo)

def ray_performence(ray):
    return [(algo["name"], algo["parameter"], number_of_steps(algo), time_per_step(algo))
            for algo in ray["algos"]]
    
def relaxed_trace(algo):
    relaxed_steps = []
    normal_steps = []
    overshot = False
    for step in algo["steps"]:
        if step["type"] == "overstep":
            overshot = step["point"]
        if overshot is False:
            relaxed_steps.append(step["point"])
        else:
            normal_steps.append(step["point"])
    return (relaxed_steps, overshot, normal_steps)

def relaxed_coherence(algo, start):
    relaxed, overshot, normal = relaxed_trace(algo)
    if not overshot is False:
        if relaxed:
            last = relaxed[0]
        else:
            last = start
        return np.linalg.norm(start - normal[0])

ray = rays[0]
start = ray["source"]
for algo in ray["algos"]:
    if algo["name"] == "relaxed":
        print(relaxed_coherence(algo, start))

