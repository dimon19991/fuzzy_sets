import matplotlib.pyplot as plt
import numpy as np

fuzzy_system = {"input_data": {}, "output_data": {}, "rules": [[], [], []], "evaluate_output": {}, "fazification": {},
                "activation": [], "accumulation": {"res": {"dom": [], "domain": []}}, "output": {}}


def domain(domain_min, domain_max, res):
    domain = np.linspace(domain_min, domain_max, res)
    dom = np.zeros(domain.shape)

    return domain, dom


def adjust_domain_val(x_val, dodain_par):
    return dodain_par[np.abs(dodain_par - x_val).argmin()]


def create_trapezoidal(name, a, b, c, d, term_data):
    this_term = fuzzy_system[term_data[0]][term_data[1]]

    t1fs = domain(this_term["min"], this_term["max"], this_term["len"])
    this_term[name] = {}
    this_term[name]["domain"] = t1fs[0]
    this_term[name]["dom"] = t1fs[1]

    a = adjust_domain_val(a, this_term[name]["domain"])
    b = adjust_domain_val(b, this_term[name]["domain"])
    c = adjust_domain_val(c, this_term[name]["domain"])
    d = adjust_domain_val(d, this_term[name]["domain"])

    this_term[name]["dom"] = np.round(
        np.minimum(
            np.maximum(np.minimum((this_term[name]["domain"] - a) / (b - a), (d - this_term[name]["domain"]) / (d - c)),
                       0), 1), 2)


def create_triangular(name, a, b, c, term_data):
    this_term = fuzzy_system[term_data[0]][term_data[1]]

    t1fs = domain(this_term["min"], this_term["max"], this_term["len"])
    this_term[name] = {}
    this_term[name]["domain"] = t1fs[0]
    this_term[name]["dom"] = t1fs[1]

    a = adjust_domain_val(a, this_term[name]["domain"])
    b = adjust_domain_val(b, this_term[name]["domain"])
    c = adjust_domain_val(c, this_term[name]["domain"])

    if b == a:
        t1fs._dom = np.round(np.maximum((c - this_term[name]["domain"]) / (c - b), 0), 2)
    elif b == c:
        t1fs._dom = np.round(np.maximum((this_term[name]["domain"] - a) / (b - a), 0), 2)
    else:
        t1fs._dom = np.round(
            np.maximum(np.minimum((this_term[name]["domain"] - a) / (b - a), (c - this_term[name]["domain"]) / (c - b)),
                       0), 2)


def create_z_linar(name, start, end, term_data):
    this_term = fuzzy_system[term_data[0]][term_data[1]]

    t1fs = domain(this_term["min"], this_term["max"], this_term["len"])
    this_term[name] = {}
    this_term[name]["domain"] = t1fs[0]
    this_term[name]["dom"] = t1fs[1]

    a = adjust_domain_val(start, this_term[name]["domain"])
    b = adjust_domain_val(end, this_term[name]["domain"])

    this_term[name]["dom"] = np.round(np.maximum(np.minimum((b - this_term[name]["domain"]) / (b - a), 1), 0), 2)


def create_s_linar(name, start, end, term_data):
    this_term = fuzzy_system[term_data[0]][term_data[1]]

    t1fs = domain(this_term["min"], this_term["max"], this_term["len"])
    this_term[name] = {}
    this_term[name]["domain"] = t1fs[0]
    this_term[name]["dom"] = t1fs[1]

    a = adjust_domain_val(start, this_term[name]["domain"])
    b = adjust_domain_val(end, this_term[name]["domain"])

    this_term[name]["dom"] = np.round(np.minimum(np.maximum((this_term[name]["domain"] - a) / (b - a), 0), 1), 2)


def InputVariable(name, min_, max_, len_):
    fuzzy_system["input_data"][name] = {}
    fuzzy_system["input_data"][name]["max"] = max_
    fuzzy_system["input_data"][name]["min"] = min_
    fuzzy_system["input_data"][name]["len"] = len_

    return "input_data", name


def OutputVariable(name, min_, max_, len_):
    fuzzy_system["output_data"][name] = {}
    fuzzy_system["output_data"][name]["max"] = max_
    fuzzy_system["output_data"][name]["min"] = min_
    fuzzy_system["output_data"][name]["len"] = len_

    return "output_data", name


def add_rule(input_rile, output_rule):
    fuzzy_system["rules"][0].append([])
    for key, value in input_rile.items():
        fuzzy_system["rules"][0][-1].append([key, value])

    fuzzy_system["rules"][1].append([])
    for key, value in output_rule.items():
        fuzzy_system["rules"][1][-1].append([key, value])


def evaluate_output(res):
    for key, value in res.items():
        fuzzy_system["evaluate_output"][key] = value


def fazification():
    for name in fuzzy_system["evaluate_output"].keys():
        fuzzy_system["fazification"][name] = {}
        x_val = fuzzy_system["evaluate_output"][name]
        kays = list(fuzzy_system["input_data"][name].keys())[3:]
        for i in kays:
            dodain_par = fuzzy_system["input_data"][name][i]["domain"]
            dod_par = fuzzy_system["input_data"][name][i]["dom"]

            nearest = np.abs(dodain_par - x_val).argmin()
            nearest_m_1 = np.abs(x_val - dodain_par[nearest - 1])
            nearest_p_1 = np.abs(x_val - dodain_par[nearest + 1])

            if nearest_m_1 < nearest_p_1:
                minn = nearest - 1
                maxx = nearest
            else:
                minn = nearest
                maxx = nearest + 1

            y = ((dod_par[maxx] - dod_par[minn]) * (x_val - dodain_par[minn])) / (dodain_par[maxx] - dodain_par[minn]) + \
                dod_par[minn]

            fuzzy_system["fazification"][name][i] = y


def aggregation():
    for i in fuzzy_system["rules"][0]:
        val_1 = fuzzy_system["fazification"][i[0][0]][i[0][1]]
        val_2 = fuzzy_system["fazification"][i[1][0]][i[1][1]]
        fuzzy_system["rules"][2].append(np.minimum(val_1, val_2))


def activation():
    k = 0
    for i in fuzzy_system["rules"][1]:
        domaine = fuzzy_system["output_data"][i[0][0]][i[0][1]]["dom"]
        agr = fuzzy_system["rules"][2][k]
        fuzzy_system["activation"].append(np.minimum(domaine, agr))
        k += 1


def accumulation():
    dom = np.zeros(list(fuzzy_system["output_data"].values())[-1]["len"])
    for i in fuzzy_system["activation"]:
        dom = np.maximum(i, dom)
    fuzzy_system["accumulation"]["res"]["dom"] = dom
    fuzzy_system["accumulation"]["res"]["domain"] = list(list(fuzzy_system["output_data"].values())[-1].values())[-1][
        "domain"]


def plot_system():
    for idx, var_name in enumerate(fuzzy_system["input_data"]):
        plot_variable(data=list(fuzzy_system["input_data"][var_name].items())[3:], name=var_name)

    for idx, var_name in enumerate(fuzzy_system["output_data"]):
        plot_variable(data=list(fuzzy_system["output_data"][var_name].items())[3:], name=var_name)

    for idx, var_name_1 in enumerate(fuzzy_system["accumulation"]):
        plot_variable(data=[*list(fuzzy_system["output_data"][var_name].items())[3:],
                            *list(fuzzy_system["accumulation"].items())], name=var_name_1)


def plot_variable(data, name, ax=None, show=True):
    if ax == None:
        ax = plt.subplot(111)

    for n, s in data:
        if len(list(filter(lambda i: i[0] == "res", data))) > 0:
            if n == "res":
                ax.plot(s["domain"], s["dom"], label=n, color="blue")
            else:
                ax.plot(s["domain"], s["dom"], label=n, color="grey")
        else:
            ax.plot(s["domain"], s["dom"], label=n)
    # Shrink current axis by 20%
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width * 0.8, pos.height])
    ax.grid(True, which='both', alpha=0.4)
    ax.set_title(name)
    ax.set(xlabel='x', ylabel='$\mu (x)$')

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if show:
        plt.show()


def defazification():
    sum_1 = 0
    sum_2 = 0
    for i in range(len(fuzzy_system["accumulation"]["res"]["dom"])):
        sum_1 += (fuzzy_system["accumulation"]["res"]["dom"][i] * fuzzy_system["accumulation"]["res"]["domain"][i])
        sum_2 += fuzzy_system["accumulation"]["res"]["dom"][i]
    fuzzy_system["output"]["centre_of_gravity"] = sum_1 / sum_2

    i = 0
    j = shape - 1
    sum_l = 0
    sum_r = 0
    while i != j:
        if sum_r > sum_l:
            sum_l += fuzzy_system["accumulation"]["res"]["dom"][j]
            j -= 1
        else:
            sum_r += fuzzy_system["accumulation"]["res"]["dom"][i]
            i += 1

    fuzzy_system["output"]["centre_of_area"] = (fuzzy_system["accumulation"]["res"]["dom"][i] *
                                                fuzzy_system["accumulation"]["res"]["domain"][i]) / \
                                               fuzzy_system["accumulation"]["res"]["dom"][i]


shape = 1000

input_1 = InputVariable('Уровень', 0, 9, shape)
create_z_linar('S', 2, 4, input_1)
create_trapezoidal('M', 2, 4, 6, 8, input_1)
create_s_linar('L', 6, 8, input_1)

input_2 = InputVariable('Расход', 0, 0.5, shape)
create_z_linar('SS', 0.2, 0.3, input_2)
create_trapezoidal('MM', 0.15, 0.25, 0.35, 0.45, input_2)
create_s_linar('LL', 0.3, 0.4, input_2)

output_1 = OutputVariable('Приток', 0, 0.5, shape)
create_z_linar('Ss', 0.20, 0.25, output_1)
create_trapezoidal('Mm', 0.20, 0.25, 0.35, 0.40, output_1)
create_s_linar('Ll', 0.35, 0.40, output_1)

add_rule(
    {'Уровень': 'S',
     'Расход': 'LL'},
    {'Приток': 'Ll'})

add_rule(
    {'Уровень': 'S',
     'Расход': 'MM'},
    {'Приток': 'Ll'})

add_rule(
    {'Уровень': 'S',
     'Расход': 'SS'},
    {'Приток': 'Mm'})

add_rule(
    {'Уровень': 'M',
     'Расход': 'LL'},
    {'Приток': 'Ll'})

add_rule(
    {'Уровень': 'M',
     'Расход': 'MM'},
    {'Приток': 'Mm'})

add_rule(
    {'Уровень': 'M',
     'Расход': 'SS'},
    {'Приток': 'Mm'})

add_rule(
    {'Уровень': 'L',
     'Расход': 'LL'},
    {'Приток': 'Mm'})

add_rule(
    {'Уровень': 'L',
     'Расход': 'MM'},
    {'Приток': 'Ss'})

add_rule(
    {'Уровень': 'L',
     'Расход': 'SS'},
    {'Приток': 'Ss'})

evaluate_output({
    'Уровень': 2.5,
    'Расход': 0.4
})

fazification()

aggregation()

activation()

accumulation()

defazification()

plot_system()

print()
