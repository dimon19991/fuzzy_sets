import numpy as np

fuzzy_system = {"input_data": {}, "output_data": {}, "rules": [[], []], "evaluate_output": {}, "fazification": {}}


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


def get_mu(name):
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

        y = ((dod_par[maxx]-dod_par[minn])*(x_val-dodain_par[minn]))/(dodain_par[maxx]-dodain_par[minn]) + dod_par[minn]

        fuzzy_system["fazification"][name][i] = y


def fazification():
    for i in fuzzy_system["evaluate_output"].keys():
        get_mu(i)


input_1 = InputVariable('Уровень', 0, 9, 1000)
create_z_linar('S', 2, 4, input_1)
create_trapezoidal('M', 2, 4, 6, 8, input_1)
create_s_linar('L', 6, 8, input_1)

input_2 = InputVariable('Расход', 0, 0.5, 1000)
create_z_linar('SS', 0.2, 0.3, input_2)
create_trapezoidal('MM', 0.15, 0.25, 0.35, 0.45, input_2)
create_s_linar('LL', 0.3, 0.4, input_2)

output_1 = OutputVariable('Speed', 0, 10, 1000)
create_z_linar('Slow', 1, 5, output_1)
create_trapezoidal('Moderate', 1, 4, 6, 9, output_1)
create_s_linar('Fast', 5, 9, output_1)

add_rule(
    {'Уровень': 'S',
     'Расход': 'LL'},
    {'Speed': 'Slow'})

add_rule(
    {'Уровень': 'S',
     'Расход': 'MM'},
    {'Speed': 'Slow'})

add_rule(
    {'Уровень': 'S',
     'Расход': 'SS'},
    {'Speed': 'Slow'})

add_rule(
    {'Уровень': 'M',
     'Расход': 'LL'},
    {'Speed': 'Moderate'})

add_rule(
    {'Уровень': 'M',
     'Расход': 'MM'},
    {'Speed': 'Moderate'})

add_rule(
    {'Уровень': 'M',
     'Расход': 'SS'},
    {'Speed': 'Moderate'})

add_rule(
    {'Уровень': 'L',
     'Расход': 'LL'},
    {'Speed': 'Fast'})

add_rule(
    {'Уровень': 'L',
     'Расход': 'MM'},
    {'Speed': 'Fast'})

add_rule(
    {'Уровень': 'L',
     'Расход': 'SS'},
    {'Speed': 'Fast'})

evaluate_output({
    'Уровень': 2.5,
    'Расход': 0.4
})

fazification()

print()
