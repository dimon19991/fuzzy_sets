import numpy as np

fuzzy_system = {"input_data": {}, "output_data": {}}

def domain(domain_min, domain_max, res):
    domain = np.linspace(domain_min, domain_max, res)
    dom = np.zeros(domain.shape)

def adjust_domain_val(x_val):
    return domain[np.abs(domain - x_val).argmin()]


def create_trapezoidal(name, domain_min, domain_max, res, a, b, c, d):
    t1fs = domain(domain_min, domain_max, res)

    a = adjust_domain_val(a)
    b = adjust_domain_val(b)
    c = adjust_domain_val(c)
    d = adjust_domain_val(d)

    t1fs._dom = np.round(
        np.minimum(np.maximum(np.minimum((t1fs._domain - a) / (b - a), (d - t1fs._domain) / (d - c)), 0), 1),
        t1fs._precision)
    return t1fs

def create_triangular(cls, name, domain_min, domain_max, res, a, b, c):
    t1fs = cls(name, domain_min, domain_max, res)

    a = t1fs._adjust_domain_val(a)
    b = t1fs._adjust_domain_val(b)
    c = t1fs._adjust_domain_val(c)

    if b == a:
        t1fs._dom = np.round(np.maximum((c - t1fs._domain) / (c - b), 0), t1fs._precision)
    elif b == c:
        t1fs._dom = np.round(np.maximum((t1fs._domain - a) / (b - a), 0), t1fs._precision)
    else:
        t1fs._dom = np.round(np.maximum(np.minimum((t1fs._domain - a) / (b - a), (c - t1fs._domain) / (c - b)), 0),
                             t1fs._precision)

    return t1fs

def create_z_linar(name, name_term, start, end):
    t1fs = cls(name, domain_min, domain_max, res)

    a = t1fs._adjust_domain_val(start)
    b = t1fs._adjust_domain_val(end)

    t1fs._dom = np.maximum(np.minimum((b - t1fs._domain) / (b - a), 1), 0)

    return t1fs


def create_s_linar(cls, name, domain_min, domain_max, res, start, end):
    t1fs = cls(name, domain_min, domain_max, res)

    a = t1fs._adjust_domain_val(start)
    b = t1fs._adjust_domain_val(end)

    t1fs._dom = np.minimum(np.maximum((t1fs._domain - a) / (b - a), 0), 1)

    return t1fs

def InputVariable(name, min_, max_, len_):
    fuzzy_system["input_data"][name] = {}
    fuzzy_system["input_data"][name]["max"] = max_
    fuzzy_system["input_data"][name]["min"] = min_
    fuzzy_system["input_data"][name]["len"] = len_




temp = InputVariable('Уровень', 0, 9, 100)
create_z_linar('S', 'Уровень', 2, 4)
# create_trapezoidal('M', 2, 4, 6, 8)
# creatr_s_linar('L', 6, 8)

# humidity = InputVariable('Расход', 0, 0.5, 100)
# humidity.add_z_linar('SS', 0.2, 0.3)
# humidity.add_trapezoidal('MM', 0.15, 0.25, 0.35, 0.45)
# humidity.add_s_linar('LL', 0.3, 0.4)
#
# motor_speed = OutputVariable('Speed', 0, 10, 100)
# motor_speed.add_z_linar('Slow', 1, 5)
# motor_speed.add_trapezoidal('Moderate', 1, 4, 6, 9)
# motor_speed.add_s_linar('Fast', 5, 9)
#
# # system = FuzzySystem()
# # system.add_input_variable(temp)
# # system.add_input_variable(humidity)
# # system.add_output_variable(motor_speed)
#
# system.add_rule(
#     {'Уровень': 'S',
#      'Расход': 'LL'},
#     {'Speed': 'Slow'})
#
# system.add_rule(
#     {'Уровень': 'S',
#      'Расход': 'MM'},
#     {'Speed': 'Slow'})
#
# system.add_rule(
#     {'Уровень': 'S',
#      'Расход': 'SS'},
#     {'Speed': 'Slow'})
#
# system.add_rule(
#     {'Уровень': 'M',
#      'Расход': 'LL'},
#     {'Speed': 'Moderate'})
#
# system.add_rule(
#     {'Уровень': 'M',
#      'Расход': 'MM'},
#     {'Speed': 'Moderate'})
#
# system.add_rule(
#     {'Уровень': 'M',
#      'Расход': 'SS'},
#     {'Speed': 'Moderate'})
#
# system.add_rule(
#     {'Уровень': 'L',
#      'Расход': 'LL'},
#     {'Speed': 'Fast'})
#
# system.add_rule(
#     {'Уровень': 'L',
#      'Расход': 'MM'},
#     {'Speed': 'Fast'})
#
# system.add_rule(
#     {'Уровень': 'L',
#      'Расход': 'SS'},
#     {'Speed': 'Fast'})
#
# output = system.evaluate_output({
#     'Уровень': 9,
#     'Расход': 100
# })
