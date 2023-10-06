import json
def load_one_report_file(reportfile):
    with open(reportfile) as fh:
        report = json.load(fh)
    return report
import glob

def filter_reports_by_kem(reports, kem):
    return list(filter(lambda report: report['kem'] == kem, reports))

def load_many_report_files(kem, reportglob):
    files = glob.glob(reportglob)
    reports = [load_one_report_file(fname) for fname in files]
    reports = filter_reports_by_kem(reports, kem)
    return reports

from collections import defaultdict
def group_reports_by_error_rate(reports):
    all_error_rates = set([report["error"] for report in reports])
    by_error_rate = defaultdict(list)
    for report in reports:
        by_error_rate[float(report["error"])].append(report)
    return by_error_rate

import numpy as np
def get_succes_rate_by_error_rate(reports):
    by_error_rate = group_reports_by_error_rate(reports)
    params = reports[0]['params']
    
    success_rate_by_error_rate = dict()
    for error_rate in sorted(by_error_rate.keys()):
        successes = [1 if report["success"] else 0 for report in by_error_rate[error_rate]]
        # Add the systematic error from the diagonal to the leak error.
        estimated_error_rate = error_rate # + (0.5-error_rate)/(params["m"]*params["t"])
        success_rate_by_error_rate[estimated_error_rate] = N(np.mean(successes))
    return success_rate_by_error_rate

def get_columns_recovered_by_error_rate(reports):
    by_error_rate = group_reports_by_error_rate(reports)
    params = reports[0]['params']

    columns_recovered_by_error_rate = dict()
    for error_rate in sorted(by_error_rate.keys()):
        columns = [report["recovered_columns"]/report["params"]["n"] for report in by_error_rate[error_rate]]
        # Add the systematic error from the diagonal to the leak error.
        estimated_error_rate = error_rate #+ (0.5-error_rate)/(params["m"]*params["t"])
        columns_recovered_by_error_rate[estimated_error_rate] = N(np.mean(columns))


import functools

@functools.cache
def calculate_V_list(m,t):
    V_list = [1]
    for d in range(1, t*m+1):
        V = V_list[-1] + binomial(m*t,d)
        V_list.append(V)
    return V_list

@functools.cache
def success_probability_from_error(m, t, error):
    V_list = calculate_V_list(m, t)
    S_1 = 0
    S_2 = 0
    for d in range(t*m+1):
        q_1 = RR( binomial(m*t,d) * error^d * (1-error)^(m*t-d) )
        q_2 = RR( (1 - V_list[d] / 2^(t*m))^(2^m) )               
        S_1 += (q_2)^(2^m) * q_1
        S_2 += (q_2) * q_1
    return S_1^(t+1) * S_2^(m*t-(t+1))

def generate_graph(kem):
    print("")
    print(kem)
    reportglob = kem + '/*.json'
    reports = load_many_report_files(kem, reportglob)
    params = reports[0]['params']
    m = params['m']
    t = params['t']
    success_rate_by_error_rate = get_succes_rate_by_error_rate(reports)

    errors = [n(x/100) for x in range(51)]
    predicted_m = [(n(x), n(success_probability_from_error(m, t, x))) for x in errors]
    
    p_mod_m = plot(spline(predicted_m), x, 0, 0.5, color="red", legend_label="model")
    G = p_mod_m
    p_exp = list_plot(list(success_rate_by_error_rate.items()), legend_label="measured")
    G = G + p_exp
    
    # Export the data.
    rates = sorted(success_rate_by_error_rate.items())
    with open("success-probability-report-{}.csv".format(kem), "w") as fh:
        fh.write("error rate;success probability\n")
        for rate, succ in rates:
            fh.write("{};{}\n".format(rate, succ))

    with open("success-probability-estimated-{}.csv".format(kem), "w") as fh:
        fh.write("error rate;success probability\n")
        for rate, succ in predicted_m:
            fh.write("{};{}\n".format(rate, succ))

    # Search for the success rate 0.5.
    crossover = None
    for i in range(len(rates)-1):
        rate1, succ1 = rates[i]
        rate2, succ2 = rates[i+1]
        if succ1 >= 0.5 and succ2 <= 0.5:
            # linear interpolation
            crossover = rate1 + (0.5 - succ1)*(rate2-rate1)/(succ2-succ1)
            # Now round down to nearest three decimal places.
            crossover = n(floor(crossover * 1000) / 1000)
    if crossover is not None:
        relative_entropy = - crossover * log(crossover, 2) - (1-crossover)*log(1-crossover, 2)
        leaked_bits = (1 - relative_entropy) * m*t*m*t
        print("Threshold error:", crossover)
        print("Relative Entropy:", relative_entropy)
        print("Information leak:", 1 - relative_entropy)
        print("Leaked Bits per column:", (1-relative_entropy)*m*t, "of", m*t)
        print("Leaked Bits total:", (1-relative_entropy)*m*t*m*t, "of", m*t*m*t)

        crossover_label = "error = " + str(round(crossover, 3)) + ", prob. = 0.5"
        crossover_label += "\n(rel. entropy: " + str(round(relative_entropy, 3)) + ")"
        p_crossover = line([(crossover, 0), (crossover, 1)], color="green", legend_label=crossover_label)
        G = G + p_crossover

    times = [report['elapsed'] for report in reports if report['success']]
    print("Successful experiments:", len(times), "/", len(reports))
    print("Average time:", np.mean(times), "+/-", np.std(times))
    
    graph_filename = "success-probability-graph-{}.png".format(kem)
    G.save(graph_filename)

generate_graph("mceliece51220")
generate_graph("botan51220")
generate_graph("mceliece102450")
generate_graph("botan102450")
generate_graph("mceliece348864")
generate_graph("botan348864")
generate_graph("mceliece460896")
generate_graph("botan460896")
generate_graph("mceliece6960119")
generate_graph("botan6960119")
generate_graph("mceliece6688128")
generate_graph("botan6688128")
generate_graph("mceliece8192128")
generate_graph("botan8192128")
