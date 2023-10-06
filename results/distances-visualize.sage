import csv
def load_distances(distances_report_fname):
    distances = []
    with open(distances_report_fname) as fh:
        reader = csv.reader(fh, delimiter=";")
        for row in reader:
            if row[0] == "round":
                continue
            distances.append([int(c) for d,c in enumerate(row[1:])])
    return distances

from sage.plot.colors import Colormaps
viridis = Colormaps()["viridis"]

def plot_distances(m,
                   t, n, distances_report_fname):
    distances=load_distances(distances_report_fname)
    nr_plots = len(distances)
    plots = []
    plot_joined=True
    for n in range(nr_plots):
        color = viridis.colors[int( viridis.N * n / nr_plots)]
        extra = {}
        if n % 5 != 0:
            extra["linestyle"] = "dashed"
            extra["thickness"] = 0.5
        else:
            extra["thickness"] = 1.1
        graph = list_plot_semilogy(distances[n], color=color, plotjoined=plot_joined, **extra)
        plots.append(graph)
    return sum(plots)

p = plot_distances(9,20,512, "distances-report-mceliece51220.csv")
p = p + text("i = 0", (40,3*10**4), background_color="white")
p = p + text("i = 5", (40,8*10**2), background_color="white")
p = p + text("i = 10", (40,2*10**1), background_color="white")
p = p + text("i = " + str(20+1), (9*20/2,2*10**4), background_color="white")
p.save_image("distances-report-mceliece51220.png")

p = plot_distances(10,50,1024, "distances-report-mceliece102450.csv")
p = p + text("i = 0", (160,6*10**5), background_color="white")
p = p + text("i = 5", (160,3*10**4), background_color="white")
p = p + text("i = 10", (160,8*10**2), background_color="white")
p = p + text("i = " + str(50+1), (10*50/2,4.2*10**4), background_color="white")
p.save_image("distances-report-mceliece102450.png")
