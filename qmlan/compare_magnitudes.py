import argparse
from obspy.core.event import read_events
from obspy.core.inventory.inventory import read_inventory
from qmlan.utils import get_magnitude, get_origin, get_hypocentral_distance
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats


def extract_station_inventories(file_name):
    inventory = read_inventory(file_name, format='StationXML')
    station_inventories = {}
    for network in inventory:
        for station in network:
            station_inventories[station.code] = station
    return station_inventories


def compare_magnitude(catalog, type1, type2, method_id1=None, method_id2=None):
    magnitudes1 = []
    magnitudes2 = []
    for event in catalog:
        magnitude1 = get_magnitude(event, magnitude_type=type1, method_id=method_id1)
        magnitude2 = get_magnitude(event, magnitude_type=type2, method_id=method_id2)
        magnitudes1.append(magnitude1)
        magnitudes2.append(magnitude2)
    return magnitudes1, magnitudes2

def  ML2Mw(ML, method='Grunthal'):
    ML2 = ML ** 2
    if method == 'Grunthal':
        Mw = 0.53 + 0.646 * ML + 0.0376 * ML2 # Grunthal 2009
    else:
        Mw = 1.2 + 0.28 * ML + 0.06 * ML2
    return Mw


def func2 (x, a, b):
    return a + b * x

def fit2(xx, yy):
    params, pcov = curve_fit(func2, xx, yy)
    a, b= params
    print(f"Fit parameters: {params}")
    print(f"Fit pcov: {pcov}")
    return a, b

def func3 (x, a, b, c):
    return a + b * x + c * x ** 2

def fit3(xx, yy):
    params, pcov = curve_fit(func3, xx, yy)
    a, b, c = params
    print(f"Fit parameters: {params}")
    print(f"Fit pcov: {pcov}")
    return a, b, c

def main():
    parser = argparse.ArgumentParser(
        prog='compare_magnitudes',
        description='Program plot two type magnitudes relation')
    parser.add_argument('-c', '--catalog', metavar='catalog_name.xml', help='Catalog')
    parser.add_argument('-s', '--station_inventory', metavar='station_inventory.xml', help='Station inventory')
    parser.add_argument('-o', '--output', metavar='output_catalog.xml', help='Output catalog')
    parser.add_argument('-1', '--type1', metavar='Mw', help='Type of first magnitudes')
    parser.add_argument('-2', '--type2', metavar='ML', help='Type of second magnitudes')
    parser.add_argument('-3', '--method1', metavar='', help='Method of first magnitudes')
    parser.add_argument('-4', '--method2', metavar='', help='method of second magnitudes')
    parser.add_argument('-i', '--input_format', default='QUAKEML', help='Input catalog format')
    parser.add_argument('-m', '--marker_size', default='3', help='Marker size')
    parser.add_argument('-e', '--error_limit', default='0.4', help='Error limit')
    parser.add_argument('-x', '--max_correct', default='3.2', help='Maximum corrected amplitudes')
    args = parser.parse_args()
    print(f"Start reading {args.catalog}")
    stations = {}
    if args.station_inventory is not None:
        stations = extract_station_inventories(args.station_inventory)
    catalog = read_events(args.catalog, format=args.input_format)
    print(f"Read {len(catalog.events)} events")
    if args.station_inventory is not None:
        distances = []
        for idx, event in enumerate(catalog.events):
            origin = get_origin(event)
            if origin is None:
                continue
            for station, station_inventory in stations.items():
                distance = get_hypocentral_distance(origin, station_inventory)
                distances.append(distance[0] / 1000.0)
        fig0, ax0 = plt.subplots()
        counts, bins = np.histogram(distances, bins=100)
        ax0.plot(bins[:-1], counts, 'b-', linewidth=2)
        plt.show(block=True)
    magnitudes1, magnitudes2 = compare_magnitude(catalog, args.type1, args.type2,
                                                 method_id1=args.method1, method_id2=args.method2)
    x = []
    y = []
    for magnitude1, magnitude2 in zip(magnitudes1, magnitudes2):
        if magnitude1 is not None and magnitude2 is not None:
            x.append(magnitude1.mag)
            y.append(magnitude2.mag)
    print(f"Found {len(x)} magnitude pairs")
    marker_size = int(args.marker_size)
    error_limit = float(args.error_limit)
    max_correct = float(args.max_correct)
    fig1, ax1 = plt.subplots()
    if args.type1 == 'E' or args.type1 == 'Energy':
        type1 = 'E'
        ax1.semilogx(x, y, 'r.', markersize=marker_size, label="Mw vs Energy")
        ax1.set_xlabel("E [J]")
        ax1.set_ylabel(args.type2 if args.method2 is None else f"{args.type2} ({args.method2})")
        xxx = np.log10(np.array(x))
        yyy = np.array(y)
        # result = stats.linregress(xxx, yyy, 'two-sided')
        # xx = np.sort(np.array(x))
        # plt.semilogx(xx, result.intercept + result.slope * np.log10(xx), 'g--',
        #              label=r'$Mw={:.2}+{:.2} \cdot {}{} [J]$'.format(result.intercept, result.slope,
        #                                                        r"log_{10}", type1))
        pars = fit2(yyy, xxx)
        yy = np.sort(yyy)
        # ax1.semilogx(np.power(10.0, func2(yy, *pars)), yy, 'b--', linewidth=1,
        #              label=r"${}={:.2}+{:.2}\cdot Mw+{:.2}\cdot Mw^2$".format("log_{10}E", pars[0], pars[1])
        #              )
        ax1.semilogx(np.power(10.0, func2(yy, *pars)), yy, 'b--', linewidth=1,
                     label=r"${}={:.2}+{:.2}\cdot Mw$".format("log_{10}E", pars[0], pars[1])
                     )
        ax1.legend()
    elif args.type2 == 'E' or args.type2 == 'Energy':
        ax1.semilogy(x, y, 'r.', markersize=marker_size)
        ax1.set_xlabel(args.type1 if args.method1 is None else f"{args.type1} ({args.method1})")
        ax1.set_ylabel(args.type2 if args.method2 is None else f"{args.type2} ({args.method2})")
    else:
        min_mag, max_mag = min(x), max(x)
        ax1.plot([min_mag, max_mag], [min_mag, max_mag], 'k:', linewidth=0.5, label=r"$Mw=ML$")
        ax1.plot(x, y, 'r.', markersize=marker_size, label="Estimated magnitudes")
        xx = np.sort(np.array(x))
        yy = ML2Mw(xx)
        ax1.plot(xx, yy, 'b--', linewidth=0.5, label="Gruntal (2009)")
        xxx = np.array(x)
        yyy = np.array(y)
        pars = fit3(xxx, yyy)
        ax1.plot(xx, func3(xx, *pars), 'g--', linewidth=1,
                 label=r"$Mw={:.2}+{:.2}\cdot ML+{:.2}\cdot ML^2$".format(pars[0], pars[1], pars[2])
                 )
        bx = []
        by = []
        e = []
        for magnitude1, magnitude2, event in zip(magnitudes1, magnitudes2, catalog.events):
            if magnitude1 is not None and magnitude2 is not None:
                if magnitude2.mag < max_correct:
                    continue
                if magnitude2.mag - func3(magnitude1.mag, *pars) < error_limit:
                    continue
                bx.append(magnitude1.mag)
                by.append(magnitude2.mag)
                e.append(event)
        print(f"{len(e)} blue events")
        if args.output:
            catalog.events = e
            catalog.write(args.output, format=args.input_format)
        ax1.plot(bx, by, 'b+', markersize=marker_size + 2, label="Outlier Mw values")
        ax1.legend()
        ax1.set_xlabel(args.type1 if args.method1 is None else f"{args.type1} ({args.method1})")
        ax1.set_ylabel(args.type2 if args.method2 is None else f"{args.type2} ({args.method2})")
        fig2, ax2 = plt.subplots()
        mag = [m.mag for m in magnitudes1 if m is not None]
        min_mag, max_mag = min(mag), max(mag)
        pre_bins = np.arange(np.floor(min_mag - 0.1), np.ceil(max_mag + 0.2), 0.1)
        counts, bins = np.histogram(mag, bins=pre_bins)
        ax2.semilogy(bins[:-1], counts, 'bo', markersize=marker_size, label="No. events")
        ax2.semilogy(bins[:-1], np.cumsum(counts[::-1])[::-1], 'b-', label="Cumulated no. events")
        # ax2.hist([m.mag for m in magnitudes1 if m is not None], bins=100)
        ax2.set_xlabel(args.type1 if args.method1 is None else f"{args.type1} ({args.method1})")
        ax2.set_ylabel("No events in 0.1 bins")
        ax2.legend()
        fig3, ax3 = plt.subplots()
        mag = [m.mag for m in magnitudes2 if m is not None]
        min_mag, max_mag = min(mag), max(mag)
        pre_bins = np.arange(np.floor(min_mag - 0.1), np.ceil(max_mag + 0.2), 0.1)
        counts, bins = np.histogram(mag, bins=pre_bins)
        ax3.semilogy(bins[:-1], counts, 'bo', markersize=marker_size, label="No. events")
        ax3.semilogy(bins[:-1], np.cumsum(counts[::-1])[::-1], 'b-', label="Cumulated no. events")
        # ax3.hist([m.mag for m in magnitudes2 if m is not None], bins=100)
        ax3.set_xlabel(args.type2 if args.method2 is None else f"{args.type2} ({args.method2})")
        ax3.set_ylabel("No. events in 0.1 bins")
        ax3.legend()
        # Tests
        fig4, ax4 = plt.subplots()
        mag = [m.mag for m,m0 in zip(magnitudes2, magnitudes1)
               if m is not None and (m.mag < max_correct or m0 is None or m.mag < ML2Mw(m0.mag) + error_limit)]
        min_mag, max_mag = min(mag), max(mag)
        pre_bins = np.arange(np.floor(min_mag - 0.1), np.ceil(max_mag + 0.2), 0.1)
        counts, bins = np.histogram(mag, bins=pre_bins)
        ax4.semilogy(bins[:-1], counts, 'bo', markersize=marker_size, label="No. events")
        ax4.semilogy(bins[:-1], np.cumsum(counts[::-1])[::-1], 'b-', label="Cumulated no. events")
        ax4.set_xlabel(args.type2 if args.method2 is None else f"{args.type2} ({args.method2})")
        ax4.set_ylabel("No. events in 0.1 bins")
        ax4.legend()
        # Tests
    plt.show(block=True)


if __name__ == '__main__':
    main()
