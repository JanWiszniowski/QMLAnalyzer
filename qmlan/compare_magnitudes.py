import argparse
from obspy.core.event import read_events
from qmlan.utils import get_magnitude
import matplotlib.pyplot as plt
import numpy as np


def compare_magnitude(catalog, type1, type2, method_id1=None, method_id2=None):
    magnitudes1 = []
    magnitudes2 = []
    for event in catalog:
        magnitude1 = get_magnitude(event, magnitude_type=type1, method_id=method_id1)
        magnitude2 = get_magnitude(event, magnitude_type=type2, method_id=method_id2)
        magnitudes1.append(magnitude1)
        magnitudes2.append(magnitude2)
    return magnitudes1, magnitudes2


def main():
    parser = argparse.ArgumentParser(
        prog='compare_magnitudes',
        description='Program plot two type magnitudes relation')
    parser.add_argument('-c', '--catalog', metavar='catalog_name.xml', help='Catalog')
    parser.add_argument('-o', '--output', metavar='output_catalog.xml', help='Output catalog')
    parser.add_argument('-1', '--type1', metavar='Mw', help='Type of first magnitudes')
    parser.add_argument('-2', '--type2', metavar='ML', help='Type of second magnitudes')
    parser.add_argument('-3', '--method1', metavar='', help='Method of first magnitudes')
    parser.add_argument('-4', '--method2', metavar='', help='method of second magnitudes')
    parser.add_argument('-i', '--input_format', default='QUAKEML', help='Input catalog format')
    parser.add_argument('-m', '--marker_size', default='3', help='Marker size')
    args = parser.parse_args()
    print(f"Start reading {args.catalog}")
    catalog = read_events(args.catalog, format=args.input_format)
    print(f"Read {len(catalog.events)} events")
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
    fig1, ax1 = plt.subplots()
    if args.type1 == 'E':
        ax1.semilogx(x, y, 'r.', markersize=marker_size)
    elif args.type2 == 'E':
        ax1.semilogy(x, y, 'r.', markersize=marker_size)
    else:
        min_mag, max_mag = min(x), max(x)
        ax1.plot([min_mag, max_mag], [min_mag, max_mag], 'k:', linewidth=0.5)
        ax1.plot(x, y, 'r.', markersize=marker_size)
        bx = []
        by = []
        e = []
        for magnitude1, magnitude2, event in zip(magnitudes1, magnitudes2, catalog.events):
            if magnitude1 is not None and magnitude2 is not None:
                if magnitude2.mag < 2.5:
                    continue
                if magnitude2.mag - magnitude1.mag < 1:
                    continue
                bx.append(magnitude1.mag)
                by.append(magnitude2.mag)
                e.append(event)
        print(f"{len(e)} blue events")
        if args.output:
            catalog.events = e
            catalog.write(args.output, format=args.input_format)
        ax1.plot(bx, by, 'b+', markersize=marker_size + 2)
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
    plt.show(block=True)


if __name__ == '__main__':
    main()
