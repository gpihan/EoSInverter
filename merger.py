import itertools
from statistics import median
import sys

ACCURACY = 1e-6

class DataRow:
    def __init__(self, e, nB, T, muB, P, S):
        self.e = e
        self.nB = nB
        self.T = T
        self.muB = muB
        self.P = P
        self.S = S

    def __lt__(self, other):
        if abs(self.e - other.e) > ACCURACY: return self.e < other.e
        if abs(self.nB - other.nB) > ACCURACY: return self.nB < other.nB
        if abs(self.T - other.T) > ACCURACY: return self.T < other.T
        if abs(self.muB - other.muB) > ACCURACY: return self.muB < other.muB
        if abs(self.P - other.P) > ACCURACY: return self.P < other.P
        return self.S < other.S

    def __eq__(self, other):
        return (abs(self.e - other.e) < ACCURACY and
                abs(self.nB - other.nB) < ACCURACY and
                abs(self.T - other.T) < ACCURACY and
                abs(self.muB - other.muB) < ACCURACY and
                abs(self.P - other.P) < ACCURACY and
                abs(self.S - other.S) < ACCURACY)

    def as_tuple(self):
        return (self.e, self.nB, self.T, self.muB, self.P, self.S)


def read_data_from_file(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            try:
                e, nB, T, muB, P, S = map(float, line.split())
                data.append(DataRow(e, nB, T, muB, P, S))
            except Exception:
                print(f"Line ignored due to format in {filename}")
    print(f"{len(data)} lines read from {filename}")
    return data


def write_data_to_file(filename, data):
    with open(filename, 'w') as f:
        for row in data:
            f.write("  ".join(f"{x:.17g}" for x in row.as_tuple()) + "\n")
    print(f"{len(data)} lines written to {filename}")


def relabs(x):
    return abs(x)


def median_from_unique_diffs(values, tol):
    if len(values) < 2:
        return 0.0
    values = sorted(values)
    uniq = []
    for v in values:
        if not uniq or abs(v - uniq[-1]) > tol:
            uniq.append(v)
    if len(uniq) < 2:
        return 0.0
    diffs = [abs(uniq[i+1] - uniq[i]) for i in range(len(uniq)-1) if abs(uniq[i+1]-uniq[i]) > 0]
    if not diffs:
        return 0.0
    return median(diffs)


def auto_tolerances(data):
    E = [r.e for r in data]
    NB = [r.nB for r in data]
    med_de = median_from_unique_diffs(E, ACCURACY)
    med_dnB = median_from_unique_diffs(NB, ACCURACY)
    if not (med_de > 0):
        med_de = (max(E)-min(E))*1e-3
    if not (med_dnB > 0):
        med_dnB = (max(NB)-min(NB))*1e-3
    return {
        'group_tol_e': 0.5*med_de,
        'group_tol_nB': 0.5*med_dnB,
        'step_min_e': 0.10*med_de,
        'step_min_nB': 0.10*med_dnB
    }


def mark_flat_along(data, coord, group_tol, step_min, val_abs, val_rel):
    N = len(data)
    idx = list(range(N))
    if coord == 0:
        idx.sort(key=lambda i: (data[i].nB, data[i].e))
    else:
        idx.sort(key=lambda i: (data[i].e, data[i].nB))
    flat = [False]*N
    start = 0
    while start < N:
        ref = data[idx[start]].nB if coord==0 else data[idx[start]].e
        end = start
        while end+1 < N:
            vnext = data[idx[end+1]].nB if coord==0 else data[idx[end+1]].e
            if relabs(vnext-ref) <= group_tol:
                end += 1
            else:
                break
        g = idx[start:end+1]
        g.sort(key=lambda i: data[i].e if coord==0 else data[i].nB)
        for k in range(1,len(g)-1):
            ip, ic, inext = g[k-1], g[k], g[k+1]
            step = (data[inext].e-data[ip].e) if coord==0 else (data[inext].nB-data[ip].nB)
            if relabs(step) < step_min:
                continue
            def small_change(a,b):
                thr = val_abs + val_rel*max(relabs(a), relabs(b))
                return relabs(a-b) <= thr
            flat_here = (small_change(data[inext].T, data[ip].T) and
                         small_change(data[inext].muB, data[ip].muB) and
                         small_change(data[inext].P, data[ip].P) and
                         small_change(data[inext].S, data[ip].S))
            if flat_here:
                flat[ic] = True
        start = end+1
    return flat


def filter_const_regions(data, group_tol_e, group_tol_nB, step_min_e, step_min_nB, val_abs=1e-9, val_rel=1e-7, require_both=False, debug=True):
    N = len(data)
    if N<3:
        return data
    flat_in_e = mark_flat_along(data,0,group_tol_nB,step_min_e,val_abs,val_rel)
    flat_in_nB = mark_flat_along(data,1,group_tol_e,step_min_nB,val_abs,val_rel)
    out = []
    for i,row in enumerate(data):
        kill = (flat_in_e[i] and flat_in_nB[i]) if require_both else (flat_in_e[i] or flat_in_nB[i])
        if not kill:
            out.append(row)
        elif debug:
            print(f"Deleting: e={row.e} nB={row.nB} T={row.T} muB={row.muB} P={row.P} S={row.S}")
    return out


def main():
    if len(sys.argv) < 3:
        print("Usage: % python merger.py EoS_Cross.dat EoS_below.dat EoS_above.dat EoS_combined.dat")
        sys.exit(1)

    inputfiles = sys.argv[1:-1]
    outputFile = sys.argv[-1]

    combinedData = []
    for fname in inputfiles:
        combinedData.extend(read_data_from_file(fname))

    print("\nTotal lines before deleting duplicates", len(combinedData))

    combinedData.sort()
    # remove duplicates
    uniqueData = []
    for k, g in itertools.groupby(combinedData):
        uniqueData.append(next(g))

    print("Total lines after deleting duplicates", len(uniqueData))

    steps = auto_tolerances(uniqueData)
    finalData = filter_const_regions(uniqueData, steps['group_tol_e'], steps['group_tol_nB'],
                                     steps['step_min_e'], steps['step_min_nB'],
                                     val_abs=1e-5, val_rel=1e-4, require_both=False, debug=True)

    print("Total lines after filtering:", len(finalData))

    write_data_to_file(outputFile, finalData)
    print("\nMerging successful!")


if __name__ == "__main__":
    main()
