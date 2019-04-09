#!/usr/bin/env python3

"""
Graph dssp .xpm map generated with "gmx do_dssp -o ss.xpm"
"""

from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import itertools

matplotlib.rcParams.update({'font.size':9})
plt.figure(figsize=(6, 6))


filename = "ss_ma10.xpm"


# The parsing should be rewritten with a parser.
lines = open(filename).readlines()
lines = [l.strip() for l in lines]  # strip the new lines characters

# parse x-axis labels
x_timepoints_str = [l.split()[2:-1] for l in lines if l.startswith("/* x-axis")]
x_timepoints = list(map(float, list(itertools.chain(*x_timepoints_str))))
print("Together timepoints:", len(x_timepoints))

# extract y-axis labels
y_residues_labels = [l.split()[2:-1] for l in lines if l.startswith("/* y-axis")]
y_residues = list(map(int, list(itertools.chain(*y_residues_labels))))
print("Together residues: ", len(y_residues))

# lines starting with " without the legend (lines that contain '/*')
data = [l for l in lines if l.startswith("\"") and '/*' not in l]
x_timepoints_number, y_residues_number, _, _ = [int(num) for num in data[0].strip(",").strip("\"").split()]
assert x_timepoints_number == len(x_timepoints) and y_residues_number == len(y_residues)

# remove the header
dssp = [l.rstrip(',').strip("\"") for l in data[1:]]
# do we have data for each residue?
assert len(dssp) == y_residues_number
# do we have data for each timepoint?
assert len(dssp[0]) == x_timepoints_number, \
    'the length of the first time point is %d but should be %d' % (len(dssp[0]), x_timepoints_number)


# construct letter to color mapping
coil_char = "~"
bsheet_char = "E"
bbridge_char = "B"
bend_char = "S"
turn_char = "T"
a_helix = 'H'
three_helix_char = "G"

# list the used colors in the growing order
sorted_colors = OrderedDict(sorted([
    (coil_char, "#FFFFFF"),
    (bsheet_char, "#FF0000"),
    (bbridge_char, "#000000"),
    (bend_char, "#008000"),
    (turn_char, "#FFFF00"),
    (a_helix, "#0000FF"),
    (three_helix_char, "#808080")
    ],
    key=lambda x:int(x[1].split('#')[1], 16)
    )
)

# asc_colors = OrderedDict(sorted(colors.items(), )
mycmap = ListedColormap(sorted_colors.values())

# the ss values have to be spaced uniformly for the colors to be applied correctly
colors_uniform = {key:i for i, key in enumerate(sorted_colors.keys())}

dssp_uniform = list(map(lambda x:[colors_uniform[letter] for letter in x], dssp))
# turn upside down, place the resid 1 at the bottom
dssp_uniform = dssp_uniform[::-1]

# graph using pcolormesh
plt.pcolormesh(dssp_uniform, cmap=mycmap)

# construct the patches for the legend
coil = mpatches.Patch(color=sorted_colors[coil_char], alpha=1, linewidth=0)
bsheet = mpatches.Patch(color=sorted_colors[bsheet_char], alpha=1, linewidth=0)
bbridge = mpatches.Patch(color=sorted_colors[bbridge_char], alpha=1, linewidth=0)
bend = mpatches.Patch(color=sorted_colors[bend_char], alpha=1, linewidth=0)
turn = mpatches.Patch(color=sorted_colors[turn_char], alpha=1, linewidth=0)
helix3 = mpatches.Patch(color=sorted_colors[three_helix_char], alpha=1, linewidth=0)

# graph the legend
plt.legend((coil, bsheet, bbridge, bend, turn, helix3), ('Coil', 'B-Sheet', 'B-Bridge', 'Bend', 'Turn', '3-Helix'),
           bbox_to_anchor=(0.5, 1.10),  # fixed o the legend (x, y)
           loc='upper center',
           ncol=6)

plt.xlabel('Time (ns)')
plt.ylabel('Residue')
plt.tight_layout()
plt.savefig('dssp_map_single.png', dpi=300)
plt.show()