__author__ = 'befulton'

import os
import sys
import subprocess

def get_times():
    d = dict()
    with open("global.time") as f:
        for line in f:
            s = line.split()
            d[s[0]] = s[1:]
    return d

times = get_times()
date = times['start'][0]
start = times['start'][1]
end = times['end'][1]
tics = int(end.split(':')[0]) + 1

prettycolors = False
if len(sys.argv) > 1:
    name = sys.argv[1]
else:
    name = os.path.split(os.getcwd())[1]

if len(sys.argv) > 2:
    cpu_count = int(sys.argv[2])
else:
    cpu_count = 64


def build_plot(files, stat):
    d = []
    for i, file in enumerate(sorted(files, key=lambda fn: int(fn.split('.')[0]))):
        title = file.split('.')[2]
        d.append("'%s' using 1:(%s) title \"%s\"  ls %s" % (file, stat, title, i+1))
    return d


def write_files(name):
    files = [file for file in os.listdir(".") if file.endswith(".sum")]

    with open("defs.gnu", 'w') as f:
        f.write("#generated gnuplot file\n")
        f.write("set xrange ['%s %s':'%s %s']\n" % (date, start, date, end))
        f.write("set xtics (")
        f.write(", ".join(["\"{0}\" '{1} {0:02d}:00:00'".format(i, date) for i in range(tics)]))
        f.write(" )\n")
        f.write("ncpu=%s" % cpu_count)
        f.write("\n")

    with open("common.gnu", 'w') as f:
        f.write(common_text.format(name))
        if prettycolors:
            colorset = color_sets[len(files)-1]
            for i, c in enumerate(colorset.split()):
                f.write("set style line %s lt 1 lc rgb \"%s\" lw mylw1 pt mypt1 ps myps1\n" % i, c)
        else:
            for i in range(len(files)):
                f.write("set style line %s lw mylw1 pt mypt1 ps myps1\n" % (i+1))

    with open('ram.gnu', 'w') as f:
        f.write(gnu_head % ('ram', end, "RAM usage GiB", ""))
        d = build_plot(files, "fg*$11")
        f.write(",\\\n     ".join(d))

    with open('cpu.gnu', 'w') as f:
        f.write(gnu_head % ('cpu', end, "Core Utilization", "ncpu"))
        d = build_plot(files, "$19/100")
        f.write(",\\\n     ".join(d))

    with open('io.gnu', 'w') as f:
        f.write(gnu_io_head % end)
        d = build_plot(files, "fm*($23+$24)")
        f.write(",\\\n     ".join(d))

    for gnu in ['ram.gnu', 'cpu.gnu', 'io.gnu']:
        subprocess.Popen(['gnuplot', gnu]).wait()

common_text = """
set terminal postscript color eps "Times" 14
#set termoption enhanced

set style data points

# ram and cpu
mypt1=7
mylw1=0
#myps1=0.3
myps1=0.6

set key below
set timefmt "%Y%m%d %H:%M:%S"
set xdata time
set format x "%k"
set xrange [*:*]
unset grid
set grid x
unset title
#set title "{0}" noenhanced
set title "{0}"

fm=1./(1024.0)
fg=1./(1024.0)/(1024.0)

#set multiplot layout 3,1
#set tmargin 0
#set tmargin 0.8
#set bmargin 0.8

#set tics nomirror scale 0.66
#set xtics scale 0

#myvertoffset=0.02

"""

gnu_head = """load 'common.gnu'
load 'defs.gnu'

set out '%s.eps'

set xlabel "Runtime [h] %s"
set ylabel "%s"
set yrange [0:%s]

plot """

gnu_io_head = """load 'common.gnu'
load 'defs.gnu'

set out 'io.eps'

set xlabel "Runtime [h] %s"
set ylabel "I/O MiB/s"
set yrange [0.005:2000]
set logscale y
set ytics ("10^{-3}" 0.001, "10^{-2}" 0.01, "10^{-1}" 0.1, "10^{0}" 1, "10^{1}" 10, "10^{2}" 100, "10^{3}" 1000)

plot """

color_sets = [
    '#f800ff',
    '#ff00e3 #f100ff',
    "#ffa900 #00ffdc #0096ff",
    "#ff005d #19ff00 #00ffc8 #c200ff",
    "#ff00c7 #9cff00 #00ff61 #0042ff #5700ff",
    "#ff00cb #ff1500 #59ff00 #00ffb5 #00b4ff #4f00ff",
    "#ff00f2 #ff000f #c0ff00 #0bff00 #00ffdd #0007ff #7800ff",
    "#ff004d #ff0100 #ffe800 #64ff00 #00ff94 #00daff #0070ff #c800ff",
    "#ff0081 #ff0031 #ffa500 #6aff00 #23ff00 #00ffb5 #0058ff #1300ff #ab00ff",
    "#ff00af #ff003b #ff4b00 #fff900 #6aff00 #00ff98 #00feff #00b6ff #3c00ff #8700ff",
    "#ff00de #ff1100 #ff7d00 #fff600 #59ff00 #1bff00 #00ff64 #00e7ff #009cff #7100ff #bc00ff",
    "#ff00b7 #ff0062 #ff5200 #ffdb00 #baff00 #54ff00 #00ff5e #00ff87 #0087ff #0083ff #7800ff #cc00ff",
    "#ff00b9 #ff0037 #ff2600 #ffa500 #bbff00 #84ff00 #00ff03 #00ff64 #00ffb5 #00cfff #0048ff #2b00ff #ba00ff",
    '#ff00d8 #ff0047 #ff0007 #ff7f00 #ffee00 #88ff00 #32ff00 #00ff67 #00ff75 #00c9ff #0070ff #0005ff #6c00ff #8f00ff',
    "#ff00db #ff0066 #ff2d00 #ff3400 #fff600 #adff00 #7eff00 #00ff03 #00ff72 #00ffc5 #00c6ff #0052ff #002fff #9500ff #bb00ff",
    "#ff00c8 #ff009b #ff0004 #ff7600 #ffbe00 #ffef00 #98ff00 #54ff00 #00ff2c #00ff6a #00ffbf #0088ff #004fff #0010ff #8900ff #b300ff",
    "#ff00ed #ff0053 #ff0005 #ff4200 #ffc200 #f9ff00 #aaff00 #52ff00 #14ff00 #00ff74 #00ffb5 #00ffe6 #00b4ff #0025ff #0003ff #6300ff #be00ff",
    "#ff00a9 #ff006c #ff004a #ff4000 #ff5a00 #ffae00 #e8ff00 #9bff00 #36ff00 #00ff3a #00ff7a #00fff4 #00bdff #007cff #0017ff #2a00ff #5d00ff #be00ff",
    "#ff00cb #ff007f #ff0054 #ff0b00 #ff7d00 #ffa700 #ffe400 #aaff00 #7cff00 #18ff00 #00ff41 #00ffb1 #00ffe1 #009cff #0071ff #001eff #1900ff #8a00ff #b400ff",
    "#ff00e1 #ff0086 #ff0037 #ff0400 #ff6f00 #ff8400 #fff100 #aeff00 #71ff00 #05ff00 #00ff22 #00ff56 #00ffd2 #00ecff #00c4ff #003eff #001cff #3300ff #8000ff #c700ff",
    "#ff00f0 #ff009c #ff004e #ff0004 #ff3300 #ff7400 #ffcc00 #f4ff00 #9fff00 #53ff00 #13ff00 #00ff3a #00ff9c #00ffb5 #00f1ff #00b4ff #0078ff #0007ff #3700ff #8a00ff #eb00ff",
    "#ff00d4 #ff00a1 #ff0064 #ff001f #ff3100 #ff9d00 #ffa200 #ffeb00 #a5ff00 #80ff00 #43ff00 #02ff00 #00ff49 #00ff88 #00ffe4 #00e6ff #008fff #004fff #001dff #4300ff #7900ff #b000ff",
    "#ff00f6 #ff00af #ff0047 #ff0012 #ff2e00 #ff6300 #ffc700 #fffc00 #bcff00 #6bff00 #63ff00 #0fff00 #00ff29 #00ff8b #00ffc5 #00fff1 #00c5ff #0072ff #0053ff #0200ff #3600ff #7700ff #bd00ff",
    "#ff00d7 #ff009e #ff0061 #ff0037 #ff3000 #ff7200 #ffaa00 #ffec00 #cfff00 #94ff00 #4eff00 #10ff00 #08ff00 #00ff59 #00ff86 #00ffc1 #00dcff #00b5ff #006cff #0041ff #1d00ff #5100ff #9800ff #d100ff",
    "#ff00e4 #ff0094 #ff007e #ff003b #ff0000 #ff4e00 #ff9100 #ffb000 #ffec00 #ccff00 #76ff00 #28ff00 #00ff17 #00ff50 #00ff5c #00ffa9 #00fff1 #00e6ff #008bff #005bff #0005ff #1b00ff #7000ff #8a00ff #c500ff",
    "#ff00ed #ff00a9 #ff0059 #ff0034 #ff0009 #ff4d00 #ff6200 #ffab00 #fffe00 #e3ff00 #8bff00 #7aff00 #41ff00 #00ff06 #00ff45 #00ff7c #00ffb5 #00fff7 #00b2ff #00a4ff #0067ff #0033ff #2e00ff #5a00ff #9d00ff #d500ff",
    "#ff00e2 #ff00af #ff0059 #ff0033 #ff1000 #ff2c00 #ff7500 #ffb100 #ffe600 #faff00 #c9ff00 #94ff00 #50ff00 #17ff00 #00ff39 #00ff74 #00ffad #00ffe6 #00fff9 #00c8ff #0089ff #004cff #0019ff #1700ff #6e00ff #9400ff #d600ff",
    "#ff00cd #ff00af #ff008a #ff003a #ff0020 #ff2f00 #ff6c00 #ff8300 #ffc200 #ffea00 #b5ff00 #82ff00 #71ff00 #12ff00 #00ff16 #00ff58 #00ff69 #00ffb6 #00ffd5 #00d0ff #00bcff #0077ff #002cff #0020ff #3c00ff #5900ff #9e00ff #be00ff",
    "#ff00d6 #ff00b4 #ff007e #ff005d #ff0013 #ff0f00 #ff4400 #ff8a00 #ffb300 #fff100 #e1ff00 #beff00 #6eff00 #31ff00 #03ff00 #00ff28 #00ff68 #00ff83 #00ffad #00ffe7 #00deff #0099ff #0053ff #0038ff #0011ff #4100ff #6d00ff #8800ff #e800ff",
    "#ff00d8 #ff00ae #ff0078 #ff0036 #ff0003 #ff0400 #ff5500 #ff8600 #ffb100 #ffd100 #cfff00 #b4ff00 #7dff00 #65ff00 #1cff00 #00ff0c #00ff4b #00ff86 #00ffb9 #00ffe5 #00f5ff #00a2ff #0086ff #0064ff #0014ff #2d00ff #4a00ff #8200ff #b100ff #de00ff",
]
if __name__ == '__main__':
    write_files(name)