__author__ = 'befulton'
import os
import datetime

def to_date(s):
    return datetime.datetime.strptime(s[0] + ' ' + s[1], "%Y%m%d %H:%M:%S")

files = [file for file in os.listdir(".") if file.endswith(".sum")]

files = [file for file in os.listdir(".") if file.endswith(".sum")]

interval = 5
runtimes = dict()
for file in files:
    app_name = file.split('.')[2]
    linecount = 0
    with open(file) as f:
        line = f.readline().split()
        if line:
            linecount = 1
            app_start = to_date(line)
            for line in f:
                linecount += 1
            if linecount > 1:
                line = line.split()
            app_end = to_date(line)
    runtimes[app_name] = linecount * interval

with open('runtime.csv', 'w') as f:
    f.write("application runtime[s]\n")
    for k, v in runtimes.iteritems():
        f.write("%s %s\n" % (k, v))

