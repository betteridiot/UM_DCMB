#Psuedocode
#glob all of the SNP files
#iteratively open and read in as tab delimited
#skip lines with pound signs
#find the ones that have both rna-seq and ribosome profiling data

from __future__ import print_function, division
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import math
import random
import string
import csv

path_name = '/home/mdsherm/Project/SNuPer_results/pythonTest/100ktest/'
snp_files = glob.glob(path_name+"*.snp")


class CommentedFile:
    def __init__(self, f, commentstring="#"):
        self.f = f
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring):
            line = self.f.next()
        return line
    def __iter__(self):
        return self


class AnnoteFinder(object):
    """callback for matplotlib to display an annotation when points are
    clicked on.  The point which is closest to the click and within
    xtol and ytol is identified.

    Register this function like this:

    scatter(xdata, ydata)
    af = AnnoteFinder(xdata, ydata, annotes)
    connect('button_press_event', af)
    """

    def __init__(self, xdata, ydata, annotes, ax=None, xtol=None, ytol=None):
        self.data = list(zip(xdata, ydata, annotes))
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.drawnAnnotations = {}
        self.links = []

    def distance(self, x1, x2, y1, y2):
        """
        return the distance between two points
        """
        return(math.sqrt((x1 - x2)**2 + (y1 - y2)**2))

    def __call__(self, event):

        if event.inaxes:

            clickX = event.xdata
            clickY = event.ydata
            if (self.ax is None) or (self.ax is event.inaxes):
                annotes = []
                # print(event.xdata, event.ydata)
                for x, y, a in self.data:
                    # print(x, y, a)
                    if ((clickX-self.xtol < x < clickX+self.xtol) and
                            (clickY-self.ytol < y < clickY+self.ytol)):
                        annotes.append(
                            (self.distance(x, clickX, y, clickY), x, y, a))
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0]
                    self.drawAnnote(event.inaxes, x, y, annote)
                    for l in self.links:
                        l.drawSpecificAnnote(annote)

    def drawAnnote(self, ax, x, y, annote):
        """
        Draw the annotation on the plot
        """
        if (x, y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x, y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.ax.figure.canvas.draw_idle()
        else:
            t = ax.text(x, y, " - %s" % (annote),)
            m = ax.scatter([x], [y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x, y)] = (t, m)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawAnnote(self.ax, x, y, a)

# Make a unique filename for each SNP written
def id_generator(size=9, chars=string.ascii_uppercase + string.digits):
    """Creates a random alphanumeric id each time it is called.

    Args:
        size (int): what is the desired size of the ID string
        chars (str): What values the function can choose from

    Returns:
        A pseudo-random alphanumeric ID
    """
    return ''.join(random.choice(chars) for _ in range(size))

sampleGroup = [[row for row in csv.reader(CommentedFile(open(file, "rb")), delimiter='\t')] for file in snp_files]
lister = []
for snp in sampleGroup:
    stepx = [(float(sample[2]), float(sample[3])) for sample in snp]
    lister.append((sum(1 for i in stepx if i[0] > 0.0)/len(stepx), sum(1 for j in stepx if j[1] > 0.0)/len(stepx)))
axis = [(snp[0], snp[1]) for snp in lister if snp[0] > .8 and snp[1] > .5]
x = [x[0] for x in axis]
y = [y[1] for y in axis]

annotes = [id_generator() for _ in range(len(sampleGroup))]
fig, ax = plt.subplots()
ax.scatter(x,y, color='orange', linewidths=0.1, edgecolors='black')
ax.set_title("Chr22 100k Test data")
ax.set_xlabel('%RNA-seq > 0.0')
ax.set_ylabel('%Ribosome Profiling > 0.0')
af =  AnnoteFinder(x,y, annotes, ax=ax)
fig.canvas.mpl_connect('button_press_event', af)
plt.show()

