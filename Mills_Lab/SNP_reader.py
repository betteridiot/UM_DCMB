#!/usr/bin/env python
from __future__ import print_function, division
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import math
import numpy as np
import os
import sys
import csv

# TODO if getting unable to connect, exit status -1: go to Run config and set DISPLAY to either 14 or 10
# path_name = '/home/mdsherm/Project/SNuPer_results/pythonTest/100ktest/'
# snp_files = glob.glob(path_name+"*.snp")

# TODO this is for globbing all SNPs through all subparts


def file_globber():
    try:
        path_name = sys.argv[1]
    except IndexError:
        path_name = '/home/mdsherm/Project/SNuPer_results/chr22'
    snp_files = []
    for dirs, _, files in os.walk(path_name):
        if "SNPs" in dirs:
            snp_files.extend(glob.glob(os.path.join(dirs, "*")))
    return path_name, snp_files


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
    def __init__(self, xdata, ydata, annotes, ax=None, xtol=5, ytol=5):
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
        return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

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

    def sub_boxplot(self, string):
        idx = qLook.get(string)
        if idx is not None:
            hetrna = np.asarray([float(row[2]) for row in sampleGroup[idx][1] if '1|0" or "0|1' in row[1]])
            hetribo = np.asarray([float(row[3]) for row in sampleGroup[idx][1] if '1|0" or "0|1' in row[1]])
            homorna = np.asarray([float(row[2]) for row in sampleGroup[idx][1] if '1|1' in row[1]])
            homoribo = np.asarray([float(row[3]) for row in sampleGroup[idx][1] if '1|1' in row[1]])
            refrna = np.asarray([float(row[2]) for row in sampleGroup[idx][1] if '0|0' in row[1]])
            refribo = np.asarray([float(row[3]) for row in sampleGroup[idx][1] if '0|0' in row[1]])
            rna = [refrna, hetrna, homorna]
            ribo = [refribo, hetribo, homoribo]
            ticks = ["0|0", "0|1", "1|1"]
            xlab = "Genotypes"
            ylab = "FPKM"
            title = string
            figrna = plt.figure()
            axrna = figrna.add_subplot(111)
            axrna.boxplot(rna, labels=ticks)
            axrna.set_title(title + "RNA-seq")
            axrna.set_ylabel(ylab)
            axrna.set_xlabel(xlab)
            figrna.show()
            figribo = plt.figure()
            axribo = figribo.add_subplot(111)
            axribo.boxplot(ribo, labels=ticks)
            axribo.set_title(title + "Ribosome Profiling")
            axribo.set_ylabel(ylab)
            axribo.set_xlabel(xlab)
            figribo.show()

    def drawAnnote(self, ax, x, y, annote, event):
        """
        Draw the annotation on the plot
        """
        global qLook
        global sampleGroup
        if (x, y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x, y)]
            for m in markers:
                m.set_visible(not m.get_visible())
            self.ax.figure.canvas.draw_idle()
        else:
            t = ax.text(x, y, " - %s" % annote,)
            m = ax.scatter([x], [y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x, y)] = (t, m)
            if event.dblclick:
                self.sub_boxplot(annote)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawAnnote(self.ax, x, y, a)


def meta_list(LIST):
    try:
        output = (sum(1 for rna in LIST if rna[0] > 0.0)/len(LIST), sum(1 for ribo in LIST if ribo[1] > 5)/len(LIST))
        return output
    except ZeroDivisionError:
        return 0


def main():
    global qLook
    global sampleGroup
    path_name, snp_files = file_globber()
    sampleGroup = [(snp.rpartition("SNPs/")[-1], [row for row in csv.reader(
        CommentedFile(open(snp, "rb")), delimiter='\t')]) for snp in snp_files]
    qLook = {entry[0]:i for (i, entry) in enumerate(sampleGroup)}
    SNP_IDs = [snp[0] for snp in sampleGroup]
    percents = []
    for snp in sampleGroup:
        step = [(float(samples[2]), float(samples[3])) for samples in snp[1]]
        percents.append((sum(1 for rna in step if rna[0] > 0.0)/len(step),
                         sum(1 for ribo in step if ribo[1] > 0.0)/len(step)))
    SNP_list = zip(SNP_IDs, percents)
    # Gives me all SNPs that have %RNA-seq >.8 and %Ribo >.5
    axis = [(snp[1][0], snp[1][1]) for snp in SNP_list]
    annotes = [snp[0] for snp in SNP_list]
    x = [x[0] for x in axis]
    y = [y[1] for y in axis]
    t = [(vector[0]*vector[1]) for vector in axis]

    # Plots the points above, and can be used to tie in individual SNP IDs
    fig, ax = plt.subplots()
    ax.scatter(x, y, color=t, cmap=cm.YlOrRd, s=10, linewidths=0.1, edgecolors='black', alpha=0.5)
    ax.set_title("Chr22")
    ax.set_xlabel('%RNA-seq > 0.0')
    ax.set_ylabel('%Ribosome Profiling > 0.0')
    af = AnnoteFinder(x, y, annotes, ax=ax)
    fig.canvas.mpl_connect('button_press_event', af)
    plt.show()


if __name__ == "__main__":
    main()
