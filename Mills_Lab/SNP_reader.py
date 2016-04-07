#!/usr/bin/env python
from __future__ import print_function, division
import matplotlib as mpl
import matplotlib.pyplot as plt
from shutil import rmtree as rm
# import matplotlib.cm as cm
import glob
import pickle
import math
import numpy as np
import os
import sys
import csv

"""If getting unable to connect, exit status -1:
go to Run config and set DISPLAY to either 14 or 10 if running debug mode.
Furthermore, you need to be SSH in MobaXterm at the same time. Otherwise, this
script works from from the command line.
"""


def file_globber(pathname):
    snp_files = []
    for dirs, _, files in os.walk(pathname):
        if "SNPs" in dirs:
            snp_files.extend(glob.glob(os.path.join(dirs, "*")))
    return snp_files


class CommentedFile:

    def __init__(self, f, commentstring="#"):
        self.f = f
        self.commentstring = commentstring
        self.length = None
        self.item = 0

    def next(self):
        line = self.f.next()
        if line.startswith("##ID"):
            line = self.f.next()
        if line.startswith("##"):
            self.length = int(line.split()[4])-int(line.split()[3])
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
        global qLook
        global sampleGroup
        global SNPs
        idx = qLook.get(string)
        if idx is not None:
            hetrna = np.asarray([float(row[2]) for row in sampleGroup[idx][2]
                                 if ('1|0' or '0|1') in row[1]])
            hetribo = np.asarray([float(row[3]) for row in sampleGroup[idx][2]
                                  if ('1|0' or '0|1') in row[1]])
            homorna = np.asarray([float(row[2]) for row in sampleGroup[idx][2]
                                  if '1|1' in row[1]])
            homoribo = np.asarray([float(row[3]) for row in sampleGroup[idx][2]
                                   if '1|1' in row[1]])
            refrna = np.asarray([float(row[2]) for row in sampleGroup[idx][2]
                                 if '0|0' in row[1]])
            refribo = np.asarray([float(row[3]) for row in sampleGroup[idx][2]
                                  if '0|0' in row[1]])
            rna = [refrna, hetrna, homorna]
            ribo = [refribo, hetribo, homoribo]
            ticks = ["0|0 (n=%d)" % len(refrna),
                     "0|1 (n=%d)" % len(hetrna),
                     "1|1 (n=%d)" % len(homorna)]
            xlab = "Genotypes"
            ylab = "FPKM"
            title = string
            figrna = plt.figure()
            axrna = figrna.add_subplot(111)
            axrna.boxplot(rna, labels=ticks)
            axrna.set_title(title + ": RNA-seq (N=%d)"
                            % sum((len(hetrna), len(homorna), len(refrna))))
            axrna.set_ylabel(ylab)
            axrna.set_xlabel(xlab)
            figrna.show()
            figribo = plt.figure()
            axribo = figribo.add_subplot(111)
            axribo.boxplot(ribo, labels=ticks)
            axribo.set_title(title + ": Ribosome Profiling (N=%d)"
                            % sum((len(hetrna), len(homorna), len(refrna))))
            axribo.set_ylabel(ylab)
            axribo.set_xlabel(xlab)
            figribo.show()

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
            t = ax.text(x, y, " - %s" % annote)
            m = ax.scatter([x], [y], marker='d', c='r', zorder=100)
            self.drawnAnnotations[(x, y)] = (t, m)
            self.sub_boxplot(annote)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawAnnote(self.ax, x, y, a)


def meta_list(LIST):
    try:
        output = (sum(1 for rna in LIST if rna[0] > 0.0)/len(LIST),
                  sum(1 for ribo in LIST if ribo[1] > 0)/len(LIST))
        return output
    except ZeroDivisionError:
        return 0


def main():
    global qLook
    global sampleGroup
    global SNPs

    try:
        if not "-n" in sys.argv[1] :
            path_name = sys.argv[1]
        elif "-n" in sys.argv[1]:
            try:
                rm(sys.argv[2] + '/pkl')
            except OSError:
                path_name = sys.argv[2]
    except IndexError:
        path_name = '/home/mdsherm/Project/SNuPer_results/chr22'
    if not os.path.isdir(path_name + '/pkl'):
        os.mkdir(path_name + '/pkl')
    if not os.path.isfile(path_name + '/pkl/plotzip.pkl'):
        snp_files = file_globber(path_name)
        sampleGroup = [(snp.rpartition("SNPs/")[-1], [row for row in csv.reader(
            open(snp, "rb"), delimiter='\t')]) for snp in snp_files]
        sampleGroup = [(snp[0], int(snp[1][1][4])-int(snp[1][1][3]),
                        snp[1][3:]) for snp in sampleGroup]
        genos = []
        for snp in sampleGroup:
            ID = snp[0]
            Lengths = snp[1]
            ref = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2] if "0|0" in sample[1]])
            het = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2] if ("1|0" or "0|1") in sample[1]])
            alt = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2] if "1|1" in sample[1]])
            raw = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2]])
            try:
                if np.mean(ref, axis=0)[1] < np.mean(
                        het, axis=0)[1] < np.mean(alt, axis=0)[1]:
                    genos.append([ID, Lengths, ref, het, alt, raw])
                else:
                    pass
            except IndexError:
                pass

        SNPs = genos[:]
        qLook = {entry[0]: i for (i, entry) in enumerate(sampleGroup)}
        pickle.dump(genos,
                    open(path_name + "/pkl/genos.pkl", "wb"))
        pickle.dump(sampleGroup,
                    open(path_name + "/pkl/SG.pkl", "wb"))
        pickle.dump(qLook,
                    open(path_name + "/pkl/qLook.pkl", "wb"))
        SNP_IDs = [snp[0] for snp in SNPs]
        SNP_len = [snp[1] for snp in SNPs]
        SNP_ratio_step = [np.mean(np.array(snp[4]), axis=0)[1] /
                          np.mean(np.array(snp[2]), axis=0)[1] for snp in SNPs]
        SNP_ratio = [np.log2(step) if step > 0 else 0 for step in SNP_ratio_step]
        # SNP_ratio = [math.log(np.mean(SNPs[4], axis=0), 2)
        #              /math.log(np.mean(SNPs[2], axis=0), 2)]
        percents = []
        for snp in range(len(SNPs)):
            step = [(sample[0], sample[1]) for sample in SNPs[snp][5]]
            percents.append(
                (float(sum(1 for rna in step if rna[0] > 0.0)/float(len(step))),
                 float(sum(1 for ribo in step if ribo[1] > 0.0)/float(len(step)))))
        pkl = zip(SNP_IDs, SNP_len, SNP_ratio, percents)
        pickle.dump(pkl,
                    open(path_name + "/pkl/plotzip.pkl", "wb"))
    else:
        pkl = pickle.load(
            open(path_name + "/pkl/plotzip.pkl", "rb"))
        sampleGroup = pickle.load(
            open(path_name + "/pkl/SG.pkl", "rb"))
        print(len(sampleGroup))
        qLook = pickle.load(
            open(path_name + "/pkl/qLook.pkl", "rb"))
        SNP_IDs, SNP_len, SNP_ratio, percents = [], [], [], []
        for row in pkl:
            SNP_IDs.append(row[0])
            SNP_len.append(row[1])
            SNP_ratio.append(row[2])
            percents.append(row[3])
    annotes = SNP_IDs
    sizes = (SNP_len / np.mean(SNP_len)) * 10
    colors = SNP_ratio
    x = [snp[0] for snp in percents]
    y = [snp[1] for snp in percents]

    # Plots the points above, and can be used to tie in individual SNP IDs
    fig, ax = plt.subplots()
    ax.scatter(x, y, color=colors, cmap=plt.get_cmap('YlOrRd'),
               s=sizes, linewidths=0.2, edgecolors='black', alpha=0.8)
    ax.set_title("Chr22")
    ax.set_xlabel('%RNA-seq > 0.0')
    ax.set_ylabel('%Ribosome Profiling > 0.0')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3", alpha=0.35)
    af = AnnoteFinder(x, y, annotes, ax=ax)
    fig.canvas.mpl_connect('button_press_event', af)
    plt.show()


if __name__ == "__main__":
    main()
