#!/usr/bin/env python
from __future__ import print_function, division
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
import math
import numpy as np
import os
import csv
import fnmatch
import argparse
from operator import itemgetter
from mpl_toolkits.mplot3d import Axes3D

"""If getting unable to connect, exit status -1:
go to Run config and set DISPLAY to either 14 or 10 if running debug mode.
Furthermore, you need to be SSH in MobaXterm at the same time. Otherwise, this
script works from from the command line.
"""

def snp_set(lst, threshold):
    """Makes a set list of SNPs that % of RNA & ribosome FPKM
    is greater than 5. As it is a set, it only has unique entries
    for each SNP, and therefore allows for control of duplicate
    SNPs later.

    Args:
        lst (list): list of metadata files from given pathname
        threshold (float): threshold for metadata cutoff

    Returns:
        setter (set): a set of unique SNPs that meet thresholds
    """
    setter = set()
    for meta in lst:
        for row in csv.reader(open(meta, 'rb'), delimiter='\t'):
            if row[0].startswith("#"):
                continue
            elif all((row[5], row[6])) >= threshold:
                setter.add(row[0])
            else:
                continue
    print("{} unique SNPs".format(len(setter)))
    return setter


def meta_catcher(ROOT, pattern):
    """Collects the paths to all metadata files created from
    ORFSNuPer.py

    Args:
        ROOT (str): root directory to walk through
        pattern (str): how selection of files are filtered

    Returns:
        metas (list): list of metadata files and their paths
    """
    metas = []
    for root, sub, files in os.walk(ROOT):
        metadata = fnmatch.filter(files, pattern)
        metas.extend(os.path.join(root, f) for f in metadata)
    return metas


def file_globber(pathname, setList):
    """Walks through path to find files with prefixes present
    in the snp set list created by snp_set()

    Args:
        pathname (str): root directory to walk through
        setList (set): set of SNPs that meet basic criteria

    Returns:
        snp_files (list): a list of SNP files and their paths
    """
    snp_files = []
    for dirs, _, files in os.walk(pathname):
        if "SNPs" in dirs:
            snp_files.extend(['{}/{}'.format(dirs, f) for f in files if f.rpartition(':')[0] in setList])
    print('{} SNP files to be processed'.format(len(snp_files)))
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
            # hopefully skips over SNPs that are too small for Spectre to analyze
            if self.length < 30:
                pass
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
        # global sampleGroup
        global SNPs
        idx = qLook.get(string)
        if idx is not None:
            hetrna = np.asarray([sample[0] for sample in SNPs[idx][3]], dtype=np.float64)
            hetribo = np.asarray([sample[1] for sample in SNPs[idx][3]], dtype=np.float64)
            normhet = np.asarray(
                [sample[1] / np.mean(hetrna) for sample in SNPs[idx][3]])
            homorna = np.asarray([sample[0] for sample in SNPs[idx][4]], dtype=np.float64)
            homoribo = np.asarray([sample[1] for sample in SNPs[idx][4]], dtype=np.float64)
            normalt = np.asarray(
                [sample[1] / np.mean(homorna) for sample in SNPs[idx][4]])
            refrna = np.asarray([sample[0] for sample in SNPs[idx][2]], dtype=np.float64)
            refribo = np.asarray([sample[1] for sample in SNPs[idx][2]], dtype=np.float64)
            normref = np.asarray(
                [sample[1]/ np.mean(refrna) for sample in SNPs[idx][2]], dtype=np.float64)
            rna = [refrna, hetrna, homorna]
            rna_median = [np.median(rna[0]), np.median(rna[1]),
                          np.median(rna[2])]
            ribo = [refribo, hetribo, homoribo]
            ribo_median = [np.median(ribo[0]), np.median(ribo[1]),
                          np.median(ribo[2])]
            norm = [normref, normhet, normalt]
            norm_masker = [~np.isinf(norm[0]), ~np.isinf(norm[1]), ~np.isinf(norm[2])]
            masked_norm = [norm[0][norm_masker[0]], norm[1][norm_masker[1]], norm[2][norm_masker[2]]]
            norm_median = [np.median(masked_norm[0]), np.median(masked_norm[1]),
                          np.median(masked_norm[2])]

            ticks = ["0|0 (n=%d)" % len(refrna),
                     "0|1 (n=%d)" % len(hetrna),
                     "1|1 (n=%d)" % len(homorna)]
            xlab = "Genotypes"
            ylab = "FPKM"
            title = string + ' log2[alt/ref] = %f' % np.log2(np.mean(homoribo)/np.mean(refribo))
            figmix, (axrna, axribo, axnorm) = plt.subplots(1,3)
            axrna.boxplot(rna, labels=ticks,  whis=[5, 95], notch=True, showmeans=True,
                          usermedians=rna_median)
            axrna.set_title("RNA-seq (N=%d)"
                            % sum((len(hetrna), len(homorna), len(refrna))),
                            fontsize=8)
            axrna.set_ylabel(ylab, fontsize=8)
            axrna.set_xlabel(xlab, fontsize=8)
            axribo.boxplot(ribo, labels=ticks, whis=[5, 95], notch=True, showmeans=True,
                           usermedians=ribo_median)
            axribo.set_title("Ribosome Profiling (N=%d)"
                            % sum((len(hetrna), len(homorna), len(refrna))),
                             fontsize=8)
                             # + \
                             #  '\n' + 'log2[alt/ref] = %f'
                             # % np.log2(np.mean(homoribo)/np.mean(refribo)),
                             # fontsize=8)
            axribo.set_ylabel(ylab, fontsize=8)
            axribo.set_xlabel(xlab, fontsize=8)
            axnorm.boxplot(norm, labels=ticks, whis=[5, 95], notch=True, showmeans=True,
                           usermedians=norm_median)
            axnorm.set_title("Normalized Ribo (N=%d)"
                            % sum((len(hetrna), len(homorna), len(refrna))),
                            fontsize=8)
            axnorm.set_ylabel(ylab, fontsize=8)
            axnorm.set_xlabel(xlab, fontsize=8)
            # plt.suptitle('log2[alt/ref] = %f'
            #                  % np.log2(np.mean(homoribo)/np.mean(refribo)),
            #                  fontsize=8)
            figmix.suptitle(title, fontsize=8)
            plt.tick_params(axis='both', labelsize=8)
            plt.tight_layout()
            figmix.show()

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
    parser = argparse.ArgumentParser(description='Plots relevant SNPs to interactive scatterplot')
    parser.add_argument('-d', action='store', dest='dir', help='/path/to/root/dir',
                        metavar="DIR", default='/home/mdsherm/Project/SNuPer_results')
    parser.add_argument('--rna', action='store', dest='rna', type=float, metavar="FLOAT",
                        help='threshold for %% of RNA-seq FPKM', default=.5)
    parser.add_argument('--ribo', action='store', dest='ribo', type=float, metavar="FLOAT",
                        help='threshold for %% of ribosomal profiling FPKM', default=.2)
    parser.add_argument('-t', action='store', dest='top', type=int, metavar='INT',
                        help='The number of top results based on log2 ratio of'
                             'homozygous alt vs homozygous ref', default=1000)
    parser.add_argument('--plot-only', action='store_true', dest='plot', default=False,
                        help='Use only with pre-compiled lists. Uses the root directory as path')
    parser.add_argument('-a', action='store', dest='meta', type=float, metavar="FLOAT",
                        help="Threshold for metadata cutoff", default=0.5)
    args = parser.parse_args()
    path_name = args.dir
    rna_thresh = round(args.rna,1)
    ribo_thresh = round(args.ribo,1)
    ignored_states = np.seterr(divide='ignore', invalid='ignore')
    if not args.plot:
        set_list = snp_set(meta_catcher(path_name, 'metadata'), args.meta)
        snp_files = file_globber(path_name, set_list)
        # pickle.dump(snp_files, open(path_name + "/pkl/snp_files.pkl", "wb"))
        sampleGroup = [(snp.rpartition("SNPs/")[-1], [row for row in csv.reader(
            open(snp, "rb"), delimiter='\t')]) for snp in snp_files]
        #
        sampleGroup = [(snp[0], int(snp[1][1][4])-int(snp[1][1][3]) + 1,
                        snp[1][3:]) for snp in sampleGroup if
                       int(snp[1][1][4])-int(snp[1][1][3]) + 1 >= 30]
        SNPs = []
        for snp in sampleGroup:
            ID = snp[0]
            Lengths = snp[1]
            ref = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2] if "0|0" in sample[1]])
            het = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2] if "1|0" in sample[1] or "0|1" in sample[1]])
            alt = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2] if "1|1" in sample[1]])
            raw = np.array([(float(sample[2]), float(sample[3])) for sample
                   in snp[2]])
            try:
                if np.mean(ref, axis=0)[1] < np.mean(
                        het, axis=0)[1] < np.mean(alt, axis=0)[1]:
                    SNPs.append([ID, Lengths, ref, het, alt, raw])
                else:
                    pass
            except IndexError:
                pass
        SNPs = [snp for snp in SNPs if snp[1] >= 30]
        SNPs = [snp for snp in SNPs if min(len(snp[2]), len(snp[3]), len(snp[3])) >2]

        # SNPs = genos[:]
        qLook = {entry[0].split(".snp")[0]: i for (i, entry) in enumerate(SNPs)}
        pickle.dump(SNPs,
                    open(path_name + "/pkl/SNPs.pkl", "wb"))
        pickle.dump(sampleGroup,
                    open(path_name + "/pkl/SG.pkl", "wb"))
        pickle.dump(qLook,
                    open(path_name + "/pkl/qLook.pkl", "wb"))
        SNP_IDs = [snp[0] for snp in SNPs]
        SNP_len = [snp[1] for snp in SNPs]
        SNP_ratio_step = [np.mean(np.array(snp[4]), axis=0)[1] /
                          np.mean(np.array(snp[2]), axis=0)[1]
                          if np.mean(np.array(snp[2]), axis=0)[1] > 0
                             and np.mean(np.array(snp[4]), axis=0)[1] > 0
                          else np.mean(np.array(snp[4]), axis=0)[1]
                          if np.mean(np.array(snp[4]), axis=0)[1] > 0
                          else 0 for snp in SNPs]
        SNP_ratio = [np.log2(step) if step > 0 else 0 for step in SNP_ratio_step]
        percents = []
        for snp in range(len(SNPs)):
            step = [(sample[0], sample[1]) for sample in SNPs[snp][5]]
            percents.append(
                (float(sum(1 for rna in step if rna[0] > rna_thresh)/float(len(step))),
                 float(sum(1 for ribo in step if ribo[1] > ribo_thresh)/float(len(step)))))
        pkl = zip(SNP_IDs, SNP_len, SNP_ratio, percents)
        pickle.dump(pkl, open(path_name + "/pkl/plotzip.pkl", "wb"))
        top = sorted(pkl, key=itemgetter(2), reverse=True)[:args.top]
        pickle.dump(top, open(path_name + "/pkl/top_%d.pkl" % args.top, "wb"))
    elif args.plot and not os.path.isfile(path_name + '/pkl/top_%d.pkl' % args.top):
        print("No precompiled list present!")
        print("Exiting")
    elif args.plot and os.path.isfile(path_name + '/pkl/top_%d.pkl' % args.top):
            top = pickle.load(open(path_name + '/pkl/top_%d.pkl' % args.top))
            # top = pickle.load(open(path_name + '/pkl/plotzip.pkl'))
            # top = [snp for snp in top if snp[2] > 1][:100]
            SNPs = pickle.load(open(path_name + "/pkl/SNPs.pkl", "rb"))
            # sampleGroup = pickle.load(open(path_name + "/pkl/SG.pkl", "rb"))
            qLook = {entry[0].split(".snp")[0]: i for (i, entry) in enumerate(SNPs)}
            SNP_len = [length[1] for length in top]
            sizes = (SNP_len / np.mean(SNP_len)) * 10
            annotes = [snp_IDs[0].split(".snp")[0] for snp_IDs in top]
            colors = [snp_ratio[2] for snp_ratio in top]
            percents = [p100[3] for p100 in top]
            print('{} minimum log2 ratio, {} maximum log2 ratio'.format(round(
                min(colors),4), round(max(colors),4)))
            x = [snp[0] for snp in percents]
            y = [snp[1] for snp in percents]

            # Plots the points above, and can be used to tie in individual SNP IDs
            fig, ax = plt.subplots()
            a = ax.scatter(x, y, color=colors, cmap=plt.get_cmap('YlOrRd'), vmin=min(colors),
                       vmax=max(colors), s=sizes, linewidths=.2, edgecolors='black', alpha=0.8)
            fig.colorbar(a, ticks=None, use_gridspec=False, shrink=0.3,
                         anchor=(0.0, 0.0), pad=0.01, drawedges=False,
                         label='log2(alt/ref)')
            if "chr" not in path_name:
                title = "Whole Genome"
            else:
                title = path_name.split('results/')[1]
            ax.set_title(title)
            ax.set_xlabel('%%RNA-seq FPKM > %f cutoff' % rna_thresh)
            ax.set_ylabel('%%Ribosome Profiling FPKM > %f cutoff' % ribo_thresh)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect('equal')
            ax.text(0, 1.1, "Selection of SNPs based on total samples above > 5 FPKM")
            ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3", alpha=0.35)
            af = AnnoteFinder(x, y, annotes, ax=ax)
            fig.canvas.mpl_connect('button_press_event', af)
            plt.show()
    else:
        print("Unforeseen error")
        print("Exiting")

if __name__ == "__main__":
    main()
