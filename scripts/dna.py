#!/usr/bin/env python3

import numpy as np
import scipy, scipy.stats
import argparse
import os
import logging
import subprocess
import json
from tempfile import NamedTemporaryFile
from sklearn.preprocessing import scale

DNA_VERSION='0.3'

class FileInput(object):
    """Simple data object. Class to store and validate files"""

    def __init__(
        self, first, second, first_headers, second_headers, pw_map, ortho=None
    ):
        super(FileInput, self).__init__()
        self.first = first
        self.second = second
        self.first_headers = first_headers
        self.second_headers = second_headers
        self.pw_map = pw_map
        self.ortho = ortho

    def validate(self):
        if not os.path.exists(self.first):
            raise FileNotFoundError("File {} does not exist".format(self.first))
        if not os.path.exists(self.second):
            raise FileNotFoundError(
                "File {} does not exist".format(self.second)
            )
        if not os.path.exists(self.first_headers):
            raise FileNotFoundError(
                "File {} does not exist".format(self.first_headers)
            )
        if not os.path.exists(self.second_headers):
            raise FileNotFoundError(
                "File {} does not exist".format(self.second_headers)
            )
        if not os.path.exists(self.pw_map):
            raise FileNotFoundError(
                "File {} does not exist".format(self.pw_map)
            )
        if self.ortho and not os.path.exists(self.ortho):
            raise FileNotFoundError("File {} does not exist".format(self.ortho))


class ExprMatrix(object):
    """Simple data object. Set up and store expression matrices"""

    def __init__(self, in_file, in_headers):
        super(ExprMatrix, self).__init__()
        self.data = np.genfromtxt(in_file, delimiter="\t")
        self.data = scale(self.data, axis=0, with_std=False)
        self.data = scale(self.data, axis=1, with_std=False)
        self.headers = []
        with open(in_headers, "r") as inh:
            for line in inh:
                for field in line.strip().split("\t"):
                    self.headers.append(field)


class PathwayMap(object):
    """Simple data object. Class for Pathway->Gene mappings"""

    def __init__(self, in_file, min_membership, max_membership):
        super(PathwayMap, self).__init__()
        self.map = {}
        with open(in_file, "r") as handle:
            for line in handle:
                pw, members = line.strip().split("\t")
                members = members.split("|")
                len_filter = False
                if min_membership != 0:
                    len_filter = len(members) >= min_membership
                if max_membership != 0:
                    len_filter = len(members) <= max_membership
                if len_filter:
                    self.map[pw] = members

        for k in self.map:
            self.map[k] = list(set(self.map[k]))


class JoinedOrthoMatrix(object):
    """docstring for JoinedOrthoMatrix"""

    def __init__(self, X1, X2, H1, H2, ortho):
        super(JoinedOrthoMatrix, self).__init__()
        self.X1 = X1
        self.X2 = X2
        self.H1 = H1
        self.H2 = H2
        self.H1M = {}
        self.H2M = {}
        self.ortho_file = ortho
        self.ortho = {}
        self.row_index = np.arange(self.X1.shape[0] + self.X2.shape[0])

        with open(self.ortho_file, "r") as of_handle:
            for line in of_handle:
                fields = line.strip().split()
                if not fields[0] in self.ortho:
                    self.ortho[fields[0]] = set()
                if not fields[1] in self.ortho:
                    self.ortho[fields[1]] = set()

                self.ortho[fields[0]].add(fields[1])
                self.ortho[fields[1]].add(fields[0])

        for i, g in enumerate(H1):
            self.H1M[g] = i
        for i, g in enumerate(H2):
            self.H2M[g] = i

    def permute(self):
        self.row_index = np.random.permutation(self.row_index)

    def get_subset(self, genes):
        no_ortho = 0
        no_expr = 0
        index_1 = []
        index_2 = []
        # Here we get all genes in a pathway, check which species they belong to
        # and also check that they have at least one ortholog and that they are
        # expressed. If that is satisfied, we add the column indices of the ortholog
        # and as many copies of the gene into the column indices which we then use
        # to subset
        for gene in genes:
            # Sanity checks
            if gene not in self.ortho:
                no_ortho += 1
                continue
            if gene not in self.H1M and gene not in self.H2M:
                no_expr += 1
                continue
            # Set up sets that explain the orthology
            first = gene in self.H1M
            second = gene in self.H2M
            if first:
                for orth in self.ortho[gene]:
                    index_1.append(self.H1M[gene])
                    index_2.append(self.H2M[orth])
            elif second:
                for orth in self.ortho[gene]:
                    index_2.append(self.H2M[gene])
                    index_1.append(self.H1M[orth])
            else:
                raise ValueError(
                    "Gene {} was neither in first nor second set. "
                    + "This probably shouldn't have happened.".format(gene)
                )
        X1_ = self.X1[:, index_1]
        X2_ = self.X2[:, index_2]
        X12 = np.concatenate([X1_, X2_])
        # Apply current permutation
        X12 = X12[self.row_index, :]
        # Set up a mostly useless header file to submit to Seidr. If we need to
        # compare networks at the gene level, we need to come up with something more
        # sophisticated than that
        header = ["G{}".format(x) for x in range(X12.shape[1])]
        s1 = X1_.shape[0]
        s2 = X2_.shape[1]

        return (X12[0:s1, :], X12[s1 : (s1 + s2), :], header)


class JoinedMatrix(object):
    """Class to operate on a joined matrix of count data"""

    def __init__(self, X1, X2, H1, H2):
        super(JoinedMatrix, self).__init__()
        logging.debug("Joining matrices...")
        self.h_isec = set()
        self.headers = []
        self.shapes = []

        logging.debug("Finding gene overlap...")
        logging.debug("Original shapes: {} || {}".format(X1.shape, X2.shape))
        self.shapes = [X1.shape, X2.shape]
        # Make a set to subset genes that are in both samples
        h1_set = set()
        for n in H1:
            h1_set.add(n)
        # Get gene set intersect
        for n in H2:
            if n in h1_set:
                self.h_isec.add(n)
        # Set up header data structures
        self.headers = list(self.h_isec)
        self.h_dict = {}
        for i, n in enumerate(self.headers):
            self.h_dict[n] = i

        index_1 = []
        for n in H1:
            if n in self.h_isec:
                index_1.append(self.h_dict[n])

        index_2 = []
        for n in H2:
            if n in self.h_isec:
                index_2.append(self.h_dict[n])

        X1_ = X1[:, index_1]
        X2_ = X2[:, index_2]

        logging.debug("Filtered shapes: {} || {}".format(X1_.shape, X2_.shape))
        self.data = np.concatenate([X1_, X2_])

        logging.debug("Concat shape: {}".format(self.data.shape))
        self.row_index = np.arange(self.data.shape[0])

    def permute(self):
        self.row_index = np.random.permutation(self.row_index)
        self.data = self.data[self.row_index, :]

    def get_subset(self, genes):
        index = []
        header = []
        for g in genes:
            if g in self.h_dict:
                index.append(self.h_dict[g])
                header.append(g)

        s1 = self.shapes[0][0]
        s2 = self.shapes[1][0]

        X_ = self.data[:, index]

        return (X_[0:s1, :], X_[s1 : (s1 + s2), :], header)


def dE(mat1, mat2, p):
    X = np.absolute(mat1 - mat2)
    X = np.power(X, p)
    s = np.sum(X)
    ng = mat1.shape[0]
    e = (ng * (ng - 1)) / 2

    # Entire pathway statistic
    d_all = np.power(s / e, 1 / p)

    # Per gene statistics
    d_sub = []
    for i in range(ng):
        X = np.absolute(mat1[:, i], mat2[:, i])
        X = np.power(X, p)
        s = np.sum(X)
        d_sub.append(np.power(s / e, 1 / p))

    return {"All": d_all, "Sub": d_sub}


def main(args):

    files = FileInput(
        args.first_expr,
        args.second_expr,
        args.first_header,
        args.second_header,
        args.map,
    )
    logging.info("Validating files...")
    files.validate()

    logging.info("Parsing pathway map {}".format(files.pw_map))
    pw_map = PathwayMap(
        files.pw_map, args.min_pathway_membership, args.max_pathway_membership
    )
    logging.info("Have {} pathways".format(len(pw_map.map)))

    logging.info(
        "Reading expression/headers of {} / {}".format(
            files.first, files.first_headers
        )
    )
    X1 = ExprMatrix(files.first, files.first_headers)

    logging.info(
        "Reading expression/headers of {} / {}".format(
            files.second, files.second_headers
        )
    )
    X2 = ExprMatrix(files.second, files.second_headers)

    if args.ortho:
        X = JoinedOrthoMatrix(
            X1.data, X2.data, X1.headers, X2.headers, args.ortho
        )
    else:
        X = JoinedMatrix(X1.data, X2.data, X1.headers, X2.headers)

    pathways = []
    if not args.pathway:
        pathways = [p for p in pw_map.map]
    else:
        pathways = args.pathway.split(",")

    print(
        "Pathway",
        "Gene",
        "P",
        "RefMean",
        "RefVar",
        "PermMean",
        "PermVar",
        "Diff",
        sep="\t",
    )

    out_data = {}

    for pw in pathways:

        logging.info("Starting pathway {}".format(pw))
        logging.info("{} genes in this pathway".format(len(pw_map.map[pw])))
        X1_, X2_, header = X.get_subset(pw_map.map[pw])
        logging.info("{} genes with expression data".format(len(header)))

        if len(header) < args.min_pathway_membership:
            logging.warning(
                "Only {} genes with expression data in pathway {}. Skipping...".format(
                    len(header), pw
                )
            )
            continue

        handle_g = NamedTemporaryFile(mode="w", delete=False)
        for g in header:
            handle_g.write("{}\n".format(g))
        handle_g.close()

        nboot = args.permutations
        nref = args.reference_runs
        p = args.p_norm
        outcomes_ref = []
        outcomes = []

        for i in range(nboot + nref):
            if i < nref:
                logging.info("Reference run {}/{}".format(i, nref))
            else:
                logging.info("Permutation run {}/{}".format(i - nref, nboot))
            if i >= nref:
                X.permute()
                X1_, X2_, header = X.get_subset(pw_map.map[pw])
            logging.debug(
                "Shapes after subset: {} || {}".format(X1_.shape, X2_.shape)
            )
            handle_1 = NamedTemporaryFile(mode="w", delete=False)
            handle_2 = NamedTemporaryFile(mode="w", delete=False)

            np.savetxt(handle_1, X1_, delimiter="\t")
            np.savetxt(handle_2, X2_, delimiter="\t")

            handle_1.close()
            handle_2.close()

            wrapper_1 = SeidrWrapper(
                handle_1.name,
                handle_g.name,
                str(args.threads),
                X1_.shape,
                args.fast,
                fail_on_error=args.strict,
            )
            wrapper_2 = SeidrWrapper(
                handle_2.name,
                handle_g.name,
                str(args.threads),
                X2_.shape,
                args.fast,
                fail_on_error=args.strict,
            )

            wrapper_1.pearson()
            wrapper_1.spearman()
            wrapper_1.pcor()
            if not args.fastest:
                wrapper_1.mi()
                wrapper_1.narromi()
                if not args.faster:
                    wrapper_1.plsnet()
                    wrapper_1.llr()
                    wrapper_1.tigress()
                    wrapper_1.genie3()
                    wrapper_1.elnet()
            wrapper_1.aggregate()
            wrapper_1.adjacency()

            wrapper_2.pearson()
            wrapper_2.spearman()
            wrapper_2.pcor()
            if not args.fastest:
                wrapper_2.mi()
                wrapper_2.narromi()
                if not args.fastest:
                    wrapper_2.plsnet()
                    wrapper_2.llr()
                    wrapper_2.tigress()
                    wrapper_2.genie3()
                    wrapper_2.elnet()
            wrapper_2.aggregate()
            wrapper_2.adjacency()

            dX1 = np.genfromtxt(wrapper_1.adj, delimiter="\t")
            dX2 = np.genfromtxt(wrapper_2.adj, delimiter="\t")

            d = dE(dX1, dX2, p)

            if i < nref:
                outcomes_ref.append(d)
            else:
                outcomes.append(d)

            os.remove(handle_1.name)
            os.remove(handle_2.name)
            wrapper_1.clean()
            wrapper_2.clean()

        os.remove(handle_g.name)
        # # Calculate P values as in the paper
        # d0 = outcomes[0]
        # I = 0
        # for i in outcomes[1:]:
        #   if d0 <= i:
        #     I += 1
        # p = (I + 1) / (nboot)
        desc1 = scipy.stats.describe([x["All"] for x in outcomes_ref])
        desc2 = scipy.stats.describe([x["All"] for x in outcomes])
        pc = scipy.stats.mannwhitneyu(
            [np.mean([x["All"] for x in outcomes_ref])],
            [x["All"] for x in outcomes],
            alternative="greater",
        )
        print(
            pw,
            "All",
            pc.pvalue,
            desc1.mean,
            desc1.variance,
            desc2.mean,
            desc2.variance,
            desc1.mean - desc2.mean,
            sep="\t",
        )
        for i in range(len(header)):
            desc1 = scipy.stats.describe([x["Sub"][i] for x in outcomes_ref])
            desc2 = scipy.stats.describe([x["Sub"][i] for x in outcomes])
            pc = scipy.stats.mannwhitneyu(
                [np.mean([x["Sub"][i] for x in outcomes_ref])],
                [x["Sub"][i] for x in outcomes],
                alternative="greater",
            )
            print(
                pw,
                header[i],
                pc.pvalue,
                desc1.mean,
                desc1.variance,
                desc2.mean,
                desc2.variance,
                desc1.mean - desc2.mean,
                sep="\t",
            )
        out_data[pw] = {"ref": outcomes_ref, "perm": outcomes}

    with open("out_data.json", "w") as od:
        json.dump(out_data, od)


class SeidrWrapper(object):
    """docstring for SeidrWrapper"""

    def __init__(self, expr, genes, threads, shape, fast, fail_on_error=False):
        super(SeidrWrapper, self).__init__()
        self.genes = genes
        self.expr = expr
        self.threads = threads
        self.results = []
        self.shape = shape
        self.adj = ""
        self.ensemble = "100" if fast else "1000"
        self.fail_on_error = fail_on_error

    def run_seidr_proc(self, cmd):
        p = subprocess.run(cmd, capture_output=True)
        if p.returncode != 0:
            logging.debug(p.stderr)
            if self.fail_on_error:
                raise RuntimeError(f"Command ${p.args} failed")

    def pearson(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "correlation",
            "-m",
            "pearson",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-o",
            handle.name,
        ]
        logging.debug("Inferring Pearson")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-A",
            "-r",
            "-u",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "Pearson",
            "-o",
            handle_i.name,
            "-F",
            "lm",
        ]
        logging.debug("Importing Pearson")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def spearman(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "correlation",
            "-m",
            "spearman",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-o",
            handle.name,
        ]
        logging.debug("Inferring Spearman")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-A",
            "-r",
            "-u",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "Spearman",
            "-o",
            handle_i.name,
            "-F",
            "lm",
        ]
        logging.debug("Importing Spearman")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def pcor(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "pcor",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-o",
            handle.name,
        ]
        logging.debug("Inferring PCor")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-A",
            "-r",
            "-u",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "PCor",
            "-o",
            handle_i.name,
            "-F",
            "lm",
        ]
        logging.debug("Importing PCor")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def mi(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "mi",
            "-m",
            "RAW",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-b",
            "10",
            "-f",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
        ]
        logging.debug("Inferring MI")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-r",
            "-u",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "MI",
            "-o",
            handle_i.name,
            "-F",
            "lm",
        ]
        logging.debug("Importing MI")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        mi_raw = handle.name

        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "mi",
            "-m",
            "CLR",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-f",
            "-M",
            mi_raw,
            "-o",
            handle.name,
        ]
        logging.debug("Inferring CLR")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_j = NamedTemporaryFile(mode="w", delete=False)
        handle_j.close()
        cmd = [
            "seidr",
            "import",
            "-r",
            "-u",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "CLR",
            "-o",
            handle_j.name,
            "-F",
            "lm",
        ]
        logging.debug("Importing CLR")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_j.name)
        os.remove(handle.name)

        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "mi",
            "-m",
            "ARACNE",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-f",
            "-M",
            mi_raw,
            "-o",
            handle.name,
        ]
        logging.debug("Inferring ARACNE")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_k = NamedTemporaryFile(mode="w", delete=False)
        handle_k.close()
        cmd = [
            "seidr",
            "import",
            "-r",
            "-u",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "ARACNE",
            "-o",
            handle_k.name,
            "-F",
            "lm",
        ]
        logging.debug("Importing ARACNE")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_k.name)
        os.remove(handle.name)
        os.remove(mi_raw)

    def plsnet(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "plsnet",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-e",
            self.ensemble,
            "-f",
            "-s",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
            "-p",
            str(self.shape[1] - 1),
        ]
        logging.debug("Inferring PLSNet")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-z",
            "-r",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "PLSNET",
            "-o",
            handle_i.name,
            "-F",
            "m",
        ]
        logging.debug("Importing PLSNet")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def narromi(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "narromi",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-f",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
        ]
        logging.debug("Inferring Narromi")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-r",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "Narromi",
            "-o",
            handle_i.name,
            "-F",
            "m",
        ]
        logging.debug("Importing Narromi")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def llr(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "llr-ensemble",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
            "-p",
            str(self.shape[1] - 3),
            "-P",
            str(self.shape[1] - 2),
            "-e",
            self.ensemble,
            "-x",
            str(self.shape[0]),
            "-X",
            str(self.shape[0]),
        ]
        logging.debug("Inferring LLR")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-z",
            "-r",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "Narromi",
            "-o",
            handle_i.name,
            "-F",
            "m",
        ]
        logging.debug("Importing LLR")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def elnet(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "el-ensemble",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
            "-p",
            str(self.shape[1] - 3),
            "-P",
            str(self.shape[1] - 2),
            "-e",
            self.ensemble,
            "-x",
            str(self.shape[0]),
            "-X",
            str(self.shape[0]),
        ]
        logging.debug("Inferring ElNet")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-z",
            "-r",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "ElNet",
            "-o",
            handle_i.name,
            "-F",
            "m",
        ]
        logging.debug("Importing ElNet")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def tigress(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "tigress",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
            "-e",
            self.ensemble,
        ]
        logging.debug("Inferring Tigress")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-z",
            "-r",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "tigress",
            "-o",
            handle_i.name,
            "-F",
            "m",
        ]
        logging.debug("Importing Tigress")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def genie3(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = [
            "genie3",
            "-i",
            self.expr,
            "-g",
            self.genes,
            "-s",
            "-f",
            "-O",
            self.threads,
            "-B",
            str(self.shape[1]),
            "-o",
            handle.name,
            "-m",
            str(self.shape[1] - 1),
            "-n",
            self.ensemble,
        ]
        logging.debug("Inferring GENIE3")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        handle_i = NamedTemporaryFile(mode="w", delete=False)
        handle_i.close()
        cmd = [
            "seidr",
            "import",
            "-z",
            "-r",
            "-i",
            handle.name,
            "-g",
            self.genes,
            "-f",
            "-n",
            "GENIE3",
            "-o",
            handle_i.name,
            "-F",
            "m",
        ]
        logging.debug("Importing GENIE3")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.results.append(handle_i.name)
        os.remove(handle.name)

    def aggregate(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = ["seidr", "aggregate", "-f", "-m", "irp", "-o", handle.name]
        for r in self.results:
            cmd.append(r)
        logging.debug("Aggregating data")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.aggregate = handle.name

    def adjacency(self):
        handle = NamedTemporaryFile(mode="w", delete=False)
        handle.close()
        cmd = ["seidr", "adjacency", "-f", "-o", handle.name, self.aggregate]
        logging.debug("Getting graph adjacency")
        logging.debug("Running cmd: {}".format(cmd))
        self.run_seidr_proc(cmd)
        self.adj = handle.name

    def clean(self):
        for r in self.results:
            os.remove(r)
        os.remove(self.aggregate)
        os.remove(self.adj)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"Differential network analysis wrapper: ${DNA_VERSION}"
    )
    parser.add_argument(
        "-1",
        "--first-expr",
        required=True,
        help="First network expression file (without headers)",
    )
    parser.add_argument(
        "-a",
        "--first-header",
        required=True,
        help="First network column headers",
    )
    parser.add_argument(
        "-2",
        "--second-expr",
        required=True,
        help="Second network expression file (without headers)",
    )
    parser.add_argument(
        "-b",
        "--second-header",
        required=True,
        help="Second network column headers",
    )
    parser.add_argument(
        "-m", "--map", required=True, help="Gene-pathway map file"
    )
    parser.add_argument("-O", "--ortho", help="Cross species orthology file")
    parser.add_argument("-p", "--pathway", help="Pathway to analyze")
    parser.add_argument(
        "-P",
        "--permutations",
        default=100,
        type=int,
        help="Number of permutation samples",
    )
    parser.add_argument(
        "-r",
        "--reference-runs",
        default=25,
        type=int,
        help="Number of runs for the reference networks",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help="Number of threads for Seidr GRN inference",
    )
    parser.add_argument(
        "-n",
        "--p-norm",
        default=2,
        type=int,
        help="Value of p for matrix p-norm",
    )
    parser.add_argument(
        "--min-pathway-membership",
        type=int,
        help="Minimum number of genes in a pathway",
        default=10,
    )
    parser.add_argument(
        "--max-pathway-membership",
        type=int,
        help="Maximum number of genes in a pathway",
        default=0,
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="count",
        help="Be more quiet with logging",
        default=0,
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Use faster, but less accurate setting in Seidr",
    )
    parser.add_argument(
        "--faster",
        action="store_true",
        help="Only run fast and intermediate Seidr algorithms",
    )
    parser.add_argument(
        "--fastest",
        action="store_true",
        help="Only run fastest Seidr algorithms",
    )
    parser.add_argument(
        "--strict",
        action='store_true',
        help="Fail if any subprocess fails",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG)

    main(args)
