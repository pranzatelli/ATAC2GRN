#!/usr/bin/env python
"""HMMtool: segment Tn5-insertion signal into chromatin states.

Fits a Poisson-emission Hidden Markov Model (built on hmmlearn's
`_BaseHMM`) to Tn5 insertion counts extracted from an ATAC2GRN `Pipes`
run's aligned/deduplicated BAM files (`Output/{project}/{project}-picard.bam`)
over a single genomic region, across one or more samples/projects. The
fitted model segments that region into discrete chromatin states, which are
written out as a bedgraph.

Optionally (only if the RNA-seq CLI arguments below are supplied), per-state
mean ATAC signal is correlated against a gene/transcript's RNA-seq expression
across the same samples. This step is entirely optional: if you omit
`--expression-quant-dir`/`--sample-metadata`/`--genes`, only the ATAC-only
HMM segmentation runs.

This tool must be run *after* an ATAC2GRN `Pipes` run has produced
`Output/{project}/{project}-picard.bam` files (see `../Pipes/README.md`).
See `README.md` in this folder for full usage and output-file documentation.
"""

import argparse
import os
import warnings
from glob import glob

import numpy as np
import pandas as pd
import pysam
from hmmlearn.base import _BaseHMM
from hmmlearn.utils import normalize
from scipy.stats import poisson

warnings.filterwarnings("ignore")


class PoissonHMM(_BaseHMM):
    """Poisson-emission HMM. This is the validated scientific core of
    HMMtool and is intentionally left algorithmically unchanged from the
    original script, other than fixing the `_generate_sample_from_state`
    bug noted below."""

    def __init__(self, n_components=1,
                 startprob_prior=1.0, transmat_prior=1.0,
                 lambda_prior=1.0, identity_weight=0.0,
                 algorithm="viterbi", random_state=None,
                 n_iter=10, tol=float('-inf'), verbose=True,
                 params="stl", init_params="stl"):
        _BaseHMM.__init__(self, n_components,
                           startprob_prior=startprob_prior,
                           transmat_prior=transmat_prior,
                           algorithm=algorithm,
                           random_state=random_state,
                           n_iter=n_iter, tol=tol, verbose=verbose,
                           params=params, init_params=init_params)
        self.lambda_prior = lambda_prior
        self.identity_weight = identity_weight

    def _init(self, X, lengths=None):
        # hmmlearn's _BaseHMM._init gained a `lengths` parameter after the
        # hmmlearn version this script was originally written against; this
        # override is pinned/tested against hmmlearn==0.3.3 (see
        # atac2grn.def) and accepts+forwards it for compatibility.
        self.n_features = getattr(self, "n_features", X.shape[1])
        super()._init(X, lengths)
        if 'l' in self.init_params:
            self.lambdas_ = np.random.random((self.n_components, self.n_features))

    def _get_n_fit_scalars_per_param(self):
        nc = self.n_components
        nf = self.n_features
        return {
            "s": nc - 1,
            "t": nc * (nc - 1),
            "l": nc * nf
        }

    def _compute_log_likelihood(self, X):
        log_likelihoods = np.zeros((X.shape[0], self.n_components))
        for i in range(self.n_components):
            d = {}
            for j in range(X.shape[1]):
                l = self.lambdas_[i, j]
                d[j] = {}
                for k in range(X[:, j].max() + 1):
                    d[j][k] = poisson.logpmf(k, l)
            log_features = np.zeros(X.shape)
            for feature, x in enumerate(X.T):
                log_features[:, feature] = np.vectorize(d[feature].get)(x)
            log_likelihoods[:, i] = log_features.sum(axis=1)
        return log_likelihoods

    def _generate_sample_from_state(self, state, random_state=None):
        # Bug fix: the original referenced an undefined variable `x` here
        # (a leftover from _compute_log_likelihood's loop variable of the
        # same name). This would raise NameError the moment `.sample()` was
        # ever called on a fitted model. Use this state's lambdas instead.
        state_lambdas = self.lambdas_[state, :]
        X = poisson.rvs(state_lambdas, random_state=random_state)
        return np.array(X)

    def _initialize_sufficient_statistics(self):
        stats = super()._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
        stats['obs'] = np.zeros((self.n_components, self.n_features))
        return stats

    def _accumulate_sufficient_statistics(self, stats, X, framelogprob,
                                           posteriors, fwdlattice, bwdlattice):
        super()._accumulate_sufficient_statistics(
            stats, X, framelogprob, posteriors, fwdlattice, bwdlattice)
        if 'l' in self.params:
            stats['post'] += posteriors.sum(axis=0)
            stats['obs'] += np.dot(posteriors.T, X)

    def _do_mstep(self, stats):
        super()._do_mstep(stats)
        if 't' in self.params:
            for c in range(self.n_components):
                if self.transmat_[c, c] < (1 - 1.0 / (self.identity_weight + 1)):
                    self.transmat_[c, c] = 0
                    normalize(self.transmat_, axis=1)
                    self.transmat_[c, c] = self.identity_weight
            normalize(self.transmat_, axis=1)
        if 'l' in self.params:
            denom = stats['post'][:, None]
            self.lambdas_ = stats['obs'] / denom


# --- ATAC signal extraction (Tn5-offset BAM-to-array logic, unchanged) ------

def construct_atac_array(chrom, start, end, sample_dir):
    """Build a per-base Tn5-insertion count array for one sample's BAM over
    [start, end) on chrom. `sample_dir` is a project output directory
    expected to contain exactly one `*picard.bam` file (as produced by the
    ATAC2GRN Snakefile, e.g. Output/{project}/{project}-picard.bam)."""
    length = end - start
    bamfile = glob(os.path.join(sample_dir, '*picard.bam'))
    Y = np.zeros(length, dtype=np.uint32)
    if len(bamfile) != 1:
        print(bamfile)
        raise RuntimeError("A non-standard number of bamfiles found in " + sample_dir + ".")
    af = pysam.AlignmentFile(bamfile[0], 'rb')
    for read in af.fetch(chrom, start, end):
        if read.is_reverse:
            Tn5_binding = read.reference_end - start - 5
        else:
            Tn5_binding = read.reference_start - start + 4
        if 0 <= Tn5_binding < length:
            Y[Tn5_binding] += 1
    return Y


def construct_all_atac(chrom, start, end, atac_output_dir, exclude_samples=None):
    """Build a {project: per-base Tn5 count array} DataFrame for every
    project subfolder under an ATAC2GRN Pipes `Output/` directory."""
    experiments = sorted(glob(os.path.join(atac_output_dir, '*')))
    experiments = [e for e in experiments if os.path.isdir(e)]
    DF_able = {}
    for experiment in experiments:
        label = os.path.basename(experiment.rstrip('/\\'))
        DF_able[label] = construct_atac_array(chrom, start, end, experiment)
    DF = pd.DataFrame(DF_able)
    if exclude_samples:
        DF = DF.drop(columns=[c for c in exclude_samples if c in DF.columns])
    return DF


# --- Optional RNA-seq expression correlation --------------------------------

def read_salm(filepath):
    """Read a salmon quant.sf file, returning {transcript_id: TPM}."""
    i = 3  # TPM column in salmon's quant.sf
    returnD = {}
    with open(filepath) as openfile:
        lines = openfile.read().split('\n')
    for line in lines[1:-1]:
        tabs = line.split('\t')
        returnD[tabs[0]] = float(tabs[i])
    return returnD


def read_sample_metadata(filepath, rna_col, atac_col):
    """Read a user-supplied sample-metadata file mapping RNA-seq sample IDs
    to ATAC-seq sample/project IDs. Must contain at least `rna_col` and
    `atac_col` columns (TSV or CSV, auto-detected by file extension/sniffing
    the first line for tabs). Returns {rna_sample_id: atac_sample_id}."""
    sep = '\t' if filepath.lower().endswith(('.tsv', '.txt')) else None
    # encoding='utf-8-sig' transparently strips a leading byte-order-mark,
    # which Excel/Windows-authored CSV/TSV files commonly include.
    DF = pd.read_csv(filepath, sep=sep, engine='python' if sep is None else 'c',
                      encoding='utf-8-sig')
    if rna_col not in DF.columns or atac_col not in DF.columns:
        raise ValueError(
            "Sample metadata file %s must contain columns %r and %r; found %r"
            % (filepath, rna_col, atac_col, list(DF.columns))
        )
    return dict(zip(DF[rna_col].astype(str), DF[atac_col].astype(str)))


def assemble_expression_df(quant_dir, sample_metadata_path, rna_col, atac_col,
                            exclude_samples=None):
    """Build a {atac_sample_id: {transcript_id: TPM}} DataFrame from a
    directory of per-sample salmon quant.sf files plus a sample-metadata
    file mapping RNA-seq sample IDs (subfolder names under quant_dir) to
    ATAC-seq sample/project IDs (subfolder names under --atac-output-dir)."""
    DTFrame = {}
    experiments = sorted(glob(os.path.join(quant_dir, '*')))
    for exp in experiments:
        if not os.path.isdir(exp):
            continue
        experiment = os.path.basename(exp.rstrip('/\\'))
        filepath = os.path.join(exp, 'quant.sf')
        if os.path.exists(filepath):
            DTFrame[experiment] = read_salm(filepath)
    DF = pd.DataFrame(DTFrame).fillna(0)
    DF = DF[(DF.T != 0).any()]
    DF.index = [value.split('.')[0] for value in DF.index]
    id_map = read_sample_metadata(sample_metadata_path, rna_col, atac_col)
    DF.columns = DF.columns.map(lambda x: id_map[x] if x in id_map else x)
    if exclude_samples:
        DF = DF.drop(columns=[c for c in exclude_samples if c in DF.columns])
    return DF


# --- Output writers ----------------------------------------------------------

def write_bed_from_scores(scores, chrom, index, output_path):
    with open(output_path, 'w') as openfile:
        start = index[0]
        length = 1
        state = scores[0]
        for i in range(len(scores)):
            if scores[i] != state:
                openfile.write(chrom + '\t' + str(start) + '\t' + str(start + length) + '\t' + str(state) + '\n')
                state = scores[i]
                start = index[i]
                length = 0
            length += 1
        openfile.write(chrom + '\t' + str(start) + '\t' + str(start + length) + '\t' + str(state) + '\n')


def size_d(states):
    """Print the mean run-length (in bins) of each non-initial state."""
    d = {}
    current = 'droppable'
    c = 0
    for state in states:
        if state != current:
            d.setdefault(current, []).append(c)
            current = state
            c = 1
        else:
            c += 1
    d.setdefault(current, []).append(c)
    del d['droppable']
    for state in d:
        print(str(state) + '\t' + str(np.mean(d[state])))


# --- CLI ---------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Segment Tn5-insertion signal from ATAC2GRN Pipes output "
                    "into chromatin states via a Poisson HMM, with optional "
                    "RNA-seq expression correlation.")
    region = parser.add_argument_group("region")
    region.add_argument("--chrom", required=True, help="Chromosome, e.g. chr11")
    region.add_argument("--start", required=True, type=int, help="Region start (0-based)")
    region.add_argument("--end", required=True, type=int, help="Region end (exclusive)")
    region.add_argument("--label", required=True,
                         help="Run label; used as the output subfolder name and "
                              "output-file prefix, e.g. 'sgonly_FOLH1'")

    atac = parser.add_argument_group("ATAC input")
    atac.add_argument("--atac-output-dir", required=True,
                       help="An ATAC2GRN Pipes run's Output/ directory "
                            "(contains one subfolder per project, each with "
                            "a {project}-picard.bam file).")
    atac.add_argument("--exclude-samples", default=None,
                       help="Comma-separated list of sample/project names to "
                            "drop from the ATAC (and, if used, expression) "
                            "data before fitting, e.g. for known-bad samples. "
                            "Names not present are silently ignored.")

    expr = parser.add_argument_group(
        "RNA-seq expression correlation (optional; omit all of these to skip)")
    expr.add_argument("--genes", default=None,
                       help="Comma-separated gene/transcript ID(s) (matching "
                            "the RNA-seq quantification file's transcript "
                            "IDs, version suffix optional) whose summed "
                            "expression is correlated against per-bin ATAC "
                            "signal. Required (with the two options below) "
                            "to run expression correlation.")
    expr.add_argument("--expression-quant-dir", default=None,
                       help="Directory containing one subfolder per RNA-seq "
                            "sample, each with a salmon quant.sf file.")
    expr.add_argument("--sample-metadata", default=None,
                       help="TSV/CSV file mapping RNA-seq sample IDs "
                            "(--expression-quant-dir subfolder names) to "
                            "ATAC-seq sample/project IDs (--atac-output-dir "
                            "subfolder names). See README.md for the "
                            "expected columns.")
    expr.add_argument("--metadata-rna-col", default="rna_sample",
                       help="Column in --sample-metadata holding RNA-seq "
                            "sample IDs (default: rna_sample).")
    expr.add_argument("--metadata-atac-col", default="atac_sample",
                       help="Column in --sample-metadata holding ATAC-seq "
                            "sample/project IDs (default: atac_sample).")

    model = parser.add_argument_group("HMM / preprocessing parameters")
    model.add_argument("--output-dir", default="Output",
                        help="Base output directory; results are written to "
                             "{output-dir}/{label}/ (default: Output)")
    model.add_argument("--n-states", type=int, default=25,
                        help="Number of HMM hidden states (default: 25)")
    model.add_argument("--n-iter", type=int, default=50,
                        help="Number of Baum-Welch training iterations (default: 50)")
    model.add_argument("--rolling-window", type=int, default=101,
                        help="Triangular rolling-window size (in bins) applied "
                             "to the raw per-base Tn5 counts before HMM "
                             "fitting (default: 101)")
    model.add_argument("--random-seed", type=int, default=None,
                        help="Random seed for the HMM (default: none/nondeterministic)")

    args = parser.parse_args(argv)

    expr_flags = [args.genes, args.expression_quant_dir, args.sample_metadata]
    if any(expr_flags) and not all(expr_flags):
        parser.error(
            "--genes, --expression-quant-dir, and --sample-metadata must all "
            "be supplied together to run expression correlation, or all "
            "omitted to skip it."
        )
    return args


def main(argv=None):
    args = parse_args(argv)
    exclude_samples = (
        [s.strip() for s in args.exclude_samples.split(',') if s.strip()]
        if args.exclude_samples else None
    )

    out_dir = os.path.join(args.output_dir, args.label)
    os.makedirs(out_dir, exist_ok=True)

    aDF = construct_all_atac(args.chrom, args.start, args.end,
                              args.atac_output_dir, exclude_samples=exclude_samples)
    aDF.index = aDF.index + args.start

    do_correlation = args.genes is not None
    gene = None
    if do_correlation:
        genes = [g.strip() for g in args.genes.split(',') if g.strip()]
        gDF = assemble_expression_df(
            args.expression_quant_dir, args.sample_metadata,
            args.metadata_rna_col, args.metadata_atac_col,
            exclude_samples=exclude_samples,
        )
        cols = gDF.columns.intersection(aDF.columns)
        if len(cols) == 0:
            raise RuntimeError(
                "No overlap between ATAC sample/project IDs and RNA-seq "
                "sample IDs after applying --sample-metadata; check the "
                "metadata file's column values."
            )
        aDF = aDF[cols]
        gDF = gDF[cols]
        gene = gDF.loc[genes].sum()

    aDFr = (aDF.rolling(args.rolling_window, win_type='triang', center=True, closed='both')
                .sum().dropna() * 100).astype('uint32')

    if do_correlation:
        corrs = aDFr.corrwith(gene, axis=1)
        aDFcorr = corrs.fillna(0)
        np.save(os.path.join(out_dir, 'corr.npy'), aDFcorr)
        aDFcorr = pd.DataFrame(aDFcorr)

    lr = PoissonHMM(n_components=args.n_states, n_iter=args.n_iter,
                    random_state=args.random_seed)
    lr.fit(aDFr)
    scores = lr.predict(aDFr)
    np.save(os.path.join(out_dir, 'scores.npy'), scores)
    np.save(os.path.join(out_dir, 'lambdas.npy'), lr.lambdas_)
    np.save(os.path.join(out_dir, 'transmat.npy'), lr.transmat_)

    if do_correlation:
        gene.to_csv(os.path.join(out_dir, args.label + '.csv'))

    # Note: the original script wrote the bedgraph using the pre-rolling
    # aDF.index while scores are predicted on the post-rolling/dropna aDFr
    # (which can be shorter at the edges depending on rolling-window size),
    # a latent length mismatch. Using aDFr.index here (which scores is
    # actually aligned to) is a correctness fix, not an algorithm change.
    write_bed_from_scores(scores, args.chrom, aDFr.index,
                           os.path.join(out_dir, args.label + '.bedgraph'))

    size_d(scores)
    print(pd.Series(scores).value_counts())
    if do_correlation:
        aDFcorr['State'] = scores
        print(aDFcorr.groupby('State').mean())


if __name__ == '__main__':
    main()
