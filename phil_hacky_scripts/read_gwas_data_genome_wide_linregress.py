from phil import *
import virscan_util
import numpy as np
from scipy.stats import linregress, mannwhitneyu, pearsonr
import parse_tsv

basedir = '/home/pbradley/csdat/virscan/IDS065/martin_gwas/'

ALL_FEATURE_TAGS = 'ALL_FEATURE_TAGS'

with Parser(locals()) as p:
    p.str('snp_matrix_file') # the matrix of snps, 398 sample rows by 10,000 (roughly) SNP columns
    p.str('feature_tsvfile')
    p.str('feature_tag')
    p.multiword('feature_tags').cast(lambda x:x.split())
    p.str('hipid_tag').default('hipid')
    p.str('grep')
    p.int('min_maf').default(3)
    p.int('force_snp_offset')
    p.float('pval_threshold').default(1e-3)
    p.flag('print_values')
    p.flag('snp_correlations')
    p.flag('setup')       # --flag_arg  (no argument passed)


if setup:

    exe = '/home/pbradley/code/tcr_scripts/read_gwas_data_genome_wide_linregress.py'
    assert exists(exe)

    xargs = '' # initialize

    # setup for running on the cluster
    runtag = 'run1'
    features_file = basedir+'tcrb_features_2020-02-22.tsv'
    feature_tags = 'n_inserts mhci_frac mait_frac'.split()
    hipid_tag = 'localID'

    runtag = 'run2'
    features_file = basedir+'emerson_vgene_out_frame_frequencies.tsv'
    feature_tags = 'TCRBV04-03 TCRBV19-01 total_tcrs'.split()
    hipid_tag = 'hipid'

    if 0:
        #runtag = 'run3'
        #runtag = 'run32'
        runtag = 'run33'

        features_file = basedir+'emerson_Jgene_frequencies_for_GWAS.tsv'
        lines = parse_tsv.parse_tsv_file(features_file)
        tags = lines[0].keys()
        max_badcount = 100
        feature_tags = []
        for qtag in [ x for x in tags if x.endswith('_Q') ]:
            vals = [ x[qtag] for x in lines]
            badcount = vals.count('None')
            if badcount <= max_badcount:
                feature_tags.extend( [ qtag, qtag[:-2]+'_IF', qtag[:-2]+'_OF'] )
                print 'good gene:', badcount, qtag[:-2]
            else:
                print 'bad  gene:', badcount, qtag[:-2]
        for f in feature_tags:
            assert f in tags
        hipid_tag = 'hipid'
        assert hipid_tag in tags

    else:

        runtag = 'run4'
        features_file = basedir+'emerson_Jgene_frequencies_for_GWAS.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run5'
        features_file = basedir+'trim_PCs.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run6'
        features_file = basedir+'trim_PCs_OF.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run7'
        features_file = basedir+'tmp.total_trims.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run8'
        features_file = basedir+'tmp.total_inserts.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run9' # clonality stats and also age, cmv
        features_file = basedir+'tmp.emerson_clonality_stats.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run10' # clonality stats and also age, cmv NOW LINREGRESS
        features_file = basedir+'tmp.emerson_clonality_stats.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        ## above this is old/not the latest min_mwu_pval setup ########################################
        runtag = 'run11' # clonality stats and also age, cmv NOW LINREGRESS and None for bad f2_fracs
        features_file = basedir+'tmp.emerson_clonality_stats.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run12' # cdr3 vdjtools features
        features_file = basedir+'tmp.all_cdr3_features_top10_v1.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run13' # cdr3 PC of aa frequencies
        features_file = basedir+'tmp.cdr3_aa_pcs.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run14' # fraction hla associated tcrs
        features_file = basedir+'tmp.hla_assoc_fraction.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run15' # fraction cmv associated tcrs
        features_file = basedir+'tmp.cmv_assoc_fraction.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run16'
        features_file = basedir+'tcrb_features_2020-02-22.tsv'
        feature_tags = 'n_inserts mhci_frac mait_frac'.split()
        hipid_tag = 'localID'

        runtag = 'run17' # cdr3 raw aa frequencies
        features_file = basedir+'tmp.cdr3_aa_raw_freqs.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run18' # cdr3 raw aa frequencies
        features_file = basedir+'tmp.cdr3_aa_raw_freqs_linreg.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run19' # cdr3 raw aa frequencies
        features_file = basedir+'tmp.all_cdr3_features_top10_v1_linreg.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run20' # cluster counts, normalized by subject bias
        features_file = basedir+'tmp.all_codist4dbs_cluster_counts.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run21' # cluster counts, normalized by subject bias -- only HLA+ subjects
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.all_codist4dbs_cluster_counts_hla_restricted.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run22' # cdr3 aa freqs and lenbin freqs in and out of frame and Qs
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.cdr3_aa_freqs_and_Qs_filt_250k.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run23'
        features_file = basedir+'tmp.total_trims.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run24'
        features_file = basedir+'tmp.total_inserts.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run25' # shannon, simpson etc from rcount and also linregressed
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.emerson_clonality_stats_rcount.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run26' # cdr3 aa freqs and lenbin freqs in and out of frame and Qs
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.cdr3_aa_freqs_lens_scores_and_Qs_filt_250k.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run27' # hla assoc fractions
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.hla_assoc_fraction_v2.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run28' # gene aa freqs split by expanded, and 'Q' for those
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.freqs_by_expanded.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run29' # pnuc fractions
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.all_pnuc_stats.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run30' # new way of getting aa freqs, all frames
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.cdr3_aa_freqs_new.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]

        runtag = 'run31' # new way of getting aa freqs, all frames
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.cdr3_aa_freqs_new2.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]
        xargs = ' --grep new_ '

        runtag = 'run34' # d gene freqs
        features_file = '/home/pbradley/csdat/martin_gwas/tmp.emerson_d_freqs.tsv'
        hipid_tag = 'hipid'
        feature_tags = [ALL_FEATURE_TAGS]
        xargs = ' '


    #################################################################
    assert exists(features_file)
    rundir = basedir+'slurm/'+runtag+'/'
    mkdir(rundir)
    cmdsfile = '{}{}_commands.txt'.format(rundir, runtag)
    out = open(cmdsfile, 'w')

    lines = parse_tsv.parse_tsv_file(features_file) ## sanity check
    for ftag in feature_tags+[hipid_tag]:
        assert ftag == ALL_FEATURE_TAGS or ftag in lines[0].keys()

    snp_matrix_files = glob('{}snp_matrix_files/snp_matrix*csv'.format(basedir))
    for ii,matfile in enumerate(snp_matrix_files):
        if ii%250:
            print ii, matfile
        matdir = rundir+matfile.split('/')[-1]+'_run/'
        mkdir(matdir)
        for ftag in feature_tags:
            outfile = '{}{}_{}.txt'.format(matdir, matfile.split('/')[-1], ftag)

            cmd = 'python {} {} --snp_matrix_file {} --feature_tsvfile {}  --feature_tag {} --hipid_tag {} --pval_threshold 1e-3 > {} 2> {}.err'\
                  .format( exe, xargs, matfile, features_file, ftag, hipid_tag, outfile, outfile )
            out.write(cmd+'\n')
    out.close()
    print 'made:', cmdsfile

    exit()


if feature_tags is None:
    if feature_tag is None:
        feature_tags = [ALL_FEATURE_TAGS]
    else:
        feature_tags = [feature_tag]
else:
    assert feature_tag is None

hip2upn, upn2hip = virscan_util.load_hip2upn()

all_hips = sorted( hip2upn.keys())
assert len(all_hips) == 666

###########################################################################
# read the mapping from GWAS "scan" id to upn/hip
sample_file = basedir+'sample_annotations_31OctSep2018.txt'
header = None
id2hip = {}
hip2id = {}

for line in open(sample_file,'rU'):
    l = line[:-1].split('\t')
    if not header:
        header = l
    else:
        id = l[ header.index('scanID') ]
        upn = int( l[ header.index('upn') ] ) # virscan_util uses ints
        assert l[ header.index('pt.dnr') ] == 'DNR'
        hip = l[ header.index('localID') ]

        hip2 = upn2hip.get(upn,None)
        upn2 = hip2upn.get(hip,None)

        if hip != hip2:
            #print 'upn2hip:',upn,hip2,hip
            assert hip.replace('R','P') == hip2
        if upn != upn2:
            pass
            #print 'hip2upn:',hip,upn2,upn

        id2hip[id] = hip2
        assert hip2 not in hip2id
        hip2id[hip2] = id


# read the scanids for the 398 samples
gwas_hips = []
scanid_file = basedir+'big_sample_data.csv'
assert exists(scanid_file)
for line in open(scanid_file,'r'):
    l = line[:-1].split(',')
    if l[1] == '"x"':
        continue
    else:
        scanid = l[1]
        gwas_hips.append( id2hip[scanid])
gwas_hips = np.array(gwas_hips)
print 'num_gwas_hips:', len(gwas_hips)

# read the snp matrix file
## read snp genotypes
## this file is basically a matrix with nSNPS+1 columns and nHIPs+1 rows, where 1st row and column are headers
num_snps = 0
genotype_matrix = []
for line in open(snp_matrix_file,'r'):
    l = line[:-1].split(',')
    if not num_snps:
        num_snps = len(l)-1
    else:
        row = [int(x)  for x in l[1:]]
        row = [-1 if x==3 else x for x in row] # switch the NA genotype to -1
        assert len(row) == num_snps
        genotype_matrix.append(row)
        assert len(genotype_matrix) == int(l[0][1:-1])
genotype_matrix = np.array(genotype_matrix).transpose()

assert genotype_matrix.shape == ( num_snps, len(gwas_hips))
print 'num_snps:', num_snps

feature_tsvlines = parse_tsv.parse_tsv_file(feature_tsvfile)

if feature_tags == [ALL_FEATURE_TAGS]:
    feature_tags = [ x for x in feature_tsvlines[0] if x != hipid_tag ]
    if grep:
        feature_tags = [ x for x in feature_tags if grep in x ]
    print ALL_FEATURE_TAGS, '=', feature_tags

good_snps = set()

for feature_tag in feature_tags:

    # read the feature for correlation with snps
    hip2feature = {}
    for l in feature_tsvlines:
        hipid = l[hipid_tag]
        val = l[feature_tag]
        hip2feature[ hipid ] = np.nan if val=='None' else float(val)

    feature_vals = np.array( [ hip2feature.get(x, np.nan) for x in gwas_hips ] )
    feature_valid = ~np.isnan( feature_vals )

    assert feature_vals.shape[0] == genotype_matrix.shape[1]

    all_snp_valid_masks = []
    all_snps = []

    for ii in range(num_snps):
        is_valid = ( genotype_matrix[ii,:] != -1 ) & feature_valid
        num_valid = np.sum(is_valid)
        if num_valid==0:
            continue
        counts = Counter( genotype_matrix[ii,:][is_valid] )
        top_count = counts.most_common(1)[0][1]
        if num_valid - top_count >= min_maf:
            all_snps.append( ii )
            all_snp_valid_masks.append(is_valid)

    print 'num_valid_snps:', len(all_snps), feature_tag


    for isnp, snp in enumerate(all_snps):
        if force_snp_offset is not None and force_snp_offset != snp:
            continue
        snp_valid_mask = all_snp_valid_masks[isnp]

        gts = genotype_matrix[snp,:][snp_valid_mask].astype(float)
        vals = feature_vals[snp_valid_mask]
        assert gts.shape == vals.shape

        reg = linregress( gts, vals )
        if reg.pvalue < pval_threshold:
            good_snps.add(snp)
            gt_ints = genotype_matrix[snp,:][snp_valid_mask]
            gt_counts = Counter(gt_ints)
            if np.sum(gt_ints==0) >= min_maf:
                _,p0 = mannwhitneyu( vals[ gt_ints==0 ], vals[ gt_ints != 0 ] )
            else:
                p0 = 1.
            if np.sum(gt_ints==2) >= min_maf:
                _,p2 = mannwhitneyu( vals[ gt_ints==2 ], vals[ gt_ints != 2 ] )
            else:
                p2 = 1.
            outline = 'pval: {:9.2e} R: {:7.3f} slope: {:7.3f} min_mwu_pval: {:9.2e} gt_counts {:3d} {:3d} {:3d} {:3d} {} {}'\
                .format( reg.pvalue, reg.rvalue, reg.slope, min(p0, p2), gt_counts[0], gt_counts[1], gt_counts[2],
                         sum(gt_counts.values()), snp, feature_tag )
            print outline

            if print_values:
                print 'gt_means: {:.6f} {:.6f} {:.6f} gt_counts {:3d} {:3d} {:3d} {:3d} {} {}'\
                    .format( np.mean(vals[ gt_ints==0 ] ), np.mean(vals[ gt_ints==1 ] ), np.mean(vals[ gt_ints==2 ] ),
                             gt_counts[0], gt_counts[1], gt_counts[2], sum(gt_counts.values()), snp, feature_tag )

            sys.stdout.flush()

if snp_correlations:
    good_snps = sorted(good_snps)
    for snp1 in good_snps:
        gts1 = genotype_matrix[snp1,:]
        for snp2 in good_snps:
            gts2 = genotype_matrix[snp2,:]
            mask = (gts1 != -1) & (gts2 != -1)
            reg1 = linregress(gts1[mask], gts2[mask])
            reg2 = linregress(2-gts1[mask], gts2[mask])
            r1, r2 = reg1.rvalue, reg2.rvalue
            if r1<0:
                assert r2>0
                rval = r2
                gts2 = 2-gts2
            else:
                assert r2<0
                rval = r1
            counts = Counter( zip( gts1[mask], gts2[mask] ) )

            print 'corr {:4d} {:4d} {:9.3f} {}'\
                .format(snp1, snp2, rval, ' '.join('{:4d}'.format(counts[(i,j)]) for i in range(3) for j in range(3)))


print 'DONE'



