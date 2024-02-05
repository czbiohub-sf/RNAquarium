import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import RANSACRegressor
from sys import argv


# plot_combined_trace.py combined_trace.txt

full = pd.read_csv(argv[0], header=0, delimiter="\t")
processes = ['prefetch', 'fastq_dump', 'filter_barcodes', 'fastp', 'priceseqfilter',
             'hisat2', 'star_counts', 'sort_bam', 'htseq_count', 'star', 'bowtie2',
             'process_bowtie2_sam', 'bowtie2_filter_by_names', 'dedup', 'gsnap',
             'process_gsnap_sam', 'stats_csv', 'gsnap_filter_by_names']
levels, categories = pd.factorize(full['NodeType'])
nodetype_handles = [matplotlib.patches.Patch(color=plt.cm.tab10(i), label=c)
             for i, c in enumerate(categories)]
for proc in processes:
    x_vars = ["reads","size","median"]
    y_vars = ['ElapsedRaw','CPUTimeRAW','%cpu','peak_rss','MaxRSS','MaxVMSize']
    fig, ax = plt.subplots(len(y_vars), len(x_vars), figsize=(12,22), sharey='row')
    ax[0, len(x_vars)-1].legend(handles=nodetype_handles, title='node type')
    ax[0, 1].set_title(proc, fontdict={'fontsize': 16})
    plt.style.context('tableau-colorblind10')
    data = full[full['process']==proc]
    for col, y_var in enumerate(y_vars):
        for row, x_var in enumerate(x_vars):
            x = data[data[y_var].notna()][x_var]
            X = np.array(x)[:,np.newaxis]
            y = data[data[y_var].notna()][y_var]
            levels, categories = pd.factorize(data[data[y_var].notna()]['NodeType'])
            colors = [plt.cm.tab10(i) for i in levels]

            estimator = RANSACRegressor(random_state=121)
            estimator.fit(X, y)
            lin_pred = estimator.predict(np.reshape(X, (-1, 1)))
            residuals = y - lin_pred
            max_residual = residuals.idxmax()
            coef = estimator.estimator_.coef_
            ax[col, row].plot(x, lin_pred, linewidth=1.5, c=plt.cm.tab10(len(levels)))
            ax[col, row].axline(xy1=(x[max_residual], y[max_residual]),
                                slope=coef, lw=1.5, c=plt.cm.tab10(len(levels)))

            ax[col, row].scatter(x, y, s=5, marker="o", c=colors)

            ax[col, row].set_xlim(0, x.max())
            #ax[col, row].set_ylim(0, y.max())
            ax[col, row].set_ylim(auto=True)

            # draw text, get bbox, erase text, draw text in better location
            ann = ax[col, row].annotate(f"slope: {coef}\ny-int: {estimator.predict([[0]])+residuals[max_residual]}\nscore: {estimator.score(X,y):.3f}",
                                  (x[max_residual], y[max_residual]),
                                  horizontalalignment='center',
                                  verticalalignment='top')
            bb = ann.get_window_extent()
            bb_datacoords = bb.transformed(ax[col, row].transData.inverted())
            (cur_x, cur_y) = ann.xy
            if bb_datacoords.x0 < ax[col, row].get_xlim()[0]:
                new_x = cur_x + (ax[col, row].get_xlim()[0] - bb_datacoords.x0)
            elif bb_datacoords.x1 > ax[col, row].get_xlim()[1]:
                new_x = cur_x - (bb_datacoords.x1 - ax[col, row].get_xlim()[1])
            else:
                new_x = cur_x

            if bb_datacoords.y0 > ax[col, row].get_ylim()[1]:
                new_y = cur_y - (bb_datacoords.y0 - ax[col, row].get_ylim()[1])
            else:
                new_y = cur_y
            ann.remove()
            ax[col, row].annotate(f"slope: {coef}\ny-int: {estimator.predict([[0]])+residuals[max_residual]}\nscore: {estimator.score(X,y):.3f}",
                                  (x[max_residual], y[max_residual]),
                                  xytext=(new_x, new_y),
                                  horizontalalignment='center',
                                  verticalalignment='top')
            
            ax[col, row].set_xlabel(x_var)
            ax[col, row].set_ylabel(y_var)
            ax[col, row].grid(True)
    fig.tight_layout()
    plt.savefig(f"reg_{proc}_2.png")
