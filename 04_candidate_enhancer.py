##### Inferring candidate enhancer regions #####
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img
from pycisTopic.utils import region_names_to_coordinates
import pickle
# topic-region (binarization of region-topic probabilites.)
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu',plot=True, num_columns=5)
region_bin_topics_top_3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3_000)
os.makedirs(os.path.join("scenic/L2_level/pycistopic/", "Topics_otsu"), exist_ok = True)
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(region_bin_topics_otsu[topic].index).sort_values(["Chromosome", "Start", "End"]
    ).to_csv(os.path.join("scenic/L2_level/pycistopic/", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",header = False, index = False)
os.makedirs(os.path.join("scenic/L2_level/pycistopic", "Topics_top_3k"), exist_ok = True)
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(region_bin_topics_top_3k[topic].index).sort_values(["Chromosome", "Start", "End"]
    ).to_csv(os.path.join("scenic/L2_level/pycistopic/", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",header = False, index = False)
# cell-topic
binarized_cell_topic = binarize_topics(cistopic_obj,target='cell',method='li',plot=True,num_columns=5, nbins=100)
# pickle.dump(binarized_cell_topic, open(os.path.join('scenic/pycistopic/binarized_cell_topic.pkl'), 'wb'))
topic_annot = topic_annotation(cistopic_obj,annot_var='L2_annotation',binarized_cell_topic=binarized_cell_topic,general_topic_thr = 0.2)
topic_qc_metrics = compute_topic_metrics(cistopic_obj)
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)
# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.show()
# Differentially Accessible Regions (DARs)
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
# variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj,
    min_disp = 0.05,min_mean = 0.0125,max_mean = 3,max_disp = np.inf,n_bins=20,n_top_features=None,plot=True)
markers_dict= find_diff_features(cistopic_obj,imputed_acc_obj,variable='L2_annotation',var_features=variable_regions,
    contrasts=None,adjpval_thr=0.05,log2fc_thr=np.log2(1.5),n_cpu=10,split_pattern = '-'
    )
print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")
# save
pickle.dump(markers_dict, open(os.path.join('scenic/L2_level/pycistopic/markers_dict.pkl'), 'wb'))
os.makedirs(os.path.join("scenic/L2_level/pycistopic", "DARs_celltype"), exist_ok = True)
for cell_type in markers_dict:
    region_names_to_coordinates(markers_dict[cell_type].index
    ).sort_values(["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join('scenic/L2_level/pycistopic/DARs_celltype/', f"{cell_type}.bed"),
        sep = "\t",header = False, index = False)