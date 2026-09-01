from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import run_cgs_models_mallet
cistopic_obj = create_cistopic_object(fragment_matrix=atac_metacell.X.T.tocsr(),
                                     cell_names = atac_metacell.obs_names.to_list(),
                                      region_names= atac_metacell.var_names.to_list(),
                                      tag_cells = False,
                                      split_pattern = ' ')
cistopic_obj.add_cell_data(atac_metacell.obs)
# pickle.dump(cistopic_obj,open('scenic/L2_level/scATAC_Peak_Metacell_cisTopic_obj.pkl', 'wb'))
cistopic_obj = pickle.load(open('scenic/L2_level/scATAC_Peak_Metacell_cisTopic_obj.pkl','rb'))
os.environ['MALLET_MEMORY'] = '250G'
!export _JAVA_OPTIONS=-Xmx250g
models=run_cgs_models_mallet(mallet_path='anaconda3/envs/scenicplus/bin/mallet',
                    cistopic_obj=cistopic_obj,
                    n_topics=[5, 20, 40, 60, 80, 100, 130, 150,180],
                    n_cpu=20,
                    n_iter=150,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    tmp_path='scenic/pycistopic/', #Use SCRATCH if many models or big data set
                    save_path='scenic/L2_level/topics/',
                    reuse_corpus=True)
pickle.dump(models,open(os.path.join('scenic/L2_level/', "models.pkl"), "wb"))
#
models = pickle.load(open('ATAC/data_process_all/scenic/pycistopic/models.pkl',"rb"))
import glob
topic_model_files = glob.glob(os.path.join('ATAC/data_process_all/scenic/L2_level/topics/', '*pkl'))
models = []
for file in topic_model_files:
    topic = pickle.load(open(file, 'rb'))
    models.append(topic)
    print(file)
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 40,
    return_model = True
)
cistopic_obj.add_LDA_model(model)
# pickle.dump(cistopic_obj,open('scenic/L2_level/scATAC_Peak_Metacell_cisTopic_obj.pkl', 'wb'))