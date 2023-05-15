import numpy as np
from ydai.src.utils.prepare_data import prepare_data
from ydai.src.utils.popphy_io import get_config


def get_data():
    config = get_config()
    dataset = config.get('Evaluation', 'DataSet')
    print("\nStarting PopPhy-CNN on %s..." % (dataset))
    path = "ydai/data/" + dataset

    my_maps, _, _, _, tree_features, labels, label_set, g, feature_df = prepare_data(path, config)
    
    num_class = len(np.unique(labels))
    if num_class == 2:
        metric = "AUC"
    else:
        metric = "MCC"

    seed = np.random.randint(100)
    np.random.seed(seed)
    np.random.shuffle(my_maps)
    np.random.seed(seed)
    np.random.shuffle(labels)

    n_values = np.max(labels) + 1
    labels_oh = np.eye(n_values)[labels]
        
    print("There are %d classes...%s" % (num_class, ", ".join(label_set)))
    print(my_maps.shape)


if __name__ == "__main__":
	get_data()