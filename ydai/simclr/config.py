########################################################################
# Data settings
########################################################################

[Evaluation]

# NumberTestSplits	Number of partitions (k) in k-fold cross-validation [integer]
# NumberRuns		Number of iterations to run cross-validation [integer]
# Dataset		Dataset located in ../data
# FilterThresh		Filter OTUs based on porportion found in samples [float]

NumberTestSplits = 6
NumberRuns = 10
DataSet = Alga
FilterThresh = 0.1


########################################################################
# PopPhy-CNN models
########################################################################

[PopPhy]

# LearningRate		Learning rate for PopPhy-CNN models [float]
# BatchSize		Batch size for PopPhy-CNN models [integer]

LearningRate = 0.001
BatchSize = 128
Patience = 40
NumberKernel = 32
KernelWidth = 5
KernelHeight = 3
NumFCNodes = 32
NumConvLayers = 1
NumFCLayers = 1
L2Lambda = 0.001
Dropout = 0.3


########################################################################
# SimCLR settings
########################################################################

[SimCLR]

# Dataset hyperparameters
unlabeled_dataset_size = 100000
labeled_dataset_size = 5000
image_size = 96
image_channels = 3

# Algorithm hyperparameters
num_epochs = 20
batch_size = 525  # Corresponds to 200 steps per epoch
width = 128
temperature = 0.1
# Stronger augmentations for contrastive, weaker ones for supervised training
contrastive_augmentation = {"brightness": 0.6, "jitter": 0.2}
classification_augmentation = {"brightness": 0.3, "jitter": 0.1}