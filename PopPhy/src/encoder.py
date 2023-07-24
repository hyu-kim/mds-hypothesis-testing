import tensorflow as tf
from utils.popphy_io import get_config

def get_encoder():
    config = get_config()
    num_col = config.getint('SimCLR', 'image_width')
    num_row = config.getint('SimCLR', 'image_height')
    num_kernel = int(config.get('PopPhy', 'NumberKernel'))
    kernel_height = int(config.get('PopPhy', 'KernelHeight'))
    kernel_width = int(config.get('PopPhy', 'KernelWidth'))
    num_fc_nodes = int(config.get('PopPhy', 'NumFCNodes'))
    num_cnn_layers = int(config.get('PopPhy', 'NumConvLayers'))
    num_fc_layers = int(config.get('PopPhy', 'NumFCLayers'))
    lamb = float(config.get('PopPhy', 'L2Lambda'))
    drop = float(config.get('PopPhy', 'Dropout'))

    reg = tf.keras.regularizers.l2(lamb)
    model = tf.keras.Sequential()

    model.add(tf.keras.layers.GaussianNoise(0.01, input_shape=(num_row, num_col, 1)))

    for i in range(0, num_cnn_layers):
        model.add(tf.keras.layers.Conv2D(
            filters=num_kernel, 
            kernel_size=(kernel_height, kernel_width), 
            activation='relu', 
            bias_regularizer=reg, 
            kernel_regularizer=reg, 
            name="conv_"+str(i)
            ))
        # model.add(tf.keras.layers.MaxPooling2D(pool_size=2))

    model.add(tf.keras.layers.Flatten())
    # model.add(tf.keras.layers.Dropout(drop))		

    for i in range(0, num_fc_layers):
        model.add(tf.keras.layers.Dense(
            num_fc_nodes, 
            activation='relu', 
            kernel_regularizer=reg, 
            bias_regularizer=reg, 
            name="fc_"+str(i)
            ))
        # model.add(tf.keras.layers.Dropout(drop))

    # # Output layer
    # model.add(tf.keras.layers.Dense(
    #     num_class, 
    #     activation='softmax', 
    #     kernel_regularizer=reg, 
    #     bias_regularizer=reg, 
    #     name="output"))
    
    return model