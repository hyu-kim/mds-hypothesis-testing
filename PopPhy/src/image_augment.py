import math
import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow_datasets as tfds

from tensorflow import keras
from tensorflow.keras import layers

from utils.popphy_io import get_config


# Distorts the color distibutions of images
class RandomColorAffine(layers.Layer):
    def __init__(self, brightness=0, jitter=0, **kwargs):
        super().__init__(**kwargs)

        self.brightness = brightness
        self.jitter = jitter

    def get_config(self):
        config = super().get_config()
        config.update({"brightness": self.brightness, "jitter": self.jitter})
        return config

    def call(self, images, training=True):
        # print(images.shape)
        if training:
            batch_size = tf.shape(images)[0]
            
            # brightness
            brightness_scales = tf.random.uniform(
                (batch_size, 1, 1, 1), minval=-self.brightness, maxval=self.brightness
            )
            images = tf.add(images, brightness_scales)
            
            # contrast
            contrast_factor = tf.random.uniform(
                (batch_size, 1, 1, 1), minval=-self.jitter, maxval=self.jitter
            )
            images_mean = tf.reduce_mean(images, axis=[1, 2, 3], keepdims=True)
            images = tf.add(images, contrast_factor * (images - images_mean))
            
            images = tf.clip_by_value(images, 0, 1)
            
        return images

# Image augmentation module
def get_augmenter(brightness, jitter):
    config = get_config()
    image_width = config.getint('SimCLR', 'image_width')
    image_height = config.getint('SimCLR', 'image_height')
    image_channels = config.getint('SimCLR', 'image_channels')
    return keras.Sequential(
        [
            keras.Input(shape=(image_height, image_width, image_channels)),
            RandomColorAffine(brightness, jitter),
        ]
    )
