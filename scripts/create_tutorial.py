from tensorflow.keras import models
import numpy as np
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))

# Load the Keras model for the DVS gesture categorization benchmark
model = models.load_model(os.path.join(PROJECT_DIR, "etc", "dvs_challenge.h5"))
inputs = np.loadtxt(os.path.join(PROJECT_DIR, "etc", "inputs.csv"))
inputs = inputs[0:1024]

dense_layer = np.reshape(model.layers[5].get_weights()[0], (9, 9, 11, 11))
dense_layer = np.swapaxes(dense_layer, 0, 1)
dense_layer = np.swapaxes(dense_layer, 2, 1)
dense_layer = np.reshape(dense_layer, (891, 11), order='C')
print(dense_layer)
np.savez("dvs_challenge.npz",
         conv1=np.rint(model.layers[0].get_weights()[0] * 420.05236577257483),
         conv2=np.rint(model.layers[1].get_weights()[0] *  351.1046444780251),
         conv3=np.rint(model.layers[2].get_weights()[0] * 276.6147837631879),
         conv4=np.rint(model.layers[3].get_weights()[0] * 371.60317670987195),
         dense1=np.rint(dense_layer * 341.41679600239286),
         inputs=np.rint(inputs),
         thresholds=np.array((255, 420, 351, 276, 371, 341)))
print(model.summary())