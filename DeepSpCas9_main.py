import os, sys
from os.path import exists
from os import system

import tensorflow as tf
tf_version = tf.__version__
print(f"TensorFlow version: {tf_version}")

if tf_version.startswith('1.'):
    import tensorflow as tf
    print("Running in TensorFlow 1.x compatibility mode.")
elif tf_version.startswith('2.'):
    import tensorflow.compat.v1 as tf
    try:
        tf.compat.v1.enable_eager_execution()
        print("Eager Execution Enabled:", tf.executing_eagerly())
    except AttributeError:
        print("Eager Execution is already enabled in TensorFlow 2.x.")
else:
    raise RuntimeError("Unsupported TensorFlow version.")

import numpy as np
import scipy.misc
import scipy.stats
from scipy.stats import rankdata

np.set_printoptions(threshold=5)

best_model_path_list = ['/cluster2/huanglab/liquan/pycode/dual/20250306_demo/DeepCas9_Final']

import logging
# set log level as WARNING
#tf.logging.set_verbosity(tf.logging.WARN)

# Model
length = 30

class DeepCas9(object):
    def __init__(self, filter_size, filter_num, node_1 = 80, node_2 = 60, l_rate = 0.005):
        self.inputs         = tf.placeholder(tf.float32, [None, 1, length, 4])
        self.targets        = tf.placeholder(tf.float32, [None, 1])
        self.is_training    = tf.placeholder(tf.bool)
        
        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, pool_shape, name):
            # setup the filter input shape for tf.nn.conv_2d
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                              num_filters]

            # initialise weights and bias for the filter
            weights   = tf.Variable(tf.truncated_normal(conv_filt_shape, stddev=0.03),
                                              name=name+'_W')
            bias      = tf.Variable(tf.truncated_normal([num_filters]), name=name+'_b')

            # setup the convolutional layer operation
            out_layer = tf.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding='VALID')

            # add the bias
            out_layer += bias

            # apply a ReLU non-linear activation
            out_layer = tf.layers.dropout(tf.nn.relu(out_layer), 0.3, self.is_training)

            # now perform max pooling
            ksize     = [1, pool_shape[0], pool_shape[1], 1]
            strides   = [1, 1, 2, 1]
            out_layer = tf.nn.avg_pool(out_layer, ksize=ksize, strides=strides, 
                                       padding='SAME')
            return out_layer
        #def end: create_new_conv_layer

        L_pool_0 = create_new_conv_layer(self.inputs, 4, filter_num[0], [1, filter_size[0]], [1, 2], name='conv1')
        L_pool_1 = create_new_conv_layer(self.inputs, 4, filter_num[1], [1, filter_size[1]], [1, 2], name='conv2')
        L_pool_2 = create_new_conv_layer(self.inputs, 4, filter_num[2], [1, filter_size[2]], [1, 2], name='conv3')
        with tf.variable_scope('Fully_Connected_Layer1'):
            layer_node_0 = int((length-filter_size[0])/2)+1
            node_num_0   = layer_node_0*filter_num[0]
            layer_node_1 = int((length-filter_size[1])/2)+1
            node_num_1   = layer_node_1*filter_num[1]
            layer_node_2 = int((length-filter_size[2])/2)+1
            node_num_2   = layer_node_2*filter_num[2]
            L_flatten_0  = tf.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1  = tf.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2  = tf.reshape(L_pool_2, [-1, node_num_2])
            L_flatten    = tf.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')
            node_num     = node_num_0 + node_num_1 + node_num_2
            W_fcl1       = tf.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1       = tf.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.nn.bias_add(tf.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.layers.dropout(L_fcl1, 0.3, self.is_training)

        with tf.variable_scope('Fully_Connected_Layer2'):
            W_fcl2       = tf.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2       = tf.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre   = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2       = tf.nn.relu(L_fcl2_pre)
            L_fcl2_drop  = tf.layers.dropout(L_fcl2, 0.3, self.is_training)
            
        with tf.variable_scope('Output_Layer'):
            W_out        = tf.get_variable("W_out", shape=[node_2, 1])#, initializer=tf.contrib.layers.xavier_initializer())
            B_out        = tf.get_variable("B_out", shape=[1])#, initializer=tf.contrib.layers.xavier_initializer())
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss    = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer   = tf.train.AdamOptimizer(l_rate).minimize(self.obj_loss)
    #def end: def __init__
#class end: DeepCas9

def Model_Finaltest(sess, TEST_X, filter_size, filter_num, if3d, model, args, load_episode, best_model_path):
    test_batch      = 500
    test_spearman   = 0.0
    optimizer       = model.optimizer
    TEST_Z          = np.zeros((TEST_X.shape[0], 1), dtype=float)
    
    for i in range(int(np.ceil(float(TEST_X.shape[0])/float(test_batch)))):
        Dict = {model.inputs: TEST_X[i*test_batch:(i+1)*test_batch], model.is_training: False}
        TEST_Z[i*test_batch:(i+1)*test_batch] = sess.run([model.outputs], feed_dict=Dict)[0]
    
    #OUT = open("RANK_final_{}.txt".format(best_model_path.split('/')[1]), "a")
    OUT = open("RANK_final_{}.txt".format(best_model_path.split('.txt')[0]), "a")
    OUT.write("Testing final \n {} ".format(tuple(TEST_Z.reshape([np.shape(TEST_Z)[0]]))))
    OUT.write("\n")
    OUT.close()
    return
#def end: Model_Finaltest


def preprocess_seq(data):
    #print("Start preprocessing the sequence done 2d")
    length  = 30
    
    DATA_X = np.zeros((len(data),1,length,4), dtype=int)
    #print(np.shape(data), len(data), length)
    for l in range(len(data)):
        for i in range(length):

            try: data[l][i]
            except: print(data[l], i, length, len(data))

            if data[l][i]in "Aa":    DATA_X[l, 0, i, 0] = 1
            elif data[l][i] in "Cc": DATA_X[l, 0, i, 1] = 1
            elif data[l][i] in "Gg": DATA_X[l, 0, i, 2] = 1
            elif data[l][i] in "Tt": DATA_X[l, 0, i, 3] = 1
            elif data[l][i] in "Nn": DATA_X[l, 0, i, :] = 0
            else:
                print ("Non-ATGC character " + data[l])
                print (i)
                print (data[l][i])
                sys.exit()
        #loop end: i
    #loop end: l
    #print("Preprocessing the sequence done")
    return DATA_X
#def end: preprocess_seq


def getseq(filenum):
    param   = parameters['%s' % filenum]
    FILE    = open(path+param, "r")
    data    = FILE.readlines()
    data_n  = len(data) - 1
    seq     = []

    for l in range(1, data_n+1):
        try:
            data_split = data[l].split()
            seq.append(data_split[1])
        except:
            print (data[l])
            seq.append(data[l])
    #loop end: l
    FILE.close()
    processed_full_seq = preprocess_seq(seq)

    return processed_full_seq, seq   


def predict_sequence(sequences, model_path='./DeepCas9_Final/'):
    
    # Input validation
    if isinstance(sequences, str):
        sequences = [sequences]
    
    for seq in sequences:
        if len(seq) != 30:
            raise ValueError("All sequences must be 30bp long. Found sequence of length {}".format(len(seq)))
        if not all(n in 'ATCGNatcgn' for n in seq):
            raise ValueError("Sequences can only contain A,T,C,G. Found invalid sequence: {}".format(seq))

    # Preprocess sequences
    processed_seqs = preprocess_seq(sequences)
    
    # Set up TensorFlow
    conf = tf.ConfigProto()
    conf.device_count['GPU'] = 0 # 禁用GPU防止崩溃
    #conf.gpu_options.allow_growth = True
    #conf.gpu_options.visible_device_list = "1"
    
    # Model parameters 
    filter_size = [3, 5, 7]
    filter_num = [100, 70, 40] 
    node_1 = 80
    node_2 = 60
    learning_rate = 0.005
    
    # Find latest model checkpoint
    model_files = [f for f in os.listdir(model_path) if f.endswith('.meta')]
    if not model_files:
        raise ValueError("No model checkpoint files found in {}".format(model_path))
    latest_model = model_files[0][:-5]  # Remove .meta extension
    
    # Initialize model and run predictions
    tf.reset_default_graph()
    with tf.Session(config=conf) as sess:
        sess.run(tf.global_variables_initializer())
        model = DeepCas9(filter_size, filter_num, node_1, node_2, learning_rate)
        
        saver = tf.train.Saver()
        saver.restore(sess, os.path.join(model_path, latest_model))
        
        # Run prediction
        feed_dict = {
            model.inputs: processed_seqs,
            model.is_training: False
        }
        predictions = sess.run(model.outputs, feed_dict=feed_dict)
    
    return predictions.flatten()


def return_session(filter_size= [3, 5, 7], filter_num= [100, 70, 40], node_1 = 80, node_2 = 60, learning_rate = 0.005,
                   model_path='./DeepCas9_Final/'):
    
    model_files = [f for f in os.listdir(model_path) if f.endswith('.meta')]
    if not model_files:
        raise ValueError("No model checkpoint files found in {}".format(model_path))
    latest_model = model_files[0][:-5]
    
    tf.reset_default_graph()
    sess = tf.Session()
    
    sess.run(tf.global_variables_initializer())
    model = DeepCas9(filter_size, filter_num, node_1, node_2, learning_rate)
    
    saver = tf.train.Saver()
    saver.restore(sess, os.path.join(model_path, latest_model))
    
    return sess


class DeepSpCas9_tf2(tf.keras.Model):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        super(DeepSpCas9_tf2, self).__init__()
        self.filter_size = filter_size
        self.filter_num = filter_num
        self.node_1 = node_1
        self.node_2 = node_2
        self.l_rate = l_rate
        self.length = 30

        # 定义卷积层
        self.conv1 = tf.keras.layers.Conv2D(
            filters=filter_num[0],
            kernel_size=(1, filter_size[0]),
            padding='valid',
            activation='relu'
        )
        self.pool1 = tf.keras.layers.AveragePooling2D(pool_size=(1, 2), strides=(1, 2), padding='same')
        self.dropout1 = tf.keras.layers.Dropout(0.3)

        self.conv2 = tf.keras.layers.Conv2D(
            filters=filter_num[1],
            kernel_size=(1, filter_size[1]),
            padding='valid',
            activation='relu'
        )
        self.pool2 = tf.keras.layers.AveragePooling2D(pool_size=(1, 2), strides=(1, 2), padding='same')
        self.dropout2 = tf.keras.layers.Dropout(0.3)

        self.conv3 = tf.keras.layers.Conv2D(
            filters=filter_num[2],
            kernel_size=(1, filter_size[2]),
            padding='valid',
            activation='relu'
        )
        self.pool3 = tf.keras.layers.AveragePooling2D(pool_size=(1, 2), strides=(1, 2), padding='same')
        self.dropout3 = tf.keras.layers.Dropout(0.3)

        # 计算全连接层的输入大小
        layer_node_0 = int((self.length - filter_size[0]) / 2) + 1
        node_num_0 = layer_node_0 * filter_num[0]
        layer_node_1 = int((self.length - filter_size[1]) / 2) + 1
        node_num_1 = layer_node_1 * filter_num[1]
        layer_node_2 = int((self.length - filter_size[2]) / 2) + 1
        node_num_2 = layer_node_2 * filter_num[2]
        self.node_num = node_num_0 + node_num_1 + node_num_2

        # 定义全连接层
        self.flatten = tf.keras.layers.Flatten()
        self.dense1 = tf.keras.layers.Dense(node_1, activation='relu')
        self.dropout4 = tf.keras.layers.Dropout(0.3)
        self.dense2 = tf.keras.layers.Dense(node_2, activation='relu')
        self.dropout5 = tf.keras.layers.Dropout(0.3)
        self.output_layer = tf.keras.layers.Dense(1)

        # 定义损失函数和优化器
        self.loss_fn = tf.keras.losses.MeanSquaredError()
        self.optimizer = tf.keras.optimizers.Adam(learning_rate=l_rate)

    def call(self, inputs, training=False):
        # 假设输入形状为 (batch_size, 1, length, 4)
        # 分别通过三个卷积层
        conv_out1 = self.conv1(inputs)
        pool_out1 = self.pool1(conv_out1)
        dropout_out1 = self.dropout1(pool_out1, training=training)

        conv_out2 = self.conv2(inputs)
        pool_out2 = self.pool2(conv_out2)
        dropout_out2 = self.dropout2(pool_out2, training=training)

        conv_out3 = self.conv3(inputs)
        pool_out3 = self.pool3(conv_out3)
        dropout_out3 = self.dropout3(pool_out3, training=training)

        # 扁平化并拼接
        flatten1 = self.flatten(dropout_out1)
        flatten2 = self.flatten(dropout_out2)
        flatten3 = self.flatten(dropout_out3)
        concat = tf.concat([flatten1, flatten2, flatten3], axis=1)

        # 全连接层
        dense1_out = self.dense1(concat)
        dropout4_out = self.dropout4(dense1_out, training=training)
        dense2_out = self.dense2(dropout4_out)
        dropout5_out = self.dropout5(dense2_out, training=training)
        outputs = self.output_layer(dropout5_out)
        
        return outputs


def build_source_model(model_path='./DeepCas9_Final/'):
    model = DeepSpCas9_tf2(filter_size=[3, 5, 7], filter_num=[100, 70, 40])

    if not model_path is None:
        model.load_weights(model_path)
    else:
        raise ValueError("No model checkpoint files found")
        
    return model