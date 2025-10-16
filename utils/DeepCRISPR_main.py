import os

import tensorflow as tf
import sonnet as snt
#import tensorflow.contrib.slim as slim
import tf_slim as slim

import pandas as pd
import numpy as np
import warnings 
warnings.filterwarnings("ignore")

import DeepCRISPR_Seq as dc
from DeepCRISPR_Seq import DCModelOntar

def process(seq):
    dna_dict={
        'A':[1,0,0,0],
        'C':[0,1,0,0],
        'G':[0,0,1,0],
        'T':[0,0,0,1],
        'N':[0,0,0,0]
    }

    def onehot(x,expand=1):
        x = np.stack( [np.vstack([ dna_dict[b] for b in seq ]) for seq in x  ] )
        x = x.transpose([0,2,1])
        if not expand is None:
            x = np.expand_dims(x, axis=expand) 
        return x
    
    return onehot(list(seq),expand=2)

def process_2(seq,expand=1):
    dna_dict={
        'A':[1,0,0,0],
        'C':[0,1,0,0],
        'G':[0,0,1,0],
        'T':[0,0,0,1],
        'N':[0,0,0,0]
    }

    def onehot(x,expand=1):
        x = np.stack( [np.vstack([ dna_dict[b] for b in seq ]) for seq in x  ] )
        if not expand is None:
            x = np.expand_dims(x, axis=expand) 
        return x
    
    return onehot(list(seq),expand=expand)

def predict_sequence(sequences, model_path='/cluster2/huanglab/liquan/pycode/dual/20250306_demo/DeepCRISPR_Seq/'):
    tf.reset_default_graph()
    sess = tf.InteractiveSession()

    ### 规定模型输入路径
    on_target_model_dir = model_path
    ### 载入模型
    dcmodel = dc.DCModelOntar(sess, on_target_model_dir, is_reg=True, seq_feature_only=True)
    ### 预测
    sequences = process(sequences)
    score = dcmodel.ontar_predict(sequences)
    
    return score

def build_source_model(model_path='/cluster2/huanglab/liquan/pycode/dual/20250306_demo/DeepCRISPR_Seq/',return_session=False):
    tf.reset_default_graph()
    sess = tf.InteractiveSession()
    tf.disable_v2_behavior()

    ### 规定模型输入路径
    on_target_model_dir = model_path
    ### 载入模型
    dcmodel = dc.DCModelOntar(sess, on_target_model_dir, is_reg=True, seq_feature_only=True)
    
    if return_session:
        return dcmodel.sess
    else:
        return dcmodel