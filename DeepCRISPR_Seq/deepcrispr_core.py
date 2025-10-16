import tensorflow as tf
from tensorflow.keras.layers import Conv2D, BatchNormalization, ReLU, Softmax
from tensorflow.keras import Model, Input

class DeepCRISPR_tf2(tf.keras.Model):
    def __init__(self,channel_size=[8, 32, 64, 64, 256, 256],cls_channel_size = [512, 512, 1024, 1]):
        super(DeepCRISPR_tf2, self).__init__()
        
        self.channel_size = channel_size
        self.cls_channel_size = cls_channel_size
        self.betas = [None] + [tf.Variable(tf.zeros(channel_size[i]), name=f'beta_{i}') for i in range(1, len(channel_size))]
        
        # Encoder layers
        self.encoder = [
            None,
            Conv2D(self.channel_size[1], kernel_size=(1, 3), padding='same', name='e_1'),
            Conv2D(self.channel_size[2], kernel_size=(1, 3), padding='same', strides=2, name='e_2'),
            Conv2D(self.channel_size[3], kernel_size=(1, 3), padding='same', name='e_3'),
            Conv2D(self.channel_size[4], kernel_size=(1, 3), padding='same', strides=2, name='e_4'),
            Conv2D(self.channel_size[5], kernel_size=(1, 3), padding='same', name='e_5')
        ]

        self.encoder_bn = [
            None,
            BatchNormalization(name='ebn_1u'),
            BatchNormalization(name='ebn_2u'),
            BatchNormalization(name='ebn_3u'),
            BatchNormalization(name='ebn_4u'),
            BatchNormalization(name='ebn_5u')
        ]
        
        # Classifier 
        self.classifier = [
            None,
            Conv2D(self.cls_channel_size[0], kernel_size=(1, 3), padding='same', strides=2, name='e_6'),
            Conv2D(self.cls_channel_size[1], kernel_size=(1, 3), padding='same', name='e_7'),
            Conv2D(self.cls_channel_size[2], kernel_size=(1, 3), padding='valid', name='e_8'),
            Conv2D(self.cls_channel_size[3], kernel_size=(1, 1), padding='valid', name='e_9')
        ]

        self.classifier_bn = [
            None,
            BatchNormalization(name='ebn_6l'),
            BatchNormalization(name='ebn_7l'),
            BatchNormalization(name='ebn_8l')
        ]
        
        self.build( (None,1,23,4) )
        
    def call(self,inputs_sg):
        
        # Encoder forward pass
        hu = inputs_sg
        for i in range(1, len(self.channel_size)):
            pre_u = self.encoder[i](hu)
            u = self.encoder_bn[i](pre_u, training=True)
            hu = ReLU()(u + self.betas[i])
            #print(f'Encoder {i}: ',hu.shape)
        
        # Classifier forward pass
        hl = hu
        for i in range(1, len(self.cls_channel_size)):
            pre_l = self.classifier[i](hl)
            l = self.classifier_bn[i](pre_l, training=True)
            hl = ReLU()(l)
            #print(f'Classifier {i}: ',hl.shape)

        l_last = self.classifier[-1](hl)
        #hl_last = Softmax(axis=-1)(l_last)
        #print(l_last.shape)

        logits_l = tf.squeeze(l_last, axis=[2, 3])
        #print(logits_l.shape)
        
        return logits_l
    
