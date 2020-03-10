import numpy as np
import tensorflow as tf

def weight_through_time(x, num_units):
    # shape of x: batch, time, num_units)
    s = tf.reduce_max(x, axis=1)
    s = tf.layers.Dense(num_units, use_bias = False)(s)
    s = tf.nn.softmax(s)
    s = tf.tile(s[:, tf.newaxis, :], [1, tf.shape(x)[1], 1])
    return tf.multiply(x, s)
    
# def CapsuleLayer(x, 

def positional_encoding(output_dim):
    position_enc = np.array([[pos / np.power(10000, 2. * i / output_dim)
                              for i in range(output_dim)]
                             for pos in range(200)])

    position_enc[:, 0::2] = np.sin(position_enc[:, 0::2])  # dim 2i
    position_enc[:, 1::2] = np.cos(position_enc[:, 1::2])  # dim 2i+1
    # return tf.Variable(tf.convert_to_tensor(position_enc, dtype=tf.float32), trainable=False, name="pos_encode")
    return tf.convert_to_tensor(position_enc, dtype=tf.float32)
    
def get_pos_embedding(position_code, time_step, batch_size):
    src_pos = tf.tile(tf.reshape(tf.range(time_step), [1, -1]), [batch_size, 1])
    pos_embedding = tf.nn.embedding_lookup(position_code, src_pos)
    return pos_embedding # batch, time_step, output_dim
    
def multihead_self_attention(input, output_dim):
    # input: batch_size, time_step, input_dim
    # https://zhuanlan.zhihu.com/p/47282410
    # https://github.com/tensorflow/tensor2tensor/blob/master/tensor2tensor/layers/common_attention.py
    q, k, v = compute_qkv(input, output_dim)
    x = dot_product_attention(q, k, v)
    x = tf.layers.Dense(output_dim, use_bias=False)(x)
    return x
    
def compute_qkv(input, output_dim):
    """Computes query, key and value.
    Args:
    input: a Tensor with shape [batch, length_q, input_shape]
    Returns:
    q, k, v : [batch, length, depth] tensors
    """
    q = compute_attention_component(input, output_dim)
    k = compute_attention_component(input, output_dim)
    v = compute_attention_component(input, output_dim)
    return q, k, v

def compute_attention_component(x, output_dim):
    """Computes attention compoenent (query, key or value).
    Args:
    x: a Tensor with shape [batch, time_step, input_dim]
    Returns:
    c : [batch, time_step, output_dim] tensor
    """
    return tf.layers.Dense(output_dim, use_bias=False)(x)
  
def dot_product_attention(q, k, v):
    """Dot-product attention.
    Args:
    q: Tensor with shape [batch, time_step, output_dim].
    k: Tensor with shape [batch, time_step, output_dim]. Leading dimensions must
      match with q.
    v: Tensor with shape [batch, time_step, output_dim] Leading dimensions must
      match with q.
    Returns:
    Tensor with shape [batch , time_step, output_dim].
    """
    logits = tf.matmul(q, k, transpose_b=True)
    weights = tf.nn.softmax(logits, name="attention_weights")
    return tf.matmul(weights, v)