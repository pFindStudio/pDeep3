import os
import time

import tensorflow as tf

from .bucket import *

np.random.seed(1337)  # for reproducibility
# tf.compat.v1.set_random_seed(1337)
tf.set_random_seed(1337)

# pdeep_lstm_cell = tf.keras.layers.LSTMCell
pdeep_lstm_cell = tf.nn.rnn_cell.LSTMCell


# pdeep_lstm_cell = tf.contrib.cudnn_rnn.CudnnCompatibleLSTMCell

class pDeepRTModel:
    def __init__(self, conf):
        self.batch_size = 1024
        self.layer_size = 256
        self.attention_size = 32
        self.epochs = 100
        self.config = conf
        self.sess = None
        self.learning_rate = 0.001
        self.num_threads = 2
        self.dropout = 0.2
        self.optim_name = "Adam"
        self.graph = tf.Graph()
        self.RT_norm = 3600

    def use_cuda(self, _use_cuda=True):
        pass

    def restart_session(self):
        self.close()
        self.init_session()

    def init_session(self):
        if self.sess is None:
            # config = tf.compat.v1.ConfigProto()
            config = tf.ConfigProto()
            config.gpu_options.allow_growth = True
            config.intra_op_parallelism_threads = self.num_threads
            # self.sess = tf.compat.v1.Session(config=config)
            self.sess = tf.Session(config=config, graph=self.graph)

    def close(self):
        if self.sess is not None:
            self.sess.close()
            self.sess = None

    def GetVariableList(self):
        # self.var_list = [v for v in tf.compat.v1.global_variables() if not self.optim_name in v.name and not "_power" in v.name]
        self.var_list = [v for v in tf.global_variables() if not self.optim_name in v.name and not "_power" in v.name]
        return self.var_list

    def GetTensorList(self):
        # self.ten_list = [n for n in tf.compat.v1.get_default_graph().as_graph_def().node]
        self.ten_list = [n for n in self.graph.as_graph_def().node]
        return self.ten_list

    def SaveVariableName(self, save_as):
        with open(save_as + ".varname", "w") as f:
            for var in self.var_list: f.write(var.name + "\n")

    def SaveTensorName(self, save_as):
        with open(save_as + ".tenname", "w") as f:
            for ten in self.ten_list: f.write(ten.name + "\n")

    def GetTensorByName(self, name):
        # return tf.compat.v1.get_default_graph().get_tensor_by_name(name)
        return self.graph.get_tensor_by_name(name)

    def GetVariableByName(self, name):
        # l = [v for v in tf.compat.v1.global_variables() if v.name == name]
        l = [v for v in tf.global_variables() if v.name == name]
        if len(l) == 0:
            return None
        else:
            return l[0]

    def GetVariablesBySubstr(self, substr):
        # return [v for v in tf.compat.v1.global_variables() if substr in v.name]
        return [v for v in tf.global_variables() if substr in v.name]

    def GetVariableListByScope(self, scope):
        return tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES, scope=scope)

    def BuildModel(self, aa_size, mod_size, output_size, nlayers=1):
        with self.graph.as_default():
            print("BuildModel ... ")

            def _input():
                self._aa_x = tf.placeholder("float", [None, None, aa_size], name="input_aa_x")
                self._mod_x = tf.placeholder("float", [None, None, mod_size], name="input_mod_x")
                self._y = tf.placeholder("float", [None, ], name="input_y")
                self._time_step = tf.placeholder(tf.int32, [None, ], name="input_time_step")
                self._dropout = tf.placeholder_with_default(0.0, shape=(), name="dropout")
                
            def _attention(x):
                attention_weight = tf.Variable(tf.random_uniform([self.config.time_step]), trainable=True, name="input_attention")
                attention_weight = tf.tile(attention_weight[tf.newaxis, tf.newaxis, ...,], [tf.shape(x)[0], 1, int(self.config.time_step/self.attention_len)+self.attention_len])
                attention_weight = tf.slice(attention_weight, [0, 0, 0], [tf.shape(x)[0], 1, self._time_step[0]])
                x = tf.matmul(attention_weight, x) 

            def MultiLayerRNN(x):
                def BiLSTM(x, id):
                    lstm_fw_cell = pdeep_lstm_cell(self.layer_size)
                    lstm_bw_cell = pdeep_lstm_cell(self.layer_size)
                    x, _ = tf.nn.bidirectional_dynamic_rnn(lstm_fw_cell, lstm_bw_cell, x, sequence_length=self._time_step,
                                                           time_major=False, dtype=tf.float32, scope="BiLSTM_%d" % id)
                    x = tf.concat(x, axis=2)
                    x = tf.nn.dropout(x, rate=self._dropout)
                    return x

                for id in range(nlayers):
                    x = BiLSTM(x, id)
                    x = tf.identity(x, name="BiLSTM_output_%d" % id)
                return x

            def _output(x):
                def _reduce_attention_dim(x):
                    attention_weight = tf.Variable(tf.random_uniform([self.attention_len]), trainable=True, name="attention_weight")
                    attention_weight = tf.tile(attention_weight[tf.newaxis, tf.newaxis, ...,], [tf.shape(x)[0], 1, int(self.config.time_step/self.attention_len)+self.attention_len])
                    attention_weight = tf.slice(attention_weight, [0, 0, 0], [tf.shape(x)[0], 1, self._time_step[0]])
                    x = tf.matmul(attention_weight, x) 
                    return tf.layers.Flatten()(x) #batch, self.layer_size*2
                def _reduce_time_step(x):
                    x1 = tf.reduce_max(x, axis=1)[..., tf.newaxis]
                    x2 = tf.reduce_min(x, axis=1)[..., tf.newaxis]
                    x3 = tf.reduce_mean(x, axis=1)[..., tf.newaxis]
                    w = tf.Variable(tf.random_uniform([3]), trainable=True, name="reduce_weight")
                    w = tf.tile(w[tf.newaxis, ..., tf.newaxis], [tf.shape(x)[0], 1, 1])
                    x = tf.matmul(tf.concat((x1, x2, x3), axis=2), w)
                    return tf.layers.Flatten()(x) #batch, self.layer_size*2
                x = tf.reduce_max(x, axis=1)
                # x = _reduce_time_step(x)
                # x = _reduce_attention_dim(x)
                with tf.variable_scope("output_nn"):
                    w = tf.Variable(tf.random_uniform([self.layer_size*2, 128]), trainable=True, name="FC1_w")
                    b = tf.Variable(tf.random_uniform([128]), trainable=True, name="FC1_b")
                    x = tf.matmul(x, w) + b
                    x = tf.contrib.layers.layer_norm(x)
                    x = tf.nn.dropout(x, rate=self._dropout)
                    
                    w = tf.Variable(tf.random_uniform([128, 64]), trainable=True, name="FC2_w")
                    b = tf.Variable(tf.random_uniform([64]), trainable=True, name="FC2_b")
                    x = tf.matmul(x, w) + b
                    x = tf.contrib.layers.layer_norm(x)
                    x = tf.nn.dropout(x, rate=self._dropout)
                    
                    w = tf.Variable(tf.random_uniform([64, 1]), trainable=True, name="FC3_w")
                    b = tf.Variable(tf.random_uniform([1]), trainable=True, name="FC3_b")
                    x = tf.matmul(x, w) + b
                    
                with tf.name_scope("output_scope"):
                    self._prediction = tf.identity(x, name="output")

            _input()
            x = tf.concat((self._aa_x, self._mod_x), axis=2)
            x = MultiLayerRNN(x)
            _output(x)

            self.GetVariableList()
            self.GetTensorList()

            self.restart_session()
            self.Optimizer()
            init_op = tf.global_variables_initializer()
            self.sess.run(init_op)

    def Optimizer(self):
        with self.graph.as_default():
            self._loss = tf.sqrt(tf.reduce_mean(tf.square(self._prediction - self._y)))

            self._optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate, name=self.optim_name)

            self._mininize = self._optimizer.minimize(self._loss)

    def BuildTransferModel(self, model_file):
        with self.graph.as_default():
            print("Fine-tuning pDeep model ...")
            self.LoadModel(model_file)

            transfer_vars = []
            transfer_vars += tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "output_nn")
            transfer_vars += tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "BiLSTM_0")
            # transfer_vars += tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, self.instrument_ce_scope)

            # print(transfer_vars)

            self._loss = tf.sqrt(tf.reduce_mean(tf.square(self._prediction - self._y)))

            self._optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate, name="transfer_" + self.optim_name)

            self._mininize = self._optimizer.minimize(self._loss, var_list=transfer_vars)

            # self.merged_summary_op = tf.summary.merge_all()

            adam_vars = self.GetVariablesBySubstr("_power") + self.GetVariablesBySubstr(self.optim_name)
            init_op = tf.variables_initializer(var_list=adam_vars, name="transfer_init")
            # init_op = tf.global_variables_initializer()
            self.sess.run(init_op)

    def TrainModel(self, buckets, save_as=None):
        with self.graph.as_default():
            if self.sess is None:
                print("[Error] no session for training tensorflow!")
                sys.exit(-1)

            bbatch = Bucket_Batch(buckets, batch_size=self.batch_size, shuffle=True)

            mean_costs = []

            for epoch in range(self.epochs):
                batch_cost = []
                batch_time_cost = []
                ith_batch = 0
                for batch in bbatch.generate_batch():
                    start_time = time.perf_counter()
                    ith_batch += 1
                    peplen = bbatch.get_data_from_batch(batch, "peplen")
                    x = np.float32(bbatch.get_data_from_batch(batch, "x"))
                    mod_x = np.float32(bbatch.get_data_from_batch(batch, "mod_x"))
                    y = bbatch.get_data_from_batch(batch, "y")/self.RT_norm

                    feed_dict = {
                        self._aa_x: x,
                        self._mod_x: mod_x,
                        self._time_step: peplen - 1,
                        self._y: y,
                        self._dropout: self.dropout
                    }

                    pred, cost, _ = self.sess.run([self._prediction, self._loss, self._mininize], feed_dict=feed_dict)
                    # print(pred)
                    end_time = time.perf_counter()
                    batch_time_cost.append(end_time - start_time)
                    batch_cost.append(cost)

                    print(
                        "Epoch = {:3d}, peplen = {:3d}, Batch={:5d}, size = {:4d}, cost = {:.5f}, time = {:.2f}s\r".format(
                            epoch + 1, peplen[0], ith_batch, len(x), cost, end_time - start_time), end="")
                mean_costs.append("Epoch = {:3d}, mean_cost = {:.5f}, time = {:.2f}s".format(epoch + 1, np.mean(batch_cost),
                                                                                             np.sum(batch_time_cost)))
                print("\n" + mean_costs[-1])

            print("")
            for l in mean_costs:
                print(l)

            if save_as is not None:
                self.SaveModel(save_as)

    def Predict(self, buckets):
        with self.graph.as_default():
            bbatch = Bucket_Batch(buckets, batch_size=self.batch_size, shuffle=False)
            output_buckets = {}
            for batch in bbatch.generate_batch():
                peplen = bbatch.get_data_from_batch(batch, "peplen")
                x = np.float32(bbatch.get_data_from_batch(batch, "x"))
                mod_x = np.float32(bbatch.get_data_from_batch(batch, "mod_x"))

                feed_dict = {
                    self._aa_x: x,
                    self._mod_x: mod_x,
                    self._time_step: peplen - 1,
                }
                predictions = self.sess.run(self._prediction, feed_dict=feed_dict)
                # print(predictions)
                _buckets = {peplen[0]: (predictions*self.RT_norm,)}
                output_buckets = merge_buckets(output_buckets, _buckets)
            return output_buckets

    def SaveModel(self, model_file):
        dir = os.path.dirname(model_file)
        if not os.path.exists(dir): os.makedirs(dir)
        # print("enter tf.train.Saver")
        saver = tf.train.Saver(var_list=self.var_list)
        # print("enter saver.save")
        save_path = saver.save(self.sess, model_file)
        print("Model save as %s" % save_path)
        self.SaveVariableName(model_file)
        self.SaveTensorName(model_file)

    # model
    def LoadModel(self, model_file):
        # run multiple instances:
        # g = tf.Graph()
        # with g.as_default():
        self.restart_session()
        with self.graph.as_default():
        # saver = tf.compat.v1.train.import_meta_graph(model_file+".meta")
            saver = tf.train.import_meta_graph(model_file + ".meta")
            saver.restore(self.sess, model_file)
            # graph = tf.compat.v1.get_default_graph()
            # graph = tf.get_default_graph()
            self.GetVariableList()
            self.GetTensorList()
            self._aa_x = self.graph.get_tensor_by_name("input_aa_x:0")
            self._mod_x = self.graph.get_tensor_by_name("input_mod_x:0")
            self._y = self.graph.get_tensor_by_name("input_y:0")
            self._time_step = self.graph.get_tensor_by_name("input_time_step:0")
            self._prediction = self.graph.get_tensor_by_name("output_scope/output:0")
            self._dropout = self.graph.get_tensor_by_name("dropout:0")
