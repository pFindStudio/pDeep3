import os
import time

import tensorflow as tf

from .bucket import *
from . import tf_ops
from .featurize import to_model_feature

np.random.seed(1337)  # for reproducibility
tf.compat.v1.set_random_seed(1337)

print("tensorflow version = %s"%tf.__version__)
if tf.__version__ < '2.0.0':
    use_tf2 = False
else:
    use_tf2 = True


class pDeepModel:
    def __init__(self, conf):
        self.batch_size = 1024
        self.layer_size = 256
        self.epochs = 100
        self.config = conf
        self.sess = None
        self.learning_rate = 0.001
        self.num_threads = 2
        self.dropout = 0.2
        self.optim_name = "Adam"
        self.graph = tf.Graph()

    def use_cuda(self, _use_cuda=True):
        pass

    def BuildModel(self, aa_size, mod_size, output_size, nlayers=2):
        if use_tf2: tf.compat.v1.disable_eager_execution()
        with self.graph.as_default():
            print("BuildModel ... ")

            def _input():
                self._aa_x = tf.compat.v1.placeholder("float", [None, None, aa_size], name="input_aa_x")
                self._mod_x = tf.compat.v1.placeholder("float", [None, None, mod_size], name="input_mod_x")
                self._y = tf.compat.v1.placeholder("float", [None, None, output_size], name="input_y")
                self._time_step = tf.compat.v1.placeholder(tf.int32, [None, ], name="input_time_step")
                self._charge = tf.compat.v1.placeholder("float", [None, ], name="input_charge")
                self._nce = tf.compat.v1.placeholder("float", [None, ], name="input_nce")
                self._instrument = tf.compat.v1.placeholder("float", [None, self.config.max_instrument_num], name="input_instrument")
                self._rnn_dropout = tf.compat.v1.placeholder_with_default(0.0, shape=(), name="rnn_dropout")
                self._dropout = tf.compat.v1.placeholder_with_default(0.0, shape=(), name="dropout")
                self.rnn_kp = None
                self.output_kp = None

            def ExpandTimeStep1D(x):
                return tf.tile(x[..., tf.newaxis, tf.newaxis], [1, self._time_step[0], 1])

            def ExpandTimeStep2D(x):
                return tf.tile(tf.expand_dims(x, 1), [1, self._time_step[0], 1])

            def LinearEmbed2D(x, in_size, out_size, scope, name):
                with tf.compat.v1.variable_scope(scope):
                    w = tf.Variable(tf.compat.v1.random_normal([in_size, out_size]), trainable=True, name=name)
                    outputs = tf.matmul(x, w)
                return outputs

            def LinearEmbed3D(x, in_size, out_size, scope, name):
                pass

            def ConcatTimeFeatures(*args):
                return tf.concat(args, axis=2)

            def MultiLayerRNN(x, ch, ins_nce):
                def BiLSTM(x, id):
                    def tf_v1(x, id):
                        lstm_fw_cell = tf.nn.rnn_cell.LSTMCell(self.layer_size)
                        lstm_bw_cell = tf.nn.rnn_cell.LSTMCell(self.layer_size)
                        # lstm_fw_cell = tf.contrib.rnn.DropoutWrapper(lstm_fw_cell, input_keep_prob=1, output_keep_prob=self.rnn_kp, state_keep_prob=self.rnn_kp, variational_recurrent=False, dtype=tf.float32)
                        # lstm_bw_cell = tf.contrib.rnn.DropoutWrapper(lstm_bw_cell, input_keep_prob=1, output_keep_prob=self.rnn_kp, state_keep_prob=self.rnn_kp, variational_recurrent=False, dtype=tf.float32)
                        x, _ = tf.nn.bidirectional_dynamic_rnn(lstm_fw_cell, lstm_bw_cell, x, sequence_length=self._time_step, time_major=False, dtype=tf.float32, scope="BiLSTM_%d" % id)
                        x = tf.concat(x, axis=2)
                        x = tf.nn.dropout(x, rate=self._dropout)
                        return x
                    def tf_v2(x, id):
                        with tf.compat.v1.variable_scope("BiLSTM_%d"%id):
                            rnn = tf.keras.layers.LSTM(self.layer_size, return_sequences=True)
                            x = tf.keras.layers.Bidirectional(rnn)(x)
                            x = tf.nn.dropout(x, rate=self._dropout)
                            return x
                    if use_tf2: return tf_v2(x, id)
                    else: return tf_v1(x, id)

                for id in range(nlayers):
                    x = BiLSTM(x, id)
                    x = tf.identity(x, name="BiLSTM_output_%d" % id)

                    x = ConcatTimeFeatures(x, ch, ins_nce)
                return x

            def _output(x):
                def OutputRNN(x):
                    def tf_v1(x):
                        cell_fw = tf.nn.rnn_cell.LSTMCell(output_size)
                        outputs, _ = tf.nn.dynamic_rnn(cell_fw, x, sequence_length=self._time_step, time_major=False, dtype=tf.float32, scope="output_nn")
                        return outputs
                    def tf_v2(x):
                        with tf.compat.v1.variable_scope("output_nn"):
                            outputs = tf.keras.layers.LSTM(output_size, return_sequences=True)(x)
                            return outputs
                    if use_tf2: return tf_v2(x)
                    else: return tf_v1(x)

                # def OutputBiRNN(x):
                    # cell_fw = tf.nn.rnn_cell.LSTMCell(output_size)
                    # cell_bw = tf.nn.rnn_cell.LSTMCell(output_size)
                    # outputs, _ = tf.nn.bidirectional_dynamic_rnn(cell_fw, cell_bw, x, sequence_length=self._time_step,
                                                                # time_major=False, dtype=tf.float32,
                                                                # scope="output_nn")  # this scope may have to be tuned in transfer learning
                    
                    # with tf.compat.v1.variable_scope("output_nn"):
                        # rnn = tf.keras.layers.LSTM(output_size, return_sequences=True)
                        # outputs = tf.keras.layers.Bidirectional(rnn)(x)
                    # outputs = tf.add(outputs[0], outputs[1])

                outputs = OutputRNN(x)
                with tf.name_scope("output_scope"):
                    self._prediction = tf.identity(outputs, name="output")

            _input()

            ins_dim = 3
            ins_nce = tf.concat((self._instrument, self._nce[..., tf.newaxis]), axis=1)
            ins_nce = LinearEmbed2D(ins_nce, self.config.max_instrument_num + 1, ins_dim, scope="ins_nce", name="ins_nce_w")

            ch = ExpandTimeStep1D(self._charge)
            ins_nce = ExpandTimeStep2D(ins_nce)

            x = ConcatTimeFeatures(self._aa_x, self._mod_x, ch, ins_nce)

            x = MultiLayerRNN(x, ch, ins_nce)
            # x = tf_ops.attention_through_time(x, self.layer_size*2+4)
            _output(x)
            
            print("[pDeep Info] # trainable parameters = %s"%tf_ops.count_parameters())

            self.GetVariableList()
            self.GetTensorList()

            self.restart_session()
            self.Optimizer()
            init_op = tf.compat.v1.global_variables_initializer()
            self.sess.run(init_op)

    def Optimizer(self):
        with self.graph.as_default():
            self._loss = tf.reduce_mean(tf.abs(self._prediction - self._y))

            self._optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=self.learning_rate, name=self.optim_name)

            self._mininize = self._optimizer.minimize(self._loss)

    def BuildTransferModel(self, model_file):
        with self.graph.as_default():
            print("Fine-tuning pDeep model ...")
            self.LoadModel(model_file)

            transfer_vars = []
            transfer_vars += tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES, "output_nn")
            transfer_vars += tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES, "BiLSTM_0")
            transfer_vars += tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES, "BiLSTM_1")
            # transfer_vars += tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES, self.instrument_ce_scope)

            # print(transfer_vars)

            self._loss = tf.reduce_mean(tf.abs(self._prediction - self._y))
            
            # RMSE:
            # RMSE prefers to predict budding peaks (very low inten peaks)
            # self._loss = tf.sqrt(tf.reduce_mean(tf.square(self._prediction - self._y)))
            
            ten_names = [n for n in tf.compat.v1.get_default_graph().as_graph_def().node]
            def name_exist(name):
                for n in ten_names:
                    if name in n.name: return True
                return False
            
            for i in range(1000):
                if not name_exist("transfer_%d_%s"%(i, self.optim_name)): 
                    transfer_scope = "transfer_%d_%s"%(i, self.optim_name)
                    break

            self._optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=self.learning_rate, name=transfer_scope)

            self._mininize = self._optimizer.minimize(self._loss, var_list=transfer_vars)

            # self.merged_summary_op = tf.compat.v1.summary.merge_all()

            adam_vars = self.GetVariablesBySubstr("_power") + self.GetVariablesBySubstr(self.optim_name)
            init_op = tf.compat.v1.variables_initializer(var_list=adam_vars, name="transfer_init")
            # init_op = tf.compat.v1.global_variables_initializer()
            self.sess.run(init_op)

    def restart_session(self):
        self.close()
        self.init_session()

    def init_session(self):
        if self.sess is None:
            config = tf.compat.v1.ConfigProto()
            config.gpu_options.allow_growth = True
            config.intra_op_parallelism_threads = self.num_threads
            self.sess = tf.compat.v1.Session(config=config, graph=self.graph)

    def close(self):
        if self.sess is not None:
            self.sess.close()
            self.sess = None

    def GetVariableList(self):
        # self.var_list = [v for v in tf.compat.v1.global_variables() if not self.optim_name in v.name and not "_power" in v.name]
        self.var_list = [v for v in tf.compat.v1.global_variables() if not self.optim_name in v.name and not "_power" in v.name]
        return self.var_list

    def GetTensorList(self):
        # self.ten_list = [n for n in tf.compat.v1.get_default_graph().as_graph_def().node]
        self.ten_list = [n for n in tf.compat.v1.get_default_graph().as_graph_def().node]
        return self.ten_list

    def SaveVariableName(self, save_as):
        with open(save_as + ".varname", "w") as f:
            for var in self.var_list: f.write(var.name + "\n")

    def SaveTensorName(self, save_as):
        with open(save_as + ".tenname", "w") as f:
            for ten in self.ten_list: f.write(ten.name + "\n")

    def GetTensorByName(self, name):
        return tf.compat.v1.get_default_graph().get_tensor_by_name(name)

    def GetVariableByName(self, name):
        l = [v for v in tf.compat.v1.global_variables() if v.name == name]
        if len(l) == 0:
            return None
        else:
            return l[0]

    def GetVariablesBySubstr(self, substr):
        return [v for v in tf.compat.v1.global_variables() if substr in v.name]

    def GetVariableListByScope(self, scope):
        return tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.GLOBAL_VARIABLES, scope=scope)

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
                    ch = np.float32(bbatch.get_data_from_batch(batch, "charge"))
                    x = np.float32(bbatch.get_data_from_batch(batch, "x"))
                    x = to_model_feature(x)
                    mod_x = np.float32(bbatch.get_data_from_batch(batch, "mod_x"))
                    instrument = np.float32(bbatch.get_data_from_batch(batch, "instrument"))
                    nce = np.float32(bbatch.get_data_from_batch(batch, "nce"))
                    y = bbatch.get_data_from_batch(batch, "y")

                    feed_dict = {
                        self._aa_x: x,
                        self._mod_x: mod_x,
                        self._charge: ch,
                        self._time_step: peplen - 1,
                        self._nce: nce,
                        self._instrument: instrument,
                        self._y: y,
                    }
                    if self.rnn_kp is not None: feed_dict[self.rnn_kp] = 1 - self.dropout
                    else: feed_dict[self._rnn_dropout] = self.dropout
                    if self.output_kp is not None: feed_dict[self.output_kp] = 1 - self.dropout
                    else: feed_dict[self._dropout] = self.dropout

                    cost, _ = self.sess.run([self._loss, self._mininize], feed_dict=feed_dict)
                    end_time = time.perf_counter()
                    batch_time_cost.append(end_time - start_time)
                    batch_cost.append(cost)

                    print("Epoch={:3d}, peplen={:3d}, Batch={:4d}, size={:4d}, cost={:.4f}, time={:.2f}s".format(
                            epoch + 1, peplen[0], ith_batch, len(x), cost, end_time - start_time), end="\r")
                mean_costs.append("Epoch={:3d}, mean_cost={:.4f}, time={:.2f}s".format(epoch + 1, np.mean(batch_cost),
                                np.sum(batch_time_cost)))
                print("\n" + mean_costs[-1])

            print("")
            for l in mean_costs:
                print(l)

            if save_as is not None:
                self.SaveModel(save_as)

    def Predict(self, buckets):
        with self.graph.as_default():
            start = time.perf_counter()
            bbatch = Bucket_Batch(buckets, batch_size=self.batch_size, shuffle=False)
            output_buckets = {}
            count = 0
            batch_count = 0
            for batch in bbatch.generate_batch():
                peplen = bbatch.get_data_from_batch(batch, "peplen")
                ch = np.float32(bbatch.get_data_from_batch(batch, "charge"))
                x = np.float32(bbatch.get_data_from_batch(batch, "x"))
                x = to_model_feature(x)
                mod_x = np.float32(bbatch.get_data_from_batch(batch, "mod_x"))
                instrument = np.float32(bbatch.get_data_from_batch(batch, "instrument"))
                nce = np.float32(bbatch.get_data_from_batch(batch, "nce"))
                count += peplen.shape[0]
                batch_count += 1

                feed_dict = {
                    self._aa_x: x,
                    self._mod_x: mod_x,
                    self._charge: ch,
                    self._time_step: peplen - 1,
                    self._nce: nce,
                    self._instrument: instrument
                }
                predictions = self.sess.run(self._prediction, feed_dict=feed_dict)
                predictions[predictions > 1] = 1
                predictions[predictions < 0] = 0
                _buckets = {peplen[0]: (predictions,)}
                output_buckets = merge_buckets(output_buckets, _buckets)
                if (batch_count % 10) == 0: print("[pDeep Info] predicted %d peptide precursors using %.3f seconds"%(count, time.perf_counter()-start), end='\r')
            print("[pDeep Info] predicted %d peptide precursors using %.3f seconds"%(count, time.perf_counter()-start))
            return output_buckets

    def SaveModel(self, model_file):
        dir = os.path.dirname(model_file)
        if not os.path.exists(dir): os.makedirs(dir)
        # print("enter tf.train.Saver")
        saver = tf.compat.v1.train.Saver(var_list=self.var_list)
        # print("enter saver.save")
        save_path = saver.save(self.sess, model_file)
        print("Model save as %s" % save_path)
        self.SaveVariableName(model_file)
        self.SaveTensorName(model_file)

    # model
    def LoadModel(self, model_file):
        print('[pDeep Info] model = {}'.format(model_file))
        if use_tf2: tf.compat.v1.disable_eager_execution()
        # run multiple instances:
        # g = tf.Graph()
        # with g.as_default():
        self.restart_session()

        with self.graph.as_default():
            saver = tf.compat.v1.train.import_meta_graph(model_file + ".meta")
            saver.restore(self.sess, model_file)
            graph = tf.compat.v1.get_default_graph()
            
            # print("[pDeep Info] # trainable parameters = %s"%tf_ops.count_parameters())
            
            self.GetVariableList()
            self.GetTensorList()
            self._aa_x = graph.get_tensor_by_name("input_aa_x:0")
            self._mod_x = graph.get_tensor_by_name("input_mod_x:0")
            self._y = graph.get_tensor_by_name("input_y:0")
            self._time_step = graph.get_tensor_by_name("input_time_step:0")
            self._charge = graph.get_tensor_by_name("input_charge:0")
            self._instrument = graph.get_tensor_by_name("input_instrument:0")
            self._nce = graph.get_tensor_by_name("input_nce:0")
            self._prediction = graph.get_tensor_by_name("output_scope/output:0")
            try:
                self._rnn_dropout = graph.get_tensor_by_name("rnn_dropout:0")
                self._dropout = graph.get_tensor_by_name("dropout:0")
                self.rnn_kp = None
                self.output_kp = None
            except:
                self.rnn_kp = graph.get_tensor_by_name("rnn_keep_prob:0")
                self.output_kp = graph.get_tensor_by_name("output_keep_prob:0")
