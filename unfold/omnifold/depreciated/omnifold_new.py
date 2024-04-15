import numpy as np
import tensorflow as tf
import tensorflow.keras.backend as K
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping


def reweight(events,model,batch_size=10000):
    f = model.predict(events, batch_size=batch_size)
    weights = f / (1. - f)
    return np.squeeze(np.nan_to_num(weights))

# Binary crossentropy for classifying two samples with weights
# Weights are "hidden" by zipping in y_true (the labels)

def weighted_binary_crossentropy(y_true, y_pred):
    weights = tf.gather(y_true, [1], axis=1) # event weights
    y_true = tf.gather(y_true, [0], axis=1) # actual y_true for loss

    # Clip the prediction value to prevent NaN's and Inf's
    epsilon = K.epsilon()
    y_pred = K.clip(y_pred, epsilon, 1. - epsilon)

    t_loss = -weights * ((y_true) * K.log(y_pred) +
                         (1 - y_true) * K.log(1 - y_pred))

    return K.mean(t_loss)

def omnifold(trainGen, trainReco, testReco, trainWts, testWts, iterations,model,verbose=0):

    weights = np.empty(shape=(iterations, 2, len(trainGen)))
    # shape = (iteration, step, event)

    X_det = np.concatenate((trainReco, testReco))
    Y_det = np.concatenate((np.zeros(len(trainReco)), np.ones(len(testReco))))

    X_model = np.concatenate((trainGen, trainGen))
    Y_model = np.concatenate((np.zeros(len(trainGen)), np.ones(len(trainGen))))

    # initial iterative weights are ones
    weights_pull = trainWts
    weights_push = trainWts

    # Early stopping criteria to avoid overtraining
    earlystopping = EarlyStopping(patience=4,
                              verbose=1,
                              restore_best_weights=True)

    for i in range(iterations):

        if (verbose>0):
            print("\nITERATION: {}\n".format(i + 1))
            pass

        # STEP 1: classify Sim. (which is reweighted by weights_push) to Data
        # weights reweighted Sim. --> Data

        if (verbose>0):
            print("STEP 1\n")
            pass

        weights_det = np.concatenate((weights_push, testWts))

        X_train_det, X_test_det, Y_train_det, Y_test_det, w_train_det, w_test_det = train_test_split(X_det, Y_det, weights_det)

        # zip ("hide") the weights with the labels
        Y_train_det = np.stack((Y_train_det, w_train_det), axis=1)
        Y_test_det = np.stack((Y_test_det, w_test_det), axis=1)   

        model.compile(loss=weighted_binary_crossentropy,
                      optimizer='Adam',
                      metrics=['accuracy'])

        model.fit(X_train_det,
                  Y_train_det,
                  epochs=100,
                  batch_size=5000,
                  validation_data=(X_test_det, Y_test_det),
                  callbacks=[earlystopping],
                  verbose=verbose)

        weights_pull = weights_push * reweight(trainReco,model)
        weights[i, :1, :] = weights_pull

        # STEP 2: classify Gen. to reweighted Gen. (which is reweighted by weights_pull)
        # weights Gen. --> reweighted Gen.

        if (verbose>0):
            print("\nSTEP 2\n")
            pass

        weights_model = np.concatenate((trainWts, weights_pull))
        # ones for Gen. (not MC weights), actual weights for (reweighted) Gen.

        X_train_model, X_test_model, Y_train_model, Y_test_model, w_train_model, w_test_model = train_test_split(X_model, Y_model, weights_model)

        # zip ("hide") the weights with the labels
        Y_train_model = np.stack((Y_train_model, w_train_model), axis=1)
        Y_test_model = np.stack((Y_test_model, w_test_model), axis=1)   

        model.compile(loss=weighted_binary_crossentropy,
                      optimizer='Adam',
                      metrics=['accuracy'])
        model.fit(X_train_model,
                  Y_train_model,
                  epochs=100,
                  batch_size=10000,
                  validation_data=(X_test_model, Y_test_model),
                  callbacks=[earlystopping],
                  verbose=verbose)

        weights_push = reweight(trainGen,model)
        weights[i, 1:2, :] = weights_push
        pass

    return weights
