from keras import models, Model
from keras.layers import Input, Embedding, Conv1D, Add, Dense, SpatialDropout1D

def _conv_block(x, filters, kernel_size):
    conv = Conv1D(filters, kernel_size, activation='relu', padding='same')(x)
    conv = Conv1D(filters, kernel_size, activation='relu', padding='same')(conv)
    return conv


def _resblock(x, filters, kernel_size):
    conv = _conv_block(x, filters, kernel_size)
    projection = Conv1D(filters, 1, padding='same')(x)
    return Add()([conv, projection])

### main model ###
def predict_bindsite(vocab_size,
                    char_embedding_size,
                    base_filters,
                    doc_embedding_size,
                    dropout):
    text = Input(shape=(None,))
    embedding = Embedding(vocab_size, char_embedding_size)(text)

    conv_1 = _resblock(embedding, base_filters, 3)
    conv_1 = SpatialDropout1D(dropout)(conv_1)
    conv_2 = _resblock(conv_1, base_filters * 2, 3)
    conv_2 = SpatialDropout1D(dropout)(conv_2)
    conv_3 = _resblock(conv_2, base_filters * 4, 3)
    conv_3 = SpatialDropout1D(dropout)(conv_3)
    conv_4 = _resblock(conv_3, base_filters * 8, 3)
    conv_4 = SpatialDropout1D(dropout)(conv_4)

    attention = AttentionWithContext()(conv_4)

    fc_1 = Dense(doc_embedding_size, activation='relu')(attention)
    fc_2 = Dense(doc_embedding_size, activation='relu')(fc_1)
    # prediction = Dense(1, activation='sigmoid')(fc_2)
    # prediction = Dense(29, activation='sigmoid')(fc_2)
    prediction = Dense(30, activation='sigmoid')(fc_2) # Originally sigmoid

    model = Model(text, prediction)
    # model.compile('adam', 'binary_crossentropy', metrics=['acc'])
    model.compile(optimizer='rmsprop', loss='categorical_crossentropy', metrics=['accuracy'])


    return model