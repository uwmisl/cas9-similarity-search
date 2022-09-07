FROM tensorflow/tensorflow:2.5.0-gpu-jupyter

WORKDIR /tf/

RUN apt-get update && apt-get install -y git wget

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip -q awscliv2.zip && \
    ./aws/install

ENV HOME=/tf/
RUN python -c "from tensorflow.keras.applications.vgg16 import VGG16; VGG16()"
RUN chmod -R go+w /tf/.keras

RUN pip install \
    requests \
    pillow \
    tqdm \
    h5py \
    tables \
    unireedsolomon \
    seaborn \
    tensorflow_probability==0.13.0 \
    biopython \
    fpdf

