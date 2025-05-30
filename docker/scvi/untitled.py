

RUN ln -s /opt/Python-${PYTHON3_VERSION}/python /opt/Python-${PYTHON3_VERSION}/python3

COPY ./requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt

# CUDA-enabled jaxlib is needed
RUN python3 -m pip install --upgrade "jax[cuda12_pip]==${JAX_VERSION}" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

pip install git+https://github.com/AllenInstitute/cell_type_mapper.git

conda create -n test3 python=3.12 pip 

pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

pip install --upgrade "jax[cuda11]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install scvi-tools scikit-image scikit-learn scikit-misc seaborn anndata muon tables pyarrow fastparquet igraph scib-metrics harmonypy
pip install git+https://github.com/AllenInstitute/cell_type_mapper.git

python3 -m pip install --upgrade "jax[cuda11_pip]==0.4.9" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# conda install -c conda-forge faiss-gpu # this could be the memory leak...

# replace leidenalg with igraph


# pip install --upgrade "jax[cuda12_pip]==${JAX_VERSION}" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

pip install git+https://github.com/AllenInstitute/cell_type_mapper.git
