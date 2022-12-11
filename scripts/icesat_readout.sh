### author: xin luo
### create: 2022.10.12
### des: read in and write out the icesat data.

cd /Users/luo/OneDrive/GitHub/icesat-tool

### --- Dianchi
## --- icesat2/alt13 dataset
## -- readout
python utils_main/read_atl13.py data/icesat2/dianchi_ATL13_down/*ATL13*.h5 -o data/icesat2/dianchi_ATL13_readout -n 4
## -- merge
python utils_main/merge_files.py data/icesat2/dianchi_ATL13_readout/*ATL13*_readout.h5 -o data/icesat2/dianchi_ATL13_readout/dianchi_ATL13_2019_1q.h5
## -- subset
python utils_main/subset_icesat.py data/icesat2/dianchi_ATL13_readout/dianchi_ATL13_2019_1q.h5 -m ./data/rs-img/dianchi_s2_20200511_wat_wgs84.tif 

# ## --- icesat2/alt03 dataset
# ## -- readout
# python utils_main/read_atl03.py data/icesat2/dianchi_ATL03_down/*ATL03*.h5 -o data/icesat2/dianchi_ATL03_readout -n 4
# ## -- merge
# python utils_main/merge_files.py data/icesat2/dianchi_ATL03_readout/*ATL03*_readout.h5 -o data/icesat2/dianchi_ATL03_readout/dianchi_ATL03_2019_1q.h5
# ## -- subset
# python utils_main/subset_icesat.py data/icesat2/dianchi_ATL03_readout/dianchi_ATL03_2019_1q.h5 -m ./data/rs-img/dianchi_s2_20200511_wat_wgs84.tif 

