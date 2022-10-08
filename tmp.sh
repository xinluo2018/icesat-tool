inp=data/rs-img/dianchi_s2_20200511_wgs84.tif
out=data/rs-img/dianchi_s2_20200511_wgs84_100m.tif
gdal_translate -outsize 20% 20% -r average -co COMPRESS=LZW $inp $out 

