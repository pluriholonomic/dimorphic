#pragma once

#include <glog/logging.h>
#include <hdf5.h>
#include <cassert>
#include <cstdio>
#include <string>
#include "blosc/blosc_filter.h"

#define HSIZE_ALLOC(ptr_name, sz) \
    ptr_name = (hsize_t *)malloc(sz * sizeof(hsize_t))

namespace Dimers {

//
// Writer is an abstract class to serve as an interface for writers of lattice
// trajectories. It mainly handles HDF5 boilerplate, whereas child classes will
// handle the precise writing of a frame.
//
template <typename T>
struct Writer {
    Writer(std::string outfile, int *dims, int ndims, bool use_blosc = false,
           bool use_int_direction = false);
    Writer(const Writer &) = delete;

    virtual ~Writer();
    virtual void write(T *lattice, int *positions) = 0;

    std::string _outfile;

    int *_dims;
    int _ndims;
    size_t _nwritten;

    herr_t status;
    hid_t file, prop;
    hid_t pos_dataspace, pos_dataset;
    hid_t dir_dataspace, dir_dataset;
    hsize_t *current_dims, *offset, *chunk_dims, *max_dims, *frame_size;
    unsigned int _cd_vals[7];  // N.B. For blosc
    char *version, *date;
    bool use_blosc;
    bool use_int_direction;
};

template <typename T>
Writer<T>::Writer(std::string outfile, int *dims, int ndims, bool use_blosc,
                  bool use_int_direction)
    : _outfile(outfile),
      _ndims(ndims),
      _nwritten(0),
      use_blosc(use_blosc),
      use_int_direction(use_int_direction) {
    // N.B. dims contains all non-time dimensions, so we really are writing an
    // ndim+1 rank tensor
    if (ndims < 2) {
        LOG(FATAL)
            << "Attempting to write a 1-dimensional time series, failing";
    }

    file = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    printf("Created HDF5 file %s with descriptor %i\n", outfile.c_str(),
           (int)file);

    // Store dimensions locally
    _dims = (int *)malloc(ndims * sizeof(int));
    for (int i = 0; i < _ndims; i++) {
        _dims[i] = dims[i];
    }

    // Allocate hsize_t arrays
    HSIZE_ALLOC(current_dims, ndims + 1);
    HSIZE_ALLOC(offset, ndims + 1);
    HSIZE_ALLOC(chunk_dims, ndims + 1);
    HSIZE_ALLOC(max_dims, ndims + 1);
    HSIZE_ALLOC(frame_size, ndims + 1);

    for (int i = 0; i < ndims + 1; i++) {
        max_dims[i] = H5S_UNLIMITED;
        offset[i] = 0;
    }
    max_dims[0] = max_dims[1] = max_dims[2] = H5S_UNLIMITED;
    offset[0] = offset[1] = offset[2] = 0;

    // FIXME: Coarse guess, probably should be defined as a function of dims,
    // ndims
    frame_size[0] = 1;
    current_dims[0] = chunk_dims[0] = 256;

    for (int i = 1; i < ndims + 1; i++) {
        current_dims[i] = chunk_dims[i] = frame_size[i] = _dims[i - 1];
    }

    // Add chunking, shuffle, compression
    prop = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(prop, H5D_CHUNKED);
    status = H5Pset_chunk(prop, ndims + 1, chunk_dims);

    // Setup blosc
    if (use_blosc) {
        status = register_blosc(&version, &date);
        if (status) {
            printf("blosc filter registered; version (%s), date (%s)\n",
                   version, date);
        }

        _cd_vals[4] = 5;  // blosc level
        _cd_vals[5] = 1;  // use shuffle
        _cd_vals[6] = BLOSC_BLOSCLZ;
        status =
            H5Pset_filter(prop, FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, _cd_vals);
        if (status < 0) {
            printf("Blosc failed\n");
        }
    }

    // Setup dataspace, sets
    pos_dataspace = H5Screate_simple(ndims + 1, current_dims, max_dims);
    pos_dataset = H5Dcreate2(file, "position", H5T_NATIVE_INT, pos_dataspace,
                             H5P_DEFAULT, prop, H5P_DEFAULT);

    dir_dataspace = H5Screate_simple(ndims + 1, current_dims, max_dims);
    if (use_int_direction) {
        dir_dataset = H5Dcreate2(file, "direction", H5T_NATIVE_INT,
                                 dir_dataspace, H5P_DEFAULT, prop, H5P_DEFAULT);
    } else {
        dir_dataset = H5Dcreate2(file, "direction", H5T_NATIVE_DOUBLE,
                                 dir_dataspace, H5P_DEFAULT, prop, H5P_DEFAULT);
    }

    assert(pos_dataset >= 0);
    assert(dir_dataset >= 0);
}

template <typename T>
Writer<T>::~Writer() {
    if (_nwritten < current_dims[0]) {
        current_dims[0] = _nwritten;
        status = H5Dset_extent(pos_dataset, current_dims);
        status = H5Dset_extent(dir_dataset, current_dims);
    }

    status = H5Pclose(prop);

    status = H5Sclose(pos_dataspace);
    status = H5Sclose(dir_dataspace);

    status = H5Dclose(pos_dataset);
    status = H5Dclose(dir_dataset);

    status = H5Fclose(file);
}

}  // namespace Dimers
