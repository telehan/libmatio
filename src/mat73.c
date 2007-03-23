/** @file mat5.c
 * Matlab MAT version 5 file functions
 * @ingroup MAT
 */
/*
 * Copyright (C) 2005-2006   Christopher C. Hulbert
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "hdf5.h"
#include "matio.h"
#include "mat73.h"
#include "matio_private.h"

static hsize_t perm_dims[10];
static hsize_t dims1[2] = {1,1};

static const char *Mat_class_names[] = {
    "",
    "cell",
    "struct",
    "object",
    "char",
    "sparse",
    "double",
    "single",
    "int8",
    "uint8",
    "int16",
    "uint16",
    "int32",
    "uint32",
    "int64",
    "uint64",
    "function"
};

/*===========================================================================
 *  Private functions
 *===========================================================================
 */
static int Mat_class_str_to_id(const char *name);

static int
Mat_class_str_to_id(const char *name)
{
    int id = 0;
    if ( NULL != name ) {
        int k;
        for ( k = 1; k < 17; k++ ) {
            if ( !strcmp(name,Mat_class_names[k]) ) {
                id = k;
                break;
            }
        }
    }
    return id;
}

static hid_t
Mat_class_type_to_hid_t(enum matio_classes class_type)
{
    switch ( class_type ) {
        case MAT_C_DOUBLE:
            return H5T_NATIVE_DOUBLE;
        case MAT_C_SINGLE:
            return H5T_NATIVE_FLOAT;
        case MAT_C_INT64:
#       if CHAR_BIT*SIZEOF_SHORT == 64
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 64
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 64
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 64
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_C_INT32:
#       if CHAR_BIT == 32
            return H5T_NATIVE_SCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 32
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 32
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 32
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 32
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_C_UINT32:
#       if CHAR_BIT == 32
            return H5T_NATIVE_UCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 32
            return H5T_NATIVE_USHORT;
#       elif CHAR_BIT*SIZEOF_INT == 32
            return H5T_NATIVE_UINT;
#       elif CHAR_BIT*SIZEOF_LONG == 32
            return H5T_NATIVE_ULONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 32
            return H5T_NATIVE_ULLONG;
#       endif
        case MAT_C_INT16:
#       if CHAR_BIT == 16
            return H5T_NATIVE_SCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 16
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 16
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 16
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 16
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_C_UINT16:
#       if CHAR_BIT == 16
            return H5T_NATIVE_UCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 16
            return H5T_NATIVE_USHORT;
#       elif CHAR_BIT*SIZEOF_INT == 16
            return H5T_NATIVE_UINT;
#       elif CHAR_BIT*SIZEOF_LONG == 16
            return H5T_NATIVE_ULONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 16
            return H5T_NATIVE_ULLONG;
#       endif
        case MAT_C_INT8:
#       if CHAR_BIT == 8
            return H5T_NATIVE_SCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 8
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 8
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 8
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 8
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_C_UINT8:
#       if CHAR_BIT == 8
            return H5T_NATIVE_UCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 8
            return H5T_NATIVE_USHORT;
#       elif CHAR_BIT*SIZEOF_INT == 8
            return H5T_NATIVE_UINT;
#       elif CHAR_BIT*SIZEOF_LONG == 8
            return H5T_NATIVE_ULONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 8
            return H5T_NATIVE_ULLONG;
#       endif
       default:
           return -1;
    }
}

static hid_t
Mat_data_type_to_hid_t(enum matio_types data_type)
{
    switch ( data_type ) {
        case MAT_T_DOUBLE:
            return H5T_NATIVE_DOUBLE;
        case MAT_T_SINGLE:
            return H5T_NATIVE_FLOAT;
        case MAT_T_INT64:
#       if CHAR_BIT*SIZEOF_SHORT == 64
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 64
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 64
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 64
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_T_INT32:
#       if CHAR_BIT == 32
            return H5T_NATIVE_SCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 32
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 32
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 32
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 32
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_T_UINT32:
#       if CHAR_BIT == 32
            return H5T_NATIVE_UCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 32
            return H5T_NATIVE_USHORT;
#       elif CHAR_BIT*SIZEOF_INT == 32
            return H5T_NATIVE_UINT;
#       elif CHAR_BIT*SIZEOF_LONG == 32
            return H5T_NATIVE_ULONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 32
            return H5T_NATIVE_ULLONG;
#       endif
        case MAT_T_INT16:
#       if CHAR_BIT == 16
            return H5T_NATIVE_SCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 16
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 16
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 16
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 16
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_T_UINT16:
#       if CHAR_BIT == 16
            return H5T_NATIVE_UCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 16
            return H5T_NATIVE_USHORT;
#       elif CHAR_BIT*SIZEOF_INT == 16
            return H5T_NATIVE_UINT;
#       elif CHAR_BIT*SIZEOF_LONG == 16
            return H5T_NATIVE_ULONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 16
            return H5T_NATIVE_ULLONG;
#       endif
        case MAT_T_INT8:
#       if CHAR_BIT == 8
            return H5T_NATIVE_SCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 8
            return H5T_NATIVE_SHORT;
#       elif CHAR_BIT*SIZEOF_INT == 8
            return H5T_NATIVE_INT;
#       elif CHAR_BIT*SIZEOF_LONG == 8
            return H5T_NATIVE_LONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 8
            return H5T_NATIVE_LLONG;
#       endif
        case MAT_T_UINT8:
#       if CHAR_BIT == 8
            return H5T_NATIVE_UCHAR;
#       elif CHAR_BIT*SIZEOF_SHORT == 8
            return H5T_NATIVE_USHORT;
#       elif CHAR_BIT*SIZEOF_INT == 8
            return H5T_NATIVE_UINT;
#       elif CHAR_BIT*SIZEOF_LONG == 8
            return H5T_NATIVE_ULONG;
#       elif CHAR_BIT*SIZEOF_LONG_LONG == 8
            return H5T_NATIVE_ULLONG;
#       endif
       default:
           return -1;
    }
}

static int
Mat_WriteNextStructField73(hid_t id,matvar_t *matvar,const char *name)
{
    unsigned long k,numel;
    hid_t mspace_id,dset_id,attr_type_id,attr_id,aspace_id;

    if ( NULL == matvar )
        return -1;

    switch ( matvar->class_type ) {
        case MAT_C_DOUBLE:
        case MAT_C_SINGLE:
        case MAT_C_INT32:
        case MAT_C_UINT32:
        case MAT_C_INT16:
        case MAT_C_UINT16:
        case MAT_C_INT8:
        case MAT_C_UINT8:
            numel = 1;
            for ( k = 0; k < matvar->rank; k++ ) {
                perm_dims[k] = matvar->dims[matvar->rank-k-1];
                numel *= perm_dims[k];
            }

            if ( matvar->isComplex ) {
                hid_t h5_complex,h5_complex_base;
                void *buf;

                h5_complex_base = Mat_class_type_to_hid_t(matvar->class_type);
                h5_complex      = H5Tcreate(H5T_COMPOUND,
                                      2*H5Tget_size(h5_complex_base));
                H5Tinsert(h5_complex,"real",0,h5_complex_base);
                H5Tinsert(h5_complex,"imag",H5Tget_size(h5_complex_base),
                          h5_complex_base);

                /* Not very memory efficient! */
                buf = malloc(2*numel*H5Tget_size(h5_complex_base));
                if ( NULL != buf ) {
                    switch ( matvar->class_type ) {
                        case MAT_C_DOUBLE:
                        {
                            double *dst = buf,
                                   *r=((struct ComplexSplit*)matvar->data)->Re,
                                   *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_SINGLE:
                        {
                            float *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_INT32:
                        {
                            mat_int32_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_UINT32:
                        {
                            mat_uint32_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_INT16:
                        {
                            mat_int16_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_UINT16:
                        {
                            mat_uint16_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_INT8:
                        {
                            mat_int8_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_UINT8:
                        {
                            mat_uint8_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                    }
                    mspace_id = H5Screate_simple(matvar->rank,perm_dims,NULL);
                    dset_id = H5Dcreate(id,name,h5_complex,mspace_id,
                                        H5P_DEFAULT);
                    attr_type_id = H5Tcopy(H5T_C_S1);
                    H5Tset_size(attr_type_id,
                                strlen(Mat_class_names[matvar->class_type])+1);
                    aspace_id = H5Screate(H5S_SCALAR);
                    attr_id = H5Acreate(dset_id,"MATLAB_class",attr_type_id,
                                        aspace_id,H5P_DEFAULT);
                    H5Awrite(attr_id,attr_type_id,
                             Mat_class_names[matvar->class_type]);
                    H5Sclose(aspace_id);
                    H5Aclose(attr_id);
                    H5Tclose(attr_type_id);
                    H5Dwrite(dset_id,h5_complex,H5S_ALL,H5S_ALL,H5P_DEFAULT,
                             buf);
                    H5Dclose(dset_id);
                    H5Sclose(mspace_id);
                    free(buf);
                }
                /* h5_complex_base is not a copy, so don't release it */
                H5Tclose(h5_complex);
            } else {
                mspace_id = H5Screate_simple(matvar->rank,perm_dims,NULL);
                dset_id = H5Dcreate(id,name,
                    Mat_class_type_to_hid_t(matvar->class_type),mspace_id,
                    H5P_DEFAULT);
                attr_type_id = H5Tcopy(H5T_C_S1);
                H5Tset_size(attr_type_id,
                            strlen(Mat_class_names[matvar->class_type])+1);
                aspace_id = H5Screate(H5S_SCALAR);
                attr_id = H5Acreate(dset_id,"MATLAB_class",attr_type_id,
                                    aspace_id,H5P_DEFAULT);
                H5Awrite(attr_id,attr_type_id,
                         Mat_class_names[matvar->class_type]);
                H5Sclose(aspace_id);
                H5Aclose(attr_id);
                H5Tclose(attr_type_id);
                H5Dwrite(dset_id,Mat_data_type_to_hid_t(matvar->data_type),
                    H5S_ALL,H5S_ALL,H5P_DEFAULT,matvar->data);
                H5Dclose(dset_id);
                H5Sclose(mspace_id);
            }
            break;
        case MAT_C_CHAR:
        {
            int matlab_int_decode = 2;
            for ( k = 0; k < matvar->rank; k++ )
                perm_dims[k] = matvar->dims[matvar->rank-k-1];

            mspace_id = H5Screate_simple(matvar->rank,perm_dims,NULL);
            switch ( matvar->data_type ) {
                case MAT_T_UTF32:
                case MAT_T_INT32:
                case MAT_T_UINT32:
                    /* Not sure matlab will actually handle this */
                    dset_id = H5Dcreate(id,name,
                        Mat_class_type_to_hid_t(MAT_C_UINT32),mspace_id,
                        H5P_DEFAULT);
                    break;
                case MAT_T_UTF16:
                case MAT_T_UTF8:
                case MAT_T_INT16:
                case MAT_T_UINT16:
                case MAT_T_INT8:
                case MAT_T_UINT8:
                    dset_id = H5Dcreate(id,name,
                        Mat_class_type_to_hid_t(MAT_C_UINT16),mspace_id,
                        H5P_DEFAULT);
                    break;
            }
            attr_type_id = H5Tcopy(H5T_C_S1);
            H5Tset_size(attr_type_id,
                        strlen(Mat_class_names[matvar->class_type])+1);
            aspace_id = H5Screate(H5S_SCALAR);
            attr_id = H5Acreate(dset_id,"MATLAB_class",attr_type_id,
                                aspace_id,H5P_DEFAULT);
            H5Awrite(attr_id,attr_type_id,Mat_class_names[matvar->class_type]);
            H5Aclose(attr_id);
            H5Tclose(attr_type_id);

            attr_type_id = H5Tcopy(H5T_NATIVE_INT);
            attr_id = H5Acreate(dset_id,"MATLAB_int_decode",attr_type_id,
                                aspace_id,H5P_DEFAULT);
            H5Awrite(attr_id,attr_type_id,&matlab_int_decode);
            H5Tclose(attr_type_id);
            H5Sclose(aspace_id);

            H5Dwrite(dset_id,Mat_data_type_to_hid_t(matvar->data_type),
                H5S_ALL,H5S_ALL,H5P_DEFAULT,matvar->data);
            H5Dclose(dset_id);
            H5Sclose(mspace_id);
            break;
        }
        case MAT_C_STRUCT:
        {
            hid_t struct_id,str_type_id,fieldnames_id;
            hsize_t    nfields,nmemb;
            matvar_t **fields;
            hvl_t     *fieldnames;
            char       id_name[128] = {'\0',};
            int        is_ref;

            (void)H5Iget_name(id,id_name,127);
            is_ref = !strcmp(id_name,"/#refs#");
            struct_id = H5Gcreate(id,name,0);
            if ( struct_id < 0 ) {
                Mat_Critical("Error creating group for struct %s",name);
            } else {
                str_type_id = H5Tcopy(H5T_C_S1);
                H5Tset_size(str_type_id,7);
                aspace_id = H5Screate(H5S_SCALAR);
                attr_id = H5Acreate(struct_id,"MATLAB_class",str_type_id,
                                    aspace_id,H5P_DEFAULT);
                H5Awrite(attr_id,str_type_id,"struct");
                H5Aclose(attr_id);

                nmemb = matvar->dims[0];
                for ( k = 1; k < matvar->rank; k++ )
                    nmemb *= matvar->dims[k];
                nfields = matvar->nbytes / (nmemb*matvar->data_size);

                fieldnames = malloc(nfields*sizeof(*fieldnames));
                fields     = matvar->data;
                for ( k = 0; k < nfields; k++ ) {
                    fieldnames[k].len = strlen(fields[k]->name);
                    fieldnames[k].p   = fields[k]->name;
                }
                H5Tset_size(str_type_id,1);
                fieldnames_id = H5Tvlen_create(str_type_id);
                aspace_id     = H5Screate_simple(1,&nfields,NULL);
                attr_id = H5Acreate(struct_id,"MATLAB_fields",fieldnames_id,
                                    aspace_id,H5P_DEFAULT);
                H5Awrite(attr_id,fieldnames_id,fieldnames);
                H5Aclose(attr_id);
                H5Sclose(aspace_id);
                H5Tclose(fieldnames_id);
                H5Tclose(str_type_id);
                free(fieldnames);

                if ( 1 == nmemb ) {
                    for ( k = 0; k < nmemb*nfields; k++ )
                        Mat_WriteNextStructField73(struct_id,fields[k],
                            fields[k]->name);
                } else {
                    hid_t refs_id;

                    if (is_ref) {
                        refs_id = id;
                    } else {
                        if ((refs_id=H5Gopen(id,"/#refs#") < 0 ))
                            refs_id = H5Gcreate(id,"/#refs#",0);
                    }
                    if ( refs_id > -1 ) {
                        char name[64];
                        hobj_ref_t **refs;
                        hsize_t      num_obj;
                        int l;

                        refs = malloc(nfields*sizeof(*refs));
                        for ( l = 0; l < nfields; l++ )
                            refs[l] = malloc(nmemb*sizeof(*refs[l]));

                        for ( k = 0; k < nmemb; k++ ) {
                            for ( l = 0; l < nfields; l++ ) {
                                (void)H5Gget_num_objs(refs_id,&num_obj);
                                sprintf(name,"%lu",num_obj);
                                Mat_WriteNextStructField73(refs_id,
                                    fields[k*nfields+l],name);
                                sprintf(name,"/#refs#/%lu",num_obj);
                                H5Rcreate(refs[l]+k,id,name,
                                          H5R_OBJECT,-1);
                            }
                        }
                        for ( k = 0; k < matvar->rank; k++ )
                            perm_dims[k] = matvar->dims[matvar->rank-k-1];

                        mspace_id=H5Screate_simple(matvar->rank,perm_dims,NULL);
                        for ( l = 0; l < nfields; l++ ) {
                            dset_id = H5Dcreate(struct_id,
                                fields[l]->name,H5T_STD_REF_OBJ,mspace_id,
                                H5P_DEFAULT);
                            H5Dwrite(dset_id,H5T_STD_REF_OBJ,H5S_ALL,H5S_ALL,
                                H5P_DEFAULT,refs[l]);
                            H5Dclose(dset_id);
                            free(refs[l]);
                        }
                        free(refs);
                        H5Sclose(mspace_id);
                        if ( !is_ref )
                            H5Gclose(refs_id);
                    }
                }
                H5Gclose(struct_id);
            }
            break;
        }
    }
    return 0;
}

/** @brief Creates a new Matlab MAT version 5 file
 *
 * Tries to create a new Matlab MAT file with the given name and optional
 * header string.  If no header string is given, the default string
 * is used containing the software, version, and date in it.  If a header
 * string is given, at most the first 116 characters is written to the file.
 * The given header string need not be the full 116 characters, but MUST be
 * NULL terminated.
 * @ingroup MAT
 * @param matname Name of MAT file to create
 * @param hdr_str Optional header string, NULL to use default
 * @return A pointer to the MAT file or NULL if it failed.  This is not a
 * simple FILE * and should not be used as one.
 */
mat_t *
Mat_Create73(const char *matname,const char *hdr_str)
{
    FILE *fp = NULL;
    mat_int16_t endian = 0, version;
    mat_t *mat = NULL;
    size_t err;
    time_t t;
    hid_t plist_id,fid;

    plist_id = H5Pcreate(H5P_FILE_CREATE);
    H5Pset_userblock(plist_id,512);
    fid = H5Fcreate(matname,H5F_ACC_TRUNC,plist_id,H5P_DEFAULT);
    H5Fclose(fid);
    H5Pclose(plist_id);

    fp = fopen(matname,"r+b");
    if ( !fp )
        return NULL;

    fseek(fp,0,SEEK_SET);

    mat = malloc(sizeof(*mat));
    if ( !mat ) {
        fclose(fp);
        return NULL;
    }

    mat->fp               = NULL;
    mat->header           = NULL;
    mat->subsys_offset    = NULL;
    mat->filename         = NULL;
    mat->version          = 0;
    mat->byteswap         = 0;
    mat->mode             = 0;
    mat->bof              = 0;
    mat->next_index       = 0;

    t = time(NULL);
    mat->filename = strdup_printf("%s",matname);
    mat->mode     = MAT_ACC_RDWR;
    mat->byteswap = 0;
    mat->header   = calloc(1,128);
    mat->subsys_offset = calloc(1,16);
    memset(mat->header,' ',128);
    if ( hdr_str == NULL ) {
        err = mat_snprintf(mat->header,116,"MATLAB 7.0 MAT-file, Platform: %s,"
                "Created by libmatio v%d.%d.%d on %s HDF5 schema 0.5",
                MATIO_PLATFORM,MATIO_MAJOR_VERSION,MATIO_MINOR_VERSION,
                MATIO_RELEASE_LEVEL,ctime(&t));
        mat->header[115] = '\0';    /* Just to make sure it's NULL terminated */    } else {
        err = mat_snprintf(mat->header,116,"%s",hdr_str);
    }
    mat->header[err] = ' ';
    mat_snprintf(mat->subsys_offset,15,"            ");
    mat->version = (int)0x0200;
    endian = 0x4d49;

    version = 0x0200;

    err = fwrite(mat->header,1,116,fp);
    err = fwrite(mat->subsys_offset,1,8,fp);
    err = fwrite(&version,2,1,fp);
    err = fwrite(&endian,2,1,fp);

    fclose(fp);

    fid = H5Fopen(matname,H5F_ACC_RDWR,H5P_DEFAULT);

    mat->fp = malloc(sizeof(hid_t));
    *(hid_t*)mat->fp = fid;

    return mat;
}

/** @brief Reads the header information for the next MAT variable
 *
 * @ingroup mat_internal
 * @param mat MAT file pointer
 * @retuen pointer to the MAT variable or NULL
 */
matvar_t *
Mat_VarReadNextInfo73( mat_t *mat )
{
    hid_t       fid,gid;
    hsize_t     num_objs;
    H5E_auto_t  efunc;
    void       *client_data;
    matvar_t   *matvar;

    if( mat == NULL )
        return NULL;

    fid = *(hid_t*)mat->fp;
    H5Gget_num_objs(fid,&num_objs);
    /* FIXME: follow symlinks, datatypes? */
    while ( H5G_DATASET != H5Gget_objtype_by_idx(fid,mat->next_index) &&
            mat->next_index < num_objs ) {
        mat->next_index++;
    }

    if ( mat->next_index >= num_objs )
        return NULL;

    matvar = Mat_VarCalloc();
    if ( NULL != matvar ) {
        ssize_t  name_len;
        /* FIXME */
        hsize_t  dims[10];
        hid_t   attr_id,type_id,dset_id,space_id;

        matvar->fp = mat;
        name_len = H5Gget_objname_by_idx(fid,mat->next_index,NULL,0);
        matvar->name = malloc(1+name_len);
        if ( matvar->name ) {
            name_len = H5Gget_objname_by_idx(fid,mat->next_index,matvar->name,
                                             1+name_len);
            matvar->name[name_len] = '\0';
        }
        dset_id = H5Dopen(fid,matvar->name);
        space_id = H5Dget_space(dset_id);
        matvar->rank = H5Sget_simple_extent_ndims(space_id);
        matvar->dims = malloc(matvar->rank*sizeof(*matvar->dims));
        if ( NULL != matvar->dims ) {
            int k;
            H5Sget_simple_extent_dims(space_id,dims,NULL);
            for ( k = 0; k < matvar->rank; k++ )
                matvar->dims[k] = dims[k];
        }
        H5Sclose(space_id);

        attr_id = H5Aopen_name(dset_id,"MATLAB_class");
        type_id  = H5Aget_type(attr_id);
        if ( H5T_STRING == H5Tget_class(type_id) ) {
            char *class_str = malloc(H5Tget_size(type_id));
            if ( NULL != class_str ) {
                hid_t class_id = H5Tcopy(H5T_C_S1);
                H5Tset_size(class_id,H5Tget_size(type_id));
                H5Aread(attr_id,class_id,class_str);
                H5Tclose(class_id);
                matvar->class_type = Mat_class_str_to_id(class_str);
                free(class_str);
            }
        }
        H5Tclose(type_id);
        H5Aclose(attr_id);

        /* Turn off error printing so testing for attributes doesn't print
         * error stacks
         */
        H5Eget_auto(&efunc,&client_data);
        H5Eset_auto((H5E_auto_t)0,NULL);

        attr_id = H5Aopen_name(dset_id,"MATLAB_global");
        /* FIXME: Check that dataspace is scalar */
        if ( -1 < attr_id ) {
            H5Aread(attr_id,H5T_NATIVE_INT,&matvar->isGlobal);
            H5Aclose(attr_id);
        }

        H5Eset_auto(efunc,client_data);
        H5Dclose(dset_id);
        mat->next_index++;
    }
    return matvar;
}

/** @brief Writes a matlab variable to a version 7.3 matlab file
 *
 * @ingroup mat_internal
 * @param mat MAT file pointer
 * @param matvar pointer to the mat variable
 * @param compress option to compress the variable
 *                 (only works for numeric types)
 * @retval 0 on success
 */
int
Mat_VarWrite73(mat_t *mat,matvar_t *matvar,int compress)
{
    unsigned long k,numel;
    hid_t mspace_id,dset_id,attr_type_id,attr_id,aspace_id;

    if ( NULL == mat || NULL == matvar )
        return -1;

    switch ( matvar->class_type ) {
        case MAT_C_DOUBLE:
        case MAT_C_SINGLE:
        case MAT_C_INT32:
        case MAT_C_UINT32:
        case MAT_C_INT16:
        case MAT_C_UINT16:
        case MAT_C_INT8:
        case MAT_C_UINT8:
            numel = 1;
            for ( k = 0; k < matvar->rank; k++ ) {
                perm_dims[k] = matvar->dims[matvar->rank-k-1];
                numel *= perm_dims[k];
            }

            if ( matvar->isComplex ) {
                hid_t h5_complex,h5_complex_base;
                void *buf;

                h5_complex_base = Mat_class_type_to_hid_t(matvar->class_type);
                h5_complex      = H5Tcreate(H5T_COMPOUND,
                                      2*H5Tget_size(h5_complex_base));
                H5Tinsert(h5_complex,"real",0,h5_complex_base);
                H5Tinsert(h5_complex,"imag",H5Tget_size(h5_complex_base),
                          h5_complex_base);

                /* Not very memory efficient! */
                buf = malloc(2*numel*H5Tget_size(h5_complex_base));
                if ( NULL != buf ) {
                    switch ( matvar->class_type ) {
                        case MAT_C_DOUBLE:
                        {
                            double *dst = buf,
                                   *r=((struct ComplexSplit*)matvar->data)->Re,
                                   *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_SINGLE:
                        {
                            float *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_INT32:
                        {
                            mat_int32_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_UINT32:
                        {
                            mat_uint32_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_INT16:
                        {
                            mat_int16_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_UINT16:
                        {
                            mat_uint16_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_INT8:
                        {
                            mat_int8_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                        case MAT_C_UINT8:
                        {
                            mat_uint8_t *dst = buf,
                                  *r=((struct ComplexSplit*)matvar->data)->Re,
                                  *i=((struct ComplexSplit*)matvar->data)->Im;
                            for ( k = numel; k--; ) {
                                *dst++ = *r++;
                                *dst++ = *i++;
                            }
                            break;
                        }
                    }
                    mspace_id = H5Screate_simple(matvar->rank,perm_dims,NULL);
                    dset_id = H5Dcreate(*(hid_t*)mat->fp,matvar->name,
                        h5_complex,mspace_id,H5P_DEFAULT);
                    attr_type_id = H5Tcopy(H5T_C_S1);
                    H5Tset_size(attr_type_id,
                                strlen(Mat_class_names[matvar->class_type])+1);
                    aspace_id = H5Screate(H5S_SCALAR);
                    attr_id = H5Acreate(dset_id,"MATLAB_class",attr_type_id,
                                        aspace_id,H5P_DEFAULT);
                    H5Awrite(attr_id,attr_type_id,
                             Mat_class_names[matvar->class_type]);
                    H5Sclose(aspace_id);
                    H5Aclose(attr_id);
                    H5Tclose(attr_type_id);
                    H5Dwrite(dset_id,h5_complex,H5S_ALL,H5S_ALL,H5P_DEFAULT,
                             buf);
                    H5Dclose(dset_id);
                    H5Sclose(mspace_id);
                    free(buf);
                }
                /* h5_complex_base is not a copy, so don't release it */
                H5Tclose(h5_complex);
            } else {

            mspace_id = H5Screate_simple(matvar->rank,perm_dims,NULL);
            dset_id = H5Dcreate(*(hid_t*)mat->fp,matvar->name,
                Mat_class_type_to_hid_t(matvar->class_type),mspace_id,
                H5P_DEFAULT);
            attr_type_id = H5Tcopy(H5T_C_S1);
            H5Tset_size(attr_type_id,
                        strlen(Mat_class_names[matvar->class_type])+1);
            aspace_id = H5Screate(H5S_SCALAR);
            attr_id = H5Acreate(dset_id,"MATLAB_class",attr_type_id,
                                aspace_id,H5P_DEFAULT);
            H5Awrite(attr_id,attr_type_id,Mat_class_names[matvar->class_type]);
            H5Sclose(aspace_id);
            H5Aclose(attr_id);
            H5Tclose(attr_type_id);
            H5Dwrite(dset_id,Mat_data_type_to_hid_t(matvar->data_type),
                H5S_ALL,H5S_ALL,H5P_DEFAULT,matvar->data);
            H5Dclose(dset_id);
            H5Sclose(mspace_id);
            }
            break;
        case MAT_C_CHAR:
        {
            int matlab_int_decode = 2;
            for ( k = 0; k < matvar->rank; k++ )
                perm_dims[k] = matvar->dims[matvar->rank-k-1];

            mspace_id = H5Screate_simple(matvar->rank,perm_dims,NULL);
            switch ( matvar->data_type ) {
                case MAT_T_UTF32:
                case MAT_T_INT32:
                case MAT_T_UINT32:
                    /* Not sure matlab will actually handle this */
                    dset_id = H5Dcreate(*(hid_t*)mat->fp,matvar->name,
                        Mat_class_type_to_hid_t(MAT_C_UINT32),mspace_id,
                        H5P_DEFAULT);
                    break;
                case MAT_T_UTF16:
                case MAT_T_UTF8:
                case MAT_T_INT16:
                case MAT_T_UINT16:
                case MAT_T_INT8:
                case MAT_T_UINT8:
                    dset_id = H5Dcreate(*(hid_t*)mat->fp,matvar->name,
                        Mat_class_type_to_hid_t(MAT_C_UINT16),mspace_id,
                        H5P_DEFAULT);
                    break;
            }
            attr_type_id = H5Tcopy(H5T_C_S1);
            H5Tset_size(attr_type_id,
                        strlen(Mat_class_names[matvar->class_type])+1);
            aspace_id = H5Screate(H5S_SCALAR);
            attr_id = H5Acreate(dset_id,"MATLAB_class",attr_type_id,
                                aspace_id,H5P_DEFAULT);
            H5Awrite(attr_id,attr_type_id,Mat_class_names[matvar->class_type]);
            H5Aclose(attr_id);
            H5Tclose(attr_type_id);

            attr_type_id = H5Tcopy(H5T_NATIVE_INT);
            attr_id = H5Acreate(dset_id,"MATLAB_int_decode",attr_type_id,
                                aspace_id,H5P_DEFAULT);
            H5Awrite(attr_id,attr_type_id,&matlab_int_decode);
            H5Tclose(attr_type_id);
            H5Sclose(aspace_id);

            H5Dwrite(dset_id,Mat_data_type_to_hid_t(matvar->data_type),
                H5S_ALL,H5S_ALL,H5P_DEFAULT,matvar->data);
            H5Dclose(dset_id);
            H5Sclose(mspace_id);
            break;
        }
        case MAT_C_STRUCT:
        {
            hid_t struct_id,str_type_id,fieldnames_id;
            hsize_t    nfields,nmemb;
            matvar_t **fields;
            hvl_t     *fieldnames;

            struct_id = H5Gcreate(*(hid_t*)mat->fp,matvar->name,0);
            if ( struct_id < 0 ) {
                Mat_Critical("Error creating group for struct %s",matvar->name);
            } else {
                str_type_id = H5Tcopy(H5T_C_S1);
                H5Tset_size(str_type_id,7);
                aspace_id = H5Screate(H5S_SCALAR);
                attr_id = H5Acreate(struct_id,"MATLAB_class",str_type_id,
                                    aspace_id,H5P_DEFAULT);
                H5Awrite(attr_id,str_type_id,"struct");
                H5Aclose(attr_id);

                nmemb = matvar->dims[0];
                for ( k = 1; k < matvar->rank; k++ )
                    nmemb *= matvar->dims[k];
                nfields = matvar->nbytes / (nmemb*matvar->data_size);

                fieldnames = malloc(nfields*sizeof(*fieldnames));
                fields     = matvar->data;
                for ( k = 0; k < nfields; k++ ) {
                    fieldnames[k].len = strlen(fields[k]->name);
                    fieldnames[k].p   = fields[k]->name;
                }
                H5Tset_size(str_type_id,1);
                fieldnames_id = H5Tvlen_create(str_type_id);
                aspace_id     = H5Screate_simple(1,&nfields,NULL);
                attr_id = H5Acreate(struct_id,"MATLAB_fields",fieldnames_id,
                                    aspace_id,H5P_DEFAULT);
                H5Awrite(attr_id,fieldnames_id,fieldnames);
                H5Aclose(attr_id);
                H5Sclose(aspace_id);
                H5Tclose(fieldnames_id);
                H5Tclose(str_type_id);
                free(fieldnames);

                if ( 1 == nmemb ) {
                    for ( k = 0; k < nmemb*nfields; k++ )
                        Mat_WriteNextStructField73(struct_id,fields[k],
                            fields[k]->name);
                } else {
                    hid_t refs_id;
                    if ((refs_id=H5Gopen(*(hid_t*)mat->fp,"/#refs#") < 0 )) {
                        refs_id = H5Gcreate(*(hid_t*)mat->fp,"/#refs#",0);
                    }
                    
                    if ( refs_id > -1 ) {
                        char name[64];
                        hobj_ref_t **refs;
                        hsize_t      num_obj;
                        int l;

                        refs = malloc(nfields*sizeof(*refs));
                        for ( l = 0; l < nfields; l++ )
                            refs[l] = malloc(nmemb*sizeof(*refs[l]));

                        for ( k = 0; k < nmemb; k++ ) {
                            for ( l = 0; l < nfields; l++ ) {
                                (void)H5Gget_num_objs(refs_id,&num_obj);
                                sprintf(name,"%lu",num_obj);
                                Mat_WriteNextStructField73(refs_id,
                                    fields[k*nfields+l],name);
                                sprintf(name,"/#refs#/%lu",num_obj);
                                H5Rcreate(refs[l]+k,*(hid_t*)mat->fp,name,
                                          H5R_OBJECT,-1);
                            }
                        }
                        for ( k = 0; k < matvar->rank; k++ )
                            perm_dims[k] = matvar->dims[matvar->rank-k-1];

                        mspace_id=H5Screate_simple(matvar->rank,perm_dims,NULL);
                        for ( l = 0; l < nfields; l++ ) {
                            dset_id = H5Dcreate(struct_id,
                                fields[l]->name,H5T_STD_REF_OBJ,mspace_id,
                                H5P_DEFAULT);
                            H5Dwrite(dset_id,H5T_STD_REF_OBJ,H5S_ALL,H5S_ALL,
                                H5P_DEFAULT,refs[l]);
                            H5Dclose(dset_id);
                            free(refs[l]);
                        }
                        free(refs);
                        H5Sclose(mspace_id);
                        H5Gclose(refs_id);
                    }
                }
                H5Gclose(struct_id);
            }
            break;
        }
    }
    return 0;
}
