/*
 * File:   minisat_ext.h
 * Author: mikolas
 *
 * Created on November 29, 2010, 5:40 PM
 */
#ifndef MINISATEXT_HH
#define	MINISATEXT_HH
#define LOGIPASIR(code)                                                        \
    do {                                                                       \
        /* code */                                                             \
    } while (0)

#if USE_IPASIR
#include "ipasir_wrap.h"
#else
#include "minisat_ext_minisat.h"
#endif

#endif	/* MINISATEXT_HH */

