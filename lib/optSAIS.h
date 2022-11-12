/*
 * Optimal SAIS (optSAIS) implementation to compute the generalized suffix array 
 * and the SAP-array of an string collection.
 * 
 * This code is adapted from https://github.com/davidecenzato/cais/blob/main/lib/cais.cpp
 * which is an implementation of the Conjugate array Induced sorting algorithm (cais), and 
 * from https://github.com/felipelouza/gsa-is/blob/master/gsacak.c which is an implementation
 * of the GSACA-K algorithm.
 *
 * void optsais_sa_sap(s, SA, SAP, n, ns, K); // compute the generalized conjugate array and the SAP-array of a string collection.
 * void optsais_sa(s, SA, n, ns, K)           // compute the generalized suffix array of a string collection.
 *
 */

#ifndef OPTSAIS_H
#define OPTSAIS_H

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <memory.h>

#include <iostream>
#include <ctime>
#include <vector>
#include <cassert>

#ifndef M64
  #define M64 0
#endif

#if M64
    typedef uint64_t uint_t;
    typedef int64_t int_t;
    #define U_MAX UINT64_MAX
    #define I_MAX INT64_MAX
#else
    typedef uint32_t uint_t;
    typedef int32_t int_t;
    #define U_MAX UINT32_MAX
    #define I_MAX INT32_MAX
#endif

// empty value
const int_t EMPTY = I_MAX;
// empty uint64_t value
const uint_t UEMPTY = U_MAX;
// function to access bitvector
#define iget(i) ( (SAP[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define iset(i, b) SAP[(i)/8]=(b) ? (mask[(i)%8]|SAP[(i)/8]) : ((~mask[(i)%8])&SAP[(i)/8])
// mask vector
const unsigned char mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};  
// get the character function
#define chr(i) (cs==sizeof(int_t)?((int_t *)s)[i]:((char *)s)[i])
// get the character function
#define uchr(i) (cs==sizeof(uint_t)?((uint_t *)s)[i]:((unsigned char *)s)[i])

/** @brief Computes the generalized suffix array SA, and the SAP-array of 
 *  the strings concatenated using (implicitly) different dollars
 *
 *  @param s	 input: strings concatenated using separators either 1 or 0
 *  @param SA	 generalized suffix array
 *  @param SAP   SAP-array (Same as previous array) - n/8 bits
 *  @param n	 string length
 *  @param ns	 number of strings
 *  @param K	 alphabet size
 *
 *  @return void.
 */
void optsais_sa_sap(char *s, uint_t *SA, unsigned char *SAP, uint_t n, uint_t ns, uint_t K);

/** @brief Computes the generalized suffix array SA of the strings concatenated
 *   using (implicitly) different dollars
 *
 *  @param s	 input: strings concatenated using separators either 1 or 0
 *  @param SA    generalized suffix array
 *  @param n	 string length
 *  @param ns	 number of strings
 *  @param K	 alphabet size
 *
 *  @return void.
 */
void optsais_sa(char *s, uint_t *SA, uint_t n, uint_t ns, uint_t K);

#endif