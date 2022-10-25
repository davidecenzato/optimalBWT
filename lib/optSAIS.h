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
#define iget(i) ( (I[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define iset(i, b) I[(i)/8]=(b) ? (mask[(i)%8]|I[(i)/8]) : ((~mask[(i)%8])&I[(i)/8])
// mask vector
const unsigned char mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};  
// get the character function
#define chr(i) (cs==sizeof(int_t)?((int_t *)s)[i]:((char *)s)[i])
// put the character function negative
#define neg(i) (cs==sizeof(int_t)?((int_t *)s)[i]*=-1:((char *)s)[i]*=-1)
// get the character function
#define uchr(i) (cs==sizeof(uint_t)?((uint_t *)s)[i]:((unsigned char *)s)[i])

/** @brief Computes the suffix array SA, and SAP-array I of the strings concatenated
 *   using (implicitly) different dollars using the optSAIS algorithm
 *
 *  @param s		input concatenated string, using separators s[i]=0 and with s[n-1]=0
 *  @param SA		suffix array 
 *  @param I    SAP-array (Same as previous array) - n/8 bits
 *  @param n		string length
 *  @param ns		number of strings
 *  @param K		alphabet size
 *
 *  @return void.
 */
void optsais(char *s, int_t *SA, unsigned char *I, int_t n, int_t ns, int K);

/** @brief Computes the suffix array SA, and SAP-array I of the strings concatenated
 *   using (implicitly) different dollars
 *
 *  @param s		input concatenated string, using separators s[i]=0 and with s[n-1]=0
 *  @param SA		Suffix array 
 *  @param I        SAP-array (Same as previous array) - n/8 bits
 *  @param n		string length
 *  @param ns		number of strings
 *  @param K		alphabet size
 *
 *  @return void.
 */
void optsais_unsigned_sa_sap(char *s, uint_t *SA, unsigned char *I, uint_t n, uint_t ns, uint_t K, int m);

/** @brief Computes the suffix array SA of the strings concatenated
 *   using (implicitly) different dollars
 *
 *  @param s		input concatenated string, using separators s[i]=0 and with s[n-1]=0
 *  @param SA		Suffix array 
 *  @param n		string length
 *  @param ns		number of strings
 *  @param K		alphabet size
 *
 *  @return void.
 */
void optsais_unsigned_sa(char *s, uint_t *SA, uint_t n, uint_t ns, uint_t K, int m);

#endif