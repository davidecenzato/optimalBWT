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

#include "optSAIS.h"

// compute the head or end of each bucket
void getBuckets(uint_t *s, uint_t *bkt, uint_t n, uint_t K, size_t cs, bool end) { 
  uint_t i, sum=0;
  for(i=0; i<K; i++) bkt[i]=0; // clear all buckets
  for(i=0; i<n; i++) bkt[uchr(i)]++;  // compute the size of each bucket
  for(i=0; i<K; i++) { sum+=bkt[i]; bkt[i]= end ? sum-1 : sum-bkt[i]; }
}

// induce the L types with a scan left to right
void induceL(uint_t *SA, uint_t *s, uint_t *bkt, uint_t n, 
             uint_t K, int cs, uint_t ns) { 
  uint_t i, j;
  // find heads of buckets
  getBuckets((uint_t*)s, bkt, n, K, cs, false);
  // scan the suffix array left to right
  // scan the dollar characters
  for(i=0; i<ns; ++i){
    j = SA[i];
    if(j != UEMPTY){
      SA[bkt[uchr(j-1)]++] = j-1;
    }
  }
  // scan all other positions
  for(i=ns; i<n; ++i){
    j = SA[i];
    // if the position is not UEMPTY
    if( j != UEMPTY ){
      // if it is not the first position
      // of the text
      if(j > 0){
        if( uchr(j-1) >= uchr(j) ){
          SA[bkt[uchr(j-1)]++] = j-1;
          SA[i] = UEMPTY;
        } 
      }
    }
  }
}

// induce the S types with a scan right to left
void induceS(uint_t *SA, uint_t *s, uint_t *bkt, uint_t n, 
             uint_t K, int cs, unsigned char sep, uint_t ns) { 
  uint_t i, j, y, x;
  // find tails of buckets
  getBuckets((uint_t*)s, bkt, n, K, cs, true);
  // scan the suffix array right to left
  for(i=0; i<n-ns; ++i){
    j = n-1-i;
    y = SA[j];
    // if the position is not UEMPTY
    if( y != UEMPTY ){
      // if it is not the first position
      // of the text
      if(y > 0){
        x = uchr(y-1);
        if( x <= uchr(y) && x != sep ){
          assert( bkt[x]<j );
          SA[bkt[x]--] = y-1;
          SA[j] = UEMPTY;
        }
        else if( x == sep ){ SA[j] = UEMPTY; } 
      }
      else{ SA[j] = UEMPTY; }
    }
  }
}

// induce the L types and the SAP-array with a scan left to right
void induceL_s2(uint_t *SA, uint_t *s, uint_t *bkt, uint_t n, uint_t K,
                int cs, unsigned char sep) {
  uint_t i, j, y;
  // find heads of buckets
  getBuckets((uint_t*)s, bkt, n, K, cs, false);
  // scan the suffix array
  for(i=0; i<n-1; ++i){
    // check if we are in a new interval
    j = SA[i];
    // if the SA value is not UEMPTY
    if( j != UEMPTY && j > 0 ) {
      y = uchr(j-1);
      if( (y >= uchr(j)) && y != sep ){
        SA[bkt[y]++] = j-1;
      }
    }
  }
}

// induce the L types and the SAP-array with a scan right to left
void induceS_s2(uint_t *SA, uint_t *s, uint_t *bkt, uint_t n, uint_t K,
                         int cs, unsigned char sep, uint_t ns) {
  uint_t i, j, c, x;
  // find tails of buckets
  getBuckets((uint_t*)s, bkt, n, K, cs, true);
  // scan the suffix array right to left
  for(i=0; i < n-ns; ++i){
    j = n-1-i;
    // if the SA value is not UEMPTY
    x = SA[j];
    if( x != UEMPTY && x > 0 ) {
      c = uchr(x-1);
      if( ( c <= uchr(x) && bkt[c] < j ) && c != sep ){
        SA[bkt[c]--]=x-1;
      } 
    }
  }
}

// induce the L types and the SAP-array with a scan left to right
void induceL_SAP(uint_t *SA, unsigned char *SAP, uint_t *s, uint_t *bkt, 
                 uint_t n, uint_t K, int cs, unsigned char sep) {
  uint_t i, j, y;
  uint_t id=0;
  bool in=false, x = iget(0), h;
  // find heads of buckets
  getBuckets((uint_t*)s, bkt, n, K, cs, false);
  // buckets storing the intervals ids
  uint_t *bkt_int = (uint_t *)malloc(sizeof(uint_t)*K);
  for(i=0; i<K; ++i){ bkt_int[i] = UEMPTY; }
  // scan the suffix array
  for(i=0; i<n-1; ++i){
    // check if we are in a new interval
    h = iget(i+1);
    if( !x ){  if(h){ in=true; } else { in=false; } id++; }
    j = SA[i];
    // if the SA value is not UEMPTY
    if( j != UEMPTY && j > 0 ) {
      y = uchr(j-1);
      if( (y >= uchr(j)) && y != sep){
        SA[bkt[y]] = j-1;
        // if we are inside and interval
        if(in){ 
          if( bkt_int[y] == id ){ iset(bkt[y],1); }
          else{ bkt_int[y] = id; iset(bkt[y],0); }
        }
        // else just put a 0
        else { iset(bkt[y],0); }
        bkt[y]++;
      }
    }
    // update the current type
    x = h;
  }
}

// induce the L types and the SAP-array with a scan right to left
void induceS_SAP(uint_t *SA, unsigned char *SAP, uint_t *s, uint_t *bkt, 
                 uint_t n, uint_t K, int cs, unsigned char sep, uint_t ns) {
  uint_t i, j, c, x;
  uint_t id = 0;
  bool in, y = iget(n-1), h;
  in = y?true:false;
  // find tails of buckets
  getBuckets((uint_t*)s, bkt, n, K, cs, true);
  // buckets for interval index
  uint_t *bkt_int = (uint_t *)malloc(sizeof(uint_t)*K);
  for(i=0; i<K; ++i){ bkt_int[i] = UEMPTY; }
  // scan the suffix array right to left
  for(i=0; i < n-ns; ++i){
    j = n-1-i;
    h = iget(j-1);
    // if the SA value is not UEMPTY
    x = SA[j];
    if( x != UEMPTY && x > 0 ) {
      c = uchr(x-1);
      if( ( c <= uchr(x) && bkt[c] < j ) && c != sep ){
        SA[bkt[c]]=x-1;
        if(in){ 
          if( bkt_int[c] == id ){ iset(bkt[c],0); iset(bkt[c]+1,1); }
          else{ bkt_int[c] = id ; iset(bkt[c],0); }
        }
        else{ iset(bkt[c],0); }
        bkt[c]--;
      } 
    }
    if ( !y ){  if(h){ in=true; } else { in=false; } id++;  }
    // update the current type
    y = h;
  }
}

// compute length of a LMS substring
uint_t getLengthOfLMS(uint_t *s, uint_t n, uint_t x, unsigned char sep, int cs) {
  // return the length of the LMS substring
  if(x==n-1 || chr(x)==sep) return 1;  
  // S positions
  uint_t dist=0, i=1;  
  while(1) {                       
    if( uchr(x+i-1) > uchr(x+i) ) break; 
    i++;
  }  
  // L positions
  while(1) {
    if(x+i >  n-1 || uchr(x+i-1) < uchr(x+i) ) break;
    if(x+i == n-1 || uchr(x+i-1) > uchr(x+i) ) dist=i;
    i++;
  }
  
  return dist+1;
}


void optSAIS(uint_t *s, uint_t *SA, unsigned char *SAP, uint_t n, uint_t ns,
             uint_t K, int cs, unsigned char sep) {
  // initialize variables
  uint_t i, j, x, o;
  uint_t g = 0; // singletons counter
  bool y, t;

  // stage 1: reduce the problem by at least 1/2
  uint_t *bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
  // initialize SA values to -1
  for(i=0; i<n; ++i) SA[i]=UEMPTY; 
  // find ends of buckets 
  getBuckets(s, bkt, n, K, cs, true); 
  // place S* suffixes in their buckets
  // process the last character
  x = uchr(n-1); o = uchr(n-2); t = true;
  if(o > x){ SA[bkt[x]--]=n-1; }
  else{ g++; }
  // process all other characters
  for(i=1; i<=n-2; ++i){
    j = n-1-i; 
    // compute the current type
    t = ( ( o<x || (o==x && t) )?1:0 );
    x = o;
    o = uchr(j-1);
    // if the current type is S
    if(t){
      if( o > x ) { SA[bkt[x]--]=j; }
      else if( x==sep ) { g++; }
    }
  }
  // if the first character is a separator it is a singleton
  if( uchr(0)==sep ) g++;

  // sort LMS-substrings
  induceL(SA, s, bkt, n, K, cs, ns); 
  induceS(SA, s, bkt, n, K, cs, sep, ns); 

  // free bucket vector
  free(bkt); 
  // compact all the sorted substrings into the first n1 items of s
  uint_t n1=0;
  for(i=0; i<n; ++i){
      uint_t pos = SA[i];
      if(pos != UEMPTY){ SA[n1++] = pos; }
  }
  // Init the name array buffer
  for(i=n1; i<n; i++) SA[i]=UEMPTY;

  // find the lexicographic names of all substrings
  uint_t name=1, name2=1;
  if(n1 > 0){
    // insert the first LMS substring
    uint_t pos = SA[0];
    SA[n1+((pos%2==0)?pos/2:(pos-1)/2)] = name-1; 
    iset(n1+((pos%2==0)?pos/2:(pos-1)/2),0);
    // compute first LMS length and pos
    uint_t len = getLengthOfLMS(s, n, pos, sep, cs);
    uint_t prev = pos, pre_len = len;
    uint_t ppos, pprev;
    for(i=1; i<n1; i++) {
      pos=SA[i]; bool diff=false, same=false;
      len = getLengthOfLMS(s, n, pos, sep, cs);
      // if the LMS length are different skip and increase name counter
      if(len != pre_len){ diff = true; }
      else{
        // if same length compare all the characters till you find a mismatch
        ppos = pos, pprev = prev;
        for( j=0; j<len; ++j ) { 
          if( uchr(pprev) != uchr(ppos) ){ diff=true; break; }
          else{
            if( uchr(pprev) == sep && uchr(ppos) == sep ) { same=true; break; }
            pprev++; ppos++;
          }
        }
      }
      // for name2 all identical dollars, for name unique dollars
      if(diff){
        name2++; name++;
        prev=pos; pre_len=len;
      }
      else if (same) name2++;
    // insert new name
    pos=(pos%2==0)?pos/2:(pos-1)/2;
    SA[n1+pos]=name-1; 
    if(same) iset((n1+pos),1); 
    }
  }
  // compact SA and I vectors
  for(i=n-1, j=n-1; i>=n1; i--){
    if(SA[i]!=UEMPTY){
      SA[j]=SA[i];
      iset(j,iget(i));
      j--;
    }
  }
  // s1 is done now
  uint_t *SA1=SA, *s1=SA+n-n1;
  // set the new offset for SAP
  o = (n-n1);

  // stage 2: solve the reduced problem
  // recurse if names are not yet unique
  if(name2<n1) {
    optSAIS((uint_t*)s1, SA1, SAP, n1, (ns-g), name, sizeof(uint_t), 0);
  } else { // generate the suffix array of s1 directly
    bkt = (uint_t *)malloc(sizeof(uint_t)*name);
    getBuckets((uint_t*)s1, bkt, n1, name, sizeof(uint_t), false); // find beginning of buckets
    // compute the SA and SAP-array of s1 and place it at the beginning of SA
    for(i=0; i<n1; ++i){
      iset(bkt[s1[i]],iget(i+o));
      SA1[bkt[s1[i]]++] = i;
    }
    free(bkt);
  }

  // stage 3: induce the result for the original problem

  // initialize the buckets for s
  bkt = (uint_t *)malloc(sizeof(uint_t)*K); 
  // compute the end of each s bucket 
  getBuckets((uint_t*)s, bkt, n, K, cs, true);
  // compute the original LMS positions of s and put them in s1
  uint_t c = 0;
  x = uchr(n-1); o = uchr(n-2); t = true;
  if(o > x){ s1[n1-1-(c++)]=n-1; } 
  // process all other characters
  for(i=1; i<=n-2; ++i){
    j = n-1-i; // n-1 -n+2 = +1
    // compute the current type
    t = ( ( o<x || (o==x && t) )?1:0 );
    x = o;
    o = uchr(j-1);
    // if the current type is S
      if( t && o > x ) { s1[n1-1-(c++)]=j; }
  }
  assert(c==n1);
  // map back procedure
  for(i=0; i<n1; ++i){ SA1[i]=s1[SA1[i]]; } 
  // clean the suffix array
  for(i=n1; i<n; ++i){ SA[i]=UEMPTY; iset(i,0); } 
  // insert sorted LMS in the correct suffix array buckets
  for(i=0; i < n1-(ns-g); ++i) {
    j = n1-1-i;
      // take SA[j] and I[j] values
    x=SA[j]; SA[j]=UEMPTY;
    y=iget(j); iset(j,0);
    // put values in the correct positions
    iset(bkt[uchr(x)],y);
    SA[bkt[uchr(x)]--]=x;
  }
  // handle the singletons if we have some
  if(g > 0){
      j=0;
      // insert the superators at the beginning of the sa
      for(i=0;i<n;++i){ if(uchr(i)== sep){ SA[j]=i; iset(j,1); j++; } }
      iset(0,0);
  }

  // sort LMS-substrings and compute SAP-array   
  induceL_SAP(SA, SAP, s, bkt, n, K, cs, sep);
  induceS_SAP(SA, SAP, s, bkt, n, K, cs, sep, ns);

  free(bkt); 
}

void optSAIS_SA(uint_t *s, uint_t *SA, uint_t n, uint_t ns,
                uint_t K, int cs, unsigned char sep) {
  // initialize variables
  uint_t i, j, x, o;
  uint_t g = 0; // singletons counter
  bool t;

  // stage 1: reduce the problem by at least 1/2
  uint_t *bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
  // initialize SA values to -1
  for(i=0; i<n; ++i) SA[i]=UEMPTY; 
  // find ends of buckets 
  getBuckets(s, bkt, n, K, cs, true); 
  // place S* suffixes in their buckets
  // process the last character
  x = uchr(n-1); o = uchr(n-2); t = true;
  if(o > x){ SA[bkt[x]--]=n-1; }
  else{ g++; }
  // process all other characters
  for(i=1; i<=n-2; ++i){
    j = n-1-i; 
    // compute the current type
    t = ( ( o<x || (o==x && t) )?1:0 );
    x = o;
    o = uchr(j-1);
    // if the current type is S
    if(t){
      if( o > x ) { SA[bkt[x]--]=j; }
      else if( x==sep ) { g++; }
    }
  }
  // if the first character is a separator it is a singleton
  if( uchr(0)==sep ) g++;

  // sort LMS-substrings
  induceL(SA, s, bkt, n, K, cs, ns); 
  induceS(SA, s, bkt, n, K, cs, sep, ns); 

  // free bucket vector
  free(bkt); 
  // compact all the sorted substrings into the first n1 items of s
  uint_t n1=0;
  for(i=0; i<n; ++i){
      uint_t pos = SA[i];
      if(pos != UEMPTY){ SA[n1++] = pos; }
  }
  // Init the name array buffer
  for(i=n1; i<n; i++) SA[i]=UEMPTY;

  // find the lexicographic names of all substrings
  uint_t name=1, name2=1;
  if(n1 > 0){
    // insert the first LMS substring
    uint_t pos = SA[0];
    SA[n1+((pos%2==0)?pos/2:(pos-1)/2)] = name-1; 
    // compute first LMS length and pos
    uint_t len = getLengthOfLMS(s, n, pos, sep, cs);
    uint_t prev = pos, pre_len = len;
    uint_t ppos, pprev;
    for(i=1; i<n1; i++) {
      pos=SA[i]; bool diff=false, same=false;
      len = getLengthOfLMS(s, n, pos, sep, cs);
      // if the LMS length are different skip and increase name counter
      if(len != pre_len){ diff = true; }
      else{
        // if same length compare all the characters till you find a mismatch
        ppos = pos, pprev = prev;
        for( j=0; j<len; ++j ) { 
          if( uchr(pprev) != uchr(ppos) ){ diff=true; break; }
          else{
            if( uchr(pprev) == sep && uchr(ppos) == sep ) { same=true; break; }
            pprev++; ppos++;
          }
        }
      }
      // for name2 all identical dollars, for name unique dollars
      if(diff){
        name2++; name++;
        prev=pos; pre_len=len;
      }
      else if (same) name2++;
    // insert new name
    pos=(pos%2==0)?pos/2:(pos-1)/2;
    SA[n1+pos]=name-1; 
    }
  }
  // compact SA and I vectors
  for(i=n-1, j=n-1; i>=n1; i--){
    if(SA[i]!=UEMPTY){
      SA[j--]=SA[i];
    }
  }
  // s1 is done now
  uint_t *SA1=SA, *s1=SA+n-n1;

  // stage 2: solve the reduced problem
  // recurse if names are not yet unique
  if(name2<n1) {
    optSAIS_SA((uint_t*)s1, SA1, n1, (ns-g), name, sizeof(uint_t), 0);
  } else { // generate the suffix array of s1 directly
    bkt = (uint_t *)malloc(sizeof(uint_t)*name);
    getBuckets((uint_t*)s1, bkt, n1, name, sizeof(uint_t), false); // find beginning of buckets
    // compute the SA and SAP-array of s1 and place it at the beginning of SA
    for(i=0; i<n1; ++i){
      SA1[bkt[s1[i]]++] = i;
    }
    free(bkt);
  }

  // stage 3: induce the result for the original problem

  // initialize the buckets for s
  bkt = (uint_t *)malloc(sizeof(uint_t)*K); 
  // compute the end of each s bucket 
  getBuckets((uint_t*)s, bkt, n, K, cs, true);
  // compute the original LMS positions of s and put them in s1
  uint_t c = 0;
  x = uchr(n-1); o = uchr(n-2); t = true;
  if(o > x){ s1[n1-1-(c++)]=n-1; } 
  // process all other characters
  for(i=1; i<=n-2; ++i){
    j = n-1-i; // n-1 -n+2 = +1
    // compute the current type
    t = ( ( o<x || (o==x && t) )?1:0 );
    x = o;
    o = uchr(j-1);
    // if the current type is S
      if( t && o > x ) { s1[n1-1-(c++)]=j; }
  }
  assert(c==n1);
  // map back procedure
  for(i=0; i<n1; ++i){ SA1[i]=s1[SA1[i]]; }
  // clean the suffix array
  for(i=n1; i<n; ++i){ SA[i]=UEMPTY; } 
  // insert sorted LMS in the correct suffix array buckets
  for(i=0; i < n1-(ns-g); ++i) {
    j = n1-1-i;
      // take SA[j] and I[j] values
    x=SA[j]; SA[j]=UEMPTY;
    // put values in the correct positions
    SA[bkt[uchr(x)]--]=x;
  }
  // handle the singletons if we have some
  if(g > 0){
      j=0;
      // insert the superators at the beginning of the sa
      for(i=0;i<n;++i){ if(uchr(i)== sep){ SA[j++]=i; } }
  }
  // sort LMS-substrings and compute SAP-array   
  induceL_s2(SA, s, bkt, n, K, cs, sep);
  induceS_s2(SA, s, bkt, n, K, cs, sep, ns);

  free(bkt); 
}

/** @brief Computes the generalized suffix array SA, and the SAP-array of 
 *  the strings concatenated using (implicitly) different dollars
 *
 *  @param s     input: strings concatenated using separators either 1 or 0
 *  @param SA    generalized suffix array
 *  @param SAP   SAP-array (Same as previous array) - n/8 bits
 *  @param n     string length
 *  @param ns    number of strings
 *  @param K     alphabet size
 *
 *  @return void.
 */
void optsais_sa_sap(char *s, uint_t *SA, unsigned char *SAP, uint_t n, uint_t ns, uint_t K){
    if((s == nullptr) || (SA == nullptr) || (n == 0)) { std::cerr << "Empty input given." << std::endl; exit(-1); }
    optSAIS((uint_t *)s, (uint_t *)SA, (unsigned char *)SAP, n, ns, K, sizeof(char), 1);
}

/** @brief Computes the generalized suffix array SA of the strings concatenated
 *   using (implicitly) different dollars
 *
 *  @param s     input: strings concatenated using separators either 1 or 0
 *  @param SA    generalized suffix array
 *  @param n     string length
 *  @param ns    number of strings
 *  @param K     alphabet size
 *
 *  @return void.
 */
void optsais_sa(char *s, uint_t *SA, uint_t n, uint_t ns, uint_t K){
    if((s == nullptr) || (SA == nullptr) || (n == 0)) { std::cerr << "Empty input given." << std::endl; exit(-1); }
    optSAIS_SA((uint_t *)s, (uint_t *)SA, n, ns, K, sizeof(char), 1);
}