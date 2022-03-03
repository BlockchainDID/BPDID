#include <ctime>
#include <iostream>
#include <vector>

#define MR_PAIRING_BN 
#define SECURITY 128

#define SECURE_PARAM_T 256

#include "miracl/pairing_3.h"

typedef struct {
    G1 g1;
    G2 g2;
    G1 g;
    G1 h;
    G1 h1;
    G1 h2;
    G1 h3; 
}Pubparams;

typedef struct {
    Big d;
    Big w1;
    Big w2;
    Big w3;
    Big w4;
}APPriparam;

typedef struct {
    G1 gt;
    G2 Yap;
    G1 delta;
}APPubparam;

typedef struct {
    Big x;
    Big y;
}IdpPrikey;

typedef struct {
    G2 Xt;
    G2 Yt;
}IdpPubkey;

typedef struct{
    G1 R;
    G1 S;
    G1 T;
    G1 W;
    G1 J;
    G1 K;
    Big h;
    Big s;
    Big nV;
    Big nT;
}Signature;