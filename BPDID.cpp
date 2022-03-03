#include "BPDID.h"
#include <time.h>

#define ITERATIONS 100

void PS_setup(PFC *pfc, G1 *g1, G2 *g2){
    cout << "PS signature setup  "  << endl;
    time_t seed;
    time(&seed);
    irand((long)seed);

    pfc->random(*g1);
    pfc->random(*g2);
}


void PS_keygen(PFC *pfc, G1 *g1, G2 *g2, Big *x, Big *y, G2 *X2, G2 *Y2){
    cout << "PS signature keygen"  << endl;
    pfc->random(*x);
    pfc->random(*y);
    *x = *x % *pfc->ord;
    *y = *y % *pfc->ord;

    *X2 = pfc->mult(*g2, *x);
    *Y2 = pfc->mult(*g2, *y);  
}

void PS_sign(PFC *pfc, G1 *g1, G2 *g2, Big *x, Big *y, Big m, G1 *A, G1 *B){
    //cout << "PS signature sign  "  << endl;
    Big u, tmp;
    G2 tmp2;

    pfc->random(u);
    u = u % *pfc->ord;

    *A = pfc->mult(*g1, u);

    tmp = *x + modmult(*y, m, *pfc->ord);
    tmp = tmp % *pfc->ord;
    //tmp = (*x + (*y) * (m));

    *B = pfc->mult(*A, tmp);
}

int PS_verify(PFC *pfc, G1 *g1, G2 *g2, G2 *X2, G2 *Y2, Big m, G1 *A, G1 *B){
    cout << "PS signature verify  "  << endl;
    G2 tmp;
    GT T1,T2;

    tmp = (*X2) + pfc->mult(*Y2, m);

    T1 = pfc->pairing(tmp, *A);
    T2 = pfc->pairing(*g2, *B);
    if(T1!= T2){
        return 0;
    }
    return 1;
}

void rangeproof_step1(PFC *pfc, G1 *g1, G2 *g2, G1 *g, G1 *h, G2 *Y2, G1 *A, G1 *B, G1 *At, G1 *Bt, G1 *D, GT *E, Big *s, Big *t, Big *rho1, Big *rho2, Big *rho3){
    cout << "PS range proof step1  "  << endl;
    //Big s, t, rho1, rho2, rho3;
    G1 tmp;
    GT tmp2, tmp3;
    
    pfc->random(*s);
    pfc->random(*t);
    pfc->random(*rho1);
    pfc->random(*rho2);
    pfc->random(*rho3);

    *At = pfc->mult(*A, *t);
    tmp = (*B) + pfc->mult(*A, *s);
    *Bt = pfc->mult(tmp, *t);

    *D = pfc->mult(*g, *rho1) + pfc->mult(*h, *rho2);

    tmp2 = pfc->pairing(*Y2, *At);
    tmp3 = pfc->pairing(*g2, *At);
    *E = pfc->power(tmp2, *rho1) * pfc->power(tmp3, *rho3);
}

void rangeproof_step2(PFC *pfc, Big *r, Big *m, Big *s, Big *t, Big *rho1, Big *rho2, Big *rho3, Big *c, Big *s1, Big *s2, Big *s3){
    cout << "PS range proof step2  "  << endl;
    *s1 = *rho1 - modmult(*c , *m, *pfc->ord);
    *s1 = *s1 % *pfc->ord;

    *s2 = *rho2 - modmult(*c , *r, *pfc->ord);
    *s2 = *s2 % *pfc->ord;

    *s3 = *rho3 - modmult(*c , *s, *pfc->ord);
    *s3 = *s3 % *pfc->ord;
}

int rangeproof_step3(PFC *pfc, G2 *g2, G1 *g, G1 *h, G2 *X2, G2 *Y2, G1 *At, G1 *Bt, Big *c, Big *s1, Big *s2, Big *s3, G1 *C, G1 *D, GT *E){
    cout << "PS range proof step3  "  << endl;
    GT Et, tmp1, tmp2, tmp3;
    G1 Dt;

    Dt = pfc->mult(*g, *s1) + pfc->mult(*h, *s2) + pfc->mult(*C, *c);
    if (Dt != *D){
        cout << "D verify failed!" << endl;
        return 0;
    }

    tmp1 = pfc->pairing(*Y2, *At);
    tmp2 = pfc->pairing(*g2, *At);
    tmp3 = pfc->pairing(*g2, *Bt) / pfc->pairing(*X2, *At); 

    Et = pfc->power(tmp1, *s1) * pfc->power(tmp2, *s3) * pfc->power(tmp3, *c);

    if(Et != *E){
        cout << "E verify failed!" << endl;
        return 0;
    }
    return 1;
}

void Setup(PFC *pfc, Pubparams *params){
    //cout << "Setup"  << endl;
    time_t seed;
    time(&seed);
    irand((long)seed);


    pfc->random(params->g1);
    pfc->random(params->g2);
    pfc->random(params->g);
    pfc->random(params->h);
    pfc->random(params->h1);
    pfc->random(params->h2);
    pfc->random(params->h3);
}

void AP_keygen(PFC *pfc, Pubparams *Pubparams, APPriparam *apPriparam, APPubparam *apPubparam){
    //cout << "AP keygen"  << endl;

    Big tmp;

    pfc->random(apPriparam->d);
    apPriparam->d = apPriparam->d % *pfc->ord;
    pfc->random(apPriparam->w1);
    pfc->random(apPriparam->w2);
    pfc->random(apPriparam->w3);
    pfc->random(apPriparam->w4);
    pfc->random(apPubparam->gt);

    apPubparam->Yap = pfc->mult(Pubparams->g2, apPriparam->d);

    tmp = modmult((apPriparam->d + apPriparam->w1), (apPriparam->d + apPriparam->w2), *pfc->ord);
    tmp = modmult(tmp, (apPriparam->d + apPriparam->w3), *pfc->ord);
    tmp = modmult(tmp, (apPriparam->d + apPriparam->w4), *pfc->ord);

    apPubparam->delta = pfc->mult(apPubparam->gt, tmp);
}

void IdP_keygen(PFC *pfc, Pubparams *Pubparams, IdpPrikey *Idpprikey, IdpPubkey *Idppubkey){
    //cout << "IdP keygen"  << endl;
    
    pfc->random(Idpprikey->x);
    pfc->random(Idpprikey->y);
    Idpprikey->x = Idpprikey->x % *pfc->ord;
    Idpprikey->y = Idpprikey->y % *pfc->ord;
    
    Idppubkey->Xt = pfc->mult(Pubparams->g2,Idpprikey->x);
    Idppubkey->Yt = pfc->mult(Pubparams->g2, Idpprikey->y);
}

void ZkpProof1_step1(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, G1 *D, G2 *E, GT *F, Big *rho1, Big *rho2, Big *rho3, Big wu, Big a1, G1 *Wu, G1 *T, G2 *a, G1 *b, Big *rw, Big *k){

    //cout << "Zkp proof 1 step1  "  << endl;
    GT tmp, tmp2, tmp3;
    
    pfc->random(*rho1);
    pfc->random(*rho2);
    pfc->random(*rho3);
    pfc->random(*rw);

    *k = modmult(*rw, wu, *pfc->ord);
    *a = pfc->mult(apPubparam->Yap, *rw);
    *b = *Wu + pfc->mult(Pubparams->g1, *rw);

    *D = pfc->mult(*T, *rho3) + pfc->mult(Pubparams->h1, a1);
    *E = pfc->mult(apPubparam->Yap, *rho1);

    tmp = pfc->pairing(*a,Pubparams->g1);
    tmp2 = pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *rho2);
    tmp3 = pfc->power(pfc->pairing(Pubparams->g2, *b), *rho3);
    *F = tmp * tmp2 / tmp3;
}

void ZkpProof1N_step1(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, G1 *D, G2 *E, GT *F, Big *rho1, Big *rho2, Big *rho3, Big wu, Big a1, Big a2, G1 *Wu, G1 *T, G2 *a, G1 *b, Big *rw, Big *k){

    //cout << "Zkp proof 1 step1  "  << endl;
    GT tmp, tmp2, tmp3;
    
    pfc->random(*rho1);
    pfc->random(*rho2);
    pfc->random(*rho3);
    pfc->random(*rw);

    *k = modmult(*rw, wu, *pfc->ord);
    *a = pfc->mult(apPubparam->Yap, *rw);
    *b = *Wu + pfc->mult(Pubparams->g1, *rw);

    *D = pfc->mult(*T, *rho3) + pfc->mult(Pubparams->h1, a1) + pfc->mult(Pubparams->h2, a2);
    *E = pfc->mult(apPubparam->Yap, *rho1);

    tmp = pfc->pairing(*a,Pubparams->g1);
    tmp2 = pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *rho2);
    tmp3 = pfc->power(pfc->pairing(Pubparams->g2, *b), *rho3);
    *F = tmp * tmp2 / tmp3;
}

void ZkpProof1_step2(PFC *pfc,Big rw, Big k, Big wu, Big c, Big rho1, Big rho2, Big rho3, Big *s1, Big *s2, Big *s3){
    //cout << "Zkp proof 1 step2  "  << endl;
    *s1 = rho1 - modmult(c , rw, *pfc->ord);
    *s1 = *s1 % *pfc->ord;

    *s2 = rho2 - modmult(c , k, *pfc->ord);
    *s2 = *s2 % *pfc->ord;

    *s3 = rho3 - modmult(c , wu, *pfc->ord);
    *s3 = *s3 % *pfc->ord;
}

int ZkpProof1_step3(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, G1 D, G2 E, GT F, Big c, Big s1, Big s2, Big s3,G1 T, G1 tau,G2 a, G1 b, Big a1){
    //cout << "Zkp proof 1 step3  "  << endl;
    G1 Dt;
    G2 Et;
    GT Ft;
    Big tmp;

    tmp = a1 - modmult(a1 , c, *pfc->ord);
    tmp = tmp % *pfc->ord;
    Dt = pfc->mult(T, s3) + pfc->mult(Pubparams->h1, tmp) + pfc->mult(tau, c);
    if (Dt != D){
        cout << "D verify failed!" << endl;
        return 0;
    }

    Et = pfc->mult(apPubparam->Yap, s1) + pfc->mult(a, c);
    if(Et != E){
        cout << "E verify failed!" << endl;
        return 0;
    }

    tmp = (1 - c) % *pfc->ord;
    Ft = pfc->power(pfc->pairing(a, Pubparams->g1), tmp) * pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), s2) * pfc->power(pfc->pairing(apPubparam->Yap, b) / pfc->pairing(Pubparams->g2, apPubparam->delta), c) / pfc->power(pfc->pairing(Pubparams->g2, b), s3);
    if(Ft != F){
        cout << "F verify failed!" << endl;
        return 0;
    }
    return 1;
}

int ZkpProof1N_step3(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, G1 *D, G2 *E, GT *F, Big *c, Big *s1, Big *s2, Big *s3,G1 *T, G1 *tau,G2 *a, G1 *b, Big *a1, Big *a2){
    //cout << "Zkp proof 1 step3  "  << endl;
    G1 Dt;
    G2 Et;
    GT Ft;
    Big tmp, tmp2;

    tmp = *a1 - modmult(*a1 , *c, *pfc->ord);
    tmp = tmp % *pfc->ord;
    tmp2 = *a2 - modmult(*a2 , *c, *pfc->ord);
    tmp2 = tmp2 % *pfc->ord;
    Dt = pfc->mult(*T, *s3) + pfc->mult(Pubparams->h1, tmp) +pfc->mult(Pubparams->h2, tmp2) + pfc->mult(*tau, *c);
    if (Dt != *D){
        cout << "D verify failed!" << endl;
        return 0;
    }

    Et = pfc->mult(apPubparam->Yap, *s1) + pfc->mult(*a, *c);
    if(Et != *E){
        cout << "E verify failed!" << endl;
        return 0;
    }

    tmp = (1 - *c) % *pfc->ord;
    Ft = pfc->power(pfc->pairing(*a, Pubparams->g1), tmp) * pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *s2) * pfc->power(pfc->pairing(apPubparam->Yap, *b) / pfc->pairing(Pubparams->g2, apPubparam->delta), *c) / pfc->power(pfc->pairing(Pubparams->g2, *b), *s3);
    if(Ft != *F){
        cout << "F verify failed!" << endl;
        return 0;
    }
    return 1;
}

void Register(PFC *pfc, APPriparam *apPriparam, APPubparam *apPubparam, Big *wu, G1 *Wu){
    //cout << "Register"  << endl;
    Big tmp;
    *wu = apPriparam->w1;
    tmp = inverse(*wu + apPriparam->d, *pfc->ord);

    *Wu = pfc->mult(apPubparam->delta, tmp);
}

void RequestCred1(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, Big wu, G1 *Wu, Big a1, G1 *T, G1 *tau){
    //cout << "Request Cred 1"  << endl;
    Big t;
    Big rho1, rho2, rho3, s1, s2, s3;
    Big c,k, rw;
    G2 a, E;
    G1 b, D;
    GT F;
    int result;

    pfc->random(t);
    *T = pfc->mult(Pubparams->g, t);
    *tau = pfc->mult(*T, wu) + pfc->mult(Pubparams->h1, a1);

    ZkpProof1_step1(pfc, Pubparams, apPubparam, &D, &E, &F, &rho1, &rho2, &rho3, wu, a1, Wu, T, &a, &b, &rw, &k );
    pfc->random(c);
    ZkpProof1_step2(pfc, rw, k, wu, c, rho1, rho2, rho3, &s1, &s2, &s3);

    result = ZkpProof1_step3(pfc, Pubparams, apPubparam, D, E, F, c, s1, s2, s3, *T, *tau, a, b, a1);

    //cout << "RequestCred1 proof verify result = " << result << endl;
}

void RequestCredN(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, Big wu, G1 *Wu, Big a1, Big a2, G1 *T, G1 *tau){
    //cout << "Request Cred N"  << endl;
    Big t;
    Big rho1, rho2, rho3, s1, s2, s3;
    Big c,k, rw;
    G2 a, E;
    G1 b, D;
    GT F;
    int result;

    pfc->random(t);
    *T = pfc->mult(Pubparams->g, t);
    *tau = pfc->mult(*T, wu) + pfc->mult(Pubparams->h1, a1) + pfc->mult(Pubparams->h2, a2);

    ZkpProof1N_step1(pfc, Pubparams, apPubparam, &D, &E, &F, &rho1, &rho2, &rho3, wu, a1, a2, Wu, T, &a, &b, &rw, &k );
    pfc->random(c);
    ZkpProof1_step2(pfc, rw, k, wu, c, rho1, rho2, rho3, &s1, &s2, &s3);

    result = ZkpProof1N_step3(pfc, Pubparams, apPubparam, &D, &E, &F, &c, &s1, &s2, &s3, T, tau, &a, &b, &a1, &a2);

    //cout << "RequestCredN proof verify result = " << result << endl;
}

void IssueCred1(PFC *pfc, Pubparams *pubparams, IdpPrikey *idpprikey, Big a1, G1 *A, G1 *B){
    //cout << "Issue Cred 1"  << endl;
    PS_sign(pfc, &pubparams->g1, &pubparams->g2, &idpprikey->x, &idpprikey->y, a1, A, B);
}

void ZkpProof2_step1(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, G1 *A, G1 *B, G1 *At, G1 *Bt, G1 *D, G2 *E, GT *F,GT *M, Big *rho1, Big *rho2, Big *rho3, Big *rho4, Big *rho5, Big wu, Big a1, G1 *Wu, G1 *T, G2 *a, G1 *b, Big *rw, Big *k, Big *s){
    //cout << "Zkp proof 2 step1  "  << endl;
    GT tmp, tmp2, tmp3;
    Big t;
    G1 tmpp;
    GT tmpp2, tmpp3;
    
    pfc->random(*rho1);
    pfc->random(*rho2);
    pfc->random(*rho3);
    pfc->random(*rho4);
    pfc->random(*rho5);
    pfc->random(*rw);
    pfc->random(*s);
    pfc->random(t);

    *k = modmult(*rw, wu, *pfc->ord);
    *a = pfc->mult(apPubparam->Yap, *rw);
    *b = *Wu + pfc->mult(Pubparams->g1, *rw);

    *D = pfc->mult(*T, *rho3) + pfc->mult(Pubparams->h1, *rho4);
    *E = pfc->mult(apPubparam->Yap, *rho1);

    tmp = pfc->pairing(*a,Pubparams->g1);
    tmp2 = pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *rho2);
    tmp3 = pfc->power(pfc->pairing(Pubparams->g2, *b), *rho3);
    *F = tmp * tmp2 / tmp3;

    *At = pfc->mult(*A, t);
    tmpp = (*B) + pfc->mult(*A, *s);
    *Bt = pfc->mult(tmpp, t);

    tmpp2 = pfc->pairing(idppubkey->Yt, *At);
    tmpp3 = pfc->pairing(Pubparams->g2, *At);
    *M = pfc->power(tmpp2, *rho4) * pfc->power(tmpp3, *rho5);
}


void ZkpProof2_step2(PFC *pfc,Big *rw, Big *k, Big *wu, Big *a1, Big *s, Big *c, Big *rho1, Big *rho2, Big *rho3,Big *rho4, Big *rho5, Big *s1, Big *s2, Big *s3, Big *s4, Big *s5){
    //cout << "Zkp proof 2 step2  "  << endl;
    *s1 = *rho1 - modmult(*c , *rw, *pfc->ord);
    *s1 = *s1 % *pfc->ord;

    *s2 = *rho2 - modmult(*c , *k, *pfc->ord);
    *s2 = *s2 % *pfc->ord;

    *s3 = *rho3 - modmult(*c , *wu, *pfc->ord);
    *s3 = *s3 % *pfc->ord;

    *s4 = *rho4 - modmult(*c , *a1, *pfc->ord);
    *s4 = *s4 % *pfc->ord;

    *s5 = *rho5 - modmult(*c , *s, *pfc->ord);
    *s5 = *s5 % *pfc->ord;
}

int ZkpProof2_step3(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, G1 *D, G2 *E, GT *F, GT *M, Big *c, Big *s1, Big *s2, Big *s3, Big *s4, Big *s5, G1 *T, G1 *tau,G2 *a, G1 *b, G1 *At, G1 *Bt){
    //cout << "Zkp proof 2 step3  "  << endl;
    G1 Dt;
    G2 Et;
    GT Ft,Mt;
    Big tmp;

    Dt = pfc->mult(*T, *s3) + pfc->mult(Pubparams->h1, *s4) + pfc->mult(*tau, *c);
    if (Dt != *D){
        cout << "D verify failed!" << endl;
        //return 0;
    }

    Et = pfc->mult(apPubparam->Yap, *s1) + pfc->mult(*a, *c);
    if(Et != *E){
        cout << "E verify failed!" << endl;
        return 0;
    }

    tmp = (1 - *c) % *pfc->ord;
    Ft = pfc->power(pfc->pairing(*a, Pubparams->g1), tmp) * pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *s2) * pfc->power(pfc->pairing(apPubparam->Yap, *b) / pfc->pairing(Pubparams->g2, apPubparam->delta), *c) / pfc->power(pfc->pairing(Pubparams->g2, *b), *s3);
    if(Ft != *F){
        cout << "F verify failed!" << endl;
        return 0;
    }

    Mt = pfc->power(pfc->pairing(idppubkey->Yt, *At), *s4) * pfc->power(pfc->pairing(Pubparams->g2, *At), *s5) * pfc->power(pfc->pairing(Pubparams->g2, *Bt) / pfc->pairing(idppubkey->Xt, *At), *c);
    if(Mt != *M){
        cout << "M verify failed!" << endl;
        return 0;
    }

    return 1;
}


void Show(PFC *pfc, Pubparams *pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, Big wu, G1 *Wu, Big a1, G1 *T, G1 *tau, G1 *A, G1 *B){
    //cout << "Show Cred 1"  << endl;

    Big t;
    Big rho1, rho2, rho3, rho4, rho5, s1, s2, s3, s4, s5;
    Big c,k, rw, s;
    G2 a, E;
    G1 b, D, At, Bt;
    GT F, M;
    int result;

    ZkpProof2_step1(pfc, pubparams, apPubparam, idppubkey, A, B, &At, &Bt, &D, &E, &F, &M, &rho1, &rho2, &rho3, &rho4, &rho5, wu, a1, Wu, T, &a, &b, &rw, &k, &s);

    pfc->random(c);
    ZkpProof2_step2(pfc, &rw, &k, &wu, &a1, &s, &c, &rho1, &rho2, &rho3, &rho4, &rho5, &s1, &s2, &s3, &s4, &s5);

    result = ZkpProof2_step3(pfc, pubparams, apPubparam, idppubkey, &D, &E, &F, &M, &c, &s1, &s2, &s3, &s4, &s5, T, tau, &a, &b, &At, &Bt);

    //cout << "Show proof verify result = " << result << endl;
}

void ZkpProof2N_step1(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, G1 *A, G1 *B, G1 *At, G1 *Bt, G1 *A2t, G1 *B2t, G1 *D, G2 *E, GT *F, GT *M, GT *N, Big *rho1, Big *rho2, Big *rho3, Big *rho4, Big *rho5, Big *rho6, Big wu, G1 *Wu, G1 *T, G2 *a, G1 *b, Big *rw, Big *k, Big *s){
    //cout << "Zkp proof 2 step1  "  << endl;
    GT tmp, tmp2, tmp3;
    Big t;
    G1 tmpp;
    GT tmpp2, tmpp3;
    
    pfc->random(*rho1);
    pfc->random(*rho2);
    pfc->random(*rho3);
    pfc->random(*rho4);
    pfc->random(*rho5);
    pfc->random(*rho6);
    pfc->random(*rw);
    pfc->random(*s);
    pfc->random(t);

    *k = modmult(*rw, wu, *pfc->ord);
    *a = pfc->mult(apPubparam->Yap, *rw);
    *b = *Wu + pfc->mult(Pubparams->g1, *rw);

    *D = pfc->mult(*T, *rho3) + pfc->mult(Pubparams->h1, *rho4) + pfc->mult(Pubparams->h2, *rho6);
    *E = pfc->mult(apPubparam->Yap, *rho1);

    tmp = pfc->pairing(*a,Pubparams->g1);
    tmp2 = pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *rho2);
    tmp3 = pfc->power(pfc->pairing(Pubparams->g2, *b), *rho3);
    *F = tmp * tmp2 / tmp3;

    *At = pfc->mult(*A, t);
    tmpp = (*B) + pfc->mult(*A, *s);
    *Bt = pfc->mult(tmpp, t);

    tmpp2 = pfc->pairing(idppubkey->Yt, *At);
    tmpp3 = pfc->pairing(Pubparams->g2, *At);
    *M = pfc->power(tmpp2, *rho4) * pfc->power(tmpp3, *rho5);

    tmpp2 = pfc->pairing(idppubkey->Yt, *A2t);
    tmpp3 = pfc->pairing(Pubparams->g2, *A2t);
    *N = pfc->power(tmpp2, *rho6) * pfc->power(tmpp3, *rho5);


}

void ZkpProof2N_step2(PFC *pfc,Big *rw, Big *k, Big *wu, Big *a1, Big *a2, Big *s, Big *c, Big *rho1, Big *rho2, Big *rho3,Big *rho4, Big *rho5, Big *rho6, Big *s1, Big *s2, Big *s3, Big *s4, Big *s5, Big *s6){
    //cout << "Zkp proof 2 step2  "  << endl;
    *s1 = *rho1 - modmult(*c , *rw, *pfc->ord);
    *s1 = *s1 % *pfc->ord;

    *s2 = *rho2 - modmult(*c , *k, *pfc->ord);
    *s2 = *s2 % *pfc->ord;

    *s3 = *rho3 - modmult(*c , *wu, *pfc->ord);
    *s3 = *s3 % *pfc->ord;

    *s4 = *rho4 - modmult(*c , *a1, *pfc->ord);
    *s4 = *s4 % *pfc->ord;

    *s5 = *rho5 - modmult(*c , *s, *pfc->ord);
    *s5 = *s5 % *pfc->ord;

    *s6 = *rho6 - modmult(*c , *a2, *pfc->ord);
    *s6 = *s6 % *pfc->ord;
}

int ZkpProof2N_step3(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, G1 *D, G2 *E, GT *F, GT *M, GT *N, Big *c, Big *s1, Big *s2, Big *s3, Big *s4, Big *s5, Big *s6, G1 *T, G1 *tau,G2 *a, G1 *b, G1 *At, G1 *Bt, G1 *A2t, G1 *B2t){
    //cout << "Zkp proof 2 step3  "  << endl;
    G1 Dt;
    G2 Et;
    GT Ft,Mt,Nt;
    Big tmp;

    Dt = pfc->mult(*T, *s3) + pfc->mult(Pubparams->h1, *s4) + pfc->mult(Pubparams->h2, *s6)+ pfc->mult(*tau, *c);
    if (Dt != *D){
        cout << "D verify failed!" << endl;
        //return 0;
    }

    Et = pfc->mult(apPubparam->Yap, *s1) + pfc->mult(*a, *c);
    if(Et != *E){
        cout << "E verify failed!" << endl;
        return 0;
    }

    tmp = (1 - *c) % *pfc->ord;
    Ft = pfc->power(pfc->pairing(*a, Pubparams->g1), tmp) * pfc->power(pfc->pairing(Pubparams->g2, Pubparams->g1), *s2) * pfc->power(pfc->pairing(apPubparam->Yap, *b) / pfc->pairing(Pubparams->g2, apPubparam->delta), *c) / pfc->power(pfc->pairing(Pubparams->g2, *b), *s3);
    if(Ft != *F){
        cout << "F verify failed!" << endl;
        return 0;
    }

    Mt = pfc->power(pfc->pairing(idppubkey->Yt, *At), *s4) * pfc->power(pfc->pairing(Pubparams->g2, *At), *s5) * pfc->power(pfc->pairing(Pubparams->g2, *Bt) / pfc->pairing(idppubkey->Xt, *At), *c);
    if(Mt != *M){
        cout << "M verify failed!" << endl;
        return 0;
    }

    Nt = pfc->power(pfc->pairing(idppubkey->Yt, *A2t), *s6) * pfc->power(pfc->pairing(Pubparams->g2, *A2t), *s5) * pfc->power(pfc->pairing(Pubparams->g2, *B2t) / pfc->pairing(idppubkey->Xt, *A2t), *c);
    if(Nt != *N){
        cout << "N verify failed!" << endl;
        return 0;
    }

    return 1;
}

void ShowN(PFC *pfc, Pubparams *pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, Big wu, G1 *Wu, Big a1, Big a2, G1 *T, G1 *tau, G1 *A, G1 *B, G1 *A2, G1 *B2){
    cout << "Show Cred N"  << endl;

    Big t;
    Big rho1, rho2, rho3, rho4, rho5, rho6, s1, s2, s3, s4, s5, s6;
    Big c,k, rw, s;
    G2 a, E;
    G1 b, D, At, Bt, A2t, B2t;
    GT F, M, N;
    int result;

    ZkpProof2N_step1(pfc, pubparams, apPubparam, idppubkey, A, B, &At, &Bt,&A2t, &B2t, &D, &E, &F, &M, &N, &rho1, &rho2, &rho3, &rho4, &rho5, &rho6, wu, Wu, T, &a, &b, &rw, &k, &s);

    pfc->random(c);
    ZkpProof2N_step2(pfc, &rw, &k, &wu, &a1, &a2, &s, &c, &rho1, &rho2, &rho3, &rho4, &rho5, &rho6, &s1, &s2, &s3, &s4, &s5, &s6);

    result = ZkpProof2N_step3(pfc, pubparams, apPubparam, idppubkey, &D, &E, &F, &M, &N, &c, &s1, &s2, &s3, &s4, &s5, &s6, T, tau, &a, &b, &At, &Bt, &A2t, &B2t);

    //cout << "Show proof verify result = " << result << endl;
}

void RevokeUser(PFC *pfc, APPubparam *apPubparam, APPriparam *apPriparam, Big w){
    Big tmp;
    tmp = inverse(w + apPriparam->d, *pfc->ord);
    apPubparam->delta = pfc->mult(apPubparam->delta, tmp);
}


void Algorithm_test(){
    PFC pfc(SECURITY);

    Pubparams pubparams;
    APPriparam apPriparam;
    APPubparam apPubparam;
    IdpPrikey idpprikey;
    IdpPubkey idppubkey;

    Big wu, a1, a2;
    G1 Wu, T, tau;
    G1 A, B, A2, B2;

    clock_t start;
    double elapsed;
    int iterations;



    Setup(&pfc, &pubparams);

    AP_keygen(&pfc, &pubparams, &apPriparam, &apPubparam);

    IdP_keygen(&pfc, &pubparams, &idpprikey, &idppubkey);

    Register(&pfc, &apPriparam, &apPubparam, &wu, &Wu);

    pfc.random(a1);
    a1 = a1 % *pfc.ord;


    RequestCred1(&pfc, &pubparams, &apPubparam, wu, &Wu, a1, &T, &tau);

    IssueCred1(&pfc, &pubparams, &idpprikey, a1, &A, &B);

    Show(&pfc, &pubparams, &apPubparam, &idppubkey, wu, &Wu, a1, &T, &tau, &A, &B);

    pfc.random(a2);
    a2 = a2 % *pfc.ord;

    RequestCredN(&pfc, &pubparams, &apPubparam, wu, &Wu, a1, a2, &T, &tau);

    IssueCred1(&pfc, &pubparams, &idpprikey, a2, &A2, &B2);

    ShowN(&pfc, &pubparams, &apPubparam, &idppubkey, wu, &Wu, a1, a2, &T, &tau, &A, &B, &A2, &B2);


}

void RequestCred1_bmark(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, Big wu, G1 *Wu, Big a1, G1 *T, G1 *tau, clock_t *request_time, clock_t *issue_time, clock_t *start){
    //cout << "Request Cred 1"  << endl;
    Big t;
    Big rho1, rho2, rho3, s1, s2, s3;
    Big c,k, rw;
    G2 a, E;
    G1 b, D;
    GT F;
    int result;
    clock_t requestp;

    pfc->random(t);
    *T = pfc->mult(Pubparams->g, t);
    *tau = pfc->mult(*T, wu) + pfc->mult(Pubparams->h1, a1);

    ZkpProof1_step1(pfc, Pubparams, apPubparam, &D, &E, &F, &rho1, &rho2, &rho3, wu, a1, Wu, T, &a, &b, &rw, &k );
    pfc->random(c);
    ZkpProof1_step2(pfc, rw, k, wu, c, rho1, rho2, rho3, &s1, &s2, &s3);
    
    *request_time += clock() - *start;
    requestp = clock();
    result = ZkpProof1_step3(pfc, Pubparams, apPubparam, D, E, F, c, s1, s2, s3, *T, *tau, a, b, a1);
    *issue_time += clock()- requestp;
    //cout << "RequestCred1 proof verify result = " << result << endl;
}

void RequestCredN_bmark(PFC *pfc, Pubparams *Pubparams, APPubparam *apPubparam, Big wu, G1 *Wu, Big a1, Big a2, G1 *T, G1 *tau, clock_t *request_time, clock_t *issue_time, clock_t *start){
    //cout << "Request Cred N"  << endl;
    Big t;
    Big rho1, rho2, rho3, s1, s2, s3;
    Big c,k, rw;
    G2 a, E;
    G1 b, D;
    GT F;
    int result;
    clock_t requestp;

    pfc->random(t);
    *T = pfc->mult(Pubparams->g, t);
    *tau = pfc->mult(*T, wu) + pfc->mult(Pubparams->h1, a1) + pfc->mult(Pubparams->h2, a2);

    ZkpProof1N_step1(pfc, Pubparams, apPubparam, &D, &E, &F, &rho1, &rho2, &rho3, wu, a1, a2, Wu, T, &a, &b, &rw, &k );
    pfc->random(c);
    ZkpProof1_step2(pfc, rw, k, wu, c, rho1, rho2, rho3, &s1, &s2, &s3);

    *request_time += clock() - *start;
    requestp = clock();
    result = ZkpProof1N_step3(pfc, Pubparams, apPubparam, &D, &E, &F, &c, &s1, &s2, &s3, T, tau, &a, &b, &a1, &a2);
    *issue_time += clock()- requestp;
    //cout << "RequestCredN proof verify result = " << result << endl;
}

void Show_bmark(PFC *pfc, Pubparams *pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, Big wu, G1 *Wu, Big a1, G1 *T, G1 *tau, G1 *A, G1 *B, clock_t *show_time, clock_t *verify_time, clock_t *start){
    //cout << "Show Cred 1"  << endl;

    Big t;
    Big rho1, rho2, rho3, rho4, rho5, s1, s2, s3, s4, s5;
    Big c,k, rw, s;
    G2 a, E;
    G1 b, D, At, Bt;
    GT F, M;
    int result;
    clock_t showp;

    ZkpProof2_step1(pfc, pubparams, apPubparam, idppubkey, A, B, &At, &Bt, &D, &E, &F, &M, &rho1, &rho2, &rho3, &rho4, &rho5, wu, a1, Wu, T, &a, &b, &rw, &k, &s);

    pfc->random(c);
    ZkpProof2_step2(pfc, &rw, &k, &wu, &a1, &s, &c, &rho1, &rho2, &rho3, &rho4, &rho5, &s1, &s2, &s3, &s4, &s5);

    *show_time += clock()-*start;
    showp = clock();
    result = ZkpProof2_step3(pfc, pubparams, apPubparam, idppubkey, &D, &E, &F, &M, &c, &s1, &s2, &s3, &s4, &s5, T, tau, &a, &b, &At, &Bt);
    *verify_time += clock()-showp;
    //cout << "Show proof verify result = " << result << endl;
}

void ShowN_bmark(PFC *pfc, Pubparams *pubparams, APPubparam *apPubparam, IdpPubkey *idppubkey, Big wu, G1 *Wu, Big a1, Big a2, G1 *T, G1 *tau, G1 *A, G1 *B, G1 *A2, G1 *B2, clock_t *show_time, clock_t *verify_time, clock_t *start){
    //cout << "Show Cred N"  << endl;

    Big t;
    Big rho1, rho2, rho3, rho4, rho5, rho6, s1, s2, s3, s4, s5, s6;
    Big c,k, rw, s;
    G2 a, E;
    G1 b, D, At, Bt, A2t, B2t;
    GT F, M, N;
    int result;
    clock_t showp;

    ZkpProof2N_step1(pfc, pubparams, apPubparam, idppubkey, A, B, &At, &Bt,&A2t, &B2t, &D, &E, &F, &M, &N, &rho1, &rho2, &rho3, &rho4, &rho5, &rho6, wu, Wu, T, &a, &b, &rw, &k, &s);

    pfc->random(c);
    ZkpProof2N_step2(pfc, &rw, &k, &wu, &a1, &a2, &s, &c, &rho1, &rho2, &rho3, &rho4, &rho5, &rho6, &s1, &s2, &s3, &s4, &s5, &s6);

    *show_time += clock()-*start;
    showp = clock();
    result = ZkpProof2N_step3(pfc, pubparams, apPubparam, idppubkey, &D, &E, &F, &M, &N, &c, &s1, &s2, &s3, &s4, &s5, &s6, T, tau, &a, &b, &At, &Bt, &A2t, &B2t);
    *verify_time += clock()-showp;

    //cout << "Show proof verify result = " << result << endl;
}

void Algorithm_bmark(){
    PFC pfc(SECURITY);

    Pubparams pubparams;
    APPriparam apPriparam;
    APPubparam apPubparam;
    IdpPrikey idpprikey;
    IdpPubkey idppubkey;

    Big wu, a1, a2;
    G1 Wu, T, tau;
    G1 A, B, A2, B2;

    clock_t start, show_time, verify_time, request_time, issue_time;
    double elapsed;
    int iterations;

    iterations = 0;
    start=clock();
    do{
        Setup(&pfc, &pubparams);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("Setup cost %8.2lf ms\n",elapsed);

    iterations = 0;
    start=clock();
    do{
        AP_keygen(&pfc, &pubparams, &apPriparam, &apPubparam);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("AP_keygen cost %8.2lf ms\n",elapsed);

    iterations = 0;
    start=clock();
    do{
        IdP_keygen(&pfc, &pubparams, &idpprikey, &idppubkey);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("IdP_keygen cost %8.2lf ms\n",elapsed);

    iterations = 0;
    start=clock();
    do{
        Register(&pfc, &apPriparam, &apPubparam, &wu, &Wu);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("Register cost %8.2lf ms\n",elapsed);


    pfc.random(a1);
    a1 = a1 % *pfc.ord;

    iterations = 0;
    request_time = 0;
    issue_time = 0;
    do{
        start=clock();
        RequestCred1_bmark(&pfc, &pubparams, &apPubparam, wu, &Wu, a1, &T, &tau, &request_time, &issue_time, &start);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=request_time/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("RequestCred1 cost %8.2lf ms\n",elapsed);

    iterations = 0;
    start=clock();
    do{
        IssueCred1(&pfc, &pubparams, &idpprikey, a1, &A, &B);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start + issue_time)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("IssueCred1 cost %8.2lf ms\n",elapsed);

    show_time = 0;
    verify_time = 0;
    do{
        start=clock();
        Show_bmark(&pfc, &pubparams, &apPubparam, &idppubkey, wu, &Wu, a1, &T, &tau, &A, &B, &show_time, &verify_time, &start);
        iterations++;
    }while(iterations < ITERATIONS); 
    elapsed=show_time/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed;
    printf("Show cost %8.2lf ms\n",elapsed);
    elapsed=verify_time/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed;
    printf("Verify cost %8.2lf ms\n",elapsed);

    pfc.random(a2);
    a2 = a2 % *pfc.ord;

    iterations = 0;
    request_time = 0;
    issue_time = 0;
    do{
        start=clock();
        RequestCredN_bmark(&pfc, &pubparams, &apPubparam, wu, &Wu, a1, a2, &T, &tau, &request_time, &issue_time, &start);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(request_time)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("RequestCredN cost %8.2lf ms\n",elapsed);

    iterations = 0;
    start=clock();
    do{
        IssueCred1(&pfc, &pubparams, &idpprikey, a2, &A2, &B2);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start + issue_time)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("IssueCredN cost %8.2lf ms\n",elapsed);

    show_time = 0;
    verify_time = 0;
    do{
        start=clock();
        ShowN_bmark(&pfc, &pubparams, &apPubparam, &idppubkey, wu, &Wu, a1, a2, &T, &tau, &A, &B, &A2, &B2, &show_time, &verify_time, &start);        
        iterations++;
    }while(iterations < ITERATIONS); 
    elapsed=show_time/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed;
    printf("ShowN cost %8.2lf ms\n",elapsed);
    elapsed=verify_time/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed;
    printf("VerifyN cost %8.2lf ms\n",elapsed);

    iterations = 0;
    start=clock();
    do{
        RevokeUser(&pfc, &apPubparam, &apPriparam, apPriparam.w1);
        iterations++;
    }while(iterations < ITERATIONS);
    elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    elapsed=1000.0*elapsed/ITERATIONS;
    printf("RevokeUser cost %8.2lf ms\n",elapsed);

}


int PS_test(){
    PFC pfc(SECURITY);

    G1 g1, A, B;
    G2 g2, X2, Y2;
    Big x, y, m;

    G1 g, h, At, Bt, C, D;
    GT E;
    Big c, r, s, t, rho1, rho2, rho3; 
    Big s1, s2, s3;


    int result;

    PS_setup(&pfc, &g1, &g2);

    PS_keygen(&pfc, &g1, &g2, &x, &y, &X2, &Y2);

    pfc.random(m);

    m = m % *pfc.ord;

    PS_sign(&pfc, &g1, &g2, &x, &y, m, &A, &B);

    result = PS_verify(&pfc, &g1, &g2, &X2, &Y2, m, &A, &B);

    cout << "PS signature verify result = " << result << endl;
    
    pfc.random(g);
    pfc.random(h);
    pfc.random(r);

    C = pfc.mult(g, m) + pfc.mult(h, r);

    rangeproof_step1(&pfc, &g1, &g2, &g, &h, &Y2, &A, &B, &At, &Bt, &D, &E, &s, &t, &rho1, &rho2, &rho3);

    pfc.random(c);
    c = c % *pfc.ord;

    rangeproof_step2(&pfc, &r, &m, &s, &t, &rho1, &rho2, &rho3, &c, &s1, &s2, &s3);

    result = rangeproof_step3(&pfc, &g2, &g, &h, &X2, &Y2, &At, &Bt, &c, &s1, &s2, &s3, &C, &D, &E);
    
    cout << "PS rangeproof verify result = " << result << endl;
    
    return 0;
}

int main() {
    //Algorithm_test();
    Algorithm_bmark();
    return 0;
}
